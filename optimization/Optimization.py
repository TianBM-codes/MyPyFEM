from FEMDataBase import FEMDataBase
from femdb.Domain import Domain
from femdb.GlobalEnum import *
import pypardiso
import numpy as np

# import matplotlib.pyplot as plt

DOF_MAP = {
    "ux": 0,
    "uy": 1,
    "uz": 2,
    "rx": 3,
    "ry": 4,
    "rz": 5,
}


def plot_history(hist, save_path=None, show=True, use_logy=True):
    """
    hist: list of tuples (loss, loss_misfit, loss_reg, grad_norm)
    save_path: e.g. "opt_history.png" or "opt_history.pdf"
    show: whether plt.show()
    use_logy: use log-scale for y-axis (recommended)
    """
    hist = np.asarray(hist, dtype=np.float64)
    if hist.ndim != 2 or hist.shape[1] != 4:
        raise ValueError("hist must be an array-like of shape (iters, 4)")

    loss = hist[:, 0]
    misfit = hist[:, 1]
    reg = hist[:, 2]
    grad_norm = hist[:, 3]

    iters = np.arange(1, len(loss) + 1)

    plt.figure(figsize=(10, 6))
    plt.plot(iters, loss, label="Total Loss")
    # plt.plot(iters, misfit, label="Misfit Loss")
    # plt.plot(iters, reg, label="Regularization Loss")
    # plt.plot(iters, grad_norm, label="Grad Norm")

    plt.xlabel("Iteration")
    plt.ylabel("Value")
    plt.title("Optimization History")
    plt.grid(True, which="both", linestyle="--", alpha=0.3)
    plt.legend()

    if use_logy:
        # 避免 log(0) 崩掉
        plt.yscale("log")

    plt.tight_layout()

    if save_path is not None:
        plt.savefig(save_path, dpi=200)

    if show:
        plt.show()

    plt.close()


class Optimization(object):
    def __init__(self):
        self.domain = Domain(False)
        self.femdb = None

    def build_element_dof_indices(self, ele):
        """
        返回该单元对应的全局结构 DOF 索引数组（长度 24）
        按 ele.search_node_ids 的 4 个角点顺序，每节点 6 DOF。
        """
        dofs = np.empty(4 * ele.node_dof_count, dtype=np.int64)  # 24
        ptr = 0
        for nidx in ele.search_node_ids:  # nidx 是 node_list 下标（你原代码就是这么用的）
            start = self.femdb.node_list[nidx].start_eq_num
            for k in range(ele.node_dof_count):
                dofs[ptr] = start + k
                ptr += 1
        return dofs

    # ----------------------------
    # 2) 预计算：每个单元 f_base（alpha=1）以及单元 DOF 映射
    # ----------------------------
    def precompute_element_f_base(self, target_element_filter=None):
        """
        预计算每个参与优化的单元：
        - f_base_e: alpha=1 时的热等效载荷 (24,)
        - dofs_e: 该单元对应的全局 DOF 索引 (24,)

        target_element_filter:
            None -> 默认所有 femdb.elements 都参与
            or callable(ele)->bool 只对某些单元参与优化（例如只壳单元）
        """
        # 确保温度场已存在
        if not hasattr(self.femdb, "temperature_res") or self.femdb.temperature_res is None:
            raise RuntimeError("femdb.temperature_res not found. Call CalculateSteadyTemperature() first.")
        if not hasattr(self.femdb, "thermal_reference"):
            raise RuntimeError("femdb.thermal_reference not found.")

        T_all = self.femdb.temperature_res
        T0 = self.femdb.thermal_reference

        f_base_list = []
        dof_map_list = []
        opt_ele_list = []

        for ele in self.femdb.elements:
            if target_element_filter is not None and (not target_element_filter(ele)):
                continue

            # 取 4 角点温度（用 temp_eq_num 的体系）
            temp_idx = [self.femdb.node_list[ii].temp_eq_num for ii in ele.search_node_ids]
            T4 = T_all[temp_idx]

            # 暂存原 alpha
            alpha_old = ele.cha_dict.get(MaterialKey.Expansion, None)

            # 设置 alpha=1，计算基热载荷
            ele.cha_dict[MaterialKey.Expansion] = 1.0
            f_base = ele.ElementThermalLoadVector(T4, T0).astype(np.float64)  # (24,)

            # 还原
            if alpha_old is None:
                del ele.cha_dict[MaterialKey.Expansion]
            else:
                ele.cha_dict[MaterialKey.Expansion] = alpha_old

            dofs = self.build_element_dof_indices(ele)

            opt_ele_list.append(ele)
            f_base_list.append(f_base)
            dof_map_list.append(dofs)

        return opt_ele_list, f_base_list, dof_map_list

    def build_obs_dofs_from_node_ids(self, node_id_list, dof_names):
        """
        输入：节点 id 列表（注意：是 node_id，不是 node_list 下标）
             dof_names: ["ux","uy","uz"] 等
        输出：I_obs (m,)
        """
        I_obs = []
        for node_id in node_id_list:
            node_index = self.femdb.node_hash[node_id]  # 你数据库里 node_hash: node_id -> index
            node = self.femdb.node_list[node_index]
            start = node.start_eq_num
            for dn in dof_names:
                I_obs.append(start + DOF_MAP[dn])
        return np.array(I_obs, dtype=np.int64)

    # ----------------------------
    # 3) 组装全局热载荷：F = sum(alpha_e * f_base_e)
    # ----------------------------
    def assemble_global_thermal_load(self, alpha_vec, f_base_list, dof_map_list, n_dof_u):
        F = np.zeros(n_dof_u, dtype=np.float64)
        for a, fbase, dofs in zip(alpha_vec, f_base_list, dof_map_list):
            F[dofs] += a * fbase
        return F

    def loss_and_grad_alpha(
            self,
            alpha_vec,
            alpha_ref,
            gamma,
            I_obs,
            d_obs,
            f_base_list,
            dof_map_list,
    ):
        """
        输入：
          alpha_vec: (n_ele,) 当前设计变量
          alpha_ref: (n_ele,) 参考值
          gamma: 正则系数
          I_obs: (m,) 观测 DOF 全局索引（结构位移向量 u 的索引）
          d_obs: (m,) 目标位移
        输出：
          loss, grad, u
        """
        K = self.femdb.global_stiff_matrix
        n_dof_u = K.shape[0]

        # 1) 组装载荷、求解位移
        F = self.assemble_global_thermal_load(alpha_vec, f_base_list, dof_map_list, n_dof_u)
        u = pypardiso.spsolve(K, F)  # (n_dof_u,)

        # 2) 位移匹配项
        r = u[I_obs] - d_obs
        loss_misfit = 0.5 * float(np.dot(r, r))

        # 3) 正则项（不要偏离太远）
        diff = alpha_vec - alpha_ref
        loss_reg = 0.5 * float(gamma * np.dot(diff, diff))

        loss = loss_misfit + loss_reg

        # 4) 构造 dL/du
        dL_du = np.zeros_like(u, dtype=np.float64)
        dL_du[I_obs] = r  # S^T(Su - d)

        # 5) 伴随解：K^T lam = dL/du （K 对称则同 K）
        lam = pypardiso.spsolve(K, dL_du)

        # 6) 梯度：grad_e = lam_e^T f_base_e + gamma*(alpha-alpha_ref)
        grad = np.zeros_like(alpha_vec, dtype=np.float64)
        for e, (fbase, dofs) in enumerate(zip(f_base_list, dof_map_list)):
            grad[e] = float(np.dot(lam[dofs], fbase))

        grad += gamma * (alpha_vec - alpha_ref)

        return loss, grad, u, loss_misfit, loss_reg

    def optimize_alpha_adam(
            self,
            alpha0,
            alpha_ref,
            gamma,
            I_obs,
            d_obs,
            f_base_list,
            dof_map_list,
            alpha_min=None,
            alpha_max=None,
            iters=200,
            lr=1e-2,
            beta1=0.9,
            beta2=0.999,
            eps=1e-8,
            verbose_every=10,

            # ---- early stopping controls ----
            tol_loss_rel=1e-6,  # loss relative change tolerance
            tol_alpha_rel=1e-6,  # alpha relative change tolerance
            tol_grad=1e-8,  # gradient norm tolerance
            patience=10,  # number of consecutive steps to trigger stop
    ):
        """
        返回：alpha_opt, history

        Early stopping:
          - |loss_{k}-loss_{k-1}|/max(1,|loss_{k-1}|) < tol_loss_rel
          - ||alpha_k - alpha_{k-1}|| / max(1, ||alpha_{k-1}||) < tol_alpha_rel
          - ||grad|| < tol_grad
          以上条件满足的连续次数达到 patience 时停止
        """
        alpha = alpha0.astype(np.float64).copy()
        m = np.zeros_like(alpha)
        v = np.zeros_like(alpha)

        history = []

        # early stop state
        stable_count = 0
        prev_loss = None
        prev_alpha = None

        for t in range(1, iters + 1):
            loss, grad, u, loss_misfit, loss_reg = self.loss_and_grad_alpha(
                alpha, alpha_ref, gamma, I_obs, d_obs, f_base_list, dof_map_list
            )

            grad_norm = float(np.linalg.norm(grad))

            # ---- Adam update ----
            m = beta1 * m + (1 - beta1) * grad
            v = beta2 * v + (1 - beta2) * (grad * grad)
            m_hat = m / (1 - beta1 ** t)
            v_hat = v / (1 - beta2 ** t)

            alpha_new = alpha - lr * m_hat / (np.sqrt(v_hat) + eps)

            # 物理范围裁剪
            if alpha_min is not None:
                alpha_new = np.maximum(alpha_new, alpha_min)
            if alpha_max is not None:
                alpha_new = np.minimum(alpha_new, alpha_max)

            # ---- convergence checks ----
            # loss relative change
            if prev_loss is None:
                loss_rel_change = np.inf
            else:
                loss_rel_change = abs(loss - prev_loss) / max(1.0, abs(prev_loss))

            # alpha relative change
            if prev_alpha is None:
                alpha_rel_change = np.inf
            else:
                da = alpha_new - prev_alpha
                alpha_rel_change = float(np.linalg.norm(da) / max(1.0, np.linalg.norm(prev_alpha)))

            # update stable counter
            cond_loss = (loss_rel_change < tol_loss_rel)
            cond_alpha = (alpha_rel_change < tol_alpha_rel)
            cond_grad = (grad_norm < tol_grad)

            # 你可以选择：
            # 1) 三个都满足才算稳定
            # stable = cond_loss and cond_alpha and cond_grad
            #
            # 2) 更“松”的：loss+alpha 满足即可（一般够用）
            stable = (cond_loss and cond_alpha) or cond_grad

            if stable:
                stable_count += 1
            else:
                stable_count = 0

            # 提交更新
            alpha = alpha_new
            prev_loss = loss
            prev_alpha = alpha.copy()

            history.append((loss, loss_misfit, loss_reg, grad_norm))

            # print
            if verbose_every and (t % verbose_every == 0 or t == 1):
                print(
                    f"iter {t:04d}  loss={loss:.6e}  misfit={loss_misfit:.6e}  reg={loss_reg:.6e}  "
                    f"|grad|={grad_norm:.3e}  alpha_mean={alpha.mean():.6e}  "
                    f"dLoss_rel={loss_rel_change:.3e}  dAlpha_rel={alpha_rel_change:.3e}  stable={stable_count}/{patience}"
                )

            # early stopping
            if stable_count >= patience:
                print(
                    f"[Early Stop] Converged at iter {t}: "
                    f"dLoss_rel={loss_rel_change:.3e}, dAlpha_rel={alpha_rel_change:.3e}, |grad|={grad_norm:.3e}"
                )
                break

        return alpha, history

    def OptExpansionByDis(self, obs_node_ids,
                          obs_dof_names,
                          d_obs,  # (m,) 目标位移
                          alpha_ref_value,  # 标量或 (n_ele,) 数组
                          gamma=1e6,
                          iters=200,
                          lr=1e-2,
                          alpha_min=0.0,
                          alpha_max=None,
                          target_element_filter=None):
        """
        优化热膨胀系数
        - obs_node_ids: [101, 203, ...] 观测节点
        - obs_dof_names: ["ux","uy","uz"] 观测哪些 DOF
        - d_obs: 与 obs_node_ids x obs_dof_names 对应的目标位移
        - alpha_ref_value: 参考热膨胀（比如 1.2e-5）
        - gamma: 正则强度（越大越不让 alpha 偏离参考）
        @return:
        """
        self.femdb = self.domain.femdb

        # 1) 构造观测 DOF
        I_obs = self.build_obs_dofs_from_node_ids(obs_node_ids, obs_dof_names)
        d_obs = np.asarray(d_obs, dtype=np.float64).reshape(-1)
        assert len(I_obs) == len(d_obs), "I_obs and d_obs length mismatch"

        # 2) 预计算单元基载荷（只对参与优化的单元）
        opt_ele_list, f_base_list, dof_map_list = self.precompute_element_f_base(target_element_filter)

        n_ele = len(opt_ele_list)
        print(f"Optimization elements: {n_ele}")

        # 3) 初始化 alpha / alpha_ref
        if np.isscalar(alpha_ref_value):
            alpha_ref = np.full(n_ele, float(alpha_ref_value), dtype=np.float64)
        else:
            alpha_ref = np.asarray(alpha_ref_value, dtype=np.float64).copy()
            assert alpha_ref.shape[0] == n_ele

        alpha0 = alpha_ref.copy()  # 从参考值出发最稳（也可以小扰动）

        # 4) 优化
        alpha_opt, hist = self.optimize_alpha_adam(
            alpha0=alpha0,
            alpha_ref=alpha_ref,
            gamma=gamma,
            I_obs=I_obs,
            d_obs=d_obs,
            f_base_list=f_base_list,
            dof_map_list=dof_map_list,
            alpha_min=alpha_min,
            alpha_max=alpha_max,
            iters=iters,
            lr=lr,
            verbose_every=10,
        )
        plot_history(hist)

        return alpha_opt, alpha_ref, opt_ele_list, hist
