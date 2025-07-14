#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
from copy import deepcopy
from pathlib import Path

import numpy as np
from femdb.FEMDataBase import *  # 有限元数据库
from femdb.GlobalFEMVariant import ModelInfo  # 模型全局信息
from scipy.sparse.linalg import factorized, spsolve  # 稀疏矩阵求解器
from femdb.GlobalEnum import *  # 全局枚举
import pypardiso  # 高性能稀疏矩阵求解库

"""
**关于稀疏矩阵的注意事项:**
   1. 要有效地构造矩阵, 请使用dok_matrix或lil_matrix。lil_matrix类支持基本切片和花式索引,
      其语法与NumPy Array类似。lil_matrix形式是基于行的，因此能够很高效地转为csr格式。
   2. 不建议直接使用NumPy函数运算稀疏矩阵。如果想将NumPy函数应用于这些矩阵，首先要检查SciPy
      是否有自己的给定稀疏矩阵类的实现，或者先将稀疏矩阵转换为NumPy数组(使用toarray()方法)。
   3. 要执行乘法或转置等操作，首先将矩阵转换为CSC或CSR格式，效率更高。CSR格式特别适用于快速矩阵向量乘法。
   4. CSR，CSC和COO格式之间的转换都是线性复杂度。
   5. 对于已知是正定对称矩阵的情况下，如何用scipy快速求解逆矩阵:
      https://stackoverflow.com/questions/40703042/...
"""


# TODO: 尝试其他稀疏矩阵求解库
# 列出了多个备选的稀疏矩阵求解库和优化方案

class Domain(object):
    """
    有限元问题域类，定义问题域并执行有限元分析
    采用单例模式，确保只有一个实例
    """

    def __init__(self, check_model):
        """
        初始化有限元问题域
        :param check_model: 是否进行模型检查的标志
        """
        self.femdb = FEMDataBase()  # 有限元数据库实例
        self.eq_count = None  # 总刚度矩阵维度
        self.free_dof_count = None  # 自由度的个数
        self.bound_dof_count = None  # 被约束的自由度个数
        self.Ub = []  # 约束指定位移
        self.Ua = None  # 未被约束的自由度
        self.Ra = None  # 未被约束自由度上的力或力矩
        self.right_hand = None  # 右端项(载荷向量)
        # 刚度矩阵相关
        self.mass_list = []  # 质量矩阵列表
        self.check_model = check_model  # 模型检查标志

    def CalAllElementStiffness(self):
        """
        计算所有单元的刚度矩阵(多线程)
        对所有的单元组进行循环计算单元刚度矩阵
        可以方便地查看各步骤运行时间
        """
        calculated_eles_count = 0
        for iter_ele in self.femdb.elements:
            iter_ele.CalculateBasic()  # 计算单元基本属性
            stiff = iter_ele.ElementStiffness()  # 计算单元刚度矩阵
            calculated_eles_count += 1

            # 模型检查相关
            if self.check_model:
                if iter_ele.id == 786:  # 调试特定单元
                    print("")
                # 检查刚度矩阵是否有全零行
                has_zero_row = (stiff >= 1e-8).all(axis=1).any()
                if has_zero_row:
                    print(f"Element {iter_ele.id} has zero row")

            self.femdb.stiff_list.append(stiff)  # 存储刚度矩阵

    def CheckRedundantNodes(self):
        """
        检查模型中是否存在不属于任何单元的冗余节点
        """
        no_dup_nodes = np.zeros(len(self.femdb.node_list), dtype=np.uint32)
        for ii, iter_node in enumerate(self.femdb.node_list):
            no_dup_nodes[ii] = iter_node.id
        nodes_set = set(no_dup_nodes)

        element_nodes = []
        for iter_ele in self.femdb.elements:
            element_nodes.extend(iter_ele.node_ids.tolist())
        ele_nodes_set = set(element_nodes)

        redundant_set = nodes_set - ele_nodes_set
        print("Redundant Nodes:", list(redundant_set))

    def AssembleStiffnessMatrixByPenalty(self, from_origin=False):
        """
        使用罚函数法组装总体刚度矩阵
        与消元法不同，罚函数法通过添加大数来施加约束，而不是消除自由度

        :param from_origin: 是否从原始数据重新计算
        """
        # 计算约束方程增加的方程数量
        ce_count = len(self.femdb.equation_constrain_couple)
        ce_add_equations = 0
        for ii in range(ce_count):
            iter_ce = self.femdb.equation_constrain_idx[ii]
            ce_add_equations += len(iter_ce)

        # 准备COO格式的三元组数据(行、列、值)
        iter_loc = 0
        rows = np.zeros(self.femdb.matrix_num_count + ce_add_equations * 4, dtype=np.uint32)
        cols = np.zeros(self.femdb.matrix_num_count + ce_add_equations * 4, dtype=np.uint32)
        datas = np.zeros(self.femdb.matrix_num_count + ce_add_equations * 4, dtype=np.float64)

        # 组装单元刚度矩阵
        for kk, ele in enumerate(self.femdb.elements):
            ele_nodes = ele.search_node_ids
            if from_origin:
                stiff_array = ele.ReCalculateElementStiffness()
            else:
                stiff_array = self.femdb.stiff_list[kk]

            n_nodes = len(ele_nodes)
            block_size = n_nodes * ModelInfo.PER_NODE_DOF
            el_dofs = np.array([x * ModelInfo.PER_NODE_DOF + np.arange(ModelInfo.PER_NODE_DOF)
                                for x in ele_nodes]).flatten()

            # 生成单元自由度对应的全局行列索引
            rows_block, cols_block = np.meshgrid(el_dofs, el_dofs)
            entries_count = block_size ** 2

            # 填充三元组数据
            rows[iter_loc:iter_loc + entries_count] = rows_block.ravel()
            cols[iter_loc:iter_loc + entries_count] = cols_block.ravel()
            datas[iter_loc:iter_loc + entries_count] = stiff_array.ravel()
            iter_loc += entries_count

        # 添加RBE2约束(多点约束)
        GlobalNodeHash = self.femdb.node_hash
        matrix_dimension = len(self.femdb.node_list) * ModelInfo.PER_NODE_DOF
        row_cursor = 0
        insert_offset = 0

        for ii in range(ce_count):
            m_node, s_node = self.femdb.equation_constrain_couple[ii]
            iter_ce_idx = self.femdb.equation_constrain_idx[ii]
            m_node_dof = GlobalNodeHash[m_node] * ModelInfo.PER_NODE_DOF
            s_node_dof = GlobalNodeHash[s_node] * ModelInfo.PER_NODE_DOF

            for jj in iter_ce_idx:
                insert_begin = self.femdb.matrix_num_count + insert_offset
                insert_end = insert_begin + 4

                # 添加约束方程对应的矩阵项
                rows[insert_begin: insert_end] = [matrix_dimension + row_cursor,
                                                  matrix_dimension + row_cursor,
                                                  s_node_dof + jj,
                                                  m_node_dof + jj]
                cols[insert_begin:insert_end] = [s_node_dof + jj,
                                                 m_node_dof + jj,
                                                 matrix_dimension + row_cursor,
                                                 matrix_dimension + row_cursor]
                datas[insert_begin: insert_end] = [1, -1, 1, -1]
                insert_offset += 4
                row_cursor += 1

        matrix_dimension += row_cursor

        # 组装总体刚度矩阵(COO格式转为CSC格式)
        self.femdb.global_stiff_matrix = sparse.coo_matrix(
            (datas, (rows, cols)),
            shape=(matrix_dimension, matrix_dimension)
        ).tocsc()

        # 模型检查
        if self.check_model:
            self.CheckRedundantNodes()  # 检查冗余节点
            print("GlobalStiffMatrix contains NaN:", np.isnan(self.femdb.global_stiff_matrix.data).any())

            # 计算矩阵条件数
            cond = np.linalg.cond(self.femdb.global_stiff_matrix.todense())
            print("Condition number:", cond)

            # 检查矩阵是否对称正定
            from scipy.linalg import eigh
            eigenvalues = eigh(self.femdb.global_stiff_matrix.todense())[0]
            print("Min eigenvalue:", np.min(eigenvalues))

        # 绘制刚度矩阵稀疏模式图
        if GlobalInfor[GlobalVariant.PlotGlobalStiffness]:
            plt.spy(self.femdb.global_stiff_matrix, markersize=1)
            plt.title("GlobalStiffnessMatrix")
            plt.savefig("GlobalStiffness.jpg")

    def CalAllElementMassMatrix(self):
        """
        计算所有单元的质量矩阵
        """
        for iter_ele in self.femdb.elements:
            self.mass_list.append(iter_ele.ElementMass())

    def AssembleMassMatrixByPerturbation(self):
        """
        组装总体质量矩阵并添加摄动项避免奇异
        """
        rows = []
        cols = []
        datas = []
        GlobalNodeHash = self.femdb.node_hash

        # 组装单元质量矩阵
        for kk, ele in enumerate(self.femdb.elements):
            nodes = ele.node_ids
            mass_array = self.mass_list[kk]

            for ii, nodeA in enumerate(nodes):
                equA = GlobalNodeHash[nodeA] * ModelInfo.PER_NODE_DOF
                for jj, nodeB in enumerate(nodes):
                    equB = GlobalNodeHash[nodeB] * ModelInfo.PER_NODE_DOF
                    for m in range(ModelInfo.PER_NODE_DOF):
                        for n in range(ModelInfo.PER_NODE_DOF):
                            TolRow = equA + m
                            TolCol = equB + n
                            eRow = ii * ModelInfo.PER_NODE_DOF + m
                            eCol = jj * ModelInfo.PER_NODE_DOF + n
                            if mass_array[eRow, eCol] != 0:
                                rows.append(TolRow)
                                cols.append(TolCol)
                                datas.append(mass_array[eRow, eCol])

        matrix_size = len(self.femdb.node_list) * ModelInfo.PER_NODE_DOF

        # 添加摄动项避免奇异
        for ii in range(matrix_size):
            rows.append(ii)
            cols.append(ii)
            datas.append(1e-12)

        # 组装总体质量矩阵
        self.femdb.global_mass_matrix = sparse.coo_matrix(
            (datas, (rows, cols)),
            shape=(matrix_size, matrix_size)
        ).tocsc()

    def CalculateGreenFunction(self, output_dir):
        """
        计算格林函数并保存
        :param output_dir: 输出目录路径
        """
        c_load = self.femdb.load_case.c_loads
        force_count = len(c_load)
        force_nodes = [ii[0] for ii in c_load]
        directories = [ii[1] for ii in c_load]
        right_hand_back_up = deepcopy(self.right_hand)

        solver = pypardiso.PyPardisoSolver()
        solver.factorize(self.femdb.global_stiff_matrix)

        for ii in range(force_count):
            # 施加单位载荷
            nd = force_nodes[ii]
            eqa = self.femdb.node_hash[nd] * ModelInfo.PER_NODE_DOF
            xyz = directories[ii]
            self.right_hand = deepcopy(right_hand_back_up)
            self.right_hand[eqa + xyz] = 1

            # 求解位移响应
            # temp_results = pypardiso.spsolve(self.femdb.global_stiff_matrix, self.right_hand)
            temp_results = solver.solve(self.femdb.global_stiff_matrix, self.right_hand)
            result_range = len(self.femdb.node_list) * ModelInfo.PER_NODE_DOF
            iter_linear_u = temp_results[:result_range]

            # 计算应力响应
            iter_sigma_xx = np.zeros(len(self.femdb.node_list))
            iter_sigma_yy = np.zeros(len(self.femdb.node_list))
            iter_sigma_zz = np.zeros(len(self.femdb.node_list))
            iter_tau_xy = np.zeros(len(self.femdb.node_list))
            iter_tau_yz = np.zeros(len(self.femdb.node_list))
            iter_tau_xz = np.zeros(len(self.femdb.node_list))

            for ele in self.femdb.elements:
                search_idx = []
                for jj in ele.search_node_ids:
                    start = jj * ModelInfo.PER_NODE_DOF
                    end = (jj + 1) * ModelInfo.PER_NODE_DOF
                    search_idx.extend(np.arange(start, end, 1).tolist())

                iter_u = iter_linear_u[search_idx].flatten()
                sigma_x, sigma_y, sigma_z, tau_yz, tau_xz, tau_xy = ele.CalculateElementStress(iter_u)

                # 节点应力平均
                for jj, n_search_id in enumerate(ele.search_node_ids):
                    iter_sigma_xx[n_search_id] += sigma_x[jj] / self.femdb.node_connected_element_count[n_search_id]
                    iter_sigma_yy[n_search_id] += sigma_y[jj] / self.femdb.node_connected_element_count[n_search_id]
                    iter_sigma_zz[n_search_id] += sigma_z[jj] / self.femdb.node_connected_element_count[n_search_id]
                    iter_tau_xy[n_search_id] += tau_xy[jj] / self.femdb.node_connected_element_count[n_search_id]
                    iter_tau_yz[n_search_id] += tau_yz[jj] / self.femdb.node_connected_element_count[n_search_id]
                    iter_tau_xz[n_search_id] += tau_xz[jj] / self.femdb.node_connected_element_count[n_search_id]

            # 保存应力结果
            iter_linear_u = np.reshape(iter_linear_u, (-1, ModelInfo.PER_NODE_DOF))
            np.save(output_dir / f"dis_node_{nd}_dir_{xyz}_x.npy", iter_linear_u[:, 0])
            np.save(output_dir / f"dis_node_{nd}_dir_{xyz}_y.npy", iter_linear_u[:, 1])
            np.save(output_dir / f"dis_node_{nd}_dir_{xyz}_z.npy", iter_linear_u[:, 2])
            np.save(output_dir / f"sigma_xx_node_{nd}_dir_{xyz}.npy", iter_sigma_xx)
            np.save(output_dir / f"sigma_yy_node_{nd}_dir_{xyz}.npy", iter_sigma_yy)
            np.save(output_dir / f"sigma_zz_node_{nd}_dir_{xyz}.npy", iter_sigma_zz)
            np.save(output_dir / f"sigma_xy_node_{nd}_dir_{xyz}.npy", iter_tau_xy)
            np.save(output_dir / f"sigma_yz_node_{nd}_dir_{xyz}.npy", iter_tau_yz)
            np.save(output_dir / f"sigma_xz_node_{nd}_dir_{xyz}.npy", iter_tau_xz)

    def CalculateResultsWithGreenFunction(self, gree_dir):
        """
        通过格林函数计算结果
        :return:
        """
        c_load = self.femdb.load_case.c_loads
        force_nodes = [ii[0] for ii in c_load]
        directories = [ii[1] for ii in c_load]
        amps = [ii[2] for ii in c_load]
        sigma_xx = np.zeros(len(self.femdb.node_list))
        sigma_yy = np.zeros(len(self.femdb.node_list))
        sigma_zz = np.zeros(len(self.femdb.node_list))
        tau_xy = np.zeros(len(self.femdb.node_list))
        tau_yz = np.zeros(len(self.femdb.node_list))
        tau_xz = np.zeros(len(self.femdb.node_list))
        displacement_x = np.zeros(len(self.femdb.node_list))
        displacement_y = np.zeros(len(self.femdb.node_list))
        displacement_z = np.zeros(len(self.femdb.node_list))
        for ii, nd in enumerate(force_nodes):
            xyz = directories[ii]
            sigma_xx_file_path = Path(gree_dir / f"sigma_xx_node_{nd}_dir_{xyz}.npy")
            sigma_yy_file_path = Path(gree_dir / f"sigma_yy_node_{nd}_dir_{xyz}.npy")
            sigma_zz_file_path = Path(gree_dir / f"sigma_zz_node_{nd}_dir_{xyz}.npy")
            sigma_xy_file_path = Path(gree_dir / f"sigma_xy_node_{nd}_dir_{xyz}.npy")
            sigma_yz_file_path = Path(gree_dir / f"sigma_yz_node_{nd}_dir_{xyz}.npy")
            sigma_xz_file_path = Path(gree_dir / f"sigma_xz_node_{nd}_dir_{xyz}.npy")
            dis_file_path_x = Path(gree_dir / f"dis_node_{nd}_dir_{xyz}_x.npy")
            dis_file_path_y = Path(gree_dir / f"dis_node_{nd}_dir_{xyz}_y.npy")
            dis_file_path_z = Path(gree_dir / f"dis_node_{nd}_dir_{xyz}_z.npy")
            if not sigma_xx_file_path.exists():
                raise FileNotFoundError(f"File Don't Exists: {sigma_xx_file_path}")
            if not sigma_yy_file_path.exists():
                raise FileNotFoundError(f"File Don't Exists: {sigma_yy_file_path}")
            if not sigma_zz_file_path.exists():
                raise FileNotFoundError(f"File Don't Exists: {sigma_zz_file_path}")
            if not sigma_xy_file_path.exists():
                raise FileNotFoundError(f"File Don't Exists: {sigma_xy_file_path}")
            if not sigma_yz_file_path.exists():
                raise FileNotFoundError(f"File Don't Exists: {sigma_yz_file_path}")
            if not sigma_xz_file_path.exists():
                raise FileNotFoundError(f"File Don't Exists: {sigma_xz_file_path}")
            if not dis_file_path_x.exists():
                raise FileNotFoundError(f"File Don't Exists: {dis_file_path_x}")
            if not dis_file_path_y.exists():
                raise FileNotFoundError(f"File Don't Exists: {dis_file_path_y}")
            if not dis_file_path_z.exists():
                raise FileNotFoundError(f"File Don't Exists: {dis_file_path_z}")
            iter_sigma_xx = np.load(sigma_xx_file_path)
            iter_sigma_yy = np.load(sigma_yy_file_path)
            iter_sigma_zz = np.load(sigma_zz_file_path)
            iter_tau_xy = np.load(sigma_xy_file_path)
            iter_tau_yz = np.load(sigma_yz_file_path)
            iter_tau_xz = np.load(sigma_xz_file_path)
            iter_dis_x = np.load(dis_file_path_x)
            iter_dis_y = np.load(dis_file_path_y)
            iter_dis_z = np.load(dis_file_path_z)
            sigma_xx += amps[ii] * iter_sigma_xx
            sigma_yy += amps[ii] * iter_sigma_yy
            sigma_zz += amps[ii] * iter_sigma_zz
            tau_xy += amps[ii] * iter_tau_xy
            tau_yz += amps[ii] * iter_tau_yz
            tau_xz += amps[ii] * iter_tau_xz
            displacement_x += amps[ii] * iter_dis_x
            displacement_y += amps[ii] * iter_dis_y
            displacement_z += amps[ii] * iter_dis_z

        term1 = (sigma_xx - sigma_yy) ** 2
        term2 = (sigma_yy - sigma_zz) ** 2
        term3 = (sigma_zz - sigma_xx) ** 2
        term4 = 6 * (tau_xy ** 2 + tau_yz ** 2 + tau_xz ** 2)
        self.femdb.linear_mises = np.sqrt(0.5 * (term1 + term2 + term3 + term4))
        self.femdb.linear_u = np.sqrt(displacement_x ** 2 + displacement_y ** 2 + displacement_z ** 2)

    def SolveDisplacement(self):
        """
        求解节点位移
        使用罚函数法处理边界条件，求解线性方程组
        """
        c_load = self.femdb.load_case.c_loads
        force_nodes = [ii[0] for ii in c_load]
        directories = [ii[1] for ii in c_load]
        amps = [ii[2] for ii in c_load]

        # 组装载荷向量
        for ii, nd in enumerate(force_nodes):
            eqa = self.femdb.node_hash[nd] * ModelInfo.PER_NODE_DOF
            xyz = directories[ii]
            self.right_hand[eqa + xyz] = amps[ii]

        try:
            # 使用高性能求解器求解
            temp_results = pypardiso.spsolve(self.femdb.global_stiff_matrix, self.right_hand)
            result_range = len(self.femdb.node_list) * ModelInfo.PER_NODE_DOF
            self.femdb.linear_u = temp_results[:result_range]
        except ValueError as e:
            # 处理奇异矩阵情况
            matrix_csr = self.femdb.global_stiff_matrix.tocsr()
            zero_rows = []
            for i in range(matrix_csr.shape[0]):
                row = matrix_csr.getrow(i)
                if row.nnz == 0:
                    zero_rows.append(i)

            self.CheckRedundantNodes()  # 检查冗余节点
            sys.exit(1)

    def SolveStress(self):
        """
        求解节点应力
        基于位移结果计算单元应力，然后平均到节点
        """
        self.femdb.linear_mises = np.zeros(len(self.femdb.node_list))
        self.femdb.sigma_xx = np.zeros(len(self.femdb.node_list))
        self.femdb.sigma_yy = np.zeros(len(self.femdb.node_list))
        self.femdb.sigma_zz = np.zeros(len(self.femdb.node_list))
        self.femdb.tau_xy = np.zeros(len(self.femdb.node_list))
        self.femdb.tau_xz = np.zeros(len(self.femdb.node_list))
        self.femdb.tau_yz = np.zeros(len(self.femdb.node_list))

        for ele in self.femdb.elements:
            search_idx = []
            for ii in ele.search_node_ids:
                start = ii * ModelInfo.PER_NODE_DOF
                end = (ii + 1) * ModelInfo.PER_NODE_DOF
                search_idx.extend(np.arange(start, end, 1).tolist())

            u = self.femdb.linear_u[search_idx].flatten()
            sigma_x, sigma_y, sigma_z, tau_yz, tau_xz, tau_xy = ele.CalculateElementStress(u)

            # 计算Mises应力
            term1 = (sigma_x - sigma_y) ** 2
            term2 = (sigma_y - sigma_z) ** 2
            term3 = (sigma_z - sigma_x) ** 2
            term4 = 6 * (tau_xy ** 2 + tau_yz ** 2 + tau_xz ** 2)
            von_mises = np.sqrt(0.5 * (term1 + term2 + term3 + term4))

            # 节点应力平均
            for ii, n_search_id in enumerate(ele.search_node_ids):
                self.femdb.linear_mises[n_search_id] += von_mises[ii] / self.femdb.node_connected_element_count[n_search_id]
                self.femdb.sigma_xx[n_search_id] += sigma_x[ii] / self.femdb.node_connected_element_count[n_search_id]
                self.femdb.sigma_yy[n_search_id] += sigma_y[ii] / self.femdb.node_connected_element_count[n_search_id]
                self.femdb.sigma_zz[n_search_id] += sigma_z[ii] / self.femdb.node_connected_element_count[n_search_id]
                self.femdb.tau_xy[n_search_id] += tau_xy[ii] / self.femdb.node_connected_element_count[n_search_id]
                self.femdb.tau_yz[n_search_id] += tau_yz[ii] / self.femdb.node_connected_element_count[n_search_id]
                self.femdb.tau_xz[n_search_id] += tau_xz[ii] / self.femdb.node_connected_element_count[n_search_id]

    def AddBoundaryByPenalty(self):
        """
        使用罚函数法施加边界条件
        通过添加大数到对角线元素来实现约束
        """
        node_count = len(self.femdb.node_list)
        ce_count = len(self.femdb.equation_constrain_couple)
        ce_add_equations = 0
        for ii in range(ce_count):
            iter_ce = self.femdb.equation_constrain_idx[ii]
            ce_add_equations += len(iter_ce)

        self.right_hand = np.zeros(node_count * ModelInfo.PER_NODE_DOF + ce_add_equations)

        # 处理边界条件
        bds = self.femdb.load_case.GetBoundaries()
        for bd in bds:
            for row in bd:
                nd = int(row[0])  # 节点ID
                dir_ = int(row[1])  # 方向
                val = float(row[2])  # 约束值
                nd_idx = self.femdb.node_hash[nd]
                base_idx = nd_idx * ModelInfo.PER_NODE_DOF

                # 添加大数施加约束
                self.femdb.global_stiff_matrix[base_idx + dir_, base_idx + dir_] += 2e21
                self.right_hand[base_idx + dir_] = val

    def NewMarkExplict(self):
        """
        显式Newmark方法求解动力学问题
        参考:《结构动力学基础》张亚辉、林家浩 P94
        """
        his_load = self.femdb.load_case.history_loads
        force_nodes = [ii[0] for ii in his_load]
        directories = [ii[1] for ii in his_load]
        scale = [ii[2] for ii in his_load]

        # 初始化参数
        steps_count = his_load[0][-1].shape[1]
        delta_t = his_load[0][-1][0, 1] - his_load[0][-1][0, 0]
        node_count = len(self.femdb.node_list)
        fem_all_dofs = ModelInfo.PER_NODE_DOF * node_count
        his_vals = np.zeros((fem_all_dofs, steps_count), dtype=float)

        # 组装时变载荷
        for ii, nd in enumerate(force_nodes):
            eqa = self.femdb.node_hash[nd] * ModelInfo.PER_NODE_DOF
            xyz = directories[ii]
            his_vals[eqa + xyz] = his_load[ii][-1][1, :] * scale[ii]

        # 初始加速度
        acc_0 = pypardiso.spsolve(self.femdb.global_mass_matrix, his_vals[:, 0])

        # Newmark参数
        alpha = 0.25
        delta = 0.5
        a0 = 1 / alpha / delta_t ** 2
        a1 = delta / alpha / delta_t
        a2 = 1 / alpha / delta_t
        a3 = 0.5 / alpha - 1
        a4 = delta / alpha - 1
        a5 = delta_t / 2 * (delta / alpha - 2)
        a6 = delta_t * (1 - delta)
        a7 = delta * delta_t

        # 预分解矩阵
        K_hat = self.femdb.global_stiff_matrix + a0 * self.femdb.global_mass_matrix
        solver = factorized(K_hat.tocsc())

        # 初始化历史变量
        u = [np.zeros(fem_all_dofs, dtype=float)]
        v = [np.zeros(fem_all_dofs, dtype=float)]
        a = [acc_0]

        # 时间步进求解
        for ii in range(steps_count - 1):
            f_hat = his_vals[:, ii + 1] + self.femdb.global_mass_matrix @ (a0 * u[ii] + a2 * v[ii] + a3 * a[ii])
            u_next = solver(f_hat)  # 求解位移
            a_next = a0 * (u_next - u[ii]) - a2 * v[ii] - a3 * a[ii]  # 计算加速度
            v_next = v[ii] + a6 * a[ii] + a7 * a_next  # 计算速度

            # 存储结果
            u.append(u_next)
            v.append(v_next)
            a.append(a_next)

        # 保存结果
        self.femdb.history_u = u
        self.femdb.history_v = v
        self.femdb.history_a = a
        self.femdb.history_step_count = len(u)


if __name__ == "__main__":
    import time

    n = 1300000
    a = np.random.rand(n)
    b = np.random.rand(n)

    start_time = time.time()
    for ii in range(600):
        c = 0.2 * a + 0.3 * b

    elapsed_time = time.time() - start_time
    print(f"Elapsed time: {elapsed_time:.6f}")

