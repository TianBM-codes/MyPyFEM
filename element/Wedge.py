from abc import ABC

from element.ElementBase import *


class C3D6(ElementBaseClass, ABC):
    """
    Wedge Element class, also known as "Pentahedral"
    TODO: 调研https://github.com/febiosoftware
    """

    def __init__(self, eid=None):
        super().__init__(eid)
        self.nodes_count = 6  # Each element has 6 nodes
        self.K = np.zeros([18, 18], dtype=float)  # 刚度矩阵
        self.vtu_type = "wedge"
        self.unv_code = 60600
        self.gs_count = 6  # 高斯积分点个数
        self.block_size = 324
        self.node_dof_count = 3

    def CalElementDMatrix(self, an_type=None):
        """
        计算本构矩阵, 弹性模量和泊松比, Bathe 上册P184
        """
        e = self.cha_dict[MaterialKey.E]
        niu = self.cha_dict[MaterialKey.Niu]
        a = e / ((1 + niu) * (1 - 2 * niu))
        b = (1 - niu) * a
        c = niu * a
        d = 0.5 * (1 - 2 * niu) * a

        self.D = np.zeros((6, 6), dtype=float)
        self.D[0, 0] = b
        self.D[1, 1] = b
        self.D[2, 2] = b
        self.D[0, 1] = c
        self.D[0, 2] = c
        self.D[1, 0] = c
        self.D[1, 2] = c
        self.D[2, 0] = c
        self.D[2, 1] = c
        self.D[3, 3] = d
        self.D[4, 4] = d
        self.D[5, 5] = d

    def ElementStiffness(self, from_origin=False):
        """
        Reference:
        1. https://www.help.febio.org/FEBio/FEBio_tm_2_7/FEBio_tm_2-7-Subsection-4.1.2.html#:~:text=Pentahedral%20elements%20%28also%20knows%20as%20%E2%80%9Cwedge%E2%80%9D%20elements%29%20consist,s%20and%20t%20and%20are%20given%20as%20follows.
        2. https://github.com/febiosoftware
        """
        assert self.node_coords.shape == (6, 3)
        self.CalElementDMatrix()

        dNdrs, weights = AllEleTypeDNDrAtGaussianPoint.C3D6
        # 重置刚度矩阵
        self.K.fill(0.0)

        # 在6个高斯点上积分 - 使用Numba加速
        contributions = np.zeros((self.gs_count, 18, 18), dtype=float)
        # 在6个高斯点上积分 - 使用Numba加速
        for ii in numba.prange(self.gs_count):  # 并行循环
            # 计算单个高斯点的贡献
            contributions[ii] = self._compute_gauss_point_contribution(
                dNdrs[ii], weights[ii], node_coords=self.node_coords, D=self.D
            )

        # 串行累加所有贡献
        for ii in range(self.gs_count):
            self.K += contributions[ii]

        return self.K

    @staticmethod
    @numba.jit(nopython=True, cache=True, fastmath=True)
    def _compute_gauss_point_contribution(dNdr, weight, node_coords, D):
        """
        计算单个高斯点对刚度矩阵的贡献
        使用JIT加速
        """
        # Jacobi 3 * 3
        J = np.dot(dNdr, node_coords)

        # 计算雅可比行列式
        det_J = np.linalg.det(J)

        # 避免求逆运算 - 使用直接解法
        J_inv = np.linalg.inv(J)  # 对于3x3矩阵，求逆通常足够高效

        # B_pre = J_inv @ dNdr
        B_pre = np.zeros((3, 6))
        for i in range(3):
            for j in range(6):
                for k in range(3):
                    B_pre[i, j] += J_inv[i, k] * dNdr[k, j]

        # 构建B矩阵 (6x18)
        B = np.zeros((6, 18))
        for i in range(6):  # 遍历6个节点
            # 节点i对应的自由度位置
            idx = i * 3

            # 填充B矩阵
            # εxx
            B[0, idx] = B_pre[0, i]
            # εyy
            B[1, idx + 1] = B_pre[1, i]
            # εzz
            B[2, idx + 2] = B_pre[2, i]
            # γxy
            B[3, idx] = B_pre[1, i]
            B[3, idx + 1] = B_pre[0, i]
            # γyz
            B[4, idx + 1] = B_pre[2, i]
            B[4, idx + 2] = B_pre[1, i]
            # γxz
            B[5, idx] = B_pre[2, i]
            B[5, idx + 2] = B_pre[0, i]

        BTDB = B.T @ (D @ B)
        return BTDB * det_J * weight

    def CalculateElementStress(self, displacement):
        """
        Calculate element stress
        """
        return np.zeros(6), np.zeros(6), np.zeros(6), np.zeros(6), np.zeros(6), np.zeros(6)

    def ElementMass(self):
        pass

    def CalculateBasic(self):
        pass

    def ReCalculateElementStiffness(self):
        pass

    def CalculateThermalStress(self, U, T, T0):
        pass

    def ElementThermalLoadVector(self, T, T0):
        pass

    def ElementThermalMatrix(self):
        pass
