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
        self.vtp_type = "wedge"
        self.unv_code = 60600
        self.gs_count = 6  # 高斯积分点个数

    def CalElementDMatrix(self, an_type=None):
        """
        计算本构矩阵, 弹性模量和泊松比, Bathe 上册P184
        """
        e = self.cha_dict[MaterialKey.E]
        niu = self.cha_dict[MaterialKey.Niu]
        a = e / ((1 + niu) * (1 - 2 * niu))
        self.D = a * np.array([[1 - niu, niu, niu, 0, 0, 0],
                               [niu, 1 - niu, niu, 0, 0, 0],
                               [niu, niu, 1 - niu, 0, 0, 0],
                               [0, 0, 0, (1 - 2 * niu) / 2., 0, 0],
                               [0, 0, 0, 0, (1 - 2 * niu) / 2., 0],
                               [0, 0, 0, 0, 0, (1 - 2 * niu) / 2.]], dtype=float)

    def ElementStiffness(self):
        """
        Reference:
        1. https://www.help.febio.org/FEBio/FEBio_tm_2_7/FEBio_tm_2-7-Subsection-4.1.2.html#:~:text=Pentahedral%20elements%20%28also%20knows%20as%20%E2%80%9Cwedge%E2%80%9D%20elements%29%20consist,s%20and%20t%20and%20are%20given%20as%20follows.
        2. https://github.com/febiosoftware
        """
        assert self.node_coords.shape == (6, 3)

        dNdrs, weights = AllEleTypeDNDr.C3D6
        # 在6个高斯点上积分
        for ii in range(self.gs_count):
            # Jacobi 3*3 & B Matrix 8*24
            J = np.matmul(dNdrs[ii], self.node_coords)
            det_J = np.linalg.det(J)
            J_inv = np.linalg.inv(J)
            B_pre = np.matmul(J_inv, dNdrs[ii])
            B = np.array([[B_pre[0, 0], 0, 0, B_pre[0, 1], 0, 0, B_pre[0, 2], 0, 0, B_pre[0, 3], 0, 0, B_pre[0, 4], 0, 0, B_pre[0, 5], 0, 0],
                          [0, B_pre[1, 0], 0, 0, B_pre[1, 1], 0, 0, B_pre[1, 2], 0, 0, B_pre[1, 3], 0, 0, B_pre[1, 4], 0, 0, B_pre[1, 5], 0],
                          [0, 0, B_pre[2, 0], 0, 0, B_pre[2, 1], 0, 0, B_pre[2, 2], 0, 0, B_pre[2, 3], 0, 0, B_pre[2, 4], 0, 0, B_pre[2, 5]],
                          [B_pre[1, 0], B_pre[0, 0], 0, B_pre[1, 1], B_pre[0, 1], 0, B_pre[1, 2], B_pre[0, 2], 0, B_pre[1, 3], B_pre[0, 3], 0, B_pre[1, 4], B_pre[0, 4], 0, B_pre[1, 5], B_pre[0, 5], 0],
                          [0, B_pre[2, 0], B_pre[1, 0], 0, B_pre[2, 1], B_pre[1, 1], 0, B_pre[2, 2], B_pre[1, 2], 0, B_pre[2, 3], B_pre[1, 3], 0, B_pre[2, 4], B_pre[1, 4], 0, B_pre[2, 5], B_pre[1, 5]],
                          [B_pre[2, 0], 0, B_pre[0, 0], B_pre[2, 1], 0, B_pre[0, 1], B_pre[2, 2], 0, B_pre[0, 2], B_pre[2, 3], 0, B_pre[0, 3], B_pre[2, 4], 0, B_pre[0, 4], B_pre[2, 5], 0, B_pre[0, 5]]], dtype=float)

            self.K = self.K + np.matmul(np.matmul(B.T, self.D), B) * det_J * weights[ii]

        return self.K

    def ElementStress(self, displacement):
        """
        Calculate element stress
        """
