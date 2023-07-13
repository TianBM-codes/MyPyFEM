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

        # Shape Function:
        N1 = 0.5 * (1 - r - s) * (1 - t)
        N2 = 0.5 * r * (1 - t)
        N3 = 0.5 * s * (1 - t)
        N4 = 0.5 * (1 - r - s) * (1 + t)
        N5 = 0.5 * r * (1 + t)
        N6 = 0.5 * s * (1 + t)

        # Partial
        pN1pr, pN1ps, pN1pt = 0.5 * (t - 1), 0.5 * (t - 1), 0.5 * (r + s - 1)
        pN2pr, pN2ps, pN2pt = 0.5 * (1 - t), 0, -0.5 * r
        pN3pr, pN3ps, pN3pt = 0, 0.5 * (1 - t), -0.5 * s
        pN4pr, pN4ps, pN4pt = -0.5 * (1 + t), -0.5 * (1 + t), 0.5 * (1 - r - s)
        pN5pr, pN5ps, pN5pt = 0.5 * (1 + t), 0, 0.5 * r
        pN6pr, pN6ps, pN6pt = 0, 0.5 * (1 + t), 0.5 * s
        """
        assert self.node_coords.shape == (6, 3)

        # Gaussian Weight
        sample_r = [0.166666667, 0.666666667, 0.166666667, 0.166666667, 0.666666667, 0.166666667]
        sample_s = [0.166666667, 0.166666667, 0.666666667, 0.166666667, 0.166666667, 0.666666667]
        sample_t = [-0.577350269, -0.577350269, -0.577350269, 0.577350269, 0.577350269, 0.577350269]
        weight = 0.166666667

        # 在6个高斯点上积分
        for ii in range(6):
            r, s, t = sample_r[ii], sample_s[ii], sample_t[ii]
            dNdr = np.array([[0.5 * (t - 1), 0.5 * (t - 1), 0.5 * (r + s - 1)],
                             [0.5 * (1 - t), 0, -0.5 * r],
                             [0, 0.5 * (1 - t), -0.5 * s],
                             [-0.5 * (1 + t), -0.5 * (1 + t), 0.5 * (1 - r - s)],
                             [0.5 * (1 + t), 0, 0.5 * r],
                             [0, 0.5 * (1 + t), 0.5 * s]], dtype=float).T

            # Jacobi 3*3 & B Matrix 8*24
            J = np.matmul(dNdr, self.node_coords)
            det_J = np.linalg.det(J)
            J_inv = np.linalg.inv(J)
            B_pre = np.matmul(J_inv, dNdr)
            B = np.array([[B_pre[0, 0], 0, 0, B_pre[0, 1], 0, 0, B_pre[0, 2], 0, 0, B_pre[0, 3], 0, 0, B_pre[0, 4], 0, 0, B_pre[0, 5], 0, 0],
                          [0, B_pre[1, 0], 0, 0, B_pre[1, 1], 0, 0, B_pre[1, 2], 0, 0, B_pre[1, 3], 0, 0, B_pre[1, 4], 0, 0, B_pre[1, 5], 0],
                          [0, 0, B_pre[2, 0], 0, 0, B_pre[2, 1], 0, 0, B_pre[2, 2], 0, 0, B_pre[2, 3], 0, 0, B_pre[2, 4], 0, 0, B_pre[2, 5]],
                          [B_pre[1, 0], B_pre[0, 0], 0, B_pre[1, 1], B_pre[0, 1], 0, B_pre[1, 2], B_pre[0, 2], 0, B_pre[1, 3], B_pre[0, 3], 0, B_pre[1, 4], B_pre[0, 4], 0, B_pre[1, 5], B_pre[0, 5], 0],
                          [0, B_pre[2, 0], B_pre[1, 0], 0, B_pre[2, 1], B_pre[1, 1], 0, B_pre[2, 2], B_pre[1, 2], 0, B_pre[2, 3], B_pre[1, 3], 0, B_pre[2, 4], B_pre[1, 4], 0, B_pre[2, 5], B_pre[1, 5]],
                          [B_pre[2, 0], 0, B_pre[0, 0], B_pre[2, 1], 0, B_pre[0, 1], B_pre[2, 2], 0, B_pre[0, 2], B_pre[2, 3], 0, B_pre[0, 3], B_pre[2, 4], 0, B_pre[0, 4], B_pre[2, 5], 0, B_pre[0, 5]]], dtype=float)

            self.K = self.K + np.matmul(np.matmul(B.T, self.D), B) * det_J * weight

        return self.K

    def ElementStress(self, displacement):
        """
        Calculate element stress
        """
