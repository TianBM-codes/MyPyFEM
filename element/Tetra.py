from element.ElementBase import *


class C3D4(ElementBaseClass):
    """ Tetra Element class """

    def __init__(self, eid=None):
        super().__init__(eid)
        self.nodes_count = 4  # Each element has 4 nodes
        self.K = np.zeros([12, 12], dtype=float)  # 刚度矩阵
        self.vtp_type = "tetra"
        self.unv_code = 40800

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
                               [0, 0, 0, 0, 0, (1 - 2 * niu) / 2.]])

    def ElementStiffness(self):
        """
        Reference:
        1. https://www.help.febio.org/FEBio/FEBio_tm_2_7/FEBio_tm_2-7-Subsection-4.1.3.html#prev

        # Shape Function:
        N1 = 1 - r - s - t
        N2 = r
        N3 = s
        N4 = t

        # Partial
        pN1pr, pN1ps, pN1pt = -1, -1, -1
        pN2pr, pN2ps, pN2pt = 1, 0, 0
        pN3pr, pN3ps, pN3pt = 0, 1, 0
        pN4pr, pN4ps, pN4pt = 0, 0, 1
        """
        assert self.node_coords.shape == (4, 3)

        # Gaussian Weight
        weight = 0.166666667
        r, s, t = 0.25, 0.25, 0.25
        dNdr = np.array([[-1, -1, -1],
                         [1, 0, 0],
                         [0, 1, 0],
                         [0, 0, 1]], dtype=float).T

        # Jacobi 3*3 & B Matrix 6*12
        J = np.matmul(dNdr, self.node_coords)
        det_J = np.linalg.det(J)
        J_inv = np.linalg.inv(J)
        B_pre = np.matmul(J_inv, dNdr)
        B = np.array([[B_pre[0, 0], 0, 0, B_pre[0, 1], 0, 0, B_pre[0, 2], 0, 0, B_pre[0, 3], 0, 0],
                      [0, B_pre[1, 0], 0, 0, B_pre[1, 1], 0, 0, B_pre[1, 2], 0, 0, B_pre[1, 3], 0],
                      [0, 0, B_pre[2, 0], 0, 0, B_pre[2, 1], 0, 0, B_pre[2, 2], 0, 0, B_pre[2, 3]],
                      [B_pre[1, 0], B_pre[0, 0], 0, B_pre[1, 1], B_pre[0, 1], 0, B_pre[1, 2], B_pre[0, 2], 0, B_pre[1, 3], B_pre[0, 3], 0],
                      [0, B_pre[2, 0], B_pre[1, 0], 0, B_pre[2, 1], B_pre[1, 1], 0, B_pre[2, 2], B_pre[1, 2], 0, B_pre[2, 3], B_pre[1, 3]],
                      [B_pre[2, 0], 0, B_pre[0, 0], B_pre[2, 1], 0, B_pre[0, 1], B_pre[2, 2], 0, B_pre[0, 2], B_pre[2, 3], 0, B_pre[0, 3]]], dtype=float)

        self.K = self.K + np.matmul(np.matmul(B.T, self.D), B) * det_J * weight  # 常数还需要乘以weight？

        return self.K

    def ElementStress(self, displacement):
        """
        Calculate element stress
        """
