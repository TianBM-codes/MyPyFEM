#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from element.ElementBase import *
import numpy as np
from abc import ABC


class CPM6(ElementBaseClass, ABC):
    """
    6-node quadratic plane strain triangle, 在这里用作Cook膜单元

    Reference:
    1. A COMPATIBLE TRIANGULAR ELEMENT INCLUDING VERTEX ROTATIONS FOR PLANE ELASTICITY ANALYSIS   D.J.ALLMAN
    2. TECHNICAL NOTE ON THE ALLMAN TRIANGLE AND A RELATED QUADRILATERAL ELEMENT  Robert D.Cook
    3. 有限元理论、格式与求解方法 上 K.J Bathe P352
    """

    def __init__(self, eid=None):
        super().__init__(eid)
        self.nodes_count = 6  # Each element has 6 nodes
        self.K = np.zeros([12, 12], dtype=float)  # 刚度矩阵
        self.vtp_type = "triangle"
        self.thickness = None
        self.T_matrix = None  # 整体坐标转到局部坐标的矩阵, 是转换位移的

    def CalElementDMatrix(self, an_type=None):
        """
        计算本构矩阵, 弹性模量和泊松比, Bathe 上册P184
        """
        e = self.cha_dict[MaterialKey.E]
        niu = self.cha_dict[MaterialKey.Niu]
        if an_type == MaterialMatrixType.PlaneStree or an_type is None:
            a = e / (1 - niu ** 2)
            self.D = a * np.array([[1, niu, 0],
                                   [niu, 1, 0],
                                   [0, 0, 0.5 * (1 - niu)]], dtype=float)
        elif an_type == MaterialMatrixType.PlaneStrain:
            a = e * (1 - niu) / (1 + niu) / (1 - 2 * niu)
            self.D = a * np.array([[1, niu / (1 - niu), 0],
                                   [niu(1 - niu), 1, 0],
                                   [0, 0, 0.5 * (1 - 2 * niu) / (1 - niu)]], dtype=float)
        else:
            mlogger.fatal("Unknown an_dimension")
            sys.exit(1)

    def ElementStiffness(self):
        """
        p代表偏导: partial, ph1pr 代表偏h1偏r
        """
        assert self.node_coords.shape == (6, 2)  # 6节点, 2个坐标分量
        points, weights = GaussIntegrationPoint.GetTrianglePointAndWeight(3)
        for ii in range(len(points)):
            r, s = points[ii]
            w = weights[ii]

            ph1pr, ph2pr, ph3pr = -3 + 4 * (r + s), -1 + 4 * r, 0
            ph4pr, ph5pr, ph6pr = 4 * (1 - 2 * r - s), 4 * s, -4 * s
            ph1ps, ph2ps, ph3ps = ph1pr, 0, -1 + 4 * s
            ph4ps, ph5ps, ph6ps = -4 * r * (1 - r), 4 * s, r * (1 - r - 2 * s)

            phpr = np.array([ph1pr, ph2pr, ph3pr, ph4pr], dtype=float)
            phps = np.array([ph1ps, ph2ps, ph3ps, ph4ps], dtype=float)

            # Jacobi 3 * 3, J_ij代表当前积分点的雅可比矩阵, 描述了几何变形
            J_ij = np.asarray([[np.matmul(phpr, self.node_coords[:, 0]), np.matmul(phpr, self.node_coords[:, 1])],
                               [np.matmul(phps, self.node_coords[:, 0]), np.matmul(phps, self.node_coords[:, 1])]], dtype=float)

            det_J = np.linalg.det(J_ij)
            J_inv = np.asarray([[J_ij[1, 1], -J_ij[0, 1]],
                                [-J_ij[1, 0], J_ij[0, 0]]], dtype=float) / det_J

            # pupxy = [pupx, pupy].T,  pvpxy = [pvpx, pvpy].T
            pupxy = np.asarray([[ph1pr, 0, ph2pr, 0, ph3pr, 0, ph4pr, 0, ph5pr, 0, ph6pr, 0],
                                [ph1ps, 0, ph2ps, 0, ph3ps, 0, ph4ps, 0, ph5ps, 0, ph6ps, 0]], dtype=float)
            pvpxy = np.asarray([[0, ph1pr, 0, ph2pr, 0, ph3pr, 0, ph4pr, 0, ph5pr, 0, ph6pr],
                                [0, ph1ps, 0, ph2ps, 0, ph3ps, 0, ph4ps, 0, ph5ps, 0, ph6ps]], dtype=float)
            B1 = np.matmul(J_inv, pupxy)
            B2 = np.matmul(J_inv, pvpxy)

            # 组装B阵, [pupx, pvpy, pupy+pvpx], B Matrix 3*12
            B = np.insert(B1, 1, B2[1, :], axis=0)
            B[2, :] += B2[0, :]
            self.K += np.matmul(np.matmul(B.T, self.D), B) * w * det_J # TODO: ??? 这里不会约分掉det_J???

        # 以上是平面单元的刚度阵, 以下转换为膜单元刚度阵, 参考Reference2
        a1 = (self.node_coords[2, 0] - self.node_coords[1, 0]) * 0.125
        a2 = (self.node_coords[0, 0] - self.node_coords[2, 0]) * 0.125
        a3 = (self.node_coords[1, 0] - self.node_coords[0, 0]) * 0.125
        b1 = (self.node_coords[1, 1] - self.node_coords[2, 1]) * 0.125
        b2 = (self.node_coords[2, 1] - self.node_coords[0, 1]) * 0.125
        b3 = (self.node_coords[0, 1] - self.node_coords[1, 1]) * 0.125
        T = np.asarray([[1, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 1, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 1, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 1, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 1, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 1, 0],
                        [0.5, 0, b3, 0.5, 0, -b3, 0, 0, 0],
                        [0, 0.5, a3, 0, 0.5, -a3, 0, 0, 0],
                        [0, 0, 0, 0.5, 0, b1, 0.5, 0, -b1],
                        [0, 0, 0, 0, 0.5, a1, 0, 0.5, -a1],
                        [0.5, 0, -b2, 0, 0, 0, 0.5, 0, b2],
                        [0, 0.5, -a2, 0, 0, 0, 0, 0.5, a2]], dtype=float)

        # self.K *= self.cha_dict[PropertyKey.ThicknessOrArea]
        self.K *= self.cha_dict["RealConst"][0] # 只适用于等厚度的壳,
        return np.matmul(self.T_matrix, np.matmul(np.matmul(np.matmul(T.T, self.K), T), self.T_matrix.T))

    def ElementStress(self, displacement):
        """
        Calculate element stress
        """


class CPM8(ElementBaseClass, ABC):
    """
    8-node quadratic plane strain quadrangle, 在这里用作Cook膜单元

    Reference:
    1. A COMPATIBLE TRIANGULAR ELEMENT INCLUDING VERTEX ROTATIONS FOR PLANE ELASTICITY ANALYSIS   D.J.ALLMAN
    2. TECHNICAL NOTE ON THE ALLMAN TRIANGLE AND A RELATED QUADRILATERAL ELEMENT  Robert D.Cook
    3. A Refined Four-Noded Membrane Element With Rotational Degrees Of Freedom
    4. 有限元法 理论、格式与求解方法 上 K.J Bath P322
    """

    def __init__(self, eid=None):
        super().__init__(eid)
        self.nodes_count = 8  # Each element has 6 nodes
        self.K = np.zeros([16, 16], dtype=float)  # 刚度矩阵
        self.vtp_type = "triangle"
        self.thickness = None
        self.T_matrix = None  # 整体坐标转到局部坐标的矩阵, 是转换位移的

    def CalElementDMatrix(self, an_type=None):
        """
        计算本构矩阵, 弹性模量和泊松比, Bathe 上册P184
        """
        e = self.cha_dict[MaterialKey.E]
        niu = self.cha_dict[MaterialKey.Niu]
        if an_type == MaterialMatrixType.PlaneStree or an_type is None:
            a = e / (1 - niu ** 2)
            self.D = a * np.array([[1, niu, 0],
                                   [niu, 1, 0],
                                   [0, 0, 0.5 * (1 - niu)]], dtype=float)
        elif an_type == MaterialMatrixType.PlaneStrain:
            a = e * (1 - niu) / (1 + niu) / (1 - 2 * niu)
            self.D = a * np.array([[1, niu / (1 - niu), 0],
                                   [niu(1 - niu), 1, 0],
                                   [0, 0, 0.5 * (1 - 2 * niu) / (1 - niu)]], dtype=float)
        else:
            mlogger.fatal("Unknown an_dimension")
            sys.exit(1)

    def ElementStiffness(self):
        """
        p代表偏导: partial, phpr 代表偏hi偏r求和
        """
        assert self.node_coords.shape == (8, 2)  # 8节点, 2个坐标分量
        points, weights = GaussIntegrationPoint.GetTrianglePointAndWeight(3)
        for ii in range(len(points)):
            r, s = points[ii]
            w = weights[ii]
            ph1pr = 0.25 * (s ** 2 + s) + 0.5 * (1 + s) * r
            ph2pr = -ph1pr
            ph3pr = -0.25 * (1 - s) + 0.25 * (1 - s ** 2) + 0.5 * r * (1 - s)
            ph4pr = 0.25 * (s ** 2 - s) + 0.5 * r*(1 - s)
            ph5pr = -r * (1 + s)
            ph6pr = 0.5 * (s ** 2 - 1)
            ph7pr = r * (s - 1)
            ph8pr = 0.5 * (1 - s ** 2)

            ph1ps = 0.25 * (r ** 2 + r) + 0.5 * s * (1 + r)
            ph2ps = 0.25 * (r ** 2 - r) + 0.5 * s * (1 - r)
            ph3ps = 0.25 * (r - r ** 2) + 0.5 * s * (1 - r)
            ph4ps = -0.25 * (r ** 2 + r) + 0.5 * s * (1 + r)
            ph5ps = 0.5 * (1 - r ** 2)
            ph6ps = s * (1 - r)
            ph7ps = 0.5 * (r ** 2 - 1)
            ph8ps = -s * (1 + r)

            phpr = np.array([ph1pr, ph2pr, ph3pr, ph4pr, ph5pr, ph6pr, ph7pr, ph8pr], dtype=float)
            phps = np.array([ph1ps, ph2ps, ph3ps, ph4ps, ph5ps, ph6ps, ph7ps, ph8ps], dtype=float)

            # Jacobi 3 * 3, J_ij代表当前积分点的雅可比矩阵, 描述了几何变形
            J_ij = np.asarray([[np.matmul(phpr, self.node_coords[:, 0]), np.matmul(phpr, self.node_coords[:, 1])],
                               [np.matmul(phps, self.node_coords[:, 0]), np.matmul(phps, self.node_coords[:, 1])]], dtype=float)

            det_J = np.linalg.det(J_ij)
            J_inv = np.asarray([[J_ij[1, 1], -J_ij[0, 1]],
                                [-J_ij[1, 0], J_ij[0, 0]]], dtype=float) / det_J

            # pupxy = [pupx, pupy].T,  pvpxy = [pvpx, pvpy].T
            pupxy = np.asarray([[ph1pr, 0, ph2pr, 0, ph3pr, 0, ph4pr, 0, ph5pr, 0, ph6pr, 0, ph7pr, 0, ph8pr, 0],
                                [ph1ps, 0, ph2ps, 0, ph3ps, 0, ph4ps, 0, ph5ps, 0, ph6ps, 0, ph7ps, 0, ph8ps, 0]], dtype=float)
            pvpxy = np.asarray([[0, ph1pr, 0, ph2pr, 0, ph3pr, 0, ph4pr, 0, ph5pr, 0, ph6pr, 0, ph7pr, 0, ph8pr],
                                [0, ph1ps, 0, ph2ps, 0, ph3ps, 0, ph4ps, 0, ph5ps, 0, ph6ps, 0, ph7ps, 0, ph8ps]], dtype=float)
            B1 = np.matmul(J_inv, pupxy)
            B2 = np.matmul(J_inv, pvpxy)

            # 组装B阵, [pupx, pvpy, pupy+pvpx], B Matrix 3*12
            B = np.insert(B1, 1, B2[1, :], axis=0)
            B[2, :] += B2[0, :]
            self.K += np.matmul(np.matmul(B.T, self.D), B) * w * det_J # TODO: ??? 这里不会约分掉det_J???

        # 以上是平面单元的刚度阵, 以下转换为膜单元刚度阵, 参考Reference2
        e = 10e-8
        a12 = (self.node_coords[1, 0] - self.node_coords[0, 0]) * 0.125
        a23 = (self.node_coords[2, 0] - self.node_coords[1, 0]) * 0.125
        a34 = (self.node_coords[3, 0] - self.node_coords[2, 0]) * 0.125
        a41 = (self.node_coords[0, 0] - self.node_coords[3, 0]) * 0.125

        b12 = (self.node_coords[1, 1] - self.node_coords[0, 1]) * 0.125
        b23 = (self.node_coords[2, 1] - self.node_coords[1, 1]) * 0.125
        b34 = (self.node_coords[3, 1] - self.node_coords[2, 1]) * 0.125
        b41 = (self.node_coords[0, 1] - self.node_coords[3, 1]) * 0.125

        T = np.asarray([[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
                        [0.5, 0, b12, 0.5, 0, -b12, 0, 0, 0, 0, 0, 0],
                        [0, 0.5, a12, 0, 0.5, -a12, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0.5, 0, b23, 0.5, 0, -b23, 0, 0, 0],
                        [0, 0, 0, 0, 0.5, a23, 0, 0.5, -a23, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0.5, 0, b34, 0.5, 0, -b34],
                        [0, 0, 0, 0, 0, 0, 0, 0.5, a34, 0, 0.5, -a34],
                        [0.5, 0, -b41, 0, 0, 0, 0, 0, 0, 0.5, 0, b41],
                        [0, 0.5, -a41, 0, 0, 0, 0, 0, 0, 0, 0.5, a41]], dtype=float)

        local_k = np.matmul(np.matmul(T.T, self.K), T) * self.cha_dict["RealConst"][0]  # 只适用于等厚度的壳
        return np.matmul(np.matmul(self.T_matrix, local_k), self.T_matrix.T)

    def ElementStress(self, displacement):
        """
        Calculate element stress
        """


if __name__ == "__main__":
    ele = CPM6()
