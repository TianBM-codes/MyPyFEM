#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from abc import ABC
from typing import Tuple

from femdb.ShapeFunctionsAndInteg import *
from element.ElementBase import *


class CSTDrill(ElementBaseClass, ABC):
    """
    6-node quadratic plane strain triangle, 在这里用作Cook膜单元

    Reference:
    1. A COMPATIBLE TRIANGULAR ELEMENT INCLUDING VERTEX ROTATIONS FOR PLANE ELASTICITY ANALYSIS   D.J.AllMan
    2. TECHNICAL NOTE ON THE AllMan TRIANGLE AND A RELATED QUADRILATERAL ELEMENT  Robert D.Cook
    3. 有限元理论、格式与求解方法 上 K.J Bathe P352
    """

    def __init__(self, eid=None):
        super().__init__(eid)
        self.nodes_count = 6  # Each element has 6 nodes
        self.K = np.zeros([12, 12], dtype=float)  # 刚度矩阵
        self.vtu_type = "triangle"
        self.thickness = None
        self.T_matrix = None  # 整体坐标转到局部坐标的矩阵, 是转换位移的
        self.h = None

    def CalElementDMatrix(self, an_type=None):
        """
        计算本构矩阵, 弹性模量和泊松比, Bathe 上册P184
        """
        e = self.cha_dict[MaterialKey.E]
        if self.sec_id:
            self.h = self.cha_dict[self.sec_id]
        elif self.cha_dict.__contains__('RealConst'):
            self.h = self.cha_dict["RealConst"]
        else:
            raise KeyError("Don't Contain RealConst and sec_id")
        niu = self.cha_dict[MaterialKey.Niu]
        a = e * self.h / (1 - niu ** 2)
        self.D = a * np.array([[1, niu, 0],
                               [niu, 1, 0],
                               [0, 0, 0.5 * (1 - niu)]], dtype=float)

    def ElementStiffness(self):
        """
        p代表偏导: partial, ph1pr 代表偏h1偏r
        """
        assert self.node_coords.shape == (3, 2)  # 3节点, 2个坐标分量
        x2 = self.node_coords[1, 0]
        x3 = self.node_coords[2, 0]
        y3 = self.node_coords[2, 1]
        A2 = y3 * x2
        B = np.array([
            [-1 / x2, 0, 0, 1 / x2, 0, 0, 0, 0, 0],
            [0, (x3 - x2) / (x2 * y3), 0, 0, -x3 / (x2 * y3), 0, 0, 1 / y3, 0],
            [(x3 - x2) / (x2 * y3), -1 / x2, 0, -x3 / (x2 * y3), 1 / x2, 0, 1 / y3, 0, 0]
        ])

        # 计算基础刚度矩阵
        G = self.cha_dict[MaterialKey.G]
        Ke = B.T @ self.D @ B * (A2 / 2)  # 矩阵乘法运算符@

        # Drilling自由度修正部分
        Q = np.array([[
            -(x2 - x3) / (2 * A2), y3 / (2 * A2), 1 / 3,
            -x3 / (2 * A2), -y3 / (2 * A2), 1 / 3,
            x2 / (2 * A2), 0, 1 / 3
        ]])

        delta2 = 0.01
        V = (A2 / 2) * self.h
        Sr = (delta2 * V * G) * (Q.T @ Q)  # 刚度修正项

        # 合并刚度矩阵
        Ke += Sr

        # 对角项增强（注意Python索引从0开始）
        Ke[2, 2] += G * V * 1e-8
        Ke[5, 5] += G * V * 1e-8
        Ke[8, 8] += G * V * 1e-8

        return Ke

    def ElementStress(self, displacement):
        """
        Calculate element stress
        """


class Q4Mem(ElementBaseClass, ABC):
    """
    王欢Q4Mem
    """

    def __init__(self, eid=None):
        super().__init__(eid)
        self.nodes_count = 4
        self.K = np.zeros([12, 12], dtype=float)  # 刚度矩阵
        self.vtu_type = "triangle"
        self.thickness = None
        self.T_matrix = None  # 整体坐标转到局部坐标的矩阵, 是转换位移的

        self.integ = IntegForm2D2P()
        self.shpFunc = shpFunc2D4Node(self.integ)
        # 用于旋转自由度
        self.shpFunc8 = shpFunc2D8Node(self.integ)
        # 用于计算不协调应变插值矩阵[G]
        self.gauss_centor = IntegForm2D1P()
        self.shpFunc_centor = shpFunc2D4Node(self.gauss_centor)
        self.shpFunc_Incom = shpFunc4NodeIncom(self.integ)

        # 单元局部坐标
        self.node_coords = None
        self.h = None  # 厚度
        self.G = None  # 剪切模量

    def CalElementDMatrix(self, an_type=None):
        """
        计算本构矩阵
        @param an_type:
        @return:
        """
        e = self.cha_dict[MaterialKey.E]
        if self.sec_id:
            self.h = self.cha_dict[self.sec_id]
        elif self.cha_dict.__contains__('RealConst'):
            self.h = self.cha_dict["RealConst"]
        else:
            raise KeyError("Don't Contain RealConst and sec_id")
        niu = self.cha_dict[MaterialKey.Niu]
        a = e * self.h / (1 - niu ** 2)
        self.D = a * np.array([[1, niu, 0],
                               [niu, 1, 0],
                               [0, 0, 0.5 * (1 - niu)]], dtype=float)
        self.G = e / 2 / (1 + niu)

    def strain_matrix(self, igauss: int) -> Tuple[np.ndarray, float]:
        """计算应变矩阵B和雅可比行列式"""
        J = self.shpFunc[igauss].DNDxi @ self.node_coords[:, :2]
        DNDx = np.linalg.solve(J, self.shpFunc[igauss].DNDxi)

        B = np.zeros((3, 12))
        for iNode in range(4):
            B[0, (iNode * 3)] = DNDx[0, iNode]
            B[1, (iNode * 3) + 1] = DNDx[1, iNode]
            B[2, (iNode * 3)] = DNDx[1, iNode]
            B[2, (iNode * 3) + 1] = DNDx[0, iNode]

        j = np.linalg.det(J)
        return B, j

    def incom_strain_matrix(self, igauss: int, j: float, J0: np.ndarray) -> np.ndarray:
        """计算不协调应变矩阵G"""
        scale = np.linalg.det(J0) / j
        DMbarDx = np.linalg.solve(scale * J0, self.shpFunc_Incom[igauss].DNDxi)

        G = np.zeros((3, 4))
        G[0, 0] = DMbarDx[0, 0]
        G[1, 1] = DMbarDx[1, 0]
        G[2, 0] = DMbarDx[1, 0]
        G[2, 1] = DMbarDx[0, 0]
        G[0, 2] = DMbarDx[0, 1]
        G[1, 3] = DMbarDx[1, 1]
        G[2, 2] = DMbarDx[1, 1]
        G[2, 3] = DMbarDx[0, 1]

        return G

    def ElementStiffness(self):
        """计算刚度矩阵"""
        # 膜部分
        Ke = np.zeros((12, 12))
        Kmud = np.zeros((12, 4))
        Kmdd = np.zeros((4, 4))

        J0 = self.shpFunc_centor[0].DNDxi @ self.node_coords[:, :2]
        V = 0.0

        for igauss in range(4):
            B, j = self.strain_matrix(igauss)
            G = self.incom_strain_matrix(igauss, j, J0)

            wgt = 1.0
            Ke += B.T @ self.D @ B * j * wgt
            Kmud += B.T @ self.D @ G * j * wgt
            Kmdd += G.T @ self.D @ G * j * wgt
            V += j * wgt * self.h

            # 凝聚
        Ke -= Kmud @ np.linalg.inv(Kmdd) @ Kmud.T

        # 钻孔部分
        delta2 = 0.01
        DN0Dx = np.linalg.solve(J0, self.shpFunc_centor[0].DNDxi)
        Q = 0.5 * np.array([
            [-DN0Dx[1, 0], DN0Dx[0, 0], -0.5, -DN0Dx[1, 1], DN0Dx[0, 1], -0.5,
             -DN0Dx[1, 2], DN0Dx[0, 2], -0.5, -DN0Dx[1, 3], DN0Dx[0, 3], -0.5]
        ])
        Sr = (delta2 * V * self.G) * (Q.T @ Q)
        Ke += Sr

        # 为旋转自由度添加小刚度
        lambda_ = 1e-4
        Ke[2, 2] += self.G * V * lambda_
        Ke[5, 5] += self.G * V * lambda_
        Ke[8, 8] += self.G * V * lambda_
        Ke[11, 11] += self.G * V * lambda_

        return Ke

    def ElementStress(self, displacement: np.array):
        """"""
        pass


class CPM6(ElementBaseClass, ABC):
    """
    6-node quadratic plane strain triangle, 在这里用作Cook膜单元

    Reference:
    1. A COMPATIBLE TRIANGULAR ELEMENT INCLUDING VERTEX ROTATIONS FOR PLANE ELASTICITY ANALYSIS   D.J.AllMan
    2. TECHNICAL NOTE ON THE AllMan TRIANGLE AND A RELATED QUADRILATERAL ELEMENT  Robert D.Cook
    3. 有限元理论、格式与求解方法 上 K.J Bathe P352
    """

    def __init__(self, eid=None):
        super().__init__(eid)
        self.nodes_count = 6  # Each element has 6 nodes
        self.K = np.zeros([12, 12], dtype=float)  # 刚度矩阵
        self.vtu_type = "triangle"
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

            ph1pr, ph2pr, ph3pr = -3 + 4 * (r + s), 4 * r - 1, 0
            ph4pr, ph5pr, ph6pr = 4 * (1 - 2 * r - s), 4 * s, -4 * s
            ph1ps, ph2ps, ph3ps = ph1pr, 0, -1 + 4 * s
            ph4ps, ph5ps, ph6ps = -4 * r, 4 * r, 4 * (1 - r - 2 * s)

            phpr = np.array([ph1pr, ph2pr, ph3pr, ph4pr, ph5pr, ph6pr], dtype=float)
            phps = np.array([ph1ps, ph2ps, ph3ps, ph4ps, ph5ps, ph6ps], dtype=float)

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
            self.K += np.matmul(np.matmul(B.T, self.D), B) * w * det_J  # TODO: ??? 这里不会约分掉det_J???

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

        return T.T @ self.K @ T * self.cha_dict[self.sec_id]  # 只适用于等厚度的壳

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
        self.vtu_type = "triangle"
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

        # 八节点四边形单元用3阶高斯积分公式
        points, weights = GaussIntegrationPoint.GetSamplePointAndWeight(3)
        for ri in range(3):
            for si in range(3):
                r, s = points[ri], points[si]
                w_r, w_s = weights[ri], weights[si]
                ph1pr = 0.25 * (s ** 2 + s) + 0.5 * (1 + s) * r
                ph2pr = -0.25 * (s ** 2 + s) + 0.5 * (1 + s) * r
                ph3pr = 0.25 * (s - s ** 2) + 0.5 * r * (1 - s)
                ph4pr = 0.25 * (s ** 2 - s) + 0.5 * r * (1 - s)
                ph5pr = -r * (1 + s)
                ph6pr = 0.5 * (s ** 2 - 1)
                ph7pr = r * (s - 1)
                ph8pr = 0.5 * (1 - s ** 2)

                ph1ps = 0.25 * (r ** 2 + r) + 0.5 * s * (1 + r)
                ph2ps = 0.25 * (r ** 2 - r) + 0.5 * s * (1 - r)
                ph3ps = 0.25 * (r - r ** 2) + 0.5 * s * (1 - r)
                ph4ps = -0.25 * (r ** 2 + r) + 0.5 * s * (1 + r)
                ph5ps = 0.5 * (1 - r ** 2)
                ph6ps = s * (r - 1)
                ph7ps = 0.5 * (r ** 2 - 1)
                ph8ps = -s * (1 + r)

                phpr = np.array([ph1pr, ph2pr, ph3pr, ph4pr, ph5pr, ph6pr, ph7pr, ph8pr], dtype=float)
                phps = np.array([ph1ps, ph2ps, ph3ps, ph4ps, ph5ps, ph6ps, ph7ps, ph8ps], dtype=float)

                # Jacobi 3 * 3, J_ij代表当前积分点的雅可比矩阵(因为在积分过程中Jacobi矩阵是变化的), 描述了几何变形
                J_ij = np.asarray([[np.matmul(phpr[:], self.node_coords[:, 0]), np.matmul(phpr[:], self.node_coords[:, 1])],
                                   [np.matmul(phps[:], self.node_coords[:, 0]), np.matmul(phps[:], self.node_coords[:, 1])]], dtype=float)

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
                self.K += np.matmul(np.matmul(B.T, self.D), B) * w_r * w_s * det_J  # TODO: ??? 这里不会约分掉det_J???

        """
        以上是平面单元的刚度阵, 以下转换为膜单元刚度阵, 参考Reference2
        """
        # e = 10e-8
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

        return T.T @ self.K @ T * self.cha_dict[self.sec_id]  # 只适用于等厚度的壳

    def ElementStress(self, displacement):
        """
        Calculate element stress
        """


if __name__ == "__main__":
    t_ele = Q4Mem()
    t_ele.sec_id = 10001
    t_ele.cha_dict = {MaterialKey.Niu: 0.3, MaterialKey.E: 2e11, 10001: 0.01}
    t_ele.node_coords = np.array([
        [0, 0, 0],
        [1, 0, 0],
        [1, 1, 0],
        [0.5, 1, 0]
    ], dtype=float)
    t_ele.CalElementDMatrix()
    Ke1 = t_ele.ElementStiffness()
    print(Ke1)
