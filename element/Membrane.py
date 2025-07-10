#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from abc import ABC
from femdb.ShapeFunctionsAndInteg import *
from element.ElementBase import *


class CPM3(ElementBaseClass, ABC):
    """
    3节点二次平面应变三角形单元(用于Cook膜单元)

    参考文献:
    1. A COMPATIBLE TRIANGULAR ELEMENT INCLUDING VERTEX ROTATIONS FOR PLANE ELASTICITY ANALYSIS - D.J.AllMan
    2. TECHNICAL NOTE ON THE AllMan TRIANGLE AND A RELATED QUADRILATERAL ELEMENT - Robert D.Cook
    3. 《有限元理论、格式与求解方法》上 - K.J Bathe P352

    特性:
    - 包含顶点旋转自由度的兼容三角形单元
    - 适用于平面弹性分析
    - 采用3点高斯积分
    """

    def __init__(self, eid=None):
        """
        初始化3节点三角形单元

        Args:
            eid (int, optional): 单元ID. Defaults to None.
        """
        super().__init__(eid)
        self.nodes_count = 3  # 单元节点数
        self.vtu_type = "triangle"  # VTK可视化类型
        self.T_matrix = None  # 全局坐标到局部坐标的转换矩阵(用于位移转换)
        self.B = []  # 应变-位移矩阵列表
        self.integ = []  # 积分权重列表

    def CalElementDMatrix(self, an_type=None):
        """
        计算平面应变问题的本构矩阵

        参考:
        《有限元方法》Bathe 上册P184

        Args:
            an_type (optional): 分析类型. Defaults to None.
        """
        e = self.cha_dict[MaterialKey.E]  # 弹性模量
        niu = self.cha_dict[MaterialKey.Niu]  # 泊松比

        # 获取单元厚度
        if self.cha_dict.__contains__('RealConst') and len(self.cha_dict["RealConst"]) != 0:
            h = self.cha_dict["RealConst"][0]
        elif self.cha_dict.__contains__(MaterialKey.Thickness):
            h = self.cha_dict[MaterialKey.Thickness]
        else:
            raise KeyError("未找到厚度参数(RealConst或Thickness)")

        # 计算平面应变问题的弹性矩阵
        a = e / (1 - niu ** 2) * h
        self.D = a * np.array([
            [1, niu, 0],
            [niu, 1, 0],
            [0, 0, 0.5 * (1 - niu)]
        ], dtype=float)

    def ElementStiffness(self, from_origin=False):
        """
        计算单元刚度矩阵

        采用3点高斯积分计算单元刚度矩阵
        包含将平面单元刚度矩阵转换为膜单元刚度矩阵的转换过程

        Args:
            from_origin (bool, optional): 是否从原始数据重新计算. Defaults to False.

        Returns:
            np.ndarray: 9x9的单元刚度矩阵
        """
        # 获取三角形积分点和权重(3点积分)
        points, weights = GaussIntegrationPoint.GetTrianglePointAndWeight(3)
        B_local = []

        # 在每个高斯点上计算刚度矩阵贡献
        for ii in range(len(points)):
            r, s = points[ii]
            w = weights[ii]

            # 计算形函数对自然坐标的导数
            ph1pr, ph2pr, ph3pr = -3 + 4 * (r + s), 4 * r - 1, 0
            ph4pr, ph5pr, ph6pr = 4 * (1 - 2 * r - s), 4 * s, -4 * s
            ph1ps, ph2ps, ph3ps = ph1pr, 0, -1 + 4 * s
            ph4ps, ph5ps, ph6ps = -4 * r, 4 * r, 4 * (1 - r - 2 * s)

            phpr = np.array([ph1pr, ph2pr, ph3pr, ph4pr, ph5pr, ph6pr], dtype=float)
            phps = np.array([ph1ps, ph2ps, ph3ps, ph4ps, ph5ps, ph6ps], dtype=float)

            # 计算雅可比矩阵(描述几何变形)
            J_ij = np.array([
                [np.matmul(phpr, self.node_coords[:, 0]), np.matmul(phpr, self.node_coords[:, 1])],
                [np.matmul(phps, self.node_coords[:, 0]), np.matmul(phps, self.node_coords[:, 1])]
            ], dtype=float)

            det_J = np.linalg.det(J_ij)  # 雅可比行列式
            J_inv = np.array([
                [J_ij[1, 1], -J_ij[0, 1]],
                [-J_ij[1, 0], J_ij[0, 0]]
            ], dtype=float) / det_J  # 雅可比逆矩阵

            # 计算位移导数矩阵
            pupxy = np.array([
                [ph1pr, 0, ph2pr, 0, ph3pr, 0, ph4pr, 0, ph5pr, 0, ph6pr, 0],
                [ph1ps, 0, ph2ps, 0, ph3ps, 0, ph4ps, 0, ph5ps, 0, ph6ps, 0]
            ], dtype=float)

            pvpxy = np.array([
                [0, ph1pr, 0, ph2pr, 0, ph3pr, 0, ph4pr, 0, ph5pr, 0, ph6pr],
                [0, ph1ps, 0, ph2ps, 0, ph3ps, 0, ph4ps, 0, ph5ps, 0, ph6ps]
            ], dtype=float)

            # 计算应变-位移矩阵B
            B1 = np.matmul(J_inv, pupxy)
            B2 = np.matmul(J_inv, pvpxy)

            # 组装完整的B矩阵(3x12)
            B = np.insert(B1, 1, B2[1, :], axis=0)
            B[2, :] += B2[0, :]  # 添加剪切应变项

            # 计算当前积分点对刚度矩阵的贡献
            if self.K is None:
                self.K = B.T @ self.D @ B * w * det_J
            else:
                self.K += B.T @ self.D @ B * w * det_J

            B_local.append(B)
            self.integ.append(w * det_J)

        # 将平面单元刚度矩阵转换为膜单元刚度矩阵(参考Cook的论文)
        a1 = (self.node_coords[2, 0] - self.node_coords[1, 0]) * 0.125
        a2 = (self.node_coords[0, 0] - self.node_coords[2, 0]) * 0.125
        a3 = (self.node_coords[1, 0] - self.node_coords[0, 0]) * 0.125
        b1 = (self.node_coords[1, 1] - self.node_coords[2, 1]) * 0.125
        b2 = (self.node_coords[2, 1] - self.node_coords[0, 1]) * 0.125
        b3 = (self.node_coords[0, 1] - self.node_coords[1, 1]) * 0.125

        # 构造转换矩阵T
        T = np.array([
            [1, 0, 0, 0, 0, 0, 0, 0, 0],
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
            [0, 0.5, -a2, 0, 0, 0, 0, 0.5, a2]
        ], dtype=float)

        self.B_global = B_local @ T  # 全局坐标系下的B矩阵
        return T.T @ self.K @ T  # 转换后的刚度矩阵

    def CalculateElementStress(self, displacement):
        """
        计算单元应力

        Args:
            displacement (np.ndarray): 节点位移向量

        Returns:
            np.ndarray: 节点应力数组
        """
        # 计算高斯点应力
        gauss_stress = self.D @ self.B_global @ displacement

        # 将高斯点应力外推到节点
        node_stress = ExtrapolateMatrix3to3() @ gauss_stress
        return node_stress

    def ElementMass(self):
        """计算单元质量矩阵(待实现)"""
        pass

    def CalculateBasic(self):
        """计算基本参数(待实现)"""
        pass

    def ReCalculateElementStiffness(self):
        """
        重新计算单元刚度矩阵

        Returns:
            np.ndarray: 更新后的单元刚度矩阵
        """
        K = np.zeros((9, 9), dtype=float)
        for iii in range(len(self.integ)):
            K += self.B_global[iii].T @ self.D @ self.B_global[iii] * self.integ[iii]
        return K


class CPM4(ElementBaseClass, ABC):
    """
    4节点二次平面应变四边形单元(用于Cook膜单元)

    参考文献:
    1. A COMPATIBLE TRIANGULAR ELEMENT INCLUDING VERTEX ROTATIONS FOR PLANE ELASTICITY ANALYSIS - D.J.ALLMAN
    2. TECHNICAL NOTE ON THE ALLMAN TRIANGLE AND A RELATED QUADRILATERAL ELEMENT - Robert D.Cook
    3. A Refined Four-Noded Membrane Element With Rotational Degrees Of Freedom
    4. 《有限元法 理论、格式与求解方法》上 - K.J Bath P322

    特性:
    - 包含旋转自由度的四边形膜单元
    - 适用于平面弹性分析
    - 采用3x3高斯积分
    """

    def __init__(self, eid=None):
        """
        初始化4节点四边形单元

        Args:
            eid (int, optional): 单元ID. Defaults to None.
        """
        super().__init__(eid)
        self.nodes_count = 4  # 单元节点数
        self.vtu_type = "quad"  # VTK可视化类型
        self.T_matrix = None  # 全局坐标到局部坐标的转换矩阵(用于位移转换)
        self.integ = []  # 积分权重列表

    def CalElementDMatrix(self, an_type=None):
        """
        计算平面应变问题的本构矩阵

        参考:
        《有限元方法》Bathe 上册P184

        Args:
            an_type (optional): 分析类型. Defaults to None.
        """
        e = self.cha_dict[MaterialKey.E]  # 弹性模量
        niu = self.cha_dict[MaterialKey.Niu]  # 泊松比

        # 获取单元厚度
        if self.cha_dict.__contains__('RealConst') and len(self.cha_dict["RealConst"]) != 0:
            h = self.cha_dict["RealConst"][0]
        elif self.cha_dict.__contains__(MaterialKey.Thickness):
            h = self.cha_dict[MaterialKey.Thickness]
        else:
            raise KeyError("未找到厚度参数(RealConst或Thickness)")

        # 计算平面应变问题的弹性矩阵
        a = e / (1 - niu ** 2) * h
        self.D = a * np.array([
            [1, niu, 0],
            [niu, 1, 0],
            [0, 0, 0.5 * (1 - niu)]
        ], dtype=float)

    def ElementStiffness(self, from_origin=False):
        """
        计算单元刚度矩阵

        采用3x3高斯积分计算单元刚度矩阵
        包含将平面单元刚度矩阵转换为膜单元刚度矩阵的转换过程

        Args:
            from_origin (bool, optional): 是否从原始数据重新计算. Defaults to False.

        Returns:
            np.ndarray: 12x12的单元刚度矩阵
        """
        # 获取四边形积分点和权重(3x3积分)
        points, weights = GaussIntegrationPoint.GetSamplePointAndWeight(3)
        B_local = []

        # 在每个高斯点上计算刚度矩阵贡献
        for ri in range(3):
            for si in range(3):
                r, s = points[ri], points[si]
                w_r, w_s = weights[ri], weights[si]

                # 计算形函数对自然坐标的导数
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

                # 计算雅可比矩阵(描述几何变形)
                J_ij = np.array([
                    [np.matmul(phpr[:], self.node_coords[:, 0]), np.matmul(phpr[:], self.node_coords[:, 1])],
                    [np.matmul(phps[:], self.node_coords[:, 0]), np.matmul(phps[:], self.node_coords[:, 1])]
                ], dtype=float)

                det_J = np.linalg.det(J_ij)  # 雅可比行列式
                J_inv = np.array([
                    [J_ij[1, 1], -J_ij[0, 1]],
                    [-J_ij[1, 0], J_ij[0, 0]]
                ], dtype=float) / det_J  # 雅可比逆矩阵

                # 计算位移导数矩阵
                pupxy = np.array([
                    [ph1pr, 0, ph2pr, 0, ph3pr, 0, ph4pr, 0, ph5pr, 0, ph6pr, 0, ph7pr, 0, ph8pr, 0],
                    [ph1ps, 0, ph2ps, 0, ph3ps, 0, ph4ps, 0, ph5ps, 0, ph6ps, 0, ph7ps, 0, ph8ps, 0]
                ], dtype=float)

                pvpxy = np.array([
                    [0, ph1pr, 0, ph2pr, 0, ph3pr, 0, ph4pr, 0, ph5pr, 0, ph6pr, 0, ph7pr, 0, ph8pr],
                    [0, ph1ps, 0, ph2ps, 0, ph3ps, 0, ph4ps, 0, ph5ps, 0, ph6ps, 0, ph7ps, 0, ph8ps]
                ], dtype=float)

                # 计算应变-位移矩阵B
                B1 = np.matmul(J_inv, pupxy)
                B2 = np.matmul(J_inv, pvpxy)

                # 组装完整的B矩阵(3x16)
                B = np.insert(B1, 1, B2[1, :], axis=0)
                B[2, :] += B2[0, :]  # 添加剪切应变项

                B_local.append(B)
                self.integ.append(w_r * w_s * det_J)

                # 计算当前积分点对刚度矩阵的贡献
                if self.K is None:
                    self.K = B.T @ self.D @ B * w_r * w_s * det_J
                else:
                    self.K += B.T @ self.D @ B * w_r * w_s * det_J

        # 将平面单元刚度矩阵转换为膜单元刚度矩阵(参考Cook的论文)
        a12 = (self.node_coords[1, 0] - self.node_coords[0, 0]) * 0.125
        a23 = (self.node_coords[2, 0] - self.node_coords[1, 0]) * 0.125
        a34 = (self.node_coords[3, 0] - self.node_coords[2, 0]) * 0.125
        a41 = (self.node_coords[0, 0] - self.node_coords[3, 0]) * 0.125

        b12 = (self.node_coords[1, 1] - self.node_coords[0, 1]) * 0.125
        b23 = (self.node_coords[2, 1] - self.node_coords[1, 1]) * 0.125
        b34 = (self.node_coords[3, 1] - self.node_coords[2, 1]) * 0.125
        b41 = (self.node_coords[0, 1] - self.node_coords[3, 1]) * 0.125

        # 构造转换矩阵T
        T = np.array([
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
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
            [0, 0.5, -a41, 0, 0, 0, 0, 0, 0, 0, 0.5, a41]
        ], dtype=float)

        self.B_global = B_local @ T  # 全局坐标系下的B矩阵
        return T.T @ self.K @ T  # 转换后的刚度矩阵

    def CalculateElementStress(self, displacement):
        """
        计算单元应力

        Args:
            displacement (np.ndarray): 节点位移向量

        Returns:
            np.ndarray: 节点应力数组
        """
        # 获取高斯积分点坐标
        points, _ = GaussIntegrationPoint.GetSamplePointAndWeight(3)
        gauss_coords = [(points[ri], points[si]) for ri in range(3) for si in range(3)]

        # 计算形函数矩阵
        A = np.array([Quad4NodeShapeFunction(r, s) for r, s in gauss_coords])

        # 计算高斯点应力
        gauss_stress = self.D @ self.B_global @ displacement

        # 将高斯点应力外推到节点(最小二乘法)
        node_stress = np.linalg.solve(A.T @ A, A.T @ gauss_stress)
        return node_stress

    def ElementMass(self):
        """计算单元质量矩阵(待实现)"""
        pass

    def CalculateBasic(self):
        """计算基本参数(待实现)"""
        pass

    def ReCalculateElementStiffness(self):
        """
        重新计算单元刚度矩阵

        Returns:
            np.ndarray: 更新后的单元刚度矩阵
        """
        K = np.zeros((12, 12), dtype=float)
        for iii in range(len(self.integ)):
            K += self.B_global[iii].T @ self.D @ self.B_global[iii] * self.integ[iii]
        return K