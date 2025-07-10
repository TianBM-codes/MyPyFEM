#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys

from element.ElementBase import *
from femdb.ShapeFunctionsAndInteg import Quad4NodeShapeFunction, ExtrapolateMatrix4to4, ExtrapolateMatrix3to3
from abc import ABC


class DKTPlate(ElementBaseClass, ABC):
    """
    DKT (Discrete Kirchhoff Triangle) 三角形板单元类
    基于离散Kirchhoff理论的三节点三角形板弯曲单元

    参考文献:
    1. 有限单元法 王勖成 P364
    2. A_Study_of_Three-Node_Triangular_Plate_Bending_Elements.pdf
    3. Note And Explanation of Formulation for DKT Shell Element and its Implementation through Open Source Finite Software,Elmer.pdf
    4. Code_Aster Elements of plate: modelings DKT, DST, DKTG and Q4G

    TODO: 实现三角形积分，研究面积坐标和r,s参数坐标的区别(面积坐标见王勖成P366)
    TODO: 实现应力杂交单元
    """

    def __init__(self, eid=None):
        """
        初始化DKT板单元
        :param eid: 单元ID
        """
        super().__init__(eid)
        self.nodes_count = 3  # 三节点单元
        self.K = np.zeros([9, 9], dtype=float)  # 单元刚度矩阵(每个节点3个自由度)
        self.vtu_type = "triangle"  # VTK可视化类型
        self.T_matrix = None  # 全局坐标到局部坐标的转换矩阵
        self.B = []  # 应变-位移矩阵列表(每个积分点一个)
        self.integ = []  # 积分权重列表

    def CalElementDMatrix(self, an_type=None):
        """
        计算本构矩阵(弹性矩阵)
        根据材料属性和厚度计算板单元的弹性矩阵
        参考Bathe《有限元方法》上册P184

        :param an_type: 分析类型(未使用)
        """
        e = self.cha_dict[MaterialKey.E]  # 弹性模量
        niu = self.cha_dict[MaterialKey.Niu]  # 泊松比

        # 获取板厚度
        if self.cha_dict.__contains__('RealConst') and len(self.cha_dict["RealConst"]) != 0:
            h = self.cha_dict["RealConst"][0]
        elif self.cha_dict.__contains__(MaterialKey.Thickness):
            h = self.cha_dict[MaterialKey.Thickness]
        else:
            raise KeyError("未找到厚度参数(RealConst或Thickness)")

        # 计算板弯曲刚度
        a = e * h ** 3 / 12 / (1 - niu ** 2)

        # 构造弹性矩阵
        self.D = a * np.array([[1, niu, 0],
                               [niu, 1, 0],
                               [0, 0, 0.5 * (1 - niu)]], dtype=float)

    def ElementStiffness(self, from_origin=False):
        """
        计算单元刚度矩阵
        基于离散Kirchhoff理论的三节点三角形板单元刚度矩阵计算

        :param from_origin: 是否从原始数据重新计算
        :return: 单元刚度矩阵(9x9)
        """
        assert self.node_coords.shape == (3, 2)  # 确保节点坐标是3x2数组

        # 计算单元几何参数
        x12 = self.node_coords[0, 0] - self.node_coords[1, 0]
        x23 = self.node_coords[1, 0] - self.node_coords[2, 0]
        x31 = self.node_coords[2, 0] - self.node_coords[0, 0]

        y12 = self.node_coords[0, 1] - self.node_coords[1, 1]
        y23 = self.node_coords[1, 1] - self.node_coords[2, 1]
        y31 = self.node_coords[2, 1] - self.node_coords[0, 1]

        # 计算边长平方
        L4_square = x12 ** 2 + y12 ** 2
        L5_square = x23 ** 2 + y23 ** 2
        L6_square = x31 ** 2 + y31 ** 2

        # 计算形函数系数
        a4 = - x12 / L4_square
        a5 = - x23 / L5_square
        a6 = - x31 / L6_square

        b4 = 0.75 * x12 * y12 / L4_square
        b5 = 0.75 * x23 * y23 / L5_square
        b6 = 0.75 * x31 * y31 / L6_square

        c4 = (0.25 * x12 ** 2 - 0.5 * y12 ** 2) / L4_square
        c5 = (0.25 * x23 ** 2 - 0.5 * y23 ** 2) / L5_square
        c6 = (0.25 * x31 ** 2 - 0.5 * y31 ** 2) / L6_square

        d4 = - y12 / L4_square
        d5 = - y23 / L5_square
        d6 = - y31 / L6_square

        e4 = (0.25 * y12 ** 2 - 0.5 * x12 ** 2) / L4_square
        e5 = (0.25 * y23 ** 2 - 0.5 * x23 ** 2) / L5_square
        e6 = (0.25 * y31 ** 2 - 0.5 * x31 ** 2) / L6_square

        # 获取三角形积分点和权重(3点积分)
        sample_pt, weight = GaussIntegrationPoint.GetTrianglePointAndWeight(3)

        # 在3个高斯点上积分计算刚度矩阵
        for ii in range(len(sample_pt)):
            r, s = sample_pt[ii]

            # 计算形函数对自然坐标的导数
            pN1pr, pN2pr, pN3pr = -3 + 4 * (r + s), -1 + 4 * r, 0
            pN4pr, pN5pr, pN6pr = 4 * (1 - 2 * r - s), 4 * s, -4 * s
            pN1ps, pN2ps, pN3ps = pN1pr, 0, -1 + 4 * s
            pN4ps, pN5ps, pN6ps = -4 * r * (1 - r), 4 * s, r * (1 - r - 2 * s)

            # 计算弯曲应变矩阵B的分量
            pHxpr = np.asarray([
                1.5 * (a4 * pN4pr - a6 * pN6pr), b4 * pN4pr + b6 * pN6pr, pN1pr - c4 * pN4pr - c6 * pN6pr,
                1.5 * (a5 * pN5pr - a4 * pN4pr), b5 * pN5pr + b4 * pN4pr, pN2pr - c5 * pN5pr - c4 * pN4pr,
                1.5 * (a6 * pN6pr - a5 * pN5pr), b6 * pN6pr + b5 * pN5pr, pN3pr - c6 * pN6pr - c5 * pN5pr
            ], dtype=float)

            pHxps = np.asarray([
                1.5 * (a4 * pN4ps - a6 * pN6ps), b4 * pN4ps + b6 * pN6ps, pN1ps - c4 * pN4ps - c6 * pN6ps,
                1.5 * (a5 * pN5ps - a4 * pN4ps), b5 * pN5ps + b4 * pN4ps, pN2ps - c5 * pN5ps - c4 * pN4ps,
                1.5 * (a6 * pN6ps - a5 * pN5ps), b6 * pN6ps + b5 * pN5ps, pN3ps - c6 * pN6ps - c5 * pN5ps
            ], dtype=float)

            pHypr = np.asarray([
                1.5 * (d4 * pN4pr - d6 * pN6pr), -pN1pr + e4 * pN4pr + e6 * pN6pr, -b4 * pN4pr - b6 * pN6pr,
                1.5 * (d5 * pN5pr - d4 * pN4pr), -pN2pr + e5 * pN5pr + e4 * pN4pr, -b5 * pN5pr - b4 * pN4pr,
                1.5 * (d6 * pN6pr - d5 * pN5pr), -pN3pr + e6 * pN6pr + e5 * pN5pr, -b6 * pN6pr - b5 * pN5pr
            ], dtype=float)

            pHyps = np.asarray([
                1.5 * (d4 * pN4ps - d6 * pN6ps), -pN1ps + e4 * pN4ps + e6 * pN6ps, -pHxps[1],
                1.5 * (d5 * pN5ps - d4 * pN4ps), -pN2ps + e5 * pN5ps + e4 * pN4ps, -pHxps[4],
                1.5 * (d6 * pN6ps - d5 * pN5ps), -pN3ps + e6 * pN6ps + e5 * pN5ps, -pHxps[7]
            ], dtype=float)

            # 计算雅可比矩阵和行列式
            detJ = x31 * y12 - x12 * y31  # 三角形面积的两倍

            j11 = y31
            j12 = y12
            j21 = -x31
            j22 = -x12

            # 组装应变-位移矩阵B
            B = np.asarray([
                j11 * pHxpr + j12 * pHxps,
                j21 * pHypr + j22 * pHyps,
                j11 * pHypr + j12 * pHyps + j21 * pHxpr + j22 * pHxpr
            ], dtype=float)

            # 保存B矩阵和积分权重
            self.B.append(B)
            self.integ.append(weight[ii] / detJ)

            # 计算当前积分点对刚度矩阵的贡献
            self.K += B.T @ self.D @ B * weight[ii] / detJ

        return self.K

    def CalculateElementStress(self, displacement):
        """
        计算单元应力
        基于位移结果计算单元应力

        :param displacement: 节点位移向量
        :return: 节点应力数组
        """
        # 计算高斯点应力
        gauss_stress = self.D @ self.B @ displacement

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
        使用保存的B矩阵和积分权重重新计算刚度矩阵

        :return: 更新后的单元刚度矩阵
        """
        K = np.zeros((9, 9), dtype=float)
        for iii in range(len(self.integ)):
            K += self.B[iii].T @ self.D @ self.B[iii] * self.integ[iii]
        return K


class DKQPlate(ElementBaseClass, ABC):
    """
    DKQ (Discrete Kirchhoff Quadrilateral) 四边形板单元类
    基于离散Kirchhoff理论的四节点四边形板弯曲单元

    参考文献:
    1. Evaluation of a new quadrilateral thin plate bending element.pdf JEAN-LOUIS BATOZ
    """

    def __init__(self, eid=None):
        """
        初始化DKQ板单元
        :param eid: 单元ID
        """
        super().__init__(eid)
        self.nodes_count = 4  # 四节点单元
        self.K = np.zeros([12, 12], dtype=float)  # 单元刚度矩阵(每个节点3个自由度)
        self.vtu_type = "quad"  # VTK可视化类型
        self.T_matrix = None  # 全局坐标到局部坐标的转换矩阵
        self.B = []  # 应变-位移矩阵列表(每个积分点一个)
        self.integ = []  # 积分权重列表

    def CalElementDMatrix(self, an_type=None):
        """
        计算本构矩阵(弹性矩阵)
        根据材料属性和厚度计算板单元的弹性矩阵
        参考Bathe《有限元方法》上册P184

        :param an_type: 分析类型(未使用)
        """
        e = self.cha_dict[MaterialKey.E]  # 弹性模量
        niu = self.cha_dict[MaterialKey.Niu]  # 泊松比

        # 获取板厚度
        if self.cha_dict.__contains__('RealConst') and len(self.cha_dict["RealConst"]) != 0:
            h = self.cha_dict["RealConst"][0]
        elif self.cha_dict.__contains__(MaterialKey.Thickness):
            h = self.cha_dict[MaterialKey.Thickness]
        else:
            raise KeyError("未找到厚度参数(RealConst或Thickness)")

        # 计算板弯曲刚度
        a = e * h ** 3 / 12 / (1 - niu ** 2)

        # 构造弹性矩阵
        self.D = a * np.array([[1, niu, 0],
                               [niu, 1, 0],
                               [0, 0, 0.5 * (1 - niu)]], dtype=float)

    def ElementStiffness(self, from_origin=False):
        """
        计算单元刚度矩阵
        基于离散Kirchhoff理论的四节点四边形板单元刚度矩阵计算

        :param from_origin: 是否从原始数据重新计算
        :return: 单元刚度矩阵(12x12)
        """
        assert self.node_coords.shape == (4, 2)  # 确保节点坐标是4x2数组

        # 计算单元几何参数
        x12 = self.node_coords[0, 0] - self.node_coords[1, 0]
        x23 = self.node_coords[1, 0] - self.node_coords[2, 0]
        x34 = self.node_coords[2, 0] - self.node_coords[3, 0]
        x41 = self.node_coords[3, 0] - self.node_coords[0, 0]
        x13 = self.node_coords[0, 0] - self.node_coords[2, 0]
        x24 = self.node_coords[1, 0] - self.node_coords[3, 0]

        y12 = self.node_coords[0, 1] - self.node_coords[1, 1]
        y23 = self.node_coords[1, 1] - self.node_coords[2, 1]
        y34 = self.node_coords[2, 1] - self.node_coords[3, 1]
        y41 = self.node_coords[3, 1] - self.node_coords[0, 1]
        y13 = self.node_coords[0, 1] - self.node_coords[2, 1]
        y24 = self.node_coords[1, 1] - self.node_coords[3, 1]

        # 计算边长平方
        L5_square = x12 ** 2 + y12 ** 2
        L6_square = x23 ** 2 + y23 ** 2
        L7_square = x34 ** 2 + y34 ** 2
        L8_square = x41 ** 2 + y41 ** 2

        # 计算形函数系数
        a5 = - x12 / L5_square
        a6 = - x23 / L6_square
        a7 = - x34 / L7_square
        a8 = - x41 / L8_square

        b5 = 0.75 * x12 * y12 / L5_square
        b6 = 0.75 * x23 * y23 / L6_square
        b7 = 0.75 * x34 * y34 / L7_square
        b8 = 0.75 * x41 * y41 / L8_square

        c5 = (0.25 * x12 ** 2 - 0.5 * y12 ** 2) / L5_square
        c6 = (0.25 * x23 ** 2 - 0.5 * y23 ** 2) / L6_square
        c7 = (0.25 * x34 ** 2 - 0.5 * y34 ** 2) / L7_square
        c8 = (0.25 * x41 ** 2 - 0.5 * y41 ** 2) / L8_square

        d5 = - y12 / L5_square
        d6 = - y23 / L6_square
        d7 = - y34 / L7_square
        d8 = - y41 / L8_square

        e5 = (0.25 * y12 ** 2 - 0.5 * x12 ** 2) / L5_square
        e6 = (0.25 * y23 ** 2 - 0.5 * x23 ** 2) / L6_square
        e7 = (0.25 * y34 ** 2 - 0.5 * x34 ** 2) / L7_square
        e8 = (0.25 * y41 ** 2 - 0.5 * x41 ** 2) / L8_square

        # 定义4点高斯积分(2x2)
        sample_pt_r = (-0.577350269189626, 0.577350269189626, 0.577350269189626, -0.577350269189626)
        sample_pt_s = (0.577350269189626, 0.577350269189626, -0.577350269189626, -0.577350269189626)
        weight = (1, 1, 1, 1)

        # 在4个高斯点上积分计算刚度矩阵
        for ii in range(4):
            r, s = sample_pt_r[ii], sample_pt_s[ii]

            # 计算形函数对自然坐标的导数
            pN1pr = 0.25 * (s ** 2 + s) + 0.5 * (1 + s) * r
            pN2pr = -0.25 * (s ** 2 + s) + 0.5 * (1 + s) * r
            pN3pr = 0.25 * (s - s ** 2) + 0.5 * r * (1 - s)
            pN4pr = 0.25 * (s ** 2 - s) + 0.5 * r * (1 - s)
            pN5pr = -r * (1 + s)
            pN6pr = 0.5 * (s ** 2 - 1)
            pN7pr = r * (s - 1)
            pN8pr = 0.5 * (1 - s ** 2)

            pN1ps = 0.25 * (r ** 2 + r) + 0.5 * s * (1 + r)
            pN2ps = 0.25 * (r ** 2 - r) + 0.5 * s * (1 - r)
            pN3ps = 0.25 * (r - r ** 2) + 0.5 * s * (1 - r)
            pN4ps = -0.25 * (r ** 2 + r) + 0.5 * s * (1 + r)
            pN5ps = 0.5 * (1 - r ** 2)
            pN6ps = s * (r - 1)
            pN7ps = 0.5 * (r ** 2 - 1)
            pN8ps = -s * (1 + r)

            # 计算弯曲应变矩阵B的分量
            pHxpr = np.asarray([
                1.5 * (a5 * pN5pr - a8 * pN8pr), b5 * pN5pr + b8 * pN8pr, pN1pr - c5 * pN5pr - c8 * pN8pr,
                1.5 * (a6 * pN6pr - a5 * pN5pr), b6 * pN6pr + b5 * pN5pr, pN2pr - c6 * pN6pr - c5 * pN5pr,
                1.5 * (a7 * pN7pr - a6 * pN6pr), b7 * pN7pr + b6 * pN6pr, pN3pr - c7 * pN7pr - c6 * pN6pr,
                1.5 * (a8 * pN8pr - a7 * pN7pr), b8 * pN8pr + b7 * pN7pr, pN4pr - c8 * pN8pr - c7 * pN7pr
            ], dtype=float)

            pHxps = np.asarray([
                1.5 * (a5 * pN5ps - a8 * pN8ps), b5 * pN5ps + b8 * pN8ps, pN1ps - c5 * pN5ps - c8 * pN8ps,
                1.5 * (a6 * pN6ps - a5 * pN5ps), b6 * pN6ps + b5 * pN5ps, pN2ps - c6 * pN6ps - c5 * pN5ps,
                1.5 * (a7 * pN7ps - a6 * pN6ps), b7 * pN7ps + b6 * pN6ps, pN3ps - c7 * pN7ps - c6 * pN6ps,
                1.5 * (a8 * pN8ps - a7 * pN7ps), b8 * pN8ps + b7 * pN7ps, pN4ps - c8 * pN8ps - c7 * pN7ps
            ], dtype=float)

            pHypr = np.asarray([
                1.5 * (d5 * pN5pr - d8 * pN8pr), -pN1pr + e5 * pN5pr + e8 * pN8pr, -b5 * pN5pr - b8 * pN8pr,
                1.5 * (d6 * pN6pr - d5 * pN5pr), -pN2pr + e6 * pN6pr + e5 * pN5pr, -b6 * pN6pr - b5 * pN5pr,
                1.5 * (d7 * pN7pr - d6 * pN6pr), -pN3pr + e7 * pN7pr + e6 * pN6pr, -b7 * pN7pr - b6 * pN6pr,
                1.5 * (d8 * pN8pr - d7 * pN7pr), -pN4pr + e8 * pN8pr + e7 * pN7pr, -b8 * pN8pr - b7 * pN7pr
            ], dtype=float)

            pHyps = np.asarray([
                1.5 * (d5 * pN5ps - d8 * pN8ps), -pN1ps + e5 * pN5ps + e8 * pN8ps, -pHxps[1],
                1.5 * (d6 * pN6ps - d5 * pN5ps), -pN2ps + e6 * pN6ps + e5 * pN5ps, -pHxps[4],
                1.5 * (d7 * pN7ps - d6 * pN6ps), -pN3ps + e7 * pN7ps + e6 * pN6ps, -pHxps[7],
                1.5 * (d8 * pN8ps - d7 * pN7ps), -pN4ps + e8 * pN8ps + e7 * pN7ps, -pHxps[10]
            ], dtype=float)

            # 计算雅可比矩阵和行列式
            J11 = 0.25 * (x12 - x34 + r * (x12 + x34))
            J12 = 0.25 * (y12 - y34 + r * (y12 + y34))
            J21 = 0.25 * (x23 - x41 + s * (x12 + x34))
            J22 = 0.25 * (y23 - y41 + s * (y12 + y34))
            detJ = 0.125 * (x13 * y24 - x24 * y13 - r * (x12 * y34 + x34 * y12) + s * (x41 * y23 - x23 * y41))

            j11 = J22
            j12 = -J12
            j21 = -J21
            j22 = J11

            # 组装应变-位移矩阵B
            B = np.asarray([
                j11 * pHxpr + j12 * pHxps,
                j21 * pHypr + j22 * pHyps,
                j11 * pHypr + j12 * pHyps + j21 * pHxpr + j22 * pHxps
            ], dtype=float)

            # 保存B矩阵和积分权重
            self.B.append(B)
            self.integ.append(weight[ii] / detJ)

            # 计算当前积分点对刚度矩阵的贡献
            self.K += B.T @ self.D @ B * weight[ii] / detJ

        return self.K

    def CalculateElementStress(self, displacement):
        """
        计算单元应力
        基于位移结果计算单元应力

        :param displacement: 节点位移向量
        :return: 节点应力数组
        """
        # 计算高斯点应力
        gauss_stress = self.D @ self.B @ displacement

        # 将高斯点应力外推到节点
        node_stress = ExtrapolateMatrix4to4() @ gauss_stress

        return node_stress

    def ElementMass(self):
        """计算单元质量矩阵"""
        pass

    def CalculateBasic(self):
        """计算基本参数"""
        pass

    def ReCalculateElementStiffness(self):
        """
        重新计算单元刚度矩阵
        使用保存的B矩阵和积分权重重新计算刚度矩阵

        :return: 更新后的单元刚度矩阵
        """
        K = np.zeros((12, 12), dtype=float)
        for iii in range(len(self.integ)):
            K += self.B[iii].T @ self.D @ self.B[iii] * self.integ[iii]
        return K