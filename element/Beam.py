#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys

from element.ElementBase import *
from typing import List
import numpy as np
from abc import ABC


class BeamCalculator:
    """
    计算梁单元截面属性

    Reference:
    1. https://www.omnicalculator.com/physics/torsional-constant
    2. https://www.jpe-innovations.com/precision-point/beam-theory-torsion/
    3. https://en.wikipedia.org/wiki/List_of_second_moments_of_area
    4. https://amesweb.info/section/i-beam-moment-of-inertia-calculator.aspx
    5. https://www.structuralbasics.com/moment-of-inertia-formulas/
    6. http://manual.midasuser.com/EN_Common/Gen/790/Start/04_Model/05_Properties/Section.htm

    Attention:
    工字梁只支持上下边等长并且上下厚度一致的, 即梁截面是中心对称的, equal flanges
    """

    @staticmethod
    def CalculateMomentOfInertiaOfArea(sec_type: BeamSectionType, sec_data: List[float]) -> dict:
        """
        计算梁横截面惯性矩, 抗扭刚度
        """
        if sec_type == BeamSectionType.I:
            assert len(sec_data) == 6
            w, h, tf, tw = sec_data[0], sec_data[2], sec_data[3], sec_data[5]
            It = w * pow(h, 3) / 12 - (w - tw) * pow((h - 2 * tf), 3) / 12
            Is = (h - 2 * tf) * pow(tw, 3) / 12 + tf * pow(w, 3) / 6

            a, b, c, d = sec_data[0], sec_data[3], sec_data[2] - 2 * sec_data[3], sec_data[5]
            K1 = a * pow(b, 3) / 3 - 0.21 * pow(b, 4) - 0.0175 * pow(b, 8) / pow(a, 4)
            K2 = c * pow(d, 3) / 3
            D = (pow(b, 2) + d * d * 0.25) / b
            if b < d:
                alpha = b / d * 0.15
            else:
                alpha = d / b * 0.15
            Torsional = 2 * K1 + K2 + 2 * alpha * pow(D, 4)
            return {SectionKey.It: It, SectionKey.Is: Is, SectionKey.Tor: Torsional}

        elif sec_type == BeamSectionType.Rectangle:
            assert len(sec_data) == 2
            a = max(sec_data)
            b = min(sec_data)
            It = b ** 3 * a / 12
            Is = a ** 3 * b / 12
            if a == b:
                Torsional = 9 * pow(a, 4) / 64
            else:
                Torsional = a * pow(b, 3) * (16 / 3 - 3.36 * b / a * (1 - pow(b / a, 4) / 12)) / 16
            return {SectionKey.It: It, SectionKey.Is: Is, SectionKey.Tor: Torsional}

        elif sec_type == BeamSectionType.CircleSolid:
            """
            sec_data: R, N, T
            R = Radius
            N = Number of divisions around the circumference; 8  N  120 (where a greater value improves accuracy slightly); default = 8
            T = Number of divisions through the radius; default = 2
            """
            I = 0.25 * np.pi * pow(sec_data[0], 4)
            Torsional = I * 2
            return {SectionKey.It: I, SectionKey.Is: I, SectionKey.Tor: Torsional}

        elif sec_type == BeamSectionType.UserInput:
            """
            A, Tor, It, Is
            """
            return {SectionKey.Area: sec_data[0], SectionKey.Tor: sec_data[1], SectionKey.It: sec_data[2], SectionKey.Is: sec_data[3]}

        else:
            mlogger.fatal("UnSupport Beam Section Type:{}".format(sec_type))
            sys.exit(1)

    @staticmethod
    def CalEffectiveShearArea(sec_type: BeamSectionType, sec_data: List[float]) -> dict:
        """
        计算截面的面积属性, 包括面积、两个方向的抗剪等效面积, 圆的输入是半径
        """
        if sec_type == BeamSectionType.Rectangle:
            assert len(sec_data) == 2
            Area = sec_data[0] * sec_data[1]
            E_Area = 5 / 6 * sec_data[0] * sec_data[1]
            return {SectionKey.Area: Area, SectionKey.At: E_Area, SectionKey.As: E_Area}

        elif sec_type == BeamSectionType.CircleSolid:
            Area = np.pi * pow(sec_data[0], 2)
            E_Area = 0.9 * np.pi * sec_data[0] ** 2
            return {SectionKey.Area: Area, SectionKey.At: E_Area, SectionKey.As: E_Area}

        elif sec_type == BeamSectionType.I:
            assert len(sec_data) == 6
            h = sec_data[2] - 2 * sec_data[3]
            Area = 2 * sec_data[0] * sec_data[3] + h * sec_data[5]
            E_Area_t = sec_data[0] * sec_data[5]
            E_Area_s = 5 / 6 * (2 * sec_data[0] * sec_data[3])
            return {SectionKey.Area: Area, SectionKey.As: E_Area_t, SectionKey.At: E_Area_s}

        else:
            mlogger.fatal("UnSupport Section Type: {}".format(sec_type))
            sys.exit(1)


class Beam188_old(ElementBaseClass, ABC):
    """
    一般直梁, 两节点混合插值单元, 横向剪应变被假设为常量(一阶高斯点处的剪应变)
    Reference:
    1. 《有限元法 理论、格式与求解方法》上 Bathe P382
    """

    def __init__(self, eid):
        super().__init__(eid)
        self.nodes_count = 2  # 每个单元包含2个节点, 表示方向的辅助节点不计算在内
        self.vtu_type = "line"
        self.stiffness = None
        self.stress = None
        self.block_size = 144
        self.node_dof_count = 6

        # 截面相关, 面积、有效面积
        self.It, self.Is, self.Tor = None, None, None
        self.A, self.effect_as, self.effect_at = None, None, None

    def CalElementDMatrix(self, an_type=None):
        pass

    def ElementStiffness(self, from_origin=False):
        """
        TODO: 有整理的pdf
        Reference:
        """
        delta = np.asarray(self.node_coords[1, :] - self.node_coords[0, :])
        E = self.cha_dict[MaterialKey.E]
        A = self.cha_dict[SectionKey.Area]
        At = self.cha_dict[SectionKey.At]
        As = self.cha_dict[SectionKey.As]
        G = self.cha_dict[MaterialKey.G]
        L = np.linalg.norm(delta)
        It = self.cha_dict[SectionKey.It]
        Is = self.cha_dict[SectionKey.Is]
        Tor = self.cha_dict[SectionKey.Tor]

        # 计算弯曲和剪切的刚度阵, 并组装成矩阵
        EA_L = E * A / L
        GAt_L = G * At / L
        GAs_L = G * As / L
        GAt_2 = G * At * 0.5
        GAs_2 = G * As * 0.5
        GAtL_4 = G * At * L * 0.25
        GAsL_4 = G * As * L * 0.25
        EIt_L = E * It / L
        EIs_L = E * Is / L
        K = np.zeros((12, 12), dtype=float)

        # 轴力因素
        K[0, 0], K[0, 6] = EA_L, -EA_L
        K[6, 0], K[6, 6] = -EA_L, EA_L

        # 扭转因素
        K[3, 3], K[3, 9] = Tor, -Tor
        K[9, 3], K[9, 9] = -Tor, Tor

        # 剪切因素
        K[1, 1], K[1, 7], K[7, 1], K[7, 7] = GAt_L, -GAt_L, -GAt_L, GAt_L
        K[1, 5], K[1, 11], K[7, 5], K[7, 11] = GAt_2, GAt_2, -GAt_2, -GAt_2
        K[5, 1], K[5, 7], K[11, 1], K[11, 7] = GAt_2, -GAt_2, GAt_2, -GAt_2
        K[5, 5], K[5, 11], K[11, 5], K[11, 11] = GAtL_4, GAtL_4, GAtL_4, GAtL_4

        K[2, 2], K[2, 8], K[8, 2], K[8, 8] = GAs_L, -GAs_L, -GAs_L, GAs_L
        K[2, 4], K[2, 10], K[8, 4], K[8, 10] = GAs_2, GAs_2, -GAs_2, -GAs_2
        K[4, 2], K[4, 8], K[10, 2], K[10, 8] = GAs_2, -GAs_2, GAs_2, -GAs_2
        K[4, 4], K[4, 10], K[10, 4], K[10, 10] = GAsL_4, GAsL_4, GAsL_4, GAsL_4

        # 弯曲因素
        K[4, 4] += EIs_L
        K[5, 5] += EIt_L
        K[10, 10] += EIs_L
        K[11, 11] += EIt_L

        K[4, 10] -= EIs_L
        K[5, 11] -= EIt_L
        K[10, 4] -= EIs_L
        K[11, 5] -= EIt_L

        """
        笛卡尔坐标变换为自然坐标, 梁主轴方向r, 梁方向点方向s, 王勖成P331
        如果没有第三个方向点, 使用默认参考方向, 如果ref与x_local共线, 则选择其他向量
        """
        x_local = delta / L
        if self.node_coords.shape[0] < 3:
            ref = np.array([0, 0, 1])
            if np.allclose(np.cross(x_local, ref), 0):
                ref = np.array([0, 1, 0])
            normal_direct = ref
        else:
            normal_direct = np.asarray(self.node_coords[2, :] - (self.node_coords[0, :] + self.node_coords[1, :]) / 2)

        z_local = np.cross(x_local, normal_direct)
        z_local /= np.linalg.norm(z_local)
        y_local = np.cross(z_local, x_local)
        y_local /= np.linalg.norm(y_local)

        """
        3x3 方向余弦矩阵, 对应每个自由度组 (0~3, 3~6, 6~9, 9~12)
        """
        trans_mat = np.zeros((12, 12))
        rot_mat = np.vstack([x_local, y_local, z_local])

        for i in range(4):
            trans_mat[i * 3:(i + 1) * 3, i * 3:(i + 1) * 3] = rot_mat

        return trans_mat.T @ K @ trans_mat

    def CalculateElementStress(self, displacement):
        """
        Calculate element stress
        """
        # self.stress = self.e / self.rod_length * (np.dot(np.asarray(x2), self.cos_angel) - np.dot(np.asarray(x1), self.cos_angel))
        return np.zeros((2, 1)), np.zeros((2, 1)), np.zeros((2, 1)), np.zeros((2, 1)), np.zeros((2, 1)), np.zeros((2, 1))

    def ElementMass(self):
        pass

    def CalculateBasic(self):
        pass

    def ReCalculateElementStiffness(self):
        pass


class Beam189(ElementBaseClass, ABC):
    """
    Beam189 Element class
    cdb节点格式：
    起始节点编号--终点节点编号--中间节点编号--方向节点编号
    方向节点编号位于中间节点编号的正上方
    """

    def __init__(self, eid):
        super().__init__(eid)
        self.nodes_count = 3  # Each element has 3 nodes
        self.vtp_type = "line3"
        self.stiffness = None
        self.stress = None
        self.I = None  # 惯性矩
        self.sec_type = None  # 截面类型
        self.sec_data = None  # 截面参数
        self.block_size = 324
        self.node_dof_count = 6

    def CalElementDMatrix(self, an_type=None):
        pass

    def ElementStiffness(self, from_origin=False):
        """
        TODO: 有整理的pdf
        Reference:
        """
        # 4个节点, 起始点, 终点, 中间点, 方向点
        assert self.node_coords.shape == (4, 3)

        # 单元参数
        delta = np.asarray(np.diff(self.node_coords, axis=0))[0]
        E = self.cha_dict[MaterialKey.E]
        A = self.cha_dict[MaterialKey.Area]
        G = self.cha_dict[MaterialKey.G]
        L = np.sqrt(np.dot(delta.T, delta))
        I = BeamCalculator.CalculateMomentOfInertiaOfArea(self.sec_type, self.sec_data)
        k = BeamCalculator.GetEquivalentCoff(self.sec_type)

        # 计算弯曲和剪切的刚度部分, 并组装成矩阵
        spt, weight = GaussIntegrationPoint.GetSamplePointAndWeight(2)
        K = np.zeros((6, 6), dtype=float)
        for i in range(2):
            Ks1 = 2 * spt[i] / L - 1 / L
            Ks2 = spt[i] - 1 / 12
            Ks3 = -4 * spt[i] / L
            Ks4 = -2 / 3
            Ks5 = 1 / L + 2 * spt[i] / L
            Ks6 = -1 / 6 - spt[i] * 0.5
            Ks = np.array([[Ks1 * Ks1, Ks1 * Ks2, Ks1 * Ks3, Ks1 * Ks4, Ks1 * Ks5, Ks1 * Ks6],
                           [Ks2 * Ks1, Ks2 * Ks2, Ks2 * Ks3, Ks2 * Ks4, Ks2 * Ks5, Ks2 * Ks6],
                           [Ks3 * Ks1, Ks3 * Ks2, Ks3 * Ks3, Ks3 * Ks4, Ks3 * Ks5, Ks3 * Ks6],
                           [Ks4 * Ks1, Ks4 * Ks2, Ks4 * Ks3, Ks4 * Ks4, Ks4 * Ks5, Ks4 * Ks6],
                           [Ks5 * Ks1, Ks5 * Ks2, Ks5 * Ks3, Ks5 * Ks4, Ks5 * Ks5, Ks5 * Ks6],
                           [Ks6 * Ks1, Ks6 * Ks2, Ks6 * Ks3, Ks6 * Ks4, Ks6 * Ks5, Ks6 * Ks6]])
            K += weight[i] * Ks

        # 几何关系, 笛卡尔坐标变换为自然坐标
        Cx, Cy, Cz = delta / L
        Cx2, Cy2, Cz2 = Cx ** 2, Cy ** 2, Cz ** 2
        CxCy, CxCz, CyCz = Cx * Cy, Cx * Cz, Cy * Cz

        return E * A / L * np.array([[Cx2, CxCy, CxCz, -Cx2, -CxCy, -CxCz],
                                     [CxCy, Cy2, CyCz, -CxCy, -Cy2, -CyCz],
                                     [CxCz, CyCz, Cz2, -CxCz, -CyCz, -Cz2],
                                     [-Cx2, -CxCy, -CxCz, Cx2, CxCy, CxCz],
                                     [-CxCy, -Cy2, -CyCz, CxCy, Cy2, CyCz],
                                     [-CxCz, -CyCz, -Cz2, CxCz, CyCz, Cz2]])

    def CalculateElementStress(self, displacement):
        """
        Calculate element stress
        """
        # fem_database = Domain()
        # x1 = fem_database.GetDisplacement(self.search_node_id[0])
        # x2 = fem_database.GetDisplacement(self.search_node_id[1])
        # self.stress = self.e / self.rod_length * (np.dot(np.asarray(x2), self.cos_angel) - np.dot(np.asarray(x1), self.cos_angel))

    def ElementMass(self):
        pass

    def CalculateBasic(self):
        pass

    def ReCalculateElementStiffness(self):
        pass


class Beam188(ElementBaseClass, ABC):
    """
    3D线性梁单元，支持端部释放条件
    """

    def __init__(self, eid):
        super().__init__(eid)
        self.nodes_count = 2  # 每个单元包含2个节点, 表示方向的辅助节点不计算在内
        self.vtu_type = "line"
        self.stiffness = None
        self.stress = None
        self.block_size = 144
        self.node_dof_count = 6

    def _calculate_length(self):
        """计算单元长度"""
        delta = self.node_coords[1] - self.node_coords[0]
        return np.linalg.norm(delta)

    def set_material_property(self, key, value):
        """设置材料属性"""
        self.cha_dict[key] = value

    def set_section_property(self, key, value):
        """设置截面属性"""
        self.cha_dict[key] = value

    def ElementStiffness(self, from_origin=False):
        """计算单元刚度矩阵"""
        E = self.cha_dict[MaterialKey.E]
        G = self.cha_dict[MaterialKey.G]
        A = self.cha_dict[SectionKey.Area]
        Jx = self.cha_dict[SectionKey.Tor]
        Iy = self.cha_dict[SectionKey.It]
        Iz = self.cha_dict[SectionKey.Is]
        L = self.length = self._calculate_length()
        self._calculate_rotation_matrix()

        oneOverL = 1.0 / L
        EoverL = E * oneOverL
        kb = np.zeros((6, 6))
        kb[0, 0] = A * EoverL
        kb[5, 5] = G * Jx * oneOverL
        EIzoverL4 = 4.0 * E * Iz / L
        EIzoverL2 = 2.0 * E * Iz / L
        kb[1, 1] = kb[2, 2] = EIzoverL4
        kb[1, 2] = kb[2, 1] = EIzoverL2

        EIyoverL4 = 4.0 * E * Iy / L
        EIyoverL2 = 2.0 * E * Iy / L
        kb[3, 3] = kb[4, 4] = EIyoverL4
        kb[3, 4] = kb[4, 3] = EIyoverL2

        # 第一步矩阵乘法: tmp = kb * T_bl
        tmp = np.zeros((12, 12))
        for i in range(6):
            tmp[i, 0] = -kb[i, 0]
            tmp[i, 1] = oneOverL * (kb[i, 1] + kb[i, 2])
            tmp[i, 2] = -oneOverL * (kb[i, 3] + kb[i, 4])
            tmp[i, 3] = -kb[i, 5]
            tmp[i, 4] = kb[i, 3]
            tmp[i, 5] = kb[i, 1]
            tmp[i, 6] = kb[i, 0]
            tmp[i, 7] = -tmp[i, 1]
            tmp[i, 8] = -tmp[i, 2]
            tmp[i, 9] = kb[i, 5]
            tmp[i, 10] = kb[i, 4]
            tmp[i, 11] = kb[i, 2]

        # 第二步矩阵乘法: kl = T_bl^T * tmp
        kl = np.zeros((12, 12))
        for i in range(12):
            kl[0, i] = -tmp[0, i]
            kl[1, i] = oneOverL * (tmp[1, i] + tmp[2, i])
            kl[2, i] = -oneOverL * (tmp[3, i] + tmp[4, i])
            kl[3, i] = -tmp[5, i]
            kl[4, i] = tmp[3, i]
            kl[5, i] = tmp[1, i]
            kl[6, i] = tmp[0, i]
            kl[7, i] = -kl[1, i]
            kl[8, i] = -kl[2, i]
            kl[9, i] = tmp[5, i]
            kl[10, i] = tmp[4, i]
            kl[11, i] = tmp[2, i]

        # 3. 将局部刚度矩阵转换为全局刚度矩阵
        # 第一步: tmp2 = kl * T_lg
        tmp2 = np.zeros((12, 12))
        for m in range(12):
            # 前3个自由度 (平移)
            tmp2[m, 0] = kl[m, 0] * self.R[0, 0] + kl[m, 1] * self.R[1, 0] + kl[m, 2] * self.R[2, 0]
            tmp2[m, 1] = kl[m, 0] * self.R[0, 1] + kl[m, 1] * self.R[1, 1] + kl[m, 2] * self.R[2, 1]
            tmp2[m, 2] = kl[m, 0] * self.R[0, 2] + kl[m, 1] * self.R[1, 2] + kl[m, 2] * self.R[2, 2]

            # 中间3个自由度 (节点I的旋转)
            tmp2[m, 3] = kl[m, 3] * self.R[0, 0] + kl[m, 4] * self.R[1, 0] + kl[m, 5] * self.R[2, 0]
            tmp2[m, 4] = kl[m, 3] * self.R[0, 1] + kl[m, 4] * self.R[1, 1] + kl[m, 5] * self.R[2, 1]
            tmp2[m, 5] = kl[m, 3] * self.R[0, 2] + kl[m, 4] * self.R[1, 2] + kl[m, 5] * self.R[2, 2]

            # 后3个自由度 (节点J的平移)
            tmp2[m, 6] = kl[m, 6] * self.R[0, 0] + kl[m, 7] * self.R[1, 0] + kl[m, 8] * self.R[2, 0]
            tmp2[m, 7] = kl[m, 6] * self.R[0, 1] + kl[m, 7] * self.R[1, 1] + kl[m, 8] * self.R[2, 1]
            tmp2[m, 8] = kl[m, 6] * self.R[0, 2] + kl[m, 7] * self.R[1, 2] + kl[m, 8] * self.R[2, 2]

            # 最后3个自由度 (节点J的旋转)
            tmp2[m, 9] = kl[m, 9] * self.R[0, 0] + kl[m, 10] * self.R[1, 0] + kl[m, 11] * self.R[2, 0]
            tmp2[m, 10] = kl[m, 9] * self.R[0, 1] + kl[m, 10] * self.R[1, 1] + kl[m, 11] * self.R[2, 1]
            tmp2[m, 11] = kl[m, 9] * self.R[0, 2] + kl[m, 10] * self.R[1, 2] + kl[m, 11] * self.R[2, 2]

        # 第二步: kg = T_lg^T * tmp2
        kg = np.zeros((12, 12))
        for m in range(12):
            # 前3个自由度 (平移)
            kg[0, m] = self.R[0, 0] * tmp2[0, m] + self.R[1, 0] * tmp2[1, m] + self.R[2, 0] * tmp2[2, m]
            kg[1, m] = self.R[0, 1] * tmp2[0, m] + self.R[1, 1] * tmp2[1, m] + self.R[2, 1] * tmp2[2, m]
            kg[2, m] = self.R[0, 2] * tmp2[0, m] + self.R[1, 2] * tmp2[1, m] + self.R[2, 2] * tmp2[2, m]

            # 中间3个自由度 (节点I的旋转)
            kg[3, m] = self.R[0, 0] * tmp2[3, m] + self.R[1, 0] * tmp2[4, m] + self.R[2, 0] * tmp2[5, m]
            kg[4, m] = self.R[0, 1] * tmp2[3, m] + self.R[1, 1] * tmp2[4, m] + self.R[2, 1] * tmp2[5, m]
            kg[5, m] = self.R[0, 2] * tmp2[3, m] + self.R[1, 2] * tmp2[4, m] + self.R[2, 2] * tmp2[5, m]

            # 后3个自由度 (节点J的平移)
            kg[6, m] = self.R[0, 0] * tmp2[6, m] + self.R[1, 0] * tmp2[7, m] + self.R[2, 0] * tmp2[8, m]
            kg[7, m] = self.R[0, 1] * tmp2[6, m] + self.R[1, 1] * tmp2[7, m] + self.R[2, 1] * tmp2[8, m]
            kg[8, m] = self.R[0, 2] * tmp2[6, m] + self.R[1, 2] * tmp2[7, m] + self.R[2, 2] * tmp2[8, m]

            # 最后3个自由度 (节点J的旋转)
            kg[9, m] = self.R[0, 0] * tmp2[9, m] + self.R[1, 0] * tmp2[10, m] + self.R[2, 0] * tmp2[11, m]
            kg[10, m] = self.R[0, 1] * tmp2[9, m] + self.R[1, 1] * tmp2[10, m] + self.R[2, 1] * tmp2[11, m]
            kg[11, m] = self.R[0, 2] * tmp2[9, m] + self.R[1, 2] * tmp2[10, m] + self.R[2, 2] * tmp2[11, m]

        return kg

    def _calculate_rotation_matrix(self):
        """计算旋转矩阵"""
        L = self.length
        delta = self.node_coords[1] - self.node_coords[0]
        x_local = delta / L  # 梁轴向方向

        # 确定参考方向
        if self.node_coords.shape[0] < 3:
            # 如果没有第三个节点，使用默认参考方向
            ref = np.array([0, 0, 1])
            if np.allclose(np.cross(x_local, ref), 0):
                ref = np.array([0, 1, 0])
            normal_direct = ref
        else:
            # 如果有第三个节点，使用它确定参考方向
            mid_point = (self.node_coords[0] + self.node_coords[1]) / 2
            normal_direct = self.node_coords[2] - mid_point

        # 构建局部坐标系
        y_local = np.cross(normal_direct, x_local)
        y_local /= np.linalg.norm(y_local)
        z_local = np.cross(x_local, y_local)
        z_local /= np.linalg.norm(z_local)

        # 构建旋转矩阵
        self.R = np.vstack([x_local, y_local, z_local])

    def ElementMass(self):
        """计算单元质量矩阵"""
        cMass = True
        rho = self.cha_dict.get(MaterialKey.Rho, 0.0)
        if rho <= 0.0:
            self.mass = np.zeros((12, 12))
            return self.mass

        L = self.length
        A = self.cha_dict[SectionKey.Area]
        Jx = self.cha_dict[SectionKey.Tor]
        Iy = self.cha_dict[SectionKey.It]
        Iz = self.cha_dict[SectionKey.Is]

        if not cMass:
            # 集中质量矩阵
            m = 0.5 * rho * L
            mass_matrix = np.zeros((12, 12))
            for i in range(3):
                mass_matrix[i, i] = m
                mass_matrix[i + 6, i + 6] = m
            self.mass = mass_matrix
        else:
            # 一致质量矩阵
            m = rho * L / 420.0
            ml = np.zeros((12, 12))

            # 轴向质量
            ml[0, 0] = ml[6, 6] = m * 140.0
            ml[0, 6] = ml[6, 0] = m * 70.0

            # 扭转质量
            ml[3, 3] = ml[9, 9] = m * (Jx / A) * 140.0
            ml[3, 9] = ml[9, 3] = m * (Jx / A) * 70.0

            # z方向弯曲质量
            ml[2, 2] = ml[8, 8] = m * 156.0
            ml[2, 8] = ml[8, 2] = m * 54.0
            ml[4, 4] = ml[10, 10] = m * 4.0 * L * L
            ml[4, 10] = ml[10, 4] = -m * 3.0 * L * L
            ml[2, 4] = ml[4, 2] = -m * 22.0 * L
            ml[8, 10] = ml[10, 8] = -ml[2, 4]
            ml[2, 10] = ml[10, 2] = m * 13.0 * L
            ml[4, 8] = ml[8, 4] = -ml[2, 10]

            # y方向弯曲质量
            ml[1, 1] = ml[7, 7] = m * 156.0
            ml[1, 7] = ml[7, 1] = m * 54.0
            ml[5, 5] = ml[11, 11] = m * 4.0 * L * L
            ml[5, 11] = ml[11, 5] = -m * 3.0 * L * L
            ml[1, 5] = ml[5, 1] = m * 22.0 * L
            ml[7, 11] = ml[11, 7] = -ml[1, 5]
            ml[1, 11] = ml[11, 1] = -m * 13.0 * L
            ml[5, 7] = ml[7, 5] = -ml[1, 11]

            self.mass = ml

        return self.mass

    def CalculateElementStress(self, displacement):
        """计算单元内力 (基本内力)"""
        return np.zeros(2), np.zeros(2), np.zeros(2), np.zeros(2), np.zeros(2), np.zeros(2)
        E = self.cha_dict[MaterialKey.E]
        G = self.cha_dict[MaterialKey.G]
        A = self.cha_dict[SectionKey.Area]
        Jx = self.cha_dict[SectionKey.Tor]
        Iy = self.cha_dict[SectionKey.It]
        Iz = self.cha_dict[SectionKey.Is]
        L = self.length

        # 提取基本位移
        v = np.array([
            displacement[6] - displacement[0],  # 轴向变形
            displacement[5],  # i端绕z轴转角
            displacement[11],  # j端绕z轴转角
            displacement[4],  # i端绕y轴转角
            displacement[10],  # j端绕y轴转角
            displacement[9] - displacement[3]  # 扭转角
        ])

        oneOverL = 1.0 / L
        EoverL = E * oneOverL
        q = np.zeros(6)

        # 轴向力
        q[0] = A * EoverL * v[0]

        # 扭矩
        q[5] = G * Jx * oneOverL * v[5]

        # z方向弯曲 (x-y平面)
        if self.releasez == 0:  # 固接
            EIzoverL4 = 4.0 * E * Iz / L
            EIzoverL2 = 2.0 * E * Iz / L
            q[1] = EIzoverL4 * v[1] + EIzoverL2 * v[2]
            q[2] = EIzoverL2 * v[1] + EIzoverL4 * v[2]
        elif self.releasez == 1:  # 释放i端
            q[1] = 0.0
            q[2] = 3.0 * E * Iz / L * v[2]
        elif self.releasez == 2:  # 释放j端
            q[1] = 3.0 * E * Iz / L * v[1]
            q[2] = 0.0
        elif self.releasez == 3:  # 两端释放
            q[1] = 0.0
            q[2] = 0.0

        # y方向弯曲 (x-z平面)
        if self.releasey == 0:  # 固接
            EIyoverL4 = 4.0 * E * Iy / L
            EIyoverL2 = 2.0 * E * Iy / L
            q[3] = EIyoverL4 * v[3] + EIyoverL2 * v[4]
            q[4] = EIyoverL2 * v[3] + EIyoverL4 * v[4]
        elif self.releasey == 1:  # 释放i端
            q[3] = 0.0
            q[4] = 3.0 * E * Iy / L * v[4]
        elif self.releasey == 2:  # 释放j端
            q[3] = 3.0 * E * Iy / L * v[3]
            q[4] = 0.0
        elif self.releasey == 3:  # 两端释放
            q[3] = 0.0
            q[4] = 0.0

        self.stress = q
        return self.stress

    def calculate_section_forces(self, displacement, x):
        """
        计算指定位置截面内力
        :param displacement: 单元位移向量 (12维)
        :param x: 沿梁长度位置 (0 <= x <= L)
        :return: 截面内力向量 (6维: N, Vy, Vz, T, My, Mz)
        """
        # 获取基本内力
        q = self.calculate_stress(displacement)
        L = self.length

        # 位置比例
        xL = x / L
        s = np.zeros(6)

        # 轴力 (沿长度不变)
        s[0] = q[0]

        # 扭矩 (沿长度不变)
        s[3] = q[5]

        # z方向弯矩和剪力
        if self.releasez == 0:  # 固接
            s[5] = q[1] * (1 - xL) + q[2] * xL  # Mz
            s[1] = (q[1] + q[2]) / L  # Vy
        elif self.releasez == 1:  # 释放i端
            s[5] = q[2] * xL  # Mz
            s[1] = q[2] / L  # Vy
        elif self.releasez == 2:  # 释放j端
            s[5] = q[1] * (1 - xL)  # Mz
            s[1] = q[1] / L  # Vy
        else:  # 两端释放
            s[5] = 0.0
            s[1] = 0.0

        # y方向弯矩和剪力
        if self.releasey == 0:  # 固接
            s[4] = q[3] * (1 - xL) + q[4] * xL  # My
            s[2] = (q[3] + q[4]) / L  # Vz
        elif self.releasey == 1:  # 释放i端
            s[4] = q[4] * xL  # My
            s[2] = q[4] / L  # Vz
        elif self.releasey == 2:  # 释放j端
            s[4] = q[3] * (1 - xL)  # My
            s[2] = q[3] / L  # Vz
        else:  # 两端释放
            s[4] = 0.0
            s[2] = 0.0

        return s

    def CalculateBasic(self):
        pass

    def ReCalculateElementStiffness(self):
        pass

    def CalElementDMatrix(self, an_type=None):
        pass


if __name__ == "__main__":
    # ele = Beam188(-1)
    # ele.cha_dict = {MaterialKey.E: 1, PropertyKey.ThicknessOrArea: np.sqrt(3)}
    # ele.node_coords = np.array([[0, 0, 0],
    #                             [1, 1, 1]], dtype=float)
    # print(ele.ElementStiffness())
    # mlogger.debug("finish")
    # print(BeamCalculator.CalculateMomentOfInertiaOfArea(BeamSectionType.Rectangle, (10, 12)))
    # print(BeamCalculator.CalculateMomentOfInertiaOfArea(BeamSectionType.CircleSolid, (10,)))
    # print(BeamCalculator.CalculateMomentOfInertiaOfArea(BeamSectionType.I, (6,6,8,0.3,0.5,0.5)))
    aa = BeamCalculator.CalEffectiveShearArea(BeamSectionType.I, (6.0, 6.0, 8.0, 0.5, 0.5, 0.5))
