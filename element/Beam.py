#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys

from element.ElementBase import *
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
    def CalculateMomentOfInertiaOfArea(sec_type: str, sec_data: list[float]) -> dict:
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
            assert len(sec_data) == 1
            I = 0.25 * np.pi * pow(sec_data[0], 4)
            Torsional = I * 2
            return {SectionKey.It: I, SectionKey.Is: I, SectionKey.Tor: Torsional}

        else:
            mlogger.fatal("UnSupport Beam Section Type:{}".format(sec_type))
            sys.exit(1)

    @staticmethod
    def CalEffectiveShearArea(sec_type: BeamSectionType, sec_data: tuple) -> dict:
        """
        计算截面的面积属性, 包括面积、两个方向的抗剪等效面积, 圆的输入是半径
        """
        if sec_type == BeamSectionType.Rectangle:
            assert len(sec_data) == 2
            Area = sec_data[0] * sec_data[1]
            E_Area = 5 / 6 * sec_data[0] * sec_data[1]
            return {SectionKey.Area: Area, SectionKey.At: E_Area, SectionKey.As: E_Area}

        elif sec_type == BeamSectionType.CircleSolid:
            assert len(sec_data) == 1
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


class Beam188(ElementBaseClass, ABC):
    """
    一般直梁, 两节点混合插值单元, 横向剪应变被假设为常量(一阶高斯点处的剪应变)
    Reference:
    1. 《有限元法 理论、格式与求解方法》上 Bathe P382
    """

    def __init__(self, eid):
        super().__init__(eid)
        self.nodes_count = 2  # 每个单元包含2个节点, 表示方向的辅助节点不计算在内
        self._vtp_type = "line"
        self.stiffness = None
        self.stress = None

        # 截面相关, 面积、有效面积
        self.It, self.Is, self.Tor = None, None, None
        self.A, self.effect_as, self.effect_at = None, None, None

    def CalElementDMatrix(self, an_type=None):
        pass

    def ElementStiffness(self):
        """
        TODO: 有整理的pdf
        Reference:
        """
        # 节点顺序为起始点, 终点, 方向点  TODO: 如何求解默认的, 也就是没定义梁方向的截面法向
        assert self.node_coords.shape == (3, 3)

        # 单元参数
        delta = np.asarray(self.node_coords[1, :] - self.node_coords[0, :])
        E = self.cha_dict[MaterialKey.E]
        A = self.cha_dict[SectionKey.Area]
        At = self.cha_dict[SectionKey.At]
        As = self.cha_dict[SectionKey.As]
        G = E / (1 + 2 * self.cha_dict[MaterialKey.Niu])
        L = np.sqrt(np.dot(delta, delta.T))
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

        # 笛卡尔坐标变换为自然坐标, 梁主轴方向r, 梁方向点方向s, 王勖成P331
        normal_direct = np.asarray(self.node_coords[2, :] - (self.node_coords[0, :] + self.node_coords[1, :]) / 2)
        cross_direct = np.cross(delta, normal_direct)
        lxx, lxy, lxz = (delta / L)[0, :]
        lyx, lyy, lyz = (normal_direct / np.linalg.norm(normal_direct))[0, :]
        lzx, lzy, lzz = (cross_direct / np.linalg.norm(cross_direct))[0, :]
        trans_mat = np.asarray([[lxx, lxy, lxz, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                [lyx, lyy, lyz, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                [lzx, lzy, lzz, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, lxx, lxy, lxz, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, lyx, lyy, lyz, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, lzx, lzy, lzz, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, lxx, lxy, lxz, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, lyx, lyy, lyz, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, lzx, lzy, lzz, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, lxx, lxy, lxz],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, lyx, lyy, lyz],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, lzx, lzy, lzz]], dtype=float)

        return np.matmul(np.matmul(trans_mat.T, K), trans_mat)

    def ElementStress(self, displacement):
        """
        Calculate element stress
        """
        # fem_database = Domain()
        # x1 = fem_database.GetDisplacement(self.search_node_id[0])
        # x2 = fem_database.GetDisplacement(self.search_node_id[1])
        # self.stress = self.e / self.rod_length * (np.dot(np.asarray(x2), self.cos_angel) - np.dot(np.asarray(x1), self.cos_angel))


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
        self._vtp_type = "line3"
        self.stiffness = None
        self.stress = None
        self.I = None  # 惯性矩
        self.sec_type = None  # 截面类型
        self.sec_data = None  # 截面参数

    def CalElementDMatrix(self, an_type=None):
        pass

    def ElementStiffness(self):
        """
        TODO: 有整理的pdf
        Reference:
        """
        # 4个节点, 起始点, 终点, 中间点, 方向点
        assert self.node_coords.shape == (4, 3)

        # 单元参数
        delta = np.asarray(np.diff(self.node_coords, axis=0))[0]
        E = self.cha_dict[MaterialKey.E]
        A = self.cha_dict[PropertyKey.ThicknessOrArea]
        G = self.cha_dict[MaterialKey.G]
        L = np.sqrt(np.dot(delta.T, delta))
        I = BeamCalculator.CalculateMomentOfInertiaOfArea(self.sec_type, self.sec_data)
        k = BeamCalculator.GetEquivalentCoff(self.sec_type)

        # 计算弯曲和剪切的刚度部分, 并组装成矩阵
        spt, weight = GaussIntegrationPoint.GetSamplePointAndWeight(2)
        K = np.np.zeros((6, 6), dtype=float)
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

    def ElementStress(self, displacement):
        """
        Calculate element stress
        """
        # fem_database = Domain()
        # x1 = fem_database.GetDisplacement(self.search_node_id[0])
        # x2 = fem_database.GetDisplacement(self.search_node_id[1])
        # self.stress = self.e / self.rod_length * (np.dot(np.asarray(x2), self.cos_angel) - np.dot(np.asarray(x1), self.cos_angel))


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
