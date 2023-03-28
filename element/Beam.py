#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from element.ElementBase import *
import numpy as np
from abc import ABC

class BeamCalculate:
    """
    计算梁单元截面属性
    """
    @staticmethod
    def CalculateRotationalInertia(sec_type, sec_data):
        """
        计算梁横截面的转动惯量
        """
        if sec_type == BeamSectionType.GongZiGang:
            pass

    @staticmethod
    def CalculateArea(sec_type, sec_data):
        """
        计算截面的面积
        """
        if sec_type == BeamSectionType.GongZiGang:
            pass

class B3D2(ElementBaseClass, ABC):
    """
    Truss Element class
    """

    def __init__(self, eid):
        super().__init__(eid)
        self.nodes_count = 2  # Each element has 2 nodes
        self._vtp_type = "line"
        self.stiffness = None
        self.stress = None

    def CalElementDMatrix(self, an_type=None):
        """
        TODO: 现在只能处理平面应力, 平面应变如何？
        计算本构矩阵, 弹性模量和泊松比, Bathe 上册P184
        """
        pass

    def ElementStiffness(self):
        """
        Reference:
        1. 王勖成P302
        2. 曾攀《有限元分析及应用_曾攀》 P110
        """
        assert self.node_coords.shape == (2, 3)

        # 单元参数
        delta = np.asarray(np.diff(self.node_coords, axis=0))[0]
        E = self.ele_mat_dict[MaterialKey.E]
        A = self.ele_mat_dict[PropertyKey.ThicknessOrArea]
        L = np.sqrt(np.dot(delta.T, delta))

        # 几何关系
        Cx, Cy, Cz = delta / L
        Cx2, Cy2, Cz2 = Cx ** 2, Cy ** 2, Cz ** 2
        CxCy, CxCz, CyCz = Cx * Cy, Cx * Cz, Cy * Cz

        return E * A / L * np.mat(np.array([[Cx2, CxCy, CxCz, -Cx2, -CxCy, -CxCz],
                                            [CxCy, Cy2, CyCz, -CxCy, -Cy2, -CyCz],
                                            [CxCz, CyCz, Cz2, -CxCz, -CyCz, -Cz2],
                                            [-Cx2, -CxCy, -CxCz, Cx2, CxCy, CxCz],
                                            [-CxCy, -Cy2, -CyCz, CxCy, Cy2, CyCz],
                                            [-CxCz, -CyCz, -Cz2, CxCz, CyCz, Cz2]]))

    def ElementStress(self, displacement):
        """
        Calculate element stress
        """
        # fem_database = Domain()
        # x1 = fem_database.GetDisplacement(self.search_node_id[0])
        # x2 = fem_database.GetDisplacement(self.search_node_id[1])
        # self.stress = self.e / self.rod_length * (np.dot(np.asarray(x2), self.cos_angel) - np.dot(np.asarray(x1), self.cos_angel))


if __name__ == "__main__":
    ele = B3D2(-1)
    ele.ele_mat_dict = {MaterialKey.E:1, PropertyKey.ThicknessOrArea:np.sqrt(3)}
    ele.node_coords = np.mat(np.array([[0,0,0],
                                       [1,1,1]],dtype=float))
    print(ele.ElementStiffness())
    mlogger.debug("finish")
