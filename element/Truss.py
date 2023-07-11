#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from element.ElementBase import *
import numpy as np
from abc import ABC


class T3D2(ElementBaseClass, ABC):
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
        桁架单元无需计算D阵
        """
        pass

    def ElementStiffness(self):
        """
        TODO: 刚度阵中加eps的这种处理方式是否合理？
        TODO: 值不对的问题可以写一个二维杆端元, 然后与曾攀《有限元分析及应用》P123页对比结果
        TODO: 对载荷也乘以全局坐标和局部坐标的转换矩阵, 在程序中如何体现
        如果不在主元上施加小量, 那么对于一些情况会刚度阵奇异, 而在对角线上添加小量, 相当于加一个刚度很小的弹簧,
        横向的位移就为0, 连接其他单元也没关系, 一个典型的情况就是杆方向为x轴方向, 那么只有这一根杆, 则刚度矩
        阵奇异, 典型计算题是/truss/trussbridge/truss_static.inp
        Reference:
        1. 王勖成P302
        2. 曾攀《有限元分析及应用_曾攀》 P110
        """
        assert self.node_coords.shape == (2, 3)

        # 单元参数
        delta = np.asarray(np.diff(self.node_coords, axis=0))[0]
        E = self.cha_dict[MaterialKey.E]
        A = self.cha_dict[PropertyKey.ThicknessOrArea]
        L = np.sqrt(np.dot(delta.T, delta))

        # 几何关系
        Cx, Cy, Cz = delta / L
        Cx2, Cy2, Cz2 = Cx ** 2, Cy ** 2, Cz ** 2
        CxCy, CxCz, CyCz = Cx * Cy, Cx * Cz, Cy * Cz

        eps = 1e-7
        return E * A / L * np.array([[Cx2 + eps, CxCy, CxCz, -Cx2, -CxCy, -CxCz],
                                     [CxCy, Cy2 + eps, CyCz, -CxCy, -Cy2, -CyCz],
                                     [CxCz, CyCz, Cz2 + eps, -CxCz, -CyCz, -Cz2],
                                     [-Cx2, -CxCy, -CxCz, Cx2 + eps, CxCy, CxCz],
                                     [-CxCy, -Cy2, -CyCz, CxCy, Cy2 + eps, CyCz],
                                     [-CxCz, -CyCz, -Cz2, CxCz, CyCz, Cz2 + eps]])

    def ElementStress(self, displacement):
        """
        Calculate element stress
        """
        # fem_database = Domain()
        # x1 = fem_database.GetDisplacement(self.search_node_id[0])
        # x2 = fem_database.GetDisplacement(self.search_node_id[1])
        # self.stress = self.e / self.rod_length * (np.dot(np.asarray(x2), self.cos_angel) - np.dot(np.asarray(x1), self.cos_angel))


if __name__ == "__main__":
    t_ele = T3D2(-1)
    t_ele.ele_mat_dict = {MaterialKey.E: 1, PropertyKey.ThicknessOrArea: np.sqrt(3)}
    t_ele.node_coords = np.array([[0, 0, 0],
                                [1, 1, 1]], dtype=float)
    print(t_ele.ElementStiffness())
    mlogger.debug("finish")
