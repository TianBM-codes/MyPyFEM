#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys

from element.ElementBase import *
import numpy as np
from abc import ABC


class MITC4(ElementBaseClass, ABC):
    """
    MITC4 Element class
    Reference:
    1.《有限元法、理论、格式与求解方法》上册Bathe P395
    2. 王欢
    """
    def __init__(self, eid=None):
        super().__init__(eid)
        self.nodes_count = 4  # Each element has 8 nodes
        self.K = np.mat(np.zeros([8, 8], dtype=float))  # 刚度矩阵
        self.vtp_type = "quad"
        self.thickness = None

    def CalElementDMatrix(self, an_type=None):
        """
        TODO: 现在只能处理平面应力, 平面应变如何？
        计算本构矩阵, 弹性模量和泊松比, Bathe 上册P184
        """
        e = self.ele_mat_dict[MaterialKey.E]
        niu = self.ele_mat_dict[MaterialKey.Niu]
        if an_type == MaterialMatrixType.PlaneStree or an_type is None:
            a = e / (1 - niu ** 2)
            self.D = a * np.mat(np.array([[1, niu, 0],
                                          [niu, 1, 0],
                                          [0, 0, 0.5 * (1 - niu)]]), dtype=float)
        elif an_type == MaterialMatrixType.PlaneStrain:
            a = e * (1 - niu) / (1 + niu) / (1 - 2 * niu)
            self.D = a * np.mat(np.array([[1, niu / (1 - niu), 0],
                                          [niu(1 - niu), 1, 0],
                                          [0, 0, 0.5 * (1 - 2 * niu) / (1 - niu)]]), dtype=float)
        else:
            mlogger.fatal("Unknown an_dimension")
            sys.exit(1)

    def ElementStiffness(self):
        """
        TODO: Wilson协调元, 王勖成P211, 剪切锁死
        Bathe 上册 P323
        dimension: 8*3, [[x1,y1,z1],[x2,y2,z2],...[x8,y8,z8]], type:np.mat, dtype:float

        # Shape Function:
        N1 = 0.25 * (1 + r) * (1 + s)
        N2 = 0.25 * (1 - r) * (1 + s)
        N3 = 0.25 * (1 - r) * (1 - s)
        N4 = 0.25 * (1 + r) * (1 - s)

        # Partial
        dN1dr, dN1ds =  0.25 * (1 + s),  0.25 * (1 + r)
        dN2dr, dN2ds = -0.25 * (1 + s),  0.25 * (1 - r)
        dN3dr, dN3ds =  0.25 * (s - 1),  0.25 * (r - 1)
        dN4dr, dN4ds =  0.25 * (1 - s), -0.25 * (1 + r)
        """
        assert self.node_coords.shape == (4,2)

        # Gaussian Weight
        sample_pt, weight = GaussIntegrationPoint.GetSamplePointAndWeight(2)

        # 在4个高斯点上积分
        for ri in range(2):
            for si in range(2):
                r, s = sample_pt[ri], sample_pt[si]
                dNdr = np.mat(np.array([[0.25 * (1 + s), -0.25 * (1 + s), 0.25 * (s - 1), 0.25 * (1 - s)],
                                        [0.25 * (1 + r), 0.25 * (1 - r), 0.25 * (r - 1), -0.25 * (1 + r)]]), dtype=float)
                g_weight = weight[ri] * weight[si]

                # Jacobi 2*2 & B Matrix 3*8
                J = np.matmul(dNdr, self.node_coords)
                det_J = np.linalg.det(J)
                J_inv = np.linalg.inv(J)
                B_pre = np.matmul(J_inv, dNdr)
                B = np.mat(np.array([[B_pre[0, 0], 0, B_pre[0, 1], 0, B_pre[0, 2], 0, B_pre[0, 3], 0],
                                     [0, B_pre[1, 0], 0, B_pre[1, 1], 0, B_pre[1, 2], 0, B_pre[1, 3]],
                                     [B_pre[1, 0], B_pre[0, 0], B_pre[1, 1], B_pre[0, 1], B_pre[1, 2], B_pre[0, 2], B_pre[1, 3], B_pre[0, 3]]]), dtype=float)

                self.K = self.K + g_weight * B.T * self.D * B * det_J * self.ele_prop_dict[PropertyKey.ThicknessOrArea]

        return self.K

    def ElementStress(self, displacement):
        """
        Calculate element stress
        """


class CPS3(ElementBaseClass, ABC):
    """ plane2D 3node Element class """

    def __init__(self, eid=None):
        super().__init__(eid)
        self.nodes_count = 3  # Each element has 3 nodes
        self.K = np.mat(np.zeros([6, 6], dtype=float))  # 刚度矩阵
        self.vtp_type = "triangle"
        self.thickness = None

    def CalElementDMatrix(self, an_type=None):
        """
        计算本构矩阵, 弹性模量和泊松比, Bathe 上册P184
        """
        e = self.ele_mat_dict[MaterialKey.E]
        niu = self.ele_mat_dict[MaterialKey.Niu]
        if an_type == MaterialMatrixType.PlaneStree or an_type is None:
            a = e / (1 - niu ** 2)
            self.D = a * np.mat(np.array([[1, niu, 0],
                                          [niu, 1, 0],
                                          [0, 0, 0.5 * (1 - niu)]]), dtype=float)
        elif an_type == MaterialMatrixType.PlaneStrain:
            a = e * (1 - niu) / (1 + niu) / (1 - 2 * niu)
            self.D = a * np.mat(np.array([[1, niu / (1 - niu), 0],
                                          [niu(1 - niu), 1, 0],
                                          [0, 0, 0.5 * (1 - 2 * niu) / (1 - niu)]]), dtype=float)
        else:
            mlogger.fatal("Unknown an_dimension")
            sys.exit(1)

    def ElementStiffness(self):
        """
        TODO: 积分过程是否正确?
        Bathe 上册 P349, 转化到参数坐标下的面积积分后, 在积分域内为常数, 所以积分等于面积 0.5
        dimension: 2*2, [[x1,y1],[x2,y2]], type:np.mat, dtype:float

        # Shape Function:
        N1 = 1 - r - s
        N2 = r
        N2 = s

        # Partial
        dN1dr, dN1ds = -1, -1
        dN2dr, dN2ds =  1,  0
        dN3dr, dN3ds =  0,  1
        """
        assert self.node_coords.shape == (3,2)

        dNdr = np.mat(np.array([[-1, 1, 0],
                                [-1, 0, 1]]), dtype=float)

        # Jacobi 2*2 & B Matrix 3*8
        J = np.matmul(dNdr, self.node_coords)
        det_J = np.linalg.det(J)
        J_inv = np.linalg.inv(J)
        B_pre = np.matmul(J_inv, dNdr)
        B = np.mat(np.array([[B_pre[0, 0], 0, B_pre[0, 1], 0, B_pre[0, 2], 0],
                             [0, B_pre[1, 0], 0, B_pre[1, 1], 0, B_pre[1, 2]],
                             [B_pre[1, 0], B_pre[0, 0], B_pre[1, 1], B_pre[0, 1], B_pre[1, 2], B_pre[0, 2]]]), dtype=float)

        return B.T * self.D * B * det_J * 0.5 * self.ele_prop_dict[PropertyKey.ThicknessOrArea]

    def ElementStress(self, displacement):
        """
        Calculate element stress
        """


if __name__ == "__main__":
    ele = CPS3(-1)
    ele.ele_mat_dict = {MaterialKey.Niu: 0.3, MaterialKey.E: 2e9}
    ele.node_coords = np.mat(np.array([[0, 0],
                                       [4, 0],
                                       [1, 3]]), dtype=float)
    ele.CalElementDMatrix(MaterialMatrixType.PlaneStree)
    ele.ElementStiffness()
    mlogger.debug("finish")
