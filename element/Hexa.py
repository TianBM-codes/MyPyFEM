#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from element.ElementBase import *
import numpy as np
from abc import ABC


class C3D8(ElementBaseClass, ABC):
    """ hexa Element class """

    def __init__(self, eid=None):
        super().__init__(eid)
        self.nodes_count = 8  # Each element has 8 nodes
        self.K = np.zeros([24, 24], dtype=float)  # 刚度矩阵
        self.vtp_type = "hexahedron"

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
        TODO: https://www.bilibili.com/video/BV19y4y1z76E/?vd_source=f964a6ab226be6b0cd5d082ed4135949 C3D20 还有二维单元的
        Bathe 上册 P323
        dimension: 8*3, [[x1,y1,z1],[x2,y2,z2],...[x8,y8,z8]], type:np.ndarray, dtype:float

        # Shape Function:
        N1 = (1 - r) * (1 - s) * (1 + t) / 8
        N2 = (1 - r) * (1 - s) * (1 - t) / 8
        N3 = (1 - r) * (1 + s) * (1 - t) / 8
        N4 = (1 - r) * (1 + s) * (1 + t) / 8
        N5 = (1 + r) * (1 - s) * (1 + t) / 8
        N6 = (1 + r) * (1 - s) * (1 - t) / 8
        N7 = (1 + r) * (1 + s) * (1 - t) / 8
        N8 = (1 + r) * (1 + s) * (1 + t) / 8

        # Partial
        dN1dr, dN1ds, dN1dt = (s - 1) * (1 + t) / 8, (r - 1) * (1 + t) / 8, (1 - r) * (1 - s) / 8
        dN2dr, dN2ds, dN2dt = (s - 1) * (1 - t) / 8, (r - 1) * (1 - t) / 8, (r - 1) * (1 - s) / 8
        dN3dr, dN3ds, dN3dt = (s + 1) * (t - 1) / 8, (1 - r) * (1 - t) / 8, (r - 1) * (1 + s) / 8
        dN4dr, dN4ds, dN4dt = (s + 1) * (-1 - t) / 8, (1 - r) * (1 + t) / 8, (1 - r) * (1 + s) / 8
        dN5dr, dN5ds, dN5dt = (1 - s) * (1 + t) / 8, -(1 + r) * (1 + t) / 8, (1 + r) * (1 - s) / 8
        dN6dr, dN6ds, dN6dt = (1 - s) * (1 - t) / 8, (1 + r) * (t - 1) / 8, (1 + r) * (s - 1) / 8
        dN7dr, dN7ds, dN7dt = (1 + s) * (1 - t) / 8, (1 + r) * (1 - t) / 8, -(1 + r) * (1 + s) / 8
        dN8dr, dN8ds, dN8dt = (1 + s) * (1 + t) / 8, (1 + r) * (1 + t) / 8, (1 + r) * (1 + s) / 8
        """
        assert len(self.node_coords) == 8

        # Gaussian Weight
        sample_pt, weight = GaussIntegrationPoint.GetSamplePointAndWeight(2)

        # 在8个高斯点上积分
        for ri in range(2):
            for si in range(2):
                for ti in range(2):
                    r, s, t = sample_pt[ri], sample_pt[si], sample_pt[ti]
                    dNdr = np.array([(s + 1) * (-1 - t) / 8, (1 - r) * (1 + t) / 8, (1 - r) * (1 + s) / 8],
                                    [(s - 1) * (1 + t) / 8, (r - 1) * (1 + t) / 8, (1 - r) * (1 - s) / 8],
                                    [(s - 1) * (1 - t) / 8, (r - 1) * (1 - t) / 8, (r - 1) * (1 - s) / 8],
                                    [(s + 1) * (t - 1) / 8, (1 - r) * (1 - t) / 8, (r - 1) * (1 + s) / 8],
                                    [(1 + s) * (1 + t) / 8, (1 + r) * (1 + t) / 8, (1 + r) * (1 + s) / 8],
                                    [(1 - s) * (1 + t) / 8, -(1 + r) * (1 + t) / 8, (1 + r) * (1 - s) / 8],
                                    [(1 - s) * (1 - t) / 8, (1 + r) * (t - 1) / 8, (1 + r) * (s - 1) / 8],
                                    [(1 + s) * (1 - t) / 8, (1 + r) * (1 - t) / 8, -(1 + r) * (1 + s) / 8]).T
                    g_weight = weight[ri] * weight[si] * weight[ti]

                    # Jacobi 3*3 & B Matrix 8*24
                    J = np.matmul(dNdr, self.node_coords)
                    det_J = np.linalg.det(J)
                    J_inv = np.linalg.inv(J)
                    B_pre = np.matmul(J_inv, dNdr)
                    B = np.array([[B_pre[0, 0], 0, 0, B_pre[0, 1], 0, 0, B_pre[0, 2], 0, 0, B_pre[0, 3], 0, 0, B_pre[0, 4], 0, 0, B_pre[0, 5], 0, 0, B_pre[0, 6], 0, 0, B_pre[0, 7], 0, 0],
                                  [0, B_pre[1, 0], 0, 0, B_pre[1, 1], 0, 0, B_pre[1, 2], 0, 0, B_pre[1, 3], 0, 0, B_pre[1, 4], 0, 0, B_pre[1, 5], 0, 0, B_pre[1, 6], 0, 0, B_pre[1, 7], 0],
                                  [0, 0, B_pre[2, 0], 0, 0, B_pre[2, 1], 0, 0, B_pre[2, 2], 0, 0, B_pre[2, 3], 0, 0, B_pre[2, 4], 0, 0, B_pre[2, 5], 0, 0, B_pre[2, 6], 0, 0, B_pre[2, 7]],
                                  [B_pre[1, 0], B_pre[0, 0], 0, B_pre[1, 1], B_pre[0, 1], 0, B_pre[1, 2], B_pre[0, 2], 0, B_pre[1, 3], B_pre[0, 3], 0, B_pre[1, 4], B_pre[0, 4], 0, B_pre[1, 5], B_pre[0, 5], 0, B_pre[1, 6], B_pre[0, 6], 0, B_pre[1, 7], B_pre[0, 7], 0],
                                  [0, B_pre[2, 0], B_pre[1, 0], 0, B_pre[2, 1], B_pre[1, 1], 0, B_pre[2, 2], B_pre[1, 2], 0, B_pre[2, 3], B_pre[1, 3], 0, B_pre[2, 4], B_pre[1, 4], 0, B_pre[2, 5], B_pre[1, 5], 0, B_pre[2, 6], B_pre[1, 6], 0, B_pre[2, 7], B_pre[1, 7]],
                                  [B_pre[2, 0], 0, B_pre[0, 0], B_pre[2, 1], 0, B_pre[0, 1], B_pre[2, 2], 0, B_pre[0, 2], B_pre[2, 3], 0, B_pre[0, 3], B_pre[2, 4], 0, B_pre[0, 4], B_pre[2, 5], 0, B_pre[0, 5], B_pre[2, 6], 0, B_pre[0, 6], B_pre[2, 7], 0, B_pre[0, 7]]])

                    self.K = self.K + g_weight * B.T * self.D * B * det_J

        return self.K

    def ElementStress(self, displacement):
        """
        Calculate element stress
        """
