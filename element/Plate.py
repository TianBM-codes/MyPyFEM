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
        self.K = np.zeros([8, 8], dtype=float)  # 刚度矩阵
        self.vtp_type = "quad"
        self.thickness = None

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
        dimension: 4*3, [[x1,y1,z1],[x2,y2,z2],...[x4,y4,z4]], type:np.ndarray, dtype:float

        # Shape Function:
        N1 = 0.25 * (1 + r) * (1 + s)
        N2 = 0.25 * (1 - r) * (1 + s)
        N3 = 0.25 * (1 - r) * (1 - s)
        N4 = 0.25 * (1 + r) * (1 - s)

        # Partial
        """
        assert self.node_coords.shape == (4, 3)

        # Gaussian Weight
        sample_pt, weight = GaussIntegrationPoint.GetSamplePointAndWeight(2)

        """
        Bending Part
        """

        """
        Shear Stain Part
        """
        tt = np.array([[1, -1, -1, 1], [1, -1, 1, -1], [1, 1, -1, -1]], dtype=float)
        Ax, Bx, Cx = tt * self.node_coords[:, 0]
        Ay, By, Cy = tt * self.node_coords[:, 1]
        r, s = 0, 0

        # Partial Derivative
        pN1pr, pN1ps = 0.25 * (1 + s), 0.25 * (1 + r)
        pN2pr, pN2ps = -0.25 * (1 + s), 0.25 * (1 - r)
        pN3pr, pN3ps = 0.25 * (s - 1), 0.25 * (r - 1)
        pN4pr, pN4ps = 0.25 * (1 - s), -0.25 * (1 + r)
        pXpr = np.matmul(np.array([pN1pr, pN2pr, pN3pr, pN4pr]), self.node_coords[:, 0])
        pYpr = np.matmul(np.array([pN1pr, pN2pr, pN3pr, pN4pr]), self.node_coords[:, 1])
        pXps = np.matmul(np.array([pN1ps, pN2ps, pN3ps, pN4ps]), self.node_coords[:, 0])
        pYps = np.matmul(np.array([pN1ps, pN2ps, pN3ps, pN4ps]), self.node_coords[:, 1])
        alpha = np.arctan(pYpr / pXpr)
        beta = np.arctan(pYps / pXps)
        detJ = pXpr * pYps - pYpr * pXps

        coeff_gama_rz = 0.125 * np.sqrt((Cx + r * Bx) ** 2 + (Cy + r * By) ** 2) / detJ
        coeff_gama_sz = 0.125 * np.sqrt((Ax + s * Bx) ** 2 + (Ay + s * By) ** 2) / detJ

        gama_rz = np.zeros((12, 1), dtype=float)
        # w1, w2, w3, w4
        gama_rz[0, 0] = 0.5 * coeff_gama_rz * (1 + s)
        gama_rz[3, 0] = -gama_rz[0, 0]
        gama_rz[6, 0] = -0.5 * coeff_gama_rz * (1 - s)
        gama_rz[9, 0] = -gama_rz[6, 0]

        # theta_1x, theta_2x, theta_3x, theta_4x
        gama_rz[1, 0] = -gama_rz[0, 0] * 0.5 * (self.node_coords[0, 1] - self.node_coords[1, 1])
        gama_rz[4, 0] = gama_rz[1, 0]
        gama_rz[7, 0] = 0.5 * gama_rz[6, 0] * (self.node_coords[2, 1] - self.node_coords[3, 1])
        gama_rz[10, 0] = gama_rz[7, 0]

        # theta_1y, theta_2y, theta_3y, theta_4y
        gama_rz[2, 0] = gama_rz[0, 0] * 0.5 * (self.node_coords[0, 0] - self.node_coords[1, 0])
        gama_rz[5, 0] = gama_rz[2, 0]
        gama_rz[8, 0] = -gama_rz[6, 0] * 0.5 * (self.node_coords[3, 0] - self.node_coords[2, 0])
        gama_rz[11, 0] = gama_rz[8, 0]

        gama_sz = np.zeros((12, 1), dtype=float)
        # w1, w2, w3, w4
        gama_sz[0, 0] = 0.5 * coeff_gama_sz * (1 + r)
        gama_sz[3, 0] = 0.5 * coeff_gama_sz * (1 - r)
        gama_sz[6, 0] = -gama_sz[3, 0]
        gama_sz[9, 0] = -gama_sz[0, 0]

        # theta_1x, theta_2x, theta_3x, theta_4x
        gama_sz[1, 0] = -gama_sz[0, 0] * 0.5 * (self.node_coords[0, 1] - self.node_coords[3, 1])
        gama_sz[4, 0] = -0.5 * gama_sz[3, 0] * (self.node_coords[1, 1] - self.node_coords[2, 1])
        gama_sz[7, 0] = gama_sz[4, 0]
        gama_sz[10, 0] = gama_sz[1, 0]

        # theta_1y, theta_2y, theta_3y, theta_4y
        gama_sz[2, 0] = gama_sz[0, 0] * 0.5 * (self.node_coords[0, 0] - self.node_coords[3, 0])
        gama_sz[5, 0] = gama_sz[1, 0] * 0.5 * (self.node_coords[1, 0] - self.node_coords[2, 0])
        gama_sz[8, 0] = gama_sz[5, 0]
        gama_sz[11, 0] = gama_sz[2, 0]

        # local coord to cartesian coord
        gama_xz = gama_rz * np.sin(beta) - gama_sz * np.sin(alpha)
        gama_yz = -gama_rz * np.cos(beta) + gama_sz * np.cos(alpha)

        # 在4个高斯点上积分
        for ri in range(2):
            for si in range(2):
                r, s = sample_pt[ri], sample_pt[si]
                g_weight = weight[ri] * weight[si]

                self.K = self.K + g_weight * B.T * self.D * B * det_J * self.cha_dict[PropertyKey.ThicknessOrArea]

        return self.K

    def ElementStress(self, displacement):
        """
        Calculate element stress
        """


class MITC3(ElementBaseClass, ABC):
    """ plate 3node Element class """

    def __init__(self, eid=None):
        super().__init__(eid)
        self.nodes_count = 3  # Each element has 3 nodes
        self.K = np.zeros([6, 6], dtype=float)  # 刚度矩阵
        self.vtp_type = "triangle"
        self.thickness = None

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
        TODO: 积分过程是否正确?
        Bathe 上册 P349, 转化到参数坐标下的面积积分后, 在积分域内为常数, 所以积分等于面积 0.5
        dimension: 2*2, [[x1,y1],[x2,y2]], type:np.ndarray, dtype:float

        # Shape Function:
        N1 = 1 - r - s
        N2 = r
        N2 = s

        # Partial
        dN1dr, dN1ds = -1, -1
        dN2dr, dN2ds =  1,  0
        dN3dr, dN3ds =  0,  1
        """
        assert self.node_coords.shape == (3, 2)

        dNdr = np.array([[-1, 1, 0],
                         [-1, 0, 1]], dtype=float)

        # Jacobi 2*2 & B Matrix 3*8
        J = np.matmul(dNdr, self.node_coords)
        det_J = np.linalg.det(J)
        J_inv = np.linalg.inv(J)
        B_pre = np.matmul(J_inv, dNdr)
        B = np.array([[B_pre[0, 0], 0, B_pre[0, 1], 0, B_pre[0, 2], 0],
                      [0, B_pre[1, 0], 0, B_pre[1, 1], 0, B_pre[1, 2]],
                      [B_pre[1, 0], B_pre[0, 0], B_pre[1, 1], B_pre[0, 1], B_pre[1, 2], B_pre[0, 2]]], dtype=float)

        return B.T * self.D * B * det_J * 0.5 * self.cha_dict[PropertyKey.ThicknessOrArea]

    def ElementStress(self, displacement):
        """
        Calculate element stress
        """


if __name__ == "__main__":
    ele = MITC3(-1)
    ele.cha_dict = {MaterialKey.Niu: 0.3, MaterialKey.E: 2e9}
    ele.node_coords = np.array([[0, 0],
                                [4, 0],
                                [1, 3]], dtype=float)
    ele.CalElementDMatrix(MaterialMatrixType.PlaneStree)
    ele.ElementStiffness()
    mlogger.debug("finish")
