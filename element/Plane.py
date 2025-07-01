#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from abc import ABC

import numpy as np

from element.ElementBase import *


class CPS4(ElementBaseClass, ABC):
    """ CPS4 Element class """

    def __init__(self, eid=None):
        super().__init__(eid)
        self.nodes_count = 4  # Each element has 8 nodes
        self.K = np.zeros([8, 8], dtype=float)  # 刚度矩阵
        self.vtu_type = "quad"
        self.B = []

    def CalElementDMatrix(self, an_type=None):
        """
        TODO: 现在只能处理平面应力, 平面应变如何？
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

    def ElementStiffness(self, from_origin=False):
        """
        TODO: Wilson协调元, 王勖成P211, 剪切锁死
        Bathe 上册 P323
        dimension: 8*3, [[x1,y1,z1],[x2,y2,z2],...[x8,y8,z8]], type:np.ndarray, dtype:float

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
        assert self.node_coords.shape == (4, 2)
        self.CalElementDMatrix()

        # Gaussian Weight
        sample_pt, weight = GaussIntegrationPoint.GetSamplePointAndWeight(2)

        # 在4个高斯点上积分
        for ri in range(2):
            for si in range(2):
                r, s = sample_pt[ri], sample_pt[si]
                dNdr = np.array([[0.25 * (1 + s), -0.25 * (1 + s), 0.25 * (s - 1), 0.25 * (1 - s)],
                                 [0.25 * (1 + r), 0.25 * (1 - r), 0.25 * (r - 1), -0.25 * (1 + r)]], dtype=float)
                g_weight = weight[ri] * weight[si]

                # Jacobi 2*2 & B Matrix 3*8
                J = dNdr @ self.node_coords
                det_J = np.linalg.det(J)
                J_inv = np.linalg.inv(J)
                B_pre = np.matmul(J_inv, dNdr)
                B = np.array([[B_pre[0, 0], 0, B_pre[0, 1], 0, B_pre[0, 2], 0, B_pre[0, 3], 0],
                              [0, B_pre[1, 0], 0, B_pre[1, 1], 0, B_pre[1, 2], 0, B_pre[1, 3]],
                              [B_pre[1, 0], B_pre[0, 0], B_pre[1, 1], B_pre[0, 1], B_pre[1, 2], B_pre[0, 2], B_pre[1, 3], B_pre[0, 3]]], dtype=float)

                self.B.append(B)

                self.K = self.K + g_weight * B.T @ self.D @ B * det_J * self.cha_dict[MaterialKey.Thickness]

        return self.K

    def CalculateElementStress(self, displacement):
        """
        Calculate element stress
        """
        gp_stresses = np.zeros((4, 3))
        for i, B in enumerate(self.B):
            strain = B @ displacement  # ε = [ε_x, ε_y, γ_xy]^T
            gp_stresses[i] = self.D @ strain  # σ = [σ_x, σ_y, τ_xy]^T

        A = np.zeros((4, 4))
        gp_nat_coords = [
            (-0.577350269189626, -0.577350269189626),  # gp0
            (-0.577350269189626, 0.577350269189626),  # gp3
            (0.577350269189626, -0.577350269189626),  # gp1
            (0.577350269189626, 0.577350269189626)  # gp2
        ]

        node_nat_coords = [(-1, -1), (1, -1), (1, 1), (-1, 1)]

        def shape_func(i, r, s):
            ri, si = node_nat_coords[i]
            return 0.25 * (1 + r * ri) * (1 + s * si)

        for i in range(4):  # 节点循环
            for j in range(4):  # 高斯点循环
                r, s = gp_nat_coords[j]
                A[i, j] = shape_func(i, r, s)

        A_inv = np.linalg.inv(A)
        node_stresses = np.zeros((4, 3), dtype=np.float32)
        for stress_idx in range(3):
            gp_stress_component = gp_stresses[:, stress_idx]
            node_stress_component = A_inv @ gp_stress_component
            node_stresses[:, stress_idx] = node_stress_component

        return node_stresses[:, 0], node_stresses[:, 1], np.zeros(4), node_stresses[:, 2], np.zeros(4), np.zeros(4)

    def ElementMass(self):
        pass

    def CalculateBasic(self):
        pass

    def ReCalculateElementStiffness(self):
        pass


class CPS3(ElementBaseClass, ABC):
    """ plane2D 3node Element class """

    def __init__(self, eid=None):
        super().__init__(eid)
        self.nodes_count = 3  # Each element has 3 nodes
        self.K = np.zeros([6, 6], dtype=float)  # 刚度矩阵
        self.vtu_type = "triangle"
        self.B = None

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

    def ElementStiffness(self, from_origin=False):
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
        self.CalElementDMatrix()

        dNdr = np.array([[-1, 1, 0],
                         [-1, 0, 1]], dtype=float)

        # Jacobi 2*2 & B Matrix 3*8
        J = np.matmul(dNdr, self.node_coords)
        det_J = np.linalg.det(J)
        J_inv = np.linalg.inv(J)
        B_pre = np.matmul(J_inv, dNdr)
        self.B = np.array([[B_pre[0, 0], 0, B_pre[0, 1], 0, B_pre[0, 2], 0],
                           [0, B_pre[1, 0], 0, B_pre[1, 1], 0, B_pre[1, 2]],
                           [B_pre[1, 0], B_pre[0, 0], B_pre[1, 1], B_pre[0, 1], B_pre[1, 2], B_pre[0, 2]]], dtype=float)

        return self.B.T @ self.D @ self.B * det_J * 0.5 * self.cha_dict[MaterialKey.Thickness]

    def CalculateElementStress(self, displacement):
        """
        Reference:
        1. 《有限单元法》王勖成 P175
        """
        sigma_xx, sigma_yy, tau_xy = self.B @ displacement
        return np.array([sigma_xx] * 3), np.array([sigma_yy] * 3), np.zeros(3), np.array([tau_xy] * 3), np.zeros(3), np.zeros(3)

    def ElementMass(self):
        pass

    def CalculateBasic(self):
        pass

    def ReCalculateElementStiffness(self):
        pass


if __name__ == "__main__":
    t_ele = CPS3(-1)
    t_ele.cha_dict = {MaterialKey.Niu: 0.3, MaterialKey.E: 2e9}
    t_ele.node_coords = np.array([[0, 0],
                                  [4, 0],
                                  [1, 3]], dtype=float)
    t_ele.CalElementDMatrix(MaterialMatrixType.PlaneStree)
    t_ele.ElementStiffness()
    mlogger.debug("finish")
