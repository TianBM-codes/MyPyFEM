#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pypardiso
import scipy.sparse.linalg

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
        self.unv_code = 80600
        self.gs_count = 8
        self.Gaussian_B = []  # 高斯积分点处的应变矩阵
        self.n_dofs_elem = 24
        self.n_face_dofs_elem = 12

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

        # data = [1 - niu, niu, niu, niu, 1 - niu, niu, niu, niu, 1 - niu, 0.5 * (1 - 2 * niu), 0.5 * (1 - 2 * niu), 0.5 * (1 - 2 * niu)]
        # rows = [0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 4, 5]
        # cols = [0, 1, 2, 0, 1, 2, 0, 1, 2, 3, 4, 5]
        # self.D = a * sparse.csc_matrix((data, (rows, cols)), shape=(6, 6))

    def ElementStiffness(self):
        """
        TODO: https://www.bilibili.com/video/BV19y4y1z76E/?vd_source=f964a6ab226be6b0cd5d082ed4135949 C3D20 还有二维单元的
        Bathe 上册 P323
        dimension: 8*3, [[x1,y1,z1],[x2,y2,z2],...[x8,y8,z8]], type:np.ndarray, dtype:float
        """
        assert self.node_coords.shape == (8, 3)

        # 在8个高斯点上积分, 这里还有么有再优化的空间?
        # dNdrs, weights = AllEleTypeDNDrAtGaussianPoint.C3D8
        dNdrs = []
        weights = []
        for ii in range(self.gs_count):
            J = np.matmul(dNdrs[ii], self.node_coords)
            det_J = np.linalg.det(J)
            # J_inv = np.linalg.inv(J)
            # B_pre = np.matmul(J_inv, dNdrs[ii])
            # B_pre = pypardiso.spsolve(sparse.csc_matrix(J), dNdrs[ii])
            B_pre = np.linalg.solve(J, dNdrs[ii])

            B_at_gs_pt = np.asarray([[B_pre[0, 0], 0, 0, B_pre[0, 1], 0, 0, B_pre[0, 2], 0, 0, B_pre[0, 3], 0, 0, B_pre[0, 4], 0, 0, B_pre[0, 5], 0, 0, B_pre[0, 6], 0, 0, B_pre[0, 7], 0, 0],
                                     [0, B_pre[1, 0], 0, 0, B_pre[1, 1], 0, 0, B_pre[1, 2], 0, 0, B_pre[1, 3], 0, 0, B_pre[1, 4], 0, 0, B_pre[1, 5], 0, 0, B_pre[1, 6], 0, 0, B_pre[1, 7], 0],
                                     [0, 0, B_pre[2, 0], 0, 0, B_pre[2, 1], 0, 0, B_pre[2, 2], 0, 0, B_pre[2, 3], 0, 0, B_pre[2, 4], 0, 0, B_pre[2, 5], 0, 0, B_pre[2, 6], 0, 0, B_pre[2, 7]],
                                     [B_pre[1, 0], B_pre[0, 0], 0, B_pre[1, 1], B_pre[0, 1], 0, B_pre[1, 2], B_pre[0, 2], 0, B_pre[1, 3], B_pre[0, 3], 0, B_pre[1, 4], B_pre[0, 4], 0, B_pre[1, 5],
                                      B_pre[0, 5], 0, B_pre[1, 6], B_pre[0, 6], 0, B_pre[1, 7], B_pre[0, 7], 0],
                                     [0, B_pre[2, 0], B_pre[1, 0], 0, B_pre[2, 1], B_pre[1, 1], 0, B_pre[2, 2], B_pre[1, 2], 0, B_pre[2, 3], B_pre[1, 3], 0, B_pre[2, 4], B_pre[1, 4], 0, B_pre[2, 5],
                                      B_pre[1, 5], 0, B_pre[2, 6], B_pre[1, 6], 0, B_pre[2, 7], B_pre[1, 7]],
                                     [B_pre[2, 0], 0, B_pre[0, 0], B_pre[2, 1], 0, B_pre[0, 1], B_pre[2, 2], 0, B_pre[0, 2], B_pre[2, 3], 0, B_pre[0, 3], B_pre[2, 4], 0, B_pre[0, 4], B_pre[2, 5], 0,
                                      B_pre[0, 5], B_pre[2, 6], 0, B_pre[0, 6], B_pre[2, 7], 0, B_pre[0, 7]]], dtype=float)

            self.Gaussian_B.append(B_at_gs_pt)
            self.K = self.K + np.matmul(np.matmul(B_at_gs_pt.T, self.D), B_at_gs_pt) * det_J * weights[ii]

        return self.K

    def ElementStress(self, displacement):
        """
        计算节点应力(已知高斯点处的应变矩阵, 以及节点位移)
        Sigma = Gaussian2Global * B * Global2Gaussian * u
        1. 计算积分点位移
        2. 计算积分点应力
        3. 外推至节点
        4. 节点平均(数据库中的方法)

        Reference:
        1. <<有限单元法>> 王勖成 P168-176
        """
        # 计算高斯点的位移(矩阵中都是正数)
        # Node2Gaussian: Nodes Displacement ==> Gaussian Points Displacement
        a = 0.49056261216234404  # 0.125*(1+1/np.sqrt(3))**3
        b = 0.13144585576580214  # 0.125*(1+1/np.sqrt(3))*2/3
        c = 0.03522081090086451  # 0.125*(1-1/np.sqrt(3))*2/3
        d = 0.00943738783765593  # 0.125*(1-1/np.sqrt(3))**3
        Global2Gaussian = np.asarray([[a, b, c, b, b, c, d, c],
                                      [b, a, b, c, c, b, c, d],
                                      [c, b, a, b, d, c, b, c],
                                      [b, c, b, a, c, d, c, b],
                                      [b, c, d, c, a, b, c, b],
                                      [c, b, c, d, b, a, b, c],
                                      [d, c, b, c, c, b, a, b],
                                      [c, d, c, b, b, c, b, a]], dtype=float)

        # 计算高斯点的应力
        gs_stress = []
        gs_dis = np.matmul(Global2Gaussian, displacement.reshape((8, 3))).reshape((24,))
        for ii in range(self.gs_count):
            gs_stress.append(np.matmul(self.D, np.matmul(self.Gaussian_B[ii], gs_dis)))
        gs_stress = np.asarray(gs_stress)

        """
        高斯点应力外推至节点
        """
        a = 2.549038105676658  # 0.25 * (5 + 3 * np.sqrt(3))
        b = -0.68301270189222  # -0.25 * (np.sqrt(3) + 1)
        c = 0.183012701892219  # 0.25 * (np.sqrt(3) - 1)
        d = -0.04903810567666  # 0.25 * (5 - 3 * np.sqrt(3))

        # Gaussian2Node: Gaussian Points Stress ==> Node Stress
        Gaussian2Global = np.asarray([[a, b, c, b, b, c, d, c],
                                      [b, a, b, c, c, b, c, d],
                                      [c, b, a, b, d, c, b, c],
                                      [b, c, b, a, c, d, c, b],
                                      [b, c, d, c, a, b, c, b],
                                      [c, b, c, d, b, a, b, c],
                                      [d, c, b, c, c, b, a, b],
                                      [c, d, c, b, b, c, b, a]], dtype=float)

        node_stress = np.matmul(Gaussian2Global, gs_stress)

        return node_stress


"""
    稀疏矩阵存储应变阵, 发现还没全矩阵计算快
        rows = [0, 0, 0, 0, 0, 0, 0, 0,
                1, 1, 1, 1, 1, 1, 1, 1,
                2, 2, 2, 2, 2, 2, 2, 2,
                3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5]
        cols = [0, 3, 6, 9, 12, 15, 18, 21,
                1, 4, 7, 10, 13, 16, 19, 22,
                2, 5, 8, 11, 14, 17, 20, 23,
                0, 1, 3, 4, 6, 7, 9, 10, 12, 13, 15, 16, 18, 19, 21, 22,
                1, 2, 4, 5, 7, 8, 10, 11, 13, 14, 16, 17, 19, 20, 22, 23,
                0, 2, 3, 5, 6, 8, 9, 11, 12, 14, 15, 17, 18, 20, 21, 23]
        data = list(B_pre.flatten())
        data.extend([B_pre[1, 0], B_pre[0, 0], B_pre[1, 1], B_pre[0, 1], B_pre[1, 2], B_pre[0, 2], B_pre[1, 3], B_pre[0, 3], B_pre[1, 4], B_pre[0, 4], B_pre[1, 5], B_pre[0, 5], B_pre[1, 6], B_pre[0, 6], B_pre[1, 7], B_pre[0, 7]])
        data.extend([B_pre[2, 0], B_pre[1, 0], B_pre[2, 1], B_pre[1, 1], B_pre[2, 2], B_pre[1, 2], B_pre[2, 3], B_pre[1, 3], B_pre[2, 4], B_pre[1, 4], B_pre[2, 5], B_pre[1, 5], B_pre[2, 6], B_pre[1, 6], B_pre[2, 7], B_pre[1, 7]])
        data.extend([B_pre[2, 0], B_pre[0, 0], B_pre[2, 1], B_pre[0, 1], B_pre[2, 2], B_pre[0, 2], B_pre[2, 3], B_pre[0, 3], B_pre[2, 4], B_pre[0, 4], B_pre[2, 5], B_pre[0, 5], B_pre[2, 6], B_pre[0, 6], B_pre[2, 7], B_pre[0, 7]])
        B = sparse.coo_matrix((data, (rows, cols)), shape=(6, 24))
"""
