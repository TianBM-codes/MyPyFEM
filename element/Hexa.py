#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pypardiso
import scipy.sparse.linalg

from element.ElementBase import *
import numpy as np
from abc import ABC


class stdC3D8(ElementBaseClass, ABC):
    """ hexa Element class """

    def __init__(self, eid=None):
        super().__init__(eid)
        self.nodes_count = 8  # Each element has 8 nodes
        self.K = np.zeros([24, 24], dtype=float)  # 刚度矩阵
        self.vtu_type = "hexahedron"
        self.unv_code = 80600
        self.gs_count = 8
        self.Gaussian_B = []  # 高斯积分点处的应变矩阵
        self.block_size = 576
        self.node_dof_count = 3

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

    def ElementStiffness(self, from_origin=False):
        """
        TODO: https://www.bilibili.com/video/BV19y4y1z76E/?vd_source=f964a6ab226be6b0cd5d082ed4135949 C3D20 还有二维单元的
        Bathe 上册 P323
        dimension: 8*3, [[x1,y1,z1],[x2,y2,z2],...[x8,y8,z8]], type:np.ndarray, dtype:float
        """
        assert self.node_coords.shape == (8, 3)
        self.CalElementDMatrix()

        # 在8个高斯点上积分, 这里还有么有再优化的空间?
        dNdrs, weights = AllEleTypeDNDrAtGaussianPoint.C3D8
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
            self.K = self.K + B_at_gs_pt.T @ self.D @ B_at_gs_pt * det_J * weights[ii]

        return self.K

    def CalculateElementStress(self, displacement):
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

        return node_stress[:, 0], node_stress[:, 1], node_stress[:, 2], node_stress[:, 3], node_stress[:, 4], node_stress[:, 5]

    def ElementMass(self):
        pass

    def CalculateBasic(self):
        pass

    def ReCalculateElementStiffness(self):
        pass


class C3D8(ElementBaseClass, ABC):
    """ hexa Element class """

    def __init__(self, eid=None):
        super().__init__(eid)
        self.nodes_count = 8  # Each element has 8 nodes
        self.K = np.zeros([24, 24], dtype=float)  # 刚度矩阵
        self.vtu_type = "hexahedron"
        self.unv_code = 80600
        self.gs_count = 8
        self.Gaussian_B = []  # 高斯积分点处的应变矩阵
        self.block_size = 576
        self.node_dof_count = 3

        self.ndm = 3  # Spatial dimension
        self.ndf = 3  # Degrees of freedom per node
        self.nstress = 6  # Stress/strain components
        self.nodes_count = 8  # Nodes per element
        self.numberGauss = 8  # Gauss points
        self.nShape = 4  # Shape function components

        # Integration points and weights (2x2x2 Gauss integration)
        self.sg = np.array([-1 / np.sqrt(3), 1 / np.sqrt(3)])  # Reduced integration points
        self.wg = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])  # Uniform weights

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

    def ElementStiffness(self, from_origin=False):
        """
        TODO: https://www.bilibili.com/video/BV19y4y1z76E/?vd_source=f964a6ab226be6b0cd5d082ed4135949 C3D20 还有二维单元的
        Bathe 上册 P323
        dimension: 8*3, [[x1,y1,z1],[x2,y2,z2],...[x8,y8,z8]], type:np.ndarray, dtype:float
        """
        assert self.node_coords.shape == (8, 3)
        self.CalElementDMatrix()

        volume = 0.0
        Shape = np.zeros((self.nShape, self.nodes_count, self.numberGauss))
        dvol = np.zeros(self.numberGauss)
        shpBar = np.zeros((self.nShape, self.nodes_count))
        gaussPoint = np.zeros(3)

        dNdrs, weights = AllEleTypeDNDrAtGaussianPoint.C3D8

        # Gauss loop to compute and save shape functions
        count = 0
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    gaussPoint = np.array([self.sg[i], self.sg[j], self.sg[k]])
                    xsj, shp = self.shp3d(gaussPoint)

                    # Save shape functions
                    for p in range(self.nShape):
                        for q in range(self.nodes_count):
                            Shape[p, q, count] = shp[p, q]

                    # Volume element
                    dvol[count] = self.wg[count] * xsj
                    volume += dvol[count]

                    # Accumulate mean shape functions
                    for p in range(self.nShape):
                        for q in range(self.nodes_count):
                            shpBar[p, q] += dvol[count] * shp[p, q]

                    count += 1

        # Normalize mean shape functions
        for p in range(self.nShape):
            for q in range(self.nodes_count):
                shpBar[p, q] /= volume

        # Second Gauss loop for residual and tangent
        for i in range(self.numberGauss):
            # Retrieve shape functions
            shp = Shape[:, :, i]

            # Compute strain
            strain = np.zeros(self.nstress)
            # for j in range(self.nodes_count):
            #     BJ = self.computeBbar(j, shp, shpBar)
            # ul = self.nodePointers[j].getTrialDisp()  # Assume returns [ux, uy, uz]
            # strain += BJ @ ul

            # Get material response
            # success = self.materialPointers[i].setTrialStrain(strain)
            # stress = self.materialPointers[i].getStress()  # Returns (6,) array
            # stress *= dvol[i]  # Scale by volume

            # Get tangent if needed
            dd = self.D * dvol[i]

            # Residual and tangent calculations
            jj = 0
            for j in range(self.nodes_count):
                BJ = self.computeBbar(j, shp, shpBar)
                BJtran = BJ.T

                # Node residual
                # residJ = BJtran @ stress
                # for p in range(self.ndf):
                #     self.resid[jj + p] += residJ[p]
                #     if not self.applyLoad:
                #         self.resid[jj + p] -= dvol[i] * self.b[p] * shp[3, j]
                #     else:
                #         self.resid[jj + p] -= dvol[i] * self.appliedB[p] * shp[3, j]

                # Tangent matrix
                BJtranD = BJtran @ dd
                kk = 0
                for k in range(self.nodes_count):
                    BK = self.computeBbar(k, shp, shpBar)
                    stiffJK = BJtranD @ BK
                    for p in range(self.ndf):
                        for q in range(self.ndf):
                            self.K[jj + p, kk + q] += stiffJK[p, q]
                    kk += self.ndf

                jj += self.ndf

        return self.K

    def CalculateElementStress(self, displacement):
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

        return node_stress[:, 0], node_stress[:, 1], node_stress[:, 2], node_stress[:, 3], node_stress[:, 4], node_stress[:, 5]

    def ElementMass(self):
        pass

    def CalculateBasic(self):
        pass

    def ReCalculateElementStiffness(self):
        pass

    def shp3d(self, ss):
        # 输入 ss: 自然坐标 [ξ, η, ζ]
        # 返回 xsj: Jacobian 行列式
        # 返回 shp: 形状函数和导数数组 [4x8]

        # 计算辅助变量
        ap1 = 1.0 + ss[0]
        am1 = 1.0 - ss[0]
        ap2 = 1.0 + ss[1]
        am2 = 1.0 - ss[1]
        ap3 = 1.0 + ss[2]
        am3 = 1.0 - ss[2]

        # 初始化 shp 数组
        shp = np.zeros((4, 8))

        # (-,-) 区域节点 (0,3,4)
        c1 = 0.125 * am1 * am2
        c2 = 0.125 * am2 * am3
        c3 = 0.125 * am1 * am3

        shp[0, 0] = -c2
        shp[0, 1] = c2
        shp[1, 0] = -c3
        shp[1, 3] = c3
        shp[2, 0] = -c1
        shp[2, 4] = c1
        shp[3, 0] = c1 * am3
        shp[3, 4] = c1 * ap3

        # (+,+) 区域节点 (2,6,7)
        c1 = 0.125 * ap1 * ap2
        c2 = 0.125 * ap2 * ap3
        c3 = 0.125 * ap1 * ap3

        shp[0, 7] = -c2
        shp[0, 6] = c2
        shp[1, 5] = -c3
        shp[1, 6] = c3
        shp[2, 2] = -c1
        shp[2, 6] = c1
        shp[3, 2] = c1 * am3
        shp[3, 6] = c1 * ap3

        # (-,+) 区域节点 (3,4,7)
        c1 = 0.125 * am1 * ap2
        c2 = 0.125 * am2 * ap3
        c3 = 0.125 * am1 * ap3

        shp[0, 4] = -c2
        shp[0, 5] = c2
        shp[1, 4] = -c3
        shp[1, 7] = c3
        shp[2, 3] = -c1
        shp[2, 7] = c1
        shp[3, 3] = c1 * am3
        shp[3, 7] = c1 * ap3

        # (+,-) 区域节点 (1,2,5)
        c1 = 0.125 * ap1 * am2
        c2 = 0.125 * ap2 * am3
        c3 = 0.125 * ap1 * am3

        shp[0, 3] = -c2
        shp[0, 2] = c2
        shp[1, 1] = -c3
        shp[1, 2] = c3
        shp[2, 1] = -c1
        shp[2, 5] = c1
        shp[3, 1] = c1 * am3
        shp[3, 5] = c1 * ap3

        # 计算 Jacobian 矩阵
        xs = np.zeros((3, 3))
        xl = self.node_coords.T  # 节点坐标数组 [3x8]

        # dx/dξ, dx/dη, dx/dζ
        xs[0, 0] = (xl[0, 1] - xl[0, 0]) * shp[0, 1] \
                   + (xl[0, 2] - xl[0, 3]) * shp[0, 2] \
                   + (xl[0, 5] - xl[0, 4]) * shp[0, 5] \
                   + (xl[0, 6] - xl[0, 7]) * shp[0, 6]

        xs[1, 0] = (xl[1, 1] - xl[1, 0]) * shp[0, 1] \
                   + (xl[1, 2] - xl[1, 3]) * shp[0, 2] \
                   + (xl[1, 5] - xl[1, 4]) * shp[0, 5] \
                   + (xl[1, 6] - xl[1, 7]) * shp[0, 6]

        xs[2, 0] = (xl[2, 1] - xl[2, 0]) * shp[0, 1] \
                   + (xl[2, 2] - xl[2, 3]) * shp[0, 2] \
                   + (xl[2, 5] - xl[2, 4]) * shp[0, 5] \
                   + (xl[2, 6] - xl[2, 7]) * shp[0, 6]

        # dy/dξ, dy/dη, dy/dζ
        xs[0, 1] = (xl[0, 2] - xl[0, 1]) * shp[1, 2] \
                   + (xl[0, 3] - xl[0, 0]) * shp[1, 3] \
                   + (xl[0, 6] - xl[0, 5]) * shp[1, 6] \
                   + (xl[0, 7] - xl[0, 4]) * shp[1, 7]

        xs[1, 1] = (xl[1, 2] - xl[1, 1]) * shp[1, 2] \
                   + (xl[1, 3] - xl[1, 0]) * shp[1, 3] \
                   + (xl[1, 6] - xl[1, 5]) * shp[1, 6] \
                   + (xl[1, 7] - xl[1, 4]) * shp[1, 7]

        xs[2, 1] = (xl[2, 2] - xl[2, 1]) * shp[1, 2] \
                   + (xl[2, 3] - xl[2, 0]) * shp[1, 3] \
                   + (xl[2, 6] - xl[2, 5]) * shp[1, 6] \
                   + (xl[2, 7] - xl[2, 4]) * shp[1, 7]

        # dz/dξ, dz/dη, dz/dζ
        xs[0, 2] = (xl[0, 4] - xl[0, 0]) * shp[2, 4] \
                   + (xl[0, 5] - xl[0, 1]) * shp[2, 5] \
                   + (xl[0, 6] - xl[0, 2]) * shp[2, 6] \
                   + (xl[0, 7] - xl[0, 3]) * shp[2, 7]

        xs[1, 2] = (xl[1, 4] - xl[1, 0]) * shp[2, 4] \
                   + (xl[1, 5] - xl[1, 1]) * shp[2, 5] \
                   + (xl[1, 6] - xl[1, 2]) * shp[2, 6] \
                   + (xl[1, 7] - xl[1, 3]) * shp[2, 7]

        xs[2, 2] = (xl[2, 4] - xl[2, 0]) * shp[2, 4] \
                   + (xl[2, 5] - xl[2, 1]) * shp[2, 5] \
                   + (xl[2, 6] - xl[2, 2]) * shp[2, 6] \
                   + (xl[2, 7] - xl[2, 3]) * shp[2, 7]

        # 计算伴随矩阵 (ad) = cofactor(xs)
        ad = np.zeros((3, 3))
        ad[0, 0] = xs[1, 1] * xs[2, 2] - xs[1, 2] * xs[2, 1]
        ad[0, 1] = xs[2, 1] * xs[0, 2] - xs[2, 2] * xs[0, 1]
        ad[0, 2] = xs[0, 1] * xs[1, 2] - xs[0, 2] * xs[1, 1]

        ad[1, 0] = xs[1, 2] * xs[2, 0] - xs[1, 0] * xs[2, 2]
        ad[1, 1] = xs[2, 2] * xs[0, 0] - xs[2, 0] * xs[0, 2]
        ad[1, 2] = xs[0, 2] * xs[1, 0] - xs[0, 0] * xs[1, 2]

        ad[2, 0] = xs[1, 0] * xs[2, 1] - xs[1, 1] * xs[2, 0]
        ad[2, 1] = xs[2, 0] * xs[0, 1] - xs[2, 1] * xs[0, 0]
        ad[2, 2] = xs[0, 0] * xs[1, 1] - xs[0, 1] * xs[1, 0]

        # 计算 Jacobian 行列式
        xsj = xs[0, 0] * ad[0, 0] + xs[0, 1] * ad[1, 0] + xs[0, 2] * ad[2, 0]
        if abs(xsj) < 1e-15:
            raise ValueError("Zero or near-zero Jacobian determinant")

        rxsj = 1.0 / xsj

        # 计算逆 Jacobian 矩阵
        xs_inv = np.zeros((3, 3))
        for j in range(3):
            for i in range(3):
                xs_inv[i, j] = ad[i, j] * rxsj

        # 将自然导数转换为全局导数
        for k in range(8):
            # 临时变量
            c1 = shp[0, k]
            c2 = shp[1, k]
            c3 = shp[2, k]

            # 全局坐标导数
            shp[0, k] = c1 * xs_inv[0, 0] + c2 * xs_inv[1, 0] + c3 * xs_inv[2, 0]
            shp[1, k] = c1 * xs_inv[0, 1] + c2 * xs_inv[1, 1] + c3 * xs_inv[2, 1]
            shp[2, k] = c1 * xs_inv[0, 2] + c2 * xs_inv[1, 2] + c3 * xs_inv[2, 2]

        return xsj, shp

    def computeBbar(self, node, shp, shpBar):
        """计算 B-bar 矩阵（体积锁定修正）"""
        one3 = 1.0 / 3.0
        Bbar = np.zeros((6, 3))

        # 提取当前节点的导数
        dNdx = shp[0, node]
        dNdy = shp[1, node]
        dNdz = shp[2, node]

        # 提取当前节点的平均导数
        dNBar_dx = shpBar[0, node]
        dNBar_dy = shpBar[1, node]
        dNBar_dz = shpBar[2, node]

        # 偏差部分 (deviatoric)
        Bdev = np.array([
            [2.0 * dNdx, -dNdy, -dNdz],
            [-dNdx, 2.0 * dNdy, -dNdz],
            [-dNdx, -dNdy, 2.0 * dNdz]
        ])

        # 体积部分 (volumetric)
        BbarVol = np.array([
            [dNBar_dx, dNBar_dy, dNBar_dz],
            [dNBar_dx, dNBar_dy, dNBar_dz],
            [dNBar_dx, dNBar_dy, dNBar_dz]
        ])

        # 组合法向项（前3行）
        for i in range(3):
            for j in range(3):
                Bbar[i, j] = one3 * (Bdev[i, j] + BbarVol[i, j])

        # 剪切项（后3行）
        Bbar[3, 0] = dNdy  # gamma_xy: ε12
        Bbar[3, 1] = dNdx

        Bbar[4, 1] = dNdz  # gamma_yz: ε23
        Bbar[4, 2] = dNdy

        Bbar[5, 0] = dNdz  # gamma_zx: ε31
        Bbar[5, 2] = dNdx

        return Bbar


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
