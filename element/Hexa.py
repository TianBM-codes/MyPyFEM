#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pypardiso
import scipy.sparse.linalg

from element.ElementBase import *
import numpy as np
from abc import ABC


class C3D8(ElementBaseClass, ABC):
    """
    8节点六面体单元类(C3D8)
    实现8节点六面体单元的计算功能，包括：
    - 单元刚度矩阵计算
    - 单元应力计算
    - 本构矩阵计算

    参考：
    《有限单元法》王勖成 P323
    """

    def __init__(self, eid=None):
        """
        初始化8节点六面体单元

        Args:
            eid (int, optional): 单元ID. Defaults to None.
        """
        super().__init__(eid)
        self.nodes_count = 8  # 单元节点数
        self.K = np.zeros([24, 24], dtype=float)  # 单元刚度矩阵(24自由度)
        self.vtu_type = "hexahedron"  # VTK可视化类型
        self.gs_count = 8  # 高斯积分点数量
        self.Gaussian_B = []  # 高斯积分点处的应变矩阵列表

    def CalElementDMatrix(self, an_type=None):
        """
        计算本构矩阵(弹性矩阵)
        根据材料属性(弹性模量和泊松比)计算单元本构矩阵

        参考：
        《有限单元法》王勖成 P184

        Args:
            an_type (optional): 分析类型. Defaults to None.
        """
        e = self.cha_dict[MaterialKey.E]  # 弹性模量
        niu = self.cha_dict[MaterialKey.Niu]  # 泊松比

        # 计算弹性矩阵系数
        a = e / ((1 + niu) * (1 - 2 * niu))

        # 构造各向同性弹性矩阵
        self.D = a * np.array([
            [1 - niu, niu, niu, 0, 0, 0],
            [niu, 1 - niu, niu, 0, 0, 0],
            [niu, niu, 1 - niu, 0, 0, 0],
            [0, 0, 0, (1 - 2 * niu) / 2., 0, 0],
            [0, 0, 0, 0, (1 - 2 * niu) / 2., 0],
            [0, 0, 0, 0, 0, (1 - 2 * niu) / 2.]
        ])

    def ElementStiffness(self, from_origin=False):
        """
        计算单元刚度矩阵
        使用8点高斯积分计算单元刚度矩阵

        TODO:
        - 参考视频 https://www.bilibili.com/video/BV19y4y1z76E/ 实现C3D20单元
        - 优化计算效率

        参考：
        《有限单元法》王勖成 P323

        Args:
            from_origin (bool, optional): 是否从原始数据重新计算. Defaults to False.

        Returns:
            np.ndarray: 24x24的单元刚度矩阵
        """
        assert self.node_coords.shape == (8, 3)  # 确保节点坐标是8x3数组
        self.CalElementDMatrix()  # 计算本构矩阵

        # 获取8个高斯积分点的形函数导数和权重
        dNdrs, weights = AllEleTypeDNDrAtGaussianPoint.C3D8

        # 在8个高斯点上积分计算刚度矩阵
        for ii in range(self.gs_count):
            # 计算雅可比矩阵
            J = np.matmul(dNdrs[ii], self.node_coords)
            det_J = np.linalg.det(J)  # 雅可比行列式

            # 计算B矩阵(应变-位移矩阵)
            B_pre = np.linalg.solve(J, dNdrs[ii])  # 解线性方程组计算B矩阵

            # 组装完整的B矩阵(6x24)
            B_at_gs_pt = np.asarray([
                [B_pre[0, 0], 0, 0, B_pre[0, 1], 0, 0, B_pre[0, 2], 0, 0, B_pre[0, 3], 0, 0, B_pre[0, 4], 0, 0, B_pre[0, 5], 0, 0, B_pre[0, 6], 0, 0, B_pre[0, 7], 0, 0],
                [0, B_pre[1, 0], 0, 0, B_pre[1, 1], 0, 0, B_pre[1, 2], 0, 0, B_pre[1, 3], 0, 0, B_pre[1, 4], 0, 0, B_pre[1, 5], 0, 0, B_pre[1, 6], 0, 0, B_pre[1, 7], 0],
                [0, 0, B_pre[2, 0], 0, 0, B_pre[2, 1], 0, 0, B_pre[2, 2], 0, 0, B_pre[2, 3], 0, 0, B_pre[2, 4], 0, 0, B_pre[2, 5], 0, 0, B_pre[2, 6], 0, 0, B_pre[2, 7]],
                [B_pre[1, 0], B_pre[0, 0], 0, B_pre[1, 1], B_pre[0, 1], 0, B_pre[1, 2], B_pre[0, 2], 0, B_pre[1, 3], B_pre[0, 3], 0, B_pre[1, 4], B_pre[0, 4], 0, B_pre[1, 5], B_pre[0, 5], 0, B_pre[1, 6], B_pre[0, 6], 0, B_pre[1, 7], B_pre[0, 7], 0],
                [0, B_pre[2, 0], B_pre[1, 0], 0, B_pre[2, 1], B_pre[1, 1], 0, B_pre[2, 2], B_pre[1, 2], 0, B_pre[2, 3], B_pre[1, 3], 0, B_pre[2, 4], B_pre[1, 4], 0, B_pre[2, 5], B_pre[1, 5], 0, B_pre[2, 6], B_pre[1, 6], 0, B_pre[2, 7], B_pre[1, 7]],
                [B_pre[2, 0], 0, B_pre[0, 0], B_pre[2, 1], 0, B_pre[0, 1], B_pre[2, 2], 0, B_pre[0, 2], B_pre[2, 3], 0, B_pre[0, 3], B_pre[2, 4], 0, B_pre[0, 4], B_pre[2, 5], 0, B_pre[0, 5], B_pre[2, 6], 0, B_pre[0, 6], B_pre[2, 7], 0, B_pre[0, 7]]
            ], dtype=float)

            # 保存当前高斯点的B矩阵
            self.Gaussian_B.append(B_at_gs_pt)

            # 计算当前高斯点对刚度矩阵的贡献
            self.K += B_at_gs_pt.T @ self.D @ B_at_gs_pt * det_J * weights[ii]

        return self.K

    def CalculateElementStress(self, displacement):
        """
        计算单元应力
        基于位移结果计算单元应力

        计算步骤：
        1. 计算高斯点位移
        2. 计算高斯点应力
        3. 将应力外推至节点
        4. 节点平均(在数据库方法中实现)

        参考：
        《有限单元法》王勖成 P168-176

        Args:
            displacement (np.ndarray): 节点位移向量(24x1)

        Returns:
            tuple: 6个应力分量(sigma_xx, sigma_yy, sigma_zz, tau_xy, tau_yz, tau_xz)
        """
        # 定义高斯点位移插值矩阵
        a = 0.49056261216234404  # 0.125*(1+1/np.sqrt(3))**3
        b = 0.13144585576580214  # 0.125*(1+1/np.sqrt(3))*2/3
        c = 0.03522081090086451  # 0.125*(1-1/np.sqrt(3))*2/3
        d = 0.00943738783765593  # 0.125*(1-1/np.sqrt(3))**3

        Global2Gaussian = np.asarray([
            [a, b, c, b, b, c, d, c],
            [b, a, b, c, c, b, c, d],
            [c, b, a, b, d, c, b, c],
            [b, c, b, a, c, d, c, b],
            [b, c, d, c, a, b, c, b],
            [c, b, c, d, b, a, b, c],
            [d, c, b, c, c, b, a, b],
            [c, d, c, b, b, c, b, a]
        ], dtype=float)

        # 计算高斯点位移(24x1)
        gs_dis = np.matmul(Global2Gaussian, displacement.reshape((8, 3))).reshape((24,))

        # 计算高斯点应力
        gs_stress = []
        for ii in range(self.gs_count):
            stress = np.matmul(self.D, np.matmul(self.Gaussian_B[ii], gs_dis))
            gs_stress.append(stress)
        gs_stress = np.asarray(gs_stress)

        # 定义应力外推矩阵(将高斯点应力外推到节点)
        a = 2.549038105676658  # 0.25 * (5 + 3 * np.sqrt(3))
        b = -0.68301270189222  # -0.25 * (np.sqrt(3) + 1)
        c = 0.183012701892219  # 0.25 * (np.sqrt(3) - 1)
        d = -0.04903810567666  # 0.25 * (5 - 3 * np.sqrt(3))

        Gaussian2Global = np.asarray([
            [a, b, c, b, b, c, d, c],
            [b, a, b, c, c, b, c, d],
            [c, b, a, b, d, c, b, c],
            [b, c, b, a, c, d, c, b],
            [b, c, d, c, a, b, c, b],
            [c, b, c, d, b, a, b, c],
            [d, c, b, c, c, b, a, b],
            [c, d, c, b, b, c, b, a]
        ], dtype=float)

        # 计算节点应力
        node_stress = np.matmul(Gaussian2Global, gs_stress)

        # 返回6个应力分量
        return (node_stress[:, 0], node_stress[:, 1], node_stress[:, 2],
                node_stress[:, 3], node_stress[:, 4], node_stress[:, 5])

    def ElementMass(self):
        """计算单元质量矩阵(待实现)"""
        pass

    def CalculateBasic(self):
        """计算基本参数(待实现)"""
        pass

    def ReCalculateElementStiffness(self):
        """重新计算单元刚度矩阵(待实现)"""
        pass