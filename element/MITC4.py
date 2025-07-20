#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from abc import ABC
import numpy as np
import math
from enum import Enum


class MaterialKey(Enum):
    """ 材料参数关键字集合 """
    E = 1
    Density = 2
    Niu = 3
    G = 4
    Thickness = 5
    Conductivity = 6  # 热传导系数
    Expansion = 7
    SpecificHeat = 8
    Area = 9


def shape2d(ss, tt, x):
    """
    计算形函数导数、Jacobian行列式、雅可比逆矩阵
    输入：
        ss, tt : 当前高斯点在局部坐标中
        x : 节点坐标，形如 2x4 的 ndarray
    返回：
        shp : 3x4 ndarray，其中 shp[0] 为 dN/dx，shp[1] 为 dN/dy，shp[2] 为 N
        xsj : Jacobian 行列式
    """
    s = np.array([-0.5, 0.5, 0.5, -0.5])
    t = np.array([-0.5, -0.5, 0.5, 0.5])

    shp = np.zeros((3, 4))

    # 临时雅可比矩阵 xs: 2x2
    xs = np.zeros((2, 2))

    for i in range(4):
        shp[2, i] = (0.5 + s[i] * ss) * (0.5 + t[i] * tt)
        shp[0, i] = s[i] * (0.5 + t[i] * tt)
        shp[1, i] = t[i] * (0.5 + s[i] * ss)

    for i in range(2):
        for j in range(2):
            for k in range(4):
                xs[i, j] += x[i, k] * shp[j, k]

    xsj = xs[0, 0] * xs[1, 1] - xs[0, 1] * xs[1, 0]

    jinv = 1.0 / xsj
    sx = np.array([
        [xs[1, 1], -xs[0, 1]],
        [-xs[1, 0], xs[0, 0]]
    ]) * jinv

    # 转换到全局导数
    for i in range(4):
        temp = shp[0, i] * sx[0, 0] + shp[1, i] * sx[1, 0]
        shp[1, i] = shp[0, i] * sx[0, 1] + shp[1, i] * sx[1, 1]
        shp[0, i] = temp

    return shp, xsj


def computeBmembrane(node, shp):
    """
    构建膜作用 B 矩阵，输出维度为 3x2
    """
    B = np.zeros((3, 2))
    B[0, 0] = shp[0, node]
    B[1, 1] = shp[1, node]
    B[2, 0] = shp[1, node]
    B[2, 1] = shp[0, node]
    return B


def computeBbend(node, shp):
    """
    构建弯曲 B 矩阵，输出维度为 3x2
    """
    B = np.zeros((3, 2))
    B[0, 1] = -shp[0, node]
    B[1, 0] = shp[1, node]
    B[2, 0] = shp[0, node]
    B[2, 1] = -shp[1, node]
    return B


def computeBdrill(node, shp, g1, g2, g3):
    """
    钻转自由度 B 矩阵计算
    输出：长度为6的ndarray
    """
    B1 = -0.5 * shp[1, node]
    B2 = 0.5 * shp[0, node]
    B6 = -shp[2, node]

    Bdrill = np.zeros(6)

    Bdrill[0] = B1 * g1[0] + B2 * g2[0]
    Bdrill[1] = B1 * g1[1] + B2 * g2[1]
    Bdrill[2] = B1 * g1[2] + B2 * g2[2]

    Bdrill[3] = B6 * g3[0]
    Bdrill[4] = B6 * g3[1]
    Bdrill[5] = B6 * g3[2]

    return Bdrill


def assembleB(Bmembrane, Bbend, Bshear, g1, g2, g3):
    """
    组装总B矩阵（8x6）: 包含膜、弯曲、剪切作用
    """
    B = np.zeros((8, 6))

    Gmem = np.array([g1, g2])  # 2x3
    Gshear = np.zeros((3, 6))
    Gshear[0, 0:3] = g3
    Gshear[1, 3:6] = g1
    Gshear[2, 3:6] = g2

    BmembraneShell = Bmembrane @ Gmem  # (3x2) x (2x3) = 3x3
    BbendShell = Bbend @ Gmem  # (3x2) x (2x3) = 3x3
    BshearShell = Bshear @ Gshear  # (2x3) x (3x6) = 2x6

    # 填入总B矩阵
    B[0:3, 0:3] = BmembraneShell
    B[3:6, 3:6] = BbendShell
    B[6:8, 0:6] = BshearShell

    return B


class ElasticMembranePlateSection:
    five6 = 5.0 / 6.0

    def __init__(self, Em, nu, h, rho=0.0, Ep_modifier=1.0):
        self.Em = Em
        self.nu = nu
        self.h = h
        self.rhoH = rho * h
        self.Ep = Em * Ep_modifier
        self.strain = np.zeros(8)
        self.stress = np.zeros(8)
        self.tangent = np.zeros((8, 8))

    def setTrialSectionDeformation(self, strain):
        self.strain = np.array(strain)
        return 0

    def getSectionDeformation(self):
        return self.strain

    def getRho(self):
        return self.rhoH

    def getStressResultant(self):
        Em = self.Em
        nu = self.nu
        h = self.h
        Ep = self.Ep

        M = Em / (1.0 - nu * nu) * h
        G = 0.5 * Em / (1.0 + nu) * h

        # Membrane
        s = self.stress
        e = self.strain
        s[0] = M * e[0] + nu * M * e[1]
        s[1] = nu * M * e[0] + M * e[1]
        s[2] = G * e[2]

        # Bending
        G *= self.five6 * (Ep / Em)
        D = Ep * h ** 3 / 12.0 / (1.0 - nu * nu)

        s[3] = -(D * e[3] + nu * D * e[4])
        s[4] = -(nu * D * e[3] + D * e[4])
        s[5] = -0.5 * D * (1.0 - nu) * e[5]
        s[6] = G * e[6]
        s[7] = G * e[7]

        return s.copy()

    def getSectionTangent(self):
        Em = self.Em
        nu = self.nu
        h = self.h
        Ep = self.Ep

        M = Em / (1.0 - nu * nu) * h
        G = 0.5 * Em / (1.0 + nu) * h
        D = Ep * h ** 3 / 12.0 / (1.0 - nu * nu)

        Gs = G * self.five6 * (Ep / Em)

        t = np.zeros((8, 8))
        t[0, 0] = t[1, 1] = M
        t[0, 1] = t[1, 0] = nu * M
        t[2, 2] = G

        t[3, 3] = t[4, 4] = -D
        t[3, 4] = t[4, 3] = -nu * D
        t[5, 5] = -0.5 * D * (1.0 - nu)

        t[6, 6] = Gs
        t[7, 7] = Gs

        return t.copy()

    def getInitTangent(self):
        return self.getSectionTangent()


class MITC4():
    """
    MITC4 Element class.
    Reference:
    1. Opensees
    """

    def __init__(self, eid):
        self.nodes_count = 4
        self.vtu_type = "quad"
        self.stiffness = None
        self.stress = None
        self.block_size = 576

    def CalElementDMatrix(self, an_type=None):
        """
        质量单元无需计算D阵
        """
        pass

    def ElementStiffness(self, from_origin=False):
        """
        Reference:
        """
        coor0, coor1 = self.node_coords[:, 0], self.node_coords[:, 1]
        coor2, coor3 = self.node_coords[:, 2], self.node_coords[:, 3]

        v1 = 0.5 * (coor2 + coor1 - coor3 - coor0)
        v2 = 0.5 * (coor3 + coor2 - coor1 - coor0)

        v1 = v1 / np.linalg.norm(v1)
        alpha = np.dot(v2, v1)
        v2_proj = alpha * v1
        v2_new = v2 - v2_proj
        v2 = v2_new / np.linalg.norm(v2_new)
        v3 = np.cross(v1, v2)
        xl = np.zeros((2, 4), dtype=float)
        for ii in range(4):
            xl[0][ii] = self.node_coords[:, ii] @ v1
            xl[1][ii] = self.node_coords[:, ii] @ v2

        # 计算节点坐标间的差值
        dx34 = xl[0][2] - xl[0][3]
        dy34 = xl[1][2] - xl[1][3]
        dx21 = xl[0][1] - xl[0][0]
        dy21 = xl[1][1] - xl[1][0]
        dx32 = xl[0][2] - xl[0][1]
        dy32 = xl[1][2] - xl[1][1]
        dx41 = xl[0][3] - xl[0][0]
        dy41 = xl[1][3] - xl[1][0]
        G = np.zeros((4, 12), dtype=float)
        one_over_four = 0.25
        G[0][0] = -0.5
        G[0][1] = -dy41 * one_over_four
        G[0][2] = dx41 * one_over_four
        G[0][9] = 0.5
        G[0][10] = -dy41 * one_over_four
        G[0][11] = dx41 * one_over_four
        G[1][0] = -0.5
        G[1][1] = -dy21 * one_over_four
        G[1][2] = dx21 * one_over_four
        G[1][3] = 0.5
        G[1][4] = -dy21 * one_over_four
        G[1][5] = dx21 * one_over_four
        G[2][3] = -0.5
        G[2][4] = -dy32 * one_over_four
        G[2][5] = dx32 * one_over_four
        G[2][6] = 0.5
        G[2][7] = -dy32 * one_over_four
        G[2][8] = dx32 * one_over_four
        G[3][6] = 0.5
        G[3][7] = -dy34 * one_over_four
        G[3][8] = dx34 * one_over_four
        G[3][9] = -0.5
        G[3][10] = -dy34 * one_over_four
        G[3][11] = dx34 * one_over_four

        # 提取节点坐标
        x = xl[0]  # x坐标数组 [x0, x1, x2, x3]
        y = xl[1]  # y坐标数组 [y0, y1, y2, y3]
        # 计算中间变量
        Ax = -x[0] + x[1] + x[2] - x[3]
        Bx = x[0] - x[1] + x[2] - x[3]
        Cx = -x[0] - x[1] + x[2] + x[3]
        Ay = -y[0] + y[1] + y[2] - y[3]
        By = y[0] - y[1] + y[2] - y[3]
        Cy = -y[0] - y[1] + y[2] + y[3]
        # 使用atan2避免除零错误
        alph = math.atan2(Ay, Ax)  # 计算α角度
        beta = math.pi / 2 - math.atan2(Cx, Cy)  # 计算β角度
        # 创建旋转矩阵
        Rot = np.zeros((2, 2))
        Rot[0, 0] = math.sin(beta)
        Rot[0, 1] = -math.sin(alph)
        Rot[1, 0] = -math.cos(beta)
        Rot[1, 1] = math.cos(alph)

        Ms = np.zeros((2, 4), dtype=float)
        Bsv = np.zeros((2, 12), dtype=float)
        Bs = np.zeros((2, 12), dtype=float)
        # 初始化变量
        r1 = 0.0
        r2 = 0.0
        r3 = 0.0

        # 高斯积分点
        sg = np.zeros(4, dtype=float)
        tg = np.zeros(4, dtype=float)
        wg = np.ones(4, dtype=float)
        one_over_root3 = 1 / np.sqrt(3)
        sg[0] = -one_over_root3
        sg[1] = one_over_root3
        sg[2] = one_over_root3
        sg[3] = -one_over_root3

        tg[0] = -one_over_root3
        tg[1] = -one_over_root3
        tg[2] = one_over_root3
        tg[3] = one_over_root3

        """
        Gauss Loop
        """
        ngauss = 4
        volume = 0
        nstress = 8
        numnodes = 4
        ndf = 6
        saveB = np.zeros((nstress, ndf, numnodes), dtype=float)
        Bshear = np.zeros((2, 3), dtype=float)
        dvol = np.zeros(ngauss)
        stiff = np.zeros((numnodes * ndf, numnodes * ndf))
        resid = np.zeros(numnodes * ndf)
        for igauss in range(ngauss):
            r1 = Cx + sg[igauss] * Bx
            r3 = Cy + sg[igauss] * By
            r1 = np.sqrt(r1 ** 2 + r3 ** 2)

            r2 = Ax + tg[igauss] * Bx
            r3 = Ay + tg[igauss] * By
            r2 = np.sqrt(r2 ** 2 + r3 ** 2)

            # 计算形函数和 Jacobian
            shp, xsj = shape2d(sg[igauss], tg[igauss], xl)

            # 计算体积微元和总体积
            dvol[igauss] = wg[igauss] * xsj
            volume += dvol[igauss]

            # 构建 Ms 矩阵 (2x4)
            Ms = np.array([
                [1 - sg[igauss], 0, 1 + sg[igauss], 0],
                [0, 1 - tg[igauss], 0, 1 + tg[igauss]]
            ])

            Bsv = Ms @ G  # 矩阵乘法

            # 缩放 Bsv 的每列
            for j in range(12):
                Bsv[0, j] *= r1 / (8 * xsj)
                Bsv[1, j] *= r2 / (8 * xsj)

            Bs = Rot @ Bsv
            epsDrill = 0.0
            strain = np.zeros(nstress)

            for j in range(numnodes):
                Bmembrane = computeBmembrane(j, shp)
                Bbend = computeBbend(j, shp)
                for p in range(3):
                    Bshear[0, p] = Bs[0, j * 3 + p]
                    Bshear[1, p] = Bs[1, j * 3 + p]
                BJ = assembleB(Bmembrane, Bbend, Bshear, v1, v2, v3)
                for p in range(nstress):
                    for q in range(ndf):
                        saveB[p][q][j] = BJ[p, q]
                # ul = self.node_disps[j] - self.init_disp[j] if hasattr(self, 'init_disp') else self.node_disps[j]

                # 计算应变
                # strain += BJ @ ul

                # 计算钻孔应变
                BdrillJ = computeBdrill(j, shp, v1, v2, v3)
                # epsDrill += np.dot(BdrillJ, ul)

            # 获取材料响应
            # success = self.materialPointers[igauss].setTrialSectionDeformation(strain)
            # stress = self.materialPointers[igauss].getStressResultant()

            # 钻孔应力
            # tauDrill = self.Ktt * epsDrill

            # 乘以体积微元
            # stress *= dvol[igauss]
            # tauDrill *= dvol[igauss]

            elastic = ElasticMembranePlateSection(self.cha_dict[MaterialKey.E],
                                                  self.cha_dict[MaterialKey.Niu],
                                                  self.cha_dict[MaterialKey.Thickness])
            dd = elastic.getSectionTangent()
            tangent = elastic.getInitTangent()[:3, :3]
            Ktt = np.min(np.linalg.eigvals(tangent))

            jj = 0  # 当前节点在刚度矩阵中的起始行索引
            for j in range(numnodes):
                # 提取当前节点的B矩阵 (8x6)
                BJ = saveB[:, :, j]

                # 计算转置 B^T (6x8)
                BJtran = BJ.T

                # 计算 B^T D (6x8)
                BJtranD = BJtran @ dd

                # 计算钻孔B矩阵 (6x1)
                BdrillJ = computeBdrill(j, shp, v1, v2, v3)

                # 缩放钻孔矩阵
                BdrillJ_scaled = BdrillJ * (Ktt * dvol[igauss])

                # 内层节点循环
                kk = 0  # 当前k节点在刚度矩阵中的起始列索引
                for k in range(numnodes):
                    # 提取k节点的B矩阵 (8x6)
                    BK = saveB[:, :, k]

                    # 计算钻孔B矩阵 (6x1)
                    BdrillK = computeBdrill(k, shp, v1, v2, v3)

                    # 计算刚度子矩阵: B_j^T D B_k (6x6)
                    stiffJK = BJtranD @ BK  # (6x8) @ (8x6) = 6x6

                    # 添加钻孔贡献: (B_drillJ)^T B_drillK
                    drill_contribution = np.outer(BdrillJ_scaled, BdrillK)

                    # 组装到全局刚度矩阵
                    stiff[jj:jj + ndf, kk:kk + ndf] += stiffJK + drill_contribution

                    kk += ndf  # 移动到下一节点的列索引

                jj += ndf  # 移动到下一节点的行索引

        print("ZZ")

    def CalculateElementStress(self, displacement):
        """
        Calculate element stress, 刚体没有形变, 所以没有应力
        """
        return np.zeros(1), np.zeros(1), np.zeros(1), np.zeros(1), np.zeros(1), np.zeros(1)

    def ElementMass(self):
        pass

    def CalculateBasic(self):
        pass

    def ReCalculateElementStiffness(self):
        pass


if __name__ == "__main__":
    t_ele = MITC4(-1)
    t_ele.cha_dict = {MaterialKey.Niu: 0.3,
                      MaterialKey.E: 1e7,
                      MaterialKey.Thickness: 1,
                      MaterialKey.G: 2e11 / 2 / (1 + 0.3)
                      }
    t_ele.node_coords = np.array([[0, 0, 0],
                                  [1, 0, 0],
                                  [1, 1.6, 0],
                                  [0, 1, 0]], dtype=float).T
    print(t_ele.ElementStiffness())
