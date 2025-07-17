#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from element.ElementBase import *
from abc import ABC
import numpy as np


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


class MITC4(ElementBaseClass, ABC):
    """
    MITC4 Element class.
    Reference:
    1. Opensees
    """

    def __init__(self, eid):
        super().__init__(eid)
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
        assert self.node_coords.shape == (1, 3)

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
    t_ele.ele_mat_dict = {MaterialKey.E: 1, MaterialKey.Area: np.sqrt(3)}
    t_ele.node_coords = np.array([[0, 0, 0],
                                  [1, 1, 1]], dtype=float)
    print(t_ele.ElementStiffness())
