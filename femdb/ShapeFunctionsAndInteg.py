import numpy as np
from typing import Dict


class IntegForm2D2P:
    """二维2点高斯积分格式"""

    def __init__(self):
        tmp1 = np.sqrt(3) / 3
        tmp2 = -np.sqrt(3) / 3
        self.point = np.array([
            [tmp2, tmp2],
            [tmp1, tmp2],
            [tmp1, tmp1],
            [tmp2, tmp1]
        ])
        self.weight = np.array([1, 1, 1, 1])


class IntegForm2D1P:
    """二维1点高斯积分格式"""

    def __init__(self):
        self.point = np.array([[0, 0]])
        self.weight = np.array([4])


class ShapeFunction:
    """形状函数基类"""

    def __init__(self):
        self.DNDxi = None  # 导数矩阵
        self.N = None  # 形函数值
        self.wgt = None  # 权重


def shpFunc4NodeIncom(integForm: IntegForm2D2P) -> Dict[int, ShapeFunction]:
    """4节点不协调单元的形状函数"""
    shpFunc = {}
    ngauss = integForm.point.shape[0]

    for i in range(ngauss):
        sf = ShapeFunction()
        xi, eta = integForm.point[i]

        sf.DNDxi = np.zeros((2, 2))
        sf.DNDxi[0, 0] = -2 * xi
        sf.DNDxi[1, 1] = -2 * eta

        shpFunc[i] = sf

    return shpFunc


def shpFunc2D4Node(integForm) -> Dict[int, ShapeFunction]:
    """二维4节点单元的形状函数"""
    shpFunc = {}
    ngauss = integForm.point.shape[0]

    for i in range(ngauss):
        sf = ShapeFunction()
        xi, eta = integForm.point[i]

        sf.DNDxi = np.zeros((2, 4))
        sf.N = np.zeros(4)

        # 导数计算
        sf.DNDxi[0, 0] = -0.25 * (1 - eta)
        sf.DNDxi[1, 0] = -0.25 * (1 - xi)

        sf.DNDxi[0, 1] = 0.25 * (1 - eta)
        sf.DNDxi[1, 1] = -0.25 * (1 + xi)

        sf.DNDxi[0, 2] = 0.25 * (1 + eta)
        sf.DNDxi[1, 2] = 0.25 * (1 + xi)

        sf.DNDxi[0, 3] = -0.25 * (1 + eta)
        sf.DNDxi[1, 3] = 0.25 * (1 - xi)

        # 形函数值
        sf.N[0] = 0.25 * (1 - xi) * (1 - eta)
        sf.N[1] = 0.25 * (1 + xi) * (1 - eta)
        sf.N[2] = 0.25 * (1 + xi) * (1 + eta)
        sf.N[3] = 0.25 * (1 - xi) * (1 + eta)

        sf.wgt = integForm.weight[i]
        shpFunc[i] = sf

    return shpFunc


def quad8NodeFunctions(r, s):
    """8节点四边形单元形函数"""
    N = np.zeros(8)
    # 角点形函数
    N[0] = 0.25 * (1 - r) * (1 - s) * (-r - s - 1)  # 节点1
    N[1] = 0.25 * (1 + r) * (1 - s) * (r - s - 1)  # 节点2
    N[2] = 0.25 * (1 + r) * (1 + s) * (r + s - 1)  # 节点3
    N[3] = 0.25 * (1 - r) * (1 + s) * (-r + s - 1)  # 节点4
    # 边中点形函数
    N[4] = 0.5 * (1 - r ** 2) * (1 - s)  # 节点5（底边）
    N[5] = 0.5 * (1 + s) * (1 - r ** 2)  # 节点6（右边）
    N[6] = 0.5 * (1 - r ** 2) * (1 + s)  # 节点7（顶边）
    N[7] = 0.5 * (1 - s) * (1 - r ** 2)  # 节点8（左边）
    return N


def shpFunc2D8Node(integForm: IntegForm2D2P) -> Dict[int, ShapeFunction]:
    """二维8节点Serendipity单元的形状函数"""
    shpFunc = {}
    ngauss = integForm.point.shape[0]

    for i in range(ngauss):
        sf = ShapeFunction()
        xi, eta = integForm.point[i]

        sf.DNDxi = np.zeros((2, 8))
        sf.N = np.zeros(8)

        # 中间节点导数 (5-8节点)
        sf.DNDxi[0, 4] = -1 * (1 - eta) * xi  # 5节点
        sf.DNDxi[1, 4] = -0.5 * (1 - xi ** 2)

        sf.DNDxi[0, 5] = 0.5 * (1 - eta ** 2)  # 6节点
        sf.DNDxi[1, 5] = -1 * (1 + xi) * eta

        sf.DNDxi[0, 6] = -1 * (1 + eta) * xi  # 7节点
        sf.DNDxi[1, 6] = 0.5 * (1 - xi ** 2)

        sf.DNDxi[0, 7] = -0.5 * (1 - eta ** 2)  # 8节点
        sf.DNDxi[1, 7] = -1 * (1 - xi) * eta

        # 角节点导数 (1-4节点)
        sf.DNDxi[0, 0] = -0.25 * (1 - eta) - 0.5 * (sf.DNDxi[0, 4] + sf.DNDxi[0, 7])
        sf.DNDxi[1, 0] = -0.25 * (1 - xi) - 0.5 * (sf.DNDxi[1, 4] + sf.DNDxi[1, 7])

        sf.DNDxi[0, 1] = 0.25 * (1 - eta) - 0.5 * (sf.DNDxi[0, 4] + sf.DNDxi[0, 5])
        sf.DNDxi[1, 1] = -0.25 * (1 + xi) - 0.5 * (sf.DNDxi[1, 4] + sf.DNDxi[1, 5])

        sf.DNDxi[0, 2] = 0.25 * (1 + eta) - 0.5 * (sf.DNDxi[0, 5] + sf.DNDxi[0, 6])
        sf.DNDxi[1, 2] = 0.25 * (1 + xi) - 0.5 * (sf.DNDxi[1, 5] + sf.DNDxi[1, 6])

        sf.DNDxi[0, 3] = -0.25 * (1 + eta) - 0.5 * (sf.DNDxi[0, 6] + sf.DNDxi[0, 7])
        sf.DNDxi[1, 3] = 0.25 * (1 - xi) - 0.5 * (sf.DNDxi[1, 6] + sf.DNDxi[1, 7])

        # 形函数值
        sf.N[4] = 0.5 * (1 - eta) * (1 - xi ** 2)  # 5节点
        sf.N[5] = 0.5 * (1 + xi) * (1 - eta ** 2)  # 6节点
        sf.N[6] = 0.5 * (1 + eta) * (1 - xi ** 2)  # 7节点
        sf.N[7] = 0.5 * (1 - xi) * (1 - eta ** 2)  # 8节点

        # 角节点形函数 (1-4节点)
        sf.N[0] = 0.25 * (1 - xi) * (1 - eta) - 0.5 * (sf.N[4] + sf.N[7])
        sf.N[1] = 0.25 * (1 + xi) * (1 - eta) - 0.5 * (sf.N[4] + sf.N[5])
        sf.N[2] = 0.25 * (1 + xi) * (1 + eta) - 0.5 * (sf.N[5] + sf.N[6])
        sf.N[3] = 0.25 * (1 - xi) * (1 + eta) - 0.5 * (sf.N[6] + sf.N[7])

        sf.wgt = integForm.weight[i]
        shpFunc[i] = sf

    return shpFunc
