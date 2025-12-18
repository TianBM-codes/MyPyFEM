import numpy as np
from typing import Dict
from femdb.Integration import GaussIntegrationPoint


def tri3_shape_grad_xy(coords):
    """
    Tri3: 3节点线性三角形
    coords: (3,2) in global x,y
    return:
        dNdx: (3,)  dNdy: (3,)
        area: scalar
    """
    x1, y1 = coords[0]
    x2, y2 = coords[1]
    x3, y3 = coords[2]

    # 面积 2A
    twoA = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)
    area = 0.5 * twoA
    if area <= 0:
        raise ValueError(f"Tri3 area <= 0 (check node ordering), area={area}")

    # 线性三角形形函数导数常数
    # N1 = (a1 + b1 x + c1 y) / (2A) etc
    b1 = y2 - y3
    c1 = x3 - x2
    b2 = y3 - y1
    c2 = x1 - x3
    b3 = y1 - y2
    c3 = x2 - x1

    dNdx = np.array([b1, b2, b3], dtype=float) / (2 * area)
    dNdy = np.array([c1, c2, c3], dtype=float) / (2 * area)
    return dNdx, dNdy, area


def quad4_shape_and_grad(xi, eta):
    """
    Quad4 bilinear shape functions and derivatives in parent coords (xi,eta)
    return:
        N: (4,)
        dN_dxi: (4,)
        dN_deta: (4,)
    node order assumed: (1)(2)(3)(4) = (-,-), (+,-), (+,+), (-,+)
    """
    N = 0.25 * np.array([
        (1 - xi) * (1 - eta),
        (1 + xi) * (1 - eta),
        (1 + xi) * (1 + eta),
        (1 - xi) * (1 + eta)
    ], dtype=float)

    dN_dxi = 0.25 * np.array([
        -(1 - eta),
        +(1 - eta),
        +(1 + eta),
        -(1 + eta)
    ], dtype=float)

    dN_deta = 0.25 * np.array([
        -(1 - xi),
        -(1 + xi),
        +(1 + xi),
        +(1 - xi)
    ], dtype=float)
    return N, dN_dxi, dN_deta


def quad4_grad_xy(coords, xi, eta):
    """
    coords: (4,2) global x,y
    return:
        N: (4,)
        dNdx: (4,)
        dNdy: (4,)
        detJ: scalar
    """
    N, dN_dxi, dN_deta = quad4_shape_and_grad(xi, eta)

    x = coords[:, 0]
    y = coords[:, 1]

    J = np.array([
        [np.dot(dN_dxi, x), np.dot(dN_deta, x)],
        [np.dot(dN_dxi, y), np.dot(dN_deta, y)]
    ], dtype=float)

    detJ = np.linalg.det(J)
    if detJ <= 0:
        raise ValueError(f"Quad4 detJ <= 0 at (xi,eta)=({xi},{eta}), detJ={detJ}")

    invJ = np.linalg.inv(J)
    # [dNdx; dNdy] = invJ * [dN_dxi; dN_deta]
    grads = invJ @ np.vstack([dN_dxi, dN_deta])  # (2,4)
    dNdx = grads[0, :]
    dNdy = grads[1, :]
    return N, dNdx, dNdy, detJ


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


def Quad8NodeShapeFunction(r, s):
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


def Quad4NodeShapeFunction(r, s):
    """
    4节点4边形形函数
    :param r:
    :param s:
    :return:
    """
    N = np.zeros(4)
    # N[0] = 0.25 * (1 - r) * (1 - s)  # 角点1 (r=-1,s=-1)
    # N[1] = 0.25 * (1 + r) * (1 - s)  # 角点2 (r=1,s=-1)
    # N[2] = 0.25 * (1 + r) * (1 + s)  # 角点3 (r=1,s=1)
    # N[3] = 0.25 * (1 - r) * (1 + s)  # 角点4 (r=-1,s=1)

    # N[0] = 0.25 * (1 - r) * (1 - s)  # 角点1 (r=-1,s=-1)
    # N[1] = 0.25 * (1 - r) * (1 + s)  # 角点4 (r=-1,s=1)
    # N[2] = 0.25 * (1 + r) * (1 - s)  # 角点2 (r=1,s=-1)
    # N[3] = 0.25 * (1 + r) * (1 + s)  # 角点3 (r=1,s=1)

    N[0] = 0.25 * (1 - r) * (1 + s)  # 角点1 (r=-1,s=-1)
    N[1] = 0.25 * (1 + r) * (1 + s)  # 角点2 (r=1,s=-1)
    N[2] = 0.25 * (1 + r) * (1 - s)  # 角点3 (r=1,s=1)
    N[3] = 0.25 * (1 - r) * (1 - s)  # 角点4 (r=-1,s=1)
    return N


_EXTRAPOLATE_4TO4_CACHE = None


def ExtrapolateMatrix4to4():
    """
    外推矩阵, 从4个高斯点应力外推4个节点应力
    :return:
    """
    global _EXTRAPOLATE_4TO4_CACHE
    if _EXTRAPOLATE_4TO4_CACHE is None:
        # sample_pt, _ = GaussIntegrationPoint.GetSamplePointAndWeight(2)
        sample_pt_r = (-0.577350269189626, 0.577350269189626, 0.577350269189626, -0.577350269189626)
        sample_pt_s = (0.577350269189626, 0.577350269189626, -0.577350269189626, -0.577350269189626)
        gauss_coords = [(sample_pt_r[ii], sample_pt_s[ii]) for ii in range(4)]
        A = np.array([Quad4NodeShapeFunction(r, s) for r, s in gauss_coords])
        _EXTRAPOLATE_4TO4_CACHE = np.linalg.inv(A)
    return _EXTRAPOLATE_4TO4_CACHE


_EXTRAPOLATE_3TO3_CACHE = None


def ExtrapolateMatrix3to3():
    """
    获取外推矩阵（带缓存功能）, 从3个高斯点应力外推3个节点应力
    """
    global _EXTRAPOLATE_3TO3_CACHE
    if _EXTRAPOLATE_3TO3_CACHE is None:
        # 获取三角形积分点（3阶）
        gauss_pts, _ = GaussIntegrationPoint.GetTrianglePointAndWeight(3)
        # 构造形函数矩阵
        N_matrix = []
        for r, s in gauss_pts:
            t = 1 - r - s
            N = [
                t * (2 * t - 1),  # N1
                r * (2 * r - 1),  # N2
                s * (2 * s - 1),  # N3
                4 * r * t,  # N4
                4 * r * s,  # N5
                4 * s * t  # N6
            ]
            N_matrix.append(N)
            # 计算伪逆并缓存
        _EXTRAPOLATE_3TO3_CACHE = np.linalg.pinv(np.array(N_matrix))

    return _EXTRAPOLATE_3TO3_CACHE


if __name__ == "__main__":
    t_sample_pt, _ = GaussIntegrationPoint.GetSamplePointAndWeight(2)
    t_gauss_coords = [(t_sample_pt[ri], t_sample_pt[si]) for ri in range(2) for si in range(2)]
    t_A = np.array([Quad4NodeShapeFunction(r, s) for r, s in t_gauss_coords])
    t_invA = np.linalg.inv(t_A)
    print(ExtrapolateMatrix4to4() @ t_A)
