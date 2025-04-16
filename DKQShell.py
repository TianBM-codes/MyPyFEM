import numpy as np
from typing import Tuple, Dict, Any


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


def shpFunc2D4Node(integForm: IntegForm2D2P) -> Dict[int, ShapeFunction]:
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


class Q4Mem:
    """
    膜单元，包含不协调应变方法和钻孔自由度(小刚度系数)
    """
    numNodeOfEle = 4
    NDIM = 3
    ngauss = 4
    numAddDof = 4

    def __init__(self):
        self.integ = IntegForm2D2P()
        self.shpFunc = shpFunc2D4Node(self.integ)

        # 用于钻孔自由度
        self.shpFunc8 = shpFunc2D8Node(self.integ)

        # 用于计算不协调应变插值矩阵[G]
        self.gauss_centor = IntegForm2D1P()
        self.shpFunc_centor = shpFunc2D4Node(self.gauss_centor)
        self.shpFunc_Incom = shpFunc4NodeIncom(self.integ)

        # 单元局部坐标
        self.LocCoord = None
        self.Dm = None  # 材料矩阵
        self.h = None  # 厚度
        self.G = None  # 剪切模量

    def set_coord(self, LoX0: np.ndarray) -> None:
        """设置单元局部坐标"""
        self.LocCoord = LoX0

    def strain_matrix(self, igauss: int) -> Tuple[np.ndarray, float]:
        """计算应变矩阵B和雅可比行列式"""
        J = self.shpFunc[igauss].DNDxi @ self.LocCoord[:, :2]
        DNDx = np.linalg.solve(J, self.shpFunc[igauss].DNDxi)

        B = np.zeros((3, 12))
        for iNode in range(4):
            B[0, (iNode * 3)] = DNDx[0, iNode]
            B[1, (iNode * 3) + 1] = DNDx[1, iNode]
            B[2, (iNode * 3)] = DNDx[1, iNode]
            B[2, (iNode * 3) + 1] = DNDx[0, iNode]

        j = np.linalg.det(J)
        return B, j

    def incom_strain_matrix(self, igauss: int, j: float, J0: np.ndarray) -> np.ndarray:
        """计算不协调应变矩阵G"""
        DMbarDx = (np.linalg.det(J0) / j) * np.linalg.solve(J0, self.shpFunc_Incom[igauss].DNDxi)

        G = np.zeros((3, 4))
        G[0, 0] = DMbarDx[0, 0]
        G[1, 1] = DMbarDx[1, 0]
        G[2, 0] = DMbarDx[1, 0]
        G[2, 1] = DMbarDx[0, 0]
        G[0, 2] = DMbarDx[0, 1]
        G[1, 3] = DMbarDx[1, 1]
        G[2, 2] = DMbarDx[1, 1]
        G[2, 3] = DMbarDx[0, 1]

        return G

    def set_constid(self, Prop: Any) -> None:
        """设置材料属性"""
        self.Dm = Prop.A
        self.h = Prop.h
        self.G = Prop.G

    def stiff_calc(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """计算刚度矩阵"""
        # 膜部分
        Ke = np.zeros((12, 12))
        Kmud = np.zeros((12, 4))
        Kmdd = np.zeros((4, 4))

        J0 = self.shpFunc_centor[0].DNDxi @ self.LocCoord[:, :2]
        V = 0.0

        for igauss in range(4):
            B, j = self.strain_matrix(igauss)
            G = self.incom_strain_matrix(igauss, j, J0)

            wgt = 1.0
            Ke += B.T @ self.Dm @ B * j * wgt
            Kmud += B.T @ self.Dm @ G * j * wgt
            Kmdd += G.T @ self.Dm @ G * j * wgt
            V += j * wgt * self.h

            # 凝聚
        Ke -= Kmud @ np.linalg.inv(Kmdd) @ Kmud.T

        # 钻孔部分
        delta2 = 0.01
        DN0Dx = np.linalg.solve(J0, self.shpFunc_centor[0].DNDxi)
        Q = 0.5 * np.array([
            [-DN0Dx[1, 0], DN0Dx[0, 0], -0.5, -DN0Dx[1, 1], DN0Dx[0, 1], -0.5,
             -DN0Dx[1, 2], DN0Dx[0, 2], -0.5, -DN0Dx[1, 3], DN0Dx[0, 3], -0.5]
        ])
        Sr = (delta2 * V * self.G) * (Q.T @ Q)
        Ke += Sr

        # 为旋转自由度添加小刚度
        lambda_ = 1e-4
        Ke[2, 2] += self.G * V * lambda_
        Ke[5, 5] += self.G * V * lambda_
        Ke[8, 8] += self.G * V * lambda_
        Ke[11, 11] += self.G * V * lambda_

        return Ke, Kmud, Kmdd

    def strain_calc(self, Loc_U: np.ndarray) -> np.ndarray:
        """计算应变"""
        strain = np.zeros((6, 4))
        J0 = self.shpFunc_centor[0].DNDxi @ self.LocCoord[:, :2]

        for iGauss in range(self.ngauss):
            B, j = self.strain_matrix(iGauss)
            G = self.incom_strain_matrix(iGauss, j, J0)
            strain[[0, 1, 3], iGauss] = B @ Loc_U

        return strain


if __name__ == "__main__":
    # 创建材料属性
    class MaterialProps:
        def __init__(self):
            self.A = np.array([
                [2.197802197802198e+09, 6.593406593406594e+08, 0],
                [6.593406593406594e+08, 2.197802197802198e+09, 0],
                [0, 0, 7.692307692307693e+08]
            ])
            self.h = 0.01  # 厚度
            self.G = 7.692307692307692e+10  # 剪切模量


    # 创建单元
    element = Q4Mem()

    # 设置单元坐标 (示例坐标)
    coords = np.array([
        [0, 0, 0],
        [1, 0, 0],
        [1, 1, 0],
        [0, 1, 0]
    ])
    element.set_coord(coords)

    # 设置材料属性
    props = MaterialProps()
    element.set_constid(props)

    # 计算刚度矩阵
    Ke, Kmud, Kmdd = element.stiff_calc()
    print("单元刚度矩阵形状:", Ke.shape)
