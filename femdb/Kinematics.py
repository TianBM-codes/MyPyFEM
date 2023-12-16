# !/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np


class Kinematics(object):
    """
    运动学相关变量存储
    """

    def __init__(self):
        self.DN_Dx = None
        self.Jx_chi = None
        self.F = None
        self.J = None
        self.b = None
        self.Ib = None
        self.lambda_ = None

    def ComputeGradients(self,
                         xlocal: np.ndarray,
                         Xlocal: np.ndarray,
                         DN_Dchi: np.ndarray):
        """
        对于有n个节点的单元,每个节点有K个自由度，p个形函数，q个参数坐标
        1. xlocal和Xlocal的shape为K*N
        2. DN_Dchi的shape为p*q, 等参单元情况下shape为N*K

        Reference:
        Javier Bonet P236 公式(9.8)
        """
        """
        Derivative of shape functions with respect to initial coordinates
        """
        DX_Dchi = np.matmul(Xlocal, DN_Dchi)
        DN_DX = np.matmul(DN_Dchi, np.linalg.inv(DX_Dchi))
        """
        current coordinates
        """
        Dx_Dchi = np.matmul(xlocal, DN_Dchi)
        DN_Dx = np.matmul(DN_Dchi, np.linalg.inv(Dx_Dchi))
        """
        Compute various strain measures
        """
        F = np.matmul(xlocal, DN_DX)
        J = np.linalg.det(F)
        b = np.matmul(F, F.T)
        """
        Storage of variables
        """
        self.DN_Dx = DN_Dx
        self.F = F
        self.J = J
        self.b = b

    def PrintVariables(self):
        print(f"F({self.F.shape}):\n", self.F)
        print(f"DN_Dx({self.DN_Dx.shape}\n:", self.DN_Dx)
        print(f"b({self.b.shape}):\n", self.b)
        print(f"J:\n", self.J)


if __name__ == "__main__":
    """
    测试运动学的算例取自Javier书EXAMPLE9.1 9.2
    """
    X = np.asarray([[0, 4, 0], [0, 0, 3]], dtype=float)
    x = np.asarray([[2, 10, 10], [3, 3, 9]], dtype=float)
    kine = Kinematics()
    kine.ComputeGradients(xlocal=x, Xlocal=X, DN_Dchi=np.asarray([[-1, -1], [1, 0], [0, 1]], dtype=float))
    kine.PrintVariables()
