# !/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np


class Kinematics(object):
    def __init__(self):
        """
        DN_x: Spatial gradient of the shape functions.
        Jx_chi: Jacobin of the mapping between spatial and iso parametric domains.
        F: Deformation gradient.
        J: Jacobin of the deformation gradient.
        b: Left Cauchy-Green strain tensor (b).
        Ib: First invariant of b.
        lambda: Principal stretches.
        n: Spatial principal directions.
        """
        self.DN_Dx = None
        self.Jx_chi = None
        self.F = None
        self.J = None
        self.b = None
        self.Ib = None
        self.lambda_ = None
        self.n = None

        self.ngauss = None
        self.ndim = None
        self.n_nodes_elem = None

    def Init(self, ndim, n_nodes_elem, ngauss):
        """
        @param ndim: number of dimension
        @param n_nodes_elem: number of element nodes
        @param ngauss: number of gaussian integration
        @return:
        """
        self.DN_Dx = np.zeros((ndim, n_nodes_elem, ngauss))
        self.Jx_chi = np.zeros((ngauss, 1))
        self.F = np.zeros((ndim, ndim, ngauss))
        self.J = np.zeros((ngauss, 1))
        self.b = np.zeros((ndim, ndim, ngauss))
        self.Ib = np.zeros((ngauss, 1))
        self.lambda_ = np.zeros((ndim, ngauss))
        self.n = np.zeros((ndim, ndim, ngauss))
        self.ngauss = ngauss
        self.ndim = ndim
        self.n_nodes_elem = n_nodes_elem

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
        for ii in range(self.ngauss):
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
            V, D = np.linalg.eig(b)
            """
            Storage of variables
            """
            self.DN_Dx[:, :, ii] = DN_Dx
            self.Jx_chi[ii] = np.abs(np.linalg.det(Dx_Dchi))
            self.F[:, ii] = F
            self.J[:, ii] = J
            self.b[:, ii] = b
            self.Ib[ii] = np.trace(b)
            self.lambda_[:, ii] = np.sqrt(np.diag(D))
            self.n[:, :, ii] = V

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
