#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import abc
from femdb.NLFEMDataBase import NLFEMDataBase
from utils.GlobalEnum import *


def flagshyp_boundary_codes(code, dim):
    """
    Returns the degrees of freedom based on the boundary code and dimension.
    """
    if code == 0:
        dof = np.array([0, 0, 0])
    elif code == 1:
        dof = np.array([1, 0, 0])
    elif code == 2:
        dof = np.array([0, 1, 0])
    elif code == 3:
        dof = np.array([1, 1, 0])
    elif code == 4:
        dof = np.array([0, 0, 1])
    elif code == 5:
        dof = np.array([1, 0, 1])
    elif code == 6:
        dof = np.array([0, 1, 1])
    elif code == 7:
        dof = np.array([1, 1, 1])
    else:
        raise ValueError("Invalid boundary code provided.")

    # Adjust the degrees of freedom array based on the dimension
    if dim == 2:
        dof = dof[:2]  # Remove the third element for 2D cases

    return dof


class BoundaryBase:
    def __init__(self):
        self.icode = None
        self.free_dof = None
        self.fixed_dof = None

    @abc.abstractmethod
    def FindFixedAndFreeDofs(self):
        pass


class FlagSHyPBoundary(BoundaryBase):
    def __init__(self):
        super().__init__()

    def FindFixedAndFreeDofs(self):
        fem_db = NLFEMDataBase()
        n_dofs = fem_db.Mesh.n_dofs
        dim = GlobalInfor[GlobalVariant.Dimension]
        npoin = fem_db.Geom.node_count

        self.free_dof = np.arange(1, n_dofs + 1)  # 使用从 1 开始的索引
        self.fixed_dof = np.zeros(n_dofs, dtype=int)

        for inode in range(npoin):
            fixed_dofs = flagshyp_boundary_codes(self.icode[inode], dim)
            self.fixed_dof[(inode * dim):(inode * dim + dim)] = fixed_dofs

        self.fixed_dof = self.fixed_dof * self.free_dof
        self.fixed_dof = self.fixed_dof[self.fixed_dof > 0]
        self.free_dof = np.delete(self.free_dof, self.fixed_dof - 1)
