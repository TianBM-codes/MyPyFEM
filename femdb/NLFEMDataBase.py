#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from femdb.LoadCase import LoadCase
from femdb.Material import Material
from femdb.ElementGroup import *
from femdb.Mesh import Mesh
from femdb.Geom import Geom
from typing import Dict
from femdb.Boundary import BoundaryBase
from femdb.SolveControl import SolveControl
from femdb.Kinematics import Kinematics
from femdb.LoadCase import RightHandItem


class IdentityTensor(object):
    def __init__(self):
        """
        Obtain entities which will be constant and only computed once.
        components of fourth order isotropic tensors
        c1 = delta(i,j)*delta(k,l)
        c2 = delta(i,k)*delta(j,l) + delta(i,l)*delta(j,k)
        (see textbook example 2.8)
        """
        self.I = None
        self.c1 = None
        self.c2 = None

    def InitVariant(self, dimension):
        self.I = np.eye(dimension)
        self.c1 = np.zeros((dimension, dimension, dimension, dimension))
        self.c2 = np.zeros((dimension, dimension, dimension, dimension))
        for l in range(dimension):
            for k in range(dimension):
                for j in range(dimension):
                    for i in range(dimension):
                        self.c1[i, j, k, l] = self.c1[i, j, k, l] + self.I[i, j] * self.I[k, l]
                        self.c2[i, j, k, l] = (self.c2[i, j, k, l] +
                                               self.I[i, k] * self.I[j, l] +
                                               self.I[i, l] * self.I[j, k])


class AuxVariant(object):
    def __init__(self):
        self.ngauss = None
        self.n_dofs_elem = None
        self.weight = None
        self.DN_Dchi = None
        self.n_nodes_element = None
        self.n_face_dofs_elem = None
        self.boundary_ngauss = None


class GlobalK(object):
    def __init__(self):
        self.indexi = None
        self.indexj = None
        self.counter = None
        self.stiffness = None


class NLFEMDataBase(object):
    """
    有限元数据库, 实例化的都存储在这里
    """

    _instance = None  # 类变量用于存储唯一的实例
    _initialized = False

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(NLFEMDataBase, cls).__new__(cls)
        return cls._instance

    def __init__(self):
        if not self._initialized:
            self.file_path = None
            self.title = None
            self.Mesh = Mesh()
            self.Geom = Geom()
            self.LoadCase = LoadCase()
            self.BC = BoundaryBase()
            self.Material = Material()
            self.Kinematics = Kinematics()
            self.Dimension = AnalyseDimension.NoAssign
            self.SolveControl = SolveControl()
            self.ElementGroupHash: Dict[int, ElementGroup] = {}
            self.identity_tensor = IdentityTensor()

            self.right_hand_item = RightHandItem()
            self.kinematics = Kinematics()
            self.global_k = GlobalK()

            self._initialized = True

    def SetDimensionVariant(self, dimension: int):
        if dimension in [2, 3]:
            self.identity_tensor.InitVariant(dimension)
        else:
            raise NoSupportDimension(dimension)
