# !/usr/bin/env python3
# -*- coding: utf-8 -*-

from femdb.NLFEMDataBase import NLFEMDataBase
from utils.GlobalEnum import *

"""
Single instance mode, convenient for programming, Connect Database
"""
fem_db = NLFEMDataBase()
KINEMATICS = fem_db.kinematics
dim = GlobalInfor[GlobalVariant.Dimension]
element_indexi = fem_db.global_k.indexi
element_indexj = fem_db.global_k.indexj
element_stiffness = fem_db.global_k.stiffness
IDENTITY_TENSOR = fem_db.identity_tensor
T_int = fem_db.right_hand_item.T_int

fem_db = NLFEMDataBase()
MAT = fem_db.Material.Mat
MESH = fem_db.Mesh


def NewtonRaphsonAlgorithm():
    pass
