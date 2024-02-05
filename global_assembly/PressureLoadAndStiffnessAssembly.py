# !/usr/bin/env python3
# -*- coding: utf-8 -*-

from femdb.NLFEMDataBase import NLFEMDataBase
from element.ElementBase import ElementBaseClass
from GlobalEnum import *
from femdb.NLDomain import NLDomain
from femdb.Material import *
import numpy as np

"""
Single instance mode, convenient for programming, Connect Database
"""
nl_domain = NLDomain()
PLAST = nl_domain.plastics
KINEMATICS = nl_domain.kinematics
dim = GlobalInfor[GlobalVariant.Dimension]
element_indexi = nl_domain.global_k.indexi
element_indexj = nl_domain.global_k.indexj
element_stiffness = nl_domain.global_k.stiffness
AUX = nl_domain.aux_variant
IDENTITY_TENSOR = nl_domain.identity_tensor
T_int = nl_domain.right_hand_item.T_int
RightHand = nl_domain.right_hand_item

fem_db = NLFEMDataBase()
MAT = fem_db.Material.Mat
MESH = fem_db.Mesh
CON = fem_db.SolveControl
LOAD_CASE = fem_db.LoadCase


def PressureLoadAndStiffnessAssembly():
    pass
