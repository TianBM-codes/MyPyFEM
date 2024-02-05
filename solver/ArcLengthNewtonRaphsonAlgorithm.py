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


def ArcLengthNewtonRaphsonAlgorithm():
    afail = CON.afail
    arcln = CON.arcln
    while CON.xlmax < CON.xlmax and CON.incrm < CON.nincr:
        CON.incrm += 1
        """
        Update the load factor. The radius is adjusted to achieve target 
        iterations per increment arbitrarily dampened by a factor of 0.7.
        """
        if not afail and not CON.fracl:
            CON.arcln = CON.arcln * (CON.itarget / CON.iterold) ** 0.7
            arcln = CON.arcln

        """
        Update nodal forces (excluding pressure) and gravity. 
        """
        RightHand.residual = (RightHand.residual -
                              CON.dlamb * RightHand.nominal_external_load)
        RightHand.external_load = (RightHand.external_load +
                                   CON.dlamb * RightHand.nominal_external_load)

        """
        Update nodal forces and stiffness matrix due to external pressure 
        boundary face (line) contributions.
        """
        if LOAD_CASE.n_pressure_loads > 0:
            from global_assembly.PressureLoadAndStiffnessAssembly import PressureLoadAndStiffnessAssembly
            PressureLoadAndStiffnessAssembly()