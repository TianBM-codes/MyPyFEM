# !/usr/bin/env python3
# -*- coding: utf-8 -*-

from femdb.NLFEMDataBase import NLFEMDataBase
from utils.GlobalEnum import *
from femdb.NLDomain import NLDomain
import numpy as np
from scipy.sparse import coo_matrix

"""
Single instance mode, convenient for programming, Connect Database
"""
nl_domain = NLDomain()
PLAST = nl_domain.plastics
KINEMATICS = nl_domain.kinematics
dim = GlobalInfor[GlobalVariant.Dimension]
element_indexi = nl_domain.global_k.indexi
element_indexj = nl_domain.global_k.indexj
GLOBAL_K = nl_domain.global_k
AUX = nl_domain.aux_variant
T_int = nl_domain.right_hand_item.T_int
RightHand = nl_domain.right_hand_item

fem_db = NLFEMDataBase()
MAT = fem_db.Material.Mat
MESH = fem_db.Mesh
CON = fem_db.SolveControl
LOAD_CASE = fem_db.LoadCase


def PressureLoadAndStiffnessAssembly():
    """
    Update nodal forces and stiffness matrix due to external pressure
    boundary face (line) contributions.
    """
    RightHand.nominal_press = np.zeros((MESH.n_dofs, 1))
    RightHand.R_pressure = np.zeros((MESH.n_dofs, 1))
    RightHand.K_pressure = coo_matrix((MESH.n_dofs, MESH.n_dofs))

    """
    Pre-allocation of memory for subsequent sparse assembly.
    Number of components of the vectors indexi, indexj and data.
    Initialise counter for storing sparse information into the tangent stiffness matrix.
    """
    n_components = AUX.n_face_dofs_elem ** 2 * AUX.ngauss * LOAD_CASE.n_pressure_loads
    indexi = np.zeros(n_components)
    indexj = np.zeros(n_components)
    stiffness = np.zeros(n_components)
    counter = 1

    """
    Loop over all boundary (pressure load) elements.
    """
    for ipressure in range(LOAD_CASE.n_pressure_loads):
        """
        Intermediate variables associated to a particular element (ipressure).
        """
        element_id = LOAD_CASE.p_loads[ipressure].ele_id
        global_nodes = LOAD_CASE.p_loads[ipressure].face_node
        cur_group = None
        for _, grp in fem_db.ElementGroupHash.items():
            if grp.IsElementInGroup(element_id):
                cur_group = grp
                break

        """
        Compute boundary (UNIT pressure load) force vector and stiffness matrix 
        contribution for a boundary (UNIT pressure load) element.
        """
        counter0 = counter
        from element_calculation.PressureElementLoadAndStiffness import PressureElementLoadAndStiffness
        R_pressure_0 = PressureElementLoadAndStiffness(cur_group, global_nodes,
                                                       indexi, indexj, stiffness, counter)

        """
        Compute boundary (NOMINAL pressure load) force vector contribution 
        for a boundary (NOMINAL pressure load) element.
        """
        nominal_press = R_pressure_0 * LOAD_CASE.p_loads[ipressure]

        """
        Assemble boundary (NOMINAL pressure load) element force vector 
        contribution into global force vector. 
        """
        RightHand.nominal_press += nominal_press

        """
        Compute boundary (CURRENT pressure load) element stiffness matrix 
        contribution.
        """
        stiffness[counter0:counter - 1] = LOAD_CASE.p_loads[ipressure] * CON.xlmax * stiffness[counter0:counter - 1]

    """
    Assemble boundary (CURRENT pressure load) stiffness matrix.
    """
    RightHand.K_pressure = coo_matrix((stiffness, (indexi, indexj)), shape=(MESH.n_dofs, MESH.n_dofs))

    """
    Add boundary (CURRENT pressure load) global force vector contribution 
    into Residual.
    """
    RightHand.residual += CON.dlamb * RightHand.nominal_press

    """
    Subtract (opposite to outward normal) boundary (CURRENT pressure load) 
    stiffness matrix into overall stiffness matrix.
    """
    GLOBAL_K.stiffness -= RightHand.K_pressure
