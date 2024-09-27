# !/usr/bin/env python3
# -*- coding: utf-8 -*-

from femdb.NLFEMDataBase import NLFEMDataBase
from utils.GlobalEnum import *
from femdb.ElementGroup import ElementGroup
import numpy as np
from scipy.sparse import coo_matrix

"""
Single instance mode, convenient for programming, Connect Database
"""
fem_db = NLFEMDataBase()
KINEMATICS = fem_db.kinematics
dim = GlobalInfor[GlobalVariant.Dimension]
element_indexi = fem_db.global_k.indexi
element_indexj = fem_db.global_k.indexj
GLOBAL_K = fem_db.global_k
T_int = fem_db.right_hand_item.T_int
RightHand = fem_db.right_hand_item

fem_db = NLFEMDataBase()
MAT = fem_db.Material.Mat
MESH = fem_db.Mesh
CON = fem_db.SolveControl
LOAD_CASE = fem_db.LoadCase


def PressureLoadAndStiffnessAssembly(grp: ElementGroup):
    """
    Update nodal forces and stiffness matrix due to external pressure
    boundary face (line) contributions.
    Javier Bonet P248
    @param grp:
    @return:
    """
    RightHand.nominal_pressure = np.zeros((MESH.n_dofs, 1), dtype=float)
    RightHand.R_pressure = np.zeros((MESH.n_dofs, 1), dtype=float)
    RightHand.K_pressure = coo_matrix((MESH.n_dofs, MESH.n_dofs), dtype=float)

    """
    Pre-allocation of memory for subsequent sparse assembly.
    Number of components of the vectors indexi, indexj and data.
    Initialise counter for storing sparse information into the tangent stiffness matrix.
    """
    n_components = grp.element_info.n_face_dofs_elem ** 2 * grp.quadrature.boundary_ngauss * LOAD_CASE.n_pressure_loads
    indexi = np.zeros(n_components, dtype=np.uint32)
    indexj = np.zeros(n_components, dtype=np.uint32)
    stiffness = np.zeros(n_components, dtype=float)
    counter = 0

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
        nominal_pressure = R_pressure_0 * LOAD_CASE.p_loads[ipressure]

        """
        Assemble boundary (NOMINAL pressure load) element force vector 
        contribution into global force vector. 
        """
        RightHand.nominal_pressure += nominal_pressure

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
    RightHand.residual += CON.dlamb * RightHand.nominal_pressure

    """
    Subtract (opposite to outward normal) boundary (CURRENT pressure load) 
    stiffness matrix into overall stiffness matrix.
    """
    GLOBAL_K.stiffness -= RightHand.K_pressure
