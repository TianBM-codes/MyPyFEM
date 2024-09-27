# !/usr/bin/env python3
# -*- coding: utf-8 -*-

from femdb.NLFEMDataBase import NLFEMDataBase
from utils.GlobalEnum import *
from femdb.ElementGroup import ElementGroup
from element.ElementBase import ElementBaseClass

"""
Single instance mode, convenient for programming, Connect Database
"""
fem_db = NLFEMDataBase()
MAT = fem_db.Material.Mat
MESH = fem_db.Mesh
KINEMATICS = fem_db.kinematics
dim = GetDomainDimension()
element_indexi = fem_db.global_k.indexi
element_indexj = fem_db.global_k.indexj
element_stiffness = fem_db.global_k.stiffness


def ConstitutiveMatrix(grp: ElementGroup,
                       ele: ElementBaseClass,
                       i_gauss: int,
                       c,
                       JW):
    """
    Compute contribution (and extract relevant information for subsequent
    assembly) of the constitutive term of the stiffness matrix.
    EQUATION (9.35)  Indicial Form
    @param grp:
    @param ele:
    @param i_gauss:
    @param c:
    @param JW:
    @return:
    """
    for bnode in range(grp.element_info.n_nodes_elem):
        for anode in range(grp.element_info.n_nodes_elem):
            for j in range(dim):
                indexj = MESH.dof_nodes[j, ele.search_node_ids[bnode]]
                for i in range(dim):
                    indexi = MESH.dof_nodes[i, ele.search_node_ids[anode]]
                    sum_ = 0
                    for k in range(dim):
                        for l in range(dim):
                            """
                            Constitutive stiffness matrix contribution.
                            """
                            sum_ += (KINEMATICS.DN_Dx[:, :, i_gauss][k, anode] *
                                     c[i, k, j, l] *
                                     KINEMATICS.DN_Dx[:, :, i_gauss][l, bnode] *
                                     JW)

                    element_indexi[fem_db.global_k.counter] = indexi
                    element_indexj[fem_db.global_k.counter] = indexj
                    element_stiffness[fem_db.global_k.counter] = sum_
                    fem_db.global_k.counter += 1
