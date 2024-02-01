# !/usr/bin/env python3
# -*- coding: utf-8 -*-

from femdb.NLFEMDataBase import NLFEMDataBase
from GlobalEnum import *
from femdb.NLDomain import NLDomain
from element.ElementBase import ElementBaseClass

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

fem_db = NLFEMDataBase()
MAT = fem_db.Material.Mat
MESH = fem_db.Mesh

"""
Compute contribution (and extract relevant information for subsequent
assembly) of the constitutive term of the stiffness matrix.
EQUATION (9.35)
"""


def ConstitutiveMatrix(ele: ElementBaseClass,
                       i_gauss: int,
                       c,
                       JW):
    """

    @param ele:
    @param i_gauss:
    @param c:
    @param JW:
    @return:
    """
    for bnode in range(nl_domain.aux_variant.n_dofs_elem):
        for anode in range(nl_domain.aux_variant.n_dofs_elem):
            for j in range(dim):
                indexj = MESH.dof_nodes[j, ele.node_ids[bnode]]
                for i in range(dim):
                    indexi = MESH.dof_nodes[i, ele.node_ids[anode]]
                    sum_ = 0
                    for k in range(dim):
                        for l in range(dim):
                            """
                            Constitutive stiffness matrix contribution.
                            """
                            sum_ += (KINEMATICS.DN_Dx[i_gauss][anode] *
                                     c[i, k, j, l] *
                                     KINEMATICS.DN_Dx[i_gauss][bnode] *
                                     JW)

                    counter = nl_domain.global_k.counter
                    element_indexi[counter] = indexi
                    element_indexj[counter] = indexj
                    element_stiffness[counter] = sum_
                    counter += 1
