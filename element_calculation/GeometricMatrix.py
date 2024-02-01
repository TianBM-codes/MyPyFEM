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
AUX = nl_domain.aux_variant
dim = GlobalInfor[GlobalVariant.Dimension]
indexi = nl_domain.global_k.indexi
indexj = nl_domain.global_k.indexj
global_stiffness = nl_domain.global_k.stiffness

fem_db = NLFEMDataBase()
MAT = fem_db.Material.Mat
MESH = fem_db.Mesh

"""
Compute contribution (and extract relevant information for subsequent
assembly) of the geometric term of the stiffness matrix.
"""


def GeometricMatrix(ele: ElementBaseClass, DN_sigma_DN, JW):
    counter = nl_domain.global_k.counter
    for bnode in range(AUX.n_nodes_element):
        for anode in range(AUX.n_nodes_element):
            # Geometric stiffness matrix contribution.
            DNa_sigma_DNb = DN_sigma_DN[anode, bnode]
            # Index for row identification.
            indexi[counter:counter + dim] = AUX.dof_nodes[:, ele.node_ids[anode] - 1]
            # Index for column identification.
            indexj[counter:counter + dim] = AUX.dof_nodes[:, ele.node_ids[bnode] - 1]
            global_stiffness[counter:counter + dim] = DNa_sigma_DNb * JW
            counter += dim
