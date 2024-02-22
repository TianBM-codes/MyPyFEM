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
dim = GetDomainDimension()
global_k = fem_db.global_k
indexi = global_k.indexi
indexj = global_k.indexj
global_stiffness = fem_db.global_k.stiffness
MESH = fem_db.Mesh

"""
Compute contribution (and extract relevant information for subsequent
assembly) of the geometric term of the stiffness matrix.
"""


def GeometricMatrix(grp: ElementGroup, ele: ElementBaseClass, DN_sigma_DN, JW):
    for bnode in range(grp.element_info.n_nodes_elem):
        for anode in range(grp.element_info.n_nodes_elem):
            # Geometric stiffness matrix contribution.
            DNa_sigma_DNb = DN_sigma_DN[anode, bnode]
            # Index for row identification.
            indexi[global_k.counter:global_k.counter + dim] = MESH.dof_nodes[:, ele.search_node_ids[anode]]
            # Index for column identification.
            indexj[global_k.counter:global_k.counter + dim] = MESH.dof_nodes[:, ele.search_node_ids[bnode]]
            global_stiffness[global_k.counter:global_k.counter + dim] = DNa_sigma_DNb * JW
            global_k.counter += dim
