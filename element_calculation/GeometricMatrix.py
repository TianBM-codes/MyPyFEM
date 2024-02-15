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
KINEMATICS = fem_db.kinematics
dim = GetDomainDimension()
indexi = fem_db.global_k.indexi
indexj = fem_db.global_k.indexj
global_stiffness = fem_db.global_k.stiffness

fem_db = NLFEMDataBase()
MAT = fem_db.Material.Mat
MESH = fem_db.Mesh

"""
Compute contribution (and extract relevant information for subsequent
assembly) of the geometric term of the stiffness matrix.
"""


def GeometricMatrix(grp:ElementGroup, ele: ElementBaseClass, DN_sigma_DN, JW):
    counter = fem_db.global_k.counter
    for bnode in range(grp.element_info.n_nodes_elem):
        for anode in range(grp.element_info.n_nodes_elem):
            # Geometric stiffness matrix contribution.
            DNa_sigma_DNb = DN_sigma_DN[anode, bnode]
            # Index for row identification.
            indexi[counter:counter + dim, 0] = MESH.dof_nodes[:, ele.search_node_ids[anode] - 1]
            # Index for column identification.
            indexj[counter:counter + dim, 0] = MESH.dof_nodes[:, ele.search_node_ids[bnode] - 1]
            global_stiffness[counter:counter + dim] = DNa_sigma_DNb * JW
            counter += dim
