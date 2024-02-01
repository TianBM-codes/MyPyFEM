# !/usr/bin/env python3
# -*- coding: utf-8 -*-

from femdb.ElementGroup import ElementGroup
from typing import Dict


class Mesh(object):
    """
    Mesh Information
    """

    def __init__(self):
        """
        n_nodes_elem: every elem nodes count
        """
        self.element_type = None
        self.nelem = -1
        self.n_dofs = -1
        self.n_nodes_elem = -1
        self.ele_grp_hash: Dict[int, ElementGroup] = {}
        self.dof_nodes = None
