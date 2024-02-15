# !/usr/bin/env python3
# -*- coding: utf-8 -*-


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
        self.dof_nodes = None
