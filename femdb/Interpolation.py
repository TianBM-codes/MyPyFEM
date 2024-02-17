#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from femdb.Quadrature import Quadrature
from femdb.ShapeFunctions import shape_functions_library_boundary, shape_functions_library
from femdb.ElementGroup import ElementInfo
from utils.GlobalEnum import GetDomainDimension


class Interpolation:
    """
    Compute the shape functions and their derivatives
    (with respect to the iso parametric domain) inside the element and on
    faces (3D) or edges (2D).
    """

    def __init__(self, e_type, ele_info: ElementInfo):
        dim = GetDomainDimension()
        self.ele_info = ele_info
        self.e_type = e_type
        self.element_N = np.zeros((ele_info.n_nodes_elem, ele_info.ngauss))
        self.element_DN_chi = np.zeros((dim, ele_info.n_nodes_elem, ele_info.ngauss))
        self.boundary_N = np.zeros((ele_info.n_face_nodes_elem, ele_info.boundary_ngauss))
        self.boundary_DN_chi = np.zeros((dim - 1, ele_info.n_face_nodes_elem, ele_info.boundary_ngauss))
        self.quadrature = Quadrature(e_type)
        self.init_variant()

    def init_variant(self):
        for ii in range(self.ele_info.ngauss):
            N, DN_chi = shape_functions_library(self.quadrature.element_chi[ii], self.e_type)
            self.element_N[:, ii] = N
            self.element_DN_chi[:, :, ii] = DN_chi

        for ii in range(self.ele_info.boundary_ngauss):
            N, DN_chi = shape_functions_library_boundary(self.quadrature.boundary_chi[ii], self.e_type)
            self.boundary_N[:, ii] = N
            self.boundary_DN_chi[:, :, ii] = DN_chi
