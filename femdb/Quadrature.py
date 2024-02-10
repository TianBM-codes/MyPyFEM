#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from utils.CustomException import *


class Quadrature:
    def __init__(self, e_type):
        self.element_chi:np.array = None
        self.element_w:np.array = None
        self.element_ngauss = None
        self.boundary_chi = None
        self.boundary_w = None
        self.boundary_ngauss = None
        self.set_quadrature_rules(e_type)
        self.set_edge_boundary_rules(e_type)

    def set_quadrature_rules(self, element_type):
        if element_type == 'tria3':
            self.element_chi = np.array([[1 / 3, 1 / 3]])
            self.element_w = np.array([0.5])

        elif element_type == 'quad4':
            self.element_chi = np.array([
                [-0.577350269189626, -0.577350269189626],
                [0.577350269189626, -0.577350269189626],
                [0.577350269189626, 0.577350269189626],
                [-0.577350269189626, 0.577350269189626]
            ])
            self.element_w = np.array([1, 1, 1, 1])

        elif element_type == 'tria6':
            self.element_chi = np.array([[0.5, 0], [0.5, 0.5], [0, 0.5]])
            self.element_w = np.array([1 / 6, 1 / 6, 1 / 6])

        elif element_type == 'tetr4':
            self.element_chi = np.array([[0.25, 0.25, 0.25]])
            self.element_w = np.array([1 / 6])

        elif element_type == 'tetr10':
            self.element_chi = np.array([
                [0.25, 0.25, 0.25], [1 / 6, 1 / 6, 1 / 6], [0.5, 1 / 6, 1 / 6],
                [1 / 6, 0.5, 1 / 6], [1 / 6, 1 / 6, 0.5]
            ])
            self.element_w = np.array([-0.1333333333333333, 0.075, 0.075, 0.075, 0.075])

        elif element_type == 'hexa8':
            self.element_chi = np.array([
                [-0.577350269189626, -0.577350269189626, -0.577350269189626],
                [0.577350269189626, -0.577350269189626, -0.577350269189626],
                [0.577350269189626, 0.577350269189626, -0.577350269189626],
                [-0.577350269189626, 0.577350269189626, -0.577350269189626],
                [-0.577350269189626, -0.577350269189626, 0.577350269189626],
                [0.577350269189626, -0.577350269189626, 0.577350269189626],
                [0.577350269189626, 0.577350269189626, 0.577350269189626],
                [-0.577350269189626, 0.577350269189626, 0.577350269189626]
            ])
            self.element_w = np.array([1, 1, 1, 1, 1, 1, 1, 1])

        # Set the number of integration points.
        self.element_ngauss = len(self.element_w)

    def set_edge_boundary_rules(self, e_type):
        """
        Obtain quadrature (numerical integration) rules for edge/surface elements.
        """
        if e_type in ['quad4', 'tria3']:
            self.boundary_chi = np.array([[-0.577350269189626], [0.577350269189626]])
            self.boundary_w = np.array([1, 1])
        elif e_type == 'tria6':
            self.boundary_chi = np.array([[0.774596669241483], [0], [-0.774596669241483]])
            self.boundary_w = np.array([0.555555555555554, 0.888888888888889, 0.555555555555554])
        elif e_type == 'tetr4':
            self.boundary_chi = np.array([[1 / 3, 1 / 3]])
            self.boundary_w = np.array([0.5])
        elif e_type == 'tetr10':
            self.boundary_chi = np.array([[0.5, 0], [0.5, 0.5], [0, 0.5]])
            self.boundary_w = np.array([1 / 6, 1 / 6, 1 / 6])
        elif e_type == 'hexa8':
            self.boundary_chi = np.array([
                [-0.577350269189626, -0.577350269189626],
                [0.577350269189626, -0.577350269189626],
                [0.577350269189626, 0.577350269189626],
                [-0.577350269189626, 0.577350269189626]
            ])
            self.boundary_w = np.array([1, 1, 1, 1])

        # Set the number of integration points for boundary.
        self.boundary_ngauss = len(self.boundary_w)
