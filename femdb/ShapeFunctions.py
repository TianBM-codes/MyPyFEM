#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from utils.CustomException import *


def shape_functions_library(Chi, e_type):
    """
    This function computes the shape functions and derivative of the shape
    functions for the specific element types considered in the code and for
    a given Gauss point.
    """
    chi = Chi[0]
    eta = Chi[1]

    if e_type == 'tria3':
        N = [1 - eta - chi, chi, eta]
        DN_chi = np.array([[-1, 1, 0],
                           [-1, 0, 1]])

    elif e_type == 'tria6':
        N_ = [(chi + eta - 0.5) * (2 * chi + 2 * eta - 2), -4 * chi * (chi + eta - 1), 2 * chi * (chi - 0.5), -4 * eta * (chi + eta - 1), 4 * chi * eta, 2 * eta * (eta - 0.5)]
        DN_chi_ = np.array([[4 * chi + 4 * eta - 3, 4 - 4 * eta - 8 * chi, 4 * chi - 1, -4 * eta, 4 * eta, 0],
                            [4 * chi + 4 * eta - 3, -4 * chi, 0, 4 - 8 * eta - 4 * chi, 4 * chi, 4 * eta - 1]])
        N = [N_[i] for i in [0, 1, 2, 4, 5, 3]]
        DN_chi = DN_chi_[:, [0, 1, 2, 4, 5, 3]]

    elif e_type == 'quad4':
        N_ = [((chi - 1) * (eta - 1)) / 4, -((chi + 1) * (eta - 1)) / 4, ((chi + 1) * (eta + 1)) / 4, -((chi - 1) * (eta + 1)) / 4]
        DN_chi_ = np.array([[eta / 4 - 0.25, 0.25 - eta / 4, -eta / 4 - 0.25, eta / 4 + 0.25],
                            [chi / 4 - 0.25, -chi / 4 - 0.25, chi / 4 + 0.25, 0.25 - chi / 4]])
        N = [N_[i] for i in [0, 1, 3, 2]]
        DN_chi = DN_chi_[:, [0, 1, 3, 2]]

    elif e_type == 'tetr4':
        iota = Chi[2]
        N = [1 - eta - iota - chi, chi, eta, iota]
        DN_chi = np.array([[-1, 1, 0, 0],
                           [-1, 0, 1, 0],
                           [-1, 0, 0, 1]])

    elif e_type == 'tetr10':
        iota = Chi[2]
        N = [(chi + eta + iota - 1 / 2) * (2 * chi + 2 * eta + 2 * iota - 2), -4 * chi * (chi + eta + iota - 1),
             2 * chi * (chi - 1 / 2), -4 * eta * (chi + eta + iota - 1), 4 * chi * eta, 2 * eta * (eta - 1 / 2),
             -4 * iota * (chi + eta + iota - 1), 4 * chi * iota, 4 * eta * iota, 2 * iota * (iota - 1 / 2)]
        DN_chi = np.array([[4 * chi + 4 * eta + 4 * iota - 3, 4 - 4 * eta - 4 * iota - 8 * chi, 4 * chi - 1,
                            -4 * eta, 4 * eta, 0, -4 * iota, 4 * iota, 0, 0],
                           [4 * chi + 4 * eta + 4 * iota - 3, -4 * chi, 0, 4 - 8 * eta - 4 * iota - 4 * chi,
                            4 * chi, 4 * eta - 1, -4 * iota, 0, 4 * iota, 0],
                           [4 * chi + 4 * eta + 4 * iota - 3, -4 * chi, 0, -4 * eta, 0, 0,
                            4 - 4 * eta - 8 * iota - 4 * chi, 4 * chi, 4 * eta, 4 * iota - 1]])

    elif e_type == 'hexa8':
        iota = Chi[2]
        N_ = np.array([
            -((chi - 1) * (eta - 1) * (iota - 1)) / 8,
            ((chi + 1) * (eta - 1) * (iota - 1)) / 8,
            ((chi - 1) * (eta + 1) * (iota - 1)) / 8,
            -((chi + 1) * (eta + 1) * (iota - 1)) / 8,
            ((chi - 1) * (eta - 1) * (iota + 1)) / 8,
            -((chi + 1) * (eta - 1) * (iota + 1)) / 8,
            -((chi - 1) * (eta + 1) * (iota + 1)) / 8,
            ((chi + 1) * (eta + 1) * (iota + 1)) / 8,
        ])

        # Derivatives of the shape functions
        DN_chi_ = np.array([
            [-(eta - 1) * (iota - 1) / 8, (eta - 1) * (iota - 1) / 8, (eta + 1) * (iota - 1) / 8, -(eta + 1) * (iota - 1) / 8,
             (eta - 1) * (iota + 1) / 8, -(eta - 1) * (iota + 1) / 8, -(eta + 1) * (iota + 1) / 8, (eta + 1) * (iota + 1) / 8],
            [-(chi - 1) * (iota - 1) / 8, (chi + 1) * (iota - 1) / 8, (chi - 1) * (iota - 1) / 8, -(chi + 1) * (iota - 1) / 8,
             (chi - 1) * (iota + 1) / 8, -(chi + 1) * (iota + 1) / 8, -(chi - 1) * (iota + 1) / 8, (chi + 1) * (iota + 1) / 8],
            [-(chi - 1) * (eta - 1) / 8, (chi + 1) * (eta - 1) / 8, (chi - 1) * (eta + 1) / 8, -(chi + 1) * (eta + 1) / 8,
             (chi - 1) * (eta - 1) / 8, -(chi + 1) * (eta - 1) / 8, -(chi - 1) * (eta + 1) / 8, (chi + 1) * (eta + 1) / 8]
        ])

        N = N_[[0, 1, 3, 2, 4, 5, 7, 6]]
        DN_chi = DN_chi_[:, [0, 1, 3, 2, 4, 5, 7, 6]]

    else:
        raise NoImplSuchShapeFunction(e_type)

    return N, DN_chi


def shape_functions_library_boundary(Chi, e_type):
    """
    Computes the shape functions and derivative of the shape
    functions for the specific element types considered in the code and for
    a given Gauss point at the edges (2D elements) or faces (3D elements) of
    those elements.
    """
    chi, eta, iota = Chi[0], Chi[1], Chi[2] if len(Chi) > 2 else None

    if e_type in ['tria3', 'quad4']:
        boundary_N = [0.5 * (1 - chi), 0.5 * (1 + chi)]
        boundary_DN_chi = 0.5 * np.array([-1, 1])

    elif e_type == 'tria6':
        N = [(chi + eta - 1 / 2) * (2 * chi + 2 * eta - 2), -4 * chi * (chi + eta - 1), 2 * chi * (chi - 1 / 2),
             -4 * eta * (chi + eta - 1), 4 * chi * eta, 2 * eta * (eta - 1 / 2)]
        DN_chi = np.array([[4 * chi + 4 * eta - 3, 4 - 4 * eta - 8 * chi, 4 * chi - 1,
                            -4 * eta, 4 * eta, 0],
                           [4 * chi + 4 * eta - 3, -4 * chi, 0, 4 - 8 * eta - 4 * chi,
                            4 * chi, 4 * eta - 1]])
        boundary_N = [N[0], N[1], N[2], N[4], N[5], N[3]]
        boundary_DN_chi = DN_chi[:, [0, 1, 2, 4, 5, 3]]

    elif e_type == 'tetr4':
        boundary_N, boundary_DN_chi = shape_functions_library(Chi, 'tria3')

    elif e_type == 'tetr10':
        boundary_N, boundary_DN_chi = shape_functions_library(Chi, 'tria6')

    elif e_type == 'hexa8':
        boundary_N, boundary_DN_chi = shape_functions_library(Chi, 'quad4')

    else:
        raise NoImplSuchShapeFunction

    return boundary_N, boundary_DN_chi
