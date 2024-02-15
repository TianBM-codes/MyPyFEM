# !/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys

from femdb.NLFEMDataBase import NLFEMDataBase
from femdb.ElementGroup import ElementGroup
from utils.GlobalEnum import *
from femdb.NLDomain import NLDomain
import numpy as np

"""
Single instance mode, convenient for programming, Connect Database
"""
nl_domain = NLDomain()
dim = GlobalInfor[GlobalVariant.Dimension]

fem_db = NLFEMDataBase()
MAT = fem_db.Material.Mat
MESH = fem_db.Mesh
CON = fem_db.SolveControl
LOAD_CASE = fem_db.LoadCase


def Levi_civita_contraction_vector(vector, dim_p):
    if dim_p == 2:
        matrix = np.array([[0, vector[2]], [-vector[2], 0]])
    elif dim_p == 3:
        matrix = np.array([
            [0, vector[2], -vector[1]],
            [-vector[2], 0, vector[0]],
            [vector[1], -vector[0], 0]
        ])
    else:
        mlogger.fatal(f"Wrong Dim:{dim_p}")
        sys.exit(7)
    return matrix


def PressureElementLoadAndStiffness(grp: ElementGroup, global_nodes,
                                    ele_indexi, ele_indexj, ele_stiffness, counter):
    R_pressure = np.zeros((grp.element_info.n_face_dofs_elem, 1))
    for ii in range(grp.quadrature.boundary_ngauss):
        N = grp.interpolation.boundary_N[:, ii]
        DN_chi = grp.interpolation.boundary_DN_chi[:, :, ii]
        xlocal_boundary = fem_db.Geom.node_list[global_nodes].coord

        """
        Determine the outward unit normal vector at the Gauss point.
        """
        if dim == 2:
            # For 2D elements, calculate the tangent vector and use it to compute the normal vector.
            dx_eta = np.dot(xlocal_boundary, DN_chi.T)
            dx_eta = np.append(dx_eta, 0)  # Append 0 to make it a 3D vector for cross product.
            k = np.array([0, 0, 1])
            normal_vector = np.cross(dx_eta, k)
            normal_vector = normal_vector[:-1]  # Remove the last component to make it a 2D vector.
        elif dim == 3:
            # For 3D elements, directly calculate the normal vector using the cross product of two tangent vectors.
            dx_chi = np.dot(xlocal_boundary, DN_chi.T)
            normal_vector = np.cross(dx_chi[:, 0], dx_chi[:, 1])
        else:
            mlogger.fatal(f"Wrong Dim:{dim}")
            sys.exit(7)

        W = grp.quadrature.boundary_w[ii]

        """
        Compute boundary (UNIT pressure load) force vector contribution.
        """
        T = np.matmul(normal_vector, grp.interpolation.boundary_N[:, ii].T)
        R_pressure += T * W

        """
        Compute boundary (UNIT pressure load) stiffness matrix contribution.
        Compute the stiffness matrix pressure load vector. We compute 
        here a matrix form of the expressions:
        -For 2D problems:
        0.5*k*(DN_etaa*Nb - Na*DN_etab).
        -For 3D problems:
        0.5*(dx_chi*(DN_etaa*Nb-Na*DN_etab)-dx_eta*(DN_chia*Nb-Na*DN_chib). 
        """
        for bnode in range(grp.element_info.n_face_dofs_elem):
            for anode in range(grp.element_info.n_face_dofs_elem):
                if dim == 2:
                    DN_eta = DN_chi
                    k = np.array([0, 0, 1])
                    scalar = DN_eta[:, anode] * N[bnode] - N[anode] * DN_eta[:, bnode]
                    # Obtain final expression for load vector T.
                    kpab = 0.5 * k * scalar

                elif dim == 3:
                    DN_chi = DN_chi[0, :]
                    DN_eta = DN_chi[1, :]
                    dx_chi = np.dot(xlocal_boundary, DN_chi.T)
                    dx_eta = np.dot(xlocal_boundary, DN_eta.T)
                    # Obtain final expression for load vector T.
                    kpab = np.outer(dx_chi, N[anode] * DN_eta[bnode]) - np.outer(dx_eta, N[anode] * DN_chi[bnode])

                else:
                    mlogger.fatal(f"Wrong Dim:{dim}")
                    sys.exit(7)
                Levi_civita_kpab = Levi_civita_contraction_vector(kpab, dim)

                """
                Compute and store indices for subsequent sparse assembly.
                """
                indexi = MESH.dof_nodes[:, global_nodes[anode]]
                indexi = np.tile(indexi, (1, dim))
                indexj = MESH.dof_nodes[:, global_nodes[bnode]]
                indexj = np.tile(indexj, (1, dim)).T
                ele_indexi[counter:counter + dim ** 2 - 1] = indexi
                ele_indexj[counter:counter + dim ** 2 - 1] = indexj
                ele_stiffness[counter:counter + dim ** 2 - 1] = Levi_civita_kpab * W
                counter += dim ** 2

    return R_pressure
