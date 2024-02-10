# !/usr/bin/env python3
# -*- coding: utf-8 -*-

from femdb.NLFEMDataBase import NLFEMDataBase
from element.ElementBase import ElementBaseClass
from utils.GlobalEnum import *
from femdb.NLDomain import NLDomain
from femdb.Material import *
import numpy as np

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
AUX = nl_domain.aux_variant
IDENTITY_TENSOR = nl_domain.identity_tensor
T_int = nl_domain.right_hand_item.T_int

fem_db = NLFEMDataBase()
MAT = fem_db.Material.Mat
MESH = fem_db.Mesh


def ElementForceAndStiffness(xlocal, Xlocal, mat_id, Ve, ele: ElementBaseClass):
    """
    Computes the element vector of global internal forces and the tangent
    stiffness matrix.
    @return:
    """
    T_internal = np.zeros((AUX.n_dofs_elem, 1))
    KINEMATICS.ComputeGradients(xlocal, Xlocal, AUX.DN_Dchi)
    mat = fem_db.Material.Mat[mat_id]

    """
    Computes element mean dilatation kinematics, pressure and bulk modulus.
    """
    if isinstance(mat, HyperElasticPlasticInPrincipal):
        DN_x_mean = np.zeros(GlobalInfor[GlobalVariant.Dimension], AUX.n_nodes_element)
        ve = 0
        for ii in range(AUX.ngauss):
            JW = KINEMATICS.Jx_chi[ii] * AUX.weight[ii]
            """
            Gauss contribution to the elemental deformed volume.Elemental averaged shape functions.
            """
            ve += JW
            DN_x_mean += AUX.DN_Dchi[:, :, ii] * JW
        DN_x_mean = DN_x_mean / ve
        Jbar = ve / Ve

        """
        Computes element mean dilatation pressure and bulk modulus.
        """
        if isinstance(mat, HyperElasticPlasticInPrincipal):
            kappa = mat.value_dict[MaterialKey.Kappa]
            press = kappa * np.log(Jbar) / Jbar
            kappa_bar = kappa / Jbar - press
        else:
            kappa_bar = 0
            press = 0

    else:
        raise NoImplSuchMaterial(mat.GetName())

    """
    Gauss quadrature integration loop. Extract kinematics at the particular Gauss point.
    Obtain stresses (for incompressible or nearly incompressible, 
    only deviatoric component) and internal variables in plasticity.
    Obtain elasticity tensor (for incompressible or nearly incompressible, 
    only deviatoric component).
    """
    from constitutive_laws.CauchyTypeSelection import CauchyTypeSelection
    from constitutive_laws.ElasticityModulus import ElasticityModulusSelection
    for jj in range(AUX.ngauss):
        Cauchy = CauchyTypeSelection(jj, mat)
        c = ElasticityModulusSelection(mat_id)
        """
        Add pressure contribution to stresses and elasticity tensor.
        """
        Cauchy = Cauchy + press * IDENTITY_TENSOR.I
        c = c + press * (IDENTITY_TENSOR.c1 - IDENTITY_TENSOR.c2)
        """
        Compute numerical integration multipliers, Computes the thickness in the deformed 
        configuration for plane stress problems.
        """
        JW = KINEMATICS.Jx_chi[jj] * AUX.weight[jj]

        """
        Compute equivalent (internal) force vector.
        """
        T = Cauchy * KINEMATICS.DN_Dx[jj]
        T_internal += T.T.reshape(T.size, 1)

        """
        Compute contribution (and extract relevant information for subsequent
        assembly) of the constitutive term of the stiffness matrix.
        """
        from element_calculation.ConstitutiveMatrix import ConstitutiveMatrix
        ConstitutiveMatrix(ele, jj, c, JW)

        """
        Compute contribution (and extract relevant information for subsequent
        assembly) of the geometric term of the stiffness matrix.
        """
        DN_sigma_DN = np.matmul(KINEMATICS.DN_Dx[jj].T, np.matmul(Cauchy, KINEMATICS.DN_Dx[jj]))

        from element_calculation.GeometricMatrix import GeometricMatrix
        GeometricMatrix(ele, DN_sigma_DN, JW)

    """
    Compute contribution (and extract relevant information for subsequent
    assembly) of the mean dilatation term (Kk) of the stiffness matrix.
    """
    counter = nl_domain.global_k.counter
    if isinstance(mat, HyperElasticPlasticInPrincipal):
        for bnode in range(AUX.n_nodes_elem):
            for anode in range(AUX.n_dofs_elem):
                DN_x_meana_DN_x_meanb = np.outer(DN_x_mean[:, anode], DN_x_mean[:, bnode])

                indexi = MESH.dof_nodes[:, ele.node_ids[anode] - 1]
                indexi = np.tile(indexi, (1, dim))

                indexj = MESH.dof_nodes[:, ele.node_ids[bnode] - 1]
                indexj = np.tile(indexj, (1, dim))
                indexj = indexj.T

                # Index for row identification.
                element_indexi[counter:counter + dim ** 2] = indexi.flatten()

                # Index for column identification.
                element_indexj[counter:counter + dim ** 2] = indexj.flatten()

                # Mean dilatation stiffness matrix contribution.
                element_stiffness[counter:counter + dim ** 2] = kappa_bar * ve * DN_x_meana_DN_x_meanb.flatten()

                counter += dim ** 2
    """
    Utility function for the assembly of element force vectors into the 
    global force vector. 
    """
    global_dofs = MESH.dof_nodes[:, ele.node_ids]
    T_int[global_dofs, 1] += T_internal


