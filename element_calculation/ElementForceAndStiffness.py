# !/usr/bin/env python3
# -*- coding: utf-8 -*-

from femdb.NLFEMDataBase import NLFEMDataBase
from element.ElementBase import ElementBaseClass
from femdb.ElementGroup import ElementGroup
from femdb.Plasticity import PlasticDeformationState
from utils.GlobalEnum import *
from femdb.Material import *
import numpy as np

"""
Single instance mode, convenient for programming, Connect Database
"""


def ElementForceAndStiffness(xlocal, Xlocal, mat_id, Ve,
                             ele: ElementBaseClass,
                             grp: ElementGroup,
                             ele_idx: int):
    """
    #TODO: 所谓的Voigt格式在哪里体现
    Computes the element vector of global internal forces and the tangent
    stiffness matrix.
    @return:
    """
    fem_db = NLFEMDataBase()
    KINEMATICS = fem_db.kinematics
    dim = GetDomainDimension()
    global_k = fem_db.global_k
    element_indexi = global_k.indexi
    element_indexj = global_k.indexj
    element_stiffness = global_k.stiffness
    IDENTITY_TENSOR = fem_db.IdentityTensor
    MESH = fem_db.Mesh
    grp_ele_info = grp.element_info
    T_internal = np.zeros((grp_ele_info.n_dofs_elem, 1), dtype=float)
    KINEMATICS.ComputeGradients(xlocal, Xlocal, grp.interpolation.element_DN_chi)
    mat = fem_db.Material.Mat[mat_id]

    """
    Computes element mean dilatation kinematics, pressure and bulk modulus.
    """
    if isinstance(mat, HyperElasticPlasticInPrincipal):
        DN_x_mean = np.zeros((GetDomainDimension(), grp_ele_info.n_nodes_elem), dtype=float)
        ve = 0
        for igauss in range(grp_ele_info.ngauss):
            """
            Gauss contribution to the elemental deformed volume.
            Elemental averaged shape functions.
            """
            JW = KINEMATICS.Jx_chi[igauss] * grp.quadrature.element_w[igauss]
            ve += JW
            DN_x_mean += KINEMATICS.DN_Dx[:, :, igauss] * JW
        DN_x_mean /= ve
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
    if isinstance(mat, HyperElasticPlasticInPrincipal):
        PLAST_element = PlasticDeformationState()
        PLAST_element.epbar = grp.global_plasticity.epbar[:, ele_idx]
        PLAST_element.invCp = grp.global_plasticity.invCp[:, :, :, ele_idx]
    else:
        raise NoImplSuchMaterial(mat.GetName())

    from constitutive_laws.CauchyTypeSelection import CauchyTypeSelection
    from constitutive_laws.ElasticityModulus import ElasticityModulusSelection
    for igauss in range(grp_ele_info.ngauss):
        Cauchy, PLAST_gauss = CauchyTypeSelection(PLAST_element, igauss, mat)
        c = ElasticityModulusSelection(PLAST_element, PLAST_gauss, igauss, mat_id)

        """
        Add pressure contribution to stresses and elasticity tensor. 
        P251 (9.54) P252 (9.56) P166 (6.39)
        """
        Cauchy = Cauchy + press * IDENTITY_TENSOR.I
        c = c + press * (IDENTITY_TENSOR.c1 - IDENTITY_TENSOR.c2)

        """
        Compute numerical integration multipliers, Computes the thickness in the deformed 
        configuration for plane stress problems.
        """
        JW = KINEMATICS.Jx_chi[igauss] * grp.interpolation.quadrature.element_w[igauss]

        """
        Compute equivalent (internal) force vector.
        """
        T = np.matmul(Cauchy, KINEMATICS.DN_Dx[:, :, igauss])
        T_internal += T.reshape((T.size, 1), order='F') * JW

        """
        Compute contribution (and extract relevant information for subsequent
        assembly) of the constitutive term of the stiffness matrix.
        """
        from element_calculation.ConstitutiveMatrix import ConstitutiveMatrix
        ConstitutiveMatrix(grp, ele, igauss, c, JW)

        """
        Compute contribution (and extract relevant information for subsequent
        assembly) of the geometric term of the stiffness matrix.
        """
        DN_sigma_DN = np.matmul(KINEMATICS.DN_Dx[:, :, igauss].T, np.matmul(Cauchy, KINEMATICS.DN_Dx[:, :, igauss]))

        from element_calculation.GeometricMatrix import GeometricMatrix
        GeometricMatrix(grp, ele, DN_sigma_DN, JW)

    """
    Compute contribution (and extract relevant information for subsequent
    assembly) of the mean dilatation term (Kk) of the stiffness matrix.
    """
    if isinstance(mat, HyperElasticPlasticInPrincipal):
        for bnode in range(grp_ele_info.n_nodes_elem):
            for anode in range(grp_ele_info.n_nodes_elem):
                DN_x_meana_DN_x_meanb = np.outer(DN_x_mean[:, anode], DN_x_mean[:, bnode])

                indexi = MESH.dof_nodes[:, ele.search_node_ids[anode]]
                indexi = np.tile(indexi, (dim, 1))

                indexj = MESH.dof_nodes[:, ele.search_node_ids[bnode]]
                indexj = np.tile(indexj, (dim, 1))

                # Index for row identification.
                element_indexi[global_k.counter:global_k.counter + dim ** 2] = indexi.flatten('F')

                # Index for column identification.
                element_indexj[global_k.counter:global_k.counter + dim ** 2] = indexj.flatten('F')

                # Mean dilatation stiffness matrix contribution.
                element_stiffness[global_k.counter:global_k.counter + dim ** 2] = kappa_bar * ve * DN_x_meana_DN_x_meanb.flatten('F')

                global_k.counter += dim ** 2

    return T_internal, PLAST_element


def ElementForceAndStiffnessTruss(xlocal, Xlocal, mat_id,
                                  ele: ElementBaseClass,
                                  grp: ElementGroup,
                                  ele_idx: int):
    """

    @param xlocal:
    @param Xlocal:
    @param mat_id:
    @param ele:
    @param grp:
    @param ele_idx:
    @return:
    """
    fem_db = NLFEMDataBase()
    global_k = fem_db.global_k
    MESH = fem_db.Mesh
    mat = fem_db.Material.Mat[mat_id]
    area = mat.value_dict[MaterialKey.Area]
    E = mat.value_dict[MaterialKey.E]
    ty0 = mat.value_dict[MaterialKey.TauY]
    H = mat.value_dict[MaterialKey.Harden]
    dim = GetDomainDimension()
    if isinstance(mat, StretchBasedHyperelasticPlastic):
        PLAST_element = PlasticDeformationState()
        PLAST_element.epbar = grp.global_plasticity.epbar[ele_idx]
        PLAST_element.ep = grp.global_plasticity.ep[ele_idx]

    else:
        raise WrongMaterialType(mat.name)

    """
    Temporary variables
    """

    L = np.linalg.norm(Xlocal[:, 1] - Xlocal[:, 0])
    dx = xlocal[:, 1] - xlocal[:, 0]
    l = np.linalg.norm(dx)
    n = dx / l
    V = area * L
    lambda_ = l / L
    epsilon = np.log(lambda_)

    """
    Trial stage
    """
    epsilon_trial = epsilon - PLAST_element.ep
    tau_trial = E * epsilon_trial

    """
    Check yield criterion
    """
    f = abs(tau_trial) - (ty0 + H * PLAST_element.epbar)

    """
    Return mapping algorithm
    """
    if f > 0:
        Dgamma = f / (E + H)
        E_computational = E * H / (E + H)
    else:
        Dgamma = 0
        E_computational = E

    Dep = Dgamma * np.sign(tau_trial)
    tau = tau_trial - E * Dep
    PLAST_element.ep += Dep
    PLAST_element.epbar += Dgamma

    """
    Computation of the internal force vector
    """
    T = tau * V / l
    Tb = T * n
    T_internal = np.reshape(np.concatenate((-Tb, Tb), axis=None), (Tb.size*2, 1))

    """
    Computation of the stiffness matrix
    """
    k = (V / l ** 2) * (E_computational - 2 * tau)
    Kbb = np.array(k * np.outer(n, n) + (T / l) * np.eye(3), dtype=float)

    top_row = np.hstack([Kbb, -Kbb])
    bottom_row = np.hstack([-Kbb, Kbb])
    K_internal = np.vstack([top_row, bottom_row])

    """
    Sparse assembly of the stiffness matrix.
    """
    tmp = MESH.dof_nodes[:, ele.search_node_ids]
    element_indexj = np.tile(tmp.reshape((tmp.size,), order='F'), (tmp.size, 1))
    element_indexi = element_indexj.T
    n_dofs_elem = (grp.element_info.n_nodes_elem * dim) ** 2
    aux_vector = np.asarray(range(global_k.counter, global_k.counter + n_dofs_elem))
    global_k.indexi[aux_vector] = element_indexi.flatten(order='F')
    global_k.indexj[aux_vector] = element_indexj.flatten(order='F')
    global_k.stiffness[aux_vector] = K_internal.flatten(order='F')

    """
    Global index for sparse assembly.
    """
    global_k.counter += n_dofs_elem

    """
    Storage of stress for postprocessing purposes
    """
    # J = lambda_ ** (1 - 2 * mat.value_dict[MaterialKey.PoissonRatio])
    # Cauchy = tau / J

    return T_internal, PLAST_element
