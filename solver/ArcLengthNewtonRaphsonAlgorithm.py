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
RightHand = nl_domain.right_hand_item


def linear_solver(K, F, fixed_dof):
    """
    Solve the linear system Ku = F for iterative displacements,
    excluding the degrees of freedom with Dirichlet boundary conditions.

    Parameters:
    K -- Stiffness matrix.
    F -- Force vector.
    fixed_dof -- Indices of degrees of freedom with Dirichlet boundary conditions.

    Returns:
    u -- Solution vector for displacements.
    rtu -- Dot product of the solution vector and the force vector.
    """
    # 将 fixed_dof 从 1 基索引转换为 0 基索引（如果从 MATLAB 转换而来）
    fixed_dof = np.array(fixed_dof) - 1

    # 从 K 和 F 中移除 Dirichlet 边界条件的自由度
    K = np.delete(K, fixed_dof, axis=0)
    K = np.delete(K, fixed_dof, axis=1)
    F = np.delete(F, fixed_dof, axis=0)

    # 解线性系统 Ku = F
    u = np.linalg.solve(K, F)

    # 计算 u 和 F 的点积
    rtu = np.dot(u.T, F)

    return u, rtu


def arclen(displ, dispf):
    fem_db = NLFEMDataBase()
    CON = fem_db.SolveControl
    CON.Arclen.afail = 0

    # Initialisation
    urdx = 0
    dxdx = 0

    # First obtains all the dot products required.
    ufuf = np.dot(dispf.T, dispf)
    urur = np.dot(displ.T, displ)
    ufur = np.dot(dispf.T, displ)

    # For the first iteration set up
    if CON.Arclen.xincr is None:
        CON.Arclen.xincr = np.zeros(displ.shape)
    ufdx = np.dot(dispf.T, CON.Arclen.xincr)
    if CON.niter == 1:
        CON.Arclen.xincr = np.zeros(displ.shape)
        displ = np.zeros(displ.shape)
    else:
        urdx = np.dot(displ.T, CON.Arclen.xincr)
        dxdx = np.dot(CON.Arclen.xincr.T, CON.Arclen.xincr)

    # Set up arcln for first increment, first iteration.
    if not CON.Arclen.farcl and CON.incrm * CON.niter == 1:
        CON.Arclen.arcln = np.sqrt(urur)

    # Finding first iteration gamma and the sign of gamma.
    if CON.niter == 1:
        gamma = np.abs(CON.Arclen.arcln) / np.sqrt(ufuf)
        if ufdx != 0:
            gamma = gamma * np.sign(ufdx)
    else:
        # Solves quadratic equation.
        a1 = ufuf
        a2 = 2 * (ufdx + ufur)
        a3 = dxdx + 2 * urdx + urur - CON.Arclen.arcln ** 2
        discr = a2 ** 2 - 4 * a1 * a3
        if discr < 0:
            CON.Arclen.afail = 1
            return displ, CON
        discr = np.sqrt(discr)
        if a2 < 0:
            discr = -discr
        discr = -(a2 + discr) / 2
        gamma1 = discr / a1
        gamma2 = a3 / discr

        gamma = gamma1
        cos1 = urdx + gamma1 * ufdx
        cos2 = urdx + gamma2 * ufdx
        if cos2 > cos1:
            gamma = gamma2

    displ += gamma * dispf
    CON.Arclen.xincr += displ

    CON.dlamb += gamma
    CON.xlamb += gamma

    if np.isnan(CON.xlamb):
        CON.Arclen.afail = 1

    return displ

def CheckResidualNorm():
    """
    Check for equilibrium convergence.
    """
    """
    Obtain the reactions
    """
    fem_db = NLFEMDataBase()
    BC = fem_db.BC
    RightHand.reactions = RightHand.residual[BC.fixed_dof] + \
                          RightHand.external_load[BC.fixed_dof]
    CON = fem_db.SolveControl

    r_norm = np.dot(RightHand.residual[BC.free_dof].T, RightHand.residual[BC.free_dof])

    # 计算外部载荷范数
    f_norm = np.dot((RightHand.nominal_external_load[BC.free_dof] - RightHand.nominal_pressure[BC.free_dof]).T,
                   (RightHand.nominal_external_load[BC.free_dof] - RightHand.nominal_pressure[BC.free_dof]))
    f_norm = f_norm * CON.xlamb ** 2

    # 计算反应力范数
    e_norm = np.dot(RightHand.Reactions.T, RightHand.Reactions)

    # 计算最终范数
    r_norm = np.sqrt(r_norm / (f_norm + e_norm))
    return r_norm

def ArcLengthNewtonRaphsonAlgorithm():
    from global_assembly.ResidualAndStiffnessAssembly import ResidualAndStiffnessAssembly
    from global_assembly.PressureLoadAndStiffnessAssembly import PressureLoadAndStiffnessAssembly
    fem_db = NLFEMDataBase()
    CON = fem_db.SolveControl
    BC = fem_db.BC
    GEOM = fem_db.Geom
    afail = CON.Arclen.afail
    arcln = CON.Arclen.arcln
    LOAD_CASE = fem_db.LoadCase

    while CON.xlmax < CON.xlmax and CON.incrm < CON.nincr:
        CON.incrm += 1
        """
        Update the load factor. The radius is adjusted to achieve target 
        iterations per increment arbitrarily dampened by a factor of 0.7.
        """
        if not afail and not CON.fracl:
            CON.arcln = CON.arcln * (CON.itarget / CON.iterold) ** 0.7
            arcln = CON.arcln

        """
        Update nodal forces (excluding pressure) and gravity. 
        """
        RightHand.residual = (RightHand.residual -
                              CON.dlamb * RightHand.nominal_external_load)
        RightHand.external_load = (RightHand.external_load +
                                   CON.dlamb * RightHand.nominal_external_load)

        """
        Update nodal forces and stiffness matrix due to external pressure 
        boundary face (line) contributions.
        """
        if LOAD_CASE.n_pressure_loads > 0:
            from global_assembly.PressureLoadAndStiffnessAssembly import PressureLoadAndStiffnessAssembly
            PressureLoadAndStiffnessAssembly()

        """
        For the case of prescribed geometry update coordinates.
        Recompute equivalent nodal forces and assembles residual force, 
        excluding pressure contributions.
        Recompute and assembles tangent stiffness matrix components, 
        excluding pressure contributions.
        """
        if LOAD_CASE.n_prescribed_displacements > 0:
            pass

        """
        Newton-Raphson iteration.
        """
        CON.niter = 0
        rnorm = 2 * CON.cnorm
        while rnorm > CON.cnorm and CON.niter < CON.miter:
            CON.niter += 1
            # Solve for iterative displacements.
            displ, _ = linear_solver(nl_domain.global_k.stiffness, -RightHand.residual, BC.fixdof)
            # Solve for displacements (for a nominal load).
            dispf, _ = linear_solver(nl_domain.global_k.stiffness,
                                  RightHand.nominal_external_load - RightHand.nominal_pressure,
                                  BC.fixdof)
            displ = arclen(displ, dispf)
            # Update coordinates.
            GEOM.x[BC.free_dof] += displ
            ResidualAndStiffnessAssembly()

            # Update nodal forces due to pressure.
            if LOAD_CASE.n_pressure_loads:
                PressureLoadAndStiffnessAssembly()

            # Check for equilibrium convergence.
            rnorm, GLOBAL = CheckResidualNorm()
            # Break iteration before residual gets unrealistic (e.g., NaN).
            if abs(rnorm) > 1e7 or np.isnan(rnorm):
                CON.ARCLEN.afail = 1
                break
            # Update value of the old iteration number.
            CON.ARCLEN.iter_old = CON.niter

        # If convergence not achieved opt to restart from previous results or terminate.
        afail = CON.ARCLEN.afail
        if CON.niter >= CON.miter or CON.ARCLEN.afail:
            terminate = input('Solution not converged. Do you want to terminate the program (y/n) ?: ')
            if terminate.lower() == 'n':
                # Restart from previous step by decreasing the fixed arclength radius to half its initially fixed value.
                print('Restart from previous step by decreasing the fixed arclength radius to half its initially fixed value.')
                CON.ARCLEN.arcln /= 2
            else:
                afail = 1
        else:
            # If convergence was achieved, update and save internal variables.
            # PLAST = save_output(updated_PLAST, PRO, FEM, GEOM, QUADRATURE, BC, MAT, LOAD, CON, CONSTANT, GLOBAL, PLAST, KINEMATICS)
            pass

        # If a failure occurs in the arc-length function, reduce the value of the radius in the arc-length to a 10th of the previous value and automatically restart.
        if afail:
            CON.ARCLEN.arcln /= 10
