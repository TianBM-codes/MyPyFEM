# !/usr/bin/env python3
# -*- coding: utf-8 -*-

from femdb.NLFEMDataBase import NLFEMDataBase
from utils.GlobalEnum import *
from scipy.sparse.linalg import spsolve
import numpy as np

"""
Single instance mode, convenient for programming, Connect Database
"""
fem_db = NLFEMDataBase()
KINEMATICS = fem_db.kinematics
dim = GlobalInfor[GlobalVariant.Dimension]
element_indexi = fem_db.global_k.indexi
element_indexj = fem_db.global_k.indexj
element_stiffness = fem_db.global_k.stiffness
IDENTITY_TENSOR = fem_db.identity_tensor
T_int = fem_db.right_hand_item.T_int
RightHand = fem_db.right_hand_item


def linear_solver(K, F, fixed_dof):
    """
    Solve for iterative displacements.
    # TODO: 需要注意的是，这些操作虽然有效，但可能不是最高效的，尤其是对于非常大的稀疏矩阵。每次这样的操作都可能涉及到底层数据的复制。如果可能的话，最好在稀疏矩阵最终形成之前，或者在转换为稠密矩阵后进行行列的删除操作。
    """
    fixed_dof = np.array(fixed_dof)

    """
    Remove from K and Residual the rows (and columns for K) associated to
    degrees of freedom with Dirichlet boundary conditions.
    """
    all_rows = np.arange(K.shape[0])
    rows_to_keep = np.isin(all_rows, fixed_dof, invert=True)
    K = K[rows_to_keep]
    K_csc = K.tocsc()
    all_cols = np.arange(K_csc.shape[1])
    cols_to_keep = np.isin(all_cols, fixed_dof, invert=True)
    K_filtered_cols = K_csc[:, cols_to_keep]
    K = K_filtered_cols.tocsr()

    F = np.delete(F, fixed_dof, axis=0)

    u = spsolve(K, F)
    rtu = np.dot(u.T, F)

    return u, rtu


def arclen(displ, dispf):
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
        CON.Arclen.xincr = np.zeros(displ.shape, dtype=float)
    ufdx = np.dot(dispf.T, CON.Arclen.xincr)
    if CON.niter == 1:
        CON.Arclen.xincr = np.zeros(displ.shape, dtype=float)
        displ = np.zeros(displ.shape, dtype=float)
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
    BC = fem_db.BC
    RightHand.reactions = (RightHand.residual[BC.fixed_dof] +
                           RightHand.external_load[BC.fixed_dof])
    CON = fem_db.SolveControl

    r_norm = np.dot(RightHand.residual[BC.free_dof].T, RightHand.residual[BC.free_dof])

    # 计算外部载荷范数
    f_norm = np.dot((RightHand.nominal_external_load[BC.free_dof] - RightHand.nominal_pressure[BC.free_dof]).T,
                    (RightHand.nominal_external_load[BC.free_dof] - RightHand.nominal_pressure[BC.free_dof]))
    f_norm = f_norm * CON.xlamb ** 2

    # 计算反应力范数
    e_norm = np.dot(RightHand.reactions.T, RightHand.reactions)

    # 计算最终范数
    r_norm = np.sqrt(r_norm / (f_norm + e_norm))
    return r_norm


def ArcLengthNewtonRaphsonAlgorithm():
    from global_assembly.ResidualAndStiffnessAssembly import ResidualAndStiffnessAssembly
    from global_assembly.PressureLoadAndStiffnessAssembly import PressureLoadAndStiffnessAssembly
    CON = fem_db.SolveControl
    BC = fem_db.BC
    GEOM = fem_db.Geom
    afail = CON.Arclen.afail
    arcln = CON.Arclen.arcln
    LOAD_CASE = fem_db.LoadCase

    while CON.xlamb < CON.xlmax and CON.incrm < CON.nincr:
        CON.incrm += 1
        """
        Update the load factor. The radius is adjusted to achieve target 
        iterations per increment arbitrarily dampened by a factor of 0.7.
        """
        if not afail and not CON.Arclen.farcl:
            arcln = arcln * (CON.Arclen.itarget / CON.Arclen.iterold) ** 0.7

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
            PressureLoadAndStiffnessAssembly(fem_db.ElementGroupHash[0])

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
            # Solve for iterative displacements.
            CON.niter += 1
            displ, _ = linear_solver(fem_db.global_k.stiffness, -RightHand.residual, BC.fixed_dof)

            # Solve for displacements (for a nominal load).
            dispf, _ = linear_solver(fem_db.global_k.stiffness,
                                     RightHand.nominal_external_load - RightHand.nominal_pressure,
                                     BC.fixed_dof)
            displ = arclen(displ, dispf)

            """
            %----------------------------------------------------------------
            % Update coordinates.
            % -Recompute equivalent nodal forces and assembles residual force, 
            %  excluding pressure contributions.
            % -Recompute and assembles tangent stiffness matrix components, 
            %  excluding pressure contributions.
            %----------------------------------------------------------------
            """
            GEOM.x.flatten()[BC.free_dof] += displ
            ResidualAndStiffnessAssembly(fem_db.ElementGroupHash[0])

            # Update nodal forces due to pressure.
            if LOAD_CASE.n_pressure_loads:
                PressureLoadAndStiffnessAssembly(fem_db.ElementGroupHash[0])

            # Check for equilibrium convergence.
            rnorm = CheckResidualNorm()

            # Break iteration before residual gets unrealistic (e.g., NaN).
            if abs(rnorm) > 1e7 or np.isnan(rnorm):
                CON.Arclen.afail = 1
                break

            # Update value of the old iteration number.
            CON.Arclen.iter_old = CON.niter

        # If convergence not achieved opt to restart from previous results or terminate.
        afail = CON.Arclen.afail
        if CON.niter >= CON.miter or CON.Arclen.afail:
            terminate = input('Solution not converged. Do you want to terminate the program (y/n) ?: ')
            if terminate.lower() == 'n':
                # Restart from previous step by decreasing the fixed arclength radius to half its initially fixed value.
                print('Restart from previous step by decreasing the fixed arclength radius to half its initially fixed value.')
                CON.Arclen.arcln /= 2
            else:
                afail = 1
        else:
            # If convergence was achieved, update and save internal variables.
            # PLAST = save_output(updated_PLAST, PRO, FEM, GEOM, QUADRATURE, BC, MAT, LOAD, CON, CONSTANT, GLOBAL, PLAST, KINEMATICS)
            pass

        # If a failure occurs in the arc-length function, reduce the value of the radius in the arc-length to a 10th of the previous value and automatically restart.
        if afail:
            CON.Arclen.arcln /= 10
