# !/usr/bin/env python3
# -*- coding: utf-8 -*-

from femdb.NLFEMDataBase import NLFEMDataBase
from utils.GlobalEnum import *
from femdb.Plasticity import *
from femdb.Material import HyperElasticPlasticInPrincipal, MaterialBase
from utils.CustomException import *
from scipy import linalg
import numpy as np

"""
Single instance mode, convenient for programming
"""
fem_db = NLFEMDataBase()
kinematics = fem_db.kinematics


def CauchyTypeSelection(PLAST_element: PlasticDeformationState, gauss_index: int, mat: MaterialBase):
    if isinstance(mat, HyperElasticPlasticInPrincipal):
        return HyperElasticPlasticInPrincipalStress(PLAST_element, gauss_index, mat)
    else:
        mlogger.fatal("Unknown Material Stress")
        raise NoImplSuchMaterialStress(mat.GetName())


def HyperElasticPlasticInPrincipalStress(PLAST_element: PlasticDeformationState, igauss: int, mat: MaterialBase):
    """
    Box 7.1 Algorithm for Rate-Independent Von Mises Plasticity with Isotropic Hardening
    Reference:
    1. scipy.linalg.eig return complex:
        https://stackoverflow.com/questions/8765310/scipy-linalg-eig-return-complex-eigenvalues-for-covariance-matrix
    """
    # TODO: prove invCp is symmetric
    dim = GetDomainDimension()

    PLAST_gauss = PlasticityGauss()
    PLAST_gauss.old.epbar = PLAST_element.epbar[igauss]
    PLAST_gauss.old.invCp = PLAST_element.invCp[:, :, igauss]
    epbar = PLAST_gauss.old.epbar
    invCp = PLAST_gauss.old.invCp

    """
    Trial stage
    """
    be_trial = np.matmul(kinematics.F[:, :, igauss],
                         np.matmul(invCp, kinematics.F[:, :, igauss].T))
    D, V = linalg.eigh(be_trial)
    lambdae_trial = np.sqrt(D)
    na_trial = V
    mu = mat.value_dict[MaterialKey.Niu]
    tauaa_trial = (2 * mu) * np.log(lambdae_trial) - \
                  (2 * mu / 3) * np.log(np.linalg.det(kinematics.F[:, :, igauss]))
    tau_trial = np.zeros(dim, dtype=float)
    for ii in range(dim):
        tau_trial = tau_trial + tauaa_trial[ii] * np.outer(na_trial[:, ii], na_trial[:, ii])

    """
    Checking for yielding
    """
    H = mat.value_dict[MaterialKey.Harden]
    ty0 = mat.value_dict[MaterialKey.TauY]
    f = np.sqrt(3 / 2) * np.linalg.norm(tau_trial, 'fro') - (
            ty0 + H * epbar)

    if f > 0:
        """
        Radial return algorithm. Return Dgamma (increment in the plastic
        multiplier) and the direction vector nua.
        """
        denominator = np.sqrt(2 / 3) * np.linalg.norm(tau_trial, 'fro')
        nu_a = tauaa_trial / denominator
        Dgamma = f / (3 * mu + H)

        lambdae_trial = np.array([0.5, 0.8, 1.2])  # 替换为实际的 lambdae_trial 数组
        norm_tauaa_trial = np.linalg.norm(tauaa_trial, 'fro')  # 计算 tauaa_trial 的 Frobenius 范数
        lambdae = np.exp(np.log(lambdae_trial) - Dgamma * nu_a)
        tauaa = (1 - 2 * mu * Dgamma / (np.sqrt(2 / 3) * norm_tauaa_trial)) * tauaa_trial
        tau = np.zeros((dim, dim), dtype=float)
        be = np.zeros((dim, dim), dtype=float)
        for ii in range(dim):
            tau += tauaa[ii] * np.outer(na_trial[:dim, ii], na_trial[:dim, ii])
            be += lambdae[ii] ** 2 * np.outer(na_trial[:dim, ii], na_trial[:dim, ii])

    else:
        tau = tau_trial[:dim, :dim]
        tauaa = tauaa_trial[:dim]
        Dgamma = 0
        nu_a = np.zeros(dim, dtype=float)
        be = be_trial[:dim, :dim]

    """
    Obtain the Cauchy stress tensor
    """
    PLAST_gauss.stress.Cauchy = np.asarray(tau / kinematics.J[igauss], dtype=float)
    PLAST_gauss.stress.Cauchyaa = np.asarray(tauaa / kinematics.J[igauss], dtype=float)

    """
    Update plasticity variables. 
    """
    invF = np.linalg.inv(kinematics.F[:, :, igauss])
    PLAST_gauss.update.invCp = np.dot(invF, np.dot(be, invF.T))
    PLAST_gauss.update.epbar = epbar + Dgamma

    """
    Store trial values (necessary for evaluation of the elasticity tensor).
    """
    PLAST_gauss.trial.lambdae = lambdae_trial[:dim]
    PLAST_gauss.tau_trial_dim = tau_trial[:dim, :dim]
    PLAST_gauss.trial.n = na_trial[:dim, :dim]

    """
    Store the yield criterion (necessary for evaluation of the elasticity tensor).
    """
    PLAST_gauss.yield_info.f = f
    PLAST_gauss.yield_info.Dgamma = Dgamma
    PLAST_gauss.yield_info.nu_a = nu_a[:dim]

    """
    Update the value of the internal variables at a Gauss point 
    (not for trusses).
    """
    PLAST_element.invCp[:, :, igauss] = PLAST_gauss.update.invCp
    PLAST_element.epbar[igauss] = PLAST_gauss.update.epbar

    return PLAST_gauss.stress.Cauchy, PLAST_gauss
