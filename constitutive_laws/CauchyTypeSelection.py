# !/usr/bin/env python3
# -*- coding: utf-8 -*-

from femdb.NLFEMDataBase import NLFEMDataBase
from GlobalEnum import *
from Kinematics import Kinematics
from femdb.Plasticity import *
from Material import HyperElasticPlasticInPrincipal
from element.ElementBase import AllEleTypeDNDrAtGaussianPoint
from utils.CustomException import *
from femdb.NLDomain import NLDomain

"""
Single instance mode, convenient for programming
"""
nl_domain = NLDomain()
plastics = nl_domain.plastics
kinematics = nl_domain.kinematics


def CauchyTypeSelection(gauss_index, mat):
    if isinstance(mat, HyperElasticPlasticInPrincipal):
        return HyperElasticPlasticInPrincipalStress(gauss_index, mat)
    else:
        mlogger.fatal("Unknown Material Stress")
        raise NoImplSuchMaterialStress(mat.GetName())


def HyperElasticPlasticInPrincipalStress(gauss_index, mat):
    dim = GlobalInfor[GlobalVariant.Dimension]
    """
    Box 7.1 Algorithm for Rate-Independent Von Mises Plasticity with Isotropic Hardening
    """
    plastics.plastic_deformation_state.oldEpbar = plastics.plastic_deformation_state.epbar
    plastics.plastic_deformation_state.oldInvCp = plastics.plastic_deformation_state.invCp

    """
    Trial stage
    """
    be_trial = np.matmul(kinematics.F,
                         np.matmul(plastics.plastic_deformation_state.oldInvCp, kinematics.F.T))
    V, D = np.linalg.eig(be_trial)
    lambdae_trial = np.sqrt(np.diag(D))
    na_trial = V
    mu = mat.value_dict[MaterialKey.Niu]
    tauaa_trial = (2 * mu) * np.log(lambdae_trial) - (2 * mu / 3) * np.log(np.linalg.det(kinematics.F))
    tau_trial = np.zeros(dim)
    for ii in range(dim):
        tau_trial += tauaa_trial[ii] * np.outer(na_trial[:, ii], na_trial[:, ii])

    """
    Checking for yielding
    """
    H = mat.value_dict[MaterialKey.Harden]
    ty0 = mat.value_dict[MaterialKey.TauY]
    f = np.sqrt(3 / 2) * np.linalg.norm(tau_trial, 'fro') - (
            ty0 + H * plastics.plastic_deformation_state.oldEpbar)

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
        tau = np.zeros((dim, dim))
        be = np.zeros((dim, dim))
        for ii in range(dim):
            tau += tauaa[ii] * np.outer(na_trial[:dim, ii], na_trial[:dim, ii])
            be += lambdae[ii] ** 2 * np.outer(na_trial[:dim, ii], na_trial[:dim, ii])
    else:
        tau = tau_trial[:dim, :dim]
        tauaa = tauaa_trial[:dim]
        Dgamma = 0
        nu_a = np.zeros(dim)
        be = be_trial[:dim, :dim]

    """
    Obtain the Cauchy stress tensor
    """
    Cauchy = tau / kinematics.J
    Cauchyaa = tauaa / kinematics.J
    invF = np.linalg.inv(kinematics.F)
    invCp_updated = np.dot(invF, np.dot(be, invF.T))
    epbar_updated = plastics.plastic_deformation_state.oldEpbar + Dgamma
    lambdae_trial_dim = lambdae_trial[:dim]
    tau_trial_dim = tau_trial[:dim, :dim]
    na_trial_dim = na_trial[:dim, :dim]
    yield_object = YieldCondition(f, Dgamma, nu_a[:dim])

    plastics.InitCalculate(Stress(Cauchy, Cauchyaa),
                           Updated(invCp_updated, epbar_updated),
                           Trial(lambdae_trial_dim, tau_trial_dim, na_trial_dim),
                           yield_object)

    return Cauchy
