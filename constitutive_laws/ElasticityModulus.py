# !/usr/bin/env python3
# -*- coding: utf-8 -*-

from femdb.NLFEMDataBase import NLFEMDataBase
from GlobalEnum import *
from femdb.Plasticity import *
from Material import HyperElasticPlasticInPrincipal
from utils.CustomException import *
from femdb.NLDomain import NLDomain

"""
Single instance mode, convenient for programming
"""
nl_domain = NLDomain()
PLAST = nl_domain.plastics
kinematics = nl_domain.kinematics
dim = GlobalInfor[GlobalVariant.Dimension]

fem_db = NLFEMDataBase()
MAT = fem_db.Material.Mat


def muab_choice(lambda_alpha, lambda_beta, sigma_alpha, sigma_beta, J, mu):
    if abs(lambda_alpha - lambda_beta) < 1e-5:
        muab = mu / J - sigma_alpha
    else:
        muab = (sigma_alpha * lambda_beta ** 2 - sigma_beta * lambda_alpha ** 2) / (
                lambda_alpha ** 2 - lambda_beta ** 2)
    return muab


def ElasticityModulusSelection(mat_id):
    mat = MAT[mat_id]
    if isinstance(mat, HyperElasticPlasticInPrincipal):
        return HyperElasticPlasticInPrincipalEModulus(mat_id, kinematics)
    else:
        raise NoImplSuchElasticityModulus(mat.GetName())


def HyperElasticPlasticInPrincipalEModulus(mat_id, kine):
    """
    P212 BOX 7.2: Tangent Modulus
    @param mat_id: material id
    @param kine:
    @return:
    """
    mu = MAT[mat_id].value_dict[MaterialKey.Niu]
    H = MAT[mat_id].value_dict[MaterialKey.Harden]
    J = kine['J']

    lambdae_trial = PLAST['trial']['lambdae']
    tau_trial = PLAST['trial']['tau']
    na_trial = PLAST['trial']['n']
    T = na_trial

    Cauchyaa = PLAST['stress']['Cauchyaa']

    Dgamma = PLAST['yield']['Dgamma']
    nua = PLAST['yield']['nu_a']

    if PLAST['yield']['f'] > 0:
        denominator = np.sqrt(2 / 3) * np.linalg.norm(tau_trial, 'fro')
        c_alphabeta = (
                (1 - 2 * mu * Dgamma / denominator) *
                (2 * mu * np.eye(dim) - 2 / 3 * mu * np.ones((dim, dim))) -
                (2 * mu) ** 2 * np.outer(nua, nua) *
                (
                        1 / (3 * mu + H) -
                        np.sqrt(2 / 3) * Dgamma /
                        np.linalg.norm(tau_trial, 'fro')
                )
        )
    else:
        c_alphabeta = 2 * mu * np.eye(dim) - 2 / 3 * mu * np.ones((dim, dim))

    c = np.zeros((dim, dim, dim, dim))

    for l in range(dim):
        for k in range(dim):
            for j in range(dim):
                for i in range(dim):
                    sum_ = 0
                    for alpha in range(dim):
                        sum_ -= 2 * Cauchyaa[alpha] * T[i, alpha] * T[j, alpha] * T[k, alpha] * T[l, alpha]
                        for beta in range(dim):
                            sum_ += (c_alphabeta[alpha, beta] / J) * \
                                    (T[i, alpha] * T[j, alpha] * T[k, beta] * T[l, beta])
                            lambda_a = lambdae_trial[alpha]
                            lambda_b = lambdae_trial[beta]
                            sigma_a = Cauchyaa[alpha]
                            sigma_b = Cauchyaa[beta]
                            if alpha != beta:
                                sum_ += muab_choice(lambda_a, lambda_b, sigma_a, sigma_b, J, mu) * \
                                        (T[i, alpha] * T[j, beta] * (
                                                T[k, alpha] * T[l, beta] + T[k, beta] * T[l, alpha]))
                    c[i, j, k, l] += sum_

    return c
