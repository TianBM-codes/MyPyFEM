# !/usr/bin/env python3
# -*- coding: utf-8 -*-

from femdb.NLFEMDataBase import NLFEMDataBase
from utils.GlobalEnum import *
from femdb.Plasticity import PlasticDeformationState, PlasticityGauss
from femdb.Material import HyperElasticPlasticInPrincipal
from utils.CustomException import *
import numpy as np

"""
Single instance mode, convenient for programming
"""
dim = GetDomainDimension()
fem_db = NLFEMDataBase()
MAT = fem_db.Material.Mat
kinematics = fem_db.kinematics


def muab_choice(lambda_alpha, lambda_beta, sigma_alpha, sigma_beta, J, mu):
    if abs(lambda_alpha - lambda_beta) < 1e-5:
        muab = mu / J - sigma_alpha
    else:
        muab = (sigma_alpha * lambda_beta ** 2 - sigma_beta * lambda_alpha ** 2) / (
                lambda_alpha ** 2 - lambda_beta ** 2)
    return muab


def ElasticityModulusSelection(PLAST_element: PlasticDeformationState,
                               PLAST_gauss: PlasticityGauss,
                               igauss: int,
                               mat_id: int):
    """
    Obtain elasticity tensor (for incompressible or nearly incompressible,
    only deviatoric component).
    :param PLAST_element:
    :param PLAST_gauss:
    :param igauss:
    :param mat_id:
    :return:
    """
    mat = MAT[mat_id]
    if isinstance(mat, HyperElasticPlasticInPrincipal):
        return HyperElasticPlasticInPrincipalEModulus(PLAST_element, PLAST_gauss, igauss, mat_id)
    else:
        raise NoImplSuchElasticityModulus(mat.GetName())


def HyperElasticPlasticInPrincipalEModulus(PLAST_element: PlasticDeformationState,
                                           PLAST_gauss: PlasticityGauss,
                                           igauss: int,
                                           mat_id: int):
    """
    P212 BOX 7.2: Tangent Modulus
    """
    mu = MAT[mat_id].value_dict[MaterialKey.Niu]
    H = MAT[mat_id].value_dict[MaterialKey.Harden]
    J = fem_db.kinematics.J[igauss]

    PLAST_gauss.old.invCp = PLAST_element.invCp[:, :, igauss]
    PLAST_gauss.old.epbar = PLAST_element.epbar[igauss]

    lambdae_trial = PLAST_gauss.trial.lambdae
    tau_trial = PLAST_gauss.trial.tau
    na_trial = PLAST_gauss.trial.n
    T = na_trial

    """
    Cauchy stress tensor and its eigenvalues.
    """
    Cauchyaa = PLAST_gauss.stress.Cauchyaa

    Dgamma = PLAST_gauss.yield_info.Dgamma
    nua = PLAST_gauss.yield_info.nu_a

    if PLAST_gauss.yield_info.f > 0:
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

    c = np.zeros((dim, dim, dim, dim), dtype=float)

    # T = np.asarray([[0.538280250409507, 0.648466455599791, -0.538280250409507],
    #                 [0.458535028126617, -0.761243230486711, -0.458535028126617],
    #                 [0.707106781186548, 0, 0.707106781186548]], dtype=float)

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
                    # try: if abs(sum_.imag) > 1e13:
                    c[i, j, k, l] += sum_
                    # if isinstance(sum_, complex):
                    #     c[i, j, k, l] += sum_.real
                    # else:
                    #     c[i, j, k, l] += sum_

    return c
