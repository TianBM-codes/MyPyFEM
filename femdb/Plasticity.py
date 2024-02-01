#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np


class Stress:
    def __init__(self, Cauchy=None, Cauchyaa=None):
        self.Cauchy = Cauchy
        self.Cauchyaa = Cauchyaa


class Updated:
    def __init__(self, invCp=None, epbar=None):
        self.invCp = invCp
        self.epbar = epbar


class Trial:
    def __init__(self, lambdae=None, tau=None, n=None):
        self.lambdae = lambdae
        self.tau = tau
        self.n = n


class YieldCondition:
    def __init__(self, f=None, Dgamma=None, nu_a=None):
        self.f = f
        self.Dgamma = Dgamma
        self.nu_a = nu_a


class PlasticDeformationState:
    def __init__(self):
        self.epbar = None
        self.invCp = None
        self.oldEpbar = None
        self.oldInvCp = None

    def InitVariant(self, ngauss, nelem, ndim):
        self.epbar = np.zeros((ngauss, nelem))
        self.invCp = np.reshape(np.repeat(np.eye(ndim), ngauss * nelem),
                                (ndim, ndim, ngauss, nelem))


class Plasticity(object):
    def __init__(self):
        self.plastic_deformation_state = PlasticDeformationState()
        self.yield_condition = YieldCondition()
        self.trial = Trial()
        self.update = Updated()
        self.stress = Stress()

    def InitPDeformationState(self, ngauss, nelem, ndim):
        self.plastic_deformation_state.InitVariant(ngauss, nelem, ndim)

    def InitCalculate(self, stress: Stress, updated: Updated, trial: Trial, yield_condition: YieldCondition):
        self.stress = stress
        self.update = updated
        self.trial = trial
        self.yield_condition = yield_condition


if __name__ == "__main__":
    ppp = Plasticity()
