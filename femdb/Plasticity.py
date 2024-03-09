#!/usr/bin/env python3
# -*- coding: utf-8 -*-


class Stress:
    def __init__(self, Cauchy=None, Cauchyaa=None):
        self.Cauchy = Cauchy
        self.Cauchyaa = Cauchyaa


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
        self.ep = None


class PlasticityGauss(object):
    def __init__(self):
        self.plastic_deformation_state = PlasticDeformationState()
        self.yield_info = YieldCondition()
        self.trial = Trial()
        self.update = PlasticDeformationState()
        self.old = PlasticDeformationState()
        self.stress = Stress()


if __name__ == "__main__":
    ppp = PlasticityGauss()
