# !/usr/bin/env python3
# -*- coding: utf-8 -*-

from femdb.NLFEMDataBase import NLFEMDataBase
from utils.GlobalEnum import *


class NLDomain(object):
    _instance = None  # 类变量用于存储唯一的实例
    _initialized = False

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(NLDomain, cls).__new__(cls)
        return cls._instance

    def __init__(self):
        if not self._initialized:
            self.fem_db = NLFEMDataBase()
            self._initialized = True

    def Initialisation(self):
        """
        初始化变量, 为求解做准备
        @return:
        """
        """
        Initialisation of internal variables for plasticity.
        """
        check = False
        for _, mat in self.fem_db.Material.Mat.items():
            if mat.GetName() in [2, 17]:
                check = True
        if check:
            # FlagSHyP only support one type element
            group = self.fem_db.ElementGroupHash[0]
            group.InitGlobalPlasticity()

        """
        Initialisation of kinematics.
        """
        for _, grp in self.fem_db.ElementGroupHash.items():
            grp_ele_info = grp.element_info
            self.fem_db.kinematics.Init(GetDomainDimension(),
                                        grp_ele_info.n_nodes_elem,
                                        grp_ele_info.ngauss
                                        )

        """
        Calculate initial volume for data checking. 
        Additionally, essential for mean dilation algorithm.
        """
        self.fem_db.Geom.InitialVolume()

        """
        Computes and assembles the initial tangent matrix and the initial  
        residual vector due to the internal contributions 
        (external contributions will be added later on). 
        """
        from global_assembly.ResidualAndStiffnessAssembly import ResidualAndStiffnessAssembly
        for _, grp in self.fem_db.ElementGroupHash.items():
            ResidualAndStiffnessAssembly(grp)

    def ChooseIncrementalAlgorithm(self):
        from solver.NewtonRaphsonAlgorithm import NewtonRaphsonAlgorithm
        from solver.LineSearchNewtonRaphsonAlgorithm import LineSearchNewtonRaphsonAlgorithm
        from solver.ArcLengthNewtonRaphsonAlgorithm import ArcLengthNewtonRaphsonAlgorithm
        CON = self.fem_db.SolveControl
        if abs(CON.Arclen.arcln) == 0:
            if not CON.searc:
                NewtonRaphsonAlgorithm()
            else:
                LineSearchNewtonRaphsonAlgorithm()
        else:
            ArcLengthNewtonRaphsonAlgorithm()


if __name__ == "__main__":
    GlobalInfor[GlobalVariant.Dimension] = 2
    dm = NLDomain()
