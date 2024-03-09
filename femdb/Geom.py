# !/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from utils.GlobalEnum import MaterialKey


class Geom(object):
    """
    存储几何信息
    """

    def __init__(self):
        self.node_count = -1
        self.NdId2ListIdx = {}
        self.ListIdx2NdId = {}
        self.x0: np.array = None
        self.x: np.array = None
        self.V_ele: np.array = None
        self.V_total = -1

    def InitialVolume(self):
        """
        Calculate initial volume for data checking.
        Additionally, essential for mean dilation algorithm.
        @return:
        """
        from femdb.NLFEMDataBase import NLFEMDataBase
        fem_db = NLFEMDataBase()
        self.V_ele = np.zeros((fem_db.Mesh.nelem, 1), dtype=float)
        self.V_total = 0
        for _, ele_grp in fem_db.ElementGroupHash.items():
            eles = ele_grp.eles
            if ele_grp.e_type in ['truss2']:
                for jj in range(len(eles)):
                    ele = eles[jj]
                    x_local = fem_db.Geom.x[:, ele.search_node_ids]
                    area = fem_db.Material[ele.mat_id].value_dict[MaterialKey.Area]
                    L = np.linalg.norm(x_local[:,1]-x_local[:,0])
                    self.V_ele[jj] = area * L
                    self.V_total += self.V_ele[jj]
            else:
                for jj in range(len(eles)):
                    ele = eles[jj]
                    x_local = fem_db.Geom.x[:, ele.search_node_ids]
                    weight = ele_grp.quadrature.element_w
                    DN_Dchi = ele_grp.interpolation.element_DN_chi
                    fem_db.Kinematics.ComputeGradients(x_local, x_local, DN_Dchi)
                    for ii in range(ele_grp.quadrature.element_ngauss):
                        jw = fem_db.Kinematics.Jx_chi[ii] * weight[ii]
                        self.V_ele[jj] += jw
                    self.V_total += self.V_ele[jj]
