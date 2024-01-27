# !/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys

from femdb.NLFEMDataBase import NLFEMDataBase
from GlobalEnum import *
from Kinematics import Kinematics
from Material import Plastics
from Material import HyperElasticPlasticInPrincipal
from LoadCase import RightHandItem
from element.ElementBase import AllEleTypeDNDrAtGaussianPoint
from utils.CustomException import *
import numpy as np


class IdentityTensor(object):
    def __init__(self, dimension):
        """
        Obtain entities which will be constant and only computed once.
        components of fourth order isotropic tensors
        c1 = delta(i,j)*delta(k,l)
        c2 = delta(i,k)*delta(j,l) + delta(i,l)*delta(j,k)
        (see textbook example 2.8)
        """
        I = np.eye(dimension)
        c1 = np.zeros((dimension, dimension, dimension, dimension))
        c2 = np.zeros((dimension, dimension, dimension, dimension))

        for l in range(dimension):
            for k in range(dimension):
                for j in range(dimension):
                    for i in range(dimension):
                        c1[i, j, k, l] = c1[i, j, k, l] + I[i, j] * I[k, l]
                        c2[i, j, k, l] = c2[i, j, k, l] + I[i, k] * I[j, l] + I[i, l] * I[j, k]


class AuxVariant(object):
    def __init__(self):
        self.ngauss = None
        self.n_dofs_elem = None
        self.weight = None
        self.DN_Dchi = None
        self.n_nodes_element = None


class NLDomain(object):
    def __init__(self):
        self.fem_db = NLFEMDataBase()
        self.identity_tensor = IdentityTensor(GlobalInfor[GlobalVariant.Dimension])
        self.right_hand_item = RightHandItem()
        self.kinematics = Kinematics()
        self.plastics = Plastics()
        self.aux_variant = AuxVariant()

    def initialisation(self):
        """
        初始化变量, 为求解做准备
        TODO: 线性的也变成这样？
        @return:
        """
        """
        Initialise undeformed geometry and initial residual and external forces. 
        Initialise external force vector contribution due to pressure
        (nominal value prior to load increment).
        """
        self.fem_db.Geom.InitX()
        self.right_hand_item.Init(self.fem_db.Mesh.n_dofs)

        """
        Initialisation of kinematics.
        """
        if GlobalInfor[GlobalVariant.InputFileSuffix] == InputFileType.FlagSHyP:
            for _, grp in self.fem_db.ElementGroupHash.items():
                self.aux_variant.ngauss = grp.eles[0].ngauss
                self.aux_variant.n_dofs_elem = grp.eles[0].n_dofs_elem
                self.aux_variant.weight, self.aux_variant.DN_Dchi = (
                    AllEleTypeDNDrAtGaussianPoint.GetElementDNDchi(grp.eles[0].e_type))
                self.aux_variant.n_nodes_element = grp.eles[0].nodes_count
                break

        self.kinematics.Init(GlobalInfor[GlobalVariant.Dimension],
                             self.fem_db.Mesh.n_nodes_elem,
                             self.aux_variant.ngauss
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
        self.ResidualAndStiffnessAssembly()

    def ResidualAndStiffnessAssembly(self):
        right_hand = self.right_hand_item
        n_dofs_elem = self.aux_variant.n_dofs_elem
        ngauss = self.aux_variant.ngauss
        ndim = GlobalInfor[GlobalVariant.Dimension]

        right_hand.external_load = self.fem_db.SolveControl.xlmax * right_hand.nominal_external_load

        """
        Pre-allocation memory to indexi, indexj and data for sparse assembly of
        the stiffness matrix.
        """
        n_components_mean_dilatation = self.fem_db.Material.n_nearly_incompressible * np.square(n_dofs_elem, 2)
        n_components_displacement_based = (self.fem_db.Mesh.nelem * np.square(n_dofs_elem, 2)
                                           + self.fem_db.Mesh.nelem * np.square(n_dofs_elem, 2) * ndim * ngauss)
        n_components = n_components_mean_dilatation + n_components_displacement_based

        """
        Initialise counter for storing sparse information into 
        global tangent stiffness matrix.
        """
        counter = 1
        indexi = np.zeros((n_components, 1))
        indexj = np.zeros((n_components, 1))
        global_stiffness = np.zeros((n_components, 1))

        """
        Main element loop
        """
        for _, grp in self.fem_db.ElementGroupHash.items():
            for ele in grp.eles:
                node_ids = ele.GetNodes()
                mat_id = ele.mat_id
                mat = self.fem_db.Material.Mat[mat_id]
                ele_id = ele.id_key
                epbar = self.plastics.epbar[:, ele_id]
                invCp = self.plastics.invCp[:, :, :, ele_id]
                # TODO 这里需要节点排好，对应关系另外存储，不要每次都查询dict, 只有在读取输入文件和写结果的时候初始化各种dict
                xlocal = self.fem_db.Geom.x[node_ids, :]
                x0local = self.fem_db.Geom.x0[node_ids, :]
                Ve = self.fem_db.Geom.V_ele[node_ids, :]
                self.ElementForceAndStiffness(xlocal, x0local, mat, Ve)

    def ElementForceAndStiffness(self, xlocal, Xlocal, mat, Ve):
        """
        Computes the element vector of global internal forces and the tangent
        stiffness matrix.
        @return:
        """
        T_internal = np.zeros((self.aux_variant.n_dofs_elem, 1))
        self.kinematics.ComputeGradients(xlocal, Xlocal, self.aux_variant.DN_Dchi)

        """
        Computes element mean dilatation kinematics, pressure and bulk modulus.
        """
        if isinstance(mat, HyperElasticPlasticInPrincipal):
            DN_x_mean = np.zeros(GlobalInfor[GlobalVariant.Dimension], self.aux_variant.n_nodes_element)
            ve = 0
            for ii in range(self.aux_variant.ngauss):
                JW = self.kinematics.Jx_chi[ii] * self.aux_variant.weight[ii]
                """
                Gauss contribution to the elemental deformed volume.Elemental averaged shape functions.
                """
                ve += JW
                DN_x_mean += self.aux_variant.DN_Dchi[:, :, ii] * JW
            DN_x_mean = DN_x_mean / ve
            Jbar = ve / Ve

            """
            Computes element mean dilatation pressure and bulk modulus.
            """
            if isinstance(mat, HyperElasticPlasticInPrincipal):
                kappa = mat.value_dict[MaterialKey.Kappa]
                press = kappa * np.log(Jbar) / Jbar
                kappa_bar = kappa / Jbar - press

        else:
            raise NoImplSuchMaterial(mat.GetName())

        """
        Gauss quadrature integration loop. Extract kinematics at the particular Gauss point.
        Obtain stresses (for incompressible or nearly incompressible, 
        only deviatoric component) and internal variables in plasticity.
        Obtain elasticity tensor (for incompressible or nearly incompressible, 
        only deviatoric component).
        """
        for jj in range(self.aux_variant.ngauss):
            self.CauchyTypeSelection(jj, mat)

    def CauchyTypeSelection(self, gauss_index, mat):
        if isinstance(mat, HyperElasticPlasticInPrincipal):
            """
            Box 7.1 Algorithm for Rate-Independent Von Mises Plasticity with Isotropic Hardening
            """
            self.plastics.oldEpbar = self.plastics.epbar
            self.plastics.invCp = self.plastics.oldInvCp

        else:
            mlogger.fatal("Unknown Material Stress")
            raise NoImplSuchMaterialStress(mat.GetName())


if __name__ == "__main__":
    GlobalInfor[GlobalVariant.Dimension] = 2
    dm = NLDomain()
