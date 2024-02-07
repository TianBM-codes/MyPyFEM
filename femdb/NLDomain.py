# !/usr/bin/env python3
# -*- coding: utf-8 -*-

from femdb.NLFEMDataBase import NLFEMDataBase
from GlobalEnum import *
from Kinematics import Kinematics
from femdb.Plasticity import *
from LoadCase import RightHandItem
from element.ElementBase import AllEleTypeDNDrAtGaussianPoint
from scipy.sparse import coo_matrix, csr_matrix
from utils.CustomException import *


class IdentityTensor(object):
    def __init__(self, dimension):
        """
        Obtain entities which will be constant and only computed once.
        components of fourth order isotropic tensors
        c1 = delta(i,j)*delta(k,l)
        c2 = delta(i,k)*delta(j,l) + delta(i,l)*delta(j,k)
        (see textbook example 2.8)
        """
        self.I = np.eye(dimension)
        self.c1 = np.zeros((dimension, dimension, dimension, dimension))
        self.c2 = np.zeros((dimension, dimension, dimension, dimension))

        for l in range(dimension):
            for k in range(dimension):
                for j in range(dimension):
                    for i in range(dimension):
                        self.c1[i, j, k, l] = self.c1[i, j, k, l] + self.I[i, j] * self.I[k, l]
                        self.c2[i, j, k, l] = (self.c2[i, j, k, l] +
                                               self.I[i, k] * self.I[j, l] +
                                               self.I[i, l] * self.I[j, k])


class AuxVariant(object):
    def __init__(self):
        self.ngauss = None
        self.n_dofs_elem = None
        self.weight = None
        self.DN_Dchi = None
        self.n_nodes_element = None
        self.n_face_dofs_elem = None


class GlobalK(object):
    def __init__(self):
        self.indexi = None
        self.indexj = None
        self.counter = None
        self.stiffness = None


class NLDomain(object):
    _instance = None  # 类变量用于存储唯一的实例

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(NLDomain, cls).__new__(cls)
        return cls._instance

    def __init__(self):
        self.fem_db = NLFEMDataBase()
        self.identity_tensor = IdentityTensor(GlobalInfor[GlobalVariant.Dimension])
        self.right_hand_item = RightHandItem()
        self.kinematics = Kinematics()
        self.plastics = Plasticity()
        self.aux_variant = AuxVariant()
        self.global_k = GlobalK()

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
        
        Welcomes the user and determines whether the problem is being
        restarted or a data file is to be read. 
        Reads all necessary input data.
        Initialises kinematic variables and internal variables.
        Compute initial tangent matrix and equivalent force vector, excluding
        pressure component.
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
                self.aux_variant.n_face_dofs_elem = grp.eles[0].n_face_dofs_elem
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
        self.global_k.counter = 1
        self.global_k.indexi = np.zeros((n_components, 1))
        self.global_k.indexj = np.zeros((n_components, 1))
        self.global_k.stiffness = np.zeros((n_components, 1))

        """
        Main element loop
        """
        for _, grp in self.fem_db.ElementGroupHash.items():
            for ele in grp.eles:
                node_ids = ele.GetNodes()
                mat_id = ele.mat_id
                mat = self.fem_db.Material.Mat[mat_id]
                ele_id = ele.id_key
                epbar = self.plastics.plastic_deformation_state.epbar[:, ele_id]
                invCp = self.plastics.plastic_deformation_state.invCp[:, :, :, ele_id]
                # TODO 这里需要节点排好，对应关系另外存储，不要每次都查询dict, 只有在读取输入文件和写结果的时候初始化各种dict
                xlocal = self.fem_db.Geom.x[node_ids, :]
                x0local = self.fem_db.Geom.x0[node_ids, :]
                Ve = self.fem_db.Geom.V_ele[node_ids, :]
                from element_calculation.ElementForceAndStiffness import ElementForceAndStiffness
                ElementForceAndStiffness(xlocal, x0local, mat_id, Ve, ele)

                """
                Storage of updated value of the internal variables. 
                """
                # TODO: 首先将节点排好序后，再储存 self.plasticity

        """
        Global tangent stiffness matrix sparse assembly except pressure contributions. 
        """
        S = coo_matrix(self.global_k.indexi, self.global_k.indexj, self.global_k.stiffness)
        self.global_k.stiffness = csr_matrix(S)
        """
        Compute global residual force vector except pressure contributions.
        """
        self.right_hand_item.residual = self.right_hand_item.T_int - self.right_hand_item.external_load

    def ChooseIncrementalAlgorithm(self):
        from solver.NewtonRaphsonAlgorithm import NewtonRaphsonAlgorithm
        from solver.LineSearchNewtonRaphsonAlgorithm import LineSearchNewtonRaphsonAlgorithm
        from solver.ArcLengthNewtonRaphsonAlgorithm import ArcLengthNewtonRaphsonAlgorithm
        CON = self.fem_db.SolveControl
        if abs(CON.arcln) == 0:
            if not CON.searc:
                NewtonRaphsonAlgorithm()
            else:
                LineSearchNewtonRaphsonAlgorithm()
        else:
            ArcLengthNewtonRaphsonAlgorithm()


if __name__ == "__main__":
    GlobalInfor[GlobalVariant.Dimension] = 2
    dm = NLDomain()
