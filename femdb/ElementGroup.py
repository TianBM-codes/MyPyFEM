#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from femdb.ElementFactory import *
from typing import List
from element.ElementBase import ElementBaseClass
from femdb.Plasticity import PlasticDeformationState


def GetElementNodeDofCount(e_type):
    """
    返回组成单元的节点所需的自由度个数
    :param e_type: 对于Abaqus是字符串, 对于ANSYS是int
    :return: 单元的自由度个数
    """
    # 1D Element
    if e_type in ["T3D2"]:
        return 3
    elif e_type in ["B31", 188]:
        return 6

    # 2D Element
    elif e_type in ["S3"]:
        return 6
    elif e_type in ["S4", "S4R", "S4RT"]:
        return 6
    elif e_type in ["CPS3", "CPS4"]:
        return 2

    # 3D Element
    elif e_type in ["C3D8", 185, "C3D6", "C3D4"]:
        return 3
    elif e_type in ["C3D8R"]:
        return 3
    elif e_type in ["C3D20R"]:
        return 3
    elif e_type in [181]:
        return 6

    raise NoImplSuchElement(e_type)


def GetElementNGauss(e_type):
    if e_type in ["hexa8"]:
        return 8
    raise NotImplementedError(e_type)


def GetNGauss(e_type):
    if e_type in ["hexa8"]:
        return 8

    raise NoImplSuchElement(e_type)


def GetNDofsElem(e_type):
    if e_type in ["hexa8"]:
        return 3 * 8
    raise NoImplSuchElement(e_type)


def GetNNodesElem(e_type):
    if e_type in ["hexa8"]:
        return 8
    raise NoImplSuchElement(e_type)


def GetNFaceDofsElem(e_type):
    if e_type in ["hexa8"]:
        return 12
    raise NoImplSuchElement(e_type)


def GetNodesCount(e_type):
    if e_type in ["hexa8"]:
        return 8
    raise NoImplSuchElement(e_type)


def GetBoundaryNGauss(e_type):
    if e_type in ["hexa8"]:
        return 4
    raise NoImplSuchElement(e_type)


def GetNFaceNodesElem(e_type):
    if e_type in ["hexa8"]:
        return 4
    raise NoImplSuchElement(e_type)


class ElementInfo:
    def __init__(self, e_type):
        self.ngauss = GetNGauss(e_type)
        self.n_dofs_elem = GetNDofsElem(e_type)
        self.n_nodes_elem = GetNNodesElem(e_type)
        self.n_face_dofs_elem = GetNFaceDofsElem(e_type)
        self.n_face_nodes_elem = GetNFaceNodesElem(e_type)
        self.nodes_count = GetNodesCount(e_type)
        self.boundary_ngauss = GetBoundaryNGauss(e_type)


class ElementGroup:
    """
    单元组, 每个单元组中存储的都是相同类型的单元, 类提供同一类型单元的一些特性
    """

    def __init__(self, e_type):
        from femdb.Interpolation import Interpolation
        self.e_type = e_type
        # 其中存储的是同一单元类型的单元, 存储了单元的全部信息
        self.eles: List[ElementBaseClass] = []
        # 单元组中的所有真实ID, 类型是set, 可以快速查找对应ID是否属于该组
        self.ele_id_sets = set()

        """
        Compute the shape functions and their derivatives
        (with respect to the iso parametric domain) inside the element and on
        faces (3D) or edges (2D).
        """
        self.element_info = ElementInfo(e_type)
        self.interpolation = Interpolation(self.e_type, self.element_info)
        self.quadrature = self.interpolation.quadrature

        """
        Plasticity variant
        """
        self.global_plasticity = PlasticDeformationState()

    def AppendElement(self, r_ele):
        """
        添加单元至单元组
        :param r_ele: 实际的单元
        """
        self.eles.append(r_ele)

    def SetEleIdSet(self, id_set):
        """
        设置单元组内包含的id, 因为需要查找某些单元是否属于该组, 所以
        设置为set类型, set是一种树结构, 比list列表结构快很多
        :param id_set: 节点id的set
        """
        self.ele_id_sets = id_set

    def Elements(self):
        return self.eles

    def IsElementInGroup(self, eid: int) -> bool:
        """
        通过单元id判断单元是否属于本Group, 从而判断单元的类型
        :param eid: 单元的真实ID
        """
        if eid in self.ele_id_sets:
            return True
        else:
            return False

    def GetElementsCurrentCount(self):
        """
        用于返回当前单元组已经添加了多少单元, 用于获取Index, 单元真实Id对应在ElementGroup中的eles中的idx
        """
        return len(self.eles)

    def GetZeroElementStiffMatrix(self):
        """
        返回一个空的Matrix, 在集成总刚的时候避免重复申请单刚
        """
        # dof_count = ElementFactory.GetElementNodeDofCount(self.e_type)
        # return np.zeros(dof_count, dof_count)

    def CalElementsEquationNum(self):
        """
        计算该单元组中所有单元包含节点对应的方程号, 为集成总刚做准备
        """
        # for e in self.eles:
        #     e.CalEquationNumber()

    def InitGlobalPlasticity(self):
        """
        Global Plasticity
        """
        self.global_plasticity.epbar = np.zeros((self.element_info.ngauss,
                                                 len(self.eles)), dtype=float)
        self.global_plasticity.invCp = np.reshape(np.repeat(np.eye(GetDomainDimension()), self.element_info.ngauss * len(self.eles)),
                                                  (GetDomainDimension(), GetDomainDimension(), self.element_info.ngauss, len(self.eles)))
