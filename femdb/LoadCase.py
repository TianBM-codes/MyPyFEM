#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys

import numpy as np
from typing import List
from utils.GlobalEnum import *


class AbaqusBoundary:
    def __init__(self, set_name, direction=None, value=None, b_type=None):
        # ABAQUS mode
        self.set_name = set_name
        self.direction = direction
        self.value = value
        self.b_type = b_type

    def __str__(self):
        return "AbaqusBoundary name: {}, direction: {}, value:{}, AbaqusBoundary Type:{}".format(
            self.set_name, self.direction, self.value, self.b_type)

    def GetSetName(self):
        """
        获取集合名字, 是ABAQUS的模式, 对本集合内的所有节点进行施加约束
        """
        return self.set_name

    def GetBoundaryType(self):
        """
        ABAQUS的约束类型包括: XSYMM YASYMM等
        """
        return self.b_type


class AnsysBoundary:
    def __init__(self):
        # ANSYS mode
        self.node_list = None
        self.directs = None
        self.con_values = None

    def SetConstraintInfor(self, node_list: list[int], directs: list[str], values: list[float]):
        """
        ANSYS格式的约束施加方式
        """
        assert len(node_list) == len(directs) == len(values)
        self.node_list = node_list
        self.directs = directs
        self.con_values = values

    def GetConstrainInfor(self):
        """
        暂时支持ANSYS的施加约束方式是每一行一个自由度
        """
        return self.node_list, self.directs, self.con_values


class NastranBoundary:
    # TODO: 还未添加测试
    def __init__(self):
        # Nastran mode
        self.node_list = None
        self.directs = None
        self.con_values = None

    def SetConstraintInfor(self, node_list: list[int], directs: list[str], values: list[float]):
        """
        Nastran格式的约束施加方式
        """
        assert len(node_list) == len(directs) == len(values)
        self.node_list = node_list
        self.directs = directs
        self.con_values = values

    def GetConstrainInfor(self):
        """
        暂时支持ANSYS的施加约束方式是每一行一个自由度
        """
        return self.node_list, self.directs, self.con_values


class ConcentratedLoad(object):
    def __init__(self):
        pass


class InpConcentratedLoad(ConcentratedLoad):
    """ ABAQUS集中力类 """

    def __init__(self, set_name=None, direction=None, value=None):
        """
        :param set_name: 集中力施加的节点集(string)
        :param direction: 力的方向(int)
        :param value: 力的大小(float)
        """
        # Abaqus inp mode
        super().__init__()
        self.set_name = set_name
        self.direction = direction
        self.value = value


class CdbConcentratedLoad(ConcentratedLoad):
    """ ANSYS集中力类 """

    def __init__(self):
        super().__init__()
        self.direction, self.node, self.value = None, None, None
        self.force_type = {"FX": 0, "FY": 1, "FZ": 2}

    def SetCForce(self, node: int, direction: str, value: float):
        """
        TODO: 添加施加弯矩
        ANSYS CDB样式的添加集中力
        :param node: 节点的真实ID
        :param direction: 方向
        :param value: 大小
        """
        self.node = node
        if self.force_type.__contains__(direction):
            self.direction = self.force_type[direction]
        else:
            mlogger.fatal("UnSupport CForce Type:{}".format(direction))
            sys.exit(1)
        self.value = value


class FlagSHyPCLoad(ConcentratedLoad):
    """ FlagSHyP Concentrate Load"""

    def __init__(self, line):
        super().__init__()
        line_split = line.split()
        self.node_id = int(line_split[0])
        force_x = float(line_split[1])
        force_y = float(line_split[2])
        if len(line_split) == 3:
            if GlobalInfor[GlobalVariant.Dimension] == AnalyseDimension.ThreeDimension:
                self.c_force = [force_x, force_y, 0]
            else:
                self.c_force = [force_x, force_y]
        elif len(line_split) == 4:
            if GlobalInfor[GlobalVariant.Dimension] == AnalyseDimension.TwoDimension:
                self.c_force = [force_x, force_y]
        else:
            mlogger.fatal(f"Wrong CLoad FlagSHyP Format{line}")
            sys.exit(1)


class PressLoad(object):
    def __init__(self):
        self.ele_id = None
        self.face_node = None
        self.p_value = None


class FlagSHyPPressLoad(PressLoad):
    """ FlagSHyP Press Load"""

    def __init__(self, line):
        super().__init__()
        line_split = line.split()
        self.ele_id = line_split[0]
        self.face_node = [int(ii) for ii in line_split[1:-1]]
        self.p_value = float(line_split[-1])


class LoadCase(object):
    """ Class LoadCase is used to store load data """

    def __init__(self):
        self.c_loads: List[ConcentratedLoad] = []
        self.boundaries = []
        self.p_loads: List[PressLoad] = []
        self.gravity = None
        self.n_pressure_loads = None
        self.n_prescribed_displacements = None

    def __str__(self):
        self.case_ = "\n  Here is LoadCase:\n"
        desc = self.case_
        for bd in self.boundaries:
            desc += "  {}\n".format(bd)
        for cld in self.c_loads:
            desc += "  set name:{}, direction:{}, value:{}".format(cld.set_name, cld.direction, cld.value)
        return desc

    def AddBoundary(self, boundary):
        """
        无论是ANSYS、ABAQUS还是NASTRAN格式的约束都存在这里
        :param boundary: AnsysBoundary、AbaqusBoundary或者NastranBoundary类型
        :return: None
        """
        self.boundaries.append(boundary)

    def AddAbaqusCLoad(self, set_name, direction, value):
        """
        ABAQUS类型的添加集中力
        """
        self.c_loads.append(InpConcentratedLoad(set_name, direction, value))

    def AddAnsysCLoad(self, node: int, direction: str, value: float):
        """
        ANSYS类型的添加集中力
        :param node: 节点的真实ID
        :param direction: 集中力的方向
        :param value: 集中力的大小
        """
        cf = CdbConcentratedLoad()
        cf.SetCForce(node, direction, value)
        self.c_loads.append(cf)

    def AddFlagSHyPCLoad(self, line):
        """
        @return:
        """

    def AddFlagSHyPPressLoad(self, p_load: FlagSHyPPressLoad):
        self.p_loads.append(p_load)

    def GetBoundaries(self):
        return self.boundaries

    def GetConcentratedLoads(self):
        return self.c_loads

    def AddGravity(self, gravity):
        assert len(gravity) == 3
        self.gravity = gravity


class RightHandItem(object):
    def __init__(self):
        self.nominal_external_load = None
        self.residual = None
        self.external_load = None
        self.nominal_press = None
        self.T_int = None
        self.reactions = None

    def Init(self, n_dof):
        self.nominal_external_load = np.zeros((n_dof, 1))
        self.residual = np.zeros((n_dof, 1))
        self.external_load = np.zeros((n_dof, 1))
        self.nominal_press = np.zeros((n_dof, 1))
        self.T_int = np.zeros((n_dof, 1))
