#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys

from femdb.GlobalEnum import *


class Boundary(object):
    def __init__(self, set_name, direction=None, value=None, b_type=None):
        # ABAQUS mode
        self.set_name = set_name
        self.direction = direction
        self.value = value
        self.b_type = b_type

        # ANSYS mode
        self.node_list = None
        self.directs = None
        self.values = None

    def __str__(self):
        return "Boundary name: {}, direction: {}, value:{}, Boundary Type:{}".format(
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

    def SetConstraintInfor(self, node_list: list, directs: list, values: list):
        """
        ANSYS格式的约束施加方式
        """
        assert len(node_list) == len(directs) == len(values)
        self.node_list = node_list
        self.directs = directs
        self.values = values


class ConcentratedLoad(object):
    """ 集中力类 """

    def __init__(self, set_name=None, direction=None, value=None):
        """
        :param set_name: 集中力施加的节点集(string)
        :param direction: 力的方向(int)
        :param value: 力的大小(float)
        """
        # Abaqus inp mode
        self.set_name = set_name
        self.direction = direction
        self.value = value

        # Ansys cdb mode
        self.direction, self.node, value = None, None, None

    def SetAnsysCForce(self, node: int, direction: str, value: float):
        """
        TODO: 添加施加弯矩
        ANSYS CDB样式的添加集中力
        :param node: 节点的真实ID
        :param direction: 方向
        :param value: 大小
        """
        self.node = node
        if direction.startswith("FX"):
            self.direction = 0
        elif direction.startswith("FY"):
            self.direction = 1
        elif direction.startswith("FZ"):
            self.direction = 2
        else:
            mlogger.fatal("UnSupport CForce Type:{}".format(direction))
            sys.exit(1)
        self.value = value


class LoadCase(object):
    """ Class LoadCase is used to store load data """

    def __init__(self):
        self.c_loads = []
        self.boundaries = []

    def __str__(self):
        self.case_ = "\n  Here is LoadCase:\n"
        desc = self.case_
        for bd in self.boundaries:
            desc += "  {}\n".format(bd)
        for cld in self.c_loads:
            desc += "  set name:{}, direction:{}, value:{}".format(cld.set_name, cld.direction, cld.value)
        return desc

    def AddBoundary(self, boundary):
        self.boundaries.append(boundary)

    def AddAbaqusCLoad(self, set_name, direction, value):
        """
        ABAQUS类型的添加集中力
        """
        self.c_loads.append(ConcentratedLoad(set_name, direction, value))

    def AddAnsysCLoad(self, node: int, direction: str, value: float):
        """
        ANSYS类型的添加集中力
        :param node: 节点的真实ID
        :param direction: 集中力的方向
        :param value: 集中力的大小
        """
        cf = ConcentratedLoad()
        cf.SetAnsysCForce(node,direction, value)
        self.c_loads.append(cf)

    def GetBoundaries(self):
        return self.boundaries

    def GetConcentratedLoad(self):
        return self.c_loads
