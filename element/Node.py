#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from femdb.GlobalEnum import *
import numpy as np
import sys


class Node(object):
    def __init__(self, nid, x=None, y=None, z=None):
        """
        # After call FEMDB.CalculateEquationNumber(), bcode stores the global equation number corresponding to each degree of freedom of the node
        """
        super().__init__()
        self.id = nid  # 节点在导入文件中的编号
        self.id_key = None  # 在Domain中, 如果节点号不连续, 那么存入Hash中的id号, self.id为value

        # 对于2D分析的支持
        if z is None:
            self.coord = np.asarray([x, y], dtype=float)
            self.origin_coord = np.asarray([x, y], dtype=float)
        else:
            self.coord = np.asarray([x, y, z], dtype=float)
            self.origin_coord = np.asarray([x, y, z], dtype=float)

        # 属性设置
        self.vtk_type = "vertex"
        self.is_assist_node = False  # 定义梁方向的节点为辅助节点, 在计算总刚维度的时候不予考虑

        # 结果保存
        self.displacement = None  # 节点位移, 计算方法为 np.sqrt(np.square(dx,dy,dz))
        self.stress = []
        self.average_stress = None

        # 根据单元自由度改变的量, 默认节点有3个自由度
        self.dof_disp = np.asarray([None] * 3, dtype=float)
        self.start_eq_num = None
        self.end_eq_num = None
        self.dof_count = None

    def __lt__(self, other):
        return self.id < other.id

    def __eq__(self, other):
        return self.id == other.id

    def GetId(self):
        return self.id

    def GetNodeCoord(self):
        return self.coord

    def GetDofCount(self):
        return self.dof_count

    def CalNodeMagnitudeDisplacement(self):
        """
        计算节点空间位移
        """
        if (self.dof_count == 3) or (self.dof_count == 6):
            self.displacement = np.sqrt(np.sum(np.square([self.dof_disp[0], self.dof_disp[1], self.dof_disp[2]])))
        elif self.dof_count == 2:
            self.displacement = np.sqrt(np.sum(np.square([self.dof_disp[0], self.dof_disp[1]])))
        else:
            mlogger.fatal("Un support dof count: {}".format(self.dof_count))
            sys.exit(1)

    def GetDisplacement(self):
        return self.displacement

    def AppendStressResult(self, stress: np.array):
        """
        给节点分配由单元计算来的应力结果, 然后平均即得到节点应力
        """
        self.stress.append(stress)

    def AverageStress(self):
        """
        平均节点的应力结果, 作为节点的应力结果输出
        """
        self.average_stress = np.sum(np.asarray(self.stress), axis=0)

    def SetDof(self, node_dof):
        """
        设置节点的自由度大小
        :param node_dof:
        :return:
        """
        if node_dof == 3:
            if self.dof_count != 6:
                self.dof_count = 3
                return
        else:
            self.dof_count = node_dof
