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
        else:
            self.coord = np.asarray([x, y, z], dtype=float)

        # 属性设置
        self.vtk_type = ""
        self.is_boundary_node = False  # 节点是否为边界节点, 即自由度是否有被约束
        self.is_assist_node = False  # 定义梁方向的节点为辅助节点, 在计算总刚维度的时候不予考虑

        # 结果保存
        self.displacement = None  # 节点位移, 计算方法为 np.sqrt(np.square(dx,dy,dz))
        self.stress = []
        self.average_stress = None

        # 根据单元自由度改变的量, 默认节点有3个自由度
        self.dof_disp = np.asarray([None] * 3, dtype=float)
        self.eq_num = np.asarray([0] * 3, dtype=np.uint32)
        self.b_code = [False] * 3  # 如果自由度被约束, 那么为True
        self.dof_count = 3

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

    def ChangeDofCount(self, dof_count):
        """
        更改节点自由度的个数, 结构单元为6, 连续介质模型为3, 会影响到位移的长度、方程号的长度
        """
        self.dof_count = dof_count
        self.dof_disp = np.asarray([None] * dof_count, dtype=float)
        self.eq_num = np.asarray([0] * dof_count, dtype=np.uint32)
        self.b_code = [False] * dof_count

    def SetEquationNumber(self, idx, eq_num):
        """
        节点每个自由度对应的方程号
        """
        self.eq_num[idx] = eq_num

    def SetAllDofEqNum(self, eq_num):
        """
        如果节点是内部节点, 那么节点的所有自由度未被约束
        :return: 节点自由度个数
        """
        for i in range(len(self.eq_num)):
            self.eq_num[i] = eq_num + i
        return len(self.eq_num)

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

    def GetEquationNumbers(self):
        """
        返回节点对应的方程号
        """
        return self.eq_num

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
