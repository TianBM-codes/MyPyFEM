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
        self.strain, self.stress = np.zeros(6), np.zeros(6)
        self.displacement = None  # 节点位移, 计算方法为 np.sqrt(np.square(dx,dy,dz))

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

    def SetBoundaryWithCDBType(self, direct_str:str, value:float):
        """
        ANSYS CDB文件施加约束的方式, 与INP文件的不同
        :param direct_str: 约束的方向
        :param value: 指定位移值
        """
        if direct_str.startswith("UX"):
            self.dof_disp[0] = value
            self.b_code[0] = True
        elif direct_str.startswith("UY"):
            self.dof_disp[1] = value
            self.b_code[1] = True
        elif direct_str.startswith("UZ"):
            self.dof_disp[2] = value
            self.b_code[2] = True
        elif direct_str.startswith("ROTX"):
            self.dof_disp[3] = value
            self.b_code[3] = True
        elif direct_str.startswith("ROTY"):
            self.dof_disp[4] = value
            self.b_code[4] = True
        elif direct_str.startswith("ROTZ"):
            self.dof_disp[5] = value
            self.b_code[5] = True
        else:
            mlogger.fatal("UnSupport boundary type:{}".format(str))
            sys.exit(1)
        self.is_boundary_node = True


    def SetBoundaryWithINPType(self, begin_idx=None, end_idx=None, value=0, b_type=None):
        """
        TODO: 对应不同自由度个数的情况补充完整
        设置边界条件, 如果idx大于6, 那么是其他的约束，比如温度, 还是壳单元自由度高？比如
        YSYMM u2=ur1=ur3=0(y轴平面应变) : 不能在Y方向上移动, 且能扫出XZ面, 法向是y, 就是没有Y轴
        XASYMM u2=u3=ru1=0 : X法向的圆形薄膜, 在x法向的压力下，一切均质的话, ur1=0, u2=u3=0
        """
        if b_type is None:
            if (begin_idx > 6) or (end_idx > 6):
                mlogger.fatal("Unknown AbaqusBoundary Index")
                sys.exit(1)
            for i in range(begin_idx, end_idx):
                self.dof_disp[i] = value
                self.b_code[i] = True
        elif b_type == "XSYMM":
            if self.dof_count == 6:
                self.dof_disp[0], self.dof_disp[3], self.dof_disp[4] = 0, 0, 0
                self.b_code[0], self.b_code[3], self.b_code[4] = True, True, True
        elif b_type == "YSYMM":
            self.dof_disp[1], self.dof_disp[3], self.dof_disp[5] = 0, 0, 0
            self.b_code[1], self.b_code[3], self.b_code[5] = True, True, True
        elif b_type == "ZSYMM":
            self.dof_disp[2], self.dof_disp[3], self.dof_disp[4] = 0, 0, 0
            self.b_code[2], self.b_code[3], self.b_code[4] = True, True, True
        elif b_type == "XASYMM":
            self.dof_disp[1], self.dof_disp[2], self.dof_disp[3] = 0, 0, 0
            self.b_code[1], self.b_code[2], self.b_code[3] = True, True, True
        elif b_type == "YASYMM":
            self.dof_disp[0], self.dof_disp[2], self.dof_disp[4] = 0, 0, 0
            self.b_code[0], self.b_code[2], self.b_code[4] = True, True, True
        elif b_type == "ZASYMM":
            self.dof_disp[0], self.dof_disp[1], self.dof_disp[5] = 0, 0, 0
            self.b_code[0], self.b_code[1], self.b_code[5] = True, True, True
        elif b_type == "PINNED":
            if self.dof_count == 3:
                self.dof_disp[0], self.dof_disp[1], self.dof_disp[2] = 0, 0, 0
                self.b_code[0], self.b_code[1], self.b_code[2] = True, True, True
        elif b_type == "ENCASTRE":
            self.dof_disp = np.asarray([0]*self.dof_count)
            self.b_code = [True] * self.dof_count

        else:
            mlogger.fatal("Unknown AbaqusBoundary Type: {}".format(b_type))
            sys.exit(1)
        self.is_boundary_node = True
