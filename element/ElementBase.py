#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import abc
import numpy as np
from femdb.GlobalEnum import *
from femdb.Integration import *
from utils.UtilsFunction import *
from scipy import sparse


class ElementBaseClass(metaclass=abc.ABCMeta):
    """
    Element base class
    All type of element classes should be derived from this base class
    单元不具有访问数据库的权限, 本构阵中的材料参数以及属性参数、方程号等由Domain传递过来计算, 否则
    会导致循环引用
    自定义单元的基类
    id : 节点的编号
    nodesCount : 包含节点个数
    pts : 包含的节点,其中为节点编号
    mat : 单元的材料
    nodeStr : 写入VTP时,需要将单元包含的节点写入,为了方便写成私有变量
          VTPType : int型,VTP规定的各个单元的编号
   新增单元需要做：
    1. 初始化id,nodesCount,VTPType
    2. 在ElementFactory类的构造函数和createElement函数中添加定义
    """

    def __init__(self, eid=None):
        """
        初始化单元实例
        """
        """
        必须要实例的成员, 因为都是流程必备的变量
        """
        self.id = eid
        self.id_key = None  # 用来重新排列后在map中的key值, self.id为value
        self.nodes_count = -1  # Number of nodes per element
        self.node_ids = np.asarray([])  # Node list of the element
        self.search_node_ids = np.asarray([])  # Domain中node_list中该节点的index
        self.cha_dict = None  # 单元的属性字典, 其中包括材料、属性、常数、惯性矩等
        self.vtp_type = None
        self.unv_code = None  # SiPESC平台显示的UNV结果, 单元代号
        self.eq_numbers = np.asarray([], dtype=np.uint32)  # 方程号, 即在求解矩阵中的第几行, 也即自由度排序后的index
        self.D = None  # 本构矩阵
        self.B = None  # 应变矩阵, 用于求解应力

        """
        包含节点的坐标, 假如有八个节点, dimension: 8 * 3,
        [[x1, y1, z1], 
         [x2, y2, z2],
          ...
         [x8, y8, z8]], type: np.ndarray, dtype: float
        """
        self.node_coords = None

        """
        不一定实例的成员, 这是由于CDB、BDF以及INP文件格式不统一造成的, 但这些变量都会转化为上述变量,
        例如在解析完成文件后会根据mat_id初始化self.cha_dict
        """
        self.mat_id = None
        self.sec_id = None
        self.real_const_id = None
        self.prop_id = None

    def SetId(self, eid: int):
        self.id = eid

    @abc.abstractmethod
    def ElementStiffness(self):
        """
        Calculate element stiffness matrix, coords of nodes
        (Upper triangular matrix, stored as an array column by colum)
        """
        pass

    @abc.abstractmethod
    def ElementStress(self, displacement):
        """ Calculate element stress """
        pass

    @abc.abstractmethod
    def CalElementDMatrix(self, an_type=None):
        """ Calculate element stress """
        pass

    def __eq__(self, other):
        return self.id == other.id

    def __lt__(self, other):
        """ 构造树的时候需要比较大小 """
        return self.id < other.id

    def SetAllCharacterAndCalD(self, cha_dict: dict):
        """ cha_dict必须包括计算单元刚度阵的所有内容 """
        self.cha_dict = cha_dict
        self.CalElementDMatrix()

    def SetEquationNumber(self, eq_nums: list[int]):
        self.eq_numbers = eq_nums

    """
    设置类相关函数
    """

    def SetNodes(self, nds: np.ndarray):
        """
        设置单元节点, 真实ID
        :param nds: 必须是列表, 内容为节点id, int型
        """
        self.node_ids = nds

    def SetNodeSearchIndex(self, idx: np.ndarray):
        """
        设置单元中所有节点对应的搜索编号
        """
        self.search_node_ids = idx

    def SetNodeCoords(self, coords: np.ndarray):
        """
        初始化单元包括的节点坐标矩阵, 包含节点的坐标, 假如有八个节点, dimension: 8 * 3,
        [[x1, y1, z1],
         [x2, y2, z2],
          ...
         [x8, y8, z8]], type: np.ndarray, dtype: float
        """
        self.node_coords = coords

    """
    获取信息相关函数
    """

    def GetNodes(self):
        """
        Return nodes of the element
        """
        return self.node_ids

    def GetNodeSearchIndex(self):
        """
        search_node_ids的含义: Domain中node_list中该节点的index
        """
        return self.search_node_ids

    def GetElementEquationNumber(self):
        """
        方程号, 即在求解矩阵中的第几行, 也即自由度排序后的index
        """
        return self.eq_numbers


class DNDrCalculator:
    """
    计算各种单元的dNdr, 也就是同一种单元涉及到的相同的地方, 这样就不需要每个单元都计算一遍,
    同一种单元计算一次即可
    """

    def __init__(self):
        # solid part
        self.C3D4 = None
        self.C3D6 = None
        self.C3D8 = None

        # shell part
        self.DKTShell = None
        self.DKQShell = None

        # dNdr: 将各种整合起来
        self.DNdr = {"C3D6": self.C3D6,
                     "C3D8": self.C3D8,
                     "DKTShell": self.DKTShell,
                     "DKQShell": self.DKQShell}

        self.CalculateAllDNDr()

    def CalculateAllDNDr(self):
        self.CalculateC3D6()
        self.CalculateC3D8()
        self.CalculateDKTShell()
        self.CalculateDKQShell()

    def GetElementDNdr(self, ele_type):
        return self.DNdr[ele_type]

    def CalculateC3D8(self):
        # Gaussian Weight
        sample_pt, weight = GaussIntegrationPoint.GetSamplePointAndWeight(2)

        # 在8个高斯点上积分
        dNdrs = []
        for ri in range(2):
            for si in range(2):
                for ti in range(2):
                    r, s, t = sample_pt[ri], sample_pt[si], sample_pt[ti]
                    dNdr = 0.125 * np.asarray([[(s + 1) * (-1 - t), (1 - r) * (1 + t), (1 - r) * (1 + s)],
                                               [(s - 1) * (1 + t), (r - 1) * (1 + t), (1 - r) * (1 - s)],
                                               [(s - 1) * (1 - t), (r - 1) * (1 - t), (r - 1) * (1 - s)],
                                               [(s + 1) * (t - 1), (1 - r) * (1 - t), (r - 1) * (1 + s)],
                                               [(1 + s) * (1 + t), (1 + r) * (1 + t), (1 + r) * (1 + s)],
                                               [(1 - s) * (1 + t), -(1 + r) * (1 + t), (1 + r) * (1 - s)],
                                               [(1 - s) * (1 - t), (1 + r) * (t - 1), (1 + r) * (s - 1)],
                                               [(1 + s) * (1 - t), (1 + r) * (1 - t), -(1 + r) * (1 + s)]]).T
                    g_weight = weight[ri] * weight[si] * weight[ti]
                    dNdrs.append(dNdr * g_weight)

        self.C3D8 = dNdrs

    def CalculateC3D6(self):
        # Gaussian Weight
        sample_r = [0.166666667, 0.666666667, 0.166666667, 0.166666667, 0.666666667, 0.166666667]
        sample_s = [0.166666667, 0.166666667, 0.666666667, 0.166666667, 0.166666667, 0.666666667]
        sample_t = [-0.577350269, -0.577350269, -0.577350269, 0.577350269, 0.577350269, 0.577350269]
        weight = 0.166666667

        # 在6个高斯点上积分
        dNdrs = []
        for ii in range(6):
            r, s, t = sample_r[ii], sample_s[ii], sample_t[ii]
            dNdrs.append(weight * np.asarray([[0.5 * (t - 1), 0.5 * (t - 1), 0.5 * (r + s - 1)],
                                              [0.5 * (1 - t), 0, -0.5 * r],
                                              [0, 0.5 * (1 - t), -0.5 * s],
                                              [-0.5 * (1 + t), -0.5 * (1 + t), 0.5 * (1 - r - s)],
                                              [0.5 * (1 + t), 0, 0.5 * r],
                                              [0, 0.5 * (1 + t), 0.5 * s]], dtype=float).T)

        self.C3D6 = dNdrs


AllEleTypeDNDr = DNDrCalculator()
