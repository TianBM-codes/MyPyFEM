#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import abc
import numpy as np
from femdb.GlobalEnum import *
from femdb.Integration import *


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
        self.nodes_count = 0  # Number of nodes per element
        self.node_ids = np.asarray([])  # Node list of the element
        self.search_node_ids = np.asarray([])  # Domain中node_list中该节点的index
        self.ele_mat_dict = None  # Pars of Material of the element
        self.ele_prop_dict = None  # 属性字典
        self.vtp_type = ""
        self.eq_numbers = np.asarray([], dtype=np.uint32)
        self.D = None  # 本构矩阵
        """
        包含节点的坐标, 假如有八个节点, dimension: 8 * 3,
        [[x1, y1, z1], [x2, y2, z2], ...[x8, y8, z8]], type: np.mat, dtype: float
        """
        self.node_coords = None
        """
        不一定实例的成员, 这是由于CDB、BDF以及INP文件格式不统一造成的, 但这些变量都会转化为上述变量,
        例如在解析完成文件后会根据mat_id初始化self.ele_mat_dict
        """
        self.mat_id = None
        self.sec_id = None
        self.real_const_id = None

    def SetId(self, eid):
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

    def SetMatAndCalBTDB(self, mat_dict, prop_dict):
        self.ele_mat_dict = mat_dict
        self.ele_prop_dict = prop_dict
        self.CalElementDMatrix()
        # self.ElementStiffness()

    def SetEquationNumber(self, eq_nums):
        self.eq_numbers = eq_nums

    """
    设置类相关函数
    """

    def SetNodes(self, nds):
        """
        设置单元节点
        :param nds: 必须是列表, 内容为节点id, int型
        """
        self.node_ids = nds

    def SetNodeSearchIndex(self, idx):
        """
        设置单元中所有节点对应的搜索编号
        """
        self.search_node_ids = idx

    def SetNodeCoords(self, coords):
        """
        初始化单元包括的节点坐标矩阵
        """
        self.node_coords = np.mat(coords)

    """
    获取信息相关函数
    """

    def GetNodes(self):
        """ Return nodes of the element """
        return self.node_ids

    def GetNodeSearchIndex(self):
        return self.search_node_ids

    def GetElementEquationNumber(self):
        return self.eq_numbers

    def GetElementMaterial(self):
        """ Return material of the element """
        return self.ele_mat_dict
