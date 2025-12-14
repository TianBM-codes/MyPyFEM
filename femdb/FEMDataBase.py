#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from utils.Singleton import Singleton
from femdb.LoadCase import *
from femdb.GlobalFEMVariant import ModelInfo
from collections import OrderedDict
from scipy import sparse


@Singleton
class FEMDataBase(object):
    """
    有限元数据库, 实例化的都存储在这里
    """

    def __init__(self):
        # 输入文件
        self.file_path = None

        # nodes & elements
        self.node_list = []  # List of all nodes in the domain, 实例化数据
        self.node_hash = {}  # 节点真实Id对应nodelist中的index的Hash表
        self.ele_hash = {}   # 单元真实Id对应elements中的index的Hash表
        self.elements = []
        self.equation_number = None
        self.node_connected_element_count = None

        # node & element sets
        self.element_set_name_hash = {}
        self.node_set_name_hash = {}
        self.element_sets = []
        self.node_sets = []

        # constrain equation
        self.equation_constrain_couple = []
        self.equation_constrain_idx = []
        self.temperature_constrain = []

        # Preprocess Fem
        self.matrix_num_count = 0
        self.stiff_list = []
        self.thermal_stiff_list = []
        self.global_stiff_matrix = None
        self.global_mass_matrix = None
        self.load_case = LoadCase()

        # Amplitudes
        self.amplitudes = {}

        # Result
        self.linear_u = None
        self.linear_mises = None

        self.sigma_xx = None
        self.sigma_yy = None
        self.sigma_zz = None
        self.tau_xy = None
        self.tau_xz = None
        self.tau_yz = None

        self.history_u = None
        self.history_v = None
        self.history_a = None
        self.history_s = None
        self.history_step_count = 0
        self.temperature_res = None

        # Plot
        self.additional_elements = []

    """ 
    以下的函数为解析文件的相关函数, 添加节点、单元、节点集、单元集、属性、材料、边界条件、LoadCase等 
    """

    def AddNode(self, node):
        """
        向有限元模型中插入节点
        :param node: 节点
        """
        self.node_list.append(node)

    def GetNodeBySearchId(self, n_id):
        return self.node_list[n_id]

    def GetModelSummary(self):
        """
        获取解析模型的信息, 包括文件位置、分析类型、单元个数、节点个数、EquationNumber(总自由度个数减去约束)
        :return:
        """
        summary_dict = OrderedDict()
        summary_dict["File Path"] = str(self.file_path)
        if not self.equation_number:
            summary_dict["Number of Equation"] = ModelInfo.PER_NODE_DOF * len(self.node_list)
        else:
            summary_dict["Number Of Equation"] = self.equation_number
        summary_dict["Number Of Node"] = len(self.node_list)  # TODO: 并不是标准的节点个数, 标准节点个数应该是由单元计算出来
        summary_dict["Number Of Element"] = len(self.elements)

        return summary_dict

    def GetNodeCoordBySearchId(self, indexes):
        """
        获取节点的坐标, 用于计算单元刚度阵
        """
        coords = []
        for idx in indexes:
            t_node = self.node_list[idx]
            coords.append((t_node.x, t_node.y, t_node.z))
        return coords

    def InitAssemblyMatrix(self, eq_count):
        """
        初始化总刚
        :param eq_count:
        """
        self.global_stiff_matrix = sparse.coo_matrix((eq_count, eq_count), dtype=float)
