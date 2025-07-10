#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 导入必要的模块
from utils.Singleton import Singleton  # 单例模式装饰器
from femdb.LoadCase import *          # 加载工况相关类
from femdb.GlobalFEMVariant import ModelInfo  # 模型全局信息
from collections import OrderedDict    # 有序字典
from scipy import sparse               # 稀疏矩阵处理

# 使用单例模式装饰器，确保只有一个数据库实例
@Singleton
class FEMDataBase(object):
    """
    有限元数据库, 实例化的都存储在这里
    用于存储和管理有限元分析中的所有数据
    """

    def __init__(self):
        # 输入文件路径
        self.file_path = None

        # 节点和单元相关数据
        self.node_list = []  # 存储所有节点的列表，包含节点实例
        self.node_hash = {}  # 节点真实ID到node_list索引的映射表
        self.ele_hash = {}   # 单元真实ID到elements列表索引的映射表
        self.elements = []   # 存储所有单元的列表
        self.equation_number = None  # 方程数量(自由度数量)
        self.node_connected_element_count = None  # 节点连接单元计数

        # 节点集和单元集
        self.element_set_name_hash = {}  # 单元集名称到索引的映射
        self.node_set_name_hash = {}    # 节点集名称到索引的映射
        self.element_sets = []          # 存储所有单元集
        self.node_sets = []             # 存储所有节点集

        # 约束方程相关
        self.equation_constrain_couple = []  # 约束方程耦合项
        self.equation_constrain_idx = []     # 约束方程索引

        # 预处理相关数据
        self.matrix_num_count = 0       # 矩阵数量计数
        self.stiff_list = []            # 刚度矩阵列表
        self.global_stiff_matrix = None  # 全局刚度矩阵
        self.global_mass_matrix = None   # 全局质量矩阵
        self.load_case = LoadCase()     # 加载工况实例

        # 幅值曲线
        self.amplitudes = {}

        # 结果数据
        self.linear_u = None       # 线性位移结果
        self.linear_mises = None   # 线性Mises应力结果

        # 应力分量
        self.sigma_xx = None  # x方向正应力
        self.sigma_yy = None  # y方向正应力
        self.sigma_zz = None  # z方向正应力
        self.tau_xy = None    # xy剪应力
        self.tau_xz = None    # xz剪应力
        self.tau_yz = None    # yz剪应力

        # 历史结果数据
        self.history_u = None  # 位移历史
        self.history_v = None  # 速度历史
        self.history_a = None  # 加速度历史
        self.history_step_count = 0  # 历史步计数

        # 绘图相关
        self.additional_elements = []  # 附加单元(用于绘图等)

    """ 
    以下是模型解析相关函数，用于添加节点、单元、节点集、单元集、属性、材料、边界条件、LoadCase等 
    """

    def AddNode(self, node):
        """
        向有限元模型中插入节点
        :param node: 要添加的节点实例
        """
        self.node_list.append(node)

    def GetNodeBySearchId(self, n_id):
        """
        通过搜索ID获取节点
        :param n_id: 节点在node_list中的索引
        :return: 对应的节点实例
        """
        return self.node_list[n_id]

    def GetModelSummary(self):
        """
        获取模型摘要信息
        包括文件位置、分析类型、单元个数、节点个数、EquationNumber等
        :return: 包含摘要信息的OrderedDict
        """
        summary_dict = OrderedDict()
        summary_dict["File Path"] = str(self.file_path)
        # 计算方程数量(自由度数量)
        if not self.equation_number:
            summary_dict["Number of Equation"] = ModelInfo.PER_NODE_DOF * len(self.node_list)
        else:
            summary_dict["Number Of Equation"] = self.equation_number
        # 节点和单元数量
        summary_dict["Number Of Node"] = len(self.node_list)  # TODO: 并不是标准的节点个数, 标准节点个数应该是由单元计算出来
        summary_dict["Number Of Element"] = len(self.elements)

        return summary_dict

    def GetNodeCoordBySearchId(self, indexes):
        """
        通过节点索引获取节点坐标
        用于计算单元刚度矩阵
        :param indexes: 节点索引列表
        :return: 包含坐标的列表，每个元素是(x,y,z)元组
        """
        coords = []
        for idx in indexes:
            t_node = self.node_list[idx]
            coords.append((t_node.x, t_node.y, t_node.z))
        return coords

    def InitAssemblyMatrix(self, eq_count):
        """
        初始化总刚度矩阵
        :param eq_count: 方程数量(矩阵维度)
        """
        # 创建空的COO格式稀疏矩阵
        self.global_stiff_matrix = sparse.coo_matrix((eq_count, eq_count), dtype=float)