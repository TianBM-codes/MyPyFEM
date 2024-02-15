#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import abc
from femdb.Integration import *
from utils.UtilsFunction import *


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
        self.node_ids = np.asarray([], dtype=np.int32)  # Node list of the element
        self.search_node_ids = np.asarray([])  # Domain中node_list中该节点的index
        self.cha_dict = None  # 单元的属性字典, 其中包括材料、属性、常数、惯性矩等
        self.vtp_type = None
        self.unv_code = None  # SiPESC平台显示的UNV结果, 单元代号
        self.eq_numbers = np.asarray([], dtype=np.uint32)  # 方程号, 即在求解矩阵中的第几行, 也即自由度排序后的index
        self.D = None  # 本构矩阵
        self.B = None  # 应变矩阵, 用于求解应力
        self.DN_Dchi = None

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
    def ElementStress(self, displacement: np.array):
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


"""
以下为计算各种单元形函数对参数坐标的导数
"""


def CalculateC3D8():
    """
    Reference:
    1. B Bar method: <<The Finite Element Method Linear Static and Dynamic Finite Element Analysis(Thomas J.R.Hughes)>> P232
    2. <<有限元法>> Bath P323

    注意节点编号规则(逆时针,右手螺旋向z正向), 即参数坐标对应具体哪个真实节点
    # Shape Function:
    N1 = (1 + r) * (1 + s) * (1 + t) / 8
    N2 = (1 - r) * (1 + s) * (1 + t) / 8
    N3 = (1 - r) * (1 - s) * (1 + t) / 8
    N4 = (1 + r) * (1 - s) * (1 + t) / 8
    N5 = (1 + r) * (1 + s) * (1 - t) / 8
    N6 = (1 - r) * (1 + s) * (1 - t) / 8
    N7 = (1 - r) * (1 - s) * (1 - t) / 8
    N8 = (1 + r) * (1 - s) * (1 - t) / 8

    # Partial
    dN1dr, dN1ds, dN1dt = (1 + s) * (1 + t) / 8, (1 + r) * (1 + t) / 8, (1 + r) * (1 + s) / 8
    dN2dr, dN2ds, dN2dt = -(1 + s) * (1 + t) / 8, (1 - r) * (1 + t) / 8, (1 - r) * (1 + s) / 8
    dN3dr, dN3ds, dN3dt = (s - 1) * (1 + t) / 8, (r - 1) * (1 + t) / 8, (1 - r) * (1 - s) / 8
    dN4dr, dN4ds, dN4dt = (1 - s) * (1 + t) / 8, -(1 + r) * (1 + t) / 8, (1 + r) * (1 - s) / 8
    dN5dr, dN5ds, dN5dt = (1 + s) * (1 - t) / 8, (1 + r) * (1 - t) / 8, -(1 + r) * (1 + s) / 8
    dN6dr, dN6ds, dN6dt = (1 + s) * (t - 1) / 8, (1 - r) * (1 - t) / 8, (r - 1) * (1 + s) / 8
    dN7dr, dN7ds, dN7dt = (s - 1) * (1 - t) / 8, (r - 1) * (1 - t) / 8, (r - 1) * (1 - s) / 8
    dN8dr, dN8ds, dN8dt = (1 - s) * (1 - t) / 8, (1 + r) * (t - 1) / 8, (1 + r) * (s - 1) / 8
    """
    # Gaussian Weight
    sample_pt, weight = GaussIntegrationPoint.GetSamplePointAndWeight(2)

    # 在8个高斯点上积分
    dNdrs = []
    weights = []
    for ri in range(2):
        for si in range(2):
            for ti in range(2):
                r, s, t = sample_pt[ri], sample_pt[si], sample_pt[ti]
                dNdr = 0.125 * np.asarray([[(1 + s) * (1 + t) / 8, (1 + r) * (1 + t) / 8, (1 + r) * (1 + s) / 8],
                                           [-(1 + s) * (1 + t) / 8, (1 - r) * (1 + t) / 8, (1 - r) * (1 + s) / 8],
                                           [(s - 1) * (1 + t) / 8, (r - 1) * (1 + t) / 8, (1 - r) * (1 - s) / 8],
                                           [(1 - s) * (1 + t) / 8, -(1 + r) * (1 + t) / 8, (1 + r) * (1 - s) / 8],
                                           [(1 + s) * (1 - t) / 8, (1 + r) * (1 - t) / 8, -(1 + r) * (1 + s) / 8],
                                           [(1 + s) * (t - 1) / 8, (1 - r) * (1 - t) / 8, (r - 1) * (1 + s) / 8],
                                           [(s - 1) * (1 - t) / 8, (r - 1) * (1 - t) / 8, (r - 1) * (1 - s) / 8],
                                           [(1 - s) * (1 - t) / 8, (1 + r) * (t - 1) / 8, (1 + r) * (s - 1) / 8]]).T
                weights.append(weight[ri] * weight[si] * weight[ti])
                dNdrs.append(dNdr)

    return dNdrs, weights


def CalculateC3D6():
    """
    Reference:
    1. https://www.help.febio.org/FEBio/FEBio_tm_2_7/FEBio_tm_2-7-Subsection-4.1.2.html#:~:text=Pentahedral%20elements%20%28also%20knows%20as%20%E2%80%9Cwedge%E2%80%9D%20elements%29%20consist,s%20and%20t%20and%20are%20given%20as%20follows.
    2. https://github.com/febiosoftware

    # Shape Function:
    N1 = 0.5 * (1 - r - s) * (1 - t)
    N2 = 0.5 * r * (1 - t)
    N3 = 0.5 * s * (1 - t)
    N4 = 0.5 * (1 - r - s) * (1 + t)
    N5 = 0.5 * r * (1 + t)
    N6 = 0.5 * s * (1 + t)

    # Partial
    pN1pr, pN1ps, pN1pt = 0.5 * (t - 1), 0.5 * (t - 1), 0.5 * (r + s - 1)
    pN2pr, pN2ps, pN2pt = 0.5 * (1 - t), 0, -0.5 * r
    pN3pr, pN3ps, pN3pt = 0, 0.5 * (1 - t), -0.5 * s
    pN4pr, pN4ps, pN4pt = -0.5 * (1 + t), -0.5 * (1 + t), 0.5 * (1 - r - s)
    pN5pr, pN5ps, pN5pt = 0.5 * (1 + t), 0, 0.5 * r
    pN6pr, pN6ps, pN6pt = 0, 0.5 * (1 + t), 0.5 * s
    """
    sample_r = [0.166666667, 0.666666667, 0.166666667, 0.166666667, 0.666666667, 0.166666667]
    sample_s = [0.166666667, 0.166666667, 0.666666667, 0.166666667, 0.166666667, 0.666666667]
    sample_t = [-0.577350269, -0.577350269, -0.577350269, 0.577350269, 0.577350269, 0.577350269]
    weight = 0.166666667

    # 在6个高斯点上积分
    dNdrs = []
    for ii in range(6):
        r, s, t = sample_r[ii], sample_s[ii], sample_t[ii]
        dNdrs.append(np.asarray([[0.5 * (t - 1), 0.5 * (t - 1), 0.5 * (r + s - 1)],
                                 [0.5 * (1 - t), 0, -0.5 * r],
                                 [0, 0.5 * (1 - t), -0.5 * s],
                                 [-0.5 * (1 + t), -0.5 * (1 + t), 0.5 * (1 - r - s)],
                                 [0.5 * (1 + t), 0, 0.5 * r],
                                 [0, 0.5 * (1 + t), 0.5 * s]], dtype=float).T)
    return dNdrs, [weight] * 6


class DNDrCalculator:
    """
    计算各种单元的dNdr, 也就是同一种单元涉及到的相同的地方, 这样就不需要每个单元都计算一遍,
    同一种单元计算一次即可
    """

    def __init__(self):
        # solid part
        self.C3D8 = CalculateC3D8()
        self.C3D6 = CalculateC3D6()

        # shell part
        self.DKTShell = None
        self.DKQShell = None

    # def GetElementDNDchi(self, ele_type):
    #     if ele_type in ["hexa"]:
    #         return self.C3D8


# AllEleTypeDNDrAtGaussianPoint = DNDrCalculator()
