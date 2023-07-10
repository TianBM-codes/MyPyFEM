#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
    SiPESC平台单元类型对应的code
    xx-yy-zz
    xx：拓扑类型，节点数
    yy：单元分类，1梁，2杆，3膜，4板，5壳，6体，7轴对称壳，8轴对称四边形，9夹层板壳，10层合板壳，11柳钉，12连接
    zz：单元子类。
    TYPE=10000,    1点三向拉压弹簧单元;
    TYPE=11100,    1点铆钉单元;
    TYPE=20100,    2点偏心梁单元;
    TYPE=20200,    2点轴力杆单元;
    TYPE=20210,    2点弹塑性轴力杆单元;
    TYPE=21200，   2点平面接触单元;
    TYPE=21201，   2点空间接触单元;
    TYPE=30100,    3点弯管单元;
    TYPE=30300,    3点平面应力膜单元;
    TYPE=30310,    3点三角形弹塑性平面应力膜单元 ;
    TYPE=30400,    3点三角形薄板单元;
    TYPE=30500,    3点三角形薄壳单元;
    TYPE=30501,    3点三角形各向异性薄壳单元;
    TYPE=30900,    3点三角形复合材料夹层板壳单元;
    TYPE=31000,    3点三角形复合材料层合板壳单元
    TYPE=40300,    4点平面应力膜单元;
    TYPE=40301,    4点平面应变膜单元;
    TYPE=40302,    4点非协调平面应力膜单元;
    TYPE=40303,    4点非协调平面应变膜单元;
    TYPE=40304,    4点矩形膜单元;
    TYPE=40310,    4点弹塑性平面应力膜单元;
    TYPE=40500,    4点任意四边形壳单元;
    TYPE=40501,    4点任意四边形各向异性壳单元;
    TYPE=40700,    4点轴对称旋转壳单元;
    TYPE=40800,    4点任意四边形4节点轴对称环体单元;
    TYPE=40900,    4点任意四边形复合材料夹层板壳单元;
    TYPE=41000,    4点任意四边形复合材料层合板壳单元;
    TYPE=50300,    5点等参平面膜单元;
    TYPE=50600,    5点金字塔单元;
    TYPE=60300,    6点高阶三角形单元;
    TYPE=60600,    6点三棱柱单元;
    TYPE=80300,    8点高阶四边形单元;
    TYPE=80600,    8点块体单元;
    TYPE=80601,    8点非协调块体单元;
    TYPE=80610,    8点弹塑性块体单元;
    TYPE=100600,   10点高阶四面体单元;
    TYPE=130600,   13点高阶金字塔单元;
    TYPE=150600,   15点高阶三棱柱单元;
    TYPE=200600,   20点高阶六面体单元;
"""
import abc
import numpy as np
from femdb.GlobalEnum import *
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
        self.nodes_count = -1  # Number of nodes per element
        self.node_ids = np.asarray([])  # Node list of the element
        self.search_node_ids = np.asarray([])  # Domain中node_list中该节点的index
        self.cha_dict = None  # 单元的属性字典, 其中包括材料、属性、常数、惯性矩等
        self.vtp_type = None
        self.unv_code = None  # SiPESC平台显示的UNV结果, 单元代号
        self.eq_numbers = np.asarray([], dtype=np.uint32)  # 方程号, 即在求解矩阵中的第几行, 也即自由度排序后的index
        self.D = None  # 本构矩阵
        self.is_degenerate_element = False

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
        # 判断单元是否为退化单元, 专门对于ANSYS来说的
        if len(list(set(nds))) != len(self.node_ids):
            self.is_degenerate_element = True

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
