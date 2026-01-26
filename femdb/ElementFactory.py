#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys

# 0D Elements
from element.ElementBase import ElementBaseClass
from element.Mass import Mass

# 1D Elements
from element.Truss import T3D2
from element.Beam import Beam188, Beam189

# 2D Elements
from element.Shell import CookQuaShell, CookTriShell
from element.MITC4 import MITC4Shell
from element.MITC3 import MITC3Shell
from element.Plane import CPS3, CPS4

# 3D Elements
from element.Tetra import C3D4
from element.Wedge import C3D6
from element.Hexa import C3D8

from femdb.GlobalEnum import *
from femdb.GlobalFEMVariant import ModelInfo

from typing import Tuple


class ElementFactory:
    """
    Reference:
    1. https://abaqus-docs.mit.edu/2017/English/SIMACAEELMRefMap/simaelm-c-shellelem.htm
    """

    @staticmethod
    def CreateElement(e_type, e_id=-1, opt=None) -> Tuple[ElementBaseClass, int, int]:
        """
        静态函数, 用于返回
        :param e_type: 单元类型，这里包含了Abaqus、Nastran和Ansys的
        :param e_id: 初始化单元需要单元ID
        :param opt: 附加参数, 比如181可能是3节点壳也可能是4节点壳, solid45可能是8节点也可能是4节点
        :return: 单元、节点个数以、单刚中包含数据个数以及单元涉及节点自由度的大小
        """
        # 0D Element
        if e_type in [21]:
            # return Mass(e_id), 1, (ModelInfo.PER_NODE_DOF + 1) * ModelInfo.PER_NODE_DOF // 2
            return Mass(e_id), 1, 3

        # 1D Element
        elif e_type in ["T3D2"]:
            return T3D2(e_id), 2, 3
        elif e_type in ["B31", 188]:
            # return Beam188(e_id), 2, 78
            return Beam188(e_id), 2, 6
        elif e_type in [189]:
            # return Beam189(e_id), 3, 171
            return Beam189(e_id), 3, 6

        # 2D Element
        elif e_type in ["CPS3"]:
            # return CPS3(e_id), 3, 21
            return CPS3(e_id), 3, 2
        elif e_type in ["CPS4"]:
            # return CPS4(e_id), 4, 36
            return CPS4(e_id), 4, 2

        # 3D Element
        elif e_type in ["S3", "S3R"]:
            # return CookTriShell(e_id), 3, 171
            return CookTriShell(e_id), 3, 6
            # return MITC3Shell(e_id), 3, 6
        elif e_type in ["S4", "S4R", "S4RT"]:
            # return CookQuaShell(e_id), 4, 300
            return CookQuaShell(e_id), 4, 6
            # return MITC4Shell(e_id), 4, 6
        elif e_type in [181, 63, 131]:
            if opt == 4:
                # return CookQuaShell(e_id), 4, 300
                return CookQuaShell(e_id), 4, 6
                # return MITC4Shell(e_id), 4, 6
            elif opt == 3:
                # return CookTriShell(e_id), 3, 171
                return CookTriShell(e_id), 3, 6
                # return MITC3Shell(e_id), 3, 6
            else:
                raise KeyError("Shell 181/63 don't support opt {}".format(opt))

        elif e_type in ["C3D8", "C3D8R"]:
            return C3D8(e_id), 8, 3
        # elif e_type in ["C3D8R"]:
        #     mlogger.fatal("No impl such element")
        #     sys.exit(1)
        elif e_type in ["C3D6"]:
            # return C3D6(e_id), 6, 171
            return C3D6(e_id), 6, 3
        elif e_type in ["C3D4"]:
            # return C3D4(e_id), 4, 78
            return C3D4(e_id), 4, 3
        elif e_type in ["C3D20R"]:
            mlogger.fatal("No impl such element")
            sys.exit(1)
        elif e_type in [185, 45, 70]:
            if opt == 8:
                # return C3D8(e_id), 8, 300
                return C3D8(e_id), 8, 3
            elif opt == 6:
                # return C3D6(e_id), 6, 171
                return C3D6(e_id), 6, 3
            elif opt == 4:
                # return C3D4(e_id), 4, 78
                return C3D4(e_id), 4, 3
            else:
                raise KeyError("Wrong opt parameter: {}".format(opt))

        raise KeyError("Fatal Error: No Such ElementType: {}".format(e_type))

    @staticmethod
    def GetElementNodeDofCount(e_type):
        """
        返回组成单元的节点所需的自由度个数
        :param e_type: 对于Abaqus是字符串, 对于ANSYS是int
        :return: 单元的自由度个数
        """
        # 0D Element
        if e_type in [21]:
            return ModelInfo.PER_NODE_DOF

        # 1D Element
        elif e_type in ["T3D2"]:
            return 3
        elif e_type in ["B31", 188]:
            return 6

        # 2D Element
        elif e_type in ["S3", "S3R"]:
            return 6
        elif e_type in ["S4", "S4R", "S4RT"]:
            return 6
        elif e_type in ["CPS3", "CPS4"]:
            return 2

        # 3D Element
        elif e_type in ["C3D8", 185, "C3D6", "C3D4", 45, 70]:
            return 3
        elif e_type in ["C3D8R"]:
            return 3
        elif e_type in ["C3D20R"]:
            return 3
        elif e_type in [181, 63, 131]:
            return 6

        raise KeyError("No Such ElementType: {}".format(e_type))
