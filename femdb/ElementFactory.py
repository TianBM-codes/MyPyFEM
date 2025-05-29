#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys

# 1D Elements
from element.ElementBase import ElementBaseClass
from element.Truss import T3D2
from element.Beam import Beam188, Beam189

# 2D Elements
from element.Shell import CookQuaShell, CookTriShell
from element.Plane import CPS3, CPS4

# 3D Elements
from element.Tetra import C3D4
from element.Wedge import C3D6
from element.Hexa import C3D8

from femdb.GlobalEnum import *


class ElementFactory:
    """
    Reference:
    1. https://abaqus-docs.mit.edu/2017/English/SIMACAEELMRefMap/simaelm-c-shellelem.htm
    """

    @staticmethod
    def CreateElement(e_type, e_id=-1, opt=None) -> tuple[ElementBaseClass, int, int]:
        """
        静态函数, 用于返回
        :param e_type: 单元类型，这里包含了Abaqus、Nastran和Ansys的
        :param e_id: 初始化单元需要单元ID
        :param opt: 附加参数, 比如181可能是3节点壳也可能是4节点壳, solid45可能是8节点也可能是4节点
        :return: 单元、节点个数以及单刚中包含数据个数
        """
        # 0D Element
        if e_type in [21]:
            return
        # 1D Element
        elif e_type in ["T3D2"]:
            return T3D2(e_id), 2, 36
        elif e_type in ["B31", 188]:
            return Beam188(e_id), 2, 144
        elif e_type in [189]:
            return Beam189(e_id), 3, 324

        # 2D Element
        elif e_type in ["CPS3"]:
            return CPS3(e_id), 3, 36
        elif e_type in ["CPS4"]:
            return CPS4(e_id), 4, 64

        # 3D Element
        elif e_type in ["S3"]:
            # return TriangleShell63(e_id), 3
            return CookTriShell(e_id), 3, 324
        elif e_type in ["S4", "S4R", "S4RT"]:
            # return QuadShell63(e_id), 4
            return CookQuaShell(e_id), 4, 576
        elif e_type in [181, 63]:
            if opt == 4:
                # return QuadShell63(e_id), 4
                return CookQuaShell(e_id), 4, 576
            elif opt == 3:
                # return TriangleShell63(e_id), 3
                return CookTriShell(e_id), 3, 324
            else:
                mlogger.fatal("Shell 181/63 don't support opt {}".format(opt))
                sys.exit(1)

        elif e_type in ["C3D8", 45]:
            return C3D8(e_id), 8, 576
        elif e_type in ["C3D8R"]:
            mlogger.fatal("No impl such element")
            sys.exit(1)
        elif e_type in ["C3D6"]:
            return C3D6(e_id), 6, 324
        elif e_type in ["C3D4"]:
            return C3D4(e_id), 4, 144
        elif e_type in ["C3D20R"]:
            mlogger.fatal("No impl such element")
            sys.exit(1)
        elif e_type == 185:
            if opt == 8:
                return C3D8(e_id), 8, 576
            elif opt == 6:
                return C3D6(e_id), 6, 324
            elif opt == 4:
                return C3D4(e_id), 4, 144
            else:
                mlogger.fatal("Wrong opt parameter: {}".format(opt))
                sys.exit(1)

        mlogger.fatal("Fatal Error: No Such ElementType: {}".format(e_type))
        sys.exit(1)

    @staticmethod
    def GetElementNodeDofCount(e_type):
        """
        返回组成单元的节点所需的自由度个数
        :param e_type: 对于Abaqus是字符串, 对于ANSYS是int
        :return: 单元的自由度个数
        """
        # 1D Element
        if e_type in ["T3D2"]:
            return 3
        elif e_type in ["B31", 188]:
            return 6

        # 2D Element
        elif e_type in ["S3"]:
            return 6
        elif e_type in ["S4", "S4R", "S4RT"]:
            return 6
        elif e_type in ["CPS3", "CPS4"]:
            return 2

        # 3D Element
        elif e_type in ["C3D8", 185, "C3D6", "C3D4"]:
            return 3
        elif e_type in ["C3D8R"]:
            return 3
        elif e_type in ["C3D20R"]:
            return 3
        elif e_type in [181, 63]:
            return 6

        raise KeyError("No Such ElementType: {}".format(e_type))
