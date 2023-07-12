#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys

# 1D Elements
from element.Truss import *
from element.Beam import *

# 2D Elements
from element.Shell import *
from element.Plane import *

# 3D Elements
from element.Plate import *
from element.Tetra import *
from element.Wedge import *
from element.Hexa import *


class ElementFactory:
    """
    Reference:
    1. https://abaqus-docs.mit.edu/2017/English/SIMACAEELMRefMap/simaelm-c-shellelem.htm
    """

    @staticmethod
    def CreateElement(e_type, e_id=-1, opt=None) -> tuple[ElementBaseClass, int]:
        """
        静态函数, 用于返回
        :param e_type: 单元类型，这里包含了Abaqus、Nastran和Ansys的
        :param e_id: 初始化单元需要单元ID
        :param opt: 附加参数, 比如181可能是3节点壳也可能是4节点壳, solid45可能是8节点也可能是4节点
        :return: 单元和节点个数
        """
        # 1D Element
        if e_type in ["T3D2"]:
            return T3D2(e_id), 2
        elif e_type in ["B31", 188]:
            return Beam188(e_id), 2
        elif e_type in [189]:
            return Beam189(e_id), 3

        # 2D Element
        elif e_type in ["CPS3"]:
            return CPS3(e_id), 3
        elif e_type in ["CPS4"]:
            return CPS4(e_id), 4

        # 3D Element
        elif e_type in ["S3"]:
            return DKTShell(e_id), 3
        elif e_type in ["S4", "S4R", "S4RT"]:
            return DKQShell(e_id), 4
        elif e_type in [181]:
            return DKQShell(e_id), 4

        elif e_type in ["C3D8", 45]:
            return C3D8(e_id), 8
        elif e_type in ["C3D8R"]:
            mlogger.fatal("No impl such element")
            sys.exit(1)
        elif e_type in ["C3D6"]:
            return C3D6(e_id), 6
        elif e_type in ["C3D4"]:
            return C3D4(e_id), 4
        elif e_type in ["C3D20R"]:
            mlogger.fatal("No impl such element")
            sys.exit(1)
        elif e_type == 185:
            if opt==8:
                return C3D8(e_id), 8
            elif opt==6:
                return C3D6(e_id), 6
            elif opt==4:
                return C3D4(e_id), 4
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
        elif e_type in [181]:
            return 6

        mlogger.fatal("No Such ElementType: {}".format(e_type))
        sys.exit(1)


