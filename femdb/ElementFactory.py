#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# 1D Elements
from element.Truss import *
from element.Beam import *

# 2D Elements
from element.Shell import *
from element.Plane import *

# 3D Elements
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
        :param opt: 附加参数, 比如181可能是3节点壳也可能是4节点壳
        :return: 单元和节点个数
        """
        # 1D Element
        if e_type in ["T3D2"]:
            return T3D2(e_id), 2
        elif e_type in ["B31", 188]:
            return B3D2(e_id), 2

        # 2D Element
        elif e_type in ["S3"]:
            return Shell3(e_id), 3
        elif e_type in ["S4", "S4R", "S4RT"]:
            return Shell4(e_id), 4
        elif e_type in [181]:
            if opt == 3:
                return MITC3(e_id), 3
            elif opt == 4:
                return MITC4(e_id), 4
        elif e_type in ["CPS3"]:
            return MITC3(e_id), 3
        elif e_type in ["CPS4"]:
            return MITC4(e_id), 4

        # 3D Element
        elif e_type in ["C3D8", 45]:
            return C3D8(e_id), 8
        elif e_type in ["C3D8R"]:
            # return C3D8R(e_id), 8
            mlogger.fatal("No impl such element")
        elif e_type in ["C3D6"]:
            return Wedge(e_id), 6
        elif e_type in ["C3D4"]:
            return TetraElement(e_id), 4
        elif e_type in ["C3D20R"]:
            mlogger.fatal("No impl such element")

        mlogger.fatal("Fatal Error: No Such ElementType: {}".format(e_type))
        sys.exit(1)

    @staticmethod
    def GetElementNodeDofCount(e_type):
        """
        返回组成单元的节点所需的自由度个数
        """
        # 1D Element
        if e_type in ["T3D2"]:
            return 3
        elif e_type in ["B31"]:
            return 6

        # 2D Element
        elif e_type in ["S3"]:
            return 6
        elif e_type in ["S4", "S4R", "S4RT"]:
            return 6
        elif e_type in ["CPS3", "CPS4"]:
            return 2

        # 3D Element
        elif e_type in ["C3D8"]:
            return 3
        elif e_type in ["C3D8R"]:
            return 3
        elif e_type in ["C3D6"]:
            return 3
        elif e_type in ["C3D4"]:
            return 3
        elif e_type in ["C3D20R"]:
            return 3

        mlogger.fatal("No Such ElementType: {}".format(e_type))
