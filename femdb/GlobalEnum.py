#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from enum import Enum
import logging
import matplotlib.pyplot as plt

# 配置logging, handler可以添加控制台, 即可实现控制台和文件同时输出
logging.basicConfig(level=logging.DEBUG)
handler = logging.FileHandler(r"D:\WorkSpace\FEM\MyPyFEM\analyse_output.log", encoding="utf-8", mode="a")
formatter = logging.Formatter('%(asctime)s - %(funcName)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
handler.setFormatter(formatter)

mlogger = logging.getLogger(__name__)
mlogger.addHandler(handler)

# 配置matplot, 关闭其logging
logging.getLogger('matplotlib.font_manager').disabled = True
plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号

mlogger.info("\n{} Analysis Calculate Begin {}".format("#" * 6, "#" * 6))

Abaqus2VTKType = {
    # trusses
    "T2D2": "line",
    "T2D2H": "line",
    "T2D3": "line3",
    "T2D3H": "line3",
    "T3D2": "line",
    "T3D2H": "line",
    "T3D3": "line3",
    "T3D3H": "line3",
    # beams
    "B21": "line",
    "B21H": "line",
    "B22": "line3",
    "B22H": "line3",
    "B31": "line",
    "B31H": "line",
    "B32": "line3",
    "B32H": "line3",
    "B33": "line3",
    "B33H": "line3",
    # surfaces
    "CPS4": "quad",
    "CPS4R": "quad",
    "S4": "quad",
    "S4R": "quad",
    "S4RS": "quad",
    "S4RSW": "quad",
    "S4R5": "quad",
    "S4RT": "quad",
    "S8R": "quad8",
    "S8R5": "quad8",
    "S9R5": "quad9",
    # "QUAD": "quad",
    # "QUAD4": "quad",
    # "QUAD5": "quad5",
    # "QUAD8": "quad8",
    # "QUAD9": "quad9",
    #
    "CPS3": "triangle",
    "STRI3": "triangle",
    "S3": "triangle",
    "S3R": "triangle",
    "S3RS": "triangle",
    "R3D3": "triangle",
    # "TRI7": "triangle7",
    # 'TRISHELL': 'triangle',
    # 'TRISHELL3': 'triangle',
    # 'TRISHELL7': 'triangle',
    #
    "STRI65": "triangle6",
    # 'TRISHELL6': 'triangle6',
    # volumes
    "C3D8": "hexahedron",
    "C3D8H": "hexahedron",
    "C3D8I": "hexahedron",
    "C3D8IH": "hexahedron",
    "C3D8R": "hexahedron",
    "C3D8RH": "hexahedron",
    # "HEX9": "hexahedron9",
    "C3D20": "hexahedron20",
    "C3D20H": "hexahedron20",
    "C3D20R": "hexahedron20",
    "C3D20RH": "hexahedron20",
    # "HEX27": "hexahedron27",
    #
    "C3D4": "tetra",
    "C3D4H": "tetra4",
    # "TETRA8": "tetra8",
    "C3D10": "tetra10",
    "C3D10H": "tetra10",
    "C3D10I": "tetra10",
    "C3D10M": "tetra10",
    "C3D10MH": "tetra10",
    # "TETRA14": "tetra14",
    #
    # "PYRAMID": "pyramid",
    "C3D6": "wedge",
    "C3D15": "wedge15",
    #
    # 4-node bilinear displacement and pore pressure
    "CAX4P": "quad",
    # 6-node quadratic
    "CPE6": "triangle6",
}

Ansys2VTKType = {
    # trusses
    # solid
    45: "hexahedron",
    # beams
    188: "line",
    # shell
    181: "quad",
}


class MaterialKey(Enum):
    """ 材料参数关键字集合 """
    E = 1
    Density = 2
    Niu = 3


class PropertyKey(Enum):
    """
    属性关键字集合, 从100开始是因为需要与MaterialKey相结合, 不能发生重复
    """
    ThicknessOrArea = 100


class FEMObject(Enum):
    """ 有限元相关类 """
    NodeSet = 1
    EleSet = 2
    Material = 3


class MaterialMatrixType(Enum):
    """
    对于平面应力和平面应变问题, D阵是不同的, 其他都是一致的
    """
    PlaneStree = 1
    PlaneStrain = 2


class AnalyseDimension(Enum):
    """
    分析类型, 暂时区分2D还是3D
    """
    TwoDimension = 1
    ThreeDimension = 2


class AnalyseType(Enum):
    """
    分析类型
    """
    LinearStatic = 1
    ModalAnalyse = 2


class InputFileType(Enum):
    """
    导入初始化数据库的文件类型
    """
    CDB = 1
    INP = 2
    BDF = 3


class BeamSectionType(Enum):
    """
    梁横截面类型
    """
    GongZiGang = 1
    Circle = 2
    Square = 3
