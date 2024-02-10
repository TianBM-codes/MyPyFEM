#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from femdb.LoadCase import LoadCase
from femdb.Material import Material
from element.Node import Node
from femdb.Sets import NodeSet, EleSet
from femdb.Property import Property
from femdb.Section import *
from femdb.ElementGroup import *
from femdb.Mesh import Mesh
from femdb.Geom import Geom
from femdb.Kinematics import Kinematics
from typing import Dict, List
from utils.CustomException import *
from femdb.Boundary import BoundaryBase
from SolveControl import SolveControl


class NLFEMDataBase(object):
    """
    有限元数据库, 实例化的都存储在这里
    """

    _instance = None  # 类变量用于存储唯一的实例

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(NLFEMDataBase, cls).__new__(cls)
        return cls._instance

    def __init__(self):
        self.file_path = None
        self.title = None
        self.Mesh = Mesh()
        self.Geom = Geom()
        self.LoadCase = LoadCase()
        self.BC = BoundaryBase()
        self.Material = Material()
        self.Kinematics = Kinematics()
        self.Dimension = AnalyseDimension.NoAssign
        self.SolveControl = SolveControl()
        self.ElementGroupHash: Dict[int, ElementGroup] = {}
