#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from element.ElementBase import *
from abc import ABC
from femdb.GlobalFEMVariant import ModelInfo
import numpy as np


class Mass(ElementBaseClass, ABC):
    """
    Mass Element class
    """

    def __init__(self, eid):
        super().__init__(eid)
        self.nodes_count = 1  # Each element has 1 node
        self.vtu_type = "vertex"
        self.stiffness = None
        self.stress = None

    def CalElementDMatrix(self, an_type=None):
        """
        质量单元无需计算D阵
        """
        pass

    def ElementStiffness(self):
        """
        Reference:
        """
        assert self.node_coords.shape == (1, 3)
        eps = 1e-7
        return eps * np.eye(ModelInfo.PER_NODE_DOF)

    def CalculateElementStress(self, displacement):
        """
        Calculate element stress, 刚体没有形变, 所以没有应力
        """
        return np.zeros(1), np.zeros(1), np.zeros(1), np.zeros(1), np.zeros(1), np.zeros(1)

    def ElementMass(self):
        pass

    def CalculateBasic(self):
        pass


if __name__ == "__main__":
    t_ele = Mass(-1)
    t_ele.ele_mat_dict = {MaterialKey.E: 1, MaterialKey.Area: np.sqrt(3)}
    t_ele.node_coords = np.array([[0, 0, 0],
                                  [1, 1, 1]], dtype=float)
    print(t_ele.ElementStiffness())
    mlogger.debug("finish")
