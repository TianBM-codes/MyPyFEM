# !/usr/bin/env python3
# -*- coding: utf-8 -*-

from element.Node import Node
from element.ElementBase import AllEleTypeDNDrAtGaussianPoint
from utils.GlobalEnum import *
from typing import List
import numpy as np


class Geom(object):
    """
    存储几何信息
    """

    def __init__(self):
        self.node_count = -1
        self.node_list: List[Node] = []
        self.node_hash = {}  # 节点真实Id对应nodelist中的index的Hash表
        self.x0 = None
        self.x = None
        self.V_ele = None
        self.V_total = -1

    def AddNode(self, node: Node, node_id: int):
        """
        向系统中添加节点，
        @param node: 文本中读取的节点
        @param node_id: 导入文件中的编号
        @return:
        """
        self.node_list.append(node)
        self.node_hash[node_id] = len(self.node_list)

    def InitX(self):
        self.x0 = np.zeros((self.node_count, GlobalInfor[GlobalVariant.Dimension]))
        for ii in range(self.node_count):
            self.x0[ii, :] = self.node_list[ii].coord

        self.x = np.zeros((self.node_count, GlobalInfor[GlobalVariant.Dimension]))
        for ii in range(self.node_count):
            self.x[ii, :] = self.node_list[ii].coord

    def InitialVolume(self):
        """
        Calculate initial volume for data checking.
        Additionally, essential for mean dilation algorithm.
        @return:
        """
        from femdb.NLFEMDataBase import NLFEMDataBase
        fem_db = NLFEMDataBase()
        self.V_ele = np.zeros((fem_db.Mesh.nelem, 1))
        self.V_total = 0
        for _, ele_grp in fem_db.Mesh.ele_grp_hash.items():
            eles = ele_grp.eles
            for jj in range(len(eles)):
                ele = eles[jj]
                nodes = ele.nodes()
                x_local = []
                for nd in nodes:
                    iter_node = self.node_list[self.node_hash[nd]]
                    x_local.append(iter_node.GetNodeCoord())
                x_local = np.asarray(x_local)
                weight, DN_Dchi = AllEleTypeDNDrAtGaussianPoint.GetElementDNDchi(ele.e_type)
                fem_db.Kinematics.ComputeGradients(x_local, x_local, DN_Dchi)
                for ii in ele.ngauss:
                    jw = fem_db.Kinematics.Jx_chi[ii] * weight[ii]
                    self.V_ele[jj] += jw
                self.V_total += self.V_ele[jj]
