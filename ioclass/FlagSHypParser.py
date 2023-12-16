#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
from femdb.NLFEMDataBase import NLFEMDataBase
from element.Node import Node
from CustomException import *
from ElementFactory import ElementFactory


class FlagSHypReader(object):
    """
    读取FlagSHy的录入文件
    """

    def __init__(self, input_path):
        """
        """
        self.fem_database = NLFEMDataBase()
        self.bbb = NLFEMDataBase()
        # self.fem_database.G
        self.fem_database.Geom = None
        self.dat_path = input_path
        self.iter_line = None
        self.et_hash = {}
        self.ele_group_hash = {}
        self.ele_count = 0

    def ParseFileAndInitFEMDB(self):
        """
        文件格式：
        Javier Book P263
        :return:
        """
        fem_db = self.fem_database
        with open(self.dat_path, 'r', encoding='utf-8') as dat_file:
            """
            读取项目基础信息, 包括项目名称, 以及算例中涉及的单元
            """
            fem_db.title = dat_file.readline()
            ele_type = dat_file.readline()

            """
            读取节点信息
            """
            fem_db.Geom.npoin = int(dat_file.readline())
            for ii in range(fem_db.Geom.npoin):
                node_line = dat_file.readline().strip().split(" ")

                # 三维问题，每个节点有xyz三个坐标
                if len(node_line) == 5:
                    Node(int(node_line[0]))
                elif len(node_line) == 4:
                    Node(int(node_line[0]))
                else:
                    raise InputTextFormat(node_line)

            """
            读取单元信息
            """
            fem_db.Mesh.element_count = int(dat_file.readline())
            for ii in range(fem_db.Mesh.element_count):
                ele, ele_node_count = ElementFactory.CreateElement(ele_type)
                ele_line = dat_file.readline().strip().split()
                # fem_db.Mat.


if __name__ == "__main__":
    kkt = FlagSHypReader("../NumericalCases/TwistingColumn/twisting_column.dat")
    kkt.ParseFileAndInitFEMDB()
