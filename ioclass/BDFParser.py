#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from femdb.FEMDataBase import *


class Element:
    def __init__(self, eid, e_type, nds: list):
        self.eid = eid
        self.e_type = e_type
        self.nds = nds


class BDFParser(object):
    """
    功能为解析BDF文件, 初始化数据库单元节点及相关信息, 只负责解析和传递, 其他功能均在其他类中,
    同样的, 其他类也不会包含任何解析输入文件的功能或者函数
    解析Nastran软件的BDF文件, 主要功能为将pyNastran库读取的结果转化为自身程序适用的, 初始化数据库
    Reference:
    1. https://pynastran-git.readthedocs.io/en/1.3/quick_start/bdf_overview.html#example-1-read-write
    """

    def __init__(self, input_path):
        self.iter_line = None
        self.bdf_path = input_path
        self.nodes = []
        self.elements = []

    def ParseFile(self):
        with open(self.bdf_path, 'r') as bdf_f:
            self.iter_line = bdf_f.readline()
            while True:
                if self.iter_line.startswith("GRID"):
                    pt_id = int(self.iter_line[4:17])
                    x = float(self.iter_line[24:32])
                    y = float(self.iter_line[32:40])
                    z = float(self.iter_line[40:48])
                    self.nodes.append(Node(pt_id, x, y, z))
                    self.iter_line = bdf_f.readline()

                elif self.iter_line.startswith("CTRIA3"):
                    eid = int(self.iter_line[8:16])
                    nodes = [int(self.iter_line[24:32]), int(self.iter_line[32:40]), int(self.iter_line[40:48])]
                    self.elements.append(Element(eid, 30500, nodes))
                    self.iter_line = bdf_f.readline()

                elif self.iter_line.startswith("CQUAD4"):
                    eid = int(self.iter_line[8:16])
                    nodes = [int(self.iter_line[24:32]), int(self.iter_line[32:40]),
                             int(self.iter_line[40:48]), int(self.iter_line[48:56])]
                    self.elements.append(Element(eid, 40500, nodes))
                    self.iter_line = bdf_f.readline()

                elif self.iter_line.startswith("CBAR"):
                    eid = int(self.iter_line[8:16])
                    nodes = [int(self.iter_line[24:32]), int(self.iter_line[32:40])]
                    self.iter_line = bdf_f.readline()
                else:
                    # 文件底部
                    if not self.iter_line:
                        break
                    # 不支持的关键字或注释直接读取下一行, 不可以strip(), 否则空行会被认作程序结束
                    self.iter_line = bdf_f.readline()


if __name__ == "__main__":
    filename = "../testcases/Nastran/BoatShowPart.bdf"
    rd = BDFParser(filename)
    rd.ParseFile()
