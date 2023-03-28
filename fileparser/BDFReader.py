#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pyNastran.bdf.bdf import BDF
from femdb.FEMDataBase import *


class BDFReader(object):
    """
    功能为解析BDF文件, 初始化数据库单元节点及相关信息, 只负责解析和传递, 其他功能均在其他类中,
    同样的, 其他类也不会包含任何解析输入文件的功能或者函数
    """

    def __init__(self, input_path):
        self.fem_data = FEMDataBase()
        self.cdb_path = input_path
        self.archive = BDF()

    def ParseFile(self):
        """
        解析Nastran软件的BDF文件, 主要功能为将pyNastran库读取的结果转化为自身程序适用的, 初始化数据库
        Reference:
        1. https://pynastran-git.readthedocs.io/en/1.3/quick_start/bdf_overview.html#example-1-read-write
        """
        mlogger.debug("Parsing CDB file: {}".format(self.cdb_path))
        bdf_filename = "./bdf"
        self.archive.read_bdf(bdf_filename)


if __name__ == "__main__":
    filename = "./.bdf"
    rd = BDFReader(filename)
