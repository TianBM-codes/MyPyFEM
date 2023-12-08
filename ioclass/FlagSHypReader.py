#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from femdb.FEMDataBase import *

class FlagSHypReader(object):
    """
    读取FlagSHy的录入文件
    """

    def __init__(self, input_path):
        self.femdb = FEMDataBase()
        self.dat_path = input_path
        self.iter_line = None
        self.et_hash = {}
        self.ele_group_hash = {}
        self.ele_count = 0

    def ParseFileAndInitFEMDB(self):
        """
        文件格式：
        
        :return:
        """