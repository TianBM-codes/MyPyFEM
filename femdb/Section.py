#!/usr/bin/env python3
# -*- coding: utf-8 -*-


class Section(object):
    """
    保存截面信息, 包括梁的还有壳的
    """

    def __init__(self, sec_name, sec_type, sec_data):
        self.name = sec_name
        self.sec_type = sec_type
        self.sec_data = sec_data
