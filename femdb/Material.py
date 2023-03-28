#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from femdb.GlobalEnum import *

class ISOMaterial(object):
    """ 各向同性材料 """

    def __init__(self, name, value_dict=None):
        self.name = name
        self.value_dict = value_dict

    def __str__(self):
        return "ISOMaterial: {},{}".format(self.name, self.value_dict)

    def GetName(self):
        return self.name

    def GetValueDict(self):
        return self.value_dict

    def SetValueDict(self, value_dict):
        self.value_dict = value_dict