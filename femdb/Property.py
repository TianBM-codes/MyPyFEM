#!/usr/bin/env python3
# -*- coding: utf-8 -*-

class Property(object):
    """
    相当于Nastran的Property, Abaqus的Section
    """

    def __init__(self, e_set_name, m_name, pars):
        self.mat_name = m_name
        self.ele_set_name = e_set_name
        self.pars = pars  # 属性参数, 字典类型

    def __str__(self):
        return "section eles: {}, material name:{}, section pars:{}".format(self.ele_set_name, self.mat_name, self.pars)

    def GetEleSetName(self):
        return self.ele_set_name

    def GetMatName(self):
        return self.mat_name

    def GetPropertyPars(self):
        """ 属性参数,暂时仅支持厚度 """
        return self.pars
