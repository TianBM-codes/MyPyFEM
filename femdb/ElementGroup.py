#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from femdb.ElementFactory import *
from typing import List
from element.ElementBase import ElementBaseClass


class ElementGroup:
    """
    单元组, 每个单元组中存储的都是相同类型的单元, 类提供同一类型单元的一些特性
    """

    def __init__(self, e_type):
        self.e_type = e_type
        # 其中存储的是同一单元类型的单元, 存储了单元的全部信息
        self.eles: List[ElementBaseClass] = []
        # 单元组中的所有真实ID, 类型是set, 可以快速查找对应ID是否属于该组
        self.ele_id_sets = set()

    def AppendElement(self, r_ele):
        """
        添加单元至单元组
        :param r_ele: 实际的单元
        """
        self.eles.append(r_ele)

    def SetEleIdSet(self, id_set):
        """
        设置单元组内包含的id, 因为需要查找某些单元是否属于该组, 所以
        设置为set类型, set是一种树结构, 比list列表结构快很多
        :param id_set: 节点id的set
        """
        self.ele_id_sets = id_set

    def Elements(self):
        return self.eles

    def IsElementInGroup(self, eid: int) -> bool:
        """
        通过单元id判断单元是否属于本Group, 从而判断单元的类型
        :param eid: 单元的真实ID
        """
        if eid in self.ele_id_sets:
            return True
        else:
            return False

    def GetElementsCurrentCount(self):
        """
        用于返回当前单元组已经添加了多少单元, 用于获取Index, 单元真实Id对应在ElementGroup中的eles中的idx
        """
        return len(self.eles)

    def GetZeroElementStiffMatrix(self):
        """
        返回一个空的Matrix, 在集成总刚的时候避免重复申请单刚
        """
        dof_count = ElementFactory.GetElementNodeDofCount(self.e_type)
        return np.zeros(dof_count, dof_count)

    def CalElementsEquationNum(self):
        """
        计算该单元组中所有单元包含节点对应的方程号, 为集成总刚做准备
        """
        for e in self.eles:
            e.CalEquationNumber()
