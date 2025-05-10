#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy


class LoadCase(object):
    """ Class LoadCase is used to store load data """

    def __init__(self):
        self.c_loads = []
        self.boundaries = []
        self.history_loads = []

    def __str__(self):
        self.case_ = "\n  Here is LoadCase:\n"
        desc = self.case_
        for bd in self.boundaries:
            desc += "  {}\n".format(bd)
        for cld in self.c_loads:
            desc += "  set name:{}, direction:{}, value:{}".format(cld.set_name, cld.direction, cld.value)
        return desc

    def AddBoundary(self, boundary: numpy.ndarray):
        """
        无论是ANSYS、ABAQUS还是NASTRAN格式的约束都存在这里
        :param boundary:
        :return: None
        """
        self.boundaries.append(boundary)

    def GetBoundaries(self):
        return self.boundaries

    def GetConcentratedLoads(self):
        return self.c_loads

    def AddHistoryLoad(self, node, directory, scale, amps):
        """
        添加时间载荷
        @param node:
        @param directory:
        @param scale:
        @param amps:
        @return:
        """
        self.history_loads.append((node, directory, scale, amps))

    def AddConcentratedLoad(self, node, directory, amp):
        """
        添加集中力载荷
        """
        self.c_loads.append((node, directory, amp))
