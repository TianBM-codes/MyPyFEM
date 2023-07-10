#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys

import numpy as np

from femdb.GlobalEnum import *


class GaussIntegrationPoint:
    @staticmethod
    def GetSamplePointAndWeight(n):
        """
        返回指定阶数的采样点以及权重, Bathe P436, 对于起始为(a,b)的,
        采样点为: (a+b)/2 + (b-a)/2 * r_i
        权重为: (b-a)/2 * alpha_i
        :return: (sample p1, sample p2,...)  (weight1, weight2, ...)
        """
        if n == 1:
            sample_point = (0.0,)
            weight = (0.0,)
            return sample_point, weight
        elif n == 2:
            sample_point = (-0.577350269189626, 0.577350269189626)
            weight = (1.0, 1.0)
            return sample_point, weight
        elif n == 3:
            sample_point = (-0.774596669241483, 0.0, 0.774596669241483)
            weight = (0.555555555555555, 0.888888888888888, 0.555555555555555)
            return sample_point, weight
        else:
            mlogger.fatal("Wrong Sample Point Count")
            sys.exit(1)

    @staticmethod
    def GetTrianglePointAndWeight(n):
        """
        三角形和四面体的积分点以及权重不同, 涉及到斜线积分, n代表积分点个数
        Reference:
        1. 《有限元分析的概念与应用》 Cook R.D pdf P230
        """
        if n == 1:
            sample_point = [(1 / 3, 1 / 3)]
            weight = (1.0,)
            return sample_point, weight
        elif n == 3:
            sample_point = [(0.5, 0), (0, 0.5), (0.5, 0.5)]
            weight = (1 / 3, 1 / 3, 1 / 3)
            return sample_point, weight
        elif n == 4:
            sample_point = [(1 / 3, 1 / 3), (3 / 5, 1 / 5), (1 / 5, 1 / 5), (1 / 5, 3 / 5)]
            weight = (-27 / 48, 25 / 48, 25 / 48, 25 / 48)
            return sample_point, weight
        else:
            mlogger.fatal("Wrong Sample Point Count")
            sys.exit(1)

    @staticmethod
    def GetTetraPointAndWeight(n):
        """
        四面体的积分点以及权重与一般的不同, 涉及到斜面积分, n代表积分点个数
        Reference:
        1. 《有限元分析的概念与应用》 Cook R.D pdf P231
        """
        if n == 1:
            sample_point = [(0.25, 0.25, 0.25)]
            weight = 1.0
            return sample_point, weight
        elif n == 4:
            a = (5 + 3 * np.sqrt(5)) / 20
            b = (5 - np.sqrt(5)) / 20
            sample_point = [(a, b, b), (b, b, b), (b, b, a), (b, a, b)]
            weight = (0.25, 0.25, 0.25, 0.25)
            return sample_point, weight
        elif n == 5:
            sample_point = [(0.25, 0.25, 0.25), (0.5, 1 / 6, 1 / 6), (1 / 6, 1 / 6, 1 / 6), (1 / 6, 1 / 6, 0.5), (1 / 6, 0.5, 1 / 6)]
            weight = (-0.8, 9 / 20, 9 / 20, 9 / 20, 9 / 20)
            return sample_point, weight
        else:
            mlogger.fatal("Wrong Sample Point Count")
            sys.exit(1)
