#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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
