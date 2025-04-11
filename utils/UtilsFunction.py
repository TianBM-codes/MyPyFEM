#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np


def GetShellGlobal2LocalTransMatrix(nodes: np.array):
    """
    第一个切向量 (沿1-2边方向), 第二个切向量 (沿1-4边方向，假设四边形单元)
    :param nodes: shape: n*3  n>=2 全局笛卡尔坐标, 单元的节点坐标
    :return 转换矩阵
    """
    origin = nodes[0]

    vec1 = nodes[1] - origin
    vec1 = vec1 / np.linalg.norm(vec1)

    vec2_initial = nodes[2] - origin
    vec2_initial = vec2_initial / np.linalg.norm(vec2_initial)

    normal = np.cross(vec1, vec2_initial)
    normal = normal / np.linalg.norm(normal)

    """
    修正第二个切向量使其正交
    """
    vec2 = np.cross(normal, vec1)
    vec2 = vec2 / np.linalg.norm(vec2)

    return np.column_stack((vec1, vec2, normal)), origin


if __name__ == "__main__":
    element_nodes = np.array([
        [20, 5, 0],
        [20, 0, 0],
        [24.33, 0, -2.5],
        [24.33, 5, -2.5]
    ])

    trans_matrix, origin0 = GetShellGlobal2LocalTransMatrix(element_nodes)

    offset_nodes = element_nodes.T - origin0[:, np.newaxis]
    print(offset_nodes.T @ trans_matrix)
