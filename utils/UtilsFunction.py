#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np

import numpy as np


def GetShellGlobal2LocalTransMatrix(nodes: np.ndarray):
    """
    鲁棒的壳单元全局-局部坐标转换矩阵生成器，适用于任意节点排序或空间方向。

    Args:
        nodes: 节点坐标数组，形状为(n,3)，n >= 3

    Returns:
        trans_matrix: 3x3 正交旋转矩阵（列向量为局部基）
        origin: 参考原点坐标（一般取节点 i）
    """
    epsilon = 1e-8

    if nodes.shape[1] != 3 or nodes.shape[0] < 3:
        raise ValueError("节点数组应为 (n,3) 且 n>=3")

    # Step 1: 动态选择两个非共线的向量构建 vec1 和 vec2_initial
    found = False
    for i in range(len(nodes)):
        for j in range(i + 1, len(nodes)):
            v1 = nodes[j] - nodes[i]
            if np.linalg.norm(v1) < epsilon:
                continue
            for k in range(len(nodes)):
                if k == i or k == j:
                    continue
                v2 = nodes[k] - nodes[i]
                if np.linalg.norm(v2) < epsilon:
                    continue
                # 判断非共线
                angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
                if abs(angle) < np.cos(np.deg2rad(5)):  # 超过5度
                    origin = nodes[i]
                    vec1 = v1 / np.linalg.norm(v1)
                    vec2_initial = v2 / np.linalg.norm(v2)
                    found = True
                    break
            if found:
                break
        if found:
            break

    if not found:
        print(nodes)
        raise np.linalg.LinAlgError("无法找到非共线节点构建局部坐标系, 向量之间夹角小于5°")

    # Step 2: 构造法向量 normal，并对 vec2 正交化
    normal = np.cross(vec1, vec2_initial)
    normal_norm = np.linalg.norm(normal)
    if normal_norm < epsilon:
        raise np.linalg.LinAlgError("法向量构造失败，可能节点共面或共线")
    normal /= normal_norm

    vec2 = np.cross(normal, vec1)
    vec2 /= np.linalg.norm(vec2)

    # Step 3: 构造正交旋转矩阵
    trans_matrix = np.column_stack((vec1, vec2, normal))

    # 验证正交性
    if not np.allclose(trans_matrix.T @ trans_matrix, np.eye(3), atol=1e-6):
        raise np.linalg.LinAlgError("局部坐标系不是正交矩阵")

    return trans_matrix, origin


if __name__ == "__main__":
    """
    test case 1
    """
    element_nodes = np.array([
        [20, 5, 0],
        [20, 0, 0],
        [24.33, 0, -2.5],
        [24.33, 5, -2.5]
    ])

    trans_matrix, origin0 = GetShellGlobal2LocalTransMatrix(element_nodes)

    offset_nodes = element_nodes.T - origin0[:, np.newaxis]
    # print(offset_nodes.T @ trans_matrix)
    """
    test case 2
    """
    nodes = np.array([
        [3558.8896116807, -95, 1610.1617154358],
        [3529.457330131, -95, 1610.1569699788],
        [3528.5926865481, -95, 1585],
        [3557.1853691219, -95, 1585]
    ])
    trans_matrix, origin0 = GetShellGlobal2LocalTransMatrix(nodes)

    offset_nodes = nodes.T - origin0[:, np.newaxis]
    # print(offset_nodes.T @ trans_matrix)
    """
    test case 3
    """
    nodes = np.array([
        [-1880.45190016, 8749.84332561, -386.],
        [-1978.26426999, 8757.2743042, -386.],
        [-1928.73000267, 8761.7603147, -386.]
    ])
    trans_matrix, origin0 = GetShellGlobal2LocalTransMatrix(nodes)

    offset_nodes = nodes.T - origin0[:, np.newaxis]
    print(offset_nodes.T @ trans_matrix)
