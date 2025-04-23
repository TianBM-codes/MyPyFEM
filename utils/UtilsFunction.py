#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np

import numpy as np


def GetShellGlobal2LocalTransMatrix2(nodes: np.ndarray):
    """
    鲁棒的壳体单元全局-局部坐标系转换矩阵生成器
    特征：
    1. 动态选择非共线节点构建正交基
    2. 处理平面退化等极端情况
    3. 增强数值稳定性

    Args:
        nodes: 节点坐标数组，形状为(n,3)，n>=2

    Returns:
        trans_matrix: 3x3转换矩阵 [vec1, vec2, normal]
        origin: 局部坐标系原点

    Raises:
        ValueError: 输入节点不符合维度要求
        np.linalg.LinAlgError:  无法生成正交基
    """
    # 输入验证
    if nodes.shape[1] != 3 or len(nodes) < 2:
        raise ValueError("节点数组形状应为(n,3)且n>=2")

    origin = nodes[0].astype(float)
    epsilon = 1e-8  # 浮点误差阈值

    # 第一基准向量（固定取节点0到节点1方向）
    vec1 = nodes[1] - origin
    vec1_norm = np.linalg.norm(vec1)
    if vec1_norm < epsilon:
        raise np.linalg.LinAlgError(" 前两个节点重合")
    vec1 = vec1 / vec1_norm

    # 动态选择第二基准向量（遍历后续节点）
    vec2_initial = None
    for i in range(2, len(nodes)):
        candidate = nodes[i] - origin
        candidate_norm = np.linalg.norm(candidate)
        if candidate_norm < epsilon:
            continue  # 跳过重合节点

        candidate_normalized = candidate / candidate_norm
        cos_theta = np.dot(vec1, candidate_normalized)

        # 选择与vec1夹角>10度的向量（cosθ<0.9848）
        if abs(cos_theta) < np.cos(np.deg2rad(10)):
            vec2_initial = candidate_normalized
            break

            # 极端情况处理：所有节点共线
    if vec2_initial is None:
        # 尝试构造垂直于vec1的向量
        if np.linalg.norm(vec1[:2]) > epsilon:  # 非纯Z向量的情况
            vec2_initial = np.array([-vec1[1], vec1[0], 0.0])
        else:  # 纯Z向量的特殊情况
            vec2_initial = np.array([1.0, 0.0, 0.0])
        vec2_initial /= np.linalg.norm(vec2_initial)

        # 计算法向量并正交化
    normal = np.cross(vec1, vec2_initial)
    normal_norm = np.linalg.norm(normal)

    if normal_norm < epsilon:  # 检测平面退化
        # 采用全局Z轴作为备选法向量
        global_z = np.array([0.0, 0.0, 1.0])
        cos_theta = np.dot(vec1, global_z)

        if abs(cos_theta) > 0.98:  # 避免与global_z接近共线
            global_ref = np.array([1.0, 0.0, 0.0])
        else:
            global_ref = global_z

        normal = np.cross(vec1, global_ref)
        normal_norm = np.linalg.norm(normal)
        if normal_norm < epsilon:
            raise np.linalg.LinAlgError(" 无法生成法向量")

    normal = normal / normal_norm

    # 正交化第二切向量
    vec2 = np.cross(normal, vec1)
    vec2 /= np.linalg.norm(vec2)

    # 构建转换矩阵
    trans_matrix = np.column_stack((vec1, vec2, normal))

    # 验证正交性
    if not np.allclose(trans_matrix.T @ trans_matrix, np.eye(3), atol=1e-6):
        raise np.linalg.LinAlgError(" 生成的非正交矩阵")

    return trans_matrix, origin


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
                if abs(angle) < np.cos(np.deg2rad(10)):  # 超过10度
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
        raise np.linalg.LinAlgError("无法找到非共线节点构建局部坐标系")

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
    print(offset_nodes.T @ trans_matrix)
    print("-------")
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
    print(offset_nodes.T @ trans_matrix)
