#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np


def GetGlobal2LocalTransMatrix(global_coord: np.array):
    """
    注意节点编号规则, 需是逆时针编号
    :param global_coord: shape: n*3  n>=2 全局笛卡尔坐标
    :return 局部坐标, x_local是节点1至节点2的单位向量, y_local是面外法向
    """
    x_local = global_coord[1,] - global_coord[0,]
    y_temp = global_coord[2,] - global_coord[1,]
    z_local = np.cross(x_local, y_temp)
    y_local = np.cross(z_local, x_local)

    return np.array([x_local[0:] / np.linalg.norm(x_local),
                     y_local[0:] / np.linalg.norm(y_local),
                     z_local[0:] / np.linalg.norm(z_local)], dtype=float)


if __name__ == "__main__":
    # t_coord = np.asarray([[3, 0, 0], [0, 3, 0], [0, 0, 3]],dtype=float)
    # t = np.asarray([[1,2,3],[4,5,6],[7,8,9]])
    # print(t.T*t_coord)
    # print(np.matmul(t,t_coord)+t_coord)
    t = np.matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    t_coord = np.array([[3, 0, 0], [0, 3, 0], [0, 0, 3]], dtype=float)
    print(GetGlobal2LocalTransMatrix(t_coord))
    print(GetGlobal2LocalTransMatrix(t_coord).shape)
    print(t.shape)
    print(t_coord.shape)
    # print(t*t_coord)
