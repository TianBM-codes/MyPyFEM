# -*- coding: utf-8 -*-

import random
import sys
from elapsedtimer import *
import os
import matplotlib.pyplot as plt

import numpy as np
import time

from scipy import sparse
from scipy.sparse.linalg import inv
from scipy.sparse.linalg import splu
from scipy.io import savemat
from scipy.io import loadmat
from scipy.sparse.linalg import norm
from matplotlib import pyplot as plt

plt.rcParams['font.sans-serif'] = ['Microsoft YaHei']  # 用来正常显示中文标签
plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号

time_format = r"{:>20s} elapsed time: {:<.3f}s"

"""
   1. 要有效地构造矩阵, 请使用dok_matrix或lil_matrix, lil_matrix类支持基本切片和花式索引, 其语法与NumPy Array类似; lil_matrix形式是基于row的
因此能够很高效的转为csr, 但是转为csc效率相对较低.
   2. 强烈建议不要直接使用NumPy函数运算稀疏矩阵如果你想将NumPy函数应用于这些矩阵，首先要检查SciPy是否有自己的给定稀疏矩阵类的实现, 或者首先将稀疏矩
阵转换为NumPy数组(使用类的toarray()方法).
   3. 要执行乘法或转置等操作，首先将矩阵转换为CSC或CSR格式，效率高CSR格式特别适用于快速矩阵矢量产品
   4. CSR，CSC和COO格式之间的所有转换都是线性复杂度.
   # https://stackoverflow.com/questions/40703042/more-efficient-way-to-invert-a-matrix-knowing-it-is-symmetric-and-positive-semi#:~:text=%3E%3E%3E%3E%20M%20%3D%20np.random.rand%20%2810%2C10%29%20%3E%3E%3E%3E%20M%20%3D,inv_M%20%3D%20np.triu%20%28inv_M%29%20%2B%20np.triu%20%28inv_M%2C%20k%3D1%29.T
"""


def GS_arnoldi(A, r, m):
    """
    Gram-Schmidt正交化方法
    """
    print(type(r))
    print(type(A))
    print(type(m))
    print("enter GS_arnoldi function {}".format(time.asctime(time.localtime(time.time()))))
    n = A.shape[0]
    # V = np.zeros((n, m + 1), dtype=float)
    V = sparse.csc_matrix((n, m + 1), dtype=float)
    H = sparse.lil_matrix((m + 1, m), dtype=float)
    # H = np.zeros((m + 1, m), dtype=float)
    # H = sparse.csc_matrix((m + 1, m), dtype=float)
    with ElapsedTimer('r/np.linalg.norm(r)'):
        V[:, 0] = r / norm(r)

    # do arnoldi
    mul2 = []
    sub = []
    total = []

    for j in range(m):
        iter_begin = time.time()
        with ElapsedTimer("A*V[:,j]:{}".format(j)):
            w = A * V[:, j]

        iter_begin2 = time.time()
        with ElapsedTimer("w.T * V[:, i] {}".format(j)):
            H[:j + 1, j] = sparse.csc_matrix.dot(w.T, V[:, 0:j + 1])[0, :j + 1].T
            # H[:j+1, j] = sparse.csc_matrix.dot(V[:, 0:j + 1], w)[:j+1,0]
            # temp = sparse.csc_matrix.dot(w.T, V[:, 0:j + 1])[0, :]
            # for i in range(j + 1):
            #     H[i, j] = sparse.csc_matrix.dot(w.T, V[:, i])[0, 0]

        iter_begin3 = time.time()
        mul2.append("{:.3f}".format(iter_begin3 - iter_begin2))

        with ElapsedTimer("w-H[i, j]*V[:, i] {}".format(j)):
            w = w - V[:, 0:j + 1] * H[0:j + 1, j]
            # for i in range(j + 1):
            #     w = w - H[i, j] * V[:, i]
            # kk = 2

        iter_begin4 = time.time()
        sub.append("{:.3f}".format(iter_begin4 - iter_begin3))
        with ElapsedTimer("norm {}".format(j)):
            H[j + 1, j] = norm(w)

        if H[j + 1, j] == 0:
            break
        else:
            with ElapsedTimer("w/H[j+1,j]"):
                V[:, j + 1] = w / H[j + 1, j]

        iter_end = time.time()
        total.append("{:.3f}".format(iter_end - iter_begin))
        print(time_format.format("iter arnoldi {}".format(j), iter_end - iter_begin))
        # if j == 80:
        #     with open('./krylov.log', 'w') as wf:
        #         for i in range(len(mul2)):
        #             wf.write(mul2[i])
        #             wf.write(' ')
        #             wf.write(sub[i])
        #             wf.write(' ')
        #             wf.write(total[i])
        #             wf.write('\n')
        #     sys.exit(10)

    return V, H


def FOM_krylov(A, x_0, m, b):
    """
    K_Matrix, x0, dimension , Q_Vector
    type(A): csc_matrix
    m: 降阶维数
    FOM方法进行求解, m为降阶维数
    """
    # r = b - A * x_0.T
    tt = A * x_0
    r = (b.T - tt).tocsc()
    print(type(r))
    deter = norm(r)
    V, H = GS_arnoldi(A, r, m)
    sparse.save_npz('big_V.npz', V)
    sparse.save_npz('big_H.npz', H.tocsc())
    print("finish GS_arnoldi {}".format(time.asctime(time.localtime(time.time()))))
    # e1 = np.eye(m, 1)

    # y = deter * np.matmul(np.linalg.inv(H[0:m, :]), e1)
    # y = deter * inv(H[0:m, :]).diagonal()
    y = deter * inv(H[0:m, :])[:, 0].toarray().flatten()
    # y = deter * np.matmul(np.linalg.inv(H[0:m, :]), e1)
    print("finish inverse matrix {}".format(time.asctime(time.localtime(time.time()))))
    error = H[m, m - 1] * np.abs(y[m - 1])
    x = x_0.toarray().flatten() + V[:, :len(y)] * y
    return x, error


def steady_from_np(K_path, Q_path, dimension):
    """
    从numpy的压缩文件中读取数据
    """
    r_begin = time.time()
    with ElapsedTimer("load npz"):
        K_Matrix = sparse.load_npz(K_path)
    with ElapsedTimer("to csc matrix"):
        K_Matrix = K_Matrix.tocsc()
    Q_Vector = np.load(Q_path)
    Q_Vector = sparse.csc_matrix(Q_Vector)
    K_row_num = K_Matrix.shape[0]
    # print("current K_Matrix's shape is:{}, len(Q_Vector) is:{}".format(K_Matrix.shape, len(Q_Vector)))
    # print("current Krylov dimension is:{}".format(dimension))
    # print("")
    print(type(K_Matrix))
    print(K_row_num)

    """
    FOM_krylov
    """
    c_begin = time.time()
    # x0 = np.zeros(K_row_num, dtype=float)
    x0 = sparse.csc_matrix((K_row_num, 1), dtype=float)
    x, e = FOM_krylov(K_Matrix, x0, dimension, Q_Vector)

    """
    print Summary
    """
    k_end = time.time()
    print(time_format.format("Read Sparse Matrix", c_begin - r_begin))
    print(time_format.format("FOM_krylov", k_end - c_begin))
    print(f"x[0]-x[9]:{x[:10]}\nerror is:{e}")


def plot_summary():
    count = 1000
    print("step 1000: elapsed time: {} hours".format(41.9 / 6400 * count ** 2 / 60 / 60))
    sum_hour = 0
    for ii in range(count):
        sum_hour += 41.9 / 6400 * ii ** 2 / 3600
    print("sum time: {:.3f} hours".format(sum_hour))

    with open("./krylov.log") as kf:
        h1 = []
        h2 = []
        t = []
        lines = kf.readlines()
        for iter_line in lines:
            iter_line = iter_line.split()
            h1.append(float(iter_line[0]))
            h2.append(float(iter_line[1]))
            t.append(float(iter_line[2]))
        x = np.arange(len(h1))
        # y = [9 / 640 * xi ** 2 for xi in x]
        y = [41.9 / 6400 * xi ** 2 for xi in x]
        plt.plot(x, h1, label='子空间迭代1')
        plt.plot(x, h2, label='子空间迭代2')
        plt.plot(x, t, label='单步耗时')
        plt.plot(x, y, label='y=9/640 x**2')
        plt.xlabel('迭代步')
        plt.ylabel('时间(s)')
        plt.legend()
        plt.show()


def matlab_time():
    with open('./matlab.txt', 'r') as mf:
        lines = mf.readlines()
        data = []
        for line in lines:
            data.append(float(line.strip()))

    x = np.arange(1, len(data) + 1)
    plt.plot(x, data)
    plt.show()


if __name__ == "__main__":
    krylov_dimension = 1500
    steady_from_np("./big_steady/steady_big_K_matrix.npz",
                   "./big_steady/steady_big_Q_Vector.npy", krylov_dimension)
    # steady_from_np("./less_dof_steady_example/steady_less_K_matrix.npz",
    #                "./less_dof_steady_example/steady_less_Q_Vector.npy", krylov_dimension)
    # steady_from_np("./K_matrix684.npz",
    #                "./Q_Vector684.npy", krylov_dimension)
    # plot_summary()
    # matlab_time()
