#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pylab as plt
import numpy as np
# import ansys.math.core.math as pymath
from ansys.mapdl.core import launch_mapdl


# def GetStiffnessMatrixFromFullFile(file_path):
#     """
#     通过full文件获取ANSYS的刚度矩阵
#     Reference:
#     1. https://math.docs.pyansys.com/version/stable/examples/eigen_solve.html
#     2. http://feacat.com/article1/ANSYS%20%E4%B8%AD%E6%95%B4%E4%BD%93%E5%88%9A%E5%BA%A6%E7%9F%A9%E9%98%B5%E7%9A%84%E8%BE%93%E5%87%BA%EF%BC%88%E4%B8%80%EF%BC%89%EF%BC%9A%E5%AD%90%E7%BB%93%E6%9E%84%E5%88%86%E6%9E%90.html
#     3. https://www.zhihu.com/question/374768464
#     4. http://feacat.com/article1/ANSYS%20%E4%B8%AD%E6%95%B4%E4%BD%93%E5%88%9A%E5%BA%A6%E7%9F%A9%E9%98%B5%E7%9A%84%E8%BE%93%E5%87%BA%EF%BC%88%E4%BA%8C%EF%BC%89%EF%BC%9AHBMAT%E5%91%BD%E4%BB%A4.html
#     5. https://blog.csdn.net/m0_61038364/article/details/127498110
#     6. https://www.leanwind.com/archives/4148.html
#     7. https://www.jishulink.com/post/1796144
#     """
#
#     # Start PyAnsys Math as a service.
#     mm = pymath.AnsMath()
#
#     # fullfile = mm._mapdl.jobname + ".full"
#     fullfile = "./fourshell0.full"
#     k = mm.stiff(fname=fullfile)
#     m = mm.mass(fname=fullfile)
#     print(m.shape)
#     print(k.shape)


def GetStiffnessMatrixFromFullFile2():
    mapdl = launch_mapdl()
    print(mapdl)
    # mm = mapdl.math
    # print(mapdl.input("./fourshell0.full"))
    # mapdl.eplot()


GetStiffnessMatrixFromFullFile2()

# -*- coding: utf-8 -*-
import random

import numpy as np
import time
from scipy import sparse
from matplotlib import pyplot as plt

plt.rcParams['font.sans-serif'] = ['Microsoft YaHei']  # 用来正常显示中文标签
plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号



def jump_out_loop(line: str):
    """
    判断是否还是保存在同一个矩阵中, 是否跳出while循环,
    """
    if line.__contains__("ROW"):
        return True
    if line.__contains__("LOAD"):
        return True
    if line.startswith("\n"):
        return True
    if line is None:
        return True
    if line == "":
        return True

    return False


def Transient():
    """
    瞬态方程求解
    """
    """
    数据准备, 初始化
    """
    time_format = r"{:>20s} elapsed time: {:<.3f}s"
    p_begin = time.time()
    nn = 180  # 矩阵维数
    K_Matrix = []
    C_Matrix = []
    Q_Vector = []
    with open("./transient_example/test9.txt", 'r') as rf:
        iter_line = rf.readline()
        while iter_line:
            if iter_line.__contains__("ROW"):
                # 读取行数据
                current_matrix = iter_line.split()[-1]
                iter_line = rf.readline()
                iter_row = []
                while not jump_out_loop(iter_line):
                    iter_row.extend([float(iv) for iv in iter_line.split()])
                    iter_line = rf.readline()

                # 保存
                if current_matrix == "1":
                    K_Matrix.append(iter_row)
                elif current_matrix == "2":
                    C_Matrix.append(iter_row)

            elif iter_line.__contains__("LOAD"):
                iter_line = rf.readline()
                iter_row = []
                while not jump_out_loop(iter_line):
                    iter_row.extend(float(iv) for iv in iter_line.split())
                    iter_line = rf.readline()
                Q_Vector.extend(iter_row)

            else:
                iter_line = rf.readline()

    K_Matrix = np.asarray(K_Matrix, dtype=np.float64)
    C_Matrix = np.asarray(C_Matrix, dtype=np.float64)
    Q_Vector = np.asarray(Q_Vector, dtype=np.float64)
    c_begin = time.time()
    print(time_format.format("Parse File", c_begin - p_begin))

    """
    降阶计算
    """
    T = 100
    N = 10000
    dt = T / N
    t = [dt * i for i in range(N)]
    # A = -np.matmul(np.linalg.inv(C_Matrix), K_Matrix)
    # b = np.matmul(np.linalg.inv(C_Matrix), Q_Vector)
    r0 = np.mat([random.randint(21, 22) for _ in range(nn)]).T
    AA, BB, V, xx = UP_GS_arnoldi(C_Matrix, K_Matrix, Q_Vector, 20, r0)

    # 降阶维数, 初始解
    y2 = [xx]
    for n in range(len(t)):
        iter_y = y2[-1] + (np.matmul(AA, y2[-1]) + BB) * dt
        y2.append(iter_y)
    y2re = np.matmul(V, np.asarray(y2))
    j_end = time.time()
    print(time_format.format("Reduced Calculate", j_end - c_begin))

    """
    未降阶计算
    """
    r1 = np.mat([random.randint(20, 22) for _ in range(nn)]).T
    y = np.mat(np.zeros((len(r1), len(t) + 1)))
    y[:, 0] = r0

    # 欧拉中心差分法
    inv_C = np.linalg.inv(C_Matrix)
    for n in range(len(t)):
        iter_y = y[:, n] + np.matmul(inv_C, (-np.matmul(K_Matrix, y[:, n]) + np.mat(Q_Vector).T)) * dt
        y[:, n + 1] = iter_y
    y = np.mat(y)
    w_end = time.time()
    print(time_format.format("Euler central diff", w_end - j_end))

    """
    Ansys参考手册计算方法
    """
    y1 = np.mat(np.zeros((len(r1), len(t) + 1)))
    y1[:, 0] = r0
    for i in range(len(t)):
        iter_ly = np.linalg.inv(C_Matrix + dt * K_Matrix)
        iter_ry = np.matmul(C_Matrix, y1[:, i]) + dt * np.mat(Q_Vector).T
        y1[:, i + 1] = np.matmul(iter_ly, iter_ry)
    a_end = time.time()
    print(time_format.format("Ansys Method", a_end - w_end))

    """
    结果绘制
    """
    t.append(T)  # matlab与python有差别的地方
    plt.plot(t, y1[22, :].T, label="Ansys参考手册")
    plt.plot(t, y[22, :].T, label="未降阶计算")
    plt.plot(t, y2re[:, 22], label="降阶计算")
    plt.xlabel(r"时间/s")
    plt.ylabel(r"温度/℃")
    plt.title("23号节点温度变化图")
    plt.legend()
    plt.show()


def steady():
    """
    稳态分析
    """
    K_Matrix = []
    Q_Vector = []
    with open("./less_dof_steady_example/testM.txt", 'r') as rf:
    # with open("./steady_example/test1205.txt", 'r') as rf:
        iter_line = rf.readline()
        while iter_line:
            if iter_line.__contains__("ROW"):
                # 读取行数据
                current_matrix = iter_line.split()[-1]
                iter_line = rf.readline()
                iter_row = []
                while not jump_out_loop(iter_line):
                    iter_row.extend([float(iv) for iv in iter_line.split()])
                    iter_line = rf.readline()

                # 保存
                if current_matrix == "1":
                    K_Matrix.append(iter_row)

            elif iter_line.__contains__("LOAD"):
                iter_line = rf.readline()
                iter_row = []
                while not jump_out_loop(iter_line):
                    iter_row.extend(float(iv) for iv in iter_line.split())
                    iter_line = rf.readline()
                Q_Vector.extend(iter_row)

            else:
                iter_line = rf.readline()

    K_Matrix = np.asarray(K_Matrix)
    Q_Vector = np.asarray(Q_Vector)

    """
    FOM_krylov
    """
    x0 = np.zeros(18, dtype=float)
    x = FOM_krylov(K_Matrix, x0, 6, Q_Vector)
    print(x[:10])


if __name__ == "__main__":
    # Transient()
    steady()
