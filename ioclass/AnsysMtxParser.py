#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from HarwellBoeingParser import HarwellBoeingMatrix
import numpy as np


def ReadANSYSStiffness(f_name, print_log=False):
    """
    ANSYS命令流如下
    /AUX2
    FILE,'Model_1',full
    HBMAT, 'Stiffness_mat', dat, , ASCII, STIFF, YES, YES
    HBMAT, 'Mass_mat', dat, ,ASCII, MASS, YES, YES
    FINISH

    Note:
    1. 保存的矩阵是Harwell-Boeing格式的, scipy只支持非对称的矩阵, 以下程序是在PyOrder中截取的
    2. 使用pyANSYS中的ansysMath, 对python版本由限制，并且环境不好配置
    3. 如果用python这个版本出现问题, 那么可以参考一个C库: D:\WorkSpace\qtcreator\HarwellBoeing,
       参考网站是 https://people.math.sc.edu/Burkardt/cpp_src/hb_io/hb_io.html
    4. Harwell-Boeing格式说明: https://zhuanlan.zhihu.com/p/587855571
    5. 还有一个较为复杂的c库: http://bebop.cs.berkeley.edu/smc/installation.html

    Reference:
    1. https://math.docs.pyansys.com/version/stable/examples/eigen_solve.html
    2. http://feacat.com/article1/ANSYS%20%E4%B8%AD%E6%95%B4%E4%BD%93%E5%88%9A%E5%BA%A6%E7%9F%A9%E9%98%B5%E7%9A%84%E8%BE%93%E5%87%BA%EF%BC%88%E4%B8%80%EF%BC%89%EF%BC%9A%E5%AD%90%E7%BB%93%E6%9E%84%E5%88%86%E6%9E%90.html
    3. https://www.zhihu.com/question/374768464
    4. http://feacat.com/article1/ANSYS%20%E4%B8%AD%E6%95%B4%E4%BD%93%E5%88%9A%E5%BA%A6%E7%9F%A9%E9%98%B5%E7%9A%84%E8%BE%93%E5%87%BA%EF%BC%88%E4%BA%8C%EF%BC%89%EF%BC%9AHBMAT%E5%91%BD%E4%BB%A4.html
    5. https://blog.csdn.net/m0_61038364/article/details/127498110
    6. https://www.leanwind.com/archives/4148.html
    7. https://www.jishulink.com/post/1796144
    8. https://zhuanlan.zhihu.com/p/198212529
    9. https://blog.csdn.net/m0_61038364/article/details/127498110
    """
    K_hb = HarwellBoeingMatrix(f_name, patternOnly=False, readRhs=True)
    (val, row, col) = K_hb.find()
    K = np.zeros((K_hb.nrow, K_hb.ncol), dtype=float)
    for ii in range(len(val)):
        K[row[ii], col[ii]] = val[ii]

    # 现在读取的矩阵是下三角矩阵, 需要使其对称化
    K = K + K.T
    for ii in range(K_hb.nrow):
        K[ii, ii] = K[ii][ii] * 0.5

    # 右端项q
    q = K_hb.rhs
    if print_log:
        print("ANSYS Stiffness.shape is:({},{})".format(K_hb.nrow, K_hb.ncol))
        print("ANSYS Displacement is :\n{}".format(np.matmul(np.linalg.inv(K), q).flatten()))


if __name__ == "__main__":
    f_path = "../../testcases/ANSYS/shell/shell/stiffness_mat.dat"
    ReadANSYSStiffness(f_path, True)
