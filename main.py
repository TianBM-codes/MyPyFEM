#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pathlib
import time
import sys
from utils.GlobalEnum import *
from ioclass.INPParser import InpParser
from ioclass.CDBParser import CDBParser
from ioclass.BDFParser import BDFParser
from ioclass.FlagSHyPParser import FlagSHyPParser
from ioclass.ResultsWriter import ResultsWriter
from femdb.Domain import Domain
from femdb.NLDomain import NLDomain

sys.path.insert(0, "./femdb")
sys.path.insert(0, "../NumericalCases")


class MyPyFEM:
    """
    TODO: 可以将单元的面编号，对面做一个可识别的ID，用于区分，ID = str(sorted(nodeIds)), 参考MySTAP C++
    """

    def __init__(self, file_path, open_paraview=False, check_model=False, plot_stiff=False):
        if isinstance(file_path, str):
            file_path = pathlib.Path(file_path)

        if not os.path.isfile(file_path):
            logging.fatal("Input File is doesn't exist! {}".format(file_path))
            sys.exit(1)

        # 程序开始时间以及解析文件完成时间
        self.program_begin = None
        self.parsed_time = None

        # 是否绘制总刚度阵
        if plot_stiff:
            GlobalInfor[GlobalVariant.PlotGlobalStiffness] = True

        self._fem_data = None
        self.input_file_path = file_path
        self.output_files = [file_path.with_suffix(".vtu"), file_path.with_suffix(".unv")]

        self.FEMAnalyseFlow(check_model)

        # 结果查看, Paraview显示, 注意要将paraview的路径加入至环境变量
        if open_paraview and not check_model:
            os.popen("paraview " + str(self.output_files[0].absolute()))

    def InitReader(self):
        """
        根据文件类型初始化不同的文件解析器, 然后读入文件初始化数据库
        """
        suffix = self.input_file_path.suffix
        if suffix == ".inp":
            GlobalInfor[GlobalVariant.InputFileSuffix] = InputFileType.INP
            return InpParser(self.input_file_path)
        elif suffix == ".cdb":
            GlobalInfor[GlobalVariant.InputFileSuffix] = InputFileType.CDB
            return CDBParser(self.input_file_path)
        elif suffix == ".bdf":
            GlobalInfor[GlobalVariant.InputFileSuffix] = InputFileType.BDF
            return BDFParser(self.input_file_path)
        elif suffix == ".dat":
            GlobalInfor[GlobalVariant.InputFileSuffix] = InputFileType.FlagSHyP
            return FlagSHyPParser(self.input_file_path)
        else:
            mlogger.fatal("UnSupport File Suffix:{}".format(suffix))
            sys.exit(1)

    def FEMAnalyseFlow(self, check_model):
        """
        TODO: 标准流程, 完成注释, 重写mlogger的debug信息, 将有限元模型的信息输出, 比如单元类型及相应个数, 自由度个数
        求解文件, 步骤如下所示, 该函数中不应包含对不同文件类型的分类, 即判断文件类型的bdf cdb等应在其他函数中完成
        """
        self.program_begin = time.time()
        mlogger.debug("{} Analysis Calculate Begin {}".format("#" * 6, "#" * 6))
        reader = self.InitReader()
        reader.ParseFileAndInitFEMDB()
        self.parsed_time = time.time()

        if check_model:
            reader.CheckModel()
            c_end = time.time()
            mlogger.debug("Elapsed time: {:.3f} seconds\n".format(c_end - self.program_begin))

        if GlobalInfor[GlobalVariant.AnaType] == AnalyseType.LinearStatic:
            """
            求解线弹性问题, 输出节点位移以及应力
            """
            domain = Domain()
            domain.AssignElementCharacter()
            time_2 = time.time()

            domain.CalBoundaryEffect()
            domain.CalculateEquationNumber()
            time_3 = time.time()

            domain.CalAllElementStiffness()
            time_4 = time.time()

            domain.AssembleStiffnessMatrix()
            time_5 = time.time()

            domain.SolveDisplacement()
            time_6 = time.time()

            domain.CalculateNodeStress()
            time_7 = time.time()

            writer = ResultsWriter()
            # writer.WriteVTPFile(output_paths[0])
            writer.WriteUNVFile(self.output_files[1])
            p_end = time.time()

            # Print FEMDB Information
            summary = domain.femdb.GetModelSummary()
            mlogger.debug(" " + "-" * 40)
            summary_format = r"{:>25s} --> {:<}"
            mlogger.debug(" Model Summary:")
            for key, value in summary.items():
                mlogger.debug(summary_format.format(key, value))
            mlogger.debug(" " + "-" * 40)

            # Define Output Format And Print Each Step Time Elapsed
            time_format = r"{:>25s} --> {:<.3f} seconds"
            last_line_format = "{:>25s} --> {:<.3f} seconds"

            mlogger.debug(" Elapsed Time Summary:")
            mlogger.debug(time_format.format("Parse File", self.parsed_time - self.program_begin))
            mlogger.debug(time_format.format("Calculate D", time_2 - self.parsed_time))
            mlogger.debug(time_format.format("Cal Equation Num", time_3 - time_2))
            mlogger.debug(time_format.format("Cal All Stiffness", time_4 - time_3))
            mlogger.debug(time_format.format("Assemble Stiffness", time_5 - time_4))
            mlogger.debug(time_format.format("Solve Displacement", time_6 - time_5))
            mlogger.debug(time_format.format("Calculate Stress", time_7 - time_6))
            mlogger.debug(time_format.format("Write Output File", p_end - time_6))
            mlogger.debug(last_line_format.format("Total Elapsed Time", p_end - self.program_begin))
            mlogger.debug(" " + "-" * 40)
            mlogger.debug(" Finish Analysis\n")

        elif GlobalInfor[GlobalVariant.AnaType] == AnalyseType.NLStatic:
            nl_domain = NLDomain()
            nl_domain.Initialisation()
            nl_domain.ChooseIncrementalAlgorithm()

        else:
            mlogger.fatal(f"UnSupport Analyse Type:{GlobalInfor[GlobalVariant.AnaType]}")
            sys.exit(1)


def check_symmetry(matrix, percentage_tol=0.1):
    # 将百分比容差转换为小数形式
    tol = percentage_tol / 100.0

    # 标记是否对称
    is_symmetric = True

    # 创建字典存储 (i, j) -> value
    matrix_dict = {(i, j): v for i, j, v in zip(matrix.row, matrix.col, matrix.data)}

    # 遍历所有非零元素
    for i, j, v in zip(matrix.row, matrix.col, matrix.data):
        # 检查对应的 (j, i) 元素
        if (j, i) in matrix_dict:
            corresponding_value = matrix_dict[(j, i)]
            # 计算相对误差
            relative_error = abs(v - corresponding_value) / max(abs(v), abs(corresponding_value), 1e-10)
            absolute_error = abs(v - corresponding_value)

            # 检查相对误差和绝对误差
            if relative_error > tol and absolute_error > 1e-4:
                print(f"Matrix element ({i}, {j}) = {v} is not approximately equal to ({j}, {i}) = {corresponding_value} "
                      f"with a relative error of {relative_error:.5f} and an absolute error of {absolute_error:.5e}")
                is_symmetric = False
        else:
            print(f"Matrix element ({i}, {j}) = {v} has no corresponding symmetric element at ({j}, {i})")
            is_symmetric = False

    if is_symmetric:
        print("The matrix is symmetric within the given tolerance.")
    else:
        print("The matrix is not symmetric within the given tolerance.")

    return is_symmetric

if __name__ == "__main__":
    # if len(sys.argv) == 1:
    #     mlogger.fatal(" No Input File Assign")
    #     sys.exit(1)
    # input_file = sys.argv[1]
    # MyPyFEM(pathlib.Path(input_file))
    import numpy as np
    import matplotlib.pyplot as plt
    import scipy.sparse

    # 加载CSV文件
    # 假设CSV文件格式是三列：[row_index, col_index, value]
    data = np.loadtxt('GLOBAL_K_2.csv', delimiter=',')

    # 获取行索引、列索引和值
    rows = data[:, 0].astype(int) - 1  # 将索引转换为0-based
    cols = data[:, 1].astype(int) - 1  # 将索引转换为0-based
    values = data[:, 2]

    # 获取矩阵的大小
    matrix_size = max(max(rows), max(cols)) + 1

    # 创建一个稀疏矩阵
    sparse_matrix = scipy.sparse.coo_matrix((values, (rows, cols)), shape=(matrix_size, matrix_size))

    # 将稀疏矩阵转换为稠密格式
    dense_matrix = sparse_matrix.toarray()
    print("Is symmetric:", check_symmetry(sparse_matrix))

    # 可视化矩阵
    plt.figure(figsize=(10, 10))
    plt.imshow(dense_matrix, cmap='viridis', interpolation='none')
    plt.colorbar(label='Matrix Values')
    plt.title('Visualization of GLOBAL.K Matrix')
    plt.xlabel('Column Index')
    plt.ylabel('Row Index')
    plt.show()
