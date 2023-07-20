#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pathlib
import time
import sys
from femdb.GlobalEnum import *
from ioclass.INPParser import InpParser
from ioclass.CDBParser import CDBParser
from ioclass.BDFParser import BDFParser
from ioclass.ResultsWriter import ResultsWriter
from femdb.Domain import Domain


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

        # 是否绘制总刚度阵
        if plot_stiff:
            GlobalInfor[GlobalVariant.PlotGlobalStiffness] = True

        self._fem_data = None
        self.input_file_path = file_path
        output_files = [file_path.with_suffix(".vtu"), file_path.with_suffix(".unv")]

        self.RunAnalyseFlow(output_files, check_model)

        # 结果查看, Paraview显示, 注意要将paraview的路径加入至环境变量
        if open_paraview and not check_model:
            os.popen("paraview " + str(output_files[0].absolute()))

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
        else:
            mlogger.fatal("UnSupport File Suffix:{}".format(suffix))
            sys.exit(1)

    def RunAnalyseFlow(self, output_paths, check_model):
        """
        TODO: 标准流程, 完成注释, 重写mlogger的debug信息, 将有限元模型的信息输出, 比如单元类型及相应个数, 自由度个数
        求解文件, 步骤如下所示, 该函数中不应包含对不同文件类型的分类, 即判断文件类型的bdf cdb等应在其他函数中完成
        """
        mlogger.debug("{} Analysis Calculate Begin {}".format("#" * 6, "#" * 6))
        # logging.debug("Step 0: Parse File And Define the Problem FEMDB")
        p_begin = time.time()
        reader = self.InitReader()
        reader.ParseFileAndInitFEMDB()
        time_1 = time.time()

        if check_model:
            reader.CheckModel()
            c_end = time.time()
            mlogger.debug("Elapsed time: {:.3f} seconds\n".format(c_end - p_begin))
            return

        domain = Domain()
        if GlobalInfor[GlobalVariant.AnaType] == AnalyseType.LinearStatic:
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
        else:
            mlogger.fatal("UnSupport Analyse Type")
            sys.exit(1)

        writer = ResultsWriter()
        # writer.WriteVTPFile(output_paths[0])
        writer.WriteUNVFile(output_paths[1])
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
        mlogger.debug(time_format.format("Parse File", time_1 - p_begin))
        mlogger.debug(time_format.format("Calculate D", time_2 - time_1))
        mlogger.debug(time_format.format("Cal Equation Num", time_3 - time_2))
        mlogger.debug(time_format.format("Cal All Stiffness", time_4 - time_3))
        mlogger.debug(time_format.format("Assemble Stiffness", time_5 - time_4))
        mlogger.debug(time_format.format("Solve Displacement", time_6 - time_5))
        mlogger.debug(time_format.format("Write Output File", p_end - time_6))
        mlogger.debug(last_line_format.format("Total Elapsed Time", p_end - p_begin))
        mlogger.debug(" " + "-" * 40)
        mlogger.debug(" Finish Analysis\n")


if __name__ == "__main__":
    if len(sys.argv) == 1:
        mlogger.fatal(" No Input File Assign")
        sys.exit(1)
    input_file = sys.argv[1]
    MyPyFEM(pathlib.Path(input_file))
