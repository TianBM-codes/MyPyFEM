#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import pathlib
import time
import sys
from femdb.GlobalEnum import *
from ioclass.INPParser import InpReader
from ioclass.CDBParser import CDBReader
from ioclass.BDFParser import BDFReader
from ioclass.ResultsWriter import ResultsWriter
from femdb.Domain import Domain


class MyPyFEM:
    """
    TODO: 可以将单元的面编号，对面做一个可识别的ID，用于区分，ID = str(sorted(nodeIds)), 参考MySTAP C++
    """

    def __init__(self, file_path, open_paraview=False):
        if isinstance(file_path, str):
            file_path = pathlib.Path(file_path)

        if not os.path.isfile(file_path):
            logging.fatal("Input File is doesn't exist! {}".format(file_path))
            sys.exit(1)

        self._fem_data = None
        self.input_file_path = file_path
        output_files = [file_path.with_suffix(".vtu"), file_path.with_suffix(".unv")]

        self.RunAnalyseFlow(output_files)

        # 结果查看, Paraview显示, 注意要将paraview的路径加入至环境变量
        if open_paraview:
            os.popen("paraview " + str(output_files[0].absolute()))

    def InitReader(self):
        """
        根据文件类型初始化不同的文件解析器, 然后读入文件初始化数据库
        """
        suffix = self.input_file_path.suffix
        if suffix == ".inp":
            GlobalInfor[GlobalVariant.InputFileSuffix] = InputFileType.INP
            return InpReader(self.input_file_path)
        elif suffix == ".cdb":
            GlobalInfor[GlobalVariant.InputFileSuffix] = InputFileType.CDB
            return CDBReader(self.input_file_path)
        elif suffix == ".bdf":
            GlobalInfor[GlobalVariant.InputFileSuffix] = InputFileType.BDF
            return BDFReader(self.input_file_path)
        else:
            mlogger.fatal("UnSupport File Suffix:{}".format(suffix))
            sys.exit(1)

    def RunAnalyseFlow(self, output_paths):
        """
        TODO: 标准流程, 完成注释
        求解文件, 步骤如下所示, 该函数中不应包含对不同文件类型的分类, 即判断文件类型的bdf cdb等应在其他函数中完成
        """
        mlogger.debug("{} Analysis Calculate Begin {}".format("#" * 6, "#" * 6))
        logging.debug("Step 0: Parse File And Define the Problem FEMDB")
        p_begin = time.time()
        reader = self.InitReader()
        reader.ParseFileAndInitFEMDB()

        domain = Domain()
        # domain.PrintFEMDBSummary()
        if GlobalInfor[GlobalVariant.AnaType] == AnalyseType.LinearStatic:
            mlogger.debug("Analyse Type: Linear Analyse")
            domain.AssignElementCharacter()
            domain.CalBoundaryEffect()
            domain.CalculateEquationNumber()
            p_end = time.time()

            mlogger.debug("Step 1: Calculate Element StiffnessMatrix & Force Vector")
            domain.AssembleStiffnessMatrix()

            mlogger.debug("Step 2: Solve the Equation")
            domain.SolveDisplacement()
            d_time = time.time()
            # Step 4: BackSubstitution
            logging.debug("Step 5: Calculate & Output of All Element")
        else:
            mlogger.fatal("UnSupport Analyse Type")
            sys.exit(1)

        writer = ResultsWriter()
        writer.WriteVTPFile(output_paths[0])
        # writer.WriteUNVFile(output_paths[1])

        # Output Each Step Time Elapsed
        w_time = time.time()
        time_format = r"{:>20s} --> {:<.3f} seconds"
        mlogger.debug("Elapsed Time Summary:")
        mlogger.debug(time_format.format("Parse File", p_end - p_begin))
        mlogger.debug(time_format.format("Solved Displacement", d_time - p_end))
        last_line_format = "{:>20s} --> {:<.3f} seconds\n"
        mlogger.debug(last_line_format.format("Write Output File", w_time - d_time))


if __name__ == "__main__":
    cal = MyPyFEM(pathlib.Path("./numerical example/ANSYS/shell/FourShell.cdb"))
