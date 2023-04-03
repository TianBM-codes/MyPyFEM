#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import pathlib
import sys
from femdb.GlobalEnum import *
from fileparser.INPReader import InpReader
from fileparser.CDBReader import CDBReader
from fileparser.BDFReader import BDFReader
from femdb.Domain import Domain


class MyPyFEM:
    """
    TODO: 可以将单元的面编号，对面做一个可识别的ID，用于区分，ID = str(sorted(nodeIds)), 参考MySTAP C++
    """

    def __init__(self, file_path):
        if not os.path.isfile(file_path):
            logging.fatal("Input File is doesn't exist! {}".format(file_path))
            sys.exit(1)

        self._fem_data = None
        self.input_file_path = file_path

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

    def RunAnalyseFlow(self, output_path):
        """
        TODO: 标准流程, 完成注释
        求解文件, 步骤如下所示, 该函数中不应包含对不同文件类型的分类, 即判断文件类型的bdf cdb等应在其他函数中完成
        """
        logging.debug("Step 0: Parse File And Define the Problem FEMDB")
        reader = self.InitReader()
        reader.ParseFileAndInitFEMDB()

        domain = Domain()
        # domain.PrintFEMDBSummary()
        if GlobalInfor[GlobalVariant.AnaType] == AnalyseType.LinearStatic:
            logging.debug("Analyse Type: Linear Analyse")
            logging.debug("Step0: Prepare For Calculate Stiffness")
            domain.AssignElementCharacter()
            domain.CalBoundaryEffect()
            domain.CalculateEquationNumber()

            logging.debug("Step 1: Calculate Element StiffnessMatrix & Force Vector")
            domain.AssembleStiffnessMatrix()

            logging.debug("Step 2: Solve the Equation")
            domain.SolveDisplacement()
            # Step 4: BackSubstitution
            logging.debug("Step 5: Calculate & Output of All Element\n")
        else:
            mlogger.fatal("UnSupport Analyse Type")
            sys.exit(1)

        domain.WriteOutPutFile(output_path)


if __name__ == "__main__":
    """
    测试的文件
    cal = MyPyFEM("./tests/static/linear/truss/Job-1.inp")
    cal = MyPyFEM("./tests/write_vtk/data/Job-guagou.inp")
    cal = MyPyFEM("./tests/write_vtk/data/shellandsolid_santong.inp")
    cal = MyPyFEM("./tests/static/linear/solid/hexaC3D8/cube/Job-1.inp")
    input_file = pathlib.Path("./tests/static/linear/Plane/Job-1.inp")
    input_file = pathlib.Path("./tests/static/linear/truss/trussbridge/truss_static.inp")
    input_file = pathlib.Path("./tests/static/linear/truss/trusstower/tower.inp")
    input_file = pathlib.Path("./tests/ANSYS/bridge/bridge.cdb")
    input_file = pathlib.Path("./tests/ANSYS/linyi/yihe_bridge.cdb")
    input_file = pathlib.Path("./tests/ANSYS/beam/beam1882/beam1882.cdb")
    """

    input_file = pathlib.Path("./tests/ANSYS/beam/beam1881/one_beam.cdb")
    cal = MyPyFEM(input_file)
    output_file = input_file.with_suffix(".vtu")
    cal.RunAnalyseFlow(output_file)

    # 结果查看, Paraview显示, 注意要将paraview的路径加入至环境变量
    os.popen("paraview " + str(output_file.absolute()))
