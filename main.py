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
    # TODO:
    1. 可以将单元的面编号，对面做一个可识别的ID，用于区分，ID = str(sorted(nodeIds)), 参考MySTAP C++
    2. SigleTon.py研究
    3. meshio中包含很多格式解析以及互相转换的功能, 包括ansys和nastran等等，需要测试
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
            return InpReader(self.input_file_path)
        elif suffix == ".cdb":
            return CDBReader(self.input_file_path)
        elif suffix == ".bdf":
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
        # domain.Prepare4Calculate()

        logging.debug("Step 1: Calculate Element StiffnessMatrix & Force Vector")
        # domain.AssembleStiffnessMatrix()

        logging.debug("Step 2: Solve the Equation")
        # domain.SolveDisplacement()
        # Step 4: BackSubstitution
        # Step 5: Calculate & Output Stress of All Elements
        logging.debug("Step 5: Calculate & Output Stress of All Element\n\n")
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
    """

    input_file = pathlib.Path("./tests/ANSYS/linyi/yihe_bridge.cdb")
    cal = MyPyFEM(input_file)
    output_file = input_file.with_suffix(".vtu")
    cal.RunAnalyseFlow(output_file)

    # 结果查看, Paraview显示, 注意要将paraview的路径加入至环境变量
    os.popen("paraview " + str(output_file.absolute()))