#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pathlib
import time
import sys
import numpy as np
from femdb.GlobalEnum import *
from ioclass.INPParser import InpParser
from ioclass.CDBParser import CDBParser
from ioclass.BDFParser import BDFParser
from ioclass.ResultsWriter import ResultsWriter
from femdb.Domain import Domain
from projects.qizhongji.zitai import MQ1330
from PySide2.QtGui import QVector3D
from PySide2.QtGui import QMatrix4x4

from flask import Flask, jsonify
app = Flask(__name__)


class MyPyFEM:
    """
    TODO: 可以将单元的面编号，对面做一个可识别的ID，用于区分，ID = str(sorted(nodeIds)), 参考MySTAP C++
    """

    def __init__(self, file_path, open_paraview=False, check_model=False, plot_stiff=False, is_server=False):
        if isinstance(file_path, str):
            file_path = pathlib.Path(file_path)

        if not os.path.isfile(file_path):
            logging.fatal("Input File is doesn't exist! {}".format(file_path))
            sys.exit(1)

        # 程序开始时间以及解析文件完成时间
        self.program_begin = None
        self.parsed_time = None
        self.check_model = check_model
        self.is_server = is_server
        self.domain = None

        # 是否绘制总刚度阵
        if plot_stiff:
            GlobalInfor[GlobalVariant.PlotGlobalStiffness] = True

        self._fem_data = None
        self.input_file_path = file_path
        self.output_files = [file_path.with_suffix(".vtu"),
                             file_path.with_suffix(".unv"),
                             file_path.with_suffix(".dat")]
        self.output_dir = file_path.parent
        self.output_name = file_path.stem

        self.FEMAnalyseFlow()

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
            return InpParser(self.input_file_path, self.check_model)
        elif suffix == ".cdb":
            GlobalInfor[GlobalVariant.InputFileSuffix] = InputFileType.CDB
            return CDBParser(self.input_file_path, self.check_model)
        elif suffix == ".bdf":
            GlobalInfor[GlobalVariant.InputFileSuffix] = InputFileType.BDF
            return BDFParser(self.input_file_path, self.check_model)
        else:
            mlogger.fatal("UnSupport File Suffix:{}".format(suffix))
            sys.exit(1)

    def FEMAnalyseFlow(self):
        """
        TODO: 标准流程, 完成注释, 重写mlogger的debug信息, 将有限元模型的信息输出, 比如单元类型及相应个数, 自由度个数
        求解文件, 步骤如下所示, 该函数中不应包含对不同文件类型的分类, 即判断文件类型的bdf cdb等应在其他函数中完成
        """
        self.program_begin = time.time()
        mlogger.debug("{} Analysis Calculate Begin {}".format("#" * 6, "#" * 6))
        reader = self.InitReader()
        reader.ParseFileAndInitFEMDB()
        if self.is_server:
            GlobalInfor[GlobalVariant.AnaType] = AnalyseType.AsServer
        self.parsed_time = time.time()

        """
        Print FEMDB Information
        """
        self.domain = Domain(self.check_model)
        summary = self.domain.femdb.GetModelSummary()
        mlogger.debug(" " + "-" * 40)
        summary_format = r"{:>25s} --> {:<}"
        mlogger.debug(" Model Summary:")
        for key, value in summary.items():
            mlogger.debug(summary_format.format(key, value))
        mlogger.debug(" " + "-" * 40)

        """
        Define Output Format And Print Each Step Time Elapsed
        """
        time_format = r"{:>25s} --> {:<.3f} seconds"
        last_line_format = "{:>25s} --> {:<.3f} seconds"

        mlogger.debug(" Elapsed Time Summary:")
        mlogger.debug(time_format.format("Parse File", self.parsed_time - self.program_begin))

        if GlobalInfor[GlobalVariant.AnaType] == AnalyseType.LinearStatic:
            self.domain.CalAllElementStiffness()
            time_1 = time.time()
            mlogger.debug(time_format.format("Calculate All Stiff", time_1 - self.parsed_time))

            self.domain.AssembleStiffnessMatrixByPenalty()
            time_2 = time.time()
            mlogger.debug(time_format.format("Assemble Global Stiff", time_2 - time_1))

            self.domain.AddBoundaryByPenalty()
            time_3 = time.time()
            mlogger.debug(time_format.format("Add Boundary Effect", time_3 - time_2))

            self.domain.SolveDisplacement()
            time_4 = time.time()
            mlogger.debug(time_format.format("Solve Displacement", time_4 - time_3))

            self.domain.SolveStress()
            time_5 = time.time()
            mlogger.debug(time_format.format("Solve Node Stress", time_5 - time_4))

            writer = ResultsWriter()
            writer.WriteStaticAnalysisVTUFile(self.output_files[0])
            # writer.WriteModel2DatFileWithoutRes(self.output_files[2])
            # writer.WriteStaticResult2DatFile(self.output_files[2])
            # writer.WriteMises2DatFile(self.output_files[2])
            p_end = time.time()
            mlogger.debug(time_format.format("Write Output", p_end - time_5))

        elif GlobalInfor[GlobalVariant.AnaType] == AnalyseType.Transient:
            """
            求解线弹性问题, 输出节点位移以及应力
            """
            self.domain.CalAllElementStiffness()
            time_1 = time.time()
            mlogger.debug(time_format.format("Calculate All Stiff", time_1 - self.parsed_time))

            self.domain.AssembleStiffnessMatrixByPenalty()
            time_2 = time.time()
            mlogger.debug(time_format.format("Assemble Global Stiff", time_2 - time_1))

            self.domain.CalAllElementMassMatrix()
            time_3 = time.time()
            mlogger.debug(time_format.format("Calculate All Mass", time_3 - time_2))

            self.domain.AssembleMassMatrixByPerturbation()
            time_4 = time.time()
            mlogger.debug(time_format.format("Assemble Global Mass", time_4 - time_3))

            self.domain.AddBoundaryByPenalty()
            time_5 = time.time()
            mlogger.debug(time_format.format("Add Boundary Effect", time_5 - time_4))

            self.domain.NewMarkExplict()
            time_6 = time.time()
            mlogger.debug(time_format.format("NewMark Analysis", time_6 - time_5))

            writer = ResultsWriter()
            writer.WriteSeriesResult(str(self.output_dir), str(self.output_name))
            p_end = time.time()
            mlogger.debug(time_format.format("Write Output", p_end - time_6))

        elif GlobalInfor[GlobalVariant.AnaType] == AnalyseType.AsServer:
            """
            有限元程序作为服务, 动态返回结果
            """
            # self.domain.CalAllElementStiffness()
            time_1 = time.time()
            mlogger.debug(time_format.format("Calculate All Stiff", time_1 - self.parsed_time))
            p_end = time.time()

        else:
            mlogger.fatal("UnSupport Analyse Type")
            sys.exit(1)

        mlogger.debug(last_line_format.format("Total Elapsed Time", p_end - self.program_begin))
        mlogger.debug(" " + "-" * 40)
        mlogger.debug(" Finish Analysis\n")

    def RotateModelAndSave2File(self, theta):
        """
        旋转模型, 然后保存
        :param theta: 旋转角度
        :return:
        """
        Abj = 43.224
        mq1330 = MQ1330()
        data = mq1330.getPartTranslate(Abj, Abj + theta)

        """
        先处理臂架、大拉杆、平衡梁的旋转
        """
        node_set_idx_hash = self.domain.femdb.node_set_name_hash
        node_sets = self.domain.femdb.node_sets

        bijia_nodes = node_sets[node_set_idx_hash["bijia_node"]]
        bijia_trans, bijia_theta, bijia_origin = data["bijia"]
        dalagan_nodes = node_sets[node_set_idx_hash["dalagan_node"]]
        dalagan_trans, dalagan_theta, dalagan_origin = data["dalagan"]
        phl_left_nodes = node_sets[node_set_idx_hash["phl_left_node"]]
        phl_right_nodes = node_sets[node_set_idx_hash["phl_right_node"]]
        phl_nodes = phl_left_nodes.extend(phl_right_nodes)
        phl_trans, phl_theta, phl_origin = data["pinghengliang"]


        """
        
        """
        for k, j in data.items():
            print(f"Name: {k} -> [{j[0]}, {j[1] * 180 / np.pi:.2f}, {j[2]}]")


    def RotateModel(self):
        """
        """
        femdb = self.domain.femdb
        print(len(femdb.node_list))
        return {
            "status": "success",
            "node_count": len(femdb.node_list)
        }


@app.route('/api/rotate_model', methods=['POST'])
def rotate_model_endpoint():
    """处理模型旋转请求的API端点"""
    global my_fem
    if not my_fem:
        return jsonify({"error": "FEM model not loaded"}), 400

    try:
        # 调用实例方法处理旋转请求
        result = my_fem.RotateModel()
        return jsonify(result), 200
    except Exception as e:
        logging.error(f"Rotation  failed: {str(e)}")
        return jsonify({"error": str(e)}), 500


if __name__ == "__main__":
    input_file = "./NumericalCases/Projects/qizhongji/last/MQ1330_remesh.cdb"
    my_fem = MyPyFEM(pathlib.Path(input_file), is_server=True)
    app.run(host='0.0.0.0', port=5000, debug=True, use_reloader=False)
    # app.run(host='0.0.0.0', port=5000, debug=True, reloader_type='watchdog')
