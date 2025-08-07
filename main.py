#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pathlib
import time
import pickle
import sys
import numpy as np

from femdb.GlobalEnum import *
from ioclass.INPParser import InpParser
from ioclass.CDBParser import CDBParser
from ioclass.BDFParser import BDFParser
from ioclass.ResultsWriter import ResultsWriter
from femdb.Domain import Domain
from femdb.GlobalFEMVariant import ModelInfo
from projects.qizhongji.QiZhongJiZiTai import MQ1330Wrapper
import scipy.sparse as sparse
from fatigue_cal import fatigue_analysis

from flask import Flask, jsonify, request
import meshio

app = Flask(__name__)

"""
Define Output Format And Print Each Step Time Elapsed
"""
time_format = r"{:>25s} --> {:<.3f} seconds"
last_line_format = "{:>25s} --> {:<.3f} seconds"


class MyPyFEM:
    """
    TODO: 可以将单元的面编号，对面做一个可识别的ID，用于区分，ID = str(sorted(nodeIds)), 参考MySTAP C++
    """

    def __init__(self, file_path, open_paraview=False, check_model=False, plot_stiff=False, AnaType=None, output_path=None):
        if isinstance(file_path, str):
            file_path = pathlib.Path(file_path)

        if not os.path.isfile(file_path):
            logging.fatal("Input File is doesn't exist! {}".format(file_path))
            sys.exit(1)

        # 程序开始时间以及解析文件完成时间
        self.program_begin = None
        self.parsed_time = None
        self.check_model = check_model
        self.ana_type = AnaType
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
        self.output_path = output_path

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
        """
        指定分析类型
        """
        if self.ana_type:
            GlobalInfor[GlobalVariant.AnaType] = self.ana_type
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
            if self.output_path is not None:
                writer.WriteStaticAnalysisVTUFile(self.output_path)
            else:
                writer.WriteStaticAnalysisVTUFile(self.output_files[0])
            # writer.WriteModel2DatFileWithoutRes(self.output_files[2])

            # wrapper = MQ1330Wrapper(30)
            # writer.WriteStaticResult2DatFile2(self.output_files[2], wrapper)
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
            self.domain.CalAllElementStiffness()
            time_1 = time.time()
            mlogger.debug(time_format.format("Calculate All Stiff", time_1 - self.parsed_time))
            p_end = time.time()

        elif GlobalInfor[GlobalVariant.AnaType] == AnalyseType.ReSortModel:
            """
            只是重排单元, 生成Dat文件
            """
            writer = ResultsWriter()
            writer.WriteResortModel2DatFile(self.output_files[2])
            p_end = time.time()

        elif GlobalInfor[GlobalVariant.AnaType] == AnalyseType.GenerateGeLinFunction:
            """
            生成格林函数文件
            """
            self.domain.CalAllElementStiffness()
            time_1 = time.time()
            mlogger.debug(time_format.format("Calculate All Stiff", time_1 - self.parsed_time))

            self.domain.AssembleStiffnessMatrixByPenalty()
            time_2 = time.time()
            mlogger.debug(time_format.format("Assemble Global Stiff", time_2 - time_1))

            self.domain.AddBoundaryByPenalty()
            time_3 = time.time()
            mlogger.debug(time_format.format("Add Boundary Effect", time_3 - time_2))

            """
            创建目录, 然后保存格林函数文件
            """
            green_dir = pathlib.Path(self.output_dir / "green")
            if not green_dir.exists():
                pathlib.Path(self.output_dir / "green").mkdir()
            self.domain.CalculateGreenFunction(green_dir)
            p_end = time.time()
            mlogger.debug(time_format.format("Generate Green Function", p_end - time_3))

        elif GlobalInfor[GlobalVariant.AnaType] == AnalyseType.CalculateByGeLin:
            """
            通过格林函数计算模型结果
            """
            green_dir = pathlib.Path(self.output_dir / "green")
            if not green_dir.exists():
                raise ImportError("./green directory Don't exists")
            self.domain.CalculateResultsWithGreenFunction(green_dir)
            writer = ResultsWriter()
            writer.WriteStaticAnalysisVTUFile(self.output_files[0])
            p_end = time.time()

        elif GlobalInfor[GlobalVariant.AnaType] == AnalyseType.GenerateDynGeLinFunction:
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

            green_dir = pathlib.Path(self.output_dir / "green")
            if not green_dir.exists():
                pathlib.Path(self.output_dir / "green").mkdir()
            self.domain.GenerateNewMarkGeLinFile(green_dir)

            p_end = time.time()
            mlogger.debug(time_format.format("Write Output", p_end - time_5))

        elif GlobalInfor[GlobalVariant.AnaType] == AnalyseType.CalculateDynByGeLin:
            """
            通过格林函数计算动力学模型结果
            """
            green_dir = pathlib.Path(self.output_dir / "green")
            if not green_dir.exists():
                raise ImportError("./green directory Don't exists")
            self.domain.CalculateDynResultsWithGreenFunction(green_dir)
            p_end = time.time()

        elif GlobalInfor[GlobalVariant.AnaType] == AnalyseType.FatigueAnalysis:
            """
            计算疲劳云图
            """
            Su = 300  # 材料抗拉强度 (MPa)
            material_params = {
                'S_endurance': 17,  # 疲劳极限 (MPa)
                'a': 1e12,  # Basquin方程参数a
                'b': -3  # Basquin方程参数b
            }
            mises_path = pathlib.Path(self.output_files[0])
            all_mises = []
            for ii in range(100):
                iter_path = str(mises_path.parent) + "/" + str(mises_path.stem) + f"_{ii}.vtu"
                iter_vtu = meshio.read(iter_path)
                iter_mises = iter_vtu.point_data['mises']
                all_mises.append(iter_mises)

            all_mises_array = np.array(all_mises) * 100
            all_damage = []
            for ii in range(all_mises_array.shape[1]):
                iter_result = fatigue_analysis(all_mises_array[:, ii], Su, material_params)
                iter_damage = iter_result['total_damage']
                all_damage.append(iter_damage)

            self.domain.femdb.damage_factor = np.array(all_damage)
            input_file_stem = pathlib.Path(self.output_files[0]).stem
            fatigure_path = pathlib.Path(self.output_files[0]).with_name(input_file_stem+"_fat.vtu")
            writer = ResultsWriter()
            writer.WriteFatigueResult(fatigure_path)
            p_end = time.time()

        else:
            mlogger.fatal("UnSupport Analyse Type")
            sys.exit(1)

        mlogger.debug(last_line_format.format("Total Elapsed Time", p_end - self.program_begin))
        mlogger.debug(" " + "-" * 40)
        mlogger.debug(" Finish Analysis\n")

    def ReCalculateFEMModel(self, e_value):
        """
        更新弹性模型, 重新计算有限元模型
        :param e_value:
        :return:
        """
        time1 = time.time()
        femdb = self.domain.femdb
        for iter_ele in femdb.elements:
            iter_ele.cha_dict[MaterialKey.E] = e_value
        self.domain.AssembleStiffnessMatrixByPenalty(from_origin=True)
        time2 = time.time()
        mlogger.debug(time_format.format("ReCalculate All Stiffness", time2 - time1))
        self.domain.AddBoundaryByPenalty()
        time3 = time.time()
        mlogger.debug(time_format.format("Add Boundary Effect", time3 - time2))

        self.domain.SolveDisplacement()
        time4 = time.time()
        mlogger.debug(time_format.format("Solve Displacement", time4 - time3))

        writer = ResultsWriter()
        writer.WriteStaticAnalysisVTUFile(self.output_files[0])

    def RotateModel(self, theta, vtu_path):
        """
        旋转模型, 以臂架为基准
        :param theta:
        :param vtu_path:
        :return:
        """
        time1 = time.time()
        femdb = self.domain.femdb
        wrapper = MQ1330Wrapper(theta)
        for nd in femdb.node_list:
            new_coord = np.array(wrapper.calculate_new_xyz(nd.id, *nd.origin_coord))
            nd.coord = new_coord

        time2 = time.time()
        mlogger.debug(time_format.format("Rotate XYZ", time2 - time1))

        for ii, stiff in enumerate(femdb.stiff_list):
            iter_ele = femdb.elements[ii]
            ele_id = iter_ele.id
            node_count = iter_ele.nodes_count
            if not wrapper.is_fixed_element(ele_id):
                if node_count == 3 or node_count == 4:
                    new_stiff = wrapper.calculate_new_stiff(ele_id, node_count, stiff)
                    femdb.stiff_list[ii] = new_stiff

        time3 = time.time()
        mlogger.debug(time_format.format("Calculate New Stiff", time3 - time2))

        self.domain.AssembleStiffnessMatrixByPenalty()
        time4 = time.time()
        mlogger.debug(time_format.format("Assemble Stiffness", time4 - time3))

        self.domain.AddBoundaryByPenalty()
        time5 = time.time()
        mlogger.debug(time_format.format("Add Boundary Effect", time5 - time4))

        self.domain.SolveDisplacement()
        time6 = time.time()
        mlogger.debug(time_format.format("Solve Displacement", time6 - time5))

        self.domain.SolveStress()
        time7 = time.time()
        mlogger.debug(time_format.format("Solve Stress", time7 - time6))

        writer = ResultsWriter(use_mysql=True)
        src = pathlib.Path(vtu_path)
        # writer.WriteStaticAnalysisVTUFile(src)
        dat_path = src.with_suffix(".dat")
        dat_path = dat_path.with_stem(dat_path.stem + f"{int(time.time())}")
        writer.WriteStaticResult2DatFile2(dat_path, wrapper)
        # writer.WriteMises2DatFile(dat_path.with_stem(dat_path.stem + "_mises"))
        time_end = time.time()
        total_time_elapsed = time_end - time1
        mlogger.debug(time_format.format("Write Output", time_end - time7))
        mlogger.debug(time_format.format("Total Time Elapsed", total_time_elapsed))
        return {
            "status": "success",
            "time_elapsed": f"{total_time_elapsed:.2f}"
        }


@app.route('/api/rotate_model', methods=['POST'])
def rotate_model_endpoint():
    """
    处理模型旋转请求的API端点
    :return:
    """
    global my_fem
    if not my_fem:
        return jsonify({"error": "FEM model not loaded"}), 400

    try:
        theta = request.args.get('rotate_theta', type=float)
        save_path = pathlib.Path(request.args.get('save_path'))
        current = 1
        step = 1
        while True:
            if current >= 33:
                step = -2
            elif current <= 1:
                step = 2
            current += step
            iter_path = save_path.with_stem(f"theta{current}")
            my_fem.RotateModel(current, iter_path)

        result = {"status": "success"}
        return jsonify(result), 200

    except Exception as e:
        logging.error(f"Rotation  failed: {str(e)}")
        return jsonify({"error": str(e)}), 500


@app.route('/api/re_calculate_stiff', methods=['POST'])
def re_calculate_element_stiff():
    global my_fem
    if not my_fem:
        return jsonify({"error": "FEM model not loaded"}), 400

    try:
        e_value = request.args.get('e_value', type=float)
        my_fem.ReCalculateFEMModel(e_value)
        result = {"status": "success"}
        return jsonify(result), 200
    except Exception as e:
        return jsonify({"error": str(e)}), 500


if __name__ == "__main__":
    # input_file = "./NumericalCases/Projects/qizhongji/last/MQ1330_remesh.cdb"
    # input_file = r"D:\WorkSpace\FEM\MyPyFEM\numerical example\ANSYS\zhijiaRenumber.cdb"
    # input_file = r"D:\WorkSpace\FEM\MyPyFEM\numerical example\ANSYS\triAndQuaCylinder.cdb"
    # input_file = r"D:\WorkSpace\FEM\MyPyFEM\numerical example\ANSYS\allTriCylinder.cdb"
    # input_file = r"D:\WorkSpace\FEM\testcases\ANSYS\shell\singleInclineQuaShell.cdb"
    input_file = r"D:\WorkSpace\WebThreeJS\PyModelToJson\model\ansys\cdb\MQ1330_remesh_resort.cdb"
    # my_fem = MyPyFEM(pathlib.Path(input_file), AnaType=AnalyseType.AsServer)
    my_fem = MyPyFEM(pathlib.Path(input_file), AnaType=AnalyseType.ReSortModel)
    # app.run(host='0.0.0.0', port=5000, debug=True, use_reloader=False)
    # app.run(host='0.0.0.0', port=5000, debug=True, reloader_type='watchdog')

    # my_fem.ReCalculateFEMModel(2.1e11)
