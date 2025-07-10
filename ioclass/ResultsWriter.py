#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pathlib
import meshio
from scipy.ndimage import distance_transform_bf

from element.MeshElementFactory import *
from femdb.FEMDataBase import *
from ioclass.MySqlMathFunction import MySQLMathFunction


class ResultsWriter(object):
    """
    计算结果导出类, 现支持UNV导出以及VTP导出
    """

    def __init__(self, use_mysql=False):
        self.femdb = FEMDataBase()
        self.use_mysql = use_mysql
        if use_mysql:
            self.mysql_db = MySQLMathFunction()

    def WriteStaticAnalysisVTUFile(self, path, reshape=True):
        """
        将结果写至vtu文件
        Reference:
        1. https://github.com/nschloe/meshio
        """
        # 模型部分
        coords = np.asarray([node.coord for node in self.femdb.node_list])
        all_eles = {}

        for iter_ele in self.femdb.elements:
            iter_relation = iter_ele.GetNodeSearchIndex().tolist()
            ele_type = iter_ele.vtu_type
            if all_eles.__contains__(ele_type):
                all_eles[ele_type].append(iter_relation)
            else:
                all_eles[ele_type] = [iter_relation]

        for iter_ele_info in self.femdb.additional_elements:
            iter_relation = [iter_ele_info[0], iter_ele_info[1]]
            ele_type = iter_ele_info[2]
            if all_eles.__contains__(ele_type):
                all_eles[ele_type].append(iter_relation)
            else:
                all_eles[ele_type] = [iter_relation]

        # 位移结果
        if reshape:
            dis_value = np.reshape(self.femdb.linear_u, (-1, ModelInfo.PER_NODE_DOF))[:, :3]
        else:
            dis_value = self.femdb.linear_u
        node_res = {"displacement": dis_value, "mises": self.femdb.linear_mises}
        meshio.write_points_cells(
            filename=path,
            points=coords,
            cells=all_eles,
            point_data=node_res,
        )

    def WriteSeriesResult(self, directory, pro_name):
        """
        将时序结果写入文件
        :param directory: 存储路径
        :param pro_name: 项目名称
        :return:
        """
        coords = np.asarray([node.coord for node in self.femdb.node_list])
        all_eles = {}

        for iter_ele in self.femdb.elements:
            iter_relation = iter_ele.GetNodeSearchIndex().tolist()
            ele_type = iter_ele.vtu_type
            if all_eles.__contains__(ele_type):
                all_eles[ele_type].append(iter_relation)
            else:
                all_eles[ele_type] = [iter_relation]

        for ii in range(self.femdb.history_step_count):
            path = directory + "/" + pro_name + f"_{ii}.vtu"
            u = self.femdb.history_u[ii]
            v = self.femdb.history_v[ii]
            a = self.femdb.history_a[ii]
            time_result = {"displacement": np.reshape(u, (-1, ModelInfo.PER_NODE_DOF))[:, :3],
                           "velocity": np.reshape(v, (-1, ModelInfo.PER_NODE_DOF))[:, :3],
                           "acceleration": np.reshape(a, (-1, ModelInfo.PER_NODE_DOF))[:, :3]}
            meshio.write_points_cells(
                filename=path,
                points=coords,
                cells=all_eles,
                point_data=time_result,
            )


if __name__ == "__main__":
    file_path = "D:/WorkSpace/FEM/MyPyFEM/NumericalCases/Projects/qizhongji/last/MQ1330_remesh.cdb"
    file_lib = pathlib.Path(file_path)
    new_name = file_lib.stem + "_mises"
    new_path = file_lib.with_stem(new_name)
    print(new_path)
