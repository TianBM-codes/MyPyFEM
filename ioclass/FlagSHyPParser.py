#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from femdb.NLFEMDataBase import NLFEMDataBase
from utils.GlobalEnum import *
from femdb.ElementFactory import ElementFactory, SetAnalyseDimension
from femdb.ElementGroup import ElementGroup
from femdb.Material import MaterialFactory
from femdb.LoadCase import FlagSHyPPressLoad, FlagSHyPCLoad
from femdb.SolveControl import FlagSHyPSolveControl


class FlagSHyPParser(object):
    """
    # TODO 这里需要节点排好，对应关系另外存储，不要每次都查询dict, 只有在读取输入文件和写结果的时候初始化各种dict
    """
    def __init__(self, input_path):
        self.fem_database = NLFEMDataBase()
        self.dat_path = input_path
        self.iter_line = None
        self.et_hash = {}
        self.ele_group_hash = {}
        self.ele_count = 0

    def ParseFileAndInitFEMDB(self):
        """
        文件格式：
        Javier Book P263
        :return:
        """
        GlobalInfor[GlobalVariant.AnaType] = AnalyseType.NLStatic
        fem_db = self.fem_database
        with open(self.dat_path, 'r', encoding='utf-8') as dat_file:
            """
            读取项目基础信息, 包括项目名称, 以及算例中涉及的单元, 整个项目是几维问题
            """
            fem_db.title = dat_file.readline()
            ele_type = dat_file.readline().strip()
            SetAnalyseDimension(ele_type)
            dim = GetDomainDimension()
            fem_db.SetDimensionVariant(dim)
            """
            读取节点个数、节点坐标以及边界条件, 三维问题, 每个节点有xyz三个坐标, 节点重排序
            """
            fem_db.Geom.node_count = int(dat_file.readline())
            fem_db.Mesh.n_dofs = dim * fem_db.Geom.node_count
            fem_db.BC.icode = np.zeros(fem_db.Geom.node_count, dtype=np.uint8)
            fem_db.Geom.x = np.zeros((GetDomainDimension(), fem_db.Geom.node_count), dtype=float)
            fem_db.Geom.x0 = np.zeros((GetDomainDimension(), fem_db.Geom.node_count), dtype=float)
            for ii in range(fem_db.Geom.node_count):
                node_line = dat_file.readline().strip().split()
                nodeId = int(node_line[0])
                fem_db.Geom.NdId2ListIdx[nodeId] = ii
                fem_db.Geom.ListIdx2NdId[ii] = nodeId
                fem_db.BC.icode[ii] = int(node_line[1])

                if len(node_line) == 5:
                    fem_db.Geom.x[:, ii] = [node_line[2], node_line[3], node_line[4]]
                    fem_db.Geom.x0[:, ii] = [node_line[2], node_line[3], node_line[4]]
                elif len(node_line) == 4:
                    fem_db.Geom.x0[:, ii] = [node_line[2], node_line[3]]
                else:
                    raise InputTextFormat(node_line)

            tmp = np.arange(dim * fem_db.Geom.node_count)
            fem_db.Mesh.dof_nodes = tmp.reshape((dim,fem_db.Geom.node_count), order='F')

            """
            读取单元信息, 依次是单元编号、材料编号以及包含节点ID(connectivity)
            """
            fem_db.Mesh.nelem = int(dat_file.readline())
            group = ElementGroup(ele_type)
            for ii in range(fem_db.Mesh.nelem):
                ele, ele_node_count = ElementFactory.CreateElement(ele_type)
                ele_line = dat_file.readline().strip().split()
                ele.id = int(ele_line[0])
                ele.mat_id = int(ele_line[1])
                ele.node_ids = np.asarray(ele_line[2:], dtype=np.uint32)
                ele.search_node_ids = np.asarray([fem_db.Geom.NdId2ListIdx[nd]
                                                  for nd in ele.node_ids], dtype=np.uint32)
                ele.e_type = ele_type
                ele.search_ele_idx = ii
                group.AppendElement(ele)

            """
            FlagSHyP只支持一种单元类型
            """
            fem_db.ElementGroupHash[0] = group

            """
            Obtain fixed and free degree of freedom numbers (dofs).
            """
            fem_db.BC.FindFixedAndFreeDofs()

            """
            Read the number of materials and material properties. 
            """
            mat_count = int(dat_file.readline())
            for ii in range(mat_count):
                mat_num, mat_type = dat_file.readline().split()
                mat = MaterialFactory.CreateMaterial(int(mat_type))
                mat.InitByFlagSHyPFormat(dat_file.readline())
                fem_db.Material.InsertMaterial(int(mat_num), mat)

            """
            Compute the number of nearly incompressible elements in the domain.
            """
            n_nearly_incompressible = 0
            for _, grp in fem_db.ElementGroupHash.items():
                for ele in grp.eles:
                    if fem_db.Material[ele.mat_id].name in [5, 7, 17]:
                        n_nearly_incompressible += 1

            fem_db.Material.n_nearly_incompressible = n_nearly_incompressible

            """
            Read nodal point loads, prescribed displacements, surface pressure loads
            and gravity (details in textbook).
            """
            load_split = dat_file.readline().split()
            concentrate_line_count = int(load_split[0])
            prescribed_dis_count = int(load_split[1])
            fem_db.LoadCase.n_pressure_loads = int(load_split[2])
            gravity = [float(load_split[3]), float(load_split[4]), float(load_split[5])]
            fem_db.LoadCase.AddGravity(gravity)

            for ii in range(concentrate_line_count):
                c_load = FlagSHyPCLoad(dat_file.readline())
                fem_db.LoadCase.AddFlagSHyPCLoad(c_load)

            for ii in range(prescribed_dis_count):
                pass

            for ii in range(fem_db.LoadCase.n_pressure_loads):
                p_load = FlagSHyPPressLoad(dat_file.readline())
                fem_db.LoadCase.AddFlagSHyPPressLoad(p_load)

            """
            Read control parameters.
            """
            solve_control = FlagSHyPSolveControl()
            solve_control.InitWithTextLine(dat_file.readline())
            fem_db.SolveControl = solve_control


if __name__ == "__main__":
    kkt = FlagSHyPParser("../NumericalCases/TwistingColumn/twisting_column.dat")
    kkt.ParseFileAndInitFEMDB()
