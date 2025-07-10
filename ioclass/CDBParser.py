#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from femdb.FEMDataBase import *
from collections import OrderedDict
from element.Node import Node
from femdb.ElementFactory import *
from element.Beam import BeamCalculator
from femdb.GlobalFEMVariant import ModelInfo

import numpy as np
import time
import re
import fortranformat as ff


class CDBParser(object):
    """
    ANSYS CDB文件解析器

    功能:
    - 解析ANSYS CDB格式的有限元模型文件
    - 初始化FEM数据库(节点、单元、材料、边界条件等)
    - 支持多种ANSYS命令和格式

    特性:
    - 仅负责文件解析和数据传递，不包含计算功能
    - 支持节点、单元、材料、截面、边界条件等解析
    - 处理多种ANSYS数据格式(Fortran格式)

    参考:
    - ANSYS Help - Mechanical APDL Element Reference
    - ANSYS Programmer's Manual
    """

    def __init__(self, input_path, check_model):
        """
        初始化CDB解析器

        Args:
            input_path (str): CDB文件路径
            check_model (bool): 是否检查模型完整性
        """
        self.femdb = FEMDataBase()  # 有限元数据库实例
        self.cdb_path = input_path  # 输入文件路径
        self.iter_line = None  # 当前读取的行
        self.et_hash = {}  # 单元类型哈希表
        self.ele_count = 0  # 单元计数器
        self.check_model = check_model  # 模型检查标志
        self.real_constant_hash = {}  # 实常数哈希表
        self.material_map = {}  # 材料属性映射表
        self.section_map = {}  # 截面属性映射表
        self.node_search_ids_list = []  # 节点搜索ID列表

    def ParseFileAndInitFEMDB(self):
        """
        解析CDB文件并初始化有限元数据库

        处理流程:
        1. 打开文件并逐行读取
        2. 根据行首关键字调用相应解析方法
        3. 将解析结果存入FEM数据库
        4. 设置单元属性和材料参数

        文件格式说明:
        1. 字母开头: ANSYS命令关键字
        2. 括号开头: Fortran格式说明符
        3. *开头: 其他信息(暂不解析)
        4. /开头: 注释
        """
        self.femdb.file_path = self.cdb_path
        GlobalInfor[GlobalVariant.AnaType] = AnalyseType.LinearStatic  # 默认静力分析

        with open(self.cdb_path, 'r') as cdb_f:
            self.iter_line = cdb_f.readline()
            while True:
                # 解析分析类型
                if self.iter_line.startswith("ANTYPE,"):
                    if self.iter_line.split(",")[1].strip() == "0":
                        GlobalInfor[GlobalVariant.AnaType] = AnalyseType.LinearStatic
                    else:
                        mlogger.fatal("UnSupport Analyse Type:{}".format(self.iter_line))
                        sys.exit(1)
                    self.iter_line = cdb_f.readline()

                # 解析单元类型
                elif self.iter_line.startswith("ET,"):
                    splits = self.iter_line.split(",")
                    assert len(splits) == 3
                    e_type = int(splits[2].strip())
                    self.et_hash[int(splits[1].strip())] = e_type
                    self.iter_line = cdb_f.readline()
                    dof_count = ElementFactory.GetElementNodeDofCount(e_type)
                    if dof_count is not None:
                        ModelInfo.PER_NODE_DOF = dof_count
                    if dof_count == 2:
                        self.femdb.an_dimension = AnalyseDimension.TwoDimension

                # 解析实常数
                elif self.iter_line.startswith("RLBLOCK,"):
                    """
                    RLBLOCK命令格式参考:
                    https://ansyshelp.ansys.com/public/account/secured?returnurl=//////Views/Secured/corp/v242/en/ans_prog/Hlp_P_INT3_3.html%23eLN4r40lcd

                    格式说明:
                    - 不以字母开头的行表示数据定义未结束
                    - 两种固定格式: (2i8,6g16.9)和(7g16.9)
                    - 当前版本不解析(7g16.9)格式的第二行
                    """
                    self.iter_line = cdb_f.readline()
                    format_list = []
                    while not self.iter_line[0].isalpha():
                        if self.iter_line[0] == "(":
                            format_list.append(ff.FortranRecordReader(self.iter_line.strip()))
                            self.iter_line = cdb_f.readline()
                        elif self.iter_line[0] == " ":
                            if len(self.iter_line.split()) == 8:
                                rl_data = format_list[0].read(self.iter_line)
                                self.real_constant_hash[rl_data[0]] = rl_data[2:]
                            self.iter_line = cdb_f.readline()
                        elif self.iter_line[0] == "!":
                            # 兼容hypermesh生成的cdb
                            self.iter_line = cdb_f.readline()

                # 解析节点块
                elif self.iter_line.startswith("NBLOCK,"):
                    node_index = 0  # 节点计数器，对应数据库中的索引
                    fortran_format = cdb_f.readline().strip()  # 节点行格式
                    if fortran_format not in ['(3i8,6e16.9)', '(3i9,6e21.13e3)']:
                        raise ValueError(f"UnSupport node format:{fortran_format}")
                    self.iter_line = cdb_f.readline()  # 读取第一行节点数据

                    while self.iter_line.startswith(" "):
                        if fortran_format == '(3i8,6e16.9)':
                            # 解析8位整数和16位浮点数格式
                            int_part = self.iter_line[:24]
                            integers = [int(int_part[i * 8:(i + 1) * 8]) for i in range(3) if int_part[i * 8:(i + 1) * 8].strip()]
                            float_part = self.iter_line[24:].rstrip()
                            floats = [
                                float(float_part[i * 16:(i + 1) * 16])
                                for i in range(min(6, len(float_part) // 16))  # 自动计算有效浮点数
                                if float_part[i * 16:(i + 1) * 16].strip()
                            ]

                            n_id = integers[0]
                            x, y, z = floats[0], floats[1] if len(floats) > 1 else 0, floats[2] if len(floats) > 2 else 0

                        elif fortran_format == '(3i9,6e21.13e3)':
                            # 解析9位整数和21位浮点数格式
                            def parse_fortran_float(chunk):
                                return float(re.sub(r'e([+-]\d{2})\d', r'e\1', chunk))

                            integers = [int(self.iter_line[i * 9:(i + 1) * 9]) for i in range(3) if self.iter_line[i * 9:(i + 1) * 9].strip()]
                            float_part = self.iter_line[27:].rstrip()
                            floats = []
                            for i in range(0, len(float_part), 21):
                                chunk = float_part[i:i + 21]
                                if chunk.strip():  # 非空校验
                                    try:
                                        floats.append(parse_fortran_float(chunk))
                                    except ValueError:
                                        print(f"警告：跳过无效浮点字段 '{chunk}'")

                            n_id = integers[0]
                            x = floats[0] if floats else 0
                            y = floats[1] if len(floats) > 1 else 0
                            z = floats[2] if len(floats) > 2 else 0

                        else:
                            raise ValueError(f"Unsupported Node Format: {fortran_format}")

                        # 添加节点到数据库
                        self.femdb.AddNode(Node(n_id, x, y, z))
                        self.femdb.node_hash[n_id] = node_index
                        node_index += 1
                        self.iter_line = cdb_f.readline()

                # 解析单元块
                elif self.iter_line.startswith("EBLOCK,"):
                    self.ReadEBlock(cdb_f)

                # 解析组件块
                elif self.iter_line.startswith("CMBLOCK"):
                    splits = self.iter_line.strip().split(",")
                    comp_name = splits[1]
                    comp_type = splits[2]
                    if "!" in splits[3]:
                        set_data_count = int(splits[3].split("!")[0])
                    else:
                        set_data_count = int(splits[3])
                    fortran_format = cdb_f.readline().strip()  # 组件数据格式
                    set_data = []
                    if fortran_format == '' or set_data_count == 0:
                        self.iter_line = cdb_f.readline()
                    else:
                        f_reader = ff.FortranRecordReader(fortran_format)
                        set_data_line = f_reader.read(cdb_f.readline())
                        input_times = 0
                        while True:
                            if len(set_data_line) != 0:
                                iter_v = set_data_line.pop(0)

                            if iter_v is not None:
                                if iter_v > 0:
                                    set_data.append(iter_v)
                                    input_times += 1
                                else:
                                    set_data.extend(list(range(set_data[-1] + 1, -iter_v + 1)))
                                    input_times += 1

                            if len(set_data_line) == 0 and input_times < set_data_count:
                                set_data_line = f_reader.read(cdb_f.readline())

                            if iter_v is None:
                                if input_times < set_data_count:
                                    self.iter_line = cdb_f.readline()
                                    set_data_line = f_reader.read(self.iter_line)
                                else:
                                    self.iter_line = cdb_f.readline()
                                    break

                    # 将组件添加到相应集合
                    if comp_type == 'NODE':
                        self.femdb.node_set_name_hash[comp_name] = len(self.femdb.node_sets)
                        self.femdb.node_sets.append(set_data)
                    elif comp_type == 'ELEMENT' or comp_type == 'ELEM':
                        self.femdb.element_set_name_hash[comp_name] = len(self.femdb.element_sets)
                        self.femdb.element_sets.append(set_data)
                    else:
                        raise KeyError(f"CMBLOCK Key: {comp_type}")

                # 解析材料属性
                elif self.iter_line.startswith("MPDATA,"):
                    self.ReadMaterial(cdb_f, "MPDATA,")

                elif self.iter_line.startswith("MP,"):
                    self.ReadMaterial(cdb_f, "MP,")

                # 解析截面属性
                elif self.iter_line.startswith("SECTYPE,"):
                    self.ReadSection(cdb_f)

                # 解析加速度
                elif self.iter_line.startswith("ACEL,"):
                    # TODO: 处理加速度, 重力
                    self.iter_line = cdb_f.readline()

                # 解析刚性约束
                elif self.iter_line.startswith("CERIG"):
                    splits = self.iter_line.strip().split(",")
                    m_node = int(splits[1])
                    s_node = int(splits[2])
                    self.femdb.equation_constrain_couple.append((m_node, s_node))
                    self.femdb.additional_elements.append((self.femdb.node_hash[m_node], self.femdb.node_hash[s_node], "line"))
                    constrain_type = splits[3].lstrip()
                    if constrain_type == 'UXYZ':
                        self.femdb.equation_constrain_idx.append([0, 1, 2, 3, 4, 5])
                    elif constrain_type == "ALL":
                        self.femdb.equation_constrain_idx.append([0, 1, 2, 3, 4, 5])
                    else:
                        raise KeyError(f"{constrain_type}")
                    self.iter_line = cdb_f.readline()

                # 解析位移约束
                elif self.iter_line.startswith("D,"):
                    # 可能包含多个自由度的约束
                    d_nodes, directs, values = [], [], []
                    while self.iter_line.startswith("D,"):
                        splits = self.iter_line.split(",")
                        d_node = [int(splits[1])]
                        d_val = [float(splits[3].strip())]
                        d_dir = splits[2].strip()
                        if "UX" in d_dir:
                            dir_idx = [0]
                        elif "UY" in d_dir:
                            dir_idx = [1]
                        elif "UZ" in d_dir:
                            dir_idx = [2]
                        elif "ROTX" in d_dir:
                            dir_idx = [3]
                        elif "ROTY" in d_dir:
                            dir_idx = [4]
                        elif "ROTZ" in d_dir:
                            dir_idx = [5]
                        elif "ALL" in d_dir:
                            if ModelInfo.PER_NODE_DOF == 2:
                                dir_idx = [0, 1]
                                d_node = [int(splits[1])] * 2
                                d_val = [float(splits[3].strip())] * 2
                            elif ModelInfo.PER_NODE_DOF == 3:
                                dir_idx = [0, 1, 2]
                                d_node = [int(splits[1])] * 3
                                d_val = [float(splits[3].strip())] * 3
                            else:
                                dir_idx = [0, 1, 2, 3, 4, 5]
                                d_node = [int(splits[1])] * 6
                                d_val = [float(splits[3].strip())] * 6
                        else:
                            raise KeyError(f"Boundary Type:{d_dir}")

                        d_nodes.extend(d_node)
                        directs.extend(dir_idx)
                        values.extend(d_val)
                        self.iter_line = cdb_f.readline()

                    # 添加边界条件到载荷工况
                    bd = np.column_stack([d_nodes, directs, values])
                    self.femdb.load_case.AddBoundary(bd)

                # 解析集中力
                elif self.iter_line.startswith("F,"):
                    # 可能包含多个集中力
                    while self.iter_line.startswith("F,"):
                        splits = self.iter_line.split(",")
                        if "FX" in splits[2]:
                            idx = 0
                        elif "FY" in splits[2]:
                            idx = 1
                        elif "FZ" in splits[2]:
                            idx = 2
                        else:
                            raise KeyError(f"Don't Support Force Key:{splits[2]}")
                        self.femdb.load_case.AddConcentratedLoad(int(splits[1]), idx, float(splits[3].strip()))
                        self.iter_line = cdb_f.readline()

                else:
                    # 文件结束
                    if not self.iter_line:
                        break
                    # 跳过不支持的关键字或注释
                    self.iter_line = cdb_f.readline()

        """
        设置单元属性和材料参数
        在文件解析完成后统一设置，确保所有计算单元矩阵的参数均已准备就绪
        """
        for iter_ele in self.femdb.elements:
            mat_dict = self.material_map[iter_ele.mat_id]
            real_const_id = iter_ele.real_const_id
            real_const = {"RealConst": self.real_constant_hash.get(real_const_id, {})}
            sec_vals = self.section_map.get(iter_ele.sec_id, {})
            iter_ele.SetAllCharacterAndCalD({**mat_dict, **real_const, **sec_vals})

        """
        计算每个节点连接的单元数量
        用于应力平均计算
        """
        _, node_connected_element_count = np.unique(self.node_search_ids_list, return_counts=True)
        self.femdb.node_connected_element_count = node_connected_element_count

    def ReadEBlock(self, f_handle):
        """
        解析单元块(EBLOCK)数据

        SOLID格式的单元块格式说明:
        - 字段1: 材料号
        - 字段2: 单元类型号
        - 字段3: 实常数号
        - 字段4: 截面ID属性(梁截面)
        - 字段5: 单元坐标系号
        - 字段6: 生死标志
        - 字段7: 实体模型参考号
        - 字段8: 单元形状标志
        - 字段9: 定义该单元的节点数(SOLID关键字时为实际值，否则为0)
        - 字段10: 未使用
        - 字段11: 单元号
        - 字段12-19: 节点号(超过8个节点时在下一行继续)
        """
        solid_type = True
        if not self.iter_line.__contains__("SOLID"):
            solid_type = False

        fortran_format = f_handle.readline().strip()  # 单元行格式
        if fortran_format == '(19i8)':
            chunk_size = 8
        elif fortran_format == '(19i9)':
            chunk_size = 9
        elif fortran_format == '(19i10)':
            chunk_size = 10
        else:
            raise ValueError(f"UnSupport Element Format: {fortran_format}")
        self.iter_line = f_handle.readline()

        if solid_type:
            # TODO: 根据单元类型判断是否为二维问题
            # TODO: 处理高阶单元(节点在第二行的情况)
            while self.iter_line.startswith(" ") and not self.iter_line.__contains__("-1"):
                # 解析单元数据
                e_data = [int(self.iter_line[i:i + chunk_size])
                          for i in range(0, min(len(self.iter_line), 19 * chunk_size), chunk_size)
                          if self.iter_line[i:i + chunk_size].strip()]
                mat_num = e_data[0]
                e_type = e_data[1]
                real_constant_num = e_data[2]
                sec_id = e_data[3]
                parsed_nodes_count = e_data[8]
                ele_num = e_data[10]

                # 处理节点ID
                node_ids = np.zeros(parsed_nodes_count, dtype=np.uint32)
                search_ids = np.zeros(parsed_nodes_count, dtype=np.uint32)
                for idx in range(parsed_nodes_count):
                    node_ids[idx] = e_data[idx + 11]
                    search_ids[idx] = self.femdb.node_hash[node_ids[idx]]

                # 创建单元并设置属性
                ele_node_list = list(OrderedDict.fromkeys(search_ids))
                self.node_search_ids_list.extend(ele_node_list)
                iter_ele, e_node_count, ele_matrix_size = ElementFactory.CreateElement(
                    e_type=self.et_hash[e_type],
                    opt=len(ele_node_list))

                iter_ele.SetNodeSearchIndex(np.asarray(ele_node_list))
                iter_ele.SetId(ele_num)
                iter_ele.SetNodes(np.asarray(list(OrderedDict.fromkeys(node_ids))))
                iter_ele.mat_id = mat_num
                iter_ele.sec_id = sec_id
                iter_ele.real_const_id = real_constant_num
                self.femdb.matrix_num_count += ele_matrix_size

                # 设置单元节点坐标
                coords = []
                for nid in ele_node_list:
                    n_coord = self.femdb.node_list[nid].GetNodeCoord()
                    coords.append(n_coord)
                iter_ele.SetNodeCoords(np.asarray(coords))

                # 添加到数据库
                self.femdb.ele_hash[ele_num] = len(self.femdb.elements)
                self.femdb.elements.append(iter_ele)

                self.ele_count += 1
                self.iter_line = f_handle.readline()
        else:
            """
            非SOLID格式的单元块格式说明:
            - 字段1: 单元号
            - 字段2: 截面类型ID
            - 字段3: 实常数号
            - 字段4: 材料号
            - 字段5: 单元坐标系号
            - 字段6-15: 节点号(超过10个节点时在下一行继续)
            """
            mlogger.fatal("UnSupport without SOLID keyword")
            sys.exit(1)

    def ReadMaterial(self, f_handle, key):
        """
        解析材料属性

        Args:
            f_handle: 文件句柄
            key (str): 材料关键字("MPDATA,"或"MP,")
        """
        jump_out = False
        while not jump_out:
            still_same_mat = True
            if self.iter_line.startswith("MPTEMP"):
                self.iter_line = f_handle.readline()
            if key == "MPDATA,":
                cur_mat_id = int(self.iter_line.split(",")[4].strip())
            elif key == "MP,":
                cur_mat_id = int(self.iter_line.split(",")[2].strip())
            else:
                raise KeyError(key)

            value_dict = {}
            while still_same_mat:
                # 跳过温度定义
                if self.iter_line.startswith("MPTEMP,"):
                    self.iter_line = f_handle.readline()
                elif self.iter_line.startswith("MPDATA,"):
                    splits = self.iter_line.split(",")
                    iter_mat_id = int(splits[4])
                    if iter_mat_id != cur_mat_id:
                        # 当前材料结束，处理下一材料或退出
                        value_dict[MaterialKey.G] = value_dict[MaterialKey.E] / 2 / (1 + value_dict[MaterialKey.Niu])
                        self.material_map[cur_mat_id] = value_dict
                        self.iter_line = f_handle.readline()
                        break
                    if splits[3].startswith("EX"):
                        value_dict[MaterialKey.E] = float(splits[6])
                    elif splits[3].startswith("DENS"):
                        value_dict[MaterialKey.Density] = float(splits[6])
                    elif splits[3].startswith("NUXY"):
                        value_dict[MaterialKey.Niu] = float(splits[6])
                    self.iter_line = f_handle.readline()
                elif self.iter_line.startswith("MP,"):
                    splits = self.iter_line.split(",")
                    iter_mat_id = int(splits[2])
                    if iter_mat_id != cur_mat_id:
                        value_dict[MaterialKey.G] = value_dict[MaterialKey.E] / 2 / (1 + value_dict[MaterialKey.Niu])
                        self.material_map[cur_mat_id] = value_dict
                        self.iter_line = f_handle.readline()
                        break
                    if splits[1].startswith("EX"):
                        value_dict[MaterialKey.E] = float(splits[3])
                    elif splits[1].startswith("DENS"):
                        value_dict[MaterialKey.Density] = float(splits[3])
                    elif splits[1].startswith("NUXY"):
                        value_dict[MaterialKey.Niu] = float(splits[3])
                    self.iter_line = f_handle.readline()
                else:
                    jump_out = not (self.iter_line.startswith("MPDATA,") or self.iter_line.startswith("MPTEMP"))
                    if jump_out:
                        value_dict[MaterialKey.G] = value_dict[MaterialKey.E] / 2 / (1 + value_dict[MaterialKey.Niu])
                        self.material_map[cur_mat_id] = value_dict
                        break

    def ReadSection(self, f_handle):
        """
        解析截面属性

        支持的截面类型:
        - BEAM: 梁截面(RECT, CSOLID)
        - SHELL: 壳截面
        """
        while True:
            splits = self.iter_line.strip().split(",")
            sec_num = int(splits[1])
            if splits[2] == "BEAM":
                beam_type = splits[3].lstrip()
                msec_data = f_handle.readline().strip().split(",")
                sec_data = [float(msec_data[idx]) for idx in range(1, len(msec_data)) if msec_data[idx]]
                f_handle.readline()  # section offset
                f_handle.readline()  # section control
                if beam_type == "RECT":
                    inertia_character = BeamCalculator.CalculateMomentOfInertiaOfArea(BeamSectionType.Rectangle, sec_data)
                    area_character = BeamCalculator.CalEffectiveShearArea(BeamSectionType.Rectangle, sec_data)
                    self.section_map[int(sec_num)] = {**inertia_character, **area_character}
                elif beam_type == "CSOLID":
                    inertia_character = BeamCalculator.CalculateMomentOfInertiaOfArea(BeamSectionType.CircleSolid, sec_data)
                    area_character = BeamCalculator.CalEffectiveShearArea(BeamSectionType.CircleSolid, sec_data)
                    self.section_map[int(sec_num)] = {**inertia_character, **area_character}
                else:
                    mlogger.fatal("UnSupport Beam Type:{}".format(beam_type))
                    sys.exit(1)
                self.iter_line = f_handle.readline()

            elif splits[2].startswith("SHELL"):
                f_handle.readline()  # secoffset
                self.iter_line = f_handle.readline()  # sec block
                if self.iter_line.startswith("SECBLOCK"):
                    self.iter_line = f_handle.readline()  # sec block data
                    splits = self.iter_line.split(",")
                    thickness = splits[0]
                    f_handle.readline()  # sec control
                    self.section_map[int(sec_num)] = {MaterialKey.Thickness: float(thickness)}
                elif self.iter_line.startswith("SECDATA"):
                    splits = self.iter_line.split(",")
                    thickness = splits[1]
                    f_handle.readline()  # sec control
                    self.section_map[int(sec_num)] = {MaterialKey.Thickness: float(thickness)}
                else:
                    raise KeyError(f"Failed to parse shell section:{splits[2]}")

            else:
                mlogger.fatal(f"UnSupport Section Type:{splits[2]}")
                sys.exit(1)

            # 跳过注释和空行
            while self.iter_line.startswith("!") or self.iter_line.startswith(" "):
                self.iter_line = f_handle.readline()
            if not self.iter_line.startswith("SECTYPE,"):
                break

    def CheckModel(self):
        """
        检查模型完整性

        检查内容:
        - 是否存在未连接到任何单元的节点
        """
        node_ids = set(nd.id for nd in self.femdb.node_list)
        ele_nodes = []
        for _, group in self.femdb.ele_grp_hash.items():
            eles = group.Elements()
            for ele in eles:
                ele_nodes.extend(ele.node_ids)
        ele_nodes_set = set(ele_nodes)

        unnessary_nds = node_ids - ele_nodes_set
        if len(unnessary_nds) != 0:
            mlogger.debug("Unnecessary Node Ids:{}".format(unnessary_nds))
        else:
            mlogger.debug("No Unnecessary Node")

        mlogger.debug("Finish Check Model")


if __name__ == "__main__":
    b_time = time.time()
    # path = "../numerical example/ANSYS/zhijiaRenumber.cdb"
    # path = "../NumericalCases/Projects/qizhongji/last/MQ1330.cdb"
    path = "../NumericalCases/Projects/qizhongji/last/MQ1330_remesh.cdb"
    cps = CDBParser(path, False)
    cps.ParseFileAndInitFEMDB()
    e_time = time.time()
    print(f"Elapsed time: {e_time - b_time:.2f} seconds")