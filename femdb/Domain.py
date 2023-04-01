#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import meshio
from femdb.FEMDataBase import *


class Domain(object):
    """
    Domain class : Define the problem domain
    Only a single instance of Domain class can be created
    Reference:
    1. https://github.com/nschloe/meshio  # write vtp
    TODO 去掉一些没用过的函数
    TODO 提高程序效率, 各个类的组织以及函数应该属于哪个类, 在哪里计算
    TODO 优化import
    """

    def __init__(self):
        """
          具体注意事项需要结合InpReader::ParseFile中的注释. 不同单元集合中的属性可能不一样, 有的集合也可能创建但没有被用到, 再计算单元
        本构阵的时候需要对ele_sets进行循环.
        """
        self.femdb = FEMDataBase()
        # mlogger.debug("In Domain, femdb id is:{}".format(id(self.femdb)))
        self.eq_count = None  # 总刚维度
        self.free_dof_idx = None  # 自由度的个数, 将总刚矩阵分为Kaa, Kab, Kba, Kbb
        self.bound_dof_count = None  # 自由度被约束的个数
        self.Ub = []  # 约束指定位移
        self.Ua = None  # 未被约束的自由度
        self.Ra = None  # 未被约束的自由度上的力或力矩

    def PrintFEMDBSummary(self):
        """
        计算前可以首先打印出数据库信息
        """
        self.femdb.PrintParseSummary()

    def AssignElementCharacter(self):
        """
        不同类型的输入文件会调用不同的函数
        """
        if self.femdb.input_file_type == InputFileType.CDB:
            self.AssignElementCharacterAnsys()
        elif self.femdb.input_file_type == InputFileType.INP:
            self.AssignElementCharacterAbaqus()
        elif self.femdb.input_file_type == InputFileType.BDF:
            self.AssignElementCharacterNastran()

    def AssignElementCharacterAnsys(self):
        """
        对于cdb文件的准备计算函数
        1. 计算各个截面属性, 赋值给对应单元
        2. 计算单元刚度阵
        3. 计算总刚矩阵维度, 初始化总刚(type: np.mat)
        """
        self.femdb.AssignElementProperty()

    def AssignElementCharacterNastran(self):
        """
        对于bdf文件的准备计算函数
        """

    def AssignElementCharacterAbaqus(self):
        """
        为计算做准备, 此时文件已经解析完成, 所有set和单元材料均已解析, 是准备开始计算的第一步
        1. 将Section属性赋予每个单元, 包括材料、几何尺寸等, 执行此操作要对section进行循环,
        2. 计算单元刚度矩阵
        3. 计算总刚矩阵维度, 初始化总刚(type: np.mat)
        """
        # 检查各个节点和单元Set是否有重名, 如果有重名, 那么是建模问题, 可以
        n_set_names = []
        e_set_names = []
        for n_set in self.femdb.node_sets:
            n_set_names.append(n_set.GetName())
        for e_set in self.femdb.ele_sets:
            e_set_names.append(e_set.GetName())
        if len(set(n_set_names)) != len(n_set_names):
            mlogger.fatal(r"包含重名的节点集合")
            sys.exit(1)
        if len(set(e_set_names)) != len(e_set_names):
            mlogger.fatal(r"包含重名的单元集合")
            sys.exit(1)

        # 分配属性参数, 并计算单元刚度阵
        self.femdb.AssignElementProperty()

    def CallBoundaryEffect(self):
        # 计算每个节点有多少自由度
        six_dof_nodes = []
        two_dof_nodes = []
        for key, ele_grp in self.femdb.ele_grp_hash.items():
            if ElementFactory.GetElementNodeDofCount(key) == 6:
                for elem in ele_grp.Elements():
                    six_dof_nodes.extend(elem.GetNodeSearchId())
            elif ElementFactory.GetElementNodeDofCount(key) == 2:
                for elem in ele_grp.Elements():
                    two_dof_nodes.extend(elem.GetNodeSearchIndex())

        # 去除重复节点, 在转换为list, set无法遍历
        six_dof_nodes = list(set(six_dof_nodes))
        two_dof_nodes = list(set(two_dof_nodes))

        # 将自由度为6(2)的节点的位移长度以及方程号长度改为6(2), 其他默认为3
        for nd_idx in six_dof_nodes:
            self.femdb.node_list[nd_idx].ChangeDofCount(6)
        for nd_idx in two_dof_nodes:
            self.femdb.node_list[nd_idx].ChangeDofCount(2)

        # 指定约束位移, 现在的位移只支持关键字约束
        for bd in self.femdb.load_case.GetBoundaries():
            node_set_name = bd.GetSetName()
            node_set = self.femdb.GetSpecificFEMObject(FEMObject.NodeSet, node_set_name)
            node_ids = node_set.GetNodeIds()
            bd_type = bd.GetBoundaryType()
            for nd in node_ids:
                nd_idx = self.femdb.node_hash[nd]
                self.femdb.node_list[nd_idx].SetBoundary(b_type=bd_type)


    def CalculateEquationNumber(self):
        """
        计算各自由度在总刚中的位置, 即方程号. K.J Bathe P178(上册), 为了组装成Kaa, 需要两次对自由度的排序. 将含有未知位移向量的
        自由度放上面, 将包含外界强制位移的放在后面
        """
        # 第一遍排序, 将未约束的自由度排在前面
        self.eq_count = 0
        for node in self.femdb.node_list:
            if not node.is_boundary_node:
                dof_count = node.SetAllDofEqNum(self.eq_count)
                self.eq_count += dof_count
            else:
                for ii in range(node.GetDofCount()):
                    if not node.b_code[ii]:
                        node.SetEquationNumber(ii, self.eq_count)
                        self.eq_count += 1

        self.free_dof_idx = self.eq_count

        # 第二遍排序, 排序约束自由度的方程号, 顺便初始化Ub
        for node in self.femdb.node_list:
            if node.is_boundary_node:
                for ii in range(node.GetDofCount()):
                    if node.b_code[ii]:
                        node.SetEquationNumber(ii, self.eq_count)
                        self.Ub.append(node.dof_disp[ii])
                        self.eq_count += 1

        # 自由度号已经计算完成, 对自由度相关变量进行初始化
        self.bound_dof_count = self.eq_count - self.free_dof_idx

        # 排序完成现在每个节点已经计算好了对应的方程号, 现在设置每个单元所包含的节点的方程号,
        # 计算每个单元节点包含节点对应的方程号, 是一个列表, 单刚是一个len(eq_numbers) * len(eq_numbers)的矩阵, 组装的
        # 时候K的第ij项加到对应总刚的eq_numbers[i]行 eq_numbers[j]列
        for _, ele_group in self.femdb.ele_grp_hash.items():
            eles = ele_group.Elements()
            for iter_ele in eles:
                eq_numbers = np.asarray([], dtype=np.uint32)
                for nid in iter_ele.search_node_ids:
                    node = self.femdb.GetNodeBySearchId(nid)
                    eq_numbers = np.concatenate((eq_numbers, node.eq_num))
                iter_ele.SetEquationNumber(eq_numbers)

        # 初始化总刚
        self.femdb.InitAssemblyMatrix(self.eq_count)

    def AssembleStiffnessMatrix(self):
        """
        Assemble the banded global stiffness matrix, STAPPy中的del Matrix是否会减少内存分配, 或提高运算速度
        """
        # Loop over for all element groups
        for key, ele_group in self.femdb.ele_grp_hash.items():
            for iter_ele in ele_group.eles:
                eq_nums = iter_ele.GetElementEquationNumber()
                stiff_mat = iter_ele.ElementStiffness()
                for row in range(len(eq_nums)):
                    for column in range(len(eq_nums)):
                        self.femdb.total_stiff_matrix[eq_nums[row]][eq_nums[column]] += stiff_mat[row, column]
        mlogger.debug("Finish Assemble Stiffness Matrix")

    def SolveDisplacement(self):
        """
        将边界条件添加至节点, 有很多方法施加, 各自对不同的情况有利, 见Reference 1
        用本质边界条件修正刚度阵, 从而先求出未知自由度位移Ua, 结合已经给定自由度的位移即求解了所有自由度的位移
        约束分为显示约束和隐式约束
        Maa * D2Ua/Dt2 + Kaa * Ua = Ra - Kab * Ub - Mab * D2Ub/Dt2,
        对于静力问题, 方程简化为下式:
        Kaa * Ua = Ra - Kab * Ub

        Reference:
        1. 《有限元分析的概念与应用》-第四版 (Robert D.Cook) P36 P421
        2. 《有限元法 理论、格式与求解方法》 (Bathe) P138 P178
        """

        # 组装Ra, 自然边界条件都是在右端项的Ua上, 不可以施加在Ub上, 现在只能处理集中载荷
        self.Ra = np.asarray([0] * self.free_dof_idx).T
        for c_load in self.femdb.load_case.GetConcentratedLoad():
            node_set = c_load.set_name
            f_dir = c_load.direction
            f_value = c_load.value
            nodes = self.femdb.GetSpecificFEMObject(FEMObject.NodeSet, node_set).GetNodeIds()
            for nd in nodes:
                cnode = self.femdb.node_list[self.femdb.node_hash[nd]]
                f_eq_num = cnode.GetEquationNumbers()[f_dir]
                self.Ra[f_eq_num] = f_value

        # 施加位移约束, 见visio文档, 求解Boundary矩阵和V矩阵, 暂未实现隐式约束, 对于显示不用求解Boundary矩阵,
        # 因为这种情况罚函数影响的只是Kbb的对角元素, 与未知位移求解没关系. 求解支反力也不需要加入罚函数的值

        # 求解
        Kaa = self.femdb.total_stiff_matrix[:self.free_dof_idx, :self.free_dof_idx]
        Kab = self.femdb.total_stiff_matrix[:self.free_dof_idx, self.free_dof_idx:]
        self.Ub = np.asarray(self.Ub).T
        self.Ua = np.matmul(np.linalg.inv(Kaa), self.Ra - np.matmul(Kab, self.Ub))
        U = np.append(self.Ua, self.Ub)

        # 按自由度分配给所有节点
        for nd in self.femdb.node_list:
            for ii in range(nd.GetDofCount()):
                nd.dof_disp[ii] = U[nd.eq_num[ii]]
            nd.CalNodeSpaceDisplacement()

    """ 
    以下的函数为计算输出部分, 后处理部分
    """

    def GetDisplacementBySearchId(self, nd_id):
        return self.femdb.node_list[nd_id].displacement

    def WriteOutPutFile(self, path):
        """
        将结果写至vtp文件
        """
        # 模型部分
        coords = np.asarray([node.coord for node in self.femdb.node_list])
        all_eles = {}
        for key, ele_grp in self.femdb.ele_grp_hash.items():
            eles = ele_grp.Elements()
            ele_count = len(eles)
            ele_node_count = len(eles[0].GetNodes())
            relations = np.zeros([ele_count, ele_node_count], dtype=int)
            for i in range(ele_count):
                relations[i] = eles[i].GetNodeSearchIndex()

            if self.femdb.input_file_type == InputFileType.CDB:
                key = self.femdb.et_hash[key]
                all_eles[Ansys2VTKType[key]] = relations
            elif self.femdb.inpu_file_type == InputFileType.INP:
                all_eles[Abaqus2VTKType[key]] = relations

        # 位移结果
        dis_value = np.asarray([node.displacement for node in self.femdb.node_list])
        displacement = {"displacement": dis_value}
        # mlogger.debug(dis_value)
        meshio.write_points_cells(
            filename=path,
            points=coords,
            cells=all_eles,
            # point_data=displacement,
            # cell_data=cell_data,
            # field_data=field_data
        )
        mlogger.debug("Finish Write Output File")
