#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import scipy.sparse

from femdb.FEMDataBase import *
from scipy.sparse.linalg import factorized
import pypardiso

"""
**关于稀疏矩阵:**
   1. 要有效地构造矩阵, 请使用dok_matrix或lil_matrix, lil_matrix类支持基本切片和花式索引, 其语法与NumPy Array类似; lil_matrix形式是基于row的
因此能够很高效的转为csr, 但是转为csc效率相对较低.
   2. 强烈建议不要直接使用NumPy函数运算稀疏矩阵如果你想将NumPy函数应用于这些矩阵，首先要检查SciPy是否有自己的给定稀疏矩阵类的实现, 或者首先将稀疏矩
阵转换为NumPy数组(使用类的toarray()方法).
   3. 要执行乘法或转置等操作，首先将矩阵转换为CSC或CSR格式，效率高CSR格式特别适用于快速矩阵矢量产品
   4. CSR，CSC和COO格式之间的所有转换都是线性复杂度.
   5. 对于已知是正定对称矩阵的情况下, 如何用scipy快速求解逆矩阵:
   https://stackoverflow.com/questions/40703042/more-efficient-way-to-invert-a-matrix-knowing-it-is-symmetric-and-positive-semi#:~:text=%3E%3E%3E%3E%20M%20%3D%20np.random.rand%20%2810%2C10%29%20%3E%3E%3E%3E%20M%20%3D,inv_M%20%3D%20np.triu%20%28inv_M%29%20%2B%20np.triu%20%28inv_M%2C%20k%3D1%29.T
"""


# TODO: 试试 https://github.com/scikit-sparse/scikit-sparse, 其是官方指定的python接口
# TODO: https://scikit-sparse.readthedocs.io/en/latest/overview.html#introduction
# https://www.zhihu.com/question/40769339
# https://github.com/BeanLiu1994/solver_speed_test
# https://pypi.org/project/pypardiso/
# https://github.com/haasad/PyPardisoProject
# https://blog.csdn.net/hu_yuhang/article/details/126674357
# https://blog.csdn.net/wzj_sxpi/article/details/116232656
# https://github.com/xmlyqing00/Cholmod-Scikit-Sparse-Windows
# Numpy-intel: https://pypi.org/project/intel-numpy/
# Numpy+Mkl
# cholesky 分解求解对称正定问题
# https://blog.csdn.net/weixin_38285131/article/details/81288338
# https://stackoverflow.com/questions/15573557/call-c-using-eigen-library-function-in-python
# https://pybind11.readthedocs.io/en/stable/advanced/cast/eigen.html
# https://kun-liu.com/2012/07/27/linear-sparse-solvers/
# https://petsc.org/release/


class Domain(object):
    """
    Domain class : Define the problem domain
    Only a single instance of Domain class can be created
    TODO 去掉一些没用过的函数
    """

    def __init__(self, check_model):
        """
          具体注意事项需要结合InpReader::ParseFile中的注释. 不同单元集合中的属性可能不一样, 有的集合也可能创建但没有被用到, 再计算单元
        本构阵的时候需要对ele_sets进行循环.
        """
        self.femdb = FEMDataBase()
        self.eq_count = None  # 总刚维度
        self.free_dof_count = None  # 自由度的个数, 将总刚矩阵分为Kaa, Kab, Kba, Kbb
        self.bound_dof_count = None  # 自由度被约束的个数
        self.Ub = []  # 约束指定位移
        self.Ua = None  # 未被约束的自由度
        self.Ra = None  # 未被约束的自由度上的力或力矩
        self.right_hand = None  # 右端项, 长度等于node_count * per_node_dof
        # 以下为刚度阵相关
        self.stiff_list = []
        self.mass_list = []
        self.eq_nums = []
        self.check_model = check_model

    def AssignElementCharacter(self):
        """
        不同类型的输入文件会调用不同的函数
        对于cdb文件的准备计算函数
        1. 计算各个截面属性, 赋值给对应单元
        2. 计算单元刚度阵
        3. 计算总刚矩阵维度, 初始化总刚(type: np.ndarray)
        """
        input_type = GlobalInfor[GlobalVariant.InputFileSuffix]
        if input_type == InputFileType.CDB:
            self.femdb.AssignElementProperty()
        elif input_type == InputFileType.INP:
            self.AssignElementCharacterAbaqus()
        elif input_type == InputFileType.BDF:
            pass

    def AssignElementCharacterAbaqus(self):
        """
        为计算做准备, 此时文件已经解析完成, 所有set和单元材料均已解析, 是准备开始计算的第一步
        1. 将Section属性赋予每个单元, 包括材料、几何尺寸等, 执行此操作要对section进行循环,
        2. 计算单元刚度矩阵
        3. 计算总刚矩阵维度, 初始化总刚(type: np.ndarray)
        """
        # 检查各个节点和单元Set是否有重名, 如果有重名, 那么是建模问题, 可以
        n_set_names = []
        e_set_names = []
        # for n_set in self.femdb.node_sets:
        #     n_set_names.append(n_set.GetName())
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

    def CalBoundaryEffect(self):
        """
        计算每个节点有多少自由度
        """
        six_dof_nodes = []
        two_dof_nodes = []
        all_ele_nodes = []
        for elem in self.femdb.elements:
            e_type = elem.e_type
            # 即使单元是退化单元也可以
            if ElementFactory.GetElementNodeDofCount(e_type) == 6:
                six_dof_nodes.extend(elem.GetNodeSearchIndex())
            elif ElementFactory.GetElementNodeDofCount(e_type) == 2:
                two_dof_nodes.extend(elem.GetNodeSearchIndex())

            if self.check_model:
                all_ele_nodes.extend(elem.node_ids)

        if self.check_model:
            if len(set(all_ele_nodes)) != len(self.femdb.node_list):
                ele_node_set = set(all_ele_nodes)
                all_node_set = set()
                for inode in self.femdb.node_list:
                    all_node_set.add(inode.id)
                diff_list = list(all_node_set - ele_node_set)
                for ii in diff_list:
                    print(f"Useless Node: {ii}")
        """
        去除重复节点, 在转换为list, set无法遍历
        """
        six_dof_nodes = list(set(six_dof_nodes))
        two_dof_nodes = list(set(two_dof_nodes))
        if len(six_dof_nodes) != 0:
            self.femdb.per_node_dof = 6
        elif len(two_dof_nodes) != 0:
            self.femdb.per_node_dof = 2
        else:
            self.femdb.per_node_dof = 3

        """
        将自由度为6(2)的节点的位移长度以及方程号长度改为6(2), 其他默认为3
        """
        for nd_idx in six_dof_nodes:
            self.femdb.node_list[nd_idx].ChangeDofCount(6)
        for nd_idx in two_dof_nodes:
            self.femdb.node_list[nd_idx].ChangeDofCount(2)

        # 指定约束位移, 现在的位移只支持关键字约束
        suffix = GlobalInfor[GlobalVariant.InputFileSuffix]
        if suffix == InputFileType.INP:
            for bd in self.femdb.load_case.GetBoundaries():
                node_set_name = bd.GetSetName()
                node_set = self.femdb.GetSpecificFEMObject(FEMObject.NodeSet, node_set_name)
                node_ids = node_set.GetNodeIds()
                bd_type = bd.GetBoundaryType()
                for nd in node_ids:
                    nd_idx = self.femdb.node_hash[nd]
                    self.femdb.node_list[nd_idx].SetBoundaryWithINPType(b_type=bd_type)

        elif suffix == InputFileType.CDB:
            for bd in self.femdb.load_case.GetBoundaries():
                node_ids, directs, con_values = bd.GetConstrainInfor()
                d_count = len(node_ids)
                for idx in range(d_count):
                    nd_idx = self.femdb.node_hash[node_ids[idx]]
                    direct = directs[idx]
                    con_value = con_values[idx]
                    self.femdb.node_list[nd_idx].SetBoundaryWithCDBType(direct, con_value)

        else:
            mlogger.fatal("UnSupport Boundary")
            sys.exit(1)

    def CalculateEquationNumber(self):
        """
        计算各自由度在总刚中的位置, 即方程号. K.J Bathe P178(上册), 为了组装成Kaa, 需要两次对自由度的排序. 将含有未知位移向量的
        自由度放上面, 将包含外界强制位移的放在后面
        """
        # 第一遍排序, 将未约束的自由度排在前面
        # TODO: 当cdb文件中包含未被引用的节点, 那么整体刚度矩阵会奇异, 所以流程应该是将单元集合然后去重复得到唯一的节点, 但是对于RBE3这种自定义的节点也要注意
        self.eq_count = 0
        for node in self.femdb.node_list:
            # 如果是辅助节点, 那么他的自由度不予考虑, 开始计算下一个节点
            if node.is_assist_node:
                continue

            # 如果不是边界节点, 是内部节点, 那么节点的所有自由度都是一个未知变量
            if not node.is_boundary_node:
                dof_count = node.SetAllDofEqNum(self.eq_count)
                self.eq_count += dof_count

            # 对于边界节点, 视约束自由度的个数添加未知数
            else:
                for ii in range(node.GetDofCount()):
                    if not node.b_code[ii]:
                        node.SetEquationNumber(ii, self.eq_count)
                        self.eq_count += 1

        self.free_dof_count = self.eq_count
        self.femdb.equation_number = self.eq_count

        """
        第二遍排序, 排序约束自由度的方程号, 顺便初始化Ub
        """
        for node in self.femdb.node_list:
            if node.is_boundary_node:
                for ii in range(node.GetDofCount()):
                    if node.b_code[ii]:
                        node.SetEquationNumber(ii, self.eq_count)
                        self.Ub.append(node.dof_disp[ii])
                        self.eq_count += 1

        """
        自由度号已经计算完成, 对自由度相关变量进行初始化
        """
        self.bound_dof_count = self.eq_count - self.free_dof_count

        """
        排序完成现在每个节点已经计算好了对应的方程号, 现在设置每个单元所包含的节点的方程号, 以用来组装总刚
        计算每个单元节点包含节点对应的方程号, 是一个列表, 单刚是一个len(eq_numbers) * len(eq_numbers)的矩阵, 组装的
        时候K的第ij项加到对应总刚的eq_numbers[i]行 eq_numbers[j]列
        """
        for iter_ele in self.femdb.elements:
            # 不能对所有的search_node_ids进行循环, 因为这其中包括了辅助节点
            eq_numbers = np.asarray([], dtype=np.uint32)
            search_node_ids = iter_ele.search_node_ids

            for idx in range(iter_ele.nodes_count):
                nid = search_node_ids[idx]
                node = self.femdb.GetNodeBySearchId(nid)
                eq_numbers = np.concatenate((eq_numbers, node.eq_num))
            iter_ele.SetEquationNumber(eq_numbers)

        """
        初始化总刚
        """
        self.femdb.InitAssemblyMatrix(self.eq_count)

    def CalAllElementStiffness(self):
        """
        调用多线程同时计算刚度矩阵
        计算所有单元的刚度阵, 对所有的单元组进行循环
        方便的查看各步骤运行时间: '%Y-%m-%d %H:%M:%S.%f')[:-3]
        """
        for iter_ele in self.femdb.elements:
            self.eq_nums.append(iter_ele.GetElementEquationNumber())
            iter_ele.CalculateBasic()
            stiff = iter_ele.ElementStiffness()
            if self.check_model:
                if iter_ele.id == 786:
                    print("")
                has_zero_row = (stiff >= 1e-8).all(axis=1).any()
                if has_zero_row:
                    # raise ValueError(f"Element {iter_ele.id} has zero row")
                    print(f"Element {iter_ele.id} has zero row")
            self.stiff_list.append(stiff)

    def AssembleStiffnessMatrixByElimination(self):
        """
        Assemble the banded global stiffness matrix, STAPPy中的del Matrix是否会减少内存分配, 或提高运算速度
        之前是lil_matrix, 但是速度很慢, 大概是现在方法的4倍左右, 原因是如下行程序所示, 需要__getitem__然后__setitem__
        self.femdb.global_stiff_matrix[eq_nums[row], eq_nums[column]] += stiff_mat[row, column]
        Reference:
        1. https://stackoverflow.com/questions/59460230/instantiate-large-sparse-matrices-for-assignment-operation
        2. https://stackoverflow.com/questions/27770906/why-are-lil-matrix-and-dok-matrix-so-slow-compared-to-common-dict-of-dicts
        """
        rows = []
        cols = []
        datas = []
        for ii in range(len(self.eq_nums)):
            eq_nums = self.eq_nums[ii]
            stiff_mat = self.stiff_list[ii]
            for row in range(len(eq_nums)):
                for column in range(len(eq_nums)):
                    if stiff_mat[row, column] != 0:
                        rows.append(eq_nums[row])
                        cols.append(eq_nums[column])
                        datas.append(stiff_mat[row, column])
        self.femdb.global_stiff_matrix = sparse.coo_matrix((datas, (rows, cols)), shape=(self.eq_count, self.eq_count))

        if GlobalInfor[GlobalVariant.PlotGlobalStiffness]:
            plt.spy(self.femdb.global_stiff_matrix, markersize=1)
            plt.title("GlobalStiffnessMatrix")
            plt.savefig("GlobalStiffness.jpg")

    def AssembleStiffnessMatrixByPenalty(self):
        """
        与AssembleStiffnessMatrixByElimination不同的是，前者是将约束的自由度去掉, 而本方法是通过罚函数或者乘大数法来实现
        本方法限制所有节点的自由度是同一个数, 即不可以让有的节点自由度是3, 有的节点自由度是6
        @return:
        """
        rows = []
        cols = []
        datas = []
        GlobalNodeHash = self.femdb.node_hash
        per_node_dof = self.femdb.per_node_dof
        for kk, ele in enumerate(self.femdb.elements):
            ele_nodes = ele.node_ids
            stiff_array = self.stiff_list[kk]
            for ii, nodeA in enumerate(ele_nodes):
                equA = GlobalNodeHash[nodeA] * per_node_dof
                for jj, nodeB in enumerate(ele_nodes):
                    equB = GlobalNodeHash[nodeB] * per_node_dof
                    for m in range(per_node_dof):
                        for n in range(per_node_dof):
                            TolRow = equA + m
                            TolCol = equB + n
                            eRow = ii * per_node_dof + m
                            eClo = jj * per_node_dof + n
                            rows.append(TolRow)
                            cols.append(TolCol)
                            datas.append(stiff_array[eRow, eClo])
        matrix_size = len(self.femdb.node_list) * per_node_dof
        self.femdb.global_stiff_matrix = sparse.coo_matrix((datas, (rows, cols)),
                                                           shape=(matrix_size, matrix_size))

    def CalAllElementMassMatrix(self):
        """
        计算所有单元的质量矩阵
        :return:
        """
        for iter_ele in self.femdb.elements:
            self.mass_list.append(iter_ele.ElementMass())

    def AssembleMassMatrixByPerturbation(self):
        """
        计算全局质量阵, 对角线增加摄动项
        :return:
        """
        rows = []
        cols = []
        datas = []
        GlobalNodeHash = self.femdb.node_hash
        per_node_dof = self.femdb.per_node_dof
        for kk, ele in enumerate(self.femdb.elements):
            nodes = ele.node_ids
            mass_array = self.mass_list[kk]
            for ii, nodeA in enumerate(nodes):
                equA = GlobalNodeHash[nodeA] * per_node_dof
                for jj, nodeB in enumerate(nodes):
                    equB = GlobalNodeHash[nodeB] * per_node_dof
                    for m in range(per_node_dof):
                        for n in range(per_node_dof):
                            TolRow = equA + m
                            TolCol = equB + n
                            eRow = ii * per_node_dof + m
                            eCol = jj * per_node_dof + n
                            if mass_array[eRow, eCol] != 0:
                                rows.append(TolRow)
                                cols.append(TolCol)
                                datas.append(mass_array[eRow, eCol])
        matrix_size = len(self.femdb.node_list) * per_node_dof

        """
        给质量阵添加扰动项，避免计算奇异
        """
        for ii in range(matrix_size):
            rows.append(ii)
            cols.append(ii)
            datas.append(1e-12)
        self.femdb.global_mass_matrix = sparse.coo_matrix((datas, (rows, cols)),
                                                          shape=(matrix_size, matrix_size)).tocsc()

    def SolveDisplacement(self):
        """
        将边界条件添加至节点, 有很多方法施加, 各自对不同的情况有利, 见Reference 1
        用本质边界条件修正刚度阵, 从而先求出未知自由度位移Ua, 结合已经给定自由度的位移即求解了所有自由度的位移
        约束分为显示约束和隐式约束
        Maa * D2Ua/Dt2 + Kaa * Ua = Ra - Kab * Ub - Mab * D2Ub/Dt2,
        对于静力问题, 方程简化为下式:
        Kaa * Ua = Ra - Kab * Ub
        历史解法(都因为过慢而淘汰):
        # method 1
        # self.Ua = spsolve(Kaa, self.Ra - Kab * self.Ub)
        # method 2
        # B = splu(Kaa)
        # self.Ua = Kaa.dot(B.solve(self.Ra - Kab*self.Ub))

        Reference:
        1. 《有限元分析的概念与应用》-第四版 (Robert D.Cook) P36 P421
        2. 《有限元法 理论、格式与求解方法》 (Bathe) P138 P178
        """
        # 组装Ra, 自然边界条件都是在右端项的Ua上, 不可以施加在Ub上, 现在只能处理集中载荷
        self.Ra = np.zeros((self.free_dof_count,), dtype=float)
        suffix = GlobalInfor[GlobalVariant.InputFileSuffix]
        if suffix == InputFileType.INP:
            """
            Abaqus格式的集中力是一个集合一个集合添加的
            """
            for c_load in self.femdb.load_case.GetConcentratedLoads():
                node_set = c_load.set_name
                f_dir = c_load.direction
                f_value = c_load.value
                # nodes = self.femdb.GetSpecificFEMObject(FEMObject.NodeSet, node_set).GetNodeIds()
                nodes = self.femdb.node_sets[node_set]
                for nd in nodes:
                    cnode = self.femdb.node_list[self.femdb.node_hash[nd]]
                    f_eq_num = cnode.GetEquationNumbers()[f_dir]
                    self.Ra[f_eq_num] = f_value

        elif suffix == InputFileType.CDB:
            """
            Ansys格式的集中力是一个自由度一个自由度添加的
            """
            for c_load in self.femdb.load_case.GetConcentratedLoads():
                cnode = self.femdb.node_list[self.femdb.node_hash[c_load.node]]
                f_eq_num = cnode.GetEquationNumbers()[c_load.direction]
                self.Ra[f_eq_num] = c_load.value

        else:
            mlogger.fatal("UnSupport File Type In SolveDisplacement")
            sys.exit(1)

        """
        施加位移约束, 见visio文档, 求解Boundary矩阵和V矩阵, 暂未实现隐式约束, 对于显示不用求解Boundary矩阵,
        因为这种情况罚函数影响的只是Kbb的对角元素, 与未知位移求解没关系. 求解支反力也不需要加入罚函数的值
        """
        self.femdb.global_stiff_matrix = self.femdb.global_stiff_matrix.tocsc()
        Kaa = self.femdb.global_stiff_matrix[:self.free_dof_count, :self.free_dof_count]
        Kab = self.femdb.global_stiff_matrix[:self.free_dof_count, self.free_dof_count:]
        self.Ub = np.asarray(self.Ub, dtype=float)

        # 求解self.Ua TODO: 没有利用Kaa是正定对称矩阵的性质, 另外Assemble对应的稀疏矩阵优化, 考虑用其他库的稀疏矩阵, 还有就是单刚的计算了
        self.Ua = pypardiso.spsolve(Kaa, self.Ra - Kab * self.Ub)

        # 所有自由度的解
        U = np.append(self.Ua, self.Ub)

        # 按自由度分配给所有节点
        for nd in self.femdb.node_list:
            for ii in range(nd.GetDofCount()):
                nd.dof_disp[ii] = U[nd.eq_num[ii]]
            nd.CalNodeMagnitudeDisplacement()

    def AddBoundaryByPenalty(self):
        """
        罚函数的方法施加约束
        @return:
        """
        # 指定约束位移, 现在的位移只支持关键字约束
        node_dof_count = self.femdb.per_node_dof
        # node_count = len(self.femdb.node_list)
        # self.right_hand = np.zeros(node_count * node_dof_count)
        for bd in self.femdb.load_case.GetBoundaries():
            node_set_name = bd.GetSetName()
            node_set = self.femdb.GetSpecificFEMObject(FEMObject.NodeSet, node_set_name)
            node_ids = node_set.GetNodeIds()
            bd_type = bd.GetBoundaryType()
            for nd in node_ids:
                nd_idx = self.femdb.node_hash[nd]
                if bd_type in ["ENCASTRE"]:
                    base_idx = nd_idx * node_dof_count
                    for ii in node_dof_count:
                        self.femdb.global_stiff_matrix[base_idx + ii, base_idx + ii] += 2e21
                        self.right_hand[base_idx + ii] = 0
                else:
                    raise KeyError(f"don't support boundary type: {bd_type}")

    def NewMarkExplict(self):
        """
        显式的动力学求解，使用NewMark方法
        Reference:
           《结构动力学基础》 张亚辉、林家浩 P94
        :return:
        """
        """
        初始化载荷, 计算初始加速度
        """
        his_load = self.femdb.load_case.history_loads
        force_nodes = [ii[0] for ii in his_load]
        directories = [ii[1] for ii in his_load]
        scale = [ii[2] for ii in his_load]
        steps_count = his_load[0][-1].shape[1]
        delta_t = his_load[0][-1][0, 1] - his_load[0][-1][0, 0]
        node_dof_count = self.femdb.per_node_dof
        node_count = len(self.femdb.node_list)
        fem_all_dofs = node_dof_count * node_count
        his_vals = np.zeros((fem_all_dofs, steps_count), dtype=float)

        for ii, nd in enumerate(force_nodes):
            eqa = self.femdb.node_hash[nd] * node_dof_count
            xyz = directories[ii]
            his_vals[eqa + xyz] = his_load[ii][-1][1, :] * scale[ii]

        acc_0 = pypardiso.spsolve(self.femdb.global_mass_matrix, his_vals[:, 0])

        """
        计算初始参数
        """
        alpha = 0.25
        delta = 0.5
        a0 = 1 / alpha / delta_t ** 2
        a1 = delta / alpha / delta_t
        a2 = 1 / alpha / delta_t
        a3 = 0.5 / alpha - 1
        a4 = delta / alpha - 1
        a5 = delta_t / 2 * (delta / alpha - 2)
        a6 = delta_t * (1 - delta)
        a7 = delta * delta_t
        K_hat = self.femdb.global_stiff_matrix + a0 * self.femdb.global_mass_matrix
        solve = factorized(K_hat.tocsc())

        """
        开始计算各个时间步的值
        """
        u = [np.zeros(fem_all_dofs, dtype=float)]
        v = [np.zeros(fem_all_dofs, dtype=float)]
        a = [acc_0]
        for ii in range(steps_count - 1):
            f_hat = his_vals[:, ii + 1] + self.femdb.global_mass_matrix @ (a0 * u[ii] + a2 * v[ii] + a3 * a[ii])
            u_next = solve(f_hat)
            a_next = a0 * (u_next - u[ii]) - a2 * v[ii] - a3 * a[ii]
            v_next = v[ii] + a6 * a[ii] + a7 * a_next
            u.append(u_next)
            v.append(v_next)
            a.append(a_next)

        self.femdb.history_u = u
        self.femdb.history_v = v
        self.femdb.history_a = a
        self.femdb.history_step_count = len(u)

    def CalculateNodeStress(self):
        """
        计算节点的应力
        """
        for _, ele_group in self.femdb.ele_grp_hash.items():
            eles = ele_group.Elements()
            for iter_ele in eles:
                # 计算各个单元的节点位移
                search_node_ids = iter_ele.search_node_ids
                ele_displacement = []
                for idx in range(iter_ele.nodes_count):
                    nid = search_node_ids[idx]
                    node = self.femdb.GetNodeBySearchId(nid)
                    ele_displacement.append(node.dof_disp)

                # 根据节点位移, 以单元为计算单位, 计算单元节点应力
                # stress = iter_ele.ElementStress(np.asarray(ele_displacement))
                # for idx in range(iter_ele.nodes_count):
                #     nid = search_node_ids[idx]
                #     node = self.femdb.GetNodeBySearchId(nid)
                #     node.AppendStressResult(stress[idx,:])

        # for node in self.femdb.node_list:
        #     node.AverageStress()

    """ 
    以下的函数为计算输出部分, 后处理部分
    """

    def GetDisplacementBySearchId(self, nd_id):
        return self.femdb.node_list[nd_id].displacement


if __name__ == "__main__":
    from scipy.linalg import cho_factor, cho_solve

    K = np.array([[4, 1], [1, 3]])
    b_list = [np.array([1, 2]), np.array([3, 4])]
    c, lower = cho_factor(K)
    x_list = [cho_solve((c, lower), b) for b in b_list]

    # for i, x in enumerate(x_list):
    #     print(f"x[{i}] = {x}")

    print(K @ np.array([0.09090909, 0.63636364]))
    print(K @ np.array([0.45454545, 1.18181818]))

    K1 = sparse.csc_matrix(K)
    solve = factorized(K1)
    x_list2 = [solve(b) for b in b_list]
    for i, x in enumerate(x_list2):
        print(f"x[{i}] = {x}")
