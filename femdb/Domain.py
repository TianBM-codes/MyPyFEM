#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
from copy import deepcopy

import numpy as np
from femdb.FEMDataBase import *
from femdb.GlobalFEMVariant import ModelInfo
from scipy.sparse.linalg import factorized, spsolve
from pathlib import Path
from femdb.GlobalEnum import *
import pypardiso

"""
**关于稀疏矩阵:**
   1. 要有效地构造矩阵, 请使用dok_matrix或lil_matrix, lil_matrix类支持基本切片和花式索引, 其语法与NumPy Array类似; lil_matrix形式是基于row的
因此能够很高效的转为csr, 但是转为csc效率相对较低.
   2. 强烈建议不要直接使用NumPy函数运算稀疏矩阵如果你想将NumPy函数应用于这些矩阵，首先要检查SciPy是否有自己的给定稀疏矩阵类的实现, 或者首先将稀疏矩
阵转换为NumPy数组(使用类的toarray()方法).
   3. 要执行乘法或转置等操作, 首先将矩阵转换为CSC或CSR格式, 效率高. CSR格式特别适用于快速矩阵矢量
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
        self.right_hand = None  # 右端项, 长度等于node_count * ModelInfo.PER_NODE_DOF
        # 以下为刚度阵相关
        self.mass_list = []
        self.check_model = check_model

    def CalAllElementStiffness(self):
        """
        调用多线程同时计算刚度矩阵
        计算所有单元的刚度阵, 对所有的单元组进行循环
        方便的查看各步骤运行时间: '%Y-%m-%d %H:%M:%S.%f')[:-3]
        """
        calculated_eles_count = 0
        for iter_ele in self.femdb.elements:
            iter_ele.CalculateBasic()
            stiff = iter_ele.ElementStiffness()
            calculated_eles_count += 1
            if self.check_model:
                if iter_ele.id == 786:
                    print("")
                has_zero_row = (stiff >= 1e-8).all(axis=1).any()
                if has_zero_row:
                    # raise ValueError(f"Element {iter_ele.id} has zero row")
                    print(f"Element {iter_ele.id} has zero row")
                # if calculated_eles_count % 2000 == 0:
                #     print(f"calculated ele's stiff count: {calculated_eles_count}")
            self.femdb.stiff_list.append(stiff)


    def CalAllElementThermalMatrixAndAssemble(self):
        """
        计算所有单元的稳态导热矩阵
        :return:
        """
        num_elems = len(self.femdb.elements)

        # 线性四面体：单元温度矩阵 4x4，有 16 个非零
        nnz_per_elem = 16
        total_nnz = num_elems * nnz_per_elem

        rows = np.zeros(total_nnz, dtype=np.uint32)
        cols = np.zeros(total_nnz, dtype=np.uint32)
        datas = np.zeros(total_nnz, dtype=np.float64)

        iter_loc = 0

        for kk, ele in enumerate(self.femdb.elements):
            search_node_ids = ele.search_node_ids
            Kt = ele.ElementThermalMatrix()
            try:
                ele_dofs = np.array(
                    [self.femdb.node_list[x].temp_eq_num for x in search_node_ids],
                    dtype=np.uint32
                )
            except AttributeError as e:
                print(e)
                raise TypeError(f"Node temp_eq_num not set correctly in Element: {ele.id}")

            rows_block, cols_block = np.meshgrid(ele_dofs, ele_dofs)
            entries_count = Kt.size

            rows[iter_loc:iter_loc + entries_count] = rows_block.ravel()
            cols[iter_loc:iter_loc + entries_count] = cols_block.ravel()
            datas[iter_loc:iter_loc + entries_count] = Kt.ravel()

            iter_loc += entries_count

        matrix_dimension = len(self.femdb.node_list)
        self.femdb.global_thermal_matrix = sparse.coo_matrix(
            (datas, (rows, cols)),
            shape=(matrix_dimension, matrix_dimension)
        ).tocsc()


    def CheckRedundantNodes(self):
        no_dup_nodes = np.zeros(len(self.femdb.node_list), dtype=np.uint32)
        for ii, iter_node in enumerate(self.femdb.node_list):
            no_dup_nodes[ii] = iter_node.id
        nodes_set = set(no_dup_nodes)

        element_nodes = []
        for iter_ele in self.femdb.elements:
            element_nodes.extend(iter_ele.node_ids.tolist())
        ele_nodes_set = set(element_nodes)

        redundant_set = nodes_set - ele_nodes_set
        print("Redundant Nodes:", list(redundant_set))

    def AssembleStiffnessMatrixByPenalty(self, from_origin=False):
        """
        与AssembleStiffnessMatrixByElimination不同的是，前者是将约束的自由度去掉, 而本方法是通过罚函数或者乘大数法来实现
        本方法限制所有节点的自由度是同一个数, 即不可以让有的节点自由度是3, 有的节点自由度是6

        Assemble the banded global stiffness matrix, STAPPy中的del Matrix是否会减少内存分配, 或提高运算速度
        之前是lil_matrix, 但是速度很慢, 大概是现在方法的4倍左右, 原因是如下行程序所示, 需要__getitem__然后__setitem__
        self.femdb.global_stiff_matrix[eq_nums[row], eq_nums[column]] += stiff_mat[row, column]
        Reference:
        1. https://stackoverflow.com/questions/59460230/instantiate-large-sparse-matrices-for-assignment-operation
        2. https://stackoverflow.com/questions/27770906/why-are-lil-matrix-and-dok-matrix-so-slow-compared-to-common-dict-of-dicts
        @return:
        """
        """
        组装总体刚度阵, 首先考虑系数矩阵一共有多少个元素, 对于约束方程来说, 施加一个自由度约束需要添加4个量
        """
        ce_count = len(self.femdb.equation_constrain_couple)
        ce_add_equations = 0
        for ii in range(ce_count):
            iter_ce = self.femdb.equation_constrain_idx[ii]
            ce_add_equations += len(iter_ce)

        iter_loc = 0
        rows = np.zeros(self.femdb.matrix_num_count + ce_add_equations * 4, dtype=np.uint32)
        cols = np.zeros(self.femdb.matrix_num_count + ce_add_equations * 4, dtype=np.uint32)
        datas = np.zeros(self.femdb.matrix_num_count + ce_add_equations * 4, dtype=np.float64)
        for kk, ele in enumerate(self.femdb.elements):
            search_node_ids = ele.search_node_ids
            if from_origin:
                stiff_array = ele.ReCalculateElementStiffness()
            else:
                stiff_array = self.femdb.stiff_list[kk]

            try:
                ele_dofs = np.array([np.arange(self.femdb.node_list[x].start_eq_num,
                                               self.femdb.node_list[x].start_eq_num + ele.node_dof_count)
                                     for x in search_node_ids])
            except TypeError as e:
                print(e)
                raise TypeError(f"Element Node Dof Count Wrong in Element:{ele.id}")

            rows_block, cols_block = np.meshgrid(ele_dofs, ele_dofs)
            entries_count = ele.block_size
            rows[iter_loc:iter_loc + entries_count] = rows_block.ravel()
            cols[iter_loc:iter_loc + entries_count] = cols_block.ravel()
            datas[iter_loc:iter_loc + entries_count] = stiff_array.ravel()

            iter_loc += entries_count

        """
        添加RBE2约束, 拉格朗日乘子法实现约束方程, 需要增加总刚维度
        """
        GlobalNodeHash = self.femdb.node_hash
        matrix_dimension = len(self.femdb.node_list) * ModelInfo.PER_NODE_DOF
        row_cursor = 0
        insert_offset = 0
        for ii in range(ce_count):
            m_node, s_node = self.femdb.equation_constrain_couple[ii]
            iter_ce_idx = self.femdb.equation_constrain_idx[ii]
            m_node_dof = GlobalNodeHash[m_node] * ModelInfo.PER_NODE_DOF
            s_node_dof = GlobalNodeHash[s_node] * ModelInfo.PER_NODE_DOF
            for jj in iter_ce_idx:
                insert_begin = self.femdb.matrix_num_count + insert_offset
                insert_end = insert_begin + 4
                rows[insert_begin: insert_end] = [matrix_dimension + row_cursor,
                                                  matrix_dimension + row_cursor,
                                                  s_node_dof + jj,
                                                  m_node_dof + jj]

                cols[insert_begin:insert_end] = [s_node_dof + jj,
                                                 m_node_dof + jj,
                                                 matrix_dimension + row_cursor,
                                                 matrix_dimension + row_cursor]
                datas[insert_begin: insert_end] = [1, -1, 1, -1]
                insert_offset += 4

                row_cursor += 1

        matrix_dimension += row_cursor

        """
        完成所有三向量数据的准备, 开始组装总刚
        """
        self.femdb.global_stiff_matrix = sparse.coo_matrix((datas, (rows, cols)),
                                                           shape=(matrix_dimension, matrix_dimension)).tocsc()

        if self.check_model:
            """
            0. 检查模型中是否存在不属于任何单元的节点
            """
            self.CheckRedundantNodes()

            """
            1. 输入矩阵或向量包含Nan或inf, 求解器无法处理非法数值
            """
            print("GlobalStiffMatrix contains NaN:", np.isnan(self.femdb.global_stiff_matrix.data).any())

            """
            2. 若A的行列式为零(或条件数极大), 求解器无法找到唯一解, 数值计算失败
            """
            cond = np.linalg.cond(self.femdb.global_stiff_matrix.todense())
            print("Condition number:", cond)
            """

            3. 检查矩阵是否对称正定
            """
            from scipy.linalg import eigh
            eigenvalues = eigh(self.femdb.global_stiff_matrix.todense())[0]
            print("Min eigenvalue:", np.min(eigenvalues))

        if GlobalInfor[GlobalVariant.PlotGlobalStiffness]:
            plt.spy(self.femdb.global_stiff_matrix, markersize=1)
            plt.title("GlobalStiffnessMatrix")
            plt.savefig("GlobalStiffness.jpg")

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
        for kk, ele in enumerate(self.femdb.elements):
            nodes = ele.node_ids
            mass_array = self.mass_list[kk]
            for ii, nodeA in enumerate(nodes):
                equA = GlobalNodeHash[nodeA] * ModelInfo.PER_NODE_DOF
                for jj, nodeB in enumerate(nodes):
                    equB = GlobalNodeHash[nodeB] * ModelInfo.PER_NODE_DOF
                    for m in range(ModelInfo.PER_NODE_DOF):
                        for n in range(ModelInfo.PER_NODE_DOF):
                            TolRow = equA + m
                            TolCol = equB + n
                            eRow = ii * ModelInfo.PER_NODE_DOF + m
                            eCol = jj * ModelInfo.PER_NODE_DOF + n
                            if mass_array[eRow, eCol] != 0:
                                rows.append(TolRow)
                                cols.append(TolCol)
                                datas.append(mass_array[eRow, eCol])
        matrix_size = len(self.femdb.node_list) * ModelInfo.PER_NODE_DOF

        """
        给质量阵添加扰动项，避免计算奇异
        """
        for ii in range(matrix_size):
            rows.append(ii)
            cols.append(ii)
            datas.append(1e-12)
        self.femdb.global_mass_matrix = sparse.coo_matrix((datas, (rows, cols)),
                                                          shape=(matrix_size, matrix_size)).tocsc()

    def CalculateGreenFunction(self, output_dir):
        """
        计算格林函数并保存
        :param output_dir: 保存路径
        :return:
        """
        c_load = self.femdb.load_case.c_loads
        force_nodes = [ii[0] for ii in c_load]
        directories = [ii[1] for ii in c_load]
        right_hand_back_up = deepcopy(self.right_hand)
        for ii, nd in enumerate(force_nodes):
            """
            保存位移结果
            """
            eqa = self.femdb.node_hash[nd] * ModelInfo.PER_NODE_DOF
            xyz = directories[ii]
            self.right_hand = deepcopy(right_hand_back_up)
            self.right_hand[eqa + xyz] = 1
            temp_results = pypardiso.spsolve(self.femdb.global_stiff_matrix, self.right_hand)
            result_range = len(self.femdb.node_list) * ModelInfo.PER_NODE_DOF
            iter_linear_u = temp_results[:result_range]
            np.save(output_dir / f"dis_node_{nd}.npy", iter_linear_u)
            """
            求解应力结果并保存
            """
            iter_sigma_xx = np.zeros(len(self.femdb.node_list))
            iter_sigma_yy = np.zeros(len(self.femdb.node_list))
            iter_sigma_zz = np.zeros(len(self.femdb.node_list))
            iter_tau_xy = np.zeros(len(self.femdb.node_list))
            iter_tau_xz = np.zeros(len(self.femdb.node_list))
            iter_tau_yz = np.zeros(len(self.femdb.node_list))
            for ele in self.femdb.elements:
                search_idx = []
                for ii in ele.search_node_ids:
                    start = ii * ModelInfo.PER_NODE_DOF
                    end = (ii + 1) * ModelInfo.PER_NODE_DOF
                    search_idx.extend(np.arange(start, end, 1).tolist())

                iter_u = iter_linear_u[search_idx].flatten()
                sigma_x, sigma_y, sigma_z, tau_yz, tau_xz, tau_xy = ele.CalculateElementStress(iter_u)

                for ii, n_search_id in enumerate(ele.search_node_ids):
                    iter_sigma_xx[n_search_id] += sigma_x[ii] / self.femdb.node_connected_element_count[n_search_id]
                    iter_sigma_yy[n_search_id] += sigma_y[ii] / self.femdb.node_connected_element_count[n_search_id]
                    iter_sigma_zz[n_search_id] += sigma_z[ii] / self.femdb.node_connected_element_count[n_search_id]
                    iter_tau_xy[n_search_id] += tau_xy[ii] / self.femdb.node_connected_element_count[n_search_id]
                    iter_tau_yz[n_search_id] += tau_yz[ii] / self.femdb.node_connected_element_count[n_search_id]
                    iter_tau_xz[n_search_id] += tau_xz[ii] / self.femdb.node_connected_element_count[n_search_id]

            np.save(output_dir / f"sigma_xx_node_{nd}.npy", iter_sigma_xx)
            np.save(output_dir / f"sigma_yy_node_{nd}.npy", iter_sigma_yy)
            np.save(output_dir / f"sigma_zz_node_{nd}.npy", iter_sigma_zz)
            np.save(output_dir / f"sigma_xy_node_{nd}.npy", iter_tau_xy)
            np.save(output_dir / f"sigma_yz_node_{nd}.npy", iter_tau_yz)
            np.save(output_dir / f"sigma_xz_node_{nd}.npy", iter_tau_xz)

    def CalculateResultsWithGreenFunction(self, green_dir):
        """
        通过格林函数计算结果
        :return:
        """
        c_load = self.femdb.load_case.c_loads
        force_nodes = [ii[0] for ii in c_load]
        directories = [ii[1] for ii in c_load]
        amps = [ii[2] for ii in c_load]
        sigma_xx = np.zeros(len(self.femdb.node_list))
        sigma_yy = np.zeros(len(self.femdb.node_list))
        sigma_zz = np.zeros(len(self.femdb.node_list))
        tau_xy = np.zeros(len(self.femdb.node_list))
        tau_yz = np.zeros(len(self.femdb.node_list))
        tau_xz = np.zeros(len(self.femdb.node_list))
        displacement_x = np.zeros(len(self.femdb.node_list))
        displacement_y = np.zeros(len(self.femdb.node_list))
        displacement_z = np.zeros(len(self.femdb.node_list))
        for ii, nd in enumerate(force_nodes):
            xyz = directories[ii]
            sigma_xx_file_path = Path(green_dir / f"sigma_xx_node_{nd}_dir_{xyz}.npy")
            sigma_yy_file_path = Path(green_dir / f"sigma_yy_node_{nd}_dir_{xyz}.npy")
            sigma_zz_file_path = Path(green_dir / f"sigma_zz_node_{nd}_dir_{xyz}.npy")
            sigma_xy_file_path = Path(green_dir / f"sigma_xy_node_{nd}_dir_{xyz}.npy")
            sigma_yz_file_path = Path(green_dir / f"sigma_yz_node_{nd}_dir_{xyz}.npy")
            sigma_xz_file_path = Path(green_dir / f"sigma_xz_node_{nd}_dir_{xyz}.npy")
            dis_file_path_x = Path(green_dir / f"dis_node_{nd}_dir_{xyz}_x.npy")
            dis_file_path_y = Path(green_dir / f"dis_node_{nd}_dir_{xyz}_y.npy")
            dis_file_path_z = Path(green_dir / f"dis_node_{nd}_dir_{xyz}_z.npy")
            if not sigma_xx_file_path.exists():
                raise FileNotFoundError(f"File Don't Exists: {sigma_xx_file_path}")
            if not sigma_yy_file_path.exists():
                raise FileNotFoundError(f"File Don't Exists: {sigma_yy_file_path}")
            if not sigma_zz_file_path.exists():
                raise FileNotFoundError(f"File Don't Exists: {sigma_zz_file_path}")
            if not sigma_xy_file_path.exists():
                raise FileNotFoundError(f"File Don't Exists: {sigma_xy_file_path}")
            if not sigma_yz_file_path.exists():
                raise FileNotFoundError(f"File Don't Exists: {sigma_yz_file_path}")
            if not sigma_xz_file_path.exists():
                raise FileNotFoundError(f"File Don't Exists: {sigma_xz_file_path}")
            if not dis_file_path_x.exists():
                raise FileNotFoundError(f"File Don't Exists: {dis_file_path_x}")
            if not dis_file_path_y.exists():
                raise FileNotFoundError(f"File Don't Exists: {dis_file_path_y}")
            if not dis_file_path_z.exists():
                raise FileNotFoundError(f"File Don't Exists: {dis_file_path_z}")
            iter_sigma_xx = np.load(sigma_xx_file_path)
            iter_sigma_yy = np.load(sigma_yy_file_path)
            iter_sigma_zz = np.load(sigma_zz_file_path)
            iter_tau_xy = np.load(sigma_xy_file_path)
            iter_tau_yz = np.load(sigma_yz_file_path)
            iter_tau_xz = np.load(sigma_xz_file_path)
            iter_dis_x = np.load(dis_file_path_x)
            iter_dis_y = np.load(dis_file_path_y)
            iter_dis_z = np.load(dis_file_path_z)
            sigma_xx += amps[ii] * iter_sigma_xx
            sigma_yy += amps[ii] * iter_sigma_yy
            sigma_zz += amps[ii] * iter_sigma_zz
            tau_xy += amps[ii] * iter_tau_xy
            tau_yz += amps[ii] * iter_tau_yz
            tau_xz += amps[ii] * iter_tau_xz
            displacement_x += amps[ii] * iter_dis_x
            displacement_y += amps[ii] * iter_dis_y
            displacement_z += amps[ii] * iter_dis_z

        term1 = (sigma_xx - sigma_yy) ** 2
        term2 = (sigma_yy - sigma_zz) ** 2
        term3 = (sigma_zz - sigma_xx) ** 2
        term4 = 6 * (tau_xy ** 2 + tau_yz ** 2 + tau_xz ** 2)
        self.femdb.linear_mises = np.sqrt(0.5 * (term1 + term2 + term3 + term4))
        self.femdb.linear_u = np.sqrt(displacement_x ** 2 + displacement_y ** 2 + displacement_z ** 2)

    def SolveDisplacement(self):
        """
        求解节点位移

        施加位移约束, 见visio文档, 求解Boundary矩阵和V矩阵, 暂未实现隐式约束, 对于显示不用求解Boundary矩阵,
        因为这种情况罚函数影响的只是Kbb的对角元素, 与未知位移求解没关系. 求解支反力也不需要加入罚函数的值

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
        c_load = self.femdb.load_case.c_loads
        force_nodes = [ii[0] for ii in c_load]
        directories = [ii[1] for ii in c_load]
        amps = [ii[2] for ii in c_load]

        for ii, nd in enumerate(force_nodes):
            eqa = self.femdb.node_hash[nd] * ModelInfo.PER_NODE_DOF
            xyz = directories[ii]
            self.right_hand[eqa + xyz] = amps[ii]

        # TODO: 没有利用Kaa是正定对称矩阵的性质, 另外Assemble对应的稀疏矩阵优化, 考虑用其他库的稀疏矩阵, 还有就是单刚的计算了
        try:
            temp_results = pypardiso.spsolve(self.femdb.global_stiff_matrix, self.right_hand)
            result_range = len(self.femdb.node_list) * ModelInfo.PER_NODE_DOF
            self.femdb.linear_u = temp_results[:result_range]
            # print(f"Displacement: {self.femdb.linear_u}")
        except ValueError as e:
            matrix_csr = self.femdb.global_stiff_matrix.tocsr()
            zero_rows = []
            for i in range(matrix_csr.shape[0]):
                row = matrix_csr.getrow(i)
                if row.nnz == 0:
                    zero_rows.append(i)
            print(f"len(zero_rows):{len(zero_rows)}\nzero_rows:{zero_rows}")
            """
            找出没有单元关系的
            """
            self.CheckRedundantNodes()
            sys.exit(1)

    def SolveStress(self, out_u=None):
        """
        求解模型节点的应力
        :return:
        """
        self.femdb.linear_mises = np.zeros(len(self.femdb.node_list))
        # self.femdb.sigma_xx = np.zeros(len(self.femdb.node_list))
        # self.femdb.sigma_yy = np.zeros(len(self.femdb.node_list))
        # self.femdb.sigma_zz = np.zeros(len(self.femdb.node_list))
        # self.femdb.tau_xy = np.zeros(len(self.femdb.node_list))
        # self.femdb.tau_xz = np.zeros(len(self.femdb.node_list))
        # self.femdb.tau_yz = np.zeros(len(self.femdb.node_list))
        for ele in self.femdb.elements:
            search_idx = []
            for ii in ele.search_node_ids:
                # start = ii * ModelInfo.PER_NODE_DOF
                # end = (ii + 1) * ModelInfo.PER_NODE_DOF
                start = self.femdb.node_list[ii].start_eq_num
                end = start + ele.node_dof_count
                search_idx.extend(np.arange(start, end, 1).tolist())

            if out_u is None:
                u = self.femdb.linear_u[search_idx].flatten()
            else:
                u = out_u[search_idx].flatten()

            """
            对于壳单元, 返回的直接是mises, 因为返回的其实是上下表面mises的最大值, 而对于一般的单元, 返回的是应力分量
            """
            stress_res = ele.CalculateElementStress(u)
            if len(stress_res) == 6:
                sigma_x, sigma_y, sigma_z, tau_yz, tau_xz, tau_xy = stress_res
                term1 = (sigma_x - sigma_y) ** 2
                term2 = (sigma_y - sigma_z) ** 2
                term3 = (sigma_z - sigma_x) ** 2
                term4 = 6 * (tau_xy ** 2 + tau_yz ** 2 + tau_xz ** 2)
                von_mises = np.sqrt(0.5 * (term1 + term2 + term3 + term4))
            elif len(stress_res) == 4 or len(stress_res) == 3:
                von_mises = stress_res
            else:
                raise ValueError(f"mises calculate error in ele's id: {ele.id}")

            for ii, n_search_id in enumerate(ele.search_node_ids):
                self.femdb.linear_mises[n_search_id] += von_mises[ii] / self.femdb.node_connected_element_count[n_search_id]
                # self.femdb.sigma_xx[n_search_id] += sigma_x[ii] / self.femdb.node_connected_element_count[n_search_id]
                # self.femdb.sigma_yy[n_search_id] += sigma_y[ii] / self.femdb.node_connected_element_count[n_search_id]
                # self.femdb.sigma_zz[n_search_id] += sigma_z[ii] / self.femdb.node_connected_element_count[n_search_id]
                # self.femdb.tau_xy[n_search_id] += tau_xy[ii] / self.femdb.node_connected_element_count[n_search_id]
                # self.femdb.tau_yz[n_search_id] += tau_yz[ii] / self.femdb.node_connected_element_count[n_search_id]
                # self.femdb.tau_xz[n_search_id] += tau_xz[ii] / self.femdb.node_connected_element_count[n_search_id]

        if out_u is not None:
            return deepcopy(self.femdb.linear_mises)
        else:
            return None

    def CalculateDynResultsWithGreenFunction(self, green_dir):
        """
        通过格林函数计算结果
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
        node_count = len(self.femdb.node_list)
        fem_all_dofs = ModelInfo.PER_NODE_DOF * node_count
        his_vals = np.zeros((fem_all_dofs, steps_count), dtype=float)

        for ii, nd in enumerate(force_nodes):
            eqa = self.femdb.node_hash[nd] * ModelInfo.PER_NODE_DOF
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
        solver = factorized(K_hat.tocsc())

        """
        开始计算各个时间步的值
        """
        u = [np.zeros(fem_all_dofs, dtype=float)]
        v = [np.zeros(fem_all_dofs, dtype=float)]
        a = [acc_0]

        """
        读取格林函数
        """
        for ii in range(steps_count - 1):
            """
            结果置零
            """
            sigma_xx = np.zeros(len(self.femdb.node_list))
            sigma_yy = np.zeros(len(self.femdb.node_list))
            sigma_zz = np.zeros(len(self.femdb.node_list))
            tau_xy = np.zeros(len(self.femdb.node_list))
            tau_yz = np.zeros(len(self.femdb.node_list))
            tau_xz = np.zeros(len(self.femdb.node_list))
            displacement_x = np.zeros(len(self.femdb.node_list))
            displacement_y = np.zeros(len(self.femdb.node_list))
            displacement_z = np.zeros(len(self.femdb.node_list))

            f_hat = his_vals[:, ii + 1] + self.femdb.global_mass_matrix @ (a0 * u[ii] + a2 * v[ii] + a3 * a[ii])
            u_next = solver(f_hat)
            a_next = a0 * (u_next - u[ii]) - a2 * v[ii] - a3 * a[ii]
            v_next = v[ii] + a6 * a[ii] + a7 * a_next
            u.append(u_next)
            v.append(v_next)
            a.append(a_next)

        term1 = (sigma_xx - sigma_yy) ** 2
        term2 = (sigma_yy - sigma_zz) ** 2
        term3 = (sigma_zz - sigma_xx) ** 2
        term4 = 6 * (tau_xy ** 2 + tau_yz ** 2 + tau_xz ** 2)
        self.femdb.linear_mises = np.sqrt(0.5 * (term1 + term2 + term3 + term4))
        self.femdb.linear_u = np.sqrt(displacement_x ** 2 + displacement_y ** 2 + displacement_z ** 2)

    def AddBoundaryByPenalty(self):
        """
        罚函数的方法施加约束
        @return:
        """
        """
        注意要增加约束方程部分
        """
        node_count = len(self.femdb.node_list)
        ce_count = len(self.femdb.equation_constrain_couple)
        ce_add_equations = 0
        for ii in range(ce_count):
            iter_ce = self.femdb.equation_constrain_idx[ii]
            ce_add_equations += len(iter_ce)
        self.right_hand = np.zeros(node_count * ModelInfo.PER_NODE_DOF + ce_add_equations)

        bds = self.femdb.load_case.GetBoundaries()
        for bd in bds:
            for row in bd:
                nd = int(row[0])
                dir_ = int(row[1])
                val = float(row[2])
                nd_idx = self.femdb.node_hash[nd]
                base_idx = nd_idx * ModelInfo.PER_NODE_DOF
                self.femdb.global_stiff_matrix[base_idx + dir_, base_idx + dir_] += 2e21
                self.right_hand[base_idx + dir_] = val

    def NewMarkExplict(self):
        """
        显式的动力学求解，使用NewMark方法
        # TODO: 使用Bathe书上算例测试该算法
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
        node_count = len(self.femdb.node_list)
        fem_all_dofs = ModelInfo.PER_NODE_DOF * node_count
        his_vals = np.zeros((fem_all_dofs, steps_count), dtype=float)

        for ii, nd in enumerate(force_nodes):
            eqa = self.femdb.node_hash[nd] * ModelInfo.PER_NODE_DOF
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
        solver = factorized(K_hat.tocsc())

        """
        开始计算各个时间步的值
        """
        u = [np.zeros(fem_all_dofs, dtype=float)]
        v = [np.zeros(fem_all_dofs, dtype=float)]
        s = [np.zeros(len(self.femdb.node_list), dtype=float)]
        a = [acc_0]
        for ii in range(steps_count - 1):
            f_hat = his_vals[:, ii + 1] + self.femdb.global_mass_matrix @ (a0 * u[ii] + a2 * v[ii] + a3 * a[ii])
            u_next = solver(f_hat)
            a_next = a0 * (u_next - u[ii]) - a2 * v[ii] - a3 * a[ii]
            v_next = v[ii] + a6 * a[ii] + a7 * a_next
            u.append(u_next)
            v.append(v_next)
            a.append(a_next)
            s.append(self.SolveStress(u_next))

        self.femdb.history_u = u
        self.femdb.history_v = v
        self.femdb.history_a = a
        self.femdb.history_s = s
        self.femdb.history_step_count = len(u)

    def GenerateNewMarkGeLinFile(self, output_dir):
        """
        显式的动力学求解，使用NewMark方法
        # TODO: 使用Bathe书上算例测试该算法
        Reference:
           《结构动力学基础》 张亚辉、林家浩 P94
        :return:
        """
        """
        初始化载荷, 计算初始加速度
        """
        his_load = self.femdb.load_case.history_loads
        delta_t = his_load[0][-1][0, 1] - his_load[0][-1][0, 0]
        alpha = 0.25
        a0 = 1 / alpha / delta_t ** 2
        K_hat = self.femdb.global_stiff_matrix + a0 * self.femdb.global_mass_matrix

        force_nodes = [ii[0] for ii in his_load]
        directories = [ii[1] for ii in his_load]
        right_hand_back_up = deepcopy(self.right_hand)
        for ii, nd in enumerate(force_nodes):
            """
            保存位移结果
            """
            eqa = self.femdb.node_hash[nd] * ModelInfo.PER_NODE_DOF
            xyz = directories[ii]
            self.right_hand = deepcopy(right_hand_back_up)
            self.right_hand[eqa + xyz] = 1
            temp_results = pypardiso.spsolve(K_hat, self.right_hand)
            result_range = len(self.femdb.node_list) * ModelInfo.PER_NODE_DOF
            iter_linear_u = temp_results[:result_range]
            np.save(output_dir / f"dis_node_{nd}.npy", iter_linear_u)
            """
            求解应力结果并保存
            """
            iter_sigma_xx = np.zeros(len(self.femdb.node_list))
            iter_sigma_yy = np.zeros(len(self.femdb.node_list))
            iter_sigma_zz = np.zeros(len(self.femdb.node_list))
            iter_tau_xy = np.zeros(len(self.femdb.node_list))
            iter_tau_xz = np.zeros(len(self.femdb.node_list))
            iter_tau_yz = np.zeros(len(self.femdb.node_list))
            for ele in self.femdb.elements:
                search_idx = []
                for ii in ele.search_node_ids:
                    start = ii * ModelInfo.PER_NODE_DOF
                    end = (ii + 1) * ModelInfo.PER_NODE_DOF
                    search_idx.extend(np.arange(start, end, 1).tolist())

                iter_u = iter_linear_u[search_idx].flatten()
                sigma_x, sigma_y, sigma_z, tau_yz, tau_xz, tau_xy = ele.CalculateElementStress(iter_u)

                for ii, n_search_id in enumerate(ele.search_node_ids):
                    iter_sigma_xx[n_search_id] += sigma_x[ii] / self.femdb.node_connected_element_count[n_search_id]
                    iter_sigma_yy[n_search_id] += sigma_y[ii] / self.femdb.node_connected_element_count[n_search_id]
                    iter_sigma_zz[n_search_id] += sigma_z[ii] / self.femdb.node_connected_element_count[n_search_id]
                    iter_tau_xy[n_search_id] += tau_xy[ii] / self.femdb.node_connected_element_count[n_search_id]
                    iter_tau_yz[n_search_id] += tau_yz[ii] / self.femdb.node_connected_element_count[n_search_id]
                    iter_tau_xz[n_search_id] += tau_xz[ii] / self.femdb.node_connected_element_count[n_search_id]

            np.save(output_dir / f"sigma_xx_node_{nd}.npy", iter_sigma_xx)
            np.save(output_dir / f"sigma_yy_node_{nd}.npy", iter_sigma_yy)
            np.save(output_dir / f"sigma_zz_node_{nd}.npy", iter_sigma_zz)
            np.save(output_dir / f"sigma_xy_node_{nd}.npy", iter_tau_xy)
            np.save(output_dir / f"sigma_yz_node_{nd}.npy", iter_tau_yz)
            np.save(output_dir / f"sigma_xz_node_{nd}.npy", iter_tau_xz)


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
