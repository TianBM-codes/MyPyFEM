#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from element.Plate import *  # 导入板单元相关类
from element.Membrane import *  # 导入膜单元相关类


class CookTriShell(ElementBaseClass, ABC):
    """
    DKT三角形壳单元类
    结合了DKT板单元和CPM3膜单元的特性
    """

    def __init__(self, eid=None):
        """
        初始化三角形壳单元
        :param eid: 单元ID
        """
        super().__init__(eid)
        self.nodes_count = 3  # 三角形单元有3个节点
        self.vtu_type = "triangle"  # VTK可视化类型
        self.K = np.zeros((18, 18), dtype=float)  # 单元刚度矩阵(18自由度)
        self.local2global_matrix = np.zeros((18, 18))  # 局部到全局坐标转换矩阵
        self.local_coord = None  # 局部坐标系下的节点坐标
        self.T_matrix = None  # 全局到局部坐标转换矩阵
        self.plate = DKTPlate(-1)  # DKT板单元实例
        self.membrane = CPM3(-1)  # CPM3膜单元实例
        self.last_e = None

    def CalElementDMatrix(self, an_type=None):
        """
        计算本构矩阵(弹性模量和泊松比)
        参考Bathe《有限元方法》上册P184
        """
        pass

    def ElementStiffness(self, from_origin=False):
        """
        计算单元刚度矩阵
        :param from_origin: 是否从原始数据重新计算
        :return: 单元刚度矩阵
        """
        # 提取局部坐标中的x,y坐标用于膜单元
        m_coords = self.local_coord[:, :2]
        # 计算边中点坐标
        mid_node = np.asarray([(self.local_coord[0, :] + self.local_coord[1, :]) * 0.5,
                               (self.local_coord[1, :] + self.local_coord[2, :]) * 0.5,
                               (self.local_coord[2, :] + self.local_coord[0, :]) * 0.5], dtype=float)[:, :2]

        # 设置膜单元和板单元的节点坐标和转换矩阵
        self.membrane.node_coords = np.append(m_coords, mid_node, axis=0)
        self.membrane.T_matrix = self.T_matrix
        self.plate.node_coords = m_coords
        self.plate.T_matrix = self.T_matrix

        # 传递材料和截面属性
        self.membrane.cha_dict = self.cha_dict
        self.membrane.sec_id = self.sec_id
        self.plate.cha_dict = self.cha_dict
        self.plate.sec_id = self.sec_id

        # 计算本构矩阵
        self.membrane.CalElementDMatrix()
        self.plate.CalElementDMatrix()

        # 组装刚度矩阵: 膜单元处理u,v,θz; 板单元处理ω,θx,θy
        if from_origin:
            k_mtx_m = self.membrane.ReCalculateElementStiffness()
            k_mtx_p = self.plate.ReCalculateElementStiffness()
            self.K = np.zeros((18, 18), dtype=float)
        else:
            k_mtx_m = self.membrane.ElementStiffness()
            k_mtx_p = self.plate.ElementStiffness()

        # 定义自由度映射关系
        index_m_g = [(0, 1), (5, 7), (11, 13)]  # 膜单元在全局矩阵中的位置
        index_m = [(0, 1), (2, 4), (5, 7)]  # 膜单元在局部矩阵中的位置
        index_p_g = [(2, 4), (8, 10), (14, 16)]  # 板单元在全局矩阵中的位置
        index_p = [(0, 2), (3, 5), (6, 8)]  # 板单元在局部矩阵中的位置

        # 组装刚度矩阵
        for ii in range(3):
            # 膜单元部分
            m_row_s, m_row_e = index_m[ii]
            m_row_g_s, m_row_g_e = index_m_g[ii]
            # 板单元部分
            p_row_s, p_row_e = index_p[ii]
            p_row_g_s, p_row_g_e = index_p_g[ii]

            for jj in range(3):
                # 组装膜单元刚度
                m_col_s, m_col_e = index_m[jj]
                m_col_g_s, m_col_g_e = index_m_g[jj]
                self.K[m_row_g_s:m_row_g_e + 1, m_col_g_s:m_col_g_e + 1] = k_mtx_m[m_row_s:m_row_e + 1, m_col_s:m_col_e + 1]

                # 组装板单元刚度
                p_col_s, p_col_e = index_p[jj]
                p_col_g_s, p_col_g_e = index_p_g[jj]
                self.K[p_row_g_s:p_row_g_e + 1, p_col_g_s:p_col_g_e + 1] = k_mtx_p[p_row_s:p_row_e + 1, p_col_s:p_col_e + 1]

            # 处理载荷项
            self.K[m_row_g_s:m_row_g_e + 1, -1] = k_mtx_m[m_row_s:m_row_e + 1, -1]
            self.K[-1, m_row_g_s:m_row_g_e + 1] = k_mtx_m[-1, m_row_s:m_row_e + 1]

        self.K[-1, -1] = k_mtx_m[-1, -1]

        # 转换到全局坐标系
        self.K = self.local2global_matrix.T @ self.K @ self.local2global_matrix

        self.last_e = self.cha_dict[MaterialKey.E]
        return self.K

    def CalculateElementStress(self, displacement):
        """
        计算单元应力
        :param displacement: 节点位移向量
        :return: 应力分量数组(sigma_xx, sigma_yy, sigma_zz, tau_yz, tau_xz, tau_xy)
        """
        # 转换到位移局部坐标系
        local_dis = self.local2global_matrix @ displacement

        # 提取膜单元和板单元对应的位移分量
        membrane_indices = [i * 6 + j for i in range(3) for j in [0, 1, 5]]  # u,v,θz
        plate_indices = [i * 6 + j for i in range(3) for j in [2, 3, 4]]  # ω,θx,θy

        # 分别计算膜应力和板应力
        membrane_stress = self.membrane.CalculateElementStress(local_dis[membrane_indices])
        plate_stress = self.plate.CalculateElementStress(local_dis[plate_indices])

        # 提取各应力分量
        sigma_xx, sigma_yy, tau_xy = membrane_stress[:, 0], membrane_stress[:, 1], membrane_stress[:, 2]
        sigma_zz, tau_yz, tau_xz = plate_stress[:, 0], plate_stress[:, 1], plate_stress[:, 2]

        return np.array([sigma_xx, sigma_yy, sigma_zz, tau_yz, tau_xz, tau_xy])

    def ElementMass(self):
        """计算单元质量矩阵(待实现)"""
        pass

    def CalculateBasic(self):
        """计算基本参数: 坐标转换矩阵和局部坐标"""
        # 获取全局到局部的转换矩阵和原点
        self.T_matrix, origin = GetShellGlobal2LocalTransMatrix(self.node_coords)
        R_matrix = self.T_matrix.T  # 局部到全局的转换矩阵

        # 组装完整的转换矩阵(18x18)
        for i in range(6):
            self.local2global_matrix[i * 3:(i + 1) * 3, i * 3:(i + 1) * 3] = R_matrix

        # 计算局部坐标
        self.local_coord = (self.node_coords.T - origin[:, np.newaxis]).T @ self.T_matrix

    def ReCalculateElementStiffness(self):
        """
        重新计算单元刚度矩阵(更新材料属性后)
        :return: 更新后的单元刚度矩阵
        """
        scale = self.cha_dict[MaterialKey.E] / self.last_e
        self.last_e = self.cha_dict[MaterialKey.E]
        return self.K * scale


class CookQuaShell(ElementBaseClass, ABC):
    """
    四边形壳单元类
    结合了DKQ板单元和CPM4膜单元的特性
    """

    def __init__(self, eid=None):
        """
        初始化四边形壳单元
        :param eid: 单元ID
        """
        super().__init__(eid)
        self.nodes_count = 4  # 四边形单元有4个节点
        self.vtu_type = "quad"  # VTK可视化类型
        self.K = np.zeros((24, 24))  # 单元刚度矩阵(24自由度)
        self.local_coord = None  # 局部坐标系下的节点坐标
        self.local2global_matrix = np.zeros((24, 24))  # 局部到全局坐标转换矩阵
        self.plate = DKQPlate(self.id)  # DKQ板单元实例
        self.membrane = CPM4(self.id)  # CPM4膜单元实例
        self.T_matrix = None  # 全局到局部坐标转换矩阵
        self.last_e = None  # 上次计算单元的弹性模量

    def CalElementDMatrix(self, an_type=None):
        """
        计算本构矩阵(弹性模量和泊松比)
        参考Bathe《有限元方法》上册P184
        """
        pass

    def ElementStiffness(self, from_origin=False):
        """
        计算单元刚度矩阵
        :param from_origin: 是否从原始数据重新计算
        :return: 单元刚度矩阵
        """
        # 提取局部坐标中的x,y坐标用于膜单元
        m_coords = self.local_coord[:, :2]
        # 计算边中点坐标
        mid_node = np.asarray([(self.local_coord[0, :] + self.local_coord[1, :]) * 0.5,
                               (self.local_coord[1, :] + self.local_coord[2, :]) * 0.5,
                               (self.local_coord[2, :] + self.local_coord[3, :]) * 0.5,
                               (self.local_coord[3, :] + self.local_coord[0, :]) * 0.5], dtype=float)[:, :2]

        # 设置膜单元和板单元的节点坐标和转换矩阵
        self.membrane.node_coords = np.append(m_coords, mid_node, axis=0)
        self.membrane.T_matrix = self.T_matrix
        self.plate.node_coords = m_coords
        self.plate.T_matrix = self.T_matrix

        # 传递材料和截面属性
        self.membrane.cha_dict = self.cha_dict
        self.membrane.sec_id = self.sec_id
        self.plate.cha_dict = self.cha_dict
        self.plate.sec_id = self.sec_id

        # 计算本构矩阵
        self.membrane.CalElementDMatrix()
        self.plate.CalElementDMatrix()

        # 组装刚度矩阵: 膜单元处理u,v,θz; 板单元处理ω,θx,θy
        if from_origin:
            k_mtx_m = self.membrane.ReCalculateElementStiffness()
            k_mtx_p = self.plate.ReCalculateElementStiffness()
        else:
            k_mtx_m = self.membrane.ElementStiffness()
            k_mtx_p = self.plate.ElementStiffness()

        # 定义自由度映射关系
        index_m_g = [(0, 1), (5, 7), (11, 13), (17, 19)]  # 膜单元在全局矩阵中的位置
        index_m = [(0, 1), (2, 4), (5, 7), (8, 10)]  # 膜单元在局部矩阵中的位置
        index_p_g = [(2, 4), (8, 10), (14, 16), (20, 22)]  # 板单元在全局矩阵中的位置
        index_p = [(0, 2), (3, 5), (6, 8), (9, 11)]  # 板单元在局部矩阵中的位置

        # 组装刚度矩阵
        for ii in range(4):
            # 膜单元部分
            m_row_s, m_row_e = index_m[ii]
            m_row_g_s, m_row_g_e = index_m_g[ii]
            # 板单元部分
            p_row_s, p_row_e = index_p[ii]
            p_row_g_s, p_row_g_e = index_p_g[ii]

            for jj in range(4):
                # 组装膜单元刚度
                m_col_s, m_col_e = index_m[jj]
                m_col_g_s, m_col_g_e = index_m_g[jj]
                self.K[m_row_g_s:m_row_g_e + 1, m_col_g_s:m_col_g_e + 1] = k_mtx_m[m_row_s:m_row_e + 1, m_col_s:m_col_e + 1]

                # 组装板单元刚度
                p_col_s, p_col_e = index_p[jj]
                p_col_g_s, p_col_g_e = index_p_g[jj]
                self.K[p_row_g_s:p_row_g_e + 1, p_col_g_s:p_col_g_e + 1] = k_mtx_p[p_row_s:p_row_e + 1, p_col_s:p_col_e + 1]

            # 处理载荷项
            self.K[m_row_g_s:m_row_g_e + 1, -1] = k_mtx_m[m_row_s:m_row_e + 1, -1]
            self.K[-1, m_row_g_s:m_row_g_e + 1] = k_mtx_m[-1, m_row_s:m_row_e + 1]

        self.K[-1, -1] = k_mtx_m[-1, -1]

        # 转换到全局坐标系
        self.K = self.local2global_matrix.T @ self.K @ self.local2global_matrix
        self.last_e = self.cha_dict[MaterialKey.E]
        return self.K

    def CalculateElementStress(self, displacement):
        """
        计算单元应力
        :param displacement: 节点位移向量
        :return: 应力分量数组(sigma_xx, sigma_yy, sigma_zz, tau_yz, tau_xz, tau_xy)
        """
        # 转换到位移局部坐标系
        local_dis = self.local2global_matrix @ displacement

        # 提取膜单元和板单元对应的位移分量
        membrane_indices = [i * 6 + j for i in range(4) for j in [0, 1, 5]]  # u,v,θz
        plate_indices = [i * 6 + j for i in range(4) for j in [2, 3, 4]]  # ω,θx,θy

        # 分别计算膜应力和板应力
        membrane_stress = self.membrane.CalculateElementStress(local_dis[membrane_indices])
        plate_stress = self.plate.CalculateElementStress(local_dis[plate_indices])

        # 提取各应力分量
        sigma_xx, sigma_yy, tau_xy = membrane_stress[:, 0], membrane_stress[:, 1], membrane_stress[:, 2]
        sigma_zz, tau_yz, tau_xz = plate_stress[:, 0], plate_stress[:, 1], plate_stress[:, 2]

        return np.array([sigma_xx, sigma_yy, sigma_zz, tau_yz, tau_xz, tau_xy])

    def ElementMass(self):
        """
        计算单元的质量矩阵(协调质量矩阵)
        :return: 单元质量矩阵
        """
        # 提取局部坐标中的x,y坐标
        m_coords = self.local_coord[:, :2]
        self.M = np.zeros((24, 24), dtype=float)

        # 获取材料属性
        rho = self.cha_dict[MaterialKey.Density]  # 密度
        thickness = self.cha_dict[MaterialKey.Thickness]  # 厚度

        # 获取高斯积分点和权重
        sample_pt, weight = GaussIntegrationPoint.GetSamplePointAndWeight(2)

        # 高斯积分计算质量矩阵
        for ri in range(2):
            for si in range(2):
                r, s = sample_pt[ri], sample_pt[si]
                # 形函数
                N1 = 0.25 * (1 - r) * (1 - s)
                N2 = 0.25 * (1 + r) * (1 - s)
                N3 = 0.25 * (1 + r) * (1 + s)
                N4 = 0.25 * (1 - r) * (1 + s)

                # 形函数矩阵
                H = np.zeros((3, 12))
                H[0, 0] = N1
                H[1, 1] = N1
                H[2, 2] = N1
                H[0, 3] = N2
                H[1, 4] = N2
                H[2, 5] = N2
                H[0, 6] = N3
                H[1, 7] = N3
                H[2, 8] = N3
                H[0, 9] = N4
                H[1, 10] = N4
                H[2, 11] = N4

                # 形函数导数
                dNdr = np.array([[0.25 * (1 + s), -0.25 * (1 + s), 0.25 * (s - 1), 0.25 * (1 - s)],
                                 [0.25 * (1 + r), 0.25 * (1 - r), 0.25 * (r - 1), -0.25 * (1 + r)]], dtype=float)
                g_weight = weight[ri] * weight[si]  # 积分权重

                # Jacobian行列式
                det_J = np.linalg.det(dNdr @ m_coords)

                # 计算当前积分点的质量矩阵贡献
                iter_mass = H.T @ H * rho * det_J * g_weight * thickness  # 12 * 12

                # 组装到全局质量矩阵
                for jj in range(4):
                    for kk in range(4):
                        row_idx = jj * 6
                        col_idx = kk * 6
                        self.M[row_idx:row_idx + 3, col_idx:col_idx + 3] += \
                            iter_mass[jj * 3:(jj + 1) * 3, kk * 3:(kk + 1) * 3]

        # 转换到全局坐标系
        return self.local2global_matrix.T @ self.M @ self.local2global_matrix

    def CalculateBasic(self):
        """
        计算基本参数: 坐标转换矩阵和局部坐标
        """
        # 获取全局到局部的转换矩阵和原点
        self.T_matrix, origin = GetShellGlobal2LocalTransMatrix(self.node_coords)
        R_matrix = self.T_matrix.T  # 局部到全局的转换矩阵

        # 组装完整的转换矩阵(24x24)
        for i in range(8):
            self.local2global_matrix[i * 3:(i + 1) * 3, i * 3:(i + 1) * 3] = R_matrix

        # 计算局部坐标
        self.local_coord = (self.node_coords.T - origin[:, np.newaxis]).T @ self.T_matrix

    def ReCalculateElementStiffness(self):
        """
        重新计算单元刚度矩阵(更新材料属性后)
        :return: 更新后的单元刚度矩阵
        """
        scale = self.cha_dict[MaterialKey.E] / self.last_e
        self.last_e = self.cha_dict[MaterialKey.E]
        return self.K * scale
