#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from element.Plate import *
from element.Membrane import *


class DKTShell(ElementBaseClass, ABC):
    """ DKTShell Element class """

    def __init__(self, eid=None):
        super().__init__(eid)
        self.nodes_count = 3  # Each element has 3 nodes
        self._nodes = [None for _ in range(self.nodes_count)]
        self._vtp_type = "triangle"
        self.K = np.zeros((18,18), dtype=float)
        self.unv_code = 30500

    def CalElementDMatrix(self, an_type=None):
        """
        计算本构矩阵, 弹性模量和泊松比, Bathe 上册P184
        """
        e = self.cha_dict[MaterialKey.E]
        niu = self.cha_dict[MaterialKey.Niu]
        a = e / ((1 + niu) * (1 - 2 * niu))
        self.D = a * np.array([[1 - niu, niu, niu, 0, 0, 0],
                               [niu, 1 - niu, niu, 0, 0, 0],
                               [niu, niu, 1 - niu, 0, 0, 0],
                               [0, 0, 0, (1 - 2 * niu) / 2., 0, 0],
                               [0, 0, 0, 0, (1 - 2 * niu) / 2., 0],
                               [0, 0, 0, 0, 0, (1 - 2 * niu) / 2.]])

    def ElementStiffness(self):
        """
        TODO: 转轴要不要加小量
        """
        """
        由膜单元和板单元构成
        """
        plate = DKTPlate(-1)
        membrane = CPM6(-1)

        """
        先转换到局部坐标
        """
        T_matrix = GetGlobal2LocalTransMatrix(self.node_coords)
        local_coord = np.matmul(self.node_coords, T_matrix)
        m_coords = local_coord[:, :2]
        mid_node = np.asarray([(local_coord[0, :] + local_coord[1, :]) * 0.5,
                               (local_coord[1, :] + local_coord[2, :]) * 0.5,
                               (local_coord[2, :] + local_coord[0, :]) * 0.5], dtype=float)[:, :2]

        """
        设置膜单元和板单元的节点坐标, 以及全局和局部坐标系的转换矩阵
        """
        membrane.node_coords = np.append(m_coords, mid_node, axis=0)
        membrane.T_matrix = T_matrix
        plate.node_coords = m_coords
        plate.T_matrix = T_matrix

        """
        设置膜单元和版单元的材料
        """
        membrane.cha_dict = self.cha_dict
        plate.cha_dict = self.cha_dict
        membrane.CalElementDMatrix()
        plate.CalElementDMatrix()

        """
        Assembly Stiffness Matrix, membrane: u, v, theta_z, plate: omega, theta_x, theta_y
        """
        # e = 10e-8
        k_mtx_m = membrane.ElementStiffness()
        k_mtx_p = plate.ElementStiffness()

        index_m_g = [(0, 1), (5, 7), (11, 13)]
        index_m = [(0, 1), (2, 4), (5, 7)]
        index_p_g = [(2, 4), (8, 10), (14, 16)]
        index_p = [(0, 2), (3, 5), (6, 8)]
        for ii in range(3):
            m_row_s = index_m[ii][0]
            m_row_e = index_m[ii][1] + 1
            m_row_g_s = index_m_g[ii][0]
            m_row_g_e = index_m_g[ii][1] + 1

            p_row_s = index_p[ii][0]
            p_row_e = index_p[ii][1] + 1
            p_row_g_s = index_p_g[ii][0]
            p_row_g_e = index_p_g[ii][1] + 1
            for jj in range(3):
                m_col_s = index_m[jj][0]
                m_col_e = index_m[jj][1] + 1

                m_col_g_s = index_m_g[jj][0]
                m_col_g_e = index_m_g[jj][1] + 1
                self.K[m_row_g_s:m_row_g_e, m_col_g_s:m_col_g_e] = k_mtx_m[m_row_s:m_row_e, m_col_s:m_col_e]

                p_col_s = index_p[jj][0]
                p_col_e = index_p[jj][1] + 1

                p_col_g_s = index_p_g[jj][0]
                p_col_g_e = index_p_g[jj][1] + 1

                self.K[p_row_g_s:p_row_g_e, p_col_g_s:p_col_g_e] = k_mtx_p[p_row_s:p_row_e, p_col_s:p_col_e]

            self.K[m_row_g_s:m_row_g_e, -1] = k_mtx_m[m_row_s:m_row_e, -1]
            self.K[-1, m_row_g_s:m_row_g_e] = k_mtx_m[-1, m_row_s:m_row_e]

        self.K[-1,-1] = k_mtx_m[-1,-1]

        global_t_matrix = np.zeros((18, 18))
        global_t_matrix[0:3, 0:3] = T_matrix
        global_t_matrix[3:6, 3:6] = T_matrix
        global_t_matrix[6:9, 6:9] = T_matrix
        global_t_matrix[9:12, 9:12] = T_matrix
        global_t_matrix[12:15, 12:15] = T_matrix
        global_t_matrix[15:18, 15:18] = T_matrix

        self.K = global_t_matrix.T @ self.K @ global_t_matrix

        return self.K

    def ElementStress(self, displacement):
        """
        Calculate element stress
        """


class DKQShell(ElementBaseClass, ABC):
    """ DKQShell Element class """

    def __init__(self, eid=None):
        super().__init__(eid)
        self.nodes_count = 4  # Each element has 4 nodes
        self._vtp_type = "quad"
        self.K = np.zeros((24, 24))
        self._nodes = [None for _ in range(self.nodes_count)]
        self.unv_code = 40500

    def CalElementDMatrix(self, an_type=None):
        """
        计算本构矩阵, 弹性模量和泊松比, Bathe 上册P184
        """
        e = self.cha_dict[MaterialKey.E]
        niu = self.cha_dict[MaterialKey.Niu]
        a = e / ((1 + niu) * (1 - 2 * niu))
        self.D = a * np.array([[1 - niu, niu, niu, 0, 0, 0],
                               [niu, 1 - niu, niu, 0, 0, 0],
                               [niu, niu, 1 - niu, 0, 0, 0],
                               [0, 0, 0, (1 - 2 * niu) / 2., 0, 0],
                               [0, 0, 0, 0, (1 - 2 * niu) / 2., 0],
                               [0, 0, 0, 0, 0, (1 - 2 * niu) / 2.]])

    def ElementStiffness(self):
        """
        壳的刚度阵由膜单元和板单元构成
        """
        plate = DKQPlate(self.id)
        membrane = CPM8(self.id)

        """
        先转换到局部坐标
        """
        T_matrix = GetGlobal2LocalTransMatrix(self.node_coords)
        local_coord = np.matmul(self.node_coords, T_matrix)
        m_coords = local_coord[:, :2]
        mid_node = np.asarray([(local_coord[0, :] + local_coord[1, :]) * 0.5,
                               (local_coord[1, :] + local_coord[2, :]) * 0.5,
                               (local_coord[2, :] + local_coord[3, :]) * 0.5,
                               (local_coord[3, :] + local_coord[0, :]) * 0.5], dtype=float)[:, :2]

        """
        设置膜单元和板单元的节点坐标, 以及全局和局部坐标系的转换矩阵
        """
        membrane.node_coords = np.append(m_coords, mid_node, axis=0)
        membrane.T_matrix = T_matrix
        plate.node_coords = m_coords
        plate.T_matrix = T_matrix

        """
        设置膜单元和版单元的材料
        """
        membrane.cha_dict = self.cha_dict
        plate.cha_dict = self.cha_dict
        membrane.CalElementDMatrix()
        plate.CalElementDMatrix()

        """
        Assembly Stiffness Matrix, membrane: u,v,theta_z, plate: omega, theta_x, theta_y
        """
        # e = 10e-8
        k_mtx_m = membrane.ElementStiffness()
        k_mtx_p = plate.ElementStiffness()

        index_m_g = [(0, 1), (5, 7), (11, 13), (17, 19)]
        index_m = [(0, 1), (2, 4), (5, 7), (8, 10)]
        index_p_g = [(2, 4), (8, 10), (14, 16), (20, 22)]
        index_p = [(0, 2), (3, 5), (6, 8), (9, 11)]

        for ii in range(4):
            m_row_s = index_m[ii][0]
            m_row_e = index_m[ii][1] + 1
            m_row_g_s = index_m_g[ii][0]
            m_row_g_e = index_m_g[ii][1] + 1

            p_row_s = index_p[ii][0]
            p_row_e = index_p[ii][1] + 1
            p_row_g_s = index_p_g[ii][0]
            p_row_g_e = index_p_g[ii][1] + 1
            for jj in range(4):
                m_col_s = index_m[jj][0]
                m_col_e = index_m[jj][1] + 1

                m_col_g_s = index_m_g[jj][0]
                m_col_g_e = index_m_g[jj][1] + 1
                self.K[m_row_g_s:m_row_g_e, m_col_g_s:m_col_g_e] = k_mtx_m[m_row_s:m_row_e, m_col_s:m_col_e]

                p_col_s = index_p[jj][0]
                p_col_e = index_p[jj][1] + 1

                p_col_g_s = index_p_g[jj][0]
                p_col_g_e = index_p_g[jj][1] + 1

                self.K[p_row_g_s:p_row_g_e, p_col_g_s:p_col_g_e] = k_mtx_p[p_row_s:p_row_e, p_col_s:p_col_e]

            self.K[m_row_g_s:m_row_g_e, -1] = k_mtx_m[m_row_s:m_row_e, -1]
            self.K[-1, m_row_g_s:m_row_g_e] = k_mtx_m[-1, m_row_s:m_row_e]

        self.K[-1,-1] = k_mtx_m[-1,-1]

        # try:
        #     np.linalg.inv(k_mtx_p_g)
        # except np.linalg.LinAlgError:
        #     print(self.id)
        #
        # try:
        #     np.linalg.inv(k_mtx_m_g)
        # except np.linalg.LinAlgError:
        #     print("mem ", self.id)
        #
        global_t_matrix = np.zeros((24, 24))
        global_t_matrix[0:3, 0:3] = T_matrix
        global_t_matrix[3:6, 3:6] = T_matrix
        global_t_matrix[6:9, 6:9] = T_matrix
        global_t_matrix[9:12, 9:12] = T_matrix
        global_t_matrix[12:15, 12:15] = T_matrix
        global_t_matrix[15:18, 15:18] = T_matrix
        global_t_matrix[18:21, 18:21] = T_matrix
        global_t_matrix[21:24, 21:24] = T_matrix

        self.K = global_t_matrix.T @ self.K @ global_t_matrix

        return self.K

    def ElementStress(self, displacement):
            """
            Calculate element stress
            """



