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
        self.K = None
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

        k_mtx_m_g = np.insert(k_mtx_m, 2, [[0] * 9] * 3, axis=0)
        k_mtx_m_g = np.insert(k_mtx_m_g, 8, [[0] * 9] * 3, axis=0)
        k_mtx_m_g = np.insert(k_mtx_m_g, 14, [[0] * 9] * 3, axis=0)

        k_mtx_m_g = np.insert(k_mtx_m_g, 2, [[0] * 18] * 3, axis=1)
        k_mtx_m_g = np.insert(k_mtx_m_g, 8, [[0] * 18] * 3, axis=1)
        k_mtx_m_g = np.insert(k_mtx_m_g, 14, [[0] * 18] * 3, axis=1)

        k_mtx_p_g = np.insert(k_mtx_p, 0, [[0] * 9] * 2, axis=0)
        k_mtx_p_g = np.insert(k_mtx_p_g, 5, [[0] * 9] * 3, axis=0)
        k_mtx_p_g = np.insert(k_mtx_p_g, 11, [[0] * 9] * 3, axis=0)
        k_mtx_p_g = np.insert(k_mtx_p_g, 17, [0] * 9, axis=0)

        k_mtx_p_g = np.insert(k_mtx_p_g, 0, [[0] * 18] * 2, axis=1)
        k_mtx_p_g = np.insert(k_mtx_p_g, 5, [[0] * 18] * 3, axis=1)
        k_mtx_p_g = np.insert(k_mtx_p_g, 11, [[0] * 18] * 3, axis=1)
        k_mtx_p_g = np.insert(k_mtx_p_g, 17, [0] * 18, axis=1)

        self.K = k_mtx_m_g + k_mtx_p_g

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

        k_mtx_m_g = np.insert(k_mtx_m, 2, [[0] * 12] * 3, axis=0)
        k_mtx_m_g = np.insert(k_mtx_m_g, 8, [[0] * 12] * 3, axis=0)
        k_mtx_m_g = np.insert(k_mtx_m_g, 14, [[0] * 12] * 3, axis=0)
        k_mtx_m_g = np.insert(k_mtx_m_g, 20, [[0] * 12] * 3, axis=0)

        k_mtx_m_g = np.insert(k_mtx_m_g, 2, [[0] * 24] * 3, axis=1)
        k_mtx_m_g = np.insert(k_mtx_m_g, 8, [[0] * 24] * 3, axis=1)
        k_mtx_m_g = np.insert(k_mtx_m_g, 14, [[0] * 24] * 3, axis=1)
        k_mtx_m_g = np.insert(k_mtx_m_g, 20, [[0] * 24] * 3, axis=1)

        k_mtx_p_g = np.insert(k_mtx_p, 0, [[0] * 12] * 2, axis=0)
        k_mtx_p_g = np.insert(k_mtx_p_g, 5, [[0] * 12] * 3, axis=0)
        k_mtx_p_g = np.insert(k_mtx_p_g, 11, [[0] * 12] * 3, axis=0)
        k_mtx_p_g = np.insert(k_mtx_p_g, 17, [[0] * 12] * 3, axis=0)
        k_mtx_p_g = np.insert(k_mtx_p_g, 23, [0] * 12, axis=0)

        k_mtx_p_g = np.insert(k_mtx_p_g, 0, [[0] * 24] * 2, axis=1)
        k_mtx_p_g = np.insert(k_mtx_p_g, 5, [[0] * 24] * 3, axis=1)
        k_mtx_p_g = np.insert(k_mtx_p_g, 11, [[0] * 24] * 3, axis=1)
        k_mtx_p_g = np.insert(k_mtx_p_g, 17, [[0] * 24] * 3, axis=1)
        k_mtx_p_g = np.insert(k_mtx_p_g, 23, [0] * 24, axis=1)

        self.K = k_mtx_m_g + k_mtx_p_g
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
