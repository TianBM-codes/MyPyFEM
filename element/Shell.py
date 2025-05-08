#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import time

from element.Plate import *
from element.Membrane import *


class DKTShell(ElementBaseClass, ABC):
    """
    DKTShell Element class
    """

    def __init__(self, eid=None):
        super().__init__(eid)
        self.nodes_count = 3  # Each element has 3 nodes
        self._nodes = [None for _ in range(self.nodes_count)]
        self.vtu_type = "triangle"
        self.e_type = 181
        self.K = np.zeros((18, 18), dtype=float)
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
        plate = KirchhoffTrianglePlate(-1)
        membrane = CSTDrill(-1)

        """
        先转换到局部坐标
        """
        T_matrix, origin = GetShellGlobal2LocalTransMatrix(self.node_coords)
        local_coord = (self.node_coords.T - origin[:, np.newaxis]).T @ T_matrix
        m_coords = local_coord[:, :2]

        """
        设置膜单元和板单元的节点坐标, 以及全局和局部坐标系的转换矩阵
        """
        membrane.node_coords = m_coords
        membrane.T_matrix = T_matrix
        plate.node_coords = m_coords
        plate.T_matrix = T_matrix

        """
        设置膜单元和版单元的材料
        """
        membrane.cha_dict = self.cha_dict
        membrane.sec_id = self.sec_id
        plate.cha_dict = self.cha_dict
        plate.sec_id = self.sec_id
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

        self.K[-1, -1] = k_mtx_m[-1, -1]

        """
        这里的T_matrix是全局==>局部，转职就是局部==>全局
        """
        R_matrix = T_matrix.T
        global_t_matrix = np.zeros((18, 18))
        global_t_matrix[0:3, 0:3] = R_matrix
        global_t_matrix[3:6, 3:6] = R_matrix
        global_t_matrix[6:9, 6:9] = R_matrix
        global_t_matrix[9:12, 9:12] = R_matrix
        global_t_matrix[12:15, 12:15] = R_matrix
        global_t_matrix[15:18, 15:18] = R_matrix

        self.K = global_t_matrix.T @ self.K @ global_t_matrix

        return self.K

    def ElementStress(self, displacement: np.array):
        """"""
        pass

    def ElementMass(self):
        """
        计算单元质量阵
        :return:
        """
        pass

    def CalculateBasic(self):
        pass


class DKQShell(ElementBaseClass, ABC):
    """
    DKQShell Element class
    """

    def __init__(self, eid=None):
        super().__init__(eid)
        self.nodes_count = 4  # Each element has 4 nodes
        self.vtu_type = "quad"
        self.e_type = 181
        self.K = np.zeros((24, 24))
        self._nodes = [None for _ in range(self.nodes_count)]
        self.unv_code = 40500
        self.local_coord = None
        self.global_t_matrix = np.zeros((24, 24))

    def CalMassMatrix(self):
        """
        计算单元的质量矩阵
        """

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
        """
        先转换到局部坐标, 初始化板单元和膜单元
        """
        plate = KirchhoffQuaPlate(-1)
        plate.sec_id = self.sec_id
        plate.node_coords = self.node_coords
        plate.cha_dict = self.cha_dict
        plate.CalElementDMatrix()
        K1 = plate.ElementStiffness()

        q4 = Q4Mem(-1)
        q4.sec_id = self.sec_id
        q4.node_coords = self.local_coord
        q4.cha_dict = self.cha_dict
        q4.CalElementDMatrix()
        K2 = q4.ElementStiffness()

        KShell = np.zeros((24, 24), dtype=float)
        index_m = [0, 1, 5, 6, 7, 11, 12, 13, 17, 18, 19, 23]
        KShell[np.ix_(index_m, index_m)] = K2

        KShell = self.global_t_matrix.T @ KShell @ self.global_t_matrix
        return K1 + KShell

    def ElementStress(self, displacement):
        """
        Calculate element stress
        """

    def ElementMass(self):
        """
        计算单元的质量阵(协调质量阵)
        :return:
        """
        """
        先转换到局部坐标
        """
        m_coords = self.local_coord[:, :2]
        self.M = np.zeros((24, 24), dtype=float)
        rho = self.cha_dict[MaterialKey.Density]
        thickness = self.cha_dict[MaterialKey.Thickness]
        sample_pt, weight = GaussIntegrationPoint.GetSamplePointAndWeight(2)
        for ri in range(2):
            for si in range(2):
                r, s = sample_pt[ri], sample_pt[si]
                N1 = 0.25 * (1 - r) * (1 - s)
                N2 = 0.25 * (1 + r) * (1 - s)
                N3 = 0.25 * (1 + r) * (1 + s)
                N4 = 0.25 * (1 - r) * (1 + s)
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

                dNdr = np.array([[0.25 * (1 + s), -0.25 * (1 + s), 0.25 * (s - 1), 0.25 * (1 - s)],
                                 [0.25 * (1 + r), 0.25 * (1 - r), 0.25 * (r - 1), -0.25 * (1 + r)]], dtype=float)
                g_weight = weight[ri] * weight[si]

                det_J = np.linalg.det(dNdr @ m_coords)
                iter_mass = H.T @ H * rho * det_J * g_weight * thickness  # 12*12
                for jj in range(4):
                    for kk in range(4):
                        row_idx = jj * 6
                        col_idx = kk * 6
                        self.M[row_idx:row_idx + 3, col_idx:col_idx + 3] += \
                            iter_mass[jj:jj + 3, kk:kk + 3]

        return self.global_t_matrix.T @ self.M @ self.global_t_matrix

    def CalculateBasic(self):
        """
        计算基本量
        :return:
        """
        T_matrix, origin = GetShellGlobal2LocalTransMatrix(self.node_coords)
        R_matrix = T_matrix.T
        self.global_t_matrix[0:3, 0:3] = R_matrix
        self.global_t_matrix[3:6, 3:6] = R_matrix
        self.global_t_matrix[6:9, 6:9] = R_matrix
        self.global_t_matrix[9:12, 9:12] = R_matrix
        self.global_t_matrix[12:15, 12:15] = R_matrix
        self.global_t_matrix[15:18, 15:18] = R_matrix
        self.global_t_matrix[18:21, 18:21] = R_matrix
        self.global_t_matrix[21:24, 21:24] = R_matrix
        self.local_coord = (self.node_coords.T - origin[:, np.newaxis]).T @ T_matrix


class CookTriShell(ElementBaseClass, ABC):
    """
    DKTShell Element class
    """

    def __init__(self, eid=None):
        super().__init__(eid)
        self.nodes_count = 3  # Each element has 3 nodes
        self._nodes = [None for _ in range(self.nodes_count)]
        self.vtu_type = "triangle"
        self.e_type = 181
        self.K = np.zeros((18, 18), dtype=float)
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
        T_matrix, origin = GetShellGlobal2LocalTransMatrix(self.node_coords)
        local_coord = (self.node_coords.T - origin[:, np.newaxis]).T @ T_matrix
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
        membrane.sec_id = self.sec_id
        plate.cha_dict = self.cha_dict
        plate.sec_id = self.sec_id
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

        self.K[-1, -1] = k_mtx_m[-1, -1]

        """
        这里的T_matrix是全局==>局部，转职就是局部==>全局
        """
        R_matrix = T_matrix.T
        global_t_matrix = np.zeros((18, 18))
        global_t_matrix[0:3, 0:3] = R_matrix
        global_t_matrix[3:6, 3:6] = R_matrix
        global_t_matrix[6:9, 6:9] = R_matrix
        global_t_matrix[9:12, 9:12] = R_matrix
        global_t_matrix[12:15, 12:15] = R_matrix
        global_t_matrix[15:18, 15:18] = R_matrix

        self.K = global_t_matrix.T @ self.K @ global_t_matrix

        return self.K

    def ElementStress(self, displacement):
        """
        Calculate element stress
        """

    def ElementMass(self):
        pass

    def CalculateBasic(self):
        pass


class CookQuaShell(ElementBaseClass, ABC):
    """
    CookShell Element class
    """

    def __init__(self, eid=None):
        super().__init__(eid)
        self.nodes_count = 4  # Each element has 4 nodes
        self.vtu_type = "quad"
        self.e_type = 181
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
        T_matrix, origin = GetShellGlobal2LocalTransMatrix(self.node_coords)
        local_coord = (self.node_coords.T - origin[:, np.newaxis]).T @ T_matrix
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
        membrane.sec_id = self.sec_id
        plate.cha_dict = self.cha_dict
        plate.sec_id = self.sec_id
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

        self.K[-1, -1] = k_mtx_m[-1, -1]

        """
        这里的T_matrix是全局==>局部，转职就是局部==>全局
        """
        R_matrix = T_matrix.T
        global_t_matrix = np.zeros((24, 24))
        global_t_matrix[0:3, 0:3] = R_matrix
        global_t_matrix[3:6, 3:6] = R_matrix
        global_t_matrix[6:9, 6:9] = R_matrix
        global_t_matrix[9:12, 9:12] = R_matrix
        global_t_matrix[12:15, 12:15] = R_matrix
        global_t_matrix[15:18, 15:18] = R_matrix
        global_t_matrix[18:21, 18:21] = R_matrix
        global_t_matrix[21:24, 21:24] = R_matrix

        self.K = global_t_matrix.T @ self.K @ global_t_matrix

        return self.K

    def ElementStress(self, displacement):
        """
        Calculate element stress
        """

    def ElementMass(self):
        pass

    def CalculateBasic(self):
        pass


if __name__ == "__main__":
    """
    测试四边形壳单元刚度阵
    """
    t_ele = DKQShell()
    t_ele.sec_id = 10001
    t_ele.cha_dict = {MaterialKey.Niu: 0.3, MaterialKey.E: 2e11, 10001: 0.01}
    t_ele.node_coords = np.array([
        [0, 0, 0],
        [1, 0, 0],
        [0.8, 0.8, 0],
        [0.5, 1, 0]
    ], dtype=float)
    t_ele.CalculateBasic()
    t_ele.CalElementDMatrix()
    Ke = t_ele.ElementStiffness()
    # np.savetxt('ke',Ke)
    # print(Ke)
    """
    测试四边形壳单元的质量阵
    """
    time1 = time.time()
    t_ele = DKQShell()
    t_ele.cha_dict = {MaterialKey.Niu: 0.3,
                      MaterialKey.E: 2e11,
                      MaterialKey.Thickness: 0.1,
                      MaterialKey.Density: 7850}
    t_ele.node_coords = np.array([
        [-1.0, -1.0, 0.0],  # 节点1
        [1.0, -1.0, 0.0],  # 节点2
        [1.0, 1.0, 0.0],  # 节点3
        [-1.0, 1.0, 0.0]  # 节点4
    ], dtype=float)
    t_ele.CalculateBasic()
    M = t_ele.ElementMass()
    print("calculated mass:", np.sum(np.diag(M)))
    print(f"theory mass:{4 * t_ele.cha_dict[MaterialKey.Density] * t_ele.cha_dict[MaterialKey.Thickness]}")
    time2 = time.time()
    print(" {:<.5f} seconds".format(time2-time1))

    """
    测试不同壳单元的刚度阵为什么差这么多
    """
    t_ele = CookQuaShell()
    t_ele.sec_id = 10001
    t_ele.cha_dict = {MaterialKey.Niu: 0.3, MaterialKey.E: 2e11, 10001: 0.01}
    t_ele.node_coords = np.array([
        [0, 0, 0],
        [1, 0, 0],
        [0.8, 0.8, 0],
        [0.5, 1, 0]
    ], dtype=float)
    t_ele.CalculateBasic()
    t_ele.CalElementDMatrix()
    Ke2 = t_ele.ElementStiffness()

    print("finish")
