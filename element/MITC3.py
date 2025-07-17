import numpy as np
from enum import Enum


class MITC3Element:
    """
    MITC3三角形壳单元类 (仅返回Mises应力最大值)
    """

    def __init__(self, eid=None):
        self.id = eid
        self.nodes_count = 3
        self.node_dof_count = 6
        self.E = None
        self.nu = None
        self.thickness = None
        self.node_coords = None

        # 几何量
        self._T_small = None
        self._T = None
        self.dN_dx = None
        self.dN_dy = None
        self.A2 = None

    def CalculateBasic(self):
        """计算基本几何量和转换矩阵"""
        dN_dxi = np.array([-1, 1, 0])
        dN_deta = np.array([-1, 0, 1])
        self._a1 = np.dot(dN_dxi, self.node_coords)
        self._a2 = np.dot(dN_deta, self.node_coords)
        self._a3 = np.cross(self._a1, self._a2)

        # 局部坐标系
        tol = 1e-10
        e1 = self._a1 / (np.linalg.norm(self._a1) + tol)
        e2 = self._a2 - (self._a2.dot(e1)) * e1
        e2 = e2 / (np.linalg.norm(e2) + tol)
        e3 = np.cross(e1, e2)
        self._T_small = np.vstack((e1, e2, e3)).T

        # 18×18转换矩阵
        self._T = np.zeros((18, 18))
        for i in range(3):
            self._T[6 * i:6 * i + 3, 6 * i:6 * i + 3] = self._T_small
            self._T[6 * i + 3:6 * i + 6, 6 * i + 3:6 * i + 6] = self._T_small

        # 单元几何参数
        translated_nodes = self.node_coords - self.node_coords[0]
        loc_nodes = translated_nodes @ self._T_small

        self.x21 = loc_nodes[1][0] - loc_nodes[0][0]
        self.y21 = loc_nodes[1][1] - loc_nodes[0][1]
        self.x32 = loc_nodes[2][0] - loc_nodes[1][0]
        self.y32 = loc_nodes[2][1] - loc_nodes[1][1]
        self.x13 = loc_nodes[0][0] - loc_nodes[2][0]
        self.y13 = loc_nodes[0][1] - loc_nodes[2][1]

        self.A2 = self.x21 * self.y32 - self.x32 * self.y21
        b = np.array([self.y32, self.y13, self.y21])
        c = np.array([-self.x32, -self.x13, -self.x21])
        self.dN_dx = -b / self.A2
        self.dN_dy = -c / self.A2

    def CalElementDMatrix(self):
        """计算本构矩阵"""
        # 薄膜本构矩阵
        factor = self.E / (1 - self.nu ** 2)
        self.D_m = factor * np.array([
            [1, self.nu, 0],
            [self.nu, 1, 0],
            [0, 0, (1 - self.nu) / 2]
        ])

        # 弯曲本构矩阵
        factor = self.E * self.thickness ** 3 / (12 * (1 - self.nu ** 2))
        self.D_b = factor * np.array([
            [1, self.nu, 0],
            [self.nu, 1, 0],
            [0, 0, (1 - self.nu) / 2]
        ])

        # 剪切本构矩阵
        k_shear = 5 / 6
        G = self.E / (2 * (1 + self.nu))
        self.D_s = k_shear * G * self.thickness * np.eye(2)

    def _calculate_von_mises(self, stress):
        """计算Mises应力"""
        sxx, syy, sxy, sxz, syz, szz = stress
        return np.sqrt(
            (sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2 +
            6 * (sxy ** 2 + sxz ** 2 + syz ** 2)
        ) / np.sqrt(2)

    def CalculateElementStress(self, displacement):
        """
        计算单元应力并返回上下表面Mises应力最大值
        返回: max_mises (上下面Mises应力最大值)
        """
        # 转换位移到局部坐标系
        U_local = self._T @ displacement

        # 薄膜应力
        active_dofs_m = [0, 1, 6, 7, 12, 13]
        U_membrane = U_local[active_dofs_m]
        membrane_strains = np.zeros(3)
        for i in range(3):
            membrane_strains[0] += self.dN_dx[i] * U_membrane[2 * i]
            membrane_strains[1] += self.dN_dy[i] * U_membrane[2 * i + 1]
            membrane_strains[2] += self.dN_dy[i] * U_membrane[2 * i] + self.dN_dx[i] * U_membrane[2 * i + 1]
        membrane_stresses = self.D_m @ membrane_strains

        # 弯曲应力
        active_dofs_b = [3, 4, 9, 10, 15, 16]
        U_bending = U_local[active_dofs_b]
        bending_curvatures = np.zeros(3)
        for i in range(3):
            bending_curvatures[0] -= self.dN_dx[i] * U_bending[2 * i + 1]
            bending_curvatures[1] += self.dN_dy[i] * U_bending[2 * i]
            bending_curvatures[2] += self.dN_dx[i] * U_bending[2 * i] - self.dN_dy[i] * U_bending[2 * i + 1]
        bending_moments = self.D_b @ bending_curvatures
        bending_stresses = (6.0 / self.thickness ** 2) * bending_moments[:3]

        # 剪切应力
        active_dofs_s = [2, 3, 4, 8, 9, 10, 14, 15, 16]
        U_shear = U_local[active_dofs_s]
        shear_strains = np.zeros(2)
        for i in range(3):
            shear_strains[0] += (-1.0 if i == 0 else 1.0 if i == 1 else 0.0) * U_shear[3 * i]
            shear_strains[1] += (-1.0 if i == 0 else 0.0 if i == 1 else 1.0) * U_shear[3 * i + 1]
        shear_stresses = self.D_s @ shear_strains / self.thickness

        # 组合上下表面应力
        upper_stress = np.array([
            membrane_stresses[0] + bending_stresses[0],  # σxx
            membrane_stresses[1] + bending_stresses[1],  # σyy
            membrane_stresses[2] + bending_stresses[2],  # σxy
            shear_stresses[0],  # σxz
            shear_stresses[1],  # σyz
            0.0  # σzz
        ])

        lower_stress = np.array([
            membrane_stresses[0] - bending_stresses[0],
            membrane_stresses[1] - bending_stresses[1],
            membrane_stresses[2] - bending_stresses[2],
            shear_stresses[0],
            shear_stresses[1],
            0.0
        ])

        # 计算Mises应力
        upper_mises = self._calculate_von_mises(upper_stress)
        lower_mises = self._calculate_von_mises(lower_stress)

        # 返回最大值
        return max(upper_mises, lower_mises)


if __name__ == "__main__":
    """测试算例：悬臂三角形壳单元"""
    # 创建单元
    elem = MITC3Element(eid=1)

    # 材料属性 (钢)
    elem.E = 2.0e11  # Pa
    elem.nu = 0.3
    elem.thickness = 0.1  # m

    # 节点坐标 (直角三角形)
    elem.node_coords = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.2, 0.0],
        [0.2, 1.5, 0.0]
    ])

    # 计算基本几何量
    elem.CalculateBasic()
    elem.CalElementDMatrix()

    # 设置位移 (自由端受载情况)
    disp = np.zeros(18)
    # disp[6 * 1 + 2] = -0.01  # 节点2 z向位移 -10mm
    # disp[6 * 2 + 2] = -0.01  # 节点3 z向位移 -10mm
    disp[6 * 1 + 2] = -0.00329218
    disp[6 * 1 + 3] = 0.00097242
    disp[6 * 1 + 4] = 0.00664678

    # 计算应力
    max_mises = elem.CalculateElementStress(disp)
    print(f"最大Mises应力: {max_mises / 1e6:.2f} MPa")
