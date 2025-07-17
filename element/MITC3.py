#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from element.ElementBase import ElementBaseClass
from abc import ABC

from femdb.GlobalEnum import MaterialKey
from femdb.GlobalFEMVariant import ModelInfo
import numpy as np


def calculate_von_mises(stress):
    sxx, syy, sxy, sxz, syz, szz = stress
    return np.sqrt(
        (sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2 +
        6 * (sxy ** 2 + sxz ** 2 + syz ** 2)
    ) / np.sqrt(2)


class MITC3(ElementBaseClass, ABC):
    """
    Reference:
    1. https://github.com/luchete80/MITC_python
    """

    def __init__(self, eid):
        super().__init__(eid)
        self.nodes_count = 3
        self.vtu_type = "triangle"
        self.block_size = 18

        self.thickness = None
        self.material = None

        self.D_m = None
        self.D_b = None
        self.D_s = None
        self.D_drill = None

        self.B_m = None
        self.B_b = None
        self.B_s = None

        self.area = None
        self.dN_dx = None
        self.dN_dy = None

        self.x21 = None
        self.y21 = None
        self.x32 = None
        self.y32 = None
        self.x13 = None
        self.y13 = None

        self._a1 = None
        self._a2 = None
        self._a3 = None

        self.stresses = {
            'membrane': [],
            'bending': [],
            'shear': [],
            'upper': [],
            'lower': [],
            'max': [],
            'min': []
        }

        self._T_small = None
        self._T = None
        self.K = np.zeros((18, 18), dtype=float)

    def CalElementDMatrix(self, an_type=None):
        E = self.cha_dict[MaterialKey.E]
        niu = self.cha_dict[MaterialKey.Niu]
        G = self.cha_dict[MaterialKey.G]
        if self.cha_dict.__contains__('RealConst') and len(self.cha_dict["RealConst"]) != 0:
            self.thickness = self.cha_dict["RealConst"][0]
        elif self.cha_dict.__contains__(MaterialKey.Thickness):
            self.thickness = self.cha_dict[MaterialKey.Thickness]
        else:
            raise KeyError("Don't Contain RealConst and Thickness")
        k_shear = 5 / 6

        self.D_m = E / (1 - niu ** 2) * np.array([
            [1, niu, 0],
            [niu, 1, 0],
            [0, 0, (1 - niu) / 2]
        ])
        self.D_b = E * self.thickness ** 3 / (12 * (1 - niu ** 2)) * np.array([
            [1, niu, 0],
            [niu, 1, 0],
            [0, 0, (1 - niu) / 2]
        ])
        self.D_s = k_shear * G * self.thickness * np.eye(2)

        # Drilling material term (scalar)
        alpha = 0.1
        self.D_drill = alpha * E * self.thickness ** 3 / 12

    def ElementStiffness(self, from_origin=False):
        """
        计算壳单元刚度矩阵
        :param from_origin:
        :return:
        """
        """
        Compute Geometry Part
        """
        translated_nodes = self.node_coords - self.node_coords[0]  # shift node 0 to origin
        loc_nodes = translated_nodes @ self._T_small

        # Compute edge vectors
        self.x21 = loc_nodes[1][0] - loc_nodes[0][0]
        self.y21 = loc_nodes[1][1] - loc_nodes[0][1]
        self.x32 = loc_nodes[2][0] - loc_nodes[1][0]
        self.y32 = loc_nodes[2][1] - loc_nodes[1][1]
        self.x13 = loc_nodes[0][0] - loc_nodes[2][0]
        self.y13 = loc_nodes[0][1] - loc_nodes[2][1]

        self.area = self.x21 * self.y32 - self.x32 * self.y21

        b = np.array([self.y32, self.y13, self.y21])
        c = np.array([-self.x32, -self.x13, -self.x21])
        self.dN_dx = -b / self.area
        self.dN_dy = -c / self.area

        """
        Membrane Stiffness
        """
        K_mem = np.zeros((18, 18))
        B_m = np.zeros((3, 6))
        for i in range(3):
            B_m[0, 2 * i] = self.dN_dx[i]
            B_m[1, 2 * i + 1] = self.dN_dy[i]
            B_m[2, 2 * i] = self.dN_dy[i]
            B_m[2, 2 * i + 1] = self.dN_dx[i]
        self.B_m = B_m

        active_dofs_m = [0, 1, 6, 7, 12, 13]
        K_m_local = B_m.T @ self.D_m @ B_m * self.area * self.thickness * 0.5
        for i, dof_i in enumerate(active_dofs_m):
            for j, dof_j in enumerate(active_dofs_m):
                K_mem[dof_i, dof_j] += K_m_local[i, j]

        """
        Plate Stiffness
        """
        B_bend = np.zeros((3, 18))

        for i in range(3):
            # Curvatures in local basis (κ_11, κ_22, κ_12)
            B_bend[0, 6 * i + 4] = -self.dN_dx[i]  # κ_11 = -∂θ2/∂g1
            B_bend[1, 6 * i + 3] = self.dN_dy[i]  # κ_22 = ∂θ1/∂g2
            B_bend[2, 6 * i + 3] = self.dN_dx[i]  # κ_12 terms
            B_bend[2, 6 * i + 4] = -self.dN_dy[i]

        K_bend = B_bend.T @ self.D_b @ B_bend * self.area * 0.5

        self.B_b = np.zeros((3,6))
        for i in range(3):
            # Compute Cartesian derivatives
            # dN_dx = invJ[0,0]*dN_dxi[i] + invJ[0,1]*dN_deta[i]
            # dN_dy = invJ[1,0]*dN_dxi[i] + invJ[1,1]*dN_deta[i]

            # Standard curvature definitions:
            self.B_b[0, 2 * i + 1] = -self.dN_dx[i]  # κ_xx = -∂θy/∂x
            self.B_b[1, 2 * i] = self.dN_dy[i]  # κ_yy = ∂θx/∂y
            self.B_b[2, 2 * i] = self.dN_dx[i]  # κ_xy terms
            self.B_b[2, 2 * i + 1] = -self.dN_dy[i]

        """
        Shear Stiffness
        """
        k_shear = 5 / 6  # Shear correction factor
        G = self.cha_dict[MaterialKey.G]

        # Tying points in natural coordinates (ξ, η)
        tying_points = [(0.5, 0.0), (0.5, 0.5), (0.0, 0.5)]

        K_shear = np.zeros((18, 18))
        K_shear_c = np.zeros((9, 9))

        x31 = -self.x13
        y31 = -self.y13
        phi1 = np.arctan2(self.y21, self.x21)
        phi2 = 0.5 * np.pi - np.arctan2(x31, y31)

        # Create and assign BsO matrix
        BsO = np.zeros((2, 2))
        BsO[0, 0] = np.sin(phi2)
        BsO[0, 1] = -np.sin(phi1)
        BsO[1, 0] = -np.cos(phi2)
        BsO[1, 1] = np.cos(phi1)

        # Create and assign BsC matrix
        BsC = np.zeros((2, 2))

        BsC[0, 0] = np.sqrt(self.x13 ** 2 + self.y13 ** 2) / self.area
        BsC[1, 1] = np.sqrt(self.x21 ** 2 + self.y21 ** 2) / self.area

        for xi, eta in tying_points:
            cross = True
            B_s = np.zeros((2, 9))

            ###Row 0: shear strain γ_xz [-1 1 0] SO DOF 6 is not written
            ### DIFFERENCE: CROSS DERIVATIVES ARE NOT CONSTANT
            ### CROSS TERMS COUPLE 0,1 0,4 and 0,
            B_s[0, 0] = -1.0  # θ_x at node 1
            if (cross): B_s[0, 1] = -(self.y21 + self.y32 * xi) / 2.0
            B_s[0, 2] = (self.x21 + self.x32 * xi) / 2.0
            B_s[0, 3] = 1.0  # θ_x at node 2
            if (cross): B_s[0, 4] = -(self.y21 + self.y13 * xi) / 2.0
            B_s[0, 5] = (self.x21 + self.x13 * xi) / 2.0
            if (cross): B_s[0, 7] = -(self.y21 * xi) / 2.0
            B_s[0, 8] = (self.x21 * xi) / 2.0
            # #############################
            # # Row 1: shear strain γ_yz ##[-1, 0, 1] SO DOF 3 not filled
            B_s[1, 0] = -1.0  # θ_y at node 1
            B_s[1, 1] = (self.y13 + self.y32 * eta) / 2.0
            if (cross): B_s[1, 2] = -(self.x13 + self.x32 * eta) / 2.0
            B_s[1, 4] = (self.y13 * eta) / 2.0
            if (cross): B_s[1, 5] = -(self.x13 * eta) / 2.0
            B_s[1, 6] = 1.0  # θ_y at node 3
            B_s[1, 7] = (self.y13 + self.y21 * eta) / 2.0
            if (cross): B_s[1, 8] = -(self.x13 + self.x21 * eta) / 2.0
            B_s_OCT = BsO @ BsC @ B_s

            K_shear_c += (B_s_OCT.T @ B_s_OCT) * (k_shear * G * self.thickness) * (self.area / 6)

        self.B_s = B_s_OCT
        active_dofs = [2, 3, 4, 8, 9, 10, 14, 15, 16]
        for i, dof_i in enumerate(active_dofs):
            for j, dof_j in enumerate(active_dofs):
                K_shear[dof_i, dof_j] = K_shear_c[i, j]

        """
        Compute only the drilling-related stiffness components
        """
        N = np.array([1 / 3, 1 / 3, 1 / 3])

        # Initialize drilling B-matrix (1×18)
        B_drill = np.zeros((1, 18))  # Only γ_drill strain

        # Drilling strain-displacement relationship
        for i in range(3):
            B_drill[0, 6 * i] = -0.5 * self.dN_dy[i]  # ½(∂v/∂x) component
            B_drill[0, 6 * i + 1] = 0.5 * self.dN_dx[i]  # -½(∂u/∂y) component
            B_drill[0, 6 * i + 5] = -N[i]  # -θz component

        # Pure drilling stiffness matrix
        K_drill = B_drill.T @ B_drill * self.D_drill * self.area
        """Compute metrics for adaptive hourglass control."""
        # Edge lengths
        L1 = np.linalg.norm(self.node_coords[1] - self.node_coords[0])
        L2 = np.linalg.norm(self.node_coords[2] - self.node_coords[1])
        L3 = np.linalg.norm(self.node_coords[0] - self.node_coords[2])

        # Aspect ratio (max edge / min edge)
        max_edge = max(L1, L2, L3)
        min_edge = min(L1, L2, L3)
        aspect_ratio = max_edge / min_edge

        # Skewness (deviation from ideal equilateral triangle)
        ideal_area = (np.sqrt(3) / 4) * (min_edge ** 2)
        skewness = abs(self.area - ideal_area) / ideal_area
        G = self.cha_dict[MaterialKey.G]

        # Base stabilization factor (similar to ANSYS defaults)
        beta_base = 0.01  # Default for well-shaped elements

        # Increase stabilization for distorted elements
        beta = beta_base * (1 + 0.5 * (aspect_ratio - 1) + 2.0 * skewness)

        # Limit to reasonable range (0.01 to 0.1)
        beta = np.clip(beta, 0.01, 0.1)

        # Hourglass stiffness (same formulation as before)
        Bhx = np.zeros(18)
        Bhy = np.zeros(18)
        for i in range(3):
            Bhx[6 * i + 5] = self.dN_dx[i]  # ∂θz/∂x
            Bhy[6 * i + 5] = self.dN_dy[i]  # ∂θz/∂y

        K_hg = beta * G * self.thickness * self.area * (
                np.outer(Bhx, Bhx) + np.outer(Bhy, Bhy))
        K_drill = K_hg + K_drill

        K_local = K_mem + K_bend + K_shear + K_drill
        self.K = self._T @ K_local @ self._T.T
        return self.K

    def ElementMass(self):
        mass = np.zeros((18, 18))
        density = 1.0  # 可按需设定材料密度
        lumped_mass = density * self.thickness * self.area / 3
        for i in range(3):
            idx = slice(i * 6, (i + 1) * 6)
            mass[idx, idx] = np.eye(6) * lumped_mass
        return mass

    def CalculateElementStress(self, displacement):
        displacement = displacement @ self._T
        active_dofs_m = [0, 1, 6, 7, 12, 13]  # u,v DOFs for membrane part
        U_membrane = displacement[active_dofs_m]  # Condensed displacement vector (6x1)

        membrane_strains = self.B_m @ U_membrane
        membrane_stresses = self.D_m @ membrane_strains

        active_dofs_b = [3, 4, 9, 10, 15, 16]  # u,v DOFs for membrane part
        U_bending = displacement[active_dofs_b]  # Condensed displacement vector (6x1)

        bending_curvatures = self.B_b @ U_bending  # [κ_xx, κ_yy, 2κ_xy]
        bending_moments = self.D_b @ bending_curvatures  # [M_xx, M_yy, M_xy] bending moments per unit length

        """
        Convert bending moments to stresses at faces
        Bending stress formula: σ_b = (M * z) / I, where I = t³/12
        For z = ±t/2, this becomes σ_b = ± (6M)/t²
        Total stresses at faces (membrane + bending)
        """
        bending_stresses = (6.0 / self.thickness ** 2) * bending_moments[:3]  # [σ_xx_b, σ_yy_b, σ_xy_b]
        upper_face = membrane_stresses + bending_stresses
        lower_face = membrane_stresses - bending_stresses

        active_dofs_s = [2, 3, 4, 8, 9, 10, 14, 15, 16]
        U_shear = displacement[active_dofs_s]
        shear_strains = self.B_s @ U_shear
        # shear_strains_local = J_pseudo @ shear_strains_cov  [γ_xz, γ_yz]

        ##Shear stresses
        shear_forces = self.D_s @ shear_strains  # [Q_x, Q_y] hear forces per unit length, not shear stresses yet.
        shear_stresses = shear_forces / self.thickness

        # Combine into full stress vectors for each face
        # Format: [σ_xx, σ_yy, σ_xy, σ_xz, σ_yz, σ_zz]
        upper_stress = np.array([
            upper_face[0],  # σ_xx
            upper_face[1],  # σ_yy
            upper_face[2],  # σ_xy
            shear_stresses[0],  # σ_xz
            shear_stresses[1],  # σ_yz
            0.0  # σ_zz (assumed zero)
        ])

        lower_stress = np.array([
            lower_face[0],  # σ_xx
            lower_face[1],  # σ_yy
            lower_face[2],  # σ_xy
            shear_stresses[0],  # σ_xz
            shear_stresses[1],  # σ_yz
            0.0  # σ_zz (assumed zero)
        ])

        upper_von_mises = calculate_von_mises(upper_stress)
        lower_von_mises = calculate_von_mises(lower_stress)

        # 取最大值
        max_von_mises = max(upper_von_mises, lower_von_mises)
        print("max_mises:", max_von_mises)

    def CalculateBasic(self):
        """Compute and cache covariant basis vectors."""
        dN_dxi = np.array([-1, 1, 0])
        dN_deta = np.array([-1, 0, 1])

        self._a1 = np.dot(dN_dxi, self.node_coords)
        self._a2 = np.dot(dN_deta, self.node_coords)
        self._a3 = np.cross(self._a1, self._a2)
        tol = 1e-10
        e1 = self._a1 / (np.linalg.norm(self._a1) + tol)
        e2 = self._a2 - (self._a2.dot(e1)) * e1
        e2 = e2 / (np.linalg.norm(e2) + tol)
        e3 = np.cross(e1, e2)
        self._T_small = np.vstack((e1, e2, e3)).T
        """Create 18×18 transformation matrix"""
        T = np.zeros((18, 18))
        for i in range(3):  # For each node
            # Apply to translational DOFs
            T[6 * i:6 * i + 3, 6 * i:6 * i + 3] = self._T_small
            # Apply to rotational DOFs
            T[6 * i + 3:6 * i + 6, 6 * i + 3:6 * i + 6] = self._T_small
        self._T = T

    def ReCalculateElementStiffness(self):
        return self.ElementStiffness()


if __name__ == "__main__":
    mitc3 = MITC3(-1)
    mitc3.cha_dict = {MaterialKey.Niu: 0.3,
                      MaterialKey.E: 2e11,
                      MaterialKey.Thickness: 0.1,
                      MaterialKey.G: 2e11 / 2 / (1 + 0.3)
                      }
    mitc3.node_coords = np.array([
        [0, 0, 0],
        [1, 0.2, 0],
        [0.2, 1.5, 0]
    ])
    mitc3.CalculateBasic()
    mitc3.CalElementDMatrix()
    K = mitc3.ElementStiffness()

    disp = np.zeros(18)
    disp[6 * 1 + 2] = -0.00329218
    disp[6 * 1 + 3] = 0.00097242
    disp[6 * 1 + 4] = 0.00664678
    stress_output = mitc3.CalculateElementStress(disp)
    # print("Stress upper:", stress_output['upper'])
    # print("Von Mises:", stress_output['von_mises'])
