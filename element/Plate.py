#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from element.ElementBase import *
from abc import ABC


# TODO: 中厚板 MITC4和MITC3

class MITC4(ElementBaseClass, ABC):
    """
    MITC4 Element class
    Reference:
    1.《有限元法、理论、格式与求解方法》上册Bathe P395
    2. 王欢Matlab程序
    """

    def __init__(self, eid=None):
        super().__init__(eid)
        self.nodes_count = 4  # Each element has 8 nodes
        self.K = np.zeros([8, 8], dtype=float)  # 刚度矩阵
        self.vtp_type = "quad"
        self.thickness = None

    def CalElementDMatrix(self, an_type=None):
        """
        计算本构矩阵, 弹性模量和泊松比, Bathe 上册P184
        """
        e = self.cha_dict[MaterialKey.E]
        niu = self.cha_dict[MaterialKey.Niu]
        if an_type == MaterialMatrixType.PlaneStree or an_type is None:
            a = e / (1 - niu ** 2)
            self.D = a * np.array([[1, niu, 0],
                                   [niu, 1, 0],
                                   [0, 0, 0.5 * (1 - niu)]], dtype=float)
        elif an_type == MaterialMatrixType.PlaneStrain:
            a = e * (1 - niu) / (1 + niu) / (1 - 2 * niu)
            self.D = a * np.array([[1, niu / (1 - niu), 0],
                                   [niu(1 - niu), 1, 0],
                                   [0, 0, 0.5 * (1 - 2 * niu) / (1 - niu)]], dtype=float)
        else:
            mlogger.fatal("Unknown an_dimension")
            sys.exit(1)

    def ElementStiffness(self):
        """
        dimension: 4*3, [[x1,y1,z1],[x2,y2,z2],...[x4,y4,z4]], type:np.ndarray, dtype:float

        # Shape Function:
        N1 = 0.25 * (1 + r) * (1 + s)
        N2 = 0.25 * (1 - r) * (1 + s)
        N3 = 0.25 * (1 - r) * (1 - s)
        N4 = 0.25 * (1 + r) * (1 - s)

        # Partial
        """
        assert self.node_coords.shape == (4, 3)

        # Gaussian Weight
        sample_pt, weight = GaussIntegrationPoint.GetSamplePointAndWeight(2)

        """
        Bending Part
        """

        """
        Shear Stain Part
        """
        tt = np.array([[1, -1, -1, 1], [1, -1, 1, -1], [1, 1, -1, -1]], dtype=float)
        Ax, Bx, Cx = tt * self.node_coords[:, 0]
        Ay, By, Cy = tt * self.node_coords[:, 1]
        r, s = 0, 0

        # Partial Derivative
        pN1pr, pN1ps = 0.25 * (1 + s), 0.25 * (1 + r)
        pN2pr, pN2ps = -0.25 * (1 + s), 0.25 * (1 - r)
        pN3pr, pN3ps = 0.25 * (s - 1), 0.25 * (r - 1)
        pN4pr, pN4ps = 0.25 * (1 - s), -0.25 * (1 + r)
        pXpr = np.matmul(np.array([pN1pr, pN2pr, pN3pr, pN4pr]), self.node_coords[:, 0])
        pYpr = np.matmul(np.array([pN1pr, pN2pr, pN3pr, pN4pr]), self.node_coords[:, 1])
        pXps = np.matmul(np.array([pN1ps, pN2ps, pN3ps, pN4ps]), self.node_coords[:, 0])
        pYps = np.matmul(np.array([pN1ps, pN2ps, pN3ps, pN4ps]), self.node_coords[:, 1])
        alpha = np.arctan(pYpr / pXpr)
        beta = np.arctan(pYps / pXps)
        detJ = pXpr * pYps - pYpr * pXps

        coeff_gama_rz = 0.125 * np.sqrt((Cx + r * Bx) ** 2 + (Cy + r * By) ** 2) / detJ
        coeff_gama_sz = 0.125 * np.sqrt((Ax + s * Bx) ** 2 + (Ay + s * By) ** 2) / detJ

        gama_rz = np.zeros((12, 1), dtype=float)
        # w1, w2, w3, w4
        gama_rz[0, 0] = 0.5 * coeff_gama_rz * (1 + s)
        gama_rz[3, 0] = -gama_rz[0, 0]
        gama_rz[6, 0] = -0.5 * coeff_gama_rz * (1 - s)
        gama_rz[9, 0] = -gama_rz[6, 0]

        # theta_1x, theta_2x, theta_3x, theta_4x
        gama_rz[1, 0] = -gama_rz[0, 0] * 0.5 * (self.node_coords[0, 1] - self.node_coords[1, 1])
        gama_rz[4, 0] = gama_rz[1, 0]
        gama_rz[7, 0] = 0.5 * gama_rz[6, 0] * (self.node_coords[2, 1] - self.node_coords[3, 1])
        gama_rz[10, 0] = gama_rz[7, 0]

        # theta_1y, theta_2y, theta_3y, theta_4y
        gama_rz[2, 0] = gama_rz[0, 0] * 0.5 * (self.node_coords[0, 0] - self.node_coords[1, 0])
        gama_rz[5, 0] = gama_rz[2, 0]
        gama_rz[8, 0] = -gama_rz[6, 0] * 0.5 * (self.node_coords[3, 0] - self.node_coords[2, 0])
        gama_rz[11, 0] = gama_rz[8, 0]

        gama_sz = np.zeros((12, 1), dtype=float)
        # w1, w2, w3, w4
        gama_sz[0, 0] = 0.5 * coeff_gama_sz * (1 + r)
        gama_sz[3, 0] = 0.5 * coeff_gama_sz * (1 - r)
        gama_sz[6, 0] = -gama_sz[3, 0]
        gama_sz[9, 0] = -gama_sz[0, 0]

        # theta_1x, theta_2x, theta_3x, theta_4x
        gama_sz[1, 0] = -gama_sz[0, 0] * 0.5 * (self.node_coords[0, 1] - self.node_coords[3, 1])
        gama_sz[4, 0] = -0.5 * gama_sz[3, 0] * (self.node_coords[1, 1] - self.node_coords[2, 1])
        gama_sz[7, 0] = gama_sz[4, 0]
        gama_sz[10, 0] = gama_sz[1, 0]

        # theta_1y, theta_2y, theta_3y, theta_4y
        gama_sz[2, 0] = gama_sz[0, 0] * 0.5 * (self.node_coords[0, 0] - self.node_coords[3, 0])
        gama_sz[5, 0] = gama_sz[1, 0] * 0.5 * (self.node_coords[1, 0] - self.node_coords[2, 0])
        gama_sz[8, 0] = gama_sz[5, 0]
        gama_sz[11, 0] = gama_sz[2, 0]

        # local global_coord to cartesian global_coord
        gama_xz = gama_rz * np.sin(beta) - gama_sz * np.sin(alpha)
        gama_yz = -gama_rz * np.cos(beta) + gama_sz * np.cos(alpha)

        # 在4个高斯点上积分
        for ri in range(2):
            for si in range(2):
                r, s = sample_pt[ri], sample_pt[si]
                g_weight = weight[ri] * weight[si]

                # self.K += g_weight * self.B.T * self.D * B * det_J * self.cha_dict[PropertyKey.ThicknessOrArea]

        return self.K

    def ElementStress(self, displacement):
        """
        Calculate element stress
        """


class MITC3(ElementBaseClass, ABC):
    """ plate 3node Element class """

    def __init__(self, eid=None):
        super().__init__(eid)
        self.nodes_count = 3  # Each element has 3 nodes
        self.K = np.zeros([6, 6], dtype=float)  # 刚度矩阵
        self.vtp_type = "triangle"
        self.thickness = None

    def CalElementDMatrix(self, an_type=None):
        """
        计算本构矩阵, 弹性模量和泊松比, Bathe 上册P184
        """
        e = self.cha_dict[MaterialKey.E]
        niu = self.cha_dict[MaterialKey.Niu]
        if an_type == MaterialMatrixType.PlaneStree or an_type is None:
            a = e / (1 - niu ** 2)
            self.D = a * np.array([[1, niu, 0],
                                   [niu, 1, 0],
                                   [0, 0, 0.5 * (1 - niu)]], dtype=float)
        elif an_type == MaterialMatrixType.PlaneStrain:
            a = e * (1 - niu) / (1 + niu) / (1 - 2 * niu)
            self.D = a * np.array([[1, niu / (1 - niu), 0],
                                   [niu(1 - niu), 1, 0],
                                   [0, 0, 0.5 * (1 - 2 * niu) / (1 - niu)]], dtype=float)
        else:
            mlogger.fatal("Unknown an_dimension")
            sys.exit(1)

    def ElementStiffness(self):
        """
        TODO: 积分过程是否正确?
        Bathe 上册 P349, 转化到参数坐标下的面积积分后, 在积分域内为常数, 所以积分等于面积 0.5
        dimension: 2*2, [[x1,y1],[x2,y2]], type:np.ndarray, dtype:float

        # Shape Function:
        N1 = 1 - r - s
        N2 = r
        N2 = s

        # Partial
        dN1dr, dN1ds = -1, -1
        dN2dr, dN2ds =  1,  0
        dN3dr, dN3ds =  0,  1
        """
        assert self.node_coords.shape == (3, 2)

        dNdr = np.array([[-1, 1, 0],
                         [-1, 0, 1]], dtype=float)

        # Jacobi 2*2 & B Matrix 3*8
        J = np.matmul(dNdr, self.node_coords)
        det_J = np.linalg.det(J)
        J_inv = np.linalg.inv(J)
        B_pre = np.matmul(J_inv, dNdr)
        B = np.array([[B_pre[0, 0], 0, B_pre[0, 1], 0, B_pre[0, 2], 0],
                      [0, B_pre[1, 0], 0, B_pre[1, 1], 0, B_pre[1, 2]],
                      [B_pre[1, 0], B_pre[0, 0], B_pre[1, 1], B_pre[0, 1], B_pre[1, 2], B_pre[0, 2]]], dtype=float)

        return B.T * self.D * B * det_J * 0.5 * self.cha_dict[PropertyKey.ThicknessOrArea]

    def ElementStress(self, displacement):
        """
        Calculate element stress
        """


class DKTPlate(ElementBaseClass, ABC):
    """
    DKT plate 3node Element class
    Reference:
    1. 有限单元法  王勖成  P364
    2. A_Study_of_Three-Node_Triangular_Plate_Bending_Elements.pdf
    3. Note And Explanation of Formulation for DKT Shell Element and its Implementation
       through Open Source Finite Software,Elmer.pdf
    4. Code_Aster Elements of plate: modelings DKT, DST, DKTG and Q4G
    TODO: 三角形积分, 面积坐标和r,s参数坐标的区别, 面积坐标见王勖成P366
    TODO: 应力杂交单元
    """

    def __init__(self, eid=None):
        super().__init__(eid)
        self.nodes_count = 3  # Each element has 3 nodes
        self.K = np.zeros([9, 9], dtype=float)  # 刚度矩阵
        self.vtp_type = "triangle"
        self.thickness = None
        self.T_matrix = None  # 整体坐标转到局部坐标的矩阵, 是转换几何坐标的

    def CalElementDMatrix(self, an_type=None):
        """
        计算本构矩阵, 弹性模量和泊松比, Bathe 上册P184
        """
        e = self.cha_dict[MaterialKey.E]
        niu = self.cha_dict[MaterialKey.Niu]
        if an_type == MaterialMatrixType.PlaneStree or an_type is None:
            a = e / (1 - niu ** 2)
            self.D = a * np.array([[1, niu, 0],
                                   [niu, 1, 0],
                                   [0, 0, 0.5 * (1 - niu)]], dtype=float)
        elif an_type == MaterialMatrixType.PlaneStrain:
            a = e * (1 - niu) / (1 + niu) / (1 - 2 * niu)
            self.D = a * np.array([[1, niu / (1 - niu), 0],
                                   [niu(1 - niu), 1, 0],
                                   [0, 0, 0.5 * (1 - 2 * niu) / (1 - niu)]], dtype=float)
        else:
            mlogger.fatal("Unknown an_dimension")
            sys.exit(1)

    def ElementStiffness(self):
        """
        Bathe 上册 P349, 转化到参数坐标下的面积积分后, 在积分域内为常数, 所以积分等于面积 0.5
        dimension: 2*2, [[x1,y1],[x2,y2]], type:np.ndarray, dtype:float
        单元的坐标是通过Shell单元给的, 不是在读取cdb文件时候给的
        """
        assert self.node_coords.shape == (3, 2)
        self.thickness = self.cha_dict["RealConst"][0]  # TODO 暂时支持各个点的厚度是一样的情形

        # 开始论文中的计算
        x12 = self.node_coords[0, 0] - self.node_coords[1, 0]
        x23 = self.node_coords[1, 0] - self.node_coords[2, 0]
        x31 = self.node_coords[2, 0] - self.node_coords[0, 0]

        y12 = self.node_coords[0, 1] - self.node_coords[1, 1]
        y23 = self.node_coords[1, 1] - self.node_coords[2, 1]
        y31 = self.node_coords[2, 1] - self.node_coords[0, 1]

        L4_square = x12 ** 2 + y12 ** 2
        L5_square = x23 ** 2 + y23 ** 2
        L6_square = x31 ** 2 + y31 ** 2

        a4 = - x12 / L4_square
        a5 = - x23 / L5_square
        a6 = - x31 / L6_square

        b4 = 0.75 * x12 * y12 / L4_square
        b5 = 0.75 * x23 * y23 / L5_square
        b6 = 0.75 * x31 * y31 / L6_square

        c4 = (0.25 * x12 ** 2 - 0.5 * y12 ** 2) / L4_square
        c5 = (0.25 * x23 ** 2 - 0.5 * y23 ** 2) / L5_square
        c6 = (0.25 * x31 ** 2 - 0.5 * y31 ** 2) / L6_square

        d4 = - y12 / L4_square
        d5 = - y23 / L5_square
        d6 = - y31 / L6_square

        e4 = (0.25 * y12 ** 2 - 0.5 * x12 ** 2) / L4_square
        e5 = (0.25 * y23 ** 2 - 0.5 * x23 ** 2) / L5_square
        e6 = (0.25 * y31 ** 2 - 0.5 * x31 ** 2) / L6_square

        sample_pt, weight = GaussIntegrationPoint.GetTrianglePointAndWeight(3)

        # 在3个高斯点上积分
        for ii in range(len(sample_pt)):
            r, s = sample_pt[ii]

            pN1pr, pN2pr, pN3pr = -3 + 4 * (r + s), -1 + 4 * r, 0
            pN4pr, pN5pr, pN6pr = 4 * (1 - 2 * r - s), 4 * s, -4 * s
            pN1ps, pN2ps, pN3ps = pN1pr, 0, -1 + 4 * s
            pN4ps, pN5ps, pN6ps = -4 * r * (1 - r), 4 * s, r * (1 - r - 2 * s)

            pHxpr = np.asarray([1.5 * (a4 * pN4pr - a6 * pN6pr), b4 * pN4pr + b6 * pN6pr, pN1pr - c4 * pN4pr - c6 * pN6pr,
                                1.5 * (a5 * pN5pr - a4 * pN4pr), b5 * pN5pr + b4 * pN4pr, pN2pr - c5 * pN5pr - c4 * pN4pr,
                                1.5 * (a6 * pN6pr - a5 * pN5pr), b6 * pN6pr + b5 * pN5pr, pN3pr - c6 * pN6pr - c5 * pN5pr], dtype=float)

            pHxps = np.asarray([1.5 * (a4 * pN4ps - a6 * pN6ps), b4 * pN4ps + b6 * pN6ps, pN1ps - c4 * pN4ps - c6 * pN6ps,
                                1.5 * (a5 * pN5ps - a4 * pN4ps), b5 * pN5ps + b4 * pN4ps, pN2ps - c5 * pN5ps - c4 * pN4ps,
                                1.5 * (a6 * pN6ps - a5 * pN5ps), b6 * pN6ps + b5 * pN5ps, pN3ps - c6 * pN6ps - c5 * pN5ps], dtype=float)

            pHypr = np.asarray([1.5 * (d4 * pN4pr - d6 * pN6pr), -pN1pr + e4 * pN4pr + e6 * pN6pr, -b4 * pN4pr - b6 * pN6pr,
                                1.5 * (d5 * pN5pr - d4 * pN4pr), -pN2pr + e5 * pN5pr + e4 * pN4pr, -b5 * pN5pr - b4 * pN4pr,
                                1.5 * (d6 * pN6pr - d5 * pN5pr), -pN3pr + e6 * pN6pr + e5 * pN5pr, -b6 * pN6pr - b5 * pN5pr], dtype=float)

            pHyps = np.asarray([1.5 * (d4 * pN4ps - d6 * pN6ps), -pN1ps + e4 * pN4ps + e6 * pN6ps, -pHxps[1],
                                1.5 * (d5 * pN5ps - d4 * pN4ps), -pN2ps + e5 * pN5ps + e4 * pN4ps, -pHxps[4],
                                1.5 * (d6 * pN6ps - d5 * pN5ps), -pN3ps + e6 * pN6ps + e5 * pN5ps, -pHxps[7]], dtype=float)

            # Jacobi 2*2
            detJ = x31 * y12 - x12 * y31

            j11 = y31 / detJ
            j12 = y12 / detJ
            j21 = -x31 / detJ
            j22 = -x12 / detJ

            # B Matrix
            B = np.asarray([j11 * pHxps + j12 * pHxps,
                            j21 * pHyps + j22 * pHypr,
                            j11 * pHyps + j12 * pHypr + j21 * pHxps + j22 * pHxpr], dtype=float)

            self.K += np.matmul(np.matmul(B.T, self.D), B) * weight[ii] * detJ

        # self.K = np.matmul(np.matmul(self.T_matrix.T, self.K), self.T_matrix) * self.thickness
        return self.K * self.thickness

    def ElementStress(self, displacement):
        """
        Calculate element stress
        """


class DKQPlate(ElementBaseClass, ABC):
    """
    DKQ plate 4node Element class
    Reference:
    1. Evaluation of a new quadrilateral thin plate bending element.pdf  JEAN-LOUIS BATOZ
    """

    def __init__(self, eid=None):
        super().__init__(eid)
        self.nodes_count = 4  # Each element has 3 nodes
        self.K = np.zeros([12, 12], dtype=float)  # 刚度矩阵
        self.vtp_type = "quad"
        self.thickness = None
        self.T_matrix = None  # 整体坐标转到局部坐标的矩阵, 是转换几何坐标的

    def CalElementDMatrix(self, an_type=None):
        """
        计算本构矩阵, 弹性模量和泊松比, Bathe 上册P184
        """
        e = self.cha_dict[MaterialKey.E]
        niu = self.cha_dict[MaterialKey.Niu]
        if an_type == MaterialMatrixType.PlaneStree or an_type is None:
            a = e / (1 - niu ** 2)
            self.D = a * np.array([[1, niu, 0],
                                   [niu, 1, 0],
                                   [0, 0, 0.5 * (1 - niu)]], dtype=float)
        elif an_type == MaterialMatrixType.PlaneStrain:
            a = e * (1 - niu) / (1 + niu) / (1 - 2 * niu)
            self.D = a * np.array([[1, niu / (1 - niu), 0],
                                   [niu(1 - niu), 1, 0],
                                   [0, 0, 0.5 * (1 - 2 * niu) / (1 - niu)]], dtype=float)
        else:
            mlogger.fatal("Unknown an_dimension")
            sys.exit(1)

    def ElementStiffness(self):
        """
        形函数与膜单元(CPM8)类似, 为8节点四边形单元
        """
        assert self.node_coords.shape == (4, 2)
        self.thickness = self.cha_dict["RealConst"][0]

        # 开始论文中的计算
        x12 = self.node_coords[0, 0] - self.node_coords[1, 0]
        x23 = self.node_coords[1, 0] - self.node_coords[2, 0]
        x34 = self.node_coords[2, 0] - self.node_coords[3, 0]
        x41 = self.node_coords[3, 0] - self.node_coords[0, 0]
        x13 = self.node_coords[0, 0] - self.node_coords[2, 0]
        x24 = self.node_coords[1, 0] - self.node_coords[3, 0]

        y12 = self.node_coords[0, 1] - self.node_coords[1, 1]
        y23 = self.node_coords[1, 1] - self.node_coords[2, 1]
        y34 = self.node_coords[2, 1] - self.node_coords[3, 1]
        y41 = self.node_coords[3, 1] - self.node_coords[0, 1]
        y13 = self.node_coords[0, 1] - self.node_coords[2, 1]
        y24 = self.node_coords[1, 1] - self.node_coords[3, 1]

        L5_square = x12 ** 2 + y12 ** 2
        L6_square = x23 ** 2 + y23 ** 2
        L7_square = x34 ** 2 + y34 ** 2
        L8_square = x41 ** 2 + y41 ** 2

        a5 = - x12 / L5_square
        a6 = - x23 / L6_square
        a7 = - x34 / L7_square
        a8 = - x41 / L8_square

        b5 = 0.75 * x12 * y12 / L5_square
        b6 = 0.75 * x23 * y23 / L6_square
        b7 = 0.75 * x34 * y34 / L7_square
        b8 = 0.75 * x41 * y41 / L8_square

        c5 = (0.25 * x12 ** 2 - 0.5 * y12 ** 2) / L5_square
        c6 = (0.25 * x23 ** 2 - 0.5 * y23 ** 2) / L6_square
        c7 = (0.25 * x34 ** 2 - 0.5 * y34 ** 2) / L7_square
        c8 = (0.25 * x41 ** 2 - 0.5 * y41 ** 2) / L8_square

        d5 = - y12 / L5_square
        d6 = - y23 / L6_square
        d7 = - y34 / L7_square
        d8 = - y41 / L8_square

        e5 = (0.25 * y12 ** 2 - 0.5 * x12 ** 2) / L5_square
        e6 = (0.25 * y23 ** 2 - 0.5 * x23 ** 2) / L6_square
        e7 = (0.25 * y34 ** 2 - 0.5 * x34 ** 2) / L7_square
        e8 = (0.25 * y41 ** 2 - 0.5 * x41 ** 2) / L8_square

        sample_pt, weight = GaussIntegrationPoint.GetSamplePointAndWeight(2)

        # 在4个高斯点上积分
        for ri in range(2):
            for si in range(2):
                r, s = sample_pt[ri], sample_pt[si]
                pN1pr = 0.25 * (s ** 2 + s) + 0.5 * (1 + s) * r
                pN2pr = -0.25 * (s ** 2 + s) + 0.5 * (1 + s) * r
                pN3pr = 0.25 * (s - s ** 2) + 0.5 * r * (1 - s)
                pN4pr = 0.25 * (s ** 2 - s) + 0.5 * r*(1 - s)
                pN5pr = -r * (1 + s)
                pN6pr = 0.5 * (s ** 2 - 1)
                pN7pr = r * (s - 1)
                pN8pr = 0.5 * (1 - s ** 2)

                pN1ps = 0.25 * (r ** 2 + r) + 0.5 * s * (1 + r)
                pN2ps = 0.25 * (r ** 2 - r) + 0.5 * s * (1 - r)
                pN3ps = 0.25 * (r - r ** 2) + 0.5 * s * (1 - r)
                pN4ps = -0.25 * (r ** 2 + r) + 0.5 * s * (1 + r)
                pN5ps = 0.5 * (1 - r ** 2)
                pN6ps = s * (r - 1)
                pN7ps = 0.5 * (r ** 2 - 1)
                pN8ps = -s * (1 + r)

                pHxpr = np.asarray([1.5 * (a5 * pN5pr - a8 * pN8pr), b5 * pN5pr + b8 * pN8pr, pN1pr - c5 * pN5pr - c8 * pN8pr,
                                    1.5 * (a6 * pN6pr - a5 * pN5pr), b6 * pN6pr + b5 * pN5pr, pN2pr - c6 * pN6pr - c5 * pN5pr,
                                    1.5 * (a7 * pN7pr - a6 * pN6pr), b7 * pN7pr + b6 * pN6pr, pN3pr - c7 * pN7pr - c6 * pN6pr,
                                    1.5 * (a8 * pN8pr - a7 * pN7pr), b8 * pN8pr + b7 * pN7pr, pN4pr - c8 * pN8pr - c7 * pN7pr], dtype=float)

                pHxps = np.asarray([1.5 * (a5 * pN5ps - a8 * pN8ps), b5 * pN5ps + b8 * pN8ps, pN1ps - c5 * pN5ps - c8 * pN8ps,
                                    1.5 * (a6 * pN6ps - a5 * pN5ps), b6 * pN6ps + b5 * pN5ps, pN2ps - c6 * pN6ps - c5 * pN5ps,
                                    1.5 * (a7 * pN7ps - a6 * pN6ps), b7 * pN7ps + b6 * pN6ps, pN3ps - c7 * pN7ps - c6 * pN6ps,
                                    1.5 * (a8 * pN8ps - a7 * pN7ps), b8 * pN8ps + b7 * pN7ps, pN4ps - c8 * pN8ps - c7 * pN7ps], dtype=float)

                pHypr = np.asarray([1.5 * (d5 * pN5pr - d8 * pN8pr), -pN1pr + e5 * pN5pr + e8 * pN8pr, -b5 * pN5pr - b8 * pN8pr,
                                    1.5 * (d6 * pN6pr - d5 * pN5pr), -pN2pr + e6 * pN6pr + e5 * pN5pr, -b6 * pN6pr - b5 * pN5pr,
                                    1.5 * (d7 * pN7pr - d6 * pN6pr), -pN3pr + e7 * pN7pr + e6 * pN6pr, -b7 * pN7pr - b6 * pN6pr,
                                    1.5 * (d8 * pN8pr - d7 * pN7pr), -pN4pr + e8 * pN8pr + e7 * pN7pr, -b8 * pN8pr - b7 * pN7pr], dtype=float)

                pHyps = np.asarray([1.5 * (d5 * pN5ps - d8 * pN8ps), -pN1ps + e5 * pN5ps + e8 * pN8ps, -pHxps[1],
                                    1.5 * (d6 * pN6ps - d5 * pN5ps), -pN2ps + e6 * pN6ps + e5 * pN5ps, -pHxps[4],
                                    1.5 * (d7 * pN7ps - d6 * pN6ps), -pN3ps + e7 * pN7ps + e6 * pN6ps, -pHxps[7],
                                    1.5 * (d8 * pN8ps - d7 * pN7ps), -pN4ps + e8 * pN8ps + e7 * pN7ps, -pHxps[10]], dtype=float)

                # Jacobi 2*2
                J11 = 0.25 * (x12 - x34 + r * (x12 + x34))
                J12 = 0.25 * (y12 - y34 + r * (y12 + y34))
                J21 = 0.25 * (x23 - x41 + s * (x12 + x34))
                J22 = 0.25 * (y23 - y41 + s * (y12 + y34))
                detJ = 0.125 * (x13 * y24 - x24 * y13 - r * (x12 * y34 + x34 * y12) + s * (x41 * y23 - x23 * y41))

                j11 = J22 / detJ
                j12 = -J12 / detJ
                j21 = -J21 / detJ
                j22 = J11 / detJ

                # B Matrix
                B = np.asarray([j11 * pHxpr + j12 * pHxps,
                                j21 * pHypr + j22 * pHyps,
                                j11 * pHypr + j12 * pHyps + j21 * pHxpr + j22 * pHxps], dtype=float)

                # 这里不会约分掉det_J
                self.K += np.matmul(np.matmul(B.T, self.D), B) * weight[si] * detJ

        # self.K = np.matmul(np.matmul(self.T_matrix.T, self.K), self.T_matrix) * self.thickness
        return self.K * self.thickness

    def ElementStress(self, displacement):
        """
        Calculate element stress
        """


if __name__ == "__main__":
    t_ele = MITC3(-1)
    t_ele.cha_dict = {MaterialKey.Niu: 0.3, MaterialKey.E: 2e9}
    t_ele.node_coords = np.array([[0, 0],
                                [4, 0],
                                [1, 3]], dtype=float)
    t_ele.CalElementDMatrix(MaterialMatrixType.PlaneStree)
    t_ele.ElementStiffness()
    mlogger.debug("finish")
