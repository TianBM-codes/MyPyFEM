import numpy as np
import math
from utils.UtilsFunction import RotateByAxisScipy


class MQ1330:
    sita = np.pi / 180  # 角度弧度转换系数
    L_bj = 24.06 * 1000  # 臂架铰点对铰点长度
    L_bj_xl = 6.63 * 1000  # 臂架底铰点对小拉杆铰点长度
    L_bj_ct = 6.13 * 1000  # 臂架底铰点对齿条架铰点长度
    A_bj = 0  # 臂架倾角，弧度
    A_bjx = 20.758 * sita  # 臂架轴线与小拉杆铰点的夹角
    A_bjc = 25.102 * sita  # 臂架轴线与齿条架铰点的夹角

    H_rzbj = 8.1 * 1000  # 人字架顶部铰点与臂架底铰点的垂直距离，臂架底铰点为0点
    D_rzbj = 6.3 * 1000  # 人字架顶部铰点与臂架底铰点的水平距离
    L_rz = 10.261 * 1000  # 人字架顶部铰点与臂架底铰点的距离
    A_rz = 52.125 * sita  # 人字架顶部铰点与臂架底铰点连线倾角

    L_dl = 22.46 * 1000  # 大拉杆铰点对铰点长度

    L_xl = 7.45 * 1000  # 小拉杆铰点对铰点长度

    L_ph = 2.4 * 1000  # 平衡梁铰点对铰点长度

    L_xb_a = 4.02 * 1000  # 象鼻梁铰点对铰点长度
    A_xb_a = 5.711 * sita  # 象鼻梁5铰点连线与轴线夹角
    A_xb_b = 2.314 * sita  # 象鼻梁1号杆铰点连线与轴线夹角

    H_ct = (6.69 - 0.5) * 1000  # 人字架顶部齿轮与齿条架切点与臂架底铰点的垂直距离，臂架底铰点为0点
    D_ctbj = 3.75 * 1000  # 人字架顶部齿轮与齿条架切点与臂架底铰点的水平距离

    point_bj = np.array([4.400, 0.500, 0]) * 1000  # 臂架底部铰点坐标
    point_rz = np.array([-1.8983, 8.601, 0]) * 1000  # 人字架铰点坐标

    L_xb_1 = 9.908 * 1000  # 象鼻梁1号杆铰点对铰点长度
    L_bj_o = 2.5 * 1000  # 臂架铰点到回转中心距离

    def getPartOrientation(self, Abj=43.224):
        ## 计算各个部件的姿态角
        ## Abj 臂架的姿态角， 角度制
        ## return [臂架角, 象鼻梁角, 大拉杆角, 小拉杆角, 平衡梁角, 齿条架角]  弧度制

        self.A_bj = Abj * self.sita  # 臂架倾角，弧度

        ## 计算大拉杆姿态角
        A_bjrz = np.pi - self.A_bj - self.A_rz  # 臂架与人字架张的角度
        L_rzxb = np.sqrt(self.L_rz ** 2 + self.L_bj ** 2 - 2 * self.L_rz * self.L_bj * math.cos(A_bjrz))  # 人字架铰点到象鼻梁铰点长度
        A_dl_1 = math.acos((self.L_dl ** 2 + L_rzxb ** 2 - self.L_xb_a ** 2) / (2 * self.L_dl * L_rzxb))  # 人字架铰点对象鼻梁两个铰点张的角度
        A_dl_2 = math.asin((self.L_bj * math.sin(self.A_bj) - self.H_rzbj) / L_rzxb)  # 人字架铰点到象鼻梁铰点连线的倾角
        A_dl = A_dl_1 + A_dl_2  # 大拉杆倾角
        # print(f"dalagan A: {A_dl*180/np.pi:.2f}")

        ## 计算象鼻梁姿态角
        A_dlxb = math.acos((self.L_dl ** 2 + self.L_xb_a ** 2 - L_rzxb ** 2) / (2 * self.L_dl * self.L_xb_a))  # 大拉杆与象鼻梁张的角度
        A_xb = np.pi - A_dlxb - A_dl - self.A_xb_a  # 象鼻梁的倾角
        # print(f"xiangbi A: {A_xb*180/np.pi:.2f}")

        ## 计算小拉杆姿态角
        L_rzbj = np.sqrt(self.L_rz ** 2 + self.L_bj_xl ** 2 - 2 * self.L_rz * self.L_bj_xl * math.cos(A_bjrz - self.A_bjx))  # 人字架铰点对臂架小拉杆铰点距离
        A_rzxl_1 = math.acos((L_rzbj ** 2 + self.L_xl ** 2 - self.L_ph ** 2) / (2 * L_rzbj * self.L_xl))  # 人字架铰点与臂架小拉杆铰点连线与小拉杆夹角
        A_rzxl_2 = math.asin((self.H_rzbj - self.L_bj_xl * math.sin(self.A_bj + self.A_bjx)) / L_rzbj)  # 人字架铰点与臂架小拉杆铰点连线倾角
        A_xl = A_rzxl_1 + A_rzxl_2  # 小拉杆倾角
        # print(f"xiaolagan A: {A_xl*180/np.pi:.2f}")

        ## 计算平衡梁姿态角
        # A_ph = math.asin((L_xl*math.sin(A_xl) + L_bj_xl*math.sin(A_bj + A_bjx) - H_rzbj)/L_ph)  #平衡梁铰点连线倾角
        A_ph = math.acos((self.D_rzbj + self.L_bj_xl * math.cos(self.A_bj + self.A_bjx) - self.L_xl * math.cos(A_xl)) / self.L_ph)  # 平衡梁铰点连线倾角
        # print(f"phl A: {A_ph*180/np.pi:.2f}")

        ## 计算齿条架姿态角
        A_ct = math.atan((self.H_ct - self.L_bj_ct * math.sin(self.A_bj + self.A_bjc)) / (self.D_ctbj + self.L_bj_ct * math.cos(self.A_bj + self.A_bjc)))  # 齿条架倾角
        # print(f"chitiao A: {A_ct*180/np.pi:.2f}")

        ## 底部圆
        self.A_dl = A_dl_1 + A_dl_2  # 大拉杆倾角
        self.A_xb = -(np.pi - A_dlxb - self.A_dl - self.A_xb_a)  # 象鼻梁的倾角
        self.D = self.L_bj * np.cos(self.A_bj) + self.L_xb_1 * np.cos(self.A_xb + self.A_xb_b) + self.L_bj_o
        self.H = self.L_bj * np.sin(self.A_bj) + self.L_xb_1 * np.sin(self.A_xb + self.A_xb_b) + self.point_bj[1]
        return [self.A_bj, -A_xb, A_dl, -A_xl, A_ph, -A_ct]

    def getPartCircle(self):
        return self.D, self.H

    def getPartTranslate(self, Abj1, Abj2):
        ## 计算各个部件的转移变换量（平移矢量，旋转角，旋转中心）
        ## Abj1 当前状态的臂架姿态角， 角度制
        ## Abj2 下一个状态的臂架姿态角， 角度制
        ## return {"PartName": [vector, Angle, point], ...} ，角度值是弧度制，长度单位米

        AngleData1 = self.getPartOrientation(Abj1)
        AngleData2 = self.getPartOrientation(Abj2)
        # print(f"AngleData1: {np.array(AngleData1)*180/np.pi}")
        # print(f"AngleData2: {np.array(AngleData2)*180/np.pi}")
        res = {}
        ## 臂架
        res["bijia"] = [np.array([0, 0, 0]), AngleData2[0] - AngleData1[0], self.point_bj]

        ## 象鼻梁
        p1 = np.array([self.L_bj * math.cos(AngleData1[0]), self.L_bj * math.sin(AngleData1[0]), 0])
        p2 = np.array([self.L_bj * math.cos(AngleData2[0]), self.L_bj * math.sin(AngleData2[0]), 0])
        vec = p2 - p1
        res["xiangbiliang"] = [vec, AngleData2[1] - AngleData1[1], self.point_bj + p2]

        ## 小拉杆
        p1 = np.array([self.L_bj_xl * math.cos(AngleData1[0] + self.A_bjx), self.L_bj_xl * math.sin(AngleData1[0] + self.A_bjx), 0])
        p2 = np.array([self.L_bj_xl * math.cos(AngleData2[0] + self.A_bjx), self.L_bj_xl * math.sin(AngleData2[0] + self.A_bjx), 0])
        vec = p2 - p1
        res["xiaolagan"] = [vec, AngleData2[3] - AngleData1[3], self.point_bj + p2]

        ## 齿条架
        p1 = np.array([self.L_bj_ct * math.cos(AngleData1[0] + self.A_bjc), self.L_bj_ct * math.sin(AngleData1[0] + self.A_bjc), 0])
        p2 = np.array([self.L_bj_ct * math.cos(AngleData2[0] + self.A_bjc), self.L_bj_ct * math.sin(AngleData2[0] + self.A_bjc), 0])
        vec = p2 - p1
        res["chitiaojia"] = [vec, AngleData2[5] - AngleData1[5], self.point_bj + p2]

        ## 大拉杆
        res["dalagan"] = [np.array([0, 0, 0]), AngleData2[2] - AngleData1[2], self.point_rz]

        ## 平衡梁
        res["pinghengliang"] = [np.array([0, 0, 0]), AngleData2[4] - AngleData1[4], self.point_rz]

        return res


class MQ1330Wrapper:
    def __init__(self, rotate_angle=30):
        """
        初始化MQ1330封装类, angle_range: int/float, 角度范围参数，默认为30
        :param rotate_angle:
        """
        self.mq1330 = MQ1330()
        self.angle_range = rotate_angle
        self.z_axis = np.array((0, 0, 1), dtype=float)
        self._init_transforms()

    def _init_transforms(self):
        """初始化所有变换参数"""
        Abj = 43.224
        data = self.mq1330.getPartTranslate(Abj, Abj + self.angle_range)

        # 初始化各个部件的变换参数
        self.bijia_trans, self.bijia_theta, bijia_origin = data["bijia"]
        self.bijia_o = bijia_origin

        self.dalagan_trans, self.dalagan_theta, dalagan_origin = data["dalagan"]
        self.dalagan_o = dalagan_origin

        self.phl_trans, self.phl_theta, phl_origin = data["pinghengliang"]
        self.phl_o = phl_origin

        self.xlg_trans, self.xlg_theta, xlg_origin = data["xiaolagan"]
        self.xlg_o = xlg_origin

        self.xbl_trans, self.xbl_theta, xbl_origin = data["xiangbiliang"]
        self.xbl_o = xbl_origin

        self.ctj_trans, self.ctj_theta, ctj_origin = data["chitiaojia"]
        self.ctj_o = ctj_origin

        # 1. 臂架 (bijia)
        R_bijia = np.array([
            [np.cos(self.bijia_theta), -np.sin(self.bijia_theta), 0],
            [np.sin(self.bijia_theta), np.cos(self.bijia_theta), 0],
            [0, 0, 1]
        ])
        self.R_bijia_4 = np.zeros((24, 24))
        self.R_bijia_4[0:3, 0:3] = R_bijia
        self.R_bijia_4[3:6, 3:6] = R_bijia
        self.R_bijia_4[6:9, 6:9] = R_bijia
        self.R_bijia_4[9:12, 9:12] = R_bijia
        self.R_bijia_4[12:15, 12:15] = R_bijia
        self.R_bijia_4[15:18, 15:18] = R_bijia
        self.R_bijia_4[18:21, 18:21] = R_bijia
        self.R_bijia_4[21:24, 21:24] = R_bijia

        self.R_bijia_3 = np.zeros((18, 18))
        self.R_bijia_3[0:3, 0:3] = R_bijia
        self.R_bijia_3[3:6, 3:6] = R_bijia
        self.R_bijia_3[6:9, 6:9] = R_bijia
        self.R_bijia_3[9:12, 9:12] = R_bijia
        self.R_bijia_3[12:15, 12:15] = R_bijia
        self.R_bijia_3[15:18, 15:18] = R_bijia

        # 2. 拉杆 (dalagan)
        R_dalagan = np.array([
            [np.cos(self.dalagan_theta), -np.sin(self.dalagan_theta), 0],
            [np.sin(self.dalagan_theta), np.cos(self.dalagan_theta), 0],
            [0, 0, 1]
        ])

        self.R_dalagan_4 = np.zeros((24, 24))
        self.R_dalagan_4[0:3, 0:3] = R_dalagan
        self.R_dalagan_4[3:6, 3:6] = R_dalagan
        self.R_dalagan_4[6:9, 6:9] = R_dalagan
        self.R_dalagan_4[9:12, 9:12] = R_dalagan
        self.R_dalagan_4[12:15, 12:15] = R_dalagan
        self.R_dalagan_4[15:18, 15:18] = R_dalagan
        self.R_dalagan_4[18:21, 18:21] = R_dalagan
        self.R_dalagan_4[21:24, 21:24] = R_dalagan

        self.R_dalagan_3 = np.zeros((18, 18))
        self.R_dalagan_3[0:3, 0:3] = R_dalagan
        self.R_dalagan_3[3:6, 3:6] = R_dalagan
        self.R_dalagan_3[6:9, 6:9] = R_dalagan
        self.R_dalagan_3[9:12, 9:12] = R_dalagan
        self.R_dalagan_3[12:15, 12:15] = R_dalagan
        self.R_dalagan_3[15:18, 15:18] = R_dalagan

        # 3. 平衡梁 (phl)
        R_phl = np.array([
            [np.cos(self.phl_theta), -np.sin(self.phl_theta), 0],
            [np.sin(self.phl_theta), np.cos(self.phl_theta), 0],
            [0, 0, 1]
        ])

        self.R_phl_4 = np.zeros((24, 24))
        self.R_phl_4[0:3, 0:3] = R_phl
        self.R_phl_4[3:6, 3:6] = R_phl
        self.R_phl_4[6:9, 6:9] = R_phl
        self.R_phl_4[9:12, 9:12] = R_phl
        self.R_phl_4[12:15, 12:15] = R_phl
        self.R_phl_4[15:18, 15:18] = R_phl
        self.R_phl_4[18:21, 18:21] = R_phl
        self.R_phl_4[21:24, 21:24] = R_phl

        self.R_phl_3 = np.zeros((18, 18))
        self.R_phl_3[0:3, 0:3] = R_phl
        self.R_phl_3[3:6, 3:6] = R_phl
        self.R_phl_3[6:9, 6:9] = R_phl
        self.R_phl_3[9:12, 9:12] = R_phl
        self.R_phl_3[12:15, 12:15] = R_phl
        self.R_phl_3[15:18, 15:18] = R_phl

        # 4. 小拉杆 (xlg)
        R_xlg = np.array([
            [np.cos(self.xlg_theta), -np.sin(self.xlg_theta), 0],
            [np.sin(self.xlg_theta), np.cos(self.xlg_theta), 0],
            [0, 0, 1]
        ])

        self.R_xlg_4 = np.zeros((24, 24))
        self.R_xlg_4[0:3, 0:3] = R_xlg
        self.R_xlg_4[3:6, 3:6] = R_xlg
        self.R_xlg_4[6:9, 6:9] = R_xlg
        self.R_xlg_4[9:12, 9:12] = R_xlg
        self.R_xlg_4[12:15, 12:15] = R_xlg
        self.R_xlg_4[15:18, 15:18] = R_xlg
        self.R_xlg_4[18:21, 18:21] = R_xlg
        self.R_xlg_4[21:24, 21:24] = R_xlg

        self.R_xlg_3 = np.zeros((18, 18))
        self.R_xlg_3[0:3, 0:3] = R_xlg
        self.R_xlg_3[3:6, 3:6] = R_xlg
        self.R_xlg_3[6:9, 6:9] = R_xlg
        self.R_xlg_3[9:12, 9:12] = R_xlg
        self.R_xlg_3[12:15, 12:15] = R_xlg
        self.R_xlg_3[15:18, 15:18] = R_xlg

        # 5. 像鼻梁 (xbl)
        R_xbl = np.array([
            [np.cos(self.xbl_theta), -np.sin(self.xbl_theta), 0],
            [np.sin(self.xbl_theta), np.cos(self.xbl_theta), 0],
            [0, 0, 1]
        ])

        self.R_xbl_4 = np.zeros((24, 24))
        self.R_xbl_4[0:3, 0:3] = R_xbl
        self.R_xbl_4[3:6, 3:6] = R_xbl
        self.R_xbl_4[6:9, 6:9] = R_xbl
        self.R_xbl_4[9:12, 9:12] = R_xbl
        self.R_xbl_4[12:15, 12:15] = R_xbl
        self.R_xbl_4[15:18, 15:18] = R_xbl
        self.R_xbl_4[18:21, 18:21] = R_xbl
        self.R_xbl_4[21:24, 21:24] = R_xbl

        self.R_xbl_3 = np.zeros((18, 18))
        self.R_xbl_3[0:3, 0:3] = R_xbl
        self.R_xbl_3[3:6, 3:6] = R_xbl
        self.R_xbl_3[6:9, 6:9] = R_xbl
        self.R_xbl_3[9:12, 9:12] = R_xbl
        self.R_xbl_3[12:15, 12:15] = R_xbl
        self.R_xbl_3[15:18, 15:18] = R_xbl

        # 6. 齿条架 (ctj)
        R_ctj = np.array([
            [np.cos(self.ctj_theta), -np.sin(self.ctj_theta), 0],
            [np.sin(self.ctj_theta), np.cos(self.ctj_theta), 0],
            [0, 0, 1]
        ])
        self.R_ctj_4 = np.zeros((24, 24))
        self.R_ctj_4[0:3, 0:3] = R_ctj
        self.R_ctj_4[3:6, 3:6] = R_ctj
        self.R_ctj_4[6:9, 6:9] = R_ctj
        self.R_ctj_4[9:12, 9:12] = R_ctj
        self.R_ctj_4[12:15, 12:15] = R_ctj
        self.R_ctj_4[15:18, 15:18] = R_ctj
        self.R_ctj_4[18:21, 18:21] = R_ctj
        self.R_ctj_4[21:24, 21:24] = R_ctj

        self.R_ctj_3 = np.zeros((18, 18))
        self.R_ctj_3[0:3, 0:3] = R_ctj
        self.R_ctj_3[3:6, 3:6] = R_ctj
        self.R_ctj_3[6:9, 6:9] = R_ctj
        self.R_ctj_3[9:12, 9:12] = R_ctj
        self.R_ctj_3[12:15, 12:15] = R_ctj
        self.R_ctj_3[15:18, 15:18] = R_ctj

    def is_fixed_element(self, e_id):
        """
        判断单元是否发生转动
        :param e_id:
        :return:
        """
        if e_id < 20000:
            return True
        else:
            return False

    def calculate_new_xyz(self, node_id, x, y, z):
        """
        计算节点的新坐标
        :param node_id:
        :param x:
        :param y:
        :param z:
        :return:
        """
        old_xyz = np.array((x, y, z))

        if node_id < 20000:
            """
            人字架不动
            """
            return x, y, z
        elif 20000 < node_id < 50000:
            """
            臂架只绕指定点旋转
            """
            return RotateByAxisScipy(self.z_axis, self.bijia_o, self.bijia_theta * 180 / np.pi, old_xyz)
        elif 50000 < node_id < 70000:
            """
            大拉杆只绕指定点旋转
            """
            return RotateByAxisScipy(self.z_axis, self.dalagan_o, self.dalagan_theta * 180 / np.pi, old_xyz)
        elif 70000 < node_id < 90000:
            """
            平衡梁只绕指定点旋转
            """
            return RotateByAxisScipy(self.z_axis, self.phl_o, self.phl_theta * 180 / np.pi, old_xyz)
        elif 90000 < node_id < 100000:
            """
            小拉杆
            """
            new_x, new_y, new_z = RotateByAxisScipy(self.z_axis, self.bijia_o, self.bijia_theta * 180 / np.pi, old_xyz)
            return RotateByAxisScipy(self.z_axis, self.xlg_o, (-self.bijia_theta + self.xlg_theta) * 180 / np.pi,
                                     np.array((new_x, new_y, new_z)))
        elif 100000 < node_id < 130000:
            """
            象鼻梁先随臂架旋转后，再绕臂架尾部旋转
            """
            new_x, new_y, new_z = RotateByAxisScipy(self.z_axis, self.bijia_o, self.bijia_theta * 180 / np.pi, old_xyz)
            return RotateByAxisScipy(self.z_axis, self.xbl_o, (-self.bijia_theta + self.xbl_theta) * 180 / np.pi,
                                     np.array((new_x, new_y, new_z)))
        elif 130000 < node_id:
            """
            齿条架
            """
            new_x, new_y, new_z = RotateByAxisScipy(self.z_axis, self.bijia_o, self.bijia_theta * 180 / np.pi, old_xyz)
            return RotateByAxisScipy(self.z_axis, self.ctj_o, (-self.bijia_theta + self.ctj_theta) * 180 / np.pi,
                                     np.array((new_x, new_y, new_z)))
        else:
            raise ValueError(f"Invalid node_id: {node_id}")

    def get_rotate_info(self, set_name):
        """
        计算不同set的旋转信息
        :param set_name:
        :return:
        """
        if set_name in ["bijia_ele"]:
            return [self.z_axis[0], self.z_axis[1], self.z_axis[2],
                    self.bijia_o[0], self.bijia_o[1], self.bijia_o[2],
                    self.bijia_theta * 180 / np.pi]
        elif set_name in ["chitiao_ele"]:
            return [self.z_axis[0], self.z_axis[1], self.z_axis[2],
                    self.bijia_o[0], self.bijia_o[1], self.bijia_o[2],
                    self.bijia_theta * 180 / np.pi,
                    self.z_axis[0], self.z_axis[1], self.z_axis[2],
                    self.ctj_o[0], self.ctj_o[1], self.ctj_o[2],
                    (-self.bijia_theta + self.ctj_theta) * 180 / np.pi]
        elif set_name in ["dalagan_ele"]:
            return [self.z_axis[0], self.z_axis[1], self.z_axis[2],
                    self.dalagan_o[0], self.dalagan_o[1], self.dalagan_o[2],
                    self.dalagan_theta * 180 / np.pi]
        elif set_name in ["phl_left_ele", "phl_right_ele"]:
            return [self.z_axis[0], self.z_axis[1], self.z_axis[2],
                    self.phl_o[0], self.phl_o[1], self.phl_o[2],
                    self.phl_theta * 180 / np.pi]
        elif set_name in ["renzijia_ele"]:
            return []
        elif set_name in ["xiangbiliang_ele"]:
            return [self.z_axis[0], self.z_axis[1], self.z_axis[2],
                    self.bijia_o[0], self.bijia_o[1], self.bijia_o[2],
                    self.bijia_theta * 180 / np.pi,
                    self.z_axis[0], self.z_axis[1], self.z_axis[2],
                    self.xbl_o[0], self.xbl_o[1], self.xbl_o[2],
                    (-self.bijia_theta + self.xbl_theta) * 180 / np.pi]
        elif set_name in ["xiaolagan_ele"]:
            return [self.z_axis[0], self.z_axis[1], self.z_axis[2],
                    self.bijia_o[0], self.bijia_o[1], self.bijia_o[2],
                    self.bijia_theta * 180 / np.pi,
                    self.z_axis[0], self.z_axis[1], self.z_axis[2],
                    self.xlg_o[0], self.xlg_o[1], self.xlg_o[2],
                    (-self.bijia_theta + self.xlg_theta) * 180 / np.pi]
        else:
            raise KeyError(set_name)

    def calculate_new_stiff(self, e_id, node_count, stiff):
        """
        计算新的刚度, 也就是旋转后的
        :param e_id:
        :param node_count:
        :param stiff:
        :return:
        """
        if 20000 < e_id < 50000:
            if node_count == 3:
                return self.R_bijia_3.T @ stiff @ self.R_bijia_3
            elif node_count == 4:
                return self.R_bijia_4.T @ stiff @ self.R_bijia_4
            else:
                raise ValueError(node_count)
        elif 50000 < e_id < 70000:
            if node_count == 3:
                return self.R_dalagan_3.T @ stiff @ self.R_dalagan_3
            elif node_count == 4:
                return self.R_dalagan_4.T @ stiff @ self.R_dalagan_4
            else:
                raise ValueError(node_count)
        elif 70000 < e_id < 90000:
            if node_count == 3:
                return self.R_phl_3.T @ stiff @ self.R_phl_3
            elif node_count == 4:
                return self.R_phl_4.T @ stiff @ self.R_phl_4
            else:
                raise ValueError(node_count)
        elif 90000 < e_id < 100000:
            if node_count == 3:
                return self.R_xlg_3.T @ stiff @ self.R_xlg_3
            elif node_count == 4:
                return self.R_xlg_4.T @ stiff @ self.R_xlg_4
            else:
                raise ValueError(node_count)
        elif 100000 < e_id < 130000:
            if node_count == 3:
                return self.R_xbl_3.T @ stiff @ self.R_xbl_3
            elif node_count == 4:
                return self.R_xbl_4.T @ stiff @ self.R_xbl_4
            else:
                raise ValueError(node_count)
        elif 130000 < e_id < 140000:
            if node_count == 3:
                return self.R_ctj_3.T @ stiff @ self.R_ctj_3
            elif node_count == 4:
                return self.R_ctj_4.T @ stiff @ self.R_ctj_4
            else:
                raise ValueError(node_count)
        else:
            raise ValueError()

    def get_circle_info(self):
        return self.mq1330.getPartCircle()
