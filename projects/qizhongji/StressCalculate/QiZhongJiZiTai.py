import numpy as np
import math


class MQ1330:
    sita = np.pi/180  #角度弧度转换系数
    L_bj = 24.06 * 1000  #臂架铰点对铰点长度
    L_bj_xl = 6.63 * 1000   #臂架底铰点对小拉杆铰点长度
    L_bj_ct = 6.13 * 1000   #臂架底铰点对齿条架铰点长度
    L_bj_o = 2.5 *1000   # 臂架铰点到回转中心距离
    A_bj = 0  #臂架倾角，弧度
    A_bjx = 20.758 * sita  #臂架轴线与小拉杆铰点的夹角
    A_bjc = 25.102 * sita  #臂架轴线与齿条架铰点的夹角

    H_rzbj = 8.1 * 1000   #人字架顶部铰点与臂架底铰点的垂直距离，臂架底铰点为0点
    D_rzbj = 6.3 * 1000   #人字架顶部铰点与臂架底铰点的水平距离
    L_rz = 10.261 * 1000   #人字架顶部铰点与臂架底铰点的距离
    A_rz = 52.125 * sita  #人字架顶部铰点与臂架底铰点连线倾角

    L_dl = 22.46 * 1000   #大拉杆铰点对铰点长度
    A_dl = 0.0 ## 大拉杆倾角

    L_xl = 7.45 * 1000   #小拉杆铰点对铰点长度
    A_xl = 0.0 ## 小拉杆倾角

    L_ph = 2.4 * 1000   #平衡梁铰点对铰点长度
    A_ph = 0.0 ## 平衡梁倾角

    L_xb_a = 4.02 * 1000   #象鼻梁铰点对铰点长度
    A_xb_a = 5.711 * sita  #象鼻梁铰点连线与轴线夹角
    A_xb = 0.0 ## 象鼻梁倾角
    L_xb_1 = 9.908 *1000  #象鼻梁1号杆铰点对铰点长度

    H_ct = (6.69 - 0.5) * 1000   #人字架顶部齿轮与齿条架切点与臂架底铰点的垂直距离，臂架底铰点为0点
    D_ctbj = 3.75 * 1000   #人字架顶部齿轮与齿条架切点与臂架底铰点的水平距离
    A_ct = 0.0 ## 齿条架倾角

    D = 0  # 起重机工作幅度
    H = 0  # 起重机象鼻梁吊钩滑轮铰点高度

    point_bj = np.array([4.400, 0.500, 0])  * 1000   #臂架底部铰点坐标
    point_rz = np.array([-1.8983, 8.601, 0]) * 1000   #人字架铰点坐标

    def getPartOrientation(self, Abj = 43.224):
        ## 计算各个部件的姿态角
        ## Abj 臂架的姿态角， 角度制
        ## return [臂架角, 象鼻梁角, 大拉杆角, 小拉杆角, 平衡梁角, 齿条架角]  弧度制

        self.A_bj = Abj * self.sita  #臂架倾角，弧度
        
        ## 计算大拉杆姿态角
        A_bjrz = np.pi - self.A_bj - self.A_rz  #臂架与人字架张的角度
        L_rzxb = np.sqrt(self.L_rz**2 + self.L_bj**2 - 2*self.L_rz*self.L_bj*math.cos(A_bjrz))  #人字架铰点到象鼻梁铰点长度
        A_dl_1 = math.acos((self.L_dl**2 + L_rzxb**2 - self.L_xb_a**2 )/(2*self.L_dl*L_rzxb))  #人字架铰点对象鼻梁两个铰点张的角度
        A_dl_2 = math.asin((self.L_bj*math.sin(self.A_bj) - self.H_rzbj)/L_rzxb)  #人字架铰点到象鼻梁铰点连线的倾角
        self.A_dl = A_dl_1 + A_dl_2  #大拉杆倾角
        # print(f"dalagan A: {A_dl*180/np.pi:.2f}")
        
        ## 计算象鼻梁姿态角
        A_dlxb = math.acos((self.L_dl**2 + self.L_xb_a**2 - L_rzxb**2 )/(2*self.L_dl*self.L_xb_a))  #大拉杆与象鼻梁张的角度
        self.A_xb = -(np.pi - A_dlxb - self.A_dl - self.A_xb_a)  #象鼻梁的倾角
        self.D = self.L_bj * np.cos(self.A_bj) + self.L_xb_1 * np.cos(self.A_xb+self.A_xb_a) + self.L_bj_o
        self.H = self.L_bj * np.sin(self.A_bj) + self.L_xb_1 * np.sin(self.A_xb+self.A_xb_a) + self.point_bj[1]
        # print(f"xiangbi A: {A_xb*180/np.pi:.2f}")
        
        ## 计算小拉杆姿态角
        L_rzbj = np.sqrt(self.L_rz**2 + self.L_bj_xl**2 - 2*self.L_rz*self.L_bj_xl*math.cos(A_bjrz-self.A_bjx))  #人字架铰点对臂架小拉杆铰点距离
        A_rzxl_1 = math.acos((L_rzbj**2 + self.L_xl**2 - self.L_ph**2 )/(2*L_rzbj*self.L_xl))  #人字架铰点与臂架小拉杆铰点连线与小拉杆夹角
        A_rzxl_2 = math.asin((self.H_rzbj - self.L_bj_xl*math.sin(self.A_bj+self.A_bjx))/L_rzbj)  #人字架铰点与臂架小拉杆铰点连线倾角
        self.A_xl = -(A_rzxl_1 + A_rzxl_2)  #小拉杆倾角
        # print(f"xiaolagan A: {A_xl*180/np.pi:.2f}")

        ## 计算平衡梁姿态角
        # A_ph = math.asin((L_xl*math.sin(A_xl) + L_bj_xl*math.sin(A_bj + A_bjx) - H_rzbj)/L_ph)  #平衡梁铰点连线倾角
        self.A_ph = math.acos((self.D_rzbj + self.L_bj_xl*math.cos(self.A_bj+self.A_bjx) - self.L_xl*math.cos(self.A_xl))/self.L_ph)  #平衡梁铰点连线倾角
        # print(f"phl A: {A_ph*180/np.pi:.2f}")

        ## 计算齿条架姿态角
        self.A_ct = -math.atan((self.H_ct - self.L_bj_ct*math.sin(self.A_bj + self.A_bjc))/(self.D_ctbj + self.L_bj_ct*math.cos(self.A_bj + self.A_bjc)))  #齿条架倾角
        # print(f"chitiao A: {A_ct*180/np.pi:.2f}")
        return [self.A_bj, self.A_xb, self.A_dl, self.A_xl, self.A_ph, self.A_ct]
    
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
        res["bijia"] = [np.array([0,0,0]), AngleData2[0] - AngleData1[0], self.point_bj]

        ## 象鼻梁
        p1 = np.array([self.L_bj * math.cos(AngleData1[0]),  self.L_bj * math.sin(AngleData1[0]), 0])
        p2 = np.array([self.L_bj * math.cos(AngleData2[0]),  self.L_bj * math.sin(AngleData2[0]), 0])
        vec = p2 - p1
        res["xiangbiliang"] = [vec, AngleData2[1] - AngleData1[1], self.point_bj+p2]

        ## 小拉杆
        p1 = np.array([self.L_bj_xl * math.cos(AngleData1[0]+self.A_bjx),  self.L_bj_xl * math.sin(AngleData1[0]+self.A_bjx), 0])
        p2 = np.array([self.L_bj_xl * math.cos(AngleData2[0]+self.A_bjx),  self.L_bj_xl * math.sin(AngleData2[0]+self.A_bjx), 0])
        vec = p2 - p1
        res["xiaolagan"] = [vec, AngleData2[3] - AngleData1[3], self.point_bj+p2]

        ## 齿条架
        p1 = np.array([self.L_bj_ct * math.cos(AngleData1[0]+self.A_bjc),  self.L_bj_ct * math.sin(AngleData1[0]+self.A_bjc), 0])
        p2 = np.array([self.L_bj_ct * math.cos(AngleData2[0]+self.A_bjc),  self.L_bj_ct * math.sin(AngleData2[0]+self.A_bjc), 0])
        vec = p2 - p1
        res["chitiaojia"] = [vec, AngleData2[5] - AngleData1[5], self.point_bj+p2]

        ## 大拉杆
        res["dalagan"] = [np.array([0,0,0]), AngleData2[2] - AngleData1[2], self.point_rz]
        
        ## 平衡梁
        res["pinghengliang"] = [np.array([0,0,0]), AngleData2[4] - AngleData1[4], self.point_rz]

        return res
        


if __name__ == "__main__":
    mq1330 = MQ1330()
    # Abj = 43.224
    # # d0 = mq1330.getPartOrientation(43.224)
    # # print(f"d0: {np.round(np.array(d0)*180/np.pi,2)}")
    # for i in range(30):
    #     data = mq1330.getPartTranslate(Abj, Abj + i)
    #     print(f"R: {mq1330.D/1000:.1f}, H: {mq1330.H/1000:.1f}")
    # for k, j in data.items():
    #     print(f"Name: {k} -> [{j[0]}, {j[1] * 180 / np.pi:.2f}, {j[2]}]")
    # 假设目标臂架角度
    Abj2 = 70 # 或 43.224 + 30

    # 调用获取绝对姿态角
    final_angles = mq1330.getPartOrientation(Abj2)

    # 打印（转换成度数更直观）
    angle_names = ["bijia", "xiangbiliang", "dalagan", "xiaolagan", "pinghengliang", "chitiaojia"]

    for name, angle_rad in zip(angle_names, final_angles):
        print(f"{name} 绝对角度: {angle_rad * 180 / np.pi:.2f}°")
