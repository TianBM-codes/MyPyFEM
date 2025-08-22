import numpy as np 
import textwrap
import math

class ParseForceFormBeam():
    L1 = 7.162 #8.766
    q12 = 9.75/180*np.pi  ## 1,2杆夹角
    A = [0.020368, 0.00608]
    def getBeamInternalForce(self, e, exy, I, A, E, Nxy):
        ## 根据应变计算梁的内力，包括弯矩和轴力
        ## e：截面测点应变值list [e1， e2， e3， e4, ...]
        ## exy: 截面测点位置坐标list [[e1x, e1y], [e2x,e2y], [e3x, e3y], [e4x, e4y], ...]
        ## I: 截面惯性矩list [Ix， Iy]
        ## H：截面尺寸list [hx， hy]
        ## A: 截面面积
        ## E: 材料弹性模量
        ## Nxy: 截面外轮廓角点坐标list [[n1x, n1y], [n2x, n2y], ...]
        ## return [Mx, My, Fz]
        ## 截面坐标系方向，正对截面，水平向右X正，垂直向上Y正，正对截面向外Z正

        try:
            # 定义 4 个坐标点
            points =[]
            for i in range(len(exy)):
                points.append([exy[i][0], exy[i][1], e[i]])
            # print("point: ", points)
            points = np.array(points)

            # 计算质心
            centroid = np.mean(points, axis=0)
            # print("centroid: ", centroid)
            # centroid = (points[0]+points[1])/2
            # print("centroid: ", centroid)
            # 将点平移到质心
            shifted_points = points - centroid
            # 计算协方差矩阵
            cov_matrix = np.cov(shifted_points.T)
            # 计算协方差矩阵的特征值和特征向量
            eigenvalues, eigenvectors = np.linalg.eig(cov_matrix)
            # 最小特征值对应的特征向量是平面的法向量
            normal_vector = eigenvectors[:, np.argmin(eigenvalues)]

            # 平面方程为 ax + by + cz + d = 0，其中法向量为 [a, b, c]
            a, b, c = normal_vector
            # 计算 d
            d = -np.dot(normal_vector, centroid)
            # 打印拟合的平面方程
            # print(f"拟合的平面方程为: {a/c:.8f}x + {b/c:.8f}y + {c/c:.8f}z + {d/c:.8f} = 0")
            # for i in range(len(points)):
            #     print(np.round((np.dot(normal_vector, points[i]) + d)*10**6,1))


            # z = -(ax + by + d)/c
            ## 内力计算公式
            ## Mx = E * Ix * -tan(xita_x)  xita_x 为法向量与z轴夹角在yz平面的投影角  tan(xita_x) = y/z;  y,z 为法向量的坐标
            ## My = E * Iy * tan(xita_y)  xita_y 为法向量与z轴夹角在xz平面的投影角  tan(xita_y) = x/z;  x,z 为法向量的坐标
            ## Fz = E * A * centroid[2]   为质心处应变为轴向拉伸应变

            Mx = E * I[0] * -normal_vector[1]/normal_vector[2]
            My = E * I[1] * normal_vector[0]/normal_vector[2]
            ## 计算外轮廓角点的应变值
            ez = []
            for nxy in Nxy:
                z = -(a * nxy[0] + b * nxy[1] + d)/c 
                ez.append(z)
            # print("ez: ", ez, "\n ezz: ", np.mean(ez))
            zz = np.mean(ez)
            Fz = E * A * zz

            return [Mx, My, Fz]
        except Exception as e:
            print(f"in getBeamInternalForce err: {e}")

    def getBarInternalForce(self, e, A, E):
        ## 根据应变计算杆的轴力
        ## e：截面测点应变值
        ## A: 截面面积
        ## E: 材料弹性模量
        ## return [Fx]
        ## Fx = e * A * E   方向：向截面内为-， 向截面外为+
        ## 截面坐标系方向，正对截面，水平向右X正，垂直向上Y正，正对截面向外Z正
        return e * A * E

    def getLiftingLoad(self, e ): ##, Q, L, A):
        ## e: 象鼻梁各个测点应变数据list(去除重力应变)
        ## Q：象鼻梁杆与杆之间的夹角
        ## L: 弯曲截面到载荷作用点的距离
        ## A2: 弯曲截面面积
        ## return liftingLoad
        ## 截面坐标系方向，正对截面，水平向右X正，垂直向上Y正
        # print("e: ", e)
        # exy1 = [[-0.173, 0.175], [0.173, 0.175], [0.405, 0.00], [0, -0.175], [-0.405, 0.00]]
        # exy1 = [[-0.173, 0.175], [0.173, 0.175], [0.405, 0.00], [-0.405, 0.00]]
        try:
            exy1 = [[-0.405, 0.00], [-0.173, 0.175], [0.173, 0.175], [0.405, 0.00]]     ## 顺序[左，左上，右上，右]
            II1 = [0.512155*10**-3, 1.733876*10**-3]
            # II1 = [4.76928*10**-4, 1.604*10**-3]
            # A1 = 0.01856
            A1 = self.A[0]
            E = 2.1*10**11
            Nxy1 = [[-0.404,-0.175], [0.404, -0.175], [0.405, 0.175], [-0.404, 0.175]]
            # e[:6] = [4.4,     5.7,  -107.1,  -221.5,  -114.7,   321.8]
            # ess = [4.4,     5.7,  -107.1,  -114.7]
            MF1 = self.getBeamInternalForce(np.array(e[:4])*10**-6, exy1, II1, A1, E, Nxy1)
            # MF1 = self.getBeamInternalForce(np.array(ess)*10**-6, exy1, II1, A1, E, Nxy1)
            # print("MF: ", MF1)

            ##对吊点进行力的平衡
            Fy1 = MF1[0]/self.L1     ## 外力，全局坐标系，y轴向上
            Fx1 = -MF1[2]      ## 外力，全局坐标系，x轴沿1杆向外
            # print(f"Fx1:{Fx1}, Fy1:{Fy1}")

            A2 = self.A[1]
            F2 = self.getBarInternalForce(e[4]*10**-6, A2, E)
            # print(f"F2:{F2}")
            q12 = self.q12  #Q ## 14.72/180*np.pi
            Fx2 = -F2*np.cos(q12)    ## 外力，全局坐标系，x轴沿1杆向外
            Fy2 = F2 * np.sin(q12)   ## 外力，全局坐标系，y轴向上
            # print(f"Fx2:{Fx2}, Fy2:{Fy2}")


            ## 吊点的外力
            Fx = -(Fx1 + Fx2)
            Fy = -(Fy1 + Fy2)
            # print(f"Fx:{Fx:.1f}, Fy:{Fy:.1f}")
            # print(f"sita F: {math.atan(Fy/Fx)*180/np.pi}")
            return [Fx, Fy, MF1[1]]  ## X向集中力，Y向集中力，侧摆力矩
        except Exception as e:
            print(f"in getLiftingLoad err: {e}")


if __name__ == "__main__":
    from pathlib import Path
    current_file = Path(__file__).resolve()
    print(current_file)          # /home/user/project/main.py

    # 所在目录
    current_dir = current_file.parent
    # 打开文件并读取所有行
    with open(current_dir/'gravity_res.txt', 'r', encoding='utf-8') as file:
        lines = file.readlines()

    # 打印每一行
    data_G =[]
    for line in lines:
        data_G.append([float(line[i*8:i*8+8]) for i in range(13)])  # 使用 strip() 去掉每行末尾的换行符

    with open(current_dir/'gravity_F_res.txt', 'r', encoding='utf-8') as file:
        lines = file.readlines()
    # 打印每一行
    data_GF =[]
    for line in lines:
        # print(len(line))
        data_GF.append([float(line[i*8:i*8+8]) for i in range(13)])  # 使用 strip() 去掉每行末尾的换行符


    # print(data_G)
    # print(data_GF)
    data=[]
    for i in range(10):
        data.append(np.array(data_GF[i]) - np.array(data_G[i]))
    # print("dataG: ", data_G[0])
    # print("dataGF: ", data_GF[0])
    # print(f"data[{1}] = {data[1]}")

    parseForce = ParseForceFormBeam()
    i = 1
    Fx,Fy,_ = parseForce.getLiftingLoad(data[i])
    FFx = 130000*np.sin(i*10*np.pi/180)
    FFy = 130000*np.cos(i*10*np.pi/180)
    # print(f"FFx: {FFx:.2f}, FFy: {FFy:.2f}")
    # print(f"dfx: {FFx-Fx:.1f} : {(FFx-Fx)/FFx*100:.1f}%;  dfy: {FFy+Fy:.1f} : {(FFy+Fy)/FFy*100:.1f}%")
    ee = list(data[i][9:])
    ee.append(data[i][8])
    L5 = 3.114
    q45 = 33.02/180*np.pi
    A2 = [0.01856, 0.00496]
    # Fx,Fy = getLiftingLoad(ee, q45, L5, A2)
    # print(f"ffy: {FFy*L1/L5:.1f}; Fy: {Fy:.1f};  df: {FFy*L1/L5 + Fy:.1f}")