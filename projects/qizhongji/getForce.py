import numpy as np 
## 以下两个库用于测试
import matplotlib.pyplot as plt
## 导入测试数据
import pathlib

class ForceParser:
    FBaseData = []
    GBaseData = []
    ABJ = np.arange(4322,7822,50)/100
    def __init__(self):
        file_path = pathlib.Path(__file__).parent / 'F10000N_res.txt'
        with open(file_path, 'r', encoding='utf-8') as f:
            for line in f:
                # 去掉换行符，按每 10 位切片
                row = [float(line[i:i+10]) for i in range(0, len(line.strip()), 10)]
                self.FBaseData.append(row)
        self.FBaseData = np.array(self.FBaseData)
        # 计算每一行的模
        self.FBaseData_norms = np.linalg.norm(self.FBaseData, axis=1, keepdims=True)
        # 将每一行除以其模，得到单位向量
        self.FBaseData_normalized = self.FBaseData / self.FBaseData_norms
        
        file_path = pathlib.Path(__file__).parent / 'G10_res.txt'
        with open(file_path, 'r', encoding='utf-8') as f:
            for line in f:
                # 去掉换行符，按每 10 位切片
                row = [float(line[i:i+10]) for i in range(0, len(line.strip()), 10)]
                self.GBaseData.append(row)
        self.GBaseData = np.array(self.GBaseData)

    def getForce(self, Abj, sensorData):
        ## Abj 臂架角度；
        ## sensorData 传感器数据，去除初始应变、重力应变后的纯载荷应变
        
        try:
            idx = np.abs(self.ABJ - Abj).argmin()  ## 获取ABJ中最接近于Abj的元素位置
            normV = np.linalg.norm(sensorData)
            F = 10000 * normV / self.FBaseData_norms[idx] 
            return F 

        except Exception as e:
            print(f"in ForceParser getForce err: {e}")
        
    '''
    fd = np.array([np.max(FBaseData[:,i])- np.min(FBaseData[:,i]) for i in range(25)])
    gd = np.array([np.max(GBaseData[:,i])- np.min(GBaseData[:,i]) for i in range(25)])
    for i in range(25):
        print(f"{i} --> F: {fd[i]:8.2f}  G: {gd[i]:8.2f}")
    ff = fd/np.linalg.norm(fd)
    gg = gd/np.linalg.norm(gd)
    print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
    for i in range(25):
        print(f"{i} --> F: {ff[i]:8.2f}  G: {gg[i]:8.2f}")
    flag = False
    while flag:
        id = input("input: ").strip()
        if(id == "exit"):
            flag = False
            break

        plt.plot(GBaseData[...,int(id)])
        plt.show()
    '''
# s0 = np.array([0.2,0.2,  -0.7,  -0.8,   2.7,  -2.5,  -2.5,  -2.5,  3.7 , -0.4,  -0.4 , -0.3,   0.6,   0.5,  -0.4 , -0.3 , -0.1 , -0.3])
# s1 = np.array([21.3,  23.2, -96.6,-104.1, 350.8,-325.0,-324.7,-326.9, 479.4, -44.7, -47.7, -42.7,  74.1,  68.5, -52.7,  -1.5, -52.3,-104.7])

# ss1 = s1/np.linalg.norm(s1)
# ss0 = s0/np.linalg.norm(s0)
# for i in range(len(s1)):
#     print(f"{s0[i]*130:.3f} -- {s1[i]:.3f}")
# import math
# print(math.acos(np.dot(ss0,ss1)))