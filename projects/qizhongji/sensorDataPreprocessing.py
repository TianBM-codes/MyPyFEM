import sys
import pathlib
current_file = pathlib.Path(__file__).resolve()
# 所在目录
current_dir = current_file.parent
# 假设 projects 文件夹在当前脚本的上一级目录
sys.path.append(str(current_dir))
sys.path.append(str(current_dir.parent.parent))
sys.path.append(str(current_dir/"StressCalculate"))
import mysql
from db_ops import query_data, save_data

# from MyPyFEM_fatigue.ioclass.MySqlMathFunction import MySQLMathFunction
from QiZhongJiZiTai import MQ1330
import StressCalculate.data_insert as stressCal

import parseForce  ## 材料力学内里计算方法识别载荷
import getForce   ## 投影方法识别载荷
import numpy as np
import pandas as pd
from scipy import signal   ## 用来计算脉动峰值，同时去趋势项

import threading
import time
import datetime

# mysql_db = MySQLMathFunction()
# sql = "SELECT T,RV1 FROM nbport_qzj_db.INCA37 ORDER BY T DESC LIMIT 10;"
# bjA = mysql_db.execute_sql(sql)
# print(bjA)

class DataPreProcessing():
    # qi_zhong_ji = {"host": "192.168.3.8", "port": 30053, "database": "nbport_qzj_db",
            #    "username": "nbport_user", "password": "P@jtMUS3Yx"}

    mysql_db = None
    mq1330 = None
    GBaseData = None
    ABJ = None
    parseBeamForce = None
    parseForceByNorm = None
    con = None
    def __init__(self):
        # self.con = mysql.connector.connect(host=self.qi_zhong_ji["host"],
        #                                        port=self.qi_zhong_ji["port"],
        #                                        database=self.qi_zhong_ji["database"],
        #                                        user=self.qi_zhong_ji["username"],
        #                                        password=self.qi_zhong_ji["password"]
        #                                        )     

        # self.mysql_db = MySQLMathFunction()
        self.mq1330 = MQ1330()
        self.GBaseData = []
        self.ABJ = np.arange(4322,7822,50)/100

        file_path = pathlib.Path(__file__).parent / 'G10_res.txt'
        with open(file_path, 'r', encoding='utf-8') as f:
            for line in f:
                # 去掉换行符，按每 10 位切片
                row = [float(line[i:i+10]) for i in range(0, len(line.strip()), 10)]
                self.GBaseData.append(row)
        self.GBaseData = np.array(self.GBaseData)

        self.parseBeamForce = parseForce.ParseForceFormBeam()
        self.parseForceByNorm = getForce.ForceParser()


    # def __del__(self):
    #     # ---------- 5. 关闭连接 ----------
    #     self.con.close()

    def updateWrokRangeByBJA(self):
        pass 
        ## 获取当前角度,1s数据，即10条数据
        sql = "SELECT T,RV1 FROM nbport_qzj_db.INCA37 ORDER BY T DESC LIMIT 10;"
        bjA = self.mysql_db.execute_sql(sql)
        t =[]
        Adata =[]
        workD = []
        H = 0
        for i in bjA:
            t.append(i[0].timestamp())
            Adata.append(i[1])
            self.mq1330.getPartOrientation(i[1])
            workD.append(self.mq1330.getPartCircle()[0])
            H = self.mq1330.getPartCircle()[1]
        dt = [t[i+1]-t[i] for i in range(len(t)-1)]
        dd = [workD[i+1] - workD[i] for i in range(len(t)-1)]
        workspeed = [dd[i]/dt[i] for i in range(len(dt))]
        speed = np.mean(workspeed)

        sql = (f"INSERT INTO t_work_status (T, value1, value5) "
           f"VALUES ('{bjA[-1][0]}', '{workD[-1] / 1000:.2f}', {H / 1000:.2f}) "
           "ON DUPLICATE KEY UPDATE "
           "value1= VALUES(value1), "
           "value5= VALUES(value5);"
           )
        print(sql)
        self.mysql_db.commit_sql(sql)

        sql = (f"INSERT INTO t_data_work_status_info_1001 (T, value7, value8) "
           f"VALUES ('{bjA[-1][0]}', '{workD[-1] / 1000:.2f}', {speed / 1000:.2f}) "
           "ON DUPLICATE KEY UPDATE "
           "value7= VALUES(value7), "
           "value8= VALUES(value8);"
           )
        print(sql)
        self.mysql_db.commit_sql(sql)

    def sectionStrain(self, data):
        ## 计算中轴柱的极值应变
        ## data: 中轴柱的6个应变值
        ## return max, min   即最大值和最小值
        try:
            # 定义 4 个坐标点
            points =[]
            r = 1.41   ## 中轴柱半径
            for i in range(6):
                points.append([r*np.cos(i*np.pi/3), r*np.sin(i*np.pi/3), data[i]])
            # print("point: ", points)
            points = np.array(points)
            if(np.isnan(points).any() or np.isinf(points).any()):
                print("in sectionStrain data is nan or inf")
                return 0,0
            
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

            # z = -(ax + by + d)/c
            x, y, _ = normal_vector*r   ## 法向量在圆周上的坐标
            mx = -(-a*x + -b*y + d)/c
            mi = -(a*x + b*y + d)/c
            return mx, mi
        except Exception as e:
            print(f"in sectionStrain err: {e}")


    def has_nan_inf(self, df: pd.DataFrame, num) -> bool:
        """
        返回 True 表示 DataFrame 中至少有一个 NaN 或 ±inf
        """
        cnt = df['sensor'].nunique(dropna=True)  ## 计算传感器个数
        if(cnt != num):
            return True 

        # 1. NaN 检测：对所有 dtype 都适用
        nan_mask = df.isna()
        # 2. inf 检测：只对数值列，跳过 object/datetime
        numeric_cols = df.select_dtypes(include=[np.number]).columns
        # numeric_cols = df["C_RV1"]
        inf_mask = pd.DataFrame(False, index=df.index, columns=df.columns)
        if not numeric_cols.empty:
            inf_mask[numeric_cols] = np.isinf(df[numeric_cols])
        return df.empty or (nan_mask | inf_mask).any().any()

    num = 0    
    def dataProcessing(self, tables, loadtime):
        ## loadtime: 吊重时间戳，即初始应变时间戳
        ## 从数据库读取最新数据，计算识别的载荷Fx，Fy，My（侧摆力矩），杆1竖向脉动应力，杆1横向脉动应力，杆2竖向脉动应力
        try:
            num_yb = 20  # 读取应变数据条数，跟采样率20Hz匹配
            num_qj = 10  # 读取倾角数据条数，跟采样率10Hz匹配
            num_cz = 5  # 读取称重数据条数，跟采样率1Hz匹配

            # ****************************************  构造动态 SQL ********************************************
            # 每张表：先过滤 T>t0，再按 id 倒序取 50 条
            sql_parts = [
                f"(SELECT %s AS sensor, T, C_RV1 FROM {tbl} WHERE T > %s ORDER BY T ASC LIMIT {num_yb})"
                for tbl in tables[0]
            ]
            sql_yb = " UNION ALL ".join(sql_parts)


            ## 臂架倾角，人字架倾角
            sql_parts = []
            sql_parts.append(f"(SELECT %s AS sensor, T, RV1 FROM {tables[1][0]} WHERE T > %s ORDER BY T ASC LIMIT {num_qj})")
            sql_parts.append(f"(SELECT %s AS sensor, T, RV1 FROM {tables[1][1]} WHERE T > %s ORDER BY T ASC LIMIT {num_qj})")
            sql_qj = " UNION ALL ".join(sql_parts)

            ## 称重
            sql_cz =f"(SELECT %s AS sensor, T, RV2, RV4, RV7 FROM {tables[2][0]} WHERE T > %s ORDER BY T ASC LIMIT {num_cz})"


            ## 
            # ****************************************  读初始应变 ********************************************
            # 参数排列顺序：sensor1, t0, sensor2, t0, ...
            print("in dataProcessing loadtime: ", loadtime)
            params = [p for tbl in tables[0] for p in (tbl, loadtime)]
            # self.mysql_db.cursor.execute(sql_yb,params)        
            # columns = [desc[0] for desc in self.mysql_db.cursor.description]
            # df0_yb = pd.DataFrame(self.mysql_db.cursor.fetchall(), columns=columns)
            df0_yb = query_data((sql_yb,params))

            ## 检查df0_yb是否存在nan，inf，数量足够，等
            has_bad = self.has_nan_inf(df0_yb, len(tables[0]))
            print("has_bad: ", has_bad)
            for i in range(15):
                if has_bad:
                    print(f"df0_yb has bad, sleep times: {i}")
                    time.sleep(3)  ## 当df_yb是空时，即表示数据还不齐，则延时2s，等待
                    df0_yb = query_data((sql_yb,params))
                    has_bad = self.has_nan_inf(df0_yb, len(tables[0]))
                else:
                    break
            if has_bad:  ## 15次延时后还是空，这直接返回
                print("df0_yb has bad! return! ")
                return


            params = [p for tbl in tables[1] for p in (tbl, loadtime)]
            # self.mysql_db.cursor.execute(sql_qj,params)
            # columns = [desc[0] for desc in self.mysql_db.cursor.description]
            # df0_qj = pd.DataFrame(self.mysql_db.cursor.fetchall(), columns=columns)
            df0_qj = query_data((sql_qj,params))

            ## 检查df0_qj是否存在nan，inf，数量足够，等
            has_bad = self.has_nan_inf(df0_qj, len(tables[1]))
            print("has_bad: ", has_bad)
            for i in range(15):
                if has_bad:
                    print(f"df0_qj has bad, sleep times: {i}")
                    time.sleep(3)  ## 当df_yb是空时，即表示数据还不齐，则延时2s，等待
                    df0_qj = query_data((sql_qj,params))
                    has_bad = self.has_nan_inf(df0_qj, len(tables[1]))
                else:
                    break
            if has_bad:  ## 15次延时后还是空，这直接返回
                print("df0_qj has bad! return! ")
                return


            yb0 = np.array([np.mean(df0_yb[df0_yb['sensor'] == tname][["C_RV1"]].values.ravel()[:10]) for tname in tables[0]])  ## 取初始时刻后的10个数据平均值作为应变初值
            Abj0 = df0_qj[df0_qj['sensor'] == tables[1][0]][['RV1']].values.ravel()[0]

            idx0 = np.abs(self.ABJ - Abj0).argmin()  ## 获取ABJ中最接近于Abj的元素位置
            ybG0 = self.GBaseData[idx0]/10   ## 初始状态的重力应变

            # ****************************************  读当前应变 ********************************************
            ## 读应变
            sqlT = f"SELECT T FROM {tables[0][0]} ORDER BY T DESC LIMIT {num_yb};"
            # self.mysql_db.cursor.execute(sqlT)
            # df = self.mysql_db.cursor.fetchall()
            df = query_data(sqlT)
            if df.empty:
                # raise ValueError("查询当前应变结果为空，无法获取T最新值")
                print("df is empty!!!")
                return 
            ti = df["T"].values[-1]      # 或 df["T"].values[-1]
            # print("df: ", df)
            # ti = df["T"].values.ravel()[-1]
            # ti += datetime.timedelta(seconds=self.num) ## 当前状态时刻  np.timedelta64(self.num, 's')  #
            ti = pd.to_datetime(ti).to_pydatetime()
            print("ti: ", ti)
            if(ti<loadtime):
                time.sleep(2)    ## 当前时间ti < loadtime时，延时2s，等待
                # raise ValueError("查询到的当前时间ti < loadtime，无法获取最新值") 
                print("ti < loadtime!! return!")
                return 
            self.num +=1
            # print("ti: ", ti)
            # ti -= datetime.timedelta(seconds=30)  ## 向后延30s，避免时间同步问题
            # ti = datetime.datetime(2025, 8, 21, 16, 50, 6, 375)   ## 2025-08-21 16:50:06.375000
            # print("tii: ", ti)
            # 参数排列顺序：sensor1, t0, sensor2, t0, ...
            params = [p for tbl in tables[0] for p in (tbl, ti)]
            # self.mysql_db.cursor.execute(sql_yb,params)
            # columns = [desc[0] for desc in self.mysql_db.cursor.description]
            # df_yb = pd.DataFrame(self.mysql_db.cursor.fetchall(), columns=columns)
            df_yb = query_data((sql_yb,params))

            ## 检查df_yb是否存在nan，inf，数量足够，等
            has_bad = self.has_nan_inf(df_yb, len(tables[0]))
            # cnt = df_yb.count()["sensor"]
            print("has_bad: ", has_bad)
            for i in range(15):
                if has_bad:
                    print(f"df_yb has bad, sleep times: {i}")
                    time.sleep(3)  ## 当df_yb是空时，即表示数据还不齐，则延时2s，等待
                    df_yb = query_data((sql_yb,params))
                    has_bad = self.has_nan_inf(df_yb, len(tables[0]))
                else:
                    break
            
            if has_bad:  ## 15次延时后还是空，这直接返回
                print("df_yb has bad! return! ")
                return
            # print("params: ", params)
            # print("df_yb: ", df_yb)


            params = [p for tbl in tables[1] for p in (tbl, ti)]
            # self.mysql_db.cursor.execute(sql_qj,params)
            # columns = [desc[0] for desc in self.mysql_db.cursor.description]
            # df_qj = pd.DataFrame(self.mysql_db.cursor.fetchall(), columns=columns)
            df_qj = query_data((sql_qj,params))

            ## 检查df_qj是否存在nan，inf，数量足够，等
            has_bad = self.has_nan_inf(df_qj, len(tables[1]))
            # cnt = df_qj.count()["sensor"]
            print("has_bad: ", has_bad)
            for i in range(15):
                if has_bad:
                    print(f"df_qj has bad, sleep times: {i}")
                    time.sleep(3)  ## 当df_yb是空时，即表示数据还不齐，则延时2s，等待
                    df_qj = query_data((sql_qj,params))
                    has_bad = self.has_nan_inf(df_qj, len(tables[1]))
                else:
                    break
            if has_bad:  ## 15次延时后还是空，这直接返回
                print("df_qj has bad!! return!")
                return


            # 取出 应变值
            yb = np.array([np.mean(df_yb[df_yb['sensor'] == tname][["C_RV1"]].values.ravel()) for tname in tables[0]])   ## 当前时刻ti，1s内的均值作为计算值
            Abj = np.mean(df_qj[df_qj['sensor'] == tables[1][0]][["RV1"]].values.ravel())  ## 当前时刻t0，1s内的均值作为计算值
            Arz = np.mean(df_qj[df_qj['sensor'] == tables[1][1]][["RV1"]].values.ravel())  ## 当前时刻t0，1s内的均值作为计算值

            # ****************************************  中轴柱极值应变 ********************************************
            ##  计算中轴柱极值应变
            if(len(yb)>25):
                # print("yb: ", yb[-6:] - yb0[-6:])
                yb_zz = yb[-6:] - yb0[-6:]
                if(np.isnan(yb_zz).any() or np.isinf(yb_zz).any()):
                    # print("df_yb: ", df_yb.to_string())
                    print("in dataProcessing zhong zhou zhu YB data is nan or inf;  data: ", yb_zz)
                    print(f"yb[-6:]; yb = {yb},   yb0[-6:]; yb0 = {yb0}")
                    zzMaxMin = (0,0)
                else:
                    zzMaxMin = self.sectionStrain(yb_zz)
            else:
                zzMaxMin = (0,0)

            # ****************************************  计算象鼻梁载荷和侧摆力矩 ********************************************
            try:
                idx = np.abs(self.ABJ - Abj).argmin()  ## 获取ABJ中最接近于Abj的元素位置
                ybG = self.GBaseData[idx]/10   ## 当前状态的重力应变
                dt_ybG = ybG - ybG0  ## 重力应变的变化量

                ## 传感器索引序列  顺序[左，左上，右上，右， 杆2]
                indexSensor = [2, 3, 1, 0, 4]    ### #6号机的排序

                yb_load = yb[indexSensor] - yb0[indexSensor] - dt_ybG[indexSensor]  ## 纯结构载荷应变
                # print("yb_load: ", yb_load)
                ## Fx: X向集中力， Fy：Y向集中力， My：侧摆力矩
                try:
                    Fx,Fy, My = self.parseBeamForce.getLiftingLoad(yb_load)      ## 载荷应变 = 当前应变 - 初始状态应变 - 重力应变变化量
                    Fx /=10000  ## 单位 吨
                    Fy /=10000  ## 单位 吨
                    My /=10000  ## 单位 吨米
                    # print(f"Fx: {Fx}, Fy: {Fy}, My: {My}")
                except Exception as err:
                    print(f"in dataProcessing parseBeamForce.getLiftingLoad err: {err}")
                try:
                    F = self.parseForceByNorm.getForce(Abj, yb_load)[0]  ## 投影法识别载荷
                    F /= 10000   ## 单位 吨
                    # print(f"F: {F}")
                except Exception as err:
                    print(f"in dataProcessing parseForceByNorm.getForce err: {err}")

            except Exception as ee:
                print(f"in dataProcessing get XBL load err: {ee}")

            # ****************************************  计算象鼻梁脉动应变 ********************************************

            ## 计算脉动应力
            cnt = num_yb*5  ## 取5秒的数据分析
            sql = f"(SELECT '{tables[0][indexSensor[1]]}' AS sensor, T, C_RV1 FROM {tables[0][indexSensor[1]]} ORDER BY T DESC LIMIT {cnt})  UNION ALL "\
                  f"(SELECT '{tables[0][indexSensor[0]]}' AS sensor, T, C_RV1 FROM {tables[0][indexSensor[0]]} ORDER BY T DESC LIMIT {cnt})  UNION ALL "\
                  f"(SELECT '{tables[0][indexSensor[3]]}' AS sensor, T, C_RV1 FROM {tables[0][indexSensor[3]]} ORDER BY T DESC LIMIT {cnt})  UNION ALL "\
                  f"(SELECT '{tables[0][indexSensor[4]]}' AS sensor, T, C_RV1 FROM {tables[0][indexSensor[4]]} ORDER BY T DESC LIMIT {cnt});"

            # self.mysql_db.cursor.execute(sql)
            # columns = [desc[0] for desc in self.mysql_db.cursor.description]
            # df = pd.DataFrame(self.mysql_db.cursor.fetchall(), columns=columns)
            df = query_data(sql)
            if df.empty:
                # raise ValueError("查询象鼻梁当前应变结果为空，无法获取最新值")
                return 

            yb1 = df[df['sensor'] == tables[0][indexSensor[1]]][["C_RV1"]].values.ravel()
            yb2 = df[df['sensor'] == tables[0][indexSensor[0]]][["C_RV1"]].values.ravel()
            yb3 = df[df['sensor'] == tables[0][indexSensor[3]]][["C_RV1"]].values.ravel()
            yb4 = df[df['sensor'] == tables[0][indexSensor[4]]][["C_RV1"]].values.ravel()

            ## 杆1 竖向脉动应力
            # 去除直流分量和趋势（推荐一步完成）
            yb1_clean = signal.detrend(yb1)  # 自动去均值 + 线性趋势
            yb1Ypp = yb1_clean.max() - yb1_clean.min()

            ## 杆1 横向脉动应力
            yb1_clean = signal.detrend(yb2 - yb3)  # 自动去均值 + 线性趋势
            yb1Xpp = yb1_clean.max() - yb1_clean.min()

            ## 杆2 竖向脉动应力
            yb1_clean = signal.detrend(yb4)  # 自动去均值 + 线性趋势
            yb2Ypp = yb1_clean.max() - yb1_clean.min()

            ## **************************************** 吊重和变幅速度 ********************************************
            # 参数排列顺序：sensor1, t0, sensor2, t0, ...
            tii = ti - datetime.timedelta(seconds=10)
            params = [p for tbl in tables[2] for p in (tbl, tii)]
            # self.mysql_db.cursor.execute(sql_cz,params)
            # columns = [desc[0] for desc in self.mysql_db.cursor.description]
            # df = pd.DataFrame(self.mysql_db.cursor.fetchall(), columns=columns)
            df = query_data((sql_cz,params))
            if df.empty:
                # raise ValueError("查询称重结果为空，无法获取最新值")
                return 

            ## 取出重量
            zaihe_sz = np.mean(df[df['sensor'] == tables[2][0]][["RV2"]].values.ravel())  ## 实重，当前时刻t0，1s内的均值作为计算值
            zaihe_pz = np.mean(df[df['sensor'] == tables[2][0]][["RV7"]].values.ravel())  ## 皮重，当前时刻t0，1s内的均值作为计算值
            fudu = df[df['sensor'] == tables[2][0]][["RV4"]].values.ravel()  ## 幅度
            # print(f"zaihe1: {zaihe_sz}, zaihe_pz: {zaihe_pz}, fudu: {fudu}")

            dd = [fudu[i+1] - fudu[i] for i in range(len(fudu)-1)]   ## 采样率为1Hz
            workspeed = np.mean(dd)  ## 取平均


            ## *************************************** 保存数据到数据库  ******************************************** 
            data_to_insert = (ti, float(np.round(My,2)), float(np.round(Arz,3)), 19.00, float(np.round(Abj,3)), float(np.round(fudu[-1],2)), float(np.round(workspeed,2)), float(np.round(zaihe_sz+zaihe_pz,2)), float(np.round(F,2)), tables[3])

            # 插入语句
            # insert_query = "INSERT INTO t_data_work_status_info (T, value1, value2, value3, value7, value8, value9, value10) VALUES (%s, %s, %s, %s, %s, %s, %s, %s)"
            insert_query = f"""INSERT INTO t_data_work_status_info 
                                (T, value1, value2, value3, value4, value7, value8, value9, value10, struct_id) 
                                VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s) 
                                ON DUPLICATE KEY UPDATE 
                                    value1 = VALUES(value1),
                                    value2 = VALUES(value2),
                                    value3 = VALUES(value3),
                                    value4 = VALUES(value4),
                                    value7 = VALUES(value7),
                                    value8 = VALUES(value8),
                                    value9 = VALUES(value9),
                                    value10 = VALUES(value10),
                                    struct_id = VALUES(struct_id)
                                """
            # self.mysql_db.cursor.execute(insert_query, data_to_insert)
            # self.mysql_db.con.commit()
            save_data(insert_query, data_to_insert)


            ## ************************************* 处理应力数据 ********************************************************
            config = {
                        'host': '192.168.3.8',
                        'port': 30053,
                        'user': 'nbport_user',
                        'password': 'P@jtMUS3Yx',
                        'database': 'nbport_qzj_db',
                        'charset': 'utf8mb4'
                        }
            stressCal.current_max_stress(config, tables[3], [float(np.round(zzMaxMin[0],2)), float(np.round(yb1Ypp,2)), float(np.round(yb1Xpp,2)), float(np.round(yb2Ypp,2))])

                      
        except Exception as e:
            print(f"in dataProcessing err: {e}")
        
if __name__ == "__main__":
    
    
    ## #6号机
    tables6_yb = [f'RSGB{i}' for i in range(47, 74)]   # ['rsg1', ... , 'rsg10']
    tables6_yb.remove("RSGB55")
    tables6_qj = ["INCA37", "INCA38"]    ## 臂架倾角, 人字架倾角
    tables6_cz = ["cza133"]   ## 称重

    ## #7号机
    tables7_yb = [f'RSGB{i}' for i in range(103, 130)]   # ['rsg1', ... , 'rsg10']
    tables7_qj = ["INCA40", "INCA42"]    ## 臂架倾角, 人字架倾角
    tables7_cz = ["cza134"]   ## 称重

    ## #8号机
    tables8_yb = [f'RSGB{i}' for i in range(76, 97)]   # ['rsg1', ... , 'rsg10']
    tables8_yb.remove("RSGB89")
    tables8_yb.remove("RSGB90")
    tables8_qj = ["INCA45", "INCA46"]    ## 臂架倾角, 人字架倾角
    tables8_cz = ["cza132"]   ## 称重

    param6 = [tables6_yb, tables6_qj, tables6_cz, 1001]
    param7 = [tables7_yb, tables7_qj, tables7_cz, 1002]
    param8 = [tables8_yb, tables8_qj, tables8_cz, 1003]
    tic = time.time()
    dp1 = DataPreProcessing()
    dp2 = DataPreProcessing()
    dp3 = DataPreProcessing()
    
    
    

    for i in range(1):
        t0 = datetime.datetime(2025, 8, 20, 9, 49, 20+i, 000) #- datetime.timedelta(seconds=5)
        # dp.dataProcessing([tables6_yb, tables6_qj, tables6_cz, "t_data_work_status_info_1001"], t0)
        threads = []
        # threading.Thread(target=dp1.dataProcessing, args=(param6, t0), daemon=True).start()
        # threading.Thread(target=dp2.dataProcessing, args=(param7, t0), daemon=True).start()
        # threading.Thread(target=dp3.dataProcessing, args=(param8, t0), daemon=True).start()
        t = threading.Thread(target=dp1.dataProcessing, args=(param6, t0), name=f"TaskThread-{1}")
        threads.append(t)
        t.start()

        t = threading.Thread(target=dp2.dataProcessing, args=(param7, t0), name=f"TaskThread-{2}")
        threads.append(t)
        t.start()

        t = threading.Thread(target=dp3.dataProcessing, args=(param8, t0), name=f"TaskThread-{3}")
        threads.append(t)
        t.start()
        # 等待所有线程完成
        for t in threads:
            t.join()

        time.sleep(1)
        print(f"i: {i}, time: {np.round(time.time()-tic, 3)}")

    print(f"ii: {i}, time: {np.round(time.time()-tic, 3)}")