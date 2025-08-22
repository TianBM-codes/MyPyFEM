from cal_stress import *
from Component_FEA import *
import numpy as np
import math
import pymysql
from datetime import datetime
import pytz
import logging
import log_config  # 导入日志配置模块
from component_max_stress import update_component_max_stress
# 臂架、象鼻梁、大拉杆的旋转角度为常规坐标系，即x轴正方向为0度
# 综合应力计算过程中，定义y轴-方向为0°，逆时针方向为-方向，顺时针方向为+方向
logger = logging.getLogger()


def get_latest_angel(config, struct_id):
    """
    根据 struct_id 获取最新的 value2 作为实时角度
    struct_id 与表名对应关系：
        1001 -> inca37
        1002 -> inca40
        1003 -> inca45
    """
    table_map = {
        1001: "inca37",
        1002: "inca40",
        1003: "inca45"
    }

    # 检查 struct_id 是否有效
    if struct_id not in table_map:
        return None

    table_name = table_map[struct_id]

    try:
        conn = pymysql.connect(**config)
        cursor = conn.cursor()

        sql = f"""
        SELECT RV1
        FROM {table_name}
        WHERE RV1 IS NOT NULL
        ORDER BY T DESC
        LIMIT 1
        """
        cursor.execute(sql)
        record = cursor.fetchone()

        if not record:
            return None

        return float(record[0])

    except Exception as e:
        return None
    finally:
        if 'conn' in locals() and conn.open:
            conn.close()

def get_latest_load(config, struct_id):
    """
    根据 struct_id 获取最新的 RV1 值作为实时负载
    struct_id与表名对应关系：
        1001 -> cza133
        1002 -> cza134
        1003 -> cza132
    """
    table_map = {
        1001: "cza133",
        1002: "cza134",
        1003: "cza132"
    }

    if struct_id not in table_map:
        logger.error(f"无效的 struct_id: {struct_id}")
        return None

    table_name = table_map[struct_id]

    try:
        conn = pymysql.connect(**config)
        cursor = conn.cursor()

        # 用参数化占位防止 SQL 注入
        sql = f"""
        SELECT RV2
        FROM {table_name}
        WHERE RV1 IS NOT NULL
        ORDER BY T DESC
        LIMIT 1
        """
        cursor.execute(sql)
        record = cursor.fetchone()

        if not record:
            return None

        return float(record[0])

    except Exception as e:
        logger.error(f"读取实时负载错误 (struct_id={struct_id}): {e}", exc_info=True)
        return None
    finally:
        if 'conn' in locals() and conn.open:
            conn.close()

def insert_sigma_values(config, sigma_rod1, sigma_rod2, sigma_rod3, sigma_rod4, sigma_rod5,
                        sigma_arm, sigma_rzj, sigma_big_rod, struct_id, sigma_zhongzz=30,
                        rod1_vertical=15, rod1_horizontal=10, rod2_vertical=50):
    connection = None
    try:
        # 获取当前北京时间
        tz = pytz.timezone('Asia/Shanghai')
        T_now = datetime.now(tz).strftime('%Y-%m-%d %H:%M:%S.%f')[:-3]

        # 建立数据库连接
        conn = pymysql.connect(**config)
        cursor = conn.cursor()

        # 插入 SQL
        sql = """
            INSERT INTO t_data_monitoring_info (
                T, value1, value2, value3, value4, value5, value6, 
                value7, value8, value9, value10, value11, value12, struct_id
            ) 
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
        """
        cursor.execute(sql, (
            T_now,  # T 北京时间
            sigma_zhongzz,  # value1 中轴柱
            sigma_rzj,  # value2 人字架
            sigma_big_rod,  # value3 大拉杆
            sigma_arm,  # value4 臂架
            sigma_rod1,  # value5 杆一
            sigma_rod2,  # value6 杆二
            sigma_rod3,  # value7 杆三
            sigma_rod4,  # value8 杆四
            sigma_rod5,  # value9 杆五
            rod1_vertical,  # value10 杆一竖向应力脉动幅值
            rod1_horizontal,  # value11 杆一横向应力脉动幅值
            rod2_vertical,  # value12 杆二竖向应力脉动幅值
            struct_id  # struct_id
        ))

        conn.commit()
        print("实时最大应力数据插入成功")


    except Exception as e:
        print(f"插入健康指数失败: {e}")
    finally:
        if 'conn' in locals() and conn.open:
            conn.close()

def insert_stability_analysis(config, sigma_rod1_Stability, sigma_rod3_Stability,
                               sigma_rod5_Stability, sigma_arm_Stability, struct_id):
    connection = None
    try:
        # 获取当前北京时间
        tz = pytz.timezone('Asia/Shanghai')
        T_now = datetime.now(tz).strftime('%Y-%m-%d %H:%M:%S.%f')[:-3]

        # 建立数据库连接
        conn = pymysql.connect(**config)
        cursor = conn.cursor()

        # 插入 SQL
        sql = """
            INSERT INTO t_data_stability_analysis (
                T, struct_id, value1, value2, value3, value4
            ) 
            VALUES (%s, %s, %s, %s, %s, %s)
        """
        cursor.execute(sql, (
            T_now,                    # T 北京时间
            struct_id,
            sigma_arm_Stability,
            sigma_rod1_Stability,
            sigma_rod3_Stability,
            sigma_rod5_Stability,
        ))

        conn.commit()
        print("稳定性分析计算结果插入成功")

    except Exception as e:
        print(f"插入稳定性分析计算结果失败: {e}")
    finally:
        if 'conn' in locals() and conn.open:
            conn.close()

# 为了仿真分析插入0值
def insert_sim_data_zero(config, struct_id):

    connection = None
    try:
        # 获取当前北京时间
        tz = pytz.timezone('Asia/Shanghai')
        T_now = datetime.now(tz).strftime('%Y-%m-%d %H:%M:%S.%f')[:-3]

        # 建立数据库连接
        conn = pymysql.connect(**config)
        cursor = conn.cursor()

        # 插入 SQL
        sql = """
            INSERT INTO t_data_sim_analysis (
                T, value1, value2, value3, value4, value5, value6, value7, value8, struct_id, fea_cal_time
            ) 
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
        """
        cursor.execute(sql, (
            T_now, 0, 0, 0, 0, 0, 0, 0, 0, struct_id, T_now
        ))

        conn.commit()

    except Exception as e:
        print(f"插入失败: {e}")
    finally:
        if 'conn' in locals() and conn.open:
            conn.close()

def insert_max_stress_zero(config, struct_id):
    """
    向 t_data_stress_analysis 表中插入4行应力全为0的记录
    mid 分别为: 人字架、大拉杆、臂架、象鼻梁
    T 为当前北京时间
    struct_id 由函数参数传入
    """
    try:
        # 获取当前北京时间
        tz = pytz.timezone('Asia/Shanghai')
        T_now = datetime.now(tz).strftime('%Y-%m-%d %H:%M:%S.%f')[:-3]

        # 建立数据库连接
        conn = pymysql.connect(**config)
        cursor = conn.cursor()

        # 插入 SQL 模板
        sql = """
            INSERT INTO t_data_stress_analysis (
                T, mid, compressive_stress, tensile_stress, struct_id
            ) 
            VALUES (%s, %s, %s, %s, %s)
        """

        # 需要插入的 mid 列表
        mid_list = ["人字架", "大拉杆", "臂架", "象鼻梁"]

        for mid in mid_list:
            print(f"插入零值记录: struct_id={struct_id}, mid={mid}, T={T_now}")
            cursor.execute(sql, (T_now, mid, 0, 0, struct_id))

        conn.commit()
        print(f"struct_id={struct_id} 零值应力记录插入完成")

    except Exception as e:
        print(f"插入失败: {e}")
    finally:
        if 'conn' in locals() and conn.open:
            conn.close()


# 记录每台起重机是否已经插过 0
last_zero_inserted_map = {
    1001: False,
    1002: False,
    1003: False
}
def current_max_stress(config, struct_id, data):
    """
    data:  中轴柱极值应力，杆1、杆2脉动应力
    根据 struct_id 获取数据并插入应力与稳定性信息。
    负载 < 0.02 时：只插一次 0
    负载 >= 0.02 时：计算应力并插入
    """
    global last_zero_inserted_map

    diao_G = get_latest_load(config, struct_id)
    arm_angle_deg = get_latest_angel(config, struct_id)

    # 如果数据无效，跳过
    if diao_G is None or arm_angle_deg is None:
        logger.info(f"[{struct_id}] 无法获取负载或角度数据，所有计算跳过")
        return

    # 负载过小
    if diao_G < 0.02:
        if not last_zero_inserted_map[struct_id]:
            insert_sigma_values(
                config,
                0, 0, 0, 0, 0,   # sigma_rod1..sigma_rod5
                0, 0, 0, struct_id, data[0], data[1], data[2], data[3]          # sigma_arm, sigma_rzj, sigma_big_rod
            )
            insert_stability_analysis(config, 0, 0, 0, 0, struct_id)
            insert_sim_data_zero(config, struct_id)
            insert_max_stress_zero(config, struct_id)
            last_zero_inserted_map[struct_id] = True
            logger.info(f"[{struct_id}] 负载 {diao_G:.3f} < 0.02，实时最大应力插入一组 0")
            logger.info(f"[{struct_id}] 负载 {diao_G:.3f} < 0.02，稳定性计算插入一组 0")
            logger.info(f"[{struct_id}] 负载 {diao_G:.3f} < 0.02，实时仿真最大应力计算插入一组 0")
            logger.info(f"[{struct_id}] 负载 {diao_G:.3f} < 0.02，部件极值应力计算插入一组 0")
        else:
            print(f"[{struct_id}] 负载 {diao_G:.3f} < 0.02，已插入过 0，跳过")
        return

    # 正常吊重
    (sigma_rod1, sigma_rod2, sigma_rod3, sigma_rod4, sigma_rod5, sigma_arm,
     sigma_rzj, sigma_big_rod, sigma_rod1_Stability, sigma_rod3_Stability,
     sigma_rod5_Stability, sigma_arm_Stability) = calculate_stress(diao_G, arm_angle_deg)

    insert_sigma_values(
        config,
        sigma_rod1, sigma_rod2, sigma_rod3, sigma_rod4, sigma_rod5,
        sigma_arm, sigma_rzj, sigma_big_rod, struct_id, data[0], data[1], data[2], data[3]
    )
    last_zero_inserted_map[struct_id] = False
    logger.info(f"[{struct_id}] 负载 {diao_G:.3f} 吨，实时最大应力已插入计算值")
    insert_stability_analysis(config, sigma_rod1_Stability, sigma_rod3_Stability,
                               sigma_rod5_Stability, sigma_arm_Stability, struct_id)
    logger.info(f"[{struct_id}] 负载 {diao_G:.3f} 吨，稳定性计算已插入计算值")
    insert_sim_data(config, struct_id)
    logger.info(f"[{struct_id}] 负载 {diao_G:.3f} 吨，仿真计算最大应力已插入计算值")
    update_component_max_stress(config, struct_id)
    logger.info(f"[{struct_id}] 负载 {diao_G:.3f} 吨，部件极值应力已更新最大值")




def main():
    """一次性处理三台起重机"""
    config = {
        'host': 'monitor.sipesc.net',
        'port': 30053,
        'user': 'nbport_user',
        'password': 'P@jtMUS3Yx',
        'database': 'nbport_qzj_db',
        'charset': 'utf8mb4'
    }
    for struct_id in (1001, 1002, 1003):
        current_max_stress(config, struct_id, [0,0,0,0])




if __name__ == "__main__":
    main()