import numpy as np
from fatigue_analysis import *
import pymysql
from datetime import datetime, timedelta
import logging

# 获取日志记录器
logger = logging.getLogger()
def get_recent_rt_files(config):

    # 固定的参考时间（可自行修改）
    # reference_time_str = "2025-04-07 00:00:00"
    # reference_time = datetime.strptime(reference_time_str, "%Y-%m-%d %H:%M:%S")
    # 使用当前时间作为参考时间
    reference_time = datetime.now()

    # 计算参考时间前一周的周一（开始时间）和周日（结束时间）
    # 先取参考时间的日期部分
    ref_date = reference_time.date()
    # 计算参考时间当天是星期几（星期一是0，星期天是6）
    weekday = ref_date.weekday()

    # 上一周的周一 = 参考时间日期 - (weekday + 7) 天
    last_week_start = datetime.combine(ref_date - timedelta(days=weekday + 7), datetime.min.time())
    # 上一周的周日 = 上一周的周一 + 6 天，时间到当天最后一秒
    last_week_end = datetime.combine(last_week_start.date() + timedelta(days=6), datetime.max.time())



    query = """
        SELECT file_name
        FROM t_calculate_nephogram_test
        WHERE result_type = 'f'
          AND calculate_time >= %s
          AND calculate_time <= %s;
    """

    try:
        connection = pymysql.connect(**config)
        with connection.cursor() as cursor:
            cursor.execute(query, (last_week_start, last_week_end))
            results = cursor.fetchall()

        if not results:
            logger.info("查询结果为空：没有符合条件的文件。")
            return []

        # 直接返回完整路径并替换 .dat 为 _mises.dat
        return [row[0] for row in results]

    except Exception as e:
        logger.error("数据库查询出错", exc_info=True)
        return []

    finally:
        if 'connection' in locals() and connection.open:
            connection.close()

if __name__ == "__main__":

    config = {
        'host': '192.168.3.8',
        'port': 30053,
        'user': 'nbport_user',
        'password': 'P@jtMUS3Yx',
        'database': 'nbport_qzj_db',
        'charset': 'utf8mb4'
    }
    files = get_recent_rt_files(config)
    # 1.dat文件所存放的文件夹
    folder_path = r"D:\SiPESC\project\701\input\test_data"
    # 2. 配置材料参数
    Su = 600  # 材料抗拉强度 (MPa)
    material_params = {
        'S_endurance': 30,  # 疲劳极限 (MPa)
        'a': 1e12,  # Basquin方程参数a
        'b': -3  # Basquin方程参数b
    }
    life_array = analyze_dat_fatigue(folder_path, Su, material_params)
    # 3. 将短疲劳结果转换为长疲劳结果
    # long_fatigue = short_to_long(life_array, r"D:\SiPESC\project\701\input\MQ1330_remesh.npy")
