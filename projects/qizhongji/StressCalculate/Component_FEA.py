import pymysql
import logging
import log_config  # 导入日志配置模块
from fatigue_analysis import ReaddatAsNumpy
import numpy as np
from datetime import datetime, timedelta
import pytz  # timezone 支持

# 获取日志记录器
logger = logging.getLogger()


def insert_sim_data(config, struct_id):
    """
    从 t_calculate_nephogram 表中，找到距离当前时间最近、result_type='mises' 的一条记录，
    返回其 file_name 和 calculate_time
    """
    try:
        conn = pymysql.connect(**config)
        cursor = conn.cursor()

        # 1. 查询最新的计算结果
        query = """
        SELECT file_name, calculate_time
        FROM t_calculate_nephogram
        WHERE result_type = 'mises' AND struct_id = %s
        ORDER BY calculate_time DESC
        LIMIT 1
        """
        cursor.execute(query, (struct_id,))
        record = cursor.fetchone()

        if not record:
            return None

        file_name, calculate_time = record
        logger.info(f"应力组件计算找到最近的记录，时间: {calculate_time}，file_name: {file_name}")

        # 2. 检查是否已存在相同calculate_time的记录
        check_query = """
        SELECT COUNT(*) FROM t_data_sim_analysis 
        WHERE struct_id = %s AND fea_cal_time = %s
        """
        cursor.execute(check_query, (struct_id, calculate_time))
        count = cursor.fetchone()[0]

        if count > 0:
            logger.info(f"struct_id={struct_id}的计算结果(时间:{calculate_time})已存在，跳过处理")
            return None

        # 3. 读取应力数据
        mises = ReaddatAsNumpy(file_name)

        # 计算各部件最大值
        # 人字架最大值 & 索引
        renzijia_rel_idx = np.argmax(mises[0:13912])
        renzijia_idx = renzijia_rel_idx  # 偏移量为 0
        renzijia_max = round(float(mises[renzijia_idx]), 2)

        # 臂架最大值 & 索引
        bijia_rel_idx = np.argmax(mises[13912:34457])
        bijia_idx = bijia_rel_idx + 13912
        bijia_max = round(float(mises[bijia_idx]), 2)

        # 大拉杆最大值 & 索引
        dalagan_rel_idx = np.argmax(mises[34457:48377])
        dalagan_idx = dalagan_rel_idx + 34457
        dalagan_max = round(float(mises[dalagan_idx]), 2)

        # 记录到日志
        logger.info(f"应力组件计算人字架最大值: {renzijia_max} 最大索引为 {renzijia_idx}")
        logger.info(f"应力组件计算臂架最大值: {bijia_max} 最大索引为 {bijia_idx}")
        logger.info(f"应力组件计算大拉杆最大值: {dalagan_max} 最大索引为 {dalagan_idx}")

        # 固定索引的杆件
        rod1_val = round(float(mises[80894]), 2)
        rod2_val = round(float(mises[69312]), 2)
        rod3_val = round(float(mises[69599]), 2)
        rod4_val = round(float(mises[69053]), 2)
        rod5_val = round(float(mises[67387]), 2)

        # 获取当前北京时间
        tz = pytz.timezone('Asia/Shanghai')
        T_now = datetime.now(tz).strftime('%Y-%m-%d %H:%M:%S.%f')[:-3]  # 毫秒精度

        # 插入数据到数据库
        sql = """
        INSERT INTO t_data_sim_analysis
        (T, value4, value5, value6, value7, value8, value1, value2, value3, 
         struct_id, fea_cal_time)
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
        """

        cursor.execute(sql, (
            T_now,
            rod1_val,
            rod2_val,
            rod3_val,
            rod4_val,
            rod5_val,
            renzijia_max,
            dalagan_max,
            bijia_max,
            struct_id,
            calculate_time
        ))

        conn.commit()
        logger.info(f"成功插入struct_id={struct_id}的计算结果(时间:{calculate_time})")

    except Exception as e:
        logging.error(f"数据库查询错误: {e}", exc_info=True)
        return None

    finally:
        if 'conn' in locals() and conn.open:
            conn.close()



    # return {
    #     "renzijia": (renzijia_max, renzijia_idx),
    #     "bijia": (bijia_max, bijia_idx),
    #     "dalagan": (dalagan_max, dalagan_idx)
    # }

if __name__ == '__main__':
    config = {
        'host': '192.168.3.8',
        'port': 30053,
        'user': 'nbport_user',
        'password': 'P@jtMUS3Yx',
        'database': 'nbport_qzj_db',
        'charset': 'utf8mb4'
    }
    insert_sim_data(config,1001)