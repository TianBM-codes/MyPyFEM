# db_ops.py
import threading
import mysql.connector
import pandas as pd
from typing import List, Tuple

# ==================== 数据库配置 ====================
DB_CONFIG = {
    'host': 'monitor.sipesc.net',
    "port": 30053,
    'user': 'nbport_user',
    'password': 'P@jtMUS3Yx',
    'database': 'nbport_qzj_db',
    'autocommit': True,
    'connection_timeout': 300,
}
# qi_zhong_ji = {"host": "monitor.sipesc.net", "port": 30053, "database": "nbport_qzj_db",
#                "username": "nbport_user", "password": "P@jtMUS3Yx"}
# 线程本地存储连接
_thread_local = threading.local()

def get_connection():
    if not hasattr(_thread_local, "connection"):
        # print(f"[{threading.current_thread().name}] 创建数据库连接")
        _thread_local.connection = mysql.connector.connect(**DB_CONFIG)
    return _thread_local.connection

# def query_data(start_t: str, end_t: str, table: str = "t_data_work_status_info") -> pd.DataFrame:
def query_data(sql) -> pd.DataFrame:
    """
    查询指定时间范围的数据
    """
    try:
        conn = get_connection()
        cursor = conn.cursor()

        if(type(sql)==str):
            cursor.execute(sql)
        else:
            cursor.execute(sql[0], sql[1])
        result = cursor.fetchall()
        columns = [desc[0] for desc in cursor.description]
        cursor.close()

        return pd.DataFrame(result, columns=columns)

    except Exception as e:
        print(f"[{threading.current_thread().name}] ❌ 查询失败: {e}")
        return pd.DataFrame()

# def save_data(df: pd.DataFrame, table: str = "t_processed_results"):
def save_data(insert_sql, data_to_insert):
    """
    保存处理结果到数据库
    """
    try:
        conn = get_connection()
        cursor = conn.cursor()
        # 插入数据
        cursor.execute(insert_sql, data_to_insert)
        conn.commit()
        cursor.close()

        print(f"[{threading.current_thread().name}] ✅ 成功保存 {len(data_to_insert)} 行数据")

    except Exception as e:
        print(f"[{threading.current_thread().name}] ❌ 保存失败: {e}")
        conn.rollback()