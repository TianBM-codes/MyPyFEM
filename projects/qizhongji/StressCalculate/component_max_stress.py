import pymysql
from datetime import datetime
import pytz
import logging
import log_config  # 导入日志配置模块

logger = logging.getLogger()

def update_component_max_stress(config, struct_id):
    try:
        tz = pytz.timezone('Asia/Shanghai')
        T_now = datetime.now(tz).strftime('%Y-%m-%d %H:%M:%S.%f')[:-3]

        conn = pymysql.connect(**config)
        cursor = conn.cursor(pymysql.cursors.DictCursor)

        # 读取 t_data_stress_analysis 最新4条数据（指定 struct_id）
        cursor.execute("""
            SELECT * FROM t_data_stress_analysis
            WHERE struct_id = %s
            ORDER BY T DESC
            LIMIT 4
        """, (struct_id,))
        latest_stress_records = cursor.fetchall()

        # 从 t_data_sim_analysis 读取最新一条数据（指定 struct_id）
        cursor.execute("""
            SELECT value1, value2, value3, value4, value5, value6, value7, value8
            FROM t_data_sim_analysis
            WHERE struct_id = %s
            ORDER BY T DESC
            LIMIT 1
        """, (struct_id,))
        latest_sim = cursor.fetchone()

        if not latest_sim:
            print(f"struct_id={struct_id} 的 t_data_sim_analysis 表没有数据，无法执行")
            return

        # 计算4个部件的应力值
        new_data = {
            "象鼻梁": {
                "compressive_stress": max(latest_sim["value5"], latest_sim["value7"]),
                "tensile_stress": max(latest_sim["value4"], latest_sim["value6"], latest_sim["value8"])
            },
            "人字架": {
                "compressive_stress": 0,
                "tensile_stress": latest_sim["value1"]
            },
            "大拉杆": {
                "compressive_stress": 0,
                "tensile_stress": latest_sim["value2"]
            },
            "臂架": {
                "compressive_stress": 0,
                "tensile_stress": latest_sim["value3"]
            }
        }

        # 判断是否插入（表为空或最新记录全为0）
        if not latest_stress_records or all(
            rec["compressive_stress"] == 0 and rec["tensile_stress"] == 0
            for rec in latest_stress_records
        ):
            print(f"struct_id={struct_id} 执行插入操作")
            insert_sql = """
                INSERT INTO t_data_stress_analysis (T, mid, compressive_stress, tensile_stress, struct_id)
                VALUES (%s, %s, %s, %s, %s)
            """
            for mid, stresses in new_data.items():
                logger.info(
                    f"插入记录 | struct_id={struct_id} | 部件={mid} | T={T_now} | 压应力={-stresses['compressive_stress']} | 拉应力={stresses['tensile_stress']}")
                cursor.execute(insert_sql, (T_now, mid, -stresses["compressive_stress"], stresses["tensile_stress"], struct_id))

        else:
            print(f"struct_id={struct_id} 执行更新操作")
            update_sql = """
                UPDATE t_data_stress_analysis
                SET compressive_stress = %s, 
                    tensile_stress = %s, 
                    T = %s
                WHERE mid = %s AND T = %s AND struct_id = %s
            """
            # 最新一次运行的T值，用于定位更新
            latest_T = latest_stress_records[0]["T"]
            logger.info(f"struct_id={struct_id} 最新批次 T = {latest_T}")

            for rec in latest_stress_records:
                mid = rec["mid"]
                old_c = rec["compressive_stress"]
                old_t = rec["tensile_stress"]

                # ---- 压应力比较 ----
                abs_old_c = abs(old_c)
                abs_new_c = abs(new_data[mid]["compressive_stress"])
                if abs_new_c > abs_old_c:
                    new_c = -abs_new_c
                else:
                    new_c = old_c

                # ---- 拉应力比较 ----
                new_t = max(old_t, new_data[mid]["tensile_stress"])



                if new_c != old_c or new_t != old_t:
                    logger.info(
                        f"{mid} | 原压应力={old_c} → 新压应力候选={new_c} | 原拉应力={old_t} → 新拉应力候选={new_t}")
                    logger.info(f"更新记录 | struct_id={struct_id} | mid={mid} | T 从 {latest_T} 改为 {T_now}")
                    cursor.execute(update_sql, (new_c, new_t, T_now, mid, latest_T, struct_id))
                else:
                    logger.info(f"跳过更新 | struct_id={struct_id} | mid={mid} (无变化)")

        conn.commit()
        print(f"struct_id={struct_id} 最大应力统计已更新完成")

    except Exception as e:
        print(f"执行失败: {e}")
    finally:
        if 'conn' in locals() and conn.open:
            conn.close()
            print(f"struct_id={struct_id} 数据库连接已关闭")


if __name__ == '__main__':
    config = {
        'host': '192.168.3.8',
        'port': 30053,
        'user': 'nbport_user',
        'password': 'P@jtMUS3Yx',
        'database': 'nbport_qzj_db',
        'charset': 'utf8mb4'
    }
    update_component_max_stress(config,1001)
