import zlib
import struct
from io import BytesIO
import numpy as np
import os
import logging
import log_config  # 导入日志配置模块
import pathlib
import pymysql
import rainflow as rf  # 需要安装rainflow库: pip install rainflow
from datetime import datetime, timedelta
from data_peprocess import *
from read_file_name import *
import math

# 获取日志记录器
logger = logging.getLogger()

def read_file_name(config, period):
    """
    根据period参数，取参考时间前一周、前一月或前一年的文件名列表。week、month、year

    :param config: 数据库连接配置字典
    :param period: 字符串，'week'、'month'或'year'，默认'week'
    :return: 文件名列表
    """

    # 固定参考时间（可自行修改）
    reference_time_str = "2025-06-01 00:00:00"
    reference_time = datetime.strptime(reference_time_str, "%Y-%m-%d %H:%M:%S")
    # # 使用当前时间作为参考时间
    # reference_time = datetime.now()

    ref_date = reference_time.date()

    if period == 'week':
        # 取参考时间前一周的周一和周日
        weekday = ref_date.weekday()
        start_date = ref_date - timedelta(days=weekday + 7)  # 上一周周一
        end_date = start_date + timedelta(days=6)           # 上一周周日

        start_datetime = datetime.combine(start_date, datetime.min.time())
        end_datetime = datetime.combine(end_date, datetime.max.time())

    elif period == 'month':
        # 取参考时间前一个月的第一天和最后一天
        # 计算前一个月的年份和月份
        year = ref_date.year
        month = ref_date.month

        if month == 1:
            prev_month_year = year - 1
            prev_month = 12
        else:
            prev_month_year = year
            prev_month = month - 1

        # 上个月第一天
        start_datetime = datetime(prev_month_year, prev_month, 1)

        # 计算上个月最后一天，方法：本月第一天减1天
        current_month_first_day = datetime(year, month, 1)
        end_datetime = current_month_first_day - timedelta(seconds=1)
        end_datetime = end_datetime.replace(hour=23, minute=59, second=59, microsecond=999999)

    elif period == 'year':
        # 取参考时间前一年的第一天和最后一天
        prev_year = ref_date.year - 1

        start_datetime = datetime(prev_year, 1, 1)
        end_datetime = datetime(prev_year, 12, 31, 23, 59, 59, 999999)

    else:
        logger.error(f"无效的period参数: {period}，必须是'week'、'month'或'year'")
        return []
    logger.info(f"从t_calculate_nephogram数据库中读取文件为起始时间为{start_datetime}到{end_datetime}")
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
            cursor.execute(query, (start_datetime, end_datetime))
            results = cursor.fetchall()

        if not results:
            logger.info("查询结果为空：没有符合条件的文件。")
            return []

        return [row[0] for row in results]

    except Exception as e:
        logger.error("数据库查询出错", exc_info=True)
        return []

    finally:
        if 'connection' in locals() and connection.open:
            connection.close()
def ReaddatAsNumpy(mises_path):
    """
    解析mises结果
    """
    with open(mises_path, 'rb') as fileData:
        rawData = fileData.read()
        uncompressed_data = zlib.decompress(rawData)
        buffer = BytesIO(uncompressed_data)
        res_length = struct.unpack('i', buffer.read(4))[0]
        mises = struct.unpack(f'{res_length}f', buffer.read(res_length * 4))
    return np.array(mises)

def process_dat_files(file_paths):
    logger.info(f"找到 {len(file_paths)} 个应力的.dat文件:")

    arrays = []
    for file_path in file_paths:
        try:
            arr = ReaddatAsNumpy(file_path).astype(np.float32)
            arrays.append(arr)
        except Exception as e:
            logger.error(f"读取文件 {file_path} 时出错: {str(e)}")
            continue

    if not arrays:
        logger.error("没有成功读取任何文件")
        return

    lengths = [len(arr) for arr in arrays]
    if len(set(lengths)) > 1:
        logger.error(f"错误: 数组长度不一致 - {lengths}")
        return

    combined = np.stack(arrays, axis=0)
    print(f"\n合并数组形状: {combined.shape} (文件数×数据长度)")

    return combined


def fatigue_analysis(x, Su, material_params):
    """
    完整的疲劳分析流程：
    1. 雨流计数获取循环数据
    2. Goodman修正平均应力影响
    3. Miner线性损伤累积预测疲劳寿命

    参数:
        x: 载荷时间序列 (numpy数组)
        Su: 材料抗拉强度 (MPa)
        material_params: 字典，包含材料S-N曲线参数:
            {
                'S_endurance': 疲劳极限 (MPa),
                'a': Basquin方程参数a,
                'b': Basquin方程参数b
            }
        Basquin方程：N = a * S^b

    返回:
        result: 包含分析结果的字典
    """
    # 1. 去除大于抗拉强度的异常值并插值
    x_cleaned = remove_outliers_by_strength(x, Su)
    # 1. 利用雨流计数法获取应力幅值、平均应力、应力循环次数、循环起点索引、终点索引，并转换为numpy数组，为n×5数组
    c = np.array(list(rf.extract_cycles(x_cleaned)))
    c[:, 0] /= 2  # 在ASTM E1049标准中，extract_cycles() 返回的range是峰谷之间的总跨度（即 max - min），而应力幅值是范围的一半：

    # 2. Goodman修正平均应力影响
    # --对非对称循环（平均应力不等于0）下的应力范围/应力幅进行修正，消除平均应力影响。
    idx = np.where(c[..., 1] < 0)  # 如果 c 是一个二维数组（例如 c 的形状为 (5, 3)），则 c[..., 1] 选择的是第二列
    c[..., 1][idx] = 0  # 将平均压应力统一设为 0，因为压应力一般对疲劳是有益的，但为了简化处理，统一设为0，相当于不再使用goodman模型修正了

    # 检查是否有平均应力大于Su的情况
    idx_invalid = np.where(c[..., 1] >= Su)
    if len(idx_invalid[0]) > 0:
        # 存在平均应力大于Su的情况，直接返回失效
        result = {
            'cycles_data': c,  # 修正后的循环数据
            'total_damage': 1.0,  # 总损伤设为1表示立即失效
            'fatigue_life': 0,  # 预测剩余循环次数为0
            'num_cycles': len(c),  # 识别的循环数量
            'material_params': material_params
        }
        return result

        # 正常Goodman修正
    c[..., 0] = c[..., 0] / (1 - c[..., 1] / Su)  # 采用Goodman修正模型算法对应力幅值进行修正

    # 3. 损伤计算
    total_damage = 0.0
    S_endurance = material_params['S_endurance']
    a = material_params['a']
    b = material_params['b']

    for cycle in c:
        S_amp, n_i = cycle[0], cycle[2]
        if S_amp > S_endurance:  # 只计算高于疲劳极限的损伤
            N_i = a * (S_amp ** b)
            total_damage += n_i / N_i

            # 4. 剩余寿命计算
    # if total_damage == 0:
    #     remaining_cycles = 5e6  # 无限寿命
    # elif total_damage >= 1:
    #     remaining_cycles = 0  # 已失效
    # else:
    #     # 计算还能承受多少次当前载荷谱的重复
    #     remaining_cycles =  (1 / total_damage - 1) * sum(c[:, 2])
    #     # remaining_cycles = 1  / total_damage * 20
    # 整理结果
    result = {
        'cycles_data': c,
        'total_damage': min(total_damage, 1.0),  # 损伤最大为1
        # 'fatigue_life': remaining_cycles,
        'num_cycles': len(c),
        'total_cycle': sum(c[:, 2]) # 添加总循环次数
    }
    return result


def analyze_dat_fatigue(file_paths, Su, material_params):
    """
    对 process_dat_files 的结果进行疲劳分析，找出总损伤最大的情况

    参数:
        file_paths: dat文件路径列表
        Su: 材料抗拉强度 (MPa)
        material_params: 材料参数字典

    返回:
        total_damage_array: 各列的总损伤数组
        total_cycle_array: 各列的总循环数组
        worst_case_series: 损伤最大的那一列原始应力时间序列（1D numpy 数组）
    """
    results_array = process_dat_files(file_paths)
    if results_array is None or results_array.shape[1] == 0:
        print("没有可分析的结果数据")
        return None, None, None

    n_groups = results_array.shape[1]
    fatigue_life_array = np.zeros(n_groups)
    total_damage_array = np.zeros(n_groups)
    total_cycle_array = np.zeros(n_groups)

    for i in range(n_groups):
        x = results_array[:, i]
        try:
            analysis_result = fatigue_analysis(x, Su, material_params)

            # fatigue_life_array[i] = analysis_result['fatigue_life']
            total_damage_array[i] = analysis_result['total_damage']
            total_cycle_array[i] = analysis_result['total_cycle']
            print(f"已经完成第 {i} 组数据的疲劳分析")
        except Exception as e:
            print(f"分析第 {i} 组数据失败: {str(e)}")
            fatigue_life_array[i] = np.nan
            total_damage_array[i] = np.nan
            total_cycle_array[i] = np.nan

    print(f"完成疲劳分析，fatigue_life 数组 shape: {fatigue_life_array.shape}")

    # === 找出最大总损伤对应的列 ===
    if np.all(np.isnan(total_damage_array)):
        print("所有组的损伤分析均失败。")
        return total_damage_array, total_cycle_array, None

    max_index = np.nanargmax(total_damage_array)
    worst_case_series = results_array[:, max_index]  # 返回该列的应力序列

    return total_damage_array, total_cycle_array, worst_case_series

def compute_damage_distribution(cycles_data, material_params):
    """
    统计不同应力幅值区间的循环次数和损伤贡献率

    参数:
        cycles_data: numpy 数组，包含 [应力幅值, 平均应力, 循环次数...]
        material_params: 包含疲劳极限、Basquin 参数a/b

    返回:
        result_dict: 每一档区间对应的 循环次数、损伤值、贡献率
    """
    bins = [0, 40, 60, 80, 100, 120, 140, float('inf')]
    bin_labels = ['≤40', '(40,60]', '(60,80]', '(80,100]', '(100,120]', '(120,140]', '>140']
    bin_counts = [0] * (len(bins) - 1)
    bin_damage = [0.0] * (len(bins) - 1)

    a = material_params['a']
    b = material_params['b']
    S_endurance = material_params['S_endurance']

    for cycle in cycles_data:
        S_amp = cycle[0] # 应力幅值
        n_i = cycle[2]  # 循环次数
        # 通过这个逻辑，找出 S_amp 落在哪个区间 bins[i] ~ bins[i+1]。将这个循环次数 n_i 加入该区间的循环总次数统计里。
        for i in range(len(bins) - 1):
            if bins[i] < S_amp <= bins[i+1]:
                bin_counts[i] += n_i
                if S_amp > S_endurance:
                    N_i = a * (S_amp ** b)
                    bin_damage[i] += n_i / N_i
                break

    total_damage = sum(bin_damage)
    damage_ratio = [d / total_damage if total_damage > 0 else 0 for d in bin_damage]

    result_dict = {}
    for i, label in enumerate(bin_labels):
        result_dict[label] = {
            'cycles': bin_counts[i],
            'damage': bin_damage[i],
            'damage_ratio': damage_ratio[i]
        }

    return result_dict

def insert_damage_distribution_to_db(result_dict, config, period='week'):
    """
    将 compute_damage_distribution 得到的结果插入数据库

    参数：
        result_dict: 来自 compute_damage_distribution 的结果字典
        config: 数据库连接配置
        period: 'week', 'month', 'year'，决定 reference_time 时间的推算方式
    """

    # 固定参考时间（可换为动态传入）
    reference_time_str = "2025-06-01 00:00:00"
    base_time = datetime.strptime(reference_time_str, "%Y-%m-%d %H:%M:%S")
    ref_date = base_time.date()

    # 计算 reference_time
    if period == 'week':
        weekday = ref_date.weekday()
        last_sunday = ref_date - timedelta(days=weekday + 1)
        reference_time = datetime.combine(last_sunday, datetime.max.time()).replace(microsecond=0)
    elif period == 'month':
        year, month = ref_date.year, ref_date.month
        first_day_this_month = datetime(year, month, 1)
        last_day_prev_month = first_day_this_month - timedelta(days=1)
        reference_time = datetime.combine(last_day_prev_month.date(), datetime.max.time()).replace(microsecond=0)
    elif period == 'year':
        reference_time = datetime(ref_date.year - 1, 12, 31, 23, 59, 59, 0)
    else:
        raise ValueError(f"无效的 period 参数: {period}")

    try:
        conn = pymysql.connect(**config)
        cursor = conn.cursor()

        for idx, (stress_range, stats) in enumerate(result_dict.items(), start=1):
            cycle_count_ceil = math.ceil(stats['cycles'])  # 向上取整

            # 跳过 cycle_count 为 0 的数据
            if cycle_count_ceil == 0:
                continue

            # damage_ratio * 100，并保留 2 位小数
            damage_ratio_pct = round(stats['damage_ratio'] * 100, 2)

            insert_sql = """
                INSERT INTO t_data_cyclic_frequency_distribution
                (T, frequency, show_order, stress_range, cycle_count, damage_ratio)
                VALUES (%s, %s, %s, %s, %s, %s)
                ON DUPLICATE KEY UPDATE
                    show_order = VALUES(show_order),
                    cycle_count = VALUES(cycle_count),
                    damage_ratio = VALUES(damage_ratio)
            """

            cursor.execute(insert_sql, (
                reference_time,
                period,
                idx,
                stress_range,
                cycle_count_ceil,
                damage_ratio_pct
            ))

        conn.commit()
        print("循环次数与损伤变化率已成功写入数据库。")

    except Exception as e:
        print("写入数据库失败:", e)

    finally:
        if cursor:
            cursor.close()
        if conn:
            conn.close()


# 示例使用
if __name__ == "__main__":
    # 读取文件路径
    config = {
        'host': '192.168.3.8',
        'port': 30053,
        'user': 'nbport_user',
        'password': 'P@jtMUS3Yx',
        'database': 'nbport_qzj_db',
        'charset': 'utf8mb4'
    }
    period = 'week'
    file_paths = read_file_name(config, period)

    # 2. 配置材料参数
    Su = 460  # 材料抗拉强度 (MPa)
    material_params = {
        'S_endurance': 100,  # 疲劳极限 (MPa)
        'a': 1e12,  # Basquin方程参数a
        'b': -3  # Basquin方程参数b
    }
    results_array = process_dat_files(file_paths)

