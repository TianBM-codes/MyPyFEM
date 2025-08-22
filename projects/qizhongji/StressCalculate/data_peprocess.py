import numpy as np
import pandas as pd


def remove_outliers_and_interpolate(data, sigma_threshold=4):
    """
    使用3σ原则剔除离群点，并用线性插值填补缺失值。

    参数:
        data (np.ndarray或pd.Series): 输入的一维时间序列数据。
        sigma_threshold (float): 西格玛倍数阈值（默认为3）。

    返回:
        np.ndarray: 处理后的数据（离群点替换为插值结果）。
    """
    # 转换为pandas Series以便处理缺失值
    if isinstance(data, np.ndarray):
        data = pd.Series(data)

    # 1. 计算均值和标准差
    mean = data.mean()
    std = data.std()

    # 2. 标记离群点（超出μ ± 3σ范围）
    lower_bound = mean - sigma_threshold * std
    upper_bound = mean + sigma_threshold * std
    outliers_mask = (data < lower_bound) | (data > upper_bound)

    # 3. 将离群点设为NaN
    data_cleaned = data.copy()
    data_cleaned[outliers_mask] = np.nan

    # 4. 线性插值填补NaN
    data_interpolated = data_cleaned.interpolate(method='linear')

    # 边缘处理：如果两端有NaN，用最近有效值填充
    data_interpolated = data_interpolated.ffill().bfill()

    return data_interpolated.to_numpy()  # 返回NumPy数组

def remove_outliers_by_strength(data, Su):
    """
    剔除大于抗拉强度Su的值，并用线性插值填补缺失值。

    参数:
        data (np.ndarray或pd.Series): 输入的一维时间序列数据。
        Su (float): 抗拉强度阈值，大于此值的点将被剔除。

    返回:
        np.ndarray: 处理后的数据（异常点替换为插值结果）。
    """
    # 转换为pandas Series以便处理缺失值
    if isinstance(data, np.ndarray):
        data = pd.Series(data)

    # 1. 标记大于抗拉强度的点
    outliers_mask = (data > Su)

    # 2. 将异常点设为NaN
    data_cleaned = data.copy()
    data_cleaned[outliers_mask] = np.nan

    # 3. 线性插值填补NaN
    data_interpolated = data_cleaned.interpolate(method='linear')

    # 边缘处理：如果两端有NaN，用最近有效值填充
    data_interpolated = data_interpolated.ffill().bfill()

    return data_interpolated.to_numpy()  # 返回NumPy数组



def filter_small_amplitudes(data, fatigue_limit):
    """
    滤除小幅值载荷（幅值 < 疲劳极限的40%），保留显著波动。

    参数:
        data (np.ndarray): 载荷时间序列（一维数组）。
        fatigue_limit (float): 材料的疲劳极限（单位与data一致）。

    返回:
        np.ndarray: 过滤后的载荷序列。
    """
    if len(data) < 2:
        return data  # 无需处理单点数据

    threshold = 0.4 * fatigue_limit
    filtered_data = [data[0]]  # 保留第一个点

    for i in range(1, len(data)):
        # 计算当前点与前一个点的差值（幅值）
        delta = abs(data[i] - filtered_data[-1])

        if delta >= threshold:
            filtered_data.append(data[i])  # 保留显著波动点
        # 否则跳过小幅波动（不添加到结果中）

    return np.array(filtered_data)
# data = np.array([100, 102, 105, 80, 78, 110, 108, 115, 30, 35, 120])
# fatigue_limit = 20  # 假设疲劳极限为50 MPa
#
# # 滤除小幅波动
# filtered_data = filter_small_amplitudes(data, fatigue_limit)
#
# print("原始数据:", data)
# print("过滤后数据:", filtered_data)
def zero_mean_adjustment(data):
    """
    将输入数据的均值归零。

    参数:
        data (np.ndarray或list): 输入的一维时间序列数据。

    返回:
        np.ndarray: 均值归零后的数据。
    """
    data = np.asarray(data)  # 确保输入为NumPy数组
    if len(data) == 0:
        return data  # 空数组直接返回

    mean = np.mean(data)  # 计算均值
    zero_mean_data = data - mean  # 减去均值

    return zero_mean_data