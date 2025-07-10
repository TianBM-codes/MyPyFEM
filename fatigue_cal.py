import numpy as np
import rainflow as rf

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
    # 1. 利用雨流计数法获取应力幅值、平均应力、应力循环次数、循环起点索引、终点索引，并转换为numpy数组，为n×5数组
    c = np.array(list(rf.extract_cycles(x)))
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
    if total_damage == 0:
        remaining_cycles = 5e6  # 无限寿命
    elif total_damage >= 1:
        remaining_cycles = 0  # 已失效
    else:
        # 计算还能承受多少次当前载荷谱的重复
        remaining_cycles =  (1 / total_damage - 1) * sum(c[:, 2])
        # remaining_cycles = 1  / total_damage * 20
    # 整理结果
    result = {
        'cycles_data': c,
        'total_damage': min(total_damage, 1.0),  # 损伤最大为1
        'fatigue_life': remaining_cycles,
        'num_cycles': len(c),
        'total_cycle': sum(c[:, 2]) # 添加总循环次数
    }
    return result

if __name__ == "__main__":
    # 1. 首先运行原始处理函数
    folder_path = r"./data/SAETransmission.dat"
    data = np.loadtxt(folder_path)

    # 2. 配置材料参数
    Su = 3800  # 材料抗拉强度 (MPa)
    material_params = {
        'S_endurance': 170,  # 疲劳极限 (MPa)
        'a': 1e12,  # Basquin方程参数a
        'b': -3  # Basquin方程参数b
    }
    results = fatigue_analysis(data, Su, material_params)
    print(results)