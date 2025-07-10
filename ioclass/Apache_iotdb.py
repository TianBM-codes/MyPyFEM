from iotdb.Session import Session
from iotdb.utils.IoTDBConstants import TSDataType
import time


def read_sensor_data_raw(device_path, sensor_name, start_time, end_time):
    """
    读取传感器原始数据

    参数:
        device_path (str): 设备路径 (如 "root.factory.device1")
        sensor_name (str): 传感器名称 (如 "temperature")
        start_time (int/str): 开始时间戳或时间表达式
        end_time (int/str): 结束时间戳或时间表达式

    返回:
        list: 原始数据列表，格式为[(timestamp1, value1), (timestamp2, value2), ...]
    """
    # 构造SQL查询语句
    query = f"""
    SELECT {sensor_name} 
    FROM {device_path} 
    WHERE time >= {start_time} AND time <= {end_time}
    """

    # 执行查询
    result = session.execute_query_statement(query)
    data = []

    # 遍历查询结果
    while result.has_next():
        row = result.next()
        # 获取时间戳和传感器值
        data.append((row.get_timestamp(), row.get_fields()[0].get_float_value()))

    return data


def batch_insert_single_sensor(session, device_path, sensor_name, data_type, data_points):
    """
    单测点批量插入函数

    参数:
        session: IoTDB会话对象
        device_path (str): 设备路径 (如 "root.factory.device1")
        sensor_name (str): 传感器名称 (如 "temperature")
        data_type (TSDataType): 数据类型 (如 TSDataType.FLOAT)
        data_points (list): 数据点列表，格式为[(timestamp1, value1), (timestamp2, value2),...]
    """
    # 准备批量插入参数
    device_ids = [device_path] * len(data_points)  # 设备ID列表
    timestamps = [ts for (ts, _) in data_points]  # 时间戳列表
    measurements_lst = [[sensor_name]] * len(data_points)  # 测量名称二维列表
    types_lst = [[data_type]] * len(data_points)  # 数据类型二维列表
    values_lst = [[val] for (_, val) in data_points]  # 值列表（每个值单独包装成列表）

    # 执行批量插入
    session.insert_records(
        device_ids=device_ids,
        times=timestamps,
        measurements_lst=measurements_lst,
        types_lst=types_lst,
        values_lst=values_lst
    )


def batch_insert_multi_sensors(session, device_path, measurements, data_types, data_points):
    """
    多测点批量插入函数

    参数:
        session: IoTDB会话对象
        device_path (str): 设备路径
        measurements (list): 传感器名称列表 ["sensor1", "sensor2"]
        data_types (list): 数据类型列表 [TSDataType.FLOAT, TSDataType.INT32]
        data_points (list): 数据点列表，格式为[(timestamp, [value1, value2,...]), ...]
    """
    # 准备批量插入参数
    device_ids = [device_path] * len(data_points)  # 设备ID列表
    timestamps = [ts for (ts, _) in data_points]  # 时间戳列表
    measurements_lst = [measurements] * len(data_points)  # 测量名称二维列表
    types_lst = [data_types] * len(data_points)  # 数据类型二维列表
    values_lst = [values for (_, values) in data_points]  # 值列表

    # 执行批量插入
    session.insert_records(
        device_ids=device_ids,
        times=timestamps,
        measurements_lst=measurements_lst,
        types_lst=types_lst,
        values_lst=values_lst
    )


if __name__ == "__main__":
    # 1. 连接IoTDB数据库
    # 参数说明:
    # host: IoTDB服务器地址
    # port: IoTDB服务端口(默认6667)
    # user: 用户名
    # password: 密码
    session = Session(host="127.0.0.1", port=6667, user="root", password="root")

    # 打开会话(False表示不启用RPC压缩)
    session.open(False)

    # 检查连接是否成功
    if not session.is_open():
        raise ConnectionError("Failed to connect to IoTDB")

    try:
        # 2. 单测点数据插入示例
        # 生成测试数据(当前时间往前5秒、4秒、3秒的温度数据)
        single_sensor_data = [
            (int(time.time() * 1000) - 5000, 24.8),  # (timestamp, value)
            (int(time.time() * 1000) - 4000, 25.1),
            (int(time.time() * 1000) - 3000, 25.3)
        ]

        # 调用单测点批量插入函数
        batch_insert_single_sensor(
            session=session,
            device_path="root.factory.device2",
            sensor_name="temperature",
            data_type=TSDataType.FLOAT,
            data_points=single_sensor_data
        )

        # 3. 多测点数据插入示例
        # 生成测试数据(当前时间往前2秒、1秒、当前时刻的温度和湿度数据)
        multi_sensor_data = [
            (int(time.time() * 1000) - 2000, [25.1, 45]),  # (timestamp, [temp, humid])
            (int(time.time() * 1000) - 1000, [25.3, 46]),
            (int(time.time() * 1000), [25.5, 47])
        ]

        # 调用多测点批量插入函数
        batch_insert_multi_sensors(
            session=session,
            device_path="root.factory.device2",
            measurements=["temperature", "humidity"],
            data_types=[TSDataType.FLOAT, TSDataType.FLOAT],
            data_points=multi_sensor_data
        )

        # 4. 查询验证插入的数据
        print("温度数据:")
        temp_data = read_sensor_data_raw(
            device_path="root.factory.device2",
            sensor_name="temperature",
            start_time="now() - 1h",  # 时间表达式: 当前时间往前1小时
            end_time="now()"  # 时间表达式: 当前时间
        )
        print(temp_data)

        print("\n湿度数据:")
        humid_data = read_sensor_data_raw(
            device_path="root.factory.device2",
            sensor_name="humidity",
            start_time="now() - 1h",
            end_time="now()"
        )
        print(humid_data)

    finally:
        # 5. 关闭数据库连接
        session.close()