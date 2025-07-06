from iotdb.Session import Session
from iotdb.utils.IoTDBConstants import TSDataType
import time


def read_sensor_data_raw(device_path, sensor_name, start_time, end_time):
    """
    返回原始数据列表，格式: [(timestamp1, value1), (timestamp2, value2), ...]
    """
    query = f"""
    SELECT {sensor_name} 
    FROM {device_path} 
    WHERE time >= {start_time} AND time <= {end_time}
    """

    result = session.execute_query_statement(query)
    data = []
    while result.has_next():
        row = result.next()
        data.append((row.get_timestamp(), row.get_fields()[0].get_float_value()))

    return data


def batch_insert_single_sensor(session, device_path, sensor_name, data_type, data_points):
    """
    单测点批量插入函数
    参数:
        session: IoTDB会话对象
        device_path: 设备路径 (如 "root.factory.device1")
        sensor_name: 单个传感器名称 (如 "temperature")
        data_type: 单个数据类型 (如 TSDataType.FLOAT)
        data_points: [(timestamp1, value1), (timestamp2, value2),...]
    """
    # 准备参数（注意所有列表都是二维的）
    device_ids = [device_path] * len(data_points)
    timestamps = [ts for (ts, _) in data_points]
    measurements_lst = [[sensor_name]] * len(data_points)  # 二维列表
    types_lst = [[data_type]] * len(data_points)  # 二维列表
    values_lst = [[val] for (_, val) in data_points]  # 关键：将每个值包装成列表

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
        device_path: 设备路径
        measurements: 传感器名称列表 ["sensor1", "sensor2"]
        data_types: 数据类型列表 [TSDataType.FLOAT, TSDataType.INT32]
        data_points: [(timestamp, [value1, value2,...]), ...]
    """
    # 准备参数
    device_ids = [device_path] * len(data_points)
    timestamps = [ts for (ts, _) in data_points]
    measurements_lst = [measurements] * len(data_points)
    types_lst = [data_types] * len(data_points)
    values_lst = [values for (_, values) in data_points]

    # 执行批量插入
    session.insert_records(
        device_ids=device_ids,
        times=timestamps,
        measurements_lst=measurements_lst,
        types_lst=types_lst,
        values_lst=values_lst
    )


if __name__ == "__main__":
    # 连接参数
    session = Session(host="127.0.0.1", port=6667, user="root", password="root")
    session.open(False)

    # 检查连接
    if not session.is_open():
        raise ConnectionError("Failed to connect to IoTDB")

    try:
        # 示例1：单测点插入
        single_sensor_data = [
            (int(time.time() * 1000) - 5000, 24.8),
            (int(time.time() * 1000) - 4000, 25.1),
            (int(time.time() * 1000) - 3000, 25.3)
        ]
        batch_insert_single_sensor(
            session=session,
            device_path="root.factory.device2",
            sensor_name="temperature",
            data_type=TSDataType.FLOAT,
            data_points=single_sensor_data
        )

        # 示例2：多测点插入
        multi_sensor_data = [
            (int(time.time() * 1000) - 2000, [25.1, 45]),
            (int(time.time() * 1000) - 1000, [25.3, 46]),
            (int(time.time() * 1000), [25.5, 47])
        ]
        batch_insert_multi_sensors(
            session=session,
            device_path="root.factory.device2",
            measurements=["temperature", "humidity"],
            data_types=[TSDataType.FLOAT, TSDataType.FLOAT],
            data_points=multi_sensor_data
        )

        # 查询验证
        print("温度数据:")
        temp_data = read_sensor_data_raw(
            device_path="root.factory.device2",
            sensor_name="temperature",
            start_time="now() - 1h",
            end_time="now()"
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
        session.close()
