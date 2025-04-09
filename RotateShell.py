import numpy as np
from main import *


def global_to_local_coords(global_point, element_nodes):
    """
    将全局坐标系中的点转换到壳单元的局部坐标系
    参数:
        global_point: 全局坐标系中的点 (3D向量)
        element_nodes: 壳单元的节点坐标列表 (至少3个节点，假设为四边形单元)
    返回:
        local_coords: 局部坐标系中的坐标 (3D向量)
    """
    # 将节点转换为numpy数组
    nodes = np.array(element_nodes)

    # 定义局部坐标系原点（假设为第一个节点）
    origin = nodes[0]

    # 计算局部坐标系的三个基向量
    # 第一个切向量 (沿1-2边方向)
    vec1 = nodes[1] - origin
    vec1 = vec1 / np.linalg.norm(vec1)  # 归一化

    # 第二个切向量 (沿1-4边方向，假设四边形单元)
    vec2_initial = nodes[3] - origin
    vec2_initial = vec2_initial / np.linalg.norm(vec2_initial)

    # 法向量 (通过叉乘得到)
    normal = np.cross(vec1, vec2_initial)
    normal = normal / np.linalg.norm(normal)

    # 修正第二个切向量使其正交
    vec2 = np.cross(normal, vec1)
    vec2 = vec2 / np.linalg.norm(vec2)

    # 构建变换矩阵
    transformation_matrix = np.vstack([vec1, vec2, normal]).T

    # 计算局部坐标
    relative_global = global_point - origin
    local_coords = np.dot(transformation_matrix, relative_global)

    return local_coords


# 示例使用
if __name__ == "__main__":
    # 定义壳单元节点（四边形单元，4个节点）
    element_nodes = [
        [0, 0, 0],  # 节点1
        [1, 0, 0],  # 节点2
        [1, 1, 0],  # 节点3
        [0, 1, 0]  # 节点4
    ]

    # 测试点（全局坐标系）
    test_point = np.array([0.5, 0.5, 0])

    # 转换到局部坐标系
    local = global_to_local_coords(test_point, element_nodes)

    input_file = r"D:\WorkSpace\FEM\MyPyFEM\numerical example\ANSYS\TwoQuaShellOneRotate.cdb"
    MyPyFEM(pathlib.Path(input_file), check_model=False, plot_stiff=False)
