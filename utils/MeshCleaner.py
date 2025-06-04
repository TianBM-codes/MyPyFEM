def MarkSurface(eles):
    """
    标记表面单元面片, 标记对应的面片, 修改的是引用, 原数据会被更改
    :return:
    """
    surface_tri_tree = {}
    surface_qua_tree = {}
    for ele in eles:
        for face in ele.getTriFaces():
            if not surface_tri_tree.__contains__(face.unique_key):
                surface_tri_tree[face.unique_key] = face
            else:
                surface_tri_tree[face.unique_key].is_surface = False
                surface_tri_tree.pop(face.unique_key)
                face.is_surface = False

        for face in ele.getQuaFaces():
            if not surface_qua_tree.__contains__(face.unique_key):
                surface_qua_tree[face.unique_key] = face
            else:
                surface_qua_tree[face.unique_key].is_surface = False
                face.is_surface = False
                surface_qua_tree.pop(face.unique_key)


def CalculateUniqueEdge(edges):
    """
    将重复的网格线删除
    @param edges:
    @return:
    """
    assert len(edges) % 2 == 0
    edge_count = int(len(edges) / 2)
    unique_edge_set = set()
    for ii in range(edge_count):
        sorted_nums = sorted([edges[ii * 2], edges[ii * 2 + 1]])
        iter_unique_key = f"{sorted_nums[0]},{sorted_nums[1]}"
        if iter_unique_key not in unique_edge_set:
            unique_edge_set.add(iter_unique_key)

    edge_lst = list(unique_edge_set)
    all_edge = []
    for edge in edge_lst:
        node_a, node_b = edge.split(",")
        all_edge.extend([int(node_a), int(node_b)])

    return all_edge


if __name__ == "__main__":
    test_edges = [1, 2, 3, 4, 2, 1]
    print(CalculateUniqueEdge(test_edges))
