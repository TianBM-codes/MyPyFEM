# !/usr/bin/env python3
# -*- coding: utf-8 -*-

from element.Node import Node


class Geom(object):
    """
    存储几何信息
    """

    def __init__(self):
        self.node_count = -1
        self.node_list = []  # List of all nodes in the domain, 实例化数据
        self.node_hash = {}  # 节点真实Id对应nodelist中的index的Hash表

    def AddNode(self, node: Node, node_id: int):
        """
        向系统中添加节点，
        @param node: 文本中读取的节点
        @param node_id: 导入文件中的编号
        @return:
        """
        self.node_list.append(node)
        self.node_hash[node_id] = len(self.node_list)