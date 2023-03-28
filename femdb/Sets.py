#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np


class NodeSet(object):
    def __init__(self, n_name, nodes):
        self.name = n_name
        self.nodes = np.asarray(nodes, dtype=np.int32)

    def __str__(self):
        return "name: {}, from {} to {}, count is: {}".format(
            self.name, self.nodes[0], self.nodes[-1], len(self.nodes))

    def GetName(self):
        return self.name

    def GetNodeIds(self):
        return self.nodes

    def SetBoundary(self, begin_idx, end_idx, value=0, b_type=None):
        pass


class EleSet(object):
    def __init__(self, e_name, eles_ids):
        """
        如果used是False, 那么这个集合是没有被用到过的. 如果某些单元集合被赋予过属性, 那么used为True
        """
        self.name = e_name
        self.eles_ids = np.asarray(eles_ids, dtype=np.int32)

    def __str__(self):
        return "name: {}, from {} to {}, count is: {}".format(
            self.name, self.eles_ids[0], self.eles_ids[-1], len(self.eles_ids))

    def GetName(self):
        return self.name

    def GetEleIds(self):
        return self.eles_ids
