#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from utils.CustomException import *
from typing import List
import numpy as np
import abc


def QuaFace2Triangles(nodes, id_order):
    triangles = [[nodes[id_order[0]], nodes[id_order[1]], nodes[id_order[3]]],
                 [nodes[id_order[1]], nodes[id_order[2]], nodes[id_order[4]]],
                 [nodes[id_order[1]], nodes[id_order[6]], nodes[id_order[3]]],
                 [nodes[id_order[1]], nodes[id_order[4]], nodes[id_order[6]]],
                 [nodes[id_order[3]], nodes[id_order[6]], nodes[id_order[5]]],
                 [nodes[id_order[4]], nodes[id_order[7]], nodes[id_order[6]]]]
    return triangles


class MeshElementFactory:
    """
    Reference:
    1. https://abaqus-docs.mit.edu/2017/English/SIMACAEELMRefMap/simaelm-c-shellelem.htm
    """

    @staticmethod
    def CreateElement(e_type, e_id=-1, opt=None, use_low_order=False,
                      except_ele_type=None,
                      except_ele_id=None):
        """
        静态函数, 用于返回
        @param e_type: 单元类型，这里包含了Abaqus、Nastran和Ansys的
        @param e_id: 初始化单元需要单元ID
        @param opt: 附加参数, 比如181可能是3节点壳也可能是4节点壳, solid45可能是8节点也可能是4节点
        @param use_low_order: 是否直接使用低阶单元，而不用高阶单元
        @param except_ele_type:
        @param except_ele_id:
        :return: 单元和节点个数
        """
        if except_ele_id is None:
            except_ele_id = []
        if except_ele_type is None:
            except_ele_type = []

        if e_type in except_ele_type:
            return None, None
        if e_id in except_ele_id:
            return None, None

        if e_type in [185]:
            if opt == 8:
                return MeshC3D8(e_id), True
            elif opt == 6:
                return MeshC3D6(e_id), True
            elif opt == 4:
                return MeshTetra(e_id), True
            elif opt == 5:
                return MeshC3D5(e_id), True
            else:
                raise NoImplSuchElement(e_type, e_id)
        elif e_type in ["C3D4"]:
            return MeshTetra(e_id), True

        elif e_type in [40600, 40800]:
            return MeshTetra(e_id), True
        elif e_type in [40500, 'S4', "quad", "CookQuaShell"]:
            return MeshCQUAD4(e_id), True
        elif e_type in [20100, 188, 10, 4, 39, 14, 180]:
            return MeshTruss(e_id), False
        elif e_type in [30500, "triangle", "CookTriShell"]:
            return MeshTRIA3(e_id), True
        elif e_type in [10000]:
            return None, False
        elif e_type in [100600, 187]:
            if use_low_order:
                return MeshTetra(e_id), True
            else:
                return MeshC3D10(e_id), True
        elif e_type in [150600]:
            if use_low_order:
                return MeshC3D6(e_id), True
            else:
                return MeshC3D15(e_id), True
        elif e_type in [200600, 186, 80600]:
            if opt == 15:
                if use_low_order:
                    return MeshC3D6(e_id), True
                else:
                    return MeshC3D15(e_id), True
            elif opt == 13:
                if use_low_order:
                    return MeshC3D5(e_id), True
                else:
                    return MeshC3D13(e_id), True
            elif opt == 20:
                if use_low_order:
                    return MeshC3D8(e_id), True
                else:
                    return MeshC3D20(e_id), True
            elif opt is None:
                return MeshC3D8(e_id), True
            else:
                raise NoImplSuchElement(e_type, opt)
        elif e_type in [130600]:
            if use_low_order:
                return MeshC3D5(e_id), True
            else:
                return MeshC3D13(e_id), True
        elif e_type in [174, 170, 21, "Mass"]:
            return None, None
        else:
            raise NoImplSuchElement(e_type, e_id)


class MeshElement(metaclass=abc.ABCMeta):
    def __init__(self):
        self.id = -1
        self.nodesCount = 0
        self.degradeEle = None  # Assuming an integer or similar
        self.upgradeNodeIndex = None  #
        self.degradeRvIndex = []
        self.pts = []
        self.triangles = []
        self.nodeStr = ""
        self.triFaces: List[TriFace] = []
        self.quaFaces: List[QuaFace] = []
        self.nodeIds = []
        self.type_code = None
        self.edges = []

    def setId(self, _id):
        self.id = _id

    @abc.abstractmethod
    def setFaces(self, _nodeIds):
        pass  # This method is meant to be overridden by derived classes

    def getQuaFaces(self):
        return self.quaFaces

    def getTriFaces(self):
        return self.triFaces

    @abc.abstractmethod
    def getAllTriangles(self, only_surface=False):
        pass

    def getEdges(self):
        # if len(self.edges) == 0:
        #     print(f"Inner Element: {self.id}")
        return self.edges

    def __eq__(self, other):
        return self.id == other.id

    def __lt__(self, other):
        return self.id < other.id


class TriFace:
    def __init__(self, node_list):
        self.nodes = node_list
        self.unique_key = ",".join([str(ii) for ii in np.sort(node_list)])
        self.is_surface = True
        self.nodeIds = []

    def calculateNodeIds(self, map_):
        self.nodeIds = [map_[ii] for ii in self.nodes]


class QuaFace:
    def __init__(self, node_list):
        self.nodes = node_list  # 文件中的原编号
        self.unique_key = ",".join([str(ii) for ii in np.sort(node_list)])
        self.is_surface = True
        self.nodeIds = []  # 文件中经过hash后的编号

    def calculateNodeIds(self, map_):
        self.nodeIds = [map_[ii] for ii in self.nodes]


class MeshC3D6(MeshElement):
    """
    三棱柱单元
    """

    def __init__(self, _id):
        super().__init__()
        self.setId(_id)
        self.nodesCount = 6
        self.upgradeNodeIndex = [0, 1, 2, 3, 4, 5]

    def setFaces(self, nodeIds):
        self.nodeIds = nodeIds
        f1 = QuaFace([nodeIds[3], nodeIds[5], nodeIds[2], nodeIds[0]])
        f2 = QuaFace([nodeIds[1], nodeIds[2], nodeIds[5], nodeIds[4]])
        f3 = QuaFace([nodeIds[0], nodeIds[1], nodeIds[4], nodeIds[3]])
        f4 = TriFace([nodeIds[0], nodeIds[2], nodeIds[1]])
        f5 = TriFace([nodeIds[3], nodeIds[4], nodeIds[5]])
        self.quaFaces = [f1, f2, f3]
        self.triFaces = [f4, f5]

    def getAllTriangles(self, only_surface=False):
        self.triangles.clear()
        if (only_surface and self.quaFaces[0].is_surface) or not only_surface:
            triangle1 = [self.quaFaces[0].nodes[0], self.quaFaces[0].nodes[1], self.quaFaces[0].nodes[2]]
            self.triangles.append(triangle1)
            triangle2 = [self.quaFaces[0].nodes[0], self.quaFaces[0].nodes[2], self.quaFaces[0].nodes[3]]
            self.triangles.append(triangle2)
            self.edges.extend([self.quaFaces[0].nodes[0], self.quaFaces[0].nodes[1]])
            self.edges.extend([self.quaFaces[0].nodes[1], self.quaFaces[0].nodes[2]])
            self.edges.extend([self.quaFaces[0].nodes[2], self.quaFaces[0].nodes[3]])
            self.edges.extend([self.quaFaces[0].nodes[3], self.quaFaces[0].nodes[0]])

        if (only_surface and self.quaFaces[1].is_surface) or not only_surface:
            triangle1 = [self.quaFaces[1].nodes[0], self.quaFaces[1].nodes[1], self.quaFaces[1].nodes[2]]
            self.triangles.append(triangle1)
            triangle2 = [self.quaFaces[1].nodes[0], self.quaFaces[1].nodes[2], self.quaFaces[1].nodes[3]]
            self.triangles.append(triangle2)
            self.edges.extend([self.quaFaces[1].nodes[0], self.quaFaces[1].nodes[1]])
            self.edges.extend([self.quaFaces[1].nodes[1], self.quaFaces[1].nodes[2]])
            self.edges.extend([self.quaFaces[1].nodes[2], self.quaFaces[1].nodes[3]])
            self.edges.extend([self.quaFaces[1].nodes[3], self.quaFaces[1].nodes[0]])

        if (only_surface and self.quaFaces[2].is_surface) or not only_surface:
            triangle1 = [self.quaFaces[2].nodes[0], self.quaFaces[2].nodes[1], self.quaFaces[2].nodes[2]]
            self.triangles.append(triangle1)
            triangle2 = [self.quaFaces[2].nodes[0], self.quaFaces[2].nodes[2], self.quaFaces[2].nodes[3]]
            self.triangles.append(triangle2)
            self.edges.extend([self.quaFaces[2].nodes[0], self.quaFaces[2].nodes[1]])
            self.edges.extend([self.quaFaces[2].nodes[1], self.quaFaces[2].nodes[2]])
            self.edges.extend([self.quaFaces[2].nodes[2], self.quaFaces[2].nodes[3]])
            self.edges.extend([self.quaFaces[2].nodes[3], self.quaFaces[2].nodes[0]])

        if (only_surface and self.triFaces[0].is_surface) or not only_surface:
            self.triangles.append(self.triFaces[0].nodes)
            self.edges.extend([self.triFaces[0].nodes[0], self.triFaces[0].nodes[1]])
            self.edges.extend([self.triFaces[0].nodes[1], self.triFaces[0].nodes[2]])
            self.edges.extend([self.triFaces[0].nodes[2], self.triFaces[0].nodes[0]])

        if (only_surface and self.triFaces[1].is_surface) or not only_surface:
            self.triangles.append(self.triFaces[1].nodes)
            self.edges.extend([self.triFaces[1].nodes[0], self.triFaces[1].nodes[1]])
            self.edges.extend([self.triFaces[1].nodes[1], self.triFaces[1].nodes[2]])
            self.edges.extend([self.triFaces[1].nodes[2], self.triFaces[1].nodes[0]])

        return self.triangles


class MeshTruss(MeshElement):
    """
    杆单元
    """

    def __init__(self, _id):
        super().__init__()
        self.id = _id
        self.nodesCount = 2
        self.triangles = []
        self.upgradeNodeIndex = [0, 1]

    def setFaces(self, node_ids):
        self.nodeIds = node_ids
        self.triangles = [node_ids[0], node_ids[1]]

        self.edges = [node_ids[0], node_ids[1]]

    def getAllTriangles(self, only_surface=False):
        return self.triangles


class MeshC3D13(MeshElement):
    """
    13节点高阶金字塔单元
    """

    def __init__(self, _id):
        super().__init__()
        self.id = _id
        self.nodesCount = 13
        self.triangles = []

    def setFaces(self, _nodeIds):
        self.nodeIds = _nodeIds
        self.quaFaces.append(QuaFace([_nodeIds[3], _nodeIds[0], _nodeIds[1], _nodeIds[2]]))
        self.triFaces.append(TriFace([_nodeIds[4], _nodeIds[3], _nodeIds[0]]))
        self.triFaces.append(TriFace([_nodeIds[4], _nodeIds[0], _nodeIds[1]]))
        self.triFaces.append(TriFace([_nodeIds[4], _nodeIds[1], _nodeIds[2]]))
        self.triFaces.append(TriFace([_nodeIds[4], _nodeIds[2], _nodeIds[3]]))

    def getAllTriangles(self, only_surface=False):
        self.triangles.clear()
        if (only_surface and self.quaFaces[0].is_surface) or not only_surface:
            self.triangles.append([self.nodeIds[3], self.nodeIds[8], self.nodeIds[7]])
            self.triangles.append([self.nodeIds[7], self.nodeIds[8], self.nodeIds[5]])
            self.triangles.append([self.nodeIds[7], self.nodeIds[5], self.nodeIds[6]])
            self.triangles.append([self.nodeIds[7], self.nodeIds[6], self.nodeIds[2]])
            self.triangles.append([self.nodeIds[8], self.nodeIds[0], self.nodeIds[5]])
            self.triangles.append([self.nodeIds[6], self.nodeIds[5], self.nodeIds[1]])

            self.edges.extend([self.nodeIds[0], self.nodeIds[1]])
            self.edges.extend([self.nodeIds[1], self.nodeIds[2]])
            self.edges.extend([self.nodeIds[2], self.nodeIds[3]])
            self.edges.extend([self.nodeIds[3], self.nodeIds[0]])

        if (only_surface and self.triFaces[0].is_surface) or not only_surface:
            self.triangles.append([self.nodeIds[4], self.nodeIds[12], self.nodeIds[9]])
            self.triangles.append([self.nodeIds[12], self.nodeIds[3], self.nodeIds[8]])
            self.triangles.append([self.nodeIds[12], self.nodeIds[8], self.nodeIds[9]])
            self.triangles.append([self.nodeIds[9], self.nodeIds[8], self.nodeIds[0]])

            self.edges.extend([self.nodeIds[0], self.nodeIds[3]])
            self.edges.extend([self.nodeIds[0], self.nodeIds[4]])
            self.edges.extend([self.nodeIds[3], self.nodeIds[4]])

        if (only_surface and self.triFaces[1].is_surface) or not only_surface:
            self.triangles.append([self.nodeIds[4], self.nodeIds[9], self.nodeIds[10]])
            self.triangles.append([self.nodeIds[9], self.nodeIds[0], self.nodeIds[5]])
            self.triangles.append([self.nodeIds[9], self.nodeIds[5], self.nodeIds[10]])
            self.triangles.append([self.nodeIds[10], self.nodeIds[5], self.nodeIds[1]])

            self.edges.extend([self.nodeIds[0], self.nodeIds[1]])
            self.edges.extend([self.nodeIds[4], self.nodeIds[1]])
            self.edges.extend([self.nodeIds[0], self.nodeIds[4]])

        if (only_surface and self.triFaces[2].is_surface) or not only_surface:
            self.triangles.append([self.nodeIds[4], self.nodeIds[10], self.nodeIds[11]])
            self.triangles.append([self.nodeIds[10], self.nodeIds[1], self.nodeIds[6]])
            self.triangles.append([self.nodeIds[10], self.nodeIds[6], self.nodeIds[11]])
            self.triangles.append([self.nodeIds[11], self.nodeIds[6], self.nodeIds[2]])

            self.edges.extend([self.nodeIds[4], self.nodeIds[2]])
            self.edges.extend([self.nodeIds[4], self.nodeIds[1]])
            self.edges.extend([self.nodeIds[1], self.nodeIds[2]])

        if (only_surface and self.triFaces[3].is_surface) or not only_surface:
            self.triangles.append([self.nodeIds[4], self.nodeIds[11], self.nodeIds[12]])
            self.triangles.append([self.nodeIds[11], self.nodeIds[2], self.nodeIds[7]])
            self.triangles.append([self.nodeIds[11], self.nodeIds[7], self.nodeIds[12]])
            self.triangles.append([self.nodeIds[12], self.nodeIds[7], self.nodeIds[3]])

            self.edges.extend([self.nodeIds[4], self.nodeIds[2]])
            self.edges.extend([self.nodeIds[4], self.nodeIds[3]])
            self.edges.extend([self.nodeIds[2], self.nodeIds[3]])

        return self.triangles


class MeshC3D5(MeshElement):
    """
    金字塔单元
    """

    def __init__(self, _id):
        super().__init__()
        self.id = _id
        self.nodesCount = 5
        self.triangles = []
        self.upgradeNodeIndex = [0, 1, 2, 3, 4]

    def setFaces(self, node_ids):
        self.nodeIds = node_ids
        self.triangles.clear()
        self.quaFaces.append(QuaFace([node_ids[0], node_ids[1], node_ids[2], node_ids[3]]))
        self.triFaces.append(TriFace([node_ids[0], node_ids[1], node_ids[4]]))
        self.triFaces.append(TriFace([node_ids[4], node_ids[1], node_ids[2]]))
        self.triFaces.append(TriFace([node_ids[4], node_ids[2], node_ids[3]]))
        self.triFaces.append(TriFace([node_ids[4], node_ids[3], node_ids[0]]))

    def getAllTriangles(self, only_surface=False):
        self.triangles.clear()
        if (only_surface and self.quaFaces[0].is_surface) or not only_surface:
            self.triangles.append([self.nodeIds[1], self.nodeIds[0], self.nodeIds[3]])
            self.triangles.append([self.nodeIds[1], self.nodeIds[3], self.nodeIds[2]])

            self.edges.extend([self.nodeIds[0], self.nodeIds[1]])
            self.edges.extend([self.nodeIds[1], self.nodeIds[2]])
            self.edges.extend([self.nodeIds[2], self.nodeIds[3]])
            self.edges.extend([self.nodeIds[3], self.nodeIds[0]])

        if (only_surface and self.triFaces[0].is_surface) or not only_surface:
            self.triangles.append([self.nodeIds[0], self.nodeIds[1], self.nodeIds[4]])
            self.edges.extend([self.nodeIds[0], self.nodeIds[1]])
            self.edges.extend([self.nodeIds[4], self.nodeIds[1]])
            self.edges.extend([self.nodeIds[0], self.nodeIds[4]])

        if (only_surface and self.triFaces[1].is_surface) or not only_surface:
            self.triangles.append([self.nodeIds[4], self.nodeIds[1], self.nodeIds[2]])
            self.edges.extend([self.nodeIds[4], self.nodeIds[2]])
            self.edges.extend([self.nodeIds[4], self.nodeIds[1]])
            self.edges.extend([self.nodeIds[1], self.nodeIds[2]])

        if (only_surface and self.triFaces[2].is_surface) or not only_surface:
            self.triangles.append([self.nodeIds[4], self.nodeIds[2], self.nodeIds[3]])
            self.edges.extend([self.nodeIds[4], self.nodeIds[2]])
            self.edges.extend([self.nodeIds[4], self.nodeIds[3]])
            self.edges.extend([self.nodeIds[2], self.nodeIds[3]])

        if (only_surface and self.triFaces[3].is_surface) or not only_surface:
            self.triangles.append([self.nodeIds[4], self.nodeIds[3], self.nodeIds[0]])
            self.edges.extend([self.nodeIds[0], self.nodeIds[3]])
            self.edges.extend([self.nodeIds[0], self.nodeIds[4]])
            self.edges.extend([self.nodeIds[3], self.nodeIds[4]])

        return self.triangles


class MeshC3D20(MeshElement):
    """
    20节点六面体
    """

    def __init__(self, _id):
        super().__init__()
        self.id = _id
        self.nodesCount = 20
        self.triangles = []

    def setFaces(self, _nodeIds):
        self.nodeIds = _nodeIds
        f1 = QuaFace([_nodeIds[0], _nodeIds[3], _nodeIds[2], _nodeIds[1]])
        f2 = QuaFace([_nodeIds[4], _nodeIds[5], _nodeIds[6], _nodeIds[7]])
        f3 = QuaFace([_nodeIds[5], _nodeIds[4], _nodeIds[0], _nodeIds[1]])
        f4 = QuaFace([_nodeIds[6], _nodeIds[2], _nodeIds[3], _nodeIds[7]])
        f5 = QuaFace([_nodeIds[1], _nodeIds[2], _nodeIds[6], _nodeIds[5]])
        f6 = QuaFace([_nodeIds[4], _nodeIds[7], _nodeIds[3], _nodeIds[0]])
        self.quaFaces = [f1, f2, f3, f4, f5, f6]

    def getAllTriangles(self, only_surface=False):
        self.triangles.clear()

        if (only_surface and self.quaFaces[0].is_surface) or not only_surface:
            self.triangles.extend(QuaFace2Triangles(self.nodeIds, [0, 8, 1, 11, 9, 3, 10, 2]))
            self.edges.extend([self.nodeIds[0], self.nodeIds[1]])
            self.edges.extend([self.nodeIds[1], self.nodeIds[2]])
            self.edges.extend([self.nodeIds[2], self.nodeIds[3]])
            self.edges.extend([self.nodeIds[3], self.nodeIds[0]])

        if (only_surface and self.quaFaces[1].is_surface) or not only_surface:
            self.triangles.extend(QuaFace2Triangles(self.nodeIds, [5, 16, 4, 17, 19, 6, 18, 7]))
            self.edges.extend([self.nodeIds[4], self.nodeIds[5]])
            self.edges.extend([self.nodeIds[5], self.nodeIds[6]])
            self.edges.extend([self.nodeIds[6], self.nodeIds[7]])
            self.edges.extend([self.nodeIds[7], self.nodeIds[4]])

        if (only_surface and self.quaFaces[2].is_surface) or not only_surface:
            self.triangles.extend(QuaFace2Triangles(self.nodeIds, [0, 12, 4, 8, 16, 1, 13, 5]))
            self.edges.extend([self.nodeIds[0], self.nodeIds[4]])
            self.edges.extend([self.nodeIds[4], self.nodeIds[5]])
            self.edges.extend([self.nodeIds[5], self.nodeIds[1]])
            self.edges.extend([self.nodeIds[1], self.nodeIds[0]])

        if (only_surface and self.quaFaces[3].is_surface) or not only_surface:
            self.triangles.extend(QuaFace2Triangles(self.nodeIds, [2, 14, 6, 10, 18, 3, 15, 7]))
            self.edges.extend([self.nodeIds[2], self.nodeIds[6]])
            self.edges.extend([self.nodeIds[6], self.nodeIds[7]])
            self.edges.extend([self.nodeIds[7], self.nodeIds[3]])
            self.edges.extend([self.nodeIds[3], self.nodeIds[2]])

        if (only_surface and self.quaFaces[4].is_surface) or not only_surface:
            self.triangles.extend(QuaFace2Triangles(self.nodeIds, [1, 13, 5, 9, 17, 2, 14, 6]))
            self.edges.extend([self.nodeIds[1], self.nodeIds[5]])
            self.edges.extend([self.nodeIds[5], self.nodeIds[6]])
            self.edges.extend([self.nodeIds[6], self.nodeIds[2]])
            self.edges.extend([self.nodeIds[2], self.nodeIds[1]])

        if (only_surface and self.quaFaces[5].is_surface) or not only_surface:
            self.triangles.extend(QuaFace2Triangles(self.nodeIds, [4, 12, 0, 19, 11, 7, 15, 3]))
            self.edges.extend([self.nodeIds[0], self.nodeIds[3]])
            self.edges.extend([self.nodeIds[3], self.nodeIds[7]])
            self.edges.extend([self.nodeIds[7], self.nodeIds[4]])
            self.edges.extend([self.nodeIds[4], self.nodeIds[0]])

        return self.triangles


class MeshCQUAD4(MeshElement):
    """
    四节点壳单元, 退化的单元为三节点壳单元, 默认壳单元全是表面单元
    """

    def __init__(self, _id):
        super().__init__()
        self.id = _id
        self.nodesCount = 4
        self.degrade_ele = 30500
        self.degrade_rv_index = [2]
        self.triangles = []

    def setFaces(self, node_ids):
        self.nodeIds = node_ids
        self.triangles.clear()

        triangle = [node_ids[0], node_ids[1], node_ids[2]]
        self.triangles.append(triangle)

        triangle = [node_ids[0], node_ids[2], node_ids[3]]
        self.triangles.append(triangle)

        self.edges = [node_ids[0], node_ids[1], node_ids[1], node_ids[2], node_ids[2], node_ids[3], node_ids[3], node_ids[0]]

    def getAllTriangles(self, only_surface=False):
        return self.triangles


class MeshC3D8(MeshElement):
    def __init__(self, _id):
        super().__init__()
        self.setId(_id)
        self.nodesCount = 8
        self.degradeEle = 60600
        self.degradeRvIndex = [4, 0]
        self.upgradeNodeIndex = [0, 1, 2, 3, 4, 5, 6, 7]

    def setFaces(self, _nodeIds):
        self.nodeIds = _nodeIds
        f1 = QuaFace([_nodeIds[0], _nodeIds[3], _nodeIds[2], _nodeIds[1]])
        f2 = QuaFace([_nodeIds[4], _nodeIds[5], _nodeIds[6], _nodeIds[7]])
        f3 = QuaFace([_nodeIds[5], _nodeIds[4], _nodeIds[0], _nodeIds[1]])
        f4 = QuaFace([_nodeIds[6], _nodeIds[2], _nodeIds[3], _nodeIds[7]])
        f5 = QuaFace([_nodeIds[1], _nodeIds[2], _nodeIds[6], _nodeIds[5]])
        f6 = QuaFace([_nodeIds[4], _nodeIds[7], _nodeIds[3], _nodeIds[0]])
        self.quaFaces = [f1, f2, f3, f4, f5, f6]

    def getAllTriangles(self, only_surface=False):
        self.triangles = []  # Reset triangles to an empty list
        for quaFace in self.quaFaces:
            # First triangle from the quadrilateral face
            if (only_surface and quaFace.is_surface) or not only_surface:
                triangle1 = [quaFace.nodes[0], quaFace.nodes[1], quaFace.nodes[2]]
                self.triangles.append(triangle1)
                # Second triangle from the quadrilateral face
                triangle2 = [quaFace.nodes[0], quaFace.nodes[2], quaFace.nodes[3]]
                self.triangles.append(triangle2)

                self.edges.extend([quaFace.nodes[0], quaFace.nodes[1]])
                self.edges.extend([quaFace.nodes[1], quaFace.nodes[2]])
                self.edges.extend([quaFace.nodes[2], quaFace.nodes[3]])
                self.edges.extend([quaFace.nodes[3], quaFace.nodes[0]])
        return self.triangles


class MeshC3D10(MeshElement):
    """
    ABAQUS十节点四面体, triangles是为了web端显示做的全部小面片
    """

    def __init__(self, _id):
        super().__init__()
        self.setId(_id)
        self.nodesCount = 10

    def setFaces(self, _nodeIds):
        self.nodeIds = _nodeIds
        self.triFaces.append(TriFace([_nodeIds[0], _nodeIds[1], _nodeIds[2]]))
        self.triFaces.append(TriFace([_nodeIds[0], _nodeIds[2], _nodeIds[3]]))
        self.triFaces.append(TriFace([_nodeIds[0], _nodeIds[1], _nodeIds[3]]))
        self.triFaces.append(TriFace([_nodeIds[1], _nodeIds[2], _nodeIds[3]]))

    def getAllTriangles(self, only_surface=False):
        self.triangles.clear()
        if (only_surface and self.triFaces[0].is_surface) or not only_surface:
            fe = [self.nodeIds[0], self.nodeIds[6], self.nodeIds[4]]
            ff = [self.nodeIds[6], self.nodeIds[2], self.nodeIds[5]]
            fg = [self.nodeIds[6], self.nodeIds[5], self.nodeIds[4]]
            fh = [self.nodeIds[4], self.nodeIds[5], self.nodeIds[1]]
            self.triangles.extend([fe, ff, fg, fh])

            self.edges.extend([self.nodeIds[0], self.nodeIds[1]])
            self.edges.extend([self.nodeIds[1], self.nodeIds[2]])
            self.edges.extend([self.nodeIds[2], self.nodeIds[0]])

        if (only_surface and self.triFaces[1].is_surface) or not only_surface:
            fa = [self.nodeIds[9], self.nodeIds[2], self.nodeIds[6]]
            fb = [self.nodeIds[3], self.nodeIds[9], self.nodeIds[7]]
            fc = [self.nodeIds[7], self.nodeIds[9], self.nodeIds[6]]
            fd = [self.nodeIds[0], self.nodeIds[7], self.nodeIds[6]]
            self.triangles.extend([fa, fb, fc, fd])

            self.edges.extend([self.nodeIds[0], self.nodeIds[2]])
            self.edges.extend([self.nodeIds[2], self.nodeIds[3]])
            self.edges.extend([self.nodeIds[3], self.nodeIds[0]])

        if (only_surface and self.triFaces[2].is_surface) or not only_surface:
            f1 = [self.nodeIds[0], self.nodeIds[4], self.nodeIds[7]]
            f2 = [self.nodeIds[4], self.nodeIds[1], self.nodeIds[8]]
            f3 = [self.nodeIds[4], self.nodeIds[8], self.nodeIds[7]]
            f4 = [self.nodeIds[8], self.nodeIds[3], self.nodeIds[7]]
            self.triangles.extend([f1, f2, f3, f4])

            self.edges.extend([self.nodeIds[0], self.nodeIds[1]])
            self.edges.extend([self.nodeIds[1], self.nodeIds[3]])
            self.edges.extend([self.nodeIds[3], self.nodeIds[0]])

        if (only_surface and self.triFaces[3].is_surface) or not only_surface:
            f5 = [self.nodeIds[1], self.nodeIds[5], self.nodeIds[8]]
            f6 = [self.nodeIds[5], self.nodeIds[2], self.nodeIds[9]]
            f7 = [self.nodeIds[5], self.nodeIds[9], self.nodeIds[8]]
            f8 = [self.nodeIds[9], self.nodeIds[3], self.nodeIds[8]]
            self.triangles.extend([f5, f6, f7, f8])

            self.edges.extend([self.nodeIds[1], self.nodeIds[2]])
            self.edges.extend([self.nodeIds[2], self.nodeIds[3]])
            self.edges.extend([self.nodeIds[3], self.nodeIds[1]])

        return self.triangles


class MeshTRIA3(MeshElement):
    """
    三角形壳单元, 默认壳单元全是表面单元
    """

    def __init__(self, _id):
        super().__init__()
        self.setId(_id)
        self.nodesCount = 3

    def setFaces(self, nodeIds):
        self.nodeIds = nodeIds
        self.triangles = [[nodeIds[0], nodeIds[1], nodeIds[2]]]

        self.edges = [nodeIds[0], nodeIds[1], nodeIds[1], nodeIds[2], nodeIds[2], nodeIds[0]]

    def getAllTriangles(self, only_surface=False):
        return self.triangles


class MeshTetra(MeshElement):
    """
    四面体单元
    """

    def __init__(self, _id):
        super().__init__()
        self.setId(_id)
        self.nodesCount = 4
        self.upgradeNodeIndex = [0, 1, 2, 3]

    def setFaces(self, _nodeIds):
        self.nodeIds = _nodeIds
        self.triFaces.append(TriFace([_nodeIds[0], _nodeIds[1], _nodeIds[2]]))
        self.triFaces.append(TriFace([_nodeIds[0], _nodeIds[1], _nodeIds[3]]))
        self.triFaces.append(TriFace([_nodeIds[0], _nodeIds[2], _nodeIds[3]]))
        self.triFaces.append(TriFace([_nodeIds[1], _nodeIds[2], _nodeIds[3]]))

    def getAllTriangles(self, only_surface=False):
        self.triangles.clear()
        if (only_surface and self.triFaces[0].is_surface) or not only_surface:
            self.triangles.append([self.nodeIds[0], self.nodeIds[1], self.nodeIds[2]])
            self.edges.extend([self.nodeIds[0], self.nodeIds[1]])
            self.edges.extend([self.nodeIds[1], self.nodeIds[2]])
            self.edges.extend([self.nodeIds[2], self.nodeIds[0]])

        if (only_surface and self.triFaces[1].is_surface) or not only_surface:
            self.triangles.append([self.nodeIds[0], self.nodeIds[1], self.nodeIds[3]])
            self.edges.extend([self.nodeIds[0], self.nodeIds[1]])
            self.edges.extend([self.nodeIds[1], self.nodeIds[3]])
            self.edges.extend([self.nodeIds[3], self.nodeIds[0]])

        if (only_surface and self.triFaces[2].is_surface) or not only_surface:
            self.triangles.append([self.nodeIds[0], self.nodeIds[2], self.nodeIds[3]])
            self.edges.extend([self.nodeIds[0], self.nodeIds[2]])
            self.edges.extend([self.nodeIds[2], self.nodeIds[3]])
            self.edges.extend([self.nodeIds[3], self.nodeIds[0]])

        if (only_surface and self.triFaces[3].is_surface) or not only_surface:
            self.triangles.append([self.nodeIds[1], self.nodeIds[2], self.nodeIds[3]])
            self.edges.extend([self.nodeIds[1], self.nodeIds[2]])
            self.edges.extend([self.nodeIds[2], self.nodeIds[3]])
            self.edges.extend([self.nodeIds[3], self.nodeIds[1]])

        return self.triangles


class MeshC3D15(MeshElement):
    """
    15节点三棱柱
    """

    def __init__(self, _id):
        super().__init__()
        self.setId(_id)
        self.nodesCount = 15

    def setFaces(self, _nodeIds):
        self.nodeIds = _nodeIds
        self.triFaces.append(TriFace([_nodeIds[0], _nodeIds[1], _nodeIds[2]]))
        self.triFaces.append(TriFace([_nodeIds[3], _nodeIds[4], _nodeIds[5]]))
        self.quaFaces.append(QuaFace([_nodeIds[0], _nodeIds[1], _nodeIds[4], _nodeIds[3]]))
        self.quaFaces.append(QuaFace([_nodeIds[1], _nodeIds[4], _nodeIds[5], _nodeIds[2]]))
        self.quaFaces.append(QuaFace([_nodeIds[0], _nodeIds[2], _nodeIds[5], _nodeIds[3]]))

    def getAllTriangles(self, only_surface=False):
        self.triangles.clear()
        if (only_surface and self.triFaces[0].is_surface) or not only_surface:
            f1 = [self.nodeIds[0], self.nodeIds[8], self.nodeIds[6]]
            f2 = [self.nodeIds[1], self.nodeIds[6], self.nodeIds[7]]
            f3 = [self.nodeIds[6], self.nodeIds[8], self.nodeIds[7]]
            f4 = [self.nodeIds[2], self.nodeIds[8], self.nodeIds[7]]
            self.triangles.extend([f1, f2, f3, f4])
            self.edges.extend([self.nodeIds[0], self.nodeIds[1]])
            self.edges.extend([self.nodeIds[1], self.nodeIds[2]])
            self.edges.extend([self.nodeIds[2], self.nodeIds[0]])

        if (only_surface and self.triFaces[1].is_surface) or not only_surface:
            f1 = [self.nodeIds[5], self.nodeIds[11], self.nodeIds[10]]
            f2 = [self.nodeIds[11], self.nodeIds[3], self.nodeIds[9]]
            f3 = [self.nodeIds[10], self.nodeIds[9], self.nodeIds[4]]
            f4 = [self.nodeIds[10], self.nodeIds[11], self.nodeIds[9]]
            self.triangles.extend([f1, f2, f3, f4])
            self.edges.extend([self.nodeIds[3], self.nodeIds[4]])
            self.edges.extend([self.nodeIds[4], self.nodeIds[5]])
            self.edges.extend([self.nodeIds[5], self.nodeIds[3]])

        if (only_surface and self.quaFaces[0].is_surface) or not only_surface:
            f1 = [self.nodeIds[9], self.nodeIds[3], self.nodeIds[12]]
            f2 = [self.nodeIds[4], self.nodeIds[9], self.nodeIds[13]]
            f3 = [self.nodeIds[9], self.nodeIds[12], self.nodeIds[6]]
            f4 = [self.nodeIds[9], self.nodeIds[6], self.nodeIds[13]]
            f5 = [self.nodeIds[12], self.nodeIds[0], self.nodeIds[6]]
            f6 = [self.nodeIds[13], self.nodeIds[6], self.nodeIds[1]]
            self.triangles.extend([f1, f2, f3, f4, f5, f6])
            self.edges.extend([self.nodeIds[0], self.nodeIds[1]])
            self.edges.extend([self.nodeIds[1], self.nodeIds[3]])
            self.edges.extend([self.nodeIds[3], self.nodeIds[4]])
            self.edges.extend([self.nodeIds[4], self.nodeIds[0]])

        if (only_surface and self.quaFaces[1].is_surface) or not only_surface:
            f1 = [self.nodeIds[11], self.nodeIds[5], self.nodeIds[14]]
            f2 = [self.nodeIds[3], self.nodeIds[11], self.nodeIds[12]]
            f3 = [self.nodeIds[11], self.nodeIds[14], self.nodeIds[8]]
            f4 = [self.nodeIds[11], self.nodeIds[8], self.nodeIds[12]]
            f5 = [self.nodeIds[14], self.nodeIds[2], self.nodeIds[8]]
            f6 = [self.nodeIds[12], self.nodeIds[8], self.nodeIds[0]]
            self.triangles.extend([f1, f2, f3, f4, f5, f6])
            self.edges.extend([self.nodeIds[0], self.nodeIds[2]])
            self.edges.extend([self.nodeIds[2], self.nodeIds[5]])
            self.edges.extend([self.nodeIds[5], self.nodeIds[3]])
            self.edges.extend([self.nodeIds[3], self.nodeIds[0]])

        if (only_surface and self.quaFaces[1].is_surface) or not only_surface:
            f1 = [self.nodeIds[10], self.nodeIds[4], self.nodeIds[13]]
            f2 = [self.nodeIds[5], self.nodeIds[10], self.nodeIds[14]]
            f3 = [self.nodeIds[10], self.nodeIds[13], self.nodeIds[7]]
            f4 = [self.nodeIds[10], self.nodeIds[7], self.nodeIds[14]]
            f5 = [self.nodeIds[13], self.nodeIds[1], self.nodeIds[7]]
            f6 = [self.nodeIds[14], self.nodeIds[7], self.nodeIds[2]]
            self.triangles.extend([f1, f2, f3, f4, f5, f6])
            self.edges.extend([self.nodeIds[1], self.nodeIds[2]])
            self.edges.extend([self.nodeIds[2], self.nodeIds[5]])
            self.edges.extend([self.nodeIds[5], self.nodeIds[4]])
            self.edges.extend([self.nodeIds[4], self.nodeIds[1]])

        return self.triangles
