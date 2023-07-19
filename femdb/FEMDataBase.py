#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from utils.Singleton import Singleton
from femdb.LoadCase import *
from femdb.Material import *
from element.Node import *
from femdb.Sets import *
from femdb.Property import *
from femdb.Section import *
from femdb.ElementGroup import *
from collections import OrderedDict
from scipy import sparse


@Singleton
class FEMDataBase(object):
    """
    有限元数据库, 实例化的都存储在这里
    """

    def __init__(self):
        # 输入文件
        self.file_path = None
        # nodes
        self.node_list = []  # List of all nodes in the domain, 实例化数据
        self.node_hash = {}  # 节点真实Id对应nodelist中的index的Hash表

        # elements
        # Dict of all Element in the domain, key: ele_keyword, value: ele_group 单元实际数据也存在这里
        self.ele_grp_hash = {}
        self.biggest_grp = ""  # 最大的组, 含义是哪个组内包含的单元最多
        self.et_hash = {}
        self.equation_number = None

        # 单元真实id对应group_hash中对应单元组中的index, (key:real_id) => (value: index), 所以存在多对一的情况
        self.ele_idx_hash = {}
        self.ele_count = 0  # 单元个数

        # Preprocess Sets
        self.node_sets = []
        self.ele_sets = []

        # Preprocess Fem
        self.properties = []
        self.sections = []
        self.materials = []
        self.global_stiff_matrix = None
        self.load_case = LoadCase()
        self.real_const_hash = {}

    """ 
    以下的函数为解析文件的相关函数, 添加节点、单元、节点集、单元集、属性、材料、边界条件、LoadCase等 
    """

    def AddNode(self, node):
        """
        向有限元模型中插入节点
        :param node: 节点
        """
        self.node_list.append(node)

    def GetNodeBySearchId(self, n_id):
        return self.node_list[n_id]

    def SetNodeHashTable(self, n_dict):
        self.node_hash = n_dict

    def SetGrpHash(self, ele_grp, ele_count):
        self.ele_grp_hash = ele_grp
        self.ele_count = ele_count

    def GetProperty(self):
        return self.properties

    def GetEleTypeById(self, eid):
        """
        通过单元的id获取单元类型, 从而确定一个set中对应的单元类型
        """

    def GetSpecificFEMObject(self, obj_type, obj_name):
        if obj_type == FEMObject.NodeSet:
            for ns in self.node_sets:
                if ns.GetName() == obj_name:
                    return ns
        elif obj_type == FEMObject.EleSet:
            for es in self.ele_sets:
                if es.GetName() == obj_name:
                    return es
        elif obj_type == FEMObject.Material:
            for mat in self.materials:
                if mat.GetName() == obj_name:
                    return mat
        elif obj_type == FEMObject.Section:
            for sec in self.sections:
                if sec.GetName() == obj_name:
                    return sec

        else:
            print("Fatal Error: Can't find object: {} with type: {}".format(obj_name, obj_type))
            sys.exit(1)

    def PrintParseSummary(self):
        """
        打印Inp解析结果, 至此文件解析结束, 下一步将具体的属性特征赋予至单元上
        """
        # Print Node And Element Summary
        parse_summary = "\nFinish Parse File, Summary is:\n  Nodes Count: {}\n".format(len(self.node_list))
        for key, ele_group in self.ele_grp_hash.items():
            iter_ele_summary = "  Element {} Count is: {}\n".format(key, ele_group.GetElementsCurrentCount())
            parse_summary += iter_ele_summary
        parse_summary += "  Total Element Count is: {}\n".format(self.ele_count)

        # Print NodeSet & ElementSet
        parse_summary += "\n  Here is NodeSets & EleSets:\n"
        parse_summary += "  Node Sets count is: {}\n".format(len(self.node_sets))
        for ns in self.node_sets:
            parse_summary += "    {}\n".format(ns)
        parse_summary += "  Elements Sets count is: {}\n".format(len(self.ele_sets))
        for es in self.ele_sets:
            parse_summary += "    {}\n".format(es)

        # Print Material & Section Summary
        parse_summary += "\n  Here is Material & Section:\n"
        parse_summary += "  Material's count: {}\n".format(len(self.materials))
        for mat in self.materials:
            parse_summary += "    {}\n".format(mat)
        parse_summary += "  Section's count: {}\n".format(len(self.properties))
        for section in self.properties:
            parse_summary += "    {}\n".format(section)

        # Print LoadCases Summary
        parse_summary += "  {}\n".format(self.load_case)
        mlogger.info(parse_summary)

    def GetModelSummary(self):
        """
        获取解析模型的信息, 包括文件位置、分析类型、单元个数、节点个数、EquationNumber(总自由度个数减去约束)
        :return:
        """
        summary_dict = OrderedDict()
        summary_dict["File Path"] = str(self.file_path)
        summary_dict["Number Of Equation"] = self.equation_number
        summary_dict["Number Of Node"] = len(self.node_list)  # TODO: 并不是标准的节点个数, 标准节点个数应该是由单元计算出来
        summary_dict["Number Of Element"] = self.ele_count

        return summary_dict

    def GetNodeCoordBySearchId(self, idxes):
        """
        获取节点的坐标, 用于计算单元刚度阵
        """
        coords = []
        for idx in idxes:
            t_node = self.node_list[idx]
            coords.append((t_node.x, t_node.y, t_node.z))
        return coords

    def AssignElementProperty(self):
        if GlobalInfor[GlobalVariant.InputFileSuffix] == InputFileType.CDB:
            self.AssignElementPropertyAnsys()
        elif GlobalInfor[GlobalVariant.InputFileSuffix] == InputFileType.INP:
            self.AssignElementPropertyAbaqus()
        elif GlobalInfor[GlobalVariant.InputFileSuffix] == InputFileType.BDF:
            self.AssignElementPropertyNastran()
        else:
            mlogger.fatal("UnSupport Input File Type:{}".format(self.input_file_type))
            sys.exit(1)

    def AssignElementPropertyNastran(self):
        pass

    def AssignElementPropertyAnsys(self):
        """
        将属性分配给对应的单元, 适配Ansys cdb格式
        """
        # 计算解析完输入文件后暂未计算的量
        for sec in self.sections:
            sec_type = sec.sec_type
            if sec.IsBeamSection:
                inertia_character = BeamCalculator.CalculateMomentOfInertiaOfArea(sec_type, sec.sec_data)
                area_character = BeamCalculator.CalEffectiveShearArea(sec_type, sec.sec_data)
                sec.SetSectionCharacter({**inertia_character, **area_character})
            else:
                mlogger.fatal("UnSupport Element Type")
                sys.exit(1)

        """
        对每个单元进行循环, 赋予属性, 计算单刚. 假设ANSYS的单元信息全部存储在了mat、sec、real_const三者中
        """
        for e_type, grp in self.ele_grp_hash.items():
            for iter_ele in grp.Elements():
                mat_dict = self.GetSpecificFEMObject(FEMObject.Material, iter_ele.mat_id).GetValueDict()

                # 并不是所有单元都有截面属性, 首先要判断是否为空
                sec_characters = self.GetSpecificFEMObject(FEMObject.Section, iter_ele.sec_id)
                if sec_characters is not None:
                    sec_characters = sec_characters.GetSectionCharacter()
                else:
                    sec_characters = {}

                # 同样, 也并不是所有单元都有实常数, 首先判断是否为空
                real_const = {"RealConst": self.real_const_hash[iter_ele.real_const_id]}

                # 所有计算单刚的参数均已设置完毕, 可以计算单刚
                iter_ele.SetAllCharacterAndCalD({**mat_dict, **sec_characters, **real_const})

    def AssignElementPropertyAbaqus(self):
        """
        将属性分配给对应的单元, 适配Abaqus inp格式
        """
        # 大部分情况是模型中绝大部分的单元是同一类型, 只有少量的单元是其他类型, 所以首先找出包含单元最多的grp, 再确定了
        # 目标单元ID不属于其他包含少量单元的grp之后, 即可确定属于包含单元最多的grp, 从而确定单元类型. 但是有一点需要注
        # 意的是前提假设是所有的单元都在grp里
        if self.biggest_grp == "":
            max_ele_count = 0
            for e_type, grp in self.ele_grp_hash.items():
                if grp.GetElementsCurrentCount() > max_ele_count:
                    max_ele_count = grp.GetElementsCurrentCount()
                    self.biggest_grp = e_type

        # 按属性将指定的set内的单元分配属性, set内单元可能属于不同的类型, 也就是属于不同的grp, 所以首先要确定单元类型
        for ele_property in self.properties:
            eset = self.GetSpecificFEMObject(FEMObject.EleSet, ele_property.GetEleSetName())
            mat = self.GetSpecificFEMObject(FEMObject.Material, ele_property.GetMatName())
            mat_dict = mat.GetValueDict()  # 材料参数用dict描述
            prop_dict = ele_property.GetPropertyPars()  # 单元的属性参数, 比如厚度

            for e_id in eset.GetEleIds():
                # 根据ele_type和list_idx就可以锁定单元, 其存储了单元所有信息, 而不只是单元号
                ele_type = self.GetElementTypeByID(e_id)
                list_idx = self.ele_idx_hash[e_id]
                assert ele_type != ""
                cur_ele = self.ele_grp_hash[ele_type].Elements()[list_idx]

                # 材料属性设置完成, 计算单元的D阵和B阵, 从而计算单元的刚度阵
                cur_ele.SetAllCharacterAndCalD(mat_dict, prop_dict)

    def GetElementTypeByID(self, eid):
        """
        根据单元真实ID获取单元类型, 因为set中可能包含不同类型的单元, 而不同类型的单元存储在不同的group里
        :param eid: 单元的真实ID
        :return: 单元的类型
        """
        ele_type = self.biggest_grp
        for iter_type, ele_grp in self.ele_grp_hash.items():
            if iter_type == self.biggest_grp:
                continue
            if ele_grp.IsElementInGroup(eid):
                ele_type = ele_grp.e_type
                break

        return ele_type

    def InitAssemblyMatrix(self, eq_count):
        """
        初始化总刚
        :param eq_count:
        """
        self.global_stiff_matrix = sparse.coo_matrix((eq_count, eq_count), dtype=float)
