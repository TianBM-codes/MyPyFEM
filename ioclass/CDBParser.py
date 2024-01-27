#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from femdb.FEMDataBase import *
import fortranformat as ff
import copy
from collections import OrderedDict


class CDBParser(object):
    """
    功能为解析CDB文件, 初始化数据库单元节点及相关信息, 只负责解析和传递, 其他功能均在其他类中,
    同样的, 其他类也不会包含任何解析输入文件的功能或者函数
    """

    def __init__(self, input_path):
        self.femdb = FEMDataBase()
        self.cdb_path = input_path
        self.iter_line = None
        self.et_hash = {}
        self.ele_group_hash = {}
        self.ele_count = 0
        self.real_constant_hash = {}

    def ParseFileAndInitFEMDB(self):
        """
        解析ANSYS软件的CDB文件
        readline()不用strip(), 在具体的splits中再strip()
        关于整体的格式：
        1. 字母开头的定义为关键字
        2. 括号开头是fortran格式符号, 注意可能会连续出现两行格式符号
        3. *开头的是其他相关, 暂时不予解析
        4. /开头的跟着的是注释

        Reference:
        1. D:\\software\\python\\setup394\\Lib\\site-packages\\ansys\\mapdl\\reader
        2. ANSYS Help - Mechanical APDL Element Reference
        """
        self.femdb.file_path = self.cdb_path
        GlobalInfor[GlobalVariant.AnaType] = AnalyseType.LinearStatic  # 默认静力分析
        with open(self.cdb_path, 'r') as cdb_f:
            self.iter_line = cdb_f.readline()
            while True:
                # 解析分析类型
                if self.iter_line.startswith("ANTYPE,"):
                    if self.iter_line.split(",")[1].strip() == "0":
                        GlobalInfor[GlobalVariant.AnaType] = AnalyseType.LinearStatic
                    else:
                        mlogger.fatal("UnSupport Analyse Type:{}".format(self.iter_line))
                        sys.exit(1)
                    self.iter_line = cdb_f.readline()

                # 解析单元类型 TODO: 解析KeyOption
                elif self.iter_line.startswith("ET,"):
                    splits = self.iter_line.split(",")
                    assert len(splits) == 3
                    self.et_hash[int(splits[1].strip())] = int(splits[2].strip())
                    self.iter_line = cdb_f.readline()

                # 解析实常数, 本行暂时不解析内容, 直接解析下一行开始的实常数
                elif self.iter_line.startswith("RLBLOCK,"):
                    self.iter_line = cdb_f.readline()
                    # 如果不是以字母开头的, 那么还没有跳出定义
                    format_list = []
                    while not self.iter_line[0].isalpha():
                        if self.iter_line[0] == "(":
                            format_list.append(ff.FortranRecordReader(self.iter_line.strip()))
                            self.iter_line = cdb_f.readline()
                        elif self.iter_line[0] == " ":  # RLBLOCK
                            rl_data = format_list.pop(0).read(self.iter_line)
                            self.real_constant_hash[rl_data[0]] = rl_data[2:]
                            self.iter_line = cdb_f.readline()

                # 解析节点信息, TODO:平面应变平面应力这种只有二维坐标的
                elif self.iter_line.startswith("NBLOCK,"):
                    node_index = 0  # 相当于节点个数, 也是对应数据库中node_list中的index
                    fortran_format = cdb_f.readline().strip()  # skip format of node line
                    f_reader = ff.FortranRecordReader(fortran_format)
                    self.iter_line = cdb_f.readline()  # 不需要strip, 否则前面的空格会滤掉

                    while self.iter_line.startswith(" "):
                        nd_data = f_reader.read(self.iter_line)
                        n_id = int(nd_data[0])
                        x, y, z = float(nd_data[3]), 0.0, 0.0
                        if nd_data[4] is not None:
                            y = float(nd_data[4])
                        if nd_data[5] is not None:
                            z = float(nd_data[5])

                        self.femdb.AddNode(Node(n_id, x, y, z))
                        self.femdb.node_hash[n_id] = node_index
                        node_index += 1
                        self.iter_line = cdb_f.readline()

                elif self.iter_line.startswith("EBLOCK,"):
                    self.ReadEBlock(cdb_f)

                elif self.iter_line.startswith("MPDATA,"):
                    self.ReadMaterial(cdb_f)

                elif self.iter_line.startswith("SECTYPE,"):
                    self.ReadSection(cdb_f)

                elif self.iter_line.startswith("ACEL,"):
                    # TODO: 处理加速度, 重力
                    self.iter_line = cdb_f.readline()

                elif self.iter_line.startswith("D,"):
                    # 一般会约束很多很多自由度, 所以先暂时不跳出分支
                    nodes, directs, values = [], [], []
                    bd = AnsysBoundary()
                    while self.iter_line.startswith("D,"):
                        splits = self.iter_line.split(",")
                        nodes.append(int(splits[1]))
                        directs.append(splits[2].strip())
                        values.append(float(splits[3].strip()))
                        self.iter_line = cdb_f.readline()

                    bd.SetConstraintInfor(nodes, directs, values)
                    self.femdb.load_case.AddBoundary(bd)

                elif self.iter_line.startswith("F,"):
                    # 一般F也不止施加在一个自由度上, 所以用while暂不跳出分支
                    while self.iter_line.startswith("F,"):
                        splits = self.iter_line.split(",")
                        self.femdb.load_case.AddAnsysCLoad(
                            int(splits[1]), splits[2].strip(), float(splits[3].strip())
                        )
                        self.iter_line = cdb_f.readline()

                else:
                    # 文件底部
                    if not self.iter_line:
                        break
                    # 不支持的关键字或注释直接读取下一行, 不可以strip(), 否则空行会被认作程序结束
                    self.iter_line = cdb_f.readline()

        self.femdb.et_hash = self.et_hash
        self.femdb.SetGrpHash(self.ele_group_hash, self.ele_count)
        self.femdb.real_const_hash = self.real_constant_hash

    def ReadEBlock(self, f_handle):
        """
        读取单元信息
        The format of the element "block" is as follows for the SOLID format:
        - Field 1 - The material number.
        - Field 2 - The element type number.
        - Field 3 - The real constant number.
        - Field 4 - The section ID attribute (beam section) number.
        - Field 5 - The element coordinate system number.
        - Field 6 - The birth/death flag.
        - Field 7 - The solid model reference number.
        - Field 8 - The element shape flag.
        - Field 9 - The number of nodes defining this element if Solid_key = SOLID;
        otherwise, Field 9 = 0.
        - Field 10 - Not used.
        - Field 11 - The element number.
        - Fields 12-19 - The node numbers. The next line will have the additional node
         numbers if there are more than eight.
        """
        solid_type = True
        if not self.iter_line.__contains__("SOLID"):
            solid_type = False

        fortran_format = f_handle.readline().strip()  # skip format line
        f_reader = ff.FortranRecordReader(fortran_format)
        self.iter_line = f_handle.readline()

        if solid_type:
            # TODO: 根据单元类型判断是否为二维问题, 即AnalyseDimension的值, 类比ABAQUS解析
            # TODO: 高阶单元, 节点在第二行的情况未处理
            while self.iter_line.startswith(" ") and not self.iter_line.__contains__("-1"):
                """
                解析字段内容
                """
                e_data = f_reader.read(self.iter_line)
                mat_num = e_data[0]
                e_type = e_data[1]
                real_constant_num = e_data[2]
                sec_id = e_data[3]
                parsed_nodes_count = e_data[8]
                ele_num = e_data[10]

                """
                保存至数据库, ANSYS单元的每一行都指定了材料等信息, 与ABAQUS不同. 在最后文件解析完成后再
                PrepareCalculateAnsys中分配各个单元信息
                """

                # 节点编号是无符号32位的, 也就是节点最大4294967295
                node_ids = np.zeros(parsed_nodes_count, dtype=np.uint32)
                search_ids = np.zeros(parsed_nodes_count, dtype=np.uint32)
                for idx in range(parsed_nodes_count):
                    node_ids[idx] = e_data[idx + 11]
                    search_ids[idx] = self.femdb.node_hash[node_ids[idx]]

                ele_node_list = list(OrderedDict.fromkeys(search_ids))
                iter_ele, e_node_count = ElementFactory.CreateElement(e_type=self.et_hash[e_type], opt=len(ele_node_list))
                iter_ele.e_type = e_type
                iter_ele.SetNodeSearchIndex(np.asarray(ele_node_list))

                # 除了组成单元所需的节点以外都是辅助节点, 默认e_node_count至parsed_nodes_count外的都是辅助节点
                for idx in range(e_node_count, parsed_nodes_count):
                    self.femdb.node_list[self.femdb.node_hash[node_ids[idx]]].is_assist_node = True

                iter_ele.SetId(ele_num)
                iter_ele.SetNodes(np.asarray(list(OrderedDict.fromkeys(node_ids))))
                iter_ele.mat_id = mat_num
                iter_ele.sec_id = sec_id
                iter_ele.real_const_id = real_constant_num

                # 计算单元包括的节点的坐标矩阵
                coords = []
                for nid in ele_node_list:
                    n_coord = self.femdb.node_list[nid].GetNodeCoord()
                    coords.append(n_coord)
                iter_ele.SetNodeCoords(np.asarray(coords))

                # 保存对应关系
                if self.ele_group_hash.__contains__(e_type):
                    ele_group = self.ele_group_hash[e_type]
                    self.femdb.ele_idx_hash[ele_num] = ele_group.GetElementsCurrentCount()
                    ele_group.AppendElement(copy.deepcopy(iter_ele))
                else:
                    new_ele_group = ElementGroup(e_type)
                    new_ele_group.AppendElement(copy.deepcopy(iter_ele))
                    self.femdb.ele_idx_hash[ele_num] = 0  # 第一个元素, 所以index为0
                    self.ele_group_hash[e_type] = new_ele_group

                # 更新并读取下一行
                self.ele_count += 1
                self.iter_line = f_handle.readline()
        else:
            """
            From the ANSYS programmer's manual
            The format without the SOLID keyword is:
            - Field 1 - The element number.
            - Field 2 - The type of section ID.
            - Field 3 - The real constant number.
            - Field 4 - The material number.
            - Field 5 - The element coordinate system number.
            - Fields 6-15 - The node numbers. The next line will have the
            additional node numbers if there are more than ten.
            """
            mlogger.fatal("UnSupport without SOLID keyword")
            sys.exit(1)

    def ReadMaterial(self, f_handle):
        """
        读取材料信息
        """
        # 如果读取材料参数结束, 那么jump_out为True, 跳出读取本elif分支读取其他信息
        jump_out = False
        while not jump_out:
            still_same_mat = True
            if self.iter_line.startswith("MPTEMP"):
                self.iter_line = f_handle.readline()
            cur_mat_id = int(self.iter_line.split(",")[4].strip())
            value_dict = {}
            while still_same_mat:
                # 暂时只读取"MPDATA"关键字
                if self.iter_line.startswith("MPTEMP,"):
                    self.iter_line = f_handle.readline()
                elif self.iter_line.startswith("MPDATA,"):
                    # 首先判断是否为同一种材料属性
                    splits = self.iter_line.split(",")
                    iter_mat_id = int(splits[4])
                    if iter_mat_id != cur_mat_id:
                        # 读取同一种材料结束, 读取下一种材料或者读取材料结束, 程序跳出材料分支
                        self.femdb.materials.append(ISOMaterial(cur_mat_id, value_dict))
                        self.iter_line = f_handle.readline()
                        break
                    if splits[3].startswith("EX"):
                        value_dict[MaterialKey.E] = float(splits[6])
                    elif splits[3].startswith("DENS"):
                        value_dict[MaterialKey.Density] = float(splits[6])
                    elif splits[3].startswith("NUXY"):
                        value_dict[MaterialKey.Niu] = float(splits[6])
                    self.iter_line = f_handle.readline()
                else:
                    # 当前行为其他信息, 跳出读材料分支, 读取其他
                    jump_out = not (self.iter_line.startswith("MPDATA,") or self.iter_line.startswith("MPTEMP"))
                    if jump_out:
                        self.femdb.materials.append(ISOMaterial(cur_mat_id, value_dict))
                        break

    def ReadSection(self, f_handle):
        """
        读取截面信息
        """
        while True:
            splits = self.iter_line.strip().split(",")
            sec_id = int(splits[1])
            if splits[2] == "BEAM":
                beam_type = splits[3]
                msec_data = f_handle.readline().strip().split(",")
                sec_data = [float(msec_data[idx]) for idx in range(1, len(msec_data)) if msec_data[idx]]
                f_handle.readline()  # section offset
                f_handle.readline()  # section control
                if beam_type == "RECT":
                    self.femdb.sections.append(BeamSection(sec_id, BeamSectionType.Rectangle, sec_data))
                else:
                    mlogger.fatal("UnSupport Beam Type:{}".format(beam_type))
                    sys.exit(1)

                self.iter_line = f_handle.readline()

            elif splits[2] == "SHELL":
                f_handle.readline()  # secoffset
                f_handle.readline()  # sec block
                f_handle.readline()  # sec block data
                f_handle.readline()  # sec control

            else:
                mlogger.fatal("UnSupport Section Type")
                sys.exit(1)

            # 当前行为其他信息, 跳出读材料分支, 读取其他
            self.iter_line = f_handle.readline()
            if not self.iter_line.startswith("SECTYPE,"):
                break

    def CheckModel(self):
        """
        检查模型是否有问题
        """
        node_ids = set(nd.id for nd in self.femdb.node_list)
        ele_nodes = []
        for _, group in self.femdb.ele_grp_hash.items():
            eles = group.Elements()
            for ele in eles:
                ele_nodes.extend(ele.node_ids)
        ele_nodes_set = set(ele_nodes)

        unnessary_nds = node_ids - ele_nodes_set
        if len(unnessary_nds) != 0:
            mlogger.debug("Unnecessary Node Ids:{}".format(unnessary_nds))
        else:
            mlogger.debug("No Unnecessary Node")

        mlogger.debug("Finish Check Model")


if __name__ == "__main__":
    list1 = [109, 169, 143, 143, 130, 130, 130, 130]
    print(np.asarray(list(OrderedDict.fromkeys(list1))))
