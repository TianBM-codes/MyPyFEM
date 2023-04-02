#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import copy

from femdb.FEMDataBase import *

"""
文件解析类是和Domain一样同样包含Element和数据库的类
"""


def ReadSectionLine(line):
    """
    多处用到的方法，提出来作为一个函数: 读取Inp文件的Section行, 返回一个字典
    """
    ret_diction = {}
    sps = line.split(",")
    for i in range(1, len(sps)):
        if "=" in sps[i]:
            keyval = sps[i].strip().split("=")
            ret_diction[keyval[0]] = keyval[1]
        else:
            mlogger.fatal("Fatal Error: in ReadSectionLine {}".format(line))
            sys.exit(1)
    return ret_diction


class InpReader(object):
    """
    功能为解析Inp文件, 初始化数据库单元节点及相关信息, 只负责解析和传递, 其他功能均在其他类中,
    同样的, 其他类也不会包含任何解析输入文件的功能或者函数
    """

    def __init__(self, input_path):
        self.fem_data = FEMDataBase()
        # mlogger.debug("In InpReader, femdb id is:{}".format(id(self.fem_data)))
        self.inp_path = input_path
        self.ele_count = 0
        self.ele_group_hash = {}  # key(单元名): value(ele_group), 存的实有值
        self.iter_line = ""

    def ParseFileAndInitFEMDB(self):
        """
        逐行读取数据，ReadPart ReadLoad等都是这样, 具体程序细节参见具体函数的注释
        由于对Inp文件的不熟悉, 现有如下假设或者限制：
        1. 只支持单个Part, 即忽略Assembly词条
        2. 所有相同类型的单元都在一起, 例如模型中所有的C3D8单元都在*Element, type=C3D8下, 即使被分为很多个ElSet, 每个ElSet中单元材料不同, 这样在
           Domain类中的ele_grp_hash的键就是唯一的, 带来的好处就是输出vtp文件时可以对ele_grp_hash进行循环.
        3. 不在读取Section的时候就赋予相关集合包含单元的材料属性等信息，因为有些材料可能还未定义, 等全部信息解析完后再赋予属性
        4. 默认Section定义的Eleset都在定义Section之前定义了, 默认定义的Section都对应Eleset, 即Eleset中的单元不为空, 这样就可以对属性进行循环来赋
          材料进而计算D阵

        Reference:
            1. http://wufengyun.com:888/books/usi/default.htm
            2. ABAQUS keyword browser table & Keyword support from the input file readers
            3.《ABAQUS有限元分析实例详解》P26  图2-30 数据库的结构示意图
        """
        with open(self.inp_path, 'r') as inp_f:
            self.iter_line = inp_f.readline()
            while True:
                # Part中的Nset和Eset一般都是该Part内部节点和单元的集合
                if self.iter_line.startswith("*Part,") or self.iter_line.startswith("*part") or self.iter_line.startswith("*PART"):
                    self.ReadPart(inp_f)

                elif self.iter_line.startswith("*Material,"):
                    self.ReadMaterial(inp_f)

                # 在Part外也会可有Nset和Elset, 比如设置约束或力的时候
                elif self.iter_line.startswith("*Nset,"):
                    self.ReadNset(inp_f)

                elif self.iter_line.startswith("*Elset,"):
                    self.ReadElset(inp_f)

                elif self.iter_line.startswith("*Step,"):
                    self.ReadLoadCase(inp_f)

                elif self.iter_line.startswith("*AbaqusBoundary"):
                    self.ReadBoundary(inp_f)

                else:
                    if not self.iter_line:
                        break
                    self.iter_line = inp_f.readline().strip()

        self.fem_data.SetGrpHash(self.ele_group_hash, self.ele_count)

    def ReadPart(self, f_handle):
        """
        读取part中的单元和节点, Part部分是以*End Part结束的, 先不读取Part中的*Elset, Part中包含节点、单元、节点集、单元集
        属性(Section)等
        :param f_handle: 文件句柄
        """
        self.iter_line = f_handle.readline().strip()
        while self.iter_line != "*End Part":
            if (self.iter_line == "*Node") or (self.iter_line == "*NODE"):
                node_index = 0  # Domain中的index
                self.iter_line = f_handle.readline().strip()
                while not self.iter_line.startswith("*"):
                    n_data = self.iter_line.split(",")
                    n_id = int(n_data[0])
                    x = float(n_data[1])
                    y = float(n_data[2])  # 最小是二维的, 即不会先让y=0.0
                    if len(n_data) == 4:
                        # 三维问题, 涉及的单元均为三维单元
                        z = float(n_data[3])
                        self.fem_data.AddNode(Node(n_id, x, y, z))
                    elif len(n_data) == 3:
                        # 二维问题, 涉及的单元均为二维单元
                        self.fem_data.AddNode(Node(n_id, x, y))

                    self.fem_data.node_hash[n_id] = node_index
                    node_index += 1
                    self.iter_line = f_handle.readline().strip()

            # 读取Elements的时候需要ElementGroups是因为输出vtp时不同单元不同样式, 默认每种Element都集中在一起,
            # 相同单元类型不同属性的话, 通过Section中的Set来区分
            elif self.iter_line.startswith("*Element,") or self.iter_line.startswith("*ELEMENT,"):
                # 解析单元类型关键字, 如果出现某些单元, 那么整个分析将变为2D分析
                e_type = self.iter_line.split(",")[1].split("=")[-1]
                mlogger.debug("Parsing Element Type: {}".format(e_type.strip()))
                if e_type in ["CPS3", "CPS4"]:
                    self.fem_data.an_dimension = AnalyseDimension.TwoDimension

                # 创建单元和单元组, ele_ids用来收集本组中单元的真实ID, nds是组成单个单元的真实节点号
                iter_ele, n_cnt = ElementFactory.CreateElement(e_type.strip())
                ele_group = ElementGroup(e_type)
                nds = np.zeros(n_cnt, dtype=np.uint32)
                ele_ids = []

                # 将*视为结束
                self.iter_line = f_handle.readline().strip()
                while not self.iter_line.startswith("*"):
                    sp_line = self.iter_line.split(",")
                    if sp_line[-1] == "":
                        sp_line.pop()
                    first_line_node_count = len(sp_line) - 1
                    iter_ele.SetId(int(sp_line[0]))
                    ele_ids.append(int(sp_line[0]))
                    sp_line[-1] = sp_line[-1].strip()  # 去掉\n换行符
                    for i in range(1, len(sp_line)):
                        nds[i - 1] = int(sp_line[i])

                    # 第一行的节点数不够, 第二行还有节点需要添加
                    if first_line_node_count < n_cnt:
                        self.iter_line = f_handle.readline().strip()
                        sp_line = self.iter_line.split(",")
                        for j in range(first_line_node_count, n_cnt):
                            nds[j] = int(sp_line[j - first_line_node_count])

                    # 读取数据完毕，首先设置单元包括的节点的搜索id
                    iter_ele.SetNodes(nds)
                    search_ids = np.array([], dtype=np.uint32)
                    for nd in nds:
                        search_ids = np.append(search_ids, self.fem_data.node_hash[nd])
                    iter_ele.SetNodeSearchIndex(search_ids)

                    # 计算单元包括的节点的坐标矩阵
                    coords = []
                    for nid in search_ids:
                        n_coord = self.fem_data.node_list[nid].GetNodeCoord()
                        # coords.append(np.stack((coords,n_coord)))
                        coords.append(n_coord)
                    coords = np.asarray(coords)
                    iter_ele.SetNodeCoords(coords)

                    # 保存对应关系
                    self.fem_data.ele_idx_hash[int(sp_line[0])] = ele_group.GetElementsCurrentCount()

                    # 需要进行深拷贝, 否则是一个单元重复了单元个数次
                    ele_group.AppendElement(copy.deepcopy(iter_ele))
                    self.ele_count += 1
                    self.iter_line = f_handle.readline().strip()

                ele_group.SetEleIdSet(set(ele_ids))
                self.ele_group_hash[e_type] = ele_group

            elif self.iter_line.startswith("*Nset,") or self.iter_line.startswith("*NSET,"):
                self.ReadNset(f_handle)

            elif self.iter_line.startswith("*Elset") or self.iter_line.startswith("ELSET"):
                self.ReadElset(f_handle)

            elif self.iter_line.startswith("*Solid Section"):
                """
                在当前程序解析属性的时候, 如果用到某个EleSet, 那么这个EleSet就是有用的
                """
                ret_dict = ReadSectionLine(self.iter_line)
                els_name = ret_dict["elset"]
                mat_name = ret_dict["material"]

                # Read Pars, 默认第一个参数为厚度参数, 当为杆的时候第一个参数为面积
                self.iter_line = f_handle.readline().strip()
                keywords = self.iter_line.split(",")
                pars = {}
                for par in keywords:
                    if par:
                        pars[PropertyKey.ThicknessOrArea] = float(par)
                section = Property(els_name, mat_name, pars)
                self.fem_data.properties.append(section)
                self.iter_line = f_handle.readline().strip()

            elif self.iter_line.startswith("*Beam Section"):
                """
                在当前程序解析属性的时候, 如果用到某个EleSet, 那么这个EleSet就是有用的
                Beam需要指定界面类型, 以及该类型的尺寸参数, 以及法线方向
                """
                ret_dict = ReadSectionLine(self.iter_line)
                els_name = ret_dict["elset"]
                mat_name = ret_dict["material"]
                section = ret_dict["section"]
                self.fem_data.GetSpecificFEMObject(FEMObject.EleSet, els_name).SetUsed(True)
                self.iter_line = f_handle.readline().strip()
                features = [float(iter_v) for iter_v in self.iter_line.split(",")]
                self.iter_line = f_handle.readline().strip()
                normal_dir = [float(di) for di in self.iter_line.split(",")]
                assert len(normal_dir) == 3
                self.fem_data.properties.append(Property(els_name, mat_name, [section, features, normal_dir]))
                self.iter_line = f_handle.readline().strip()

            elif self.iter_line.startswith("*Shell Section"):
                """
                在当前程序解析属性的时候, 如果用到某个EleSet, 那么这个EleSet就是有用的
                """
                # Shell
                ret_dict = ReadSectionLine(self.iter_line)
                els_name = ret_dict["elset"]
                mat_name = ret_dict["material"]
                self.fem_data.GetSpecificFEMObject(FEMObject.EleSet, els_name).SetUsed(True)
                self.iter_line = f_handle.readline().strip()
                features = [float(iter_v) for iter_v in self.iter_line.split(",")]
                self.iter_line = f_handle.readline().strip()
                self.fem_data.properties.append(Property(els_name, mat_name, features))
                self.iter_line = f_handle.readline().strip()

            else:
                # 对于暂不支持的内容直接读取下一行
                self.iter_line = f_handle.readline().strip()
                # 读取结束
                if not self.iter_line:
                    return
        self.iter_line = f_handle.readline().strip()

    def ReadMaterial(self, f_handle):
        """
        读取材料信息
        """
        mat_name = self.iter_line.split(",")[1].strip().split("=")[1]
        pars_dict = {}
        self.iter_line = f_handle.readline().strip()
        while not self.iter_line.startswith("**"):
            if self.iter_line == "*Density":
                self.iter_line = f_handle.readline().strip()
                pars_dict[MaterialKey.Density] = float(self.iter_line.split(",")[0])
                self.iter_line = f_handle.readline().strip()
            elif self.iter_line == "*Elastic":
                self.iter_line = f_handle.readline().strip()
                pars_dict[MaterialKey.E] = float(self.iter_line.split(",")[0])
                pars_dict[MaterialKey.Niu] = float(self.iter_line.split(",")[1])
                self.iter_line = f_handle.readline().strip()
            else:
                mlogger.fatal("Fatal Error: UnSupport Material Para Line {}".format(self.iter_line))
                sys.exit(1)

        self.fem_data.materials.append(ISOMaterial(mat_name, pars_dict))
        self.iter_line = f_handle.readline().strip()

    def ReadNset(self, f_handle):
        """
        读取节点集合
        """
        keywords = self.iter_line.strip().split(",")
        set_name = keywords[1].split("=")[1]
        is_generate = False
        for kw in keywords:
            if kw.strip() == "generate":
                is_generate = True
                break
        self.iter_line = f_handle.readline().strip()
        nodes = []

        # 有两种形式, 如果是generate, 那么是start, end, inc形式, 否则都按照罗列法
        if is_generate:
            begin_idx, end_idx, inc, = self.iter_line.split(",")
            for i in range(int(begin_idx), int(end_idx) + 1, int(inc)):
                nodes.append(i)
        else:
            while not self.iter_line.startswith("*"):
                for nd in self.iter_line.split(","):
                    if nd:
                        nodes.append(int(nd))
                self.iter_line = f_handle.readline().strip()
        self.fem_data.node_sets.append(NodeSet(set_name, nodes))

    def ReadElset(self, f_handle):
        """
        读取单元集合
        """
        keywords = self.iter_line.strip().split(",")
        set_name = keywords[1].split("=")[1]
        is_generate = False
        for kw in keywords:
            if kw.strip() == "generate":
                is_generate = True
                break
        self.iter_line = f_handle.readline().strip()

        # 有两种形式, 如果是generate, 那么是start, end, inc形式, 只占一行, 否则都按照罗列法
        eles = []
        if is_generate:
            begin_idx, end_idx, inc, = self.iter_line.split(",")
            for i in range(int(begin_idx), int(end_idx) + 1, int(inc)):
                eles.append(i)
            self.iter_line = f_handle.readline().strip()
        else:
            while not self.iter_line.startswith("*"):
                for ele_id in self.iter_line.split(","):
                    if ele_id:
                        eles.append(int(ele_id))
                self.iter_line = f_handle.readline().strip()
        self.fem_data.ele_sets.append(EleSet(set_name, eles))

    def ReadLoadCase(self, f_handle):
        """
        读取工况信息, 作为一次求解的信息, 包括约束、外力以及输出
        """
        self.iter_line = f_handle.readline().strip()
        while self.iter_line != "*End Step":
            if self.iter_line == "*AbaqusBoundary":
                self.ReadBoundary(f_handle)
            elif self.iter_line == "*Cload":
                self.iter_line = f_handle.readline().strip()
                keywords = self.iter_line.split(",")
                self.fem_data.load_case.AddAbaqusCLoad(keywords[0].strip(), int(keywords[1]), float(keywords[2]))
                self.iter_line = f_handle.readline().strip()
            else:
                # 其他情况先读取下一行, 直到遇到 *End Step为止
                self.iter_line = f_handle.readline().strip()
                if not self.iter_line:
                    return
        self.iter_line = f_handle.readline().strip()

    def ReadBoundary(self, f_handle):
        """
        读取边界条件
        """
        self.iter_line = f_handle.readline().strip()
        keywords = self.iter_line.split(",")
        if len(keywords) == 2:
            # 如果改行为一个逗号间隔, 那么约束方式是用字符串来标识的
            self.fem_data.load_case.AddBoundary(AbaqusBoundary(keywords[0].strip(), b_type=keywords[1].strip()))
        elif len(keywords) == 3:
            # 如果该行用两个逗号间隔, 那么是固定该自由度为零
            self.fem_data.load_case.AddBoundary(AbaqusBoundary(keywords[0].strip(), int(keywords[1]), 0.0))
        elif len(keywords) == 4:
            # 如果该行为三个逗号间隔, 那么为指定位移形式
            self.fem_data.load_case.AddBoundary(AbaqusBoundary(keywords[0].strip(), int(keywords[1]), float(keywords[3])))
        self.iter_line = f_handle.readline().strip()


if __name__ == "__main__":
    print(ReadSectionLine("*Beam Section, elset=_PickedSet8, material=Material-1, temperature=GRADIENTS, section=PIPE\n"))
