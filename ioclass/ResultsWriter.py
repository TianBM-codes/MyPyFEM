#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import meshio
from femdb.FEMDataBase import *


class ResultsWriter(object):
    """
    计算结果导出类, 现支持UNV导出以及VTP导出
    """

    def __init__(self):
        self.femdb = FEMDataBase()

    def WriteVTPFile(self, path):
        """
        将结果写至vtp文件
        Reference:
        1. https://github.com/nschloe/meshio
        """
        # 模型部分
        coords = np.asarray([node.coord for node in self.femdb.node_list])
        all_eles = {}
        for key, ele_grp in self.femdb.ele_grp_hash.items():
            eles = ele_grp.Elements()
            ele_count = len(eles)
            ele_node_count = eles[0].nodes_count
            relations = np.zeros([ele_count, ele_node_count], dtype=int)
            for i in range(ele_count):
                relations[i] = eles[i].GetNodeSearchIndex()[:ele_node_count]

            suffix = GlobalInfor[GlobalVariant.InputFileSuffix]
            if suffix == InputFileType.CDB:
                key = self.femdb.et_hash[key]
                all_eles[Ansys2VTKType[key]] = relations
            elif suffix == InputFileType.INP:
                all_eles[Abaqus2VTKType[key]] = relations

        # 位移结果
        dis_value = np.asarray([node.displacement for node in self.femdb.node_list])
        displacement = {"displacement": dis_value}
        meshio.write_points_cells(
            filename=path,
            points=coords,
            cells=all_eles,
            point_data=displacement,
            # cell_data=cell_data,
            # field_data=field_data
        )

    def WriteUNVFile(self, u_path):
        """
        将结果写入UNV文件用SiPESC平台查看
        :param u_path: unv文件路径
        :return:
        """
        with open(u_path, 'w') as uf:
            """
            写入标题开头
            """
            uf.write('{ Header\n'
                     '( " ", 3.0, 1)\n'
                     '}\n'
                     '{ Node\n')

            """
            写入节点信息
            """
            uf.write(f' ({len(self.femdb.node_list)})\n')
            for nd in self.femdb.node_list:
                uf.write("( {}, {}, {}, {}, 0)\n".format(nd.id, nd.coord[0], nd.coord[1], nd.coord[2]))

            """
            写入单元信息
            """
            uf.write('}\n'
                     '{ Element\n')
            uf.write(f'( {len(self.femdb.ele_count)})\n')
            for _, group in self.femdb.ele_grp_hash.items:
                for ele in group:
                    node_str = ""
                    for nd in ele.node_ids:
                        node_str = node_str + str(nd.id) + ", "
                    node_str = node_str[:-1]
                    ele_line = "({},{}, 1, 0, 0, {})\n".format(ele.id, ele.unv_code, node_str)
                    uf.write(ele_line)

            """
            写入位移结果
            """
            uf.write('}\n'
                     '{ ResultSet\n'
                     '( name="Static Displacement", target="NodeResult", type="Vector", varlabels="X|Y|Z", subcase=0)\n')

            # 节点只保存位移幅值, 所以这里都保存在x方向的位移上
            for ii in range(len(self.femdb.node_list)):
                displacement = self.femdb.node_list[ii].displacement
                uf.write("( {}, {}, 0, 0)\n".format(nd.id, displacement))

            uf.write('}\n')
