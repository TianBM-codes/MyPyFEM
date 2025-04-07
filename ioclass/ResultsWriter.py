#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
    SiPESC平台单元类型对应的code
    xx-yy-zz
    xx：拓扑类型，节点数
    yy：单元分类，1梁，2杆，3膜，4板，5壳，6体，7轴对称壳，8轴对称四边形，9夹层板壳，10层合板壳，11柳钉，12连接
    zz：单元子类。
    TYPE=10000,    1点三向拉压弹簧单元;
    TYPE=11100,    1点铆钉单元;
    TYPE=20100,    2点偏心梁单元;
    TYPE=20200,    2点轴力杆单元;
    TYPE=20210,    2点弹塑性轴力杆单元;
    TYPE=21200，   2点平面接触单元;
    TYPE=21201，   2点空间接触单元;
    TYPE=30100,    3点弯管单元;
    TYPE=30300,    3点平面应力膜单元;
    TYPE=30310,    3点三角形弹塑性平面应力膜单元 ;
    TYPE=30400,    3点三角形薄板单元;
    TYPE=30500,    3点三角形薄壳单元;
    TYPE=30501,    3点三角形各向异性薄壳单元;
    TYPE=30900,    3点三角形复合材料夹层板壳单元;
    TYPE=31000,    3点三角形复合材料层合板壳单元
    TYPE=40300,    4点平面应力膜单元;
    TYPE=40301,    4点平面应变膜单元;
    TYPE=40302,    4点非协调平面应力膜单元;
    TYPE=40303,    4点非协调平面应变膜单元;
    TYPE=40304,    4点矩形膜单元;
    TYPE=40310,    4点弹塑性平面应力膜单元;
    TYPE=40500,    4点任意四边形壳单元;
    TYPE=40501,    4点任意四边形各向异性壳单元;
    TYPE=40700,    4点轴对称旋转壳单元;
    TYPE=40800,    4点任意四边形4节点轴对称环体单元;
    TYPE=40900,    4点任意四边形复合材料夹层板壳单元;
    TYPE=41000,    4点任意四边形复合材料层合板壳单元;
    TYPE=50300,    5点等参平面膜单元;
    TYPE=50600,    5点金字塔单元;
    TYPE=60300,    6点高阶三角形单元;
    TYPE=60600,    6点三棱柱单元;
    TYPE=80300,    8点高阶四边形单元;
    TYPE=80600,    8点块体单元;
    TYPE=80601,    8点非协调块体单元;
    TYPE=80610,    8点弹塑性块体单元;
    TYPE=100600,   10点高阶四面体单元;
    TYPE=130600,   13点高阶金字塔单元;
    TYPE=150600,   15点高阶三棱柱单元;
    TYPE=200600,   20点高阶六面体单元;
"""

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
            for i in range(ele_count):
                ele_node_count = eles[i].nodes_count
                iter_relation = eles[i].GetNodeSearchIndex()[:ele_node_count].tolist()
                key2 = str(self.femdb.et_hash[key]) + "_" + str(ele_node_count)
                if all_eles.__contains__(Ansys2VTKType[key2]):
                    all_eles[Ansys2VTKType[key2]].append(iter_relation)
                else:
                    all_eles[Ansys2VTKType[key2]] = [iter_relation]


        # 位移结果
        dis_value = np.asarray([node.dof_disp[:3] for node in self.femdb.node_list])
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
            uf.write('{ Header;\n'
                     '( "Model database", 2.0,1;)\n'
                     '}\n'
                     '{ Node;\n')

            """
            写入节点信息
            """
            uf.write(f'(   {len(self.femdb.node_list)};)\n')
            for nd in self.femdb.node_list:
                uf.write("( {}, {}, {}, {},   3;)\n".format(nd.id, nd.coord[0], nd.coord[1], nd.coord[2]))

            """
            写入单元信息
            """
            uf.write('}\n'
                     '{ Element;\n')
            uf.write(f'(  {self.femdb.ele_count};)\n')
            for _, group in self.femdb.ele_grp_hash.items():
                eles = group.Elements()
                for ele in eles:
                    node_str = ""
                    for nd in ele.node_ids:
                        node_str = node_str + str(nd) + ", "
                    node_str = node_str[:-2] + ";"
                    ele_line = "({},{}, 1, 0, 0, {})\n".format(ele.id, ele.unv_code, node_str)
                    uf.write(ele_line)

            """
            写入位移结果
            """
            uf.write('}\n'
                     '{ ResultSet\n'
                     '( name="Static Displacement", target="NodeResult", type="Vector", varlabels="X|Y|Z", subcase=0)\n')

            # 节点无需保存幅值, SiPESC平台自动计算
            for ii in range(len(self.femdb.node_list)):
                node = self.femdb.node_list[ii]
                node_id = node.id
                displacement = node.dof_disp[:3]
                uf.write("( {}, {:.6}, {:.6}, {:.6})\n".format(node_id, displacement[0], displacement[1], displacement[2]))

            # uf.write('}\n')

            """
            写入应力结果
            """
            # uf.write('{ StaticStrs;\n')
            # uf.write('( 1, {};)\n'.format(len(self.femdb.node_list)))
            # uf.write('{ StaticStrsSet;\n')
            # uf.write('( 1, "Stress", 2;)\n')
            # for ii in range(len(self.femdb.node_list)):
            #     node = self.femdb.node_list[ii]
            #     xx, yy, zz, xy, yz, xz = node.average_stress
            #     uf.write("( {}, {:.6}, {:.6}, {:.6}, {:.6}, {:.6},{:.6})\n".format(node.id, xx, yy, zz, xy, yz, xz))
            #
            # uf.write('}\n')  # 应力结果结束
            uf.write('}\n')  # 整个文件结束
