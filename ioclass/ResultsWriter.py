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
import pathlib
import zlib
import struct
from io import BytesIO
import meshio
import numpy as np

from utils.MeshCleaner import *
from element.MeshElementFactory import *
from femdb.FEMDataBase import *
from ioclass.MySqlMathFunction import MySQLMathFunction
from projects.qizhongji.QiZhongJiZiTai import MQ1330Wrapper


def compute_normal_vector(p1, p2, p3):
    """
    计算三点法向量并归一化
    @param p1:
    @param p2:
    @param p3:
    @return: 单位法向量（numpy数组）
    """
    v1 = p2 - p1  # 向量AB
    v2 = p3 - p1  # 向量AC
    normal = np.cross(v1, v2)
    norm = np.linalg.norm(normal)

    if np.isclose(norm, 0):  # 浮点数精度容错
        raise ValueError("三点共线，无法构成有效平面")

    return normal / norm


class ResultsWriter(object):
    """
    计算结果导出类, 现支持UNV导出以及VTP导出
    """

    def __init__(self, use_mysql=False):
        self.femdb = FEMDataBase()
        self.use_mysql = use_mysql
        if use_mysql:
            self.mysql_db = MySQLMathFunction()

    def WriteStaticAnalysisVTUFile(self, path):
        """
        将结果写至vtu文件
        Reference:
        1. https://github.com/nschloe/meshio
        """
        # 模型部分
        coords = np.asarray([node.coord for node in self.femdb.node_list])
        all_eles = {}

        for iter_ele in self.femdb.elements:
            iter_relation = iter_ele.GetNodeSearchIndex().tolist()
            ele_type = iter_ele.vtu_type
            if all_eles.__contains__(ele_type):
                all_eles[ele_type].append(iter_relation)
            else:
                all_eles[ele_type] = [iter_relation]

        for iter_ele_info in self.femdb.additional_elements:
            iter_relation = [iter_ele_info[0], iter_ele_info[1]]
            ele_type = iter_ele_info[2]
            if all_eles.__contains__(ele_type):
                all_eles[ele_type].append(iter_relation)
            else:
                all_eles[ele_type] = [iter_relation]

        # 位移结果
        dis_value = np.reshape(self.femdb.linear_u, (-1, ModelInfo.PER_NODE_DOF))[:, :3]
        # node_res = {"displacement": dis_value, "mises": self.femdb.linear_mises}
        node_res = {"displacement": dis_value}
        # node_res = {"displacement": dis_value, "mises": self.femdb.linear_mises,
        #             "sigma_xx": self.femdb.sigma_xx, "sigma_yy": self.femdb.sigma_yy, "sigma_zz": self.femdb.sigma_zz,
        #             "tau_xy": self.femdb.tau_xy, "tau_xz": self.femdb.tau_xz, "tau_yz": self.femdb.tau_yz}
        meshio.write_points_cells(
            filename=path,
            points=coords,
            cells=all_eles,
            point_data=node_res,
            # cell_data=cell_data,
            # field_data=field_data
        )

    def WriteMises2DatFile(self, dat_path, struct_id):
        """
        将Mises结果写入文件
        :param dat_path:
        :param struct_id:
        :return:
        """
        struct_id += 1001
        src = pathlib.Path(dat_path)
        buffer = BytesIO()
        buffer.write(struct.pack('i', len(self.femdb.linear_mises)))
        buffer.write(struct.pack(
            f'{len(self.femdb.linear_mises)}f',
            *self.femdb.linear_mises
        ))
        raw_data = buffer.getvalue()
        compressed = zlib.compress(raw_data, level=9)

        with open(src, 'wb') as f:
            f.write(compressed)

        if self.use_mysql:
            sql = (f"INSERT INTO t_calculate_nephogram (calculate_time, file_name, result_type, struct_id) "
                   f"VALUES (NOW(), '{dat_path}', 'mises', {struct_id}) "
                   "ON DUPLICATE KEY UPDATE "
                   "calculate_time = VALUES(calculate_time), "
                   "file_name = VALUES(file_name), "
                   "result_type = VALUES(result_type);"
                   )
            self.mysql_db.commit_sql(sql)

    def WriteStaticResult2DatFile2(self, dat_path, wrp, struct_id, load):
        """
        对模型进行重新排序
        :param dat_path:
        :param wrp:
        :param struct_id:
        :param load:
        :return:
        """
        struct_id += 1001
        buffer = BytesIO()
        coords = np.asarray([node.coord for node in self.femdb.node_list])
        model_data = {"plot_type": 4,
                      "boundary_box": [np.min(coords[:, 0]), np.max(coords[:, 0]),
                                       np.min(coords[:, 1]), np.max(coords[:, 1]),
                                       np.min(coords[:, 2]), np.max(coords[:, 2])],
                      "nFrames": 2,
                      "global_min": 0,
                      "global_max": 0,
                      "part_count": 1,
                      "rotation_info": [],
                      "rot_length": []}

        """
        处理单元信息
        """
        all_eles = []
        each_part_ele_count = []
        for key, value in self.femdb.element_set_name_hash.items():
            ele_iter_ids = self.femdb.element_sets[value]
            iter_part_ele_count = 0
            for ii in ele_iter_ids:
                if ii not in self.femdb.ele_hash:
                    continue
                idx = self.femdb.ele_hash[ii]
                ele = self.femdb.elements[idx]
                ele_type = ele.__class__.__name__
                ele_nodes = ele.search_node_ids
                iter_ele, _ = MeshElementFactory.CreateElement(ele_type, ii)
                if iter_ele is not None:
                    iter_ele.setFaces(ele_nodes)
                    all_eles.append(iter_ele)
                    iter_part_ele_count += 1
            info = wrp.get_rotate_info(key)
            model_data["rotation_info"].extend(info)
            model_data["rot_length"].append(len(info))
            each_part_ele_count.append(iter_part_ele_count)
        D, H = wrp.get_circle_info()

        displacement = np.reshape(self.femdb.linear_u, (-1, ModelInfo.PER_NODE_DOF))[:, :3]
        dis_mag = np.sqrt(displacement[:, 0] ** 2 + displacement[:, 1] ** 2 + displacement[:, 2] ** 2)
        all_nodes = []
        all_res = []
        map_npy = []
        displacement_mag = []
        mises = []

        cursor = 0
        cur_part_idx = 0
        tri_face_idx = [0]
        iter_part_face_count = 0
        for ii, ele in enumerate(all_eles):
            triangles = ele.getAllTriangles(only_surface=True)
            for tri in triangles:
                all_nodes.extend([coords[row].flatten() for row in tri])
                iter_part_face_count += 1
                map_npy.extend(tri)
                displacement_mag.extend([dis_mag[row] for row in tri])
                mises.extend([self.femdb.linear_mises[row] for row in tri])

            cursor = cursor + 1
            if cursor >= each_part_ele_count[cur_part_idx]:
                cur_part_idx += 1
                tri_face_idx.append(iter_part_face_count + tri_face_idx[-1])
                cursor = 0
                iter_part_face_count = 0

        model_data["node_count"] = len(all_nodes)
        model_data["scalar_length"] = len(all_res)
        model_data["result"] = all_res

        """
        0. 保存对应关系, 每个三角面片都要保存一遍结果, 也就是结果存在冗余
        """
        # npy_file_path = pathlib.Path(dat_path).with_suffix(".npy")
        # np.save(npy_file_path, np.array(map_npy))

        """
        1. 头部数据写入
        """
        buffer.write(struct.pack('i', model_data["plot_type"]))
        buffer.write(struct.pack('6f', *model_data["boundary_box"]))
        buffer.write(struct.pack('i', model_data["nFrames"]))

        """
        2. 变换信息写入
        """
        rot_info_length = len(np.array(model_data["rot_length"]))
        buffer.write(struct.pack('i', rot_info_length))
        buffer.write(struct.pack(f'{rot_info_length}i', *model_data["rot_length"]))

        rr = len(model_data["rotation_info"])
        buffer.write(struct.pack('i', rr))
        buffer.write(struct.pack(f'{rr}f', *model_data["rotation_info"]))

        buffer.write(struct.pack('i', len(tri_face_idx)))
        buffer.write(struct.pack(f'{len(tri_face_idx)}i', *tri_face_idx))

        """
        3. 部件循环写入
        """
        # 标量结果（多帧）
        buffer.write(struct.pack('i', len(displacement_mag)))
        buffer.write(struct.pack(
            f'{len(displacement_mag)}f',
            # *model_data[f"iter_result_{iter_frame}"]
            *displacement_mag
        ))
        buffer.write(struct.pack('i', len(mises)))
        buffer.write(struct.pack(
            f'{len(mises)}f',
            # *model_data[f"iter_result_{iter_frame}"]
            *mises
        ))

        """
        4. 起重机地下圈的显示
        """
        max_dis = np.max(np.array(dis_mag))
        max_mises = np.max(np.array(self.femdb.linear_mises))

        sql = "SELECT value7 FROM nbport_qzj_db.t_data_work_status_info order by T desc limit 1;"
        D = self.mysql_db.execute_sql(sql)[0][0]

        buffer.write(struct.pack('f', D))
        buffer.write(struct.pack('f', H / 1000))
        buffer.write(struct.pack('f', max_dis))
        buffer.write(struct.pack('f', max_mises))
        buffer.write(struct.pack('f', float(load)))

        # D_H_xyz = [D + 1900, 0, 0, D + 1900, H, 0]
        # buffer.write(struct.pack(
        #     '6f',
        #     # *model_data[f"iter_result_{iter_frame}"]
        #     *D_H_xyz
        # ))

        """
        5. 压缩并写入文件
        """
        raw_data = buffer.getvalue()
        compressed = zlib.compress(raw_data, level=9)

        with open(dat_path, 'wb') as f:
            f.write(compressed)

        """
        6. 写入韩昊兵的疲劳结果
        """
        with open("/root/TianPyFem/mypyfem/model/year_calculate_1.dat", 'rb') as fileData:
            rawData = fileData.read()
            uncompressed_data = zlib.decompress(rawData)
            buffer = BytesIO(uncompressed_data)
            _ = struct.unpack('i', buffer.read(4))[0]
            res_length = struct.unpack('i', buffer.read(4))[0]
            total_damage = struct.unpack(f'{res_length}f', buffer.read(res_length * 4))
        buffer2 = BytesIO()
        buffer2.write(struct.pack('i', model_data["plot_type"]))
        buffer2.write(struct.pack('6f', *model_data["boundary_box"]))
        buffer2.write(struct.pack('i', model_data["nFrames"]))

        rot_info_length = len(np.array(model_data["rot_length"]))
        buffer2.write(struct.pack('i', rot_info_length))
        buffer2.write(struct.pack(f'{rot_info_length}i', *model_data["rot_length"]))

        rr = len(model_data["rotation_info"])
        buffer2.write(struct.pack('i', rr))
        buffer2.write(struct.pack(f'{rr}f', *model_data["rotation_info"]))

        buffer2.write(struct.pack('i', len(tri_face_idx)))
        buffer2.write(struct.pack(f'{len(tri_face_idx)}i', *tri_face_idx))

        buffer2.write(struct.pack('i', len(total_damage)))
        buffer2.write(struct.pack(
            f'{len(total_damage)}f',
            # *model_data[f"iter_result_{iter_frame}"]
            *total_damage
        ))
        raw_data = buffer2.getvalue()
        compressed = zlib.compress(raw_data, level=9)

        fatigue_path = dat_path.with_stem(dat_path.stem + "_fatigue")
        with open(fatigue_path, 'wb') as f:
            f.write(compressed)

        """
        7. 如果涉及MySQL数据库, 将结果写入数据库
        """
        if self.use_mysql:
            sql = (f"INSERT INTO t_calculate_nephogram (calculate_time, file_name, result_type, struct_id) "
                   f"VALUES (NOW(), '{dat_path}', 'rt', {struct_id}) "
                   "ON DUPLICATE KEY UPDATE "
                   "calculate_time = VALUES(calculate_time), "
                   "file_name = VALUES(file_name), "
                   "result_type = VALUES(result_type);"
                   )
            self.mysql_db.commit_sql(sql)
            sql = (f"INSERT INTO t_calculate_nephogram (calculate_time, file_name, result_type, struct_id) "
                   f"VALUES (NOW(), '{fatigue_path}', 'fatigue', {struct_id}) "
                   "ON DUPLICATE KEY UPDATE "
                   "calculate_time = VALUES(calculate_time), "
                   "file_name = VALUES(file_name), "
                   "result_type = VALUES(result_type);"
                   )
            self.mysql_db.commit_sql(sql)
            sql = (f"INSERT INTO t_work_status (T, value1, value2, value3, value4, value5, struct_id) "
                   f"VALUES (NOW(), '{D / 1000:.2f}', '13', {max_dis}, {max_mises}, {H / 1000:.2f}, {struct_id}) "
                   "ON DUPLICATE KEY UPDATE "
                   "value1= VALUES(value1), "
                   "value2= VALUES(value2), "
                   "value3= VALUES(value3), "
                   "value4= VALUES(value4), "
                   "value5= VALUES(value5);"
                   )
            self.mysql_db.commit_sql(sql)

    def WriteResortModel2DatFile(self, dat_path):
        """
        对模型按单元集合进行重新排序
        :param dat_path:
        :return:
        """
        buffer = BytesIO()
        coords = np.asarray([node.coord for node in self.femdb.node_list])
        model_data = {"plot_type": 4,
                      "boundary_box": [np.min(coords[:, 0]), np.max(coords[:, 0]),
                                       np.min(coords[:, 1]), np.max(coords[:, 1]),
                                       np.min(coords[:, 2]), np.max(coords[:, 2])],
                      "nFrames": 0,
                      "var_name": "displacement",
                      "global_min": 0,
                      "global_max": 0,
                      "part_count": 1,
                      "rotation_info": [],
                      "rot_length": []}

        """
        处理单元信息
        """
        all_eles = []
        for key, value in self.femdb.element_set_name_hash.items():
            ele_iter_ids = self.femdb.element_sets[value]
            for ii in ele_iter_ids:
                if ii not in self.femdb.ele_hash:
                    continue
                idx = self.femdb.ele_hash[ii]
                ele = self.femdb.elements[idx]
                ele_type = ele.__class__.__name__
                ele_nodes = ele.search_node_ids
                iter_ele, _ = MeshElementFactory.CreateElement(ele_type, ii)
                if iter_ele is not None:
                    iter_ele.setFaces(ele_nodes)
                    all_eles.append(iter_ele)

        # MarkSurface(all_eles)
        all_nodes = []
        all_tri_faces = []
        all_res = []
        edges = []
        all_node_normal = []
        edge_count = 0
        tri_face_count = 0
        map_npy = []

        for ii, ele in enumerate(all_eles):
            triangles = ele.getAllTriangles(only_surface=True)
            for tri in triangles:
                all_nodes.extend([coords[row].flatten() for row in tri])
                all_node_normal.extend([compute_normal_vector(*[coords[row] for row in tri]).flatten()] * 3)
                begin_idx = len(all_tri_faces)
                all_tri_faces.extend([begin_idx, begin_idx + 1, begin_idx + 2])
                edges.extend([begin_idx, begin_idx + 1, begin_idx + 1, begin_idx + 2, begin_idx + 2, begin_idx])
                edge_count += 3
                tri_face_count += 3
                map_npy.extend(tri)

        model_data["node_count"] = len(all_nodes)
        model_data["node_coords"] = list(np.array(all_nodes).flatten())
        model_data["node_normal"] = list(np.array(all_node_normal).flatten())
        model_data["scalar_length"] = len(all_res)
        model_data["result"] = all_res
        model_data["tri_face_count"] = tri_face_count
        model_data["tri_faces"] = all_tri_faces
        model_data["edges"] = edges
        model_data["edge_count"] = edge_count
        model_data["outline_count"] = 0
        model_data["outline"] = []

        """
        0. 保存对应关系, 每个三角面片都要保存一遍结果, 也就是结果存在冗余
        """
        npy_file_path = pathlib.Path(dat_path).with_suffix(".npy")
        np.save(npy_file_path, np.array(map_npy))

        """
        1. 头部数据写入
        """
        buffer.write(struct.pack('i', model_data["plot_type"]))
        buffer.write(struct.pack('6f', *model_data["boundary_box"]))
        buffer.write(struct.pack('i', model_data["nFrames"]))

        """
        2. 变量名写入（固定256字节）
        """
        var_name = model_data["var_name"].encode('utf-8')
        buffer.write(var_name.ljust(256, b'\x00')[:256])  # 确保256字节

        """
        3. 极值写入
        """
        buffer.write(struct.pack('f', model_data["global_min"]))
        buffer.write(struct.pack('f', model_data["global_max"]))

        """
        5. 部件循环写入
        """
        buffer.write(struct.pack('i', model_data["part_count"]))
        for i_part in range(model_data["part_count"]):
            buffer.write(struct.pack('i', model_data["node_count"]))
        buffer.write(struct.pack(f'{model_data["node_count"] * 3}f', *model_data["node_coords"]))
        buffer.write(struct.pack(f'{model_data["node_count"] * 3}f', *model_data["node_normal"]))

        # 标量结果（多帧）
        buffer.write(struct.pack('i', model_data["scalar_length"]))
        for iter_frame in range(model_data["nFrames"]):
            buffer.write(struct.pack(
                f'{model_data["scalar_length"]}f',
                # *model_data[f"iter_result_{iter_frame}"]
                *model_data["result"]
            ))

        # 面片数据
        buffer.write(struct.pack('i', model_data["tri_face_count"]))
        buffer.write(struct.pack(
            f'{model_data["tri_face_count"]}i',
            *model_data["tri_faces"]
        ))

        # 轮廓线（原outline字段）
        buffer.write(struct.pack('i', model_data["outline_count"]))
        buffer.write(struct.pack(
            f'{model_data["outline_count"]}i',
            *model_data["outline"]
        ))

        # 单元边（原edges字段）
        buffer.write(struct.pack('i', model_data["edge_count"] * 2))
        buffer.write(struct.pack(
            f'{model_data["edge_count"] * 2}i',
            *model_data["edges"]
        ))

        """
        6. 压缩并写入文件
        """
        raw_data = buffer.getvalue()
        compressed = zlib.compress(raw_data, level=9)

        with open(dat_path, 'wb') as f:
            f.write(compressed)

    def WriteModel2DatFileWithoutRes(self, dat_path):
        """
        只写模型信息，不写入结果
        :param dat_path:
        :return:
        """
        buffer = BytesIO()
        coords = np.asarray([node.coord for node in self.femdb.node_list])
        model_data = {"plot_type": 4,
                      "boundary_box": [np.min(coords[:, 0]), np.max(coords[:, 0]),
                                       np.min(coords[:, 1]), np.max(coords[:, 1]),
                                       np.min(coords[:, 2]), np.max(coords[:, 2])],
                      "nFrames": 0, "var_name": "displacement"}
        displacement = np.reshape(self.femdb.linear_u, (-1, ModelInfo.PER_NODE_DOF))[:, :3]
        dis_mag = np.sqrt(displacement[:, 0] ** 2 + displacement[:, 1] ** 2 + displacement[:, 2] ** 2)
        model_data["global_min"] = np.min(dis_mag)
        model_data["global_max"] = np.max(dis_mag)
        model_data["part_count"] = 1

        """
        处理单元信息
        """
        all_eles = []
        for ii, ele in enumerate(self.femdb.elements):
            ele_type = ele.__class__.__name__
            ele_nodes = ele.search_node_ids
            iter_ele, _ = MeshElementFactory.CreateElement(ele_type, ii)
            if iter_ele is not None:
                iter_ele.setFaces(ele_nodes)
                all_eles.append(iter_ele)

        MarkSurface(all_eles)
        all_nodes = []
        all_tri_faces = []
        all_res = []
        edges = []
        all_node_normal = []
        edge_count = 0
        tri_face_count = 0
        map_npy = []
        for ii, ele in enumerate(all_eles):
            triangles = ele.getAllTriangles(only_surface=True)
            for tri in triangles:
                all_nodes.extend([coords[row].flatten() for row in tri])
                all_node_normal.extend([compute_normal_vector(*[coords[row] for row in tri]).flatten()] * 3)
                begin_idx = len(all_tri_faces)
                all_tri_faces.extend([begin_idx, begin_idx + 1, begin_idx + 2])
                edges.extend([begin_idx, begin_idx + 1, begin_idx + 1, begin_idx + 2, begin_idx + 2, begin_idx])
                edge_count += 3
                tri_face_count += 3
                map_npy.extend(tri)

        model_data["node_count"] = len(all_nodes)
        model_data["node_coords"] = np.array(all_nodes).flatten()
        model_data["node_normal"] = np.array(all_node_normal).flatten()
        model_data["scalar_length"] = len(all_res)
        model_data["result"] = all_res
        model_data["tri_face_count"] = tri_face_count
        model_data["tri_faces"] = all_tri_faces
        model_data["edges"] = edges
        model_data["edge_count"] = edge_count
        model_data["outline_count"] = 0
        model_data["outline"] = []

        """
        0. 保存对应关系, 每个三角面片都要保存一遍结果, 也就是结果存在冗余
        """
        npy_file_path = pathlib.Path(dat_path).with_suffix(".npy")
        np.save(npy_file_path, np.array(map_npy))

        """
        1. 头部数据写入
        """
        buffer.write(struct.pack('i', model_data["plot_type"]))
        buffer.write(struct.pack('6f', *model_data["boundary_box"]))
        buffer.write(struct.pack('i', model_data["nFrames"]))

        """
        2. 变量名写入（固定256字节）
        """
        var_name = model_data["var_name"].encode('utf-8')
        buffer.write(var_name.ljust(256, b'\x00')[:256])  # 确保256字节

        """
        3. 极值写入
        """
        buffer.write(struct.pack('f', model_data["global_min"]))
        buffer.write(struct.pack('f', model_data["global_max"]))

        """
        4. 部件循环写入
        """
        buffer.write(struct.pack('i', model_data["part_count"]))
        for i_part in range(model_data["part_count"]):
            # 顶点数据
            buffer.write(struct.pack('i', model_data["node_count"]))
            buffer.write(struct.pack(f'{model_data["node_count"] * 3}f', *model_data["node_coords"]))
            buffer.write(struct.pack(f'{model_data["node_count"] * 3}f', *model_data["node_normal"]))

            # 标量结果（多帧）
            buffer.write(struct.pack('i', model_data["scalar_length"]))
            for iter_frame in range(model_data["nFrames"]):
                buffer.write(struct.pack(
                    f'{model_data["scalar_length"]}f',
                    # *model_data[f"iter_result_{iter_frame}"]
                    *model_data["result"]
                ))

            # 面片数据
            buffer.write(struct.pack('i', model_data["tri_face_count"]))
            buffer.write(struct.pack(
                f'{model_data["tri_face_count"]}i',
                *model_data["tri_faces"]
            ))

            # 轮廓线（原outline字段）
            buffer.write(struct.pack('i', model_data["outline_count"]))
            buffer.write(struct.pack(
                f'{model_data["outline_count"]}i',
                *model_data["outline"]
            ))

            # 单元边（原edges字段）
            buffer.write(struct.pack('i', model_data["edge_count"] * 2))
            buffer.write(struct.pack(
                f'{model_data["edge_count"] * 2}i',
                *model_data["edges"]
            ))

        """
        5. 压缩并写入文件
        """
        raw_data = buffer.getvalue()
        compressed = zlib.compress(raw_data, level=9)

        with open(dat_path, 'wb') as f:
            f.write(compressed)

    def WriteStaticResult2DatFile(self, dat_path):
        """
        将静态结果写入到Dat文件中
        :param dat_path:
        :return:
        """
        buffer = BytesIO()
        coords = np.asarray([node.coord for node in self.femdb.node_list])
        model_data = {"plot_type": 4,
                      "boundary_box": [np.min(coords[:, 0]), np.max(coords[:, 0]),
                                       np.min(coords[:, 1]), np.max(coords[:, 1]),
                                       np.min(coords[:, 2]), np.max(coords[:, 2])],
                      "nFrames": 2, "var_name": "displacement"}
        displacement = np.reshape(self.femdb.linear_u, (-1, ModelInfo.PER_NODE_DOF))[:, :3]
        dis_mag = np.sqrt(displacement[:, 0] ** 2 + displacement[:, 1] ** 2 + displacement[:, 2] ** 2)
        model_data["global_min"] = np.min(dis_mag)
        model_data["global_max"] = np.max(dis_mag)
        model_data["part_count"] = 1

        """
        处理单元信息
        """
        all_eles = []
        for ii, ele in enumerate(self.femdb.elements):
            ele_type = ele.__class__.__name__
            ele_nodes = ele.search_node_ids
            iter_ele, _ = MeshElementFactory.CreateElement(ele_type, ii)
            if iter_ele is not None:
                iter_ele.setFaces(ele_nodes)
                all_eles.append(iter_ele)

        MarkSurface(all_eles)
        all_nodes = []
        all_tri_faces = []
        displacement_mag = []
        mises = []
        edges = []
        all_node_normal = []
        edge_count = 0
        tri_face_count = 0
        for ii, ele in enumerate(all_eles):
            triangles = ele.getAllTriangles(only_surface=True)
            for tri in triangles:
                all_nodes.extend([coords[row].flatten() for row in tri])
                all_node_normal.extend([compute_normal_vector(*[coords[row] for row in tri]).flatten()] * 3)
                begin_idx = len(all_tri_faces)
                all_tri_faces.extend([begin_idx, begin_idx + 1, begin_idx + 2])
                edges.extend([begin_idx, begin_idx + 1, begin_idx + 1, begin_idx + 2, begin_idx + 2, begin_idx])
                displacement_mag.extend([dis_mag[row] for row in tri])
                mises.extend([self.femdb.linear_mises[row] for row in tri])
                edge_count += 3
                tri_face_count += 3

        model_data["node_count"] = len(all_nodes)
        model_data["node_coords"] = np.array(all_nodes).flatten()
        model_data["node_normal"] = np.array(all_node_normal).flatten()
        model_data["scalar_length"] = len(displacement_mag)
        model_data["displacement_mag"] = displacement_mag
        model_data["mises"] = mises
        model_data["tri_face_count"] = tri_face_count
        model_data["tri_faces"] = all_tri_faces
        model_data["edges"] = edges
        model_data["edge_count"] = edge_count
        model_data["outline_count"] = 0
        model_data["outline"] = []

        """
        1. 头部数据写入
        """
        buffer.write(struct.pack('i', model_data["plot_type"]))
        buffer.write(struct.pack('6f', *model_data["boundary_box"]))
        buffer.write(struct.pack('i', model_data["nFrames"]))

        """
        2. 变量名写入（固定256字节）
        """
        var_name = model_data["var_name"].encode('utf-8')
        buffer.write(var_name.ljust(256, b'\x00')[:256])  # 确保256字节

        """
        3. 极值写入
        """
        buffer.write(struct.pack('f', model_data["global_min"]))
        buffer.write(struct.pack('f', model_data["global_max"]))

        """
        4. 部件循环写入
        """
        buffer.write(struct.pack('i', model_data["part_count"]))
        for i_part in range(model_data["part_count"]):
            # 顶点数据
            buffer.write(struct.pack('i', model_data["node_count"]))
            buffer.write(struct.pack(f'{model_data["node_count"] * 3}f', *model_data["node_coords"]))
            buffer.write(struct.pack(f'{model_data["node_count"] * 3}f', *model_data["node_normal"]))

            # 标量结果（多帧）
            buffer.write(struct.pack('i', model_data["scalar_length"]))
            buffer.write(struct.pack(
                f'{model_data["scalar_length"]}f',
                # *model_data[f"iter_result_{iter_frame}"]
                *model_data["displacement_mag"]
            ))
            buffer.write(struct.pack(
                f'{model_data["scalar_length"]}f',
                # *model_data[f"iter_result_{iter_frame}"]
                *model_data["mises"]
            ))

            # 面片数据
            buffer.write(struct.pack('i', model_data["tri_face_count"]))
            buffer.write(struct.pack(
                f'{model_data["tri_face_count"]}i',
                *model_data["tri_faces"]
            ))

            # 轮廓线（原outline字段）
            buffer.write(struct.pack('i', model_data["outline_count"]))
            buffer.write(struct.pack(
                f'{model_data["outline_count"]}i',
                *model_data["outline"]
            ))

            # 单元边（原edges字段）
            buffer.write(struct.pack('i', model_data["edge_count"] * 2))
            buffer.write(struct.pack(
                f'{model_data["edge_count"] * 2}i',
                *model_data["edges"]
            ))

        """
        5. 压缩并写入文件
        """
        raw_data = buffer.getvalue()
        compressed = zlib.compress(raw_data, level=9)

        with open(dat_path, 'wb') as f:
            f.write(compressed)

    def WriteSeriesResult(self, directory, pro_name):
        """
        将时序结果写入文件
        :param directory: 存储路径
        :param pro_name: 项目名称
        :return:
        """
        coords = np.asarray([node.coord for node in self.femdb.node_list])
        all_eles = {}

        for iter_ele in self.femdb.elements:
            iter_relation = iter_ele.GetNodeSearchIndex().tolist()
            ele_type = iter_ele.vtu_type
            if all_eles.__contains__(ele_type):
                all_eles[ele_type].append(iter_relation)
            else:
                all_eles[ele_type] = [iter_relation]

        for ii in range(self.femdb.history_step_count):
            path = directory + "/" + pro_name + f"_{ii}.vtu"
            u = self.femdb.history_u[ii]
            v = self.femdb.history_v[ii]
            a = self.femdb.history_a[ii]
            s = self.femdb.history_s[ii]
            time_result = {"displacement": np.reshape(u, (-1, ModelInfo.PER_NODE_DOF))[:, :3],
                           "velocity": np.reshape(v, (-1, ModelInfo.PER_NODE_DOF))[:, :3],
                           "acceleration": np.reshape(a, (-1, ModelInfo.PER_NODE_DOF))[:, :3],
                           "mises": s
                           }
            meshio.write_points_cells(
                filename=path,
                points=coords,
                cells=all_eles,
                point_data=time_result,
            )

    def WriteFatigueResult(self, path):
        """
        将疲劳结果谢至vtu文件
        :param path:
        :return:
        """
        coords = np.asarray([node.coord for node in self.femdb.node_list])
        all_eles = {}

        for iter_ele in self.femdb.elements:
            iter_relation = iter_ele.GetNodeSearchIndex().tolist()
            ele_type = iter_ele.vtu_type
            if all_eles.__contains__(ele_type):
                all_eles[ele_type].append(iter_relation)
            else:
                all_eles[ele_type] = [iter_relation]

        for iter_ele_info in self.femdb.additional_elements:
            iter_relation = [iter_ele_info[0], iter_ele_info[1]]
            ele_type = iter_ele_info[2]
            if all_eles[ele_type].__contains__(ele_type):
                all_eles[ele_type].append(iter_relation)
            else:
                all_eles[ele_type] = [iter_relation]

        # 疲劳结果
        node_res = {"fatigue": self.femdb.damage_factor}
        meshio.write_points_cells(
            filename=path,
            points=coords,
            cells=all_eles,
            point_data=node_res
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
            uf.write(f'(  {len(self.femdb.elements)};)\n')
            eles = self.femdb.elements
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


if __name__ == "__main__":
    file_path = "D:/WorkSpace/FEM/MyPyFEM/NumericalCases/Projects/qizhongji/last/MQ1330_remesh.cdb"
    file_lib = pathlib.Path(file_path)
    new_name = file_lib.stem + "_mises"
    new_path = file_lib.with_stem(new_name)
    print(new_path)
    res = ResultsWriter(use_mysql=True)

