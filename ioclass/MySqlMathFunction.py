# -*- coding: utf-8 -*-
import copy
import datetime
import logging
import sys
import time
import os

# import _mysql_connector
import mysql.connector
import numpy as np
from dateutil.relativedelta import relativedelta
from scipy.spatial import distance_matrix

qi_zhong_ji = {"host": "192.168.3.8", "port": 30053, "database": "nbport_qzj_db",
               "username": "nbport_user", "password": "P@jtMUS3Yx"}


class MySQLMathFunction(object):
    """
    提供方便的MySQL查询函数，返回搜索值，搜索结果的样式为[(time1,value1), (time2,value2)...]
    Reference:
    1. https://www.runoob.com/python3/python-mysql-connector.html
    """

    def __init__(self):
        db_info = qi_zhong_ji
        self.con = mysql.connector.connect(host=db_info["host"],
                                           port=db_info["port"],
                                           database=db_info["database"],
                                           user=db_info["username"],
                                           password=db_info["password"]
                                           )

        self.cursor = self.con.cursor()

    def isdb_connected(self):
        return self.con is not None

    def execute_sql_batch(self, query, data):
        """
        执行批量插入操作
        """
        self.cursor = self.con.cursor()
        self.cursor.executemany(query, data)
        self.con.commit()

    def close(self):
        """
        释放连接线程
        """
        self.cursor.close()
        self.con.close()

    def execute_sql(self, sql):
        """
        直接执行sql, 只支持查询
        """
        self.cursor = self.con.cursor()
        try:
            self.cursor.execute(sql)
        except BaseException as e:
            self.cursor.close()
            sys.exit(7)
        sql_data = self.cursor.fetchall()
        self.cursor.close()

        return sql_data

    def commit_sql(self, sql):
        """
        提交sql
        """
        self.cursor = self.con.cursor()
        try:
            self.cursor.execute(sql)
            self.con.commit()
            row_count = "{} record(s) affected".format(self.cursor.rowcount)
        except BaseException as e:
            self.cursor.close()
            sys.exit(7)
        self.cursor.close()

        return row_count


if __name__ == "__main__":
    print(int(time.time()))

    sql_db = MySQLMathFunction()
    dat_path_p = "/mnt/proj/NBQZJ/tbm/zz.dat"
    sql_ = (f"INSERT INTO t_calculate_nephogram (calculate_time, file_name, result_type) "
            f"VALUES (NOW(), '{dat_path_p}', 'rt') "
            "ON DUPLICATE KEY UPDATE "
            "calculate_time = VALUES(calculate_time), "
            "file_name = VALUES(file_name), "
            "result_type = VALUES(result_type);"
            )
    # sql_db.commit_sql(sql_)

    D = 12
    H = 14
    dis_mag = [1, 23.1234654, 8, 20]
    linear_mises = [2, 8.89765, 1, 0]
    # sql = (f"INSERT INTO t_work_status (T, value1, value2, value3, value4) "
    #        f"VALUES (NOW(), '{D}', '13', '{np.max(np.array(dis_mag))}', '{np.max(np.array(linear_mises))}') "
    #        "ON DUPLICATE KEY UPDATE "
    #        "value1= VALUES(value1), "
    #        "value2= VALUES(value2), "
    #        "value3= VALUES(value3), "
    #        "value4= VALUES(value4);"
    #        )
    # sql = (f"INSERT INTO t_work_status (T, value1, value2, value3, value4) "
    #        f"VALUES (NOW(), '{D / 1000}', '13', '{np.max(np.array(dis_mag)):.3f}', '{np.max(np.array(linear_mises)):.3f}') "
    #        "ON DUPLICATE KEY UPDATE "
    #        "value1= VALUES(value1), "
    #        "value2= VALUES(value2), "
    #        "value3= VALUES(value3), "
    #        "value4= VALUES(value4);"
    #        )
    sql = (f"INSERT INTO t_work_status (T, value1, value2, value3, value4, value5) "
           f"VALUES (NOW(), '{D / 1000:.2f}', '13', '{np.max(np.array(dis_mag)):.3f}', '{np.max(np.array(linear_mises) / 1000000):.3f}', {H / 1000:.2f}) "
           "ON DUPLICATE KEY UPDATE "
           "value1= VALUES(value1), "
           "value2= VALUES(value2), "
           "value3= VALUES(value3), "
           "value4= VALUES(value4), "
           "value5= VALUES(value5);"
           )
    print(sql)
    sql_db.commit_sql(sql)
