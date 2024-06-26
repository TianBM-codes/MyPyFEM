#!/usr/bin/env python3
# -*- coding: utf-8 -*-


class NoSupportDimension(Exception):
    def __init__(self, err_msg):
        self.err_msg = err_msg

    def __str__(self):
        return f"不支持的分析维度:\n{self.err_msg}"


class NoSupportLoadCase(Exception):
    def __init__(self, err_msg):
        self.err_msg = err_msg

    def __str__(self):
        return f"不支持的载荷工况:\n{self.err_msg}"


class InputTextFormat(Exception):
    def __init__(self, err_msg):
        self.err_msg = err_msg

    def __str__(self):
        return f"该行导入文件未能解析成功:\n{self.err_msg}"


class NoImplSuchElement(Exception):
    def __init__(self, err_msg):
        self.err_msg = err_msg

    def __str__(self):
        return f"未实现该类型的单元:\n{self.err_msg}"


class ElementNoUNV(Exception):
    def __init__(self, err_msg):
        self.err_msg = err_msg

    def __str__(self):
        return f"这类单元没有UNV对应格式:\n{self.err_msg}"


class NoImplSuchMaterialStress(Exception):
    def __init__(self, err_msg):
        self.err_msg = err_msg

    def __str__(self):
        return f"未实现该类型材料的应力计算:\n{self.err_msg}"


class NoSupportOption(Exception):
    def __init__(self, ele_type, opt):
        self.ele_type = ele_type
        self.opt = opt

    def __str__(self):
        return f"单元类型{self.ele_type}不支持{self.opt}选项"


class NoImplSuchMaterial(Exception):
    def __init__(self, err_msg):
        self.err_msg = err_msg

    def __str__(self):
        return f"未实现该材料类型:\n{self.err_msg}"


class WrongMaterialType(Exception):
    def __init__(self, err_msg):
        self.err_msg = err_msg

    def __str__(self):
        return f"错误的材料属性:\n{self.err_msg}"


class NoImplSuchShapeFunction(Exception):
    def __init__(self, err_msg):
        self.err_msg = err_msg

    def __str__(self):
        return f"未实现该单元形函数:\n{self.err_msg}"


class NoImplSuchElasticityModulus(Exception):
    def __init__(self, err_msg):
        self.err_msg = err_msg

    def __str__(self):
        return f"未实现该材料本构:\n{self.err_msg}"


class OtherException(Exception):
    def __init__(self, err_msg):
        self.err_msg = err_msg

    def __str__(self):
        return self.err_msg
