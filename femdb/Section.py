#!/usr/bin/env python3
# -*- coding: utf-8 -*-


class Section(object):
    """
    保存截面信息, 包括梁的还有壳的
    """

    def __init__(self, sec_name:str, sec_type:str, sec_data:list[float]):
        self.name = sec_name
        self.sec_type = sec_type
        self.sec_data = sec_data
        self.characters = None

    def GetName(self):
        return self.name

    def SetSectionCharacter(self, characters):
        """
        设置截面属性的相关值, t,s代表自然坐标系方向
        是字典属性, 梁包括内容：It, Is, Tor, A, At, As
        """
        self.characters = characters

    def GetSectionCharacter(self):
        return self.characters