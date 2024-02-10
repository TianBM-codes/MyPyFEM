#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from utils.GlobalEnum import *


class BeamSection(object):
    """
    保存梁的截面信息
    """

    def __init__(self, sec_name: int, sec_type: BeamSectionType, sec_data: list[float]):
        self.name = sec_name
        self.sec_type = sec_type
        self.sec_data = sec_data
        self.characters = None
        self.IsBeamSection = True

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
