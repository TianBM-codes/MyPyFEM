#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from utils.CustomException import *
from abc import ABC
from typing import Dict
from utils.GlobalEnum import MaterialKey


class MaterialBase(object):
    def __init__(self, name, value_dict):
        self.name = name
        self.value_dict = value_dict

    def GetName(self):
        return self.name

    def GetValueDict(self):
        return self.value_dict

    def SetValueDict(self, value_dict):
        self.value_dict = value_dict


class Material(object):
    def __init__(self):
        self.Mat: Dict[int, MaterialBase] = {}
        self.n_nearly_incompressible = None

    def InsertMaterial(self, mat_id, mat: MaterialBase):
        self.Mat[mat_id] = mat


class ISOMaterial(MaterialBase, ABC):
    """ 各向同性材料 """

    def __init__(self, name, value_dict=None):
        super().__init__(name, value_dict)

    def __str__(self):
        return "ISOMaterial: {},{}".format(self.name, self.value_dict)


class CompressibleNeoHookean(MaterialBase, ABC):
    """
    plane stain or three-dimensional compressible neo-Hookean
    """

    def __init__(self, name, value_dict):
        super().__init__(name, value_dict)


class HyperElasticPlasticInPrincipal(MaterialBase, ABC):
    """
    plane strain or three-dimensional nearly incompressible hyperelastic
    plastic in principle directions
    """

    def __init__(self, name, value_dict):
        super().__init__(name, value_dict)

    def InitByFlagSHyPFormat(self, line):
        lineSplit = line.split(" ")
        self.value_dict[MaterialKey.Density] = float(lineSplit[0])
        self.value_dict[MaterialKey.Niu] = float(lineSplit[1])
        self.value_dict[MaterialKey.Lamda] = float(lineSplit[2])
        self.value_dict[MaterialKey.Kappa] = (self.value_dict[MaterialKey.Lamda]
                                              + 2 / 3 * self.value_dict[MaterialKey.Niu])
        self.value_dict[MaterialKey.TauY] = float(lineSplit[3])
        self.value_dict[MaterialKey.Harden] = float(lineSplit[4])


class MaterialFactory(object):

    @staticmethod
    def CreateMaterial(mat_type: int):
        """

        @param mat_type:
        @return:
        """
        if mat_type in [1]:
            return CompressibleNeoHookean(None, None)
        elif mat_type in [17]:
            return HyperElasticPlasticInPrincipal(None, None)
        else:
            raise NoImplSuchMaterial(mat_type)
