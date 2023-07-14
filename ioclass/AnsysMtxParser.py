#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pylab as plt
import numpy as np


def GetStiffnessMatrixFromFullFile(file_path):
    """
    通过full文件获取ANSYS的刚度矩阵
    Reference:
    1. https://math.docs.pyansys.com/version/stable/examples/eigen_solve.html
    """

    import ansys.math.core.math as pymath

    # Start PyAnsys Math as a service.
    mm = pymath.AnsMath()

    # fullfile = mm._mapdl.jobname + ".full"
    fullfile  = "./fourshell0.full"
    k = mm.stiff(fname=fullfile)
    m = mm.mass(fname=fullfile)
    print(m.shape)
    print(k.shape)


GetStiffnessMatrixFromFullFile("")