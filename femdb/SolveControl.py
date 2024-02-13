#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from utils.CustomException import *

class ARCLEN:
    def __init__(self):
        self.arcln = None
        self.farcl = None
        self.xincr = None
        self.afail = None
        self.itarget = None
        self.iter_old = None

class Output:
    def __init__(self):
        self.inc_out = None
        self.nwant = None
        self.iwant = None

class SolveControl(object):
    def __init__(self):
        self.nincr = None
        self.xlmax = None
        self.dlamb = None
        self.miter = None
        self.cnorm = None
        self.searc = None
        self.msearch = None
        self.incrm = 0
        self.xlamb = None
        self.niter = 0
        self.Arclen = ARCLEN()
        self.Output = Output()



class FlagSHyPSolveControl(SolveControl, object):
    """
    Control Information, BOOK P264 P267 Note5
    0.  Number of load/displacement increments
    1.  Maximum value of load-scaling parameter
    2.  Load parameter increment
    3.  Maximum number of iterations per increment
    4.  Convergence tolerance
    5.  Line Search parameter(if 0.0 not in use)
    6.  Arc Length parameter (if 0.0 not in use)
    7.  Output counter (e.g. 5 for every five increments)
    8.  Target iterations per increment
    9.  Single output node (0 if not used)
    10. Output degree of freedom (0 if not used)

    If arc length control is employed then lambda is controlled indirectly;
    A positive value of CON.ARCLEN.arcln will produce a variable arc length,
    the value of which is determined by the desired number of iterations
    per increment itarget. Irrespective of the magnitude input by the user
    (i.e. CON.ARCLEN.arcln= 1), the program will work out the most appropriate
    value, discarding that entered by the user. A negative value (simply
    as an indicator) for CON.ARCLEN.arcln will provide a constant arc length.
    In this latter case, some experimentation with values will be necessary.
    If the arc length option is not to be used, input CON.ARCLEN.arcln = 0
    and CON.ARCLEN.itarget = 0.
    """

    def __init__(self):
        super().__init__()

    def InitWithTextLine(self, line):
        split_line = line.split()
        self.nincr = int(split_line[0])
        self.xlmax = int(split_line[1])
        self.dlamb = float(split_line[2])
        self.miter = int(split_line[3])
        self.cnorm = float(split_line[4])
        self.searc = float(split_line[5])
        self.Arclen.arcln = float(split_line[6])
        self.Output.inc_out = int(split_line[7])
        self.Arclen.itarget = float(split_line[8])
        self.Arclen.iter_old = self.Arclen.itarget
        self.Output.nwant = int(split_line[9])
        self.Output.iwant = int(split_line[10])
        if abs(self.searc * abs(self.Arclen.arcln)):
            raise OtherException("Error reading solution control parameters.\n"
                                 "Line search and arc length methods cannot be both activated.\n")
