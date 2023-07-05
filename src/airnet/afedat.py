# SPDX-FileCopyrightText: 2023-present Oak Ridge National Laboratory, managed by UT-Battelle
#
# SPDX-License-Identifier: BSD-3-Clause
import math

class AirflowElement:
    pass

class Afe_Plr(AirflowElement):
    def __init__(self, init = 0.0, lam = 0.0, turb = 0.0, expt = 0.5):
        self.init = init # laminar initialization coefficient
        self.lam = lam # laminar flow coefficient
        self.turb = turb # turbulent flow coefficient
        self.expt = expt # turbulent flow exponent

    def calculate(self, link, pdrop):
        if pdrop > 0.0:
            cdm = self.lam * link.node0.dvisc
            fl = cdm * pdrop
            ft = self.turb * link.node0.sqrtd * math.pow(pdrop, self.expt)
            if fl <= ft:
                f = fl
                df = cdm
            else:
                f = ft
                df = ft * self.expt / pdrop
        elif pdrop < 0.0:
            cdm = self.lam * link.node1.dvisc
            fl = cdm * pdrop
            ft = -self.turb * link.node1.sqrtd * math.pow(-pdrop, self.expt)
            if fl >= ft:
                f = fl
                df = cdm
            else:
                f = ft
                df = ft * self.expt / pdrop
        else:
            cdm = 0.5 * self.lam * (link.node0.dvisc + link.node1.dvisc) # original code used node0
            f = fl = ft = 0.0
            df = cdm
        return f, 0.0, df, 0.0

class Afe_Dwc(AirflowElement):
    def __init__(self, length = 0.0, hdia = 0.0, area = 0.0, rough = 0.0, tdlc = 0.0,
                 lflc = 0.0, ldlc = 0.0, linit = 0.0, ed = 0.0, ld = 0.0, f = 0.0):
        self.length = length # length of the duct (m)
        self.hdia = hdia # hydraulic diameter (m)
        self.area = area # cross sectional area (m^2)
        self.rough = rough # roughness dimension (m)
        self.tdlc = tdlc # turbulent dynamic loss coefficient
        self.lflc = lflc # laminar friction loss coefficient
        self.ldlc = ldlc # laminar dynamic loss coefficient
        self.linit = linit # laminar initialization coefficient
        self.ed = ed # relative roughness (rough/hdia)
        self.ld = ld # relative length (length/hdia)
        self.f = f # Darcy friction factor

class Afe_Qfr(AirflowElement):
    def __init__(self, a = 0.0, b = 0.0):
        self.a = a # pdrop = a*f + b*f*f
        self.b = b # 

class Afe_Dor(AirflowElement):
    def __init__(self, init = 0.0, lam = 0.0, turb = 0.0, expt = 0.5, dtmin = 0.0,
                 ht = 0.0, wd = 0.0, cd = 0.0):
        self.init = init # laminar initialization coefficient
        self.lam = lam # laminar flow coefficient
        self.turb = turb # turbulent flow coefficient
        self.expt = expt # turbulent flow exponent (0.5)
        self.dtmin = dtmin # minimum temperature difference for two-way flow (C)
        self.ht = ht # height of doorway (m)
        self.wd = wd # width of doorway (m)
        self.cd = cd # discharge coefficient

class Afe_Cfr(AirflowElement):
    def __init__(self, flow = 0.0):
        self.flow = flow # flow rate (kg/s)

class Afe_Fan(AirflowElement):
    def __init__(self, init = 0.0, lam = 0.0, turb = 0.0, expt = 0.5, rdens = 0.0,
                 fdf = 0.0, sop = 0.0, off = 0.0, mf1 = 0.0, pts = None): #prl = None, fpc = None):
        self.init = init # laminar initialization coefficient
        self.lam = lam # laminar flow coefficient
        self.turb = turb # turbulent flow coefficient
        self.expt = expt # turbulent flow exponent
        self.rdens = rdens # reference fluid density (kg/m^3)
        self.fdf = fdf # free delivery flow (prise = 0) (kg/s)
        self.sop = sop # shut-off pressure (flow = 0) (Pa)
        self.off = off # fan is off if (RPM/rated RPM) < off
        #self.*prl = *prl # vector of pressure range limits [0..nfr]
        #self.**fpc = **fpc # array of fan performance coefficients [1..nfr][0..3]

class Afe_Cpf(AirflowElement):
    def __init__(self, upo = 0.0, prmin = 0.0, ftyp = 0.0):
        self.upo = upo # useful power output (W)
        self.prmin = prmin # minimum pressure rise (Pa)
        self.ftyp = ftyp # typical mass flow rate (kg/s)

class Afe_Ckv(AirflowElement):
    def __init__(self, dp0 = 0.0, coef = 0.0):
        self.dp0 = dp0 # cut-off pressure
        self.coef = coef # flow coefficient

class Afe_Prv(AirflowElement):
    def __init__(self, fpos = 0.0, cpos = 0.0, fneg = 0.0, cneg = 0.0):
        self.fpos = fpos # design flow rate - positive direction (kg/s)
        self.cpos = cpos # positive pressure coefficient
        self.fneg = fneg # design flow rate - negative direction (kg/s)
        self.cneg = cneg # negitive pressure coefficient

object_lookup = {'plr': Afe_Plr, 'dwc': Afe_Dwc, 'qfr': Afe_Qfr,
                 'dor': Afe_Dor, 'cfr': Afe_Cfr, 'fan': Afe_Fan,
                 'cpf': Afe_Cpf, 'ckv': Afe_Ckv, 'prv': Afe_Prv}