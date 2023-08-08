# SPDX-FileCopyrightText: 2023-present Oak Ridge National Laboratory, managed by UT-Battelle
#
# SPDX-License-Identifier: BSD-3-Clause
import math

class MissingArgument(Exception):
    pass

class AirflowElement:
    def get_arg(self, dictionary, name0, name1, object):
        try:
            return dictionary.get(name0)
        except KeyError:
            pass
        try:
            return dictionary.get(name1)
        except KeyError:
            raise MissingArgument('Argument "%s" or "%s" is required for %s objects' % (name0, name1, object))


class Afe_Plr(AirflowElement):
    def __init__(self, **kwargs):
        self.initializer = kwargs.get('init', None) # laminar initialization coefficient
        if self.initializer is None:
            self.initializer = kwargs.get('initializer', None)
        self.linear = self.get_arg(kwargs, 'lam', 'linear', 'Afe_Plr') # laminar/linear flow coefficient
        self.coefficient = self.get_arg(kwargs, 'turb', 'coefficient', 'Afe_Plr') # turbulent flow coefficient
        self.exponent = kwargs.get('expt', 0.5) # turbulent flow exponent
        if self.initializer is None:
            self.initializer = self.linear

    def type(self):
        return 'plr'

    def linearize(self, link):
        return 0.5 * self.initializer * (link.node0.dvisc + link.node1.dvisc) # original code used node1

    def calculate(self, link, pdrop):
        if pdrop > 0.0:
            cdm = self.linear * link.node0.dvisc
            fl = cdm * pdrop
            ft = self.coefficient * link.node0.sqrt_density * math.pow(pdrop, self.exponent)
            if fl <= ft:
                f = fl
                df = cdm
            else:
                f = ft
                df = ft * self.exponent / pdrop
        elif pdrop < 0.0:
            cdm = self.linear * link.node1.dvisc
            fl = cdm * pdrop
            ft = -self.coefficient * link.node1.sqrt_density * math.pow(-pdrop, self.exponent)
            if fl >= ft:
                f = fl
                df = cdm
            else:
                f = ft
                df = ft * self.exponent / pdrop
        else:
            cdm = 0.5 * self.linear * (link.node0.dvisc + link.node1.dvisc) # original code used node0
            f = fl = ft = 0.0
            df = cdm
        return 1, f, 0.0, df, 0.0

class Afe_Dwc:
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

    def type(self):
        return 'dwc'

class Afe_Qfr:
    def __init__(self, a = 0.0, b = 0.0):
        self.a = a # pdrop = a*f + b*f*f
        self.b = b #

    def type(self):
        return 'qfr'

class Afe_Dor(Afe_Plr):
    def __init__(self, **kwargs): #init = 0.0, lam = 0.0, turb = 0.0, expt = 0.5, dtmin = 0.0, ht = 0.0, wd = 0.0, cd = 0.0):
        super().__init__(kwargs)
        self.min_temperature_difference = self.get_arg(kwargs, 'dtmin', 'min_temperature_difference') # minimum temperature difference for two-way flow (C)
        self.height = self.get_arg(kwargs, 'ht', 'height') # height of doorway (m)
        self.width = self.get_arg(kwargs, 'wd', 'width') # width of doorway (m)
        self.discharge_coefficient = self.get_arg(kwargs, 'cd', 'discharge_coefficient') # discharge coefficient
    
    def type(self):
        return 'dor'

    def calculate(self, link, pdrop):
        f1 = 0.0 # computed flow rate
        df1 = 0.0 # partial derivative: df/dp
        f2 = 0.0 # computed flow rate
        df2 = 0.0 # partial derivative: df/dp

        drho = link.node0.density - link.node1.density
        dt = link.node0.temperature - link.node1.temperature
        gdrho = 9.8 * drho

        if abs(dt) < self.min_temperature_difference:
            return super().calculate(link, pdrop - 0.5 * self.height * gdrho)
        else:
            y = pdrop / gdrho # Possible two-way flow

            c = 1.414214 * self.width * self.discharge_coefficient
            df0 = c * math.sqrt(abs(pdrop))/abs(gdrho)
            f0 = 0.666667 * c * math.sqrt(abs(gdrho*y))*abs(y)
            dfh = c * math.sqrt(abs((self.height-y)/gdrho))
            fh = 0.666667 * dfh * abs(gdrho*(self.height-y))
            nf = 1

            if y < 0.0: # One-way flow (1 to 0)
                if drho > 0.0:
                    f1 = -link.node1.sqrt_density * abs(fh-f0)
                    df1 = link.node1.sqrt_density * abs(dfh-df0)
                else:
                    f1 =  link.node0.sqrt_density * abs(fh-f0)
                    df1 = link.node0.sqrt_density * abs(dfh-df0)
            elif y > self.ht: # One-way flow (0 to 1)
                if drho > 0.0:
                    f1 =  link.node0.sqrt_density * abs(fh-f0)
                    df1 = link.node0.sqrt_density * abs(dfh-df0)
                else:
                    f1 = -link.node1.sqrt_density * abs(fh-f0)
                    df1 = link.node1.sqrt_density * abs(dfh-df0)
            else: # Two-way flow
                nf = 2
                if drho > 0.0:
                    f1 = -link.node1.sqrt_density * fh
                    df1 = link.node1.sqrt_density * dfh
                    f2 =  link.node0.sqrt_density * f0
                    df2 = link.node0.sqrt_density * df0
                else:
                    f1 =  link.node0.sqrt_density * fh
                    df1 = link.node0.sqrt_density * dfh
                    f2 = -link.node1.sqrt_density * f0
                    df2 = link.node1.sqrt_density * df0
        return nf, f1, f2, df1, df2

class Afe_Cfr:
    def __init__(self, flow = 0.0):
        self.flow = flow # flow rate (kg/s)
    
    def type(self):
        return 'cfr'

class Afe_Fan(Afe_Plr):
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

    def type(self):
        return 'fan'

class Afe_Cpf:
    def __init__(self, upo = 0.0, prmin = 0.0, ftyp = 0.0):
        self.upo = upo # useful power output (W)
        self.prmin = prmin # minimum pressure rise (Pa)
        self.ftyp = ftyp # typical mass flow rate (kg/s)

    def type(self):
        return 'cpf'

class Afe_Ckv:
    def __init__(self, dp0 = 0.0, coef = 0.0):
        self.dp0 = dp0 # cut-off pressure
        self.coef = coef # flow coefficient

    def type(self):
        return 'ckv'

class Afe_Prv:
    def __init__(self, fpos = 0.0, cpos = 0.0, fneg = 0.0, cneg = 0.0):
        self.fpos = fpos # design flow rate - positive direction (kg/s)
        self.cpos = cpos # positive pressure coefficient
        self.fneg = fneg # design flow rate - negative direction (kg/s)
        self.cneg = cneg # negitive pressure coefficient

    def type(self):
        return 'prv'

object_lookup = {'plr': Afe_Plr, 'dwc': Afe_Dwc, 'qfr': Afe_Qfr,
                 'dor': Afe_Dor, 'cfr': Afe_Cfr, 'fan': Afe_Fan,
                 'cpf': Afe_Cpf, 'ckv': Afe_Ckv, 'prv': Afe_Prv}