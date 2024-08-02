# SPDX-FileCopyrightText: 2023-present Oak Ridge National Laboratory, managed by UT-Battelle
#
# SPDX-License-Identifier: BSD-3-Clause
import math

class Afe_Plr:
    """Power law flow element.
    
    Parameters
    ----------
    init: optional
        The linear coefficient used in initialization.
    lam: optional
        The linear flow coefficient used in simulation.
    turb: optional
        The nonlinear flow coefficient used in simulation.
    expt: optional
        The nonlinear flow exponent used in simulation.
    """
    def __init__(self, init:float=0.0, lam:float=0.0, turb:float=0.0, expt:float=0.5):
        self.init = init # laminar initialization coefficient
        self.lam = lam # laminar flow coefficient
        self.turb = turb # turbulent flow coefficient
        self.expt = expt # turbulent flow exponent

    def type(self):
        """Returns the flow element three-character type.
        
        Returns
        -------
        str:
            The three-character type of the flow element.
        """
        return 'plr'

    def linearize(self, link, multiplier:float=1.0, control:float=1.0):
        """Compute the linear flow coeffient.
        
        The function computes the linear coefficient for use in initialization.
        
        Parameters
        ----------
        link:
            The link object that the linear coefficient is needed for.
        multiplier: optional
            The multiplier determines how many elements this linkage represents. Defaults to 1.
        control: optional
            The control signal that is used to control flow for some elements. Defaults to 1.
        """
        return 0.5 * self.init * (link.node0.dvisc + link.node1.dvisc) # original code used node1

    def flow(self, link, pdrop:float, multiplier:float=1.0, control:float=1.0):
        """Compute the flow through an element.
        
        Parameters
        ----------
        link:
            The link object that the linear coefficient is needed for.
        pdrop:
            The pressure drop driving the flow.
        multiplier: optional
            The multiplier determines how many elements this linkage represents. Defaults to 1.
        control: optional
            The control signal that is used to control flow for some elements. Defaults to 1.
        """
        f = fl = ft = 0.0
        if pdrop > 0.0:
            cdm = self.lam * link.node0.dvisc
            fl = cdm * pdrop
            ft = self.turb * link.node0.sqrt_density * math.pow(pdrop, self.expt)
            if fl <= ft:
                f = fl
            else:
                f = ft
        elif pdrop < 0.0:
            cdm = self.lam * link.node1.dvisc
            fl = cdm * pdrop
            ft = -self.turb * link.node1.sqrt_density * math.pow(-pdrop, self.expt)
            if fl >= ft:
                f = fl
            else:
                f = ft
        return 1, f, 0.0

    def jacobian(self, link, pdrop:float, multiplier:float=1.0, control:float=1.0):
        """Compute the flow through an element and the associated Jacobian terms.
        
        Parameters
        ----------
        link:
            The link object that the linear coefficient is needed for.
        pdrop:
            The pressure drop driving the flow.
        multiplier: optional
            The multiplier determines how many elements this linkage represents. Defaults to 1.
        control: optional
            The control signal that is used to control flow for some elements. Defaults to 1.
        """
        if pdrop > 0.0:
            cdm = self.lam * link.node0.dvisc
            fl = cdm * pdrop
            ft = self.turb * link.node0.sqrt_density * math.pow(pdrop, self.expt)
            if fl <= ft:
                f = fl
                df = cdm
            else:
                f = ft
                df = ft * self.expt / pdrop
        elif pdrop < 0.0:
            cdm = self.lam * link.node1.dvisc
            fl = cdm * pdrop
            ft = -self.turb * link.node1.sqrt_density * math.pow(-pdrop, self.expt)
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
        return 1, f, 0.0, df, 0.0

class Afe_Dwc:
    """Darcy-Weisbach duct flow element."""
    def __init__(self, length:float=0.0, hdia:float=0.0, area:float=0.0, rough:float=0.0, tdlc:float=0.0,
                 lflc:float=0.0, ldlc:float=0.0, linit:float=0.0, ed:float=0.0, ld:float=0.0, f:float=0.0):
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
        """Returns the flow element three-character type.
        
        Returns
        -------
        str:
            The three-character type of the flow element.
        """
        return 'dwc'

class Afe_Qfr:
    """Quadratic flow element.
    
    Parameters
    ----------
    a:
        First order coefficient the quadratic flow relation.
    b:
        Second order coefficent in the quadratic flow relation.
    """
    def __init__(self, a:float=0.0, b:float=0.0):
        self.a = a # pdrop = a*f + b*f*f
        self.b = b #

    def type(self):
        """Returns the flow element three-character type.
        
        Returns
        -------
        str:
            The three-character type of the flow element.
        """
        return 'qfr'

class Afe_Dor(Afe_Plr):
    sqrt_two = 1.414214
    two_thirds = 0.666667
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
    
    def type(self):
        """Returns the flow element three-character type.
        
        Returns
        -------
        str:
            The three-character type of the flow element.
        """
        return 'dor'
    
    def one_way_flow(self, link, pdrop:float):
        return (link.node0.temperature - link.node1.temperature) < self.dtmin

    def flow(self, link, pdrop:float, multiplier:float=1.0, control:float=1.0):
        """Compute the flow through an element.
        
        Parameters
        ----------
        link:
            The link object that the linear coefficient is needed for.
        pdrop:
            The pressure drop driving the flow.
        multiplier: optional
            The multiplier determines how many elements this linkage represents. Defaults to 1.
        control: optional
            The control signal that is used to control flow for some elements. Defaults to 1.
        """
        f1 = 0.0 # computed flow rate
        f2 = 0.0 # computed flow rate

        drho = link.node0.density - link.node1.density
        gdrho = 9.8 * drho

        if self.one_way_flow(link, pdrop):
            return super().calculate(link, pdrop - 0.5 * self.ht * gdrho)
        else:
            y = pdrop / gdrho # Possible two-way flow

            c = self.sqrt_two * self.wd * self.cd
            f0 = self.two_thirds * c * math.sqrt(abs(gdrho*y))*abs(y)
            dfh = c * math.sqrt(abs((self.ht-y)/gdrho))
            fh = self.two_thirds * dfh * abs(gdrho*(self.ht-y))
            nf = 1

            if y < 0.0: # One-way flow (1 to 0)
                if drho > 0.0:
                    f1 = -link.node1.sqrt_density * abs(fh-f0)
                else:
                    f1 =  link.node0.sqrt_density * abs(fh-f0)
            elif y > self.ht: # One-way flow (0 to 1)
                if drho > 0.0:
                    f1 =  link.node0.sqrt_density * abs(fh-f0)
                else:
                    f1 = -link.node1.sqrt_density * abs(fh-f0)
            else: # Two-way flow
                nf = 2
                if drho > 0.0:
                    f1 = -link.node1.sqrt_density * fh
                    f2 =  link.node0.sqrt_density * f0
                else:
                    f1 =  link.node0.sqrt_density * fh
                    f2 = -link.node1.sqrt_density * f0
        return nf, f1, f2

    def jacobian(self, link, pdrop:float, multiplier:float=1.0, control:float=1.0):
        """Compute the flow through an element and the associated Jacobian terms.
        
        Parameters
        ----------
        link:
            The link object that the linear coefficient is needed for.
        pdrop:
            The pressure drop driving the flow.
        multiplier: optional
            The multiplier determines how many elements this linkage represents. Defaults to 1.
        control: optional
            The control signal that is used to control flow for some elements. Defaults to 1.
        """
        f1 = 0.0 # computed flow rate
        df1 = 0.0 # partial derivative: df/dp
        f2 = 0.0 # computed flow rate
        df2 = 0.0 # partial derivative: df/dp

        drho = link.node0.density - link.node1.density
        gdrho = 9.8 * drho

        if self.one_way_flow(link, pdrop):
            return super().calculate(link, pdrop - 0.5 * self.ht * gdrho)
        else:
            y = pdrop / gdrho # Possible two-way flow

            c = self.sqrt_two * self.wd * self.cd
            df0 = c * math.sqrt(abs(pdrop))/abs(gdrho)
            f0 = self.two_thirds * c * math.sqrt(abs(gdrho*y))*abs(y)
            dfh = c * math.sqrt(abs((self.ht-y)/gdrho))
            fh = self.two_thirds * dfh * abs(gdrho*(self.ht-y))
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
    def __init__(self, flow:float=0.0):
        self.flow = flow # flow rate (kg/s)
    
    def type(self):
        """Returns the flow element three-character type.
        
        Returns
        -------
        str:
            The three-character type of the flow element.
        """
        return 'cfr'

class Afe_Fan(Afe_Plr):
    def __init__(self, init:float=0.0, lam:float=0.0, turb:float=0.0, expt:float=0.5, rdens:float=0.0,
                 fdf:float=0.0, sop:float=0.0, off:float=0.0, mf1:float=0.0, pts:float=None): #prl = None, fpc = None):
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
        """Returns the flow element three-character type.
        
        Returns
        -------
        str:
            The three-character type of the flow element.
        """
        return 'fan'

class Afe_Cpf:
    def __init__(self, upo = 0.0, prmin = 0.0, ftyp = 0.0):
        self.upo = upo # useful power output (W)
        self.prmin = prmin # minimum pressure rise (Pa)
        self.ftyp = ftyp # typical mass flow rate (kg/s)

    def type(self):
        """Returns the flow element three-character type.
        
        Returns
        -------
        str:
            The three-character type of the flow element.
        """
        return 'cpf'

class Afe_Ckv:
    def __init__(self, dp0:float=0.0, coef:float=0.0):
        self.dp0 = dp0 # cut-off pressure
        self.coef = coef # flow coefficient

    def type(self):
        """Returns the flow element three-character type.
        
        Returns
        -------
        str:
            The three-character type of the flow element.
        """
        return 'ckv'

class Afe_Prv:
    def __init__(self, fpos:float=0.0, cpos:float=0.0, fneg:float=0.0, cneg:float=0.0):
        self.fpos = fpos # design flow rate - positive direction (kg/s)
        self.cpos = cpos # positive pressure coefficient
        self.fneg = fneg # design flow rate - negative direction (kg/s)
        self.cneg = cneg # negitive pressure coefficient

    def type(self):
        """Returns the flow element three-character type.
        
        Returns
        -------
        str:
            The three-character type of the flow element.
        """
        return 'prv'

object_lookup = {'plr': Afe_Plr, 'dwc': Afe_Dwc, 'qfr': Afe_Qfr,
                 'dor': Afe_Dor, 'cfr': Afe_Cfr, 'fan': Afe_Fan,
                 'cpf': Afe_Cpf, 'ckv': Afe_Ckv, 'prv': Afe_Prv}