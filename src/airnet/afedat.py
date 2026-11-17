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

    def linearize(self, link, pdrop:float, multiplier:float=1.0,
                  control:float=1.0):
        """Compute the linear flow coeffient.
        
        The function computes the linear coefficient for use in initialization.
        
        Parameters
        ----------
        link:
            The link object that the linear coefficient is needed for.
        pdrop:
            The pressure drop used to select the upwind air properties.
        multiplier: optional
            The multiplier determines how many elements this linkage represents. Defaults to 1.
        control: optional
            The control signal that is used to control flow for some elements. Defaults to 1.
        """
        if pdrop > 0.0:
            return self.init * link.node0.dvisc
        return self.init * link.node1.dvisc

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
    def __init__(self, length:float=0.0, hdia:float=0.0, area:float=0.0,
                 rough:float=0.0, tdlc:float=0.0, lflc:float=0.0,
                 ldlc:float=0.0, linit:float=0.0):
        self.length = length # length of the duct (m)
        self.hdia = hdia # hydraulic diameter (m)
        self.area = area # cross sectional area (m^2)
        self.rough = rough # roughness dimension (m)
        self.tdlc = tdlc # turbulent dynamic loss coefficient
        self.lflc = lflc # laminar friction loss coefficient
        self.ldlc = ldlc # laminar dynamic loss coefficient
        self.linit = linit # laminar initialization coefficient
        self.ed = rough / hdia # relative roughness (rough/hdia)
        self.ld = length / hdia # relative length (length/hdia)
        estimate = 1.14 - 0.868589 * math.log(self.ed)
        self.f = 1.0 / (estimate * estimate) # Darcy friction factor

    def type(self):
        """Returns the flow element three-character type.
        
        Returns
        -------
        str:
            The three-character type of the flow element.
        """
        return 'dwc'

    def linearize(self, link, pdrop:float, multiplier:float=1.0,
                  control:float=1.0):
        """Compute AIRNET's laminar initialization coefficient."""
        node = link.node0 if pdrop >= 0.0 else link.node1
        return (
            2.0 * node.density * self.area * self.hdia
            / (node.viscosity * self.linit * self.ld)
        )

    def _friction_factor(self, reynolds:float):
        """Update the Darcy friction factor using AIRNET's Colebrook solve."""
        if self.ed > 1.0e-7:
            a = self.ed / 3.7
            b = 2.51 / reynolds
            coefficient = 0.868589
            g = 1.0 / math.sqrt(self.f)
            while True:
                delta = (
                    g + coefficient * math.log(a + g * b)
                ) / (1.0 + coefficient * b / (a + g * b))
                g -= delta
                if abs(delta / g) <= 1.0e-5:
                    break
            self.f = 1.0 / (g * g)
        elif reynolds < 1.0e5:
            self.f = 0.3164 / math.pow(reynolds, 0.25)
        else:
            self.f = 0.0032 + 0.221 / math.pow(reynolds, 0.237)
        return self.f

    def jacobian(self, link, pdrop:float, multiplier:float=1.0,
                 control:float=1.0):
        """Compute duct mass flow and AIRNET's pressure derivative."""
        node = link.node0 if pdrop >= 0.0 else link.node1
        pressure = abs(pdrop)

        if self.ldlc > 0.001:
            a2 = self.ldlc / (2.0 * node.density * self.area * self.area)
            a1 = (
                node.viscosity * self.lflc * self.ld
                / (2.0 * node.density * self.area * self.hdia)
            )
            root = math.sqrt(a1 * a1 + 4.0 * a2 * pressure)
            laminar_flow = (root - a1) / (2.0 * a2)
            laminar_jacobian = 1.0 / root
        else:
            laminar_jacobian = (
                2.0 * node.density * self.area * self.hdia
                / (node.viscosity * self.lflc * self.ld)
            )
            laminar_flow = laminar_jacobian * pressure

        reynolds = (
            laminar_flow * self.hdia / (node.viscosity * self.area)
        )
        turbulent_flow = laminar_flow
        if reynolds > 10.0:
            numerator = math.sqrt(2.0 * node.density * pressure) * self.area
            turbulent_flow = numerator / math.sqrt(
                self.ld * self.f + self.tdlc
            )
            while True:
                previous = turbulent_flow
                reynolds = (
                    turbulent_flow * self.hdia
                    / (node.viscosity * self.area)
                )
                friction = self._friction_factor(reynolds)
                turbulent_flow = numerator / math.sqrt(
                    self.ld * friction + self.tdlc
                )
                if abs(previous - turbulent_flow) / turbulent_flow <= 0.001:
                    break

        if laminar_flow <= turbulent_flow:
            flow = laminar_flow
            jacobian = laminar_jacobian
        else:
            flow = turbulent_flow
            # AIRNET treats the friction factor as locally constant here.
            jacobian = 0.5 * turbulent_flow / pressure

        if pdrop < 0.0:
            flow = -flow
        return 1, flow, 0.0, jacobian, 0.0

    def flow(self, link, pdrop:float, multiplier:float=1.0,
             control:float=1.0):
        """Compute the one-way duct mass flow."""
        nf, flow0, flow1, _, _ = self.jacobian(
            link, pdrop, multiplier, control
        )
        return nf, flow0, flow1

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
        if a <= 0.0 or b <= 0.0:
            raise ValueError('Quadratic flow coefficients must be positive')
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

    def linearize(self, link, pdrop:float, multiplier:float=1.0,
                  control:float=1.0):
        """Return the zero-flow linearization coefficient."""
        return 1.0 / self.a

    def jacobian(self, link, pdrop:float, multiplier:float=1.0,
                 control:float=1.0):
        """Compute quadratic-relation mass flow and its pressure derivative."""
        pressure = abs(pdrop)
        root = math.sqrt(self.a * self.a + 4.0 * self.b * pressure)
        # This equivalent form avoids cancellation as pressure approaches zero.
        flow = 2.0 * pressure / (self.a + root)
        if pdrop < 0.0:
            flow = -flow
        derivative = 1.0 / (self.a + 2.0 * self.b * abs(flow))
        return 1, flow, 0.0, derivative, 0.0

    def flow(self, link, pdrop:float, multiplier:float=1.0,
             control:float=1.0):
        """Compute quadratic-relation mass flow without derivative terms."""
        nf, flow0, flow1, _, _ = self.jacobian(
            link, pdrop, multiplier, control
        )
        return nf, flow0, flow1

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
    
    def one_way_flow(self, link):
        temperature_difference = link.node0.temperature - link.node1.temperature
        drho = link.node0.density - link.node1.density
        return abs(temperature_difference) < self.dtmin or drho == 0.0

    def _power_law_pressure_drop(self, link, pdrop):
        drho = link.node0.density - link.node1.density
        return pdrop - 0.5 * self.ht * 9.8 * drho

    def linearize(self, link, pdrop:float, multiplier:float=1.0,
                  control:float=1.0):
        """Compute the power-law initialization coefficient."""
        adjusted_pdrop = self._power_law_pressure_drop(link, pdrop)
        return super().linearize(link, adjusted_pdrop, multiplier, control)

    def jacobian(self, link, pdrop:float, multiplier:float=1.0,
                 control:float=1.0):
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

        if self.one_way_flow(link):
            adjusted_pdrop = self._power_law_pressure_drop(link, pdrop)
            return super().jacobian(
                link, adjusted_pdrop, multiplier, control
            )

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

    def flow(self, link, pdrop:float, multiplier:float=1.0,
             control:float=1.0):
        """Compute doorway flow without returning derivative terms."""
        nf, flow0, flow1, _, _ = self.jacobian(
            link, pdrop, multiplier, control
        )
        return nf, flow0, flow1

class Afe_Cfr:
    connects_initialization_pressures = False
    uses_input_orientation = True

    def __init__(self, flow:float=0.0):
        self.flow_rate = flow # rated flow rate (kg/s)
    
    def type(self):
        """Returns the flow element three-character type.
        
        Returns
        -------
        str:
            The three-character type of the flow element.
        """
        return 'cfr'

    def initialization(self, link, pdrop:float, multiplier:float=1.0,
                       control:float=1.0):
        """Return the initialization flow and its pressure derivative."""
        return self.flow_rate * control, 0.0

    def jacobian(self, link, pdrop:float, multiplier:float=1.0,
                 control:float=1.0):
        """Return the pressure-independent flow and zero derivative."""
        return 1, self.flow_rate * control, 0.0, 0.0, 0.0

    def flow(self, link, pdrop:float, multiplier:float=1.0,
             control:float=1.0):
        """Return the pressure-independent mass flow."""
        nf, flow0, flow1, _, _ = self.jacobian(
            link, pdrop, multiplier, control
        )
        return nf, flow0, flow1

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
    connects_initialization_pressures = False
    uses_input_orientation = True

    def __init__(self, dp0:float=0.0, coef:float=0.0):
        if coef <= 0.0:
            raise ValueError('Check-valve flow coefficient must be positive')
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

    def linearize(self, link, pdrop:float, multiplier:float=1.0,
                  control:float=1.0):
        """Return AIRNET's closed-valve initialization derivative."""
        return 0.0

    def initialization_flow(self, link, pdrop:float,
                            multiplier:float=1.0, control:float=1.0):
        """Return AIRNET's zero check-valve initialization flow."""
        return 0.0

    def jacobian(self, link, pdrop:float, multiplier:float=1.0,
                 control:float=1.0):
        """Compute the positive-input-direction check-valve flow."""
        pressure_above_cutoff = pdrop - self.dp0
        if pressure_above_cutoff <= 0.0:
            return 1, 0.0, 0.0, 0.0, 0.0

        flow = math.pow(
            self.coef * pressure_above_cutoff, 1.1
        )
        derivative = 1.1 * flow / pressure_above_cutoff
        return 1, flow, 0.0, derivative, 0.0

    def flow(self, link, pdrop:float, multiplier:float=1.0,
             control:float=1.0):
        """Compute check-valve flow without derivative terms."""
        nf, flow0, flow1, _, _ = self.jacobian(
            link, pdrop, multiplier, control
        )
        return nf, flow0, flow1

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
