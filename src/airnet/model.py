# SPDX-FileCopyrightText: 2023-present Oak Ridge National Laboratory, managed by UT-Battelle
#
# SPDX-License-Identifier: BSD-3-Clause
from __future__ import annotations

import math
from os import PathLike

import numpy
import scipy

from .afedat import object_lookup
from .reader import InputType, Reader
from .weather import Weather
from .wind import read_wind_profiles


def devnull(msg: str):
    return

class Node:
    """A class representing a node in the pressure network.
    
    Parameters
    ----------
    name:
        The name of the node.
    variable:
        Boolean flag that determines whether the pressure is variable.
    ht: optional
        The elevation of the node for stack-based pressure calculation.
    temp: optional
        The temperature of the node.
    pres: optional
        The pressure of the node.
    index: optional
        The index of the node in the system of equations, only meaningful for variable nodes.
    input_c: optional
        Flag determining the input temperature units. If True, the temperature is in Celcius, otherwise Kelvin.
    """
    def __init__(self, name:str|None=None, variable:bool=True, ht:float=0.0, temp:float=293.15,
                 pres:float=0.0, index:int|None=None, input_c:bool=True,
                 ambient:bool=False):
        self.name = name
        self.variable = variable
        self.height = ht
        self.temperature = temp
        if input_c:
            self.temperature += 273.15
        self.pressure = pres
        self.index = index
        self.ambient = ambient
        self.density = 0.0
        self.viscosity = 0.0
        self.sqrt_density = 0.0
        self.dvisc = 0.0 # Density divided by viscosity
    #def copy_state(self, other):
    #    self.temperature = other.temperature
    #    self.pressure = other.pressure
    #    self.density = other.density
    #    self.viscosity = other.viscosity
    #    self.sqrt_density = other.sqrt_density
    #    self.dvisc = other.dvisc


class Link:
    def __init__(self, name=None, node0=None, ht0=0.0, node1=None, ht1=0.0, element=None,
                 wind=None, wpmod=0.0, mult=1.0, control=1.0, flipped=False):
        self.name = name
        self.node0 = node0
        self.node1 = node1
        self.ht0 = ht0
        self.ht1 = ht1
        self.element = element
        self.wind = wind
        self.wpmod = wpmod
        self.multiplier = mult
        self.control = control
        self.flipped = flipped
        self.flow_count = 0
        self.flow0 = 0.0
        self.flow1 = 0.0
        self.pdrop = 0.0

    @property
    def input_pressure_drop(self):
        """Pressure drop expressed in the original AIRNET link order."""
        return -self.pdrop if self.flipped else self.pdrop

    @property
    def input_flow0(self):
        """Primary flow expressed in the original AIRNET link order."""
        return -self.flow0 if self.flipped else self.flow0

    @property
    def input_flow1(self):
        """Secondary flow expressed in the original AIRNET link order."""
        return -self.flow1 if self.flipped else self.flow1

class BadNetwork(Exception):
    """Raised if the network is bad."""


class DuplicateNetworkName(BadNetwork):
    """Raised when a node, element, or link name is defined more than once."""


class MissingNetworkReference(BadNetwork):
    """Raised when a link references a node or element that is not defined."""


class UnsupportedNetworkFeature(BadNetwork):
    """Raised when input requests behavior that is not implemented."""


class MissingWeather(BadNetwork):
    """Raised when numerical work is requested without bound weather."""


class SolverError(RuntimeError):
    """Base class for explicit airflow-solver failures."""


class LinearSolverError(SolverError):
    """Raised when a linear solve fails or does not converge."""

    def __init__(self, message, *, info, phase, nonlinear_iteration=None):
        super().__init__(message)
        self.info = info
        self.phase = phase
        self.nonlinear_iteration = nonlinear_iteration


class InitializationError(LinearSolverError):
    """Raised when the linear initialization solve fails."""

    def __init__(self, *, info, max_iterations):
        if info > 0:
            detail = f'did not converge within {max_iterations} iterations'
        else:
            detail = 'failed'
        super().__init__(
            f'Initialization linear solve {detail} (solver info={info})',
            info=info,
            phase='initialization',
        )
        self.max_iterations = max_iterations


class NonlinearSolverError(SolverError):
    """Raised when the nonlinear airflow solve does not converge."""

    def __init__(
        self,
        *,
        max_iterations,
        maximum_absolute_residual,
        relative_residual,
        absolute_tolerance,
        relative_tolerance,
    ):
        super().__init__(
            f'Airflow solve failed to converge within {max_iterations} iterations '
            f'(maximum absolute residual={maximum_absolute_residual:.9e})'
        )
        self.max_iterations = max_iterations
        self.maximum_absolute_residual = maximum_absolute_residual
        self.relative_residual = relative_residual
        self.absolute_tolerance = absolute_tolerance
        self.relative_tolerance = relative_tolerance


class NumericalSolverError(SolverError):
    """Raised when an element or solver produces a non-finite value."""


class StaleSolverError(SolverError):
    """Raised when a model's weather changes after solver construction."""


class Model:
    """Physical network, environmental binding, and current solution state.

    AIRNET network records may be constructed directly or read with
    :meth:`read_airnet`. Weather is bound separately through :attr:`weather`.
    """
    supported_element_types = frozenset(
        {'plr', 'dwc', 'qfr', 'dor', 'cfr', 'fan', 'cpf', 'ckv'}
    )

    @staticmethod
    def read_airnet(
        input_source,
        *,
        element_lookup: dict = object_lookup,
        node_object=Node,
        link_object=Link,
        global_temperature: float | None = None,
        global_density: float | None = None,
    ):
        """Read an AIRNET network source without binding run weather.

        ``input_source`` may be a path or an open text stream. AIRNET weather
        and wind coefficients are separate run inputs, so the returned model's
        :attr:`weather` is always ``None``.
        """
        options = {
            'element_lookup': element_lookup,
            'node_object': node_object,
            'link_object': link_object,
            'global_temperature': global_temperature,
            'global_density': global_density,
        }
        if hasattr(input_source, 'read'):
            return Model(Reader(input_source), **options)
        if not isinstance(input_source, (str, bytes, PathLike)):
            raise TypeError('AIRNET input source must be a path or text stream')
        with open(input_source, encoding='utf-8') as input_stream:
            return Model(Reader(input_stream), **options)

    def __init__(self, network_input, element_lookup:dict = object_lookup, node_object=Node,
                 link_object=Link, global_temperature:float|None=None, global_density:float|None=None,
                 weather:Weather|None=None):
        self.title = ''
        self.nodes = {}
        self.links = []
        self.elements = {}
        self._weather = None
        self._weather_revision = 0
        self._global_temperature = global_temperature
        self._global_density = global_density

        # Copy every record before removing construction-only fields. Reader
        # output is also part of the public API and remains owned by the caller.
        links = []
        link_names = set()
        title_seen = False
        for input_record in network_input:
            record = dict(input_record)
            try:
                input_type = record.pop('input_type')
            except KeyError as error:
                raise BadNetwork('Network input record has no input_type') from error

            if input_type == InputType.TITLE:
                if title_seen:
                    raise BadNetwork('Network input contains more than one title')
                title_seen = True
                self.title = record.get('title', '')
            elif input_type == InputType.NODE:
                try:
                    name = record['name']
                except KeyError as error:
                    raise BadNetwork('Node record has no name') from error
                if name in self.nodes:
                    raise DuplicateNetworkName(f'Node "{name}" is defined more than once')
                node_type = record.pop('type', 'c')
                if node_type not in {'a', 'c', 'v'}:
                    raise BadNetwork(f'Node "{name}" has invalid type "{node_type}"')
                try:
                    self.nodes[name] = node_object(
                        variable=node_type == 'v', ambient=node_type == 'a', **record
                    )
                except (TypeError, ValueError) as error:
                    raise BadNetwork(f'Node "{name}" has invalid data: {error}') from error
            elif input_type == InputType.LINK:
                try:
                    name = record['name']
                except KeyError as error:
                    raise BadNetwork('Link record has no name') from error
                if name in link_names:
                    raise DuplicateNetworkName(f'Link "{name}" is defined more than once')
                link_names.add(name)
                links.append(record)
            elif input_type == InputType.ELEMENT:
                try:
                    element_type = record.pop('type')
                    name = record.pop('name')
                except KeyError as error:
                    raise BadNetwork(
                        f'Element record has no {error.args[0]}'
                    ) from error
                if name in self.elements:
                    raise DuplicateNetworkName(
                        f'Element "{name}" is defined more than once'
                    )
                if element_type not in element_lookup:
                    raise UnsupportedNetworkFeature(
                        f'Element "{name}" has unsupported type "{element_type}"'
                    )
                if (
                    element_lookup is object_lookup
                    and element_type not in self.supported_element_types
                ):
                    raise UnsupportedNetworkFeature(
                        f'Element "{name}" uses unsupported solver element '
                        f'type "{element_type}"'
                    )
                try:
                    element = element_lookup[element_type](**record)
                except (TypeError, ValueError, ZeroDivisionError) as error:
                    raise BadNetwork(
                        f'Element "{name}" has invalid data: {error}'
                    ) from error
                missing_methods = [
                    method
                    for method in ('linearize', 'flow', 'calculate')
                    if not callable(getattr(element, method, None))
                ]
                if missing_methods:
                    raise UnsupportedNetworkFeature(
                        f'Element "{name}" does not implement '
                        f'{", ".join(missing_methods)}'
                    )
                self.elements[name] = element
            else:
                raise UnsupportedNetworkFeature(
                    f'Unsupported network input type: {input_type!r}'
                )

        adjacency = {name: set() for name in self.nodes}
        pressure_adjacency = {name: set() for name in self.nodes}
        for link in links:
            name = link['name']
            wind = link.get('wind')
            wpmod = link.get('wpmod', 0.0)
            if not math.isfinite(wpmod):
                raise BadNetwork(f'Link "{name}" wind-pressure modifier must be finite')
            multiplier = link.get('mult', 1.0)
            if not math.isfinite(multiplier) or multiplier <= 0.0:
                raise BadNetwork(
                    f'Link "{name}" multiplier must be finite and positive'
                )

            try:
                node0_name = link.pop('node-1')
                node1_name = link.pop('node-2')
                element_name = link.pop('element')
                ht0 = link.pop('ht-1')
                ht1 = link.pop('ht-2')
            except KeyError as error:
                raise BadNetwork(
                    f'Link "{name}" has no {error.args[0]}'
                ) from error
            if node0_name not in self.nodes:
                raise MissingNetworkReference(
                    f'Link "{name}" references undefined node "{node0_name}"'
                )
            if node1_name not in self.nodes:
                raise MissingNetworkReference(
                    f'Link "{name}" references undefined node "{node1_name}"'
                )
            if element_name not in self.elements:
                raise MissingNetworkReference(
                    f'Link "{name}" references undefined element "{element_name}"'
                )
            if node0_name == node1_name:
                raise BadNetwork(f'Link "{name}" connects node "{node0_name}" to itself')

            node0 = self.nodes[node0_name]
            node1 = self.nodes[node1_name]
            element = self.elements[element_name]
            if 'control' in link and 'ctrl' in link:
                raise BadNetwork(
                    f'Link "{name}" defines both control and ctrl'
                )
            control = link.pop('control', link.pop('ctrl', 1.0))
            if not math.isfinite(control) or not 0.0 <= control <= 1.0:
                raise BadNetwork(
                    f'Link "{name}" control must be between 0 and 1'
                )
            if control != 1.0 and element.type() not in {'fan', 'cpf', 'cfr'}:
                raise UnsupportedNetworkFeature(
                    f'Link "{name}" controls an element that AIRNET does not control'
                )
            link['control'] = control
            adjacency[node0_name].add(node1_name)
            adjacency[node1_name].add(node0_name)
            if element.type() not in ['ckv', 'cfr']:
                pressure_adjacency[node0_name].add(node1_name)
                pressure_adjacency[node1_name].add(node0_name)
            if wind in (None, 'null'):
                link['wind'] = None
            else:
                link['wind'] = wind
            try:
                self.links.append(
                    link_object(
                        node0=node0,
                        ht0=ht0,
                        node1=node1,
                        ht1=ht1,
                        element=element,
                        **link,
                        flipped=False,
                    )
                )
            except (TypeError, ValueError) as error:
                raise BadNetwork(f'Link "{name}" has invalid data: {error}') from error

        disconnected = sorted(name for name, neighbors in adjacency.items() if not neighbors)
        if disconnected:
            raise BadNetwork(f'Disconnected nodes found: {", ".join(disconnected)}')

        visited = set()
        for start, start_node in self.nodes.items():
            if not start_node.variable:
                continue
            if start in visited:
                continue
            pending = [start]
            component = []
            while pending:
                name = pending.pop()
                if name in visited:
                    continue
                visited.add(name)
                component.append(name)
                pending.extend(pressure_adjacency[name] - visited)
            if not any(not self.nodes[name].variable for name in component):
                raise BadNetwork(
                    'Variable-pressure component has no linear-initialization '
                    'path to a fixed-pressure node: '
                    + ', '.join(sorted(component))
                )

        if weather is not None:
            self.weather = weather

    @property
    def weather(self) -> Weather | None:
        """Environmental inputs currently bound to this model."""
        return self._weather

    @weather.setter
    def weather(self, weather: Weather | None):
        if weather is None:
            self._weather = None
            self._weather_revision += 1
            return
        if not isinstance(weather, Weather):
            raise TypeError('Model weather must be a Weather object or None')

        for link in self.links:
            if link.wind is not None and link.wind not in weather.wind_profiles:
                raise MissingNetworkReference(
                    f'Link "{link.name}" references undefined wind profile '
                    f'"{link.wind}"'
                )

        property_states = self._property_states(
            self.nodes.values(),
            weather,
            global_temperature=self._global_temperature,
            global_density=self._global_density,
        )
        self._weather = weather
        self._apply_property_states(property_states)
        self._weather_revision += 1

    def _require_weather(self) -> Weather:
        if self.weather is None:
            raise MissingWeather(
                'Model has no weather; assign model.weather before solving'
            )
        return self.weather

    @staticmethod
    def _apply_property_states(property_states):
        for node, temperature, density, viscosity in property_states:
            node.temperature = temperature
            node.density = density
            node.sqrt_density = math.sqrt(density)
            node.viscosity = viscosity
            node.dvisc = density / viscosity

    @staticmethod
    def _property_states(
        nodes,
        weather,
        *,
        global_temperature=None,
        global_density=None,
    ):
        property_states = []
        for node in nodes:
            if global_temperature is not None:
                temperature = global_temperature
            elif node.ambient:
                temperature = weather.ambient_temperature + 273.15
            else:
                temperature = node.temperature

            if not math.isfinite(temperature) or temperature <= 0.0:
                raise BadNetwork(
                    f'Node "{node.name}" temperature must be finite and positive'
                )
            density = (
                global_density
                if global_density is not None
                else 0.0034838
                * (weather.barometric_pressure + node.pressure)
                / temperature
            )
            if not math.isfinite(density) or density <= 0.0:
                raise BadNetwork(
                    f'Node "{node.name}" density must be finite and positive'
                )
            viscosity = 1.71432e-5 + 4.828e-8 * (temperature - 273.15)
            if not math.isfinite(viscosity) or viscosity <= 0.0:
                raise BadNetwork(
                    f'Node "{node.name}" viscosity must be finite and positive'
                )
            property_states.append((node, temperature, density, viscosity))
        return property_states

    def summary(self):
        """Summarize the pressure network in the model.
        
        Returns
        -------
        str:
            A multiline string that describes the model.
        """
        string = 'Title: %s\n\nElements:\n=========\n' % self.title
        elements = {}
        for el in self.elements.values():
            tag = el.type()
            if tag in elements:
                elements[tag] += 1
            else:
                elements[tag] = 1
        for key, value in elements.items():
            string += '%s: %d\n' % (key, value)

        variable_count = sum(node.variable for node in self.nodes.values())
        string += '\nNodes: %s\n\nLinks: %s\n' % (len(self.nodes), len(self.links))
        string += '\nSystem size: %d x %x\n' % (variable_count, variable_count)
        return string
    
    def results_summary(self):
        """Summarize the result of simulation of the pressure network.
        
        Returns
        -------
        str:
            A multiline string that describes the results.
        """
        string = 'Title: %s\n\nNodes:\n======\n' % self.title
        for name, node in self.nodes.items():
            nr = node.index
            if nr is None:
                nr = 0
            string += '%4d %s: %e %e %e\n' % (nr, node.name, node.pressure, node.temperature, node.density)
        string += '\n'
    
        variable_count = sum(node.variable for node in self.nodes.values())
        string += '\nNodes: %s\n\nLinks: %s\n' % (len(self.nodes), len(self.links))
        string += '\nSystem size: %d x %x\n' % (variable_count, variable_count)
        return string

    def set_properties(self, nodes:list, global_temperature:float|None=None, global_density:float|None=None):
        """Loop and set the properties of the nodes.
        
        Parameters
        ----------
        nodes:
            The list of nodes to modify.
        global_temperature: optional
            If not None, use this temperature everywhere.
        global_density: optional
            If not None, use this density everywhere.
        """
        weather = self._require_weather()
        property_states = self._property_states(
            nodes,
            weather,
            global_temperature=global_temperature,
            global_density=global_density,
        )
        self._apply_property_states(property_states)

    def initialize(self, maxiter:int=100):
        """Initialize through a temporary :class:`Solver`.

        This compatibility entry point remains until the solver algorithms are
        moved in the next partitioning slice. New code should construct one
        solver and call :meth:`Solver.initialize` on it.
        """
        return Solver(self).initialize(maxiter=maxiter)

    def _initialize(self, solver, maxiter:int=100):
        """Run the linear initialization algorithm using solver-owned state.
        
        Compute flows and pressure drops using linear flow representation for all elements.
        
        Parameters
        ----------
        maxiter: optional
            The maximum number of iterations allowed in the solution.
            
        Raises
        ------
        InitializationError
            If the linear initialization solve fails or does not converge.
        NumericalSolverError
            If initialization produces a non-finite coefficient, solution,
            or flow.
        """
        self._require_weather()
        if maxiter < 1:
            raise ValueError('Maximum initialization iterations must be at least 1')

        # AIRNET includes stack and wind terms in the pressure drop passed to
        # each element during linear initialization.
        solver._compute_environmental_pressure_drops()
        solver._assemble_initialization()
        preconditioner = solver._linear_preconditioner('initialization')
        solution, info = scipy.sparse.linalg.cg(
            solver.A,
            solver.x,
            rtol=1.0e-12,
            atol=0.0,
            maxiter=maxiter,
            M=preconditioner,
        )
        if info != 0:
            raise InitializationError(info=info, max_iterations=maxiter)
        if not bool(numpy.all(numpy.isfinite(solution))):
            raise NumericalSolverError(
                'Initialization linear solve produced a non-finite pressure'
            )

        solver.x = solution
        # Update the nodal pressures
        for node in self.nodes.values():
            if node.variable:
                node.pressure = solver.x[node.index]
        # Update the flows
        solver._synchronize_pressure_drops()
        for link in self.links:
            pdrop = link.pdrop
            flow, derivative = link.element.linearize(
                link, pdrop=pdrop, control=link.control
            )
            if not math.isfinite(derivative) or not math.isfinite(flow):
                raise NumericalSolverError(
                    f'Element on link "{link.name}" produced a non-finite '
                    'initialization flow'
                )
            link.flow_count = 1
            link.flow0 = flow
            link.flow1 = 0.0
    
    def compute_pressure_drops(self):
        """Compute the pressure drop across all the links in the model."""
        weather = self._require_weather()
        for link in self.links:
            # Stack contribution
            sp0 = -9.80 * link.node0.density * link.ht0
            sp1 =  9.80 * link.node1.density * link.ht1
            dhx = (link.node0.height - link.node1.height) + (link.ht0 - link.ht1)
            spx = 0.0
            if dhx != 0.0:
                if link.flow0 > 0.0:
                    spx = 9.80 * link.node0.density * dhx
                elif link.flow0 < 0.0:
                    spx = 9.80 * link.node1.density * dhx
                else:
                    spx = 4.90 * (link.node0.density + link.node1.density) * dhx
            wind_pressure = 0.0
            if link.wind is not None:
                dynamic_pressure = (
                    0.5 * weather.ambient_density * weather.wind_speed**2
                )
                wind_pressure = (
                    link.wpmod
                    * dynamic_pressure
                    * weather.wind_profiles[link.wind].coefficient(
                        weather.wind_direction
                    )
                )
            link.pdrop = sp0 + sp1 + spx + wind_pressure

    def air_movement(self, maxiter:int=100, max_subiter:int=100,
                     tolerance:float=1.0e-8, status_function=devnull,
                     relative_tolerance:float|None=None, relaxation:float=1.0,
                     minimum_iterations:int=1):
        """Solve through a temporary :class:`Solver`.

        This compatibility entry point remains until the solver algorithms are
        moved in the next partitioning slice. New code should use
        :meth:`Solver.steady_solve`.
        """
        return Solver(self).steady_solve(
            maxiter=maxiter,
            max_subiter=max_subiter,
            tolerance=tolerance,
            status_function=status_function,
            relative_tolerance=relative_tolerance,
            relaxation=relaxation,
            minimum_iterations=minimum_iterations,
        )

    def _air_movement(self, solver, maxiter:int=100, max_subiter:int=100,
                      tolerance:float=1.0e-8, status_function=devnull,
                      relative_tolerance:float|None=None, relaxation:float=1.0,
                      minimum_iterations:int=1):
        """Run the nonlinear algorithm using solver-owned state.
        
        Parameters
        ----------
        maxiter: optional
            The maximum number of Newton iterations allowed.
        max_subiter: optional
            The maximum number of conjugate gradient iterations allowed in the linear solve.
        tolerance: optional
            The absolute convergence tolerance.
        relative_tolerance: optional
            The per-node relative airflow convergence tolerance. When omitted,
            convergence is based only on ``tolerance``.
        relaxation: optional
            The AIRNET under-relaxation coefficient. A correction is relaxed
            when the global relative residual is not falling rapidly.
        minimum_iterations: optional
            The minimum number of nonlinear residual evaluations.
        status_function: optional
            Function that handles message output, the default is to discard all messages.
        
        Returns
        -------
        int:
            The number of Newton iterations used in the solve.

        Raises
        ------
        LinearSolverError
            If a linear correction fails or does not converge.
        NonlinearSolverError
            If the airflow iteration does not converge within ``maxiter``
            iterations.
        NumericalSolverError
            If an element or linear solve produces a non-finite value.
        """
        self._require_weather()
        if maxiter < 1:
            raise ValueError('Maximum nonlinear iterations must be at least 1')
        if max_subiter < 1:
            raise ValueError('Maximum linear iterations must be at least 1')
        if tolerance < 0.0:
            raise ValueError('Absolute convergence tolerance cannot be negative')
        if relative_tolerance is not None and relative_tolerance < 0.0:
            raise ValueError('Relative convergence tolerance cannot be negative')
        if not 0.0 < relaxation <= 1.0:
            raise ValueError('Relaxation coefficient must be greater than 0 and at most 1')
        if minimum_iterations < 1:
            raise ValueError('Minimum iterations must be at least 1')

        if not solver._environmental_pressure_is_set:
            solver._compute_environmental_pressure_drops()
        solver._synchronize_pressure_drops()
        status_function('iter | Max Resid|\n==== ===============')
        previous_relative_residual = 0.0
        for iter in range(1,maxiter+1):
            sum_abs_flow = solver._assemble_nonlinear(iter)
            residual = numpy.asarray(solver.x).reshape(-1)
            if not bool(numpy.all(numpy.isfinite(residual))):
                raise NumericalSolverError(
                    f'Non-finite residual at nonlinear iteration {iter}'
                )
            absolute_residual = numpy.abs(residual)
            maxf = float(absolute_residual.max(initial=0.0))
            converged = absolute_residual < tolerance
            if relative_tolerance is not None:
                relative_by_node = numpy.full(solver.size, numpy.inf)
                numpy.divide(
                    absolute_residual,
                    sum_abs_flow,
                    out=relative_by_node,
                    where=sum_abs_flow > 0.0,
                )
                converged |= relative_by_node <= relative_tolerance

            total_abs_flow = float(sum_abs_flow.sum())
            if total_abs_flow > 0.0:
                relative_residual = float(absolute_residual.sum()) / total_abs_flow
            else:
                relative_residual = 0.0

            if bool(numpy.all(converged)) and iter >= minimum_iterations:
                status_function('%4d %15.9e < %e' % (iter, maxf, tolerance))
                return iter

            status_function('%4d %15.9e' % (iter, maxf))
            preconditioner = solver._linear_preconditioner(
                f'nonlinear iteration {iter}'
            )
            solver.correction, info = scipy.sparse.linalg.cg(
                solver.A,
                residual,
                rtol=1.0e-12,
                atol=0.0,
                maxiter=max_subiter,
                M=preconditioner,
            )
            if info != 0:
                raise LinearSolverError(
                    f'Linear correction solve failed at nonlinear iteration {iter} '
                    f'(solver info={info})',
                    info=info,
                    phase='nonlinear correction',
                    nonlinear_iteration=iter,
                )
            if not bool(numpy.all(numpy.isfinite(solver.correction))):
                raise NumericalSolverError(
                    f'Linear correction solve produced a non-finite value at '
                    f'nonlinear iteration {iter}'
                )

            correction_factor = 1.0
            if relative_residual > 0.30 * previous_relative_residual:
                correction_factor = relaxation
            previous_relative_residual = relative_residual

            # Update the nodal pressures; flows are set above.
            for node in self.nodes.values():
                if node.variable:
                    node.pressure -= (
                        correction_factor * solver.correction[node.index]
                    )
            # Update pressure drops incrementally, as AIRNET does, to avoid
            # cancellation when two large node pressures differ only slightly.
            for link in self.links:
                correction0 = (
                    solver.correction[link.node0.index]
                    if link.node0.variable
                    else 0.0
                )
                correction1 = (
                    solver.correction[link.node1.index]
                    if link.node1.variable
                    else 0.0
                )
                link.pdrop += correction_factor * (correction1 - correction0)
            solver._pressure_drop_node_pressures = {
                node: node.pressure for node in self.nodes.values()
            }

        raise NonlinearSolverError(
            max_iterations=maxiter,
            maximum_absolute_residual=maxf,
            relative_residual=relative_residual,
            absolute_tolerance=tolerance,
            relative_tolerance=relative_tolerance,
        )


class Solver:
    """Numerical interface for initializing and solving a weather-bound model.

    The solver updates pressure and flow state on its model. Rebinding the
    model's weather invalidates the solver, so each environmental scenario
    requires a newly constructed solver.
    """

    def __init__(self, model: Model):
        if not isinstance(model, Model):
            raise TypeError('Solver model must be a Model object')
        model._require_weather()
        self.model = model
        self._weather_revision = model._weather_revision
        self._environmental_pressure_is_set = False
        self._pressure_drop_node_pressures = None

        # Equation preparation intentionally mutates model state. Nodes retain
        # their row indices, and fixed-to-variable links retain their normalized
        # orientation. Link input-oriented properties reverse result signs for
        # presentation.
        for link in model.links:
            if link.flipped or link.node0.variable or not link.node1.variable:
                continue
            link.node0, link.node1 = link.node1, link.node0
            link.ht0, link.ht1 = link.ht1, link.ht0
            link.multiplier = -link.multiplier
            link.wpmod = -link.wpmod
            link.flipped = True
        variable_count = 0
        for node in model.nodes.values():
            if node.variable:
                node.index = variable_count
                variable_count += 1
            else:
                node.index = None
        self.size = variable_count

        row = []
        col = []
        data = []
        for link in model.links:
            if link.node0.variable:
                row.append(link.node0.index)
                col.append(link.node0.index)
                data.append(1.0)
                if link.node1.variable:
                    row.extend((link.node1.index, link.node0.index, link.node1.index))
                    col.extend((link.node1.index, link.node1.index, link.node0.index))
                    data.extend((1.0, 1.0, 1.0))

        matrix = scipy.sparse.coo_matrix(
            (
                numpy.array(data, dtype=numpy.double),
                (numpy.array(row), numpy.array(col)),
            ),
            shape=(self.size, self.size),
        )
        self.A = scipy.sparse.csr_matrix(matrix)
        self.x = numpy.zeros([self.size, 1], dtype=numpy.double)
        self.correction = numpy.zeros(self.size, dtype=numpy.double)

    def _require_current_weather(self):
        self.model._require_weather()
        if self.model._weather_revision != self._weather_revision:
            raise StaleSolverError(
                'Model weather changed after Solver construction; create a new Solver'
            )

    def _compute_environmental_pressure_drops(self):
        self.model.compute_pressure_drops()
        self._environmental_pressure_is_set = True
        self._pressure_drop_node_pressures = None

    def _synchronize_pressure_drops(self):
        previous = self._pressure_drop_node_pressures
        if previous is None:
            for link in self.model.links:
                link.pdrop += link.node0.pressure - link.node1.pressure
        else:
            for link in self.model.links:
                link.pdrop += (
                    link.node0.pressure
                    - previous[link.node0]
                    - link.node1.pressure
                    + previous[link.node1]
                )
        self._pressure_drop_node_pressures = {
            node: node.pressure for node in self.model.nodes.values()
        }

    def _linear_preconditioner(self, phase):
        diagonal = self.A.diagonal()
        if diagonal.size == 0:
            return None
        if not bool(numpy.all(numpy.isfinite(diagonal))) or bool(
            numpy.any(diagonal <= 0.0)
        ):
            raise NumericalSolverError(
                f'{phase.capitalize()} matrix has a non-positive or non-finite '
                'diagonal coefficient'
            )
        # Preserve the established numerical path for ordinary networks.
        # Diagonal scaling is needed for stress cases such as 133, whose
        # coefficient range spans more than thirteen orders of magnitude.
        if diagonal.max() <= 1.0e6 * diagonal.min():
            return None
        return scipy.sparse.diags(1.0 / diagonal)

    def _assemble_initialization(self):
        self.A.data.fill(0.0)
        self.x.fill(0.0)
        for link in self.model.links:
            if not link.node0.variable:
                continue
            pdrop = link.node0.pressure - link.node1.pressure + link.pdrop
            flow, derivative = link.element.linearize(
                link, pdrop=pdrop, control=link.control
            )
            right_hand_side = -flow
            if not math.isfinite(derivative) or not math.isfinite(
                right_hand_side
            ):
                raise NumericalSolverError(
                    f'Element on link "{link.name}" produced a non-finite '
                    'initialization derivative or residual flow'
                )
            self.x[link.node0.index] += right_hand_side
            self.A[link.node0.index, link.node0.index] += derivative
            if link.node1.variable:
                self.x[link.node1.index] -= right_hand_side
                self.A[link.node1.index, link.node1.index] += derivative
                self.A[link.node0.index, link.node1.index] -= derivative
                self.A[link.node1.index, link.node0.index] -= derivative

    def _assemble_nonlinear(self, iteration):
        self.A.data.fill(0.0)
        self.x.fill(0.0)
        sum_abs_flow = numpy.zeros(self.size, dtype=numpy.double)
        for link in self.model.links:
            pdrop = link.pdrop
            if link.node0.variable:
                nf, link.flow0, link.flow1, df0, df1 = link.element.calculate(
                    link, pdrop, control=link.control
                )
                if not all(
                    math.isfinite(value)
                    for value in (link.flow0, link.flow1, df0, df1)
                ):
                    raise NumericalSolverError(
                        f'Element on link "{link.name}" produced a non-finite '
                        f'flow or derivative at nonlinear iteration {iteration}'
                    )
                link.flow_count = nf
                flow = link.flow0 + link.flow1
                derivative = df0 + df1
                # AIRNET assembles the net flow and net derivative while
                # retaining both components on the link for reporting.
                self.A[link.node0.index, link.node0.index] += derivative
                self.x[link.node0.index] += flow
                sum_abs_flow[link.node0.index] += abs(flow)
                if link.node1.variable:
                    self.A[link.node1.index, link.node1.index] += derivative
                    self.x[link.node1.index] -= flow
                    sum_abs_flow[link.node1.index] += abs(flow)
                    self.A[link.node0.index, link.node1.index] -= derivative
                    self.A[link.node1.index, link.node0.index] -= derivative
            else:
                nf, link.flow0, link.flow1 = link.element.flow(
                    link, pdrop, control=link.control
                )
                if not all(
                    math.isfinite(value) for value in (link.flow0, link.flow1)
                ):
                    raise NumericalSolverError(
                        f'Element on link "{link.name}" produced a non-finite '
                        f'flow at nonlinear iteration {iteration}'
                    )
                link.flow_count = nf
        return sum_abs_flow

    def initialize(self, maxiter: int = 100):
        """Initialize model pressures using element linearizations."""
        self._require_current_weather()
        return self.model._initialize(self, maxiter=maxiter)

    def steady_solve(
        self,
        maxiter: int = 100,
        max_subiter: int = 100,
        tolerance: float = 1.0e-8,
        status_function=devnull,
        relative_tolerance: float | None = None,
        relaxation: float = 1.0,
        minimum_iterations: int = 1,
    ):
        """Solve the model's steady airflow equations."""
        self._require_current_weather()
        return self.model._air_movement(
            self,
            maxiter=maxiter,
            max_subiter=max_subiter,
            tolerance=tolerance,
            status_function=status_function,
            relative_tolerance=relative_tolerance,
            relaxation=relaxation,
            minimum_iterations=minimum_iterations,
        )

    def transient_step(self, *args, **kwargs):
        """Report that AIRNET has no transient airflow formulation."""
        raise NotImplementedError('AIRNET does not define a transient solve')

def write_results_csv(nodes, links, csv_file_name):
    """Write node and link results in airnet's CSV interchange format."""
    with open(csv_file_name, 'w', encoding='utf-8', newline='') as output:
        output.write('node header, name, time id, pressure, temperature, density\n')
        output.writelines(
            f'node, {node.name}, 0, {node.pressure:21.15e}, '
            f'{node.temperature:21.15e}, {node.density:21.15e}\n'
            for node in nodes
        )
        output.write('link header, name, time id, pressure drop, flow0, flow1\n')
        output.writelines(
            f'link, {link.name}, 0, {link.input_pressure_drop:21.15e}, '
            f'{link.input_flow0:21.15e}, {link.input_flow1:21.15e}\n'
            for link in links
        )


def run_simulate(
    input_file,
    verbose=False,
    global_temperature=None,
    global_density=None,
    *,
    output_file='airnetsim.csv',
    wind_file=None,
    ambient_temperature=20.0,
    barometric_pressure=101325.0,
    wind_speed=0.0,
    wind_direction=0.0,
    initialization_iterations=100,
    maximum_iterations=100,
    maximum_linear_iterations=100,
    absolute_tolerance=1.0e-8,
    relative_tolerance=None,
    relaxation=1.0,
):
    """Read, solve, and optionally write one AIRNET network simulation.

    This convenience API performs a single, non-interactive run. It returns
    the solved :class:`Model`; callers that need structured results can read
    the model's nodes and links directly. The output file is written only
    after successful initialization and nonlinear convergence.

    Parameters are expressed in AIRNET's SI units: degrees Celsius, pascals,
    metres per second, degrees clockwise from north, and kilograms per second.
    Set ``output_file`` to ``None`` to suppress CSV output.
    """
    if verbose:
        print(f'Opening input file "{input_file}"...')
    with open(input_file, encoding='utf-8') as input_stream:
        reader = Reader(input_stream)
        if verbose:
            print(f'Reading input file "{input_file}"...')
        items = list(reader)
    if verbose:
        print(f'Closing input file "{input_file}".')

    wind_profiles = None
    if wind_file is not None:
        if verbose:
            print(f'Reading wind profiles from "{wind_file}"...')
        with open(wind_file, encoding='utf-8') as wind_stream:
            wind_profiles = read_wind_profiles(wind_stream)

    weather = Weather(
        ambient_temperature=ambient_temperature,
        barometric_pressure=barometric_pressure,
        wind_speed=wind_speed,
        wind_direction=wind_direction,
        wind_profiles=wind_profiles or {},
    )
    model = Model(
        items,
        global_temperature=global_temperature,
        global_density=global_density,
        weather=weather,
    )
    if verbose:
        print(model.summary())

    solver = Solver(model)
    solver.initialize(maxiter=initialization_iterations)
    iterations = solver.steady_solve(
        maxiter=maximum_iterations,
        max_subiter=maximum_linear_iterations,
        tolerance=absolute_tolerance,
        relative_tolerance=relative_tolerance,
        relaxation=relaxation,
        status_function=print if verbose else devnull,
    )
    if verbose:
        print(f'iters = {iterations}')
        print(model.results_summary())

    if output_file is not None:
        write_results_csv(
            [node for node in model.nodes.values() if node.variable],
            model.links,
            output_file,
        )
    return model
