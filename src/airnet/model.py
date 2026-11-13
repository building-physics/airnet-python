# SPDX-FileCopyrightText: 2023-present Oak Ridge National Laboratory, managed by UT-Battelle
#
# SPDX-License-Identifier: BSD-3-Clause
from __future__ import annotations # Remove when dropping 3.9
import argparse
import os
import math
import scipy
import numpy
from .afedat import object_lookup
from .reader import Reader, InputType

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
                 pres:float=0.0, index:int|None=None, input_c:bool=True):
        self.name = name
        self.variable = variable
        self.height = ht
        self.temperature = temp
        if input_c:
            self.temperature += 273.15
        self.pressure = pres
        self.index = index
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
                 wind=None, wpmod=0.0, mult=1.0, flipped=False):
        self.name = name
        self.node0 = node0
        self.node1 = node1
        self.ht0 = ht0
        self.ht1 = ht1
        self.element = element
        self.wind = wind
        self.wpmod = wpmod
        self.multiplier = mult
        self.flipped = flipped
        self.flow_count = 0
        self.flow0 = 0.0
        self.flow1 = 0.0
        self.pdrop = 0.0

class BadNetwork(Exception):
    """Raised if the network is bad."""
    pass

class Model:
    """A class containing nodes, links, and other data representing a pressure network.
    """
    def __init__(self, network_input, element_lookup:dict = object_lookup, node_object=Node,
                 link_object=Link, global_temperature:float|None=None, global_density:float|None=None):
        self.title = ''
        self.nodes = {}
        self.links = []
        self.elements = {}
        # Handle the network inputs
        links = []
        for el in network_input:
            if el['input_type'] == InputType.TITLE:
                self.title = el['title']
            elif el['input_type'] == InputType.NODE:
                variable = el.pop('type', 'c') == 'v'
                el.pop('input_type')
                self.nodes[el['name']] = node_object(variable=variable, **el)
            elif el['input_type'] == InputType.LINK:
                el.pop('input_type')
                links.append(el)
            elif el['input_type'] == InputType.ELEMENT:
                el.pop('input_type')
                type = el.pop('type')
                name = el.pop('name')
                self.elements[name] = element_lookup[type](**el)
        for link in links:
            flip = False
            node0 = self.nodes[link.pop('node-1')]
            ht0 = link.pop('ht-1')
            node1 = self.nodes[link.pop('node-2')]
            ht1 = link.pop('ht-2')
            if not node0.variable:
                if node1.variable:
                    flip = True
            element = self.elements[link.pop('element')]
            if flip:
                self.links.append(link_object(node0=node1, ht0=ht1, node1=node0, ht1=ht0,
                                              element=element, **link, flipped=flip))
            else:
                self.links.append(link_object(node0=node0, ht0=ht0, node1=node1, ht1=ht1,
                                              element=element, **link, flipped=flip))
        # Figure out the size of the matrix
        self.variable_nodes = []
        count = 0
        for node in self.nodes.values():
            if node.variable:
                node.index = count
                count += 1
                self.variable_nodes.append(node)
        assert count == len(self.variable_nodes)
        self.size = count
        row = []
        col = []
        data = []
        for link in self.links:
            if link.node0.variable:
                # diagonal term
                row.append(link.node0.index)
                col.append(link.node0.index)
                data.append(1.0)
                if link.node1.variable:
                    # diagonal term
                    row.append(link.node1.index)
                    col.append(link.node1.index)
                    data.append(1.0)
                    # off diagonal terms
                    row.append(link.node0.index)
                    col.append(link.node1.index)
                    data.append(1.0)
                    # off diagonal terms
                    row.append(link.node1.index)
                    col.append(link.node0.index)
                    data.append(1.0)

        matrix = scipy.sparse.coo_matrix((numpy.array(data, dtype=numpy.double),
                                          (numpy.array(row), numpy.array(col))),
                                          shape=(count, count))
        self.A = scipy.sparse.csr_matrix(matrix)
        self.x = numpy.zeros([self.size, 1], dtype=numpy.double)

        self.set_properties(self.nodes.values(), global_temperature=global_temperature,
                            global_density=global_density)

        # Check for disconnected nodes
        problems = []
        for i,el in enumerate(self.A.diagonal()):
            if el < 1:
                for node in self.nodes.values():
                    if node.index == i:
                        problems.append(node)
                        break
        if problems:
           raise BadNetwork('Disconnected nodes found: %s' % ', '.join([el.name for el in problems]))

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

        string += '\nNodes: %s\n\nLinks: %s\n' % (len(self.nodes), len(self.links))
        string += '\nSystem size: %d x %x\n' % (len(self.variable_nodes), len(self.variable_nodes))
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
    
        string += '\nNodes: %s\n\nLinks: %s\n' % (len(self.nodes), len(self.links))
        string += '\nSystem size: %d x %x\n' % (len(self.variable_nodes), len(self.variable_nodes))
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
        if global_temperature is None:
            if global_density is None:
                # Compute the full set of properties
                for node in nodes:
                    node.density = 0.0034838*(101325.0+node.pressure)/node.temperature
                    node.sqrt_density = math.sqrt(node.density)
                    node.viscosity = 1.71432e-5 + 4.828E-8 * (node.temperature - 273.15)
                    node.dvisc = node.density / node.viscosity
            else:
                # Set the density everywhere
                for node in nodes:
                    node.density = global_density
                    node.sqrt_density = math.sqrt(node.density)
                    node.viscosity = 1.71432e-5 + 4.828E-8 * (node.temperature - 273.15)
                    node.dvisc = node.density / node.viscosity
        else:
            if global_density is None:
                # Set the temperature everywhere
                for node in nodes:
                    node.temperature = global_temperature
                    node.density = 0.0034838*(101325.0+node.pressure)/node.temperature
                    node.sqrt_density = math.sqrt(node.density)
                    node.viscosity = 1.71432e-5 + 4.828E-8 * (node.temperature - 273.15)
                    node.dvisc = node.density / node.viscosity
            else:
                # Set the temperature and density everywhere
                for node in nodes:
                    node.temperature = global_temperature
                    node.density = global_density
                    node.sqrt_density = math.sqrt(node.density)
                    node.viscosity = 1.71432e-5 + 4.828E-8 * (node.temperature - 273.15)
                    node.dvisc = node.density / node.viscosity

    def initialize(self, maxiter:int=100):
        """Initialize the flow network.
        
        Compute flows and pressure drops using linear flow representation for all elements.
        
        Parameters
        ----------
        maxiter: optional
            The maximum number of iterations allowed in the solution.
            
        Returns
        -------
        bool:
            True is returned if the initialization has converged in the allowed number of iterations, False otherwise.
        """
        self.A.data.fill(0.0)
        self.x.fill(0.0)
        for link in self.links:
            if link.node0.variable:
                pdrop = link.node0.pressure - link.node1.pressure
                c = link.element.linearize(link, pdrop=pdrop)
                # diagonal term
                self.A[link.node0.index, link.node0.index] += c
                if link.node1.variable:
                    # diagonal term
                    self.A[link.node1.index, link.node1.index] += c
                    # off diagonal terms
                    self.A[link.node0.index, link.node1.index] -= c
                    self.A[link.node1.index, link.node0.index] -= c
                else:
                    self.x[link.node0.index] += c*link.node1.pressure
        self.x, info = scipy.sparse.linalg.cg(
            self.A, self.x, rtol=1.0e-12, atol=0.0, maxiter=maxiter
        )
        if info == 0:
            # Update the nodal pressures
            for node in self.variable_nodes:
                node.pressure = self.x[node.index]
            # Update the flows:
            for link in self.links:
                pdrop = link.node0.pressure - link.node1.pressure
                c = link.element.linearize(link, pdrop=pdrop)
                link.flow_count = 1
                link.flow0 = c*(link.node0.pressure-link.node1.pressure)
                link.flow1 = 0.0
            return True
        return False
    
    def compute_pressure_drops(self):
        """Compute the pressure drop across all the links in the model."""
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
                    spx = 4.90 * (link.node0.dens + link.node1.dens) * dhx
            # Wind pressure contribution goes here
            link.pdrop =  sp0 + sp1 + spx

    def air_movement(self, maxiter:int=100, max_subiter:int=100,
                     tolerance:float=1.0e-8, status_function=devnull,
                     relative_tolerance:float|None=None, relaxation:float=1.0,
                     minimum_iterations:int=1):
        """Compute steady airflows in the model.
        
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
        RuntimeError
            If either a linear correction fails or the nonlinear iteration
            does not converge within ``maxiter`` iterations.
        """
        if tolerance < 0.0:
            raise ValueError('Absolute convergence tolerance cannot be negative')
        if relative_tolerance is not None and relative_tolerance < 0.0:
            raise ValueError('Relative convergence tolerance cannot be negative')
        if not 0.0 < relaxation <= 1.0:
            raise ValueError('Relaxation coefficient must be greater than 0 and at most 1')
        if minimum_iterations < 1:
            raise ValueError('Minimum iterations must be at least 1')

        self.compute_pressure_drops()
        status_function('iter | Max Resid|\n==== ===============')
        previous_relative_residual = 0.0
        for iter in range(1,maxiter+1):
            self.A.data.fill(0.0)
            self.x.fill(0.0)
            sum_abs_flow = numpy.zeros(self.size, dtype=numpy.double)
            for link in self.links:
                pdrop = link.node0.pressure - link.node1.pressure + link.pdrop
                if link.node0.variable:
                    nf, link.flow0, link.flow1, df0, df1 = link.element.jacobian(link, pdrop)
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
                    nf, link.flow0, link.flow1 = link.element.flow(link, pdrop)
                    link.flow_count = nf

            residual = numpy.asarray(self.x).reshape(-1)
            absolute_residual = numpy.abs(residual)
            maxf = float(absolute_residual.max(initial=0.0))
            converged = absolute_residual < tolerance
            if relative_tolerance is not None:
                relative_by_node = numpy.full(self.size, numpy.inf)
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
                # Update the pressure drops, up until now only secondary terms present
                for link in self.links:
                    link.pdrop += link.node0.pressure - link.node1.pressure
                return iter

            status_function('%4d %15.9e' % (iter, maxf))
            correction, info = scipy.sparse.linalg.cg(
                self.A, residual, rtol=1.0e-12, atol=0.0, maxiter=max_subiter
            )
            if info == 0:
                correction_factor = 1.0
                if relative_residual > 0.30 * previous_relative_residual:
                    correction_factor = relaxation
                previous_relative_residual = relative_residual

                # Update the nodal pressures; flows are set above.
                for node in self.variable_nodes:
                    node.pressure -= correction_factor * correction[node.index]
            else:
                raise RuntimeError('Newton iteration solve failed at %d iterations' % iter)
        raise RuntimeError('Airflow solve failed to converge in %d iterations' % maxiter)

def write_results_csv(nodes, links, csv_file_name: str)   :
    fp = open(csv_file_name, 'w')
    fp.write('node header, name, time id, pressure, temperature, density\n')
    for node in nodes:
        fp.write('node, %s, 0, %21.15e, %21.15e, %21.15e\n' % (node.name, node.pressure,
                                                               node.temperature, node.density))
    fp.write('link header, name, time id, pressure drop, flow0, flow1\n')
    for link in links:
        fp.write('link, %s, 0, %21.15e, %21.15e, %21.15e\n' % (link.name, link.pdrop, link.flow0, link.flow1))
    fp.close()

def run_simulate(input_file, verbose=False, global_temperature=None, global_density=None):

    if verbose:
        print('Opening input file "%s"...' % input_file)

    if not os.path.exists(input_file):
        print('airnetsim: error: the input file "%s" does not exist' % input_file)
        return 1

    fp = open(input_file, 'r')
    reader = Reader(fp)
    if verbose:
        print('Reading input file "%s"...' % input_file)
    items = []
    for item in reader:
        items.append(item)
    if verbose:
        print('Closing input file "%s".' % input_file)
    fp.close()

    model = Model(items, global_temperature=global_temperature, global_density=global_density)
    
    if verbose:
        print(model.summary())

    model.initialize()

    if verbose:
        iters = model.air_movement(status_function=print)
        print('iters = %d' % iters)
        print(model.results_summary())
    else:
        model.air_movement()

    write_results_csv([el for el in model.nodes.values() if el.index is not None], model.links, 'airnetsim.csv')

