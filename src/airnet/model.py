# SPDX-FileCopyrightText: 2023-present Oak Ridge National Laboratory, managed by UT-Battelle
#
# SPDX-License-Identifier: BSD-3-Clause
import argparse
import os
import math
import scipy
import numpy
from .afedat import object_lookup
from .reader import Reader, InputType

class Node:
    def __init__(self, name=None, variable=True, ht=0.0, temp=293.15, pres=0.0,
                 index=None):
        self.name = name
        self.variable = variable
        self.height = ht
        self.temperature = temp
        self.pressure = pres
        self.index = index
        self.density = 0.0
        self.viscosity = 0.0
        self.sqrt_density = 0.0
        self.dvisc = 0.0 # Density divided by viscosity
    def copy_state(self, other):
        self.temperature = other.temperature
        self.pressure = other.pressure
        self.density = other.density
        self.viscosity = other.viscosity
        self.sqrt_density = other.sqrt_density
        self.dvisc = other.dvisc


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
        self.flow0 = 0.0
        self.flow1 = 0.0
        self.pdrop = 0.0

class BadNetwork(Exception):
    pass

class Model:
    def __init__(self, network_input, element_lookup = object_lookup, node_object=Node,
                 link_object=Link):
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

        self.set_properties()

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

    def set_variable_properties(self):
        for node in self.variable_nodes:
            node.density = 0.0034838*(101325.0+node.pressure)/node.temperature
            node.sqrt_density = math.sqrt(node.density)
            node.viscosity = 1.71432e-5 + 4.828E-8 * (node.temperature - 273.15)
            node.dvisc = node.density / node.viscosity

    def set_properties(self):
        for node in self.nodes.values():
            node.density = 0.0034838*(101325.0+node.pressure)/node.temperature
            node.sqrt_density = math.sqrt(node.density)
            node.viscosity = 1.71432e-5 + 4.828E-8 * (node.temperature - 273.15)
            node.dvisc = node.density / node.viscosity

    def initialize(self, maxiter=100):
        self.A.data.fill(0.0)
        self.x.fill(0.0)
        for link in self.links:
            if link.node0.variable:
                c = link.element.linearize(link)
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
        self.x, info = scipy.sparse.linalg.cg(self.A, self.x, maxiter=maxiter)
        if info == 0:
            # Update the nodal pressures
            for node in self.variable_nodes:
                node.pressure = self.x[node.index]
            # Update the flows:
            for link in self.links:
                c = link.element.linearize(link)
                link.flow0 = c*(link.node0.pressure-link.node1.pressure)
                link.flow1 = 0.0
            return True
        return False
    
    def compute_pressure_drops(self):
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

    def air_movement(self, maxiter=100, max_subiter=100, tolerance=1.0e-8):
        self.compute_pressure_drops()
        for iter in range(1,maxiter+1):
            self.A.data.fill(0.0)
            self.x.fill(0.0)
            for link in self.links:
                if link.node0.variable:
                    pdrop = link.node0.pressure - link.node1.pressure + link.pdrop
                    nf, link.flow0, link.flow1, df0, df1 = link.element.jacobian(link, pdrop)
                    print(link.name, nf, link.flow0, link.flow1, df0, df1)
                    if nf == 1:
                        # diagonal term
                        self.A[link.node0.index, link.node0.index] += df0
                        self.x[link.node0.index] += link.flow0
                        if link.node1.variable:
                            # diagonal term
                            self.A[link.node1.index, link.node1.index] += df0
                            self.x[link.node1.index] -= link.flow0
                            # off diagonal terms
                            self.A[link.node0.index, link.node1.index] -= df0
                            self.A[link.node1.index, link.node0.index] -= df0
                    else:
                        pass # Do this later...
            maxf = max(self.x, key=abs)
            if maxf <= tolerance:
                return iter
            self.x, info = scipy.sparse.linalg.cg(self.A, self.x, maxiter=max_subiter)
            info = 0
            if info == 0:
                # Update the nodal pressures, flows are set above
                for node in self.variable_nodes:
                    node.pressure -= self.x[node.index]
            else:
                raise 'STOPSTOPSTOP'
        return maxiter
            

def summarize_input():
    parser = argparse.ArgumentParser(description='Summarize an AIRNET network input file.')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true',
                        default=False, help='operate verbosely')
    #parser.add_argument('-d', '--debug', dest='debug', action='store_true',
    #                    default=False, help='produce debug output')
    parser.add_argument('input', metavar='NETWORK_FILE')

    args = parser.parse_args()

    if args.verbose:
        print('Opening input file "%s"...' % args.input)

    if not os.path.exists(args.input):
        print('summarize_airnet_input: error: the input file "%s" does not exist' % args.input)
        return 1

    fp = open(args.input, 'r')
    reader = Reader(fp)
    if args.verbose:
        print('Reading input file "%s"...' % args.input)
    items = []
    for item in reader:
        items.append(item)
    if args.verbose:
        print('Closing input file "%s".' % args.input)
    fp.close()

    if reader.title:
        print('Title:', reader.title)

    print('\nElements:\n=========')
    elements = {}
    nodes = 0
    links = 0
    for el in items:
        if el['input_type'] == InputType.ELEMENT:
            if el['type'] in elements:
                elements[el['type']] += 1
            else:
                elements[el['type']] = 1
        elif el['input_type'] == InputType.NODE:
            nodes += 1
        elif el['input_type'] == InputType.LINK:
            links += 1
    for key, value in elements.items():
        print(key, value)

    print('\nNodes: %s\n\nLinks: %s' % (nodes, links))

    if args.verbose:
        print('Done.')
    return 0

def simulate():
    parser = argparse.ArgumentParser(description='Simulate airflow in an AIRNET network.')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true',
                        default=False, help='operate verbosely')
    #parser.add_argument('-d', '--debug', dest='debug', action='store_true',
    #                    default=False, help='produce debug output')
    parser.add_argument('input', metavar='NETWORK_FILE')

    args = parser.parse_args()

    if args.verbose:
        print('Opening input file "%s"...' % args.input)

    if not os.path.exists(args.input):
        print('airnetsim: error: the input file "%s" does not exist' % args.input)
        return 1

    fp = open(args.input, 'r')
    reader = Reader(fp)
    if args.verbose:
        print('Reading input file "%s"...' % args.input)
    items = []
    for item in reader:
        items.append(item)
    if args.verbose:
        print('Closing input file "%s".' % args.input)
    fp.close()

    model = Model(items)
    
    if args.verbose:
        print(model.summary())

def gui_simulate():
    pass
