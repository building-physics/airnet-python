# SPDX-FileCopyrightText: 2023-present Oak Ridge National Laboratory, managed by UT-Battelle
#
# SPDX-License-Identifier: BSD-3-Clause
import argparse
import os
import math
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


class Link:
    def __init__(self, name=None, node0=None, ht0=0.0, node1=None, ht1=0.0, element=None,
                 wind=None, wpmod=0.0, mult=1.0):
        self.name = name
        self.node0 = node0
        self.node1 = node1
        self.ht0 = ht0
        self.ht1 = ht1
        self.element = element
        self.wind = wind
        self.wpmod = wpmod
        self.mult = mult

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
                variable = el.pop('type', 'c') == 'c'
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
            node0 = self.nodes[link.pop('node-1')]
            ht0 = link.pop('ht-1')
            node1 = self.nodes[link.pop('node-2')]
            ht1 = link.pop('ht-2')
            element = self.elements[link.pop('element')]
            self.links.append(link_object(node0=node0, ht0=ht0, node1=node1, ht1=ht1, element=element, **link))
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

    def set_properties(self):
        for node in self.variable_nodes:
            node.density = 0.0034838*(101325.0+node.pressure)/node.temperature
            node.sqrt_density = math.sqrt(node.density)
            node.viscosity = 1.71432e-5 + 4.828E-8 * (node.temperature - 273.15)
            node.dvisc = node.density / node.viscosity

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
