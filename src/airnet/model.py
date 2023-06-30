# SPDX-FileCopyrightText: 2023-present Oak Ridge National Laboratory, managed by UT-Battelle
#
# SPDX-License-Identifier: BSD-3-Clause
import airnet
import argparse
import os

def summarize_input():
    parser = argparse.ArgumentParser(description='Summarize an AIRNET network input file.')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true',
                        default=False, help='operate verbosely')
    #parser.add_argument('-d', '--debug', dest='debug', action='store_true',
    #                    default=False, help='produce debug output')
    parser.add_argument('input', metavar='INPUT_FILE')

    args = parser.parse_args()

    if args.verbose:
        print('Opening input file "%s"...' % args.input)

    if not os.path.exists(args.input):
        print('summarize_airnet_input: error: the input file "%s" does not exist' % args.input)
        return 1

    fp = open(args.input, 'r')
    reader = airnet.Reader(fp)
    if args.verbose:
        print('Reading input file "%s"...' % args.input)
    items = []
    for item in reader:
        items.append(item)
    if args.verbose:
        print('Closing input file "%s"...' % args.input)
    fp.close()

    if reader.title:
        print('Title:', reader.title)

    print('\nElements:\n=========')
    elements = {}
    nodes = 0
    links = 0
    for el in items:
        if el['input_type'] == airnet.InputType.ELEMENT:
            if el['type'] in elements:
                elements[el['type']] += 1
            else:
                elements[el['type']] = 1
        elif el['input_type'] == airnet.InputType.NODE:
            nodes += 1
        elif el['input_type'] == airnet.InputType.LINK:
            links += 1
    for value in airnet.ElementType:
        if value in elements:
            print(value, elements[value])

    print('\nNodes: %s\n\nLinks: %s' % (nodes, links))

    if args.verbose:
        print('Done.')
    return 0

def simulate():
    pass

def gui_simulate():
    pass
