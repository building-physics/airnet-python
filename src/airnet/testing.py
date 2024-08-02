# SPDX-FileCopyrightText: 2023-present Oak Ridge National Laboratory, managed by UT-Battelle
#
# SPDX-License-Identifier: BSD-3-Clause
import csv
from contextlib import contextmanager
import tempfile
import os

@contextmanager
def temporary_directory():
    origin = os.getcwd()
    try:
        tmp = tempfile.TemporaryDirectory()
        os.chdir(tmp.name)
        yield
    finally:
        os.chdir(origin)
        tmp.cleanup()

def compare_nodes(line1, line2, pressure_tolerance, temperature_tolerance, density_tolerance, line_number):
    messages = []
    # Check names
    if line1[1].strip() != line2[1].strip():
        messages.append('Nodes at line %d have different names: %s, %s' % (line_number, line1[1], line2[1]))
    # Check the time id
    if int(line1[2]) != int(line2[2]):
        messages.append('Nodes at line %d have difference time index: %s, %s' % (line_number, line1[2], line2[2]))
    # Check pressure, temperature, and density
    index = [3, 4, 5]
    vars = ['pressure', 'temperature', 'density']
    tols = [pressure_tolerance, 1.0e-5, density_tolerance]
    for i,v,t in zip(index,vars,tols):
        delta = abs(float(line1[i]) - float(line2[i]))
        if delta > t:
            messages.append('Nodes at line %d have %s difference > %e: %s, %s' % (line_number, v, t, line1[i], line2[i]))

    return messages

def compare_csvs(csv1, csv2, node_pressure_tolerance=1.0e-7, node_temperature_tolerance=1.0e-5,
                 node_density_tolerance=1.0e-7, node_tolerance=None, link_tolerance=1.0e-8):
    messages = []
    if node_tolerance is not None:
        node_pressure_tolerance = node_tolerance
        node_temperature_tolerance = node_tolerance
        node_density_tolerance = node_tolerance
    with open(csv1) as fp1, open(csv2) as fp2:
        r1 = csv.reader(fp1)
        r2 = csv.reader(fp2)
        line_count = 1
        for line1, line2 in zip(r1, r2):
            # The first column is supposed to be the type of the row
            if line1[0] != line2[0]:
                messages.append('Line %d type mismatch (%s vs %s)' % (line_count, line1[0], line2[0]))
            else:
                if line1[0].endswith('header'):
                    pass
                elif line1[0] == 'node':
                    messages.extend(compare_nodes(line1, line2, node_pressure_tolerance, 
                                                  node_temperature_tolerance, node_density_tolerance, line_count))
                elif line1[0] == 'link':
                    pass
                else:
                    messages.append('Unknown row type: %s' % line1[0])
            line_count += 1
    return messages
