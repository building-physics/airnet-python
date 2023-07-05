# SPDX-FileCopyrightText: 2023-present Oak Ridge National Laboratory, managed by UT-Battelle
#
# SPDX-License-Identifier: BSD-3-Clause
import airnet
from unittest import mock
import argparse

AFDATA_PL2 = '''/*subfile:  afdata.pl2  ******************************************************/
/
title  powerlaw test #2 input file

node   node-1    c   0.0   20.0   0.0
node   node-2    v   0.0   20.0
node   node-3    v   0.0   20.0
node   node-4    c   0.0   20.0  -100.0

element  orf-0.0001 plr 8.124e-09 8.124e-09 8.48528e-05 0.5 ! orf - 0.0001 m^2
element  orf-0.001  plr 2.569e-07 2.569e-07 0.000848528 0.5 ! orf - 0.001 m^2
element  orf-0.01   plr 8.124e-06 8.124e-06 0.00848528  0.5 ! orf - 0.01 m^2
element  orf-0.1    plr 0.0002569 0.0002569 0.0848528   0.5 ! orf - 0.1 m^2
element  orf-1.0    plr 0.008124  0.008124  0.848528    0.5 ! orf - 1.0 m^2
element  orf-10.0   plr 0.2569    0.2569    8.48528     0.5 ! orf - 10.0 m^2
element  orf-100.0  plr 8.124     8.124     84.8528     0.5 ! orf - 100.0 m^2

link   link-1   node-1   0.0   node-2   0.0   orf-0.0001   null
link   link-2   node-2   0.0   node-3   0.0   orf-0.0001   null
link   link-3   node-3   0.0   node-4   0.0   orf-0.0001   null

*********'''

def test_model_creation():
    reader = airnet.Reader(iter(AFDATA_PL2.splitlines()))
    items = []

    for item in reader:
        items.append(item)

    model = airnet.Model(items)
    assert len(model.nodes) == 4
    assert len(model.links) == 3
    assert len(model.elements) == 7