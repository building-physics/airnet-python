# SPDX-FileCopyrightText: 2023-present Oak Ridge National Laboratory, managed by UT-Battelle
#
# SPDX-License-Identifier: BSD-3-Clause
import airnet
import math

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

class AltModel(airnet.Model):
    def set_variable_properties(self):
        pass

def test_model_creation():
    reader = airnet.Reader(iter(AFDATA_PL2.splitlines()))
    items = []

    for item in reader:
        items.append(item)

    model = AltModel(items)
    assert len(model.nodes) == 4
    assert len(model.variable_nodes) == 2
    assert len(model.links) == 3
    assert len(model.elements) == 7

    for el in model.elements.values():
        assert el.type() == 'plr'

    assert not model.nodes['node-1'].variable
    assert model.nodes['node-2'].variable
    assert model.nodes['node-3'].variable
    assert not model.nodes['node-4'].variable

    # Mess with the properties so that the properties are the same so 
    # we can do simple tests.
    node_4_pressure = model.nodes['node-4'].pressure
    node_1 = model.nodes['node-1']
    for node in model.nodes.values():
        node.copy_state(node_1)
    model.nodes['node-4'].pressure = node_4_pressure

    assert model.initialize(maxiter=10)

    # Compute the flow we should get (playing a bit loose with the state props involved)
    # f = C1(p1-p2) = C2(p2-p3) = C1(p3-p4) = Ce(p1-p4) = Ce((p1-p2) + (p2-p3) + (p3-p4))
    # f = Ce*f(1/C1 + 1/C2 + 1/C3)
    # 1/Ce = 1/C1 + 1/C2 + 1/C3
    Ce = model.nodes['node-1'].dvisc/(1.0/model.links[0].element.init + 1.0/model.links[1].element.init + 1.0/model.links[2].element.init)
    f = Ce * 100.0
    assert abs(model.links[0].flow0 + f) < 1.0e-12
    assert abs(model.links[0].flow0 + model.links[1].flow0) < 1.0e-12
    assert abs(model.links[1].flow0 - model.links[2].flow0) < 1.0e-12

    # Compute the pressure drops (playing a bit loose with the state props involved)
    # f = C1(p1-p2) => p2 = p1 - f/C1
    # f = C2(p2-p3) => p3 = p2 - f/C2
    assert model.nodes['node-1'].pressure == 0.0
    assert abs(model.nodes['node-2'].pressure - model.nodes['node-1'].pressure + f/(model.links[0].element.init * model.nodes['node-1'].dvisc)) < 1.0e-12
    assert abs(model.nodes['node-3'].pressure - model.nodes['node-2'].pressure + f/(model.links[1].element.init * model.nodes['node-1'].dvisc)) < 1.0e-12
    assert model.nodes['node-4'].pressure == -100.0
    assert model.links[0].flipped
    

    model.set_variable_properties()
    iters = model.air_movement(maxiter=50)
    assert iters < 5

    # Compute the flow we should get (playing a bit loose with the state props involved)
    # f = C1(p1-p2)^n = C2(p2-p3)^n = C1(p3-p4)^n = Ce(p1-p4)^n
    # (f/C1)^(1/n) = p1-p2
    # (f/C2)^(1/n) = p2-p3
    # (f/C3)^(1/n) = p3-p4
    # p1-p4 = (f/C1)^(1/n) + (f/C2)^(1/n) + (f/C3)^(1/n) = f^(1/n) ((1/C1)^(1/n) + (1/C2)^(1/n) + (1/C3)^(1/n))
    # (p1-p4)^n = f((1/C1)^(1/n) + (1/C2)^(1/n) + (1/C3)^(1/n))^n
    # f = ((1/C1)^(1/n) + (1/C2)^(1/n) + (1/C3)^(1/n))^(-n) (p1-p4)^n
    # Ce = ((1/C1)^(1/n) + (1/C2)^(1/n) + (1/C3)^(1/n))^(-n)
    Ce = 1.0/math.pow(math.pow(1.0/model.links[0].element.turb, 1.0/model.links[0].element.expt) +
                      math.pow(1.0/model.links[1].element.turb, 1.0/model.links[0].element.expt) +
                      math.pow(1.0/model.links[2].element.turb, 1.0/model.links[0].element.expt), model.links[0].element.expt)
    f = model.nodes['node-1'].sqrt_density * Ce * 10.0
    assert abs(model.links[0].flow0 + f) < 1.0e-12
    assert abs(model.links[0].flow0 + model.links[1].flow0) < 1.0e-12
    assert abs(model.links[1].flow0 - model.links[2].flow0) < 1.0e-12

    # Compute the pressure drops (playing a bit loose with the state props involved)
    # f = C1(p1-p2)^n => p2 = p1 - (f/C1)^(1/n)
    # f = C2(p2-p3)^n => p3 = p2 - (f/C2)^(1/n)
    assert model.nodes['node-1'].pressure == 0.0
    assert abs(model.nodes['node-2'].pressure - model.nodes['node-1'].pressure + math.pow(f/(model.links[0].element.turb * model.nodes['node-1'].sqrt_density), 1.0/model.links[0].element.expt)) < 1.0e-12
    assert abs(model.nodes['node-3'].pressure - model.nodes['node-2'].pressure + math.pow(f/(model.links[1].element.turb * model.nodes['node-1'].sqrt_density), 1.0/model.links[0].element.expt)) < 1.0e-12
    assert model.nodes['node-4'].pressure == -100.0