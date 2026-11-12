# SPDX-FileCopyrightText: 2023-present Oak Ridge National Laboratory, managed by UT-Battelle
#
# SPDX-License-Identifier: BSD-3-Clause
from types import SimpleNamespace

import pytest

from airnet.afedat import Afe_Dwc


def make_node(density, viscosity):
    return SimpleNamespace(density=density, viscosity=viscosity)


def test_dwc_linearization_uses_upwind_properties():
    element = Afe_Dwc(
        length=2.0,
        hdia=0.25,
        area=0.04909,
        rough=0.00015,
        tdlc=0.0,
        lflc=64.0,
        ldlc=0.0,
        linit=128.0,
    )
    link = SimpleNamespace(
        node0=make_node(density=1.2, viscosity=1.8e-5),
        node1=make_node(density=0.9, viscosity=2.0e-5),
    )

    positive = 2.0 * 1.2 * 0.04909 * 0.25 / (1.8e-5 * 128.0 * 8.0)
    negative = 2.0 * 0.9 * 0.04909 * 0.25 / (2.0e-5 * 128.0 * 8.0)

    assert element.linearize(link, 1.0) == pytest.approx(positive)
    assert element.linearize(link, -1.0) == pytest.approx(negative)


def test_dwc_laminar_jacobian_matches_finite_difference():
    element = Afe_Dwc(
        length=2.0,
        hdia=0.25,
        area=0.04909,
        rough=0.00015,
        tdlc=0.0,
        lflc=64.0,
        ldlc=0.0,
        linit=128.0,
    )
    link = SimpleNamespace(
        node0=make_node(density=1.2, viscosity=1.8e-5),
        node1=make_node(density=1.2, viscosity=1.8e-5),
    )
    pressure = 1.0e-8
    step = 1.0e-10

    _, flow, _, derivative, _ = element.jacobian(link, pressure)
    _, upper_flow, _ = element.flow(link, pressure + step)
    _, lower_flow, _ = element.flow(link, pressure - step)
    finite_difference = (upper_flow - lower_flow) / (2.0 * step)

    assert flow == pytest.approx(derivative * pressure)
    assert derivative == pytest.approx(finite_difference)


def test_dwc_turbulent_flow_matches_dw1_parallel_duct():
    element = Afe_Dwc(
        length=10.0,
        hdia=0.25,
        area=0.04909,
        rough=0.00015,
        tdlc=0.0,
        lflc=64.0,
        ldlc=0.0,
        linit=128.0,
    )
    temperature = 293.15
    density = 0.0034838 * (101325.0 + 9.0) / temperature
    viscosity = 1.71432e-5 + 4.828e-8 * 20.0
    link = SimpleNamespace(
        node0=make_node(density=density, viscosity=viscosity),
        node1=make_node(density=density, viscosity=viscosity),
    )

    # AIRNET retains the friction-factor estimate between nonlinear
    # evaluations; DW1 converges on its second reported iteration.
    element.jacobian(link, 9.0)
    _, flow, _, derivative, _ = element.jacobian(link, 9.0)

    assert flow == pytest.approx(0.2448779, rel=1.0e-6)
    assert derivative == pytest.approx(0.5 * flow / 9.0)
