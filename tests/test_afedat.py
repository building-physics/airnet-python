# SPDX-FileCopyrightText: 2023-present Oak Ridge National Laboratory, managed by UT-Battelle
#
# SPDX-License-Identifier: BSD-3-Clause
from types import SimpleNamespace

import pytest

from airnet.afedat import (
    Afe_Cfr,
    Afe_Ckv,
    Afe_Cpf,
    Afe_Dor,
    Afe_Dwc,
    Afe_Fan,
    Afe_Plr,
    Afe_Qfr,
)


def make_node(density, viscosity, temperature=293.15):
    return SimpleNamespace(
        density=density,
        viscosity=viscosity,
        temperature=temperature,
        sqrt_density=density**0.5,
        dvisc=density / viscosity,
    )


@pytest.mark.parametrize("pressure", [-1.0, 0.0, 1.0])
def test_plr_one_way_return_contract(pressure):
    element = Afe_Plr(init=0.01, lam=0.01, turb=0.1, expt=0.5)
    link = SimpleNamespace(
        node0=make_node(density=1.2, viscosity=1.8e-5),
        node1=make_node(density=0.9, viscosity=2.0e-5),
        multiplier=1.0,
    )

    count, flow0, flow1, derivative0, derivative1 = element.calculate(
        link, pressure
    )

    assert count == 1
    assert flow1 == 0.0
    assert derivative1 == 0.0
    assert element.flow(link, pressure) == pytest.approx(
        (count, flow0, flow1)
    )
    assert derivative0 > 0.0


@pytest.mark.parametrize(
    ('pressure', 'upwind_node'),
    [(1.5, 0), (-1.5, 1), (0.0, 1)],
)
def test_plr_linearize_returns_flow_and_derivative(
    pressure, upwind_node
):
    element = Afe_Plr(init=0.01)
    nodes = (
        make_node(density=1.2, viscosity=1.8e-5),
        make_node(density=0.9, viscosity=2.0e-5),
    )
    link = SimpleNamespace(node0=nodes[0], node1=nodes[1], multiplier=1.0)
    derivative = element.init * nodes[upwind_node].dvisc

    assert element.linearize(link, pressure) == pytest.approx(
        (derivative * pressure, derivative)
    )


def test_passive_element_multiplier_scales_flow_and_derivative():
    element = Afe_Plr(init=0.01, lam=0.01, turb=0.1, expt=0.5)
    link = SimpleNamespace(
        node0=make_node(density=1.2, viscosity=1.8e-5),
        node1=make_node(density=0.9, viscosity=2.0e-5),
        multiplier=1.0,
    )
    base_linear = element.linearize(link, 1.5)
    base_result = element.calculate(link, 1.5)

    link.multiplier = 2.5

    assert element.linearize(link, 1.5) == pytest.approx(
        tuple(2.5 * value for value in base_linear)
    )
    scaled_result = element.calculate(link, 1.5)
    assert scaled_result[1] == pytest.approx(2.5 * base_result[1])
    assert scaled_result[3] == pytest.approx(2.5 * base_result[3])


def test_dwc_linearize_uses_upwind_properties():
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
        multiplier=1.0,
    )

    positive = 2.0 * 1.2 * 0.04909 * 0.25 / (1.8e-5 * 128.0 * 8.0)
    negative = 2.0 * 0.9 * 0.04909 * 0.25 / (2.0e-5 * 128.0 * 8.0)

    assert element.linearize(link, 1.0) == pytest.approx(
        (positive, positive)
    )
    assert element.linearize(link, -1.0) == pytest.approx(
        (-negative, negative)
    )


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
        multiplier=1.0,
    )
    pressure = 1.0e-8
    step = 1.0e-10

    count, flow, flow1, derivative, derivative1 = element.calculate(
        link, pressure
    )
    _, upper_flow, _ = element.flow(link, pressure + step)
    _, lower_flow, _ = element.flow(link, pressure - step)
    finite_difference = (upper_flow - lower_flow) / (2.0 * step)

    assert count == 1
    assert flow1 == 0.0
    assert derivative1 == 0.0
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
        multiplier=1.0,
    )

    # AIRNET retains the friction-factor estimate between nonlinear
    # evaluations; DW1 converges on its second reported iteration.
    element.calculate(link, 9.0)
    _, flow, _, derivative, _ = element.calculate(link, 9.0)

    assert flow == pytest.approx(0.2448779, rel=1.0e-6)
    assert derivative == pytest.approx(0.5 * flow / 9.0)


@pytest.mark.parametrize('pressure', [-10.0, 0.0, 10.0])
def test_qfr_flow_satisfies_quadratic_relation(pressure):
    element = Afe_Qfr(a=0.31, b=44.8)
    link = SimpleNamespace(multiplier=1.0)

    count, flow, flow1, derivative, derivative1 = element.calculate(
        link, pressure
    )

    reconstructed_pressure = element.a * flow
    reconstructed_pressure += element.b * flow * abs(flow)
    assert count == 1
    assert flow1 == 0.0
    assert derivative1 == 0.0
    assert reconstructed_pressure == pytest.approx(pressure)
    assert derivative == pytest.approx(
        1.0 / (element.a + 2.0 * element.b * abs(flow))
    )
    assert element.flow(link, pressure) == pytest.approx(
        (count, flow, flow1)
    )


@pytest.mark.parametrize('pressure', [-10.0, -1.0e-8, 0.0, 1.0e-8, 10.0])
def test_qfr_derivative_matches_finite_difference(pressure):
    element = Afe_Qfr(a=0.31, b=44.8)
    link = SimpleNamespace(multiplier=1.0)
    step = 1.0e-10 if abs(pressure) < 1.0e-7 else 1.0e-6

    _, _, _, derivative, _ = element.calculate(link, pressure)
    _, upper, _ = element.flow(link, pressure + step)
    _, lower, _ = element.flow(link, pressure - step)

    assert derivative == pytest.approx(
        (upper - lower) / (2.0 * step), rel=1.0e-6
    )


def test_qfr_linearize_returns_flow_and_derivative():
    element = Afe_Qfr(a=0.25, b=3.0)

    assert element.linearize(
        SimpleNamespace(multiplier=1.0), -2.0
    ) == (-8.0, 4.0)


def test_cfr_flow_is_pressure_independent():
    element = Afe_Cfr(flow=1.25)
    link = SimpleNamespace(multiplier=1.0)

    assert element.linearize(link, 100.0) == (1.25, 0.0)
    assert element.linearize(link, -100.0) == (1.25, 0.0)
    assert element.calculate(link, -100.0) == (1, 1.25, 0.0, 0.0, 0.0)
    assert element.calculate(link, 100.0) == (1, 1.25, 0.0, 0.0, 0.0)
    link.multiplier = -1.0
    assert element.linearize(link, 100.0) == (
        -1.25,
        0.0,
    )
    assert element.calculate(link, 100.0) == (
        1,
        -1.25,
        0.0,
        0.0,
        0.0,
    )


def make_fan():
    return Afe_Fan(
        init=8.1216e-6,
        lam=8.1216e-6,
        turb=0.00848528,
        expt=0.5,
        rdens=1.204,
        fdf=764.429,
        sop=5.46267,
        off=0.02,
        mf1=10000.0,
        pts=[
            {
                'a1': 764.429,
                'a2': -18.2932,
                'a3': 19.4639,
                'a4': -7.63948,
                'mf': -10000.0,
            }
        ],
    )


def test_fan_curve_flow_and_derivative_match_finite_difference():
    fan = make_fan()
    link = SimpleNamespace(
        node0=make_node(1.204, 1.8e-5),
        node1=make_node(1.204, 1.8e-5),
        multiplier=1.0,
    )
    pressure = -660.3155
    # The recovered AIRNET calculation rounds its Cardano temporaries to
    # single precision, so use a pressure step large enough to rise above
    # those quantization jumps.
    step = 0.1

    count, flow, flow1, derivative, derivative1 = fan.calculate(
        link, pressure, control=1.0
    )
    _, upper, _ = fan.flow(link, pressure + step, control=1.0)
    _, lower, _ = fan.flow(link, pressure - step, control=1.0)

    assert count == 1
    assert flow == pytest.approx(3.156744, rel=1.0e-5)
    assert flow1 == 0.0
    assert derivative1 == 0.0
    assert derivative == pytest.approx(
        (upper - lower) / (2.0 * step), rel=0.02
    )


def test_fan_below_off_control_uses_power_law_relation():
    fan = make_fan()
    link = SimpleNamespace(
        node0=make_node(1.204, 1.8e-5),
        node1=make_node(1.204, 1.8e-5),
        multiplier=1.0,
    )

    assert fan.calculate(link, 4.0, control=0.01) == Afe_Plr.calculate(
        fan, link, 4.0, control=0.01
    )


def test_fan_multiplier_sign_orients_and_magnitude_scales_curve():
    fan = make_fan()
    node0 = make_node(1.204, 1.8e-5)
    node1 = make_node(1.204, 1.8e-5)
    forward_link = SimpleNamespace(
        node0=node0, node1=node1, multiplier=1.0
    )
    reversed_link = SimpleNamespace(
        node0=node1, node1=node0, multiplier=-2.0
    )

    forward = fan.calculate(forward_link, -660.0)
    reversed_scaled = fan.calculate(reversed_link, 660.0)

    assert reversed_scaled[1] == pytest.approx(-2.0 * forward[1])
    assert reversed_scaled[3] == pytest.approx(2.0 * forward[3])


def test_cpf_flow_and_derivative_match_finite_difference():
    cpf = Afe_Cpf(upo=10.0, prmin=1.0, ftyp=2.0)
    link = SimpleNamespace(
        node0=make_node(1.204, 1.8e-5),
        node1=make_node(1.204, 1.8e-5),
        multiplier=1.0,
    )
    pressure = -6.0
    step = 1.0e-5

    assert cpf.linearize(link, pressure) == pytest.approx((2.0, 1.0 / 3.01))
    _, flow, _, derivative, _ = cpf.calculate(link, pressure)
    _, upper, _ = cpf.flow(link, pressure + step)
    _, lower, _ = cpf.flow(link, pressure - step)

    assert flow == pytest.approx(2.0066666666666664)
    assert derivative == pytest.approx(
        (upper - lower) / (2.0 * step), rel=1.0e-6
    )


@pytest.mark.parametrize('pressure', [-1.0, 0.0, 0.5])
def test_ckv_is_closed_at_and_below_cutoff(pressure):
    element = Afe_Ckv(dp0=0.5, coef=10.0)
    link = SimpleNamespace(multiplier=1.0)

    assert element.calculate(link, pressure) == (1, 0.0, 0.0, 0.0, 0.0)


def test_ckv_linearize_is_closed():
    element = Afe_Ckv(dp0=0.5, coef=10.0)

    assert element.linearize(
        SimpleNamespace(multiplier=1.0), 2.0
    ) == (0.0, 0.0)


def test_ckv_open_flow_and_derivative_match_finite_difference():
    element = Afe_Ckv(dp0=0.5, coef=10.0)
    link = SimpleNamespace(multiplier=1.0)
    pressure = 2.0
    step = 1.0e-6

    count, flow, flow1, derivative, derivative1 = element.calculate(
        link, pressure
    )
    _, upper, _ = element.flow(link, pressure + step)
    _, lower, _ = element.flow(link, pressure - step)

    assert count == 1
    assert flow > 0.0
    assert flow1 == 0.0
    assert derivative1 == 0.0
    assert derivative == pytest.approx(
        (upper - lower) / (2.0 * step), rel=1.0e-6
    )


def test_ckv_multiplier_preserves_input_direction_when_link_is_reordered():
    element = Afe_Ckv(dp0=0.5, coef=10.0)
    link = SimpleNamespace(multiplier=1.0)

    forward_result = element.calculate(link, 2.0)
    link.multiplier = -1.0
    reordered_result = element.calculate(link, -2.0)

    assert reordered_result[1] == pytest.approx(-forward_result[1])
    assert reordered_result[3] == pytest.approx(forward_result[3])


def make_doorway():
    return Afe_Dor(
        init=0.015575,
        lam=0.015575,
        turb=1.76494,
        expt=0.5,
        dtmin=0.0001,
        ht=2.0,
        wd=0.8,
        cd=0.78,
    )


def make_duct():
    return Afe_Dwc(
        length=10.0,
        hdia=0.25,
        area=0.04909,
        rough=0.00015,
        tdlc=0.0,
        lflc=64.0,
        ldlc=0.0,
        linit=128.0,
    )


@pytest.mark.parametrize(
    ('element_factory', 'pressure', 'control'),
    [
        (
            lambda: Afe_Plr(
                init=0.01, lam=0.01, turb=0.1, expt=0.5
            ),
            3.0,
            1.0,
        ),
        (make_duct, 9.0, 1.0),
        (lambda: Afe_Qfr(a=0.31, b=44.8), 5.0, 1.0),
        (make_doorway, 0.98, 1.0),
        (lambda: Afe_Cfr(flow=1.25), 3.0, 0.6),
        (make_fan, -660.0, 1.0),
        (lambda: Afe_Cpf(upo=10.0, prmin=1.0, ftyp=2.0), -6.0, 1.0),
        (lambda: Afe_Ckv(dp0=0.5, coef=10.0), 2.0, 1.0),
    ],
    ids=['plr', 'dwc', 'qfr', 'dor', 'cfr', 'fan', 'cpf', 'ckv'],
)
def test_supported_element_multiplier_scales_all_flow_and_derivative_terms(
    element_factory, pressure, control
):
    multiplier = 2.5
    node0 = make_node(1.20, 1.8e-5, temperature=291.15)
    node1 = make_node(1.10, 1.9e-5, temperature=295.15)
    base_link = SimpleNamespace(
        node0=node0, node1=node1, multiplier=1.0
    )
    scaled_link = SimpleNamespace(
        node0=node0, node1=node1, multiplier=multiplier
    )

    base_linear = element_factory().linearize(
        base_link, pressure, control=control
    )
    scaled_linear = element_factory().linearize(
        scaled_link, pressure, control=control
    )
    assert scaled_linear == pytest.approx(
        tuple(multiplier * value for value in base_linear)
    )

    base_calculation = element_factory().calculate(
        base_link, pressure, control=control
    )
    scaled_calculation = element_factory().calculate(
        scaled_link, pressure, control=control
    )
    assert scaled_calculation[0] == base_calculation[0]
    assert scaled_calculation[1:] == pytest.approx(
        tuple(multiplier * value for value in base_calculation[1:])
    )

    base_flow = element_factory().flow(
        base_link, pressure, control=control
    )
    scaled_flow = element_factory().flow(
        scaled_link, pressure, control=control
    )
    assert scaled_flow[0] == base_flow[0]
    assert scaled_flow[1:] == pytest.approx(
        tuple(multiplier * value for value in base_flow[1:])
    )


def test_dor_uses_midheight_power_law_below_temperature_threshold():
    doorway = make_doorway()
    power_law = Afe_Plr(
        init=doorway.init,
        lam=doorway.lam,
        turb=doorway.turb,
        expt=doorway.expt,
    )
    link = SimpleNamespace(
        node0=make_node(1.20, 1.8e-5, temperature=293.15),
        node1=make_node(1.10, 1.9e-5, temperature=293.15005),
        multiplier=1.0,
    )
    pressure = 1.25
    adjusted_pressure = pressure - 0.5 * doorway.ht * 9.8 * 0.10

    assert doorway.linearize(link, pressure) == pytest.approx(
        power_law.linearize(link, adjusted_pressure)
    )
    assert doorway.calculate(link, pressure) == pytest.approx(
        power_law.calculate(link, adjusted_pressure)
    )


def test_dor_uses_power_law_when_density_difference_is_zero():
    doorway = make_doorway()
    power_law = Afe_Plr(
        init=doorway.init,
        lam=doorway.lam,
        turb=doorway.turb,
        expt=doorway.expt,
    )
    link = SimpleNamespace(
        node0=make_node(1.20, 1.8e-5, temperature=291.15),
        node1=make_node(1.20, 1.8e-5, temperature=295.15),
        multiplier=1.0,
    )

    assert doorway.calculate(link, 0.25) == pytest.approx(
        power_law.calculate(link, 0.25)
    )


def test_dor_two_way_component_jacobians_match_finite_difference():
    doorway = make_doorway()
    link = SimpleNamespace(
        node0=make_node(1.20, 1.8e-5, temperature=291.15),
        node1=make_node(1.10, 1.9e-5, temperature=295.15),
        multiplier=1.0,
    )
    pressure = 0.98
    step = 1.0e-6

    count, flow0, flow1, derivative0, derivative1 = doorway.calculate(
        link, pressure
    )
    _, upper0, upper1 = doorway.flow(link, pressure + step)
    _, lower0, lower1 = doorway.flow(link, pressure - step)

    assert count == 2
    assert flow0 < 0.0 < flow1
    assert derivative0 == pytest.approx(
        (upper0 - lower0) / (2.0 * step), rel=1.0e-6
    )
    assert derivative1 == pytest.approx(
        (upper1 - lower1) / (2.0 * step), rel=1.0e-6
    )
    assert doorway.flow(link, pressure) == pytest.approx(
        (count, flow0, flow1)
    )


def test_dor_two_way_flow_reverses_with_link_orientation():
    doorway = make_doorway()
    node0 = make_node(1.20, 1.8e-5, temperature=291.15)
    node1 = make_node(1.10, 1.9e-5, temperature=295.15)
    forward_link = SimpleNamespace(node0=node0, node1=node1, multiplier=1.0)
    reverse_link = SimpleNamespace(node0=node1, node1=node0, multiplier=1.0)

    forward = doorway.calculate(forward_link, 0.98)
    reverse = doorway.calculate(reverse_link, -0.98)

    assert forward[0] == reverse[0] == 2
    assert reverse[1] == pytest.approx(-forward[1])
    assert reverse[2] == pytest.approx(-forward[2])
    assert reverse[3] == pytest.approx(forward[3])
    assert reverse[4] == pytest.approx(forward[4])


@pytest.mark.parametrize(
    ("pressure", "expected_sign"),
    [(-0.49, -1.0), (2.45, 1.0)],
)
def test_dor_neutral_plane_outside_doorway_is_one_way(
    pressure, expected_sign
):
    doorway = make_doorway()
    link = SimpleNamespace(
        node0=make_node(1.20, 1.8e-5, temperature=291.15),
        node1=make_node(1.10, 1.9e-5, temperature=295.15),
        multiplier=1.0,
    )

    count, flow0, flow1, derivative0, derivative1 = doorway.calculate(
        link, pressure
    )

    assert count == 1
    assert flow0 * expected_sign > 0.0
    assert flow1 == 0.0
    assert derivative0 > 0.0
    assert derivative1 == 0.0


@pytest.mark.parametrize(
    ("neutral_height", "zero_flow_index"),
    [(0.0, 2), (2.0, 1)],
)
def test_dor_neutral_plane_boundaries_remain_two_way(
    neutral_height, zero_flow_index
):
    doorway = make_doorway()
    link = SimpleNamespace(
        node0=make_node(1.20, 1.8e-5, temperature=291.15),
        node1=make_node(1.10, 1.9e-5, temperature=295.15),
        multiplier=1.0,
    )
    pressure = neutral_height * 9.8 * (
        link.node0.density - link.node1.density
    )

    result = doorway.calculate(link, pressure)

    assert result[0] == 2
    assert result[zero_flow_index] == pytest.approx(0.0)
