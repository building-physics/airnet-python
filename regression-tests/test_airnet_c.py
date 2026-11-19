# SPDX-FileCopyrightText: 2023-present Oak Ridge National Laboratory, managed by UT-Battelle
#
# SPDX-License-Identifier: BSD-3-Clause
import json
import math
from pathlib import Path

import pytest

import airnet

FIXTURE_DIRECTORY = Path(__file__).parent / "airnet-c"


def load_case(name):
    with (FIXTURE_DIRECTORY / name).open(encoding="utf-8") as stream:
        case = json.load(stream)
    with (FIXTURE_DIRECTORY / case["network"]).open(encoding="utf-8") as stream:
        records = list(airnet.Reader(stream))
    for record in records:
        initial_pressure = case.get("initial_node_pressures", {}).get(
            record.get("name")
        )
        if record.get("input_type") == airnet.InputType.NODE and initial_pressure is not None:
            record["pres"] = initial_pressure
    wind_profiles = None
    if wind_file := case.get("wind_profiles"):
        with (FIXTURE_DIRECTORY / wind_file).open(encoding="utf-8") as stream:
            wind_profiles = airnet.read_wind_profiles(stream)
    weather = case.get(
        "weather",
        {
            "ambient_temperature_c": 20.0,
            "barometric_pressure_pa": 101325.0,
            "wind_speed_m_per_s": 0.0,
            "wind_direction_degrees": 0.0,
        },
    )
    model = airnet.Model(
        records,
        wind_profiles=wind_profiles,
        ambient_temperature=weather["ambient_temperature_c"],
        barometric_pressure=weather["barometric_pressure_pa"],
        wind_speed=weather["wind_speed_m_per_s"],
        wind_direction=weather["wind_direction_degrees"],
    )
    # AIRNET evaluates properties from the pressure state inherited from the
    # preceding run, then performs a fresh linear pressure initialization.
    # Reproduce that hidden state without carrying its pressure guess into
    # this independent single-shot solve.
    if case.get("initial_node_pressures"):
        for node in model.variable_nodes:
            node.pressure = 0.0
    links = {link.name: link for link in model.links}
    for link_name, control in case.get("link_controls", {}).items():
        links[link_name].control = control
    return case, model


def input_oriented_link_values(link):
    """Return link results in the orientation used by the input file."""
    direction = -1.0 if link.flipped else 1.0
    return {
        "flow_count": link.flow_count,
        "pressure_drop": direction * link.pdrop,
        "flow0": direction * link.flow0,
        "flow1": direction * link.flow1,
    }


def assert_finite_results(model):
    for node in model.nodes.values():
        assert math.isfinite(node.pressure), f"non-finite pressure at {node.name}"
    for link in model.links:
        assert math.isfinite(link.pdrop), f"non-finite pressure drop at {link.name}"
        assert math.isfinite(link.flow0), f"non-finite flow0 at {link.name}"
        assert math.isfinite(link.flow1), f"non-finite flow1 at {link.name}"


def variable_node_mass_balance(model):
    residuals = {node.name: 0.0 for node in model.variable_nodes}
    for link in model.links:
        if link.node0.variable:
            residuals[link.node0.name] += link.flow0 + link.flow1
        if link.node1.variable:
            residuals[link.node1.name] -= link.flow0 + link.flow1
    return residuals


@pytest.mark.parametrize(
    "fixture_name",
    [
        "pl1-run-01.json",
        "pl2-run-01.json",
        "pl3-run-01.json",
        "dw1-run-01.json",
        "dr1-01-run-01.json",
        "dr1-run-01.json",
        "st1-run-01.json",
        "st2-run-01.json",
        "wp1-run-01.json",
        "wp2-run-02.json",
        "qf1-run-01.json",
        "cf1-run-01.json",
        "cv1-run-01.json",
        "fn1-run-01.json",
        "fn1-run-02.json",
        "fn1-run-03.json",
        "fn1-run-04.json",
        "fn1-run-05.json",
        "fn2-run-01.json",
        "fn3-run-01.json",
        "fn4-run-02.json",
        "fn4-run-03.json",
        "fn4-run-04.json",
        "fn4-run-05.json",
        "cp1-run-01.json",
        "cp2-run-01.json",
    ],
)
def test_airnet_c_single_shot(fixture_name):
    case, model = load_case(fixture_name)
    tolerances = case["tolerances"]

    model.initialize(maxiter=case["run_control"]["maximum_iterations"])
    iterations = model.air_movement(
        maxiter=case["run_control"]["maximum_iterations"],
        tolerance=case["run_control"]["absolute_airflow_tolerance"],
        relative_tolerance=case["run_control"]["relative_airflow_tolerance"],
        relaxation=case["run_control"]["relaxation_coefficient"],
        minimum_iterations=2,
    )

    # Solver implementation details may reduce the nonlinear iteration count;
    # exceeding the AIRNET-C reference count is the regression concern.
    assert 1 <= iterations <= case.get(
        "maximum_python_iterations", case["expected_iterations"]
    )
    assert_finite_results(model)

    for name, expected in case["nodes"].items():
        assert model.nodes[name].pressure == pytest.approx(
            expected["pressure"],
            rel=tolerances["relative"],
            abs=tolerances["pressure_absolute"],
        )

    links = {link.name: link for link in model.links}
    for name, expected in case["links"].items():
        actual = input_oriented_link_values(links[name])
        if "flow_count" in expected:
            assert actual["flow_count"] == expected["flow_count"]
        assert actual["pressure_drop"] == pytest.approx(
            expected["pressure_drop"],
            rel=tolerances["relative"],
            abs=tolerances["pressure_absolute"],
        )
        assert actual["flow0"] == pytest.approx(
            expected["flow0"],
            rel=tolerances["relative"],
            abs=expected.get("flow_absolute", tolerances["flow_absolute"]),
        )
        assert actual["flow1"] == pytest.approx(
            expected["flow1"],
            rel=tolerances["relative"],
            abs=expected.get("flow_absolute", tolerances["flow_absolute"]),
        )

    for name, residual in variable_node_mass_balance(model).items():
        assert residual == pytest.approx(
            0.0, abs=tolerances["mass_balance_absolute"]
        ), f"mass imbalance at {name}"
