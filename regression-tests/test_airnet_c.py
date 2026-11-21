# SPDX-FileCopyrightText: 2023-present Oak Ridge National Laboratory, managed by UT-Battelle
#
# SPDX-License-Identifier: BSD-3-Clause
import json
import math
from pathlib import Path

import pytest

import airnet

FIXTURE_DIRECTORY = Path(__file__).parent / "airnet-c"

SWEEP_NODE_PRESSURES = {
    "PL1": {2: ("node-3", -0.1), 3: ("node-3", -0.01),
            4: ("node-3", -0.00161), 5: ("node-3", -0.0016076),
            6: ("node-3", -0.0000426), 7: ("node-3", -0.000042),
            8: ("node-3", -0.00004)},
    "FN4": {2: ("node-3", 200.0), 3: ("node-3", 800.0),
            4: ("node-3", -200.0), 5: ("node-3", -800.0)},
    "QF1": {
        run: ("source", pressure)
        for run, pressure in enumerate(
            (0.02, 0.04, 0.08, 0.16, 0.32, 0.64, 1.28, 2.56,
             5.12, 10.24, 20.48, 40.96),
            start=2,
        )
    },
    "CV1": {
        run: ("node-1", pressure)
        for run, pressure in enumerate((4.0, 1.0, 0.25, 0.05, 0.01, 0.0, -4.0), start=2)
    },
}

SWEEP_LINK_ELEMENTS = {
    "PL2": {
        run: ("link-2", element)
        for run, element in enumerate(
            ("orf-0.001", "orf-0.01", "orf-0.1", "orf-1.0",
             "orf-10.0", "orf-100.0"),
            start=2,
        )
    }
}


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
        weather=airnet.Weather(
            ambient_temperature=weather["ambient_temperature_c"],
            barometric_pressure=weather["barometric_pressure_pa"],
            wind_speed=weather["wind_speed_m_per_s"],
            wind_direction=weather["wind_direction_degrees"],
            wind_profiles=wind_profiles or {},
        ),
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


def load_collection_run(fixture_name, run_number):
    with (FIXTURE_DIRECTORY / fixture_name).open(encoding="utf-8") as stream:
        collection = json.load(stream)
    run = collection["runs"][run_number - 1]
    with (FIXTURE_DIRECTORY / collection["network"]).open(encoding="utf-8") as stream:
        records = list(airnet.Reader(stream))

    pressure_change = SWEEP_NODE_PRESSURES.get(collection["case"], {}).get(run_number)
    if pressure_change:
        name, pressure = pressure_change
        next(record for record in records if record.get("name") == name)["pres"] = pressure
    element_change = SWEEP_LINK_ELEMENTS.get(collection["case"], {}).get(run_number)
    if element_change:
        name, element = element_change
        next(record for record in records if record.get("name") == name)["element"] = element

    if run_number > 1:
        previous_pressures = collection["runs"][run_number - 2]["nodes"]
        for record in records:
            if (
                record.get("input_type") == airnet.InputType.NODE
                and record.get("type") == "v"
            ):
                record["pres"] = previous_pressures[record["name"]]

    wind_profiles = None
    if any(record.get("wind") not in (None, "null") for record in records):
        with (FIXTURE_DIRECTORY / "WIND").open(encoding="utf-8") as stream:
            wind_profiles = airnet.read_wind_profiles(stream)
    weather = run["weather"]
    model = airnet.Model(
        records,
        weather=airnet.Weather(
            ambient_temperature=weather["ambient_temperature_c"],
            barometric_pressure=weather["barometric_pressure_pa"],
            wind_speed=weather["wind_speed_m_per_s"],
            wind_direction=weather["wind_direction_degrees"],
            wind_profiles=wind_profiles or {},
        ),
    )
    if run_number > 1:
        for node in model.variable_nodes:
            node.pressure = 0.0
    return collection, run, model


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


COLLECTION_FIXTURES = [
    ("pl1-sweep.json", range(2, 9)),
    ("pl2-sweep.json", range(2, 8)),
    ("wp1-sweep.json", range(2, 5)),
    ("wp2-sweep.json", [1, *range(3, 14)]),
    ("fn4-sweep.json", [1]),
    ("qf1-sweep.json", range(2, 14)),
    ("cv1-sweep.json", range(2, 9)),
    ("dw2-runs.json", [1]),
    ("cf2-runs.json", [1]),
    ("qf2-runs.json", [1]),
    ("04f-runs.json", [1]),
    ("08f-runs.json", [1]),
    ("16f-runs.json", [1]),
    ("32f-runs.json", [1]),
    ("133-runs.json", [1]),
]


@pytest.mark.parametrize(
    ("fixture_name", "run_number"),
    [
        (fixture_name, run_number)
        for fixture_name, run_numbers in COLLECTION_FIXTURES
        for run_number in run_numbers
    ],
)
def test_airnet_c_collection_run(fixture_name, run_number):
    collection, run, model = load_collection_run(fixture_name, run_number)
    tolerances = collection["tolerances"]
    run_control = run["run_control"]

    linear_iterations = collection.get(
        "maximum_python_linear_iterations",
        max(run_control["maximum_iterations"], 10 * model.size),
    )
    model.initialize(maxiter=linear_iterations)
    iterations = model.air_movement(
        maxiter=run_control["maximum_iterations"],
        tolerance=run_control["absolute_airflow_tolerance"],
        relative_tolerance=run_control["relative_airflow_tolerance"],
        relaxation=run_control["relaxation_coefficient"],
        minimum_iterations=2,
        max_subiter=collection.get(
            "maximum_python_linear_iterations", max(100, 10 * model.size)
        ),
    )

    assert 1 <= iterations <= run["expected_iterations"]
    assert_finite_results(model)
    for name, expected in run["nodes"].items():
        assert model.nodes[name].pressure == pytest.approx(
            expected,
            rel=tolerances["relative"],
            abs=tolerances["pressure_absolute"],
        )

    links = {link.name: link for link in model.links}
    for name, expected in run["links"].items():
        actual = input_oriented_link_values(links[name])
        assert actual["pressure_drop"] == pytest.approx(
            expected[0],
            rel=tolerances["relative"],
            abs=tolerances["pressure_absolute"],
        )
        assert actual["flow0"] == pytest.approx(
            expected[1],
            rel=tolerances["relative"],
            abs=tolerances["flow_absolute"],
        )
        assert actual["flow1"] == pytest.approx(
            expected[2],
            rel=tolerances["relative"],
            abs=tolerances["flow_absolute"],
        )

    for name, residual in variable_node_mass_balance(model).items():
        assert residual == pytest.approx(
            0.0, abs=tolerances["mass_balance_absolute"]
        ), f"mass imbalance at {name}"
