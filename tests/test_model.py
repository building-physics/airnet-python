# SPDX-FileCopyrightText: 2023-present Oak Ridge National Laboratory, managed by UT-Battelle
#
# SPDX-License-Identifier: BSD-3-Clause
import copy
import math

import pytest

import airnet

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

    model = airnet.Model(items, weather=airnet.Weather())
    assert len(model.nodes) == 4
    assert sum(node.variable for node in model.nodes.values()) == 2
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
    node_1_temp = node_1.temperature
    node_1_dens = node_1.density
    model.set_properties(model.nodes.values(), global_temperature=node_1_temp, global_density=node_1_dens)
    model.nodes['node-4'].pressure = node_4_pressure

    solver = airnet.Solver(model)
    solver.initialize(maxiter=10)

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
    

    #model.set_variable_properties()
    iters = solver.steady_solve(maxiter=50)
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


def test_linear_initialization_uses_flow_direction_density():
    model = airnet.Model(
        airnet.Reader(iter(AFDATA_PL2.splitlines())), weather=airnet.Weather()
    )

    airnet.Solver(model).initialize(maxiter=10)
    assert model.nodes['node-2'].pressure == pytest.approx(-100.0 / 3.0)
    assert model.nodes['node-3'].pressure == pytest.approx(-200.0 / 3.0)


def test_steady_solve_reports_nonconvergence():
    model = airnet.Model(
        airnet.Reader(iter(AFDATA_PL2.splitlines())), weather=airnet.Weather()
    )
    solver = airnet.Solver(model)
    solver.initialize(maxiter=10)

    with pytest.raises(airnet.NonlinearSolverError, match='failed to converge') as error:
        solver.steady_solve(maxiter=1, tolerance=0.0)

    assert error.value.max_iterations == 1
    assert error.value.maximum_absolute_residual == 0.0
    assert error.value.relative_residual == 0.0
    assert error.value.absolute_tolerance == 0.0


def test_initialization_reports_linear_solver_failure(monkeypatch):
    model = airnet.Model(
        airnet.Reader(iter(AFDATA_PL2.splitlines())), weather=airnet.Weather()
    )

    def fail_linear_solve(matrix, right_hand_side, **kwargs):
        return right_hand_side, 7

    monkeypatch.setattr('airnet.model.scipy.sparse.linalg.cg', fail_linear_solve)

    with pytest.raises(airnet.InitializationError, match='did not converge') as error:
        airnet.Solver(model).initialize(maxiter=5)

    assert error.value.info == 7
    assert error.value.phase == 'initialization'
    assert error.value.max_iterations == 5


def test_steady_solve_reports_linear_correction_failure(monkeypatch):
    model = airnet.Model(
        airnet.Reader(iter(AFDATA_PL2.splitlines())), weather=airnet.Weather()
    )
    solver = airnet.Solver(model)
    solver.initialize(maxiter=10)
    model.nodes['node-2'].pressure += 1.0

    def fail_linear_solve(matrix, right_hand_side, **kwargs):
        return right_hand_side, -1

    monkeypatch.setattr('airnet.model.scipy.sparse.linalg.cg', fail_linear_solve)

    with pytest.raises(airnet.LinearSolverError, match='solver info=-1') as error:
        solver.steady_solve(maxiter=10)

    assert error.value.info == -1
    assert error.value.phase == 'nonlinear correction'
    assert error.value.nonlinear_iteration == 1


def test_steady_solve_rejects_non_finite_element_results(monkeypatch):
    model = airnet.Model(
        airnet.Reader(iter(AFDATA_PL2.splitlines())), weather=airnet.Weather()
    )
    solver = airnet.Solver(model)
    solver.initialize(maxiter=10)

    monkeypatch.setattr(
        model.links[0].element,
        'calculate',
        lambda link, pdrop, control=1.0: (
            1,
            float('nan'),
            0.0,
            1.0,
            0.0,
        ),
    )

    with pytest.raises(airnet.NumericalSolverError, match='link "link-1"'):
        solver.steady_solve(maxiter=10)


def test_run_simulate_does_not_write_results_after_initialization_failure(
    monkeypatch, tmp_path
):
    input_file = tmp_path / 'network.txt'
    input_file.write_text(AFDATA_PL2, encoding='utf-8')
    monkeypatch.chdir(tmp_path)

    def fail_initialization(self, maxiter=100):
        raise airnet.InitializationError(info=3, max_iterations=maxiter)

    monkeypatch.setattr(airnet.Solver, 'initialize', fail_initialization)

    with pytest.raises(airnet.InitializationError):
        airnet.run_simulate(input_file)

    assert not (tmp_path / 'airnetsim.csv').exists()


def test_run_simulate_uses_partitioned_api_and_writes_selected_output(
    monkeypatch, tmp_path
):
    input_file = tmp_path / 'network.txt'
    output_file = tmp_path / 'selected-results.csv'
    input_file.write_text(AFDATA_PL2, encoding='utf-8')
    original_read_airnet = airnet.Model.read_airnet
    read_models = []

    def record_read(input_source, **options):
        read_model = original_read_airnet(input_source, **options)
        assert read_model.weather is None
        read_models.append(read_model)
        return read_model

    monkeypatch.setattr(
        airnet.Model,
        'read_airnet',
        staticmethod(record_read),
    )

    model = airnet.run_simulate(input_file, output_file=output_file)

    assert isinstance(model, airnet.Model)
    assert read_models == [model]
    assert isinstance(model.weather, airnet.Weather)
    output = output_file.read_text(encoding='utf-8')
    assert output.startswith('node header,')
    link_1 = next(line for line in output.splitlines() if line.startswith('link, link-1,'))
    fields = [field.strip() for field in link_1.split(',')]
    assert float(fields[3]) > 0.0
    assert float(fields[4]) > 0.0


class TwoWayElement:
    def type(self):
        return 'two'

    def linearize(self, link, pdrop, control=1.0):
        return 2.0 * pdrop, 2.0

    def calculate(self, link, pdrop, control=1.0):
        return 2, pdrop + 2.0, pdrop - 1.0, 1.0, 1.0

    def flow(self, link, pdrop, control=1.0):
        nf, flow0, flow1, _, _ = self.calculate(
            link, pdrop, control=control
        )
        return nf, flow0, flow1


def two_way_network(variable_node=True):
    return [
        {'input_type': airnet.InputType.TITLE, 'title': 'two-way assembly'},
        {
            'input_type': airnet.InputType.NODE,
            'name': 'node-1',
            'type': 'v' if variable_node else 'c',
            'ht': 0.0,
            'temp': 20.0,
            'pres': 0.0,
        },
        {
            'input_type': airnet.InputType.NODE,
            'name': 'node-2',
            'type': 'c',
            'ht': 0.0,
            'temp': 20.0,
            'pres': 0.0,
        },
        {
            'input_type': airnet.InputType.ELEMENT,
            'type': 'two',
            'name': 'two-way',
        },
        {
            'input_type': airnet.InputType.LINK,
            'name': 'link-1',
            'node-1': 'node-1',
            'ht-1': 0.0,
            'node-2': 'node-2',
            'ht-2': 0.0,
            'element': 'two-way',
            'wind': None,
            'wpmod': 0.0,
        },
    ]


def test_steady_solve_assembles_two_way_flow_and_jacobian_components():
    model = airnet.Model(
        two_way_network(),
        element_lookup={'two': TwoWayElement},
        weather=airnet.Weather(),
    )
    solver = airnet.Solver(model)
    solver.initialize(maxiter=10)

    iterations = solver.steady_solve(maxiter=10, minimum_iterations=2)

    assert iterations == 2
    assert model.nodes['node-1'].pressure == pytest.approx(-0.5)
    assert model.links[0].flow_count == 2
    assert model.links[0].flow0 == pytest.approx(1.5)
    assert model.links[0].flow1 == pytest.approx(-1.5)


def test_steady_solve_retains_fixed_to_fixed_two_way_flows():
    model = airnet.Model(
        two_way_network(variable_node=False),
        element_lookup={'two': TwoWayElement},
        weather=airnet.Weather(),
    )

    iterations = airnet.Solver(model).steady_solve(maxiter=10)

    assert iterations == 1
    assert model.links[0].flow_count == 2
    assert model.links[0].flow0 == pytest.approx(2.0)
    assert model.links[0].flow1 == pytest.approx(-1.0)


def valid_network_records():
    return [
        {'input_type': airnet.InputType.TITLE, 'title': 'validation network'},
        {
            'input_type': airnet.InputType.NODE,
            'name': 'variable',
            'type': 'v',
            'ht': 0.0,
            'temp': 20.0,
        },
        {
            'input_type': airnet.InputType.NODE,
            'name': 'fixed',
            'type': 'c',
            'ht': 0.0,
            'temp': 20.0,
            'pres': 0.0,
        },
        {
            'input_type': airnet.InputType.ELEMENT,
            'name': 'opening',
            'type': 'plr',
            'init': 1.0,
            'lam': 1.0,
            'turb': 1.0,
            'expt': 0.5,
        },
        {
            'input_type': airnet.InputType.LINK,
            'name': 'connection',
            'node-1': 'variable',
            'ht-1': 0.0,
            'node-2': 'fixed',
            'ht-2': 0.0,
            'element': 'opening',
        },
    ]


def record_named(records, name):
    return next(record for record in records if record.get('name') == name)


def test_cfr_uses_negative_multiplier_when_link_is_reordered():
    records = valid_network_records()
    records.extend(
        [
            {
                'input_type': airnet.InputType.ELEMENT,
                'name': 'supply',
                'type': 'cfr',
                'flow': 1.0,
            },
            {
                'input_type': airnet.InputType.LINK,
                'name': 'supply-link',
                'node-1': 'fixed',
                'ht-1': 0.0,
                'node-2': 'variable',
                'ht-2': 0.0,
                'element': 'supply',
            },
        ]
    )
    model = airnet.Model(records, weather=airnet.Weather())

    solver = airnet.Solver(model)
    solver.initialize(maxiter=10)

    supply_link = next(link for link in model.links if link.name == 'supply-link')
    assert supply_link.flipped
    assert supply_link.multiplier == -1.0
    assert supply_link.flow0 == -1.0

    solver.steady_solve(maxiter=20)

    assert supply_link.flow0 == -1.0


def test_ckv_uses_negative_multiplier_when_link_is_reordered():
    records = valid_network_records()
    records.extend(
        [
            {
                'input_type': airnet.InputType.NODE,
                'name': 'high-pressure',
                'type': 'c',
                'ht': 0.0,
                'temp': 20.0,
                'pres': 10.0,
            },
            {
                'input_type': airnet.InputType.ELEMENT,
                'name': 'check-valve',
                'type': 'ckv',
                'dp0': 0.5,
                'coef': 0.1,
            },
            {
                'input_type': airnet.InputType.LINK,
                'name': 'check-valve-link',
                'node-1': 'high-pressure',
                'ht-1': 0.0,
                'node-2': 'variable',
                'ht-2': 0.0,
                'element': 'check-valve',
            },
        ]
    )
    model = airnet.Model(records, weather=airnet.Weather())

    solver = airnet.Solver(model)
    solver.initialize(maxiter=10)
    solver.steady_solve(maxiter=50)

    opening_link = next(link for link in model.links if link.name == 'connection')
    valve_link = next(
        link for link in model.links if link.name == 'check-valve-link'
    )
    assert valve_link.flipped
    assert valve_link.multiplier == -1.0
    assert valve_link.flow0 < 0.0
    assert opening_link.flow0 == pytest.approx(-valve_link.flow0)


def test_qfr_uses_passive_element_initialization_sign():
    records = valid_network_records()
    record_named(records, 'fixed')['pres'] = 10.0
    element = record_named(records, 'opening')
    element.clear()
    element.update(
        {
            'input_type': airnet.InputType.ELEMENT,
            'name': 'opening',
            'type': 'qfr',
            'a': 2.0,
            'b': 3.0,
        }
    )
    model = airnet.Model(records, weather=airnet.Weather())

    airnet.Solver(model).initialize(maxiter=10)

    assert model.nodes['variable'].pressure == pytest.approx(10.0)
    assert model.links[0].flow0 == pytest.approx(0.0)


def test_model_construction_does_not_mutate_input_records():
    records = valid_network_records()
    original = copy.deepcopy(records)

    airnet.Model(records)

    assert records == original


@pytest.mark.parametrize(
    ('record_name', 'exception', 'message'),
    [
        ('variable', airnet.DuplicateNetworkName, 'Node "variable"'),
        ('opening', airnet.DuplicateNetworkName, 'Element "opening"'),
        ('connection', airnet.DuplicateNetworkName, 'Link "connection"'),
    ],
)
def test_model_rejects_duplicate_names(record_name, exception, message):
    records = valid_network_records()
    records.append(copy.deepcopy(record_named(records, record_name)))

    with pytest.raises(exception, match=message):
        airnet.Model(records)


@pytest.mark.parametrize(
    ('field', 'missing_name', 'message'),
    [
        ('node-1', 'missing-node', 'undefined node "missing-node"'),
        ('element', 'missing-element', 'undefined element "missing-element"'),
    ],
)
def test_model_rejects_missing_link_references(field, missing_name, message):
    records = valid_network_records()
    record_named(records, 'connection')[field] = missing_name

    with pytest.raises(airnet.MissingNetworkReference, match=message):
        airnet.Model(records)


def test_model_rejects_disconnected_nodes():
    records = valid_network_records()
    records.append(
        {
            'input_type': airnet.InputType.NODE,
            'name': 'unconnected',
            'type': 'c',
            'ht': 0.0,
            'temp': 20.0,
            'pres': 0.0,
        }
    )

    with pytest.raises(airnet.BadNetwork, match='Disconnected nodes found: unconnected'):
        airnet.Model(records)


def test_model_rejects_variable_component_without_fixed_pressure_node():
    records = valid_network_records()
    fixed = record_named(records, 'fixed')
    fixed['type'] = 'v'
    fixed.pop('pres')

    with pytest.raises(
        airnet.BadNetwork,
        match='Variable-pressure component has no linear-initialization path',
    ):
        airnet.Model(records)


def test_constant_flow_does_not_anchor_variable_pressure_component():
    records = valid_network_records()
    element = record_named(records, 'opening')
    element.clear()
    element.update(
        {
            'input_type': airnet.InputType.ELEMENT,
            'name': 'opening',
            'type': 'cfr',
            'flow': 1.0,
        }
    )

    with pytest.raises(
        airnet.BadNetwork,
        match='no linear-initialization path to a fixed-pressure node',
    ):
        airnet.Model(records)


def test_check_valve_does_not_anchor_linear_initialization():
    records = valid_network_records()
    element = record_named(records, 'opening')
    element.clear()
    element.update(
        {
            'input_type': airnet.InputType.ELEMENT,
            'name': 'opening',
            'type': 'ckv',
            'dp0': 0.0,
            'coef': 10.0,
        }
    )

    with pytest.raises(
        airnet.BadNetwork,
        match='no linear-initialization path to a fixed-pressure node',
    ):
        airnet.Model(records)


def test_model_constructs_supported_fan_element():
    records = valid_network_records()
    element = record_named(records, 'opening')
    element.clear()
    element.update(
        {
            'input_type': airnet.InputType.ELEMENT,
            'name': 'opening',
            'type': 'fan',
            'init': 1.0,
            'lam': 1.0,
            'turb': 1.0,
            'expt': 0.5,
            'rdens': 1.2,
            'fdf': 1.0,
            'sop': 1.0,
                'off': 0.0,
                'mf1': 100.0,
                'pts': [
                    {
                        'a1': 1.0,
                        'a2': -1.0,
                        'a3': 0.0,
                        'a4': -1.0,
                        'mf': -100.0,
                    }
                ],
            }
        )

    model = airnet.Model(records)

    assert model.elements['opening'].type() == 'fan'


@pytest.mark.parametrize(
    ('updates', 'message'),
    [
        ({'control': 0.5}, 'does not control'),
    ],
)
def test_model_rejects_unsupported_link_features(updates, message):
    records = valid_network_records()
    record_named(records, 'connection').update(updates)

    with pytest.raises(airnet.UnsupportedNetworkFeature, match=message):
        airnet.Model(records)


def test_model_rejects_missing_wind_profile():
    records = valid_network_records()
    record_named(records, 'connection').update(
        {'wind': 'missing-profile', 'wpmod': 1.0}
    )

    model = airnet.Model(records)

    with pytest.raises(
        airnet.MissingNetworkReference, match='undefined wind profile'
    ):
        model.weather = airnet.Weather()

    assert model.weather is None


def test_model_applies_interpolated_wind_pressure_in_stored_orientation():
    records = valid_network_records()
    connection = record_named(records, 'connection')
    connection.update(
        {
            'node-1': 'fixed',
            'node-2': 'variable',
            'wind': 'constant',
            'wpmod': 2.0,
        }
    )
    profile = airnet.WindProfile('constant', (1.0,) * 16)

    weather = airnet.Weather(
        ambient_temperature=20.0,
        barometric_pressure=101325.0,
        wind_speed=5.0,
        wind_direction=11.25,
        wind_profiles={'constant': profile},
    )
    model = airnet.Model(records, weather=weather)
    airnet.Solver(model)
    model.compute_pressure_drops()

    expected_input_pressure = 2.0 * 0.5 * weather.ambient_density * 5.0**2
    # The fixed-to-variable input link is stored variable-to-fixed.
    assert model.links[0].flipped
    assert model.links[0].pdrop == pytest.approx(-expected_input_pressure)


def test_stack_pressure_uses_average_density_before_flow_direction_is_known():
    records = valid_network_records()
    record_named(records, 'fixed')['temp'] = 0.0
    link = record_named(records, 'connection')
    link['ht-1'] = 5.0
    link['ht-2'] = 0.0

    model = airnet.Model(records, weather=airnet.Weather())
    airnet.Solver(model)
    model.compute_pressure_drops()

    stored_link = model.links[0]
    expected = (
        -9.80 * stored_link.node0.density * stored_link.ht0
        + 9.80 * stored_link.node1.density * stored_link.ht1
        + 4.90
        * (stored_link.node0.density + stored_link.node1.density)
        * (
            stored_link.node0.height
            - stored_link.node1.height
            + stored_link.ht0
            - stored_link.ht1
        )
    )
    assert stored_link.pdrop == pytest.approx(expected)


def test_model_rejects_nonpositive_link_multiplier():
    records = valid_network_records()
    record_named(records, 'connection')['mult'] = 0.0

    with pytest.raises(
        airnet.BadNetwork, match='multiplier must be finite and positive'
    ):
        airnet.Model(records)


def test_model_accepts_positive_link_multiplier():
    records = valid_network_records()
    record_named(records, 'connection')['mult'] = 2.5

    model = airnet.Model(records)

    assert model.links[0].multiplier == 2.5


def test_model_uses_weather_temperature_for_ambient_node():
    records = valid_network_records()
    record_named(records, 'fixed')['type'] = 'a'

    model = airnet.Model(
        records, weather=airnet.Weather(ambient_temperature=-5.0)
    )

    ambient = model.nodes['fixed']
    assert ambient.ambient
    assert not ambient.variable
    assert ambient.temperature == pytest.approx(268.15)


def test_read_airnet_leaves_weather_unbound(tmp_path):
    input_file = tmp_path / 'network.txt'
    input_file.write_text(AFDATA_PL2, encoding='utf-8')

    model = airnet.Model.read_airnet(input_file)

    assert model.weather is None
    with pytest.raises(airnet.MissingWeather, match='assign model.weather'):
        airnet.Solver(model)


def test_weather_binding_is_atomic_when_wind_profile_is_missing():
    records = valid_network_records()
    record_named(records, 'connection').update(
        {'wind': 'constant', 'wpmod': 1.0}
    )
    profile = airnet.WindProfile('constant', (1.0,) * 16)
    original_weather = airnet.Weather(
        ambient_temperature=10.0,
        wind_profiles={'constant': profile},
    )
    model = airnet.Model(records, weather=original_weather)
    original_properties = {
        name: (node.temperature, node.density, node.viscosity)
        for name, node in model.nodes.items()
    }

    with pytest.raises(
        airnet.MissingNetworkReference, match='undefined wind profile'
    ):
        model.weather = airnet.Weather(ambient_temperature=-5.0)

    assert model.weather is original_weather
    assert {
        name: (node.temperature, node.density, node.viscosity)
        for name, node in model.nodes.items()
    } == original_properties


def test_solver_delegates_steady_solution_and_rejects_rebound_weather():
    model = airnet.Model(
        airnet.Reader(iter(AFDATA_PL2.splitlines())), weather=airnet.Weather()
    )
    solver = airnet.Solver(model)

    solver.initialize(maxiter=10)
    iterations = solver.steady_solve(maxiter=50)

    assert iterations < 5
    model.weather = airnet.Weather(ambient_temperature=10.0)
    with pytest.raises(airnet.StaleSolverError, match='create a new Solver'):
        solver.steady_solve()


def test_solver_owns_workspace_and_assigns_indices_on_model_nodes():
    records = valid_network_records()
    connection = record_named(records, 'connection')
    connection['node-1'], connection['node-2'] = (
        connection['node-2'],
        connection['node-1'],
    )
    model = airnet.Model(records, weather=airnet.Weather())
    link = model.links[0]

    assert model.nodes['variable'].index is None
    assert link.node0 is model.nodes['fixed']
    assert not link.flipped
    for name in (
        'A',
        'x',
        'correction',
        'size',
        'variable_nodes',
        'initialize',
        'air_movement',
    ):
        assert not hasattr(model, name)

    first = airnet.Solver(model)

    assert model.nodes['variable'].index == 0
    assert link.node0 is model.nodes['variable']
    assert link.flipped
    assert link.multiplier == -1.0
    assert not hasattr(first, 'variable_nodes')
    assert first.A.shape == (1, 1)
    assert first.x.shape == (1, 1)
    assert first.correction.shape == (1,)

    second = airnet.Solver(model)

    assert second.A is not first.A
    assert second.x is not first.x
    assert second.correction is not first.correction
    assert link.node0 is model.nodes['variable']
    assert link.multiplier == -1.0
    first.A.data.fill(7.0)
    first.x.fill(7.0)
    first.correction.fill(7.0)
    assert not (second.A.data == 7.0).all()
    assert not (second.x == 7.0).all()
    assert not (second.correction == 7.0).all()


def test_input_oriented_link_results_reverse_solver_normalization():
    records = valid_network_records()
    connection = record_named(records, 'connection')
    connection['node-1'], connection['node-2'] = (
        connection['node-2'],
        connection['node-1'],
    )
    model = airnet.Model(records, weather=airnet.Weather())
    link = model.links[0]
    airnet.Solver(model)
    link.pdrop = -3.0
    link.flow0 = -2.0
    link.flow1 = 0.5

    assert link.input_pressure_drop == 3.0
    assert link.input_flow0 == 2.0
    assert link.input_flow1 == -0.5


def test_transient_step_is_explicitly_unsupported_without_model_mutation():
    model = airnet.Model(
        airnet.Reader(iter(AFDATA_PL2.splitlines())), weather=airnet.Weather()
    )
    solver = airnet.Solver(model)
    pressures = {name: node.pressure for name, node in model.nodes.items()}
    flows = [(link.flow0, link.flow1) for link in model.links]

    with pytest.raises(NotImplementedError, match='does not define a transient'):
        solver.transient_step(60.0)

    assert {name: node.pressure for name, node in model.nodes.items()} == pressures
    assert [(link.flow0, link.flow1) for link in model.links] == flows
