# SPDX-FileCopyrightText: 2023-present Oak Ridge National Laboratory, managed by UT-Battelle
#
# SPDX-License-Identifier: BSD-3-Clause
from pathlib import Path

from click.testing import CliRunner

import airnet
from airnet.__about__ import __version__
from airnet.cli import main, simulate, summarize

TEST_DIRECTORY = Path(__file__).parent


def test_main_help_lists_commands():
    result = CliRunner().invoke(main, ['--help'])

    assert result.exit_code == 0
    assert 'simulate' in result.output
    assert 'summarize' in result.output


def test_main_version():
    result = CliRunner().invoke(main, ['--version'])

    assert result.exit_code == 0
    assert result.output.strip() == f'airnet, version {__version__}'


def test_summarize_rejects_missing_file():
    result = CliRunner().invoke(summarize, ['missing-airnet-input.txt'])

    assert result.exit_code == 2
    assert 'does not exist' in result.output


def test_summarize_input():
    input_file = TEST_DIRECTORY / 'AFDATA.FNX'

    result = CliRunner().invoke(summarize, [str(input_file)])

    assert result.exit_code == 0
    assert 'Elements:' in result.output
    assert 'Nodes: 5' in result.output
    assert 'Links: 5' in result.output


def test_summarize_input_verbose():
    input_file = TEST_DIRECTORY / 'AFDATA.FN3'

    result = CliRunner().invoke(summarize, ['--verbose', str(input_file)])

    assert result.exit_code == 0
    assert 'Opening input file' in result.output
    assert 'Reading input file' in result.output
    assert 'Closing input file' in result.output
    assert result.output.rstrip().endswith('Done.')


def test_simulate_calls_package_api(monkeypatch):
    input_file = TEST_DIRECTORY / 'AFDATA.FNX'
    calls = []

    def record_call(path, **kwargs):
        calls.append((path, kwargs))

    monkeypatch.setattr('airnet.run_simulate', record_call)

    result = CliRunner().invoke(
        simulate,
        [
            '--verbose',
            '--output', 'custom.csv',
            '--ambient-temperature', '-5.0',
            '--relative-tolerance', '0.0001',
            str(input_file),
        ],
    )

    assert result.exit_code == 0
    assert calls[0][0] == input_file
    assert calls[0][1]['verbose'] is True
    assert calls[0][1]['output_file'] == Path('custom.csv')
    assert calls[0][1]['ambient_temperature'] == -5.0
    assert calls[0][1]['relative_tolerance'] == 0.0001
    assert result.output.rstrip().endswith('Results written to "custom.csv".')


def test_simulate_reports_solver_failure_as_cli_error(monkeypatch):
    input_file = TEST_DIRECTORY / 'AFDATA.FNX'

    def fail_simulation(path, **kwargs):
        raise airnet.InitializationError(info=4, max_iterations=100)

    monkeypatch.setattr('airnet.run_simulate', fail_simulation)

    result = CliRunner().invoke(simulate, [str(input_file)])

    assert result.exit_code == 1
    assert 'Error: Initialization linear solve' in result.output
    assert 'solver info=4' in result.output


def test_simulate_reports_output_write_failure(monkeypatch):
    input_file = TEST_DIRECTORY / 'AFDATA.FNX'

    def fail_simulation(path, **kwargs):
        raise OSError('cannot write result file')

    monkeypatch.setattr('airnet.run_simulate', fail_simulation)

    result = CliRunner().invoke(simulate, [str(input_file)])

    assert result.exit_code == 1
    assert 'Error: cannot write result file' in result.output


def test_summarize_reports_bad_input_as_cli_error(tmp_path):
    input_file = tmp_path / 'bad-input.txt'
    input_file.write_text('element bad unknown 1.0\n', encoding='utf-8')

    result = CliRunner().invoke(summarize, [str(input_file)])

    assert result.exit_code == 1
    assert 'Error: Element type "unknown" not recognized' in result.output
