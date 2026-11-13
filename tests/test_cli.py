# SPDX-FileCopyrightText: 2023-present Oak Ridge National Laboratory, managed by UT-Battelle
#
# SPDX-License-Identifier: BSD-3-Clause
from pathlib import Path

from click.testing import CliRunner

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

    def record_call(path, verbose=False):
        calls.append((path, verbose))

    monkeypatch.setattr('airnet.run_simulate', record_call)

    result = CliRunner().invoke(simulate, ['--verbose', str(input_file)])

    assert result.exit_code == 0
    assert calls == [(input_file, True)]
