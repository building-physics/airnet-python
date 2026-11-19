# SPDX-FileCopyrightText: 2024-present Oak Ridge National Laboratory, managed by UT-Battelle, Alliance for Energy Innovation, LLC, and contributors
#
# SPDX-License-Identifier: BSD-3-Clause
from pathlib import Path

import click

import airnet
from airnet.__about__ import __version__


@click.command()
@click.argument(
    'input_file',
    metavar='INPUT_FILE',
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option('-v', '--verbose', is_flag=True, show_default=True, default=False, help='Operate verbosely.')
@click.option(
    '-o', '--output', 'output_file',
    type=click.Path(dir_okay=False, path_type=Path),
    default=Path('airnetsim.csv'), show_default=True,
    help='Write converged node and link results to this CSV file.',
)
@click.option(
    '--wind-file', type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help='Read AIRNET 16-direction wind-pressure profiles from this file.',
)
@click.option('--ambient-temperature', type=float, default=20.0, show_default=True)
@click.option('--barometric-pressure', type=float, default=101325.0, show_default=True)
@click.option('--wind-speed', type=click.FloatRange(min=0.0), default=0.0, show_default=True)
@click.option('--wind-direction', type=float, default=0.0, show_default=True)
@click.option(
    '--initialization-iterations', type=click.IntRange(min=1),
    default=100, show_default=True,
)
@click.option(
    '--maximum-iterations', type=click.IntRange(min=1), default=100, show_default=True,
)
@click.option(
    '--maximum-linear-iterations', type=click.IntRange(min=1),
    default=100, show_default=True,
)
@click.option(
    '--absolute-tolerance', type=click.FloatRange(min=0.0),
    default=1.0e-8, show_default=True,
)
@click.option('--relative-tolerance', type=click.FloatRange(min=0.0))
@click.option(
    '--relaxation', type=click.FloatRange(min=0.0, max=1.0, min_open=True),
    default=1.0, show_default=True,
)
def simulate(
    input_file,
    verbose,
    output_file,
    wind_file,
    ambient_temperature,
    barometric_pressure,
    wind_speed,
    wind_direction,
    initialization_iterations,
    maximum_iterations,
    maximum_linear_iterations,
    absolute_tolerance,
    relative_tolerance,
    relaxation,
):
    """Simulate the pressure network in the specified input file."""
    try:
        airnet.run_simulate(
            input_file,
            verbose=verbose,
            output_file=output_file,
            wind_file=wind_file,
            ambient_temperature=ambient_temperature,
            barometric_pressure=barometric_pressure,
            wind_speed=wind_speed,
            wind_direction=wind_direction,
            initialization_iterations=initialization_iterations,
            maximum_iterations=maximum_iterations,
            maximum_linear_iterations=maximum_linear_iterations,
            absolute_tolerance=absolute_tolerance,
            relative_tolerance=relative_tolerance,
            relaxation=relaxation,
        )
    except (
        OSError,
        airnet.BadNetworkInput,
        airnet.BadNetwork,
        airnet.SolverError,
    ) as error:
        raise click.ClickException(str(error)) from error
    click.echo(f'Results written to "{output_file}".')


@click.command()
@click.argument(
    'input_file',
    metavar='INPUT_FILE',
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option('-v', '--verbose', is_flag=True, show_default=True, default=False, help='Operate verbosely.')
def summarize(input_file, verbose):
    """Summarize the pressure network in the specified input file."""
    try:
        if verbose:
            click.echo(f'Opening input file "{input_file}"...')

        with input_file.open(encoding='utf-8') as input_stream:
            reader = airnet.Reader(input_stream)
            if verbose:
                click.echo(f'Reading input file "{input_file}"...')
            items = list(reader)
    except (OSError, airnet.BadNetworkInput) as error:
        raise click.ClickException(str(error)) from error
    if verbose:
        click.echo(f'Closing input file "{input_file}".')

    if reader.title:
        click.echo(f'Title: {reader.title}')

    click.echo('\nElements:\n=========')
    elements = {}
    nodes = 0
    links = 0
    for el in items:
        if el['input_type'] == airnet.InputType.ELEMENT:
            if el['type'] in elements:
                elements[el['type']] += 1
            else:
                elements[el['type']] = 1
        elif el['input_type'] == airnet.InputType.NODE:
            nodes += 1
        elif el['input_type'] == airnet.InputType.LINK:
            links += 1
    for key, value in elements.items():
        click.echo(f'{key} {value}')

    click.echo(f'\nNodes: {nodes}\n\nLinks: {links}')

    if verbose:
        click.echo('Done.')


@click.group(context_settings={'help_option_names': ['-h', '--help']}, invoke_without_command=False)
@click.version_option(version=__version__, prog_name='airnet')
def main():
    """Solve AIRNET pressure networks."""


main.add_command(summarize)
main.add_command(simulate)
