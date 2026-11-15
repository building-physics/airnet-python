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
def simulate(input_file, verbose):
    """Simulate the pressure network in the specified input file."""
    try:
        airnet.run_simulate(input_file, verbose=verbose)
    except (airnet.BadNetworkInput, airnet.BadNetwork, airnet.SolverError) as error:
        raise click.ClickException(str(error)) from error


@click.command()
@click.argument(
    'input_file',
    metavar='INPUT_FILE',
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option('-v', '--verbose', is_flag=True, show_default=True, default=False, help='Operate verbosely.')
def summarize(input_file, verbose):
    """Summarize the pressure network in the specified input file."""
    if verbose:
        click.echo(f'Opening input file "{input_file}"...')

    with input_file.open(encoding='utf-8') as input_stream:
        reader = airnet.Reader(input_stream)
        if verbose:
            click.echo(f'Reading input file "{input_file}"...')
        items = list(reader)
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
