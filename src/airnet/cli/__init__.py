# SPDX-FileCopyrightText: 2024-present Oak Ridge National Laboratory, managed by UT-Battelle, Alliance for Energy Innovation, LLC, and contributors
#
# SPDX-License-Identifier: BSD-3-Clause
import click
import os
import airnet


@click.command()
@click.argument('input',metavar='INPUT_FILE', type=click.Path(exists=True))
@click.option('-v', '--verbose', is_flag=True, show_default=True, default=False, help='Operate verbosely.')
def simulate(input, verbose):
    airnet.run_simulate(input, verbose=verbose)


@click.command()
@click.argument('input',metavar='INPUT_FILE', type=click.Path(exists=True))
@click.option('-v', '--verbose', is_flag=True, show_default=True, default=False, help='Operate verbosely.')
def summarize(input, verbose):
    if verbose:
        print('Opening input file "%s"...' % input)

    fp = open(input, 'r')
    reader = airnet.Reader(fp)
    if verbose:
        print('Reading input file "%s"...' % input)
    items = []
    for item in reader:
        items.append(item)
    if verbose:
        print('Closing input file "%s".' % input)
    fp.close()

    if reader.title:
        print('Title:', reader.title)

    print('\nElements:\n=========')
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
        print(key, value)

    print('\nNodes: %s\n\nLinks: %s' % (nodes, links))

    if verbose:
        print('Done.')


@click.group(context_settings={'help_option_names': ['-h', '--help']}, invoke_without_command=False)
@click.version_option(version=__version__, prog_name='airnet')
@click.pass_context
def airnetx(ctx: click.Context):
    pass

airnetx.add_command(summarize)
