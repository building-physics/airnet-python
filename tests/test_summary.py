# SPDX-FileCopyrightText: 2023-present Oak Ridge National Laboratory, managed by UT-Battelle
#
# SPDX-License-Identifier: BSD-3-Clause
from airnet import model
from unittest import mock
import argparse

@mock.patch('argparse.ArgumentParser.parse_args',
            return_value=argparse.Namespace(input='somefile.txt', verbose=False))
def test_missing_file(mock_args):
    res = model.summarize_input()
    assert res == 1

@mock.patch('argparse.ArgumentParser.parse_args',
            return_value=argparse.Namespace(input='tests/AFDATA.FN3', verbose=True))
def test_complete_verbose(mock_args):
    res = model.summarize_input()
    assert res == 0

@mock.patch('argparse.ArgumentParser.parse_args',
            return_value=argparse.Namespace(input='tests/AFDATA.FNX', verbose=False))
def test_complete(mock_args):
    res = model.summarize_input()
    assert res == 0