# SPDX-FileCopyrightText: 2023-present Oak Ridge National Laboratory, managed by UT-Battelle
#
# SPDX-License-Identifier: BSD-3-Clause
import airnet
import os

def test_law_office_simple():
    contam_csv = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'law-office-simple-3405.csv')
    net_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'examples', 'law-office-simple.txt')
    old_cwd = os.getcwd()
    with airnet.temporary_directory():
        airnet.run_simulate(net_file)
        assert os.path.exists('airnetsim.csv')
        assert airnet.compare_csvs(contam_csv, 'airnetsim.csv') == []
        os.chdir(old_cwd)
