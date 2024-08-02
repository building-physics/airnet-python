# SPDX-FileCopyrightText: 2023-present Oak Ridge National Laboratory, managed by UT-Battelle
#
# SPDX-License-Identifier: BSD-3-Clause
import airnet
import os

def test_law_office_simple():
    current_dir = os.path.dirname(os.path.abspath(__file__))
    contam_csv = os.path.join(current_dir, 'law-office-simple-3405.csv')
    airnet_csv = os.path.join(current_dir, 'law-office-simple-airnet.csv')
    net_file = os.path.join(current_dir, '..', 'examples', 'law-office-simple.txt')
    with airnet.temporary_directory():
        airnet.run_simulate(net_file, global_density=1.2040973677927915)
        assert os.path.exists('airnetsim.csv')
        assert airnet.compare_csvs(contam_csv, 'airnetsim.csv') == []
        assert airnet.compare_csvs(airnet_csv, 'airnetsim.csv', node_tolerance=1.0e-15) == []
