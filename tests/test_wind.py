# SPDX-FileCopyrightText: 2023-present Oak Ridge National Laboratory, managed by UT-Battelle
#
# SPDX-License-Identifier: BSD-3-Clause
import pytest

import airnet

WIND_DATA = '''
/
title wind coefficients
north 1.0 0.924 0.707 0.383 0.0 -0.383 -0.707 -0.924
      -1.0 -0.924 -0.707 -0.383 0.0 0.383 0.707 0.924
constant 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
'''


def test_read_wind_profiles_and_interpolate_between_sectors():
    profiles = airnet.read_wind_profiles(WIND_DATA.splitlines())

    assert set(profiles) == {'north', 'constant'}
    assert profiles['north'].coefficient(0.0) == 1.0
    assert profiles['north'].coefficient(360.0) == 1.0
    assert profiles['north'].coefficient(11.25) == pytest.approx(0.962)
    assert profiles['north'].coefficient(-11.25) == pytest.approx(0.962)


def test_wind_profile_requires_sixteen_finite_coefficients():
    with pytest.raises(ValueError, match='exactly 16'):
        airnet.WindProfile('short', (1.0,) * 15)

    with pytest.raises(ValueError, match='finite'):
        airnet.WindProfile('nonfinite', (1.0,) * 15 + (float('nan'),))


def test_read_wind_profiles_reports_incomplete_profile():
    with pytest.raises(airnet.BadNetworkInput, match='fewer than 16'):
        airnet.read_wind_profiles(['incomplete 1.0 2.0'])
