# SPDX-FileCopyrightText: 2023-present Oak Ridge National Laboratory, managed by UT-Battelle
#
# SPDX-License-Identifier: BSD-3-Clause
from dataclasses import FrozenInstanceError

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


def test_weather_is_immutable_and_owns_an_immutable_profile_mapping():
    profile = airnet.WindProfile('constant', (1.0,) * 16)
    weather = airnet.Weather(
        ambient_temperature=10.0,
        wind_profiles={'ignored-key': profile},
    )

    assert weather.ambient_density == pytest.approx(
        0.0034838 * 101325.0 / 283.15
    )
    assert weather.wind_profiles['constant'] is profile
    with pytest.raises(FrozenInstanceError):
        weather.wind_speed = 5.0
    with pytest.raises(TypeError):
        weather.wind_profiles['another'] = profile


@pytest.mark.parametrize(
    ('kwargs', 'message'),
    [
        ({'ambient_temperature': -273.15}, 'absolute zero'),
        ({'barometric_pressure': 0.0}, 'must be positive'),
        ({'wind_speed': -1.0}, 'cannot be negative'),
        ({'wind_direction': float('nan')}, 'must be finite'),
    ],
)
def test_weather_rejects_invalid_environmental_inputs(kwargs, message):
    with pytest.raises(airnet.BadWeather, match=message):
        airnet.Weather(**kwargs)
