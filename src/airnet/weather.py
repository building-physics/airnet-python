# SPDX-FileCopyrightText: 2023-present Oak Ridge National Laboratory, managed by UT-Battelle
#
# SPDX-License-Identifier: BSD-3-Clause
from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass, field
from types import MappingProxyType

from .wind import WindProfile


class BadWeather(ValueError):
    """Raised when environmental inputs cannot define valid weather."""


@dataclass(frozen=True)
class Weather:
    """Immutable environmental inputs for one airflow-network solve.

    Temperatures are expressed in degrees Celsius, pressure in pascals, wind
    speed in metres per second, and wind direction in degrees clockwise from
    north.
    """

    ambient_temperature: float = 20.0
    barometric_pressure: float = 101325.0
    wind_speed: float = 0.0
    wind_direction: float = 0.0
    wind_profiles: Mapping[str, WindProfile] = field(default_factory=dict)

    def __post_init__(self):
        values = (
            self.ambient_temperature,
            self.barometric_pressure,
            self.wind_speed,
            self.wind_direction,
        )
        if not all(math.isfinite(value) for value in values):
            raise BadWeather('Weather values must be finite')
        if self.ambient_temperature <= -273.15:
            raise BadWeather('Ambient temperature must be above absolute zero')
        if self.barometric_pressure <= 0.0:
            raise BadWeather('Barometric pressure must be positive')
        if self.wind_speed < 0.0:
            raise BadWeather('Wind speed cannot be negative')

        profiles = {}
        source = (
            self.wind_profiles.values()
            if isinstance(self.wind_profiles, Mapping)
            else self.wind_profiles
        )
        for profile in source:
            if not isinstance(profile, WindProfile):
                raise BadWeather('Wind profiles must be WindProfile objects')
            if profile.name in profiles:
                raise BadWeather(
                    f'Wind profile "{profile.name}" is defined more than once'
                )
            profiles[profile.name] = profile
        object.__setattr__(self, 'wind_profiles', MappingProxyType(profiles))

    @property
    def ambient_density(self) -> float:
        """Return ambient air density in kilograms per cubic metre."""
        return (
            0.0034838
            * self.barometric_pressure
            / (self.ambient_temperature + 273.15)
        )
