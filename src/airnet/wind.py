# SPDX-FileCopyrightText: 2023-present Oak Ridge National Laboratory, managed by UT-Battelle
#
# SPDX-License-Identifier: BSD-3-Clause
from __future__ import annotations

import math
from dataclasses import dataclass

from .reader import BadNetworkInput


@dataclass(frozen=True)
class WindProfile:
    """AIRNET wind-pressure coefficients at 22.5-degree intervals."""

    name: str
    coefficients: tuple[float, ...]

    def __post_init__(self):
        coefficients = tuple(float(value) for value in self.coefficients)
        if len(coefficients) != 16:
            raise ValueError('A wind profile must contain exactly 16 coefficients')
        if not all(math.isfinite(value) for value in coefficients):
            raise ValueError('Wind-profile coefficients must be finite')
        if not self.name or self.name == 'null':
            raise ValueError('A wind profile must have a non-reserved name')
        object.__setattr__(self, 'coefficients', coefficients)

    def coefficient(self, direction_degrees: float) -> float:
        """Linearly interpolate a coefficient using AIRNET's 16 sectors."""
        direction = direction_degrees % 360.0
        lower = int(direction / 22.5)
        upper = (lower + 1) % 16
        fraction = (direction - 22.5 * lower) / 22.5
        return (
            (1.0 - fraction) * self.coefficients[lower]
            + fraction * self.coefficients[upper]
        )


def read_wind_profiles(lines) -> dict[str, WindProfile]:
    """Read an AIRNET wind-coefficient file into named profiles."""
    profiles = {}
    tokens = []
    for line_number, line in enumerate(lines, start=1):
        content = line.split('!', 1)[0].strip()
        if not content or content.startswith('/'):
            continue
        data = content.split()
        if data[0] == 'title':
            continue
        tokens.extend((line_number, token) for token in data)

    position = 0
    while position < len(tokens):
        line_number, name = tokens[position]
        position += 1
        if position + 16 > len(tokens):
            raise BadNetworkInput(
                f'Wind profile "{name}" at line {line_number} has fewer than 16 coefficients'
            )
        values = []
        for coefficient_line, token in tokens[position:position + 16]:
            try:
                values.append(float(token))
            except ValueError as error:
                raise BadNetworkInput(
                    f'Wind profile "{name}" has invalid coefficient "{token}" '
                    f'at line {coefficient_line}'
                ) from error
        position += 16
        if name in profiles:
            raise BadNetworkInput(f'Wind profile "{name}" is defined more than once')
        try:
            profiles[name] = WindProfile(name, tuple(values))
        except ValueError as error:
            raise BadNetworkInput(str(error)) from error
    return profiles
