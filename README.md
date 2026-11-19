# airnet

[![PyPI - Version](https://img.shields.io/pypi/v/airnet.svg)](https://pypi.org/project/airnet)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/airnet.svg)](https://pypi.org/project/airnet)
![Python package](https://github.com/building-physics/airnet-python/actions/workflows/python-package.yml/badge.svg)

-----

**Table of Contents**

- [Installation](#installation)
- [Command line](#command-line)
- [Documentation](#documentation)
- [Scope](#scope)
- [License](#license)

## Installation

```console
pip install airnet
```

## Command line

Run one non-interactive simulation and select its CSV output:

```console
airnet simulate -o results.csv network.txt
```

Use `airnet simulate --help` for weather, wind-profile, convergence, and
relaxation options. `airnet summarize network.txt` reports the contents of an
input file without solving it.

## Documentation

The Sphinx documentation combines narrative material from `docs/` with API
docstrings. From CMD, build the HTML site and PDF manual with:

```bat
hatch run docs:build
```

The results are written to `docs\_build\html` and `docs\_build\pdf`.

## Scope

Version 2 targets AIRNET single-shot calculations. Its regression suite covers
all 81 historical runs using PLR, DWC, QFR, DOR, CFR, FAN, CPF, and CKV
elements. AIRNET's interactive editing is out of scope; MRR and PRV are
deferred to version 3. This package is not independently validated
building-safety software.

## License

`airnet` is distributed under the terms of the [BSD-3-Clause](https://spdx.org/licenses/BSD-3-Clause.html) license.
