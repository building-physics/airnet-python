Development and verification
============================

The compatibility reference is the adjacent ``AIRNET-C`` repository. Use
historical reports first, then the NIST equations, then reviewed executable
output. Do not regenerate expected data from Python output.

Running checks
--------------

From CMD, run focused tests first and then the full suite:

.. code-block:: bat

   hatch run +py=312 test:pytest tests\test_afedat.py -q
   hatch run +py=312 test:test
   hatch run +py=311 test:test
   hatch run lint:check

The default Hatch environment also includes pytest and pytest-cov for an
interactive development session:

.. code-block:: bat

   hatch shell
   pytest

The ``pytest`` command in this shell is provided by the Hatch environment and
the checkout is installed as an editable package.

The test matrix in ``pyproject.toml`` also declares later supported Python
versions for CI even when those interpreters are not installed locally.

Continuous integration
----------------------

GitHub Actions runs the complete unit and AIRNET regression suites on CPython
3.11 through 3.15 on Ubuntu, with an additional CPython 3.12 Windows run. A
separate Python 3.11 job exercises the minimum supported dependency versions,
and a canonical Python 3.12 job reports branch coverage. PyPy testing is
deferred pending a reproducible environment that provides a compatible SciPy
binary.

CI also runs Ruff, builds the package distributions, and builds both Sphinx
formats with warnings treated as errors. The wheel, source distribution, and
PDF manual are retained as workflow artifacts. The HTML documentation is
validated but is not uploaded as an artifact.

The release workflow publishes the validated distribution artifacts with PyPI
Trusted Publishing. Repository administrators must configure a protected
``pypi`` GitHub environment and authorize ``python-publish.yml`` in the PyPI
project's trusted-publisher settings before publishing a release.

Regression fixtures
-------------------

Structured fixtures and copied historical inputs live under
``regression-tests\airnet-c``. Each fixture records its source report, run
control, weather, expected values, and tolerances. Multi-run collections are
solved as independent single shots while reproducing AIRNET's documented
prior-pressure property state where required.

The regression harness reconstructs each variable node's mass balance from
the stored link-flow results. A node satisfies this check when either its
case-specific absolute residual threshold is met or its residual divided by
the total absolute incident flow meets the run's relative airflow tolerance.
This absolute-plus-relative check avoids making otherwise converged large-flow
networks depend on sparse-solver rounding details.

``scripts\extract_airnet_reference.py`` converts archived AIRNET reports to
JSON without consulting Python results. AFOUT.133's detailed English-unit
report is converted back to SI with tolerances that account for its printed
precision.

Documentation publication
-------------------------

Narrative reStructuredText is maintained in ``docs``. API pages are generated
from package docstrings by Sphinx autodoc, so public signature and docstring
changes should be checked with a warning-as-error documentation build:

.. code-block:: bat

   hatch run docs:build

Both HTML and PDF are release artifacts. A change is not documentation-clean
if either builder emits a warning.
