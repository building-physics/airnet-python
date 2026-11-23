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

Regression fixtures
-------------------

Structured fixtures and copied historical inputs live under
``regression-tests\airnet-c``. Each fixture records its source report, run
control, weather, expected values, and tolerances. Multi-run collections are
solved as independent single shots while reproducing AIRNET's documented
prior-pressure property state where required.

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
