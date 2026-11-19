Tested support matrix
=====================

Element support
---------------

.. list-table:: Version 2 element status
   :header-rows: 1
   :widths: 12 14 24 50

   * - Record
     - Status
     - Regression evidence
     - Notes
   * - PLR
     - Supported
     - PL1, PL2, PL3, 04F--133
     - Power-law relationship, including zero and reversed pressure drop.
   * - DWC
     - Supported
     - DW1, DW2
     - Darcy-Weisbach duct with laminar and turbulent regimes.
   * - QFR
     - Supported
     - QF1, QF2
     - Quadratic resistance; QF1 uses the corrected AIRNET-C initialization.
   * - DOR
     - Supported
     - DR1, DR1-01
     - One- and two-way doorway flow with both flow components checked.
   * - CFR
     - Supported
     - CF1, CF2, 133
     - Pressure-independent source flow, control, and orientation.
   * - FAN
     - Supported
     - FN1, FN2, FN3, FN4
     - Legacy curve-range and single-precision cubic behavior retained.
   * - CPF
     - Supported
     - CP1, CP2
     - Constant-pressure fan and link controls.
   * - CKV
     - Supported with provenance caveat
     - CV1
     - Matches recovered AIRNET-C behavior; its equation is not specified in
       the principal NIST AIRNET report.
   * - PRV
     - Deferred
     - PR1, PR2 not run
     - Incomplete internal class, no reader, rejected by the built-in model.
   * - MRR
     - Deferred
     - MR1, MR2 not run
     - Not parsed or represented; deferred to version 3.

Environmental and solver behavior
---------------------------------

.. list-table:: Cross-cutting status
   :header-rows: 1
   :widths: 30 20 50

   * - Feature
     - Status
     - Evidence
   * - Stack pressure
     - Supported
     - ST1, ST2, and full DR1.
   * - Ambient node and weather
     - Supported
     - 133 and the weather-varied regression collections.
   * - 16-direction wind profiles
     - Supported
     - All WP1 and WP2 directions, including interpolation and wraparound.
   * - Link multiplier
     - Supported programmatically
     - Direct unit coverage for every supported solver element.
   * - Link control
     - Supported programmatically
     - FN1--FN4 and CP1--CP2; allowed on CFR, FAN, and CPF.
   * - Explicit non-convergence
     - Supported
     - Initialization, correction, non-finite, and iteration-limit tests.

Historical case coverage
------------------------

The pytest regression harness represents all 81 runs from the historical
reports whose elements are in the version 2 scope:

``04F``, ``08F``, ``16F``, ``32F``, ``133``, ``CF1``, ``CF2``, ``CP1``,
``CP2``, ``CV1``, ``DR1``, ``DW1``, ``DW2``, ``FN1``, ``FN2``, ``FN3``,
``FN4``, ``PL1``, ``PL2``, ``PL3``, ``QF1``, ``QF2``, ``ST1``, ``ST2``,
``WP1``, and ``WP2``.

The additional ``DR1-01`` fixture is a deliberately reduced doorway-only
diagnostic, not another historical case. MR1, MR2, PR1, and PR2 account for
the excluded MRR/PRV reports.

Parity assertions cover pressures, pressure drops, both flow components,
flow direction, convergence bounds, finite values, and mass balance using
case-specific absolute-plus-relative tolerances. A successful parse alone is
not counted as support.
