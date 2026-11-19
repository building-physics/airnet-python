# AIRNET-C regression fixtures

These fixtures capture reviewed, single-shot scenarios from the historical
AIRNET-C tests. Python regression tests do not execute AIRNET-C, so they remain
portable and deterministic in CI.

Each JSON file identifies:

- the copied AIRNET network input;
- the immutable AIRNET-C `AFOUT.*` report used as its source;
- which run in that report is represented;
- the AIRNET run control and weather values needed to interpret the result;
- structured node and link values transcribed at the report's printed
  precision; and
- comparison tolerances selected to account for that printed precision.

Interactive changes in an AIRNET run script are represented as separate
single-shot fixtures named `<case>-run-<number>.json`, such as
`pl1-run-01.json`. Do not update expected values from Python output. New or
changed expected values must be checked against the historical report or a
reviewed run of AIRNET-C.

Complete multi-run reports are stored as `<case>-sweep.json` or
`<case>-runs.json` collections. The parameter changes made between runs are
declared in `regression-tests/test_airnet_c.py`, so each report run is still
constructed and solved independently. The version 2.0.0 matrix contains all
81 historical runs that use the documented, supported elements. MR1, MR2,
PR1, and PR2 are intentionally deferred because MRR and PRV are outside that
release scope.

Some repeated AIRNET runs establish air properties from the preceding node
pressure state before performing a new linear initialization. Their
independent fixtures record that state in `initial_node_pressures`; the test
harness uses it only to establish properties and then resets variable-node
pressure guesses. Controlled links are recorded explicitly in
`link_controls`.

`AFDATA.DR1-01` is an explicitly simplified, doorway-only derivative of the
historical input. It retains the two nodes, DOR element, and `link-11`, but
omits the ten PLR links. It is not the full historical DR1 input, and its
expected equilibrium comes from a reviewed AIRNET-C run of the reduced case.

`AFDATA.DR1` is the full historical DR1 network. The distinct `-01` suffix
keeps the focused doorway-only regression available without displacing the
canonical historical name. `WIND` is a copy of the historical 16-direction
coefficient file used by WP1 and WP2.

The development utility `scripts/extract_airnet_reference.py` converts the
immutable reports to structured JSON and never reads Python solver output.
For example, from the repository root in CMD:

```bat
python scripts\extract_airnet_reference.py ..\AIRNET-C\c\tests\reference\AFOUT.WP1 WP1 AFDATA.WP1 regression-tests\airnet-c\wp1-sweep.json
```

QF1 is the exception to the archived-report rule: its collection records a
reviewed run of the rebuilt executable after correction of the QFR linear
initialization, and all 13 runs converge in two iterations. The immutable
historical report is not modified. AFOUT.133 uses AIRNET's detailed English
output rather than its normal metric summary; the extractor converts its
0.001 in. H2O pressure and 0.1 SCFM flow fields back to SI, and its tolerances
account for that printed precision.
