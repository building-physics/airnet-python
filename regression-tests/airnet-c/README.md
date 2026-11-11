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
