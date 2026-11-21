Overview
========

An airflow network represents rooms, zones, and the outdoors as pressure
nodes joined by airflow elements. Given fixed boundary pressures, weather,
and element data, ``airnet`` solves for the pressures at variable nodes that
balance mass flow. Link flows are then evaluated from the resulting pressure
drops.

Calculation sequence
--------------------

A single run has five stages:

#. :meth:`airnet.model.Model.read_airnet` parses and validates the legacy
   network without inventing weather that is absent from the AIRNET file.
#. An immutable :class:`airnet.weather.Weather` is assigned to
   :attr:`airnet.model.Model.weather`. This validates wind-profile references
   and establishes weather-dependent physical properties.
#. :class:`airnet.model.Solver` is constructed for the weather-bound model.
   It assigns variable-node equation indices, normalizes link orientation when
   required, and allocates its own sparse numerical workspace.
#. :meth:`airnet.model.Solver.initialize` obtains an initial pressure state
   using each element's linear representation.
#. :meth:`airnet.model.Solver.steady_solve` performs the nonlinear airflow
   solve and either returns an iteration count or raises a solver exception.

The model never reports a failed or non-converged solve as success. CSV output
is written only after the nonlinear convergence criteria are satisfied. Link
results in that output retain the node ordering and sign convention of the
AIRNET input even when the solver normalized a link internally.

Relationship to AIRNET
----------------------

AIRNET-C and its NIST documentation are the compatibility reference. The
Python implementation intentionally differs in software structure:

* elements are objects with ``flow``, ``calculate``, and ``linearize``
  operations;
* network data, environmental inputs, numerical solution, and presentation
  have separate ``Model``, ``Weather``, and ``Solver`` responsibilities;
* invalid references and numerical failures raise explicit exceptions; and
* regression tests compare structured numerical results rather than report
  typography or interactive transcripts.

The version 2 regression matrix covers all 81 historical runs that use the
in-scope elements, plus a reduced doorway-only diagnostic case. See
:doc:`support-matrix` for the precise boundary.
