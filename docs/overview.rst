Overview
========

An airflow network represents rooms, zones, and the outdoors as pressure
nodes joined by airflow elements. Given fixed boundary pressures, weather,
and element data, ``airnet`` solves for the pressures at variable nodes that
balance mass flow. Link flows are then evaluated from the resulting pressure
drops.

Calculation sequence
--------------------

A single run has four stages:

#. :class:`airnet.reader.Reader` parses the legacy network records.
#. :class:`airnet.model.Model` validates references and constructs nodes, elements,
   links, environmental pressure terms, and the sparse system.
#. :meth:`airnet.model.Model.initialize` obtains an initial pressure state using
   each element's linear representation.
#. :meth:`airnet.model.Model.air_movement` performs the nonlinear airflow solve and
   either returns an iteration count or raises a solver exception.

The model never reports a failed or non-converged solve as success. CSV output
is written only after the nonlinear convergence criteria are satisfied.

Relationship to AIRNET
----------------------

AIRNET-C and its NIST documentation are the compatibility reference. The
Python implementation intentionally differs in software structure:

* elements are objects with ``flow``, ``calculate``, and ``linearize``
  operations;
* parsing, model construction, numerical solution, and presentation are
  separate;
* invalid references and numerical failures raise explicit exceptions; and
* regression tests compare structured numerical results rather than report
  typography or interactive transcripts.

The version 2 regression matrix covers all 81 historical runs that use the
in-scope elements, plus a reduced doorway-only diagnostic case. See
:doc:`support-matrix` for the precise boundary.
