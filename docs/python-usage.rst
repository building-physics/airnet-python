Python usage
============

Convenience run
---------------

:func:`airnet.model.run_simulate` provides the shortest single-shot path and returns
the solved :class:`airnet.model.Model`:

.. code-block:: python

   import airnet

   model = airnet.run_simulate(
       "network.txt",
       output_file="results.csv",
       wind_file="WIND",
       ambient_temperature=-5.0,
       wind_speed=8.0,
       wind_direction=337.5,
       relative_tolerance=1.0e-4,
   )
   pressures = {
       name: node.pressure for name, node in model.nodes.items()
   }

Pass ``output_file=None`` when no CSV is wanted. No file is written when
parsing, construction, initialization, or nonlinear solution fails.

Explicit construction and solve
-------------------------------

The lower-level API separates each stage:

.. code-block:: python

   import airnet

   with open("network.txt", encoding="utf-8") as stream:
       records = list(airnet.Reader(stream))

   model = airnet.Model(records, ambient_temperature=20.0)
   model.initialize(maxiter=100)
   iterations = model.air_movement(
       maxiter=100,
       max_subiter=100,
       tolerance=1.0e-8,
       relative_tolerance=1.0e-4,
       relaxation=0.75,
       minimum_iterations=2,
   )

Node pressures are in pascals and link flows in kg/s. ``flow0`` and ``flow1``
are both meaningful for a two-way DOR link; one-way elements return zero in
the second component. Internally, links connected from a fixed node to a
variable node may be reordered for sparse assembly. Consult ``link.flipped``
or use the original input records when an input-oriented report is required.

Failure handling
----------------

Catch :class:`airnet.reader.BadNetworkInput` for syntax problems,
:class:`airnet.model.BadNetwork` for invalid constructed models, and
:class:`airnet.model.SolverError` for numerical failure:

.. code-block:: python

   try:
       model = airnet.run_simulate("network.txt", output_file=None)
   except (airnet.BadNetworkInput, airnet.BadNetwork) as error:
       print(f"Invalid model: {error}")
   except airnet.SolverError as error:
       print(f"Simulation failed: {error}")

Non-convergence is never returned as a successful partial result.
