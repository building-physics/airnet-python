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

   model = airnet.Model.read_airnet("network.txt")
   model.weather = airnet.Weather(ambient_temperature=20.0)
   solver = airnet.Solver(model)
   solver.initialize(maxiter=100)
   iterations = solver.steady_solve(
       maxiter=100,
       max_subiter=100,
       tolerance=1.0e-8,
       relative_tolerance=1.0e-4,
       relaxation=0.75,
       minimum_iterations=2,
   )

AIRNET keeps its network, wind coefficients, and run weather separate.
Consequently, :meth:`airnet.model.Model.read_airnet` always returns a model
whose ``weather`` is ``None``. When wind profiles are needed, read the separate
file and include them in the weather before binding it:

.. code-block:: python

   with open("WIND", encoding="utf-8") as stream:
       profiles = airnet.read_wind_profiles(stream)

   model.weather = airnet.Weather(
       wind_speed=8.0,
       wind_direction=337.5,
       wind_profiles=profiles,
   )

Changing ``model.weather`` invalidates an existing solver. Construct a new
solver for each weather scenario.

Node pressures are in pascals and link flows in kg/s. ``flow0`` and ``flow1``
are both meaningful for a two-way DOR link; one-way elements return zero in
the second component. Internally, links connected from a fixed node to a
variable node may be reordered for sparse assembly. Consult ``link.flipped``
or use the original input records when an input-oriented report is required.

Failure handling
----------------

Catch :class:`airnet.reader.BadNetworkInput` for syntax problems,
:class:`airnet.weather.BadWeather` for invalid environmental inputs,
:class:`airnet.model.BadNetwork` for invalid constructed models or missing
weather, and :class:`airnet.model.SolverError` for numerical failure or a
solver made stale by rebinding weather:

.. code-block:: python

   try:
       model = airnet.run_simulate("network.txt", output_file=None)
   except (airnet.BadNetworkInput, airnet.BadWeather, airnet.BadNetwork) as error:
       print(f"Invalid model: {error}")
   except airnet.SolverError as error:
       print(f"Simulation failed: {error}")

Non-convergence is never returned as a successful partial result.
