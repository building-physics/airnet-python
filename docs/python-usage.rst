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
Internally, this convenience function follows the same stages as the explicit
API: it calls ``Model.read_airnet``, binds one ``Weather`` object, constructs
one ``Solver``, initializes, and performs the steady solve. It returns that
same solved model rather than maintaining a separate result object.

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

Initialization and steady solution are operations on ``Solver``. ``Model``
does not provide duplicate solving entry points; it holds the physical network
and the pressure and flow state updated by the solver.

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
the second component. Constructing a solver assigns equation indices to
variable nodes and may reorder a fixed-to-variable link on the model for
sparse assembly. ``link.flipped`` records that change. Stored ``pdrop``,
``flow0``, and ``flow1`` use the solver orientation; their
``input_pressure_drop``, ``input_flow0``, and ``input_flow1`` counterparts use
the original AIRNET input order. CSV output always uses the input order.

The sparse matrix, residual vector, and pressure synchronization cache are
solver state. They are not allocated while reading or constructing a model.
``Solver`` determines variable-node ordering by assigning each node's
``index``; that index remains on the model so the relationship between model
state and matrix rows is explicit without a duplicate node-reference mapping
inside the solver. Constructing another solver rebuilds independent numerical
workspace and leaves an already normalized link unchanged.

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
