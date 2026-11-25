Network input format
====================

``airnet`` reads the record-oriented text format used by AIRNET. Fields are
whitespace separated. A line beginning with ``*`` terminates the input; blank
lines and unrecognized comment records are skipped. The package preserves SI
units used by AIRNET-C.

Core records
------------

Title
~~~~~

.. code-block:: text

   title Example pressure network

Only one title record is accepted.

Nodes
~~~~~

.. code-block:: text

   node NAME v HEIGHT TEMPERATURE
   node NAME c HEIGHT TEMPERATURE PRESSURE
   node NAME a HEIGHT TEMPERATURE PRESSURE

``v`` denotes a variable-pressure node, ``c`` a fixed-pressure node, and ``a``
an ambient fixed-pressure node. Heights are metres, temperatures are degrees
Celsius, and pressures are pascals. The run's ambient temperature replaces
the temperature field of an ``a`` node.

Links
~~~~~

.. code-block:: text

   link NAME NODE1 HEIGHT1 NODE2 HEIGHT2 ELEMENT null
   link NAME NODE1 HEIGHT1 NODE2 HEIGHT2 ELEMENT WIND_PROFILE WIND_MODIFIER

Attachment heights are relative to their nodes. Link order defines the
positive input direction. The wind modifier scales the dynamic pressure and
the coefficient interpolated from the selected wind profile. The network
reader retains ``WIND_PROFILE`` as a name because AIRNET stores coefficient
definitions in a separate file. Assigning :class:`airnet.weather.Weather` to
the model resolves and validates every referenced profile.

Model construction preserves this input order. Solver construction may
reverse a fixed-to-variable link to simplify equation assembly and records the
change in ``link.flipped``. The link's ``input_pressure_drop``,
``input_flow0``, and ``input_flow1`` properties, and all CSV output, restore
the original input direction.

Elements
--------

The reader accepts these documented element records:

``plr``
   ``element NAME plr INIT LAM TURB EXPONENT``

``dwc``
   Two lines containing length, hydraulic diameter, area, roughness, then
   turbulent dynamic loss, laminar friction loss, laminar dynamic loss, and
   the initialization coefficient.

``qfr``
   ``element NAME qfr A B``

``dor``
   A PLR first line followed by minimum temperature difference, opening
   height, width, and discharge coefficient.

``cfr``
   ``element NAME cfr MASS_FLOW``

``fan``
   A PLR first line, fan reference data, and one or more five-field curve
   ranges in the historical AIRNET format.

``cpf``
   ``element NAME cpf OPERATING_PRESSURE MINIMUM_RISE FLOW_TYPE``

``ckv``
   ``element NAME ckv CUTOFF_PRESSURE COEFFICIENT``

The exact field names exposed by the parser are documented under
:class:`airnet.reader.Reader`. Examples in ``regression-tests\airnet-c`` are copied
from historical AIRNET cases and exercise every supported file element.

Programmatic-only link fields
-----------------------------

When constructing records in Python, a link may additionally contain a
positive ``mult`` value and a ``control`` value from zero through one. Control
is supported for CFR, FAN, and CPF elements. These extensions are solver
features but are not currently parsed from the legacy network-file link line.

Unsupported records
-------------------

MRR is not parsed or solved. PRV has incomplete internal code but is rejected
by the built-in model and has no reader. Both are deferred to version 3. The
interactive AIRNET command stream, online editing, repeated runs, and report
format controls are not network-file features in ``airnet``.
