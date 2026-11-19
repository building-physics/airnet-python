Command-line interface
======================

The installed entry point is a Click command group with two subcommands:

.. code-block:: bat

   airnet --help
   airnet summarize INPUT_FILE
   airnet simulate INPUT_FILE

``summarize``
-------------

``summarize`` parses a network file and reports its title and the number of
nodes, links, and elements of each type. It does not construct or solve the
model.

.. code-block:: bat

   airnet summarize regression-tests\airnet-c\AFDATA.DW1

Use ``--verbose`` to show file-processing messages. Invalid records are
reported as command-line errors with a nonzero exit status.

``simulate``
------------

``simulate`` performs one complete initialization and nonlinear solve. By
default it writes ``airnetsim.csv`` in the current directory:

.. code-block:: bat

   airnet simulate examples\law-office-simple.txt
   airnet simulate -o results.csv examples\law-office-simple.txt

Weather and wind options use AIRNET SI units:

.. code-block:: text

   airnet simulate --wind-file WIND --ambient-temperature -5 ^
     --barometric-pressure 101325 --wind-speed 8 ^
     --wind-direction 337.5 -o windy.csv AFDATA.WP2

``--wind-direction`` is degrees clockwise from north. A wind file contains
named profiles of 16 coefficients at 22.5-degree intervals.

Solver options
~~~~~~~~~~~~~~

The principal numerical options are:

``--initialization-iterations``
   Maximum conjugate-gradient iterations for linear initialization.

``--maximum-iterations``
   Maximum nonlinear airflow iterations.

``--maximum-linear-iterations``
   Maximum conjugate-gradient iterations for each nonlinear correction.

``--absolute-tolerance`` and ``--relative-tolerance``
   Per-node mass-flow convergence criteria in kg/s and as a fraction of the
   absolute flow through the node. A node satisfies either enabled criterion.

``--relaxation``
   AIRNET's conditional pressure-correction factor, greater than zero and no
   greater than one.

Use ``airnet simulate --help`` for defaults and the complete option list. A
parse, validation, initialization, convergence, or output-write failure exits
nonzero and does not leave a newly written result file.
