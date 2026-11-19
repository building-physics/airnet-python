Python API reference
====================

This section combines the narrative guide with docstrings from the installed
Python package.

Public package
--------------

.. automodule:: airnet
   :members:
   :undoc-members:
   :show-inheritance:

Model and solver
----------------

.. automodule:: airnet.model
   :members: Node, Link, Model, BadNetwork, DuplicateNetworkName, MissingNetworkReference, UnsupportedNetworkFeature, SolverError, LinearSolverError, InitializationError, NonlinearSolverError, NumericalSolverError, write_results_csv, run_simulate
   :show-inheritance:

Input reader
------------

.. automodule:: airnet.reader
   :members: BadNetworkInput, InputType, Reader
   :show-inheritance:

Airflow elements
----------------

.. automodule:: airnet.afedat
   :members: Afe_Plr, Afe_Dwc, Afe_Qfr, Afe_Dor, Afe_Cfr, Afe_Fan, Afe_Cpf, Afe_Ckv
   :show-inheritance:

Wind profiles
-------------

.. automodule:: airnet.wind
   :members: WindProfile, read_wind_profiles
   :show-inheritance:
