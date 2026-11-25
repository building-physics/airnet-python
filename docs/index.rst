airnet documentation
====================

``airnet`` is an object-oriented Python implementation of the NIST AIRNET
bulk-airflow model. Its version 2 scope is single-shot network simulation:
read a network, initialize it, solve for steady pressures and mass flows, and
return or export the converged state.

The implementation is designed for understandable experiments and regression
testing while preserving AIRNET's documented equations, units, link
orientation, and two-way-flow semantics. It does not reproduce AIRNET's
interactive editor or repeated stdin-driven session.

.. important::

   ``airnet`` is engineering calculation software, but it is not an
   independently validated building-safety application. Use the
   :doc:`support-matrix` to determine what has actually been tested.

.. toctree::
   :maxdepth: 2
   :caption: User guide

   overview
   installation
   command-line
   input-format
   python-usage
   support-matrix

.. toctree::
   :maxdepth: 2
   :caption: Reference and development

   api
   development
