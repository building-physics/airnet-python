:html_theme.sidebar_secondary.remove:

Running
=======

The package provides two programs for use on :ref:`network input files <input_file>`.

.. _airnetsim:

airnetsim
---------

.. code-block:: console

   $ airnetsim --help
   usage: airnetsim [-h] [-v] NETWORK_FILE

   Simulate airflow in an AIRNET network.

   positional arguments:
     NETWORK_FILE

   options:
     -h, --help     show this help message and exit
     -v, --verbose  operate verbosely

.. _summarize_airnet_input:

summarize_airnet_input
----------------------

.. code-block:: console

   $ summarize_airnet_input --help
   usage: summarize_airnet_input [-h] [-v] NETWORK_FILE

   Summarize an AIRNET network input file.

   positional arguments:
     NETWORK_FILE

   options:
     -h, --help     show this help message and exit
     -v, --verbose  operate verbosely
