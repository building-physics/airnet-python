:html_theme.sidebar_secondary.remove:

Development Environment
=======================

.. _setup_hatch:

Hatch Environment Setup
-----------------------

To set up a development environment using hatch, first install `hatch` with pip:

.. code-block:: console

   $ pip install hatch

Clone the `repository <https://github.com/building-physics/airnet-python/>`_ to the location of your choice and start up a terminal in the
root folder of the repo. Execute the following to generate an environment that has everything that is needed:

.. code-block:: console

   $ hatch env create

To make sure that everything has worked, run

.. code-block:: console

   $ hatch shell

to enter the environment that was created, and then execute

.. code-block:: console

   $ airnetsim --help

You should see the help output from the program.

.. _vscode:

Using Visual Studio Code
------------------------

To point Visual Studio Code at the created environment, find the environment with

.. code-block:: console

   $ hatch run python -c "import sys;print(sys.executable)"

and copy the result. In Visual Studio Code, hit `ctrl-shift-P` to bring up the command palette, select "Python: Select Interpreter", and paste in the result from above. Any warnings (yellow squiqqly underlines) in the source files should go away.