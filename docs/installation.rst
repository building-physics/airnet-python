Installation
============

Install the published package with pip:

.. code-block:: bat

   py -m pip install airnet

The package requires Python 3.11 or newer and installs the ``airnet`` command.

Development checkout
--------------------

The project uses Hatch to create isolated test, lint, and documentation
environments:

.. code-block:: bat

   py -m pip install hatch
   hatch run +py=312 test:test

Use ``312`` rather than ``3.12`` in Hatch matrix selection. Select another
installed interpreter with ``+py=311``, for example.

Building the documentation
--------------------------

Build both publication formats from the repository root:

.. code-block:: bat

   hatch run docs:build

The outputs are:

* ``docs\_build\html\index.html`` for the HTML site; and
* ``docs\_build\pdf\airnet.pdf`` for the PDF manual.

Individual formats can be built with ``hatch run docs:html`` and
``hatch run docs:pdf``. The ``docs\make.bat`` wrapper provides the same
``html``, ``pdf``, and ``all`` targets.
