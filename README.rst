
=======
pyJHTDB
=======

Python wrapper for the JHU Turbulence Database Cluster library.
More information can be found at http://turbulence.pha.jhu.edu/.

Installing pypi version (RECOMMENDED)
=======================

If you have ``pip``, you can simply do this:

.. code:: bash

    pip install pyJHTDB

If you're running unix (i.e. some MacOS or GNU/Linux variant), you will
probably need to have a ``sudo`` in front of the ``pip`` command.
If you don't have ``pip`` on your system, it is quite easy to get it
following the instructions at
http://pip.readthedocs.org/en/latest/installing.html.

Installing from source
======================

**ubuntu 14.04**

Bare-bone installation:

.. code:: bash

    git clone https://github.com/idies/pyJHTDB.git
    cd pyJHTDB
    python update_turblib.py
    pip install --upgrade --force-reinstall ./

The procedures are similar in MacOS.

Basic usage
===========

On first contact with this library, we recommend that you first run
``test_plain``. To be more specific:

.. code:: python

    from pyJHTDB import test_plain
    test_plain()

The code that is executed can be found in "pyJHTDB/test.py", and it's
the simplest example of how to access the turbulence database.

Configuration
=============

While our service is open to anyone, we would like to keep track of who
is using the service, and how. To this end, we would like each user or
site to obtain an authorization token from us:
http://turbulence.pha.jhu.edu/help/authtoken.aspx
For simple experimentation, the default token included in the package
should be valid.

