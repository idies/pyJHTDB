=======
pyJHTDB
=======

Python wrapper for the JHU Turbulence Database Cluster library.
More information can be found at http://turbulence.pha.jhu.edu/.

Installing
==========

In theory setuptools should take care of everything so that the
package is installed properly, but I noticed that it doesn't.

ubuntu 14.04
------------

Bare-bone installation::

    sudo apt-get install build-essential gfortran
    sudo apt-get install python-setuptools
    sudo apt-get install python-dev
    sudo easy_install numpy
    sudo python setup.py install

Note that doing this should, in principle, also install ``sympy`` on your
system, since it's used by ``pyJHTDB``.

Happy fun installation::

    sudo apt-get install build-essential gfortran
    sudo apt-get install python-setuptools
    sudo apt-get install python-dev
    sudo apt-get install libpng-dev libfreetype6-dev
    sudo apt-get install libhdf5-dev
    sudo easy_install numpy
    sudo easy_install h5py
    sudo easy_install matplotlib
    sudo python setup.py install

Please note that before you run ``sudo python setup.py install``, you
will need to edit the ``turblib-20140606/turblib.c`` file: on line 53,
please add a comma after the closing curly brace, before the ``//``.

I haven't tested the installation on any other system, but I think
reasonable variations on the above should work for the minimal
installation on all unix systems (i.e. for MacOS as well).
If you manage to get it working (i.e. you import test_plain like the
README says and you can run it), please let me know what steps you
needed to take for your system, so I can append the instructions to
this file.

Basic usage
===========

Although this particular Python wrapper is still a work in progress, it
is mature enough to be used in production work.

On first contact with this library, we recommend that you first run
"test_plain". To be more specific:

.. code:: python

    >>> from pyJHTDB import test_plain
    >>> test_plain()

The code that is executed can be found in "pyJHTDB/test.py", and it's
the simplest example of how to access the turbulence database.

