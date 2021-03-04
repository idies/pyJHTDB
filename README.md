# pyJHTDB

Python wrapper for the [JHU Turbulence Database Cluster](http://turbulence.pha.jhu.edu/) library.
More information can be found at [http://turbulence.pha.jhu.edu/help/python/](http://turbulence.pha.jhu.edu/help/python/).

## Use pyJHTDB through SciServer  (RECOMMENDED)
The SciServer is a cloud-based data-driven cluster, of The Institute for Data Intensive Engineering and Science (IDIES) at Johns Hopkins University. Users get the advantages of more reliable and faster data access since the SciServer is directly connected to JHTDB through a 10 Gigabit ethernet connection. SciServer provides docker containers with the pyJHTDB library pre-installed.

To use pyJHTDB through Sciserver:
```
Login to [SciServer](http://turbulence.pha.jhu.edu/) (may need to create a new account first).
Click on *Compute* and then *Create container* (You could also run jobs in batch mode, by selecting Compute Jobs).
Type in *Container name*, select *JH Turbulence DB* in *Compute Image* and then click on *Create*.
Click on the container you just created, then you could start using pyJHTDB with Python or IPython Notebook.
```
Please go to [SciServer](http://turbulence.pha.jhu.edu/) for more information on SciServer as well as the help on Sciserver.

## Use pyJHTDB on local computers

Prerequisites: numpy>=1.15.0, scipy>=1.1.0, sympy>=1.2, h5py>=2.8.0, matplotlib>=3.0.0, 
More prerequisites since 01/08/2021: wurlitzer>=2.0.2 

### Installing pypi version (Linux and MacOS)

If you have *pip*, you can simply do this:
```
pip install pyJHTDB
```
If you're running unix (i.e. some MacOS or GNU/Linux variant), you will probably need to have a `sudo` in front of the `pip` command. If you don't have `pip` on your system, it is quite easy to get it following the instructions at [http://pip.readthedocs.org/en/latest/installing.html](http://pip.readthedocs.org/en/latest/installing.html).

### Installing from source (Linux and MacOS)

```
git clone https://github.com/idies/pyJHTDB.git
cd pyJHTDB
python update_turblib.py
pip install --upgrade ./
```

### In Windows

We notice the compilation error in Windows, and provide a temperary workaround in `examples\Use_JHTDB_in_windows.ipynb`. This method use [zeep](https://python-zeep.readthedocs.io/en/master/) package. Should be as convient as pyJHTDB.

## Basic usage

On first contact with this library, we recommend that you first run
``test_plain``. To be more specific:
```
from pyJHTDB import test_plain
test_plain()
```

The code that is executed can be found in *pyJHTDB/test.py*, and it's the simplest example of how to access the turbulence database.

## Configuration

While our service is open to anyone, we would like to keep track of who is using the service, and how. To this end, we would like each user or site to obtain an authorization token from us: [http://turbulence.pha.jhu.edu/help/authtoken.aspx](http://turbulence.pha.jhu.edu/help/authtoken.aspx)

For simple experimentation, the default token included in the package should be valid.
