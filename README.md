# diskvert

## Installation

### Numerical codes

Program requires **GCC 6** (or newer), including **gfortran** or other compiler suite (such as Intel Parallel Studio).
Developement package for Lapack, OpenBLAS or other compatible library is required.

```sh
cd build
make
# install system-wide in /usr/local
sudo make install
```

If one has no root access, installation can be made in user's home directory.
Typical location is ``~/.local/bin`` and ``~/.local/lib``, which is achieved by following commands:

```sh
# install for current user only
make install prefix=~/.local
# add the directory to system path
# this line can be also added to .bashrc
export PATH="$PATH:$HOME/.local/bin"
```

Following files will be installed in the directory given by ``prefix`` parameter:

```
bin/diskvert
bin/dv-rad1
bin/dv-mag-rx
bin/dv-mag
bin/dv-alpha-rx
bin/dv-alpha
lib/libdiskvert.so
lib/diskvert/modules/gfortran/alphadisk.mod
lib/diskvert/modules/gfortran/globals.mod
lib/diskvert/modules/gfortran/heatbalance.mod
(and a couple of other .mod files...)
lib/pkgconfig/diskvert.pc
```

### Installing the python package

To install the Python package, a setup script is provided.
Is it advised (but not required) that you use a virtual environment.
Installing the package system-wide is risky and can conflict with your distribution's packages.

```sh
cd ..
virtualenv venv
. venv/bin/activate
cd python
python setup.py install
```

Alternatively, once can install the Python package for the current user only:

```sh
python setup.py install --user
```

### Special builds

During normal use, there is no need for changing the default compiler flags, which provide optimal execution speed.
However, if for some reason different compiler options are needed (such as debugging or a specific architecture), they can be overriden in a following way:
```sh
# check for array bounds and append debug info:
make FFLAGS='-g -fcheck=all'
# optimize more for the current machine only:
make FFLAGS='-O3 -march=native -funsafe-math-optimizations'
```

### Other compiler vendors

Diskvert can be compiled using **icc** and **ifort**.
They provide about 10 percent performance improvement due to more advanced optimization and faster algebra library.
One can either specify the alternative compiler from the command line, or edit the Makefile to make the change permanent.
```sh
make CC=icc FC=ifort FFLAGS='-O2 -xhost'
```
PGI compilers (**pgcc** and **pgf90**) do not work with this code (as of 2018) but hopefully they will soon.

## Reading input files

The library for reading key-value files can be obtained separately from here: [libconfort](https://github.com/gronki/libconfort).

## Usage

### Programs

The recommended way to produce models is to execute ``diskvert`` program, which refers to the most up-to-date version of the model.
(Other programs: ``dv-alpha``, ``dv-alpha-rx``, ``dv-mag`` and ``dv-rad1`` work but are not officially supported.)
After successful build an installation, execution of program does not require changing the source code.
For example, the following command will process the input file ``input.par``, and produce three output files (``model.dat``, ``model.txt`` and ``model.col``):

```sh
cat input.par | diskvert -o model
```

It may be helpful to store the files in archive file (here called ``model.tgz``) to save disk space:

```sh
tar czf model.{tgz,col,dat,txt} && rm model.{col,dat,txt}
```

#### Input files

The structure of the input file consists of key-value pairs (case-sensitive!), for example:

```
mbh 10
mdot 0.1
radius 6.5
# this is a comment
```

The following keywords are allowed, all taking numerical values (required keywords are in **bold**):

 - **``mbh``** is the black hole mass (in solar masses)
 - **``mdot``** is the accretion rate (in units of Eddington rate)
 - **``radius``** is the radius from the center of the black hole (in Schwarzschild radii)
 - **``alpha``**, ``eta`` (default = ``sqrt(alpha)``) and ``nu`` (default = 0) are magnetic parameters (refer to the paper for details)

#### Command-line parameters

 - ``-equilibrium``/``-dyfu`` (default) assumes that radiation and gas temperature are equal
 - ``-compton`` solves heating and cooling balance but only using Compton term
 - ``-corona``/``-balance`` attempts to solve the full heating/cooling balance (will fail to converge) if the instability occurs
 - ``-post-corona`` solves the full heating-cooling balance after relaxation using constant density. It is used to get the solution despite thermal instability (can be used with ``-dyfu`` or ``-compton`` and is ignored with ``-corona``).
 - ``-rel`` includes the relativistic term in Compton scattering
 - ``-quench`` enables switching off MRI (often causes divergence)

### Reading model files

```python
import diskvert as dv
import matplotlib.pyplot as plt
d,p = dv.col2python('disk.dat')
plt.plot(d['z'] / p.zdisk, d['trad'], color = '#6F6F6E')
plt.plot(d['z'] / p.zdisk, d['temp'], color = '#E94626')
plt.show()
```

### Supplementary scripts

```bash
diskvert-plot -show disk.dat
diskvert-plot -tau disk.dat -o disk.png
diskvert-random
diskvert-random | diskvert -compton -post-corona -o disk && diskvert-plot -show -xl disk.dat
```
