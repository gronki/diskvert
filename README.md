# diskvert

## Installation

### Obtain code

Remember to use ``--recursive`` option while downloading the repo:

```sh
git clone --recursive --depth 1 https://github.com/gronki/diskvert.git
```

### Docker

The easiest way is to build a Docker container (only do it once):

```sh
cd diskvert
docker build -t diskvert .
```

Then run the container with:

```sh
docker run -it -v "$PWD":/work diskvert
```

where ``$PWD`` can be replaced with path to computation directory.

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

During normal use, there is no need for changing the default compiler flags, which provide best execution speed.
If you encounter any errors, you might want to build with debug flags that will cause slow execution but will provide debug info and checks:
```sh
# check for array bounds and append debug info:
make FFLAGS='-g -fcheck=all'
```
Provided Makefile respects ``ARCH_FLAGS`` environment variable.
You might want to export it if your execution environment has some special needs, such as computing clusters and ARM (by default it will be ``-march=native``).
For example:
```sh
# pro tip: add export to rc file
# and you never have to type it again
export ARCH_FLAGS='-msse4.2 -mavx'
make
```

### Other compiler vendors

Diskvert can be compiled using **icc** and **ifort**.
They provide about 10 percent performance improvement due to more advanced optimization and faster algebra library.
One can either specify the alternative compiler from the command line, or edit the Makefile to make the change permanent.
```sh
make CC=icc FC=ifort FFLAGS='-O2 -xhost'
```
PGI compilers (**pgcc** and **pgf90**) do not work with this code (as of 2018) but hopefully they will upgrade their compilers soon.

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
alpha 0.02
# this is a comment
```

The following keywords are allowed, all taking numerical values:

 - ``mbh`` is the black hole mass (in solar masses)
 - ``mdot`` is the accretion rate (in units of Eddington rate)
 - ``radius`` is the radius from the center of the black hole (in Schwarzschild radii)
 - ``alpha``, ``eta`` and ``nu`` are magnetic parameters ($\eta$ and $\nu$ are not required, refer to the paper for details and relations)

#### Command-line parameters

Note: every switch can be negated by adding ``-no`` prefix (for example: ``-no-fluxcorr``).

 - ``-equilibrium``/``-dyfu`` (default) assumes that radiation and gas temperature are equal
 - ``-compton`` solves heating and cooling balance but only using Compton term
 - ``-corona``/``-balance`` attempts to solve the full heating/cooling balance (will fail to converge if the instability occurs)
 - ``-post-corona`` solves the full heating-cooling balance after relaxation using constant density. It is used to get the solution despite thermal instability (can be used with ``-dyfu`` or ``-compton`` and is ignored with ``-corona``).
 - ``-rel`` includes the relativistic term in Compton scattering
 - ``-quench`` enables switching off MRI (often causes divergence)
 - ``-fluxcorr`` (default) accounts for correction of flux for very strongly magnetized disks
 - ``-klnish`` enables Klein-Nishina cross section (does not converge, obviously)
 - ``-alpha-prad`` (default) includes radiation pressure in alpha prescription

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
diskvert-random | diskvert -compton -post-corona -o disk \
	&& diskvert-plot -show -xl disk.dat
```
