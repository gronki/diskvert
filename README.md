# diskvert

> **If you use this code in your research, please cite:**
>
> - Gronkiewicz D., Różańska A., 2020, *Warm and thick corona for a magnetically supported disk in galactic black hole binaries*, A&A 633, A35, [doi:10.1051/0004-6361/201935033](https://doi.org/10.1051/0004-6361/201935033)
> - Gronkiewicz D., Różańska A., Petrucci P.-O., Belmont R., 2023, *Thermal instability as a constraint for warm X-ray coronas in active galactic nuclei*, A&A 675, A198, [doi:10.1051/0004-6361/202244410](https://doi.org/10.1051/0004-6361/202244410)
>
> Machine-readable metadata is in [`CITATION.cff`](CITATION.cff).

`diskvert` computes the **1D vertical structure of an accretion disk** around a black hole at a single radius, from the equatorial plane ($z=0$) up to the top of the atmosphere.
The disk can be supported by magnetic pressure, and heated by magnetic dissipation, which allows a warm corona to form self-consistently above the disk.
It is meant for both galactic black hole binaries and AGN.

The code solves a set of coupled ODEs (hydrostatic equilibrium, radiative diffusion, magnetic field evolution, energy transport and, optionally, the local heating–cooling balance) as a boundary-value problem, using a **relaxation (Newton) method** on a grid.
All physics is described in the papers above; this README tells you how to run the code and how it is put together.

Contents:

1. [Getting started](#getting-started)
2. [Running a model](#running-a-model)
3. [Output files](#output-files)
4. [Reading results in Python](#reading-results-in-python)
5. [How the code works](#how-the-code-works)
6. [Changing the equations (`coefficients.py`)](#changing-the-equations-coefficientspy)
7. [Python bindings to the Fortran library](#python-bindings-to-the-fortran-library)
8. [Repository map](#repository-map)

---

## Getting started

There are four ways to get `diskvert`:

| | for whom |
|---|---|
| [prebuilt AppImage](#prebuilt-releases) | you only want to run models |
| [native build](#native-build) | development on your own machine: Fortran compiler + OpenBLAS + a current Python |
| [conda](#conda) | as native, but with everything (incl. the old Python tools) in one environment |
| [Docker / Apptainer](#docker-and-apptainer) | a ready Ubuntu 20.04 environment with the Python tools; Apptainer for clusters |

The Python tools for plotting and reading output (`python/`) have not been maintained since 2021 and need old package versions (Python 3.8, see `python/requirements.txt`); conda or the containers provide them.
`sympy/coefficients.py`, which generates the equations, runs on current Python (tested with Python 3.14 + sympy 1.14 and Python 3.8 + sympy 1.3).

To build from source, clone the repository **with** `--recursive`; the input-file parser [libconfort](https://github.com/gronki/libconfort) is a git submodule and the build fails without it:

```sh
git clone --recursive https://github.com/gronki/diskvert.git
cd diskvert
```

### Prebuilt releases

Each tagged version on the [Releases page](https://github.com/gronki/diskvert/releases) provides two files for Linux x86_64:

- `diskvert-x86_64.AppImage`: only the `diskvert` program, no installation needed (glibc 2.31 or newer):
  `chmod +x diskvert-x86_64.AppImage && ./diskvert-x86_64.AppImage -corona -o mydisk < input.txt`
- `diskvert-x86_64.sif`: the full Apptainer image, including the Python tools:
  `apptainer exec diskvert-x86_64.sif diskvert -corona -o mydisk < input.txt`

They are built by `.github/workflows/release.yml`; the AppImage can be built locally with `scripts/build_appimage.sh` (see the header of that script).

### Native build

Requires `make`, `gfortran` (GCC 6 or newer) and OpenBLAS (the Makefile links `-lopenblas`). To regenerate the equations you also need Python 3 with `sympy` and `numpy`.

```sh
make -C build
make -C build install                    # into /usr/local if writable, otherwise ~/.local
make -C build install prefix=/opt/somewhere
```

This installs the programs into `bin/`, `libdiskvert.so` into `lib/`, Fortran modules into `lib/diskvert/modules/` and a `pkg-config` file. Make sure `bin/` is on your `PATH` (for `~/.local`: `export PATH="$PATH:$HOME/.local/bin"`).

### conda

```sh
conda env create -f environment.yml   # once
conda activate diskvert
make -C build install prefix=$CONDA_PREFIX
pip install ./python                  # Python tools (optional)
```

The programs are installed into the environment, so `diskvert` is on the `PATH` while it is active. `mamba`/`micromamba` work the same.

### Docker and Apptainer

Both images are based on Ubuntu 20.04 and contain the compiler and the old Python tools. The repository on the host is mounted in the container, so you edit the code on the host as usual.

**Docker:**

```sh
./run_docker.sh
```

This builds the image `diskvert` and starts an interactive shell. Directories are shared with the host:

| host               | container      | purpose                                             |
|--------------------|----------------|-----------------------------------------------------|
| `.` (repository)   | `/source`      | source code                                         |
| `work/`            | `/work`        | your working directory (shell starts here)          |
| `build/dist/`      | `/opt/diskvert`| installed binaries and `libdiskvert.so`             |

On first start, `work/` is created and filled with the files from `examples/`.

**Apptainer** (e.g. on a cluster):

```sh
./build_apptainer.sh   # builds diskvert.sif (uses sudo)
./run_apptainer.sh     # starts a shell, same mounts as Docker except work/
```

Instead of building, you can download `diskvert-x86_64.sif` from the releases and save it as `diskvert.sif` in the repository root.

In both, `build/dist/` from the host is mounted over `/opt/diskvert`, which hides the binaries compiled while building the image. **On first start run `rebuild.sh`** (`make` in `/source/build` + install into `/opt/diskvert`). The installed files stay on the host in `build/dist/`, so later you only need to rebuild after changing the code.

### First test

```sh
cat examples/input.txt | diskvert -corona -o mydisk   # compute
diskvert-plot mydisk.dat                              # writes mydisk.png (needs the Python tools)
```

### Development loop

After editing `sympy/coefficients.py` (or any Fortran source), regenerate the equations and rebuild:

| environment | regenerate + rebuild |
|---|---|
| native | `(cd sympy && python3 coefficients.py) && make -C build install` |
| conda | `(cd sympy && python coefficients.py) && make -C build install prefix=$CONDA_PREFIX` |
| Docker / Apptainer shell | `coefficients.sh && rebuild.sh` |

Skip the regeneration step if you changed only Fortran. `make` recompiles only what changed, so a rebuild takes seconds. With containers, you can also run it from the host without entering the shell: `./run_docker.sh bash -c 'coefficients.sh && rebuild.sh'` (the same works with `./run_apptainer.sh`).

---

## Running a model

`diskvert` reads the parameters from **standard input** and options from the **command line**:

```sh
cat input.txt | diskvert [options] -o NAME
```

### Input file

Key–value pairs, one per line, keys are case-sensitive, `#` starts a comment:

```
mbh    10      # black hole mass [M_sun]
mdot   0.05    # accretion rate [Eddington units]
radius 10      # distance from the black hole [Schwarzschild radii]
alpha  0.1     # magnetic alpha (field generation by dynamo)
eta    0.05    # magnetic eta (buoyant rise speed of the field)
nu     1       # magnetic nu (reconnection rate)
```

| key      | required | meaning |
|----------|----------|---------|
| `mbh`    | yes | black hole mass in $M_\odot$ |
| `mdot`   | yes | accretion rate in units of $\dot M_{\rm Edd} = 4\pi G M / (c\,\kappa_{\rm es}\,\eta_{\rm acc})$, with efficiency $\eta_{\rm acc}=0.083$ (`globals.f90`) |
| `radius` | yes | radius in units of $R_{\rm schw} = 2GM/c^2$; the inner edge of the disk is fixed at $3\,R_{\rm schw}$ |
| `alpha`  | yes | $\alpha_B$: field generation (dynamo) parameter |
| `eta`    | no  | $\eta_B$: field rise velocity $v_B = \eta_B\,\Omega z$; default $0.5\sqrt{\alpha/0.5}$ |
| `nu`     | no  | $\nu_B$: reconnection parameter; default $0.1/\alpha$ |
| `output` | no  | output base name (same as `-o`) |

On start the program prints derived quantities: $\beta_0 = (2\eta + \alpha(\nu-1))/\alpha$ (expected plasma beta in the midplane) and $q_{\rm cor} = (2\eta + \alpha(\nu-1))/\eta$ (magnetic pressure falls off as $P_{\rm mag}\propto z^{-q_{\rm cor}}$ high above the disk). The program stops if $\beta_0<0$ or $q_{\rm cor}\le 0$. See the 2020 paper for the meaning of the magnetic parameters.

Chemical composition (X, Z), the mean molecular weight ($\mu=0.5$) and the opacity coefficients are hard-coded in `src/globals.f90`.

`diskvert-random` writes a random but sensible input to standard output (options `-mbh`, `-mdot`, `-radius`, `-qcor`; `-mbh agn` or `-mbh xrb` draws a typical mass):

```sh
diskvert-random -mbh agn -qcor 3 | diskvert -corona -o disk
```

### Choosing the physics: temperature treatment

This is the most important choice. Exactly one of:

| option | equations solved | |
|---|---|---|
| `-dyfu`, `-equilibrium` | gas temperature = radiation temperature (no corona) | **default** |
| `-compton` | separate gas temperature from heating–cooling balance, **Compton cooling only** | |
| `-corona`, `-balance` | separate gas temperature from the **full** heating–cooling balance (Compton + free–free) | |

The full balance (`-corona`) does not converge when the atmosphere is thermally unstable (see the 2023 paper). In that case one can use:

- `-post-corona`: compute the structure first (e.g. with `-dyfu` or `-compton`), then solve the full heating–cooling balance once at fixed density. This breaks the hydrostatic equilibrium slightly.
- `-iter-post-corona`: as above, but iterated a few times.

### Other options

Every on/off switch has a negated form with `-no-`, e.g. `-no-fluxcorr`.

| option | default | meaning |
|---|---|---|
| `-o NAME` | `disk` | base name of the output files |
| `-alpha` | off | non-magnetic $\alpha$-disk (no corona is possible in this mode) |
| `-quench` | **on** | switch off MRI (and thus magnetic field generation) where the field exceeds the maximum allowed by MRI; may cause divergence |
| `-alpha-prad` | **on** | include radiation pressure in the pressure used by the $\alpha$ prescription |
| `-fluxcorr` | **on** | add the energy flux escaping above the top of the grid (relevant for strongly magnetized disks) |
| `-rel` | **on** | relativistic correction to Compton cooling |
| `-klnish` | off | Klein–Nishina cross-section (also turns on `-rel`); usually does not converge |
| `-ff`, `-bf` | on | free–free and bound–free absorption opacity |
| `-simple-hydro` | off | in hydrostatic equilibrium, use the radiation temperature gradient instead of the gas temperature gradient (used by default in `dv.run`) |
| `-n N` | 2048 | number of grid points (300–4000, rounded to a multiple of 16; increased by ⅓ for the corona step) |
| `-htop H` | auto | top of the grid in units of $H$ (see `zscale` below); by default estimated from the parameters |
| `-log`, `-linear`, `-linlog`, `-asinh`, `-pow2` | `-log` | grid type |
| `-trim` | off | cut the top of the grid where density becomes negligible |
| `-estim-full`, `-estim-compt`, `-estim-instab`, `-no-estim` | automatic | how to guess the gas temperature before the corona step |
| `-all` | off | write every iteration to `NAME.NNN.dat` (useful to watch the relaxation) |
| `-conduction` | off | thermal conduction; experimental, marked as not working in `coefficients.py` |

### Did it work?

The terminal shows one line per iteration with the relative error and "ramps" (see [How the code works](#how-the-code-works)). A successful run ends with `SUCCESS`; otherwise you see `diverged` or `not converged`. The output files are **written anyway**, so always check `converged` in `NAME.txt`.

To run many models, use e.g. GNU `parallel` (installed in the Docker/Apptainer images); see `tests/fluxcorr/run.sh` for an example.

---

## Output files

Each run produces three text files:

| file | content |
|---|---|
| `NAME.dat` | the vertical structure: one row per grid point |
| `NAME.col` | names and types of the columns of `.dat` (read by the Python tools) |
| `NAME.txt` | scalar quantities: input, derived parameters, convergence flags, values at characteristic heights |

To save space you can pack them; the Python reader opens archives directly:

```sh
tar czf NAME.tgz NAME.{col,dat,txt} && rm NAME.{col,dat,txt}
```

All values are in **cgs units**.

### Columns of `NAME.dat`

The first three columns are the point index `i`, the height `z` [cm] and `h` = $z/H$, where $H$ = `zscale` = $\sqrt{2kT_{\rm eff}/m_{\rm H}}/\Omega$.
Most useful remaining columns:

| column | meaning |
|---|---|
| `rho` | density [g cm⁻³] |
| `temp`, `trad` | gas and radiation temperature [K] |
| `pgas`, `prad`, `pmag` | gas, radiation, magnetic pressure [dyn cm⁻²] |
| `frad`, `fmag`, `fcnd` | radiative, magnetic (Poynting, $2P_{\rm mag}v_B$), conductive flux [erg cm⁻² s⁻¹] |
| `heat` | heating rate [erg cm⁻³ s⁻¹]; `heatm` and `heatr` are its two magnetic contributions |
| `cool`, `coolb`, `coolc` | total, free–free and Compton cooling rate [erg cm⁻³ s⁻¹]; `compfr` = Compton fraction |
| `coolerr` | (heat − cool)/heat: how well the local energy balance is satisfied |
| `tau`, `taues`, `tauth` | optical depth from the top: absorption+scattering, scattering only, effective (thermalization) |
| `tavg` | absorption-weighted average gas temperature above this point [K] |
| `coldens` | column density from the midplane [g cm⁻²]; `nhtot`: hydrogen column from the top [cm⁻²] |
| `beta` | plasma beta $P_{\rm gas}/P_{\rm mag}$ |
| `qcor` | local slope $-d\ln P_{\rm mag}/d\ln z$ |
| `qmri` | MRI quenching factor (1 = MRI active) |
| `instabp` | thermal instability criterion at constant pressure; **negative = unstable** |
| `ksct`, `kabs`, `kabp` | scattering, absorption (incl. scattering), Planck opacity [cm² g⁻¹] |
| `compy` | Compton $y$ parameter integrated from the top |
| `ionxi` | $P_{\rm rad}/P_{\rm gas}$ |
| `ts_cool` | cooling time $P_{\rm gas}$/`cool` [s] |

The remaining columns are diagnostics; their definitions are in the subroutine `fillcols` in `src/prog/diskvert.f90`.

### Contents of `NAME.txt`

Lines of `key value # comment`. Important keys:

- `converged`, `has_corona`, `has_instability`, `model`, `niter`, `ngrid`;
- input and derived parameters: `mbh`, `mdot`, `radius`, `alpha`, `eta`, `nu`, `qcor`, `beta_0`, `facc` (flux produced by accretion), `teff`, `omega`, `zscale`, `rschw`;
- `zdisk`/`hdisk`: density-weighted disk thickness [cm] / [$H$]; `coldens`, `tavgr`;
- values at characteristic heights, with a suffix: `_midpl` (midplane), `_phot` ($\tau=1$), `_therm` ($\tau^\star=1$), `_tmin` (temperature minimum), `_cor` (corona base = max of the two previous), `_instabp_lo/_hi/_md` (extent of the unstable zone). For each there are e.g. `zphot`, `hphot`, `temp_phot`, `rho_phot`, `tavg_phot`, `taues_phot`, …

---

## Reading results in Python

```python
import diskvert as dv
import matplotlib.pyplot as plt

d, p = dv.col2python('mydisk.dat')   # also accepts mydisk.tgz / .tar.gz
# d: numpy structured array, columns by name (names are lower-case, e.g. 'cool_dt')
# p: object with all keys of mydisk.txt as attributes

plt.plot(d['z'] / p.zdisk, d['trad'], label='T_rad')
plt.plot(d['z'] / p.zdisk, d['temp'], label='T_gas')
plt.yscale('log')
plt.legend()
plt.show()
```

Helpers:

- `dv.run(pars, options=['-corona', '-simple-hydro'])` runs `diskvert` in a temporary directory for a `dict` of input parameters and returns `(d, p)`. Raises an exception with the log if the program fails.
- `dv.standard_pars_q(q)` returns `(alpha, eta, nu)` for a desired `qcor` (one particular parametrization, see `python/diskvert/random.py`).
- `diskvert-plot NAME.dat` draws an overview figure. The x axis is `-col` (column density, default), `-xl` (log height), `-x` (linear height), `-tau` or `-grid`; `-show` opens a window, `-o file.png` sets the output.

---

## How the code works

### Unknowns and equations

At every grid point the unknowns are (depending on the model) $\rho$, $T$, $T_{\rm rad}$, $F_{\rm rad}$, $P_{\rm mag}$ and $F_{\rm cond}$. The equations, as written in `sympy/coefficients.py`, are:

- **hydrostatic equilibrium**: $\dfrac{dP_{\rm gas}}{dz} + \dfrac{dP_{\rm mag}}{dz} = \dfrac{\kappa\rho}{c}F_{\rm rad} - \Omega^2\rho z$
- **radiative diffusion**: $16\sigma T_{\rm rad}^3\dfrac{dT_{\rm rad}}{dz} = -3\kappa\rho F_{\rm rad}$
- **magnetic field**: $\eta z\dfrac{dP_{\rm mag}}{dz} = \alpha P_{\rm gen} - (2\eta + \alpha\nu)P_{\rm mag}$, where $P_{\rm gen} = P_{\rm gas} + P_{\rm rad} + P_{\rm mag}$ (with default options)
- **energy**: $\dfrac{dF_{\rm rad}}{dz} = \mathcal{H}$, with heating $\mathcal{H} = 2(\eta + \alpha\nu)\Omega P_{\rm mag} - \alpha\Omega P_{\rm gen}$
- **heating = cooling** (only with `-compton`/`-corona`): $\mathcal{H} = 4\sigma\rho\left[\kappa_P(T^4 - T_{\rm rad}^4) + \kappa_{\rm es}T_{\rm rad}^4\dfrac{4k}{m_ec^2}(T - T_{\rm rad})\right]$, where the free–free term is dropped for `-compton`

Boundary conditions: at $z=0$, $F_{\rm rad}=0$ and the magnetic equation without the derivative; at the top, $F_{\rm rad} + F_{\rm mag} = F_{\rm acc}$ (plus the `-fluxcorr` term) and $F_{\rm rad} = 2\sigma T_{\rm rad}^4$.

In the non-magnetic mode (`-alpha`) $P_{\rm mag}=0$ and the heating is $\alpha\Omega P$. The `-quench`, `-rel`, `-alpha-prad` and `-fluxcorr` switches modify these expressions; see `coefficients.py` for the exact form.

### Model numbers

Each combination of physics is a separate "model" with its own set of generated routines. Names and numbers (printed as `model` in `NAME.txt`):

| nr | name | temperature | magnetic | conduction |
|---|---|---|---|---|
| 1 | `ad` | $T=T_{\rm rad}$ | no | no |
| 5 | `md` | $T=T_{\rm rad}$ | yes | no |
| 6 | `mw` | Compton balance | yes | no |
| 8 | `mc` | full balance | yes | no |
| 13 | `mdt` | $T=T_{\rm rad}$ | yes | yes |
| 16 | `mct` | full balance | yes | yes |

The number is $1 + [\text{Compton}] + 3[\text{full balance}] + 4[\text{magnetic}] + 8[\text{conduction}]$ (`mrx_number` in `src/relaxation.f90`).

### Relaxation scheme

1. Global quantities at the given radius ($\Omega$, $F_{\rm acc}$, $T_{\rm eff}$, $H$) are computed in `cylinder` (`globals.f90`), and a Shakura–Sunyaev solution (`ss73solution.f90`) gives the initial guess.
2. The ODEs are discretized on the grid (central differences between neighbouring points). For an approximate solution $Y$, `mrx_matrix` (`relaxation.f90`) builds the residuals $A$ and the Jacobian $M$. The linear system $M\,\delta Y = A$ is banded and is solved with LAPACK `dgbsv`.
3. The correction is applied with a damping factor ("ramp" in the terminal output): $Y \leftarrow Y + \text{ramp}\cdot\delta Y$, where ramp < 1 when the error is large.
4. For stability, the absorption opacity is switched on gradually (`R_OPA` column, `opacities_kill` in the code) and so is MRI quenching (`R_MRI`, `qmri_kill`).
5. First the model with $T=T_{\rm rad}$ is converged. If a corona was requested, the solution is copied to the model with a separate gas temperature (`mrx_transfer`), the grid is refined, the gas temperature is estimated (`-estim-*`) and relaxation continues.
6. Finally, the output columns are computed (`fillcols`) and files are written.

The main program `src/prog/diskvert.f90` drives all of this; most physics lives in the generated file `src/mrxcoeff.fi`.

---

## Changing the equations (`coefficients.py`)

The residuals and their partial derivatives (the Jacobian) are **not written by hand**. `sympy/coefficients.py` defines the equations symbolically, differentiates them with SymPy and prints Fortran code into `src/`:

| generated file | content |
|---|---|
| `mrxcoeff.fi` | for each model: `mrx_coeff1_*` (1st-order ODEs: residual `A`, Jacobians `MY` w.r.t. values and `MD` w.r.t. derivatives), `mrx_coeff0_*` (algebraic equations such as heating = cooling), `mrx_coeffbl_*` / `mrx_coeffbr_*` (boundary conditions at the midplane / top), `mrx_output_*` (output quantities) |
| `mrxdims.fi` | numbers of variables, equations and boundary conditions per model |
| `mrxptrs.fi` | selects the routines above for a model number |
| `mrxhash.fi` | position of each variable ($\rho$, $T$, $T_{\rm rad}$, $F_{\rm rad}$, $P_{\rm mag}$, $F_{\rm cond}$) in the state vector, 0 if absent |
| `mrxfout.fi` | selects the output routine |

These files are included into `src/relaxation.f90`. Things to know before editing:

- Opacities are abstract functions `fkabs`, `fkabp`, `fksct`, `fkcnd` of $(\rho, T)$. In the generated code they are replaced by the matrix `F(3,4)`: row 1 is the value, rows 2 and 3 the derivatives with respect to $\rho$ and $T$; columns are absorption, Planck, scattering, conduction. The Fortran side fills `F` from `kappabs`, `kappabp`, `kappesp`, `kappcnd` in `globals.f90`. To change opacities, edit those Fortran functions (value **and** derivatives); no regeneration is needed.
- Run-time switches (`use_prad_in_alpha`, `use_relcompt`, `use_quench_mri`, `use_flux_correction`) become `merge(...)` expressions in Fortran, so they can be changed without regenerating.
- The list of models is `choices` at the top of the loop.
- The order of the list `yout` must match the constants `c_rho`, `c_temp`, … in `relaxation.f90`, and `yvar_all` must match `cc_rho`, …. If you add an output, update `N_RELAX_OUTPUTS` too.
- The script checks that the number of equations equals the number of unknowns and that there is one boundary condition per first-order equation.

To regenerate and rebuild, see [Development loop](#development-loop). The script must be run from the `sympy/` directory because it writes to `../src/`. Commit the regenerated `.fi` files together with the change in `coefficients.py`.

Different sympy versions format the generated code differently (order of terms, line breaks), so regenerating with another version changes many lines of `mrxcoeff.fi` even if the equations did not change. The code is equivalent: models computed with files from sympy 1.3 and 1.14 were bit-identical. To keep git diffs readable, regenerate with the same sympy version as the committed files (currently 1.3, as pinned in `python/requirements.txt`).

---

## Python bindings to the Fortran library

`make` also builds `libdiskvert.so`. `python/diskvert/models.py` loads it with `ctypes` (the directory with the library must be on `LD_LIBRARY_PATH`; the containers set this) and exposes the relaxation routines, so the Jacobian can be inspected from Python:

| Python | Fortran (`relaxation_c` in `relaxation.f90`) | purpose |
|---|---|---|
| `mrx_init(mbh, mdot, r, alpha, eta, nu)` | `mrx_init_c` | set the disk parameters |
| `mrx_model(balance, magnetic, conduction)` | `mrx_number_c` | model number; `balance` is `b'D'`, `b'W'` (Compton) or `b'C'` (full) |
| `mrx_get_ny(nr)` | `mrx_get_ny_c` | number of variables per grid point |
| `mrx_matrix(nr, z, nz, Y, ny, A, na, M)` | `mrx_matrix_c` | fill residuals `A` and Jacobian `M` (Fortran order) |
| `mrx_advance(nr, z, nz, Y, dY, ny)` | `mrx_advance_c` | solve $M\,\delta Y = A$, returns the error code; `dY` receives the correction |

The state vector `Y` has length `nz*ny` and is ordered point by point: `[ρ₁, T₁, …, ρ₂, T₂, …]`; the order of variables within a point is given by `mrxhash.fi`. Run-time switches and the opacity ramp keep their module defaults (all on), i.e. the binding does not reproduce the gradual start used by `diskvert`.

`examples/draw_matrix.py` builds an initial guess and plots the matrix structure. These bindings are incomplete and serve mostly for inspection; to compute models, run the `diskvert` program.

---

## Repository map

```
src/prog/diskvert.f90     main program: options, relaxation loop, output
src/relaxation.f90        building the Jacobian, model selection, C bindings
src/mrx*.fi               generated by sympy/coefficients.py — do not edit by hand
src/globals.f90           constants, composition, opacities, Ω, F_acc, T_eff
src/balance.f90           local heating–cooling balance solver (corona estimate, -post-corona)
src/ss73solution.f90      Shakura–Sunyaev solution (initial guess)
src/util/settings.f90     reading input keys and common command-line options
sympy/coefficients.py     symbolic equations → Fortran code
python/diskvert/          reading output, running models, ctypes bindings
python/scripts/           diskvert-plot and other scripts
examples/                 example input files (copied to work/)
build/Makefile            build and install
libconfort/               submodule: key–value input parser
```

The programs `dv-alpha`, `dv-mag`, `dv-rad1`, `dv-ss73` and `dv-opactbl` are legacy tools built alongside `diskvert`; they are not maintained.
