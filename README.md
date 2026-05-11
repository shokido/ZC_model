# Zebiak-Cane Type intermediate atmosphere-ocean coupled ENSO Model

This repository contains a Fortran implementation of a Zebiak-Cane type
intermediate coupled ocean-atmosphere model for tropical Pacific ENSO
experiments. It includes standalone ocean, standalone atmosphere, and coupled
model drivers, together with example namelists, grid files, forcing files, and
postprocessing notebooks.

The code is intended for research and experimentation. Exact scientific
configuration, validation status, and recommended default experiment should be
confirmed before using the model for publication-quality results.


## Model Overview

The repository provides four main executables:

- `exec_solver_ocn_dyn.out`: standalone ocean dynamical model forced by
  prescribed winds.
- `exec_solver_ocn_full.out`: standalone ocean model including the SST equation,
  forced by prescribed winds.
- `exec_solver_atm.out`: atmosphere/Gill-type response model forced by SST.
- `exec_solver_couple_full.out`: fully coupled ocean-atmosphere model.

The coupled configuration reads an ocean grid, an atmosphere grid, a coupler
mapping file, background ocean fields, SST and wind forcing files, and parameter
settings from a Fortran namelist. The available example namelists are mostly for
equatorial Pacific configurations (`eqpac`) and climatological or
annual-cycle background states.

## Example Visualization
Simulated sea surface temperature and thermocline depth for fully-coupled simulation.
![CGCM SST anomaly animation](GALLERY/cgcm_ssta_thermocline.gif)


## Repository Structure

```text
CODES/
  Fortran source files, Makefile, and generated executables/modules after build.

INPUT/
  Grid and coupling files.
  ATM/       Atmosphere grid files and grid-generation scripts.
  OCN/       Ocean grid, masks, topography, and grid-generation scripts.
  COUPLER/   Ocean-atmosphere coupler mapping files and generation scripts.

Forcing/
  External forcing and background fields, mostly NetCDF.
  OCN/       Ocean background fields and scripts to make mean fields.
  WIND/      Wind forcing files and wind preprocessing scripts.
  SST/       SST climatology/anomaly files and preprocessing scripts.
  Tzm/       Background vertical temperature-gradient fields.

RUN/
  Example namelist files for model execution.
  OGCM/      Standalone ocean runs.
  AGCM/      Standalone atmosphere runs.
  CGCM/      Coupled model runs.
  MIROC6/    Additional experiment files; purpose to be confirmed.
  Proj/      Project-specific coupled examples; purpose to be confirmed.

OUTPUTS/
  Example or generated NetCDF outputs, organized by model component.

GALLERY/
  Plotting notebooks and example figures.
```

## Requirements

Confirmed from the repository:

- A Fortran compiler. The provided `CODES/Makefile` uses `gfortran`.
- NetCDF C and NetCDF Fortran libraries. The Makefile links with
  `-lnetcdf -lnetcdff`.
- `make`.
- Python 3 for preprocessing and analysis scripts.
- Python packages used by scripts/notebooks are to be confirmed, but likely
  include `numpy`, `xarray`, `netCDF4`, `matplotlib`, and `jupyter`.

The Makefile currently contains local Homebrew paths such as:

```make
NETCDF_INCDIR=-I/opt/homebrew/Cellar/netcdf-fortran/4.6.1/include
NETCDF_LIBDIR=-L/opt/homebrew/Cellar/netcdf-fortran/4.6.1/lib -L/opt/homebrew/Cellar/netcdf/4.9.2_2/lib
```

These paths are machine-specific. New users will usually need to edit
`NETCDF_INCDIR` and `NETCDF_LIBDIR` for their environment.

## Example Workflow
### 1. Clone and enter the repository.
```bash
git clone https://github.com/shokido/ZC_model
cd ZC_model
```
### 2. Edit CODES/Makefile if NetCDF include/library paths differ.
```bash
cd CODES
make -f Makefile
differ.
```
### 3. Run the standalone ocean spin-up example.
```bash
cd ..
mkdir -p OUTPUTS/OGCM
cd RUN/OGCM
../../CODES/exec_solver_ocn_dyn.out < do_ogcm_spn_eqpac_30_cd013H120.nml
```
### 4. Run the coupled example.
```bash
cd ../..
mkdir -p OUTPUTS/CGCM
cd RUN/CGCM
../../CODES/exec_solver_couple_full.out < do_cgcm_eqpac_30_ann_c1.4H120_dt3600_c10day.nml
```
### 5. Inspect results with notebooks in GALLERY/.
```bash
cd ../../GALLERY
jupyter lab
```

For a short test run, make a copy of the namelist first:

```bash
cd RUN/CGCM
cp do_cgcm_eqpac_30_ann_c1.4H120_dt3600_c10day.nml test_short.nml
```

Then edit `test_short.nml` to shorten the date range and change all output
filenames to avoid overwriting existing files.

## Build Instructions
From the repository root:
```bash
cd CODES
make -f Makefile
```

Expected executables:

```text
CODES/exec_solver_ocn_dyn.out
CODES/exec_solver_ocn_full.out
CODES/exec_solver_atm.out
CODES/exec_solver_couple_full.out
```

Useful Makefile targets:

```bash
make -f Makefile          # build all executables
make -f Makefile ogcm     # build standalone ocean executables
make -f Makefile agcm     # build standalone atmosphere executable
make -f Makefile clean    # remove generated objects, modules, and executables
```

## Running a Standalone Ocean Simulation

Start with a standalone ocean spin-up before running the coupled model. A
standard ocean example corresponding to the coupled example below is:

```text
RUN/OGCM/do_ogcm_spn_eqpac_30_cd013H120.nml
```

Run it from the repository root with:

```bash
mkdir -p OUTPUTS/OGCM
cd RUN/OGCM
../../CODES/exec_solver_ocn_dyn.out < do_ogcm_spn_eqpac_30_cd013H120.nml
```

This run writes an ocean average file and an ocean restart file under
`OUTPUTS/OGCM/`. See `RUN/OGCM/README.md` for the input and output files used
by this namelist.

The coupled example below is ordered after the standalone ocean example for a
standard workflow, but the current coupled namelist has `flag_ini_ocn="F"` and
does not directly read the OGCM restart file. How best to use OGCM spin-up
output when preparing new coupled experiments is to be confirmed.

The previous README described a workflow in which an OGCM spin-up is followed
by generating mean ocean background fields:

```bash
cd Forcing/OCN
python make_meanfield_ann.py   # fields without annual cycle
python make_meanfield_clm.py   # fields with annual cycle
```

The exact recommended spin-up/mean-field workflow for new experiments is to be
confirmed.

## Running a Standard Coupled Experiment

After the standalone ocean simulation, run the coupled model. A standard coupled
example available in the repository is:

```text
RUN/CGCM/do_cgcm_eqpac_30_ann_c1.4H120_dt3600_c10day.nml
```

This namelist runs an equatorial Pacific coupled configuration without annual-cycle
background/forcing files. Run it from the repository root with:

```bash
mkdir -p OUTPUTS/CGCM
cd RUN/CGCM
../../CODES/exec_solver_couple_full.out < do_cgcm_eqpac_30_ann_c1.4H120_dt3600_c10day.nml
```

The model writes progress information to standard output and NetCDF files to
the paths specified in the namelist. See `RUN/CGCM/README.md` for the input and
output files used by this namelist.

Runtime may be long because this example spans 1970-01-01 to 2020-01-01 with
`dt=3600.0`. For a quick smoke test, copy the namelist, shorten `end_yymmdd`,
and write to new output filenames so existing outputs are not overwritten.

## Standalone Atmospheric Simulation Example

A standalone atmosphere example is:

```bash
mkdir -p OUTPUTS/AGCM
cd RUN/AGCM
../../CODES/exec_solver_atm.out < do_agcm_hindcast_clm_256_ZC87.nml
```

## Required Input Files

Required files are defined by each namelist, not hard-coded globally. Before
running, inspect the selected `RUN/**/*.nml` file and verify that every listed
file exists.

Detailed input lists for the standard examples are documented with the
namelists:

- `RUN/OGCM/README.md`: standalone ocean input files.
- `RUN/CGCM/README.md`: coupled model input files.

In general, standalone ocean runs require an ocean grid file from `INPUT/OCN/`
and wind files from `Forcing/WIND/`. Coupled runs additionally require an
atmosphere grid, a coupler mapping file, background ocean fields, atmosphere SST
forcing, and ocean background/forcing files.

Some NetCDF forcing and output files in this working tree may be generated data
rather than version-controlled source files. Confirm which large files should be
shared with collaborators.

## Output Files

Output filenames are configured in each namelist. The standard standalone ocean
example writes `avg_*` and `rst_*` files under `OUTPUTS/OGCM/`; the standard
coupled example writes ocean average, ocean diagnostic, atmosphere average, and
ocean restart files under `OUTPUTS/CGCM/`.

Detailed output lists for the standard examples are documented in
`RUN/OGCM/README.md` and `RUN/CGCM/README.md`.

The meaning and units of all output variables should be confirmed from the
Fortran source and experiment documentation before scientific interpretation.

## Governing equations
### Overview
This model is a Fortran implementation of a Zebiak--Cane type intermediate
coupled ocean--atmosphere model. The ocean component consists of two parts:
a dynamical model for the upper-ocean current and thermocline-depth anomalies,
and a thermodynamic model for the sea surface temperature (SST) anomaly. The
dynamical ocean model combines a reduced-gravity shallow-water model, which
describes the low-order baroclinic response to wind stress, with a diagnostic
Ekman model representing higher vertical-mode contributions. The atmosphere
component is based on the steady Gill-type response of the first baroclinic
mode of the tropical atmosphere to diabatic heating.

### Ocean shallow-water model
The reduced-gravity ocean dynamics are governed by
```math
\frac{\partial u_{\rm SW}}{\partial t}=fv_{\rm SW}-g' \frac{\partial h}{\partial x}+\frac{\tau^x_0}{\rho_0 (H_1+H_2)}
-r u_{\rm SW}+
\nu\left(\frac{\partial^2}{\partial x^2}+\frac{\partial^2}{\partial y^2} \right) u_{\rm SW},
```
```math
\frac{\partial v_{\rm SW}}{\partial t}=-f u_{\rm SW}-g' \frac{\partial h}{\partial y}+\frac{\tau^y_0}{\rho_0 (H_1+H_2)}-r v_{\rm SW}+\nu \left(\frac{\partial^2}{\partial x^2}+\frac{\partial^2}{\partial y^2}\right)v_{\rm SW},
```
```math
\frac{\partial h}{\partial t}=-(H_1+H_2)\left(\frac{\partial u_{\rm SW}}{\partial x}+\frac{\partial v_{\rm SW}}{\partial y}\right)-rh .
```

Here, $u_{\rm SW}$ and $v_{\rm SW}$ are the zonal and meridional
components of the vertically averaged upper-ocean current associated with the
shallow-water mode, and $h$ is the thermocline-depth anomaly. The variables
$\tau^x_0$ and $\tau^y_0$ denote the zonal and meridional wind stress, $H_{1}$ is the mixed-layer depth, $H_{2}$ is the thickness between the mixed
layer and the reference thermocline depth, $g'$ is the reduced gravity,$\rho_{0}$ is the reference density, $r$ is the linear drag coefficient, and $\nu$ is the horizontal viscosity coefficient.

The shallow-water current is related to the layer currents by
```math
u_{\rm SW}=\frac{H_1 u_{o1} + H_2 u_{o2}}{H_1+H_2},
\qquad v_{\rm SW} =
\frac{H_1 v_{o1} + H_2 v_{o2}}{H_1+H_2},
```
where $(u_{o1},v_{o1})$ and $(u_{o2},v_{o2}$ denote the currents in the mixed layer and in the layer between the mixed layer and the reference
thermocline depth, respectively.

In the original Zebiak and Cane (1987) model, the meridional momentum tendency
is neglected and the ocean response is decomposed into equatorial Kelvin and
Rossby wave components. This approach filters out inertia--gravity waves and
allows a relatively long time step, but it is less convenient for complex
land--sea masks. In the present implementation, the shallow-water equations are
instead directly discretized. A horizontal diffusion term is added to suppress
the growth of high-frequency and short-wavelength noise.

The equations are discretized on an Arakawa C grid in space and integrated
using a leapfrog scheme in time. An Asselin filter is applied to suppress the
computational mode. The boundary condition is no-normal flow at the western,
eastern, northern, southern, and land boundaries. For the tangential component
along the boundary, either a no-slip condition,
```math
u_l = 0,
```
or a free-slip condition,
```math
\frac{\partial u_l}{\partial n} = 0,
```
can be selected using `slip_ind`.
The corresponding implementation is mainly found in

```text
mod_ocn_solver_zc.f90: solve_rg_vgeo_ocn(ogrd, oset, dt)
```
### Ocean Ekman model
The vertical shear between the mixed layer and the layer below is diagnosed
using a linear Ekman balance. Defining
```math
u_{\rm EK} = u_{o1} - u_{o2},
\qquad
v_{\rm EK} = v_{o1} - v_{o2},
```
The Ekman equations are
```math
\epsilon_{so} u_{\rm EK}-f v_{\rm EK}=\frac{\tau^x_0}{\rho_0 H_1},
```
```math
\epsilon_{so} v_{\rm EK}+f u_{\rm EK}=\frac{\tau^y_0}{\rho_0 H_1},
```
where $\epsilon_{so}$ is the Ekman-layer damping coefficient. Solving these
diagnostic equations gives
```math
u_{\rm EK}=\frac{
\epsilon_{so}\tau^x_0 + f \tau^y_0
}{
\rho_0 H_1 (f^2 + \epsilon_{so}^2 )
}
```
```math
v_{\rm EK}
=
\frac{
\epsilon_{so}\tau^y_0 - f \tau^x_0
}{
\rho_0 H_1 (f^2 + \epsilon_{so}^2 )
}.
```
The mixed-layer current is then obtained by combining the shallow-water and
Ekman components:
```math
u_{o1}=u_{\rm SW}+\frac{H_2}{H_1+H_2} u_{\rm EK},
```
```math
v_{o1}=v_{\rm SW}+\frac{H_2}{H_1+H_2} v_{\rm EK}.
```
The vertical velocity at the base of the mixed layer is diagnosed from
continuity as
```math
w_{o}=H_{1}
\left(\frac{\partial u_{o1}}{\partial x}+\frac{\partial v_{o1}}{\partial y}
\right).
```
These calculations are implemented in
```text
mod_ocn_solver_zc.f90: solve_ekman_ocn(ogrd, oset)
mod_ocn_solver_zc.f90: solve_totalcurrent_ocn(ogrd, oset)
```

### SST anomaly equation
The SST anomaly equation describes the evolution of the mixed-layer temperature
anomaly $T$. In the present model, it can be written schematically as
```math
\begin{aligned}
\frac{\partial T'}{\partial t}
=-u'_{o1}\frac{\partial \overline{T}}{\partial x}
-\overline{u}_{o1}\frac{\partial T'}{\partial x}
-u'_{o1}\frac{\partial T'}{\partial x}
-v'_{o1}\frac{\partial \overline{T}}{\partial y}
-\overline{v}_{o1}\frac{\partial T'}{\partial y}
-v'_{o1}\frac{\partial T'}{\partial y}
\\
&
-
M(\overline{w}_o+w'_o)
\frac{\partial \overline{T}}{\partial z}
-
M(\overline{w}_o+w'_o)
\frac{T' - T'_e}{H_1}
-
\alpha_s T' .
\end{aligned}
```
Here, overbars denote prescribed background fields and primes denote anomalies.
The first two lines represent horizontal temperature advection by the mean and
anomalous mixed-layer currents. The third line represents vertical entrainment
and upwelling effects at the base of the mixed layer, together with a linear
thermal damping term. The function $M(x)$ is the upwelling operator, usually
defined as
```math
M(x)=
\begin{cases}
x, & x > 0, \\
0, & x \le 0,
\end{cases}
```
so that only upward motion contributes to entrainment cooling. The parameter$\alpha_s$ is the SST damping coefficient, and $T_{e}^{\prime}$ denotes the
temperature anomaly entrained from below the mixed layer.
The SST equation is implemented in the ocean thermodynamic component, mainly in
```text
mod_ocn_solver_zc.f90: get_tsub_ZC(oset,h,hbar,sst,Tsub,Tent)
mod_ocn_solver_zc.f90: solve_sst_ocn_ZC(ogrd,oset,mudata,mvdata,mwdata,mhdata,mTzdata,msstdata,dt)
```

### Atmospheric Gill model
The atmospheric component represents the steady response of the first baroclinic mode of the tropical atmosphere to heating, following the Gill-type model. The governing equations are
```math
\epsilon_a u_a-\beta_0 y v_a=-\frac{\partial \phi}{\partial x},
```
```math
\epsilon_a v_a+\beta_0 y u_a=-\frac{\partial \phi}{\partial y},
```
```math
\epsilon_a \phi+c_a^2\left(\frac{\partial u_a}{\partial x}+\frac{\partial v_a}{\partial y}\right)=
-Q_0 .
```
Here, $u_{a}$ and $v_{a}$ are the zonal and meridional surface wind responses,
$\phi$ is the geopotential height perturbation, $c_{a}$ is the phase speed of
the first baroclinic mode, $\epsilon_{a}$ is the atmospheric damping coefficient, and $Q_{0}$ is the prescribed heating associated with SST
anomalies.
Following Zebiak (1982), these equations are solved in spectral space after Fourier transformation. The corresponding implementation is found in
```text
mod_atm_solver_gill.f90
```

## Known Limitations and Caveats

- The Makefile contains machine-specific NetCDF paths.
- There is no automated test suite currently documented.
- The exact validated/default experiment is to be confirmed.
- Units and scientific interpretation of several parameters and output
  variables should be confirmed from the source and/or publications.
- Reproducibility may depend on compiler, NetCDF library versions, and exact
  forcing files.

## References and Citation

Core reference:

- Zebiak, S. E., and M. A. Cane, 1987: A Model El Nino-Southern Oscillation.
  *Monthly Weather Review*, 115, 2262-2278.

If you use this repository in a publication or shared analysis, also cite this
code repository and record the commit hash, namelist file, compiler, NetCDF
versions, and input forcing files used.
