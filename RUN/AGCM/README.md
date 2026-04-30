# RUN/AGCM

This directory contains namelist files for standalone atmosphere model
experiments. The default AGCM hindcast example documented here is:

```text
do_agcm_hindcast_clm_256_ZC87.nml
```

This namelist is used with `CODES/exec_solver_atm.out`, which is built from
`CODES/zc_agcm.f90`.

## How to Run

From the repository root:

```bash
mkdir -p OUTPUTS/AGCM
cd RUN/AGCM
../../CODES/exec_solver_atm.out < do_agcm_hindcast_clm_256_ZC87.nml
```

The run writes progress information to standard output and writes NetCDF output
to the path specified by `fname_avg_atm`.

## Standard Namelist Summary

For `do_agcm_hindcast_clm_256_ZC87.nml`:

- Time step: `dt=86400.0` seconds, or 1 day.
- Start date: `20100101`.
- End date: `20160101`.
- Intended period: 2010-01 through 2015-12.
- Atmosphere grid: `../../INPUT/ATM/grid_256.nc`.
- Heating formulation: `aset%heating_type="ZC87"`.
- SST climatology is cyclic with `Lcycle_atm_sstm="T"` and
  `Tcycle_atm_sstm=365`.
- SST anomaly forcing is non-cyclic with `Lcycle_atm_ssta="F"`.

The exact scientific meaning and units of all physical parameters are to be
confirmed before publication use.

## Required Input Files

Required by `do_agcm_hindcast_clm_256_ZC87.nml`:

- `../../INPUT/ATM/grid_256.nc`
- `../../Forcing/SST/ERSST_v5_atm_256_clm.nc`
- `../../Forcing/SST/ERSST_v5_atm_256_anm_2010_2015.nc`

Variables read from these files include:

- Atmosphere-grid SST climatology: `sst`.
- Atmosphere-grid SST anomaly: `sst`.

## Output Files

The default namelist writes:

- Atmosphere average:
  `../../OUTPUTS/AGCM/avg_agcm_hindcast_256_atm_clm_ZC87.nc`

The atmosphere average file contains fields such as `q_atm`, `p_atm`, `u_atm`,
`v_atm`, `ssta_atm`, and `sstm_atm`, plus coordinate variables. The output
attributes in `mod_io_avg.f90` describe `q_atm` and `p_atm` as pressure-like
fields and `u_atm`/`v_atm` as velocity fields; exact interpretation should be
confirmed before scientific analysis.

## Namelist Variables for `zc_agcm.f90`

`zc_agcm.f90` runs the atmosphere response model forced by prescribed SST
climatology and SST anomalies. Depending on `aset%heating_type`, it uses one of
the Gill-type response routines in `CODES/mod_atm_solver_gill.f90`.

### Read Order

The namelist is read from standard input in this order:

1. `&date`
2. `&io_atm`
3. `&param_atm`
4. `&sstm_param_atm`
5. `&sstm_io_atm`
6. `&ssta_param_atm`
7. `&ssta_io_atm`

If `aset%heating_type=="ZC87_conv"`, the code then also reads:

8. `&uam_param_atm`
9. `&uam_io_atm`
10. `&vam_param_atm`
11. `&vam_io_atm`

Keep this order when creating new namelist files.

### `&date`

| Variable | Type | Meaning in `zc_agcm.f90` | Notes |
| --- | --- | --- | --- |
| `dt` | real | Model time step in seconds. | Used to compute total step count and output schedule. |
| `start_yymmdd` | integer | Simulation start date in `YYYYMMDD` format. | Example: `20100101`. |
| `start_hhmmss` | integer | Simulation start time in `HHMMSS` format. | Example: `0` for 00:00:00. |
| `end_yymmdd` | integer | Simulation end date in `YYYYMMDD` format. | `20160101` gives a 2010-01 to 2015-12 hindcast period. |
| `end_hhmmss` | integer | Simulation end time in `HHMMSS` format. | Example: `0` for 00:00:00. |

### `&io_atm`

| Variable | Type | Meaning in `zc_agcm.f90` | Notes |
| --- | --- | --- | --- |
| `fname_grd_atm` | character | Path to the atmosphere grid NetCDF file. | Read by `read_atm_grd`; coordinates are prepared by `set_coord_atm`. |
| `fname_avg_atm` | character | Path to the atmosphere average-output NetCDF file. | Created by `create_avg_atm_ZC`. |
| `out_avg_flag` | integer | Time unit used for average-output scheduling. | See output-time-unit table below. |
| `out_avg_int` | integer | Output interval in units of `out_avg_flag`. | Example: `out_avg_flag=100`, `out_avg_int=1` gives one output per calendar month. |

Output-time-unit flags used by the calendar routines:

| Flag value | Unit |
| --- | --- |
| `-10000` | seconds |
| `-100` | minutes |
| `-1` | hours |
| `1` | days |
| `100` | months |
| `10000` | years |

### `&param_atm`

`param_atm` contains one variable, `aset`, a derived type defined as
`type(atm_set)` in `CODES/run_types.f90`.

| Variable | Default in `run_types.f90` | Meaning in `zc_agcm.f90` | Notes |
| --- | --- | --- | --- |
| `aset%cp_atm` | `60.0` | Atmospheric gravity-wave speed parameter in m/s. | Source comment gives units of m/s; used in atmosphere coordinate scaling and Gill-type response routines. |
| `aset%eps_s_atm_day` | `1.106364457323311` | Atmosphere damping timescale parameter in days. | Used to form nondimensional damping in atmosphere response routines. |
| `aset%alpha_gill_atm` | `0.031` | SST-to-heating coefficient for the Gill-type atmosphere. | Source comment gives units of `m^2*s^(-3)*K^(-1)`. |
| `aset%heating_type` | `"ZC87"` | Selects atmosphere heating/response formulation. | Code branches for `"ZC87"`, `"ZC87_conv"`, and `"GJ22"`. |

Available heating branches in the current source:

- `"ZC87"`: calls `return_uvp_fromSST_ZC87`.
- `"ZC87_conv"`: calls `return_uvp_fromSST_ZC87_conv`, enables atmospheric
  convergence feedback, and requires additional atmosphere mean-wind input
  groups.
- `"GJ22"`: calls `return_uvp_fromSST_GJ22`, an alternative heating
  parameterization proposed by Geng and Jin (2022).

Relevant references:

- Zebiak, S. E., 1986: Atmospheric Convergence Feedback in a Simple Model for
  El Nino. *Monthly Weather Review*, 114, 1263-1271.
  https://doi.org/10.1175/1520-0493(1986)114<1263:ACFIAS>2.0.CO;2.
- Geng, L., and F.-F. Jin, 2022: ENSO Diversity Simulated in a Revised
  Cane-Zebiak Model. *Frontiers in Earth Science*, 10, 899323.
  https://doi.org/10.3389/feart.2022.899323.

### SST Input Groups

Both SST input groups use the same pattern:

- `nfile_*`: number of NetCDF files to read.
- `timename_*`: name of the time coordinate variable.
- `varname_*`: name of the physical variable to read.
- `Lcycle_*`: cyclic-input switch. `"T"` uses day-of-year cyclic interpolation;
  other values use absolute model-time interpolation/clamping.
- `Tcycle_*`: cycle length, usually `365` days for annual-cycle climatology.
- `fnames_*(1:nfile_*)`: paths to the NetCDF files.

| Namelist groups | Variables | Field read by the code | Role in `zc_agcm.f90` |
| --- | --- | --- | --- |
| `&sstm_param_atm`, `&sstm_io_atm` | `nfile_atm_sstm`, `timename_atm_sstm`, `varname_atm_sstm`, `Lcycle_atm_sstm`, `Tcycle_atm_sstm`, `fnames_atm_sstm(:)` | SST climatology, for example `sst`. | Assigned to `agrd%sstm_atm`. |
| `&ssta_param_atm`, `&ssta_io_atm` | `nfile_atm_ssta`, `timename_atm_ssta`, `varname_atm_ssta`, `Lcycle_atm_ssta`, `Tcycle_atm_ssta`, `fnames_atm_ssta(:)` | SST anomaly, for example `sst`. | Assigned to `agrd%ssta_atm` and used to force the atmosphere response. |

### Optional Mean-Wind Input Groups

These groups are read only when `aset%heating_type=="ZC87_conv"`:

| Namelist groups | Variables | Field read by the code | Role in `zc_agcm.f90` |
| --- | --- | --- | --- |
| `&uam_param_atm`, `&uam_io_atm` | `nfile_atm_uam`, `timename_atm_uam`, `varname_atm_uam`, `Lcycle_atm_uam`, `Tcycle_atm_uam`, `fnames_atm_uam(:)` | Atmosphere mean zonal wind. | Passed to `return_uvp_fromSST_ZC87_conv`. |
| `&vam_param_atm`, `&vam_io_atm` | `nfile_atm_vam`, `timename_atm_vam`, `varname_atm_vam`, `Lcycle_atm_vam`, `Tcycle_atm_vam`, `fnames_atm_vam(:)` | Atmosphere mean meridional wind. | Passed to `return_uvp_fromSST_ZC87_conv`. |

The default `do_agcm_hindcast_clm_256_ZC87.nml` uses `aset%heating_type="ZC87"`,
so the optional mean-wind groups are not required.
