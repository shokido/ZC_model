# RUN/OGCM

This directory contains namelist files for standalone ocean model experiments.
The standard ocean spin-up example used by the top-level README is:

```text
do_ogcm_spn_eqpac_30_cd013H120.nml
```

Run from the repository root:

```bash
mkdir -p OUTPUTS/OGCM
cd RUN/OGCM
../../CODES/exec_solver_ocn_dyn.out < do_ogcm_spn_eqpac_30_cd013H120.nml
```

## Standard Namelist Summary

For `do_ogcm_spn_eqpac_30_cd013H120.nml`:

- Time step: `dt=3600.0` seconds.
- Start date: `19700101`.
- End date: `19900101`.
- Ocean grid: `../../INPUT/OCN/grid_eqpac_30.nc`.
- Initial condition: none, because `flag_ini_ocn="F"`.
- Wind forcing is cyclic with `Tcycle_ocn_taux=365` and
  `Tcycle_ocn_tauy=365`.

The exact scientific meaning and units of all physical parameters are to be
confirmed before publication use.

## Required Input Files

Required by `do_ogcm_spn_eqpac_30_cd013H120.nml`:

- `../../INPUT/OCN/grid_eqpac_30.nc`
- `../../Forcing/WIND/jra55do_v1.3_u10_clm_ocn_eqpac_30.nc`
- `../../Forcing/WIND/jra55do_v1.3_v10_clm_ocn_eqpac_30.nc`

The wind variables are read using:

- `varname_ocn_taux="u10"`
- `varname_ocn_tauy="v10"`

## Output Files

The standard namelist writes:

- `../../OUTPUTS/OGCM/avg_spinup_clm_eqpac_30_H120_cd1.3_1_20.nc`
- `../../OUTPUTS/OGCM/rst_spinup_clm_eqpac_30_H120_cd1.3_1_20.nc`

The average file contains ocean dynamical fields such as `u_ocn_sw`,
`v_ocn_sw`, `h_ocn_sw`, `u_ocn_ek`, `v_ocn_ek`, `w_ocn_ek`, `u_ocn_1`,
`v_ocn_1`, `w_ocn_1`, `taux`, and `tauy`, plus grid and mask variables.

The restart file is written at the end of the run and can be used as an initial
condition in later ocean experiments if the corresponding namelist is configured
with `flag_ini_ocn="T"` and `fname_ini_ocn` pointing to the restart file.

## Namelist Variables for `zc_ogcm_dyn.f90`

The executable `exec_solver_ocn_dyn.out` is built from `CODES/zc_ogcm_dyn.f90`.
This driver runs the ocean dynamical part of the ZC model without the SST
equation. It reads wind forcing, computes wind stress, Ekman currents,
baroclinic/geostrophic current, thermocline-depth anomaly, total current, and
writes averaged current/thermocline fields plus a restart file.

The namelist is read from standard input in this order:

1. `&date`
2. `&io_ocn`
3. `&param_ocn`
4. `&taux_param_ocn`
5. `&taux_io_ocn`
6. `&tauy_param_ocn`
7. `&tauy_io_ocn`

Because the program reads namelist groups sequentially, keep this order when
creating new namelist files.

### `&date`

| Variable | Type | Meaning in `zc_ogcm_dyn.f90` | Notes |
| --- | --- | --- | --- |
| `dt` | real | Model time step in seconds. | The number of time steps is computed from the run length divided by `dt`. |
| `start_yymmdd` | integer | Simulation start date in `YYYYMMDD` format. | Example: `19700101`. |
| `start_hhmmss` | integer | Simulation start time in `HHMMSS` format. | Example: `0` for 00:00:00. |
| `end_yymmdd` | integer | Simulation end date in `YYYYMMDD` format. | The run length is calculated from start to end using calendar utilities. |
| `end_hhmmss` | integer | Simulation end time in `HHMMSS` format. | Example: `0` for 00:00:00. |

### `&io_ocn`

| Variable | Type | Meaning in `zc_ogcm_dyn.f90` | Notes |
| --- | --- | --- | --- |
| `fname_grd_ocn` | character | Path to the ocean grid NetCDF file. | Read by `read_ocn_dyn_grd`. The grid must match the forcing dimensions. |
| `flag_ini_ocn` | character | Switch for restart initialization. | If exactly `"T"`, the code reads `fname_ini_ocn`; otherwise it starts from initialized zero/analyzed fields. |
| `fname_ini_ocn` | character | Path to an ocean dynamical restart file. | Used only when `flag_ini_ocn="T"`. |
| `fname_avg_ocn` | character | Path to the NetCDF average-output file. | Created before the time loop by `create_avg_ocn_dyn_ZC`. Existing files may be overwritten; behavior to be confirmed. |
| `out_avg_flag` | integer | Time unit used for average-output scheduling. | Passed to calendar routines. See the output-time-unit table below. |
| `out_avg_int` | integer | Output interval in units of `out_avg_flag`. | Example: `out_avg_flag=100`, `out_avg_int=1` gives one output per calendar month. |
| `fname_rst_ocn` | character | Path to the ocean dynamical restart file written at the end of the run. | Written after the final time step by `write_restart_ocn_dyn`. |

Output-time-unit flags used by the calendar routines:

| Flag value | Unit |
| --- | --- |
| `-10000` | seconds |
| `-100` | minutes |
| `-1` | hours |
| `1` | days |
| `100` | months |
| `10000` | years |

The output `time` coordinate is written at the center of each averaging window,
while the write step is scheduled at the end of each window.

### `&param_ocn`

`param_ocn` contains one variable, `oset`, which is a Fortran derived type
defined as `type(ocn_set)` in `CODES/run_types.f90`. The table below lists the
members of `oset`.

| Variable | Default in `run_types.f90` | Meaning in `zc_ogcm_dyn.f90` | Notes |
| --- | --- | --- | --- |
| `oset%wbc_flag_p` | `"Clo"` | Western boundary condition for p-grid variables, including `h_sw`. | Accepted values in boundary routines are `Clo`, `Gra`, and `Per` (case variants also accepted). |
| `oset%ebc_flag_p` | `"Clo"` | Eastern boundary condition for p-grid variables. | `Clo` sets boundary values to zero, `Gra` uses a zero-gradient copy, `Per` applies periodic wrapping. |
| `oset%nbc_flag_p` | `"Clo"` | Northern boundary condition for p-grid variables. | Same accepted values as above. |
| `oset%sbc_flag_p` | `"Clo"` | Southern boundary condition for p-grid variables. | Same accepted values as above. |
| `oset%wbc_flag_u` | `"Clo"` | Western boundary condition for u-grid velocity. | Used by `set_bc_u`. |
| `oset%ebc_flag_u` | `"Clo"` | Eastern boundary condition for u-grid velocity. | Used by `set_bc_u`. |
| `oset%nbc_flag_u` | `"Clo"` | Northern boundary condition for u-grid velocity. | Used by `set_bc_u`; closed meridional ghost cells depend on `slip_ind`. |
| `oset%sbc_flag_u` | `"Clo"` | Southern boundary condition for u-grid velocity. | Used by `set_bc_u`; closed meridional ghost cells depend on `slip_ind`. |
| `oset%wbc_flag_v` | `"Clo"` | Western boundary condition for v-grid velocity. | Used by `set_bc_v`; closed zonal ghost cells depend on `slip_ind`. |
| `oset%ebc_flag_v` | `"Clo"` | Eastern boundary condition for v-grid velocity. | Used by `set_bc_v`; closed zonal ghost cells depend on `slip_ind`. |
| `oset%nbc_flag_v` | `"Clo"` | Northern boundary condition for v-grid velocity. | Used by `set_bc_v`. |
| `oset%sbc_flag_v` | `"Clo"` | Southern boundary condition for v-grid velocity. | Used by `set_bc_v`. |
| `oset%slip_ind` | `0.0` | Slip/no-slip index used when setting masks and closed velocity boundaries. | Code comments state `0` gives slip-like condition (`du/dx=0`) and `1` gives no-slip-like condition (`u=0`). Exact physical interpretation to be confirmed. |
| `oset%cp_ocn` | `2.9` | Baroclinic gravity-wave speed. | Comment in source gives units of m/s. Used to compute reduced gravity as `cp_ocn^2 / (H1 + H2)`. |
| `oset%r_ocn_day` | `912.5` | Inverse damping coefficient for baroclinic current, in days. | Used in `initialize_ocn_visc` to scale grid damping arrays by `1/(r_ocn_day * day_to_sec)`. |
| `oset%H1` | `50.0` | Mixed-layer depth. | Comment in source gives units of m. Used in Ekman forcing and total-current calculations. |
| `oset%H2` | `80.0` | Upper-layer/deep-layer depth parameter. | Comment in source gives units of m. Used with `H1` in dynamical equations and total-current calculation. Exact layer interpretation to be confirmed. |
| `oset%eps_s_ocn_day` | `2.0` | Inverse damping timescale for Ekman current, in days. | Used in `solve_ekman_ocn` as `1/(eps_s_ocn_day * day_to_sec)`. |
| `oset%eps_s_sst_day` | `125.0` | SST damping timescale parameter. | Present in `ocn_set`, but not used by `zc_ogcm_dyn.f90` because this driver does not solve the SST equation. |
| `oset%cd_bulk` | `0.0018` | Bulk coefficient for wind-stress calculation. | With the current Makefile, `-DWSTRESS_BULK` is enabled, so `u10`/`v10` winds are converted to stress using `rho_air * cd_bulk * speed * wind`. |
| `oset%nu` | `100000.0` | Horizontal viscosity. | Comment in source gives units of m^2/s. Used to fill `visc_2D`. |
| `oset%Tsub_T1` | `28.0` | Subsurface-temperature parameter. | Present in `ocn_set`, but not used by `zc_ogcm_dyn.f90`; relevant to SST/subsurface-temperature code paths. |
| `oset%Tsub_T2` | `-40.0` | Subsurface-temperature parameter. | Not used by `zc_ogcm_dyn.f90`. Exact meaning to be confirmed. |
| `oset%Tsub_b1` | `1/80` | Subsurface-temperature parameter. | Not used by `zc_ogcm_dyn.f90`. Exact meaning to be confirmed. |
| `oset%Tsub_b2` | `1/33` | Subsurface-temperature parameter. | Not used by `zc_ogcm_dyn.f90`. Exact meaning to be confirmed. |
| `oset%Tsub_gamma` | `0.75` | Subsurface-temperature parameter. | Not used by `zc_ogcm_dyn.f90`. Exact meaning to be confirmed. |
| `oset%kick_WWB` | `"F"` | Westerly-wind-burst switch. | Present in `ocn_set`, but not used by `zc_ogcm_dyn.f90`. |
| `oset%G0_wwb` | `0.09` | WWB occurrence/intensity parameter. | Not used by `zc_ogcm_dyn.f90`; exact meaning to be confirmed. |
| `oset%G1_wwb` | `0.0` | WWB occurrence/intensity parameter. | Not used by `zc_ogcm_dyn.f90`; exact meaning to be confirmed. |
| `oset%dur_wwb` | `20.0` | WWB duration parameter. | Not used by `zc_ogcm_dyn.f90`; source comment indicates days in WWB routine. |
| `oset%us0_wwb` | `6.5` | WWB wind-amplitude parameter. | Not used by `zc_ogcm_dyn.f90`; exact meaning to be confirmed. |
| `oset%widx_wwb` | `20.0` | WWB zonal width parameter. | Not used by `zc_ogcm_dyn.f90`; exact units to be confirmed. |
| `oset%widy_wwb` | `6.0` | WWB meridional width parameter. | Not used by `zc_ogcm_dyn.f90`; exact units to be confirmed. |
| `oset%x0_wwb` | `160.0` | WWB center longitude parameter. | Not used by `zc_ogcm_dyn.f90`. |
| `oset%y0_wwb` | `0.0` | WWB center latitude parameter. | Not used by `zc_ogcm_dyn.f90`. |

### `&taux_param_ocn`

| Variable | Type | Meaning in `zc_ogcm_dyn.f90` | Notes |
| --- | --- | --- | --- |
| `nfile_ocn_taux` | integer | Number of zonal wind/wind-stress NetCDF files to read. | The program allocates `fnames_ocn_taux(nfile_ocn_taux)` after reading this value. |
| `timename_ocn_taux` | character | Name of the time coordinate variable in the zonal forcing file(s). | Example: `"time"`. The time variable must have a CF-like `units` attribute containing `since`. |
| `varname_ocn_taux` | character | Name of the zonal forcing variable to read. | With `-DWSTRESS_BULK`, this is interpreted as zonal wind speed; otherwise as zonal wind stress. |
| `Lcycle_ocn_taux` | character(1) | Cyclic-forcing switch for zonal forcing. | If `"T"`, the code uses day-of-year cyclic interpolation; otherwise it interpolates/clamps in absolute model time. |
| `Tcycle_ocn_taux` | real | Cycle length for cyclic zonal forcing. | Usually `365` days for climatological annual-cycle forcing. |

### `&taux_io_ocn`

| Variable | Type | Meaning in `zc_ogcm_dyn.f90` | Notes |
| --- | --- | --- | --- |
| `fnames_ocn_taux(:)` | character array | Paths to zonal wind/wind-stress NetCDF files. | Number of entries must match `nfile_ocn_taux`; each file is concatenated in time. |

### `&tauy_param_ocn`

| Variable | Type | Meaning in `zc_ogcm_dyn.f90` | Notes |
| --- | --- | --- | --- |
| `nfile_ocn_tauy` | integer | Number of meridional wind/wind-stress NetCDF files to read. | The program allocates `fnames_ocn_tauy(nfile_ocn_tauy)` after reading this value. |
| `timename_ocn_tauy` | character | Name of the time coordinate variable in the meridional forcing file(s). | Example: `"time"`. |
| `varname_ocn_tauy` | character | Name of the meridional forcing variable to read. | With `-DWSTRESS_BULK`, this is interpreted as meridional wind speed; otherwise as meridional wind stress. |
| `Lcycle_ocn_tauy` | character(1) | Cyclic-forcing switch for meridional forcing. | If `"T"`, the code uses day-of-year cyclic interpolation; otherwise it interpolates/clamps in absolute model time. |
| `Tcycle_ocn_tauy` | real | Cycle length for cyclic meridional forcing. | Usually `365` days for climatological annual-cycle forcing. |

### `&tauy_io_ocn`

| Variable | Type | Meaning in `zc_ogcm_dyn.f90` | Notes |
| --- | --- | --- | --- |
| `fnames_ocn_tauy(:)` | character array | Paths to meridional wind/wind-stress NetCDF files. | Number of entries must match `nfile_ocn_tauy`; each file is concatenated in time. |

## Namelist Variables for `zc_ogcm_full.f90`

The executable `exec_solver_ocn_full.out` is built from `CODES/zc_ogcm_full.f90`.
This driver runs the OGCM anomaly model with both ocean dynamics and the SST
anomaly equation. In addition to wind anomalies, it reads mean SST, mean winds,
and mean ocean fields used by the SST equation.

A corresponding hindcast example is:

```text
do_ogcm_hindcast_eqpac_30_cd013H120.nml
```

Run from the repository root:

```bash
mkdir -p OUTPUTS/OGCM OUTPUTS/CGCM
cd RUN/OGCM
../../CODES/exec_solver_ocn_full.out < do_ogcm_hindcast_eqpac_30_cd013H120.nml
```

Note: this namelist writes `fname_diag_ocn` under `../../OUTPUTS/CGCM/`, even
though the run is launched from `RUN/OGCM`. Confirm whether this is intentional
before sharing a standard workflow.

### Standard Full-Model Inputs and Outputs

For `do_ogcm_hindcast_eqpac_30_cd013H120.nml`, the main inputs are:

- `../../INPUT/OCN/grid_eqpac_30.nc`
- `../../Forcing/WIND/jra55do_v1.3_u10_anm_ocn_eqpac_30.nc`
- `../../Forcing/WIND/jra55do_v1.3_v10_anm_ocn_eqpac_30.nc`
- `../../Forcing/OCN/basic_clm_eqpac_30_H120_cd1.3_1_20.nc`
- `../../Forcing/WIND/jra55do_v1.3_u10_clm_ocn_eqpac_30.nc`
- `../../Forcing/WIND/jra55do_v1.3_v10_clm_ocn_eqpac_30.nc`
- `../../Forcing/Tzm/basic_2D_eqpac_30_ZC87.nc`

The standard namelist writes:

- Average output:
  `../../OUTPUTS/OGCM/avg_ogcm_hindcast_eqpac_30_H120_cd1.3_1_20.nc`
- Diagnostic output:
  `../../OUTPUTS/CGCM/diag_ogcm_hindcast_eqpac_30_H120_cd1.3_1_20.nc`
- Restart output:
  `../../OUTPUTS/OGCM/rst_ogcm_hindcast_eqpac_30_H120_cd1.3_1_20.nc`

The average file includes the dynamical variables from `zc_ogcm_dyn.f90` and
also `ssta`, the simulated SST anomaly. The diagnostic file contains SST-budget
terms including `dSSTdt`, `uaTm`, `umTa`, `uaTa`, `vaTm`, `vmTa`, `vaTa`,
`waTm`, `wmTa`, `waTa`, `qh`, `Tsub`, and `Te`.

### Read Order

The namelist is read from standard input in this order:

1. `&date`
2. `&io_ocn`
3. `&param_ocn`
4. `&taux_param_ocn`
5. `&taux_io_ocn`
6. `&tauy_param_ocn`
7. `&tauy_io_ocn`
8. `&sstm_param_ocn`
9. `&sstm_io_ocn`
10. `&tauxm_param_ocn`
11. `&tauxm_io_ocn`
12. `&tauym_param_ocn`
13. `&tauym_io_ocn`
14. `&hm_param_ocn`
15. `&hm_io_ocn`
16. `&um_param_ocn`
17. `&um_io_ocn`
18. `&vm_param_ocn`
19. `&vm_io_ocn`
20. `&wm_param_ocn`
21. `&wm_io_ocn`
22. `&Tzm_param_ocn`
23. `&Tzm_io_ocn`

Keep this order when creating new namelists. Some current hindcast namelists
also include `&sstm_param_atm` and `&sstm_io_atm` after `&Tzm_io_ocn`; these
groups are not read by `zc_ogcm_full.f90` and appear to be unused here, to be
confirmed.

### Shared Groups

The following groups have the same basic meaning as in `zc_ogcm_dyn.f90`:

- `&date`: `dt`, `start_yymmdd`, `start_hhmmss`, `end_yymmdd`, `end_hhmmss`.
- Most of `&io_ocn`: `fname_grd_ocn`, `flag_ini_ocn`, `fname_ini_ocn`,
  `fname_avg_ocn`, `out_avg_flag`, `out_avg_int`, and `fname_rst_ocn`.
- `&param_ocn`: `oset`.
- `&taux_param_ocn`, `&taux_io_ocn`, `&tauy_param_ocn`, `&tauy_io_ocn`.

Important differences from `zc_ogcm_dyn.f90`:

| Variable or group | Meaning in `zc_ogcm_full.f90` | Notes |
| --- | --- | --- |
| `fname_diag_ocn` | Path to the SST diagnostic NetCDF file. | Created by `create_diag_ocn_ZC`. |
| `out_diag_flag` | Time unit used for diagnostic-output scheduling. | Uses the same flag values as `out_avg_flag`: seconds, minutes, hours, days, months, or years. |
| `out_diag_int` | Diagnostic output interval in units of `out_diag_flag`. | Example: `out_diag_flag=100`, `out_diag_int=1` gives one diagnostic output per calendar month. |
| `flag_ini_ocn` / `fname_ini_ocn` | If `flag_ini_ocn="T"`, both dynamical and SST restart fields are read. | Calls both `read_restart_ocn_dyn` and `read_restart_ocn_sst`. |
| `oset%eps_s_sst_day` | SST damping timescale parameter. | Used in the SST equation as `1/(eps_s_sst_day * day_to_sec)`. |
| `oset%Tsub_T1`, `oset%Tsub_T2`, `oset%Tsub_b1`, `oset%Tsub_b2`, `oset%Tsub_gamma` | Subsurface/entrainment-temperature parameters. | Used through `get_tsub_ZC` inside the SST equation; exact scientific interpretation to be confirmed. |

With the current Makefile, `-DWSTRESS_BULK` is enabled. In this mode,
`zc_ogcm_full.f90` reads wind anomalies plus mean winds and computes wind-stress
anomalies using `ua_to_stress_anm`. If `-DWSTRESS_BULK` is disabled, the input
variables are treated as wind stress and the code subtracts the mean stress
fields directly.

### Additional Input Groups

All additional full-model input groups use the same pattern:

- `nfile_*`: number of NetCDF files to read.
- `timename_*`: name of the time coordinate variable.
- `varname_*`: name of the physical variable to read.
- `Lcycle_*`: cyclic-input switch. `"T"` uses day-of-year cyclic interpolation;
  other values use absolute model-time interpolation/clamping.
- `Tcycle_*`: cycle length, usually `365` days for annual-cycle climatology.
- `fnames_*(1:nfile_*)`: paths to the NetCDF files.

| Namelist groups | Variables | Field read by the code | Grid reader | Role in `zc_ogcm_full.f90` |
| --- | --- | --- | --- | --- |
| `&sstm_param_ocn`, `&sstm_io_ocn` | `nfile_ocn_sstm`, `timename_ocn_sstm`, `varname_ocn_sstm`, `Lcycle_ocn_sstm`, `Tcycle_ocn_sstm`, `fnames_ocn_sstm(:)` | Mean/background SST, for example `sst_ocn`. | `read_data_TLL_p` | Used as the mean SST field for horizontal SST gradients and total SST context. |
| `&tauxm_param_ocn`, `&tauxm_io_ocn` | `nfile_ocn_tauxm`, `timename_ocn_tauxm`, `varname_ocn_tauxm`, `Lcycle_ocn_tauxm`, `Tcycle_ocn_tauxm`, `fnames_ocn_tauxm(:)` | Mean/background zonal wind, for example `u10`. | `read_data_TLL_p` | Used with anomalous zonal wind to compute wind-stress anomalies when `-DWSTRESS_BULK` is enabled. |
| `&tauym_param_ocn`, `&tauym_io_ocn` | `nfile_ocn_tauym`, `timename_ocn_tauym`, `varname_ocn_tauym`, `Lcycle_ocn_tauym`, `Tcycle_ocn_tauym`, `fnames_ocn_tauym(:)` | Mean/background meridional wind, for example `v10`. | `read_data_TLL_p` | Used with anomalous meridional wind to compute wind-stress anomalies when `-DWSTRESS_BULK` is enabled. |
| `&hm_param_ocn`, `&hm_io_ocn` | `nfile_ocn_hm`, `timename_ocn_hm`, `varname_ocn_hm`, `Lcycle_ocn_hm`, `Tcycle_ocn_hm`, `fnames_ocn_hm(:)` | Mean thermocline-depth field, for example `h_ocn`. | `read_data_TLL_p` | Passed to `solve_sst_ocn_ZC` as `mhdata`; used in subsurface-temperature/upwelling feedback calculations. |
| `&um_param_ocn`, `&um_io_ocn` | `nfile_ocn_um`, `timename_ocn_um`, `varname_ocn_um`, `Lcycle_ocn_um`, `Tcycle_ocn_um`, `fnames_ocn_um(:)` | Mean zonal current, for example `u_ocn`. | `read_data_TLL_u` | Used as climatological zonal current in SST advection terms. |
| `&vm_param_ocn`, `&vm_io_ocn` | `nfile_ocn_vm`, `timename_ocn_vm`, `varname_ocn_vm`, `Lcycle_ocn_vm`, `Tcycle_ocn_vm`, `fnames_ocn_vm(:)` | Mean meridional current, for example `v_ocn`. | `read_data_TLL_v` | Used as climatological meridional current in SST advection terms. |
| `&wm_param_ocn`, `&wm_io_ocn` | `nfile_ocn_wm`, `timename_ocn_wm`, `varname_ocn_wm`, `Lcycle_ocn_wm`, `Tcycle_ocn_wm`, `fnames_ocn_wm(:)` | Mean vertical velocity, for example `w_ocn`. | `read_data_TLL_p` | Used as climatological vertical velocity in SST/upwelling terms. |
| `&Tzm_param_ocn`, `&Tzm_io_ocn` | `nfile_ocn_Tzm`, `timename_ocn_Tzm`, `varname_ocn_Tzm`, `Lcycle_ocn_Tzm`, `Tcycle_ocn_Tzm`, `fnames_ocn_Tzm(:)` | Mean vertical temperature-gradient field, for example `t_an`. | `read_data_TLL_p` | Used as `dTbdz` in vertical-advection/upwelling feedback terms. |

### SST Diagnostic Output Variables

The SST diagnostic file is controlled by `fname_diag_ocn`, `out_diag_flag`, and
`out_diag_int`. The following variables are defined by `create_diag_ocn_ZC`:

| Variable | Meaning from source attributes/code | Units in output |
| --- | --- | --- |
| `dSSTdt` | SST tendency. | `degrees_celsius/s` |
| `uaTm` | Advection of climatological zonal SST gradient by anomalous zonal current. | `degrees_celsius/s` |
| `umTa` | Advection of anomalous zonal SST gradient by climatological zonal current. | `degrees_celsius/s` |
| `uaTa` | Advection of anomalous zonal SST gradient by anomalous zonal current. | `degrees_celsius/s` |
| `vaTm` | Advection of climatological meridional SST gradient by anomalous meridional current. | `degrees_celsius/s` |
| `vmTa` | Advection of anomalous meridional SST gradient by climatological meridional current. | `degrees_celsius/s` |
| `vaTa` | Advection of anomalous meridional SST gradient by anomalous meridional current. | `degrees_celsius/s` |
| `waTm` | Advection of climatological vertical SST gradient by anomalous vertical current. | `degrees_celsius/s` |
| `wmTa` | Advection of anomalous vertical SST gradient by climatological vertical current. | `degrees_celsius/s` |
| `waTa` | Advection of anomalous vertical SST gradient by anomalous vertical current. | `degrees_celsius/s` |
| `qh` | Thermal damping. | `degrees_celsius/s` |
| `Tsub` | Subsurface temperature. | `degrees_celsius` |
| `Te` | Entrainment temperature. | `degrees_celsius` |

The exact interpretation of each SST-budget term should be confirmed against
the model equations before scientific analysis.
