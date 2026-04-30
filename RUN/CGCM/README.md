# RUN/CGCM

This directory contains namelist files for coupled ocean-atmosphere model
experiments. The default coupled example documented here is:

```text
do_cgcm_eqpac_30_clm_c1.4H120_dt3600_c10day.nml
```

This namelist is used with `CODES/exec_solver_couple_full.out`, which is built
from `CODES/zc_cgcm_full.f90`.

## How to Run

Before running the CGCM, prepare the ocean basic-state fields from a standalone
OGCM spin-up. The standard workflow is:

1. Run the OGCM spin-up documented in `RUN/OGCM/README.md`.
2. Create the ocean mean/basic fields in `Forcing/OCN`.
3. Run the CGCM namelist from this directory.

After completing the OGCM spin-up, create either annual-mean or monthly
climatological ocean fields:

```bash
cd Forcing/OCN
python make_meanfield_ann.py   # annual-mean climatology / no annual cycle
python make_meanfield_clm.py   # monthly climatology / seasonal cycle
```

Confirm which of these outputs matches the selected CGCM namelist. The default
namelist documented here uses climatological fields with an annual cycle
(`Lcycle="T"`, `Tcycle=365`).

Then run the CGCM from the repository root:

```bash
mkdir -p OUTPUTS/CGCM
cd RUN/CGCM
../../CODES/exec_solver_couple_full.out < do_cgcm_eqpac_30_clm_c1.4H120_dt3600_c10day.nml
```

The run writes progress information to standard output and writes NetCDF files
to the paths specified in the namelist.

Runtime may be long because this default experiment spans 1970-01-01 to
2020-01-01 with `dt=3600.0`. For a quick smoke test, copy the namelist, shorten
`end_yymmdd`, and change all output filenames in the copy.

## Standard Namelist Summary

For `do_cgcm_eqpac_30_clm_c1.4H120_dt3600_c10day.nml`:

- Time step: `dt=3600.0` seconds.
- Start date: `19700101`.
- End date: `20200101`.
- Coupling interval: `time_couple=864000` seconds, or 10 days.
- Ocean grid: `../../INPUT/OCN/grid_eqpac_30.nc`.
- Atmosphere grid: `../../INPUT/ATM/grid_256.nc`.
- Coupler file: `../../INPUT/COUPLER/coupler_eqpac_30_atm_256.nc`.
- Initial ocean restart: none, because `flag_ini_ocn="F"`.
- Initial wind perturbation: enabled, because `kick_ini="T"`.
- Ocean bulk coefficient: `oset%cd_bulk=0.0014`.
- Background/forcing fields are climatological and cyclic with `Tcycle=365`.

Although this CGCM example is often run after standalone OGCM experiments, this
namelist does not directly read an OGCM restart file. If restart initialization
is desired, copy the namelist and set `flag_ini_ocn="T"` and `fname_ini_ocn` to
the desired restart file. 

The exact scientific meaning and units of all physical parameters are to be
confirmed before publication use.

## Required Input Files

Required by `do_cgcm_eqpac_30_clm_c1.4H120_dt3600_c10day.nml`:

- `../../INPUT/OCN/grid_eqpac_30.nc`
- `../../INPUT/ATM/grid_256.nc`
- `../../INPUT/COUPLER/coupler_eqpac_30_atm_256.nc`
- `../../Forcing/OCN/basic_clm_eqpac_30_H120_cd1.4_1_20.nc`
- `../../Forcing/WIND/jra55do_v1.3_u10_clm_ocn_eqpac_30.nc`
- `../../Forcing/WIND/jra55do_v1.3_v10_clm_ocn_eqpac_30.nc`
- `../../Forcing/Tzm/basic_2D_eqpac_30_ZC87.nc`
- `../../Forcing/SST/ERSST_v5_atm_256_clm.nc`

Variables read from these files include:

- Ocean background SST: `sst_ocn`.
- Ocean background winds: `u10`, `v10`.
- Ocean background fields: `h_ocn`, `u_ocn`, `v_ocn`, `w_ocn`.
- Background vertical temperature-gradient field: `t_an`.
- Atmosphere background SST: `sst`.

All listed forcing/background fields are cyclic in the default namelist with
`Tcycle=365`.

## Output Files

The default namelist writes:

- Ocean average:
  `../../OUTPUTS/CGCM/avg_cgcm_full_eqpac_30_ocn_clm_H120_cd1.4_dt3600_c10day.nc`
- Ocean diagnostics:
  `../../OUTPUTS/CGCM/diag_cgcm_full_eqpac_30_ocn_clm_H120_cd1.4_dt3600_c10day.nc`
- Atmosphere average:
  `../../OUTPUTS/CGCM/avg_cgcm_full_eqpac_30_atm_clm_H120_cd1.4_dt3600_c10day.nc`
- Ocean restart:
  `../../OUTPUTS/CGCM/rst_cgcm_full_eqpac_30_ocn_clm_H120_cd1.4_dt3600_c10day.nc`

The ocean average file contains ocean dynamical fields such as `u_ocn_sw`,
`v_ocn_sw`, `h_ocn_sw`, `u_ocn_ek`, `v_ocn_ek`, `w_ocn_ek`, `u_ocn_1`,
`v_ocn_1`, `w_ocn_1`, `taux`, `tauy`, and `ssta`, plus grid and mask
variables.

The ocean diagnostic file contains SST-budget terms such as `dSSTdt`, `uaTm`,
`umTa`, `uaTa`, `vaTm`, `vmTa`, `vaTa`, `waTm`, `wmTa`, `waTa`, `qh`, `Tsub`,
and `Te`.

The atmosphere average file contains fields such as `q_atm`, `u_atm`, `v_atm`,
`p_atm`, `ssta_atm`, and `sstm_atm`.

The meaning and units of all output variables should be confirmed from the
Fortran source and experiment documentation before scientific interpretation.

## Namelist Variables for `zc_cgcm_full.f90`

`zc_cgcm_full.f90` runs the fully coupled ocean-atmosphere model. It reads
background ocean fields, background atmosphere SST, a coupler mapping file,
computes the atmosphere response to ocean SST anomalies, maps atmosphere winds
back to the ocean, solves ocean dynamics and SST anomaly, and writes ocean and
atmosphere averages plus ocean diagnostics and restart fields.

### Read Order

The namelist is read from standard input in this order:

1. `&date`
2. `&io_ocn`
3. `&io_atm`
4. `&io_coupler`
5. `&param_ocn`
6. `&param_atm`
7. `&param_coupler`
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
24. `&sstm_param_atm`
25. `&sstm_io_atm`

If `aset%heating_type=="ZC87_conv"`, the code then also reads:

26. `&uam_param_atm`
27. `&uam_io_atm`
28. `&vam_param_atm`
29. `&vam_io_atm`

Keep this order when creating new namelists.

### `&date`

| Variable | Type | Meaning in `zc_cgcm_full.f90` | Notes |
| --- | --- | --- | --- |
| `dt` | real | Model time step in seconds. | Used for the ocean time stepping and to compute total step count. |
| `start_yymmdd` | integer | Simulation start date in `YYYYMMDD` format. | Example: `19700101`. |
| `start_hhmmss` | integer | Simulation start time in `HHMMSS` format. | Example: `0` for 00:00:00. |
| `end_yymmdd` | integer | Simulation end date in `YYYYMMDD` format. | The run length is calculated from start to end using calendar utilities. |
| `end_hhmmss` | integer | Simulation end time in `HHMMSS` format. | Example: `0` for 00:00:00. |
| `time_couple` | real | Coupling interval in seconds. | The code sets `ntime_couple=int(time_couple/dt)` and updates the atmosphere/coupler when `mod(itime, ntime_couple)==1`. |

### `&io_ocn`

| Variable | Type | Meaning in `zc_cgcm_full.f90` | Notes |
| --- | --- | --- | --- |
| `fname_grd_ocn` | character | Path to the ocean grid NetCDF file. | Read by `read_ocn_dyn_grd` and `read_ocn_sst_grd`; also used for WWB mask if enabled. |
| `flag_ini_ocn` | character | Switch for restart initialization. | If exactly `"T"`, both ocean dynamics and SST restart fields are read from `fname_ini_ocn`. |
| `fname_ini_ocn` | character | Path to an ocean restart file. | Used only when `flag_ini_ocn="T"`. |
| `fname_avg_ocn` | character | Path to the ocean average-output NetCDF file. | Contains ocean dynamics and `ssta`. |
| `out_avg_flag` | integer | Time unit used for average-output scheduling. | Shared by ocean and atmosphere average-output creation. |
| `out_avg_int` | integer | Average-output interval in units of `out_avg_flag`. | Example: `out_avg_flag=100`, `out_avg_int=1` gives one output per calendar month. |
| `fname_diag_ocn` | character | Path to the ocean SST diagnostic NetCDF file. | Created by `create_diag_ocn_ZC`. |
| `out_diag_flag` | integer | Time unit used for diagnostic-output scheduling. | Uses the same flag values as `out_avg_flag`. |
| `out_diag_int` | integer | Diagnostic-output interval in units of `out_diag_flag`. | Example: monthly diagnostics when `out_diag_flag=100`, `out_diag_int=1`. |
| `fname_rst_ocn` | character | Path to the ocean restart file written at the end of the run. | Written after the final step by `write_restart_ocn_dyn` and `write_restart_ocn_sst`. |

Output-time-unit flags used by the calendar routines:

| Flag value | Unit |
| --- | --- |
| `-10000` | seconds |
| `-100` | minutes |
| `-1` | hours |
| `1` | days |
| `100` | months |
| `10000` | years |

### `&io_atm`

| Variable | Type | Meaning in `zc_cgcm_full.f90` | Notes |
| --- | --- | --- | --- |
| `fname_grd_atm` | character | Path to the atmosphere grid NetCDF file. | Read by `read_atm_grd`; coordinates are then prepared by `set_coord_atm`. |
| `fname_avg_atm` | character | Path to the atmosphere average-output NetCDF file. | Contains atmosphere response fields and SST fields on the atmosphere grid. |
| `out_avg_flag` | integer | Average-output time unit. | Same variable name as in `&io_ocn`; because it is global in `mod_io_avg`, keep the ocean and atmosphere values consistent. |
| `out_avg_int` | integer | Average-output interval. | Same note as above. |

### `&io_coupler`

| Variable | Type | Meaning in `zc_cgcm_full.f90` | Notes |
| --- | --- | --- | --- |
| `fname_coupler` | character | Path to the ocean-atmosphere coupler mapping file. | Read by `read_coupler`; used by `exchange_OtoA` and `exchange_AtoO`. |

### `&param_ocn`

`param_ocn` contains one variable, `oset`, a derived type defined as
`type(ocn_set)` in `CODES/run_types.f90`. The coupled model uses the same ocean
parameters documented in `RUN/OGCM/README.md`, with the following especially
important members:

| Variable | Default in `run_types.f90` | Meaning in `zc_cgcm_full.f90` | Notes |
| --- | --- | --- | --- |
| `oset%cd_bulk` | `0.0018` | Bulk coefficient for wind-stress calculation. | With the current Makefile, `-DWSTRESS_BULK` is enabled and atmosphere winds are converted to stress using the bulk formula. The default namelist sets this to `0.0014`. |
| `oset%kick_WWB` | `"F"` | Westerly-wind-burst switch. | If `"T"`, `zc_cgcm_full.f90` calls `ua_to_stress_anm_WWB`; otherwise it calls `ua_to_stress_anm`. |
| `oset%G0_wwb` | `0.09` | Additive WWB occurrence-rate parameter. | The daily occurrence rate is `lambda_perday=max(0, G0_wwb + G1_wwb * nino34)`. |
| `oset%G1_wwb` | `0.0` | State-dependent WWB occurrence-rate parameter. | Controls the part of the event probability that depends on the model Nino3.4-like SST anomaly. |
| `oset%dur_wwb` | `20.0` | WWB duration parameter. | Used in the Gaussian time envelope; units are days. The wind pulse is centered at `t0_wwb + dur_wwb`. |
| `oset%us0_wwb` | `6.5` | Peak zonal wind-amplitude parameter. | Added to `ua_ocn` inside the WWB spatial-temporal envelope before bulk wind-stress conversion. |
| `oset%widx_wwb` | `20.0` | WWB zonal width parameter. | Longitude width scale in the Gaussian spatial envelope. |
| `oset%widy_wwb` | `6.0` | WWB meridional width parameter. | Latitude width scale in the Gaussian spatial envelope. |
| `oset%x0_wwb` | `160.0` | WWB center longitude. | Center longitude of the Gaussian WWB patch. |
| `oset%y0_wwb` | `0.0` | WWB center latitude. | Center latitude of the Gaussian WWB patch. |

The WWB parameterization is implemented in `CODES/mod_wstress.f90`. At each
time step, the code computes an SST-index value over `mask_wwb`, then sets
`lambda_perday=max(0, G0_wwb + G1_wwb * nino34)`. The event threshold for one
time step is `1 - exp(-lambda_perday * dt_in_days)`. If a random number is
smaller than this threshold, a new WWB event starts at the current model time.
The added zonal wind has a Gaussian envelope in time, longitude, and latitude,
with amplitude `us0_wwb`, duration scale `dur_wwb`, center
`(x0_wwb, y0_wwb)`, and width scales `(widx_wwb, widy_wwb)`.

### `&param_atm`

`param_atm` contains one variable, `aset`, a derived type defined as
`type(atm_set)` in `CODES/run_types.f90`. See `RUN/AGCM/README.md` for the
detailed table of atmospheric parameters. In the CGCM, these parameters control
the Gill-type atmosphere response used after ocean SST anomalies are mapped to
the atmosphere grid.

The default CGCM namelist leaves `aset%heating_type` unset, so the type default
`"ZC87"` is used. If `aset%heating_type=="ZC87_conv"`, additional atmosphere
mean-wind input groups are required; see the atmosphere input section below.

### `&param_coupler`

| Variable | Type | Meaning in `zc_cgcm_full.f90` | Notes |
| --- | --- | --- | --- |
| `kick_ini` | character(1) | Switch for an initial zonal-wind perturbation. | If `"T"`, the code applies a temporary `ua_ocn` perturbation over 145E-190E for the first 120 model days, following Zebiak and Cane (1987) |

### Ocean Background Input Groups

All ocean background input groups use the same pattern:

- `nfile_*`: number of NetCDF files to read.
- `timename_*`: name of the time coordinate variable.
- `varname_*`: name of the physical variable to read.
- `Lcycle_*`: cyclic-input switch. `"T"` uses day-of-year cyclic interpolation;
  other values use absolute model-time interpolation/clamping.
- `Tcycle_*`: cycle length, usually `365` days for annual-cycle climatology.
- `fnames_*(1:nfile_*)`: paths to the NetCDF files.

| Namelist groups | Variables | Field read by the code | Grid reader | Role in `zc_cgcm_full.f90` |
| --- | --- | --- | --- | --- |
| `&sstm_param_ocn`, `&sstm_io_ocn` | `nfile_ocn_sstm`, `timename_ocn_sstm`, `varname_ocn_sstm`, `Lcycle_ocn_sstm`, `Tcycle_ocn_sstm`, `fnames_ocn_sstm(:)` | Ocean mean/background SST, for example `sst_ocn`. | `read_data_TLL_p` | Used as mean SST in the ocean SST equation. |
| `&tauxm_param_ocn`, `&tauxm_io_ocn` | `nfile_ocn_tauxm`, `timename_ocn_tauxm`, `varname_ocn_tauxm`, `Lcycle_ocn_tauxm`, `Tcycle_ocn_tauxm`, `fnames_ocn_tauxm(:)` | Mean/background zonal wind, for example `u10`. | `read_data_TLL_p` | Used with atmosphere wind anomalies to compute ocean wind-stress anomalies. |
| `&tauym_param_ocn`, `&tauym_io_ocn` | `nfile_ocn_tauym`, `timename_ocn_tauym`, `varname_ocn_tauym`, `Lcycle_ocn_tauym`, `Tcycle_ocn_tauym`, `fnames_ocn_tauym(:)` | Mean/background meridional wind, for example `v10`. | `read_data_TLL_p` | Used with atmosphere wind anomalies to compute ocean wind-stress anomalies. |
| `&hm_param_ocn`, `&hm_io_ocn` | `nfile_ocn_hm`, `timename_ocn_hm`, `varname_ocn_hm`, `Lcycle_ocn_hm`, `Tcycle_ocn_hm`, `fnames_ocn_hm(:)` | Mean thermocline-depth field, for example `h_ocn`. | `read_data_TLL_p` | Used in subsurface-temperature/upwelling feedback calculations. |
| `&um_param_ocn`, `&um_io_ocn` | `nfile_ocn_um`, `timename_ocn_um`, `varname_ocn_um`, `Lcycle_ocn_um`, `Tcycle_ocn_um`, `fnames_ocn_um(:)` | Mean zonal current, for example `u_ocn`. | `read_data_TLL_u` | Used as climatological zonal current in SST advection terms. |
| `&vm_param_ocn`, `&vm_io_ocn` | `nfile_ocn_vm`, `timename_ocn_vm`, `varname_ocn_vm`, `Lcycle_ocn_vm`, `Tcycle_ocn_vm`, `fnames_ocn_vm(:)` | Mean meridional current, for example `v_ocn`. | `read_data_TLL_v` | Used as climatological meridional current in SST advection terms. |
| `&wm_param_ocn`, `&wm_io_ocn` | `nfile_ocn_wm`, `timename_ocn_wm`, `varname_ocn_wm`, `Lcycle_ocn_wm`, `Tcycle_ocn_wm`, `fnames_ocn_wm(:)` | Mean vertical velocity, for example `w_ocn`. | `read_data_TLL_p` | Used as climatological vertical velocity in SST/upwelling terms. |
| `&Tzm_param_ocn`, `&Tzm_io_ocn` | `nfile_ocn_Tzm`, `timename_ocn_Tzm`, `varname_ocn_Tzm`, `Lcycle_ocn_Tzm`, `Tcycle_ocn_Tzm`, `fnames_ocn_Tzm(:)` | Mean vertical temperature-gradient field, for example `t_an`. | `read_data_TLL_p` | Used as `dTbdz` in vertical-advection/upwelling feedback terms. |

### Atmosphere Background Input Groups

| Namelist groups | Variables | Field read by the code | Grid reader | Role in `zc_cgcm_full.f90` |
| --- | --- | --- | --- | --- |
| `&sstm_param_atm`, `&sstm_io_atm` | `nfile_atm_sstm`, `timename_atm_sstm`, `varname_atm_sstm`, `Lcycle_atm_sstm`, `Tcycle_atm_sstm`, `fnames_atm_sstm(:)` | Atmosphere-grid mean/background SST, for example `sst`. | `read_data_TLL_atm` | Used as `agrd%sstm_atm`; combined with exchanged ocean SST anomalies for the atmosphere response. |
| `&uam_param_atm`, `&uam_io_atm` | `nfile_atm_uam`, `timename_atm_uam`, `varname_atm_uam`, `Lcycle_atm_uam`, `Tcycle_atm_uam`, `fnames_atm_uam(:)` | Atmosphere-grid mean zonal wind. | `read_data_TLL_atm` | Read only if `aset%heating_type=="ZC87_conv"`; used in the convective heating formulation. |
| `&vam_param_atm`, `&vam_io_atm` | `nfile_atm_vam`, `timename_atm_vam`, `varname_atm_vam`, `Lcycle_atm_vam`, `Tcycle_atm_vam`, `fnames_atm_vam(:)` | Atmosphere-grid mean meridional wind. | `read_data_TLL_atm` | Read only if `aset%heating_type=="ZC87_conv"`; used in the convective heating formulation. |

The default `do_cgcm_eqpac_30_clm_c1.4H120_dt3600_c10day.nml` does not set
`aset%heating_type`, so the default `"ZC87"` is used and the optional
`uam_*`/`vam_*` groups are not required.

### Coupling Sequence

At each coupling interval:

1. The atmosphere-grid mean SST is read/interpolated.
2. Ocean SST anomaly is mapped to the atmosphere grid with `exchange_OtoA`.
3. The Gill-type atmosphere solver computes `u_atm`, `v_atm`, and `p_atm`.
4. Atmosphere wind response is mapped back to the ocean grid with
   `exchange_AtoO`.
5. Ocean wind stress is computed from atmosphere wind anomalies and mean winds.
6. Ocean dynamics and the SST anomaly equation are advanced every ocean time
   step.

This sequence is inferred from `zc_cgcm_full.f90`.
