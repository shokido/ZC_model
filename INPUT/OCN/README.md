### Preprocessing of ocean gridfile
#### Step 1: Define lon-lat grid
First, generate a curvilinear ocean grid template based on a specified longitude–latitude domain and horizontal resolution. It creates: (1) A NetCDF grid file containing p-, u-, and v-grid coordinates, horizontal distances (1D and 2D), Coriolis parameter, masks and damping coefficients, and (2) A grid description (.dat) file compatible with remapping tools (e.g., CDO)

- The grid name and geographic domain are defined at the beginning of the script and are used consistently to construct output filenames.

#### How to run
Edit the grid configuration parameters in the script header, then run:
```bash 
python make_grid.py
```

The following parameters define the grid name, geographic extent, resolution, and output files:
<details>
<summary>List of variable</summary>

Variable	|Example value	| Description
---|---|---
gname |	"eqpac_ref"	| Grid name identifier. Used as a suffix in output filenames.
lon_w	|110.0	|Western boundary longitude (degrees east).
lon_e	|300.0	|Eastern boundary longitude (degrees east).
lat_s	|-30.0	|Southern boundary latitude (degrees north).
lat_n	|30.0	|Northern boundary latitude (degrees north).
dlon	|2.5	|Zonal grid spacing in degrees.
dlat	|1.0	 |Meridional grid spacing in degrees.
fname_grid	|"grid_eqpac_ref.nc"	|Output NetCDF file containing the ocean grid definition.
fname_gdes	|"des_grid_eqpac_ref.dat"	|Output grid description file for remapping utilities.

</details>


-   Notes
    -   Longitudes are assumed to be in degrees east and may exceed 180° (e.g., 0–360 convention).
    -   Grid spacing is uniform in longitude–latitude space, while physical distances are computed using spherical geometry.
    -   The script currently produces a template grid; land/sea masking and damping coefficients are initialized uniformly and can be modified later.

#### Step 2: Set tentative lon-lat boundary
Generates a tentative land–sea mask on the model grid by remapping an external topography dataset onto the ocean grid and applying a depth threshold.

The workflow is as follows:
1. Remap the source topography file onto the target ocean grid using bilinear interpolation (cdo remapbil).
1. Read the remapped topography and grid coordinates.
1. Set grid points shallower than a specified depth threshold to land (mask = 0).
1. Output a NetCDF mask file containing longitude, latitude, topography, and land–sea mask fields.

The resulting mask file can be used directly in ocean model preprocessing.
edit create_mask_from_topo.py
```bash
python create_mask_from_topo.py
```
<details>
<summary>List of variable</summary>

Variable	|Example value	| Description
---|---|---
grd_name	| "eqpac_ref"	|  Grid name identifier. Used consistently in input and output filenames.
fname_gdes	| "des_grid_eqpac_ref.dat"	| Grid description file used by CDO for remapping.
fname_mask	| "mask_eqpac_ref.nc"| 	Output NetCDF file containing the land–sea mask and remapped topography.
fname_topo	| "etopo5.nc"	|  Input topography dataset (NetCDF).
fname_grid	| "grid_eqpac_ref.nc"	| Ocean grid file defining the target grid coordinates.
varname	 | "topo"	| Variable name of topography in the input topography file.
lonname	 | "lon_p"	| Longitude variable name (updated after reading the grid file).
latname	| "lat_p"	| Latitude variable name (updated after reading the grid file).
hthres	| -100.0| 	Depth threshold (meters). Grid points shallower than this value are classified as land.
</details>


-   Notes
    -   Bilinear interpolation is performed using CDO, which must be installed and available in the system path.
    -   Depth values are assumed to be negative below sea level.
    -   The mask convention is:
    mask = 1 : ocean
    mask = 0 : land
#### Step 3: Manually edit p-grid
