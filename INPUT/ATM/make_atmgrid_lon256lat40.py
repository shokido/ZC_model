import numpy as np
from netCDF4 import Dataset

out_nc = "grid_256.nc"

# -----------------------------
# 1) Create coordinate
# -----------------------------
# lat: -40, -39, ..., 39, 40  (81 points)
lat = np.arange(-40, 41, 1, dtype=np.float32)

# lon: divide 0〜360 into 256 points （Include edges）
lon = np.linspace(0.0, 360.0, 256, dtype=np.float32)
# -----------------------------
# 2) Output netCDF
# -----------------------------
with Dataset(out_nc, "w", format="NETCDF4_CLASSIC") as nc:
    # dimensions
    nc.createDimension("lon", lon.size)
    nc.createDimension("lat", lat.size)
    # variables
    vlon = nc.createVariable("lon", "f4", ("lon",))
    vlat = nc.createVariable("lat", "f4", ("lat",))
    # attributes
    vlon.units = "degrees_east"
    vlat.units = "degrees_north"
    # data
    vlon[:] = lon
    vlat[:] = lat

