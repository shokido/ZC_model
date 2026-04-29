import netCDF4 as ncdf
import numpy as np
lon1=190;lon2=240;lat1=-5;lat2=5
grd_name="eqpac_30"
fname_grid="grid_"+grd_name+".nc"
nc1=ncdf.Dataset(fname_grid,"r+")
lon_p=nc1.variables["lon_p"][:]
lat_p=nc1.variables["lat_p"][:]
nlon=len(lon_p);nlat=len(lat_p)
mask_wwb=np.zeros((nlat,nlon))
ilon=np.where((lon_p>=lon1)&(lon_p<=lon2))[0]
ilat=np.where((lat_p>=lat1)&(lat_p<=lat2))[0]
for i in ilat:
    mask_wwb[i,ilon]=1
if "mask_wwb" not in nc1.variables:
    nc1.createVariable("mask_wwb", "f4", ("y_p", "x_p"))
nc1.variables["mask_wwb"][:]=mask_wwb
nc1.close()
