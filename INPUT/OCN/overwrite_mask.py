# Python script for construct LOM grid from ROMS grid
import netCDF4 as ncdf
import numpy as np

grd_name="eqpac_v2"
grd_name="eqpac_ref"
fname_grd="grid_"+grd_name+".nc"
fname_mask="mask_"+grd_name+"_rev.nc"
fname_mask_sst="mask_"+grd_name+"_rev_sst.nc"

nc_grd=ncdf.Dataset(fname_grd,"r+")
nc_mask=ncdf.Dataset(fname_mask,"r")
nc_mask_sst=ncdf.Dataset(fname_mask_sst,"r")
mask_p=nc_mask.variables["mask"][:,:]
mask_sst=nc_mask_sst.variables["mask"][:,:]
nc_grd.variables["mask_p"][:,:]=mask_p
nc_grd.variables["mask_sst"][:,:]=mask_sst
nc_grd.close()
nc_mask.close()
nc_mask_sst.close()
