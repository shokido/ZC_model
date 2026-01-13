from ncdf_read import *
from ncdf_write import *
import netCDF4 as ncdf
import datetime as dt
import numpy as np
ref_dt=dt.datetime(1900,1,1,0,0,0)

fname_grid="../../INPUT/OCN/grid_eqpac_v1.nc"
percent=[80,90,100,110,120]
fnames_in=[]
fnames_out=[]
varnames=[]
percent_out=[]
for i in range(0,len(percent)):
    percent_out.append(percent[i])
    fnames_in.append("jra55do_v1.3_u10_clm_ocn_eqpac_v1.nc")
    fnames_out.append("jra55do_v1.3_u10_anm_clm"+str(int(percent[i]))+"per_ocn_eqpac_v1.nc")
    varnames.append("u10")
    percent_out.append(percent[i])
    fnames_in.append("jra55do_v1.3_v10_clm_ocn_eqpac_v1.nc")
    fnames_out.append("jra55do_v1.3_v10_anm_clm"+str(int(percent[i]))+"per_ocn_eqpac_v1.nc")
    varnames.append("v10")
nc_grid=ncdf.Dataset(fname_grid,"r")
lon_in=nc_grid.variables["lon_p"][:]
lat_in=nc_grid.variables["lat_p"][:]
nc_grid.close()

start_dt=dt.datetime(1970,1,1,0,0,0)
end_dt=dt.datetime(1970,12,31,0,0,0)
time_out=[]
for iy in range(start_dt.year,end_dt.year+1):
    for im in range(start_dt.month,end_dt.month+1):                    
        dt_out=dt.datetime(iy,im,1,0,0,0)
        time_out.append(1.0*(dt_out-ref_dt).days)
time_out=np.asarray(time_out)
ref_units="days since "+str(ref_dt)
ntime_out=len(time_out)
nx=len(lon_in)
ny=len(lat_in)
lon_out,lat_out=np.meshgrid(lon_in,lat_in)
for ifile in range(0,len(fnames_in)):
    nc_out=ncdf.Dataset(fnames_out[ifile],"w")
    nc_out.createDimension("x",nx)
    nc_out.createDimension("y",ny)
    nc_out.createDimension("time",ntime_out)
    nc_out.createVariable("lon","double",("y","x"))
    nc_out.createVariable("lat","double",("y","x"))
    nc_out.createVariable("time","double",("time"))
    nc_out.createVariable(varnames[ifile],"double",("time","y","x"))
    nc_out["lon"].units="degrees_east"
    nc_out["lat"].units="degrees_north"
    nc_out["time"].units=ref_units
    nc_out[varnames[ifile]].units="m/s"

    nc_out["lon"][:]=lon_out
    nc_out["lat"][:]=lat_out
    nc_out["time"][:]=time_out
    nc_in=ncdf.Dataset(fnames_in[ifile],"r")
    var_in=nc_in.variables[varnames[ifile]][:,:,:]
    nc_in.close()
    for it in range(0,ntime_out):
        it_in=it%12
        print(it_in)
        tmp=var_in[it_in,:,:]
        nc_out[varnames[ifile]][it,:,:]=tmp*(0.01*percent_out[ifile]-1)
    nc_out.close()
