import netCDF4 as ncdf
import numpy as np
import matplotlib.pyplot as plt
from stat_ncl import *
import datetime as dt
fnames_in=["ERSST_v5_eqpac_v1_clm.nc","ERSST_v5_atm_256_clm.nc"]
ref_dt=dt.datetime(1900,1,1,0,0,0)
ref_units="days since "+str(ref_dt)
target_year=1970
for i in range(0,len(fnames_in)):
    nc_in=ncdf.Dataset(fnames_in[i],"r+")
    time=nc_in.variables["time"][:]
    time_units=nc_in.variables["time"].units
    cal_in=ncdf.num2date(time,units=time_units,calendar="standard")
    nc_in.variables["time"].units=ref_units
    for itime in range(0,len(cal_in)):
        tmp=cal_in[itime]
        tmp_dt=dt.datetime(target_year,tmp.month,tmp.day,tmp.hour,tmp.minute,tmp.second)
        print(tmp_dt)
        time_out=1.0*(tmp_dt-ref_dt).days
        nc_in.variables["time"][itime]=time_out
    nc_in.close()
