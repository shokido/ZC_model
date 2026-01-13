import netCDF4 as ncdf
import numpy as np
import matplotlib.pyplot as plt
from stat_ncl import *
from ncdf_read import *
dt1=dt.datetime(1970,1,1,0,0,0)
dt2=dt.datetime(2014,12,31,0,0,0)
dt1_1=dt.datetime(dt1.year,1,3,0,0,0)
dt2_2=dt.datetime(dt2.year+1,1,1,0,0,0)
#dt2=dt.datetime(1989,12,31,0,0,0)
lon1=120;lon2=270;lat1=-2;lat2=2
fname_model="../../OUTPUTS/avg_spinup_ann_ocn_1_20.nc"
fname_dta1="data_mean_clm_ocn_120.nc";hm1=120
fname_dta2="data_mean_clm_ocn_1_20.nc";hm2=150

nc_model=ncdf.Dataset(fname_model,"r")
lon=nc_model.variables["lon_p"][:]
lat=nc_model.variables["lat_p"][:]
time=nc_model.variables["time"][:]
time_u=nc_model.variables["time"].units
cal=ncdf.num2date(time,units=time_u,calendar="standard")
ical=np.where((cal>=dt1_1)&(cal<=dt2_2))[0]
ilon=np.where((lon>=lon1)&(lon<=lon2))[0]
ilat=np.where((lat>=lat1)&(lat<=lat2))[0]
taux=nc_model.variables["taux"][ical,ilat,ilon]
h=nc_model.variables["h_ocn_sw"][ical,ilat,ilon]
x=nc_model.variables["x_p_2d"][ilat,ilon]
nc_model.close()

nc_dta1=ncdf.Dataset(fname_dta1,"r")
h1=nc_dta1.variables["h_ocn"][:,ilat,ilon]
nc_dta1.close()
nc_dta2=ncdf.Dataset(fname_dta2,"r")
h2=nc_dta2.variables["h_ocn"][:,ilat,ilon]
nc_dta2.close()
h1=np.mean(np.mean(h1,axis=1),axis=0)
h2=np.mean(np.mean(h2,axis=1),axis=0)
x=np.mean(x,axis=0)
plt.subplot(1,1,1)
plt.plot(lon[ilon],h1)
plt.plot(lon[ilon],h2)
plt.show()
