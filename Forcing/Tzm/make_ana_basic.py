import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as ncdf
import datetime as dt
fname_grid="../../INPUT/OCN/grid_eqpac_v1.nc"
fname_out="basic_2D_ZC87.nc"

nc_grid=ncdf.Dataset(fname_grid,"r")
x_p=nc_grid.variables["x_p_2d"][:]
y_p=nc_grid.variables["y_p_2d"][:]
lon_basic=nc_grid.variables["lon_p"][:]
lat_basic=nc_grid.variables["lat_p"][:]
mask_basic=nc_grid.variables["mask_p"][:]
nc_grid.close()


lon_grid=np.linspace(129.375,360-84.375,27)
nlon=len(lon_grid)
Tzm=np.zeros(nlon)
hbar=np.zeros(nlon)
for i in range(0,10):
    Tzm[i]=2.0/100
    hbar[i]=165
i=10;Tzm[i]=1.8/100;hbar[i]=170
i=11;Tzm[i]=1.5/100;hbar[i]=175
i=12;Tzm[i]=1.3/100;hbar[i]=175
i=13;Tzm[i]=2.0/100;hbar[i]=175
i=14;Tzm[i]=2.5/100;hbar[i]=150
i=15;Tzm[i]=3.0/100;hbar[i]=125
i=16;Tzm[i]=3.5/100;hbar[i]=100
i=17;Tzm[i]=4.0/100;hbar[i]=93
i=18;Tzm[i]=4.5/100;hbar[i]=86
i=19;Tzm[i]=5.0/100;hbar[i]=80
i=20;Tzm[i]=5.5/100;hbar[i]=73
i=21;Tzm[i]=5.5/100;hbar[i]=66
i=22;Tzm[i]=5.5/100;hbar[i]=60
i=23;Tzm[i]=5.5/100;hbar[i]=52
i=24;Tzm[i]=5.5/100;hbar[i]=50
i=25;Tzm[i]=5.5/100;hbar[i]=50
i=26;Tzm[i]=5.5/100;hbar[i]=50
hbar_basic=np.interp(lon_basic,lon_grid,hbar)
Tzm_basic=np.interp(lon_basic,lon_grid,Tzm)

hbar_basic_2D=np.zeros((len(lat_basic),len(lon_basic)))
Tzm_basic_2D=np.zeros((len(lat_basic),len(lon_basic)))
for ilat in range(0,len(lat_basic)):
    hbar_basic_2D[ilat,:]=hbar_basic
    Tzm_basic_2D[ilat,:]=Tzm_basic

miss=-99.999
Tzm_basic_2D=Tzm_basic_2D*mask_basic+(1-mask_basic)*miss
hbar_basic_2D=hbar_basic_2D*mask_basic+(1-mask_basic)*miss
time_out=[]
ref_dt=dt.datetime(1900,1,1,0,0,0)
ref_units="days since "+str(ref_dt)
for i in range(1,2):
    tmp=dt.datetime(1970,i,15,0,0,0)-ref_dt
    time_out.append(tmp.days)
time_out=np.asarray(time_out)                    

nc_out=ncdf.Dataset(fname_out,"w")
ndim_p=np.shape(x_p)
nc_out.createDimension("x_p",ndim_p[1])
nc_out.createDimension("y_p",ndim_p[0])
nc_out.createDimension("time",1)
nc_out.createVariable("x_p_2d",time_out.dtype,("y_p","x_p"))
nc_out["x_p_2d"][:]=x_p[:]
nc_out.createVariable("y_p_2d",time_out.dtype,("y_p","x_p"))
nc_out["y_p_2d"][:]=y_p[:]
nc_out.createVariable("time",time_out.dtype,("time"))
nc_out["time"].units=ref_units
nc_out["time"][:]=time_out[0:1]
nc_out.createVariable("h_ocn",hbar_basic_2D.dtype,("time","y_p","x_p"))
nc_out.createVariable("t_an",Tzm_basic_2D.dtype,("time","y_p","x_p"))
nc_out["h_ocn"].units="m"
nc_out["t_an"].units="degrees_celsius/m"
nc_out["h_ocn"].missing_value=miss
nc_out["t_an"].missing_value=miss
nc_out["h_ocn"][:]=hbar_basic_2D
nc_out["t_an"][:]=Tzm_basic_2D
nc_out.close()
