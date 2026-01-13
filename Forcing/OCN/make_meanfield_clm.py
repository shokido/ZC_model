import netCDF4 as ncdf
import numpy as np
import matplotlib.pyplot as plt
from stat_ncl import *
import datetime as dt
dir_io="../..//OUTPUTS/SPIN/"
fflag="bulk_H130_cd_2.0_1_20"
#fflag="bulk_H150_cd_2.1_1_20"
#fflag="bulk_H150_cd_2.25_1_20"
#fflag="bulk_H150_cd_2.5_1_20"
fflag="bulk_H120_cd_1.7_1_20"
#fflag="bulk_H150_cd_1.8_1_20"
#fflag="bulk_H120_cd_1.8_1_20"

fname_base=dir_io+"avg_spinup_clm_"+fflag+".nc"
fname_base_sst="../SST/ERSST_v5_eqpac_v1_clm.nc"
fname_out="basic_clm_"+fflag+".nc"

nc_base=ncdf.Dataset(fname_base,"r")
lon=nc_base.variables["lon_p"][:]
lat=nc_base.variables["lat_p"][:]
time=nc_base.variables["time"][:]
time_units=nc_base.variables["time"].units
x_p=nc_base.variables["x_p_2d"][:]
y_p=nc_base.variables["y_p_2d"][:]
x_u=nc_base.variables["x_u_2d"][:]
y_u=nc_base.variables["y_u_2d"][:]
x_v=nc_base.variables["x_v_2d"][:]
y_v=nc_base.variables["y_v_2d"][:]
u=nc_base.variables["u_ocn_1"][12:,:,:]
v=nc_base.variables["v_ocn_1"][12:,:,:]
w=nc_base.variables["w_ocn_1"][12:,:,:]
h=nc_base.variables["h_ocn_sw"][12:,:,:]
miss=nc_base.variables["h_ocn_sw"].missing_value
H1_base=nc_base.variables["time"].H1
H2_base=nc_base.variables["time"].H2
nc_base.close()
nc_base_sst=ncdf.Dataset(fname_base_sst,"r")
sst_base=nc_base_sst["sst"][:]
nc_base.close()

sst_base_ann=clmmon(sst_base)
time_out=[]
ref_dt=dt.datetime(1900,1,1,0,0,0)
ref_units="days since "+str(ref_dt)
for i in range(1,13):
    tmp=dt.datetime(1970,i,15,0,0,0)-ref_dt
    time_out.append(tmp.days)
time_out=np.asarray(time_out)                    

u_base_ann=clmmon(u)
v_base_ann=clmmon(v)
w_base_ann=clmmon(w)
h_base_ann=clmmon(h)
# ilon=np.where((lon>=120) & (lon<=280))[0]
# ilat=np.where(np.abs(lat)<=10)[0]
# view=h_base_ann[:,ilat,:]
# print(np.nanmean(view[:,:,ilon]))
nc_out=ncdf.Dataset(fname_out,"w")
ndim_p=np.shape(x_p)
nc_out.createDimension("x_p",ndim_p[1])
nc_out.createDimension("y_p",ndim_p[0])
ndim_u=np.shape(x_u)
nc_out.createDimension("x_u",ndim_u[1])
nc_out.createDimension("y_u",ndim_u[0])
ndim_v=np.shape(x_v)
nc_out.createDimension("x_v",ndim_v[1])
nc_out.createDimension("y_v",ndim_v[0])
nc_out.createDimension("time",12)
nc_out.createVariable("x_p_2d",time.dtype,("y_p","x_p"))
nc_out["x_p_2d"][:]=x_p[:]
nc_out.createVariable("y_p_2d",time.dtype,("y_p","x_p"))
nc_out["y_p_2d"][:]=y_p[:]
nc_out.createVariable("x_u_2d",time.dtype,("y_u","x_u"))
nc_out["x_u_2d"][:]=x_u[:]
nc_out.createVariable("y_u_2d",time.dtype,("y_u","x_u"))
nc_out["y_u_2d"][:]=y_u[:]
nc_out.createVariable("x_v_2d",time.dtype,("y_v","x_v"))
nc_out["x_v_2d"][:]=x_v[:]
nc_out.createVariable("y_v_2d",time.dtype,("y_v","x_v"))
nc_out["y_v_2d"][:]=y_v[:]
nc_out.createVariable("time",time.dtype,("time"))
nc_out["time"].units=ref_units
nc_out["time"][:]=time_out[0:12]

nc_out.createVariable("u_ocn",u.dtype,("time","y_u","x_u"))
nc_out.createVariable("v_ocn",v.dtype,("time","y_v","x_v"))
nc_out.createVariable("w_ocn",w.dtype,("time","y_p","x_p"))
nc_out.createVariable("h_ocn",h.dtype,("time","y_p","x_p"))
nc_out.createVariable("sst_ocn",h.dtype,("time","y_p","x_p"))
nc_out["u_ocn"].units="m/s"
nc_out["v_ocn"].units="m/s"
nc_out["w_ocn"].units="m/s"
nc_out["h_ocn"].units="m/s"
nc_out["sst_ocn"].units="degrees_celsius"
nc_out["u_ocn"].missing_value=miss
nc_out["v_ocn"].missing_value=miss
nc_out["w_ocn"].missing_value=miss
nc_out["h_ocn"].missing_value=miss
nc_out["sst_ocn"].missing_value=miss

nc_out["u_ocn"][:]=u_base_ann[:]
nc_out["v_ocn"][:]=v_base_ann[:]
nc_out["w_ocn"][:]=w_base_ann[:]
h_out=h_base_ann[:]+(H1_base+H2_base)
h_out[h_out<0]=10.0
nc_out["h_ocn"][:]=h_out
nc_out["sst_ocn"][:]=sst_base_ann[:]
nc_out.close()
