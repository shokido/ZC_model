import subprocess
import glob

grd_name="eqatm_256"
fname_gdes="../../INPUT/ATM/des_256.txt"
grd_name="eqpac_v1"
fname_gdes="../../INPUT/OCN/des_grid_"+grd_name+".dat"
#fnames_atm=["/Volumes/HV620S/DATA/WIND/JRA55do/jra55do_v1.3_u10_clm_monthly.nc","/Volumes/HV620S/DATA/WIND/JRA55do/jra55do_v1.3_v10_clm_monthly.nc"]

iy1=1958
iy2=2015
fnames_in=[];fnames_out=[]
dir_in="/latent/kido/ATM_REANALYSIS/JRA55do/version1.3/Output/Monthly/"
fnames_in.append(dir_in+"jra55do_v1.3_uw_monthly_clm.nc")
fnames_out.append("jra55do_v1.3_uw_clm_ocn_"+grd_name+".nc")
fnames_in.append(dir_in+"jra55do_v1.3_vw_monthly_clm.nc")
fnames_out.append("jra55do_v1.3_vw_clm_ocn_"+grd_name+".nc")
dir_in="/latent/kido/ATM_REANALYSIS/JRA55do/version1.3/Output/Monthly/"
for iy in range(iy1,iy2):
   fnames_in.append(dir_in+"jra55do_v1.3_uw_"+str(iy)+"_monthly.nc")
   fnames_out.append("jra55do_v1.3_uw_"+str(iy)+"_ocn_"+grd_name+".nc")
   fnames_in.append(dir_in+"jra55do_v1.3_vw_"+str(iy)+"_monthly.nc")
   fnames_out.append("jra55do_v1.3_vw_"+str(iy)+"_ocn_"+grd_name+".nc")

for i in range(0,len(fnames_in)):
    command="cdo remapbil,"+fname_gdes+" "+fnames_in[i]+" "+fnames_out[i]
    subprocess.call(command.split())
