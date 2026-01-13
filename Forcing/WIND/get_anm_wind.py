import subprocess
import glob

grd_name="eqpac_v1"
fname_gdes="../../INPUT/OCN/des_grid_"+grd_name+".dat"
iy1=1958
iy2=2015
fnames_in=[];
fnames_clm=[]
fnames_out=[]
dir_in="/latent/kido/ATM_REANALYSIS/JRA55do/version1.3/Output/Monthly/"
for iy in range(iy1,iy2):
    fnames_in.append("jra55do_v1.3_u10_"+str(iy)+"_ocn_"+grd_name+".nc")
    fnames_clm.append("jra55do_v1.3_u10_clm_ocn_"+grd_name+".nc")
    fnames_out.append("jra55do_v1.3_u10_"+str(iy)+"_ocn_"+grd_name+"_anm.nc")
    fnames_in.append("jra55do_v1.3_v10_"+str(iy)+"_ocn_"+grd_name+".nc")
    fnames_clm.append("jra55do_v1.3_v10_clm_ocn_"+grd_name+".nc")
    fnames_out.append("jra55do_v1.3_v10_"+str(iy)+"_ocn_"+grd_name+"_anm.nc")
for i in range(0,len(fnames_in)):
    command="cdo sub "+fnames_in[i]+" "+fnames_clm[i]+" "+fnames_out[i]
    subprocess.call(command.split())
