import subprocess
import glob

grd_name="eqpac_v1";fname_gdes="../..//INPUT/OCN/des_grid_"+grd_name+".dat"
grd_name="atm_256";fname_gdes="../../INPUT/ATM/des_256.txt"
fname_in="ERSST_v5_1958_2019.nc"
fname_clm="ERSST_v5_"+str(grd_name)+"_clm.nc"
fname_anm="ERSST_v5_"+str(grd_name)+"_anm.nc"

iy1=1958
iy2=2015
command="cdo -selyear,"+str(iy1)+"/"+str(iy2)+" "+fname_in+" a.nc"
subprocess.call(command.split())
command="cdo remapbil,"+fname_gdes+" a.nc tmp.nc"
subprocess.call(command.split())
command="mv tmp.nc a.nc"
subprocess.call(command.split())
command="cdo ymonmean a.nc "+fname_clm
subprocess.call(command.split())
command="cdo setmisstonn "+fname_clm+" b.nc"
subprocess.call(command.split())
command="mv b.nc "+fname_clm
subprocess.call(command.split())
command="cdo ymonsub a.nc "+fname_clm+" "+fname_anm
subprocess.call(command.split())
command="cdo setmisstonn "+fname_anm+" b.nc"
subprocess.call(command.split())
command="mv b.nc "+fname_anm
subprocess.call(command.split())
command="rm a.nc"
subprocess.call(command.split())
