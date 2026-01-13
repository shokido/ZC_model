program zc_ogcm_dyn
  ! Fortran code for running the OGCM part of the ZC model without SST equation
  ! Mainly used for obtaining the mean fields u,v,h,and w, which is used for subsequent SST calculation
  ! Input:
  ! 1. Wind anomalies
  ! Output
  ! 1: Horizontal and vertical current (u,v,w) and thermocline depth (h)
  use run_params;use run_types
  use calendar_sub
  use mod_io_master
  use mod_io_avg
  use mod_io_diag
  use mod_wstress
  use mod_ocn_solver_zc
  use mod_ocn_dta
  implicit none
  type(ocn_dta) :: ogrd
  type(ocn_set) :: oset
  character(len=maxlen) :: fname_grd_ocn
  character(len=maxlen) :: flag_ini_ocn,fname_ini_ocn
  character(len=maxlen) :: fname_avg_ocn
  character(len=maxlen) :: fname_diag_ocn
  character(len=maxlen) :: fname_rst_ocn
  ! Wind forcing (assumed to be surface wind)
  ! Zonal
  type(TLL_dta) :: ocn_taux_dta
  integer :: nfile_ocn_taux
  character(len=maxlen),allocatable :: fnames_ocn_taux(:)
  character(len=maxlen) :: timename_ocn_taux,varname_ocn_taux
  character(1) :: Lcycle_ocn_taux
  real(idx) :: Tcycle_ocn_taux
  ! Meridional
  type(TLL_dta) :: ocn_tauy_dta
  integer :: nfile_ocn_tauy
  character(len=maxlen),allocatable :: fnames_ocn_tauy(:)
  character(len=maxlen) :: timename_ocn_tauy,varname_ocn_tauy
  character(1) :: Lcycle_ocn_tauy
  real(idx) :: Tcycle_ocn_tauy

  integer :: itime,ntime,total_time
  real(idx) :: dt,tmp1,time_int
  integer :: start_yymmdd,start_hhmmss
  integer :: end_yymmdd,end_hhmmss
  integer :: ix,iy
  integer :: tmp_yymmdd,tmp_hhmmss

  namelist/date/dt,start_yymmdd,start_hhmmss,end_yymmdd,end_hhmmss
  namelist/io_ocn/fname_grd_ocn
  namelist/io_ocn/flag_ini_ocn,fname_ini_ocn
  namelist/io_ocn/fname_avg_ocn,out_avg_flag,out_avg_int
  namelist/io_ocn/fname_rst_ocn
  namelist/param_ocn/oset
  
  namelist/taux_param_ocn/nfile_ocn_taux,timename_ocn_taux,varname_ocn_taux,Lcycle_ocn_taux,Tcycle_ocn_taux
  namelist/taux_io_ocn/fnames_ocn_taux
  namelist/tauy_param_ocn/nfile_ocn_tauy,timename_ocn_tauy,varname_ocn_tauy,Lcycle_ocn_tauy,Tcycle_ocn_tauy
  namelist/tauy_io_ocn/fnames_ocn_tauy

  read(5,date)
  read(5,io_ocn)
  read(5,param_ocn)
  read(5,taux_param_ocn)
  allocate(fnames_ocn_taux(nfile_ocn_taux))
  read(5,taux_io_ocn)
  read(5,tauy_param_ocn)
  allocate(fnames_ocn_tauy(nfile_ocn_tauy))
  read(5,tauy_io_ocn)
  ! Parameter
  ! Time setting
  call calendar_cal_length_ymdhms(start_yymmdd,start_hhmmss,end_yymmdd,end_hhmmss,1,tmp1)
  total_time=int(tmp1);ntime=int(tmp1/(dt*sec_to_day))

  ! Read ocean grid
  call read_ocn_dyn_grd(fname_grd_ocn,ogrd)
  ! Set masking
  call set_mask_ocn(ogrd,oset%slip_ind)
  call initialize_ocn_dyn(ogrd)
  call initialize_ocn_visc(ogrd,oset)
  if (flag_ini_ocn == "T") then
     call read_restart_ocn_dyn(fname_ini_ocn,ogrd)
  end if
  ! Read forcing fields
  call read_data_TLL_p(nfile_ocn_taux,fnames_ocn_taux,timename_ocn_taux,varname_ocn_taux,&
       & ogrd,ocn_taux_dta,start_yymmdd,start_hhmmss)
  ocn_taux_dta%Lcycle=Lcycle_ocn_taux;ocn_taux_dta%Tcycle=Tcycle_ocn_taux
  call read_data_TLL_p(nfile_ocn_tauy,fnames_ocn_tauy,timename_ocn_tauy,varname_ocn_tauy,&
       & ogrd,ocn_tauy_dta,start_yymmdd,start_hhmmss)
  ocn_tauy_dta%Lcycle=Lcycle_ocn_tauy;ocn_tauy_dta%Tcycle=Tcycle_ocn_tauy
  ! Prepare averaged file
  call create_avg_ocn_dyn_ZC(fname_avg_ocn,ogrd,oset,&
       & start_yymmdd,start_hhmmss,end_yymmdd,end_hhmmss,dt,&
       & out_avg_flag,out_avg_int,missing_value,istep_avg,ntime_avg)
  call initialize_ocn_dyn_avg(ogrd)
  ! Prepare diagnostic file
  ! Write output file
  iavg=1;iavg_count=0
  do itime=1,ntime
     time_int=dt*itime*sec_to_day
     iavg_count = iavg_count + 1
     ! Get uw and vw
     call get_data_TLL_p(time_int,ogrd,ocn_taux_dta)
     call get_data_TLL_p(time_int,ogrd,ocn_tauy_dta)
#if defined WSTRESS_BULK
     ! Read surface wind
     ogrd%ua_ocn%val=ocn_taux_dta%data_now%val
     ogrd%va_ocn%val=ocn_tauy_dta%data_now%val
     call ua_to_stress_total(ogrd,oset)
#else
     ! Read surface wind stress
     ogrd%taux_ocn%val=ocn_taux_dta%data_now%val
     ogrd%tauy_ocn%val=ocn_tauy_dta%data_now%val
#endif
     ! Calculate Ekman currents
     call solve_ekman_ocn(ogrd,oset)
     ! Calculate geostrophic current (with time stepping)
     call solve_rg_vgeo_ocn(ogrd,oset,dt)
     ! Calculate total current
     call solve_totalcurrent_ocn(ogrd,oset)
     ! Add snapshot to mean fiedls
     call oper_avg_ocn_dyn(ogrd)
     ! I/O operation
     if (itime .eq. istep_avg(iavg)) then
       call calendar_cal_ymdhms_after(start_yymmdd,start_hhmmss,dt*itime*sec_to_day,1,tmp_yymmdd,tmp_hhmmss)
       write(*,*) "Step (average) =",iavg," ",tmp_yymmdd,tmp_hhmmss,iavg_count
        call write_avg_ocn_dyn(fname_avg_ocn,ogrd,iavg_count,iavg)
        write(*,*) maxval(ogrd%u_sw%val),minval(ogrd%u_sw%val)
        iavg_count=0
        iavg=iavg+1
        iavg=min(iavg,ntime_avg)
     end if
  end do
  ! Write restart file
  call write_restart_ocn_dyn(fname_rst_ocn,ogrd,missing_value)
  call deallocate_ocn_all(ogrd)
  write(*,*) "Finish run, ended at ",end_yymmdd,end_hhmmss
  deallocate(fnames_ocn_taux)
  deallocate(fnames_ocn_tauy)
end program zc_ogcm_dyn
