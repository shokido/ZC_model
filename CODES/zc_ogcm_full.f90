program zc_ogcm_full
  ! Fortran code for running the OGCM part of the ZC model (anomaly model)
  ! Input:
  ! 1. Interannualy varying wind stress (zonal and meridional)
  ! 2. Mean SST
  ! 3. Mean wind stress
  ! 4: Mean oceaninc fields (Tz, u,v,w, and h)
  ! Output
  ! 1: Anomalous Horizontal and vertical current (u,v,w) and thermocline depth (h)
  ! 2: SST anomalies
  ! 3: Diagnostic file for SST
  use run_params
  use run_types
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
  ! Wind forcing (assumed to be surface wind anomalies)
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
  type(TLL_dta) :: ocn_sstm_dta
  integer :: nfile_ocn_sstm
  character(len=maxlen),allocatable :: fnames_ocn_sstm(:)
  character(len=maxlen) :: timename_ocn_sstm,varname_ocn_sstm
  character(1) :: Lcycle_ocn_sstm
  real(idx) :: Tcycle_ocn_sstm
  ! Mean zonal wind
  type(TLL_dta) :: ocn_tauxm_dta
  integer :: nfile_ocn_tauxm
  character(len=maxlen),allocatable :: fnames_ocn_tauxm(:)
  character(len=maxlen) :: timename_ocn_tauxm,varname_ocn_tauxm
  character(1) :: Lcycle_ocn_tauxm
  real(idx) :: Tcycle_ocn_tauxm
  ! Mean meridonal wind
  type(TLL_dta) :: ocn_tauym_dta
  integer :: nfile_ocn_tauym
  character(len=maxlen),allocatable :: fnames_ocn_tauym(:)
  character(len=maxlen) :: timename_ocn_tauym,varname_ocn_tauym
  character(1) :: Lcycle_ocn_tauym
  real(idx) :: Tcycle_ocn_tauym

  ! Mean oceanic fields
  ! Themocline depth (hbar)
  type(TLL_dta) :: ocn_hm_dta
  integer :: nfile_ocn_hm
  character(len=maxlen),allocatable :: fnames_ocn_hm(:)
  character(len=maxlen) :: timename_ocn_hm,varname_ocn_hm
  character(1) :: Lcycle_ocn_hm
  real(idx) :: Tcycle_ocn_hm
  ! Surface zonal current (ubar)
  type(TLL_dta) :: ocn_um_dta
  integer :: nfile_ocn_um
  character(len=maxlen),allocatable :: fnames_ocn_um(:)
  character(len=maxlen) :: timename_ocn_um,varname_ocn_um
  character(1) :: Lcycle_ocn_um
  real(idx) :: Tcycle_ocn_um
  ! Surface meridional current (vbar)
  type(TLL_dta) :: ocn_vm_dta
  integer :: nfile_ocn_vm
  character(len=maxlen),allocatable :: fnames_ocn_vm(:)
  character(len=maxlen) :: timename_ocn_vm,varname_ocn_vm
  character(1) :: Lcycle_ocn_vm
  real(idx) :: Tcycle_ocn_vm
  ! Vertical velocity (wbar)
  type(TLL_dta) :: ocn_wm_dta
  integer :: nfile_ocn_wm
  character(len=maxlen),allocatable :: fnames_ocn_wm(:)
  character(len=maxlen) :: timename_ocn_wm,varname_ocn_wm
  character(1) :: Lcycle_ocn_wm
  real(idx) :: Tcycle_ocn_wm
  ! Mean temperature gradient at the bottom of the mixed layer (Tzbar)  
  type(TLL_dta) :: ocn_Tzm_dta
  integer :: nfile_ocn_Tzm
  character(len=maxlen),allocatable :: fnames_ocn_Tzm(:)
  character(len=maxlen) :: timename_ocn_Tzm,varname_ocn_Tzm
  character(1) :: Lcycle_ocn_Tzm
  real(idx) :: Tcycle_ocn_Tzm

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
  namelist/io_ocn/fname_diag_ocn,out_diag_flag,out_diag_int
  namelist/io_ocn/fname_rst_ocn
  namelist/param_ocn/oset
  
  namelist/taux_param_ocn/nfile_ocn_taux,timename_ocn_taux,varname_ocn_taux,Lcycle_ocn_taux,Tcycle_ocn_taux
  namelist/taux_io_ocn/fnames_ocn_taux
  namelist/tauy_param_ocn/nfile_ocn_tauy,timename_ocn_tauy,varname_ocn_tauy,Lcycle_ocn_tauy,Tcycle_ocn_tauy
  namelist/tauy_io_ocn/fnames_ocn_tauy

  namelist/sstm_param_ocn/nfile_ocn_sstm,timename_ocn_sstm,varname_ocn_sstm,Lcycle_ocn_sstm,Tcycle_ocn_sstm
  namelist/sstm_io_ocn/fnames_ocn_sstm
  namelist/tauxm_param_ocn/nfile_ocn_tauxm,timename_ocn_tauxm,varname_ocn_tauxm,Lcycle_ocn_tauxm,Tcycle_ocn_tauxm
  namelist/tauxm_io_ocn/fnames_ocn_tauxm
  namelist/tauym_param_ocn/nfile_ocn_tauym,timename_ocn_tauym,varname_ocn_tauym,Lcycle_ocn_tauym,Tcycle_ocn_tauym
  namelist/tauym_io_ocn/fnames_ocn_tauym
  namelist/hm_param_ocn/nfile_ocn_hm,timename_ocn_hm,varname_ocn_hm,Lcycle_ocn_hm,Tcycle_ocn_hm
  namelist/hm_io_ocn/fnames_ocn_hm
  namelist/um_param_ocn/nfile_ocn_um,timename_ocn_um,varname_ocn_um,Lcycle_ocn_um,Tcycle_ocn_um
  namelist/um_io_ocn/fnames_ocn_um
  namelist/vm_param_ocn/nfile_ocn_vm,timename_ocn_vm,varname_ocn_vm,Lcycle_ocn_vm,Tcycle_ocn_vm
  namelist/vm_io_ocn/fnames_ocn_vm
  namelist/wm_param_ocn/nfile_ocn_wm,timename_ocn_wm,varname_ocn_wm,Lcycle_ocn_wm,Tcycle_ocn_wm
  namelist/wm_io_ocn/fnames_ocn_wm
  namelist/Tzm_param_ocn/nfile_ocn_Tzm,timename_ocn_Tzm,varname_ocn_Tzm,Lcycle_ocn_Tzm,Tcycle_ocn_Tzm
  namelist/Tzm_io_ocn/fnames_ocn_Tzm
  read(5,date)
  read(5,io_ocn)
  read(5,param_ocn)
  read(5,taux_param_ocn)
  allocate(fnames_ocn_taux(nfile_ocn_taux))
  read(5,taux_io_ocn)
  read(5,tauy_param_ocn)
  allocate(fnames_ocn_tauy(nfile_ocn_tauy))
  read(5,tauy_io_ocn)
  read(5,sstm_param_ocn)
  allocate(fnames_ocn_sstm(nfile_ocn_sstm))
  read(5,sstm_io_ocn)
  read(5,tauxm_param_ocn)
  allocate(fnames_ocn_tauxm(nfile_ocn_tauxm))
  read(5,tauxm_io_ocn)
  read(5,tauym_param_ocn)
  allocate(fnames_ocn_tauym(nfile_ocn_tauym))
  read(5,tauym_io_ocn)
  read(5,hm_param_ocn)
  allocate(fnames_ocn_hm(nfile_ocn_hm))
  read(5,hm_io_ocn)
  read(5,um_param_ocn)
  allocate(fnames_ocn_um(nfile_ocn_um))
  read(5,um_io_ocn)
  read(5,vm_param_ocn)
  allocate(fnames_ocn_vm(nfile_ocn_vm))
  read(5,vm_io_ocn)
  read(5,wm_param_ocn)
  allocate(fnames_ocn_wm(nfile_ocn_wm))
  read(5,wm_io_ocn)
  read(5,Tzm_param_ocn)
  allocate(fnames_ocn_Tzm(nfile_ocn_Tzm))
  read(5,Tzm_io_ocn)

  ! Time setting
  call calendar_cal_length_ymdhms(start_yymmdd,start_hhmmss,end_yymmdd,end_hhmmss,1,tmp1)
  total_time=int(tmp1);ntime=int(tmp1/(dt*sec_to_day))
  write(*,*) "Start=",start_yymmdd,"End=",end_yymmdd
  write(*,*) "dt=",dt
  write(*,*) "Number of step=",ntime
  ! Read ocean grid
  call read_ocn_dyn_grd(fname_grd_ocn,ogrd)
  ! Set masking
  call set_mask_ocn(ogrd,oset%slip_ind)
  call read_ocn_sst_grd(fname_grd_ocn,ogrd)

  call initialize_ocn_dyn(ogrd)
  call initialize_ocn_sst(ogrd)
  call initialize_ocn_visc(ogrd,oset)
  if (flag_ini_ocn == "T") then
     call read_restart_ocn_dyn(fname_ini_ocn,ogrd)
     call read_restart_ocn_sst(fname_ini_ocn,ogrd)
  end if
  ! Read forcing fields
  call read_data_TLL_p(nfile_ocn_taux,fnames_ocn_taux,timename_ocn_taux,varname_ocn_taux,&
       & ogrd,ocn_taux_dta,start_yymmdd,start_hhmmss) ! Zonal wind
  ocn_taux_dta%Lcycle=Lcycle_ocn_taux;ocn_taux_dta%Tcycle=Tcycle_ocn_taux
  call read_data_TLL_p(nfile_ocn_tauy,fnames_ocn_tauy,timename_ocn_tauy,varname_ocn_tauy,&
       & ogrd,ocn_tauy_dta,start_yymmdd,start_hhmmss) ! Meridional wind
  ocn_tauy_dta%Lcycle=Lcycle_ocn_tauy;ocn_tauy_dta%Tcycle=Tcycle_ocn_tauy
  call read_data_TLL_p(nfile_ocn_sstm,fnames_ocn_sstm,timename_ocn_sstm,&
       & varname_ocn_sstm,ogrd,ocn_sstm_dta,start_yymmdd,start_hhmmss) ! Mean SST
  ocn_sstm_dta%Lcycle=Lcycle_ocn_sstm;ocn_sstm_dta%Tcycle=Tcycle_ocn_sstm
  call read_data_TLL_p(nfile_ocn_tauxm,fnames_ocn_tauxm,timename_ocn_tauxm,&
       & varname_ocn_tauxm,ogrd,ocn_tauxm_dta,start_yymmdd,start_hhmmss)
  ocn_tauxm_dta%Lcycle=Lcycle_ocn_tauxm;ocn_tauxm_dta%Tcycle=Tcycle_ocn_tauxm
  call read_data_TLL_p(nfile_ocn_tauym,fnames_ocn_tauym,timename_ocn_tauym,&
       & varname_ocn_tauym,ogrd,ocn_tauym_dta,start_yymmdd,start_hhmmss)
  ocn_tauym_dta%Lcycle=Lcycle_ocn_tauym;ocn_tauym_dta%Tcycle=Tcycle_ocn_tauym
  call read_data_TLL_u(nfile_ocn_um,fnames_ocn_um,timename_ocn_um,varname_ocn_um,ogrd,&
       & ocn_um_dta,start_yymmdd,start_hhmmss)
  ocn_um_dta%Lcycle=Lcycle_ocn_um;ocn_um_dta%Tcycle=Tcycle_ocn_um
  call read_data_TLL_v(nfile_ocn_vm,fnames_ocn_vm,timename_ocn_vm,varname_ocn_vm,ogrd,&
       & ocn_vm_dta,start_yymmdd,start_hhmmss)
  ocn_vm_dta%Lcycle=Lcycle_ocn_vm;ocn_vm_dta%Tcycle=Tcycle_ocn_vm
  call read_data_TLL_p(nfile_ocn_wm,fnames_ocn_wm,timename_ocn_wm,&
       & varname_ocn_wm,ogrd,ocn_wm_dta,start_yymmdd,start_hhmmss)
  ocn_wm_dta%Lcycle=Lcycle_ocn_wm;ocn_wm_dta%Tcycle=Tcycle_ocn_wm
  call read_data_TLL_p(nfile_ocn_hm,fnames_ocn_hm,timename_ocn_hm,&
       & varname_ocn_hm,ogrd,ocn_hm_dta,start_yymmdd,start_hhmmss)
  ocn_hm_dta%Lcycle=Lcycle_ocn_hm;ocn_hm_dta%Tcycle=Tcycle_ocn_hm
  call read_data_TLL_p(nfile_ocn_Tzm,fnames_ocn_Tzm,timename_ocn_Tzm,&
       & varname_ocn_Tzm,ogrd,ocn_Tzm_dta,start_yymmdd,start_hhmmss)
  ocn_Tzm_dta%Lcycle=Lcycle_ocn_Tzm;ocn_Tzm_dta%Tcycle=Tcycle_ocn_Tzm
  ! Prepare averaged file
  call create_avg_ocn_dyn_ZC(fname_avg_ocn,ogrd,oset,&
       & start_yymmdd,start_hhmmss,end_yymmdd,end_hhmmss,dt,&
       & out_avg_flag,out_avg_int,missing_value,istep_avg,ntime_avg)
  call initialize_ocn_dyn_avg(ogrd)

  call create_avg_ocn_sst_ZC(fname_avg_ocn,ogrd,oset,missing_value)
  call initialize_ocn_sst_avg(ogrd)
  ! Prepare diagnostic file
  call create_diag_ocn_ZC(fname_diag_ocn,ogrd,&
       & start_yymmdd,start_hhmmss,end_yymmdd,end_hhmmss,dt,&
       & out_diag_flag,out_diag_int,missing_value,istep_diag,ntime_diag)
  call initialize_ocn_diag(ogrd)
  ! Write output file
  iavg=1;iavg_count=0
  idiag=1;idiag_count=0
  do itime=1,ntime
     time_int=dt*itime*sec_to_day
     iavg_count = iavg_count + 1
     idiag_count = idiag_count + 1
     ! Obtain wind anomalies
     call get_data_TLL_p(time_int,start_yymmdd,start_hhmmss,ogrd,ocn_taux_dta)
     call get_data_TLL_p(time_int,start_yymmdd,start_hhmmss,ogrd,ocn_tauy_dta)
     ! Obtain mean SST fields
     call get_data_TLL_p(time_int,start_yymmdd,start_hhmmss,ogrd,ocn_sstm_dta)
     ! Obtain mean wind fields
     call get_data_TLL_p(time_int,start_yymmdd,start_hhmmss,ogrd,ocn_tauxm_dta)
     ! Obtain mean vw fields
     call get_data_TLL_p(time_int,start_yymmdd,start_hhmmss,ogrd,ocn_tauym_dta)
     ! Obtain mean ocean fields
     call get_data_TLL_u(time_int,start_yymmdd,start_hhmmss,ogrd,ocn_um_dta)
     call get_data_TLL_v(time_int,start_yymmdd,start_hhmmss,ogrd,ocn_vm_dta)
     call get_data_TLL_p(time_int,start_yymmdd,start_hhmmss,ogrd,ocn_wm_dta)
     call get_data_TLL_p(time_int,start_yymmdd,start_hhmmss,ogrd,ocn_hm_dta)
     call get_data_TLL_p(time_int,start_yymmdd,start_hhmmss,ogrd,ocn_Tzm_dta)

     ! Set wind stress
#if defined WSTRESS_BULK
     ! Read surface wind
     ogrd%ua_ocn%val=ocn_taux_dta%data_now%val
     ogrd%va_ocn%val=ocn_tauy_dta%data_now%val
     call ua_to_stress_anm(ogrd,oset,ocn_tauxm_dta,ocn_tauym_dta)
#else
     ! Read surface wind stress     #
     ogrd%taux_ocn%val=ocn_taux_dta%data_now%val-ocn_tauxm_dta%data_now%val
     ogrd%tauy_ocn%val=ocn_tauy_dta%data_now%val-ocn_tauym_dta%data_now%val
#endif
     ! Calculate Ekman currents
     call solve_ekman_ocn(ogrd,oset)
     ! Calculate geostrophic currents (with time stepping)
     call solve_rg_vgeo_ocn(ogrd,oset,dt)
     ! Calculate total current
     call solve_totalcurrent_ocn(ogrd,oset)
     ! Solve SST equation
     call solve_sst_ocn_ZC(ogrd,oset,ocn_um_dta,ocn_vm_dta,ocn_wm_dta,&
          & ocn_hm_dta,ocn_Tzm_dta,ocn_sstm_dta,dt)

     ! Add snapshot to mean fiedls
     call oper_avg_ocn_dyn(ogrd)
     call oper_avg_ocn_sst(ogrd)
     ! I/O operation
     if (itime .eq. istep_avg(iavg)) then
       call calendar_cal_ymdhms_after(start_yymmdd,start_hhmmss,dt*itime*sec_to_day,1,tmp_yymmdd,tmp_hhmmss)
        write(*,*) "Step (average) =",iavg," ",tmp_yymmdd,tmp_hhmmss,iavg_count
        call write_avg_ocn_dyn(fname_avg_ocn,ogrd,iavg_count,iavg)
        call write_avg_ocn_sst(fname_avg_ocn,ogrd,iavg_count,iavg)
        iavg_count=0
        iavg=iavg+1
        iavg=min(iavg,ntime_avg)
     end if
     if (itime .eq. istep_diag(idiag)) then
        write(*,*) "Step (Diagnose) =",idiag," ",itime,idiag_count
        call write_diag_ocn(fname_diag_ocn,ogrd,idiag_count,idiag)
        idiag_count=0
        idiag=idiag+1
        idiag=min(idiag,ntime_diag)
     end if
  end do
  ! Write restart file
  call write_restart_ocn_dyn(fname_rst_ocn,ogrd,missing_value)
  call write_restart_ocn_sst(fname_rst_ocn,ogrd,missing_value)
  call deallocate_ocn_all(ogrd)
  write(*,*) "Finish run, ended at ",end_yymmdd,end_hhmmss
  deallocate(fnames_ocn_taux)
  deallocate(fnames_ocn_tauy)
  deallocate(fnames_ocn_sstm)
  deallocate(fnames_ocn_tauxm)
  deallocate(fnames_ocn_tauym)
  deallocate(fnames_ocn_hm)
  deallocate(fnames_ocn_um)
  deallocate(fnames_ocn_vm)
  deallocate(fnames_ocn_wm)
  deallocate(fnames_ocn_Tzm)
end program zc_ogcm_full
