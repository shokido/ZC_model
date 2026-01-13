program zc_cgcm_full
  use run_params
  use run_types
  use calendar_sub
  use mod_io_master
  use mod_io_avg
  use mod_io_diag
  use mod_wstress
  use mod_ocn_solver_zc
  use mod_ocn_dta
  use mod_atm_dta
  use mod_atm_solver_gill
  use mod_coupler_dta
  
  implicit none
  integer ntime_couple
  type(couple_dta) :: cgrd
  type(ocn_dta) :: ogrd
  type(ocn_set) :: oset
  character(len=maxlen) :: fname_grd_ocn
  character(len=maxlen) :: flag_ini_ocn,fname_ini_ocn
  character(len=maxlen) :: fname_avg_ocn
  character(len=maxlen) :: fname_diag_ocn
  character(len=maxlen) :: fname_rst_ocn
  type(atm_dta) :: agrd
  type(atm_set) :: aset
  character(len=maxlen) :: fname_grd_atm
  character(len=maxlen) :: fname_avg_atm
  integer:: iy,ix
  integer :: itime,ntime,total_time
  real(idx) :: dt,tmp1,time_int,time_couple
  integer :: start_yymmdd,start_hhmmss
  integer :: end_yymmdd,end_hhmmss
  ! Wind stress fields
  type(TLL_dta) :: ocn_sstm_dta
  integer :: nfile_ocn_sstm
  character(len=maxlen),allocatable :: fnames_ocn_sstm(:)
  character(len=maxlen) :: timename_ocn_sstm,varname_ocn_sstm
  character(1) :: Lcycle_ocn_sstm
  real(idx) :: Tcycle_ocn_sstm
  ! Mean wind
  type(TLL_dta) :: ocn_tauxm_dta
  integer :: nfile_ocn_tauxm
  character(len=maxlen),allocatable :: fnames_ocn_tauxm(:)
  character(len=maxlen) :: timename_ocn_tauxm,varname_ocn_tauxm
  character(1) :: Lcycle_ocn_tauxm
  real(idx) :: Tcycle_ocn_tauxm
  ! Mean tauy
  type(TLL_dta) :: ocn_tauym_dta
  integer :: nfile_ocn_tauym
  character(len=maxlen),allocatable :: fnames_ocn_tauym(:)
  character(len=maxlen) :: timename_ocn_tauym,varname_ocn_tauym
  character(1) :: Lcycle_ocn_tauym
  real(idx) :: Tcycle_ocn_tauym

  ! Mean oceanic fields
  type(TLL_dta) :: ocn_hm_dta
  integer :: nfile_ocn_hm
  character(len=maxlen),allocatable :: fnames_ocn_hm(:)
  character(len=maxlen) :: timename_ocn_hm,varname_ocn_hm
  character(1) :: Lcycle_ocn_hm
  real(idx) :: Tcycle_ocn_hm

  type(TLL_dta) :: ocn_um_dta
  integer :: nfile_ocn_um
  character(len=maxlen),allocatable :: fnames_ocn_um(:)
  character(len=maxlen) :: timename_ocn_um,varname_ocn_um
  character(1) :: Lcycle_ocn_um
  real(idx) :: Tcycle_ocn_um
  type(TLL_dta) :: ocn_vm_dta
  integer :: nfile_ocn_vm
  character(len=maxlen),allocatable :: fnames_ocn_vm(:)
  character(len=maxlen) :: timename_ocn_vm,varname_ocn_vm
  character(1) :: Lcycle_ocn_vm
  real(idx) :: Tcycle_ocn_vm
  type(TLL_dta) :: ocn_wm_dta
  integer :: nfile_ocn_wm
  character(len=maxlen),allocatable :: fnames_ocn_wm(:)
  character(len=maxlen) :: timename_ocn_wm,varname_ocn_wm
  character(1) :: Lcycle_ocn_wm
  real(idx) :: Tcycle_ocn_wm
  
  type(TLL_dta) :: ocn_Tzm_dta
  integer :: nfile_ocn_Tzm
  character(len=maxlen),allocatable :: fnames_ocn_Tzm(:)
  character(len=maxlen) :: timename_ocn_Tzm,varname_ocn_Tzm
  character(1) :: Lcycle_ocn_Tzm
  real(idx) :: Tcycle_ocn_Tzm
  ! SST forcing
  type(TLL_dta) :: atm_sstm_dta
  integer :: nfile_atm_sstm
  character(len=maxlen),allocatable :: fnames_atm_sstm(:)
  character(len=maxlen) :: timename_atm_sstm,varname_atm_sstm
  character(1) :: Lcycle_atm_sstm
  real(idx) :: Tcycle_atm_sstm
#if defined CONV
  type(TLL_dta) :: atm_uam_dta
  integer :: nfile_atm_uam
  character(len=maxlen),allocatable :: fnames_atm_uam(:)
  character(len=maxlen) :: timename_atm_uam,varname_atm_uam
  character(1) :: Lcycle_atm_uam
  real(idx) :: Tcycle_atm_uam
  type(TLL_dta) :: atm_vam_dta
  integer :: nfile_atm_vam
  character(len=maxlen),allocatable :: fnames_atm_vam(:)
  character(len=maxlen) :: timename_atm_vam,varname_atm_vam
  character(1) :: Lcycle_atm_vam
  real(idx) :: Tcycle_atm_vam
#endif
  integer :: tmp_yymmdd,tmp_hhmmss

  namelist/date/dt,start_yymmdd,start_hhmmss,end_yymmdd,end_hhmmss,time_couple
  namelist/io_ocn/fname_grd_ocn
  namelist/io_ocn/flag_ini_ocn,fname_ini_ocn
  namelist/io_ocn/fname_avg_ocn,out_avg_flag,out_avg_int
  namelist/io_ocn/fname_diag_ocn,out_diag_flag,out_diag_int
  namelist/io_ocn/fname_rst_ocn
  namelist/param_ocn/oset
  namelist/io_atm/fname_grd_atm
  namelist/io_atm/fname_avg_atm,out_avg_flag,out_avg_int
  namelist/param_atm/aset
  
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
  namelist/sstm_param_atm/nfile_atm_sstm,timename_atm_sstm,varname_atm_sstm,Lcycle_atm_sstm,Tcycle_atm_sstm
  namelist/sstm_io_atm/fnames_atm_sstm
#if defined CONV
  namelist/uam_param_atm/nfile_atm_uam,timename_atm_uam,varname_atm_uam,Lcycle_atm_uam,Tcycle_atm_uam
  namelist/uam_io_atm/fnames_atm_uam
  namelist/vam_param_atm/nfile_atm_vam,timename_atm_vam,varname_atm_vam,Lcycle_atm_vam,Tcycle_atm_vam
  namelist/vam_io_atm/fnames_atm_vam
#endif
  read(5,date)
  read(5,io_ocn)
  read(5,io_atm)
  read(5,param_ocn)
  read(5,param_atm)
  ntime_couple=int(time_couple/dt)
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
  read(5,sstm_param_atm)
  allocate(fnames_atm_sstm(nfile_atm_sstm))
  read(5,sstm_io_atm)
#if defined CONV
  read(5,uam_param_atm)
  allocate(fnames_atm_uam(nfile_atm_uam))
  read(5,uam_io_atm)
  read(5,vam_param_atm)
  allocate(fnames_atm_vam(nfile_atm_vam))
  read(5,vam_io_atm)
#endif
  ! Parameter
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

  ! Atmospheric grid
  call read_atm_grd(fname_grd_atm,agrd)
  call set_coord_atm(agrd,aset)
 
  ! Set coupler
  call initialize_coupler(cgrd,ogrd,agrd)
  
  ! Read forcing field
  call read_data_TLL_p(nfile_ocn_sstm,fnames_ocn_sstm,timename_ocn_sstm,&
       & varname_ocn_sstm,ogrd,ocn_sstm_dta,start_yymmdd,start_hhmmss)
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
  ! Read forcing field
  call read_data_TLL_atm(nfile_atm_sstm,fnames_atm_sstm,timename_atm_sstm,&
       & varname_atm_sstm,agrd,atm_sstm_dta,start_yymmdd,start_hhmmss)
  atm_sstm_dta%Lcycle=Lcycle_atm_sstm;atm_sstm_dta%Tcycle=Tcycle_atm_sstm
#if defined CONV
  call read_data_TLL_atm(nfile_atm_uam,fnames_atm_uam,timename_atm_uam,&
       & varname_atm_uam,agrd,atm_uam_dta,start_yymmdd,start_hhmmss)
  atm_uam_dta%Lcycle=Lcycle_atm_uam;atm_uam_dta%Tcycle=Tcycle_atm_uam
  call read_data_TLL_atm(nfile_atm_vam,fnames_atm_vam,timename_atm_vam,&
       & varname_atm_vam,agrd,atm_vam_dta,start_yymmdd,start_hhmmss)
  atm_vam_dta%Lcycle=Lcycle_atm_vam;atm_vam_dta%Tcycle=Tcycle_atm_vam
#endif  
  ! Make averaged file
  call create_avg_ocn_dyn_ZC(fname_avg_ocn,ogrd,oset,&
       & start_yymmdd,start_hhmmss,end_yymmdd,end_hhmmss,dt,&
       & out_avg_flag,out_avg_int,missing_value,istep_avg,ntime_avg)
  call create_avg_ocn_sst_ZC(fname_avg_ocn,ogrd,oset,missing_value)
  call initialize_ocn_dyn_avg(ogrd)
  call initialize_ocn_sst_avg(ogrd)
  ! Make diagnostic file
  call create_diag_ocn_ZC(fname_diag_ocn,ogrd,&
       & start_yymmdd,start_hhmmss,end_yymmdd,end_hhmmss,dt,&
       & out_diag_flag,out_diag_int,missing_value,istep_diag,ntime_diag)
  call initialize_ocn_diag(ogrd)
  ! Make averaged file
  call create_avg_atm_ZC(fname_avg_atm,agrd,&
       & start_yymmdd,start_hhmmss,end_yymmdd,end_hhmmss,dt,&
       & out_avg_flag,out_avg_int,missing_value,istep_avg,ntime_avg)
  ! Write output file
  iavg=1;iavg_count=0
  idiag=1;idiag_count=0
  do itime=1,ntime
     time_int=dt*itime*sec_to_day
     iavg_count = iavg_count + 1
     idiag_count = idiag_count + 1
     ! Obtain mean SST fields
     call get_data_TLL_p(time_int,ogrd,ocn_sstm_dta)
     ! Obtain mean uw fields
     call get_data_TLL_p(time_int,ogrd,ocn_tauxm_dta)
     ! Obtain mean vw fields
     call get_data_TLL_p(time_int,ogrd,ocn_tauym_dta)
     ! Obtain mean ocean fields
     call get_data_TLL_u(time_int,ogrd,ocn_um_dta)
     call get_data_TLL_v(time_int,ogrd,ocn_vm_dta)
     call get_data_TLL_p(time_int,ogrd,ocn_wm_dta)
     call get_data_TLL_p(time_int,ogrd,ocn_hm_dta)
     call get_data_TLL_p(time_int,ogrd,ocn_Tzm_dta)
     if (mod(itime,ntime_couple)==0) then
       ! Get climatological SST
        call get_data_TLL_atm(time_int,agrd,atm_sstm_dta)
        agrd%sstm_atm%val=atm_sstm_dta%data_now%val
#if defined CONV
        call get_data_TLL_atm(time_int,agrd,atm_uam_dta)
        call get_data_TLL_atm(time_int,agrd,atm_vam_dta)
#endif
        ! Send SSTA to ATM.grid
        call exchange_OtoA(cgrd,ogrd,agrd)
#if defined CONV
        call return_uvp_fromSST_conv(agrd,aset,atm_uam_dta,atm_vam_dta)
#else
        call return_uvp_fromSST_ZC87(agrd,aset)
!        call return_uvp_fromSST_GJ20(agrd,aset)
#endif
        ! Exchange information
        call exchange_AtoO(cgrd,ogrd,agrd)
     end if
     ! Convert surface wind to wind stress
#if defined INI
     if (time_int <= 30*4) then
        do iy=0,ogrd%ny_p+1
           do ix=0,ogrd%nx_p+1
               if (ogrd%lon_p%val(ix)>=145.0_idx .and. &
               & ogrd%lon_p%val(ix)<=190.0_idx) then
                   ogrd%ua_ocn%val(ix,iy)=2.0_idx*&
                   & exp(-1.0_idx*((ogrd%lat_p%val(iy)-0.0)/20)**2)                    
               end if
!              ogrd%ua_ocn%val(ix,iy)=2*&
!                   & exp(-1.0_idx*((ogrd%lon_p%val(ix)-167.5)/22.5)**2)*&
!                   & exp(-1.0_idx*((ogrd%lat_p%val(iy)-0.0)/20)**2) !&
!              ogrd%va_ocn%val(ix,iy)=0.0_idx
           end do
        end do
     end if
#endif
     call ua_to_stress_anm(ogrd,oset,ocn_tauxm_dta,ocn_tauym_dta)
     ! Calculate Ekman current
     call solve_ekman_ocn(ogrd,oset)
     ! Calculate geostrophic current
     call solve_rg_vgeo_ocn(ogrd,oset,dt)
!     call solve_rg_vgeo_ocn3(ogrd,oset,dt)
     ! Calculate total current
     call solve_totalcurrent_ocn(ogrd,oset)
     call solve_sst_ocn_ZC(ogrd,oset,ocn_um_dta,ocn_vm_dta,ocn_wm_dta,&
          & ocn_hm_dta,ocn_Tzm_dta,ocn_sstm_dta,dt)
!     call limit_total_SST(ogrd,ocn_sstm_dta)
     call oper_avg_ocn_dyn(ogrd)
     call oper_avg_ocn_sst(ogrd)
     call oper_avg_atm(agrd)
     if (itime .eq. istep_avg(iavg)) then
        call calendar_cal_ymdhms_after(start_yymmdd,start_hhmmss,dt*itime*sec_to_day,1,tmp_yymmdd,tmp_hhmmss)
        write(*,*) "Step (average) =",iavg," ",tmp_yymmdd,tmp_hhmmss,iavg_count
        call write_avg_ocn_dyn(fname_avg_ocn,ogrd,iavg_count,iavg)
        call write_avg_ocn_sst(fname_avg_ocn,ogrd,iavg_count,iavg)
        call write_avg_atm(fname_avg_atm,agrd,iavg_count,iavg)
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
  call write_restart_ocn_dyn(fname_rst_ocn,ogrd,missing_value)
  call write_restart_ocn_sst(fname_rst_ocn,ogrd,missing_value)
  call deallocate_ocn_all(ogrd)
  call deallocate_atm(agrd)
  call deallocate_coupler(cgrd)
  write(*,*) "Finish run, ended at ",end_yymmdd,end_hhmmss
  deallocate(fnames_ocn_sstm)
  deallocate(fnames_ocn_tauxm)
  deallocate(fnames_ocn_tauym)
  deallocate(fnames_ocn_hm)
  deallocate(fnames_ocn_um)
  deallocate(fnames_ocn_vm)
  deallocate(fnames_ocn_wm)
  deallocate(fnames_ocn_Tzm)
end program zc_cgcm_full
