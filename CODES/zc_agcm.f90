program zc_agcm
  use run_params
  use run_types
  use calendar_sub
  use mod_io_master
  use mod_io_avg
  use mod_atm_dta
  use mod_atm_solver_gill
  implicit none
  type(atm_dta) :: agrd
  type(atm_set) :: aset
  character(len=maxlen) :: fname_grd_atm
  character(len=maxlen) :: fname_avg_atm
  integer :: itime,ntime,total_time
  real(idx) :: dt,tmp1,time_int
  integer :: start_yymmdd,start_hhmmss
  integer :: end_yymmdd,end_hhmmss
  integer :: ix,iy
  ! Climatological SST
  type(TLL_dta) :: atm_sstm_dta
  integer :: nfile_atm_sstm
  character(len=maxlen),allocatable :: fnames_atm_sstm(:)
  character(len=maxlen) :: timename_atm_sstm,varname_atm_sstm
  character(1) :: Lcycle_atm_sstm
  real(idx) :: Tcycle_atm_sstm
  ! SST anomaly
  type(TLL_dta) :: atm_ssta_dta
  integer :: nfile_atm_ssta
  character(len=maxlen),allocatable :: fnames_atm_ssta(:)
  character(len=maxlen) :: timename_atm_ssta,varname_atm_ssta
  character(1) :: Lcycle_atm_ssta
  real(idx) :: Tcycle_atm_ssta
  ! Climatological zonal wind
  type(TLL_dta) :: atm_uam_dta
  integer :: nfile_atm_uam
  character(len=maxlen),allocatable :: fnames_atm_uam(:)
  character(len=maxlen) :: timename_atm_uam,varname_atm_uam
  character(1) :: Lcycle_atm_uam
  real(idx) :: Tcycle_atm_uam
  ! Climatological zonal wind
  type(TLL_dta) :: atm_vam_dta
  integer :: nfile_atm_vam
  character(len=maxlen),allocatable :: fnames_atm_vam(:)
  character(len=maxlen) :: timename_atm_vam,varname_atm_vam
  character(1) :: Lcycle_atm_vam
  real(idx) :: Tcycle_atm_vam
  namelist/date/dt,start_yymmdd,start_hhmmss,end_yymmdd,end_hhmmss
  namelist/io_atm/fname_grd_atm
  namelist/io_atm/fname_avg_atm,out_avg_flag,out_avg_int
  namelist/param_atm/aset
  namelist/sstm_param_atm/nfile_atm_sstm,timename_atm_sstm,varname_atm_sstm,Lcycle_atm_sstm,Tcycle_atm_sstm
  namelist/sstm_io_atm/fnames_atm_sstm
  namelist/ssta_param_atm/nfile_atm_ssta,timename_atm_ssta,varname_atm_ssta,Lcycle_atm_ssta,Tcycle_atm_ssta
  namelist/ssta_io_atm/fnames_atm_ssta
  namelist/uam_param_atm/nfile_atm_uam,timename_atm_uam,varname_atm_uam,Lcycle_atm_uam,Tcycle_atm_uam
  namelist/uam_io_atm/fnames_atm_uam
  namelist/vam_param_atm/nfile_atm_vam,timename_atm_vam,varname_atm_vam,Lcycle_atm_vam,Tcycle_atm_vam
  namelist/vam_io_atm/fnames_atm_vam
  read(5,date)
  read(5,io_atm)
  read(5,param_atm)
  read(5,sstm_param_atm)
  allocate(fnames_atm_sstm(nfile_atm_sstm))
  read(5,sstm_io_atm)
  read(5,ssta_param_atm)
  allocate(fnames_atm_ssta(nfile_atm_ssta))
  read(5,ssta_io_atm)
  if (aset%heating_type=="ZC87_conv") then
    read(5,uam_param_atm)
    allocate(fnames_atm_uam(nfile_atm_uam))
    read(5,uam_io_atm)
    read(5,vam_param_atm)
    allocate(fnames_atm_vam(nfile_atm_vam))
    read(5,vam_io_atm)
  endif
  ! Time setting
  call calendar_cal_length_ymdhms(start_yymmdd,start_hhmmss,end_yymmdd,end_hhmmss,1,tmp1)
  total_time=int(tmp1);ntime=int(tmp1/(dt*sec_to_day))

  ! Read atm grid
  call read_atm_grd(fname_grd_atm,agrd)
  ! Set masking
  call set_coord_atm(agrd,aset)
  ! Read forcing field
  call read_data_TLL_atm(nfile_atm_sstm,fnames_atm_sstm,timename_atm_sstm,&
       & varname_atm_sstm,agrd,atm_sstm_dta,start_yymmdd,start_hhmmss)
  atm_sstm_dta%Lcycle=Lcycle_atm_sstm;atm_sstm_dta%Tcycle=Tcycle_atm_sstm
  call read_data_TLL_atm(nfile_atm_ssta,fnames_atm_ssta,timename_atm_ssta,&
       & varname_atm_ssta,agrd,atm_ssta_dta,start_yymmdd,start_hhmmss)
  atm_ssta_dta%Lcycle=Lcycle_atm_ssta;atm_ssta_dta%Tcycle=Tcycle_atm_ssta
if (aset%heating_type=="ZC87_conv") then
  call read_data_TLL_atm(nfile_atm_uam,fnames_atm_uam,timename_atm_uam,&
       & varname_atm_uam,agrd,atm_uam_dta,start_yymmdd,start_hhmmss)
  atm_uam_dta%Lcycle=Lcycle_atm_uam;atm_uam_dta%Tcycle=Tcycle_atm_uam
  call read_data_TLL_atm(nfile_atm_vam,fnames_atm_vam,timename_atm_vam,&
       & varname_atm_vam,agrd,atm_vam_dta,start_yymmdd,start_hhmmss)
  atm_vam_dta%Lcycle=Lcycle_atm_vam;atm_vam_dta%Tcycle=Tcycle_atm_vam
end if
  ! Make averaged file
  call create_avg_atm_ZC(fname_avg_atm,agrd,&
       & start_yymmdd,start_hhmmss,end_yymmdd,end_hhmmss,dt,&
       & out_avg_flag,out_avg_int,missing_value,istep_avg,ntime_avg)

  ! Write output file
  iavg=1;iavg_count=0
  do itime=1,ntime
     time_int=dt*itime*sec_to_day
     iavg_count = iavg_count + 1
     ! Get climatological SST
     call get_data_TLL_atm(time_int,agrd,atm_ssta_dta)
     ! Get SST anomaly
     call get_data_TLL_atm(time_int,agrd,atm_sstm_dta)
if (aset%heating_type=="ZC87_conv") then
     ! Get climatological Wind
     call get_data_TLL_atm(time_int,agrd,atm_uam_dta)
     call get_data_TLL_atm(time_int,agrd,atm_vam_dta)
end if
     agrd%ssta_atm%val=atm_ssta_dta%data_now%val
     agrd%sstm_atm%val=atm_sstm_dta%data_now%val
if (aset%heating_type=="ZC87") then
     call return_uvp_fromSST_ZC87(agrd,aset)
else if (aset%heating_type=="ZC87_conv") then
     call return_uvp_fromSST_ZC87_conv(agrd,aset,atm_uam_dta,atm_vam_dta)
else if (aset%heating_type=="GJ22") then
     call return_uvp_fromSST_GJ22(agrd,aset)
end if
!     call return_uvp_fromSST_GJ20(agrd,aset)
     call oper_avg_atm(agrd)
     if (itime .eq. istep_avg(iavg)) then
        write(*,*) "Step (average) =",iavg," ",tmp_yymmdd,tmp_hhmmss,iavg_count
        call write_avg_atm(fname_avg_atm,agrd,iavg_count,iavg)
        iavg_count=0
        iavg=iavg+1
        iavg=min(iavg,ntime_avg)
     end if
  end do
  deallocate(fnames_atm_ssta)
  deallocate(fnames_atm_sstm)
if (aset%heating_type=="ZC87_conv") then
  deallocate(fnames_atm_uam)
  deallocate(fnames_atm_vam)
end if
end program zc_agcm
