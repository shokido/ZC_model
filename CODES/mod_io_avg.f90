module mod_io_avg
  use run_params
  use run_types
  use calendar_sub
  use ncdf_write
  implicit none
  ! Time array
  integer :: ntime_avg,iavg,iavg_count
  real(idx),allocatable :: time_avg(:)
  integer,allocatable :: istep_avg(:)  
  integer :: out_avg_flag,out_avg_int
  character(maxlen) :: fname_out_avg
  ! Namelist
  namelist/output_avg/out_avg_flag,out_avg_int
  namelist/output_avg/fname_out_avg
contains
  !==========================!
  ! Creation of average file !
  !==========================!
  subroutine initialize_ocn_dyn_avg(grd)
    implicit none
    type(ocn_dta),intent(inout) :: grd
    integer :: nx_p,ny_p
    nx_p=grd%nx_p;ny_p=grd%ny_p
    allocate(grd%h_sw_avg%val(0:nx_p+1,0:ny_p+1));grd%h_sw_avg%val(0:nx_p+1,0:ny_p+1)=0.0_idx
    allocate(grd%u_sw_avg%val(1:nx_p+1,0:ny_p+1));grd%u_sw_avg%val(1:nx_p+1,0:ny_p+1)=0.0_idx
    allocate(grd%v_sw_avg%val(0:nx_p+1,1:ny_p+1));grd%v_sw_avg%val(0:nx_p+1,1:ny_p+1)=0.0_idx
    allocate(grd%w_ek_avg%val(0:nx_p+1,0:ny_p+1));grd%w_ek_avg%val(0:nx_p+1,0:ny_p+1)=0.0_idx
    allocate(grd%u_ek_avg%val(1:nx_p+1,0:ny_p+1));grd%u_ek_avg%val(1:nx_p+1,0:ny_p+1)=0.0_idx
    allocate(grd%v_ek_avg%val(0:nx_p+1,1:ny_p+1));grd%v_ek_avg%val(0:nx_p+1,1:ny_p+1)=0.0_idx
    allocate(grd%u_ocn_1_avg%val(1:nx_p+1,0:ny_p+1));grd%u_ocn_1_avg%val(1:nx_p+1,0:ny_p+1)=0.0_idx
    allocate(grd%v_ocn_1_avg%val(0:nx_p+1,1:ny_p+1));grd%v_ocn_1_avg%val(0:nx_p+1,1:ny_p+1)=0.0_idx
    allocate(grd%w_ocn_1_avg%val(0:nx_p+1,0:ny_p+1));grd%w_ocn_1_avg%val(0:nx_p+1,0:ny_p+1)=0.0_idx
    ! Average array
    allocate(grd%taux_ocn_avg%val(0:nx_p+1,0:ny_p+1));grd%taux_ocn_avg%val(0:nx_p+1,0:ny_p+1)=0.0_idx
    allocate(grd%tauy_ocn_avg%val(0:nx_p+1,0:ny_p+1));grd%tauy_ocn_avg%val(0:nx_p+1,0:ny_p+1)=0.0_idx    
  end subroutine initialize_ocn_dyn_avg
  subroutine initialize_ocn_sst_avg(grd)
    implicit none
    type(ocn_dta),intent(inout) :: grd
    integer :: nx_p,ny_p
    nx_p=grd%nx_p;ny_p=grd%ny_p
    allocate(grd%ssta_ocn_avg%val(0:nx_p+1,0:ny_p+1));grd%ssta_ocn_avg%val(0:nx_p+1,0:ny_p+1)=0.0_idx
  end subroutine initialize_ocn_sst_avg
  subroutine create_avg_atm_ZC(fname,agrd,&
       & start_yymmdd,start_hhmmss,end_yymmdd,end_hhmmss,dt,&
       & out_flag,out_int,missing_value,istep,nt)
    implicit none
    character(len=*),intent(in) :: fname
    type(atm_dta),intent(in) :: agrd
    real(idx),intent(in):: dt
    integer,intent(in) :: start_yymmdd,start_hhmmss,end_yymmdd,end_hhmmss,out_flag,out_int
    real(idx),intent(in) :: missing_value
    integer,allocatable,intent(inout) :: istep(:)
    integer,intent(inout) :: nt
    real(idx),allocatable :: time(:)
    integer :: it
    character(len=maxlen) :: ref_time
    integer :: tmp_yymmdd,tmp_hhmmss
    character(len=maxlen) :: dim_names(3),dim_types(3)
    real(idx) :: tmp
    integer :: nx,ny
    ref_time=calendar_create_time_att(start_yymmdd,start_hhmmss,1)
    call calendar_cal_length_ymdhms(start_yymmdd,start_hhmmss,end_yymmdd,end_hhmmss,out_flag,tmp)
    nt= int(tmp / out_int)
    nt=max(nt,1)
    if (allocated(istep) .eqv. .true.) then
       deallocate(istep)
    end if    
    if (allocated(time) .eqv. .true.) then
       deallocate(time)
    end if    
    allocate(time(nt)); allocate(istep(nt))
    do it=1,nt
       time(it) = (real(it))* out_int
       call calendar_cal_ymdhms_after(start_yymmdd,start_hhmmss,time(it),out_flag,tmp_yymmdd,tmp_hhmmss)
       call calendar_cal_length_ymdhms(start_yymmdd,start_hhmmss,tmp_yymmdd,tmp_hhmmss,1,tmp)
       istep(it)=int(tmp/(dt*sec_to_day))
       time(it) = (real(it-0.5))* out_int
       call calendar_cal_ymdhms_after(start_yymmdd,start_hhmmss,time(it),out_flag,tmp_yymmdd,tmp_hhmmss)
       call calendar_cal_length_ymdhms(start_yymmdd,start_hhmmss,tmp_yymmdd,tmp_hhmmss,1,tmp)
       time(it) = real(tmp)
    end do
    dim_names(1)="lon"
    dim_names(2)="lat"
    dim_names(3)="time"
    do it=1,3
       dim_types(it)="double"
    end do
    nx=agrd%nx_atm
    ny=agrd%ny_atm
    call writenet_def_dim(trim(fname),3,(/nx,ny,nt/),dim_names,dim_types)
    call add_var_att(trim(fname),"lon","long_name","longitude")
    call add_var_att(trim(fname),"lon","units","degrees_east")
    call add_var_att(trim(fname),"lat","long_name","latitude")
    call add_var_att(trim(fname),"lat","units","degrees_north")
    call add_var_att(trim(fname),"time","long_name","time")
    call add_var_att(trim(fname),"time","units",ref_time)

    ! Output p field
    dim_names(1)="lon";dim_names(2)="lat";dim_names(3)="time"
    call writenet_def_var(trim(fname),3,1,dim_names(1:3),(/"q_atm"/),(/"double"/))
    call add_var_att(trim(fname),"q_atm","long_name","pressure")
    call add_var_att(trim(fname),"q_atm","units","m^2/s^2")
    call add_var_att(trim(fname),"q_atm","missing_value",missing_value)

    call writenet_def_var(trim(fname),3,1,dim_names(1:3),(/"p_atm"/),(/"double"/))
    call add_var_att(trim(fname),"p_atm","long_name","pressure")
    call add_var_att(trim(fname),"p_atm","units","m^2/s^2")
    call add_var_att(trim(fname),"p_atm","missing_value",missing_value)

    call writenet_def_var(trim(fname),3,1,dim_names(1:3),(/"u_atm"/),(/"double"/))
    call add_var_att(trim(fname),"u_atm","long_name","velocity")
    call add_var_att(trim(fname),"u_atm","units","m^/s")
    call add_var_att(trim(fname),"u_atm","missing_value",missing_value)
    call writenet_def_var(trim(fname),3,1,dim_names(1:3),(/"v_atm"/),(/"double"/))
    call add_var_att(trim(fname),"v_atm","long_name","velocity")
    call add_var_att(trim(fname),"v_atm","units","m^/s")
    call add_var_att(trim(fname),"v_atm","missing_value",missing_value)
    call writenet_def_var(trim(fname),3,1,dim_names(1:3),(/"ssta_atm"/),(/"double"/))
    call add_var_att(trim(fname),"ssta_atm","long_name","SST anomaly")
    call add_var_att(trim(fname),"ssta_atm","units","degrees_celsius")
    call add_var_att(trim(fname),"ssta_atm","missing_value",missing_value)
    call writenet_def_var(trim(fname),3,1,dim_names(1:3),(/"sstm_atm"/),(/"double"/))
    call add_var_att(trim(fname),"sstm_atm","long_name","SST clim")
    call add_var_att(trim(fname),"sstm_atm","units","degrees_celsius")
    call add_var_att(trim(fname),"sstm_atm","missing_value",missing_value)
    ! Write coordinate variables
    call writenet_wv(trim(fname),"lon",(/1/),(/nx/),agrd%lon_atm%val)
    call writenet_wv(trim(fname),"lat",(/1/),(/ny/),agrd%lat_atm%val)
    call writenet_wv(trim(fname),"time",(/1/),(/nt/),time(1:nt))    
  end subroutine create_avg_atm_ZC
  ! Ocean
  subroutine create_avg_ocn_dyn_ZC(fname,ogrd,oset,&
       & start_yymmdd,start_hhmmss,end_yymmdd,end_hhmmss,dt,&
       & out_flag,out_int,missing_value,istep,nt)
    implicit none
    character(len=*),intent(in) :: fname
    type(ocn_dta),intent(in) :: ogrd
    type(ocn_set),intent(in) :: oset
    real(idx),intent(in):: dt
    integer,intent(in) :: start_yymmdd,start_hhmmss,end_yymmdd,end_hhmmss,out_flag,out_int
    real(idx),intent(in) :: missing_value
    integer,allocatable,intent(inout) :: istep(:)
    integer,intent(inout) :: nt
    real(idx),allocatable :: time(:)
    integer :: it
    character(len=maxlen) :: ref_time
    integer :: tmp_yymmdd,tmp_hhmmss
    character(len=maxlen) :: dim_names(7),dim_types(7)
    real(idx) :: tmp
    ref_time=calendar_create_time_att(start_yymmdd,start_hhmmss,1)
    call calendar_cal_length_ymdhms(start_yymmdd,start_hhmmss,end_yymmdd,end_hhmmss,out_flag,tmp)
    nt= int(tmp / out_int)
    nt=max(nt,1)
    if (allocated(istep) .eqv. .true.) then
       deallocate(istep)
    end if    
    if (allocated(time) .eqv. .true.) then
       deallocate(time)
    end if    
    allocate(time(nt)); allocate(istep(nt))
    do it=1,nt
       time(it) = (real(it))* out_int
       call calendar_cal_ymdhms_after(start_yymmdd,start_hhmmss,time(it),out_flag,tmp_yymmdd,tmp_hhmmss)
       call calendar_cal_length_ymdhms(start_yymmdd,start_hhmmss,tmp_yymmdd,tmp_hhmmss,1,tmp)
       istep(it)=int(tmp/(dt*sec_to_day))
       time(it) = (real(it-0.5))* out_int
       call calendar_cal_ymdhms_after(start_yymmdd,start_hhmmss,time(it),out_flag,tmp_yymmdd,tmp_hhmmss)
       call calendar_cal_length_ymdhms(start_yymmdd,start_hhmmss,tmp_yymmdd,tmp_hhmmss,1,tmp)
       time(it) = real(tmp)
    end do
    dim_names(1)="x_p"
    dim_names(2)="y_p"
    dim_names(3)="x_u"
    dim_names(4)="y_u"
    dim_names(5)="x_v"
    dim_names(6)="y_v"
    dim_names(7)="time"
    do it=1,7
       dim_types(it)="double"
    end do
    call writenet_def_dim(trim(fname),7,(/ogrd%nx_p+2,ogrd%ny_p+2,&
         & ogrd%nx_p+1,ogrd%ny_p+2,ogrd%nx_p+2,ogrd%ny_p+1&
         & ,nt/),dim_names,dim_types)

    call writenet_def_var(trim(fname),1,1,(/"x_p"/),(/"lon_p"/),(/"double"/))
    call add_var_att(trim(fname),"x_p","long_name","longitude on p-grid")
    call add_var_att(trim(fname),"x_p","units","degrees_east")
    call writenet_def_var(trim(fname),1,1,(/"y_p"/),(/"lat_p"/),(/"double"/))
    call add_var_att(trim(fname),"y_p","long_name","latitude on p-grid")
    call add_var_att(trim(fname),"y_p","units","degrees_north")
    call writenet_def_var(trim(fname),1,1,(/"x_u"/),(/"lon_u"/),(/"double"/))
    call add_var_att(trim(fname),"x_u","long_name","longitude on u-grid")
    call add_var_att(trim(fname),"x_u","units","degrees_east")
    call writenet_def_var(trim(fname),1,1,(/"y_u"/),(/"lat_u"/),(/"double"/))
    call add_var_att(trim(fname),"y_u","long_name","latitude on u-grid")
    call add_var_att(trim(fname),"y_u","units","degrees_north")
    call writenet_def_var(trim(fname),1,1,(/"x_v"/),(/"lon_v"/),(/"double"/))
    call add_var_att(trim(fname),"x_v","long_name","longitude on v-grid")
    call add_var_att(trim(fname),"x_v","units","degrees_east")
    call writenet_def_var(trim(fname),1,1,(/"y_v"/),(/"lat_v"/),(/"double"/))
    call add_var_att(trim(fname),"y_v","long_name","latitude on v-grid")
    call add_var_att(trim(fname),"y_v","units","degrees_north")
    call writenet_def_var(trim(fname),2,2,(/"x_p","y_p"/),(/"x_p_2d","y_p_2d"/),(/"double","double"/))
    call add_var_att(trim(fname),"x_p_2d","long_name","zonal distance on p-grid")
    call add_var_att(trim(fname),"x_p_2d","units","m")
    call add_var_att(trim(fname),"y_p_2d","long_name","meridional distance on p-grid")
    call add_var_att(trim(fname),"y_p_2d","units","m")
    call writenet_def_var(trim(fname),2,2,(/"x_u","y_u"/),(/"x_u_2d","y_u_2d"/),(/"double","double"/))
    call add_var_att(trim(fname),"x_u_2d","long_name","zonal distance on u-grid")
    call add_var_att(trim(fname),"x_u_2d","units","m")
    call add_var_att(trim(fname),"y_u_2d","long_name","meridional distance on u-grid")
    call add_var_att(trim(fname),"y_u_2d","units","m")
    call writenet_def_var(trim(fname),2,2,(/"x_v","y_v"/),(/"x_v_2d","y_v_2d"/),(/"double","double"/))
    call add_var_att(trim(fname),"x_v_2d","long_name","zonal distance on v-grid")
    call add_var_att(trim(fname),"x_v_2d","units","m")
    call add_var_att(trim(fname),"y_v_2d","long_name","meridional distance on v-grid")
    call add_var_att(trim(fname),"y_v_2d","units","m")
    call add_var_att(trim(fname),"time","long_name","time")
    call add_var_att(trim(fname),"time","units",ref_time)
    write(*,*)
    call add_var_att(trim(fname),"time","slip_ind",oset%slip_ind)
    call add_var_att(trim(fname),"time","cp_ocn",oset%cp_ocn)
    call add_var_att(trim(fname),"time","r_ocn_day",oset%r_ocn_day)
    call add_var_att(trim(fname),"time","H1",oset%H1)
    call add_var_att(trim(fname),"time","H2",oset%H2)
    call add_var_att(trim(fname),"time","eps_s_ocn_day",oset%eps_s_ocn_day)
    call add_var_att(trim(fname),"time","cd_bulk",oset%cd_bulk)

    ! Output u
    dim_names(1)="x_u";dim_names(2)="y_u";dim_names(3)="time"
    call writenet_def_var(trim(fname),3,1,dim_names(1:3),(/"u_ocn_sw"/),(/"double"/))
    call add_var_att(trim(fname),"u_ocn_sw","long_name","SW zonal velocity")
    call add_var_att(trim(fname),"u_ocn_sw","units","m/s")
    call add_var_att(trim(fname),"u_ocn_sw","missing_value",missing_value)
    call writenet_def_var(trim(fname),3,1,dim_names(1:3),(/"u_ocn_ek"/),(/"double"/))
    call add_var_att(trim(fname),"u_ocn_ek","long_name","Ekman zonal velocity")
    call add_var_att(trim(fname),"u_ocn_ek","units","m/s")
    call add_var_att(trim(fname),"u_ocn_ek","missing_value",missing_value)
    call writenet_def_var(trim(fname),3,1,dim_names(1:3),(/"u_ocn_1"/),(/"double"/))
    call add_var_att(trim(fname),"u_ocn_1","long_name","zonal velocity")
    call add_var_att(trim(fname),"u_ocn_1","units","m/s")
    call add_var_att(trim(fname),"u_ocn_1","missing_value",missing_value)
    call writenet_def_var(trim(fname),2,1,dim_names(1:2),(/"mask_u"/),(/"double"/))
    call add_var_att(trim(fname),"mask_u","long_name","mask of u-grid")
    call add_var_att(trim(fname),"mask_u","units","")
    call add_var_att(trim(fname),"mask_u","missing_value",missing_value)
    call writenet_wv(trim(fname),"mask_u",(/1,1/),(/ogrd%nx_p+1,ogrd%ny_p+2/),ogrd%mask_u%val(1:ogrd%nx_p+1,0:ogrd%ny_p+1))

    ! Output v field
    dim_names(1)="x_v";dim_names(2)="y_v";;dim_names(3)="time"
    call writenet_def_var(trim(fname),3,1,dim_names(1:3),(/"v_ocn_sw"/),(/"double"/))
    call add_var_att(trim(fname),"v_ocn_sw","long_name","SW meridional velocity")
    call add_var_att(trim(fname),"v_ocn_sw","units","m/s")
    call add_var_att(trim(fname),"v_ocn_sw","missing_value",missing_value)
    call writenet_def_var(trim(fname),3,1,dim_names(1:3),(/"v_ocn_ek"/),(/"double"/))
    call add_var_att(trim(fname),"v_ocn_ek","long_name","Ekman meridional velocity")
    call add_var_att(trim(fname),"v_ocn_ek","units","m/s")
    call add_var_att(trim(fname),"v_ocn_ek","missing_value",missing_value)
    call writenet_def_var(trim(fname),3,1,dim_names(1:3),(/"v_ocn_1"/),(/"double"/))
    call add_var_att(trim(fname),"v_ocn_1","long_name","meridional velocity")
    call add_var_att(trim(fname),"v_ocn_1","units","m/s")
    call add_var_att(trim(fname),"v_ocn_1","missing_value",missing_value)
    call writenet_def_var(trim(fname),2,1,dim_names(1:2),(/"mask_v"/),(/"double"/))
    call add_var_att(trim(fname),"mask_v","long_name","mask of v-grid")
    call add_var_att(trim(fname),"mask_v","units","")
    call add_var_att(trim(fname),"mask_v","missing_value",missing_value)
    call writenet_wv(trim(fname),"mask_v",(/1,1/),(/ogrd%nx_p+2,ogrd%ny_p+1/),ogrd%mask_v%val(0:ogrd%nx_p+1,1:ogrd%ny_p+1))

    ! Output p field
    dim_names(1)="x_p";dim_names(2)="y_p";dim_names(3)="time"
    call writenet_def_var(trim(fname),3,1,dim_names(1:3),(/"h_ocn_sw"/),(/"double"/))
    call add_var_att(trim(fname),"h_ocn_sw","long_name","thermocline depth")
    call add_var_att(trim(fname),"h_ocn_sw","units","m")
    call add_var_att(trim(fname),"h_ocn_sw","missing_value",missing_value)
    call writenet_def_var(trim(fname),3,1,dim_names(1:3),(/"w_ocn_ek"/),(/"double"/))
    call add_var_att(trim(fname),"w_ocn_ek","long_name","Ekman Vertical velocity")
    call add_var_att(trim(fname),"w_ocn_ek","units","m/s")
    call add_var_att(trim(fname),"w_ocn_ek","missing_value",missing_value)
    call writenet_def_var(trim(fname),3,1,dim_names(1:3),(/"w_ocn_1"/),(/"double"/))
    call add_var_att(trim(fname),"w_ocn_1","long_name","Vertical velocity")
    call add_var_att(trim(fname),"w_ocn_1","units","m/s")
    call add_var_att(trim(fname),"w_ocn_1","missing_value",missing_value)
    call writenet_def_var(trim(fname),2,1,dim_names(1:2),(/"mask_p"/),(/"double"/))
    call add_var_att(trim(fname),"mask_p","long_name","mask of p-grid")
    call add_var_att(trim(fname),"mask_p","units","")
    call add_var_att(trim(fname),"mask_p","missing_value",missing_value)
    call writenet_wv(trim(fname),"mask_p",(/1,1/),(/ogrd%nx_p+2,ogrd%ny_p+2/),ogrd%mask_p%val(0:ogrd%nx_p+1,0:ogrd%ny_p+1))

    call writenet_def_var(trim(fname),3,1,dim_names(1:3),(/"taux"/),(/"double"/))
    call add_var_att(trim(fname),"taux","long_name","Zonal wind stress")
    call add_var_att(trim(fname),"taux","units","N/m^2")
    call add_var_att(trim(fname),"taux","missing_value",missing_value)
    call writenet_def_var(trim(fname),3,1,dim_names(1:3),(/"tauy"/),(/"double"/))
    call add_var_att(trim(fname),"tauy","long_name","Meridional wind stress")
    call add_var_att(trim(fname),"tauy","units","N/m^2")
    call add_var_att(trim(fname),"tauy","missing_value",missing_value)
    ! Write coordinate variables
    call writenet_wv(trim(fname),"x_p_2d",(/1,1/),(/ogrd%nx_p+2,ogrd%ny_p+2/),ogrd%x_p%val)
    call writenet_wv(trim(fname),"y_p_2d",(/1,1/),(/ogrd%nx_p+2,ogrd%ny_p+2/),ogrd%y_p%val)
    call writenet_wv(trim(fname),"x_u_2d",(/1,1/),(/ogrd%nx_p+1,ogrd%ny_p+2/),ogrd%x_u%val)
    call writenet_wv(trim(fname),"y_u_2d",(/1,1/),(/ogrd%nx_p+1,ogrd%ny_p+2/),ogrd%y_u%val)
    call writenet_wv(trim(fname),"x_v_2d",(/1,1/),(/ogrd%nx_p+2,ogrd%ny_p+1/),ogrd%x_v%val)
    call writenet_wv(trim(fname),"y_v_2d",(/1,1/),(/ogrd%nx_p+2,ogrd%ny_p+1/),ogrd%y_v%val)
    call writenet_wv(trim(fname),"lon_p",(/1/),(/ogrd%nx_p+2/),ogrd%lon_p%val)
    call writenet_wv(trim(fname),"lat_p",(/1/),(/ogrd%ny_p+2/),ogrd%lat_p%val)
    call writenet_wv(trim(fname),"lon_u",(/1/),(/ogrd%nx_p+1/),ogrd%lon_u%val)
    call writenet_wv(trim(fname),"lat_u",(/1/),(/ogrd%ny_p+2/),ogrd%lat_u%val)
    call writenet_wv(trim(fname),"lon_v",(/1/),(/ogrd%nx_p+2/),ogrd%lon_v%val)
    call writenet_wv(trim(fname),"lat_v",(/1/),(/ogrd%ny_p+1/),ogrd%lat_v%val)
    call writenet_wv(trim(fname),"time",(/1/),(/nt/),time(1:nt))
  end subroutine create_avg_ocn_dyn_ZC
  subroutine create_avg_ocn_sst_ZC(fname,ogrd,oset,missing_value)
    implicit none
    character(len=*),intent(in) :: fname
    type(ocn_dta),intent(in) :: ogrd
    type(ocn_set),intent(in) :: oset
    real(idx),intent(in) :: missing_value
    character(len=maxlen) :: dim_names(3),dim_types(3)
    dim_names(1)="x_p";dim_names(2)="y_p";dim_names(3)="time"
    dim_types(1)="double";dim_types(2)="double";dim_types(3)="double"
    call writenet_def_var(trim(fname),3,1,dim_names(1:3),(/"ssta"/),(/"double"/))
    call add_var_att(trim(fname),"ssta","long_name","sst anomaly")
    call add_var_att(trim(fname),"ssta","units","degrees_celsius")
    call add_var_att(trim(fname),"ssta","missing_value",missing_value)    
    call writenet_def_var(trim(fname),2,1,dim_names(1:2),(/"mask_sst"/),(/"double"/))
    call add_var_att(trim(fname),"mask_sst","long_name","mask of SST grid")
    call add_var_att(trim(fname),"mask_sst","units","")
    call add_var_att(trim(fname),"mask_sst","missing_value",missing_value)
    call writenet_wv(trim(fname),"mask_sst",(/1,1/),(/ogrd%nx_p+2,ogrd%ny_p+2/),ogrd%mask_sst%val(0:ogrd%nx_p+1,0:ogrd%ny_p+1))
  end subroutine create_avg_ocn_sst_ZC

  subroutine oper_avg_ocn_dyn(ogrd)
    implicit none
    type(ocn_dta),intent(inout) :: ogrd
    ! U-sw
    ogrd%u_sw_avg%val=ogrd%u_sw_avg%val+ogrd%u_sw%val
    ogrd%v_sw_avg%val=ogrd%v_sw_avg%val+ogrd%v_sw%val
    ogrd%h_sw_avg%val=ogrd%h_sw_avg%val+ogrd%h_sw%val
    ogrd%u_ek_avg%val=ogrd%u_ek_avg%val+ogrd%u_ek%val
    ogrd%v_ek_avg%val=ogrd%v_ek_avg%val+ogrd%v_ek%val
    ogrd%w_ek_avg%val=ogrd%w_ek_avg%val+ogrd%w_ek%val
    ogrd%u_ocn_1_avg%val=ogrd%u_ocn_1_avg%val+ogrd%u_ocn_1%val
    ogrd%v_ocn_1_avg%val=ogrd%v_ocn_1_avg%val+ogrd%v_ocn_1%val       
    ogrd%w_ocn_1_avg%val=ogrd%w_ocn_1_avg%val+ogrd%w_ocn_1%val
    ogrd%taux_ocn_avg%val=ogrd%taux_ocn_avg%val+ogrd%taux_ocn%val
    ogrd%tauy_ocn_avg%val=ogrd%tauy_ocn_avg%val+ogrd%tauy_ocn%val
  end subroutine oper_avg_ocn_dyn
  subroutine oper_avg_ocn_sst(ogrd)
    implicit none
    type(ocn_dta),intent(inout) :: ogrd
    ogrd%ssta_ocn_avg%val=ogrd%ssta_ocn_avg%val+ogrd%ssta_ocn%val
  end subroutine oper_avg_ocn_sst
  subroutine write_avg_ocn_dyn(fname,ogrd,iavg_count,iavg)
    implicit none
    character(len=*),intent(in) :: fname
    type(ocn_dta),intent(inout) :: ogrd
    integer,intent(inout) :: iavg_count
    integer,intent(in) :: iavg
    real(idx), allocatable:: tmp_3d(:,:,:)

    ogrd%u_sw_avg%val=ogrd%u_sw_avg%val/iavg_count
    ogrd%u_sw_avg%val=ogrd%u_sw_avg%val*ogrd%mask_u%val+&
         & (1-ogrd%mask_u%val)*missing_value
    ogrd%u_ek_avg%val=ogrd%u_ek_avg%val/iavg_count
    ogrd%u_ek_avg%val=ogrd%u_ek_avg%val*ogrd%mask_u%val+&
         & (1-ogrd%mask_u%val)*missing_value
    ogrd%v_sw_avg%val=ogrd%v_sw_avg%val/iavg_count
    ogrd%v_sw_avg%val=ogrd%v_sw_avg%val*ogrd%mask_v%val+&
         & (1-ogrd%mask_v%val)*missing_value
    ogrd%v_ek_avg%val=ogrd%v_ek_avg%val/iavg_count
    ogrd%v_ek_avg%val=ogrd%v_ek_avg%val*ogrd%mask_v%val+&
         & (1-ogrd%mask_v%val)*missing_value

    ogrd%h_sw_avg%val=ogrd%h_sw_avg%val/iavg_count
    ogrd%h_sw_avg%val=ogrd%h_sw_avg%val*ogrd%mask_p%val+&
         & (1-ogrd%mask_p%val)*missing_value
    ogrd%w_ek_avg%val=ogrd%w_ek_avg%val/iavg_count
    ogrd%w_ek_avg%val=ogrd%w_ek_avg%val*ogrd%mask_p%val+&
         & (1-ogrd%mask_p%val)*missing_value
       
    ogrd%u_ocn_1_avg%val=ogrd%u_ocn_1_avg%val/iavg_count
    ogrd%u_ocn_1_avg%val=ogrd%u_ocn_1_avg%val*ogrd%mask_u%val+&
         & (1-ogrd%mask_u%val)*missing_value
    ogrd%v_ocn_1_avg%val=ogrd%v_ocn_1_avg%val/iavg_count
    ogrd%v_ocn_1_avg%val=ogrd%v_ocn_1_avg%val*ogrd%mask_v%val+&
         & (1-ogrd%mask_v%val)*missing_value
    ogrd%w_ocn_1_avg%val=ogrd%w_ocn_1_avg%val/iavg_count
    ogrd%w_ocn_1_avg%val=ogrd%w_ocn_1_avg%val*ogrd%mask_p%val+&
         & (1-ogrd%mask_p%val)*missing_value

    ogrd%taux_ocn_avg%val=ogrd%taux_ocn_avg%val/iavg_count
    ogrd%taux_ocn_avg%val=ogrd%taux_ocn_avg%val*ogrd%mask_p%val+&
         & (1-ogrd%mask_p%val)*missing_value
    ogrd%tauy_ocn_avg%val=ogrd%tauy_ocn_avg%val/iavg_count
    ogrd%tauy_ocn_avg%val=ogrd%tauy_ocn_avg%val*ogrd%mask_p%val+&
         & (1-ogrd%mask_p%val)*missing_value
    allocate(tmp_3d(1:ogrd%nx_p+1,0:ogrd%ny_p+1,1))
    tmp_3d(:,:,1)=ogrd%u_sw_avg%val
    call writenet_wv(trim(fname),"u_ocn_sw",(/1,1,iavg/),(/ogrd%nx_p+1,ogrd%ny_p+2,iavg/),tmp_3d)
    deallocate(tmp_3d)
    allocate(tmp_3d(1:ogrd%nx_p+1,0:ogrd%ny_p+1,1))
    tmp_3d(:,:,1)=ogrd%u_ek_avg%val
    call writenet_wv(trim(fname),"u_ocn_ek",(/1,1,iavg/),(/ogrd%nx_p+1,ogrd%ny_p+2,iavg/),tmp_3d)
    deallocate(tmp_3d)
    allocate(tmp_3d(1:ogrd%nx_p+1,0:ogrd%ny_p+1,1))
    tmp_3d(:,:,1)=ogrd%u_ocn_1_avg%val
    call writenet_wv(trim(fname),"u_ocn_1",(/1,1,iavg/),(/ogrd%nx_p+1,ogrd%ny_p+2,iavg/),tmp_3d)
    deallocate(tmp_3d)
    ! V
    allocate(tmp_3d(0:ogrd%nx_p+1,1:ogrd%ny_p+1,1))
    tmp_3d(:,:,1)=ogrd%v_sw_avg%val
    call writenet_wv(trim(fname),"v_ocn_sw",(/1,1,iavg/),(/ogrd%nx_p+2,ogrd%ny_p+1,iavg/),tmp_3d)
    deallocate(tmp_3d)
    allocate(tmp_3d(0:ogrd%nx_p+1,1:ogrd%ny_p+1,1))
    tmp_3d(:,:,1)=ogrd%v_ek_avg%val
    call writenet_wv(trim(fname),"v_ocn_ek",(/1,1,iavg/),(/ogrd%nx_p+2,ogrd%ny_p+1,iavg/),tmp_3d)
    deallocate(tmp_3d)
    allocate(tmp_3d(0:ogrd%nx_p+1,1:ogrd%ny_p+1,1))
    tmp_3d(:,:,1)=ogrd%v_ocn_1_avg%val
    call writenet_wv(trim(fname),"v_ocn_1",(/1,1,iavg/),(/ogrd%nx_p+2,ogrd%ny_p+1,iavg/),tmp_3d)
    deallocate(tmp_3d)
    ! h
    allocate(tmp_3d(0:ogrd%nx_p+1,0:ogrd%ny_p+1,1))
    tmp_3d(:,:,1)=ogrd%h_sw_avg%val
    call writenet_wv(trim(fname),"h_ocn_sw",(/1,1,iavg/),(/ogrd%nx_p+2,ogrd%ny_p+2,iavg/),tmp_3d)
    deallocate(tmp_3d)
    allocate(tmp_3d(0:ogrd%nx_p+1,0:ogrd%ny_p+1,1))
    tmp_3d(:,:,1)=ogrd%w_ek_avg%val
    call writenet_wv(trim(fname),"w_ocn_ek",(/1,1,iavg/),(/ogrd%nx_p+2,ogrd%ny_p+2,iavg/),tmp_3d)
    deallocate(tmp_3d)
    allocate(tmp_3d(0:ogrd%nx_p+1,0:ogrd%ny_p+1,1))
    tmp_3d(:,:,1)=ogrd%w_ocn_1_avg%val
    call writenet_wv(trim(fname),"w_ocn_1",(/1,1,iavg/),(/ogrd%nx_p+2,ogrd%ny_p+2,iavg/),tmp_3d)
    deallocate(tmp_3d)
    allocate(tmp_3d(0:ogrd%nx_p+1,0:ogrd%ny_p+1,1))
    tmp_3d(:,:,1)=ogrd%taux_ocn_avg%val
    call writenet_wv(trim(fname),"taux",(/1,1,iavg/),(/ogrd%nx_p+2,ogrd%ny_p+2,iavg/),tmp_3d)
    deallocate(tmp_3d)
    allocate(tmp_3d(0:ogrd%nx_p+1,0:ogrd%ny_p+1,1))
    tmp_3d(:,:,1)=ogrd%tauy_ocn_avg%val
    call writenet_wv(trim(fname),"tauy",(/1,1,iavg/),(/ogrd%nx_p+2,ogrd%ny_p+2,iavg/),tmp_3d)
    deallocate(tmp_3d)
    ogrd%u_sw_avg%val=0.0_idx;ogrd%u_ek_avg%val=0.0_idx;ogrd%u_ocn_1_avg%val=0.0_idx
    ogrd%v_sw_avg%val=0.0_idx;ogrd%v_ek_avg%val=0.0_idx;ogrd%v_ocn_1_avg%val=0.0_idx
    ogrd%h_sw_avg%val=0.0_idx;ogrd%w_ek_avg%val=0.0_idx;ogrd%w_ocn_1_avg%val=0.0_idx
    ogrd%taux_ocn_avg%val=0.0_idx;ogrd%tauy_ocn_avg%val=0.0_idx
  end subroutine write_avg_ocn_dyn
  subroutine write_avg_ocn_sst(fname,ogrd,iavg_count,iavg)
    implicit none
    character(len=*),intent(in) :: fname
    type(ocn_dta),intent(inout) :: ogrd
    integer,intent(inout) :: iavg_count
    integer,intent(in) :: iavg
    real(idx), allocatable:: tmp_3d(:,:,:)
    ogrd%ssta_ocn_avg%val=ogrd%ssta_ocn_avg%val/iavg_count
    ogrd%ssta_ocn_avg%val=ogrd%ssta_ocn_avg%val*ogrd%mask_sst%val+&
         & (1-ogrd%mask_sst%val)*missing_value
    allocate(tmp_3d(0:ogrd%nx_p+1,0:ogrd%ny_p+1,1))
    tmp_3d(:,:,1)=ogrd%ssta_ocn_avg%val
    call writenet_wv(trim(fname),"ssta",(/1,1,iavg/),(/ogrd%nx_p+2,ogrd%ny_p+2,iavg/),tmp_3d)
    deallocate(tmp_3d)
    ogrd%ssta_ocn_avg%val=0.0_idx
  end subroutine write_avg_ocn_sst

  subroutine oper_avg_atm(agrd)
    implicit none
    type(atm_dta),intent(inout) :: agrd
    agrd%qa_atm_avg%val=agrd%qa_atm_avg%val+agrd%qa_atm%val
    agrd%ua_atm_avg%val=agrd%ua_atm_avg%val+agrd%ua_atm%val
    agrd%va_atm_avg%val=agrd%va_atm_avg%val+agrd%va_atm%val
    agrd%pa_atm_avg%val=agrd%pa_atm_avg%val+agrd%pa_atm%val
    agrd%ssta_atm_avg%val=agrd%ssta_atm_avg%val+agrd%ssta_atm%val
    agrd%sstm_atm_avg%val=agrd%sstm_atm_avg%val+agrd%sstm_atm%val
  end subroutine oper_avg_atm
  subroutine write_avg_atm(fname,agrd,iavg_count,iavg)
    implicit none
    character(len=*),intent(in) :: fname
    type(atm_dta),intent(inout) :: agrd
    integer,intent(inout) :: iavg_count
    integer,intent(in) :: iavg
    real(idx), allocatable:: tmp_3d(:,:,:)
    integer :: nx,ny
    nx=agrd%nx_atm
    ny=agrd%ny_atm
    agrd%pa_atm_avg%val=agrd%pa_atm_avg%val/iavg_count
    agrd%qa_atm_avg%val=agrd%qa_atm_avg%val/iavg_count
    agrd%ua_atm_avg%val=agrd%ua_atm_avg%val/iavg_count
    agrd%va_atm_avg%val=agrd%va_atm_avg%val/iavg_count
    agrd%sstm_atm_avg%val=agrd%sstm_atm_avg%val/iavg_count
    agrd%ssta_atm_avg%val=agrd%ssta_atm_avg%val/iavg_count
    allocate(tmp_3d(1:nx,1:ny,1))
    tmp_3d(1:nx,1:ny,1)=agrd%qa_atm_avg%val
    call writenet_wv(trim(fname),"q_atm",(/1,1,iavg/),(/nx,ny,iavg/),tmp_3d)
    deallocate(tmp_3d)
    allocate(tmp_3d(1:nx,1:ny,1))
    tmp_3d(1:nx,1:ny,1)=agrd%ua_atm_avg%val
    call writenet_wv(trim(fname),"u_atm",(/1,1,iavg/),(/nx,ny,iavg/),tmp_3d)
    deallocate(tmp_3d)
    allocate(tmp_3d(1:nx,1:ny,1))
    tmp_3d(1:nx,1:ny,1)=agrd%va_atm_avg%val
    call writenet_wv(trim(fname),"v_atm",(/1,1,iavg/),(/nx,ny,iavg/),tmp_3d)
    deallocate(tmp_3d)
    allocate(tmp_3d(1:nx,1:ny,1))
    tmp_3d(1:nx,1:ny,1)=agrd%pa_atm_avg%val
    call writenet_wv(trim(fname),"p_atm",(/1,1,iavg/),(/nx,ny,iavg/),tmp_3d)
    deallocate(tmp_3d)
    allocate(tmp_3d(1:nx,1:ny,1))
    tmp_3d(1:nx,1:ny,1)=agrd%ssta_atm_avg%val
    call writenet_wv(trim(fname),"ssta_atm",(/1,1,iavg/),(/nx,ny,iavg/),tmp_3d)
    deallocate(tmp_3d)
    allocate(tmp_3d(1:nx,1:ny,1))
    tmp_3d(1:nx,1:ny,1)=agrd%sstm_atm_avg%val
    call writenet_wv(trim(fname),"sstm_atm",(/1,1,iavg/),(/nx,ny,iavg/),tmp_3d)
    deallocate(tmp_3d)
    agrd%pa_atm_avg%val=0.0_idx
    agrd%ua_atm_avg%val=0.0_idx
    agrd%va_atm_avg%val=0.0_idx
    agrd%qa_atm_avg%val=0.0_idx
    agrd%ssta_atm_avg%val=0.0_idx
    agrd%sstm_atm_avg%val=0.0_idx
  end subroutine write_avg_atm
end module mod_io_avg
