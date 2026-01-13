module mod_io_diag
  use run_params
  use run_types
  use calendar_sub
  use ncdf_write
  implicit none
  ! Time array
  integer :: ntime_diag,idiag,idiag_count
  real(idx),allocatable :: time_diag(:)
  integer,allocatable :: istep_diag(:)  
  integer :: out_diag_flag,out_diag_int
  character(maxlen) :: fname_out_diag
  ! Namelist
  namelist/output_diag/out_diag_flag,out_diag_int
  namelist/output_diag/fname_out_diag
contains
  !==========================!
  ! Creation of average file !
  !==========================!
  subroutine initialize_ocn_diag(grd)
    implicit none
    type(ocn_dta),intent(inout) :: grd
    integer :: nx_p,ny_p
    nx_p=grd%nx_p;ny_p=grd%ny_p
    allocate(grd%dTdt_ocn_diag%val(0:nx_p+1,0:ny_p+1));grd%dTdt_ocn_diag%val(0:nx_p+1,0:ny_p+1)=0.0_idx
    allocate(grd%uaTm_ocn_diag%val(0:nx_p+1,0:ny_p+1));grd%uaTm_ocn_diag%val(0:nx_p+1,0:ny_p+1)=0.0_idx
    allocate(grd%umTa_ocn_diag%val(0:nx_p+1,0:ny_p+1));grd%umTa_ocn_diag%val(0:nx_p+1,0:ny_p+1)=0.0_idx
    allocate(grd%uaTa_ocn_diag%val(0:nx_p+1,0:ny_p+1));grd%uaTa_ocn_diag%val(0:nx_p+1,0:ny_p+1)=0.0_idx
    allocate(grd%vaTm_ocn_diag%val(0:nx_p+1,0:ny_p+1));grd%vaTm_ocn_diag%val(0:nx_p+1,0:ny_p+1)=0.0_idx
    allocate(grd%vmTa_ocn_diag%val(0:nx_p+1,0:ny_p+1));grd%vmTa_ocn_diag%val(0:nx_p+1,0:ny_p+1)=0.0_idx
    allocate(grd%vaTa_ocn_diag%val(0:nx_p+1,0:ny_p+1));grd%vaTa_ocn_diag%val(0:nx_p+1,0:ny_p+1)=0.0_idx
    allocate(grd%waTm_ocn_diag%val(0:nx_p+1,0:ny_p+1));grd%waTm_ocn_diag%val(0:nx_p+1,0:ny_p+1)=0.0_idx
    allocate(grd%wmTa_ocn_diag%val(0:nx_p+1,0:ny_p+1));grd%wmTa_ocn_diag%val(0:nx_p+1,0:ny_p+1)=0.0_idx
    allocate(grd%waTa_ocn_diag%val(0:nx_p+1,0:ny_p+1));grd%waTa_ocn_diag%val(0:nx_p+1,0:ny_p+1)=0.0_idx
    allocate(grd%qh_ocn_diag%val(0:nx_p+1,0:ny_p+1)) ;grd%qh_ocn_diag%val(0:nx_p+1,0:ny_p+1)=0.0_idx
    allocate(grd%Tsub_ocn_diag%val(0:nx_p+1,0:ny_p+1)) ;grd%Tsub_ocn_diag%val(0:nx_p+1,0:ny_p+1)=0.0_idx
    allocate(grd%Te_ocn_diag%val(0:nx_p+1,0:ny_p+1)) ;grd%Te_ocn_diag%val(0:nx_p+1,0:ny_p+1)=0.0_idx
  end subroutine initialize_ocn_diag
  subroutine create_diag_ocn_ZC(fname,ogrd,&
       & start_yymmdd,start_hhmmss,end_yymmdd,end_hhmmss,dt,&
       & out_flag,out_int,missing_value,istep,nt)
    implicit none
    character(len=*),intent(in) :: fname
    type(ocn_dta),intent(in) :: ogrd
    real(idx),intent(in):: dt
    integer,intent(in) :: start_yymmdd,start_hhmmss,end_yymmdd,end_hhmmss,out_flag,out_int
    real(idx),intent(in) :: missing_value
    integer,allocatable,intent(inout) :: istep(:)
    integer,intent(inout) :: nt
    real(idx),allocatable :: time(:)
    integer :: it
    character(len=maxlen) :: ref_time,varname
    integer :: tmp_yymmdd,tmp_hhmmss
    character(len=maxlen) :: dim_names(7),dim_types(7)
    real(idx) :: tmp
    ref_time=calendar_create_time_att(start_yymmdd,start_hhmmss,1)
    call calendar_cal_length_ymdhms(start_yymmdd,start_hhmmss,end_yymmdd,end_hhmmss,out_flag,tmp)
    nt= int(tmp / out_int)
    nt=max(nt,1)
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

    ! Output p field
    dim_names(1)="x_p";dim_names(2)="y_p";dim_names(3)="time"
    varname="dSSTdt"
    call writenet_def_var(trim(fname),3,1,dim_names(1:3),(/trim(varname)/),(/"double"/))
    call add_var_att(trim(fname),trim(varname),"long_name","SST tendency")
    call add_var_att(trim(fname),trim(varname),"units","degrees_celsius/s")
    call add_var_att(trim(fname),trim(varname),"missing_value",missing_value)
    varname="uaTm"
    call writenet_def_var(trim(fname),3,1,dim_names(1:3),(/trim(varname)/),(/"double"/))
    call add_var_att(trim(fname),trim(varname),"long_name",&
         & "Advection of climatological zonal SST gradient by anomalous zonal current")
    call add_var_att(trim(fname),trim(varname),"units","degrees_celsius/s")
    call add_var_att(trim(fname),trim(varname),"missing_value",missing_value)
    varname="umTa"
    call writenet_def_var(trim(fname),3,1,dim_names(1:3),(/trim(varname)/),(/"double"/))
    call add_var_att(trim(fname),trim(varname),"long_name",&
         & "Advection of anomalous zonal SST gradient by climatological zonal current")
    call add_var_att(trim(fname),trim(varname),"units","degrees_celsius/s")
    call add_var_att(trim(fname),trim(varname),"missing_value",missing_value)
    varname="uaTa"
    call writenet_def_var(trim(fname),3,1,dim_names(1:3),(/trim(varname)/),(/"double"/))
    call add_var_att(trim(fname),trim(varname),"long_name"&
         & ,"Advection of anomalous zonal SST gradient by anomalous zonal current")
    call add_var_att(trim(fname),trim(varname),"units","degrees_celsius/s")
    call add_var_att(trim(fname),trim(varname),"missing_value",missing_value)
    varname="vaTm"
    call writenet_def_var(trim(fname),3,1,dim_names(1:3),(/trim(varname)/),(/"double"/))
    call add_var_att(trim(fname),trim(varname),"long_name",&
         & "Advection of climatological meridional SST gradient by anomalous meridional current")
    call add_var_att(trim(fname),trim(varname),"units","degrees_celsius/s")
    call add_var_att(trim(fname),trim(varname),"missing_value",missing_value)
    varname="vmTa"
    call writenet_def_var(trim(fname),3,1,dim_names(1:3),(/trim(varname)/),(/"double"/))
    call add_var_att(trim(fname),trim(varname),"long_name",&
         & "Advection of anomalous meridional SST gradient by climatological meridional current")
    call add_var_att(trim(fname),trim(varname),"units","degrees_celsius/s")
    call add_var_att(trim(fname),trim(varname),"missing_value",missing_value)
    varname="vaTa"
    call writenet_def_var(trim(fname),3,1,dim_names(1:3),(/trim(varname)/),(/"double"/))
    call add_var_att(trim(fname),trim(varname),"long_name",&
          & "Advection of anomalous meridional SST gradient by anomalous meridional current")
    call add_var_att(trim(fname),trim(varname),"units","degrees_celsius/s")
    call add_var_att(trim(fname),trim(varname),"missing_value",missing_value)
    varname="waTm"
    call writenet_def_var(trim(fname),3,1,dim_names(1:3),(/trim(varname)/),(/"double"/))
    call add_var_att(trim(fname),trim(varname),"long_name",&
         & "Advection of climatological vertical SST gradient by anomalous vertical current")
    call add_var_att(trim(fname),trim(varname),"units","degrees_celsius/s")
    call add_var_att(trim(fname),trim(varname),"missing_value",missing_value)
    varname="wmTa"
    call writenet_def_var(trim(fname),3,1,dim_names(1:3),(/trim(varname)/),(/"double"/))
    call add_var_att(trim(fname),trim(varname),"long_name",&
         & "Advection of anomalous vertical SST gradient by climatological vertical current")
    call add_var_att(trim(fname),trim(varname),"units","degrees_celsius/s")
    call add_var_att(trim(fname),trim(varname),"missing_value",missing_value)
    varname="waTa"
    call writenet_def_var(trim(fname),3,1,dim_names(1:3),(/trim(varname)/),(/"double"/))
    call add_var_att(trim(fname),trim(varname),"long_name",&
          & "Advection of anomalous vertical SST gradient by anomalous vertical current")
    call add_var_att(trim(fname),trim(varname),"units","degrees_celsius/s")
    call add_var_att(trim(fname),trim(varname),"missing_value",missing_value)
    varname="qh"
    call writenet_def_var(trim(fname),3,1,dim_names(1:3),(/trim(varname)/),(/"double"/))
    call add_var_att(trim(fname),trim(varname),"long_name",&
          & "Thermal damping")
    call add_var_att(trim(fname),trim(varname),"units","degrees_celsius/s")
    call add_var_att(trim(fname),trim(varname),"missing_value",missing_value)
    varname="Tsub"
    call writenet_def_var(trim(fname),3,1,dim_names(1:3),(/trim(varname)/),(/"double"/))
    call add_var_att(trim(fname),trim(varname),"long_name",&
          & "Subsurface temperature")
    call add_var_att(trim(fname),trim(varname),"units","degrees_celsius")
    call add_var_att(trim(fname),trim(varname),"missing_value",missing_value)
    varname="Te"
    call writenet_def_var(trim(fname),3,1,dim_names(1:3),(/trim(varname)/),(/"double"/))
    call add_var_att(trim(fname),trim(varname),"long_name",&
          & "Entrainment temperature")
    call add_var_att(trim(fname),trim(varname),"units","degrees_celsius")
    call add_var_att(trim(fname),trim(varname),"missing_value",missing_value)    
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
  end subroutine create_diag_ocn_ZC
  subroutine write_diag_ocn(fname,ogrd,idiag_count,idiag)
    implicit none
    character(len=*),intent(in) :: fname
    type(ocn_dta),intent(inout) :: ogrd
    integer,intent(inout) :: idiag_count
    integer,intent(in) :: idiag
    real(idx), allocatable:: tmp_3d(:,:,:)
    ogrd%dTdt_ocn_diag%val=ogrd%dTdt_ocn_diag%val/idiag_count
    ogrd%dTdt_ocn_diag%val=ogrd%dTdt_ocn_diag%val*ogrd%mask_p%val+&
         & (1-ogrd%mask_p%val)*missing_value
    ogrd%uaTm_ocn_diag%val=ogrd%uaTm_ocn_diag%val/idiag_count
    ogrd%uaTm_ocn_diag%val=ogrd%uaTm_ocn_diag%val*ogrd%mask_p%val+&
         & (1-ogrd%mask_p%val)*missing_value
    ogrd%umTa_ocn_diag%val=ogrd%umTa_ocn_diag%val/idiag_count
    ogrd%umTa_ocn_diag%val=ogrd%umTa_ocn_diag%val*ogrd%mask_p%val+&
         & (1-ogrd%mask_p%val)*missing_value
    ogrd%uaTa_ocn_diag%val=ogrd%uaTa_ocn_diag%val/idiag_count
    ogrd%uaTa_ocn_diag%val=ogrd%uaTa_ocn_diag%val*ogrd%mask_p%val+&
         & (1-ogrd%mask_p%val)*missing_value
    ogrd%vaTm_ocn_diag%val=ogrd%vaTm_ocn_diag%val/idiag_count
    ogrd%vaTm_ocn_diag%val=ogrd%vaTm_ocn_diag%val*ogrd%mask_p%val+&
         & (1-ogrd%mask_p%val)*missing_value
    ogrd%vmTa_ocn_diag%val=ogrd%vmTa_ocn_diag%val/idiag_count
    ogrd%vmTa_ocn_diag%val=ogrd%vmTa_ocn_diag%val*ogrd%mask_p%val+&
         & (1-ogrd%mask_p%val)*missing_value
    ogrd%vaTa_ocn_diag%val=ogrd%vaTa_ocn_diag%val/idiag_count
    ogrd%vaTa_ocn_diag%val=ogrd%vaTa_ocn_diag%val*ogrd%mask_p%val+&
         & (1-ogrd%mask_p%val)*missing_value
    ogrd%waTm_ocn_diag%val=ogrd%waTm_ocn_diag%val/idiag_count
    ogrd%waTm_ocn_diag%val=ogrd%waTm_ocn_diag%val*ogrd%mask_p%val+&
         & (1-ogrd%mask_p%val)*missing_value
    ogrd%wmTa_ocn_diag%val=ogrd%wmTa_ocn_diag%val/idiag_count
    ogrd%wmTa_ocn_diag%val=ogrd%wmTa_ocn_diag%val*ogrd%mask_p%val+&
         & (1-ogrd%mask_p%val)*missing_value
    ogrd%waTa_ocn_diag%val=ogrd%waTa_ocn_diag%val/idiag_count
    ogrd%waTa_ocn_diag%val=ogrd%waTa_ocn_diag%val*ogrd%mask_p%val+&
         & (1-ogrd%mask_p%val)*missing_value
    ogrd%qh_ocn_diag%val=ogrd%qh_ocn_diag%val/idiag_count
    ogrd%qh_ocn_diag%val=ogrd%qh_ocn_diag%val*ogrd%mask_p%val+&
         & (1-ogrd%mask_p%val)*missing_value
    ogrd%Tsub_ocn_diag%val=ogrd%Tsub_ocn_diag%val/idiag_count
    ogrd%Tsub_ocn_diag%val=ogrd%Tsub_ocn_diag%val*ogrd%mask_p%val+&
         & (1-ogrd%mask_p%val)*missing_value
    ogrd%Te_ocn_diag%val=ogrd%Te_ocn_diag%val/idiag_count
    ogrd%Te_ocn_diag%val=ogrd%Te_ocn_diag%val*ogrd%mask_p%val+&
         & (1-ogrd%mask_p%val)*missing_value
    allocate(tmp_3d(0:ogrd%nx_p+1,0:ogrd%ny_p+1,1))
    tmp_3d(:,:,1)=ogrd%dTdt_ocn_diag%val
    call writenet_wv(trim(fname),"dSSTdt",(/1,1,idiag/),(/ogrd%nx_p+2,ogrd%ny_p+2,idiag/),tmp_3d)

    tmp_3d(:,:,1)=ogrd%uaTm_ocn_diag%val
    call writenet_wv(trim(fname),"uaTm",(/1,1,idiag/),(/ogrd%nx_p+2,ogrd%ny_p+2,idiag/),tmp_3d)
    tmp_3d(:,:,1)=ogrd%umTa_ocn_diag%val
    call writenet_wv(trim(fname),"umTa",(/1,1,idiag/),(/ogrd%nx_p+2,ogrd%ny_p+2,idiag/),tmp_3d)
    tmp_3d(:,:,1)=ogrd%uaTa_ocn_diag%val
    call writenet_wv(trim(fname),"uaTa",(/1,1,idiag/),(/ogrd%nx_p+2,ogrd%ny_p+2,idiag/),tmp_3d)

    tmp_3d(:,:,1)=ogrd%vaTm_ocn_diag%val
    call writenet_wv(trim(fname),"vaTm",(/1,1,idiag/),(/ogrd%nx_p+2,ogrd%ny_p+2,idiag/),tmp_3d)
    tmp_3d(:,:,1)=ogrd%vmTa_ocn_diag%val
    call writenet_wv(trim(fname),"vmTa",(/1,1,idiag/),(/ogrd%nx_p+2,ogrd%ny_p+2,idiag/),tmp_3d)
    tmp_3d(:,:,1)=ogrd%vaTa_ocn_diag%val
    call writenet_wv(trim(fname),"vaTa",(/1,1,idiag/),(/ogrd%nx_p+2,ogrd%ny_p+2,idiag/),tmp_3d)
    tmp_3d(:,:,1)=ogrd%waTm_ocn_diag%val
    call writenet_wv(trim(fname),"waTm",(/1,1,idiag/),(/ogrd%nx_p+2,ogrd%ny_p+2,idiag/),tmp_3d)
    tmp_3d(:,:,1)=ogrd%wmTa_ocn_diag%val
    call writenet_wv(trim(fname),"wmTa",(/1,1,idiag/),(/ogrd%nx_p+2,ogrd%ny_p+2,idiag/),tmp_3d)
    tmp_3d(:,:,1)=ogrd%waTa_ocn_diag%val
    call writenet_wv(trim(fname),"waTa",(/1,1,idiag/),(/ogrd%nx_p+2,ogrd%ny_p+2,idiag/),tmp_3d)
    tmp_3d(:,:,1)=ogrd%qh_ocn_diag%val
    call writenet_wv(trim(fname),"qh",(/1,1,idiag/),(/ogrd%nx_p+2,ogrd%ny_p+2,idiag/),tmp_3d)
    tmp_3d(:,:,1)=ogrd%Tsub_ocn_diag%val
    call writenet_wv(trim(fname),"Tsub",(/1,1,idiag/),(/ogrd%nx_p+2,ogrd%ny_p+2,idiag/),tmp_3d)
    tmp_3d(:,:,1)=ogrd%Te_ocn_diag%val
    call writenet_wv(trim(fname),"Te",(/1,1,idiag/),(/ogrd%nx_p+2,ogrd%ny_p+2,idiag/),tmp_3d)
    deallocate(tmp_3d)
    ogrd%dTdt_ocn_diag%val=0.0_idx
    ogrd%uaTm_ocn_diag%val=0.0_idx
    ogrd%umTa_ocn_diag%val=0.0_idx
    ogrd%uaTa_ocn_diag%val=0.0_idx
    ogrd%vaTm_ocn_diag%val=0.0_idx
    ogrd%vmTa_ocn_diag%val=0.0_idx
    ogrd%vaTa_ocn_diag%val=0.0_idx
    ogrd%waTm_ocn_diag%val=0.0_idx
    ogrd%wmTa_ocn_diag%val=0.0_idx
    ogrd%waTa_ocn_diag%val=0.0_idx
    ogrd%qh_ocn_diag%val=0.0_idx
    ogrd%Tsub_ocn_diag%val=0.0_idx
    ogrd%Te_ocn_diag%val=0.0_idx
  end subroutine write_diag_ocn
end module mod_io_diag
