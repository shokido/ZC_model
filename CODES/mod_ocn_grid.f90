module mod_ocn_dta
  use ncdf_read
  use ncdf_write
  use run_types
  implicit none
contains
    subroutine read_ocn_dyn_grd(fname,grd)
    implicit none
    character(len=maxlen),intent(in) :: fname
    type(ocn_dta),intent(inout) :: grd
    integer :: ntmp
    real(idx),allocatable :: v_1d(:)
    real(idx),allocatable :: v_2d(:,:)
    call get_dimsize(fname,"x_p",ntmp)
    grd%nx_p=ntmp-2
    call get_dimsize(fname,"y_p",ntmp)
    grd%ny_p=ntmp-2
    allocate(grd%lon_p%val(0:grd%nx_p+1))
    call get_variable(fname,"lon_p",(/1/),(/grd%nx_p+2/),v_1d)
    grd%lon_p%val(0:grd%nx_p+1)=v_1d(1:grd%nx_p+2)
    allocate(grd%lon_u%val(1:grd%nx_p+1))
    call get_variable(fname,"lon_u",(/1/),(/grd%nx_p+1/),v_1d)
    grd%lon_u%val(1:grd%nx_p+1)=v_1d(1:grd%nx_p+1)
    allocate(grd%lon_v%val(0:grd%nx_p+1))
    call get_variable(fname,"lon_v",(/1/),(/grd%nx_p+2/),v_1d)
    grd%lon_v%val(0:grd%nx_p+1)=v_1d(1:grd%nx_p+2)
    allocate(grd%lat_p%val(0:grd%nx_p+1))
    call get_variable(fname,"lat_p",(/1/),(/grd%ny_p+2/),v_1d)
    grd%lat_p%val(0:grd%ny_p+1)=v_1d(1:grd%ny_p+2)
    allocate(grd%lat_u%val(0:grd%nx_p+1))
    call get_variable(fname,"lat_u",(/1/),(/grd%ny_p+2/),v_1d)
    grd%lat_u%val(0:grd%ny_p+1)=v_1d(1:grd%ny_p+2)
    allocate(grd%lat_v%val(0:grd%nx_p+1))
    call get_variable(fname,"lat_v",(/1/),(/grd%ny_p+1/),v_1d)
    grd%lat_u%val(1:grd%ny_p+1)=v_1d(1:grd%ny_p+1)

    call get_variable(fname,"x_p_2d",(/1,1/),(/grd%nx_p+2,grd%ny_p+2/),v_2d)
    allocate(grd%x_p%val(0:grd%nx_p+1,0:grd%ny_p+1))
    grd%x_p%val(0:grd%nx_p+1,0:grd%ny_p+1)=v_2d(1:grd%nx_p+2,1:grd%ny_p+2)
    call get_variable(fname,"y_p_2d",(/1,1/),(/grd%nx_p+2,grd%ny_p+2/),v_2d)
    allocate(grd%y_p%val(0:grd%nx_p+1,0:grd%ny_p+1))
    grd%y_p%val(0:grd%nx_p+1,0:grd%ny_p+1)=v_2d(1:grd%nx_p+2,1:grd%ny_p+2)

    call get_variable(fname,"x_u_2d",(/1,1/),(/grd%nx_p+1,grd%ny_p+2/),v_2d)
    allocate(grd%x_u%val(1:grd%nx_p+1,0:grd%ny_p+1))
    grd%x_u%val(1:grd%nx_p+1,0:grd%ny_p+1)=v_2d(1:grd%nx_p+1,1:grd%ny_p+2)
    call get_variable(fname,"y_u_2d",(/1,1/),(/grd%nx_p+1,grd%ny_p+2/),v_2d)
    allocate(grd%y_u%val(1:grd%nx_p+1,0:grd%ny_p+1))
    grd%y_u%val(1:grd%nx_p+1,0:grd%ny_p+1)=v_2d(1:grd%nx_p+1,1:grd%ny_p+2)
    call get_variable(fname,"damp_u",(/1,1/),(/grd%nx_p+1,grd%ny_p+2/),v_2d)
    allocate(grd%damp_u%val(1:grd%nx_p+1,0:grd%ny_p+1))
    grd%damp_u%val(1:grd%nx_p+1,0:grd%ny_p+1)=v_2d(1:grd%nx_p+1,1:grd%ny_p+2)

    call get_variable(fname,"x_v_2d",(/1,1/),(/grd%nx_p+2,grd%ny_p+1/),v_2d)
    allocate(grd%x_v%val(0:grd%nx_p+1,1:grd%ny_p+1))
    grd%x_v%val(0:grd%nx_p+1,1:grd%ny_p+1)=v_2d(1:grd%nx_p+2,1:grd%ny_p+1)
    call get_variable(fname,"y_v_2d",(/1,1/),(/grd%nx_p+2,grd%ny_p+1/),v_2d)
    allocate(grd%y_v%val(0:grd%nx_p+1,1:grd%ny_p+1))
    grd%y_v%val(0:grd%nx_p+1,1:grd%ny_p+1)=v_2d(1:grd%nx_p+2,1:grd%ny_p+1)
    call get_variable(fname,"damp_v",(/1,1/),(/grd%nx_p+2,grd%ny_p+1/),v_2d)
    allocate(grd%damp_v%val(0:grd%nx_p+1,1:grd%ny_p+1))
    grd%damp_v%val(0:grd%nx_p+1,1:grd%ny_p+1)=v_2d(1:grd%nx_p+2,1:grd%ny_p+1)

    call get_variable(fname,"f",(/1,1/),(/grd%nx_p+2,grd%ny_p+2/),v_2d)
    allocate(grd%f%val(0:grd%nx_p+1,0:grd%ny_p+1))
    grd%f%val(0:grd%nx_p+1,0:grd%ny_p+1)=v_2d(1:grd%nx_p+2,1:grd%ny_p+2)
    call get_variable(fname,"mask_p",(/1,1/),(/grd%nx_p+2,grd%ny_p+2/),v_2d)
    allocate(grd%mask_p%val(0:grd%nx_p+1,0:grd%ny_p+1))
    grd%mask_p%val(0:grd%nx_p+1,0:grd%ny_p+1)=v_2d(1:grd%nx_p+2,1:grd%ny_p+2)
  end subroutine read_ocn_dyn_grd
    subroutine read_ocn_sst_grd(fname,grd)
    implicit none
    character(len=maxlen),intent(in) :: fname
    type(ocn_dta),intent(inout) :: grd
    integer :: ntmp
    real(idx),allocatable :: v_2d(:,:)
    call get_variable(fname,"mask_sst",(/1,1/),(/grd%nx_p+2,grd%ny_p+2/),v_2d)
    allocate(grd%mask_sst%val(0:grd%nx_p+1,0:grd%ny_p+1))
    grd%mask_sst%val(0:grd%nx_p+1,0:grd%ny_p+1)=v_2d(1:grd%nx_p+2,1:grd%ny_p+2)
  end subroutine read_ocn_sst_grd
  subroutine initialize_ocn_dyn(grd)
    implicit none
    type(ocn_dta),intent(inout) :: grd
    integer :: nx_p,ny_p
    nx_p=grd%nx_p;ny_p=grd%ny_p
    allocate(grd%h_sw%val(0:nx_p+1,0:ny_p+1));grd%h_sw%val(0:nx_p+1,0:ny_p+1)=0.0_idx
    allocate(grd%u_sw%val(1:nx_p+1,0:ny_p+1));grd%u_sw%val(1:nx_p+1,0:ny_p+1)=0.0_idx
    allocate(grd%v_sw%val(0:nx_p+1,1:ny_p+1));grd%v_sw%val(0:nx_p+1,1:ny_p+1)=0.0_idx
    allocate(grd%h_sw_past%val(0:nx_p+1,0:ny_p+1));grd%h_sw_past%val(0:nx_p+1,0:ny_p+1)=0.0_idx
    allocate(grd%u_sw_past%val(1:nx_p+1,0:ny_p+1));grd%u_sw_past%val(1:nx_p+1,0:ny_p+1)=0.0_idx
    allocate(grd%v_sw_past%val(0:nx_p+1,1:ny_p+1));grd%v_sw_past%val(0:nx_p+1,1:ny_p+1)=0.0_idx
    allocate(grd%h_sw_next%val(0:nx_p+1,0:ny_p+1));grd%h_sw_next%val(0:nx_p+1,0:ny_p+1)=0.0_idx
    allocate(grd%u_sw_next%val(1:nx_p+1,0:ny_p+1));grd%u_sw_next%val(1:nx_p+1,0:ny_p+1)=0.0_idx
    allocate(grd%v_sw_next%val(0:nx_p+1,1:ny_p+1));grd%v_sw_next%val(0:nx_p+1,1:ny_p+1)=0.0_idx
    allocate(grd%u_ek%val(1:nx_p+1,0:ny_p+1));grd%u_ek%val(1:nx_p+1,0:ny_p+1)=0.0_idx
    allocate(grd%v_ek%val(0:nx_p+1,1:ny_p+1));grd%v_ek%val(0:nx_p+1,1:ny_p+1)=0.0_idx
    allocate(grd%w_ek%val(0:nx_p+1,0:ny_p+1));grd%w_ek%val(0:nx_p+1,0:ny_p+1)=0.0_idx
    allocate(grd%u_ocn_1%val(1:nx_p+1,0:ny_p+1));grd%u_ocn_1%val(1:nx_p+1,0:ny_p+1)=0.0_idx
    allocate(grd%v_ocn_1%val(0:nx_p+1,1:ny_p+1));grd%v_ocn_1%val(0:nx_p+1,1:ny_p+1)=0.0_idx
    allocate(grd%w_ocn_1%val(0:nx_p+1,0:ny_p+1));grd%w_ocn_1%val(0:nx_p+1,0:ny_p+1)=0.0_idx
    allocate(grd%ua_ocn%val(0:nx_p+1,0:ny_p+1));grd%ua_ocn%val(0:nx_p+1,0:ny_p+1)=0.0_idx
    allocate(grd%va_ocn%val(0:nx_p+1,0:ny_p+1));grd%va_ocn%val(0:nx_p+1,0:ny_p+1)=0.0_idx
    allocate(grd%taux_ocn%val(0:nx_p+1,0:ny_p+1));grd%taux_ocn%val(0:nx_p+1,0:ny_p+1)=0.0_idx
    allocate(grd%tauy_ocn%val(0:nx_p+1,0:ny_p+1));grd%tauy_ocn%val(0:nx_p+1,0:ny_p+1)=0.0_idx
  end subroutine initialize_ocn_dyn
  subroutine initialize_ocn_sst(grd)
    implicit none
    type(ocn_dta),intent(inout) :: grd
    integer :: nx_p,ny_p
    nx_p=grd%nx_p;ny_p=grd%ny_p
    allocate(grd%ssta_ocn%val(0:nx_p+1,0:ny_p+1));grd%ssta_ocn%val(0:nx_p+1,0:ny_p+1)=0.0_idx
    allocate(grd%ssta_ocn_next%val(0:nx_p+1,0:ny_p+1));grd%ssta_ocn_next%val(0:nx_p+1,0:ny_p+1)=0.0_idx
  end subroutine initialize_ocn_sst
  subroutine initialize_ocn_visc(grd,oset)
    implicit none
    type(ocn_dta),intent(inout) :: grd
    type(ocn_set),intent(in) :: oset
    integer :: nx_p,ny_p
    real(idx),parameter :: scale=10
    integer,parameter :: nys=10
    nx_p=grd%nx_p;ny_p=grd%ny_p
    allocate(grd%damp_p%val(0:nx_p+1,0:ny_p+1))
    grd%damp_p%val(0:nx_p+1,0:ny_p+1)=(1.0_idx/(oset%r_ocn_day*day_to_sec))*grd%damp_p%val(0:nx_p+1,0:ny_p+1)
    grd%damp_u%val(1:nx_p+1,0:ny_p+1)=(1.0_idx/(oset%r_ocn_day*day_to_sec))*grd%damp_u%val(1:nx_p+1,0:ny_p+1)
    grd%damp_v%val(0:nx_p+1,1:ny_p+1)=(1.0_idx/(oset%r_ocn_day*day_to_sec))*grd%damp_v%val(0:nx_p+1,1:ny_p+1)
    allocate(grd%visc_2D%val(0:nx_p+1,0:ny_p+1))
    grd%visc_2D%val(0:nx_p+1,1:ny_p+1)=oset%nu
  end subroutine initialize_ocn_visc
  ! Masking
  subroutine set_mask_ocn(grd,slip)
    type(ocn_dta),intent(inout) :: grd
    real(idx),intent(in) :: slip
    integer :: iy,ix,nx,ny
    nx=grd%nx_p;ny=grd%ny_p
    ! slip_0 (dun/dx=0),
    ! slip=2 (un=0-> dun/dx = (un - (-un))/dx)
    allocate(grd%mask_u%val(1:nx+1,0:ny+1))
    allocate(grd%mask_v%val(0:nx+1,1:ny+1))
    allocate(grd%mask_phi_u%val(1:nx+1,1:ny+1))
    allocate(grd%mask_phi_v%val(1:nx+1,1:ny+1))
    grd%mask_u%val(:,:)=1.0_idx
    grd%mask_v%val(:,:)=1.0_idx
    grd%mask_phi_u%val(:,:)=1.0_idx
    grd%mask_phi_v%val(:,:)=1.0_idx
    iy = 0
    do ix = 1,nx+1
       grd%mask_u%val(ix,iy)=grd%mask_p%val(ix-1,iy)*grd%mask_p%val(ix,iy)
    end do

    do iy = 1,ny+1
       ix = 0
       grd%mask_v%val(ix,iy)=grd%mask_p%val(ix,iy-1)*grd%mask_p%val(ix,iy)       
       do ix = 1,nx+1
          grd%mask_u%val(ix,iy)=grd%mask_p%val(ix-1,iy)*grd%mask_p%val(ix,iy)
          grd%mask_v%val(ix,iy)=grd%mask_p%val(ix,iy-1)*grd%mask_p%val(ix,iy)
          if (grd%mask_p%val(ix-1,iy-1) == 0.0_idx .and. grd%mask_p%val(ix,iy-1)==0.0_idx) then
             grd%mask_phi_u%val(ix,iy) = slip
          end if
          if (grd%mask_p%val(ix-1,iy) == 0.0_idx .and. grd%mask_p%val(ix,iy)==0.0_idx) then
             grd%mask_phi_u%val(ix,iy) = slip
          end if
          if (grd%mask_p%val(ix-1,iy-1) == 0.0_idx .and. grd%mask_p%val(ix-1,iy) == 0.0_idx) then
             grd%mask_phi_v%val(ix,iy) = slip
          end if
          if (grd%mask_p%val(ix,iy-1) == 0.0_idx .and. grd%mask_p%val(ix,iy) == 0.0_idx) then
             grd%mask_phi_v%val(ix,iy) = slip
          end if
       end do
    end do
  end subroutine set_mask_ocn
  subroutine write_restart_ocn_dyn(fname,ogrd,missing_value)
    implicit none
    type(ocn_dta),intent(inout) :: ogrd
    character(len=maxlen),intent(in) :: fname
    real(idx),intent(in) :: missing_value
    character(len=maxlen) :: dim_names(6),dim_types(6)
    real(idx), allocatable:: tmp_2d(:,:)
    integer :: it
    dim_names(1)="x_p"
    dim_names(2)="y_p"
    dim_names(3)="x_u"
    dim_names(4)="y_u"
    dim_names(5)="x_v"
    dim_names(6)="y_v"
    do it=1,6
       dim_types(it)="double"
    end do

    call writenet_def_dim(trim(fname),6,(/ogrd%nx_p+2,ogrd%ny_p+2,&
         & ogrd%nx_p+1,ogrd%ny_p+2,ogrd%nx_p+2,ogrd%ny_p+1/),dim_names,dim_types)
    dim_names(1)="x_u";dim_names(2)="y_u"
    call writenet_def_var(trim(fname),2,1,dim_names(1:2),(/"u_sw"/),(/"double"/))
    call add_var_att(trim(fname),"u_sw","long_name","zonal velocity")
    call add_var_att(trim(fname),"u_sw","units","m/s")
    call add_var_att(trim(fname),"u_sw","missing_value",missing_value)
    call writenet_def_var(trim(fname),2,1,dim_names(1:2),(/"u_sw_past"/),(/"double"/))
    call add_var_att(trim(fname),"u_sw_past","long_name","zonal velocity")
    call add_var_att(trim(fname),"u_sw_past","units","m/s")
    call add_var_att(trim(fname),"u_sw_past","missing_value",missing_value)
    call writenet_def_var(trim(fname),2,1,dim_names(1:2),(/"u_sw_next"/),(/"double"/))
    call add_var_att(trim(fname),"u_sw_next","long_name","zonal velocity")
    call add_var_att(trim(fname),"u_sw_next","units","m/s")
    call add_var_att(trim(fname),"u_sw_next","missing_value",missing_value)
    ! U
    allocate(tmp_2d(1:ogrd%nx_p+1,1:ogrd%ny_p+2))
    tmp_2d(1:ogrd%nx_p+1,1:ogrd%ny_p+2)=ogrd%u_sw%val
    call writenet_wv(trim(fname),"u_sw",(/1,1/),(/ogrd%nx_p+1,ogrd%ny_p+2/),tmp_2d)
    tmp_2d(1:ogrd%nx_p+1,1:ogrd%ny_p+2)=ogrd%u_sw_past%val
    call writenet_wv(trim(fname),"u_sw_past",(/1,1/),(/ogrd%nx_p+1,ogrd%ny_p+2/),tmp_2d)
    tmp_2d(1:ogrd%nx_p+1,1:ogrd%ny_p+2)=ogrd%u_sw_next%val
    call writenet_wv(trim(fname),"u_sw_next",(/1,1/),(/ogrd%nx_p+1,ogrd%ny_p+2/),tmp_2d)
    deallocate(tmp_2d)

    ! V
    dim_names(1)="x_v";dim_names(2)="y_v"
    call writenet_def_var(trim(fname),2,1,dim_names(1:2),(/"v_sw"/),(/"double"/))
    call add_var_att(trim(fname),"v_sw","long_name","Meridional velocity")
    call add_var_att(trim(fname),"v_sw","units","m/s")
    call add_var_att(trim(fname),"v_sw","missing_value",missing_value)
    call writenet_def_var(trim(fname),2,1,dim_names(1:2),(/"v_sw_past"/),(/"double"/))
    call add_var_att(trim(fname),"v_sw_past","long_name","Meridional velocity")
    call add_var_att(trim(fname),"v_sw_past","units","m/s")
    call add_var_att(trim(fname),"v_sw_past","missing_value",missing_value)
    call writenet_def_var(trim(fname),2,1,dim_names(1:2),(/"v_sw_next"/),(/"double"/))
    call add_var_att(trim(fname),"v_sw_next","long_name","Meridional velocity")
    call add_var_att(trim(fname),"v_sw_next","units","m/s")
    call add_var_att(trim(fname),"v_sw_next","missing_value",missing_value)
    allocate(tmp_2d(1:ogrd%nx_p+2,1:ogrd%ny_p+1))
    tmp_2d(1:ogrd%nx_p+2,1:ogrd%ny_p+1)=ogrd%v_sw%val
    call writenet_wv(trim(fname),"v_sw",(/1,1/),(/ogrd%nx_p+2,ogrd%ny_p+1/),tmp_2d)
    tmp_2d(1:ogrd%nx_p+2,1:ogrd%ny_p+1)=ogrd%v_sw_past%val
    call writenet_wv(trim(fname),"v_sw_past",(/1,1/),(/ogrd%nx_p+2,ogrd%ny_p+1/),tmp_2d)
    tmp_2d(1:ogrd%nx_p+2,1:ogrd%ny_p+1)=ogrd%v_sw_next%val
    call writenet_wv(trim(fname),"v_sw_next",(/1,1/),(/ogrd%nx_p+2,ogrd%ny_p+1/),tmp_2d)
    deallocate(tmp_2d)
    ! P
    dim_names(1)="x_p";dim_names(2)="y_p"
    call writenet_def_var(trim(fname),2,1,dim_names(1:2),(/"h_sw"/),(/"double"/))
    call add_var_att(trim(fname),"h_sw","long_name","Meridional velocity")
    call add_var_att(trim(fname),"h_sw","units","m/s")
    call add_var_att(trim(fname),"h_sw","missing_value",missing_value)
    call writenet_def_var(trim(fname),2,1,dim_names(1:2),(/"h_sw_past"/),(/"double"/))
    call add_var_att(trim(fname),"h_sw_past","long_name","Meridional velocity")
    call add_var_att(trim(fname),"h_sw_past","units","m/s")
    call add_var_att(trim(fname),"h_sw_past","missing_value",missing_value)
    call writenet_def_var(trim(fname),2,1,dim_names(1:2),(/"h_sw_next"/),(/"double"/))
    call add_var_att(trim(fname),"h_sw_next","long_name","Meridional velocity")
    call add_var_att(trim(fname),"h_sw_next","units","m/s")
    call add_var_att(trim(fname),"h_sw_next","missing_value",missing_value)
    allocate(tmp_2d(1:ogrd%nx_p+2,1:ogrd%ny_p+2))
    tmp_2d(1:ogrd%nx_p+2,1:ogrd%ny_p+2)=ogrd%h_sw%val(0:ogrd%nx_p+1,0:ogrd%ny_p+1)
    call writenet_wv(trim(fname),"h_sw",(/1,1/),(/ogrd%nx_p+2,ogrd%ny_p+2/),tmp_2d)
    tmp_2d(1:ogrd%nx_p+2,1:ogrd%ny_p+2)=ogrd%h_sw_past%val(0:ogrd%nx_p+1,0:ogrd%ny_p+1)
    call writenet_wv(trim(fname),"h_sw_past",(/1,1/),(/ogrd%nx_p+2,ogrd%ny_p+2/),tmp_2d)
    tmp_2d(1:ogrd%nx_p+2,1:ogrd%ny_p+2)=ogrd%h_sw_next%val(0:ogrd%nx_p+1,0:ogrd%ny_p+1)
    call writenet_wv(trim(fname),"h_sw_next",(/1,1/),(/ogrd%nx_p+2,ogrd%ny_p+2/),tmp_2d)
    deallocate(tmp_2d)
  end subroutine write_restart_ocn_dyn

  subroutine write_restart_ocn_sst(fname,ogrd,missing_value)
    implicit none
    type(ocn_dta),intent(inout) :: ogrd
    character(len=maxlen),intent(in) :: fname
    real(idx),intent(in) :: missing_value
    character(len=maxlen) :: dim_names(2)
    real(idx), allocatable:: tmp_2d(:,:)
    ! P
    dim_names(1)="x_p";dim_names(2)="y_p"
    call writenet_def_var(trim(fname),2,1,dim_names(1:2),(/"ssta"/),(/"double"/))
    call add_var_att(trim(fname),"ssta","long_name","SST anomaly")
    call add_var_att(trim(fname),"ssta","units","degrees_celsius")
    call add_var_att(trim(fname),"ssta","missing_value",missing_value)
    allocate(tmp_2d(1:ogrd%nx_p+2,1:ogrd%ny_p+2))
    tmp_2d(1:ogrd%nx_p+2,1:ogrd%ny_p+2)=ogrd%ssta_ocn%val(0:ogrd%nx_p+1,0:ogrd%ny_p+1)
    call writenet_wv(trim(fname),"ssta",(/1,1/),(/ogrd%nx_p+2,ogrd%ny_p+2/),tmp_2d)
    deallocate(tmp_2d)
  end subroutine write_restart_ocn_sst
  subroutine read_restart_ocn_dyn(fname,grd)
    implicit none
    character(len=maxlen),intent(in) :: fname
    type(ocn_dta),intent(inout) :: grd
    integer :: nx_p,ny_p
    real(idx),allocatable :: v_2d(:,:)
    write(*,*) "Ocean dynamical initial condition is from "//trim(fname)
    nx_p=grd%nx_p;ny_p=grd%ny_p
    call get_variable(fname,"h_sw",(/1,1/),(/grd%nx_p+2,grd%ny_p+2/),v_2d)
    grd%h_sw%val(0:nx_p+1,0:ny_p+1)=v_2d(1:nx_p+2,1:ny_p+2)
    call get_variable(fname,"h_sw_past",(/1,1/),(/grd%nx_p+2,grd%ny_p+2/),v_2d)
    grd%h_sw_past%val(0:nx_p+1,0:ny_p+1)=v_2d(1:nx_p+2,1:ny_p+2)
    call get_variable(fname,"h_sw_next",(/1,1/),(/grd%nx_p+2,grd%ny_p+2/),v_2d)
    grd%h_sw_next%val(0:nx_p+1,0:ny_p+1)=v_2d(1:nx_p+2,1:ny_p+2)
    call get_variable(fname,"u_sw",(/1,1/),(/grd%nx_p+1,grd%ny_p+2/),v_2d)
    grd%u_sw%val(1:nx_p+1,0:ny_p+1)=v_2d(1:nx_p+1,1:ny_p+2)
    call get_variable(fname,"u_sw_past",(/1,1/),(/grd%nx_p+1,grd%ny_p+2/),v_2d)
    grd%u_sw_past%val(1:nx_p+1,0:ny_p+1)=v_2d(1:nx_p+1,1:ny_p+2)
    call get_variable(fname,"u_sw_next",(/1,1/),(/grd%nx_p+1,grd%ny_p+2/),v_2d)
    grd%u_sw_next%val(1:nx_p+1,0:ny_p+1)=v_2d(1:nx_p+1,1:ny_p+2)
    call get_variable(fname,"v_sw",(/1,1/),(/grd%nx_p+2,grd%ny_p+1/),v_2d)
    grd%v_sw%val(0:nx_p+1,1:ny_p+1)=v_2d(1:nx_p+2,1:ny_p+1)
    call get_variable(fname,"v_sw_past",(/1,1/),(/grd%nx_p+2,grd%ny_p+1/),v_2d)
    grd%v_sw_past%val(0:nx_p+1,1:ny_p+1)=v_2d(1:nx_p+2,1:ny_p+1)
    call get_variable(fname,"v_sw_next",(/1,1/),(/grd%nx_p+2,grd%ny_p+1/),v_2d)
    grd%v_sw_next%val(0:nx_p+1,1:ny_p+1)=v_2d(1:nx_p+2,1:ny_p+1)
  end subroutine read_restart_ocn_dyn
  subroutine read_restart_ocn_sst(fname,grd)
    implicit none
    character(len=maxlen),intent(in) :: fname
    type(ocn_dta),intent(inout) :: grd
    integer :: nx_p,ny_p
    real(idx),allocatable :: v_2d(:,:)
    write(*,*) "Ocean SST initial condition is from "//trim(fname)
    nx_p=grd%nx_p;ny_p=grd%ny_p
    call get_variable(fname,"ssta",(/1,1/),(/grd%nx_p+2,grd%ny_p+2/),v_2d)
    grd%ssta_ocn%val(0:nx_p+1,0:ny_p+1)=v_2d(1:nx_p+2,1:ny_p+2)
  end subroutine read_restart_ocn_sst
  subroutine deallocate_ocn_all(grd)
    implicit none
    type(ocn_dta),intent(inout) :: grd
    deallocate(grd%lon_p%val); deallocate(grd%lat_p%val)
    deallocate(grd%lon_u%val); deallocate(grd%lat_u%val)
    deallocate(grd%lon_v%val);deallocate(grd%lat_v%val)
    deallocate(grd%x_p%val); deallocate(grd%y_p%val)
    deallocate(grd%x_u%val); deallocate(grd%y_u%val)
    deallocate(grd%x_v%val);deallocate(grd%y_v%val)
    deallocate(grd%f%val)
    deallocate(grd%mask_p%val);deallocate(grd%mask_u%val);deallocate(grd%mask_v%val)
    deallocate(grd%h_sw%val); deallocate(grd%h_sw_past%val); deallocate(grd%h_sw_next%val)
    deallocate(grd%u_sw%val);deallocate(grd%u_sw_past%val);deallocate(grd%u_sw_next%val)
    deallocate(grd%v_sw%val); deallocate(grd%v_sw_past%val);deallocate(grd%v_sw_next%val)
    deallocate(grd%u_ek%val);deallocate(grd%v_ek%val)
    deallocate(grd%w_ek%val)
    deallocate(grd%u_ocn_1%val);deallocate(grd%v_ocn_1%val);deallocate(grd%w_ocn_1%val)
    deallocate(grd%ua_ocn%val);deallocate(grd%va_ocn%val)
    deallocate(grd%taux_ocn%val);deallocate(grd%tauy_ocn%val)
  end subroutine deallocate_ocn_all
  subroutine deallocate_ocn_sst(grd)
    implicit none
    type(ocn_dta),intent(inout) :: grd
    deallocate(grd%ssta_ocn%val)
    deallocate(grd%ssta_ocn_next%val)
  end subroutine deallocate_ocn_sst
end module mod_ocn_dta
