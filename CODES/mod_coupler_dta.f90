module mod_coupler_dta
  use ncdf_read
  use run_types
  implicit none
contains
  subroutine cal_dist_spherical_min(lon1,lon2,lat1,lat2,dist)
    ! lon1,lon2,lat1,lat2 should be in radian
    real(idx),intent(in) :: lon1,lon2,lat1,lat2
    real(idx),intent(inout) :: dist
    real(idx) :: dlon,dlat,a
    real(8),parameter :: pi=4.0*atan(1.0)
    dlon=lon1-lon2
    dlon=min(dlon,2.0*pi-dlon)
    dlat=lat1-lat2
    dlat=min(dlat,2.0*pi-dlat)
    a=sin(dlat/2.0) * sin(dlat/2.0) + &
         & cos(lat1) *cos(lat2) * sin(dlon/2.0) * sin(dlon/2.0)
    dist= 2.0_idx *atan2(sqrt(a),sqrt(1-a))
  end subroutine cal_dist_spherical_min
  subroutine get_nearest(nx_in,ny_in,lon_in,lat_in,mask_in,&
       & lon_tar,lat_tar,ind_x_out,ind_y_out)
    implicit none
    integer,intent(in) :: nx_in,ny_in
    real(8),intent(in) :: lon_in(nx_in,ny_in),lat_in(nx_in,ny_in)
    real(8),intent(in) :: mask_in(nx_in,ny_in)
    real(8),intent(in) :: lon_tar,lat_tar
    integer,intent(inout) :: ind_x_out,ind_y_out
    integer :: i,j,nval_in,imask
    integer,allocatable :: idx_i_in(:),idx_j_in(:)
    real(8),allocatable :: x_1d_in(:),y_1d_in(:)
    real(8) :: x_1d_out,y_1d_out
    real(8),parameter :: pi=4.0*atan(1.0)
    real(8) :: dist_tmp,dist

    nval_in=int(sum(mask_in))
    allocate(x_1d_in(nval_in));allocate(y_1d_in(nval_in))
    allocate(idx_i_in(nval_in));allocate(idx_j_in(nval_in))
    imask=1
    do j=1,ny_in       
       do i=1,nx_in
          if (mask_in(i,j).ne.0) then
             x_1d_in(imask)=lon_in(i,j)*pi/180.0
             y_1d_in(imask)=lat_in(i,j)*pi/180.0
             idx_i_in(imask)=i
             idx_j_in(imask)=j
             imask=imask+1
          end if
       end do
    end do
    dist_tmp=1.0e20
    x_1d_out=lon_tar*pi/180.0
    y_1d_out=lat_tar*pi/180.0

    do j=1,nval_in
       call cal_dist_spherical_min(x_1d_in(j),x_1d_out,y_1d_in(j),y_1d_out,dist)
       if (dist<=dist_tmp) then
          dist_tmp=dist
          ind_x_out=idx_i_in(j)
          ind_y_out=idx_j_in(j)
       end if
    end do
    deallocate(idx_i_in);deallocate(idx_j_in)
    deallocate(x_1d_in);deallocate(y_1d_in)
  end subroutine get_nearest
  subroutine read_coupler(cinfo,ogrd,agrd,fname)
    character(len=maxlen),intent(in) :: fname
    type(couple_dta),intent(inout) :: cinfo
    type(ocn_dta),intent(in) :: ogrd
    type(atm_dta),intent(in) :: agrd
    integer,allocatable :: i_3d(:,:,:)
    real(idx),allocatable :: v_3d(:,:,:)
    integer :: nx_a,ny_a,nx_o,ny_o
    integer :: ix,iy
    integer :: ix_t,iy_t
    nx_o=ogrd%nx_p
    ny_o=ogrd%ny_p
    nx_a=agrd%nx_atm
    ny_a=agrd%ny_atm
    allocate(cinfo%ind_AtoO_x%val(0:nx_o+1,0:ny_o+1,4))
    allocate(cinfo%ind_AtoO_y%val(0:nx_o+1,0:ny_o+1,4))
    allocate(cinfo%wgt_AtoO%val(0:nx_o+1,0:ny_o+1,4))
    allocate(cinfo%ind_OtoA_x%val(1:nx_a,1:ny_a,4))
    allocate(cinfo%ind_OtoA_y%val(1:nx_a,1:ny_a,4))
    allocate(cinfo%wgt_OtoA%val(1:nx_a,1:ny_a,4))
    call get_variable(fname,"ind_AtoO_x",(/1,1,1/),(/nx_o+2,ny_o+2,4/),i_3d)
    cinfo%ind_AtoO_x%val(0:nx_o+1,0:ny_o+1,1:4)=i_3d(1:nx_o+2,1:ny_o+2,1:4)
    call get_variable(fname,"ind_AtoO_y",(/1,1,1/),(/nx_o+2,ny_o+2,4/),i_3d)
    cinfo%ind_AtoO_y%val(0:nx_o+1,0:ny_o+1,1:4)=i_3d(1:nx_o+2,1:ny_o+2,1:4)
    call get_variable(fname,"wgt_AtoO",(/1,1,1/),(/nx_o+2,ny_o+2,4/),v_3d)
    cinfo%wgt_AtoO%val(0:nx_o+1,0:ny_o+1,1:4)=v_3d(1:nx_o+2,1:ny_o+2,1:4)
    call get_variable(fname,"ind_OtoA_x",(/1,1,1/),(/nx_a,ny_a,4/),i_3d)
    cinfo%ind_OtoA_x%val(1:nx_a,0:ny_a,1:4)=i_3d(1:nx_a,1:ny_a,1:4)
    call get_variable(fname,"ind_OtoA_y",(/1,1,1/),(/nx_a,ny_a,4/),i_3d)
    cinfo%ind_OtoA_y%val(1:nx_a,1:ny_a,1:4)=i_3d(1:nx_a,1:ny_a,1:4)
    call get_variable(fname,"wgt_OtoA",(/1,1,1/),(/nx_a,ny_a,4/),v_3d)
    cinfo%wgt_OtoA%val(1:nx_a,1:ny_a,1:4)=v_3d(1:nx_a,1:ny_a,1:4)
   end subroutine read_coupler
  subroutine initialize_coupler(cinfo,ogrd,agrd)
    implicit none
    type(couple_dta),intent(inout) :: cinfo
    type(ocn_dta),intent(in) :: ogrd
    type(atm_dta),intent(in) :: agrd
    integer :: nx_a,ny_a,nx_o,ny_o
    integer :: ix,iy
    integer :: ix_t,iy_t
    real(idx),allocatable :: lon_in(:,:),lat_in(:,:),mask_in(:,:)
    nx_o=ogrd%nx_p
    ny_o=ogrd%ny_p
    nx_a=agrd%nx_atm
    ny_a=agrd%ny_atm
    allocate(cinfo%ind_AtoO_x%val(0:nx_o+1,0:ny_o+1,4))
    allocate(cinfo%ind_AtoO_y%val(0:nx_o+1,0:ny_o+1,4))
    allocate(cinfo%wgt_AtoO%val(0:nx_o+1,0:ny_o+1,4))
    allocate(cinfo%ind_OtoA_x%val(1:nx_a,1:ny_a,4))
    allocate(cinfo%ind_OtoA_y%val(1:nx_a,1:ny_a,4))
    allocate(cinfo%wgt_OtoA%val(1:nx_a,1:ny_a,4))
    ! ATM -> OCN
    allocate(lon_in(1:nx_a,1:ny_a))
    allocate(lat_in(1:nx_a,1:ny_a))
    allocate(mask_in(1:nx_a,1:ny_a))
    do iy=1,ny_a
       do ix=1,nx_a
          lon_in(ix,iy)=agrd%lon_atm%val(ix)
          lat_in(ix,iy)=agrd%lat_atm%val(iy)
          mask_in(ix,iy)=1.0_idx
       end do
    end do
    do iy=0,ny_o+1
       do ix=0,nx_o+1
          if (ogrd%mask_p%val(ix,iy)==0) then
             cinfo%ind_AtoO_x%val(ix,iy,1:4)=1
             cinfo%ind_AtoO_y%val(ix,iy,1:4)=1
             cinfo%wgt_AtoO%val(ix,iy,1:4)=0.0_idx
          else
             call get_nearest(nx_a,ny_a,lon_in,lat_in,mask_in,&
                  & ogrd%lon_p%val(ix),ogrd%lat_p%val(iy),ix_t,iy_t)
             cinfo%ind_AtoO_x%val(ix,iy,1)=ix_t
             cinfo%ind_AtoO_y%val(ix,iy,1)=iy_t
             cinfo%wgt_AtoO%val(ix,iy,1)=1.0_idx
             cinfo%ind_AtoO_x%val(ix,iy,2:4)=1
             cinfo%ind_AtoO_y%val(ix,iy,2:4)=1
             cinfo%wgt_AtoO%val(ix,iy,2:4)=0.0_idx
          end if
       end do
    end do
    deallocate(lon_in)
    deallocate(lat_in)
    deallocate(mask_in)

    ! OCN to ATM
    allocate(lon_in(1:nx_o+2,1:ny_o+2))
    allocate(lat_in(1:nx_o+2,1:ny_o+2))
    allocate(mask_in(1:nx_o+2,1:ny_o+2))
    do iy=1,ny_o+2
       do ix=1,nx_o+2
          lon_in(ix,iy)=ogrd%lon_p%val(ix-1)
          lat_in(ix,iy)=ogrd%lat_p%val(iy-1)
          mask_in(ix,iy)=ogrd%mask_p%val(ix-1,iy-1)
       end do
    end do
    do iy=1,ny_a
       do ix=1,nx_a
          call get_nearest(nx_o+2,ny_o+2,lon_in,lat_in,mask_in,&
               & agrd%lon_atm%val(ix),agrd%lat_atm%val(iy),ix_t,iy_t)
          cinfo%ind_OtoA_x%val(ix,iy,1)=ix_t-1
          cinfo%ind_OtoA_y%val(ix,iy,1)=iy_t-1
          if ((agrd%lon_atm%val(ix) .ge. minval(lon_in)) .and. &
               & (agrd%lon_atm%val(ix) .le. maxval(lon_in)) .and. &
               & (agrd%lat_atm%val(iy) .ge. minval(lat_in)) .and. &
               & (agrd%lat_atm%val(iy) .le. maxval(lat_in))) then
             cinfo%wgt_OtoA%val(ix,iy,1)=1.0_idx
          else
             cinfo%wgt_OtoA%val(ix,iy,1)=0.0_idx
          end if
          cinfo%ind_OtoA_x%val(ix,iy,2:4)=1
          cinfo%ind_OtoA_y%val(ix,iy,2:4)=1
          cinfo%wgt_OtoA%val(ix,iy,2:4)=0.0_idx
       end do
    end do
    deallocate(lon_in)
    deallocate(lat_in)
    deallocate(mask_in)

  end subroutine initialize_coupler
  subroutine exchange_AtoO(cinfo,ogrd,agrd)
    ! Import atmospheric information onto ocean grid
    implicit none
    type(couple_dta),intent(in) :: cinfo
    type(ocn_dta),intent(inout) :: ogrd
    type(atm_dta),intent(inout) :: agrd
    integer :: nx_o,ny_o
    integer :: ix,iy
    integer :: ix_1,ix_2,ix_3,ix_4
    integer :: iy_1,iy_2,iy_3,iy_4
    nx_o=ogrd%nx_p
    ny_o=ogrd%ny_p
    do iy=0,ny_o+1
       do ix=0,nx_o+1
         ! Ensure mask-p is not zero
          ix_1=cinfo%ind_AtoO_x%val(ix,iy,1);iy_1=cinfo%ind_AtoO_y%val(ix,iy,1)
          ix_2=cinfo%ind_AtoO_x%val(ix,iy,2);iy_2=cinfo%ind_AtoO_y%val(ix,iy,2)
          ix_3=cinfo%ind_AtoO_x%val(ix,iy,3);iy_3=cinfo%ind_AtoO_y%val(ix,iy,3)
          ix_4=cinfo%ind_AtoO_x%val(ix,iy,4);iy_4=cinfo%ind_AtoO_y%val(ix,iy,4)
          ogrd%ua_ocn%val(ix,iy)= &
               & (cinfo%wgt_AtoO%val(ix,iy,1)*agrd%ua_atm%val(ix_1,iy_1)+&
               & cinfo%wgt_AtoO%val(ix,iy,2)*agrd%ua_atm%val(ix_2,iy_2)+&
               & cinfo%wgt_AtoO%val(ix,iy,3)*agrd%ua_atm%val(ix_3,iy_3)+&
               & cinfo%wgt_AtoO%val(ix,iy,4)*agrd%ua_atm%val(ix_4,iy_4))
          ogrd%va_ocn%val(ix,iy)=&
               & (cinfo%wgt_AtoO%val(ix,iy,1)*agrd%va_atm%val(ix_1,iy_1)+&
               & cinfo%wgt_AtoO%val(ix,iy,2)*agrd%va_atm%val(ix_2,iy_2)+&
               & cinfo%wgt_AtoO%val(ix,iy,3)*agrd%va_atm%val(ix_3,iy_3)+&
               & cinfo%wgt_AtoO%val(ix,iy,4)*agrd%va_atm%val(ix_4,iy_4))
       end do
    end do
  end subroutine exchange_AtoO
  subroutine exchange_OtoA(cinfo,ogrd,agrd)
    ! Import atmospheric information onto ocean grid
    implicit none
    type(couple_dta),intent(in) :: cinfo
    type(ocn_dta),intent(inout) :: ogrd
    type(atm_dta),intent(inout) :: agrd
    integer :: nx_a,ny_a
    integer :: ix,iy
    integer :: ix_1,ix_2,ix_3,ix_4
    integer :: iy_1,iy_2,iy_3,iy_4
    nx_a=agrd%nx_atm
    ny_a=agrd%ny_atm
    do iy=1,ny_a
       do ix=1,nx_a
          ix_1=cinfo%ind_OtoA_x%val(ix,iy,1);iy_1=cinfo%ind_OtoA_y%val(ix,iy,1)
          ix_2=cinfo%ind_OtoA_x%val(ix,iy,2);iy_2=cinfo%ind_OtoA_y%val(ix,iy,2)
          ix_3=cinfo%ind_OtoA_x%val(ix,iy,3);iy_3=cinfo%ind_OtoA_y%val(ix,iy,3)
          ix_4=cinfo%ind_OtoA_x%val(ix,iy,4);iy_4=cinfo%ind_OtoA_y%val(ix,iy,4)
          agrd%ssta_atm%val(ix,iy)= &
               & cinfo%wgt_OtoA%val(ix,iy,1)*ogrd%ssta_ocn%val(ix_1,iy_1)+&
               & cinfo%wgt_OtoA%val(ix,iy,2)*ogrd%ssta_ocn%val(ix_2,iy_2)+&
               & cinfo%wgt_OtoA%val(ix,iy,3)*ogrd%ssta_ocn%val(ix_3,iy_3)+&
               & cinfo%wgt_OtoA%val(ix,iy,4)*ogrd%ssta_ocn%val(ix_4,iy_4)
       end do
    end do
  end subroutine exchange_OtoA
  subroutine deallocate_coupler(cinfo)
    type(couple_dta),intent(inout) :: cinfo
    deallocate(cinfo%ind_AtoO_x%val)
    deallocate(cinfo%ind_AtoO_y%val)
    deallocate(cinfo%wgt_AtoO%val)
    deallocate(cinfo%ind_OtoA_x%val)
    deallocate(cinfo%ind_OtoA_y%val)
    deallocate(cinfo%wgt_OtoA%val)
  end subroutine deallocate_coupler
end module mod_coupler_dta
