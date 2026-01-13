module mod_atm_dta
  use ncdf_read
  use ncdf_write
  use run_types
  use run_params
  implicit none
contains
  subroutine read_atm_grd(fname,grd)
    implicit none
    character(len=maxlen),intent(in) :: fname
    type(atm_dta),intent(inout) :: grd
    integer :: ntmp
    real(idx),allocatable :: v_1d(:)
    call get_dimsize(fname,"lon",ntmp)
    grd%nx_atm=ntmp
    call get_dimsize(fname,"lat",ntmp)
    grd%ny_atm=ntmp
    allocate(grd%lon_atm%val(1:grd%nx_atm))
    call get_variable(fname,"lon",(/1/),(/grd%nx_atm/),v_1d)
    grd%lon_atm%val(1:grd%nx_atm)=v_1d(1:grd%nx_atm)
    allocate(grd%lat_atm%val(1:grd%nx_atm))
    call get_variable(fname,"lat",(/1/),(/grd%ny_atm/),v_1d)
    grd%lat_atm%val(1:grd%ny_atm)=v_1d(1:grd%ny_atm)
  end subroutine read_atm_grd

  subroutine set_coord_atm(agrd,aset)
    implicit none
    type(atm_dta),intent(inout) :: agrd
    type(atm_set),intent(in) :: aset
    integer :: ix,iy
    real(idx) :: k_factor,Lstar
    allocate(agrd%x_atm%val(agrd%nx_atm))
    allocate(agrd%y_atm%val(agrd%ny_atm))
    allocate(agrd%qa_atm%val(agrd%nx_atm,agrd%ny_atm)) ;
    agrd%qa_atm%val(1:agrd%nx_atm,1:agrd%ny_atm)=0.0_idx
    
    allocate(agrd%pa_atm%val(1:agrd%nx_atm,1:agrd%ny_atm))
    agrd%pa_atm%val(1:agrd%nx_atm,1:agrd%ny_atm)=0.0_idx
    allocate(agrd%ua_atm%val(1:agrd%nx_atm,1:agrd%ny_atm))
    agrd%ua_atm%val(1:agrd%nx_atm,1:agrd%ny_atm)=0.0_idx
    allocate(agrd%va_atm%val(1:agrd%nx_atm,1:agrd%ny_atm))
    agrd%va_atm%val(1:agrd%nx_atm,1:agrd%ny_atm)=0.0_idx
    allocate(agrd%qa_atm_avg%val(1:agrd%nx_atm,1:agrd%ny_atm)) ;
    agrd%qa_atm_avg%val(1:agrd%nx_atm,1:agrd%ny_atm)=0.0_idx
    allocate(agrd%pa_atm_avg%val(1:agrd%nx_atm,1:agrd%ny_atm))
    agrd%pa_atm_avg%val(1:agrd%nx_atm,1:agrd%ny_atm)=0.0_idx
    allocate(agrd%ua_atm_avg%val(1:agrd%nx_atm,1:agrd%ny_atm))
    agrd%ua_atm_avg%val(1:agrd%nx_atm,1:agrd%ny_atm)=0.0_idx
    allocate(agrd%va_atm_avg%val(1:agrd%nx_atm,1:agrd%ny_atm))
    agrd%va_atm_avg%val(1:agrd%nx_atm,1:agrd%ny_atm)=0.0_idx
    allocate(agrd%ssta_atm_avg%val(1:agrd%nx_atm,1:agrd%ny_atm))
    agrd%ssta_atm_avg%val(1:agrd%nx_atm,1:agrd%ny_atm)=0.0_idx
    allocate(agrd%sstm_atm_avg%val(1:agrd%nx_atm,1:agrd%ny_atm))
    agrd%sstm_atm_avg%val(1:agrd%nx_atm,1:agrd%ny_atm)=0.0_idx
    allocate(agrd%k_atm%val(1:agrd%nx_atm))

    allocate(agrd%sstm_atm%val(1:agrd%nx_atm,1:agrd%ny_atm))
    agrd%sstm_atm%val(1:agrd%nx_atm,1:agrd%ny_atm)=0.0_idx
    allocate(agrd%ssta_atm%val(1:agrd%nx_atm,1:agrd%ny_atm))
    agrd%ssta_atm%val(1:agrd%nx_atm,1:agrd%ny_atm)=0.0_idx

    Lstar = sqrt(aset%cp_atm/ (2.0_idx*beta))
    do ix = 1,agrd%nx_atm
       agrd%x_atm%val(ix) = lat_to_dis*agrd%lon_atm%val(ix) / Lstar
    end do
    do iy = 1,agrd%ny_atm
       agrd%y_atm%val(iy) = lat_to_dis*agrd%lat_atm%val(iy) / Lstar
    end do
    k_factor=2.0_idx*pi/(agrd%x_atm%val(agrd%nx_atm)- agrd%x_atm%val(1))
    !k_factor * xmax = 2.0 * pi
    !==========================================
    ! Wave number setting
    !==========================================
   if (mod(agrd%nx_atm,2) .eq. 0) then
      do ix =1,agrd%nx_atm/2+1
         agrd%k_atm%val(ix) = (ix-1) *k_factor
      end do
      do ix = agrd%nx_atm/2+2,agrd%nx_atm
         agrd%k_atm%val(ix) = (ix-1-agrd%nx_atm)*k_factor
      end do
   else
      agrd%k_atm%val(1) = 0
      do ix =2,agrd%nx_atm/2+1
         agrd%k_atm%val(ix) = (ix-1) *k_factor
      end do
      do ix =agrd%nx_atm/2+2,agrd%nx_atm
         agrd%k_atm%val(ix) = (ix-1-agrd%nx_atm)*k_factor
      end do
   end if
  end subroutine set_coord_atm
  subroutine deallocate_atm(grd)
    implicit none
    type(atm_dta),intent(inout) :: grd
    deallocate(grd%x_atm%val);
    deallocate(grd%y_atm%val);
    deallocate(grd%lon_atm%val); deallocate(grd%lat_atm%val);
    deallocate(grd%k_atm%val);
    deallocate(grd%qa_atm%val);    deallocate(grd%pa_atm%val)
    deallocate(grd%ua_atm%val);    deallocate(grd%va_atm%val)
    deallocate(grd%sstm_atm%val);    deallocate(grd%ssta_atm%val)
  end subroutine deallocate_atm
end module mod_atm_dta
