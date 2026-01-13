module mod_wstress
  use ncdf_read
  use run_params
  use run_types
  use calendar_sub
  implicit none
contains
  subroutine ua_to_stress_total(ogrd,oset)
    implicit none
    type(ocn_dta),intent(inout) :: ogrd
    type(ocn_set),intent(in) :: oset
    real(idx):: ua,va,uw,vw
    real(idx),parameter :: rhoa=1.275_idx ! kg/m^3; reference air density
    integer :: ix,iy
     do ix = 0,ogrd%nx_p+1
        do iy = 0,ogrd%ny_p+1
           ua=ogrd%ua_ocn%val(ix,iy)
           va=ogrd%va_ocn%val(ix,iy)
           uw=rhoa*oset%cd_bulk*sqrt(ua*ua+va*va)*ua
           vw=rhoa*oset%cd_bulk*sqrt(ua*ua+va*va)*va
           ogrd%taux_ocn%val(ix,iy)=uw
           ogrd%tauy_ocn%val(ix,iy)=vw
        end do
     end do
   end subroutine ua_to_stress_total
  subroutine ua_to_stress_anm(ogrd,oset,uwmgrd,vwmgrd)
    implicit none
    type(ocn_dta),intent(inout) :: ogrd
    type(ocn_set),intent(in) :: oset
    type(TLL_dta),intent(inout) :: uwmgrd,vwmgrd
    real(idx):: ua,va,uw,vw
    real(idx):: ua_mean,va_mean
    real(idx) :: ut,vt
    real(idx),parameter :: rhoa=1.275_idx ! kg/m^3; reference air density
    integer :: ix,iy
     do ix = 0,ogrd%nx_p+1
        do iy = 0,ogrd%ny_p+1
           ua=ogrd%ua_ocn%val(ix,iy)
           va=ogrd%va_ocn%val(ix,iy)
           ua_mean=uwmgrd%data_now%val(ix,iy)
           va_mean=vwmgrd%data_now%val(ix,iy)
           ut=(ua+ua_mean);vt=(va+va_mean)
           uw=rhoa*oset%cd_bulk*(sqrt(ut*ut+vt*vt)*ut-&
                & sqrt(ua_mean*ua_mean+va_mean*va_mean)*ua_mean)
           vw=rhoa*oset%cd_bulk*(sqrt(ut*ut+vt*vt)*vt-&
                & sqrt(ua_mean*ua_mean+va_mean*va_mean)*va_mean)
           ogrd%taux_ocn%val(ix,iy)=uw
           ogrd%tauy_ocn%val(ix,iy)=vw
        end do
     end do
   end subroutine ua_to_stress_anm
  subroutine ua_to_stress_anm_WWB_GJ20(ogrd,oset,uwmgrd,vwmgrd,model_time,t0_wwb,wt0)
    implicit none
    type(ocn_dta),intent(inout) :: ogrd
    type(ocn_set),intent(in) :: oset
    type(TLL_dta),intent(inout) :: uwmgrd,vwmgrd
    real(idx):: ua,va,uw,vw
    real(idx):: ua_mean,va_mean
    real(idx) :: ut,vt
    real(idx),parameter :: rhoa=1.275_idx ! kg/m^3; reference air density
    integer :: ix,iy
   real(idx),intent(in) :: model_time
   real(idx),intent(inout) :: t0_wwb,wt0 ! in day
   real(idx) :: prob_wwb(1)
   real(idx) :: thres_WWB,G,nino34
   real(idx),parameter :: us0_wwb=6.5_idx,dur_wwb=20.0_idx,interval_wwb=30.0_idx
   real(idx),parameter :: widx_wwb=20.0_idx,widy_wwb=6.0_idx
   real(idx),parameter :: x0_wwb=160.0_idx,y0_wwb=0.0_idx
   G=0.0_idx
   thres_WWB= 0.015_idx+G*nino34
   prob_wwb(1)=1.0_idx
   if (t0_wwb+interval_wwb < model_time) then
      call random_number(prob_wwb(1))  
   end if   
   if (prob_wwb(1) < thres_WWB) then
      t0_wwb=model_time
      wt0=1.0
      write(*,*) t0_wwb
   end if
    do ix = 0,ogrd%nx_p+1
        do iy = 0,ogrd%ny_p+1
           ua=ogrd%ua_ocn%val(ix,iy)+&
         & wt0*us0_wwb*exp(-1.0*((model_time-t0_wwb-dur_wwb)**2)/(dur_wwb)**2) *&
         & exp(-1.0* ((ogrd%lon_p%val(ix)-x0_wwb)**2)/(widx_wwb**2)) *&
         & exp(-1.0* ((ogrd%lat_p%val(iy)-y0_wwb)**2)/(widy_wwb**2))
           va=ogrd%va_ocn%val(ix,iy)
           ua_mean=uwmgrd%data_now%val(ix,iy)
           va_mean=vwmgrd%data_now%val(ix,iy)
           ut=(ua+ua_mean);vt=(va+va_mean)
           uw=rhoa*oset%cd_bulk*(sqrt(ut*ut+vt*vt)*ut-&
                & sqrt(ua_mean*ua_mean+va_mean*va_mean)*ua_mean)
           vw=rhoa*oset%cd_bulk*(sqrt(ut*ut+vt*vt)*vt-&
                & sqrt(ua_mean*ua_mean+va_mean*va_mean)*va_mean)
           ogrd%taux_ocn%val(ix,iy)=uw
           ogrd%tauy_ocn%val(ix,iy)=vw
        end do
     end do
   end subroutine ua_to_stress_anm_WWB_GJ20
end module mod_wstress
