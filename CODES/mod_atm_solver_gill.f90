module mod_atm_solver_gill
  use run_types
  use run_params
  implicit none
contains
  ! Obtain heatsource from SST anomaly and climatology
  subroutine get_heatsource_ZC1987(nx,ny,sst_clm,sst_anm,alpha,q_d_to_nd,q)
    implicit none
    integer,intent(in) :: nx,ny
    real(idx),intent(in) :: sst_clm(nx,ny),sst_anm(nx,ny)
    real(idx),intent(in)  :: alpha,q_d_to_nd
    real(idx),intent(inout) :: q(nx,ny)
    integer :: ix,iy
    do iy = 1,ny
       do ix = 1,nx
          q(ix,iy) = q_d_to_nd*alpha *  sst_anm(ix,iy) * exp((sst_clm(ix,iy)-30.0_idx)/16.7_idx)
       end do
    end do
  end subroutine get_heatsource_ZC1987
  ! Obtain heatsource from SST anomaly and climatology
  subroutine get_heatsource_GJ2022(nx,ny,sst_clm,sst_anm,q_d_to_nd,q)
    implicit none
    integer,intent(in) :: nx,ny
    real(idx),intent(in) :: sst_clm(nx,ny),sst_anm(nx,ny)
    real(idx),intent(inout) :: q(nx,ny)
    real(idx),parameter :: alpha_Q=7.6e-8_idx
    real(idx),parameter :: b_Q=0.5_idx
    real(idx),intent(in)  :: q_d_to_nd
    integer :: ix,iy      
    do iy = 1,ny
       do ix = 1,nx
          q(ix,iy) = q_d_to_nd*alpha_Q*exp(b_Q*(sst_clm(ix,iy)))&
               & *(b_Q*sst_anm(ix,iy)&
               &+(b_Q*(sst_anm(ix,iy)**2))/2.0_idx &
               &+(b_Q*(sst_anm(ix,iy)**3))/6.0_idx &
               &)
       end do
    end do
  end subroutine get_heatsource_GJ2022
  ! Get fourier coefficients of U,V,P from that of Q
  subroutine get_uvpk(nx,ny,k,y,eps,qk,uk,vk,pk)
    implicit none
    integer,intent(in) :: nx,ny
    real(idx),intent(in) :: k(nx),y(ny)
    real(idx),intent(in) :: eps
    complex  :: qk(nx,ny)
    complex,intent(inout) :: uk(nx,ny),vk(nx,ny),pk(nx,ny)
    integer :: ik,iy
    complex  :: A(ny),B(ny),C(ny),D(ny)
    do ik = 1,nx
       A(1)=0.0_idx; B(1)=1.0_idx ; C(1)=0.0_idx; D(1)=0.0_idx
       A(ny)=0.0_idx; B(ny)=1.0_idx ; C(ny)=0.0_idx; D(ny)=0.0_idx
       do iy = 2,ny-1
          D(iy)=  (0.0_idx,1.0_idx)*y(iy)*k(ik)*qk(ik,iy) / (2.0_idx*eps) &
               & -1.0_idx * (qk(ik,(iy+1))-qk(ik,(iy-1)))/(y(iy+1)-y(iy-1))
          A(iy)= 1.0_idx / ((y(iy)-y(iy-1))*0.5_idx*(y(iy+1)-y(iy-1))) !V[iy-1]
          C(iy)= 1.0_idx / ((y(iy+1)-y(iy))*0.5_idx*(y(iy+1)-y(iy-1))) !V[iy+1]
          B(iy) = -1.0_idx / ((y(iy+1)-y(iy))*0.5_idx*(y(iy+1)-y(iy-1))) &
               & - 1.0_idx / ((y(iy)-y(iy-1))*0.5_idx*(y(iy+1)-y(iy-1)))&
               & - 1.0_idx * (y(iy)**2) / 4.0_idx + &
               &      (0.0_idx,1.0_idx) * k(ik) / (2.0_idx*eps)- k(ik)**2- eps**2
       end do
       call solve_tri_cmp(ny,A,B,C,D,Vk(ik,:))
       !       Vk(ik,:)=solve_axb_tri(ny,A,B)
       Uk(ik,1)=0.0_idx
       Uk(ik,ny)=0.0_idx
       Pk(ik,1)=0.0_idx
       Pk(ik,ny)=0.0_idx
       do iy = 2,ny-1
          Uk(ik,iy) = 0.5_idx * eps * y(iy) * Vk(ik,iy)  &
               &    +(0.0_idx,1.0_idx) * k(ik) * (Vk(ik,(iy+1))-Vk(ik,(iy-1))) / (y(iy+1)-y(iy-1)) &
               &    +(0.0_idx,1.0_idx) * k(ik) * qk(ik,iy)
          Uk(ik,iy) = Uk(ik,iy) / (k(ik)**2 + eps**2)
          Pk(ik,iy) = -1.0_idx * eps * qk(ik,iy) &
               & - 1.0_idx * eps * (Vk(ik,(iy+1))-Vk(ik,(iy-1))) / (y(iy+1)-y(iy-1)) &
               & - (0.0_idx,1.0_idx) * k(ik) * 0.5_idx * y(iy) * Vk(ik,iy)
          Pk(ik,iy) = Pk(ik,iy) / (k(ik)**2 + eps**2)
       end do
    end do
  end subroutine get_uvpk
  subroutine return_uvp_fromSST_ZC87(agrd,aset)
    implicit none
    type(atm_dta),intent(inout) :: agrd
    type(atm_set),intent(in) :: aset
!    type(TLL_dta),intent(in) :: sstm_grd,ssta_grd
!    real(idx),intent(in) :: k(1:nx),y(1:ny)
!    real(idx),intent(in) :: sst_clm(1:nx,1:ny),sst_anm(1:nx,1:ny)
!    real(idx),intent(inout) :: u_atm(1:nx,1:ny),v_atm(1:nx,1:ny),p_atm(1:nx,1:ny),q_atm(1:nx,1:ny)
    integer :: iy      
    complex,allocatable :: qk(:,:),pk(:,:),uk(:,:),vk(:,:)
    integer :: nx,ny
    real(idx) :: q_d_to_nd,alpha,eps
    real(idx),parameter :: alpha_0 = 1.0_idx / 300.0_idx ! in [1/K]
    real(idx),parameter :: m0 = 1.0 / 5.0e3 ! in [1/m]
    real(idx),parameter :: cp2=1.0e3 ! in m^2*s^(-2)*K^(-1)
    nx=agrd%nx_atm
    ny=agrd%ny_atm   
    allocate(qk(1:nx,1:ny)) ; allocate(pk(1:nx,1:ny))
    allocate(uk(1:nx,1:ny)) ; allocate(vk(1:nx,1:ny))

    alpha=aset%alpha_gill_atm
    eps = (1.0 / (aset%eps_s_atm_day * 60.0_idx*60.0_idx*24.0_idx))*&
         & (1.0_idx / sqrt(2.0_idx * beta * aset%cp_atm))
    q_d_to_nd  = (g * alpha_0) / ((aset%cp_atm * aset%cp_atm * sqrt(2.0 * beta * aset%cp_atm)) *(m0*cp2))
    ! In [s^3/m^2]
    call get_heatsource_ZC1987(nx,ny,&
         &   agrd%sstm_atm%val(1:nx,1:ny),&
         &   agrd%ssta_atm%val(1:nx,1:ny),&
         &   alpha,q_d_to_nd,agrd%qa_atm%val(1:nx,1:ny))
    !==========================================
    ! Fourier transform of heating function
    !==========================================
    do iy = 1,ny
       qk(1:nx,iy) = fft_for(nx,agrd%qa_atm%val(1:nx,iy))
    end do
    !==========================================!
    ! Obtain fourier coefficients of u,v,p     !
    !==========================================!
    call get_uvpk(nx,ny,agrd%k_atm%val,agrd%y_atm%val,&
         & eps,qk,uk,vk,pk)
    !==========================================!
    !Inverse Fourier transform of p,u,v       !
    !==========================================!
    do iy = 1,ny
       agrd%pa_atm%val(1:nx,iy)=real(fft_back(nx,pk(1:nx,iy)))*&
            & aset%cp_atm*aset%cp_atm*1.3_idx / 100.0_idx ! in [hpa]
       agrd%ua_atm%val(1:nx,iy)=real(fft_back(nx,uk(1:nx,iy)))*&
            & aset%cp_atm ! in [m/s]
       agrd%va_atm%val(1:nx,iy)=real(fft_back(nx,vk(1:nx,iy)))*&
            & aset%cp_atm ! in [m/s]
       agrd%qa_atm%val(1:nx,iy)=agrd%qa_atm%val(1:nx,iy)/q_d_to_nd
    end do
  end subroutine return_uvp_fromSST_ZC87
  subroutine return_uvp_fromSST_GJ22(agrd,aset)
   ! Heating function based on Geng and Jin 2022
   ! "ENSO Diversity Simulated in a Revised Cane-Zebiak Model"
   ! https://doi.org/10.3389/feart.2022.899323
    implicit none
    type(atm_dta),intent(inout) :: agrd
    type(atm_set),intent(in) :: aset
    integer :: iy      
    complex,allocatable :: qk(:,:),pk(:,:),uk(:,:),vk(:,:)
    integer :: nx,ny
    real(idx) :: q_d_to_nd,alpha,eps
    real(idx),parameter :: alpha_0 = 1.0_idx / 300.0_idx ! in [1/K]
    real(idx),parameter :: m0 = 1.0 / 5.0e3 ! in [1/m]
    real(idx),parameter :: cp2=1.0e3 ! in m^2*s^(-2)*K^(-1)
    nx=agrd%nx_atm
    ny=agrd%ny_atm   
    allocate(qk(1:nx,1:ny)) ; allocate(pk(1:nx,1:ny))
    allocate(uk(1:nx,1:ny)) ; allocate(vk(1:nx,1:ny))

    alpha=aset%alpha_gill_atm
    eps = (1.0 / (aset%eps_s_atm_day * 60.0_idx*60.0_idx*24.0_idx))*&
         & (1.0_idx / sqrt(2.0_idx * beta * aset%cp_atm))
    q_d_to_nd  = (g * alpha_0) / ((aset%cp_atm * aset%cp_atm * sqrt(2.0 * beta * aset%cp_atm)) *(m0*cp2))
    call get_heatsource_GJ2022(nx,ny,&
         &   agrd%sstm_atm%val(1:nx,1:ny),&
         &   agrd%ssta_atm%val(1:nx,1:ny),&
         &   q_d_to_nd,agrd%qa_atm%val(1:nx,1:ny))
    !==========================================
    ! Fourier transform of heating function
    !==========================================
    do iy = 1,ny
       qk(1:nx,iy) = fft_for(nx,agrd%qa_atm%val(1:nx,iy))
    end do
    !==========================================!
    ! Obtain fourier coefficients of u,v,p     !
    !==========================================!
    call get_uvpk(nx,ny,agrd%k_atm%val,agrd%y_atm%val,&
         & eps,qk,uk,vk,pk)
    !==========================================!
    !Inverse Fourier transform of p,u,v       !
    !==========================================!
    do iy = 1,ny
       agrd%pa_atm%val(1:nx,iy)=real(fft_back(nx,pk(1:nx,iy)))*&
            & aset%cp_atm*aset%cp_atm*1.3_idx / 100.0_idx ! in [hpa]
       agrd%ua_atm%val(1:nx,iy)=real(fft_back(nx,uk(1:nx,iy)))*&
            & aset%cp_atm ! in [m/s]
       agrd%va_atm%val(1:nx,iy)=real(fft_back(nx,vk(1:nx,iy)))*&
            & aset%cp_atm ! in [m/s]
       agrd%qa_atm%val(1:nx,iy)=agrd%qa_atm%val(1:nx,iy)/q_d_to_nd
    end do
  end subroutine return_uvp_fromSST_GJ22
  subroutine return_uvp_fromSST_ZC87_conv(agrd,aset,uam_grd,vam_grd)
    implicit none
    type(atm_dta),intent(inout) :: agrd
    type(atm_set),intent(in) :: aset
    type(TLL_dta),intent(in) :: uam_grd,vam_grd
    integer :: ix,iy,iter,ix1,ix2
    real(idx),allocatable :: utmp(:,:),vtmp(:,:)
    real(idx),allocatable :: q0(:,:),q(:,:),cnvm(:,:),cnv(:,:)
    complex ,allocatable :: qk(:,:),pk(:,:),uk(:,:),vk(:,:)
    integer :: nx,ny,niter
    real(idx) :: q_d_to_nd,alpha,eps,b1,b2
    real(idx),parameter :: alpha_0 = 1.0_idx / 300.0_idx ! in [1/K]
    real(idx),parameter :: m0 = 1.0 / 5.0e3 ! in [1/m]
    real(idx),parameter :: cp2=1.0e3 ! in m^2*s^(-2)*K^(-1)
    real(idx),parameter :: beta_cnv_atm=1.6e4 ! In [m^2/s^2]
    real(idx) :: cnv_nd_to_d
    niter=3
    nx=agrd%nx_atm
    ny=agrd%ny_atm
    allocate(qk(1:nx,1:ny)) ; allocate(pk(1:nx,1:ny))
    allocate(uk(1:nx,1:ny)) ; allocate(vk(1:nx,1:ny))
    allocate(utmp(1:nx,1:ny));allocate(vtmp(1:nx,1:ny))
    allocate(q(1:nx,1:ny))
    allocate(q0(1:nx,1:ny))
    allocate(cnv(1:nx,1:ny))
    allocate(cnvm(1:nx,1:ny))
    cnv_nd_to_d=sqrt(aset%cp_atm*2.0_idx*beta) ! (1/s); Nondimensional to dimensional
    ! Calculate convergence
    cnvm(1:nx,1:ny)=0.0_idx
    do iy=2,ny-1
       do ix=1,nx
          ix2=ix-1
          if (ix2 < 1) then
             ix2=nx
          end if
          ix1=ix+1
          if (ix1 > nx) then
             ix1=1
          end if
          cnvm(ix,iy)=cnv_nd_to_d*(-1.0_idx*(uam_grd%data_now%val(ix1,iy)-uam_grd%data_now%val(ix2,iy))/&
               & (agrd%x_atm%val(ix1)-agrd%x_atm%val(ix2))&
               & -(vam_grd%data_now%val(ix,iy+1)-vam_grd%data_now%val(ix,iy-1))/&
               & (agrd%y_atm%val(iy+1)-agrd%y_atm%val(iy-1)))/(aset%cp_atm)
       end do
    end do

    alpha=aset%alpha_gill_atm
    eps = (1.0 / (aset%eps_s_atm_day * 60.0_idx*60.0_idx*24.0_idx))*&
         & (1.0_idx / sqrt(2.0_idx * beta * aset%cp_atm))
    q_d_to_nd  = (g * alpha_0) / ((aset%cp_atm * aset%cp_atm * sqrt(2.0 * beta * aset%cp_atm)) *(m0*cp2))
    ! First guess of heat source
    call get_heatsource_ZC1987(nx,ny,&
         &   agrd%sstm_atm%val(1:nx,1:ny),&
         &   agrd%ssta_atm%val(1:nx,1:ny),&
         &   alpha,q_d_to_nd,q0(1:nx,1:ny))
    q(1:nx,1:ny)=q0(1:nx,1:ny)
    do iter =1,niter
       !==========================================
       ! Fourier transform of heating function
       !==========================================
       do iy = 1,ny
          qk(1:nx,iy) = fft_for(nx,q(1:nx,iy))
       end do
       !==========================================!
       ! Obtain fourier coefficients of u,v,p     !
       !==========================================!
       call get_uvpk(nx,ny,agrd%k_atm%val,agrd%y_atm%val,&
            & eps,qk,uk,vk,pk)
       !==========================================!
       !Inverse Fourier transform of p,u,v       !
       !==========================================!
       do iy = 1,ny
          utmp(1:nx,iy)=real(fft_back(nx,uk(1:nx,iy)))*&
               & aset%cp_atm ! in [m/s]
          vtmp(1:nx,iy)=real(fft_back(nx,vk(1:nx,iy)))*&
               & aset%cp_atm ! in [m/s]
       end do

       ! Calculate convergence
       cnv(1:nx,1:ny)=0.0_idx
       do iy=2,ny-1
          do ix=1,nx
             ix1=ix+1
             if (ix1 > nx) then
                ix1=1
             end if
             ix2=ix-1
             if (ix2 < 1) then
                ix2=nx
             end if
             cnv(ix,iy)=cnv_nd_to_d*(-1.0_idx*(utmp(ix1,iy)-utmp(ix2,iy))/&
                  & (agrd%x_atm%val(ix1)-agrd%x_atm%val(ix2))&
                  & -(vtmp(ix,iy+1)-vtmp(ix,iy-1))/&
                  & (agrd%y_atm%val(iy+1)-agrd%y_atm%val(iy-1)))/(aset%cp_atm)
          end do
       end do
       ! Calculate heating
       do iy=1,ny
          do ix=1,nx
             if (cnv(ix,iy)+cnvm(ix,iy)>0) then
                b1=cnv(ix,iy)+cnvm(ix,iy)
             else
                b1=0.0_idx
             end if
             if (cnvm(ix,iy)>0) then
                b2=cnvm(ix,iy)
             else
                b2=0.0_idx
             end if
             if (q0(ix,iy).ne. 0.0_idx) then
             q(ix,iy)=q0(ix,iy)+q_d_to_nd*beta_cnv_atm*&
                  & (b1-b2)             
             else
               q(ix,iy)=0.0_idx
             end if
          end do
       end do
    end do
    do iy = 1,ny
       agrd%pa_atm%val(1:nx,iy)=real(fft_back(nx,pk(1:nx,iy)))*&
            & aset%cp_atm*aset%cp_atm*1.3_idx / 100.0_idx ! in [hpa]
       agrd%ua_atm%val(1:nx,iy)=utmp(1:nx,iy)
       agrd%va_atm%val(1:nx,iy)=vtmp(1:nx,iy)
       agrd%qa_atm%val(1:nx,iy)=q(1:nx,iy)/q_d_to_nd
    end do
    deallocate(utmp);deallocate(vtmp)
    deallocate(qk);deallocate(pk)
    deallocate(uk);deallocate(vk)
    deallocate(q0);deallocate(q)
    deallocate(cnv);deallocate(cnvm)
  end subroutine return_uvp_fromSST_ZC87_conv

  ! Solver
  subroutine solve_tri_cmp(N,A_in,B_in,C_in,D_in,U)
    ! Solve A * u(i-1) + B * u(i) + C * u(i+1) = D
    implicit none
    integer,intent(in) :: N
    complex ,intent(in) :: A_in(N),B_in(N),C_in(N),D_in(N)
    complex ,intent(inout) :: U(N)
    complex  :: A(N),B(N),C(N),D(N)
    integer :: i
    complex  :: m
    A = A_in ; B = B_in ; C = C_in; D = D_in
    ! initialize array
    U(1:N) = 0.0_idx
    do i=2,N
       m = A(i) / B(i-1)
       B(i) = B(i) - m * C(i-1)
       D(i) = D(i) - m * D(i-1)
       U(i) = 0.0_idx
    end do
    U(N) = D(N) / B(N)
    do i = N-1,1,-1
       U(i) = (D(i)-C(i)*U(i+1)) / B(i)
    end do
  end subroutine solve_tri_cmp

  ! FFT forward routine (using FFTPACK)
  function fft_for(n,array) result(c)
    implicit none
    integer,intent(in) :: n
    real(idx),intent(in) :: array(n)
    real(idx),allocatable :: wsave(:),work(:)
    integer :: lensav,ier,inc,lenc,lenwrk
    complex :: c(n)
    integer :: i
    lensav=2*n+n+4
    allocate(wsave(lensav))
    call cfft1i(n,wsave,lensav,ier)
    do i = 1,n
       c(i)=array(i)
    end do
    inc=1
    lenc=inc*(N-1)+1
    lenwrk=3*N
    allocate(work(lenwrk))
    call cfft1f ( n, inc, c, lenc, wsave, lensav, work, lenwrk, ier)
    deallocate(wsave)
    deallocate(work)
  end function fft_for
  ! FFT backward routine (using FFTPACK)
  function fft_back(n,c) result(array)
    implicit none
    integer,intent(in) :: n
    complex,intent(in) :: c(n)
    real(idx),allocatable :: wsave(:),work(:)
    integer :: lensav,ier,inc,lenc,lenwrk
    real(idx) :: array(n)
    integer :: i
    lensav=2*N+N+4
    allocate(wsave(lensav))
    call cfft1i (n,wsave,lensav,ier)
    inc=1
    lenc=inc*(N-1)+1
    lenwrk=3*N
    allocate(work(lenwrk))
    call cfft1b(n,inc, c, lenc, wsave, lensav, work, lenwrk, ier)
    do i =1,n
       array(i)=real(c(i))
    end do
    deallocate(wsave)
    deallocate(work)
  end function fft_back

!   function fft_for(n,array) result(c)
!     implicit none
!     integer,intent(in) :: n
!     real(idx),intent(in) :: array(n)
!     real(idx),allocatable :: wsave(:),work(:)
!     integer :: lensav,ier,inc,lenc,lenwrk
!     complex  :: c(n)
!     integer :: i
!     lensav=2*N+N+4
!     allocate(wsave(lensav))
!     call cfft1i(n,wsave,lensav,ier)
!     do i = 1,n
!        c(i)=array(i)
!     end do
!     inc=1
!     lenc=inc*(N-1)+1
!     lenwrk=3*N
!     allocate(work(lenwrk))
!     call cfft1f ( n, inc, c, lenc, wsave, lensav, work, lenwrk, ier)
!     deallocate(wsave)
!     deallocate(work)
!   end function fft_for
!   ! FFT backward routine (using FFTPACK)
!   function fft_back(n,c) result(array)
!     implicit none
!     integer,intent(in) :: n
!     complex ,intent(in) :: c(n)
!     real(idx),allocatable :: wsave(:),work(:)
!     integer :: lensav,ier,inc,lenc,lenwrk
!     real(idx) :: array(n)
!     integer :: i
!     lensav=2*n+n+4
!     allocate(wsave(lensav))
!     call cfft1i (n,wsave,lensav,ier)
!     inc=1
!     lenc=inc*(n-1)+1
!     lenwrk=3*n
!     allocate(work(lenwrk))
!     call cfft1b(n,inc, c, lenc, wsave, lensav, work, lenwrk, ier)
!     do i =1,n
! !      array(i)=real(c(i))
!        array(i) = abs(c(i))*sign(1.0_idx,real(c(i),idx))
! !       array(i) = abs(c(i))*sign (1.0,real(c(i),idx))
!     end do
!     deallocate(wsave)
!     deallocate(work)
!   end function fft_back
end module mod_atm_solver_gill


