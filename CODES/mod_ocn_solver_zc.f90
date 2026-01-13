module mod_ocn_solver_zc
  use run_types
  use run_params
  implicit none
  character(len=maxlen) :: wbc_flag_p,ebc_flag_p,nbc_flag_p,sbc_flag_p
  character(len=maxlen) :: wbc_flag_u,ebc_flag_u,nbc_flag_u,sbc_flag_u
  character(len=maxlen) :: wbc_flag_v,ebc_flag_v,nbc_flag_v,sbc_flag_v
contains
  subroutine solve_ekman_ocn(ogrd,oset)
    implicit none
    type(ocn_dta),intent(inout) :: ogrd
    type(ocn_set),intent(in) :: oset
    integer :: ix,iy
    real(idx) :: uw,vw,H1,H2,eps_s_ocn,y
    real(idx),allocatable :: u_ek_p(:,:),v_ek_p(:,:)
    eps_s_ocn=1.0/(oset%eps_s_ocn_day*day_to_sec)
    H1=oset%H1
    H2=oset%H2
    allocate(u_ek_p(0:ogrd%nx_p+1,0:ogrd%ny_p+1))
    allocate(v_ek_p(0:ogrd%nx_p+1,0:ogrd%ny_p+1))
    do iy=0,ogrd%ny_p+1
       do ix =0,ogrd%nx_p+1
          y=lat_to_dis*ogrd%lat_p%val(iy)
          uw=ogrd%taux_ocn%val(ix,iy)/(rho0_o*H1)
          vw=ogrd%tauy_ocn%val(ix,iy)/(rho0_o*H1)
          u_ek_p(ix,iy)=(eps_s_ocn*uw+ogrd%f%val(ix,iy)*vw) / &
               & (((ogrd%f%val(ix,iy))**2+eps_s_ocn**2))
          v_ek_p(ix,iy)=(eps_s_ocn*vw-ogrd%f%val(ix,iy)*uw) / &
               & (((ogrd%f%val(ix,iy))**2+eps_s_ocn**2))
       end do
    end do
    do iy=1,ogrd%ny_p+1
       ix=0
       ogrd%v_ek%val(ix,iy)=0.5_idx*(v_ek_p(ix,iy)+v_ek_p(ix,iy-1))
       do ix=1,ogrd%nx_p+1
          ogrd%u_ek%val(ix,iy)=0.5_idx*(u_ek_p(ix,iy)+u_ek_p(ix-1,iy))
          ogrd%v_ek%val(ix,iy)=0.5_idx*(v_ek_p(ix,iy)+v_ek_p(ix,iy-1))
       end do
    end do
    iy=0
    do ix=1,ogrd%nx_p+1
       ogrd%u_ek%val(ix,iy)=0.5_idx*(u_ek_p(ix,iy)+u_ek_p(ix-1,iy))
    end do
    do iy=1,ogrd%ny_p
       do ix=1,ogrd%nx_p
          ogrd%w_ek%val(ix,iy)=((ogrd%u_ek%val(ix+1,iy)-ogrd%u_ek%val(ix,iy))/&
               & (ogrd%x_u%val(ix+1,iy)-ogrd%x_u%val(ix,iy))  &
               & +(ogrd%v_ek%val(ix,iy+1)-ogrd%v_ek%val(ix,iy))/&
               & (ogrd%y_v%val(ix,iy+1)-ogrd%y_v%val(ix,iy)))*H1
       end do
    end do
    deallocate(u_ek_p)
    deallocate(v_ek_p)
  end subroutine solve_ekman_ocn
  subroutine set_bc_p(nx,ny,p,wbc_flag,ebc_flag,nbc_flag,sbc_flag)
    implicit none
    integer,intent(in) :: nx,ny
    real(idx),intent(inout) :: p(0:nx+1,0:ny+1)
    character(len=*),intent(in) :: wbc_flag,ebc_flag,nbc_flag,sbc_flag
    integer :: ix,iy
    !=================================
    ! Western boundary condition
    !=================================
    select case(wbc_flag)
    case ("Clo","CLO","clo") ! Closed boundary condition
       do iy = 1,ny
          p(0,iy) = 0.0_idx
       end do
    case ("Gra","GRA","gra") ! Gradient boundary condition
       do iy = 1,ny
          p(0,iy) = p(1,iy)
       end do
    end select
    !=================================
    ! Eastern boundary condition
    !=================================    
    select case(ebc_flag)
    case ("Clo","CLO","clo") ! Closed boundary condition
       do iy = 1,ny
          p(nx+1,iy) = 0.0_idx
       end do
    case ("Gra","GRA","gra") ! Gradient boundary condition
       do iy = 1,ny
          p(nx+1,iy) = p(nx,iy)
       end do
    end select
    !=================================
    ! Southern boundary condition
    !=================================
    select case(sbc_flag)
    case ("Clo","CLO","clo") ! Closed boundary condition
       do ix = 1,nx
          p(ix,0) = 0.0_idx
       end do
    case ("Gra","GRA","gra") ! Gradient boundary condition
       do ix = 1,nx
          p(ix,0) = p(ix,1)
       end do
    end select
    !=================================
    ! Northern boundary condition
    !=================================
    select case(nbc_flag)
    case ("Clo","CLO","clo") ! Closed boundary condition
       do ix = 1,nx
          p(ix,ny+1) = 0.0_idx
       end do
    case ("Gra","GRA","gra") ! Gradient boundary condition
       do ix = 1,nx
          p(ix,ny+1) = p(ix,ny)
       end do
    end select
  end subroutine set_bc_p
  subroutine set_bc_u(nx,ny,u,wbc_flag,ebc_flag,nbc_flag,sbc_flag,slip_ind)
    implicit none
    integer,intent(in) :: nx,ny
    real(idx),intent(inout) :: u(1:nx+1,0:ny+1)
    character(len=*),intent(in) :: wbc_flag,ebc_flag,nbc_flag,sbc_flag
    real(idx),intent(in) :: slip_ind
    real(idx) :: gamma2
    integer :: ix,iy
    gamma2 = 1.0_idx - 2*slip_ind
    ! gamma2=1 for slip_ind=0 (du/dx=0),
    ! gamma2=-1 for slip_ind=1 (u=0)
    !=================================
    ! Western boundary condition
    !=================================
    select case(wbc_flag)
    case ("Clo","CLO","clo") ! Closed boundary condition
       do iy = 1,ny
          u(1,iy) = 0.0_idx
       end do
    case ("Gra","GRA","gra") ! Gradient boundary condition
       do iy = 1,ny
          u(1,iy) = u(2,iy)
       end do
    end select
    !=================================
    ! Eastern boundary condition
    !=================================    
    select case(ebc_flag)
    case ("Clo","CLO","clo") ! Closed boundary condition
       do iy = 1,ny
          u(nx+1,iy) = 0.0_idx
       end do
    case ("Gra","GRA","gra") ! Gradient boundary condition
       do iy = 1,ny
          u(nx+1,iy) = u(nx,iy)
       end do
    end select
    !=================================
    ! Southern boundary condition
    !=================================
    select case(sbc_flag)
    case ("Clo","CLO","clo") ! Closed boundary condition
       do ix = 2,nx
          u(ix,0) = gamma2 * u(ix,1)
       end do
    case ("Gra","GRA","gra") ! Gradient boundary condition
       do ix = 2,nx
          u(ix,0) = u(ix,1)
       end do
    end select
    !=================================
    ! Northern boundary condition
    !=================================
    select case(nbc_flag)
    case ("Clo","CLO","clo") ! Closed boundary condition
       do ix = 2,nx
          u(ix,ny+1) = gamma2 * u(ix,ny)
       end do
    case ("Gra","GRA","gra") ! Gradient boundary condition
       do ix = 2,nx
          u(ix,ny+1) = u(ix,ny)
       end do
    end select
  end subroutine set_bc_u
  subroutine set_bc_v(nx,ny,v_next,wbc_flag,ebc_flag,nbc_flag,sbc_flag,slip_ind)
    implicit none
    integer,intent(in) :: nx,ny
!    real(idx),intent(in) :: v(0:nx+1,1:ny+1),v_past(0:nx+1,1:ny+1)
    real(idx),intent(inout) :: v_next(0:nx+1,1:ny+1)
    character(len=*),intent(in) :: wbc_flag,ebc_flag,nbc_flag,sbc_flag
    real(idx),intent(in) :: slip_ind
    real(idx) :: gamma2
    integer :: ix,iy
    gamma2 = 1.0_idx - 2*slip_ind 
    ! gamma2=1 for slip_ind=0 (du/dx=0),
    ! gamma2=-1 for slip_ind=1 (u=0)
    !=================================
    ! Western boundary condition
    !=================================
    select case(wbc_flag)
    case ("Clo","CLO","clo") ! Closed boundary condition
       do iy = 2,ny
          v_next(0,iy) = gamma2 * v_next(1,iy)
       end do
    case ("Gra","GRA","gra") ! Gradient boundary condition
       do iy = 2,ny
          v_next(0,iy) = v_next(1,iy)
       end do
    end select
    !=================================
    ! Eastern boundary condition
    !=================================    
    select case(ebc_flag)
    case ("Clo","CLO","clo") ! Closed boundary condition
       do iy = 2,ny
          v_next(nx+1,iy) = gamma2 * v_next(nx,iy)
       end do
    case ("Gra","GRA","gra") ! Gradient boundary condition
       do iy = 2,ny
          v_next(nx+1,iy) = v_next(nx,iy)
       end do
    end select
    !=================================
    ! Southern boundary condition
    !=================================
    select case(sbc_flag)
    case ("Clo","CLO","clo") ! Closed boundary condition
       do ix = 1,nx
          v_next(ix,1) = 0.0_idx
       end do
    case ("Gra","GRA","gra") ! Gradient boundary condition
       do ix = 1,nx
          v_next(ix,1) = v_next(ix,2)
       end do
    end select
    !=================================
    ! Northern boundary condition
    !=================================
    select case(nbc_flag)
    case ("Clo","CLO","clo") ! Closed bondary condition
       do ix = 1,nx
          v_next(ix,ny+1) = 0.0_idx
       end do
    case ("Gra","GRA","gra") ! Gradient boundary condition
       do ix = 1,nx
          v_next(ix,ny+1) = v_next(ix,ny)
       end do
    end select
  end subroutine set_bc_v
 subroutine get_rhs_u(ogrd,oset,u,v,h,rhs_u)
   type(ocn_dta),intent(in) :: ogrd
    type(ocn_set),intent(in) :: oset
   real(idx),intent(in) :: u(1:ogrd%nx_p+1,0:ogrd%ny_p+1)
   real(idx),intent(in) :: v(0:ogrd%nx_p+1,1:ogrd%ny_p+1)
   real(idx),intent(in) :: h(0:ogrd%nx_p+1,0:ogrd%ny_p+1)
   real(idx),intent(inout) :: rhs_u(1:ogrd%nx_p+1,0:ogrd%ny_p+1)   
   integer :: ix,iy
   real(idx) :: dudx_e,dudx_w,dudy_n,dudy_s
   real(idx) :: drag_u,corx,pgrdx,tx,diffu
   real(idx) :: H1,H2,gprime
   real(idx) :: sigma_n,sigma_s,sigma_w,sigma_e
   rhs_u=0.0_idx
   H1=oset%H1; H2=oset%H2
   gprime=(oset%cp_ocn*oset%cp_ocn)/(H1+H2)
   do iy = 1,ogrd%ny_p
      do ix = 2,ogrd%nx_p
         ! Coriolis force
         corx= 0.125_idx * ((ogrd%f%val(ix,iy)+ogrd%f%val(ix,iy+1))*&
              & (v(ix,iy+1)+v(ix-1,iy+1)) + &
              & (ogrd%f%val(ix,iy-1)+ogrd%f%val(ix,iy)) * (v(ix,iy)+v(ix-1,iy)))
         !Pressure gradient force
         pgrdx= -1.0_idx*gprime*(h(ix,iy)-h(ix-1,iy)) &
              & / (ogrd%x_p%val(ix,iy)-ogrd%x_p%val(ix-1,iy))
         ! Wind forcing
         tx= 0.5_idx*(ogrd%taux_ocn%val(ix-1,iy) + ogrd%taux_ocn%val(ix,iy)) * &
              (1.0_idx/((H1+H2)* rho0_o))            
          ! Drag force
          drag_u = -1.0_idx * u(ix,iy) * (ogrd%damp_u%val(ix,iy))
          ! U-viscosity
          dudx_e = (u(ix+1,iy) - u(ix,iy)) / &
               & (ogrd%x_u%val(ix+1,iy)-ogrd%x_u%val(ix,iy))
          dudx_w = (u(ix,iy) - u(ix-1,iy)) / &
               & (ogrd%x_u%val(ix,iy)-ogrd%x_u%val(ix-1,iy))
          dudy_n = ogrd%mask_phi_u%val(ix,iy+1) * &
               & (u(ix,iy+1)-u(ix,iy)) / &
               & (ogrd%y_u%val(ix,iy+1)-ogrd%y_u%val(ix,iy))
          dudy_s = ogrd%mask_phi_u%val(ix,iy) * &
               & (u(ix,iy)-u(ix,iy-1)) / &
               & (ogrd%y_u%val(ix,iy)-ogrd%y_u%val(ix,iy-1))
          sigma_n = 0.25_idx * (ogrd%visc_2D%val(ix-1,iy)+ogrd%visc_2D%val(ix,iy)+&
               & ogrd%visc_2D%val(ix-1,iy+1)+ogrd%visc_2D%val(ix,iy+1)) * dudy_n
          sigma_s = 0.25_idx * (ogrd%visc_2D%val(ix-1,iy-1)+ogrd%visc_2D%val(ix,iy-1)+&
               & ogrd%visc_2D%val(ix-1,iy)+ogrd%visc_2D%val(ix,iy)) * dudy_s
          sigma_e = ogrd%visc_2D%val(ix,iy) * dudx_e
          sigma_w = ogrd%visc_2D%val(ix-1,iy) * dudx_w
          diffu = (sigma_e-sigma_w) / (ogrd%x_p%val(ix,iy)-ogrd%x_p%val(ix-1,iy))&
               & + (sigma_n - sigma_s) / (ogrd%y_v%val(ix,iy+1)-ogrd%y_v%val(ix,iy))
          rhs_u(ix,iy) = corx + pgrdx + tx+drag_u +diffu
      end do
   end do
 end subroutine get_rhs_u
  subroutine get_rhs_v(ogrd,oset,u,v,h,rhs_v)
   type(ocn_dta),intent(in) :: ogrd
    type(ocn_set),intent(in) :: oset
   real(idx),intent(in) :: u(1:ogrd%nx_p+1,0:ogrd%ny_p+1)
   real(idx),intent(in) :: v(0:ogrd%nx_p+1,1:ogrd%ny_p+1)
   real(idx),intent(in) :: h(0:ogrd%nx_p+1,0:ogrd%ny_p+1)
   real(idx),intent(inout) :: rhs_v(0:ogrd%nx_p+1,1:ogrd%ny_p+1)   
   integer :: ix,iy
   real(idx) :: drag_v,cory,pgrdy,ty,diffv
   real(idx) :: sigma_n,sigma_s,sigma_w,sigma_e
   real(idx) :: dvdx_e,dvdx_w,dvdy_n,dvdy_s
   real(idx) :: H1,H2,gprime
   rhs_v=0.0_idx
   H1=oset%H1; H2=oset%H2
   gprime=(oset%cp_ocn*oset%cp_ocn)/(H1+H2)
   
   do iy = 2,ogrd%ny_p
      do ix = 1,ogrd%nx_p
         cory = -0.25_idx *  (ogrd%f%val(ix,iy) * u(ix,iy)+ ogrd%f%val(ix,iy)*&
              & u(ix+1,iy)+ogrd%f%val(ix,iy-1) * u(ix,iy-1)+ ogrd%f%val(ix,iy-1)* u(ix+1,iy-1))! Coriolis force (ENG conserve) !
          ! Pressure gradient
          pgrdy=-1.0_idx*gprime*(h(ix,iy)-h(ix,iy-1)) / &
               & (ogrd%y_p%val(ix,iy)-ogrd%y_p%val(ix,iy-1)) 
          ! Wind forcing
          ty=0.5_idx*(ogrd%tauy_ocn%val(ix,iy-1) + ogrd%tauy_ocn%val(ix,iy)) * &
               & (1.0_idx/((H1+H2)* rho0_o))
          ! Drag
          drag_v = -1.0_idx * v(ix,iy) * (ogrd%damp_v%val(ix,iy))
          ! Viscosity
          dvdx_e = ogrd%mask_phi_v%val(ix+1,iy) * (v(ix+1,iy) -v(ix,iy)) &
               & / (ogrd%x_v%val(ix+1,iy)-ogrd%x_v%val(ix,iy))
          dvdx_w = ogrd%mask_phi_v%val(ix,iy) * (v(ix,iy) -  v(ix-1,iy)) / &
          & (ogrd%x_v%val(ix,iy)-ogrd%x_v%val(ix-1,iy))
          dvdy_n = (v(ix,iy+1) - v(ix,iy)) / &
               & (ogrd%y_v%val(ix,iy+1)-ogrd%y_v%val(ix,iy))
          dvdy_s = (v(ix,iy) - v(ix,iy-1)) / &
               & (ogrd%y_v%val(ix,iy)-ogrd%y_v%val(ix,iy-1))
          sigma_n = ogrd%visc_2D%val(ix,iy)*dvdy_n
          sigma_s = ogrd%visc_2D%val(ix,iy-1)*dvdy_s
          sigma_e = 0.25_idx * (ogrd%visc_2D%val(ix,iy-1)+ogrd%visc_2D%val(ix+1,iy-1)&
               & +ogrd%visc_2D%val(ix,iy)+ogrd%visc_2D%val(ix+1,iy))*dvdx_e
          sigma_w = 0.25_idx * (ogrd%visc_2D%val(ix-1,iy-1)+ogrd%visc_2D%val(ix,iy-1)&
               & +ogrd%visc_2D%val(ix-1,iy)+ogrd%visc_2D%val(ix,iy))*dvdx_w
          diffv = (sigma_e-sigma_w) / (ogrd%x_u%val(ix+1,iy)-ogrd%x_u%val(ix,iy)) + &
               & (sigma_n - sigma_s) / (ogrd%y_p%val(ix,iy)-ogrd%y_p%val(ix,iy-1))
          rhs_v(ix,iy)=drag_v + cory + pgrdy + ty +diffv
       end do
    end do
  end subroutine get_rhs_v
  subroutine get_rhs_p(ogrd,oset,u,v,h,rhs_p)
   type(ocn_dta),intent(in) :: ogrd
    type(ocn_set),intent(in) :: oset
   real(idx),intent(in) :: u(1:ogrd%nx_p+1,0:ogrd%ny_p+1)
   real(idx),intent(in) :: v(0:ogrd%nx_p+1,1:ogrd%ny_p+1)
   real(idx),intent(in) :: h(0:ogrd%nx_p+1,0:ogrd%ny_p+1)
   real(idx),intent(inout) :: rhs_p(0:ogrd%nx_p+1,0:ogrd%ny_p+1)   
   integer :: ix,iy
    real(idx) :: dudx,dvdy,drag_p
   real(idx) :: H1,H2
   rhs_p=0.0_idx
   H1=oset%H1; H2=oset%H2
      do iy = 1,ogrd%ny_p
       do ix = 1,ogrd%nx_p
          !=================================
          ! p equation
          !=================================
          dudx=-1.0_idx*(H1+H2)*(u(ix+1,iy)-u(ix,iy)) &
               & / (ogrd%x_u%val(ix+1,iy)-ogrd%x_u%val(ix,iy))
          dvdy=-1.0_idx*(H1+H2)*(v(ix,iy+1)-v(ix,iy)) &
               & / (ogrd%y_v%val(ix,iy+1)-ogrd%y_v%val(ix,iy))
          drag_p = -1.0_idx * h(ix,iy) * (ogrd%damp_p%val(ix,iy))
          rhs_p(ix,iy) = drag_p + dudx + dvdy
       end do
    end do
  end subroutine get_rhs_p
  subroutine solve_rg_vgeo_ocn(ogrd,oset,dt)
    implicit none
    type(ocn_dta),intent(inout) :: ogrd
    type(ocn_set),intent(in) :: oset
    real(idx),intent(in) :: dt
    integer :: ix,iy
    real(idx),allocatable :: rhs1_u(:,:),rhs1_v(:,:),rhs1_h(:,:)
    allocate(rhs1_h(0:ogrd%nx_p+1,0:ogrd%ny_p+1));
    allocate(rhs1_u(1:ogrd%nx_p+1,0:ogrd%ny_p+1));
    allocate(rhs1_v(0:ogrd%nx_p+1,1:ogrd%ny_p+1));
    call get_rhs_u(ogrd,oset,ogrd%u_sw%val,ogrd%v_sw%val,ogrd%h_sw%val,rhs1_u)
    ! ! u    x:2~nx_p # 1:WB nx+1:EB
    ! !      y:1~ny # 0:SB ny+1:NB
    do iy = 1,ogrd%ny_p
       do ix = 2,ogrd%nx_p
          ! Update u
          ogrd%u_sw_next%val(ix,iy)=ogrd%u_sw_past%val(ix,iy) +2.0_idx*dt * ogrd%mask_u%val(ix,iy) * rhs1_u(ix,iy)
       end do
    end do
    call get_rhs_v(ogrd,oset,ogrd%u_sw%val,ogrd%v_sw%val,ogrd%h_sw%val,rhs1_v)
    ! v    x:1~nx # 0:WB nx+1:EB
    !      y:2~ny # 1:SB ny+1:NB
    do iy = 2,ogrd%ny_p
       do ix = 1,ogrd%nx_p
          ogrd%v_sw_next%val(ix,iy)=ogrd%v_sw_past%val(ix,iy)+2.0_idx*dt*ogrd%mask_v%val(ix,iy)*rhs1_v(ix,iy)
       end do
    end do
    ! P
    call get_rhs_p(ogrd,oset,ogrd%u_sw%val,ogrd%v_sw%val,ogrd%h_sw%val,rhs1_h)
    do iy = 1,ogrd%ny_p
       do ix = 1,ogrd%nx_p
          ogrd%h_sw_next%val(ix,iy)=ogrd%h_sw_past%val(ix,iy)+2.0_idx*dt*ogrd%mask_p%val(ix,iy)*rhs1_h(ix,iy)
       end do
    end do
    ! Asselin filter
    do iy = 1,ogrd%ny_p
      do ix =1,ogrd%nx_p
         ogrd%h_sw%val(ix,iy) = ogrd%h_sw%val(ix,iy)  &
              & +0.5_idx * 0.2_idx * ogrd%mask_p%val(ix,iy)&
              & *(ogrd%h_sw_next%val(ix,iy) + ogrd%h_sw_past%val(ix,iy)&
              & -2.0_idx*ogrd%h_sw%val(ix,iy))
      end do
   end do
   do iy = 1,ogrd%ny_p
      do ix =2,ogrd%nx_p
         ogrd%u_sw%val(ix,iy) = ogrd%u_sw%val(ix,iy)  &
              & +0.5_idx * 0.2_idx * ogrd%mask_u%val(ix,iy)&
              & *(ogrd%u_sw_next%val(ix,iy) + ogrd%u_sw_past%val(ix,iy)&
              & -2.0_idx*ogrd%u_sw%val(ix,iy))
      end do
   end do
   do iy = 2,ogrd%ny_p
      do ix =1,ogrd%nx_p
         ogrd%v_sw%val(ix,iy) = ogrd%v_sw%val(ix,iy)  &
              & +0.5_idx * 0.2_idx * ogrd%mask_v%val(ix,iy)&
              & *(ogrd%v_sw_next%val(ix,iy) + ogrd%v_sw_past%val(ix,iy)&
              & -2.0_idx*ogrd%v_sw%val(ix,iy))
      end do
   end do
   
   call set_bc_p(ogrd%nx_p,ogrd%ny_p,&
        & ogrd%h_sw_next%val(0:ogrd%nx_p+1,0:ogrd%ny_p+1),&
        & oset%wbc_flag_p,oset%ebc_flag_p,oset%nbc_flag_p,oset%sbc_flag_p)
   call set_bc_u(ogrd%nx_p,ogrd%ny_p,&
        & ogrd%u_sw_next%val(1:ogrd%nx_p+1,0:ogrd%ny_p+1),&
        & oset%wbc_flag_u,oset%ebc_flag_u,oset%nbc_flag_u,oset%sbc_flag_u,&
        & oset%slip_ind)
   call set_bc_v(ogrd%nx_p,ogrd%ny_p,&
        & ogrd%v_sw_next%val(0:ogrd%nx_p+1,1:ogrd%ny_p+1),&
        & oset%wbc_flag_v,oset%ebc_flag_v,oset%nbc_flag_v,oset%sbc_flag_v,&
        & oset%slip_ind)
   ! Update
   ogrd%u_sw_past%val(1:ogrd%nx_p+1,0:ogrd%ny_p+1) = ogrd%u_sw%val(1:ogrd%nx_p+1,0:ogrd%ny_p+1)
   ogrd%u_sw%val(1:ogrd%nx_p+1,0:ogrd%ny_p+1) = ogrd%u_sw_next%val(1:ogrd%nx_p+1,0:ogrd%ny_p+1)
   ogrd%v_sw_past%val(0:ogrd%nx_p+1,1:ogrd%ny_p+1) = ogrd%v_sw%val(0:ogrd%nx_p+1,1:ogrd%ny_p+1)
   ogrd%v_sw%val(0:ogrd%nx_p+1,1:ogrd%ny_p+1) = ogrd%v_sw_next%val(0:ogrd%nx_p+1,1:ogrd%ny_p+1)
   ogrd%h_sw_past%val(0:ogrd%nx_p+1,0:ogrd%ny_p+1) = ogrd%h_sw%val(0:ogrd%nx_p+1,0:ogrd%ny_p+1)
   ogrd%h_sw%val(0:ogrd%nx_p+1,0:ogrd%ny_p+1) = ogrd%h_sw_next%val(0:ogrd%nx_p+1,0:ogrd%ny_p+1)
   deallocate(rhs1_h);deallocate(rhs1_u);deallocate(rhs1_v)
 end subroutine solve_rg_vgeo_ocn
  subroutine solve_rg_vgeo_ocn_rk2(ogrd,oset,dt)
    implicit none
    type(ocn_dta),intent(inout) :: ogrd
    type(ocn_set),intent(in) :: oset
    real(idx),intent(in) :: dt
    integer :: ix,iy
    real(idx),allocatable :: rhs1_u(:,:),rhs1_v(:,:),rhs1_h(:,:)
    real(idx),allocatable :: rhs2_u(:,:),rhs2_v(:,:),rhs2_h(:,:)
    real(idx),allocatable :: u_tmp(:,:),v_tmp(:,:),h_tmp(:,:)
    allocate(rhs1_h(0:ogrd%nx_p+1,0:ogrd%ny_p+1));
    allocate(rhs2_h(0:ogrd%nx_p+1,0:ogrd%ny_p+1));
    allocate(h_tmp(0:ogrd%nx_p+1,0:ogrd%ny_p+1))
    allocate(rhs1_u(1:ogrd%nx_p+1,0:ogrd%ny_p+1));
    allocate(rhs2_u(1:ogrd%nx_p+1,0:ogrd%ny_p+1));
    allocate(u_tmp(1:ogrd%nx_p+1,0:ogrd%ny_p+1))
    allocate(rhs1_v(0:ogrd%nx_p+1,1:ogrd%ny_p+1));
    allocate(rhs2_v(0:ogrd%nx_p+1,1:ogrd%ny_p+1));
    allocate(v_tmp(0:ogrd%nx_p+1,1:ogrd%ny_p+1))
    call get_rhs_u(ogrd,oset,ogrd%u_sw%val,ogrd%v_sw%val,ogrd%h_sw%val,rhs1_u)
    ! ! u    x:2~nx_p # 1:WB nx+1:EB
    ! !      y:1~ny # 0:SB ny+1:NB
    do iy = 1,ogrd%ny_p
       do ix = 2,ogrd%nx_p
          ! Update u
          u_tmp(ix,iy)=ogrd%u_sw%val(ix,iy) +dt * ogrd%mask_u%val(ix,iy) * rhs1_u(ix,iy)
       end do
    end do
    call get_rhs_v(ogrd,oset,ogrd%u_sw%val,ogrd%v_sw%val,ogrd%h_sw%val,rhs1_v)
    ! v    x:1~nx # 0:WB nx+1:EB
    !      y:2~ny # 1:SB ny+1:NB
    do iy = 2,ogrd%ny_p
       do ix = 1,ogrd%nx_p
          v_tmp(ix,iy)=ogrd%v_sw%val(ix,iy)+ dt*ogrd%mask_v%val(ix,iy)*rhs1_v(ix,iy)
       end do
    end do
    ! P
    call get_rhs_p(ogrd,oset,ogrd%u_sw%val,ogrd%v_sw%val,ogrd%h_sw%val,rhs1_h)
    do iy = 1,ogrd%ny_p
       do ix = 1,ogrd%nx_p
          h_tmp(ix,iy)=ogrd%h_sw%val(ix,iy)+dt*ogrd%mask_p%val(ix,iy)*rhs1_h(ix,iy)
       end do
    end do
   call set_bc_p(ogrd%nx_p,ogrd%ny_p,h_tmp(0:ogrd%nx_p+1,0:ogrd%ny_p+1),&
        & oset%wbc_flag_p,oset%ebc_flag_p,oset%nbc_flag_p,oset%sbc_flag_p)
   call set_bc_u(ogrd%nx_p,ogrd%ny_p,u_tmp(1:ogrd%nx_p+1,0:ogrd%ny_p+1),&
        & oset%wbc_flag_u,oset%ebc_flag_u,oset%nbc_flag_u,oset%sbc_flag_u,&
        & oset%slip_ind)
   call set_bc_v(ogrd%nx_p,ogrd%ny_p,v_tmp(0:ogrd%nx_p+1,1:ogrd%ny_p+1),&
        & oset%wbc_flag_v,oset%ebc_flag_v,oset%nbc_flag_v,oset%sbc_flag_v,&
        & oset%slip_ind)

   ! RK-2 step
    call get_rhs_u(ogrd,oset,u_tmp,v_tmp,h_tmp,rhs2_u)
    do iy = 1,ogrd%ny_p
       do ix = 2,ogrd%nx_p
          ! Update u
          ogrd%u_sw_next%val(ix,iy)=ogrd%u_sw%val(ix,iy) +0.5_idx*dt*ogrd%mask_u%val(ix,iy)*(rhs1_u(ix,iy)+rhs2_u(ix,iy))
       end do
    end do
    call get_rhs_v(ogrd,oset,u_tmp,v_tmp,h_tmp,rhs2_v)
     do iy = 2,ogrd%ny_p
       do ix = 1,ogrd%nx_p
          ogrd%v_sw_next%val(ix,iy)=ogrd%v_sw%val(ix,iy)+0.5_idx*dt*ogrd%mask_v%val(ix,iy)*(rhs1_v(ix,iy)+rhs2_v(ix,iy))
       end do
    end do
    ! P
    call get_rhs_p(ogrd,oset,u_tmp,v_tmp,h_tmp,rhs2_h)
    do iy = 1,ogrd%ny_p
       do ix = 1,ogrd%nx_p
          ogrd%h_sw_next%val(ix,iy)=ogrd%h_sw%val(ix,iy)+0.5_idx*dt*ogrd%mask_p%val(ix,iy)*(rhs1_h(ix,iy)+rhs2_h(ix,iy))
       end do
    end do
    call set_bc_p(ogrd%nx_p,ogrd%ny_p,&
        & ogrd%h_sw_next%val(0:ogrd%nx_p+1,0:ogrd%ny_p+1),&
        & oset%wbc_flag_p,oset%ebc_flag_p,oset%nbc_flag_p,oset%sbc_flag_p)
   call set_bc_u(ogrd%nx_p,ogrd%ny_p,&
        & ogrd%u_sw_next%val(1:ogrd%nx_p+1,0:ogrd%ny_p+1),&
        & oset%wbc_flag_u,oset%ebc_flag_u,oset%nbc_flag_u,oset%sbc_flag_u,&
        & oset%slip_ind)
   call set_bc_v(ogrd%nx_p,ogrd%ny_p,&
        & ogrd%v_sw_next%val(0:ogrd%nx_p+1,1:ogrd%ny_p+1),&
        & oset%wbc_flag_v,oset%ebc_flag_v,oset%nbc_flag_v,oset%sbc_flag_v,&
        & oset%slip_ind)
   ! Update
   ogrd%u_sw%val(1:ogrd%nx_p+1,0:ogrd%ny_p+1) = ogrd%u_sw_next%val(1:ogrd%nx_p+1,0:ogrd%ny_p+1)
   ogrd%v_sw%val(0:ogrd%nx_p+1,1:ogrd%ny_p+1) = ogrd%v_sw_next%val(0:ogrd%nx_p+1,1:ogrd%ny_p+1)
   ogrd%h_sw%val(0:ogrd%nx_p+1,0:ogrd%ny_p+1) = ogrd%h_sw_next%val(0:ogrd%nx_p+1,0:ogrd%ny_p+1)
   deallocate(rhs1_h);deallocate(rhs1_u);deallocate(rhs1_v)
   deallocate(rhs2_h);deallocate(rhs2_u);deallocate(rhs2_v)
   deallocate(h_tmp);deallocate(u_tmp);deallocate(v_tmp)
 end subroutine solve_rg_vgeo_ocn_rk2

  subroutine solve_rg_vgeo_ocn_rk3(ogrd,oset,dt)
    implicit none
    type(ocn_dta),intent(inout) :: ogrd
    type(ocn_set),intent(in) :: oset
    real(idx),intent(in) :: dt
    integer :: ix,iy
    real(idx),allocatable :: rhs1_u(:,:),rhs1_v(:,:),rhs1_h(:,:)
    real(idx),allocatable :: rhs2_u(:,:),rhs2_v(:,:),rhs2_h(:,:)
    real(idx),allocatable :: rhs3_u(:,:),rhs3_v(:,:),rhs3_h(:,:)
    real(idx),allocatable :: u_tmp1(:,:),v_tmp1(:,:),h_tmp1(:,:)
    real(idx),allocatable :: u_tmp2(:,:),v_tmp2(:,:),h_tmp2(:,:)
    allocate(rhs1_h(0:ogrd%nx_p+1,0:ogrd%ny_p+1));
    allocate(rhs2_h(0:ogrd%nx_p+1,0:ogrd%ny_p+1));
    allocate(rhs3_h(0:ogrd%nx_p+1,0:ogrd%ny_p+1));
    allocate(h_tmp1(0:ogrd%nx_p+1,0:ogrd%ny_p+1))
    allocate(h_tmp2(0:ogrd%nx_p+1,0:ogrd%ny_p+1))
    allocate(rhs1_u(1:ogrd%nx_p+1,0:ogrd%ny_p+1));
    allocate(rhs2_u(1:ogrd%nx_p+1,0:ogrd%ny_p+1));
    allocate(rhs3_u(1:ogrd%nx_p+1,0:ogrd%ny_p+1));
    allocate(u_tmp1(1:ogrd%nx_p+1,0:ogrd%ny_p+1))
    allocate(u_tmp2(1:ogrd%nx_p+1,0:ogrd%ny_p+1))
    allocate(rhs1_v(0:ogrd%nx_p+1,1:ogrd%ny_p+1));
    allocate(rhs2_v(0:ogrd%nx_p+1,1:ogrd%ny_p+1));
    allocate(rhs3_v(0:ogrd%nx_p+1,1:ogrd%ny_p+1));
    allocate(v_tmp1(0:ogrd%nx_p+1,1:ogrd%ny_p+1))
    allocate(v_tmp2(0:ogrd%nx_p+1,1:ogrd%ny_p+1))

    !RK1
    call get_rhs_u(ogrd,oset,ogrd%u_sw%val,ogrd%v_sw%val,ogrd%h_sw%val,rhs1_u)
    do iy = 1,ogrd%ny_p
       do ix = 2,ogrd%nx_p
          u_tmp1(ix,iy)=ogrd%u_sw%val(ix,iy) +dt * ogrd%mask_u%val(ix,iy) * rhs1_u(ix,iy)
       end do
    end do
    call get_rhs_v(ogrd,oset,ogrd%u_sw%val,ogrd%v_sw%val,ogrd%h_sw%val,rhs1_v)
    do iy = 2,ogrd%ny_p
       do ix = 1,ogrd%nx_p
          v_tmp1(ix,iy)=ogrd%v_sw%val(ix,iy)+ dt*ogrd%mask_v%val(ix,iy)*rhs1_v(ix,iy)
       end do
    end do
    call get_rhs_p(ogrd,oset,ogrd%u_sw%val,ogrd%v_sw%val,ogrd%h_sw%val,rhs1_h)
    do iy = 1,ogrd%ny_p
       do ix = 1,ogrd%nx_p
          h_tmp1(ix,iy)=ogrd%h_sw%val(ix,iy)+dt*ogrd%mask_p%val(ix,iy)*rhs1_h(ix,iy)
       end do
    end do
   call set_bc_p(ogrd%nx_p,ogrd%ny_p,h_tmp1(0:ogrd%nx_p+1,0:ogrd%ny_p+1),&
        & oset%wbc_flag_p,oset%ebc_flag_p,oset%nbc_flag_p,oset%sbc_flag_p)
   call set_bc_u(ogrd%nx_p,ogrd%ny_p,u_tmp1(1:ogrd%nx_p+1,0:ogrd%ny_p+1),&
        & oset%wbc_flag_u,oset%ebc_flag_u,oset%nbc_flag_u,oset%sbc_flag_u,&
        & oset%slip_ind)
   call set_bc_v(ogrd%nx_p,ogrd%ny_p,v_tmp1(0:ogrd%nx_p+1,1:ogrd%ny_p+1),&
        & oset%wbc_flag_v,oset%ebc_flag_v,oset%nbc_flag_v,oset%sbc_flag_v,&
        & oset%slip_ind)
   ! RK-2 step
    call get_rhs_u(ogrd,oset,u_tmp1,v_tmp1,h_tmp1,rhs2_u)
    do iy = 1,ogrd%ny_p
       do ix = 2,ogrd%nx_p
          ! Update u
          u_tmp2(ix,iy)=ogrd%u_sw%val(ix,iy) +dt*ogrd%mask_u%val(ix,iy)*0.25_idx*(rhs1_u(ix,iy)+rhs2_u(ix,iy))
       end do
    end do
    call get_rhs_v(ogrd,oset,u_tmp1,v_tmp1,h_tmp1,rhs2_v)
     do iy = 2,ogrd%ny_p
       do ix = 1,ogrd%nx_p
          v_tmp2(ix,iy)=ogrd%v_sw%val(ix,iy) +dt*ogrd%mask_v%val(ix,iy)*0.25_idx*(rhs1_v(ix,iy)+rhs2_v(ix,iy))
       end do
    end do
    ! P
    call get_rhs_p(ogrd,oset,u_tmp1,v_tmp1,h_tmp1,rhs2_h)
    do iy = 1,ogrd%ny_p
       do ix = 1,ogrd%nx_p
          h_tmp2(ix,iy)=ogrd%h_sw%val(ix,iy) +dt*ogrd%mask_p%val(ix,iy)*0.25_idx*(rhs1_h(ix,iy)+rhs2_h(ix,iy))
       end do
    end do
   call set_bc_p(ogrd%nx_p,ogrd%ny_p,h_tmp2(0:ogrd%nx_p+1,0:ogrd%ny_p+1),&
        & oset%wbc_flag_p,oset%ebc_flag_p,oset%nbc_flag_p,oset%sbc_flag_p)
   call set_bc_u(ogrd%nx_p,ogrd%ny_p,u_tmp2(1:ogrd%nx_p+1,0:ogrd%ny_p+1),&
        & oset%wbc_flag_u,oset%ebc_flag_u,oset%nbc_flag_u,oset%sbc_flag_u,&
        & oset%slip_ind)
   call set_bc_v(ogrd%nx_p,ogrd%ny_p,v_tmp2(0:ogrd%nx_p+1,1:ogrd%ny_p+1),&
        & oset%wbc_flag_v,oset%ebc_flag_v,oset%nbc_flag_v,oset%sbc_flag_v,&
        & oset%slip_ind)
   ! RK3
   ! Update
    call get_rhs_u(ogrd,oset,u_tmp2,v_tmp2,h_tmp2,rhs3_u)
    do iy = 1,ogrd%ny_p
       do ix = 2,ogrd%nx_p
          ! Update u
          ogrd%u_sw_next%val(ix,iy)=ogrd%u_sw%val(ix,iy)+dt*ogrd%mask_u%val(ix,iy)*&
               & (rhs1_u(ix,iy)/6.0_idx+rhs2_u(ix,iy)/6.0_idx+2.0_idx*rhs3_u(ix,iy)/3.0_idx)
       end do
    end do
    call get_rhs_v(ogrd,oset,u_tmp2,v_tmp2,h_tmp2,rhs3_v)
     do iy = 2,ogrd%ny_p
       do ix = 1,ogrd%nx_p
          ogrd%v_sw_next%val(ix,iy)=ogrd%v_sw%val(ix,iy)+dt*ogrd%mask_v%val(ix,iy)*&
               & (rhs1_v(ix,iy)/6.0_idx+rhs2_v(ix,iy)/6.0_idx+2.0_idx*rhs3_v(ix,iy)/3.0_idx)
       end do
    end do
    ! P
    call get_rhs_p(ogrd,oset,u_tmp2,v_tmp2,h_tmp2,rhs3_h)
    do iy = 1,ogrd%ny_p
       do ix = 1,ogrd%nx_p
          ogrd%h_sw_next%val(ix,iy)=ogrd%h_sw%val(ix,iy)+dt*ogrd%mask_p%val(ix,iy)*&
               & (rhs1_h(ix,iy)/6.0_idx+rhs2_h(ix,iy)/6.0_idx+2.0_idx*rhs3_h(ix,iy)/3.0_idx)
       end do
    end do
    call set_bc_p(ogrd%nx_p,ogrd%ny_p,&
        & ogrd%h_sw_next%val(0:ogrd%nx_p+1,0:ogrd%ny_p+1),&
        & oset%wbc_flag_p,oset%ebc_flag_p,oset%nbc_flag_p,oset%sbc_flag_p)
   call set_bc_u(ogrd%nx_p,ogrd%ny_p,&
        & ogrd%u_sw_next%val(1:ogrd%nx_p+1,0:ogrd%ny_p+1),&
        & oset%wbc_flag_u,oset%ebc_flag_u,oset%nbc_flag_u,oset%sbc_flag_u,&
        & oset%slip_ind)
   call set_bc_v(ogrd%nx_p,ogrd%ny_p,&
        & ogrd%v_sw_next%val(0:ogrd%nx_p+1,1:ogrd%ny_p+1),&
        & oset%wbc_flag_v,oset%ebc_flag_v,oset%nbc_flag_v,oset%sbc_flag_v,&
        & oset%slip_ind)

   ogrd%u_sw%val(1:ogrd%nx_p+1,0:ogrd%ny_p+1) = ogrd%u_sw_next%val(1:ogrd%nx_p+1,0:ogrd%ny_p+1)
   ogrd%v_sw%val(0:ogrd%nx_p+1,1:ogrd%ny_p+1) = ogrd%v_sw_next%val(0:ogrd%nx_p+1,1:ogrd%ny_p+1)
   ogrd%h_sw%val(0:ogrd%nx_p+1,0:ogrd%ny_p+1) = ogrd%h_sw_next%val(0:ogrd%nx_p+1,0:ogrd%ny_p+1)
   deallocate(rhs1_h);deallocate(rhs1_u);deallocate(rhs1_v)
   deallocate(rhs2_h);deallocate(rhs2_u);deallocate(rhs2_v)
   deallocate(rhs3_h);deallocate(rhs3_u);deallocate(rhs3_v)
   deallocate(h_tmp1);deallocate(u_tmp1);deallocate(v_tmp1)
   deallocate(h_tmp2);deallocate(u_tmp2);deallocate(v_tmp2)
 end subroutine solve_rg_vgeo_ocn_rk3
 ! Get total current
 subroutine solve_totalcurrent_ocn(ogrd,oset)
    implicit none
    type(ocn_dta),intent(inout) :: ogrd
    type(ocn_set),intent(in) :: oset
    integer :: ix,iy
    real(idx) :: H1,H2
    H1=oset%H1
    H2=oset%H2
    do iy=0,ogrd%ny_p+1
       do ix =1,ogrd%nx_p+1
          ogrd%u_ocn_1%val(ix,iy)= ogrd%u_sw%val(ix,iy)+ogrd%u_ek%val(ix,iy)*(H2)/(H1+H2)
       end do
    end do
    do iy=1,ogrd%ny_p+1
       do ix =0,ogrd%nx_p+1
          ogrd%v_ocn_1%val(ix,iy)= ogrd%v_sw%val(ix,iy)+ogrd%v_ek%val(ix,iy)*(H2)/(H1+H2)
       end do
    end do
    do iy=1,ogrd%ny_p
       do ix=1,ogrd%nx_p
          ogrd%w_ocn_1%val(ix,iy)=((ogrd%u_ocn_1%val(ix+1,iy)-ogrd%u_ocn_1%val(ix,iy))/&
               & (ogrd%x_u%val(ix+1,iy)-ogrd%x_u%val(ix,iy))  &
               & +(ogrd%v_ocn_1%val(ix,iy+1)-ogrd%v_ocn_1%val(ix,iy))/&
               & (ogrd%y_v%val(ix,iy+1)-ogrd%y_v%val(ix,iy)))*H1
       end do
    end do
  end subroutine solve_totalcurrent_ocn
 subroutine get_tsub_ZC(oset,h,hbar,sst,Tsub,Tent)
   implicit none
   type(ocn_set),intent(in) :: oset
   real(idx),intent(in) :: h,hbar,sst
   real(idx),intent(inout) :: Tsub,Tent
   real(idx) :: T1,T2,b1,b2,gamma
   T1=oset%Tsub_T1;T2=oset%Tsub_T2
   b1=oset%Tsub_b1;b2=oset%Tsub_b2
   gamma=oset%Tsub_gamma
   if (h>0) then
      Tsub=T1*(tanh(b1*(hbar+h))-tanh(b1*hbar))
   else
      Tsub=T2*(tanh(b2*(hbar-h))-tanh(b2*hbar))
   end if
   Tent=gamma*Tsub+(1.0_idx-gamma)*sst
 end subroutine get_tsub_ZC
 subroutine get_tsub_GJ2022(h,hbar,sst,Tsub,Tent)
   implicit none
   real(idx),intent(in) :: h,hbar,sst
   real(idx),intent(inout) :: Tsub,Tent
   real(idx) :: x1,x2,hsub,hstar,Asub,gamma
   hsub=75.0_idx
   hstar=60.0_idx
   Asub=11.0_idx
   x1=(hbar+h-hsub)/hstar
   x2=(hbar-hsub)/hstar
   Tsub=Asub*(tanh(x1)-tanh(x2))
   gamma=hsub/hbar
   Tent=gamma*Tsub+(1.0_idx-gamma)*sst
 end subroutine get_tsub_GJ2022
 subroutine solve_sst_ocn_ZC(ogrd,oset,mudata,mvdata,mwdata,mhdata,mTzdata,msstdata,dt)
    implicit none
    type(ocn_dta),intent(inout) :: ogrd
    type(ocn_set),intent(in) :: oset
    type(TLL_dta),intent(inout) :: mudata,mvdata,mwdata,mhdata,mTzdata,msstdata
    real(idx),intent(in) :: dt
    integer :: ix,iy,nx,ny
    real(idx) :: umTa,uaTm,uaTa,znadv
    real(idx) :: dTadx,dTmdx,ua,uc
    real(idx) :: vmTa,vaTm,vaTa,mnadv
    real(idx) :: dTady,dTmdy,va,vc
    real(idx) :: u1,v1,eps_s_sst
    real(idx) :: vadv,wtotal,mwtotal,wbar,mwbar
    real(idx) :: dTbdz,waTm,wmTa,waTa
    real(idx) :: rhs_sst,QH,H1,Te,hbar,Tsub
    real(idx) :: F1,F2
    H1=oset%H1
    eps_s_sst=1.0/(oset%eps_s_sst_day*day_to_sec)
    nx=ogrd%nx_p
    ny=ogrd%ny_p
    do iy=1,ny
       do ix=1,nx
         ! Zonal advection
          u1=0.5_idx*(mudata%data_now%val(ix,iy)*ogrd%mask_u%val(ix,iy)+&
               & mudata%data_now%val(ix+1,iy)*ogrd%mask_u%val(ix+1,iy))
          if (((u1.gt. 0) .or. (ogrd%mask_u%val(ix+1,iy).eq.0.0_idx)) .and. &
               &             (ogrd%mask_u%val(ix,iy).ne.0.0_idx)) then
             if (ogrd%mask_sst%val(ix,iy)==1.0_idx .and. &
                  & ogrd%mask_sst%val(ix-1,iy)==1.0_idx) then
                ! Upstream
                dTadx=(ogrd%ssta_ocn%val(ix,iy)*ogrd%mask_sst%val(ix,iy)&
                     & -ogrd%ssta_ocn%val(ix-1,iy)*ogrd%mask_sst%val(ix-1,iy))/&
                     & (ogrd%x_p%val(ix,iy)-ogrd%x_p%val(ix-1,iy))
                uc= mudata%data_now%val(ix,iy)*ogrd%mask_u%val(ix,iy)
             else
                uc=0.0_idx
                dTadx=0.0_idx
             end if
          else if (ogrd%mask_u%val(ix+1,iy) .ne.0.0_idx) then
             if (ogrd%mask_sst%val(ix,iy)==1.0_idx .and. &
                  & ogrd%mask_sst%val(ix+1,iy)==1.0_idx) then
                ! Downstream
                dTadx=(ogrd%ssta_ocn%val(ix+1,iy)*ogrd%mask_sst%val(ix+1,iy)&
                     & -ogrd%ssta_ocn%val(ix,iy)*ogrd%mask_sst%val(ix,iy))/&
                     & (ogrd%x_p%val(ix+1,iy)-ogrd%x_p%val(ix,iy))
                uc= mudata%data_now%val(ix+1,iy)*ogrd%mask_u%val(ix+1,iy)
             else
                uc=0.0_idx
                dTadx=0.0_idx
             end if
          else
             uc=0.0_idx
             dTadx=0.0_idx
          end if
          umTa=-1.0_idx*uc*dTadx

          u1=0.5_idx*(ogrd%u_ocn_1%val(ix,iy)*ogrd%mask_u%val(ix,iy)+&
               & ogrd%u_ocn_1%val(ix+1,iy)*ogrd%mask_u%val(ix+1,iy))
          if (((u1.gt. 0) .or. (ogrd%mask_u%val(ix+1,iy).eq.0.0_idx)) .and. &
               &             (ogrd%mask_u%val(ix,iy).ne.0.0_idx)) then
             if (ogrd%mask_sst%val(ix,iy)==1.0_idx .and. &
                  & ogrd%mask_sst%val(ix-1,iy)==1.0_idx) then
                ! Upstream
                dTadx=(ogrd%ssta_ocn%val(ix,iy)*ogrd%mask_sst%val(ix,iy)&
                     & -ogrd%ssta_ocn%val(ix-1,iy)*ogrd%mask_sst%val(ix-1,iy))/&
                     & (ogrd%x_p%val(ix,iy)-ogrd%x_p%val(ix-1,iy))
                dTmdx=(msstdata%data_now%val(ix,iy)*ogrd%mask_sst%val(ix,iy)-&
                     & msstdata%data_now%val(ix-1,iy)*ogrd%mask_sst%val(ix-1,iy))/&
                  & (ogrd%x_p%val(ix,iy)-ogrd%x_p%val(ix-1,iy))
                ua= ogrd%u_ocn_1%val(ix,iy)*ogrd%mask_u%val(ix,iy)
             else
                ua=0.0_idx
                dTadx=0.0_idx
                dTmdx=0.0_idx
             end if
          else if (ogrd%mask_u%val(ix+1,iy) .ne.0.0_idx) then
             if (ogrd%mask_sst%val(ix,iy)==1.0_idx .and. &
                  & ogrd%mask_sst%val(ix+1,iy)==1.0_idx) then
                ! Downstream
                dTadx=(ogrd%ssta_ocn%val(ix+1,iy)*ogrd%mask_sst%val(ix+1,iy)&
                     & -ogrd%ssta_ocn%val(ix,iy)*ogrd%mask_sst%val(ix,iy))/&
                     & (ogrd%x_p%val(ix+1,iy)-ogrd%x_p%val(ix,iy))
                dTmdx=(msstdata%data_now%val(ix+1,iy)*ogrd%mask_sst%val(ix+1,iy)&
                  & -msstdata%data_now%val(ix,iy)*ogrd%mask_sst%val(ix,iy))/&
                  & (ogrd%x_p%val(ix+1,iy)-ogrd%x_p%val(ix,iy))
                ua=ogrd%u_ocn_1%val(ix+1,iy)*ogrd%mask_u%val(ix+1,iy)
             else
                ua=0.0_idx
                dTadx=0.0_idx
                dTmdx=0.0_idx
             end if
          else
             ua=0.0_idx
             dTadx=0.0_idx
             dTmdx=0.0_idx
          end if
          uaTm=-1.0_idx*ua*dTmdx
          uaTa=-1.0_idx*ua*dTadx
          znadv=umTa+uaTm+uaTa
          
          ! Meridional advection
          v1=0.5_idx*(mvdata%data_now%val(ix,iy)*ogrd%mask_v%val(ix,iy)+&
               & mvdata%data_now%val(ix,iy+1)*ogrd%mask_v%val(ix,iy+1))
          if (((v1.gt. 0) .or. (ogrd%mask_v%val(ix,iy+1).eq.0.0_idx)) .and. &
               &             (ogrd%mask_v%val(ix,iy).ne.0.0_idx)) then
             if (ogrd%mask_sst%val(ix,iy)==1.0_idx .and. &
                  & ogrd%mask_sst%val(ix,iy-1)==1.0_idx) then
                ! Upstream
                dTady=(ogrd%ssta_ocn%val(ix,iy)*ogrd%mask_sst%val(ix,iy)&
                     & -ogrd%ssta_ocn%val(ix,iy-1)*ogrd%mask_sst%val(ix,iy-1))/&
                     & (ogrd%y_p%val(ix,iy)-ogrd%y_p%val(ix,iy-1))
                vc=mvdata%data_now%val(ix,iy)*ogrd%mask_v%val(ix,iy)
             else
                vc=0.0_idx;dTady=0.0_idx             
             end if
          else if (ogrd%mask_v%val(ix,iy) .ne.0.0_idx) then
             if (ogrd%mask_sst%val(ix,iy)==1.0_idx .and. &
                  & ogrd%mask_sst%val(ix,iy+1)==1.0_idx) then
             ! Downstream
                dTady=(ogrd%ssta_ocn%val(ix,iy+1)*ogrd%mask_sst%val(ix,iy+1)&
                     & -ogrd%ssta_ocn%val(ix,iy)*ogrd%mask_sst%val(ix,iy))/&
                     & (ogrd%y_p%val(ix,iy+1)-ogrd%y_p%val(ix,iy))
                vc=mvdata%data_now%val(ix,iy+1)*ogrd%mask_v%val(ix,iy+1)
             else
                vc=0.0_idx;dTady=0.0_idx             
             end if
          else
             vc=0.0_idx;dTady=0.0_idx             
          end if            
          vmTa=-1.0_idx*vc*dTady

          v1=0.5_idx*(ogrd%v_ocn_1%val(ix,iy)*ogrd%mask_v%val(ix,iy)+&
               & ogrd%v_ocn_1%val(ix,iy+1)*ogrd%mask_v%val(ix,iy+1))
          if (((v1.gt. 0) .or. (ogrd%mask_v%val(ix,iy+1).eq.0.0_idx)) .and. &
               &             (ogrd%mask_v%val(ix,iy).ne.0.0_idx)) then
             if (ogrd%mask_sst%val(ix,iy)==1.0_idx .and. &
                  & ogrd%mask_sst%val(ix,iy-1)==1.0_idx) then
                ! Upstream
                dTady=(ogrd%ssta_ocn%val(ix,iy)*ogrd%mask_sst%val(ix,iy)&
                     & -ogrd%ssta_ocn%val(ix,iy-1)*ogrd%mask_sst%val(ix,iy-1))/&
                     & (ogrd%y_p%val(ix,iy)-ogrd%y_p%val(ix,iy-1))
                dTmdy=(msstdata%data_now%val(ix,iy)*ogrd%mask_sst%val(ix,iy)&
                     & -msstdata%data_now%val(ix,iy-1)*ogrd%mask_sst%val(ix,iy-1))/&
                     & (ogrd%y_p%val(ix,iy)-ogrd%y_p%val(ix,iy-1))
                va=ogrd%v_ocn_1%val(ix,iy)*ogrd%mask_v%val(ix,iy)
             else
                va=0.0_idx;  dTmdy=0.0_idx; dTady=0.0_idx                
             end if
          else if (ogrd%mask_v%val(ix,iy) .ne.0.0_idx) then
             if (ogrd%mask_sst%val(ix,iy)==1.0_idx .and. &
                  & ogrd%mask_sst%val(ix,iy+1)==1.0_idx) then
                dTady=(ogrd%ssta_ocn%val(ix,iy+1)*ogrd%mask_sst%val(ix,iy+1)&
                     & -ogrd%ssta_ocn%val(ix,iy)*ogrd%mask_sst%val(ix,iy))/&
                     & (ogrd%y_p%val(ix,iy+1)-ogrd%y_p%val(ix,iy))
                dTmdy=(msstdata%data_now%val(ix,iy+1)*ogrd%mask_sst%val(ix,iy+1)&
                     & -msstdata%data_now%val(ix,iy)*ogrd%mask_sst%val(ix,iy))/&
                     & (ogrd%y_p%val(ix,iy+1)-ogrd%y_p%val(ix,iy))
                va=ogrd%v_ocn_1%val(ix,iy+1)*ogrd%mask_v%val(ix,iy+1)
             else
                va=0.0_idx;  dTmdy=0.0_idx; dTady=0.0_idx                
             end if
          else
             va=0.0_idx;  dTmdy=0.0_idx; dTady=0.0_idx
          end if
          vaTm=-1.0_idx*va*dTmdy
          vaTa=-1.0_idx*va*dTady
          mnadv=vmTa+vaTm+vaTa
          QH=-1.0_idx*ogrd%ssta_ocn%val(ix,iy)*eps_s_sst
          ! 
          ! Upwelling feedback
          hbar=mhdata%data_now%val(ix,iy)
          call get_tsub_ZC(oset,ogrd%h_sw%val(ix,iy),hbar,ogrd%ssta_ocn%val(ix,iy),Tsub,Te)
          ogrd%Tsub_ocn_diag%val(ix,iy)=ogrd%Tsub_ocn_diag%val(ix,iy)+ogrd%mask_p%val(ix,iy)*Tsub
          ogrd%Te_ocn_diag%val(ix,iy)=ogrd%Te_ocn_diag%val(ix,iy)+ogrd%mask_p%val(ix,iy)*Te         
          wtotal=(mwdata%data_now%val(ix,iy)+ogrd%w_ocn_1%val(ix,iy))*&
               &ogrd%mask_p%val(ix,iy)
          wbar=mwdata%data_now%val(ix,iy)*ogrd%mask_p%val(ix,iy)
          if (wbar >0) then
            F1=wbar
          else
            F1=0.0_idx
          end if  
          if (wbar <=0) then
            if (wtotal <=0) then
               F2=0.0
            else
               F2=wtotal
            end if
         else
            if (wtotal<=0) then
               F2=-1.0*wbar
            else
               F2=ogrd%w_ocn_1%val(ix,iy)*ogrd%mask_p%val(ix,iy)
            end if
         end if
         dTbdz=mTzdata%data_now%val(ix,iy)*ogrd%mask_p%val(ix,iy)
         waTa= -1.0_idx*F2*(ogrd%ssta_ocn%val(ix,iy)-Te)/H1
         waTm= -1.0_idx*F2*dTbdz
         wmTa= -1.0_idx*F1*(ogrd%ssta_ocn%val(ix,iy)-Te)/H1
         if (wtotal > 0) then
           mwtotal=wtotal
         else
           mwtotal=0.0_idx
         end if
         if (wbar > 0) then
          mwbar=wbar
           else
          mwbar=0.0_idx
         end if
         !  ! Ekman feedback
         waTm=-1.0_idx*(mwtotal-mwbar)*dTbdz
         wmTa=-1.0_idx*mwbar*(ogrd%ssta_ocn%val(ix,iy)-Te)/H1
         waTa=-1.0_idx*(mwtotal-mwbar)*(ogrd%ssta_ocn%val(ix,iy)-Te)/H1
         vadv=wmTa+waTm+waTa           
         rhs_sst=znadv+mnadv+vadv+qh                 
          ogrd%ssta_ocn_next%val(ix,iy)=ogrd%ssta_ocn%val(ix,iy)+&
               & dt*ogrd%mask_p%val(ix,iy)*rhs_sst
          ogrd%dTdt_ocn_diag%val(ix,iy)=ogrd%dTdt_ocn_diag%val(ix,iy)+dt*ogrd%mask_p%val(ix,iy)*rhs_sst
          ogrd%uaTm_ocn_diag%val(ix,iy)=ogrd%uaTm_ocn_diag%val(ix,iy)+dt*ogrd%mask_p%val(ix,iy)*uaTm
          ogrd%umTa_ocn_diag%val(ix,iy)=ogrd%umTa_ocn_diag%val(ix,iy)+dt*ogrd%mask_p%val(ix,iy)*umTa
          ogrd%uaTa_ocn_diag%val(ix,iy)=ogrd%uaTa_ocn_diag%val(ix,iy)+dt*ogrd%mask_p%val(ix,iy)*uaTa
          ogrd%vaTm_ocn_diag%val(ix,iy)=ogrd%vaTm_ocn_diag%val(ix,iy)+dt*ogrd%mask_p%val(ix,iy)*vaTm
          ogrd%vmTa_ocn_diag%val(ix,iy)=ogrd%vmTa_ocn_diag%val(ix,iy)+dt*ogrd%mask_p%val(ix,iy)*vmTa
          ogrd%vaTa_ocn_diag%val(ix,iy)=ogrd%vaTa_ocn_diag%val(ix,iy)+dt*ogrd%mask_p%val(ix,iy)*vaTa
          ogrd%waTm_ocn_diag%val(ix,iy)=ogrd%waTm_ocn_diag%val(ix,iy)+dt*ogrd%mask_p%val(ix,iy)*waTm
          ogrd%wmTa_ocn_diag%val(ix,iy)=ogrd%wmTa_ocn_diag%val(ix,iy)+dt*ogrd%mask_p%val(ix,iy)*wmTa
          ogrd%waTa_ocn_diag%val(ix,iy)=ogrd%waTa_ocn_diag%val(ix,iy)+dt*ogrd%mask_p%val(ix,iy)*waTa
          ogrd%qh_ocn_diag%val(ix,iy)=ogrd%qh_ocn_diag%val(ix,iy)+dt*ogrd%mask_p%val(ix,iy)*qh          
       end do
    end do
   call set_bc_p(ogrd%nx_p,ogrd%ny_p,&
        & ogrd%ssta_ocn_next%val(0:ogrd%nx_p+1,0:ogrd%ny_p+1),&
        & oset%wbc_flag_p,oset%ebc_flag_p,oset%nbc_flag_p,oset%sbc_flag_p)
   ogrd%ssta_ocn%val(0:ogrd%nx_p+1,0:ogrd%ny_p+1) =ogrd%ssta_ocn_next%val(0:ogrd%nx_p+1,0:ogrd%ny_p+1)
 end subroutine solve_sst_ocn_ZC
 subroutine solve_sst_ocn_ZC_ori(ogrd,oset,mudata,mvdata,mwdata,mhdata,mTzdata,msstdata,dt)
    implicit none
    type(ocn_dta),intent(inout) :: ogrd
    type(ocn_set),intent(in) :: oset
    type(TLL_dta),intent(inout) :: mudata,mvdata,mwdata,mhdata,mTzdata,msstdata
    real(idx),intent(in) :: dt
    integer :: ix,iy,nx,ny
    real(idx) :: umTa,uaTm,uaTa,znadv
    real(idx) :: dTadx,dTmdx,ua,uc
    real(idx) :: vmTa,vaTm,vaTa,mnadv
    real(idx) :: dTady,dTmdy,va,vc
    real(idx) :: u1,v1,eps_s_sst
    real(idx) :: vadv,wtotal,mwtotal,wbar,mwbar
    real(idx) :: dTbdz,waTm,wmTa,waTa
    real(idx) :: rhs_sst,QH,H1,Te,hbar,Tsub
    H1=oset%H1
    eps_s_sst=1.0/(oset%eps_s_sst_day*day_to_sec)
    nx=ogrd%nx_p
    ny=ogrd%ny_p
    do iy=1,ny
       do ix=1,nx
          u1=0.5_idx*(mudata%data_now%val(ix,iy)*ogrd%mask_u%val(ix,iy)+&
               & mudata%data_now%val(ix+1,iy)*ogrd%mask_u%val(ix+1,iy))
          if (((u1.gt. 0) .or. (ogrd%mask_u%val(ix+1,iy).eq.0.0_idx)) .and. &
               &             (ogrd%mask_u%val(ix,iy).ne.0.0_idx)) then
             if (ogrd%mask_sst%val(ix,iy)==1.0_idx .and. &
                  & ogrd%mask_sst%val(ix-1,iy)==1.0_idx) then
                ! Upstream
                dTadx=(ogrd%ssta_ocn%val(ix,iy)*ogrd%mask_sst%val(ix,iy)&
                     & -ogrd%ssta_ocn%val(ix-1,iy)*ogrd%mask_sst%val(ix-1,iy))/&
                     & (ogrd%x_p%val(ix,iy)-ogrd%x_p%val(ix-1,iy))
                uc= mudata%data_now%val(ix,iy)*ogrd%mask_u%val(ix,iy)
             else
                uc=0.0_idx
                dTadx=0.0_idx
             end if
          else if (ogrd%mask_u%val(ix+1,iy) .ne.0.0_idx) then
             if (ogrd%mask_sst%val(ix,iy)==1.0_idx .and. &
                  & ogrd%mask_sst%val(ix+1,iy)==1.0_idx) then
                ! Downstream
                dTadx=(ogrd%ssta_ocn%val(ix+1,iy)*ogrd%mask_sst%val(ix+1,iy)&
                     & -ogrd%ssta_ocn%val(ix,iy)*ogrd%mask_sst%val(ix,iy))/&
                     & (ogrd%x_p%val(ix+1,iy)-ogrd%x_p%val(ix,iy))
                uc= mudata%data_now%val(ix+1,iy)*ogrd%mask_u%val(ix+1,iy)
             else
                uc=0.0_idx
                dTadx=0.0_idx
             end if
          else
             uc=0.0_idx
             dTadx=0.0_idx
          end if
          umTa=-1.0_idx*uc*dTadx

          u1=0.5_idx*(ogrd%u_ocn_1%val(ix,iy)*ogrd%mask_u%val(ix,iy)+&
               & ogrd%u_ocn_1%val(ix+1,iy)*ogrd%mask_u%val(ix+1,iy))
          if (((u1.gt. 0) .or. (ogrd%mask_u%val(ix+1,iy).eq.0.0_idx)) .and. &
               &             (ogrd%mask_u%val(ix,iy).ne.0.0_idx)) then
             if (ogrd%mask_sst%val(ix,iy)==1.0_idx .and. &
                  & ogrd%mask_sst%val(ix-1,iy)==1.0_idx) then
                ! Upstream
                dTadx=(ogrd%ssta_ocn%val(ix,iy)*ogrd%mask_sst%val(ix,iy)&
                     & -ogrd%ssta_ocn%val(ix-1,iy)*ogrd%mask_sst%val(ix-1,iy))/&
                     & (ogrd%x_p%val(ix,iy)-ogrd%x_p%val(ix-1,iy))
                dTmdx=(msstdata%data_now%val(ix,iy)*ogrd%mask_sst%val(ix,iy)-&
                     & msstdata%data_now%val(ix-1,iy)*ogrd%mask_sst%val(ix-1,iy))/&
                  & (ogrd%x_p%val(ix,iy)-ogrd%x_p%val(ix-1,iy))
                ua= ogrd%u_ocn_1%val(ix,iy)*ogrd%mask_u%val(ix,iy)
             else
                ua=0.0_idx
                dTadx=0.0_idx
                dTmdx=0.0_idx
             end if
          else if (ogrd%mask_u%val(ix+1,iy) .ne.0.0_idx) then
             if (ogrd%mask_sst%val(ix,iy)==1.0_idx .and. &
                  & ogrd%mask_sst%val(ix+1,iy)==1.0_idx) then
                ! Downstream
                dTadx=(ogrd%ssta_ocn%val(ix+1,iy)*ogrd%mask_sst%val(ix+1,iy)&
                     & -ogrd%ssta_ocn%val(ix,iy)*ogrd%mask_sst%val(ix,iy))/&
                     & (ogrd%x_p%val(ix+1,iy)-ogrd%x_p%val(ix,iy))
                dTmdx=(msstdata%data_now%val(ix+1,iy)*ogrd%mask_sst%val(ix+1,iy)&
                  & -msstdata%data_now%val(ix,iy)*ogrd%mask_sst%val(ix,iy))/&
                  & (ogrd%x_p%val(ix+1,iy)-ogrd%x_p%val(ix,iy))
                ua=ogrd%u_ocn_1%val(ix+1,iy)*ogrd%mask_u%val(ix+1,iy)
             else
                ua=0.0_idx
                dTadx=0.0_idx
                dTmdx=0.0_idx
             end if
          else
             ua=0.0_idx
             dTadx=0.0_idx
             dTmdx=0.0_idx
          end if
          uaTm=-1.0_idx*ua*dTmdx
          uaTa=-1.0_idx*ua*dTadx
          uaTa=0.0_idx
          znadv=umTa+uaTm+uaTa
          
          ! Meridional advection
          v1=0.5_idx*(mvdata%data_now%val(ix,iy)*ogrd%mask_v%val(ix,iy)+&
               & mvdata%data_now%val(ix,iy+1)*ogrd%mask_v%val(ix,iy+1))
          if (((v1.gt. 0) .or. (ogrd%mask_v%val(ix,iy+1).eq.0.0_idx)) .and. &
               &             (ogrd%mask_v%val(ix,iy).ne.0.0_idx)) then
             if (ogrd%mask_sst%val(ix,iy)==1.0_idx .and. &
                  & ogrd%mask_sst%val(ix,iy-1)==1.0_idx) then
                ! Upstream
                dTady=(ogrd%ssta_ocn%val(ix,iy)*ogrd%mask_sst%val(ix,iy)&
                     & -ogrd%ssta_ocn%val(ix,iy-1)*ogrd%mask_sst%val(ix,iy-1))/&
                     & (ogrd%y_p%val(ix,iy)-ogrd%y_p%val(ix,iy-1))
                vc=mvdata%data_now%val(ix,iy)*ogrd%mask_v%val(ix,iy)
             else
                vc=0.0_idx;dTady=0.0_idx             
             end if
          else if (ogrd%mask_v%val(ix,iy) .ne.0.0_idx) then
             if (ogrd%mask_sst%val(ix,iy)==1.0_idx .and. &
                  & ogrd%mask_sst%val(ix,iy+1)==1.0_idx) then
             ! Downstream
                dTady=(ogrd%ssta_ocn%val(ix,iy+1)*ogrd%mask_sst%val(ix,iy+1)&
                     & -ogrd%ssta_ocn%val(ix,iy)*ogrd%mask_sst%val(ix,iy))/&
                     & (ogrd%y_p%val(ix,iy+1)-ogrd%y_p%val(ix,iy))
                vc=mvdata%data_now%val(ix,iy+1)*ogrd%mask_v%val(ix,iy+1)
             else
                vc=0.0_idx;dTady=0.0_idx             
             end if
          else
             vc=0.0_idx;dTady=0.0_idx             
          end if            
          vmTa=-1.0_idx*vc*dTady

          v1=0.5_idx*(ogrd%v_ocn_1%val(ix,iy)*ogrd%mask_v%val(ix,iy)+&
               & ogrd%v_ocn_1%val(ix,iy+1)*ogrd%mask_v%val(ix,iy+1))
          if (((v1.gt. 0) .or. (ogrd%mask_v%val(ix,iy+1).eq.0.0_idx)) .and. &
               &             (ogrd%mask_v%val(ix,iy).ne.0.0_idx)) then
             if (ogrd%mask_sst%val(ix,iy)==1.0_idx .and. &
                  & ogrd%mask_sst%val(ix,iy-1)==1.0_idx) then
                ! Upstream
                dTady=(ogrd%ssta_ocn%val(ix,iy)*ogrd%mask_sst%val(ix,iy)&
                     & -ogrd%ssta_ocn%val(ix,iy-1)*ogrd%mask_sst%val(ix,iy-1))/&
                     & (ogrd%y_p%val(ix,iy)-ogrd%y_p%val(ix,iy-1))
                dTmdy=(msstdata%data_now%val(ix,iy)*ogrd%mask_sst%val(ix,iy)&
                     & -msstdata%data_now%val(ix,iy-1)*ogrd%mask_sst%val(ix,iy-1))/&
                     & (ogrd%y_p%val(ix,iy)-ogrd%y_p%val(ix,iy-1))
                va=ogrd%v_ocn_1%val(ix,iy)*ogrd%mask_v%val(ix,iy)
             else
                va=0.0_idx;  dTmdy=0.0_idx; dTady=0.0_idx                
             end if
          else if (ogrd%mask_v%val(ix,iy) .ne.0.0_idx) then
             if (ogrd%mask_sst%val(ix,iy)==1.0_idx .and. &
                  & ogrd%mask_sst%val(ix,iy+1)==1.0_idx) then
                dTady=(ogrd%ssta_ocn%val(ix,iy+1)*ogrd%mask_sst%val(ix,iy+1)&
                     & -ogrd%ssta_ocn%val(ix,iy)*ogrd%mask_sst%val(ix,iy))/&
                     & (ogrd%y_p%val(ix,iy+1)-ogrd%y_p%val(ix,iy))
                dTmdy=(msstdata%data_now%val(ix,iy+1)*ogrd%mask_sst%val(ix,iy+1)&
                     & -msstdata%data_now%val(ix,iy)*ogrd%mask_sst%val(ix,iy))/&
                     & (ogrd%y_p%val(ix,iy+1)-ogrd%y_p%val(ix,iy))
                va=ogrd%v_ocn_1%val(ix,iy+1)*ogrd%mask_v%val(ix,iy+1)
             else
                va=0.0_idx;  dTmdy=0.0_idx; dTady=0.0_idx                
             end if
          else
             va=0.0_idx;  dTmdy=0.0_idx; dTady=0.0_idx
          end if
          vaTm=-1.0_idx*va*dTmdy
          vaTa=-1.0_idx*va*dTady
          vaTa=0.0_idx
          mnadv=vmTa+vaTm+vaTa
          QH=-1.0_idx*ogrd%ssta_ocn%val(ix,iy)*eps_s_sst
          ! 
          ! Upwelling feedback
          hbar=mhdata%data_now%val(ix,iy)
          call get_tsub_ZC(oset,ogrd%h_sw%val(ix,iy),hbar,ogrd%ssta_ocn%val(ix,iy),Tsub,Te)
          ogrd%Tsub_ocn_diag%val(ix,iy)=ogrd%Tsub_ocn_diag%val(ix,iy)+ogrd%mask_p%val(ix,iy)*Tsub
          ogrd%Te_ocn_diag%val(ix,iy)=ogrd%Te_ocn_diag%val(ix,iy)+ogrd%mask_p%val(ix,iy)*Te         
          wtotal=(mwdata%data_now%val(ix,iy)+ogrd%w_ocn_1%val(ix,iy))*&
               &ogrd%mask_p%val(ix,iy)
          wbar=mwdata%data_now%val(ix,iy)*ogrd%mask_p%val(ix,iy)
           if (wtotal > 0) then
              mwtotal=wtotal
           else
              mwtotal=0.0_idx
           end if
           if (wbar > 0) then
              mwbar=wbar
           else
              mwbar=0.0_idx
           end if
          dTbdz=mTzdata%data_now%val(ix,iy)*ogrd%mask_p%val(ix,iy)
          ! Ekman feedback
          waTa=-1.0_idx*(mwtotal-mwbar)*(ogrd%ssta_ocn%val(ix,iy)-Te)/H1
          waTm=-1.0_idx*dTbdz*(mwtotal-mwbar)
          wmTa=-1.0_idx*mwbar*(ogrd%ssta_ocn%val(ix,iy)-Te)/H1
          waTm=0.0_idx
          vadv=wmTa+waTm+waTa          
          rhs_sst=znadv+mnadv+vadv+qh                 
          ogrd%ssta_ocn_next%val(ix,iy)=ogrd%ssta_ocn%val(ix,iy)+&
               & dt*ogrd%mask_p%val(ix,iy)*rhs_sst
          ogrd%dTdt_ocn_diag%val(ix,iy)=ogrd%dTdt_ocn_diag%val(ix,iy)+dt*ogrd%mask_p%val(ix,iy)*rhs_sst
          ogrd%uaTm_ocn_diag%val(ix,iy)=ogrd%uaTm_ocn_diag%val(ix,iy)+dt*ogrd%mask_p%val(ix,iy)*uaTm
          ogrd%umTa_ocn_diag%val(ix,iy)=ogrd%umTa_ocn_diag%val(ix,iy)+dt*ogrd%mask_p%val(ix,iy)*umTa
          ogrd%uaTa_ocn_diag%val(ix,iy)=ogrd%uaTa_ocn_diag%val(ix,iy)+dt*ogrd%mask_p%val(ix,iy)*uaTa
          ogrd%vaTm_ocn_diag%val(ix,iy)=ogrd%vaTm_ocn_diag%val(ix,iy)+dt*ogrd%mask_p%val(ix,iy)*vaTm
          ogrd%vmTa_ocn_diag%val(ix,iy)=ogrd%vmTa_ocn_diag%val(ix,iy)+dt*ogrd%mask_p%val(ix,iy)*vmTa
          ogrd%vaTa_ocn_diag%val(ix,iy)=ogrd%vaTa_ocn_diag%val(ix,iy)+dt*ogrd%mask_p%val(ix,iy)*vaTa
          ogrd%waTm_ocn_diag%val(ix,iy)=ogrd%waTm_ocn_diag%val(ix,iy)+dt*ogrd%mask_p%val(ix,iy)*waTm
          ogrd%wmTa_ocn_diag%val(ix,iy)=ogrd%wmTa_ocn_diag%val(ix,iy)+dt*ogrd%mask_p%val(ix,iy)*wmTa
          ogrd%waTa_ocn_diag%val(ix,iy)=ogrd%waTa_ocn_diag%val(ix,iy)+dt*ogrd%mask_p%val(ix,iy)*waTa
          ogrd%qh_ocn_diag%val(ix,iy)=ogrd%qh_ocn_diag%val(ix,iy)+dt*ogrd%mask_p%val(ix,iy)*qh          
       end do
    end do
   call set_bc_p(ogrd%nx_p,ogrd%ny_p,&
        & ogrd%ssta_ocn_next%val(0:ogrd%nx_p+1,0:ogrd%ny_p+1),&
        & oset%wbc_flag_p,oset%ebc_flag_p,oset%nbc_flag_p,oset%sbc_flag_p)
   ogrd%ssta_ocn%val(0:ogrd%nx_p+1,0:ogrd%ny_p+1) =ogrd%ssta_ocn_next%val(0:ogrd%nx_p+1,0:ogrd%ny_p+1)
 end subroutine solve_sst_ocn_ZC_ori
 subroutine limit_total_SST(ogrd,msstdata)
    implicit none
    type(ocn_dta),intent(inout) :: ogrd
    type(TLL_dta),intent(inout) :: msstdata
    real(idx) :: total_sst
    real(idx) :: sst_max
    integer :: nx,ny,iy,ix
    sst_max  = 30.0
    nx=ogrd%nx_p
    ny=ogrd%ny_p
    do iy=1,ny
       do ix=1,nx
            if (ogrd%mask_sst%val(ix,iy)==1.0_idx) then
               total_sst=ogrd%ssta_ocn%val(ix,iy)+&
              &    msstdata%data_now%val(ix,iy)
               if (total_sst > sst_max) then
                  ogrd%ssta_ocn%val(ix,iy)=sst_max-msstdata%data_now%val(ix,iy)
               end if
            end if   
       end do
      end do
   end subroutine limit_total_SST
end module mod_ocn_solver_zc
