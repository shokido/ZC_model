module run_types
  use run_params
  type :: vector_1d
     real(idx),allocatable :: val(:)
  end type vector_1d
  type :: vector_2d
     real(idx),allocatable :: val(:,:)
  end type vector_2d
  type :: int_2d
     integer,allocatable :: val(:,:)
  end type int_2d
  type :: int_3d
     integer,allocatable :: val(:,:,:)
  end type int_3d
  type :: vector_3d
     real(idx),allocatable :: val(:,:,:)
  end type vector_3d
  type :: ocn_dta
     integer :: nx_p,ny_p
     integer :: nx_u,ny_u
     integer :: nx_v,ny_v
     type(vector_1d) :: lon_p, lat_p
     type(vector_1d) :: lon_u, lat_u
     type(vector_1d) :: lon_v, lat_v
     type(vector_2d) :: x_p, y_p
     type(vector_2d) :: x_u, y_u
     type(vector_2d) :: x_v, y_v
     type(vector_2d) :: f,mask_p,mask_u,mask_v
     type(vector_2d) :: mask_phi_u,mask_phi_v
     type(vector_2d) :: mask_sst
     type(vector_2d) :: mask_wwb
     ! Velocity array
     type(vector_2d) :: u_sw,v_sw,h_sw
     type(vector_2d) :: u_sw_past,v_sw_past,h_sw_past
     type(vector_2d) :: u_sw_next,v_sw_next,h_sw_next
     type(vector_2d) :: u_ek,v_ek,w_ek
     type(vector_2d) :: u_ocn_1,v_ocn_1,w_ocn_1
     type(vector_2d) :: ssta_ocn,ssta_ocn_next

     ! Velocity averaged array
     type(vector_2d) :: u_sw_avg,v_sw_avg,h_sw_avg
     type(vector_2d) :: u_ek_avg,v_ek_avg,w_ek_avg
     type(vector_2d) :: u_ocn_1_avg,v_ocn_1_avg,w_ocn_1_avg
     type(vector_2d) :: damp_u,damp_v,damp_p
     type(vector_2d) :: visc_2D
     ! SST averaged array
     type(vector_2d) :: ssta_ocn_avg

     ! Velocity averaged array
     type(vector_2d) :: ua_ocn,va_ocn
     type(vector_2d) :: taux_ocn,tauy_ocn
     type(vector_2d) :: taux_ocn_avg,tauy_ocn_avg
     ! SST tendencies
     type(vector_2d) :: dTdt_ocn_diag
     type(vector_2d) :: uaTm_ocn_diag,umTa_ocn_diag,uaTa_ocn_diag
     type(vector_2d) :: vaTm_ocn_diag,vmTa_ocn_diag,vaTa_ocn_diag
     type(vector_2d) :: waTm_ocn_diag,wmTa_ocn_diag,waTa_ocn_diag
     type(vector_2d) :: qh_ocn_diag
     type(vector_2d) :: Te_ocn_diag
     type(vector_2d) :: Tsub_ocn_diag
  end type ocn_dta
  type :: ocn_set
     character(len=maxlen) :: wbc_flag_p="Clo" ! "Clo" or "Gra"
     character(len=maxlen) :: ebc_flag_p="Clo"
     character(len=maxlen) :: nbc_flag_p="Clo"
     character(len=maxlen) :: sbc_flag_p="Clo"
     character(len=maxlen) :: wbc_flag_u="Clo"
     character(len=maxlen) :: ebc_flag_u="Clo"
     character(len=maxlen) :: nbc_flag_u="Clo"
     character(len=maxlen) :: sbc_flag_u="Clo"
     character(len=maxlen) :: wbc_flag_v="Clo"
     character(len=maxlen) :: ebc_flag_v="Clo"
     character(len=maxlen) :: nbc_flag_v="Clo"
     character(len=maxlen) :: sbc_flag_v="Clo"
     real(idx) :: slip_ind=0.0_idx ! slip_ind=0 (du/dx=0), slip_ind=1 (u=0)
     real(idx) :: cp_ocn=2.9_idx ! Speed of baroclinic gravity wave (in [m/s])
     real(idx) :: r_ocn_day=912.5_idx ! Inverse of damping coefficient for baroclinic current (in [days])
     real(idx) :: H1=50.0_idx ! Mixed layer depth (in [m])
     real(idx) :: H2=80.0_idx ! Upper layer depth (in [m])
     real(idx) :: eps_s_ocn_day=2.0_idx ! Inverse of damping coefficient for Ekman current (in [days])
     real(idx) :: eps_s_sst_day=125.0_idx ! Inverse of damping coefficient for SST (in [days])
     real(idx) :: cd_bulk=0.0018_idx ! Bulk coefficient for wind stress calculation
     real(idx) :: nu=100000.0 ! Horizontal viscocity in [m^2/s]
     real(idx) :: Tsub_T1=28.0_idx 
     real(idx) :: Tsub_T2=-40.0_idx
     real(idx) :: Tsub_b1=1.0_idx/(80.0_idx)
     real(idx) :: Tsub_b2=1.0_idx/(33.0_idx)
     real(idx) :: Tsub_gamma=0.75_idx
      ! WWB parameter 
     character(1) :: kick_WWB="F"
     real(idx) :: G0_wwb=0.09_idx,G1_wwb=0.0_idx,dur_wwb=20.0_idx
     real(idx) :: us0_wwb=6.5_idx,widx_wwb=20.0_idx,widy_wwb=6.0_idx
     real(idx) :: x0_wwb=160.0_idx,y0_wwb=0.0_idx
     end type ocn_set
  type :: atm_dta
     integer :: nx_atm,ny_atm
     type(vector_1d) :: lon_atm, lat_atm
     type(vector_1d) :: x_atm, y_atm,k_atm
     type(vector_2d) :: qa_atm, pa_atm, ua_atm,va_atm
     type(vector_2d) :: qa_atm_avg, pa_atm_avg, ua_atm_avg,va_atm_avg
     type(vector_2d) :: sstm_atm, ssta_atm
     type(vector_2d) :: sstm_atm_avg, ssta_atm_avg
  end type atm_dta
  type :: atm_set
     real(idx) :: cp_atm=60.0_idx ! Speed of atmospheric gravity wave (in [m/s])
     real(idx) :: eps_s_atm_day=1.106364457323311_idx
     real(idx) :: alpha_gill_atm=0.031_idx     ! In [m^2*s^(-3)*K^(-1)]
     character(len=maxlen) :: heating_type="ZC87" ! Heating function form; "ZC87"; "ZC87conv"; "GJ22"
  end type atm_set
  type :: TLL_dta
     integer :: ind1,ind2,ntime
     real(idx) :: wgt1,wgt2
     type(vector_1d) :: time
     type(vector_1d) :: time_cyc 
     type(vector_2d) :: data_mod
     type(vector_2d) :: data_now
     type(vector_3d) :: data
     character(1) :: Lcycle='F'
     real(idx) :: Tcycle=365.0_idx
  end type TLL_dta
  type :: couple_dta
     type(int_3d) :: ind_AtoO_x,ind_AtoO_y
     type(int_3d) :: ind_OtoA_x,ind_OtoA_y
     type(vector_3d) :: wgt_AtoO,wgt_OtoA
  end type couple_dta
  type :: stat_dta
     integer :: ntime,ind1
     type(vector_2d) :: mask_O
     type(vector_3d) :: reg_taux,reg_tauy
  end type stat_dta
end module run_types
