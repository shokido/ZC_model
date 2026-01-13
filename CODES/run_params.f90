module run_params
  implicit none
  integer,parameter :: idx = 8,maxlen=400
  real(idx) :: missing_value=-9999.9_idx
  real(idx),parameter :: rho0_o=1.024e3,g=9.8_idx
  real(idx),parameter :: pi = 4.0_idx *atan(1.0_idx),beta=2.28e-11_idx
  real(idx),parameter :: day_to_sec=60.0_idx * 60.0_idx * 24.0_idx  ! [s/day]
  real(idx),parameter :: sec_to_day=1.0_idx / (60.0_idx * 60.0_idx*24.0_idx)  ![day/s]
  real(idx),parameter :: year_to_sec=60.0_idx * 60.0_idx * 24.0_idx * 365 ! [s/year]
  real(idx),parameter :: dis_to_lat= 1.0_idx / (1.11e5),lat_to_dis=1.11e5
end module run_params
