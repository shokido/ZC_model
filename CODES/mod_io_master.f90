module mod_io_master
  use run_params
  use run_types
  use calendar_sub
  use ncdf_read
  implicit none  
contains
  function modify_time(ntime,time_in,time_units,start_yymmdd,start_hhmmss) result(time)
    implicit none
    integer,intent(in) :: ntime
    real(idx),intent(in) :: time_in(ntime)
    character(len=*),intent(in) :: time_units
    integer,intent(in) :: start_yymmdd,start_hhmmss
    real(idx) :: time(ntime)
    character :: flag_char*8,yr_char*4,mn_char*2,dy_char*2,hr_char*2,min_char*2,sec_char
    integer :: ref_year,ref_month,ref_day,ref_hour,ref_min,ref_sec
    integer :: ind1,ind2
    integer :: ref_yymmdd,ref_hhmmss
    integer :: it
    integer :: flag
    integer :: tmp_yymmdd,tmp_hhmmss
    ind1=1
    ind2=index(time_units,"since")-2
    flag_char=time_units(ind1:ind2)
    ind1=index(time_units,"since")+6
    ind2=ind1+index(time_units(ind1:),"-")-1
    yr_char=time_units(ind1:ind2-1)
    ind1=ind2+1
    ind2=ind1+index(time_units(ind1:),"-")-1
    mn_char=time_units(ind1:ind2-1)
    ind1=ind2+1
    ind2=ind1+index(time_units(ind1:)," ")-1
    dy_char=time_units(ind1:ind2-1)
    ind1=ind2+1
    ind2=ind1+index(time_units(ind1:),":")-1
    hr_char=time_units(ind1:ind2-1)
    ind1=ind2+1
    ind2=ind1+index(time_units(ind1:),":")-1
    min_char=time_units(ind1:ind2-1)
    ind1=ind2+1
    ind2=len_trim(time_units)
    sec_char=time_units(ind1:ind2)

    read(yr_char,*) ref_year ; read(mn_char,*) ref_month ; read(dy_char,*) ref_day
    ref_yymmdd=ref_year*10000+ref_month*100+ref_day
    read(hr_char,*) ref_hour ; read(min_char,*) ref_min ; read(sec_char,*) ref_sec
    ref_hhmmss=ref_hour*10000+ref_min*100+ref_sec

    select case(trim(flag_char))
    case("seconds")
       flag=-10000
    case("minitues")
       flag=-100
    case("hours")
       flag=-1
    case("days")
       flag=1
    case("months")
       flag=100
    case("years")
       flag=10000
    end select
    do it = 1,ntime
       call calendar_cal_ymdhms_after(ref_yymmdd,ref_hhmmss,time_in(it),flag,tmp_yymmdd,tmp_hhmmss)
       call calendar_cal_length_ymdhms(start_yymmdd,start_hhmmss,tmp_yymmdd,tmp_hhmmss,1,time(it))
    end do
  end function modify_time
  subroutine time_wgt(time_array,time,i1,i2,w1,w2)
    real(idx),intent(in) :: time_array(:),time
    integer,intent(out) :: i1,i2    
    integer :: ntime
    real(idx) :: w1,w2
    ntime=sum(shape(time_array))
    if (ntime .ne. 1) then
       if (time  .gt. minval(time_array) .and. time .lt. maxval(time_array)) then
          i1 = sum(maxloc(time_array,mask=(time_array<=time)))
          i2 = i1 + 1
          w1 = (time_array(i2)-time) / (time_array(i2)-time_array(i1))
          w2 = (time-time_array(i1)) / (time_array(i2)-time_array(i1))
       else if (time .le. minval(time_array)) then
          i1 = 1
          i2 = 1
          w1 = 1.0_idx
          w2 = 0.0_idx
       else
          i1 = ntime
          i2 = ntime
          w1 = 1.0_idx
          w2 = 0.0_idx
       end if
    else
       i1 = 1; i2 =1
       w1 = 1.0_idx ; w2 = 0.0_idx
    end if
  end subroutine time_wgt
  function set_data(ind1,ind2,wgt1,wgt2,data_1d) result(data_ret)
    implicit none
    integer,intent(in) :: ind1,ind2
    real(idx),intent(in) :: wgt1,wgt2
    real(idx),intent(in) :: data_1d(:)
    real(idx) :: data_ret
    data_ret= wgt1*data_1d(ind1)+wgt2*data_1d(ind2)
  end function set_data
  subroutine read_data_TLL_atm(nfile,fnames,timename,varname,agrd,wgrd,start_yymmdd,start_hhmmss)
    implicit none
    integer,intent(in) :: nfile
    character(len=maxlen),intent(in) :: fnames(nfile),timename,varname
    integer,intent(in) :: start_yymmdd,start_hhmmss
    type(atm_dta),intent(in) :: agrd
    type(TLL_dta),intent(inout) :: wgrd
    integer :: nt,nx,ny,ntime
    real(idx),allocatable :: time_tmp(:),time(:),time_out(:),data(:,:,:)
    character(len=maxlen) :: time_units
    integer :: ifile,it1,it2
    ntime=0
    do ifile=1,nfile
       call get_dimsize(fnames(ifile),timename,nt)
       ntime=ntime+nt
    end do
    nx=agrd%nx_atm
    ny=agrd%ny_atm

    it1=1;it2=0
    allocate(wgrd%time%val(1:ntime))
    allocate(wgrd%data%val(1:nx,1:ny,1:ntime))
    allocate(wgrd%data_now%val(1:nx,1:ny))
    allocate(wgrd%data_mod%val(1:nx,1:ny))

    do ifile=1,nfile
       call get_dimsize(trim(fnames(ifile)),timename,nt)
       call get_variable(trim(fnames(ifile)),timename,(/1/),(/nt/),time_tmp)
       call get_attribute(trim(fnames(ifile)),timename,"units",time_units)
       allocate(time(nt))
       time=modify_time(nt,time_tmp,time_units,start_yymmdd,start_hhmmss)
       it2=it1+nt-1
       call get_variable(trim(fnames(ifile)),varname,(/1,1,1/),(/nx,ny,nt/),data)
       wgrd%time%val(it1:it2)=time(1:nt)
       wgrd%data%val(1:nx,1:ny,it1:it2)=data(1:nx,1:ny,1:nt)
       it1=it2+1
       deallocate(time)
    end do
  end subroutine read_data_TLL_atm
  subroutine read_data_TLL_p(nfile,fnames,timename,varname,ogrd,wgrd,start_yymmdd,start_hhmmss)
    implicit none
    integer,intent(in) :: nfile
    character(len=maxlen),intent(in) :: fnames(nfile),timename,varname
    integer,intent(in) :: start_yymmdd,start_hhmmss
    type(ocn_dta),intent(in) :: ogrd
    type(TLL_dta),intent(inout) :: wgrd
    integer :: nt,nx,ny,ntime
    real(idx),allocatable :: time_tmp(:),time(:),time_out(:),data(:,:,:)
    character(len=maxlen) :: time_units
    integer :: ifile,it1,it2
    ntime=0
    do ifile=1,nfile
       call get_dimsize(fnames(ifile),timename,nt)
       ntime=ntime+nt
    end do
    nx=ogrd%nx_p
    ny=ogrd%ny_p

    it1=1;it2=0
    allocate(wgrd%time%val(1:ntime))
    allocate(wgrd%data%val(0:nx+1,0:ny+1,1:ntime))
    allocate(wgrd%data_now%val(0:nx+1,0:ny+1))
    allocate(wgrd%data_mod%val(0:nx+1,0:ny+1))

    do ifile=1,nfile
       call get_dimsize(trim(fnames(ifile)),timename,nt)
       call get_variable(trim(fnames(ifile)),timename,(/1/),(/nt/),time_tmp)
       call get_attribute(trim(fnames(ifile)),timename,"units",time_units)
       allocate(time(nt))
       time=modify_time(nt,time_tmp,time_units,start_yymmdd,start_hhmmss)
       it2=it1+nt-1
       call get_variable(trim(fnames(ifile)),varname,(/1,1,1/),(/nx+2,ny+2,nt/),data)
       wgrd%time%val(it1:it2)=time(1:nt)
       wgrd%data%val(0:nx+1,0:ny+1,it1:it2)=data(1:nx+2,1:ny+2,1:nt)
       it1=it2+1
       deallocate(time)
    end do
  end subroutine read_data_TLL_p
  subroutine read_data_TLL_u(nfile,fnames,timename,varname,ogrd,wgrd,start_yymmdd,start_hhmmss)
    implicit none
    integer,intent(in) :: nfile
    character(len=maxlen),intent(in) :: fnames(nfile),timename,varname
    integer,intent(in) :: start_yymmdd,start_hhmmss
    type(ocn_dta),intent(in) :: ogrd
    type(TLL_dta),intent(inout) :: wgrd
    integer :: nt,nx,ny,ntime
    real(idx),allocatable :: time_tmp(:),time(:),time_out(:),data(:,:,:)
    character(len=maxlen) :: time_units
    integer :: ifile,it1,it2
    ntime=0
    do ifile=1,nfile
       call get_dimsize(fnames(ifile),timename,nt)
       ntime=ntime+nt
    end do

    nx=ogrd%nx_p
    ny=ogrd%ny_p
    it1=1;it2=0
    allocate(wgrd%time%val(1:ntime))
    allocate(wgrd%data%val(1:nx+1,0:ny+1,1:ntime))
    allocate(wgrd%data_now%val(1:nx+1,0:ny+1))
    allocate(wgrd%data_mod%val(1:nx+1,0:ny+1))
    do ifile=1,nfile
       call get_dimsize(trim(fnames(ifile)),timename,nt)
       call get_variable(trim(fnames(ifile)),timename,(/1/),(/nt/),time_tmp)
       call get_attribute(trim(fnames(ifile)),timename,"units",time_units)
       allocate(time(nt))
       time=modify_time(nt,time_tmp,time_units,start_yymmdd,start_hhmmss)
       it2=it1+nt-1
       call get_variable(trim(fnames(ifile)),varname,(/1,1,1/),(/nx+1,ny+2,nt/),data)
       wgrd%time%val(it1:it2)=time(1:nt)
       wgrd%data%val(1:nx+1,0:ny+1,it1:it2)=data(1:nx+1,1:ny+2,1:nt)
       it1=it2+1
       deallocate(time)
    end do
  end subroutine read_data_TLL_u
  subroutine read_data_TLL_v(nfile,fnames,timename,varname,ogrd,wgrd,start_yymmdd,start_hhmmss)
    implicit none
    integer,intent(in) :: nfile
    character(len=maxlen),intent(in) :: fnames(nfile),timename,varname
    integer,intent(in) :: start_yymmdd,start_hhmmss
    type(ocn_dta),intent(in) :: ogrd
    type(TLL_dta),intent(inout) :: wgrd
    integer :: nt,nx,ny,ntime
    real(idx),allocatable :: time_tmp(:),time(:),time_out(:),data(:,:,:)
    character(len=maxlen) :: time_units
    integer :: ifile,it1,it2
    ntime=0
    do ifile=1,nfile
       call get_dimsize(fnames(ifile),timename,nt)
       ntime=ntime+nt
    end do
    nx=ogrd%nx_p
    ny=ogrd%ny_p

    it1=1;it2=0
    allocate(wgrd%time%val(1:ntime))
    allocate(wgrd%data%val(0:nx+1,1:ny+1,1:ntime))
    allocate(wgrd%data_now%val(0:nx+1,1:ny+1))
    allocate(wgrd%data_mod%val(0:nx+1,1:ny+1))

    do ifile=1,nfile
       call get_dimsize(trim(fnames(ifile)),timename,nt)
       call get_variable(trim(fnames(ifile)),timename,(/1/),(/nt/),time_tmp)
       call get_attribute(trim(fnames(ifile)),timename,"units",time_units)
       allocate(time(nt))
       time=modify_time(nt,time_tmp,time_units,start_yymmdd,start_hhmmss)
       it2=it1+nt-1
       call get_variable(trim(fnames(ifile)),varname,(/1,1,1/),(/nx+2,ny+1,nt/),data)
       wgrd%time%val(it1:it2)=time(1:nt)
       wgrd%data%val(0:nx+1,1:ny+1,it1:it2)=data(1:nx+2,1:ny+1,1:nt)
       it1=it2+1
       deallocate(time)
    end do
  end subroutine read_data_TLL_v
  subroutine read_coupler(fname,ogrd,agrd,sgrd)
    implicit none
    character(len=maxlen),intent(in) :: fname
    type(ocn_dta),intent(in) :: ogrd
    type(atm_dta),intent(in) :: agrd
    type(stat_dta),intent(inout) :: sgrd
    integer :: nx_a,ny_a,nx_o,ny_o,ntime
    integer :: it
    real(idx),allocatable :: data_2d(:,:),data_3d(:,:,:)
    call get_dimsize(fname,"time",ntime)
    nx_o=ogrd%nx_p;ny_o=ogrd%ny_p
    nx_a=agrd%nx_atm; ny_a=agrd%ny_atm
    allocate(sgrd%mask_O%val(0:nx_o+1,0:ny_o+1))
    allocate(sgrd%reg_taux%val(1:nx_a,1:ny_a,1:ntime))
    allocate(sgrd%reg_tauy%val(1:nx_a,1:ny_a,1:ntime))
    sgrd%ind1=1
    call get_variable(fname,"mask_O",(/1,1,1/),(/nx_o+2,ny_o+2/),data_2d)
    sgrd%mask_O%val(0:nx_o+1,0:ny_o+1)=data_2d(1:nx_o+2,1:ny_o+2)
    do it =1,ntime
       call get_variable(fname,"reg_taux",(/1,1,it/),(/nx_a,ny_a,it/),data_3d)
       sgrd%reg_taux%val(1:nx_a,1:ny_a,it)=data_3d(1:nx_a,1:ny_a,1)
       call get_variable(fname,"reg_tauy",(/1,1,it/),(/nx_a,ny_a,it/),data_3d)
       sgrd%reg_tauy%val(1:nx_a,1:ny_a,it)=data_3d(1:nx_a,1:ny_a,1)
    end do
    deallocate(data_2d)
    deallocate(data_3d)
  end subroutine read_coupler
  subroutine get_data_TLL_atm(time_now,agrd,wgrd)
    implicit none
    real(idx),intent(in) :: time_now
    type(atm_dta),intent(in) :: agrd
    type(TLL_dta),intent(inout) :: wgrd
    real(idx) :: time_int
    real(idx) :: tmp1
    integer :: ix,iy
    if (wgrd%Lcycle .eq. "T") then
       time_int=time_now-int(time_now/(wgrd%Tcycle))*wgrd%Tcycle
    else
       time_int=time_now
    end if
    call time_wgt(wgrd%time%val,time_int,wgrd%ind1,wgrd%ind2,wgrd%wgt1,wgrd%wgt2)
    do iy = 1,agrd%ny_atm
       do ix = 1,agrd%nx_atm
          tmp1=set_data(wgrd%ind1,wgrd%ind2,wgrd%wgt1,wgrd%wgt2,wgrd%data%val(ix,iy,:))
          wgrd%data_now%val(ix,iy)=tmp1
       end do
    end do
  end subroutine get_data_TLL_atm
  subroutine get_data_TLL_p(time_now,ogrd,wgrd)
    implicit none
    real(idx),intent(in) :: time_now
    type(ocn_dta),intent(in) :: ogrd
    type(TLL_dta),intent(inout) :: wgrd
    real(idx) :: time_int
    real(idx) :: tmp1
    integer :: ix,iy
    if (wgrd%Lcycle .eq. "T") then
       time_int=time_now-int(time_now/(wgrd%Tcycle))*wgrd%Tcycle
    else
       time_int=time_now
    end if
    call time_wgt(wgrd%time%val,time_int,wgrd%ind1,wgrd%ind2,wgrd%wgt1,wgrd%wgt2)
    do iy = 0,ogrd%ny_p+1
       do ix = 0,ogrd%nx_p+1
          tmp1=set_data(wgrd%ind1,wgrd%ind2,wgrd%wgt1,wgrd%wgt2,wgrd%data%val(ix,iy,:))
          wgrd%data_now%val(ix,iy)=tmp1
       end do
    end do
  end subroutine get_data_TLL_p
  subroutine get_data_TLL_u(time_now,ogrd,wgrd)
    implicit none
    real(idx),intent(in) :: time_now
    type(ocn_dta),intent(in) :: ogrd
    type(TLL_dta),intent(inout) :: wgrd
    real(idx) :: time_int
    real(idx) :: tmp1
    integer :: ix,iy
    if (wgrd%Lcycle .eq. "T") then
       time_int=time_now-int(time_now/(wgrd%Tcycle))*wgrd%Tcycle
    else
       time_int=time_now
    end if
    call time_wgt(wgrd%time%val,time_int,wgrd%ind1,wgrd%ind2,wgrd%wgt1,wgrd%wgt2)
    do iy = 0,ogrd%ny_p+1
       do ix = 1,ogrd%nx_p+1
          tmp1=set_data(wgrd%ind1,wgrd%ind2,wgrd%wgt1,wgrd%wgt2,wgrd%data%val(ix,iy,:))
          wgrd%data_now%val(ix,iy)=tmp1
       end do
    end do
  end subroutine get_data_TLL_u
  subroutine get_data_TLL_v(time_now,ogrd,wgrd)
    implicit none
    real(idx),intent(in) :: time_now
    type(ocn_dta),intent(in) :: ogrd
    type(TLL_dta),intent(inout) :: wgrd
    real(idx) :: time_int
    real(idx) :: tmp1
    integer :: ix,iy
    if (wgrd%Lcycle .eq. "T") then
       time_int=time_now-int(time_now/(wgrd%Tcycle))*wgrd%Tcycle
    else
       time_int=time_now
    end if
    call time_wgt(wgrd%time%val,time_int,wgrd%ind1,wgrd%ind2,wgrd%wgt1,wgrd%wgt2)
    do iy = 1,ogrd%ny_p+1
       do ix = 0,ogrd%nx_p+1
          tmp1=set_data(wgrd%ind1,wgrd%ind2,wgrd%wgt1,wgrd%wgt2,wgrd%data%val(ix,iy,:))
          wgrd%data_now%val(ix,iy)=tmp1
       end do
    end do
  end subroutine get_data_TLL_v  
end module mod_io_master

