## Makefile of Linearly continuosly stratified mode
#==========================
#Set compiler
FC  = gfortran
#==========================
#Set compile option
NETCDF_INCDIR=-I/opt/homebrew/Cellar/netcdf-fortran/4.6.1/include -I/opt/homebrew/Cellar/netcdf-fortran/4.6.1/include
NETCDF_LIBDIR=-L/opt/homebrew/Cellar/netcdf-fortran/4.6.1/lib -L/opt/homebrew/Cellar/netcdf/4.9.2_2/lib
#==========================
FFLAGS += $(NETCDF_INCDIR)
LDFLAGS += $(NETCDF_LIBDIR) -lnetcdf -lnetcdff
LDFLAGS += -cpp
LDFLAGS += -O3
FFLAGS += -DINI
FFLAGS += -DWSTRESS_BULK
#FFLAGS += -CB
FFLAGS_LEGACY+=-std=legacy -Wno-argument-mismatch
#==========================
# You don't have to edit the following part
#==========================
TARGET_ogcm_dyn=exec_solver_ocn_dyn.out
TARGET_ogcm_full=exec_solver_ocn_full.out
TARGET_agcm=exec_solver_atm.out
TARGET_cgcm_full=exec_solver_couple_full.out

OBJ_ncdf=ncdf_read.o \
	 ncdf_write.o
OBJ_calendar=calendar_sub.o
OBJ_fft=fftpack5.1.o
OBJ_run= run_params.o \
	run_types.o 
OBJ_io= mod_io_master.o \
	mod_io_avg.o \
	mod_io_diag.o
OBJ_atm_driver=mod_atm_grid.o \
	mod_atm_solver_gill.o
OBJ_ocn_driver=mod_ocn_grid.o \
	mod_wstress.o \
	mod_ocn_solver_zc.o 
OBJ_couple_driver=mod_coupler_dta.o 

OBJ_ocn_full = zc_ogcm_full.o
OBJ_ocn_dyn = zc_ogcm_dyn.o
OBJ_atm = zc_agcm.o
OBJ_couple_full = zc_cgcm_full.o
MODS  = ncdf_read.mod \
	ncdf_write.mod \
	calendar_sub.mod \
	run_params.mod \
	run_types.mod \
	mod_io_master.mod  \
	mod_io_avg.mod  \
	mod_io_diag.mod  \
	mod_ocn_dta.mod  \
	mod_wstress.mod  \
	mod_ocn_solver_zc.mod  \
	mod_coupler_dta.mod  \
	mod_atm_dta.mod  \
	mod_atm_solver_gill.mod 

all : $(TARGET_ogcm_dyn) $(TARGET_ogcm_full) $(TARGET_agcm) $(TARGET_cgcm_full)  
ogcm: $(TARGET_ogcm_dyn) $(TARGET_ogcm_full)
agcm: $(TARGET_agcm)
cgcm: $(TARGET_cgcm)
$(TARGET_ogcm_full): $(OBJ_ncdf) $(OBJ_calendar) $(OBJ_run) $(OBJ_io) $(OBJ_ocn_driver) $(OBJ_ocn_full)
	$(FC) -o $@ $(OBJ_ncdf) $(OBJ_calendar) $(OBJ_run) $(OBJ_io) $(OBJ_ocn_driver) $(OBJ_ocn_full) ${LDFLAGS} ${FFLAGS}
$(TARGET_ogcm_dyn): $(OBJ_ncdf) $(OBJ_calendar) $(OBJ_run) $(OBJ_io) $(OBJ_ocn_driver) $(OBJ_ocn_dyn)
	$(FC) -o $@ $(OBJ_ncdf) $(OBJ_calendar) $(OBJ_run) $(OBJ_io) $(OBJ_ocn_driver) $(OBJ_ocn_dyn) ${LDFLAGS} ${FFLAGS}

$(TARGET_agcm): $(OBJ_ncdf) $(OBJ_calendar) $(OBJ_fft) $(OBJ_run) $(OBJ_io) $(OBJ_atm_driver) $(OBJ_atm)
	$(FC) -o $@ $(OBJ_ncdf) $(OBJ_calendar) $(OBJ_fft) $(OBJ_run) $(OBJ_io) $(OBJ_atm) $(OBJ_atm_driver) ${LDFLAGS} ${FFLAGS}

$(TARGET_cgcm_full): $(OBJ_ncdf) $(OBJ_calendar) $(OBJ_fft) $(OBJ_run) $(OBJ_io) $(OBJ_atm_driver) $(OBJ_ocn_driver) $(OBJ_couple_driver) $(OBJ_couple_full)
	$(FC) -o $@ $(OBJ_ncdf) $(OBJ_calendar) $(OBJ_fft) $(OBJ_run) $(OBJ_io) $(OBJ_atm_driver) $(OBJ_ocn_driver) $(OBJ_couple_driver) $(OBJ_couple_full) ${LDFLAGS} ${FFLAGS}

#Objects--------------------------------------
fftpack5.1.o:fftpack5.1.f90
	$(FC) -c $< ${LDFLAGS} ${FFLAGS_LEGACY}
%.o: %.f90
	$(FC) -c $< ${LDFLAGS} ${FFLAGS}
clean:
	rm -f $(TARGET_ogcm_full) $(TARGET_ogcm_dyn) $(TARGET_agcm) $(TARGET_cgcm_full) $(OBJ_ncdf) $(OBJ_calendar) $(OBJ_fft) $(OBJ_run) $(OBJ_io) $(OBJ_ocn_dyn) $(OBJ_ocn_full) $(OBJ_ocn_driver) $(OBJ_atm) $(OBJ_atm_driver) $(OBJ_stcouple) $(OBJ_couple_full) $(OBJ_couple_driver)  $(MODS)
oclean:
	rm -f $(TARGET_ogcm)  $(OBJ_run) $(OBJ_io) $(OBJ_ocn_full) $(OBJ_ocn_dyn) $(OBJ_ocn_driver) $(OBJ_couple_full) 
