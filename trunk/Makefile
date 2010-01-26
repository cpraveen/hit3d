# Makefile with logic in it.  For different hosts runs different things

# Define the correct flags and compilers

MPIF90 = blah

# Franklin cluster, NERSC
ifeq ($(NERSC_HOST), franklin)
        MPIF90 = ftn

        FCFLAGS = -target=linux -i8 -r8 -O4 -c
        LDFLAGS = -target=linux -i8 -r8 -O4 -lfftw3
        FCFLAGS_F77 = -Mextend
endif

# Hopper cluster, NERSC
ifeq ($(NERSC_HOST), hopper)
        MPIF90 = ftn

        FCFLAGS = -target=linux -i8 -r8 -O4 -c
        LDFLAGS = -target=linux -i8 -r8 -O4 -lfftw3
        FCFLAGS_F77 = -Mextend
endif

# Yellowrail cluster, LANL
ifeq ($(HOSTNAME), yr-fe1.lanl.gov)
	MPIF90 = mpif90

	FCFLAGS = -i8 -r8 -O4 -c $(MPI_COMPILE_FLAGS) -I$(FFTW_INCLUDE)
	LDFLAGS = -i8 -r8 -O4 -lfftw3 $(MPI_LD_FLAGS) -L$(FFTW_HOME)/lib
endif

# Coyote cluster, LANL
ifeq ($(HOSTNAME), cy-c2.lanl.gov)
	MPIF90 = mpif90
	FCFLAGS = -i8 -r8 -O4 -c $(MPI_COMPILE_FLAGS) -I$(FFTW_INCLUDE)
	LDFLAGS = -i8 -r8 -O4 -lfftw3 $(MPI_LD_FLAGS) -L$(FFTW_HOME)/lib
endif

# WCR cluster, Stanford
ifeq ($(HOSTNAME), glacial.stanford.edu)
        FFTW_HOME = /home/chumakov/packages/fftw-3.2/pgi
        MPIF90 = mpif90
        FCFLAGS = -i8 -r8 -c $(MPI_COMPILE_FLAGS) -I$(FFTW_HOME)/include
        LDFLAGS = -i8 -r8 $(MPI_LD_FLAGS) -L$(FFTW_HOME)/lib -lfftw3 -lm
        FCFLAGS_F77 = -Mextend

endif

# sparrow.stanford.edu, Mac OS X box
ifeq ($(HOSTNAME), sparrow.stanford.edu)
        MPIF90 = mpif90
        FCFLAGS = -fdefault-real-8 -fdefault-integer-8 -finit-integer=0 -finit-real=zero -c 
        FCFLAGS_F77 = -ffixed-form -ffixed-line-length-none
        FCFLAGS_F90 = -ffree-form -ffree-line-length-none
        LDFLAGS = -fdefault-real-8 -fdefault-integer-8 -finit-integer=0 -finit-real=zero -lmpi -lmpi -lfftw3 -lm
endif

# Program name
PROG    = hit3d.x

# Modules
MODULES = m_openmpi.o\
	m_io.o\
	m_parameters.o\
	m_work.o\
	m_fields.o\
	m_timing.o\
	x_fftw.o\
	m_filter_xfftw.o\
	m_particles.o\
	m_stats.o\
	m_force.o\
	m_rand_knuth.o\
	m_les.o\
	RANDu.o

# Objects
OBJ     = main.o\
	begin_new.o\
	begin_restart.o\
	dealias_all.o\
	get_file_ext.o\
	init_velocity.o\
	init_scalars.o\
	io_write_4.o\
	my_dt.o\
	my_exit.o\
	pressure.o\
	restart_io.o\
	rhs_velocity.o\
	rhs_scalars.o\
	velocity_rescale.o\
	write_tmp4.o


# -------------------------------------------------------
# link

$(PROG):  $(MODULES) $(OBJ)
	$(MPIF90) $(MODULES) $(OBJ) -o $(PROG) $(LDFLAGS)
# -------------------------------------------------------
# compile

$(OBJ): $(MODULES) 

%.o: %.f
	$(MPIF90) $(FCFLAGS) $(FCFLAGS_F77)  $<

%.o: %.f90
	$(MPIF90) $(FCFLAGS) $<

clean:
	rm *.o *.mod $(PROG)

