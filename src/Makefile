# MOLLIE has been compiled and run with the Gnu compilers and the Intel compilers.
# Because parts of MOLLIE are quite old, it includes Fortran programming as well as C. 
# There is no problem compiling and linking both languages, but the compilers have to be
# matched to each other. In other words, C and Fortran compilers of the same type, Gnu or Intel,
# and the same version. The latest release of the Gnu compilers will always work together even if their
# version numbers are not identical. 
# Furthermore, the MPI installation and MOLLIE have to be built with the same compilers.
#
# The next two lines indicate my private MPI installation built with the Intel C and Fortran compilers
CC = /sma/data04/mpi/mpich2-build/bin/mpicc
FCC = /sma/data04/mpi/mpich2-build/bin/mpif77
FFLAGS= -O2 -w
CFLAGS= -w
FOBJECTS= fsub.o fion.o lamdbh2o.o nh3.o ch3cn.o hcn.o rate.o subcritical.o
COBJECTS= cmain.o define_model.o
all: simulate

fsub.o: fsub.f nlines_f77.h
	$(FCC)  -c ${FFLAGS} fsub.f
fion.o: fion.f nlines_f77.h
	$(FCC)  -c ${FFLAGS} fion.f
lamdbh2o.o: lamdbh2o.f nlines_f77.h
	$(FCC)  -c ${FFLAGS} lamdbh2o.f
nh3.o: nh3.f nlines_f77.h
	$(FCC)  -c ${FFLAGS} nh3.f
ch3cn.o: ch3cn.f nlines_f77.h
	$(FCC)  -c ${FFLAGS} ch3cn.f
hcn.o: hcn.f nlines_f77.h
	$(FCC)  -c ${FFLAGS} hcn.f
rate.o: rate.f nlines_f77.h
	$(FCC)  -c ${FFLAGS} rate.f
subcritical.o: subcritical.f nlines_f77.h
	$(FCC)  -c ${FFLAGS} subcritical.f



cmain.o: cmain.c setup.c global.h cmain.h nlines_C.h
	$(CC)  -c ${CFLAGS} cmain.c
define_model.o: define_model.c global.h cmain.h nlines_C.h
	$(CC)  -c ${CFLAGS} define_model.c

all: simulate

# This is the link command for the newer Gnu C and Fortran compilers. In this case we link with the Fortran
# compiler. We don't need to indicate the F77 I/O libraries but we indicate the C math library "-lm". 
# This is the one I use with Ubuntu linux on my laptop.
# simulate: ${FOBJECTS} ${COBJECTS}
#	$(FCC) -o c.x ${COBJECTS} ${FOBJECTS} -lm

# This is the link command for the Intel C and Fortran compilers that I still use. These are old versions now.
# With the older compilers we link with the C compiler and explicitly include the libraries for 
# the Fortran I/O statements. The location and names of these libraries vary with the operating system 
# and installation.
simulate: ${FOBJECTS} ${COBJECTS}
	$(CC) -o mollie ${COBJECTS} ${FOBJECTS} -lm -L/sma/data04/mpi/intel/fce/10.0.023/lib -lifcore_pic

# This was the link command for an old DEC alpha
# simulate: ${FOBJECTS} ${COBJECTS}
#	$(CC) -o c.x ${FFLAGS} ${COBJECTS} ${FOBJECTS} -lCEPCF90 -lF90 -lm
# This was the link command for an old linux with MPI built with older gcc and g77
# simulate: ${FOBJECTS} ${COBJECTS}
#	$(CC) -o c.x ${COBJECTS} ${FOBJECTS} -lm -lg2c 


clean:
	-rm cmain.o define_model.o fsub.o fion.o lamdbh2o.o nh3.o ch3cn.o hcn.o subcritical.o mollie
