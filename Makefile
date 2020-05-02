########################################################################
# Makefile for DoNOF program (Date: April 2020)
########################################################################

# Intel Fortran
F90 = ifort -i8 -r8 -fpp -static -Ofast
MPIF90 = mpiifort -DMPI -r8 -i8 -fpp -Ofast

# GNU Fortran
SFLAGS = -fdefault-integer-8 -fdefault-real-8 -cpp -O3 -ffpe-summary=none
F90g = gfortran $(SFLAGS)

########################################################################

all: serial mpi serialg 

########################################################################

serial:

	$(F90) -o donof.x donof1.f donof2.f90 

########################################################################

mpi:

	$(MPIF90) -o donofmpi.x donof1.f donof2.f90

########################################################################

serialg:

	$(F90g) -o donofgnu.x donof1.f donof2.f90 

########################################################################

