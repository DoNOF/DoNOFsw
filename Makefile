########################################################################
# Makefile for DoNOF program (Date: April 2020)
########################################################################

# Intel Fortran
F90 = ifort -i8 -r8 -fpp -static -O3 -mkl
MPIF90 = mpiifort -DMPI -r8 -i8 -fpp -O3 -mkl

# GNU Fortran
SFLAGS = -fdefault-integer-8 -fdefault-real-8 -cpp -O3 -ffpe-summary=none -llapack -lblas
F90g = gfortran $(SFLAGS)

########################################################################

all: serial serialg mpi

########################################################################

serial:

	$(F90) -o donof.x donof1.f mbpt.f donof2.f90 gauss_legendre.f90 

########################################################################

mpi:

	$(MPIF90) -o donofmpi.x donof1.f mbpt.f donof2.f90 gauss_legendre.f90

########################################################################

serialg:

	$(F90g) -o donofgnu.x donof1.f mbpt.f donof2.f90 gauss_legendre.f90

########################################################################

tar:   
	rm DoNOF.tar 
	tar -cvf DoNOF.tar *f *f90 *h Makefile 

########################################################################
