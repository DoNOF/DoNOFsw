###############################################################################################
#                      Makefile for DoNOF program (Date: January 2021)                        #
###############################################################################################
SFLAGS  = -i8 -r8 -fpp -O2 
# You can use -O3 with recent versions of the Intel Fortran compiler
F90     = ifort          $(SFLAGS)
MPIF90  = mpiifort -DMPI $(SFLAGS)
#
SFLAGSg = -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -cpp -ffpe-summary=none -O3
F90g    = gfortran     $(SFLAGSg) 
MPIF90g = mpif90 -DMPI $(SFLAGSg)
#
###############################################################################################

all: serial serialg mpi mpig

#########################################################################

serial:
	./gitversion.sh	
	$(F90) -o DoNOF.x donof.f90 mbpt.f90 gitver.f90 lapack.f

mpi:
	./gitversion.sh	
	$(MPIF90) -o DoNOFmpi.x donof.f90 mbpt.f90 gitver.f90 lapack.f
	
serialg:
	./gitversion.sh	
	$(F90g) -o DoNOFg.x donof.f90 mbpt.f90 gitver.f90 lapack.f
	
mpig:
	./gitversion.sh	
	$(MPIF90g) -o DoNOFmpig.x donof.f90 mbpt.f90 gitver.f90 lapack.f
	

#########################################################################



