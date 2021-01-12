###############################################################################################
#                      Makefile for DoNOF program (Date: January 2021)                        #
###############################################################################################
SFLAGS  = -i8 -r8 -fpp -O2
F90     = ifort          $(SFLAGS)
MPIF90  = mpiifort -DMPI $(SFLAGS)
#
SFLAGSg = -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -cpp -ffpe-summary=none -O3
F90g    = gfortran     $(SFLAGSg) 
MPIF90g = mpif90 -DMPI $(SFLAGSg)
#
###############################################################################################

all: serial mpi serialg mpig

#########################################################################

serial:
	$(F90) -o DoNOF.x donof.f90 lapack.f

mpi:
	$(MPIF90) -o DoNOFmpi.x donof.f90 lapack.f
	
serialg:
	$(F90g) -o DoNOFg.x donof.f90 lapack.f
	
mpig:
	$(MPIF90g) -o DoNOFmpig.x donof.f90 lapack.f
	

#########################################################################



