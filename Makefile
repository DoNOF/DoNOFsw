###############################################################################################
#                      Makefile for DoNOF program (Date: April 2024)                          #
###############################################################################################
PROG=$(PWD)
SOU=$(PROG)/src
EXC=$(PROG)/exe
BACKUP=$(PROG)/backup
TEST=$(PROG)/test
Cln=/bin/rm -rf
#
SFLAGS  = -i8 -r8 -fpp -O2
F90     = ifort          $(SFLAGS)
MPIF90  = mpiifort -DMPI $(SFLAGS)
#
SFLAGSg = -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -cpp -ffpe-summary=none -O1
F90g    = gfortran     $(SFLAGSg)
MPIF90g = mpif90 -DMPI $(SFLAGSg)
#for new versions of gnumpi:
MPIF90ng = mpif90 -DMPI -fallow-argument-mismatch $(SFLAGSg)
#
F90_FILES = $(SOU)/*.f90
F_FILES = $(SOU)/*.f
###############################################################################################

all: serial mpi serialg mpig

########################################################################################################################################################################################

serial:
	cd $(PROG)/src && ./gitversion.sh && $(F90) -o DoNOF.x $(F90_FILES) $(F_FILES)
	
	mv $(SOU)/DoNOF.x  $(EXC)/DoNOF.x

mpi:
	cd $(PROG)/src && ./gitversion.sh && $(MPIF90) -o DoNOFmpi.x $(F90_FILES) $(F_FILES)
	
	mv $(SOU)/DoNOFmpi.x  $(EXC)/DoNOFmpi.x
	
serialg:
	cd $(PROG)/src && ./gitversion.sh && $(F90g) -o DoNOFg.x $(F90_FILES) $(F_FILES)
	
	mv $(SOU)/DoNOFg.x  $(EXC)/DoNOFg.x
	
mpig:
	cd $(PROG)/src && ./gitversion.sh && $(MPIF90g) -o DoNOFmpig.x $(F90_FILES) $(F_FILES) 
	
	mv $(SOU)/DoNOFmpig.x $(EXC)/DoNOFmpig.x

mping:
	cd $(PROG)/src && ./gitversion.sh && $(MPIF90ng) -o DoNOFmpig.x $(F90_FILES) $(F_FILES) 
	
	mv $(SOU)/DoNOFmpig.x $(EXC)/DoNOFmpig.x
########################################################################################################################################################################################
