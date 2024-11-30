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
SFLAGSg = -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -cpp -ffpe-summary=none -O2
F90g    = gfortran     $(SFLAGSg)
MPIF90g = mpif90 -DMPI $(SFLAGSg)
#for new versions of gnumpi:
MPIF90ng = mpif90 -DMPI -fallow-argument-mismatch $(SFLAGSg)
#
F90_FILES = $(SOU)/*.f90
F_FILES = $(SOU)/*.f
###############################################################################################

all: serial mpi omp hybrid serialg mpig ompg hybridg

########################################################################################################################################################################################

serial:
	cd $(PROG)/src && $(F90) -o DoNOF.x $(F90_FILES) $(F_FILES)
	
	mv $(SOU)/DoNOF.x  $(EXC)/DoNOF.x

mpi:
	cd $(PROG)/src && $(MPIF90) -o DoNOFmpi.x $(F90_FILES) $(F_FILES)
	
	mv $(SOU)/DoNOFmpi.x  $(EXC)/DoNOFmpi.x
	
omp:
	cd $(PROG)/src && $(F90) -fopenmp -o DoNOFomp.x $(F90_FILES) $(F_FILES)
	
	mv $(SOU)/DoNOFomp.x  $(EXC)/DoNOFomp.x

hybrid:
	cd $(PROG)/src && $(MPIF90) -fopenmp -o DoNOFhybrid.x $(F90_FILES) $(F_FILES)
	
	mv $(SOU)/DoNOFhybrid.x  $(EXC)/DoNOFhybrid.x
	
serialg:
	cd $(PROG)/src && $(F90g) -o DoNOFg.x $(F90_FILES) $(F_FILES)
	
	mv $(SOU)/DoNOFg.x  $(EXC)/DoNOFg.x
	
mpig:
	cd $(PROG)/src && $(MPIF90g) -o DoNOFmpig.x $(F90_FILES) $(F_FILES) 
	
	mv $(SOU)/DoNOFmpig.x $(EXC)/DoNOFmpig.x

mping:
	cd $(PROG)/src && $(MPIF90ng) -o DoNOFmpig.x $(F90_FILES) $(F_FILES) 
	
	mv $(SOU)/DoNOFmpig.x $(EXC)/DoNOFmpig.x
ompg:
	cd $(PROG)/src && $(F90g) -fopenmp -o DoNOFompg.x $(F90_FILES) $(F_FILES)
	
	mv $(SOU)/DoNOFompg.x  $(EXC)/DoNOFompg.x
	
hybridg:
	cd $(PROG)/src && $(MPIF90ng) -fopenmp -o DoNOFhybridg.x $(F90_FILES) $(F_FILES) 
	
	mv $(SOU)/DoNOFhybridg.x $(EXC)/DoNOFhybridg.x
########################################################################################################################################################################################
