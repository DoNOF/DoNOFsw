###############################################################################################
#                      Makefile for DoNOF program (Date: October 2025)                        #
###############################################################################################
PROG=$(PWD)
SOU=$(PROG)/src
EXC=$(PROG)/exe
BACKUP=$(PROG)/backup
TEST=$(PROG)/test
Cln=/bin/rm -rf
LIBS = -lcint -llapack -lblas
#
SFLAGS  = -r8 -fpp -O2
CXX    =  icx
F90     = ifort          $(SFLAGS)
MPIF90  = mpiifort -DMPI $(SFLAGS)
#
SFLAGSg = -fdefault-real-8 -fdefault-double-8 -cpp -ffpe-summary=none -O2
CXXg    = gcc
F90g    = gfortran     $(SFLAGSg)
MPIF90g = mpif90 -DMPI $(SFLAGSg)
#for new versions of gnumpi:
MPIF90ng = mpif90 -DMPI -fallow-argument-mismatch $(SFLAGSg)
#
F90_FILES = $(SOU)/*.f90
###############################################################################################

all: serial mpi omp hybrid serialg mpig ompg hybridg

########################################################################################################################################################################################

serial:
	cd $(PROG)/src && $(CXX) -c -O1 nr_ecp.c
	
	cd $(PROG)/src && $(F90) -o DoNOF.x nr_ecp.o $(F90_FILES) $(LIBS)
	
	mv $(SOU)/DoNOF.x  $(EXC)/DoNOF.x

mpi:
	cd $(PROG)/src && $(CXX) -c -O1 nr_ecp.c
	
	cd $(PROG)/src && $(MPIF90) -o DoNOFmpi.x nr_ecp.o $(F90_FILES) $(LIBS)
	
	mv $(SOU)/DoNOFmpi.x  $(EXC)/DoNOFmpi.x
	
omp:
	cd $(PROG)/src && $(CXX) -c -O1 nr_ecp.c
	
	cd $(PROG)/src && $(F90) -fopenmp -o DoNOFomp.x nr_ecp.o $(F90_FILES) $(LIBS)
	
	mv $(SOU)/DoNOFomp.x  $(EXC)/DoNOFomp.x

hybrid:
	cd $(PROG)/src && $(CXX) -c -O1 nr_ecp.c
	
	cd $(PROG)/src && $(MPIF90) -fopenmp -o DoNOFhybrid.x nr_ecp.o $(F90_FILES) $(LIBS)
	
	mv $(SOU)/DoNOFhybrid.x  $(EXC)/DoNOFhybrid.x
	
serialg:
	cd $(PROG)/src && $(CXXg) -c nr_ecp.c
	
	cd $(PROG)/src && $(F90g) -o DoNOFg.x nr_ecp.o $(F90_FILES) $(LIBS)
	
	mv $(SOU)/DoNOFg.x  $(EXC)/DoNOFg.x
	
mpig:
	cd $(PROG)/src && $(CXXg) -c nr_ecp.c
	
	cd $(PROG)/src && $(MPIF90g) -o DoNOFmpig.x nr_ecp.o $(F90_FILES) $(LIBS)
	
	mv $(SOU)/DoNOFmpig.x $(EXC)/DoNOFmpig.x

mping:
	cd $(PROG)/src && $(CXXg) -c nr_ecp.c
	
	cd $(PROG)/src && $(MPIF90ng) -o DoNOFmpig.x nr_ecp.o $(F90_FILES) $(LIBS)
	
	mv $(SOU)/DoNOFmpig.x $(EXC)/DoNOFmpig.x
ompg:
	cd $(PROG)/src && $(CXXg) -c nr_ecp.c
	
	cd $(PROG)/src && $(F90g) -fopenmp -o DoNOFompg.x nr_ecp.o $(F90_FILES) $(LIBS)
	
	mv $(SOU)/DoNOFompg.x  $(EXC)/DoNOFompg.x
	
hybridg:
	cd $(PROG)/src && $(CXXg) -c nr_ecp.c
	
	cd $(PROG)/src && $(MPIF90ng) -fopenmp -o DoNOFhybridg.x nr_ecp.o $(F90_FILES) $(LIBS)
	
	mv $(SOU)/DoNOFhybridg.x $(EXC)/DoNOFhybridg.x
########################################################################################################################################################################################
