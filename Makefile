###############################################################################################
#                      Makefile for DoNOF program (Date: May 2026)                            #
###############################################################################################
PROG=$(PWD)
SOU=$(PROG)/src
EXC=$(PROG)/exe
TEST=$(PROG)/test
Cln=/bin/rm -rf
LIBS = -lcint -llapack -lblas
#
SFLAGSg = -fdefault-real-8 -fdefault-double-8 -cpp -ffpe-summary=none -O2
CXXg    = gcc
F90g    = gfortran     $(SFLAGSg)
MPIF90g = mpif90 -DMPI -fallow-argument-mismatch $(SFLAGSg)
#
F90_FILES = $(SOU)/*.f90
###############################################################################################

all: serialg ompg mpig

###############################################################################################

serialg:
	cd $(PROG)/src && $(CXXg) -c nr_ecp.c
	
	cd $(PROG)/src && $(F90g) -o DoNOFg.x nr_ecp.o $(F90_FILES) $(LIBS)
	
	mv $(SOU)/DoNOFg.x  $(EXC)/DoNOFg.x
	
mpig:
	cd $(PROG)/src && $(CXXg) -c nr_ecp.c
	
	cd $(PROG)/src && $(MPIF90g) -o DoNOFmpig.x nr_ecp.o $(F90_FILES) $(LIBS)
	
	mv $(SOU)/DoNOFmpig.x $(EXC)/DoNOFmpig.x

ompg:
	cd $(PROG)/src && $(CXXg) -c nr_ecp.c
	
	cd $(PROG)/src && $(F90g) -fopenmp -o DoNOFompg.x nr_ecp.o $(F90_FILES) $(LIBS)
	
	mv $(SOU)/DoNOFompg.x  $(EXC)/DoNOFompg.x
	
###############################################################################################
