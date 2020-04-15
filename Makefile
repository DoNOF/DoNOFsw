########################################################################
# Makefile for DoNOF program (Date: April 2020)
########################################################################
PROG=$(HOME)/DoNOF/
SOU=$(PROG)/sources
OBJ1=$(PROG)/objects/serial
OBJ2=$(PROG)/objects/mpi
EXC=$(PROG)/exe
BACKUP=$(PROG)/backup
SCR=$(PROG)/scripts
EXAMP=$(PROG)/examples
DOC=$(PROG)/doc/
Cln=/bin/rm -rf
#
SFLAGS   = -i8 -r8 -fpp -O0
SFLAGSO0 = -i8 -r8 -fpp -static -O0
SFLAGSw  = -i8 -r8 -fpp -static -O0 -warn all,nodeclarations
SFLAGSb  = -i8 -r8 -fpp -static -O0 -CB
SFLAGSO3 = -i8 -r8 -fpp -static -O3
SFLAGSO4 = -i8 -r8 -fpp -static -Ofast
F90      = ifort $(SFLAGSO4)
#
PFLAGSO0 = -DMPI -i8 -r8 -fpp -O0
PFLAGSw  = -DMPI -i8 -r8 -fpp -Ofast -warn all,nodeclarations
PFLAGSb  = -DMPI -i8 -r8 -fpp -Ofast -CB
PFLAGSO3 = -DMPI -i8 -r8 -fpp -O3
PFLAGSO4 = -DMPI -i8 -r8 -fpp -Ofast
MPIF90   = mpiifort $(PFLAGSO4)
#
########################################################################

all: serial mpi

########################################################################

serial:
	cd $(OBJ1)&& $(Cln) *.o *genmod*
	
	cd $(SOU) && $(F90) -c donof1.f donof2.f90 
		
	cd $(SOU) && $(F90) -o donof.x donof1.o donof2.o 
	
# move exe and object files
	
	mv $(SOU)/donof.x  $(EXC)/donof.x
	
	make mvfiles1
# clean
	cd $(SOU) && $(Cln) *genmod*

########################################################################

mpi:
	cd $(OBJ2)&& $(Cln) *.o *genmod*
	
	cd $(SOU) && $(MPIF90) -c donof1.f donof2.f90
	
	cd $(SOU) && $(MPIF90) -o donofmpi.x donof1.o donof2.o
	
# move exe and object files
	
	mv $(SOU)/donofmpi.x  $(EXC)/donofmpi.x
	
	make mvfiles2
	
# clean
	cd $(SOU) && $(Cln) *genmod*
#
mvfiles1:
	mv $(SOU)/donof1.o      $(OBJ1)/donof1.o
	mv $(SOU)/donof2.o      $(OBJ1)/donof2.o
#
mvfiles2:
	mv $(SOU)/donof1.o      $(OBJ2)/donof1.o
	mv $(SOU)/donof2.o      $(OBJ2)/donof2.o
	
###############################################################################
#
tar:    
	cd $(BACKUP)/ && tar -zPcvf DoNOF_2020.04.tar.gz               \
			            $(SOU) $(DOC) $(EXAMP) $(SCR)
#
###############################################################################
