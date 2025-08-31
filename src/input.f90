!======================================================================!
!                                                                      !
!                   I N P U T   S U B R O U T I N E S                  !
!                                                                      !
!======================================================================!
!======================================================================!

! START                                            
      SUBROUTINE START(IGTYP,NATOMS,NATmax,NBF,NBFaux,NUQMT,NELEC,NALP, &
                      NBET,NSHELL,NSHELLmax,NSHELLSaux,NPRIMI,NPRIMImax,&
                      NPRIMIaux,ZAN,Cxyz,IAN,IMIN,IMAX,ZMASS,KSTART,    &
                      KATOM,KTYPE,KNG,KLOC,KMIN,KMAX,INTIPO,ISHPIR,     &
                      ITIPO,C1,C2,EX,CS,CP,CD,CF,CG,CH,CI,SIZE_ENV,ENV, &
                      ATM,BAS) !LIBCINT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER(4) :: ERITYP,GEN 
      CHARACTER(5) :: RITYP 
      CHARACTER(8) :: UNITS
      LOGICAL EFLDL,LINEAR                                                  
!
      COMMON/INFOA/NAT,ICH,MUL,NUM,NQMT,NE,NA,NB
      COMMON/INFOB/NUMaux,NSHELLaux
      COMMON/EFLDC_1/EFLDL
      COMMON/EFLDC_2/EVEC(3) 
      COMMON/ELPROP/IEMOM      
      COMMON/CONTROL/UNITS  
      COMMON/ZMAT/LINEAR 
      COMMON/INTOPT/ISCHWZ,IECP,NECP            
      COMMON/RUNTYPE/IRUNTYP
      COMMON/WRTGCF/IWRTGCF                                                  
      LOGICAL SMCD
      COMMON/ERITYPE/IERITYP,IRITYP,IGEN,ISTAR,MIXSTATE,SMCD
!      
      INTEGER :: IGTYP
      INTEGER,DIMENSION(NATmax) :: IAN,IMIN,IMAX
      INTEGER,DIMENSION(NSHELLmax) :: KSTART,KATOM,KTYPE,KNG
      INTEGER,DIMENSION(NSHELLmax) :: KLOC,KMIN,KMAX,INTIPO
      INTEGER,DIMENSION(NPRIMImax) :: ISHPIR,ITIPO
      DOUBLE PRECISION,DIMENSION(NPRIMImax) :: C1,C2,EX,CS,CP
      DOUBLE PRECISION,DIMENSION(NPRIMImax) :: CD,CF,CG,CH,CI
      DOUBLE PRECISION,DIMENSION(NATmax) :: ZAN,ZMASS
      DOUBLE PRECISION,DIMENSION(3,NATmax) :: Cxyz
      DOUBLE PRECISION,DIMENSION(3) :: VMOI
!
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: COM

      INTEGER :: SIZE_ENV                         !LIBCINT
      DOUBLE PRECISION :: ENV(SIZE_ENV)           !LIBCINT
      INTEGER :: ATM(6,NATmax), BAS(8,NSHELLmax)  !LIBCINT
!-----------------------------------------------------------------------
!     Read the Molecule and its normal basis set
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL MOLECULE(IGTYP,NSHELL,NPRIMI,NATmax,NSHELLmax,NPRIMImax,     &
                    NPRIMIaux,C1,C2,EX,CS,CP,CD,CF,CG,CH,CI,KSTART,     &
                    KATOM,KTYPE,KNG,KLOC,KMIN,KMAX,IMIN,IMAX,INTIPO,    &
                    ISHPIR,ITIPO,ZAN,Cxyz,SIZE_ENV,ENV,ATM,BAS) !LIBCINT
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     True Nuclear Charges
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO I=1,NAT                                                      
       IAN(I) = INT(ZAN(I)+0.001D0)                             
      ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Atomic Mass Table  ( ZMASS: normal masses )
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL AMT(ZAN,ZMASS,NAT)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Check for Linear Molecule
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(COM(3,NAT))
      CALL VCLR(VMOI,1,3)                                               
      IF(NAT>0)CALL INRTIA(Cxyz,COM,ZMASS,VMOI,NAT)                    
      LINEAR = .FALSE.                                                    
      IF(VMOI(1)<1.0D-04)LINEAR = .TRUE. 
      IF(NAT==1)LINEAR = .FALSE.               
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Print Input Run Options
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      WRITE(6,1)IRUNTYP,MUL,ICH,IECP,IEMOM,UNITS
      if(EFLDL)then
       WRITE(6,'(/1X,A15,3F10.5)')'Electric Field:',(EVEC(I),I=1,3)
      end if
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Write GCF file for Restart calculations
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IWRTGCF = 1
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     ECP Input
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL ECPPAR(KTYPE,NSHELL,ZAN,Cxyz)                                                           
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Schwarz inequality on if NAT>5: Integrals < 10.0**(-9) aren't used
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(NAT>5)ISCHWZ = 1                                                    
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     ERITYP & GEN
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    
      IF(IERITYP==1)THEN
       ERITYP = 'FULL'
       WRITE(6,2)1.0D-09,ISCHWZ,ERITYP      
      ELSE IF(IERITYP==2 .OR. IERITYP==3)THEN
       IF(IERITYP==2) ERITYP = 'RI'
       IF(IERITYP==3) ERITYP = 'MIX'
       IF(IRITYP==1) RITYP = 'JKFIT'
       IF(IRITYP==2) RITYP = 'GEN'
       IF(IRITYP==1) THEN
        WRITE(6,3)1.0D-09,ISCHWZ,ERITYP,RITYP       
       ELSE IF(IRITYP==2) THEN
        IF(ISTAR==0)WRITE(GEN,'(A1,I1)')    'A',IGEN
        IF(ISTAR==1)WRITE(GEN,'(A1,I1,A1)') 'A',IGEN,'*'
        WRITE(6,4)1.0D-09,ISCHWZ,ERITYP,RITYP,GEN       
       END IF
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Abort Program if LMAXIMA > 5
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL BASCHK(LMAXIMA,KTYPE,NSHELL)
      IF(LMAXIMA>5)THEN                                                
       WRITE(6,'(35A)')'Functions with LMAX > 5 not allowed'           
       CALL ABRT                                                      
      END IF                                                            
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Passing values to DoNOF through START
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      NBF    = NUM
      NATOMS = NAT
      NUQMT  = NQMT
      NELEC  = NE
      NALP   = NA
      NBET   = NB
!     RI
      NBFaux = NUMaux
      NSHELLSaux = NSHELLaux
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DEALLOCATE(COM)
      RETURN                                                            
!-----------------------------------------------------------------------
    1 FORMAT(/1X,'Input Run Options',/,                                 &
              1X,'-----------------',/,                                 &
       1X,'IRUNTYP =',I2,2X,'MULT =',I2,2X,'ICHARG =',I3,1X,'IECP =',I2,&
       2X,'IEMOM =',I2,2X,'UNITS = ',A8)                    
    2 FORMAT(/1X,'Integral Options'/,1X,16(1H-)/1X,'CUTOFF =',1P,E8.1,  &
              3X,'ISCHWZ =',I2,3X,'ERITYP =',A5)
    3 FORMAT(/1X,'Integral Options'/,1X,16(1H-)/1X,'CUTOFF =',1P,E8.1,  &
              3X,'ISCHWZ =',I2,3X,'ERITYP =',A5,3X,'AUX =',A6)
    4 FORMAT(/1X,'Integral Options'/,1X,16(1H-)/1X,'CUTOFF =',1P,E8.1,  &
              3X,'ISCHWZ =',I2,3X,'ERITYP =',A5,3X,'AUX =',A6,           &
              3X,'GEN =',A5)
!-----------------------------------------------------------------------
      END                                                               

! MOLECULE
      SUBROUTINE MOLECULE(IGTYP,NSHELL,NPRIMI,NATmax,NSHELLmax,         &
                         NPRIMImax,NPRIMIaux,C1,C2,EX,CS,CP,CD,CF,CG,CH,&
                         CI,KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX,IMIN, &
                         IMAX,INTIPO,ISHPIR,ITIPO,ZAN,Cxyz,SIZE_ENV,ENV,&
                         ATM,BAS) !LIBCINT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/FRAME /U1,U2,U3,V1,V2,V3,WW1,WW2,WW3,X0,Y0,Z0     ! PTGRP            
      CHARACTER(80) :: TITLE, BASIS_FILE
      COMMON/TIT/TITLE
      COMMON/BASIS_FILE/BASIS_FILE
      INTEGER,DIMENSION(NSHELLmax) :: KSTART,KATOM,KTYPE,KNG
      INTEGER,DIMENSION(NSHELLmax) :: KLOC,KMIN,KMAX
      DOUBLE PRECISION,DIMENSION(NPRIMImax) :: C1,C2,EX,CS,CP
      DOUBLE PRECISION,DIMENSION(NPRIMImax) :: CD,CF,CG,CH,CI
      DOUBLE PRECISION,DIMENSION(NATmax) :: ZAN
      DOUBLE PRECISION,DIMENSION(3,NATmax) :: Cxyz
!
      INTEGER :: IGTYP
      INTEGER,DIMENSION(NATmax)  :: IMIN,IMAX
      INTEGER,DIMENSION(NSHELLmax) :: INTIPO
      INTEGER,DIMENSION(NPRIMImax) :: ISHPIR,ITIPO

      INTEGER :: SIZE_ENV                           !LIBCINT
      DOUBLE PRECISION :: ENV(SIZE_ENV)             !LIBCINT
      INTEGER :: ATM(6,NATmax), BAS(8,NSHELLmax)    !LIBCINT
!
#include "mpip.h"
!-----------------------------------------------------------------------
!     Read Atom Coordinates and Basis from the Input File ($DATA)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#ifdef MPI
      CONTINUE
#else
      REWIND(5)
#endif
      CALL FNDGRP(5,' $DATA  ',IEOF)
      IF(IEOF/=0)THEN
       WRITE(6,1)' $DATA  '
       CALL ABRT
      ENDIF
      CALL OPNCRD
      READ (5,'(A80)')TITLE
      WRITE(6,2)TITLE
      READ (5,'(A80)')BASIS_FILE
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Read Atoms
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL ATOMS(IGTYP,NSHELL,NPRIMI,NATmax,NSHELLmax,NPRIMImax,        &
                 NPRIMIaux,C1,C2,EX,CS,CP,CD,CF,CG,CH,CI,KSTART,KATOM,  &
                 KTYPE,KNG,KLOC,KMIN,KMAX,IMIN,IMAX,INTIPO,ISHPIR,ITIPO,&
                 ZAN,Cxyz,SIZE_ENV,ENV,ATM,BAS) !LIBCINT
      CALL SETLAB(KATOM,KMIN,KMAX,NSHELL,ZAN)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
!-----------------------------------------------------------------------
    1 FORMAT(1X,'**** ERROR, NO ',A8,' GROUP WAS FOUND')
    2 FORMAT(/1X,'Run Title',/1X,9(1H-)/1X,A80)
!-----------------------------------------------------------------------
      END

! ATOMS                                            
      SUBROUTINE ATOMS(IGTYP,NSHELL,NPRIMI,NATmax,NSHELLmax,NPRIMImax,  &
                       NPRIMIaux,C1PIR,C2PIR,EX,CS,CP,CD,CF,CG,CH,CI,   &
                       KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX,IMINPIR,   &
                       IMAXPIR,INTIPO,ISHPIR,ITIPO,ZAN,Cxyz,SIZE_ENV,   &
                       ENV,ATM,BAS) !LIBCINT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION,PARAMETER :: PT2953=29.53125D0
      DOUBLE PRECISION,PARAMETER :: PT1624=162.421875D0
      DOUBLE PRECISION,PARAMETER :: PT75=0.75D0
      DOUBLE PRECISION,PARAMETER :: PT187=1.875D0
      DOUBLE PRECISION,PARAMETER :: TM10=1.0D-10
      DOUBLE PRECISION,PARAMETER :: PT6562=6.5625D0
      DOUBLE PRECISION,PARAMETER :: UNIT=0.52917724924D0
      CHARACTER(8), DIMENSION(103,7) :: ABASIS
      INTEGER, DIMENSION(103,7) :: IAGAUS
      COMMON/INTNAL/NATIN 
      COMMON/INTOPT/ISCHWZ,IECP,NECP            
      COMMON/INFOA/NAT,ICH,MUL,NUM,NQMT,NE,NA,NB
      COMMON/CONV/ACURCY,EN,Etot,EHF,EHF0,DIFF,ITER,ICALP,ICBET                    
      COMMON/INFOB/NUMaux,NSHELLaux
      COMMON/MAPSHEL/MAPSHL(600,48),NT
      COMMON/TRANSF/XP,YP,ZP 
!
      COMMON/USELIBCINT/ILIBCINT
      LOGICAL SMCD
      COMMON/ERITYPE/IERITYP,IRITYP,IGEN,ISTAR,MIXSTATE,SMCD
!
      INTEGER,DIMENSION(NSHELLmax) :: KSTART,KATOM,KTYPE,KNG
      INTEGER,DIMENSION(NSHELLmax) :: KLOC,KMIN,KMAX
      DOUBLE PRECISION,DIMENSION(NPRIMImax) :: C1PIR,C2PIR
      DOUBLE PRECISION,DIMENSION(NPRIMImax) :: EX,CS,CP,CD,CF,CG,CH,CI
      DOUBLE PRECISION,DIMENSION(NATmax) :: ZAN
      DOUBLE PRECISION,DIMENSION(3,NATmax) :: Cxyz
!      
      INTEGER,DIMENSION(NATmax)  :: IMINPIR,IMAXPIR
      INTEGER,DIMENSION(NSHELLmax) :: INTIPO
      INTEGER,DIMENSION(NPRIMImax) :: ISHPIR,ITIPO       
      INTEGER,DIMENSION(NATmax,48) :: MAPCTR
      CHARACTER(8),DIMENSION(NATmax) :: ANAM
      CHARACTER(2),DIMENSION(NATmax) :: BNAM
!
      CHARACTER(8) :: BLANK
      DATA BLANK /'        '/
      CHARACTER(8) :: CBASIS                       
      CHARACTER(8),DIMENSION(8) :: LABEL
      DATA LABEL/'S       ','P       ','D       ','F       ',           &
                 'G       ','H       ','I       ','L       '/
      INTEGER,DIMENSION(8) :: NBFS,MINF,MAXF,NANGM
      DATA NBFS / 1, 3,  6, 10, 15, 21, 28, 4/                 
      DATA MINF / 1, 2,  5, 11, 21, 36, 57, 1/                 
      DATA MAXF / 1, 4, 10, 20, 35, 56, 84, 4/                 
      DATA NANGM/ 1, 2,  3,  4,  5,  6, 7,  2/                 
      CHARACTER(8) :: LETK
      DATA LETK/'K       '/ 
      CHARACTER(8) :: BASIS
      CHARACTER(10) :: ATOMNM,ENDWRD
      DATA ENDWRD /'$END      '/                                        
!
      INTEGER :: IGTYP
      LOGICAL :: SPRKLE,QMCHKA,QMCHKB,FILE_EXISTS      
      CHARACTER(80) :: BASIS_FILE,PREFIX,BASIS_FILE_1
      COMMON/BASIS_FILE/BASIS_FILE
!
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: CSINP,CPINP,CDINP
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: CFINP,CGINP,CHINP
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: CIINP,EXX,CSS
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: CPP,CDD,SCFAC
      INTEGER,ALLOCATABLE,DIMENSION(:) :: INTYP,NS,KS
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: CINP

      INTEGER :: OFF, OFF_PRIM                     !LIBCINT
      INTEGER :: SIZE_ENV                          !LIBCINT
      DOUBLE PRECISION :: ENV(SIZE_ENV)            !LIBCINT
      INTEGER :: ATM(6,NATmax), BAS(8,NSHELLmax)   !LIBCINT
      DOUBLE PRECISION,EXTERNAL :: CINTgto_norm    !LIBCINT

      COMMON/NSHELaux/KATOMaux(700),KTYPEaux(700),KLOCaux(700),         &
                      KSTARTaux(700),KNGaux(700)
      COMMON/EXCaux/EXaux(2000),Caux(2000)
!-----------------------------------------------------------------------
      ALLOCATE(INTYP(NPRIMImax),NS(NATmax),KS(NATmax))  
      ALLOCATE(EXX(6),CSS(6),CPP(6),CDD(6),SCFAC(4))
      ALLOCATE(CSINP(NPRIMImax),CPINP(NPRIMImax),CDINP(NPRIMImax))
      ALLOCATE(CFINP(NPRIMImax),CGINP(NPRIMImax))
      ALLOCATE(CHINP(NPRIMImax),CIINP(NPRIMImax))
      IF(ILIBCINT==1)ALLOCATE(CINP(NPRIMImax))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      NT = 1                         ! No symmetry point group - > C1
      PI = 2.0d0*DASIN(1.0d0)
      PI32 = PI*SQRT(PI)                                              
      CBASIS = BLANK                                                 
      BASIS  = BLANK
      IDUM   = 0                                                          
      IGAUSS = 0                                                        
      NAT    = 0                                                           
      NATIN  = 0                                                         
      NE     = 0                                                            
      NSHELL = 0                                                        
      LOC    = 0                                                           
      NGAUSS = 0                                                        
      DO I = 1,NPRIMImax                                               
       EX(I) = 0.0D0                                                   
       CS(I) = 0.0D0                                                   
       CP(I) = 0.0D0                                                   
       CD(I) = 0.0D0                                                   
       CF(I) = 0.0D0                                                   
       CG(I) = 0.0D0                                                   
       CH(I) = 0.0D0                                                   
       CI(I) = 0.0D0                                                   
       CSINP(I) = 0.0D0                                                
       CPINP(I) = 0.0D0                                                
       CDINP(I) = 0.0D0                                                
       CFINP(I) = 0.0D0                                                
       CGINP(I) = 0.0D0                                                
       CHINP(I) = 0.0D0                                                
       CIINP(I) = 0.0D0                                                
      END DO
      ZNUC = 0.0D0                                                       
      X = 0.0D0                                                          
      Y = 0.0D0                                                          
      Z = 0.0D0                                                          
      SCFAC(1) = 0.0D0                                                   
      SCFAC(2) = 0.0D0                                                   
      SCFAC(3) = 0.0D0                                                   
      SCFAC(4) = 0.0D0                                                   
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                               Read Atoms
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Open Basis Set file if exists
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(LEN_TRIM(BASIS_FILE)>0)THEN
       CALL GETENV( 'HOME', PREFIX ) 
       BASIS_FILE_1 = TRIM(PREFIX)//'/DoNOFsw/basis/'//                 &
                      TRIM(ADJUSTL(BASIS_FILE))//".bas"
       INQUIRE(FILE=BASIS_FILE_1,EXIST=FILE_EXISTS)
       if(FILE_EXISTS)then                            ! in DoNOFsw
        BASIS_FILE = BASIS_FILE_1       
!      DoNOF
       else
        BASIS_FILE_1= TRIM(PREFIX)//'/DoNOF/basis/'//                   &
                      TRIM(ADJUSTL(BASIS_FILE))//".bas"                  
        INQUIRE(FILE=BASIS_FILE_1,EXIST=FILE_EXISTS)                     
        if(FILE_EXISTS)then                           ! in DoNOF         
         BASIS_FILE = BASIS_FILE_1                                       
        else                                                             
         CALL GETENV( 'PWD', PREFIX )                                    
         BASIS_FILE_1 = TRIM(PREFIX)//"/"//                             &
                        TRIM(ADJUSTL(BASIS_FILE))//".bas"                
         INQUIRE(FILE=BASIS_FILE_1,EXIST=FILE_EXISTS)                    
         if(FILE_EXISTS)then                          ! in pwd           
          BASIS_FILE = BASIS_FILE_1                                      
         else                                                            
          BASIS_FILE_1 = TRIM(ADJUSTL(BASIS_FILE))                       
          INQUIRE(FILE=BASIS_FILE_1,EXIST=FILE_EXISTS)                   
          if(FILE_EXISTS)then                         ! in given path    
           BASIS_FILE = BASIS_FILE_1                                     
          else                             ! basis set file not found    
           WRITE(6,*)"Basis File ",TRIM(BASIS_FILE)," does not exist"    
           CALL ABRT                                                     
          endif                                                          
         endif                                                           
        endif                                                            
       endif                                                             
!      Open Basis Set File ( Unit = 50 )                                 
       OPEN(50,FILE=BASIS_FILE,STATUS='UNKNOWN',                        &
               FORM='FORMATTED',ACCESS='SEQUENTIAL')                     
       WRITE(6,'(/1X,A9,/1X,9(1H-)/1X,A80)')'Basis Set',BASIS_FILE       
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      INPATM = 0                                                                       
    1 CONTINUE                                                          
      INPATM = INPATM + 1   
      CALL READAT(ATOMNM,ZNUC,X,Y,Z)  
!
      IF(ADJUSTL(ATOMNM)/=ADJUSTL(ENDWRD))THEN
       NAT = NAT+1                                                       
       IF(NAT+1>NATmax)THEN
        WRITE(6,'(/,1X,A22,I6)')'Stop: NAT+1 > NATmax =',NAT                                 
        CALL ABRT                                                         
       ENDIF
       READ(UNIT=ATOMNM,FMT='(A8,A2)')ANAM(NAT),BNAM(NAT) 
!
       if(LEN_TRIM(BASIS_FILE)>0)then      ! read from file.bas
        REWIND(50)
        CALL FNDATMBASIS(ANAM(NAT),IEOF)   
       endif
!
       NS(NAT) = 0                                                       
       KS(NAT) = NSHELL+1  
       X = X / UNIT
       Y = Y / UNIT                                                     
       Z = Z / UNIT                                                     
       Cxyz(1,NAT) = X                                                      
       Cxyz(2,NAT) = Y                                                      
       Cxyz(3,NAT) = Z                                                      
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       ZEXTRA = ABS(ZNUC) - INT(ABS(ZNUC))                               
       SPRKLE = ZEXTRA>1.0D-05                                        
       IF(.NOT.SPRKLE)THEN                                              
        IF(ZNUC>0.0D0)THEN                                          
         ZAN(NAT) = ZNUC                                             
         NE = NE + INT(ZNUC)                                         
        ELSE                                                           
         ZAN(NAT) = 0.0D0                                             
         ZNUC = ABS(ZNUC)                                            
        END IF                                                         
       ELSE                                                              
        ZAN(NAT) = ZNUC                                                
       END IF                                                            
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       ICNT = 0                                                          
    2  CONTINUE                                                          
       ICNT = ICNT + 1                                                   
       NUCZ = INT(ZNUC)                                                  
       MPCORE = 0
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                        
!      Read Basis from the input file or basisname.bas
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       IEOF = 0                                                       
       IERR = 0 
       if(LEN_TRIM(BASIS_FILE)>0)then      ! read from file.bas
        CALL RDCARD(50,'$DATA 6U',IEOF)  
       else
        CALL RDCARD(5,'$DATA 6U',IEOF)
       endif
       KSIZE = -8
       CALL GSTRNG(CBASIS,KSIZE)
       READ(UNIT=CBASIS,FMT='(A8)')BASIS   
       IF (BASIS==LETK) BASIS=LABEL(1)                           
       IGAUSS = IFIND('NGAUSS  ',IERR)                             
       DO I=1,4                                                   
        SCFAC(I) = RFIND('SCFAC   ',IERR)                           
        IF(IERR/=0)CALL ABRT                          
       END DO
       IF(ICNT==1)THEN                            
        ABASIS(NUCZ,1) = BASIS                                         
        IAGAUS(NUCZ,1) = IGAUSS                                        
       END IF                                                            
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      A blank string means this Atom is done
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
       IF(BASIS==BLANK)THEN                                           
        XS = X                                                            
        YS = Y                                                            
        ZS = Z                                                            
        NAT0 = NAT                                                        
        QMCHKA = .FALSE.                                                  
        QMCHKB = .FALSE.                                                  
        GOTO 1                                                         
       END IF  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
       IF(IGAUSS>30)THEN                                        
        WRITE(6,'(/,1X,A17,I6)')'Stop: IGAUSS > 30',IGAUSS                                                                                
        WRITE(6,100) NAT,INPATM,ATOMNM,ZNUC,X*UNIT,Y*UNIT,Z*UNIT                      
        CALL ABRT                                                      
       END IF                                                            
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -        
       ITYP = 0                                                          
       DO I=1,8                                                         
        IF(BASIS==LABEL(I))ITYP=I                                   
       ENDDO                                                             
       IF(ITYP==0)THEN
!       Stop: Basis is unrecognized
        WRITE(6,'(/1X,A33,A8,1X,A9,I4)')                                &
        'Stop: Illegal basis function type',BASIS,'IGAUSS = ',IGAUSS
        WRITE(6,100) NAT,INPATM,ATOMNM,ZNUC,X/UNIT,Y/UNIT,Z/UNIT                         
        CALL BERROR(4)                                                    
        CALL ABRT                                                         
       ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                       
!      General Basis Set: 
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                       
!      NSHELL- IS THE TOTAL NUMBER OF SHELLS.  A P SHELL MEANS X,Y,Z,
!       D SHELL MEANS XX,YY,ZZ,XY,XY,YZ, AND SO ON FOR F,G,H,I.      
!       CARTESIAN GAUSSIAN SHELLS CONTAIN ALL X**L Y**M Z**N         
!       PRODUCTS WITH L+M+N= CONSTANT, I.E. L+M+N=3 FOR F SHELL.     
!      EX- GAUSSIAN EXPONENTS, FOR EVERY SYMMETRY UNIQUE PRIMITIVE.  
!      CS- THROUGH -CI- ARE S,P,D,F,G,H,I CONTRACTION COEFFICIENTS.  
!       NORMALLY ONLY ONE OF THE -CX- ARRAYS WILL BE NON-ZERO,       
!       FOR ANY GIVEN EXPONENT IN -EX-.  THE EXCEPTION IS "L" SHELLS,
!       WHERE BOTH -CS- AND -CP- WILL HAVE (DIFFERENT) VALUES.       
!                                                                    
!            THE VARIOUS "K"S DEFINE EACH SHELL'S CONTENTS:          
!                                                                    
!      KATOM- TELLS WHICH ATOM THE SHELL IS CENTERED ON, NORMALLY    
!       MORE THAN ONE SHELL EXISTS ON EVERY ATOM.                    
!      KLOC- GIVES THE LOCATION OF THIS SHELL IN THE TOTAL AO BASIS, 
!       PLEASE READ THE EXAMPLE.                                     
!      KSTART- IS THE LOCATION OF THE FIRST EXPONENT AND THE FIRST   
!       CONTRACTION COEFFICIENT CONTAINED IN A PARTICULAR SHELL.     
!       THUS, -KLOC- IS AN AO COUNTER, -KSTART- A PRIMITIVE COUNTER. 
!      KNG- IS THE NUMBER OF GAUSSIANS IN THIS SHELL.  THEIR DATA    
!       ARE STORED CONSECUTIVELY BEGINNING AT THE -KSTART- VALUE.    
!      KTYPE- IS 1,2,3,4,5,6,7 FOR S,P,D,F,G,H,I.  NOTE THAT THE     
!       VALUE STORED IN -KTYPE- FOR AN "L" SHELL IS A 2, SO THAT     
!       BY ITSELF, -KTYPE- CANNOT DISTINGUISH A "P" FROM A "L".      
!       THUS, KTYPE IS ONE HIGHER THAN THE TRUE ANGULAR MOMENTUM.    
!      KMIN- AND -KMAX- ARE THE STARTING AND ENDING INDICES OF THE   
!       SHELL.  THESE ARE DEFINED AS                                 
!                 S    P    D    F   G   H   I   L                   
!          KMIN   1    2    5   11  21  34  57   1                   
!          KMAX   1    4   10   20  35  56  84   4                   
!       SO YOU CAN TELL AN "L" SHELL BY ITS RUNNING FROM 1 TO 4,     
!       NAMELY S,X,Y,Z, WHEREAS A "P" SHELL RUNS 2,3,4 FOR X,Y,Z.    
!       THE TABLE ABOVE IS GENERATED BY WRITING ALL CARTESIAN        
!       PRODUCTS, "MAXIMUM POWERS FIRST", BACK TO BACK:              
!          S,  X,Y,Z,  XX,YY,ZZ,XY,XZ,YZ,                            
!          1   2 3 4    5  6  7  8  9 10                             
!          XXX,YYY,ZZZ,XXY,XXZ,YYX,YYZ,ZZX,ZZY,XYZ, ... G,H,I        
!          11  12  13,  14  15  16  17  18  19  20, ... G,H,I        
!                                                                    
!      AN EXAMPLE, TO TRY TO MAKE THIS CONCRETE, IS A 6-311G(D,P) BASIS 
!      FOR THE MOLECULE CSIH.  JUST THOSE THREE ATOMS, IN THAT ORDER:   
!              S  L  L  L  D    S  L  L  L  L  D    S  S  S  P          
!      KATOM   1  1  1  1  1    2  2  2  2  2  2    3  3  3  3          
!      KNG     6  3  1  1  1    6  6  3  1  1  1    3  1  1  1          
!      KTYPE   1  2  2  2  3    1  2  2  2  2  3    1  1  1  2          
!      KMIN    1  1  1  1  5    1  1  1  1  1  5    1  1  1  2          
!      KMAX    1  4  4  4 10    1  4  4  4  4 10    1  1  1  4          
!      KSTART  1  7 10 11 12   13 19 25 28 29 30   31 34 35 36 (SUM KNG)
!      KLOC    1  2  6 10 14   20 21 25 29 33 37   43 44 45 46          
!                                                                       
!      KLOC- HELPS POINT TO THE RIGHT AO INDEX, E.G. THE D SHELL       
!      OF THE SI ATOM CONTAINS AOS NUMBERED 37,38,39,40,41,42.          
!      KLOC(I) = KLOC(I-1) + KMAX(I) - KMIN(I) + 1                      
!      TOTAL NUMBER OF AOS (NUM IN COMMON -INFOA-) IN THIS EXAMPLE      
!      IS 48, FROM THE HYPOTHETICAL NEXT KLOC OF 46 + 4 - 2 + 1.        
!      CLEARLY -NSHELL- IS 15, THE NUMBER OF COLUMNS GIVEN HERE.        
!                                                                       
!      NOTE THAT THIS EXAMPLE SHOWS YOU HOW TO TELL A -P- FROM A -L-,   
!      EVEN THOUGH -KTYPE- IS 2 FOR EACH.                               
!                                                                       
!      NOTE THAT D SHELLS ALWAYS HAVE 6 MEMBERS, FOR SPHERICAL          
!      HARMONICS ARE NOT TAKEN CARE OF IN THE BASIS (ALWAYS A           
!      CARTESIAN GAUSSIAN BASIS IS SET UP) BUT RATHER AT THE TIME       
!      OF VARYING THE MOS       
!                                                                       
!      IF OUR MOLECULE WAS REALLY CSIH3, WITH C3V SYMMETRY, SO THAT     
!      THE INPUT GAVE ONLY ONE OF THE HYDROGENS, HOW DOES -NSHEL-       
!      CHANGE?  IT WOULD BE EXTENDED BY TWO MORE ATOMS,                 
!              S  S  S  P      S  S  S  P                               
!      KATOM   4  4  4  4      5  5  5  5                               
!      KNG     3  1  1  1      3  1  1  1                               
!      KTYPE   1  1  1  2      1  1  1  2                               
!      KMIN    1  1  1  2      1  1  1  2                               
!      KMAX    1  1  1  4      1  1  1  4                               
!      KSTART 31 34 35 36     31 34 35 36                               
!      KLOC   49 50 51 52     55 56 57 58                               
!                                                                       
!      NOTE THAT SINCE THESE ARE SYMMETRY EQUIVALENT, -KSTART-          
!      POINTS TO THE ORIGINAL GAUSSIAN DETAILS IN -EX- AND -CX-,        
!      BUT THESE ARE ADDITIONAL AOS, SO -KLOC- DOES GO UP.              
!      NSHELL- IS NOW 24, AND -NUM- IS NOW 60.                         
!                                                                       
!      A MOLECULE MAY VERY WELL HAVE MANY HYDROGENS, PERHAPS USING      
!      IDENTICAL BASIS SETS, BUT EVERY DIFFERENT SET OF EQUIVALENT      
!      HYDROGENS GETS SEPARATE STORAGE OF ITS EXPONENTS/CONTRACTION     
!      COEFFICIENTS (STORED AT DIFFERENT -KSTART- VALUES).              
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                       
       NSHELL = NSHELL + 1
       if(IECP>0.or.ILIBCINT==1)then
        IF(NSHELL>600)THEN
         WRITE(6,'(1X,A30,I6,A16)')                                      &
          'Stop: No more than NSHELLmax (',NSHELLmax,') shells allowed'
         CALL ABRT
        ENDIF
       end if
!                                                                        
       KMIN(NSHELL)  = MINF(ITYP)                                          
       KMAX(NSHELL)  = MAXF(ITYP)                                          
       KSTART(NSHELL)= NGAUSS+1                                          
       KATOM(NSHELL) = NAT                                                
       KTYPE(NSHELL) = NANGM(ITYP)                                        
       INTYP(NSHELL) = ITYP                                               
       KNG(NSHELL)   = IGAUSS                                               
       KLOC(NSHELL)  = LOC+1                                               
       NGAUSS = NGAUSS + IGAUSS                                             
       IF(NGAUSS>NPRIMImax)THEN                                               
        WRITE(6,'(1X,A30,I6,A19)')                                      &
         'Stop: No more than NPRIMImax (',NPRIMImax,') Gaussians allowed'                                   
        CALL ABRT                                                        
       ENDIF                                                             
!                                                                        
       LOC = LOC + NBFS(ITYP)                                            
       IF(LOC>8192)THEN                                                  
        WRITE(6,'(1X,A18,I6,A24)')                                      &
        'Stop: No more than',8192,' basis functions allowed'                                  
        CALL ABRT
       END IF                                                         
       K1 = KSTART(NSHELL)                                               
       K2 = K1 + KNG(NSHELL) - 1                                             
       NS(NAT) = NS(NAT) + 1
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -     
!      S, P, D, F, G, H, I, or L Basis Sets
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       IF(ITYP<=8)THEN
        IF(SCFAC(1) <= 0.0D0) SCFAC(1) = 1.0D0                             
        IDUM=0                                                            
        DO K = K1,K2                                                  
         C1 = 0.0D0                                                      
         C2 = 0.0D0                                                      
         IEOF = 0                                                    
         IERR = 0
         if(LEN_TRIM(BASIS_FILE)>0)then      ! read from file.bas
          CALL RDCARD(50,'$DATA 7U',IEOF)  
         else
          CALL RDCARD(5,'$DATA 7U',IEOF)
         endif
         IDUM = IFIND('IDUM    ',IERR)                               
         IF(IERR/=0)CALL ABRT                           
         EX(K) = RFIND('ZETA    ',IERR)                               
         IF(IERR/=0) CALL ABRT                           
         IF(EX(K)==0.0D0) THEN                            
          WRITE(6,'(1X,A10,I4,A39,A4)')'Stop: Atom',NATIN,              &
                ' has 0 exponent for basis function type',LABEL(ITYP)
          CALL ABRT                                      
         END IF                                            
         C1 = RFIND('C1      ',IERR)                                  
         IF(IERR/=0) CALL ABRT                           
         C2 = RFIND('C2      ',IERR)                                  
         IF(IERR/=0) CALL ABRT                           
         IF(C1==0.0D0.and.EX(K)==1.0D0) THEN                    
          WRITE(6,*)'Stop: Contraction Coefficient not found'      
          WRITE(6,100)NAT,INPATM,ATOMNM,ZNUC,X*UNIT,Y*UNIT,Z*UNIT                
          CALL ABRT                                                
         END IF                                                      
         IF(IGAUSS==1)C1=1.0D0                                      
         IF(IGAUSS==1)C2=1.0D0                                      
         EX(K) = EX(K)*SCFAC(1)**2                                    
         IF(ITYP==1) CSINP(K) = C1                                    
         IF(ITYP==2) CPINP(K) = C1                                    
         IF(ITYP==3) CDINP(K) = C1                                    
         IF(ITYP==4) CFINP(K) = C1                                    
         IF(ITYP==5) CGINP(K) = C1                                    
         IF(ITYP==6) CHINP(K) = C1                                    
         IF(ITYP==7) CIINP(K) = C1                                    
         IF(ITYP==8) CSINP(K) = C1                                    
         IF(ITYP==8) CPINP(K) = C2                                    
         CS(K) = CSINP(K)                                               
         CP(K) = CPINP(K)                                               
         CD(K) = CDINP(K)                                               
         CF(K) = CFINP(K)                                               
         CG(K) = CGINP(K)                                               
         CH(K) = CHINP(K)                                               
         CI(K) = CIINP(K)
!        Store unnormalized coefficients for libint
         IF(ILIBCINT==1) CINP(K) = C1
         !CINP(K) = C1 !LIBCINT
        END DO
       END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                          
!      Normalize Primitive Basis Functions
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                          
       DO IG = K1,K2                                                 
        EE = EX(IG)+EX(IG)                                             
        FACS = PI32/(EE*SQRT(EE))                                      
        FACP = 0.5D0*FACS/EE                                          
        FACD = PT75  *FACS/(EE*EE)                                     
        FACF = PT187 *FACS/(EE**3)                                     
        FACG = PT6562*FACS/(EE**4)                                     
        FACH = PT2953*FACS/(EE**5)                                     
        FACI = PT1624*FACS/(EE**6)                                     
        CS(IG) = CS(IG)/SQRT(FACS)                                     
        CP(IG) = CP(IG)/SQRT(FACP)                                     
        CD(IG) = CD(IG)/SQRT(FACD)                                     
        CF(IG) = CF(IG)/SQRT(FACF)                                     
        CG(IG) = CG(IG)/SQRT(FACG)                                     
        CH(IG) = CH(IG)/SQRT(FACH)                                     
        CI(IG) = CI(IG)/SQRT(FACI)                                     
       END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                          
!      Normalize the contracted basis functions to Unity
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                          
       FACS = 0.0D0                                                       
       FACP = 0.0D0                                                       
       FACD = 0.0D0                                                       
       FACF = 0.0D0                                                       
       FACG = 0.0D0                                                       
       FACH = 0.0D0                                                       
       FACI = 0.0D0                                                       
       DO IG = K1,K2                                                 
        DO JG = K1,IG                                              
         EE = EX(IG)+EX(JG)                                          
         FAC = EE*SQRT(EE)                                           
         DUMS = CS(IG)*CS(JG)/FAC                                    
         DUMP = 0.5D0*CP(IG)*CP(JG)/(EE*FAC)                        
         DUMD = PT75  *CD(IG)*CD(JG)/(EE*EE*FAC)                     
         DUMF = PT187 *CF(IG)*CF(JG)/(EE**3*FAC)                     
         DUMG = PT6562*CG(IG)*CG(JG)/(EE**4*FAC)                     
         DUMH = PT2953*CH(IG)*CH(JG)/(EE**5*FAC)                     
         DUMI = PT1624*CI(IG)*CI(JG)/(EE**6*FAC)                     
         IF(IG /= JG) THEN                                         
          DUMS = DUMS+DUMS                                         
          DUMP = DUMP+DUMP                                         
          DUMD = DUMD+DUMD                                         
          DUMF = DUMF+DUMF                                         
          DUMG = DUMG+DUMG                                         
          DUMH = DUMH+DUMH                                         
          DUMI = DUMI+DUMI                                         
         END IF                                                      
         FACS = FACS+DUMS                                            
         FACP = FACP+DUMP                                            
         FACD = FACD+DUMD                                            
         FACF = FACF+DUMF                                            
         FACG = FACG+DUMG                                            
         FACH = FACH+DUMH                                            
         FACI = FACI+DUMI                                            
        END DO
       END DO
       IF(FACS < TM10) THEN                                           
        FACS = 0.0D0                                                      
       ELSE                                                              
        FACS = 1.0D0/SQRT(FACS*PI32)                                     
       END IF                                                            
       IF(FACP < TM10) THEN                                           
        FACP=0.0D0                                                      
       ELSE                                                              
        FACP = 1.0D0/SQRT(FACP*PI32)                                     
       END IF                                                            
       IF(FACD < TM10) THEN                                           
        FACD=0.0D0                                                      
       ELSE                                                              
        FACD = 1.0D0/SQRT(FACD*PI32)                                     
       END IF                                                            
       IF(FACF < TM10) THEN                                           
        FACF=0.0D0                                                      
       ELSE                                                              
        FACF = 1.0D0/SQRT(FACF*PI32)                                     
       END IF                                                            
       IF(FACG < TM10) THEN                                           
        FACG=0.0D0                                                      
       ELSE                                                              
        FACG = 1.0D0/SQRT(FACG*PI32)                                     
       END IF                                                            
       IF(FACH < TM10) THEN                                           
        FACH=0.0D0                                                      
       ELSE                                                              
        FACH = 1.0D0/SQRT(FACH*PI32)                                     
       END IF                                                            
       IF(FACI < TM10) THEN                                           
        FACI=0.0D0                                                      
       ELSE                                                              
        FACI = 1.0D0/SQRT(FACI*PI32)                                     
       END IF                                                            
       DO IG = K1,K2                                                 
        CS(IG) = CS(IG) * FACS                                         
        CP(IG) = CP(IG) * FACP                                         
        CD(IG) = CD(IG) * FACD                                         
        CF(IG) = CF(IG) * FACF                                         
        CG(IG) = CG(IG) * FACG                                         
        CH(IG) = CH(IG) * FACH                                         
        CI(IG) = CI(IG) * FACI                                         
        CSINP(IG) = CSINP(IG) * FACS                                   
        CPINP(IG) = CPINP(IG) * FACP                                   
        CDINP(IG) = CDINP(IG) * FACD                                   
        CFINP(IG) = CFINP(IG) * FACF                                   
        CGINP(IG) = CGINP(IG) * FACG                                   
        CHINP(IG) = CHINP(IG) * FACH                                   
        CIINP(IG) = CIINP(IG) * FACI                                   
       END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                       
       GO TO 2                                                         
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                       
!      Generate equivalent Atoms
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                       
       XS = X                                                            
       YS = Y                                                            
       ZS = Z                                                            
       NAT0 = NAT                                                        
       QMCHKA = .FALSE.                                                  
       QMCHKB = .FALSE.                                                  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                       
       GO TO 1                                                         
      END IF      
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(NGAUSS==0.or.NSHELL==0)THEN    
       WRITE(6,'(32A)')'Stop: No basis functions defined' 
       CALL ABRT                                                      
      END IF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Fill ENV, ATM and BAS for Libcint
!     The specifications can be found in the LIBCINT documentation
!     but are briefly summarized in the following.
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     ENV
!- - - - - - - - - - - - - - - -  
!     (0) X 20          : Free for internal use (do not fill)
!     (X,Y,Z) X NAT     : Coordinates of each atom
!     (Exps) + (Coeffs) : Exponents and Coefficients of the primitives
!         X NSHELL        of a given shell
!- - - - - - - - - - - - - - - -  
!     ATM
!- - - - - - - - - - - - - - - -  
!     ATM(1, IAT)       : Nuclear Charge of I atom
!     ATM(2, IAT)       : Index of the X coordinate of I atom in ENV
!- - - - - - - - - - - - - - - -  
!     BAS
!- - - - - - - - - - - - - - - -  
!     BAS(1, IBAS)      : Atom index where the shell is centered
!     BAS(2, IBAS)      : Angular momentum of the shell
!     BAS(3, IBAS)      : Number of primitives in the shell
!     BAS(4, IBAS)      : Number of coefficientes per primitive
!     BAS(6, IBAS)      : Index of the first exponent of the shell in ENV
!     BAS(7, IBAS)      : Index of the first coeff of the shell in ENV
!- - - - - - - - - - - - - - - -  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(ILIBCINT==1) THEN
        LOC = 0
        ENV = 0.0D0
        ATM = 0
        BAS = 0
        OFF = 20
        OFF_PRIM = 20 + NAT*3
        DO IAT=1,NAT
          OFF = 20 + (IAT-1)*3
          ATM(1,IAT) = ZAN(IAT)
          ATM(2,IAT) = OFF
          ENV(OFF+1) = Cxyz(1,IAT)
          ENV(OFF+2) = Cxyz(2,IAT)
          ENV(OFF+3) = Cxyz(3,IAT)

          NS1 = KS(IAT)
          NS2 = NS1+NS(IAT)-1
          DO ISH = NS1,NS2
            BAS(1,ISH) = IAT - 1
            BAS(2,ISH) = INTYP(ISH) - 1
            BAS(3,ISH) = KNG(ISH)
            BAS(4,ISH) = 1
            BAS(6,ISH) = OFF_PRIM
            BAS(7,ISH) = OFF_PRIM + KNG(ISH)

            I1 = KSTART(ISH)
            I2 = I1+KNG(ISH)-1
            I0 = 0
            DO IG = I1,I2
               I0 = I0 + 1
               E = EX(IG)
              ENV(OFF_PRIM + I0) = EX(IG)
              ENV(OFF_PRIM + KNG(ISH) + I0) = CINP(IG) *                  &
                      CINTgto_norm(BAS(2,ISH), ENV(BAS(6,ISH)+I0))
            END DO
            OFF_PRIM = OFF_PRIM + 2*KNG(ISH)
            L = BAS(2,ISH)
            IF(IGTYP==1) THEN
              LOC = LOC + (L + 1) * (L + 2) / 2
            ELSE IF(IGTYP==2) THEN
              LOC = LOC + 2 * L + 1
            END IF
          END DO
        END DO
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!             Form Transformation Tables for Atoms and Shells
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO IAT = 1,NAT                                               
       NS1 = KS(IAT)-1                                                
       NS2 = NS(IAT)                                                  
       X = Cxyz(1,IAT)                                                   
       Y = Cxyz(2,IAT)                                                   
       Z = Cxyz(3,IAT)                                                   
       XS = X                                                         
       YS = Y                                                         
       ZS = Z                                                         
       DO IT = 1,NT                                              
        NN = 9*(IT-1)                                               
        XP = XS
        YP = YS
        ZP = ZS
        ICTR = -2**20                                               
        DO 1010 I = 1,NAT                                           
         TEST = (XP-Cxyz(1,I))**2+(YP-Cxyz(2,I))**2+(ZP-Cxyz(3,I))**2      
         IF(TEST>TM10 .OR. (QMCHKA.NEQV.QMCHKB)) GO TO 1010
         ICTR = I                                                 
         GO TO 1020                                               
 1010   CONTINUE                                                    
 1020   MAPCTR(IAT,IT) = ICTR                                       
        NS3 = KS(ICTR)-1                                            
        DO ISH = 1,NS2                                         
         MAPSHL(NS1+ISH,IT) = NS3+ISH                             
        END DO
       END DO
      END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Print atomic coordinates
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      WRITE (6,200)                                             
      DO IAT = 1,NAT                                         
       WRITE(6,300)ANAM(IAT),BNAM(IAT),ZAN(IAT),                        &
                   Cxyz(1,IAT),Cxyz(2,IAT),Cxyz(3,IAT)
      END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Print Basis Functions
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --                                                   
      WRITE(6,'(/1X,A16,/1X,16(1H-))')'Atomic Basis Set'                                                 
      WRITE(6,'(1X,A52)')                                               &
       'The contracted primitive functions were unnormalized'                                         
      WRITE(6,'(1X,A58)')                                               &
       'The contracted basis functions are now normalized to unity'       
      WRITE(6,'(/2X,A21,9X,A33)')'Shell Type  Primitive',               &
                                 'Exponent  Contraction Coefficient'                          
      MLDUDF=0                                                          
      MLDNDA=0
!-----------------------------------------------------------------------
      IMINPIR(1)=1                                                      
!-----------------------------------------------------------------------
      DO IAT = 1,NAT                                               
       DO IT = 1,NT                                              
        IF(MAPCTR(IAT,IT) > IAT) GO TO 20                      
       END DO
       WRITE(6,'(/1X,A8,A2)')ANAM(IAT),BNAM(IAT)                                  
       NS1 = KS(IAT)                                                  
       NS2 = NS1+NS(IAT)-1                                            
       MLDNDA=MLDNDA+1                                                
       DO ISH = NS1,NS2                                          
        WRITE(6,*)                                             
        I1 = KSTART(ISH)                                            
        I2 = I1+KNG(ISH)-1                                          
        ITYP = INTYP(ISH)                                           
        DO IG = I1,I2                                          
         GOTO(11,12,13,14,15,16,17,18)ITYP                   
   11    CONTINUE                                                 
         C1=CSINP(IG)                                             
!--------S--------------------------------------------------------------
         ISHPIR(IG)=ISH                                           
         ITIPO(IG)=ITYP                                         
         C1PIR(IG)=C1                                             
         NPRIMI=IG                                                
!-----------------------------------------------------------------------
         WRITE(6,400) ISH,LABEL(ITYP),IG,EX(IG),C1              
         GO TO 19                                               
   12    CONTINUE                                                 
         C1=CPINP(IG)                                             
!--------P--------------------------------------------------------------
         ISHPIR(IG)=ISH                                           
         ITIPO(IG)=ITYP                                         
         C1PIR(IG)=C1                                             
         NPRIMI=IG                                                
!-----------------------------------------------------------------------
         WRITE(6,400) ISH,LABEL(ITYP),IG,EX(IG),C1              
         GO TO 19                                               
   13    CONTINUE                                                 
         C1=CDINP(IG)                                             
!--------D--------------------------------------------------------------
         ISHPIR(IG)=ISH                                           
         ITIPO(IG)=ITYP                                         
         C1PIR(IG)=C1                                             
         NPRIMI=IG                                                
!-----------------------------------------------------------------------
         WRITE(6,400) ISH,LABEL(ITYP),IG,EX(IG),C1              
         GO TO 19                                               
   14    CONTINUE                                                 
         C1=CFINP(IG)                                             
!--------F--------------------------------------------------------------
         ISHPIR(IG)=ISH                                           
         ITIPO(IG)=ITYP                                         
         C1PIR(IG)=C1                                             
         NPRIMI=IG                                                
!-----------------------------------------------------------------------
         WRITE(6,400) ISH,LABEL(ITYP),IG,EX(IG),C1              
         GO TO 19                                               
   15    CONTINUE                                                 
         C1=CGINP(IG)                                             
!--------G--------------------------------------------------------------
         ISHPIR(IG)=ISH                                           
         ITIPO(IG)=ITYP                                         
         C1PIR(IG)=C1                                             
         NPRIMI=IG                                                
!-----------------------------------------------------------------------
         WRITE(6,400) ISH,LABEL(ITYP),IG,EX(IG),C1              
         GO TO 19                                               
   16    CONTINUE                                                 
         C1=CHINP(IG)                                             
!--------H--------------------------------------------------------------
         ISHPIR(IG)=ISH                                           
         ITIPO(IG)=ITYP                                         
         C1PIR(IG)=C1                                             
         NPRIMI=IG                                                
!-----------------------------------------------------------------------
         WRITE(6,400) ISH,LABEL(ITYP),IG,EX(IG),C1              
         GO TO 19                                               
   17    CONTINUE                                                 
         C1=CIINP(IG)                                             
!--------I--------------------------------------------------------------
         ISHPIR(IG)=ISH                                           
         ITIPO(IG)=ITYP                                         
         C1PIR(IG)=C1                                             
         NPRIMI=IG                                                
!-----------------------------------------------------------------------
         WRITE(6,400) ISH,LABEL(ITYP),IG,EX(IG),C1              
         GO TO 19                                               
   18    CONTINUE                                                 
         C1=CSINP(IG)                                             
         C2=CPINP(IG)                                             
!--------L--------------------------------------------------------------
         ISHPIR(IG)=ISH                                           
         ITIPO(IG)=ITYP                                         
         C1PIR(IG)=C1                                             
         C2PIR(IG)=C2                                             
         NPRIMI=IG                                                
!-----------------------------------------------------------------------
         WRITE (6,400) ISH,LABEL(ITYP),IG,EX(IG),C1,C2          
         GO TO 19                                               
          WRITE(6,*) 'These shells do not exist'     
         CALL ABRT                                                
   19    CONTINUE                                                    
        END DO
       END DO
!-----------------------------------------------------------------------
       IMAXPIR(IAT)=I2                                                  
       IMINPIR(IAT+1)=IMAXPIR(IAT)+1                                    
!-----------------------------------------------------------------------
   20  CONTINUE
      END DO
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --                                                          
!     NE: Number of Electrons
!     NA: Number of Alpha Electrons
!     NB: Number of Beta Electrons
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                       
      NUM  = LOC                                                        
      NQMT = LOC                                                        
      NE = NE-ICH                                                       
      NA = (NE+MUL-1)/2                                                 
      NB = (NE-MUL+1)/2 
      IF(NA+NB/=NE)THEN
       WRITE(6,'(/,1X,19A)')'Number of Electrons'      
       WRITE(6,'(/,1X,I6)')NE
       WRITE(6,'(/,1X,34A)')'Check your Charge and Multiplicity'
       WRITE(6,'(/,1X,2I3,/)')ICH,MUL 
       CALL ABRT                                                   
      END IF
      IMULNE=MUL+NE                                                    
      IF(MUL>NE+1 .or. MUL<0 .or. 2*INT(IMULNE/2)==IMULNE)THEN  
       WRITE(6,'(28A)')'Impossible Spin Multiplicity'     
       CALL ABRT                                                    
      END IF                                                         
!-----------------------------------------------------------------------
      DO I=1,NSHELL                                                     
       INTIPO(I)=INTYP(I)                                             
      ENDDO
!-----------------------------------------------------------------------
      IF(IGTYP==1) WRITE(6,500)NSHELL,NPRIMI,NUM
      IF(IGTYP==2) WRITE(6,501)NSHELL,NPRIMI,NUM
      NSHELLaux = 0
      NPRIMIaux = 0
      IF (IERITYP==2 .or. IERITYP==3) THEN
       if(ILIBCINT==0)then
        IF(IRITYP==1) THEN
          CALL AUXREAD(IGTYP,NAT,NSHELLaux,NUMaux,NATmax,ANAM)       
        ELSE IF(IRITYP==2) THEN
          CALL AUXGEN(IGTYP,NAT,NPRIMI,ITIPO,IMINPIR,IMAXPIR,NSHELLaux, &
                 NUMaux,IGEN,ISTAR,EX,ZAN)       
        END IF
       else if(ILIBCINT==1)then
        IF(IRITYP==1) THEN
          CALL AUXREADlib(IGTYP,NAT,NSHELLaux,NUMaux,NATmax,ANAM)       
        ELSE IF(IRITYP==2) THEN
          CALL AUXGENlib(IGTYP,NAT,NPRIMI,ITIPO,IMINPIR,IMAXPIR,        &
                  NSHELLaux,NUMaux,EX,ZAN,Cxyz)
        END IF
       end if
      IF(IGTYP==1) WRITE(6,550)NSHELLaux,NUMaux
      IF(IGTYP==2) WRITE(6,551)NSHELLaux,NUMaux
      END IF
      WRITE(6,600)NE,ICH,MUL,NA,NB,NAT          
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Fill ENV, and BAS for Libcint with Auxiliary basis
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(ILIBCINT==1) THEN
        I1 = 0
        DO I=1,NSHELLaux
          IBAS = NSHELL + I
          BAS(1,IBAS) = KATOMaux(I) - 1
          BAS(2,IBAS) = KTYPEaux(I) - 1
          BAS(3,IBAS) = KNGaux(I)
          BAS(4,IBAS) = 1
          BAS(6,IBAS) = OFF_PRIM
          BAS(7,IBAS) = OFF_PRIM + KNGaux(I)
          DO I0=1,KNGaux(I)
          I1 = I1 + 1
          ENV(OFF_PRIM + I0) = EXaux(I1)
          ENV(OFF_PRIM + KNGaux(I) + I0) = Caux(I1) *                    &
                         CINTgto_norm(BAS(2,IBAS), ENV(BAS(6,IBAS)+I0))
          END DO
          OFF_PRIM = OFF_PRIM + 2*KNGaux(I)
        END DO
        NPRIMIaux = I1 
      END IF 
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Nuclear Energy with true nuclear charges
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      EN = ENUC(NAT,ZAN,Cxyz)            
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     IECP /=0 -> Write a warning with respect to NE and ENUC
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(IECP/=0)WRITE(6,'(/1X,A16,F20.10/,/1X,A59,/1X,A63)')           &
      'Nuclear Energy = ',EN,                                           &
      'Note: This calculation is using an effective core potential',    &
      'so NE and ENUC will be adjusted later after removal of the Core'
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DEALLOCATE(CSINP,CPINP,CDINP,CFINP,CGINP,CHINP,CIINP)
      DEALLOCATE(INTYP,NS,KS,EXX,CSS,CPP,CDD,SCFAC)
      IF(ILIBCINT==1)DEALLOCATE(CINP)      
!-----------------------------------------------------------------------                                                                       
  100 FORMAT(1X,'ERROR OCCURED READING ATOM NO.',I4,                    &
                ', INPUT ATOM NO.',I5,', NAME=',A10/                    &
             1X,'CHARGE=',F5.1,' X,Y,Z=',3F15.8)                         
  200 FORMAT(/1X,'Atom',6X,'Charge',16X,'Coordinates (Bohr)'/           &
               27X,'x',13X,'y',13X,'z')                                  
  300 FORMAT(1X,A8,A2,F5.1,3F14.4)                                       
  400 FORMAT(1X,I6,3X,A2,I7,F22.7,2F18.12)                               
  500 FORMAT(/1X,'Total Number of Basis Set Shells             =',I5/   &
              1X,'Total Number of Primitives                   =',I5/   &
              1X,'Number of Cartesian Gaussian Basis Functions =',I5)
  501 FORMAT(/1X,'Total Number of Basis Set Shells             =',I5/   &
              1X,'Total Number of Primitives                   =',I5/   &
              1X,'Number of Spherical Gaussian Basis Functions =',I5)
  550 FORMAT( 1X,'Total Number of Auxiliary Basis Set Shells   =',I5/   &
              1X,'Number of Cart. Gaussian Auxiliary Set Func. =',I5)
  551 FORMAT( 1X,'Total Number of Auxiliary Basis Set Shells   =',I5/   &
              1X,'Number of Sph. Gaussian Auxiliary Set Func.  =',I5)
  600 FORMAT( 1X,'Number of Electrons                          =',I5/   &
              1X,'Charge of Molecule                           =',I5/   &
              1X,'Spin Multiplicity                            =',I5/   &
              1X,'Number of Occupied Orbitals (Alpha)          =',I5/   &
              1X,'Number of Occupied Orbitals (Beta )          =',I5/   &
              1X,'Total Number of Atoms                        =',I5) 
!-----------------------------------------------------------------------  
      RETURN                                                            
      END                                                               

! FNDATMBASIS
      SUBROUTINE FNDATMBASIS(ATOMNAME,IEOF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*8 PREVWORD,WORD,ATOMNAME
      IEOF = 0
      WORD = "        "
    1 CONTINUE
      PREVWORD = WORD
      READ(50,'(A2)',END=2)WORD
      IF(WORD==ATOMNAME .and. LEN(TRIM(PREVWORD))==0)GO TO 3
      CALL UPRCAS(WORD,8)
      IF(WORD==ATOMNAME .and. LEN(TRIM(PREVWORD))==0)GO TO 3
      GOTO 1
    2 CONTINUE
      IEOF = 1
    3 CONTINUE
      RETURN
      END
            
! READAT                                           
      SUBROUTINE READAT(ATOMNM,ZNUC,X,Y,Z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      COMMON/INFO  /IUNTRD            
      COMMON/INTNAL/NATIN            
      CHARACTER*10 ATOMNM,ENDWRD,BLANK10                                  
      DATA ENDWRD,BLANK10/'$END      ','          '/
!-----------------------------------------------------------------------
!     Read Cartesian Coordinates: READ(5,*) ATOMNM,ZNUC,X,Y,Z                                      
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                       
      DFACT = 1.0D0                                                       
      IF(IUNTRD==-1)DFACT = 0.52917724924D0                                     
!                                         
      IEOF = 0                                                          
      IERR = 0                                                          
      CALL RDCARD(5,'$DATA 5U',IEOF)                                      
      IF(IEOF==1) CALL ABRT                                        
!                                                                       
      KSIZE = -10                                                       
      CALL GSTRNG(ATOMNM,KSIZE)                                         
      IF (ADJUSTL(ATOMNM)==ADJUSTL(ENDWRD)) RETURN                                      
!                                                                       
      NATIN = NATIN+1                                                   
      IF (ATOMNM==BLANK10) THEN                                         
         WRITE(6,1) NATIN                               
         CALL ABRT                                                      
      END IF                                                            
!                                                                       
      ZNUC = RFIND('ZNUC    ',IERR)                                     
             IF(IERR/=0) CALL ABRT                                    
      X = DFACT*RFIND('X       ',IERR)                                  
             IF(IERR/=0) CALL ABRT                                    
      Y = DFACT*RFIND('Y       ',IERR)                                  
             IF(IERR/=0) CALL ABRT                                    
      Z = DFACT*RFIND('Z       ',IERR)                                  
             IF(IERR/=0) CALL ABRT
!-----------------------------------------------------------------------
    1 FORMAT(//1X,'*** ERROR!'/                                         &
               1X,'BLANK CARD FOUND WHILE TRYING TO READ INPUT ATOM',I5/&     
               1X,'POSSIBLE ERRORS INCLUDE:'/                           &     
               1X,'1. C1 GROUP SHOULD NOT HAVE A BLANK CARD AFTER IT.'/ &     
               1X,'2. BOTH $BASIS GROUP AND BASIS SET IN $DATA GIVEN?'/ &     
               1X,'3. EXTRANEOUS BLANK CARDS IN $DATA?')                     
!-----------------------------------------------------------------------                                    
      RETURN                                                            
      END                                                               

! SETLAB
      SUBROUTINE SETLAB(KATOM,KMIN,KMAX,NSHELL,ZAN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/INFOA/NAT,ICH,MUL,NUM,NQMT,NE,NA,NB
      DOUBLE PRECISION,DIMENSION(NAT) :: ZAN
      COMMON/RUNLAB/BFLAB(8192)
      INTEGER,DIMENSION(NSHELL) :: KATOM,KMIN,KMAX
      CHARACTER*4 :: LABEL                                               
      CHARACTER*8 :: BFL                                                 
      CHARACTER*4, DIMENSION(35) :: BFNAM1                               
      DATA BFNAM1/'  S ','  X ','  Y ','  Z ',                          &
                  ' XX ',' YY ',' ZZ ',' XY ',' XZ ',' YZ ',            &
                  ' XXX',' YYY',' ZZZ',' XXY',' XXZ',                   &
                  ' YYX',' YYZ',' ZZX',' ZZY',' XYZ',                   &
                  'XXXX','YYYY','ZZZZ','XXXY','XXXZ',                   &
                  'YYYX','YYYZ','ZZZX','ZZZY','XXYY',                   &
                  'XXZZ','YYZZ','XXYZ','YYXZ','ZZXY'/                    
      CHARACTER*6, DIMENSION(49) :: BFNAM2                               
      DATA BFNAM2/' XXXXX',' YYYYY',' ZZZZZ',' XXXXY',' XXXXZ',         &
                  ' YYYYX',' YYYYZ',' ZZZZX',' ZZZZY',' XXXYY',         &
                  ' XXXZZ',' YYYXX',' YYYZZ',' ZZZXX',' ZZZYY',         &
                  ' XXXYZ',' YYYXZ',' ZZZXY',' XXYYZ',' XXZZY',         &
                  ' YYZZX',                                             &
                  '    X6','    Y6','    Z6','   X5Y','   X5Z',         &
                  '   Y5X','   Y5Z','   Z5X','   Z5Y','  X4Y2',         &
                  '  X4Z2','  Y4X2','  Y4Z2','  Z4X2','  Z4Y2',         &
                  '  X4YZ','  Y4XZ','  Z4XY','  X3Y3','  X3Z3',         &
                  '  Y3Z3',' X3Y2Z',' X3Z2Y',' Y3X2Z',' Y3Z2X',         &
                  ' Z3X2Y',' Z3Y2X','X2Y2Z2'/                            
      CHARACTER*4 :: BONDF                                               
      DATA BONDF/' BF '/                                                 
      CHARACTER*4,DIMENSION(106) :: ATMLAB                               
      DATA ATMLAB/'H   ','HE  ','LI  ','BE  ','B   ','C   ',            &
                  'N   ','O   ','F   ','NE  ','NA  ','MG  ',            &
                  'AL  ','SI  ','P   ','S   ','CL  ','AR  ',            &
                  'K   ','CA  ','SC  ','TI  ','V   ','CR  ',            &
                  'MN  ','FE  ','CO  ','NI  ','CU  ','ZN  ',            &
                  'GA  ','GE  ','AS  ','SE  ','BR  ','KR  ',            &
                  'RB  ','SR  ','Y   ','ZR  ','NB  ','MO  ',            &
                  'TC  ','RU  ','RH  ','PD  ','AG  ','CD  ',            &
                  'IN  ','SN  ','SB  ','TE  ','I   ','XE  ',            &
                  'CS  ','BA  ','LA  ','CE  ','PR  ','ND  ',            &
                  'PM  ','SM  ','EU  ','GD  ','TB  ','DY  ',            &
                  'HO  ','ER  ','TM  ','YB  ','LU  ','HF  ',            &
                  'TA  ','W   ','RE  ','OS  ','IR  ','PT  ',            &
                  'AU  ','HG  ','TL  ','PB  ','BI  ','PO  ',            &
                  'AT  ','RN  ','FR  ','RA  ','AC  ','TH  ',            &
                  'PA  ','U   ','NP  ','PU  ','AM  ','CM  ',            &
                  'BK  ','CF  ','ES  ','FM  ','MD  ','NO  ',            &
                  'LR  ','RF  ','X   ','BQ  '/
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Basis Function Table
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      N = 0
      DO III = 1,NSHELL
       IAT = KATOM(III)
       J = INT(ZAN(IAT))
       IF(J<=0)THEN
        LABEL = BONDF
       ELSE
        IF(J>106)J = 106
        LABEL = ATMLAB(J)
       ENDIF
       MINI = KMIN(III)
       MAXI = KMAX(III)
       DO I = MINI,MAXI
        N = N+1
        IF(I<=35) THEN
         WRITE(UNIT=BFL,FMT='(A2,I2,A4)')LABEL,MOD(IAT,100),BFNAM1(I)
        ELSE
         WRITE(UNIT=BFL,FMT='(A2,A6)')LABEL,BFNAM2(I-35)
        END IF
        READ(UNIT=BFL,FMT='(A8)')BFLAB(N)
       END DO
      END DO
!-----------------------------------------------------------------------
      RETURN
      END

! AMT
      SUBROUTINE AMT(ZAN,ZMASS,NAT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DOUBLE PRECISION, DIMENSION(106) :: AMS
      DOUBLE PRECISION, DIMENSION(NAT) :: ZAN,ZMASS,AMASS
!-----------------------------------------------------------------------
!     Mass of most abundant Isotopes
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DATA (AMS(I),I=1,106)  /                                          &
         1.007825D0,4.0026D0,7.01600D0,9.01218D0,11.00931D0,            &
         12.0D0,14.00307D0,15.99491D0,18.99840D0,19.99244D0,            &
         22.9898D0,23.98504D0,26.98153D0,27.97693D0,                    &
         30.97376D0,31.97207D0,34.96885D0,39.948D0,                     &
         38.96371D0,39.96259D0,44.95592D0,47.90D0,50.9440D0,            &
         51.9405D0,54.9381D0,55.9349D0,58.9332D0,57.9353D0,             &
         62.9298D0,63.9291D0,68.9257D0,73.9219D0,74.9216D0,             &
         79.9165D0,78.9183D0,83.9115D0,                                 &
         84.9117D0,87.9056D0,89.9054D0,89.9043D0,92.9060D0,             &
         97.9055D0,97.0D0,101.9037D0,102.9048D0,105.9032D0,             &
         106.9041D0,113.9036D0,114.9041D0,119.9022D0,                   &
         120.9038D0,129.9067D0,126.9044D0,131.9042D0,                   &
         132.9054D0,137.9052D0,138.9063D0,139.9054D0,                   &
         140.9076D0,141.9077D0,144.9127D0,151.9197D0,                   &
         152.9212D0,157.9241D0,158.9253D0,163.9292D0,                   &
         164.9303D0,165.9303D0,168.9342D0,173.9389D0,                   &
         174.9408D0,179.9465D0,180.9480D0,183.9509D0,                   &
         186.9557D0,191.9615D0,192.9629D0,194.9648D0,                   &
         196.9665D0,201.9706D0,                                         &
         204.9744D0,207.9766D0,208.9804D0,208.9824D0,                   &
         209.9871D0,222.0176D0,                                         &
         223.0197D0,226.0254D0,                                         &
         227.0278D0,232.0381D0,231.0359D0,238.0508D0,                   &
         237.0482D0,244.0642D0,243.0614D0,247.0703D0,                   &
         247.0703D0,251.0796D0,252.0829D0,257.0751D0,                   &
         258.0986D0,259.1009D0,260.1053D0,261.1087D0,                   &
         2*0.0D0/
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Note: ZAN is the True Charge because this input is before ECP
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO IAT = 1,NAT
       ZNUC = ZAN(IAT)
       NUCZ = INT(ZNUC)
       IF(NUCZ>=1.and.NUCZ<=106)THEN
        AMASS(IAT) = AMS(NUCZ)
       ELSE
        AMASS(IAT) = 0.0D0
       END IF
      END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Returns Normal Masses
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO I=1,NAT
       ZMASS(I) = AMASS(I)
      END DO
!-----------------------------------------------------------------------
      RETURN
      END

! INRTIA
      SUBROUTINE INRTIA(C,COM,ZMASS,VMOI,NPART)
      IMPLICIT NONE
      DOUBLE PRECISION :: Xcm,Ycm,Zcm
      COMMON/CMCoord/Xcm,Ycm,Zcm
!     ARGUMENTS
      INTEGER,INTENT(IN) :: NPART
      DOUBLE PRECISION,DIMENSION(3,NPART),INTENT(IN) :: C
      DOUBLE PRECISION,DIMENSION(3,NPART),INTENT(OUT) :: COM
      DOUBLE PRECISION,DIMENSION(3),INTENT(OUT) :: VMOI
      DOUBLE PRECISION,DIMENSION(NPART),INTENT(IN) :: ZMASS
!     VARIABLES
      DOUBLE PRECISION,DIMENSION(3):: CMASS,WRK
      DOUBLE PRECISION,DIMENSION(3,3):: TROT,TMOIsq
      DOUBLE PRECISION,DIMENSION(6):: TMOI
      DOUBLE PRECISION :: XX,YY,ZZ,XY,XZ,YZ,WEIGHT,TOTWT,XC,YC,ZC
      INTEGER :: I,J
!-----------------------------------------------------------------------
!     C: Nuclear Coordinates
!     ZMASS: Nuclear Masses
!     COM: Center of Mass Coordinates
!     VMOI: Products of Inertia (Eigenvalues)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Center of Mass
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      TOTWT=0.0D0
      CMASS(1)=0.0D0
      CMASS(2)=0.0D0
      CMASS(3)=0.0D0
      DO I=1,NPART
       WEIGHT=ZMASS(I)
       TOTWT=TOTWT+WEIGHT
       DO J=1,3
        CMASS(J)=CMASS(J)+WEIGHT*C(J,I)
       END DO
      END DO
      DO I=1,3
       CMASS(I)=CMASS(I)/TOTWT
      END DO
!
      Xcm = CMASS(1)
      Ycm = CMASS(2)
      Zcm = CMASS(3)
!
      DO I=1,NPART
       DO J=1,3
        COM(J,I) = C(J,I)-CMASS(J)
       END DO
      END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Moment of Inertia Tensor
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      XX=0.0D0
      YY=0.0D0
      ZZ=0.0D0
      XY=0.0D0
      XZ=0.0D0
      YZ=0.0D0
      DO I=1,NPART
       WEIGHT=ZMASS(I)
       XC = COM(1,I)
       YC = COM(2,I)
       ZC = COM(3,I)
       XX = XX + WEIGHT*(YC*YC+ZC*ZC)
       YY = YY + WEIGHT*(XC*XC+ZC*ZC)
       ZZ = ZZ + WEIGHT*(XC*XC+YC*YC)
       XY = XY - WEIGHT*XC*YC
       XZ = XZ - WEIGHT*XC*ZC
       YZ = YZ - WEIGHT*YC*ZC
      END DO
      TMOI(1) = XX
      TMOI(2) = XY
      TMOI(3) = YY
      TMOI(4) = XZ
      TMOI(5) = YZ
      TMOI(6) = ZZ
!
      CALL CPYTSQ(TMOI,TMOIsq,3)            
      CALL DIAG(3,TMOIsq,TROT,VMOI,WRK)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END

! ENUC                                             
      DOUBLE PRECISION FUNCTION ENUC(N,Z,Cxyz)  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      LOGICAL EFLDL                                                  
      COMMON/EFLDC_1/EFLDL
      COMMON/EFLDC_2/EVEC(3) 
      DIMENSION Z(N),Cxyz(3,N)                                             
!-----------------------------------------------------------------------                                                        
      REPNUC = 0.0d0                                                     
!                                                        
      IF(N/=1)THEN
       DO I = 2,N                                             
        NI = I-1                                                    
        DO J = 1,NI                                         
         RR = 0.0d0                                                
         DO K = 1,3                                           
          RR = RR+(Cxyz(K,I)-Cxyz(K,J))**2                            
         END DO
         IF(RR/=0.0d0) REPNUC = REPNUC + Z(I)*Z(J)/SQRT(RR)        
        END DO
       END DO
      END IF                                                        
!
      IF(EFLDL)THEN                                                   
       ANUCF = 0.0d0                                                   
       DO J = 1,N                                                 
        DO I = 1,3                                              
         ANUCF = ANUCF - EVEC(I)*Cxyz(I,J)*Z(J)                      
        END DO
       END DO
       REPNUC = REPNUC + ANUCF                                        
      END IF 
!                                                           
      ENUC = REPNUC 
!-----------------------------------------------------------------------                                                                                                            
      RETURN                                                            
      END                                                               

! BASCHK                                           
      SUBROUTINE BASCHK(LMAXIMA,KTYPE,NSHELL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      INTEGER,DIMENSION(NSHELL) :: KTYPE
!-----------------------------------------------------------------------
!     Return the highest angular momentum in the basis
!-----------------------------------------------------------------------                                                     
      KANG = 0                                                          
      DO N=1,NSHELL                                                 
       IF(KTYPE(N)>KANG)KANG = KTYPE(N)                          
      ENDDO
      LMAXIMA = KANG-1
!-----------------------------------------------------------------------                                                     
      RETURN                                                            
      END                                                               

! BERROR
      SUBROUTINE BERROR(K)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      CHARACTER(8), DIMENSION(4) :: BASIS
      DATA BASIS /'MINIMAL ','EXTENDED','GENERAL ','        '/
!-----------------------------------------------------------------------
      WRITE(6,1)BASIS(K)
      CALL ABRT
      RETURN
!-----------------------------------------------------------------------
    1 FORMAT(/2X,'Stop: illegal ',A8,' basis function')
!-----------------------------------------------------------------------
      END
      
! FNDGRP
      SUBROUTINE FNDGRP(LUIN,GRPNAM,IEOF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      CHARACTER*8 WORD,GRPNAM
      IEOF = 0
    1 CONTINUE
      READ(LUIN,'(A8)',END=2)WORD
      CALL UPRCAS(WORD,8)
      IF(ADJUSTL(WORD)==ADJUSTL(GRPNAM))GO TO 3
      GOTO 1
    2 CONTINUE
      IEOF = 1
    3 CONTINUE
      RETURN
      END

! UPRCAS
      SUBROUTINE UPRCAS(STRING,LENSTR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      CHARACTER*(*) STRING
      CHARACTER*26 UCASE,LCASE
      DATA UCASE /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      DATA LCASE /'abcdefghijklmnopqrstuvwxyz'/
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Lower Case -> Upper Case
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO I=1,LENSTR
       IC = INDEX(LCASE,STRING(I:I))
       IF (IC>0) STRING(I:I) = UCASE(IC:IC)
      END DO
      RETURN
      END

! OPNCRD
      SUBROUTINE OPNCRD
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      CHARACTER*1 LCONT,LEOD,LEOC
      CHARACTER*80 CARD
      COMMON /FREFM1/ NCOL,LSTCOL,MAXCOL,KEOF,KERR,LUOUT,LUERR,KOLSV
      COMMON /FREFM2/ CARD,LCONT,LEOD,LEOC
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Free format: Read from 5, output to 6
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      MAXCOL= 80
      KEOF  = 0
      KERR  = 0
      LUOUT = 6
      LUERR = 6
      LCONT = '>'
      LEOD  = '!'
      LEOC  = ';'
      RETURN
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -      
      END

! RDCARD
      SUBROUTINE RDCARD(LUIN,ROUTIN,IEOF)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL FIRST
      CHARACTER*8 ROUTIN,ROUTNE,STAR,STARS
      CHARACTER*4 KEYWRD,KCOL,KSET,KEOD,KEOC
      CHARACTER*1 LCONT,LEOD,LEOC
      CHARACTER*80 CARD
      COMMON /FREFM1/ NCOL,LSTCOL,MAXCOL,KEOF,KERR,LUOUT,LUERR,KOLSV
      COMMON /FREFM2/ CARD,LCONT,LEOD,LEOC
      CHARACTER(80) :: BASIS_FILE
      COMMON/BASIS_FILE/BASIS_FILE
      DATA FIRST/.TRUE./
      SAVE FIRST,MORE,STAR,STARS,KCOL,KSET,KEOD,KEOC,LASTC
!-----------------------------------------------------------------------
      IZERO = 0
      IONE = 1
      KEOF = IEOF
      IF(FIRST) THEN
       MORE=0
       STAR='   *    '
       STARS='   **   '
       KCOL='COL '
       KSET='SET '
       KEOD='EOD '
       KEOC='EOC '
       FIRST=.FALSE.
      END IF
      IF(KEOF==1) GO TO 170
         ROUTNE=ROUTIN
         KPRNT=1
         IF(KEOF==(-1)) GO TO 160
            KPRNT=2
            IF(MORE==1) GO TO 120
  105          KPRNT=0
               NSTART=MAXCOL+1
               LSTCOL=MAXCOL
               KEOF=1
               READ(LUIN,900,END=180) CARD
!              IF(LUOUT>0)WRITE(LUOUT,901) ROUTNE,CARD
               IF(ROUTNE==STAR) GO TO 118
!
!                   RECURSIVE CALLS TO GSTRNG AND DECODN WILL BE AVOIDED
!
                  NSET=0
                  NCOL=1
  112             CONTINUE
                  KEYWRD='   *'
                  LGSTR=-4
                  CALL GSTRNG(KEYWRD,LGSTR)
                  IF(KEYWRD/=KSET) GO TO 117
                     LGSTR=-4
                     CALL GSTRNG(KEYWRD,LGSTR)
                     IF(KEYWRD/=KCOL) GO TO 113
                        IERR=0
                        MAXCOL=IFIND('MAXCOL  ',IERR)
                        NSET=NSET+1
                        GO TO 112
!
  113                IF(KEYWRD/=KEOD) GO TO 114
                        LGSTR=-4
                        CALL GSTRNG(LEOD,LGSTR)
                        NSET=NSET+1
                        GO TO 112
!
  114                IF(KEYWRD/=KEOC) GO TO 115
                        LGSTR=-4
                        CALL GSTRNG(LEOC,LGSTR)
                        NSET=NSET+1
                        GO TO 112
!
  115                IF(KEYWRD/=LCONT) GO TO 117
                        LGSTR=-4
                        CALL GSTRNG(LCONT,LGSTR)
                        NSET=NSET+1
                        GO TO 112
!
  117             IF(NSET==0) GO TO 118
                     IF(LUOUT>0) WRITE(LUOUT,904)                       &
                        NSET,MAXCOL,LEOD,LEOC,LCONT
                     ROUTNE=STARS
                     GO TO 105
!
  118          CONTINUE
               LASTC=0
  120       CONTINUE
!
!           CHECK FOR DATA FIELD TERMINATION
!
            MORE=1
            NSTART=LASTC+1
            LASTC=MAXCOL+1
            DO 130 N=NSTART,MAXCOL
               IF(CARD(N:N)==LEOD) LASTC=N
               IF(CARD(N:N)==LEOC) GO TO 150
  130       CONTINUE
            N=LASTC
            MORE=0
  150       LSTCOL=MIN(LASTC,N)-1
            LASTC=N
  160    CONTINUE
         IF(KPRNT==1 .AND. LUOUT>0)                                     &
            WRITE(LUOUT,902) ROUTNE                                      
         IF(KPRNT==2 .AND. LUOUT>0)                                     &
            WRITE(LUOUT,903) NSTART,LSTCOL,ROUTNE                        
         KEOF=0                                                          
  170 CONTINUE                                                           
      NCOL=NSTART                                                        
      KOLSV=NCOL                                                         
      IEOF = KEOF                                                        
      RETURN                                                             
  180 CONTINUE                                                           
      NCOL=NSTART                                                        
      KOLSV=NCOL                                                         
      IEOF = KEOF                                                        
      RETURN
!                                                                        
  900 FORMAT(A80)                                                        
  902 FORMAT(10X,10HREREAD AT ,A8)                                       
  903 FORMAT(10X,8HCOLUMNS ,I2,3H - ,I2,1X,8HREAD AT ,A8)                
  904 FORMAT(/6X,I2,' CARD PARAMETERS HAVE BEEN RESET --'/1X,           &
             7H  COL =,I3,3X,7H EOD = ,A4,7H EOC = ,A4,8H CONT = ,A4/)
      END

! RFIND
      FUNCTION RFIND(VARABL,IERR)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*8 VARABL
      COMMON /FREFM1/ NCOL,LSTCOL,MAXCOL,KEOF,KERR,LUOUT,LUERR,KOLSV
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      KOLSV = NCOL
      KERR = IERR
      CALL DECODN(VARABL,VALUE,ISIGN,FIXPNT,FRACT,IEXPFR,ISIGNE,IEXP)
      RFIND = VALUE
      IERR = KERR
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -      
      RETURN
      END

! DECODN
      SUBROUTINE DECODN(VARABL,VALUE,ISIGN,FIXPNT,FRACT,IEXPFR,ISIGNE,  &
                        IEXP)                                            
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                 
      LOGICAL FIRST                                                      
      CHARACTER*8 VARABL                                                 
      CHARACTER*1 BLANK,PLUS,DECPT,SLASH,COMMA,EQUALS,RPAREN,TAB,       &
                  MINUS,LPAREN,LETE,LETD,LETS,LETQ,LETR,LETT,           &
                  NUM(10),LETTER                                         
      CHARACTER*1 LCONT,LEOD,LEOC                                        
      CHARACTER*80 CARD                                                  
      COMMON /FREFM1/ NCOL,LSTCOL,MAXCOL,KEOF,KERR,LUOUT,LUERR,KOLSV     
      COMMON /FREFM2/ CARD,LCONT,LEOD,LEOC                               
      CHARACTER(80) :: BASIS_FILE                                        
      COMMON/BASIS_FILE/BASIS_FILE                                       
      DATA FIRST/.TRUE./                                                 
      SAVE FIRST,TAB,NUM,BLANK,PLUS,MINUS,DECPT,SLASH,COMMA,EQUALS,     &
           LPAREN,RPAREN,LETE,LETD,LETS,LETQ,LETR,LETT
!-----------------------------------------------------------------------
      IF(FIRST) THEN
       TAB = CHAR(9)
       NUM(1) = '0'
       NUM(2) = '1'
       NUM(3) = '2'
       NUM(4) = '3'
       NUM(5) = '4'
       NUM(6) = '5'
       NUM(7) = '6'
       NUM(8) = '7'
       NUM(9) = '8'
       NUM(10)= '9'
       BLANK = ' '
       PLUS =  '+'
       MINUS = '-'
       DECPT = '.'
       SLASH = '/'
       COMMA = ','
       EQUALS= '='
       LPAREN= '('
       RPAREN= ')'
       LETE =  'E'
       LETD  = 'D'
       LETS  = 'S'
       LETQ  = 'Q'
       LETR  = 'R'
       LETT  = 'T'
       FIRST = .FALSE.
      END IF
      TNUM=0.0D0
      ISQRT=0
  110 LSIGN=1
      ICOMMA=0
      NEWCRD=0
      X=0.0D0
      LSTERR=KERR
      KERR = 0
      NDECPT=0
      KINDST=1
      ISIGN=1
      ISIGNE=1
      VALUE=0.0D0
      FIXPNT=0.0D0
      FRACT=0.0D0
      IEXPFR=0
      IEXP=0
      LEND=0
      LAST=KOLSV
!
  120 LFIRST=0
  130 IF (NCOL>LSTCOL) GO TO 240
      LETTER=CARD(NCOL:NCOL)
      CALL UPRCAS(LETTER,1)
      IF(LETTER==PLUS  ) GO TO 170
      IF(LETTER==MINUS ) GO TO 190
      IF((LFIRST==0  .AND.  LAST<=KOLSV)  .AND.                         &
         (LETTER==LETD  .OR.  LETTER==LETE)) GO TO 290
      IF(LETTER==LETE  ) GO TO 200
      IF(LETTER==LETD  ) GO TO 200
      IF(LETTER==LCONT ) GO TO 208
      IF(LETTER==COMMA ) GO TO 210
      IF(LETTER==LPAREN) GO TO 210
      IF(LETTER==RPAREN) GO TO 210
      IF(LETTER==DECPT ) GO TO 220
      IF(LETTER==SLASH ) GO TO 230
      IF(LETTER==BLANK ) GO TO 250
      IF(LETTER==TAB  ) GO TO 250
      IF(LETTER==EQUALS) GO TO 250
      DO 140 J=1,10
      IF(LETTER==NUM(J)) GO TO 150
  140 CONTINUE
      LAST=NCOL
      IF(LETTER/=LETS) GO TO 290
      LAST=LAST+1
      IF(CARD(LAST:LAST)/=LETQ) GO TO 290
      LAST=LAST+1
      IF(CARD(LAST:LAST)/=LETR) GO TO 290
      LAST=LAST+1
      IF(CARD(LAST:LAST)/=LETT) GO TO 290
      LAST=LAST+1
      IF(CARD(LAST:LAST)/=LPAREN) GO TO 290
      NCOL=LAST+1
      ISQRT=1
      GO TO 130
!
  150 IF (LFIRST==0) LFIRST=NCOL
      IF(NCOL==MAXCOL) NEWCRD=1
  160 LEND=1
  170 NCOL=NCOL+1
      GO TO 130
!
  180 IF (LEND>0) GO TO 350
      GO TO 170
!
  190 LSIGN=-1
      GO TO 170
!
  200 IF (LFIRST/=0) GO TO 260
      IF (LEND==0) FIXPNT=1.0D0
      IF (CARD(NCOL+1:NCOL+1)==BLANK) NCOL=NCOL+1
      LSIGN=1
      KINDST=3
      GO TO 160
!
  208 NEWCRD=1
      IF(LFIRST/=0) GO TO 260
      IF(LEND>0) GO TO 350
      if(LEN_TRIM(BASIS_FILE)>0)then      ! read from file.bas
       CALL RDCARD(50,'   *    ',KEOF)  
      else
       CALL RDCARD(5,'   *    ',KEOF)
      endif
      NEWCRD=0
      GO TO 130
!
  210 ICOMMA=1
      IF(LFIRST/=0) GO TO 260
      GO TO 350
!
  220 CONTINUE
      IF(NDECPT>1) GO TO 290
      NDECPT=NDECPT+1
      IF (LFIRST/=0) GO TO 260
      ISIGN=LSIGN
      KINDST=2
      GO TO 160
!
  230 IF(LFIRST/=0) GO TO 260
      TNUM=(FIXPNT+FRACT*(10.0D0**IEXPFR))*(10.0D0**(ISIGNE*IEXP))
      IF(ISIGN<0) TNUM=-TNUM
      NCOL=NCOL+1
      GO TO 110
!
  240 IF (LEND==0) NCOL=LSTCOL+1
      LEND=1
  250 IF (LFIRST==0) GO TO 180
!                  DECODE DIGIT STRING AND STORE IN X
  260 LAST=NCOL-1
      X=0.0D0
      FAC=1.0D0
  270 CONTINUE
      DO 280 J=1,10
         IF (CARD(LAST:LAST)==NUM(J)) GO TO 300
  280 CONTINUE
!
  290 CONTINUE
      IF(LSTERR/=30 .AND. LUOUT>0) THEN
            WRITE(LUOUT,900) VARABL,LAST
            WRITE(LUOUT,910) CARD,(I,I=1,8)
      END IF
      IF(LSTERR/=30 .AND. LUERR>0 .AND. LUOUT<0)THEN
            WRITE(LUERR,900) VARABL,LAST
            WRITE(LUERR,910) CARD,(I,I=1,8)
      END IF
      KERR = 1
      RETURN
!
  300 X=X+FAC*(J-1)
      IF (LAST==LFIRST) GO TO 310
      LAST=LAST-1
      FAC = 10.0D0 * FAC
      GO TO 270
!
  310 CONTINUE
      IF(KINDST==3) GO TO 340
         IF(KINDST==2) GO TO 330
            FIXPNT=X
            ISIGN=LSIGN
            GO TO 120
!
  330    FRACT=X
         IEXPFR=LFIRST-NCOL
         GO TO 120
!
  340 IEXP=INT(X)
      ISIGNE=LSIGN
      GO TO 120
!
  350 CONTINUE
      VALUE=(FIXPNT+FRACT*(10.0D0**IEXPFR))*(10.0D0**(ISIGNE*IEXP))
      IF(TNUM/=0.0D0) VALUE=TNUM/VALUE
      IF(ISQRT==1) VALUE=SQRT(VALUE)
      IF(ISIGN<0) VALUE=-VALUE
      NCOL=NCOL+ICOMMA
      IF(NEWCRD==1)THEN
       if(LEN_TRIM(BASIS_FILE)>0)then      ! read from file.bas
        CALL RDCARD(50,'   *    ',KEOF)  
       else
        CALL RDCARD(5,'   *    ',KEOF)
       endif
      ENDIF
      RETURN
!
  900 FORMAT(1X,'**** ERROR READING VARIABLE ',A8,' CHECK COLUMN',I3)
  910 FORMAT(1X,A80/1X,8('....V....',I1))
      END

! GSTRNG
      SUBROUTINE GSTRNG(STRING,LENGTH)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL FIRST
      CHARACTER*(*) STRING
      CHARACTER*1 BLANK,TAB,QUOTE,EQUAL,COMMA,LPAREN
      CHARACTER*1 LCONT,LEOD,LEOC
      CHARACTER*80 CARD
      COMMON /FREFM1/ NCOL,LSTCOL,MAXCOL,KEOF,KERR,LUOUT,LUERR,KOLSV
      COMMON /FREFM2/ CARD,LCONT,LEOD,LEOC
      CHARACTER(80) :: BASIS_FILE
      COMMON/BASIS_FILE/BASIS_FILE
      DATA FIRST/.TRUE./
      SAVE FIRST,BLANK,QUOTE,EQUAL,COMMA,LPAREN,TAB
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      KOLSV = NCOL
      IF(FIRST) THEN
       BLANK = ' '
       QUOTE = ''''
       EQUAL = '='
       COMMA = ','
       LPAREN = '('
       TAB= CHAR(9)
       FIRST=.FALSE.
      END IF
!
  100 CONTINUE
      NOCHAR=MAX(IABS(LENGTH),4)
      NCC=0
      NBC=NOCHAR
      IF(NCOL>LSTCOL) GO TO 320
!
      LSTART=NCOL-1
      LSTOP=LSTART+NOCHAR
!
      DO 110 N=NCOL,LSTCOL
         IF(CARD(N:N)==TAB) GO TO 110
         IF(CARD(N:N)/=BLANK) GO TO 120
  110 CONTINUE
      NCOL=LSTCOL+1
      GO TO 300
!
  120 CONTINUE
      IF(CARD(N:N)/=LCONT) GO TO 130
!
      IF(STRING(1:4)=='   *') RETURN
      if(LEN_TRIM(BASIS_FILE)>0)then      ! read from file.bas
       CALL RDCARD(50,'   *    ',KEOF)  
      else
       CALL RDCARD(5,'   *    ',KEOF)
      endif
      GO TO 100
!
  130 CONTINUE
      NCOL=LSTOP+1
      IF(LENGTH<0) GO TO 200
      IF(LENGTH>0) GO TO 300
!
      NCOL=N
      IF(CARD(N:N)/=QUOTE) GO TO 300
!
      LSTART=NCOL-1
  140 CONTINUE
      NCOL=NCOL+1
!
      DO 150 N=NCOL,LSTCOL
         CARD(N-1:N-1)=CARD(N:N)
  150 CONTINUE
      CARD(LSTCOL:LSTCOL)=BLANK
!
      DO 160 N=NCOL,LSTCOL
         IF(CARD(N:N)==QUOTE) GO TO 170
  160 CONTINUE
      N=LSTCOL
  170 CONTINUE
      NCOL=N+1
!
      IF(NCOL<LSTCOL .AND. CARD(NCOL:NCOL)==QUOTE) GO TO 140
      LSTOP=N-1
      LENGTH=LSTOP-LSTART
      NOCHAR=MAX(LENGTH,4)
      GO TO 300
!
  200 CONTINUE
      NCOL=N
      LSTART=NCOL-1
!
      DO 210 N=NCOL,LSTCOL
         IF(CARD(N:N)==BLANK) GO TO 220
         IF(CARD(N:N)==EQUAL) GO TO 220
         IF(CARD(N:N)==COMMA) GO TO 220
         IF(CARD(N:N)==LPAREN) GO TO 220
  210 CONTINUE
      N=LSTCOL+1
  220 CONTINUE
      NCOL=N+1
      LSTOP=MIN(LSTART+NOCHAR,N-1)
!
  300 CONTINUE
      NCC=MIN(LSTOP,LSTCOL) - LSTART
      NBC=NOCHAR - NCC
      IF(NCC<=0) GO TO 320
         DO 310 N=1,NCC
            STRING(N:N)=CARD(LSTART+N:LSTART+N)
  310    CONTINUE
  320 CONTINUE
      IF(NBC==0) GO TO 340
!
         DO 330 N=1,NBC
            STRING(NCC+N:NCC+N)=BLANK
  330    CONTINUE
  340 CONTINUE
!
      IF(LENGTH<0) CALL UPRCAS(STRING,IABS(LENGTH))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -      
      RETURN
      END

! IFIND
      FUNCTION IFIND(VARABL,IERR)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*8 VARABL
      COMMON /FREFM1/ NCOL,LSTCOL,MAXCOL,KEOF,KERR,LUOUT,LUERR,KOLSV
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      KOLSV = NCOL
      KERR = IERR
      CALL DECODN(VARABL,VALUE,ISIGN,FIXPNT,FRACT,IEXPFR,ISIGNE,IEXP)
!
      IFRACT=INT(FRACT)
      INTT=INT(FIXPNT)
      I1=10**IEXP
      I2=IEXPFR+ISIGNE*IEXP
      I3=10**IABS(I2)
      IF (ISIGNE<0) GO TO 130
!
      IF (I2<0) THEN
         I5=IFRACT/I3
      ELSE
         I5=IFRACT*I3
      END IF
!
      IFIND=ISIGN*(INTT*I1+I5)
      IERR = KERR
      RETURN
!
  130 IFIND=ISIGN*INTT/I1
      IERR = KERR
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -      
      RETURN
      END

! ECPPAR
      SUBROUTINE ECPPAR(KTYPE,NSHELL,ZAN,Cxyz)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/INTOPT/ISCHWZ,IECP,NECP
      COMMON/INFOA/NAT,ICH,MUL,NUM,NQMT,NE,NA,NB                                    
      !JFHLewYee: Changed NATOMS allowed dimension from 100 to 1000
      COMMON/ECP2/CLP(4004),ZLP(4004),NLP(4004),KFRST(1001,6),          &
                  KLAST(1001,6),LMAX(1001),LPSKIP(1001),IZCORE(1001)
      COMMON/ECPDIM/NCOEF1,NCOEF2,J1LEN,J2LEN,LLIM,NLIM,NTLIM,J4LEN      
      CHARACTER(80) ::  BASIS_FILE
      COMMON/BASIS_FILE/BASIS_FILE
!
      INTEGER,DIMENSION(NSHELL) :: KTYPE
      DOUBLE PRECISION,DIMENSION(NAT) :: ZAN
      DOUBLE PRECISION,DIMENSION(3,NAT) :: Cxyz
      CHARACTER(8) :: PPNAME,PPTYPE,PTYPE
      CHARACTER(8), DIMENSION(100) :: TYPELP
      CHARACTER(8) :: ANONE,GEN,BLANK,ENDWRD
      DATA ANONE /'NONE    '/, GEN   /'GEN     '/
      DATA BLANK /'        '/, ENDWRD/'$END    '/
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     NECP: Number of electrons removed
!     IZCORE: Number of electrons removed by each atom
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      NECP = 0
      DO I=1,NAT
       IZCORE(I)= 0
       LMAX(I)  = 0
       LPSKIP(I)= 1
      END DO
      IF(IECP==0)RETURN                           ! No ECP is used
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      WRITE(6,'(/1X,A14,/1X,14(1H-))')'ECP Potentials'
      PPNAME = BLANK
      PPTYPE = BLANK
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Find Group $ECP in the Input File ( IECP = 1 )
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(IECP==1) THEN
#ifdef MPI
       CONTINUE
#else
       if(LEN_TRIM(BASIS_FILE)>0)then      ! read from file.bas
        REWIND(50)
       else
        REWIND(5)
       endif
#endif
       if(LEN_TRIM(BASIS_FILE)>0)then      ! read from file.bas
        CALL FNDGRP(50,'$ECP    ',JEOF)  
       else
        CALL FNDGRP(5,'$ECP    ',JEOF)
       endif
!
       IF(JEOF==1)THEN
        WRITE(6,'(/36A)')'Stop: No ECP group found, IECP=1'
        CALL ABRT
       END IF
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IADD= 0
      PTYPE= ANONE
      DO ICNTR=1,NAT
       NUCZ = INT(ZAN(ICNTR))
       LPSKIP(ICNTR) = 0
       IEOF= 0
!- - - - - - - - - - - - - - - - - - - - - - - - -
!      IECP = 1
!- - - - - - - - - - - - - - - - - - - - - - - - -
       IF(IECP==1)THEN
        IEOF= 0
        IERR= 0
        if(LEN_TRIM(BASIS_FILE)>0)then      ! read from file.bas
         CALL RDCARD(50,'$ECP   1',IEOF)  
        else
         CALL RDCARD(5,'$ECP   1',IEOF)
        endif
        KSIZE=-8
        CALL GSTRNG(PPNAME,KSIZE)
        IF(ADJUSTL(PPNAME)==ADJUSTL(ENDWRD)) THEN
         WRITE(6,'(A32)')'Stop: An unexpected $END in $ECP'
         CALL ABRT
        END IF
        READ(UNIT=PPNAME,FMT='(A8)')TYPELP(ICNTR)
        KSIZE=-8
        CALL GSTRNG(PPTYPE,KSIZE)
        READ(UNIT=PPTYPE,FMT='(A8)') PTYPE
        IZCORE(ICNTR)= IFIND('IZCORE  ',IERR)
        IF(IERR/=0) CALL ABRT
        LMAX(ICNTR)= IFIND('LMAX    ',IERR)
        IF(IERR/=0) CALL ABRT
       END IF
!- - - - - - - - - - - - - - - - - - - - - - - - -
!      No ECP Type
!- - - - - - - - - - - - - - - - - - - - - - - - -
       IF(PTYPE==ANONE) THEN
        LPSKIP(ICNTR)=1
        IZCORE(ICNTR)=0
        LMAX(ICNTR)=0
        GO TO 1
       END IF
!- - - - - - - - - - - - - - - - - - - - - - - - -
!       Type
!- - - - - - - - - - - - - - - - - - - - - - - - -
!- - - - - - - - - - - - - - - - - - - - - - - - -
!      Check if PP has occured before
!- - - - - - - - - - - - - - - - - - - - - - - - -
       IF(ICNTR/=1)THEN
        NCNTR = ICNTR-1
        DO JCNTR=1,NCNTR
         IF(TYPELP(ICNTR)/=TYPELP(JCNTR)) GO TO 3
         LMN = 1
         LMX = LMAX(JCNTR)+1
         DO  L=LMN,LMX
          KFRST(ICNTR,L)= KFRST(JCNTR,L)
          KLAST(ICNTR,L)= KLAST(JCNTR,L)
         END DO
         LMAX(ICNTR)= LMAX(JCNTR)
         IZCORE(ICNTR)= IZCORE(JCNTR)
         GO TO 2
    3    CONTINUE
        END DO
       END IF
!- - - - - - - - - - - - - - - - - - - - - - - - -
!      GEN Type
!- - - - - - - - - - - - - - - - - - - - - - - - -
       IF(PTYPE==GEN) THEN
        LMN = 1
        LMX = LMAX(ICNTR)+1
        DO L=LMN,LMX
         IEOF= 0
         IERR= 0
         if(LEN_TRIM(BASIS_FILE)>0)then      ! read from file.bas
          CALL RDCARD(50,'$ECP   2',IEOF)  
         else
          CALL RDCARD(5,'$ECP   2',IEOF)
         endif
         NGPOT = IFIND('NGPOT   ',IERR)
         KF = IADD+1
         KL = IADD+NGPOT
         IF(KL>404) THEN
          WRITE(6,*)'Stop: KL > 404'
          CALL ABRT
         END IF
         KFRST(ICNTR,L)= KF
         KLAST(ICNTR,L)= KL
         DO K=KF,KL
          IEOF= 0
          IERR= 0
          if(LEN_TRIM(BASIS_FILE)>0)then      ! read from file.bas
           CALL RDCARD(50,'$ECP   3',IEOF)  
          else
           CALL RDCARD(5,'$ECP   3',IEOF)
          endif
          CLP(K)= RFIND('ECP COEF',IERR)
          IF(IERR/=0) CALL ABRT
          NLP(K)= IFIND('POLYNO N',IERR)
          IF(IERR/=0) CALL ABRT
          ZLP(K)= RFIND('ECP EXP.',IERR)
          IF(IERR/=0) CALL ABRT
         END DO
         IADD = KL
        END DO
       END IF
!- - - - - - - - - - - - - - - - - - - - - - - - -
    2  CONTINUE
!- - - - - - - - - - - - - - - - - - - - - - - - -
!      Number of Electrons and Nuclear Charge
!- - - - - - - - - - - - - - - - - - - - - - - - -
       NE = NE-IZCORE(ICNTR)
       NA = NA-IZCORE(ICNTR)/2
       NB = NB-IZCORE(ICNTR)/2
       ZAN(ICNTR) = ZAN(ICNTR)-IZCORE(ICNTR)
    1  CONTINUE
      END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Print the Pseudo Potentials
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      NUMECP = 0
      DO ICNTR=1,NAT
       IF(LPSKIP(ICNTR)==1) GOTO 4
       NUMECP = NUMECP+1
       NCNTR  = ICNTR-1
       DO JCNTR=1,NCNTR
        IF(TYPELP(ICNTR)/=TYPELP(JCNTR)) GO TO 5
         WRITE(6,10) TYPELP(ICNTR),ICNTR,JCNTR
        GOTO 4
    5   CONTINUE
       END DO

       WRITE(6,20)TYPELP(ICNTR),ICNTR,IZCORE(ICNTR),LMAX(ICNTR)
       LMN = 1
       LMX = LMAX(ICNTR)+1
       DO L=LMN,LMX
        LM1=L-2
        IF(L==LMN) LM1=LMX-1
         WRITE(6,30) LM1
        KF = KFRST(ICNTR,L)
        KL = KLAST(ICNTR,L)
        KK= 1
        DO K=KF,KL
          WRITE(6,40) KK,CLP(K),NLP(K),ZLP(K)
         KK = KK+1
        ENDDO
       END DO
    4  CONTINUE
      END DO

      IF(NUMECP==0) THEN
       WRITE(6,*)'Stop: This run has no atoms with ECP'
       CALL ABRT
      END IF

      NECP= 0
      DO I=1,NAT
       NECP = NECP+IZCORE(I)
      END DO
      WRITE(6,50)NECP,NECP,NE,NA,NB
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Adjusted Nuclear Energy
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Vnn = ENUC(NAT,ZAN,Cxyz)
      WRITE(6,'(/A26,F20.10/)')' Adjusted Nuclear Energy =',Vnn
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Check the Basis
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL BASCHK(MAXANG,KTYPE,NSHELL)
      IF(MAXANG>=5)THEN
       WRITE(6,*)'Stop: ECP integrals only for S,P,D,F,G'
       CALL ABRT
      END IF
      NLIM = MAXANG+1
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Check the Potential Type
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      LLIM = 0
      DO IAT=1,NAT
       LLIM = MAX0(LMAX(IAT),LLIM)
      END DO
      NLIM= MAX0(NLIM,LLIM)
      IF(LLIM>=5) THEN
       WRITE(6,*)'Stop: ECP integrals only for S,P,D,F,G'
       CALL ABRT
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Test NLIM
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      NTLIM=(NLIM*(NLIM+1)*(NLIM+2))/6
      NCOEF1 = 8520
      NCOEF2 = 3424
      IF(NLIM==5) THEN
       NCOEF1 = 71660
       NCOEF2 = 10555
      ELSE IF(NLIM==6) THEN
       NCOEF1 = 280000
       NCOEF2 = 28940
      ELSE IF(NLIM==7) THEN
       NCOEF1= 892584
       NCOEF2= 60382
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      J1LEN = (NTLIM*NTLIM+NTLIM)/2+1
      J2LEN = (LLIM-1)*(LLIM+1)*NTLIM+NTLIM+1
      J4LEN = J2LEN-1
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   10 FORMAT(' Parameters for "',A8,'" on Atom',I3,                     &
             ' are the same as atom',I3)                                 
   20 FORMAT(' Parameters for "',A8,'" on Atom',I3,                     &
             ' with ZCORE',I3,' and LMAX',I2,' are'/)                    
   30 FORMAT(' for L=',I2,6X,'COEFF',4X,'N',11X,'ZETA')                  
   40 FORMAT(I5,F15.5,I5,F15.5)                                          
   50 FORMAT(/' The ECP removes',I4,' core electrons and',I4,' protons'/&
             /' Number of electrons kept in the run =',I5               &
             /' Number of occupied orbitals (Alpha) =',I5               &
             /' Number of occupied orbitals (Beta)  =',I5)
!-----------------------------------------------------------------------
      END
