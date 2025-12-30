!======================================================================!
!                                                                      !
!             O R B I T A L       O P T I M I Z A T I O N              !
!                                                                      !
!             Variables: MOs {Cvq}  (NVAR = NSQ = NBF*NBF)             !
!                                                                      !
!   OrbOpt: Minimize the energy with respect to the MOs {Cvq} = COEF   !
!                                                                      !
!======================================================================!
!                                                                      !
!   Method:                                                            !
!                                                                      !
!   1) OrbOptFMIUGr: Iterative Diagonalization using a Generalized     !
!                    Fock Matrix (FMIUG). MP, JMU, JCC 30, 2078, 2009  !
!                                                                      !
!   2) OrbOptADAM: Orbital rotaions [Cn=C*exp(Y)] through a deep       !
!                  learning-inspired adaptative momentum technique     !
!                  (ADAM). JFHLY, JMC, MP; PRL 134, 206401, 2025       !
!                                                                      !
!======================================================================!

! OrbOpt
      SUBROUTINE OrbOpt(ITCALL,ITLIM,OVERLAP,AHCORE,IJKL,XIJKL,XIJKaux, &
                        QD,COEF,RO,CJ12,CK12,ELAG,FMIUG0,DIPN,ILOOP,    &
                        IORBOPT,OCCTIME,IPRINTOPT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                        
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPFILE_NIJKL/NINTMX,NIJKL,NINTCR,NSTORE
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5      
      LOGICAL ERIACTIVATED
      COMMON/ERIACT/ERIACTIVATED,NIJKaux,NINTCRaux,NSTOREaux,IAUXDIM
      COMMON/SKIPRINT/ISKIP
!      
      INTEGER :: ITCALL,ITLIM,ILOOP,IPRINTOPT,IORBOPT
      INTEGER(8),DIMENSION(NSTORE)::IJKL
      DOUBLE PRECISION,DIMENSION(NSTORE)::XIJKL
      DOUBLE PRECISION,DIMENSION(NSTOREaux)::XIJKaux 
      DOUBLE PRECISION,DIMENSION(NBF)::FMIUG0      
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::OVERLAP,AHCORE,COEF,ELAG
      DOUBLE PRECISION,DIMENSION(NBF,NBF,NBF)::QD
      DOUBLE PRECISION,DIMENSION(NBF5)::RO
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5)::CJ12,CK12
      DOUBLE PRECISION,DIMENSION(3)::DIPN
!-----------------------------------------------------------------------
      IF(IORBOPT==1)THEN
       CALL OrbOptFMIUGr(ITCALL,ITLIM,OVERLAP,AHCORE,IJKL,XIJKL,        &
                         XIJKaux,QD,COEF,RO,CJ12,CK12,ELAG,FMIUG0,      &
                         DIPN,ILOOP,IPRINTOPT)
      ELSE IF(IORBOPT==2)THEN
       ISKIP=0
       CALL OrbOptADAM(ITCALL,OVERLAP,AHCORE,IJKL,XIJKL,XIJKaux,QD,COEF,&
                       RO,CJ12,CK12,ELAG,DIPN,ILOOP,OCCTIME,IPRINTOPT)
      ELSE IF(3<=IORBOPT .and. IORBOPT<=6)THEN
!OO36  CALL OrbOpt36(ITCALL,OVERLAP,AHCORE,IJKL,XIJKL,XIJKaux,QD,COEF,  &
!OO36                RO,CJ12,CK12,ELAG,DIPN,ILOOP,IORBOPT,OCCTIME,      &
!OO36                IPRINTOPT)
      END IF
!-----------------------------------------------------------------------
      RETURN
      END

! OrbOptFMIUGr      
      SUBROUTINE OrbOptFMIUGr(ITCALL,ITLIM,OVERLAP,AHCORE,IJKL,XIJKL,   &
                              XIJKaux,QD,COEF,RO,CJ12,CK12,ELAG,FMIUG0, &
                              DIPN,ILOOP,IPRINTOPT)
!----------------------------------------------------------------------!
!                                                                      !
!      This subroutine minimizes w.r.t. MOs using a Generalized Fock   !
!      Matrix for Spin-compensated Systems (Restricted Shells)         !
!                                                                      !
!      1) MSpin = 0 ( Restricted Closed: rc )                          !
!         Singlet States (S=0,Ms=0) and Multiplet States (S>0,Ms=0)    !
!                                                                      !
!      2) MSpin > 0 ( Restricted Open: ro )                            !
!         High-Spin States (Ms=S)                                      !
!                                                                      !
!----------------------------------------------------------------------!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      LOGICAL CONVGDELAG,EFIELDL,SMCD,RESTART,SCALING
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/CONVERGENCE/DUMEL,PCONV,CONVGDELAG
      COMMON/INP_EFIELDL/EFX,EFY,EFZ,EFIELDL
      COMMON/ERITYPE/IERITYP,IRITYP,IGEN,ISTAR,MIXSTATE,SMCD      
      COMMON/INPNOF_RSTRT/RESTART,INPUTGAMMA,INPUTC,INPUTFMIUG,INPUTCXYZ
      COMMON/INPSCALING/SCALING,NZEROS,NZEROSm,NZEROSr,ITZITER      
      COMMON/INPNOF_THRESH/THRESHL,THRESHE,THRESHEC,THRESHEN
      COMMON/INPNOF_COEFOPT/MAXLOOP
      COMMON/INPNOF_PRINT/NPRINT,IWRITEC,IMULPOP,IAIMPAC,IFCHK
      COMMON/INPNOF_PNOF/IPNOF,NTWOPAR
      COMMON/INPNOF_MOD/Imod
      LOGICAL ERIACTIVATED,HighSpin
      COMMON/ERIACT/ERIACTIVATED,NIJKaux,NINTCRaux,NSTOREaux,IAUXDIM      
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
      LOGICAL DIIS,PERDIIS,DAMPING,EXTRAP            
      COMMON/INPNOF_DIIS/DIIS,PERDIIS,NDIIS,NTHDIIS,THDIIS
      COMMON/INPNOF_DAMPEXTRAP/DAMPING,EXTRAP      
      COMMON/ACONV/RRSHFT,EXTTOL,DMPTOL      
      COMMON/INPFILE_NIJKL/NINTMX,NIJKL,NINTCR,NSTORE
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      COMMON/RUNTYPE/IRUNTYP
      COMMON/CONVERGESUM/SUMDIF,SUMDIF_OLD
      COMMON/EHFEN/EHF,EN
      COMMON/ENERGIAS/EELEC,EELEC_OLD,DIF_EELEC,EELEC_MIN
!
      INTEGER :: NSOC,NDNS,MSpin
      INTEGER :: ITCALL,ITLIM,ILOOP,IPRINTOPT
      INTEGER(8),DIMENSION(NSTORE)::IJKL
      DOUBLE PRECISION,DIMENSION(NSTORE)::XIJKL
      DOUBLE PRECISION,DIMENSION(NSTOREaux)::XIJKaux
      DOUBLE PRECISION,DIMENSION(NBF)::FMIUG0
      DOUBLE PRECISION,DIMENSION(NBF5)::RO
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::OVERLAP,AHCORE,COEF,ELAG
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5)::CJ12,CK12
      DOUBLE PRECISION,DIMENSION(3)::DIPN
      DOUBLE PRECISION,DIMENSION(NBF,NBF,NBF)::QD
!      
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::EVA,TEMP,CFM
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::DMP,DM,WRK1,WRK2,WRK3
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::FAO,FAO0,FAO1,FAO2
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::DEN
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::FMIUG,W      
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::COEFNEW,BFM,G
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:)::FK      
!-----------------------------------------------------------------------
!     First Call to OrbOptFMIUGr
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(ITCALL==1)THEN
       IF(IRUNTYP==3.or.IRUNTYP==4.or.IRUNTYP==5)NZEROS=NZEROSr
       IF(IPRINTOPT==1)WRITE(6,1)
       EELEC_OLD  = EELEC
       SUMDIF_OLD = 0.0d0
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Density Matrices (DEN)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(DEN(NBF,NBF))
      CALL DENMATr(DEN,COEF,RO,NBF,1,NBF5)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Compute Gradient & Lagrangian Multipliers (ELAG)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(G(NBF,NBF5))
      IF( (IPNOF==5.or.IPNOF==7.or.(IPNOF==8.and.Imod==0)) .and.        &
              MSpin==0 .and. IERITYP==1 )THEN
!AO          
       OVERLAP(1,1) = OVERLAP(1,1) ! avoiding warning
!AO    CALL ELAGaor(OVERLAP,AHCORE,IJKL,XIJKL,COEF,RO,DEN,ELAG)
       CALL ELAGr(AHCORE,IJKL,XIJKL,COEF,RO,DEN,ELAG,G,IPNOF)
      ELSE
       CALL ELAG1r(AHCORE,IJKL,XIJKL,XIJKaux,QD,COEF,RO,CJ12,CK12,ELAG,G)
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Calculate Electronic Energy (EELEC)       
!     Note: This Energy is equal to EELEC from MOLOCUPrc
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -      
      EELEC = TRACE(DEN,AHCORE,NBF)
      DO iq=1,NBF5
       EELEC = EELEC + ELAG(iq,iq)
      END DO
!     Include Nuclear Dipoles if electric field =/ 0
      IF(EFIELDL)EELEC = EELEC - EFX*DIPN(1)-EFY*DIPN(2)-EFZ*DIPN(3)
      Etot = EELEC + EN
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -      
!     Check for the Symmetry of Lagrangian & Energy
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL PCONVE(ELAG,DUMEL,MAXI,MAXJ,SUMDIF)         ! DUMEL
      PCONV=ABS(DIF_EELEC)                             ! PCONV
      IF(ITCALL==1)THEN
       PCONV=THRESHE+1.0d0                             ! avoiding stop        
       IF(IPRINTOPT==1)THEN
        DIF_EELEC = EELEC-EELEC_OLD
        WRITE(6,2)ITCALL,EELEC,Etot,DIF_EELEC,DUMEL
       ENDIF      
      ENDIF
      IF( DUMEL<THRESHL .and. PCONV<THRESHE )THEN
       IF(ITCALL>1.and.IPRINTOPT==1)THEN
        WRITE(6,2)ITCALL,EELEC,Etot,DIF_EELEC,DUMEL 
       END IF
       CONVGDELAG = .TRUE.
       RETURN
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Check ITCALL and SUMDIF for changing the scaling factor
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(SCALING.and.ITCALL>2.and.ITCALL>ITLIM.and.                     &
         SUMDIF>SUMDIF_OLD)THEN
       NZEROS=NZEROS+1
       ITLIM=ITCALL+ITZITER
!      NZEROS>NZEROSm -> Restart from here with NZEROSr
       IF(NZEROS>NZEROSm)THEN
        NZEROS=NZEROSr
       ENDIF
      ENDIF
      SUMDIF_OLD=SUMDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -      
!     Initialize Variables
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(FMIUG(NBF,NBF),W(NBF,NBF),EVA(NBF))
      ALLOCATE(TEMP(NBF),COEFNEW(NBF,NBF))
      IF(DIIS)THEN
       ALLOCATE(CFM(MAXLOOP+1),BFM(MAXLOOP+1,MAXLOOP+1))
       ALLOCATE(FK(MAXLOOP,NBF,NBF))
      ENDIF
      IF(DAMPING.or.EXTRAP)THEN
       ALLOCATE(DMP(NBFT),DM(NBFT),WRK1(NSQ),WRK2(NSQ),WRK3(NSQ))
       ALLOCATE(FAO(NBFT),FAO0(NBFT),FAO1(NBFT),FAO2(NBFT))
      ENDIF
!      
      ILOOP = 0
      DELE = 1.0d20
      IDIIS = 0
!0    IF(ITCALL==1.and.INPUTFMIUG==0)THEN
!0     MAXLP=1
!0    ELSE
!0     MAXLP=MAXLOOP
!0    ENDIF
!
      IF(INPUTFMIUG==0)THEN
       if(ITCALL==1)then
        MAXLP=1
       else
        MAXLP=3*ITCALL                                       
        IF(MAXLP>MAXLOOP)MAXLP=MAXLOOP                       
       endif
      ELSE
       MAXLP=MAXLOOP
      ENDIF
!      
      DAMP   = 0.0d0                                                       
      DAMP0  = 0.0d0                                                      
      IF(DAMPING)DAMP = 1.0d0
      DIFFP = 0.0d0                                                       
      DIFF  = 0.0d0 
      ITERV = 0                                                                      
      RRSHFT = 0.0d0
      EXTTOL = 1.0D-03                                                  
      DMPTOL = 1.0D-04                                                  
      ICALP  = 0
      ICBET  = 0
!-----------------------------------------------------------------------
!                       START SCF-ITERATION CYCLE
!-----------------------------------------------------------------------
      DO WHILE( ILOOP<MAXLP .and. DABS(DELE)>THRESHEC )
       ILOOP = ILOOP+1
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Generalized Fock Matrix (FMIUG)
!
!      Convergent technique:
!
!      1) SCALING:   Decrease FMIUG using a scaling factor.
!                    The scaling factor varies until the number of
!                    ZEROS (.000##) is equal for all elements Fij
!      2) DIIS       Direct Inversion in the Iterative Subspace 
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       CALL FORMFMIUG(FMIUG,ELAG,FMIUG0,ITCALL)
       IF(DUMEL<THDIIS)THEN
        IF(DAMPING.or.EXTRAP)CALL SQUARETRIAN(FMIUG,FAO,NBF,NBFT)
        IF(DIIS)THEN
         CALL FFMIUG_DIIS(NBF,FMIUG,CFM,BFM,FK,IDIIS)
         if(IDIIS>0)then
          if(DAMPING)DAMP = 1.0d0                                                   
          RRSHFT = 0.0d0                                                 
         endif                                                         
        ENDIF
       ENDIF        
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -   
!      Damping and Extrapolation
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       IF(ILOOP==2)DEAVG = ABS(DELE)
       IF(IDIIS==0.and.DUMEL<THDIIS)THEN
        IF(DAMPING.and.ILOOP>2)THEN
         DEAVG = (ABS(DELE)+ABS(DELE0)+0.2d0*DEAVG)/2.2D0          
         DAMP0 = DAMP
         CALL DAMPD(DELE,DELE0,DEAVG,DAMP,1.0D-05,DIFF,DIFFP,1.0D-02)
        END IF
        IF(DAMP<0.0d0)DAMP = 0.0d0 
        IF(DAMPING.or.EXTRAP)THEN        
         CALL EXTRAPOL(DELE,DAMP,DAMP0,FAO,WRK1,WRK2,WRK3,FAO0,FAO1,    &
                       FAO2,NBF,NBFT,ITERV,1,1,ILOOP,ICALP,ICBET)
         CALL TRIANSQUARE(FMIUG,FAO,NBF,NBFT)                       
        ENDIF                      
       END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      DIAGONALIZE SQUARE MATRIX (FMIUG) FOR REAL SYMMETRIC CASE
!      W - EIGENVECTORS, EVA - EIGENVALUES IN ALGEBRAIC DESCENDING ORDER
!      HOUSEHOLDER METHOD
!      NOTE: ONLY LOWER TRIANGLE IS USED + THIS IS DESTROYED !!!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       CALL DIAG(NBF,FMIUG,W,EVA,TEMP)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Move EVA -> FMIUG0
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       DO I=1,NBF
        FMIUG0(I) = EVA(I)
       ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      New Coefficients (COEFNEW=COEF*W), Move COEFNEW -> COEF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       CALL COEFW1(NBF,NBF,COEFNEW,COEF,W)
       COEF = COEFNEW
       CALL DENMATr(DEN,COEF,RO,NBF,1,NBF5) ! Density Matrix
!       
       IF( (DAMPING.or.EXTRAP) .and. DUMEL<THDIIS )THEN
        CALL DCOPY(NBFT,DM,1,DMP,1) ! DM -> DMP
        DIFFP = DIFF
        CALL SQUARETRIAN(DEN,DM,NBF,NBFT)
        CALL DDIFF(DMP,DM,NBFT,DIFF)
        ICALP = ICALP+1                                                   
        ICBET = ICBET+1
       ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Compute Gradient & Lagrangian Multipliers (ELAG)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       IF( (IPNOF==5.or.IPNOF==7.or.(IPNOF==8.and.Imod==0)) .and.       &
              MSpin==0 .and. IERITYP==1 )THEN
!AO     CALL ELAGaor(OVERLAP,AHCORE,IJKL,XIJKL,COEF,RO,DEN,ELAG)
        CALL ELAGr(AHCORE,IJKL,XIJKL,COEF,RO,DEN,ELAG,G,IPNOF)
       ELSE
        CALL ELAG1r(AHCORE,IJKL,XIJKL,XIJKaux,QD,COEF,RO,CJ12,CK12,ELAG,G)
       ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -              
!      Calculate Electronic Energy (EELEC)       
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
       EELEC0 = EELEC
       Etot0 = Etot
       DELE0 = DELE       
       EELEC = TRACE(DEN,AHCORE,NBF)
       DO iq=1,NBF5
        EELEC = EELEC + ELAG(iq,iq)
       END DO
!      Include Nuclear Dipoles if electric field =/ 0
       IF(EFIELDL)EELEC = EELEC - EFX*DIPN(1)-EFY*DIPN(2)-EFZ*DIPN(3)
       DELE = EELEC - EELEC0
       Etot = EELEC + EN       
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -       
!      Intermediate Output of the internal iteration (Nprint=2)                                    
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       IF(NPRINT==2.AND.IPRINTOPT==1)WRITE(6,3)ILOOP,EELEC,Etot,DELE
       CALL PCONVE(ELAG,DUMEL,MAXI,MAXJ,SUMDIF)
!-----------------------------------------------------------------------
!                       LOOP-END OF SCF-ITERATION
!-----------------------------------------------------------------------
      ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Itermediate Output of the external iteration
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DIF_EELEC = EELEC-EELEC_OLD
      EELEC_OLD = EELEC
      IF(ITCALL>1.and.IPRINTOPT==1)THEN
       if(DABS(DELE)>THRESHEC)then
        WRITE(6,2)ITCALL,EELEC,EELEC+EN,DIF_EELEC,DUMEL       
       else
        WRITE(6,21)ITCALL,EELEC,EELEC+EN,DIF_EELEC,DUMEL
       endif
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -      
!     FORMAT STATEMENTS
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -      
    1 FORMAT(//2X,'RESULTS OF OCCUPATION-COEFFICIENT S.C.F. PROCEDURE'  &
            ,/1X,'====================================================',&
        //2X,'Iter',5X,'Electronic Energy',6X,'Total Energy',           &
          3X,'Energy Convergency',4X,'Max Mul-Lag Diff',/)
    2 FORMAT(I5,'.',1X,F20.10,1X,F19.10,2X,F15.10,8X,F11.6)
   21 FORMAT(I5,'.*',F20.10,1X,F19.10,2X,F15.10,8X,F11.6)
    3 FORMAT(2X,I3,'.',3X,F17.8,4X,F15.8,6X,F11.6)
!-----------------------------------------------------------------------
      DEALLOCATE(COEFNEW,DEN,FMIUG,W,EVA,TEMP,G)
      IF(DIIS)DEALLOCATE(CFM,BFM,FK)
      IF(DAMPING.or.EXTRAP)THEN
       DEALLOCATE(DMP,DM,WRK1,WRK2,WRK3,FAO,FAO0,FAO1,FAO2)
      ENDIF
      RETURN
      END

! OrbOptADAM
      SUBROUTINE OrbOptADAM(ITCALL,OVERLAP,AHCORE,IJKL,XIJKL,XIJKaux,   &
                            QD,COEF,RO,CJ12,CK12,ELAG,DIPN,ILOOP,       &
                            OCCTIME,IPRINTOPT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL CONVGDELAG,EFIELDL,SMCD,ERIACTIVATED,AUTOLR
      COMMON/INPADAM/LR,FACT,BETA1,BETA2,AUTOLR
      COMMON/CONVERGENCE/DUMEL,PCONV,CONVGDELAG
      COMMON/INP_EFIELDL/EFX,EFY,EFZ,EFIELDL
      COMMON/ERITYPE/IERITYP,IRITYP,IGEN,ISTAR,MIXSTATE,SMCD
      COMMON/ERIACT/ERIACTIVATED,NIJKaux,NINTCRaux,NSTOREaux,IAUXDIM
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPFILE_NIJKL/NINTMX,NIJKL,NINTCR,NSTORE
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      COMMON/EHFEN/EHF,EN
      COMMON/INPNOF_THRESH/THRESHL,THRESHE,THRESHEC,THRESHEN
      COMMON/INPNOF_PRINT/NPRINT,IWRITEC,IMULPOP,IAIMPAC,IFCHK
      COMMON/INPNOF_ORBSPACE0/NO1,NDOC,NCO,NCWO,NVIR,NAC,NO0
      COMMON/ENERGIAS/EELEC,EELEC_OLD,DIF_EELEC,EELEC_MIN
      COMMON/INPNOF_COEFOPT/MAXLOOP
      COMMON/INPNOF_PNOF/IPNOF,NTWOPAR
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
      COMMON/INPNOF_MOD/Imod
      COMMON/INPNOF_HFID/KOOPMANS
      COMMON/INPNOF_RHF/CONVRHFDM,IRHF,IRHFTYP,NCONVRHF,MAXITRHF
      COMMON/SKIPRINT/ISKIP
!
      INTEGER :: NO1,NDOC,NCO,NCWO,NVIR,NAC,NO0
      INTEGER :: NSOC,NDNS,MSpin
      LOGICAL :: HighSpin
      INTEGER :: NIT
      INTEGER :: ITCALL,ILOOP,IPRINTOPT
      INTEGER(8),DIMENSION(NSTORE)::IJKL
      DOUBLE PRECISION,DIMENSION(NSTORE)::XIJKL
      DOUBLE PRECISION,DIMENSION(NSTOREaux)::XIJKaux
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::OVERLAP,AHCORE,COEF,ELAG
      DOUBLE PRECISION,DIMENSION(NBF,NBF,NBF)::QD
      DOUBLE PRECISION,DIMENSION(NBF5)::RO
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5)::CJ12,CK12
      DOUBLE PRECISION,DIMENSION(3)::DIPN
!
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::Y,GRAD
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::M,V,MHAT,VHAT,VHATMAX
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::BEST_COEF
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::G,DEN
      LOGICAL :: IMPROVED
      DOUBLE PRECISION LR, FACT, BETA1, BETA2
      INTEGER :: CR, CM, TIMESTART, TIMEFINISH
      DOUBLE PRECISION :: RATE
!-----------------------------------------------------------------------
!     Initialization for system_clock
!-----------------------------------------------------------------------
      CALL SYSTEM_CLOCK(COUNT_RATE=CR)
      CALL SYSTEM_CLOCK(COUNT_MAX=CM)
      RATE = REAL(CR)
      CALL SYSTEM_CLOCK(TIMESTART)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     First Call to OrbOptADAM
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(ITCALL==1)THEN
       if(IPRINTOPT==1)then
        if(ISKIP==0)WRITE(6,10)
        if(ISKIP==1)WRITE(6,1)
       end if
       EELEC_OLD = EELEC
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ILOOP = ILOOP
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Density Matrices (DEN)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(G(NBF,NBF5))
      ALLOCATE(DEN(NBF,NBF))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Compute Gradient & Lagrangian Multipliers (ELAG)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL DENMATr(DEN,COEF,RO,NBF,1,NBF5) ! Density Matrix
      IF( (IPNOF==5.or.IPNOF==7.or.(IPNOF==8.and.Imod==0)) .and.        &
              MSpin==0 .and. IERITYP==1 )THEN
!AO    CALL ELAGaor(OVERLAP,AHCORE,IJKL,XIJKL,COEF,RO,DEN,ELAG)
       CALL ELAGr(AHCORE,IJKL,XIJKL,COEF,RO,DEN,ELAG,G,IPNOF)
      ELSE
       CALL ELAG1r(AHCORE,IJKL,XIJKL,XIJKaux,QD,COEF,RO,CJ12,CK12,ELAG,G)
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Calculate Electronic Energy (EELEC)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      EELEC = TRACE(DEN,AHCORE,NBF)
      DO iq=1,NBF5
       EELEC = EELEC + ELAG(iq,iq)
      END DO
!     Include Nuclear Dipoles if electric field =/ 0
      IF(EFIELDL)EELEC = EELEC - EFX*DIPN(1)-EFY*DIPN(2)-EFZ*DIPN(3)
      Etot = EELEC + EN
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Check for the Symmetry of Lagrangian & Energy
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL PCONVE(ELAG,DUMEL,MAXI,MAXJ,SUMDIF)        ! DUMEL
      PCONV=ABS(DIF_EELEC)                            ! PCONV
      IF( DUMEL<THRESHL .and. PCONV<THRESHE)THEN
       CALL SYSTEM_CLOCK(TIMEFINISH)
       ORBTIME = (TIMEFINISH - TIMESTART)/RATE
       IF(IPRINTOPT==1)WRITE(6,2)ITCALL,OCCTIME,ORBTIME,ILOOP,EELEC,    &
                                 Etot,DIF_EELEC,DUMEL
       CONVGDELAG = .TRUE.
       AUTOLR=.TRUE.
       LR=0.01D0
       MAXLOOP=10
       RETURN
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Initialize variables
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      NV = NBF*(NBF-1)/2
      ALLOCATE(M(NV), V(NV))
      ALLOCATE(Y(NV), GRAD(NV), MHAT(NV), VHAT(NV), VHATMAX(NV))
      Y = 0.0D0
      GRAD = 0.0D0
      M = 0.0D0
      V = 0.0D0
      MHAT = 0.0D0
      VHAT = 0.0D0
      VHATMAX = 0.0D0
      IMPROVED = .FALSE.
      ALLOCATE(BEST_COEF(NBF,NBF))
      BEST_E = EELEC
      BEST_COEF = COEF
!-----------------------------------------------------------------------
!                       START SCF-ITERATION CYCLE
!-----------------------------------------------------------------------
      DELE = 1.0D20
      ILOOP = 0
!     DO WHILE( ILOOP<MAXLOOP .and. DABS(DELE)>THRESHE )
      DO WHILE( ILOOP<MAXLOOP.and.DUMEL>THRESHL)
       ILOOP = ILOOP+1
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Compute Gradient from Lagrangian Multipliers (ELAG)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       GRAD = 0.0D0
       K = 1
       DO I=1,NBF-1
        DO J=I+1,NBF
         GRAD(K) = 2.0d0 * ( ELAG(I,J) - ELAG(J,I) )
         K = K + 1
        END DO
       END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Rotate Oribtals with adaptative momentum step
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       DO J=1,NV
         M(J) = BETA1 * M(J) + (1.0D0-BETA1) * GRAD(J)
         V(J) = BETA2 * V(J) + (1.0D0-BETA2) * GRAD(J)**2
         MHAT(J) = M(J) / (1.0D0-BETA1**ILOOP)
         VHAT(J) = V(J) / (1.0D0-BETA2**ILOOP)
         VHATMAX(J) = MAX(VHATMAX(J), VHAT(J))
         Y(J) = -LR*MHAT(J)/(SQRT(VHATMAX(J)) + 10D-16)
       END DO
       CALL RotOrb(NV,Y,COEF)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Compute Gradient & Lagrangian Multipliers (ELAG)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       CALL DENMATr(DEN,COEF,RO,NBF,1,NBF5) ! Density Matrix
       IF( (IPNOF==5.or.IPNOF==7.or.(IPNOF==8.and.Imod==0)) .and.       &
             MSpin==0 .and. IERITYP==1 )THEN
!AO     CALL ELAGaor(OVERLAP,AHCORE,IJKL,XIJKL,COEF,RO,DEN,ELAG)
        CALL ELAGr(AHCORE,IJKL,XIJKL,COEF,RO,DEN,ELAG,G,IPNOF)
       ELSE
        CALL ELAG1r(AHCORE,IJKL,XIJKL,XIJKaux,QD,COEF,RO,CJ12,CK12,ELAG,G)
       ENDIF
       CALL PCONVE(ELAG,DUMEL,MAXI,MAXJ,SUMDIF)        ! DUMEL
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Calculate Electronic Energy (EELEC)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       EELEC0 = EELEC
       EELEC = TRACE(DEN,AHCORE,NBF)
       DO iq=1,NBF5
        EELEC = EELEC + ELAG(iq,iq)
       END DO
!      Include Nuclear Dipoles if electric field =/ 0
       IF(EFIELDL)EELEC = EELEC - EFX*DIPN(1)-EFY*DIPN(2)-EFZ*DIPN(3)
       DELE = EELEC - EELEC0
       Etot = EELEC + EN
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Store best energy and its orbitals
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       IF(EELEC < BEST_E) THEN
        IMPROVED = .TRUE.
        BEST_E = EELEC
        BEST_COEF = COEF
       END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Intermediate Output of the internal iteration (Nprint=2)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       IF(NPRINT==2.AND.IPRINTOPT==1)WRITE(6,3)ILOOP,EELEC,Etot,DELE
!-----------------------------------------------------------------------
!                       LOOP-END OF SCF-ITERATION
!-----------------------------------------------------------------------
      END DO
      DEALLOCATE(GRAD, M, V, MHAT, VHAT, VHATMAX, Y)
      IF(AUTOLR) THEN
       IF(.NOT. IMPROVED) THEN
        EELEC = BEST_E
        COEF = BEST_COEF
       ELSE
        AUTOLR = .FALSE.
       END IF
      END IF
      DEALLOCATE(BEST_COEF)
!
      IF (.NOT. IMPROVED) THEN
       LR = FACT*LR
       MAXLOOP = MAXLOOP+20
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Compute Gradient & Lagrangian Multipliers (ELAG)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL DENMATr(DEN,COEF,RO,NBF,1,NBF5) ! Density Matrix
      IF( (IPNOF==5.or.IPNOF==7.or.(IPNOF==8.and.Imod==0)) .and.        &
              MSpin==0 .and. IERITYP==1 )THEN
!AO    CALL ELAGaor(OVERLAP,AHCORE,IJKL,XIJKL,COEF,RO,DEN,ELAG)
       CALL ELAGr(AHCORE,IJKL,XIJKL,COEF,RO,DEN,ELAG,G,IPNOF)
      ELSE
       CALL ELAG1r(AHCORE,IJKL,XIJKL,XIJKaux,QD,COEF,RO,CJ12,CK12,ELAG,G)
      ENDIF
      DEALLOCATE(G,DEN)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Check convergence for the Energy
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Etot = EELEC + EN
      DIF_EELEC = EELEC - EELEC_OLD
      EELEC_OLD = EELEC
      PCONV = ABS(DIF_EELEC)
      CALL PCONVE(ELAG,DUMEL,MAXI,MAXJ,SUMDIF)
      CALL SYSTEM_CLOCK(TIMEFINISH)
      ORBTIME = (TIMEFINISH - TIMESTART)/RATE
      IF(IPRINTOPT==1)THEN
       WRITE(6,2)ITCALL,OCCTIME,ORBTIME,ILOOP,EELEC,Etot,DIF_EELEC,DUMEL
      END IF
      IF(PCONV<THRESHE .and. DUMEL<THRESHL) THEN
       CONVGDELAG=.TRUE.
       AUTOLR=.TRUE.
       LR=0.01D0
       MAXLOOP=10
      END IF
!-----------------------------------------------------------------------
!     FORMAT STATEMENTS
!-----------------------------------------------------------------------
    1 FORMAT(//2X,'HARTREE-FOCK',/1X,'===============',                 &
        //2X,'Ex It',2X,'Occ Time',2X,'Orb Time',2X,'Orb It',3X,        &
        'Electronic Energy',5X,'Total Energy',3X,'Energy Convergency',  &
        4X,'Max Mul-Lag Diff',/)
   10 FORMAT(//2X,'RESULTS OF OCCUPATION-COEFFICIENT OPTIMIZATION'      &
            ,/1X,'================================================',    &
        //2X,'Ex It',2X,'Occ Time',2X,'Orb Time',2X,'Orb It',3X,        &
        'Electronic Energy',5X,'Total Energy',3X,'Energy Convergency',  &
        4X,'Max Mul-Lag Diff',/)
    2 FORMAT(I5,3X,ES8.1,2X,ES8.1,3X,I5,1X,F18.10,1X,F19.10,2X,F15.10,  &
            8X,F11.6)
    3 FORMAT(2X,I3,'.',3X,F17.8,4X,F15.8,6X,F11.6)
!-----------------------------------------------------------------------
      RETURN
      END

!----------------------------------------------------------------------!
!                                                                      !
!     Subroutines for the Orbital Optimization using OrbOptFMIUGr      !
!                                                                      !
!   ELAG1r: Calculate the Lagrangian (ELAG)                            !
!   PCONVE: Check for the symmetry of Lagrangian [ ELAG(ij)-ELAG(ji) ] !
!   FFJMN1rc,ro: Calculate the generalized Fock matrix Fj(m,n)         !
!   ELG: Calculate the Lagrange Multipliers                            !
!   FORMFMIUG: Form & decrease gen-Fock matrix using a scaling factor  !
!   F01: Fij = Fij*0.1                                                 !
!   FFMIUG_DIIS: Direct Inversion in the Iterative Subspace Technique  !
!   TRACEFF: Calculate the trace by i,j of F(m,i,j)*F(idiis,j,i)       !
!                                                                      !
!----------------------------------------------------------------------!

! ELAG1r
      SUBROUTINE ELAG1r(AHCORE,IJKL,XIJKL,XIJKaux,QD,COEF,RO,CJ12,CK12, &
                        ELAG,G)
!-----------------------------------------------------------------------
!     Calculate the Lagrange Multipliers
!-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      LOGICAL ERIACTIVATED,HighSpin
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
      COMMON/ERIACT/ERIACTIVATED,NIJKaux,NINTCRaux,NSTOREaux,IAUXDIM
      COMMON/INPFILE_NIJKL/NINTMX,NIJKL,NINTCR,NSTORE
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      LOGICAL SMCD
      COMMON/ERITYPE/IERITYP,IRITYP,IGEN,ISTAR,MIXSTATE,SMCD
      COMMON/ENERGIAS/EELEC,EELEC_OLD,DIF_EELEC,EELEC_MIN
!
      INTEGER :: NSOC,NDNS,MSpin
      INTEGER(8),DIMENSION(NSTORE)::IJKL
      DOUBLE PRECISION,DIMENSION(NSTORE)::XIJKL
      DOUBLE PRECISION,DIMENSION(NSTOREaux)::XIJKaux
      DOUBLE PRECISION,DIMENSION(NBF5)::RO
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5)::CJ12,CK12
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::AHCORE,COEF,ELAG
      DOUBLE PRECISION,DIMENSION(NBF,NBF5)::G      
      DOUBLE PRECISION,DIMENSION(NBF,NBF,NBF)::QD
!
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::WJj,WKj,WF,AUX1
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: AUX2,AUX3
!-----------------------------------------------------------------------
!     Calculate Dj: QD(j,miu,niu), Jj(miu,niu), Kj(miu,niu) (j=1,NBF5) 
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE (WJj(NSQ,NBF5),WKj(NSQ,NBF5),AUX2(NSQ))
!
      IF(IERITYP==1 .or. (IERITYP==3 .and. MIXSTATE==2))THEN
       ALLOCATE(AUX1(NBF,NBF))
       DO j=1,NBF5
        CALL DENMATj(j,AUX1,COEF,NBF)
        QD(j,1:NBF,1:NBF) = AUX1(1:NBF,1:NBF) 
        CALL HSTARJ(AUX2,AUX1,IJKL,XIJKL)
        WJj(1:NSQ,j) = AUX2(1:NSQ)
        CALL HSTARK(AUX2,AUX1,IJKL,XIJKL)
        WKj(1:NSQ,j) = AUX2(1:NSQ)
       ENDDO
       DEALLOCATE(AUX1)
      ELSE IF(IERITYP==2 .or. (IERITYP==3 .and. MIXSTATE==1))THEN 
       ALLOCATE(AUX3(NSQ))
       DO j=1,NBF5
        CALL HSTARJKRI(AUX2,AUX3,XIJKaux,COEF(1:NBF,j))
        WJj(1:NSQ,j) = AUX2(1:NSQ)
        WKj(1:NSQ,j) = AUX3(1:NSQ)
       ENDDO
       DEALLOCATE(AUX3)
      END IF
      DEALLOCATE(AUX2)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Form F Matrix and keep it in WF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE (WF(NSQ,NBF5))
      if(MSpin==0)then
       CALL FFJMN1rc(RO,CJ12,CK12,AHCORE,WJj,WKj,WF)
      else if(MSpin>0)then
       CALL FFJMN1ro(RO,CJ12,CK12,AHCORE,WJj,WKj,WF)
      end if      
      DEALLOCATE (WJj,WKj)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Calculate G Matrix
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO iq=1,NBF5
       do i=1,nbf
        G(i,iq) = FC(i,iq,WF(1,iq),COEF)
       enddo
      ENDDO
      DEALLOCATE (WF)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Lagrangian Multipliers (ELAG)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL ELG(ELAG,COEF,G)
!-----------------------------------------------------------------------
      RETURN
      END

! PCONVE
      SUBROUTINE PCONVE(ELAG,DUM,MAXI,MAXJ,SUMDIF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::ELAG
!-----------------------------------------------------------------------
!     Check for the symmetry of Lagrangian [ ELAG(ij)-ELAG(ji) ], i.e.,
!     the convergence of the difference of Lagrange multipliers
!-----------------------------------------------------------------------
      DUM=0.0d0
      SUMDIF=0.0d0
      DO IQ=1,NBF
       DO JQ=1,NBF
        GCF=ABS(ELAG(IQ,JQ)-ELAG(JQ,IQ))
        SUMDIF=SUMDIF+GCF
        IF(GCF>DUM)THEN
         DUM=GCF
         MAXI=IQ
         MAXJ=JQ
        ENDIF
       ENDDO
      ENDDO
!     WRITE(6,'(2I4,2F10.5)')MAXI,MAXJ,ELAG(MAXI,MAXJ),ELAG(MAXJ,MAXI)
!-----------------------------------------------------------------------
      RETURN
      END

! FFJMN1rc
      SUBROUTINE FFJMN1rc(RO,CJ12,CK12,H,DJ,DK,F)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL HighSpin
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_ORBSPACE0/NO1,NDOC,NCO,NCWO,NVIR,NAC,NO0
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
#include "mpip.h"
      INTEGER :: NSOC,NDNS,MSpin
      DOUBLE PRECISION,DIMENSION(NBF5)::RO
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5)::CJ12,CK12
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::H
      DOUBLE PRECISION,DIMENSION(NSQ,NBF5)::DJ,DK,F
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::FF
      INTEGER LL,UL,EQPART,UNEQPART
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Wake up the nodes for the task
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#ifdef MPI
      DO I=1,NPROCS-1
        NOPT=6
        CALL MPI_SEND(NOPT,1,MPI_INTEGER,I,I,MPI_COMM_WORLD,IERR)
        CALL MPI_SEND(NBF,1,MPI_INTEGER,I,I,MPI_COMM_WORLD,IERR)
      ENDDO
      CALL MPI_BCAST(NBF5,1,MPI_INTEGER,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(DJ,NSQ*NBF5,MPI_REAL8,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(DK,NSQ*NBF5,MPI_REAL8,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(CJ12,NBF5*NBF5,MPI_REAL8,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(CK12,NBF5*NBF5,MPI_REAL8,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(RO,NBF5,MPI_REAL8,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(H,NBF*NBF,MPI_REAL8,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NA,1,MPI_INTEGER,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NB,1,MPI_INTEGER,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NO1,1,MPI_INTEGER,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NSOC,1,MPI_INTEGER,MASTER,MPI_COMM_WORLD,IERR)
      ALLOCATE(FF(NSQ,NBF5))
      F = 0.0D0
      FF=0.0D0
      EQPART = INT(NBF*NBF/NPROCS)
      UNEQPART = NBF*NBF - NPROCS*EQPART
      LL = 1
      UL = LL + EQPART - 1
      IF(0<UNEQPART) UL = UL + 1
#else
      LL = 1
      UL = NBF*NBF
#endif
!-----------------------------------------------------------------------
!     Calculate Fj(m,n)
!-----------------------------------------------------------------------
      IF(NO1>1)THEN
        !$OMP PARALLEL DO PRIVATE(i,k,ik,J)
        do ik=LL,UL
          i = INT((ik-1)/nbf) + 1
          k = MOD(ik-1,nbf) + 1
          F(ik,1) = H(i,k) + PRODWCWij(ik,DJ,CJ12)-PRODWCWij(ik,DK,CK12)
          DO J=NO1+1,NB
            F(ik,J) = RO(J) * ( H(i,k) + DJ(ik,J) )                       &
                  + PRODWCWijq(ik,J,DJ,CJ12) - PRODWCWijq(ik,J,DK,CK12)
          ENDDO
          if(NSOC>0)then
            DO J=NB+1,NA
              F(ik,J) = RO(J) * H(i,k)                                     &
                    + PRODWCWijq(ik,J,DJ,CJ12) - PRODWCWijq(ik,J,DK,CK12)
            ENDDO
          end if
          DO J=NA+1,NBF5
            F(ik,J) = RO(J) * ( H(i,k) + DJ(ik,J) )                       &
                  + PRODWCWijq(ik,J,DJ,CJ12) - PRODWCWijq(ik,J,DK,CK12)
          ENDDO
        enddo
        !$OMP END PARALLEL DO
        DO J=2,NO1
          F(1:NSQ,J) = F(1:NSQ,1)
        ENDDO
      ELSE
        !$OMP PARALLEL DO PRIVATE(i,k,ik,J)
        do ik=LL,UL
          i = INT((ik-1)/nbf) + 1
          k = MOD(ik-1,nbf) + 1
          DO J=1,NB
            F(ik,J) = RO(J) * ( H(i,k) + DJ(ik,J) )                       &
                  + PRODWCWijq(ik,J,DJ,CJ12) - PRODWCWijq(ik,J,DK,CK12)
          ENDDO
          if(NSOC>0)then
            DO J=NB+1,NA
              F(ik,J) = RO(J) * H(i,k)                                     &
                    + PRODWCWijq(ik,J,DJ,CJ12) - PRODWCWijq(ik,J,DK,CK12)
            ENDDO
          end if
          DO J=NA+1,NBF5
            F(ik,J) = RO(J) * ( H(i,k) + DJ(ik,J) )                       &
                  + PRODWCWijq(ik,J,DJ,CJ12) - PRODWCWijq(ik,J,DK,CK12)
          ENDDO
        enddo
        !$OMP END PARALLEL DO
      ENDIF
#ifdef MPI
      CALL MPI_REDUCE(F,FF,NSQ*NBF5,MPI_REAL8,MPI_SUM,MASTER,           &
            MPI_COMM_WORLD,IERR)
      F = FF
      DEALLOCATE(FF)
#endif
!-----------------------------------------------------------------------
      RETURN
      END
      
! FFJMN1ro
      SUBROUTINE FFJMN1ro(RO,CJ12,CK12,H,DJ,DK,F)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      COMMON/INPNOF_ORBSPACE0/NO1,NDOC,NCO,NCWO,NVIR,NAC,NO0
!
      DOUBLE PRECISION,DIMENSION(NBF5)::RO
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5)::CJ12,CK12
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::H
      DOUBLE PRECISION,DIMENSION(NSQ,NBF5)::DJ,DK,F
!-----------------------------------------------------------------------
!     Calculate Fj(m,n)
!-----------------------------------------------------------------------
      IF(NO1>1)THEN
       do i=1,nbf
        do k=1,nbf
         ik=i+(k-1)*nbf
         F(ik,1) = H(i,k)+PRODWCWij1(ik,DJ,CJ12)-PRODWCWij1(ik,DK,CK12) &
                 + SUMWij(ik,DJ) - 0.50d0*SUMWij(ik,DK)                  
         DO J=NO1+1,NB                                                   
          F(ik,J) = RO(J) * ( H(i,k) + DJ(ik,J) )                       &
                  + PRODWCWijq1(ik,J,DJ,CJ12)-PRODWCWijq1(ik,J,DK,CK12) &
                  + RO(J) * ( SUMWij(ik,DJ) - 0.50d0*SUMWij(ik,DK) )          
         ENDDO                                                           
         DO J=NB+1,NA                                                    
          F(ik,J) = 0.50d0 * ( H(i,k)+SUMWijq(ik,J,DJ)-SUMWijq(ik,J,DK) &
                      + 2.0d0*PRODWROij(ik,DJ,RO)-PRODWROij(ik,DK,RO) )  
         ENDDO                                                           
         DO J=NA+1,NBF5                                                  
          F(ik,J) = RO(J) * ( H(i,k) + DJ(ik,J) )                       &
                  + PRODWCWijq2(ik,J,DJ,CJ12)-PRODWCWijq2(ik,J,DK,CK12) &
                  + RO(J) * ( SUMWij(ik,DJ) - 0.50d0*SUMWij(ik,DK) )          
         ENDDO                                                           
        enddo                                                            
       enddo                                                             
       DO J=2,NO1                                                        
        F(1:NSQ,J) = F(1:NSQ,1)                                          
       ENDDO                                                             
      ELSE                                                               
       do i=1,nbf                                                        
        do k=1,nbf                                                       
         ik=i+(k-1)*nbf                                                  
         DO J=1,NB                                                       
          F(ik,J) = RO(J) * ( H(i,k) + DJ(ik,J) )                       &
                  + PRODWCWijq1(ik,J,DJ,CJ12)-PRODWCWijq1(ik,J,DK,CK12) &
                  + RO(J) * ( SUMWij(ik,DJ) - 0.50d0*SUMWij(ik,DK) )          
         ENDDO                                                           
         DO J=NB+1,NA                                                    
          F(ik,J) = 0.50d0 * ( H(i,k)+SUMWijq(ik,J,DJ)-SUMWijq(ik,J,DK) &
                      + 2.0d0*PRODWROij(ik,DJ,RO)-PRODWROij(ik,DK,RO) )  
         ENDDO                                                           
         DO J=NA+1,NBF5                                                  
          F(ik,J) = RO(J) * ( H(i,k) + DJ(ik,J) )                       &
                  + PRODWCWijq2(ik,J,DJ,CJ12)-PRODWCWijq2(ik,J,DK,CK12) &
                  + RO(J) * ( SUMWij(ik,DJ) - 0.50d0*SUMWij(ik,DK) )         

         ENDDO
        enddo
       enddo
      ENDIF
!-----------------------------------------------------------------------
      RETURN
      END

! ELG
      SUBROUTINE ELG(ELAG,C,G)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::ELAG,C
      DOUBLE PRECISION,DIMENSION(NBF,NBF5)::G
!-----------------------------------------------------------------------
      DO IQ=1,NBF
       DO JQ=1,NBF5
        ELAG(IQ,JQ) = 0.0d0
        do i=1,nbf
         ELAG(IQ,JQ) = ELAG(IQ,JQ) + C(i,IQ)*G(i,JQ)
        enddo
       ENDDO
       DO JQ=NBF5+1,NBF
        ELAG(IQ,JQ) = 0.0d0
       ENDDO
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! FORMFMIUG      
      SUBROUTINE FORMFMIUG(FMIUG,ELAG,FMIUG0,ITCALL) 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      LOGICAL RESTART,SCALING
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPSCALING/SCALING,NZEROS,NZEROSm,NZEROSr,ITZITER
      COMMON/INPNOF_RSTRT/RESTART,INPUTGAMMA,INPUTC,INPUTFMIUG,INPUTCXYZ
      DOUBLE PRECISION,DIMENSION(NBF)::FMIUG0
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::FMIUG,ELAG
!-----------------------------------------------------------------------
!     Generalized Fock Matrix (FMIUG)
!-----------------------------------------------------------------------
      IF(ITCALL==1.AND.INPUTFMIUG==0)THEN          ! only for itcall==1
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       DO IQ=1,NBF
        DO JQ=1,IQ-1
         FMIUG(IQ,JQ)=(ELAG(IQ,JQ)+ELAG(JQ,IQ))/2.0        ! Nondiagonal 
         FMIUG(JQ,IQ)=FMIUG(IQ,JQ)
        ENDDO
        FMIUG(IQ,IQ)=ELAG(IQ,IQ)                            ! Diagonal
       ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ELSE
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       if(SCALING)then
        DO IQ=1,NBF
         DO JQ=1,IQ-1
          FMIUG(IQ,JQ)=ELAG(IQ,JQ)-ELAG(JQ,IQ)           ! Nondiagonal 
!-.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.-
!         Decrease FMIUG using a scaling factor
!         The scaling factor varies until the number of
!         ZEROS (.000##) is equal for all elements Fij
!-.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.-
          CALL F01(NZEROS+9,FMIUG(IQ,JQ))
!-.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.-
          FMIUG(JQ,IQ)=FMIUG(IQ,JQ)                      ! Fji=Fij
         ENDDO
         FMIUG(IQ,IQ)=FMIUG0(IQ)                         ! Diagonal
        ENDDO
       else
        DO IQ=1,NBF
         DO JQ=1,IQ-1
          FMIUG(IQ,JQ)=ELAG(IQ,JQ)-ELAG(JQ,IQ)           ! Nondiagonal 
          FMIUG(JQ,IQ)=FMIUG(IQ,JQ)                      ! Fji=Fij
         ENDDO
         FMIUG(IQ,IQ)=FMIUG0(IQ)                         ! Diagonal
        ENDDO
       endif
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ENDIF
!-----------------------------------------------------------------------
      RETURN
      END
            
! F01
      SUBROUTINE F01(imax,Fij)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      do i=0,imax
       VAL=ABS(Fij)
       if(VAL>10.0d0**(9-i).and.VAL<10.0d0**(10-i))then
        Fij = Fij * 0.1d0
       endif
      enddo
      RETURN
      END

! FFMIUG_DIIS
      SUBROUTINE FFMIUG_DIIS(NUM,FMIUG,CFM,BFM,FK,IDIIS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      LOGICAL DIIS,PERDIIS
      COMMON/INPNOF_DIIS/DIIS,PERDIIS,NDIIS,NTHDIIS,THDIIS
      COMMON/INPNOF_COEFOPT/MAXLOOP
!
      DOUBLE PRECISION,DIMENSION(MAXLOOP+1)::CFM
      DOUBLE PRECISION,DIMENSION(MAXLOOP+1,MAXLOOP+1)::BFM
      DOUBLE PRECISION,DIMENSION(NUM,NUM)::FMIUG
      DOUBLE PRECISION,DIMENSION(MAXLOOP,NUM,NUM)::FK
!
      INTEGER,ALLOCATABLE,DIMENSION(:)::IPIV
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::A
!-----------------------------------------------------------------------
      IDIIS = IDIIS+1
      FK(IDIIS,1:NUM,1:NUM) = FMIUG(1:NUM,1:NUM)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Form BFM(1:idiis,idiis)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      NN = idiis+1
      do m=1,idiis
       BFM(m,idiis) = TRACEFF(MAXLOOP,NUM,m,idiis,FK)
       BFM(idiis,m) = BFM(m,idiis)
       BFM(m,NN)    = -1.0d0
       BFM(NN,m)    = -1.0d0
      enddo
      BFM(NN,NN) = 0.0d0
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(idiis>NDIIS)then
!- - - - - - - - - - - - - - - - - - - - - - -
       ALLOCATE (A(NN,NN),IPIV(NN))
       IPIV = 0
       A = BFM(1:NN,1:NN)
       CFM(1:idiis) = 0.0d0
       CFM(NN) =  -1.0d0
!- - - - - - - - - - - - - - - - - - - - - - -
       CALL DGESV(NN,1,A,NN,IPIV,CFM,NN,INFO)
       DEALLOCATE (A,IPIV)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      New Generalized Fock Matrix: F'=SUM_k [CFM(k)*F(k)]
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       do i=1,NUM
        do j=1,i-1
         FMIUG(i,j) = 0.0d0
         do k=1,idiis
          FMIUG(i,j) = FMIUG(i,j) + CFM(k) * FK(k,i,j)
         enddo
         FMIUG(j,i) = FMIUG(i,j)
        enddo
       enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      IDIIS is nullified for a periodic DIIS
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       IF(PERDIIS)IDIIS = 0
!- - - - - - - - - - - - - - - - - - - - - - -
      endif
!-----------------------------------------------------------------------
      RETURN
      END

! TRACEFF
      FUNCTION TRACEFF(maxloop,N,m,idiis,F)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DOUBLE PRECISION,DIMENSION(maxloop,N,N)::F
!-----------------------------------------------------------------------
      TRACEFF = 0.0d0
      do i=1,N
       do j=1,i-1
        TRACEFF = TRACEFF + F(m,i,j)*F(idiis,j,i)
       enddo
      enddo
!-----------------------------------------------------------------------
      RETURN
      END

!----------------------------------------------------------------------!
!                                                                      !
!     Subroutines for the Orbital Optimization using Density Matrices  !
!                        (PNOF5, PNOF7 and GNOF)                       !
!                                                                      !
!   ELAGr:   Calculate the Lagrangian (ELAG)                           !
!   DENMATg: Calculate Density Matrix for a given subspace 'l'         !
!   HSTARK3: Determine three skeleton K from atomic integrals (AUX)    !
!                                                                      !
!----------------------------------------------------------------------!

! ELAGr (MSpin=0,IERITYP=1)
      SUBROUTINE ELAGr(AHCORE,IJKL,XIJKL,COEF,RO,DEN,ELAG,G,IPNOF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/INPFILE_NIJKL/NINTMX,NIJKL,NINTCR,NSTORE      
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5      
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_ORBSPACE0/NO1,NDOC,NCO,NCWO,NVIR,NAC,NO0
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
      COMMON/INPNOF_STATIC/Ista
      COMMON/INPNOF_MOD/Imod
      INTEGER :: NSOC,NDNS,MSpin
      LOGICAL :: HighSpin
      INTEGER(8),DIMENSION(NSTORE) :: IJKL
      DOUBLE PRECISION,DIMENSION(NSTORE) :: XIJKL
      DOUBLE PRECISION,DIMENSION(NBF5) :: RO
      DOUBLE PRECISION,DIMENSION(NBF,NBF) :: AHCORE,COEF,DEN,ELAG
      DOUBLE PRECISION,DIMENSION(NBF,NBF5) :: G
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)   :: BETA,FIs,ROd,Rd
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: DM1,DM2,DM3
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: WF1,WF2,WF3
!-----------------------------------------------------------------------!
!                               P N O F 5                               !
!-----------------------------------------------------------------------!
      ALLOCATE(DM1(NBF,NBF),DM2(NBF,NBF),DM3(NBF,NBF))      
      ALLOCATE(WF1(NBF,NBF),WF2(NBF,NBF),WF3(NBF,NBF))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Fock Matrix & Hartree-Fock Component of the G Matrix
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL FORM2JK(WF1,DEN,IJKL,XIJKL)
      WF2 = AHCORE + WF1  ! Fockian
      DO iq=1,NBF5
       do i=1,nbf
        G(i,iq) = RO(iq) * FC(i,iq,WF2,COEF)
       enddo
      ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Intrapair Component of the G Matrix [ Sum_g (Eg) ]
!     Eliminating the Intrapair Component associated to RO (RHF)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(BETA(NBF5))
      BETA = - DSQRT(RO)
      DO i=1,NB
       BETA(i) = DSQRT(RO(i))
      END DO      
      DO l=1,NDOC
       ig = NO1+l
       CALL DENMATg(l,DM1,COEF,BETA)           
       CALL HSTARK(WF1,DM1,IJKL,XIJKL)  ! WK(l),beta
       CALL DENMATg(l,DM2,COEF,RO)              
       CALL FORM2JK(WF2,DM2,IJKL,XIJKL) ! 2JK(l),ro       
       do i=1,nbf
        G(i,ig) = G(i,ig) + BETA(ig)*FC(i,ig,WF1,COEF)                  &
                          - RO(ig) * FC(i,ig,WF2,COEF)
       enddo
       do j=NDNS+NDOC-l+1,NDNS+NDOC-l+1+(NCWO-1)*NDOC,NDOC
        ip = NO1+j
        do i=1,nbf
         G(i,ip) = G(i,ip) + BETA(ip)*FC(i,ip,WF1,COEF)                 &
                           - RO(ip) * FC(i,ip,WF2,COEF)
        enddo
       end do
      END DO      
      DEALLOCATE(BETA)      
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Eliminating the Component associated to the Multiplet (MSpin=0)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(NSOC>0)then
       CALL DENMATr(DM3,COEF,RO,NBF,NB+1,NA)       
       CALL HSTARK(WF1,DM3,IJKL,XIJKL) 
       DO ig=NB+1,NA
        do i=1,nbf
         G(i,ig) = G(i,ig) - RO(ig) * FC(i,ig,WF1,COEF)
        enddo
       ENDDO
      endif
!-----------------------------------------------------------------------!
!                               P N O F 7                               !
!-----------------------------------------------------------------------!
      IF(IPNOF==7)THEN
       ALLOCATE(FIs(NBF5))
       if(Ista==0)then
        FIs = DSQRT(RO*(1.0d0-RO))
       else
        FIs = 2.0d0*RO*(1.0d0-RO)
       endif
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -       
!      Additional Density & Fock & G Matrices associated to FIs
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
       CALL DENMATr(DM1,COEF,FIs,NBF,1,NBF5)
       CALL HSTARK(WF1,DM1,IJKL,XIJKL)        ! WK,fis
       DO iq=1,NBF5
        do i=1,nbf
         G(i,iq) = G(i,iq) - FIs(iq) * FC(i,iq,WF1,COEF)
        enddo
       ENDDO      
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Eliminating the Intrapair Component associated to FIs
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       DO l=1,NDOC
        ig = NO1+l
        CALL DENMATg(l,DM2,COEF,FIs)        
        CALL HSTARK(WF1,DM2,IJKL,XIJKL)    ! WK(l),fis
        do i=1,nbf
         G(i,ig) = G(i,ig) + FIs(ig) * FC(i,ig,WF1,COEF)
        enddo
        do j=NDNS+NDOC-l+1,NDNS+NDOC-l+1+(NCWO-1)*NDOC,NDOC
         ip = NO1+j
         do i=1,nbf
          G(i,ip) = G(i,ip) + FIs(ip) * FC(i,ip,WF1,COEF)         
         enddo
        end do
       END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Eliminating the Component associated to the Multiplet (MSpin=0)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       if(NSOC>0)then
        CALL DENMATr(DM3,COEF,FIs,NBF,NB+1,NA)       
        CALL HSTARK(WF1,DM3,IJKL,XIJKL) 
        DO ig=NB+1,NA
         do i=1,nbf
          G(i,ig) = G(i,ig) + FIs(ig) * FC(i,ig,WF1,COEF)
         enddo
        ENDDO
       endif
       DEALLOCATE(FIs)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -       
      ENDIF
!-----------------------------------------------------------------------!
!                                G N O F                                !
!-----------------------------------------------------------------------!
      IF(IPNOF==8 .AND. Imod==0)THEN
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Dynamic & Static Occupation Numbers (ROd,Rd,FIs)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       ALLOCATE(ROd(NBF5))
       ROd = 0.0d0
       Hcut = 0.02d0*DSQRT(2.0d0)
       DO i=1,NDOC
        in = NO1+i      !NO1+1:NB
        HOLEin = 1.0d0-RO(in)
        COC = HOLEin/Hcut       
        Fin = DEXP(-COC*COC)            
        ROd(in) = RO(in) * Fin
        IF(NCWO>1)THEN  !NA+1:NBF5
         do iw=1,ncwo                                                                   
          im = no1+(na-nb)+ndoc*(iw+1)-i+1
          ROd(im) = RO(im) * Fin   
         enddo
        ELSE      !perfect-pairing                       
         icf = na+ndoc-i+1
         ROd(icf) = RO(icf) * Fin  
        ENDIF
       ENDDO
!- - - - - - - - - - - - - - - - - -
       ALLOCATE(Rd(NBF5))
       Rd = - DSQRT(ROd)       
       DO i=1,NB
        Rd(i) = DSQRT(ROd(i))
       END DO
!- - - - - - - - - - - - - - - - - -
       ALLOCATE(FIs(NBF5))
       FIs = DSQRT(RO*(1.0d0-RO))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -       
!      Additional Density & Fock & G Matrices associated to FIs,ROd & Rd
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
       CALL DENMATr(DM1,COEF,FIs,NBF,1,NBF5)
       CALL DENMATr(DM2,COEF,ROd,NBF,1,NBF5)       
       CALL DENMATr(DM3,COEF,Rd ,NBF,1,NBF5)
       CALL HSTARK3(WF1,DM1,WF2,DM2,WF3,DM3,IJKL,XIJKL)  ! WFi
       DO iq=1,NBF5
        do i=1,nbf
         G(i,iq) = G(i,iq) - FIs(iq) * FC(i,iq,WF1,COEF)                &
                           + ROd(iq) * FC(i,iq,WF2,COEF)                &
                           + Rd(iq)  * FC(i,iq,WF3,COEF)                           
        enddo
       ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Eliminating the Intrapair Component associated to FIs,ROd & Rd
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       DO l=1,NDOC
        ig = NO1+l
        CALL DENMAT0g(l,DM1,COEF,FIs)
        CALL DENMAT0g(l,DM2,COEF,ROd)        
        CALL DENMAT0g(l,DM3,COEF,Rd)        
        CALL HSTARK3(WF1,DM1,WF2,DM2,WF3,DM3,IJKL,XIJKL) ! WFi(l)
        do i=1,nbf
         G(i,ig) = G(i,ig) + FIs(ig) * FC(i,ig,WF1,COEF)               &
                           - ROd(ig) * FC(i,ig,WF2,COEF)               &
                           - Rd(ig)  * FC(i,ig,WF3,COEF)
        enddo
        CALL DENMATg(l,DM1,COEF,FIs)        
        CALL DENMATg(l,DM2,COEF,ROd)       
        CALL DENMATg(l,DM3,COEF,Rd)
        CALL HSTARK3(WF1,DM1,WF2,DM2,WF3,DM3,IJKL,XIJKL) ! WFi(l)        
        do j=NDNS+NDOC-l+1,NDNS+NDOC-l+1+(NCWO-1)*NDOC,NDOC
         ip = NO1+j
         do i=1,nbf
          G(i,ip) = G(i,ip) + FIs(ip) * FC(i,ip,WF1,COEF)               &
                            - ROd(ip) * FC(i,ip,WF2,COEF)               &
                            - Rd(ip)  * FC(i,ip,WF3,COEF)
         enddo
        end do
       END DO       
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Eliminating the Below-Below Contribution (up to NB)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       CALL DENMATr(DM1,COEF,FIs,NBF,1,NB)       
       CALL DENMATr(DM2,COEF,ROd,NBF,1,NB)              
       CALL DENMATr(DM3,COEF,Rd ,NBF,1,NB)              
       CALL HSTARK3(WF1,DM1,WF2,DM2,WF3,DM3,IJKL,XIJKL) ! WFi
       DO iq=1,NB
        do i=1,nbf
         G(i,iq) = G(i,iq) + FIs(iq) * FC(i,iq,WF1,COEF)                &
                           - ROd(iq) * FC(i,iq,WF2,COEF)                &
                           - Rd(iq)  * FC(i,iq,WF3,COEF)                                                      
        enddo
       ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Eliminating the Component associated to the Multiplet (MSpin=0)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       if(NSOC>0)then
        DO iq=NB+1,NA
         do i=1,nbf
          G(i,iq) = G(i,iq) + 0.5d0 * FIs(iq) * FC(i,iq,WF1,COEF)
         enddo
        ENDDO
        CALL DENMATr(DM3,COEF,FIs,NBF,NB+1,NA)       
        CALL HSTARK(WF1,DM3,IJKL,XIJKL) 
        DO iq=1,NB
         do i=1,nbf
          G(i,iq) = G(i,iq) + 0.5d0 * FIs(iq) * FC(i,iq,WF1,COEF)
         enddo
        ENDDO
        DO ig=NB+1,NA
         do i=1,nbf
          G(i,ig) = G(i,ig) + FIs(ig) * FC(i,ig,WF1,COEF)
         enddo
        ENDDO
       endif
       DEALLOCATE(ROd,Rd,FIs)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ENDIF      
!-----------------------------------------------------------------------                                 
!     Lagrangian Multipliers (ELAG)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                 
      CALL ELG(ELAG,COEF,G)
      DEALLOCATE(DM1,DM2,DM3,WF1,WF2,WF3)      
!-----------------------------------------------------------------------
      RETURN                                     
      END

! DENMATg
      SUBROUTINE DENMATg(l,DMl,COEF,RO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      LOGICAL HighSpin      
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_ORBSPACE0/NO1,NDOC,NCO,NCWO,NVIR,NAC,NO0
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5      
      DOUBLE PRECISION,DIMENSION(NBF,NBF) :: DMl,COEF
      DOUBLE PRECISION,DIMENSION(NBF5) :: RO      
!-----------------------------------------------------------------------
      ig = NO1+l
      do m=1,NBF
       do n=m,NBF
        DMl(m,n) = RO(ig)*COEF(m,ig)*COEF(n,ig)
        do i=NDNS+NDOC-l+1,NDNS+NDOC-l+1+(NCWO-1)*NDOC,NDOC
         ip = NO1+i
         DMl(m,n)= DMl(m,n) + RO(ip)*COEF(m,ip)*COEF(n,ip)
        end do
        DMl(n,m) = DMl(m,n) 
       end do
      end do
      DMl = 2.0d0*DMl
!-----------------------------------------------------------------------      
      RETURN
      END
      
! DENMAT0g
      SUBROUTINE DENMAT0g(l,DMl,COEF,RO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_ORBSPACE0/NO1,NDOC,NCO,NCWO,NVIR,NAC,NO0
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5      
      INTEGER :: NSOC,NDNS,MSpin
      LOGICAL :: HighSpin
      DOUBLE PRECISION,DIMENSION(NBF,NBF) :: DMl,COEF
      DOUBLE PRECISION,DIMENSION(NBF5) :: RO      
!-----------------------------------------------------------------------
      do m=1,NBF
       do n=m,NBF
        DMl(m,n) = 0.00
        do i=NDNS+NDOC-l+1,NDNS+NDOC-l+1+(NCWO-1)*NDOC,NDOC
         ip = NO1+i
         DMl(m,n)= DMl(m,n) + RO(ip)*COEF(m,ip)*COEF(n,ip)
        end do
        DMl(n,m) = DMl(m,n) 
       end do
      end do
      DMl = 2.0d0*DMl
!-----------------------------------------------------------------------      
      RETURN
      END      

! HSTARK3
      SUBROUTINE HSTARK3(FM1,PM1,FM2,PM2,FM3,PM3,IERI,ERI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPFILE_NIJKL/NINTMX,NIJKL,NINTCR,NSTORE
      COMMON/USEHUBBARD/IHUB
#include "mpip.h"
      INTEGER(8),DIMENSION(NSTORE) :: IERI
      INTEGER(8) :: LABEL
      DOUBLE PRECISION,DIMENSION(NSTORE) :: ERI
      DOUBLE PRECISION,DIMENSION(NBF,NBF):: FM1,PM1,FM2,PM2,FM3,PM3
      ALLOCATABLE :: P1(:),P2(:),P3(:),F1(:),F2(:),F3(:)
#ifdef MPI 
      ALLOCATABLE :: FF1(:),FF2(:),FF3(:)
#endif
      ALLOCATE (P1(NBFT),F1(NBFT),P2(NBFT),F2(NBFT),P3(NBFT),F3(NBFT))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Wake up the nodes for the task
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#ifdef MPI
      ALLOCATE (FF1(NBFT),FF2(NBFT),FF3(NBFT))
      DO I=1,NPROCS-1
        NOPT=5 
        CALL MPI_SEND(NOPT,1,MPI_INTEGER,I,I,MPI_COMM_WORLD,IERR)
        CALL MPI_SEND(NBFT,1,MPI_INTEGER,I,I,MPI_COMM_WORLD,IERR)
      ENDDO
#endif
      CALL SQUARETRIAN(PM1,P1,NBF,NBFT)
      CALL SQUARETRIAN(PM2,P2,NBF,NBFT)
      CALL SQUARETRIAN(PM3,P3,NBF,NBFT)
      F1 = 0.0d0
      F2 = 0.0d0
      F3 = 0.0d0
#ifdef MPI 
      FF1 = 0.0d0
      FF2 = 0.0d0
      FF3 = 0.0d0
      CALL MPI_BCAST(NBF,1,MPI_INTEGER,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(P1,NBFT,MPI_REAL8,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(P2,NBFT,MPI_REAL8,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(P3,NBFT,MPI_REAL8,MASTER,MPI_COMM_WORLD,IERR)
#endif
      !$OMP PARALLEL DO PRIVATE(LABEL, I, J, K, L, XJ, XK, NIJ, NKL, NIK, NJL, NIL, NJK)  REDUCTION(+:F1,F2,F3)
      DO M=1,NINTCR
        LABEL = IERI(M)
        CALL LABELIJKL(LABEL,I,J,K,L)
        XJ = ERI(M)
        NIJ = I*(I-1)/2 + J
        NKL = K*(K-1)/2 + L
        XJ = 0.25*XJ
        IF(IHUB==0)CALL OTTOINTEGR(I,J,K,L,NIJ,NKL,XJ)
        XK = XJ
        NIK = I*(I-1)/2 + K 
        NJL = MAX0(J,L)*(MAX0(J,L)-1)/2 + MIN0(J,L)

        IF(I==K.OR.J==L) XK=XK+XK
        F1(NIK)=F1(NIK)+P1(NJL)*XK
        F2(NIK)=F2(NIK)+P2(NJL)*XK
        F3(NIK)=F3(NIK)+P3(NJL)*XK

        IF(NIK/=NJL)THEN
          F1(NJL)=F1(NJL)+P1(NIK)*XK
          F2(NJL)=F2(NJL)+P2(NIK)*XK
          F3(NJL)=F3(NJL)+P3(NIK)*XK
        ENDIF

        IF(I/=J.and.K/=L)THEN
          NIL = I*(I-1)/2 + L
          NJK = MAX0(J,K)*(MAX0(J,K)-1)/2 + MIN0(J,K)
          IF(I==L.OR.J==K) XJ=XJ+XJ

          F1(NIL)=F1(NIL)+P1(NJK)*XJ
          F2(NIL)=F2(NIL)+P2(NJK)*XJ
          F3(NIL)=F3(NIL)+P3(NJK)*XJ

          IF(NIL/=NJK)THEN
            F1(NJK)=F1(NJK)+P1(NIL)*XJ
            F2(NJK)=F2(NJK)+P2(NIL)*XJ
            F3(NJK)=F3(NJK)+P3(NIL)*XJ
          ENDIF
        ENDIF
      ENDDO
      !$OMP END PARALLEL DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Get the pieces from slaves
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#ifdef MPI 
      CALL MPI_REDUCE(F1,FF1,NBFT,MPI_REAL8,MPI_SUM,MASTER,             &
            MPI_COMM_WORLD,IERR)
      CALL MPI_REDUCE(F2,FF2,NBFT,MPI_REAL8,MPI_SUM,MASTER,             &
            MPI_COMM_WORLD,IERR)
      CALL MPI_REDUCE(F3,FF3,NBFT,MPI_REAL8,MPI_SUM,MASTER,             &
            MPI_COMM_WORLD,IERR)

      CALL TRIANSQUARE(FM1,FF1,NBF,NBFT)
      CALL TRIANSQUARE(FM2,FF2,NBF,NBFT)
      CALL TRIANSQUARE(FM3,FF3,NBF,NBFT)

      DEALLOCATE (P1,P2,P3,F1,F2,F3,FF1,FF2,FF3)
#else
      CALL TRIANSQUARE(FM1,F1,NBF,NBFT)
      CALL TRIANSQUARE(FM2,F2,NBF,NBFT)
      CALL TRIANSQUARE(FM3,F3,NBF,NBFT)

      DEALLOCATE (P1,P2,P3,F1,F2,F3)
#endif
!----------------------------------------------------------------------
      RETURN
      END

!----------------------------------------------------------------------!
!   ELAGaor: Calculate the Lagrangian (ELAGao) in AO Basis             !
!   FTERTRA: Interpair - Intrapair Components (- Below-Below if GNOF)  !
!   DENMATg,0g: Calculate Density Matrix for a given subspace 'l'      !
!----------------------------------------------------------------------!

! ELAGaor (MSpin=0,IERITYP=1)
      SUBROUTINE ELAGaor(OVERLAP,AHCORE,IJKL,XIJKL,COEF,RO,DEN,ELAG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL HighSpin
      COMMON/INPFILE_NIJKL/NINTMX,NIJKL,NINTCR,NSTORE
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5      
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_ORBSPACE0/NO1,NDOC,NCO,NCWO,NVIR,NAC,NO0
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
      COMMON/INPNOF_PNOF/IPNOF,NTWOPAR      
      COMMON/INPNOF_STATIC/Ista      
      COMMON/INPNOF_MOD/Imod
!
      INTEGER :: NSOC,NDNS,MSpin
      INTEGER(8),DIMENSION(NSTORE)        :: IJKL
      DOUBLE PRECISION,DIMENSION(NSTORE)  :: XIJKL
      DOUBLE PRECISION,DIMENSION(NBF5)    :: RO      
      DOUBLE PRECISION,DIMENSION(NBF,NBF) :: OVERLAP,AHCORE
      DOUBLE PRECISION,DIMENSION(NBF,NBF) :: COEF,DEN,ELAG
!
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)   :: AA,FI,ROd,Rd      
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: DM,DMS,WF,AUX 
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: ELAGao,COEFt
!-----------------------------------------------------------------------!
!                               P N O F 5                               !
!-----------------------------------------------------------------------!
      ALLOCATE(DM(NBF,NBF),DMS(NBF,NBF),WF(NBF,NBF))
      ALLOCATE(AUX(NBF,NBF),ELAGao(NBF,NBF),COEFt(NBF,NBF))
      COEFt = TRANSPOSE(COEF)     
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -      
!                   Interpair Component: Fock Matrix                    
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DM = DEN                           ! DM
      CALL FORM2JK(WF,DM,IJKL,XIJKL)     ! Skeleton of the Fockian (2JK)
      AUX = AHCORE + WF                  ! Fockian (F=H+2JK)
      DMS = MATMUL(DM,OVERLAP)           ! DM*S      
      ELAGao = MATMUL(AUX,DMS)           ! Lagrangian Multipliers in AO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Eliminating the Intrapair Component: Sum_g [ (2JK)g * DMg * S ]
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO lg=1,NDOC
       CALL DENMATg(lg,DM,COEF,RO)       ! DMg 
       CALL FORM2JK(WF,DM,IJKL,XIJKL)    ! (2JK)g
       DMS = MATMUL(DM,OVERLAP)          ! DMg*S       
       ELAGao = ELAGao - MATMUL(WF,DMS) 
      END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Eliminating the Component associated to the Multiplet (MSpin=0)      
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -      
      if(NSOC>0)then      
       CALL DENMATr(DM,COEF,RO,NBF,NB+1,NA) ! DM_op            
       CALL HSTARK(WF,DM,IJKL,XIJKL)        ! K_op
       DMS = MATMUL(DM,OVERLAP)             ! DM_op*S
       ELAGao = ELAGao - MATMUL(WF,DMS)  
      endif      
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!              Intrapair Component: Sum_g [ Ka_g * Ag * S )             
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Natural Amplitudes    
!- - - - - - - - - - - - - -
      ALLOCATE(AA(NBF5))
      AA = - DSQRT(RO)               
      DO i=1,NB
       AA(i) = DSQRT(RO(i))
      END DO
!- - - - - - - - - - - - - - - - - - - - -       
      DO lg=1,NDOC
       CALL DENMATg(lg,DM,COEF,AA)       ! Ag
       CALL HSTARK(WF,DM,IJKL,XIJKL)     ! Ka_g
       DMS = MATMUL(DM,OVERLAP)          ! Ag*S       
       ELAGao = ELAGao + MATMUL(WF,DMS)  
      END DO       
!-----------------------------------------------------------------------!
!                               P N O F 7                               !
!-----------------------------------------------------------------------!
      IF(IPNOF==7)THEN
       ALLOCATE(FI(NBF5))      
       if(Ista==0)then
        FI = DSQRT(RO*(1.0d0-RO))
       else
        FI = 2.0d0*RO*(1.0d0-RO)
       endif
!      Interpair Component: - Sum (K*FI*FI)                 
       CALL FTERTRA(-1.0,FI,COEF,DM,DMS,WF,IJKL,XIJKL,OVERLAP,ELAGao)
      END IF
!-----------------------------------------------------------------------!
!                                G N O F                                !
!-----------------------------------------------------------------------!
      IF(IPNOF==8 .AND. Imod==0)THEN
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
       ALLOCATE(FI(NBF5),ROd(NBF5),Rd(NBF5))      
!       
       FI = DSQRT(RO*(1.0d0-RO))
!       
       ROd = 0.0d0
       Hcut = 0.02d0*DSQRT(2.0d0)
       DO i=1,NDOC
        in = NO1+i      !NO1+1:NB
        HOLEin = 1.0d0-RO(in)
        COC = HOLEin/Hcut       
        Fin = DEXP(-COC*COC)            
        ROd(in) = RO(in) * Fin
        IF(NCWO>1)THEN !NA+1:NBF5
         do iw=1,ncwo                                                                   
          im = no1+(na-nb)+ndoc*(iw+1)-i+1
          ROd(im) = RO(im) * Fin   
         enddo
        ELSE     !perfect-pairing                       
         icf = na+ndoc-i+1
         ROd(icf) = RO(icf) * Fin  
        ENDIF
       ENDDO
!       
       Rd = - DSQRT(ROd)       
       DO i=1,NB
        Rd(i) = DSQRT(ROd(i))
       END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Interpair Component: - Sum (K*FI*FI)                             
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       CALL FTERTRA(-1.0,FI,COEF,DM,DMS,WF,IJKL,XIJKL,OVERLAP,ELAGao)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Interpair Dynamic Component: + Sum (K*ROd*ROd)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       CALL FTERTRA(1.0,ROd,COEF,DM,DMS,WF,IJKL,XIJKL,OVERLAP,ELAGao)
!- - - - - - - - - - - - - - - - - - - - - - -
!      Interpair Dynamic Component [ + Sum (K*Rd*Rd) ]
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       CALL FTERTRA(1.0,Rd,COEF,DM,DMS,WF,IJKL,XIJKL,OVERLAP,ELAGao)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ENDIF      
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -      
      ELAGao = ELAGao / 2.0d0              ! Since DMs are * 2.0d0
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Calculate Lagrangian Multipliers in NO basis
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      AUX = MATMUL(ELAGao,COEF)
      ELAG = MATMUL(COEFt,AUX)
!-----------------------------------------------------------------------
      DEALLOCATE(DM,DMS,WF,AUX,ELAGao,AA,COEFt)
      IF(IPNOF==7)DEALLOCATE(FI)
      IF(IPNOF==8.AND.Imod==0)DEALLOCATE(FI,ROd,Rd)
      RETURN                                     
      END

! FTERTRA
      SUBROUTINE FTERTRA(SGN,FI,COEF,DM,DMS,WF,IJKL,XIJKL,OVERLAP,ELAGao)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL HighSpin      
      COMMON/INPFILE_NIJKL/NINTMX,NIJKL,NINTCR,NSTORE 
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5 
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_ORBSPACE0/NO1,NDOC,NCO,NCWO,NVIR,NAC,NO0
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
      COMMON/INPNOF_PNOF/IPNOF,NTWOPAR      
      COMMON/INPNOF_MOD/Imod
!
      INTEGER(8),DIMENSION(NSTORE)        :: IJKL
      DOUBLE PRECISION,DIMENSION(NSTORE)  :: XIJKL
      DOUBLE PRECISION,DIMENSION(NBF5)    :: FI
      DOUBLE PRECISION,DIMENSION(NBF,NBF) :: COEF,DM,DMS
      DOUBLE PRECISION,DIMENSION(NBF,NBF) :: WF,OVERLAP,ELAGao       
!-----------------------------------------------------------------------
!     Interpair Component: - Sum (K*FI*FI)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL DENMATr(DM,COEF,FI,NBF,1,NBF5)     ! Dfi
      CALL HSTARK(WF,DM,IJKL,XIJKL)           ! Kfi
      DMS = MATMUL(DM,OVERLAP)                ! Dfi*S
      ELAGao = ELAGao + SGN*MATMUL(WF,DMS)        
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Eliminating the Intrapair Component: Sum_g [ Kfi_g * Dfi_g * S ]
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO lg=1,NDOC
       CALL DENMATg(lg,DM,COEF,FI)            ! Dfi_g
       CALL HSTARK(WF,DM,IJKL,XIJKL)          ! Kfi_g
       DMS = MATMUL(DM,OVERLAP)               ! Dfi_g*S       
       ELAGao = ELAGao - SGN*MATMUL(WF,DMS) 
      END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Eliminating the Component associated to the Multiplet (MSpin=0)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(NSOC>0.and.(IPNOF==8 .AND. Imod==0))stop ! needs change (1/2, see cjckd8) !!!
      if(NSOC>0)then      
       CALL DENMATr(DM,COEF,FI,NBF,NB+1,NA)   ! Dfi_op            
       CALL HSTARK(WF,DM,IJKL,XIJKL)          ! Kfi_op
       DMS = MATMUL(DM,OVERLAP)               ! DM_op*S
       ELAGao = ELAGao - SGN*MATMUL(WF,DMS)  
      endif
!-----------------------------------------------------------------------
      IF(IPNOF==8 .AND. Imod==0)THEN
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Eliminating the Below-Below Contribution (up to NB)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       CALL DENMATr(DM,COEF,FI,NBF,1,NB)       ! Dfi_nb
       CALL HSTARK(WF,DM,IJKL,XIJKL)           ! Kfi_nb
       DMS = MATMUL(DM,OVERLAP)                ! Dfi*S_nb
       ELAGao = ELAGao - SGN*MATMUL(WF,DMS)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Eliminating the Intrapair Component
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       DO lg=1,NDOC
        ig = NO1+lg
        CALL DENMATj(ig,DM,COEF,NBF)           ! 2CCg
        CALL HSTARK(WF,DM,IJKL,XIJKL)          ! Kg
        DMS = MATMUL(DM,OVERLAP)               ! CCg*S       
        ELAGao = ELAGao + SGN*FI(ig)*FI(ig)*MATMUL(WF,DMS) 
       END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -       
      END IF       
!-----------------------------------------------------------------------      
      RETURN
      END

!----------------------------------------------------------------------!
!                                                                      !
!     Subroutines for the Orbital Optimization using OrbOptADAM        !

!                                                                      !
!----------------------------------------------------------------------!

! RotOrb
      SUBROUTINE RotOrb(NV,Y,COEF)
!-----------------------------------------------------------------------
!     This subroutine generate an unitary matrix U = exp(YM), with YM
!     an antisymmetric matrix [YM(I,J)=-YM(J,I)]. The variable Y store
!     the upper triangular part of YM. The parameters of Y can be varied
!     in a minimization procedure to optimize the natural orbitals via
!     the rotation COEFn = COEF * U.
!-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      DOUBLE PRECISION,DIMENSION(NV) :: Y
      DOUBLE PRECISION,DIMENSION(NBF,NBF) :: COEF
!
      INTEGER,ALLOCATABLE,DIMENSION(:) :: IPIV
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: RWORK
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: U
      COMPLEX,ALLOCATABLE,DIMENSION(:) :: W,WORK
      COMPLEX,ALLOCATABLE,DIMENSION(:,:) :: YM,VL,VR,VRinv

      DOUBLE COMPLEX,ALLOCATABLE,DIMENSION(:,:) :: MAT,MAT2
      DOUBLE COMPLEX :: ALPHA,BETA
      INTEGER :: INFO
!-----------------------------------------------------------------------
!     Use Y to fill YM
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(YM(NBF,NBF))
      YM = (0.0D0, 0.0D0)
      N = 1
      DO I=1,NBF-1
       DO J=I+1,NBF
        YM(I,J) =  Y(N)
        YM(J,I) = -Y(N)
        N = N + 1
       END DO
      END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Diagonalize antisymmetric matrix YM
!     Find eigenvectors and eigenvalues of YM
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(W(NBF),VL(1,NBF),VR(NBF,NBF),VRinv(NBF,NBF))
      ALLOCATE(WORK(2*NBF),RWORK(2*NBF),IPIV(NBF),U(NBF,NBF))
      ALLOCATE(MAT(NBF,NBF),MAT2(NBF,NBF))
      CALL ZGEEV("N","V",NBF,YM,NBF,W,VL,1,VR,NBF,WORK,2*NBF,RWORK,INFO)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Invert eigenvectors VR
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      VRinv = VR
      CALL ZGETRF(NBF,NBF,VRinv,NBF,IPIV,INFO)
      CALL ZGETRI(NBF,VRinv,NBF,IPIV,WORK,NBF,INFO)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Take exponential of eigenvalues and do MAT = M*exp(W)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      W = EXP(W)
      DO I=1,NBF
       DO J=1,NBF
        MAT(I,J) = VR(I,J)*W(J)
       END DO
      END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Build the unitary transformation U = exp(Y) = M*exp(W)*M^-1
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      MAT2 = (0.0D0, 0.0D0)
!      ALPHA = (1.0D0, 0.0D0)
!      BETA = (0.0D0, 0.0D0)
!      CALL ZGEMM('N', 'C', NBF, NBF, NBF, ALPHA, MAT, NBF, VR, NBF,     &
!              BETA, MAT2, NBF)

      MAT2 = MATMUL(MAT,VRinv)

      U = REAL(MAT2)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Rotate the orbials COEFn = COEF * U
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      COEF = MATMUL(COEF,U)
!-----------------------------------------------------------------------
      DEALLOCATE(YM,W,VL,VR,VRinv,WORK,RWORK,IPIV,U,MAT,MAT2)
      RETURN
      END

!----------------------------------------------------------------------!
