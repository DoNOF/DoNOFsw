!----------------------------------------------------------------------!
!                                                                      !
! O T H E R   O R B I T A L   O P T I M I Z A T I O N   M E T H O D S  !
!                                                                      !
!----------------------------------------------------------------------!
!                                                                      !
! IORBOPT = 3  ADABelief: Adaptive with Belief in gradients            !
! IORBOPT = 4  YOGI                                                    !
! IORBOPT = 5  DEMON: Decaying Momentum                                !
! IORBOPT = 6  OrbOptSQP: Use a Sequential Quadratic Program (SQP)     !
!                                                                      !
!----------------------------------------------------------------------!

      SUBROUTINE OrbOpt36(ITCALL,OVERLAP,AHCORE,IJKL,XIJKL,XIJKaux,     &
                          QD,COEF,RO,CJ12,CK12,ELAG,DIPN,ILOOP,IORBOPT, &
                          OCCTIME,IPRINTOPT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPFILE_NIJKL/NINTMX,NIJKL,NINTCR,NSTORE
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      LOGICAL ERIACTIVATED
      COMMON/ERIACT/ERIACTIVATED,NIJKaux,NINTCRaux,NSTOREaux,IAUXDIM
!
      INTEGER :: ITCALL,ILOOP,IPRINTOPT,IORBOPT
      INTEGER(8),DIMENSION(NSTORE)::IJKL
      DOUBLE PRECISION,DIMENSION(NSTORE)::XIJKL
      DOUBLE PRECISION,DIMENSION(NSTOREaux)::XIJKaux
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::OVERLAP,AHCORE,COEF,ELAG
      DOUBLE PRECISION,DIMENSION(NBF,NBF,NBF)::QD
      DOUBLE PRECISION,DIMENSION(NBF5)::RO
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5)::CJ12,CK12
      DOUBLE PRECISION,DIMENSION(3)::DIPN
!-----------------------------------------------------------------------
      IF(IORBOPT==3)THEN
       CALL OrbOptADABelief(ITCALL,OVERLAP,AHCORE,IJKL,XIJKL,XIJKaux,   &
                            QD,COEF,RO,CJ12,CK12,ELAG,DIPN,ILOOP,       &
                            OCCTIME,IPRINTOPT)
      ELSE IF(IORBOPT==4)THEN
       CALL OrbOptYOGI(ITCALL,OVERLAP,AHCORE,IJKL,XIJKL,XIJKaux,QD,COEF,&
                       RO,CJ12,CK12,ELAG,DIPN,ILOOP,OCCTIME,IPRINTOPT)
      ELSE IF(IORBOPT==5)THEN
       CALL OrbOptDEMON(ITCALL,OVERLAP,AHCORE,IJKL,XIJKL,XIJKaux,QD,    &
                        COEF,RO,CJ12,CK12,ELAG,DIPN,ILOOP,IPRINTOPT)
      ELSE IF(IORBOPT==6)THEN
       CALL OrbOptSQP(ITCALL,OVERLAP,AHCORE,IJKL,XIJKL,XIJKaux,         &
                      QD,COEF,RO,CJ12,CK12,ELAG,DIPN,ILOOP,IPRINTOPT)
      END IF
!-----------------------------------------------------------------------
      RETURN
      END

! OrbOptADABelief
      SUBROUTINE OrbOptADABelief(ITCALL,OVERLAP,AHCORE,IJKL,XIJKL,      &
                                 XIJKaux,QD,COEF,RO,CJ12,CK12,ELAG,     &
                                 DIPN,ILOOP,OCCTIME,IPRINTOPT)
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
!
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
!     First Call to OrbOptADABelief
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(ITCALL==1)THEN
       IF(IPRINTOPT==1)WRITE(6,1)
       EELEC_OLD = EELEC
      ENDIF
!-----------------------------------------------------------------------
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
       IF(IPRINTOPT==1)WRITE(6,2)ITCALL,OCCTIME,ORBTIME,ILOOP,EELEC,   &
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
      !DO WHILE( ILOOP<MAXLOOP .and. DABS(DELE)>THRESHEC )
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
          V(J) = BETA2 * V(J) + (1.0D0-BETA2) * (GRAD(J) - M(J))**2
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
!AO      CALL ELAGaor(OVERLAP,AHCORE,IJKL,XIJKL,COEF,RO,DEN,ELAG)
         CALL ELAGr(AHCORE,IJKL,XIJKL,COEF,RO,DEN,ELAG,G,IPNOF)
        ELSE
         CALL ELAG1r(AHCORE,IJKL,XIJKL,XIJKaux,QD,COEF,RO,CJ12,CK12,ELAG,G)
        ENDIF
        CALL PCONVE(ELAG,DUMEL,MAXI,MAXJ,SUMDIF)        ! DUMEL
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Calculate Electronic Energy (EELEC)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        EELEC0 = EELEC
        EELEC = TRACE(DEN,AHCORE,NBF)
        DO iq=1,NBF5
         EELEC = EELEC + ELAG(iq,iq)
        END DO
!       Include Nuclear Dipoles if electric field =/ 0
        IF(EFIELDL)EELEC = EELEC - EFX*DIPN(1)-EFY*DIPN(2)-EFZ*DIPN(3)
        DELE = EELEC - EELEC0
        Etot = EELEC + EN
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Store best energy and its orbitals
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
      IF(IPRINTOPT==1)WRITE(6,2)ITCALL,OCCTIME,ORBTIME,ILOOP,           &
                                EELEC,Etot,DIF_EELEC,DUMEL
      IF(PCONV<THRESHE .and. DUMEL<THRESHL) THEN
       CONVGDELAG=.TRUE.
       AUTOLR=.TRUE.
       LR=0.01D0
       MAXLOOP=10
      END IF
!-----------------------------------------------------------------------
!     FORMAT STATEMENTS
!-----------------------------------------------------------------------
    1 FORMAT(//2X,'RESULTS OF OCCUPATION-COEFFICIENT OPTIMIZATION'      &
            ,/1X,'================================================',    &
        //2X,'Ex It',2X'Occ Time',2X'Orb Time',2X,'Orb It',3X,          &
        'Electronic Energy',5X,'Total Energy',3X,'Energy Convergency',  &
        4X,'Max Mul-Lag Diff',/)
    2 FORMAT(I5,3X,ES8.1,2X,ES8.1,3X,I5,1X,F18.10,1X,F19.10,2X,F15.10,  &
            8X,F11.6)
    3 FORMAT(2X,I3,'.',3X,F17.8,4X,F15.8,6X,F11.6)
!-----------------------------------------------------------------------
      RETURN
      END

! OrbOptYOGI
      SUBROUTINE OrbOptYOGI(ITCALL,OVERLAP,AHCORE,IJKL,XIJKL,XIJKaux,   &
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
!
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
!     First Call to OrbOptYOGI
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(ITCALL==1)THEN
       IF(IPRINTOPT==1)WRITE(6,1)
       EELEC_OLD = EELEC
      ENDIF
!-----------------------------------------------------------------------
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
      IF( (IPNOF==5.or.IPNOF==7.or.(IPNOF==8.and.Imod==0)) .and.       &
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
       IF(IPRINTOPT==1)WRITE(6,2)ITCALL,OCCTIME,ORBTIME,ILOOP,EELEC,   &
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
      !DO WHILE( ILOOP<MAXLOOP .and. DABS(DELE)>THRESHEC )
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
          V(J) = V(J) - (1.0D0-BETA2)*SIGN(GRAD(J)**2,V(J)-GRAD(J)**2)
          MHAT(J) = M(J) / (1.0D0-BETA1**ILOOP)
          VHAT(J) = V(J) / (1.0D0-BETA2**ILOOP)
          VHATMAX(J) = MAX(VHATMAX(J), VHAT(J))
          Y(J) = -LR*MHAT(J)/(SQRT(VHATMAX(J)) + 10D-16)
        END DO
        CALL RotOrb(NV,Y,COEF)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       Compute Gradient & Lagrangian Multipliers (ELAG)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        CALL DENMATr(DEN,COEF,RO,NBF,1,NBF5) ! Density Matrix
        IF( (IPNOF==5.or.IPNOF==7.or.(IPNOF==8.and.Imod==0)) .and.      &
              MSpin==0 .and. IERITYP==1 )THEN
!AO      CALL ELAGaor(OVERLAP,AHCORE,IJKL,XIJKL,COEF,RO,DEN,ELAG)
         CALL ELAGr(AHCORE,IJKL,XIJKL,COEF,RO,DEN,ELAG,G,IPNOF)
        ELSE
         CALL ELAG1r(AHCORE,IJKL,XIJKL,XIJKaux,QD,COEF,RO,CJ12,CK12,    &
                     ELAG,G)
        ENDIF
        CALL PCONVE(ELAG,DUMEL,MAXI,MAXJ,SUMDIF)        ! DUMEL
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       Calculate Electronic Energy (EELEC)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        EELEC0 = EELEC
        EELEC = TRACE(DEN,AHCORE,NBF)
        DO iq=1,NBF5
         EELEC = EELEC + ELAG(iq,iq)
        END DO
!       Include Nuclear Dipoles if electric field =/ 0
        IF(EFIELDL)EELEC = EELEC - EFX*DIPN(1)-EFY*DIPN(2)-EFZ*DIPN(3)
        DELE = EELEC - EELEC0
        Etot = EELEC + EN
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       Store best energy and its orbitals
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        IF(EELEC < BEST_E) THEN
          IMPROVED = .TRUE.
          BEST_E = EELEC
          BEST_COEF = COEF
        END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       Intermediate Output of the internal iteration (Nprint=2)
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
       CALL ELAG1r(AHCORE,IJKL,XIJKL,XIJKaux,QD,COEF,RO,CJ12,CK12,      &
                   ELAG,G)
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
      IF(IPRINTOPT==1)WRITE(6,2)ITCALL,OCCTIME,ORBTIME,ILOOP,EELEC,Etot,DIF_EELEC,DUMEL
      IF(PCONV<THRESHE .and. DUMEL<THRESHL) THEN
       CONVGDELAG=.TRUE.
       AUTOLR=.TRUE.
       LR=0.01D0
       MAXLOOP=10
      END IF
!-----------------------------------------------------------------------
!     FORMAT STATEMENTS
!-----------------------------------------------------------------------
    1 FORMAT(//2X,'RESULTS OF OCCUPATION-COEFFICIENT OPTIMIZATION'      &
            ,/1X,'================================================',    &
        //2X,'Ex It',2X'Occ Time',2X'Orb Time',2X,'Orb It',3X,          &
        'Electronic Energy',5X,'Total Energy',3X,'Energy Convergency',  &
        4X,'Max Mul-Lag Diff',/)
    2 FORMAT(I5,3X,ES8.1,2X,ES8.1,3X,I5,1X,F18.10,1X,F19.10,2X,F15.10,  &
            8X,F11.6)
    3 FORMAT(2X,I3,'.',3X,F17.8,4X,F15.8,6X,F11.6)
!-----------------------------------------------------------------------
      RETURN
      END

! OrbOptDEMON
      SUBROUTINE OrbOptDEMON(ITCALL,OVERLAP,AHCORE,IJKL,XIJKL,XIJKaux,  &
                             QD,COEF,RO,CJ12,CK12,ELAG,DIPN,ILOOP,      &
                             IPRINTOPT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL CONVGDELAG,SMCD,ERIACTIVATED
      COMMON/CONVERGENCE/DUMEL,PCONV,CONVGDELAG
      COMMON/ERITYPE/IERITYP,IRITYP,IGEN,ISTAR,MIXSTATE,SMCD
      COMMON/ERIACT/ERIACTIVATED,NIJKaux,NINTCRaux,NSTOREaux,IAUXDIM
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPFILE_NIJKL/NINTMX,NIJKL,NINTCR,NSTORE
      COMMON/PUNTEROSROT/M1,M2,M3,M4,M5,M6,M7,M8,M9,M10,M11,MUSE
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      COMMON/EHFEN/EHF,EN
      COMMON/INPNOF_THRESH/THRESHL,THRESHE,THRESHEC,THRESHEN
      COMMON/INPNOF_ORBSPACE0/NO1,NDOC,NCO,NCWO,NVIR,NAC,NO0
      COMMON/ENERGIAS/EELEC,EELEC_OLD,DIF_EELEC,EELEC_MIN
      COMMON/INPNOF_COEFOPT/MAXLOOP
!
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
      INTEGER(8),ALLOCATABLE,DIMENSION(:) :: IUSER
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: USER,Y,GRAD
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: M,V,VHAT,VHATMAX
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: BEST_COEF
      LOGICAL :: IMPROVED
      DOUBLE PRECISION :: LR, BETA1, BETA2, RATE, BINIT
!-----------------------------------------------------------------------
      ILOOP = ILOOP
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     First Call to OrbOptDEMON
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(ITCALL==1)THEN
       IF(IPRINTOPT==1)WRITE(6,1)
       EELEC_OLD = EELEC
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     IUSER
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(IUSER(NSTORE))
      CALL IXtoIX0(IJKL,IUSER,NSTORE)  ! IUSER = IJKL(NSTORE)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     USER
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      M1  = 1                    ! USER(M1) = COEF(NBF,NBF)
      M2  = M1  +  NSQ           ! USER(M2) = AHCORE(NBF,NBF)
      M3  = M2  +  NSQ           ! USER(M3) = XIJKL(NSTORE)
      M4  = M3  +  NSTORE        ! USER(M4) = QD(NBF,NBF,NBF)
      M5  = M4  +  NBF*NSQ       ! USER(M5) = RO(NBF5)
      M6  = M5  +  NBF5          ! USER(M6) = CJ12(NBF5,NBF5)
      M7  = M6  +  NSQ5          ! USER(M7) = CK12(NBF5,NBF5)
      M8  = M7  +  NSQ5          ! USER(M8) = ELAG(NBF,NBF)
      M9  = M8  +  NSQ           ! USER(M9) = DIPN(3)
      M10 = M9  +  3             ! USER(M10)= OVERLAP(NBF,NBF)
      if(IERITYP==1)then
       MUSE = M10 - M1 + NSQ
      else if(IERITYP==2 .or. IERITYP==3)then
       M11 = M10 +  NSQ          ! USER(M11)= XIJKaux(NSTOREaux)
       MUSE = M11 - M1 + NSTOREaux
      end if
      ALLOCATE (USER(MUSE))
!
      CALL XtoX0(   COEF,USER(M1),NSQ)
      CALL XtoX0( AHCORE,USER(M2),NSQ)
      CALL XtoX0(  XIJKL,USER(M3),NSTORE)
      CALL XtoX0(     QD,USER(M4),NBF*NSQ)
      CALL XtoX0(     RO,USER(M5),NBF5)
      CALL XtoX0(   CJ12,USER(M6),NSQ5)
      CALL XtoX0(   CK12,USER(M7),NSQ5)
      CALL XtoX0(   ELAG,USER(M8),NSQ)
      CALL XtoX0(   DIPN,USER(M9),3)
      CALL XtoX0(OVERLAP,USER(M10),NSQ)
      if(IERITYP>1)CALL XtoX0(XIJKaux,USER(M11),NSTOREaux)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Initialize variables for Y
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      NV = NBF*(NBF-1)/2
      ALLOCATE(M(NV), V(NV))
      ALLOCATE(Y(NV), GRAD(NV), VHAT(NV), VHATMAX(NV))
      Y = 0.0D0
      GRAD = 0.0D0
      M = 0.0D0
      V = 0.0D0
      VHAT = 0.0D0
      VHATMAX = 0.0D0

      LR = 0.02 !0.001
      BETA1 = 0.7 !0.9
      BETA2 = 0.9

      BEST_E = EELEC

      IMPROVED = .FALSE.
      ALLOCATE(BEST_COEF(NBF,NBF))
      BEST_COEF = COEF

      BINIT = BETA1

      NIT = 0
      DO I=1, MAXLOOP
        NIT = NIT + 1
        CALL CALCORBG(NV,Y,NF,GRAD,IUSER,USER)
        IF(NORM2(GRAD) < 10E-3 .AND. IMPROVED) EXIT

        RATE = 1.0D0*(I-1)/(MAXLOOP-1)
        BETA1 = BINIT*(1.0D0-RATE)/((1.0D0-BINIT)+BINIT*(1.0D0-RATE))
        DO J=1,NV
          M(J) = BETA1 * M(J) + GRAD(J)
          V(J) = BETA2 * V(J) + (1.0D0-BETA2) * GRAD(J)**2
          Y(J) = -LR*M(J)/(SQRT(V(J) + 10D-8))
        END DO

        CALL RotOrb(NV,Y,COEF)
        CALL XtoX0(COEF,USER(M1),NSQ)
        Y = 0.0D0
        CALL CALCORBE(NV,Y,NF,EELEC,IUSER,USER)
        IF(EELEC < BEST_E) THEN
          IMPROVED = .TRUE.
          BEST_E = EELEC
          BEST_COEF = COEF
        END IF

      END DO
      COEF = BEST_COEF
      EELEC = BEST_E
      DEALLOCATE(BEST_COEF)
      DEALLOCATE(GRAD, M, V, VHAT, VHATMAX)

      IF (.NOT. IMPROVED) THEN
       LR = LR/2
       MAXLOOP = MIN(300,MAXLOOP+30)
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Update Orbitals
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     CALL RotOrb(NV,Y,COEF)
      CALL XtoX0(COEF,USER(M1),NSQ)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Calculate Electronic Energy (EELEC) for new Orbitals
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Y = 0.0D0
      CALL CALCORBE(NV,Y,NF,EELEC,IUSER,USER)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Update QD and ELAG
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL XtoX0(USER(M4),QD,NBF*NSQ)
      CALL XtoX0(USER(M8),ELAG,NSQ)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Check convergence for the Energy
!     Check for the Symmetry of Lagrangian (ELAG)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Etot = EELEC + EN
      DIF_EELEC = EELEC - EELEC_OLD
      EELEC_OLD = EELEC
      PCONV = ABS(DIF_EELEC)
      CALL PCONVE(ELAG,DUMEL,MAXI,MAXJ,SUMDIF)
      IF(IPRINTOPT==1)WRITE(6,2)ITCALL,NIT,EELEC,Etot,DIF_EELEC,DUMEL
      IF(PCONV<THRESHE .and. DUMEL<THRESHL)CONVGDELAG=.TRUE.
!-----------------------------------------------------------------------
!     FORMAT STATEMENTS
!-----------------------------------------------------------------------
    1 FORMAT(//2X,'RESULTS OF OCCUPATION-COEFFICIENT OPTIMIZATION'      &
            ,/1X,'================================================',    &
        //2X,'Ex It',2X,'Orb It',4X,'Electronic Energy',5X,           &
        'Total Energy',3X,'Energy Convergency',4X,'Max Mul-Lag Diff',/)
    2 FORMAT(I5,3X,I5,1X,F20.10,1X,F19.10,2X,F15.10,8X,F11.6)
!-----------------------------------------------------------------------
      DEALLOCATE (IUSER,USER,Y)
      RETURN
      END

! CALCORBE
      SUBROUTINE CALCORBE(NV,Y,NF,ENERGYEL,IUSER,USER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL EFIELDL,SMCD,HighSpin
      COMMON/INPFILE_NIJKL/NINTMX,NIJKL,NINTCR,NSTORE
      COMMON/PUNTEROSROT/M1,M2,M3,M4,M5,M6,M7,M8,M9,M10,M11,MUSE
      COMMON/INPNOF_PNOF/IPNOF,NTWOPAR
      COMMON/INP_EFIELDL/EFX,EFY,EFZ,EFIELDL
      COMMON/ENERGIAS/EELEC,EELEC_OLD,DIF_EELEC,EELEC_MIN
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
      COMMON/ERITYPE/IERITYP,IRITYP,IGEN,ISTAR,MIXSTATE,SMCD
      COMMON/INPNOF_MOD/Imod
!
      INTEGER(8),DIMENSION(NSTORE) :: IUSER
      DOUBLE PRECISION,DIMENSION(NV) :: Y
      DOUBLE PRECISION,DIMENSION(MUSE) :: USER
!
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: COEF,DEN,G
!-----------------------------------------------------------------------
!     Avoiding warnings
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      NF = NF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Copy USER(M1) to COEF (Orbitals without rotation)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(COEF(NBF,NBF))
      CALL XtoX0(USER(M1),COEF,NSQ)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Rotate Orbitals by COEF*EXP(YM) -> COEF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL RotOrb(NV,Y,COEF)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Compute Lagrangian Multipliers (ELAG)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(DEN(NBF,NBF),G(NBF,NBF5))
      CALL DENMATr(DEN,COEF,USER(M5),NBF,1,NBF5) ! Density Matrix
      IF( (IPNOF==5.or.IPNOF==7.or.(IPNOF==8.and.Imod==0)) .and.        &
              MSpin==0 .and. IERITYP==1 )THEN
!AO    CALL ELAGaor(USER(M10),USER(M2),IUSER,USER(M3),COEF,             &
!AO                 USER(M5),DEN,USER(M8))
       CALL ELAGr(USER(M2),IUSER,USER(M3),COEF,USER(M5),DEN,USER(M8),   &
                  G,IPNOF)
      ELSE
       CALL ELAG1r(USER(M2),IUSER,USER(M3),USER(M11),USER(M4),          &
                   COEF,USER(M5),USER(M6),USER(M7),USER(M8),G)
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Calculate Electronic Energy (ENERGYEL)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ENERGYEL = TRACE(DEN,USER(M2),NBF)
      DO iq=1,NBF5
       iqiq = M8-1 + iq+(iq-1)*nbf  ! ELAG(iq,iq)
       ENERGYEL = ENERGYEL + USER(iqiq)
      END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Include Nuclear Dipoles if electric field =/ 0
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(EFIELDL)THEN
       ENERGYEL = ENERGYEL- EFX*USER(M9)-EFY*USER(M9+1)-EFZ*USER(M9+2)
      END IF
!-----------------------------------------------------------------------
      EELEC = ENERGYEL
      DEALLOCATE(COEF,DEN,G)
      RETURN
      END

! CALCORBG
      SUBROUTINE CALCORBG(NV,Y,NF,GRAD,IUSER,USER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/INPFILE_NIJKL/NINTMX,NIJKL,NINTCR,NSTORE
      COMMON/PUNTEROSROT/M1,M2,M3,M4,M5,M6,M7,M8,M9,M10,M11,MUSE
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
!
      INTEGER(8),DIMENSION(NSTORE) :: IUSER
      DOUBLE PRECISION,DIMENSION(NV) :: Y,GRAD
      DOUBLE PRECISION,DIMENSION(MUSE) :: USER
!
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: ELAG
!-----------------------------------------------------------------------
!     Avoiding warnings
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      NF = NF
      IUSER(1) = IUSER(1)
      Y(1) = Y(1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(ELAG(NBF,NBF))
      CALL XtoX0(USER(M8),ELAG,NSQ)
!
      GRAD = 0.0D0
      K = 1
      DO I=1,NBF-1
       DO J=I+1,NBF
        GRAD(K) = 2.0d0 * ( ELAG(I,J) - ELAG(J,I) )
        K = K + 1
       END DO
      END DO
!-----------------------------------------------------------------------
      DEALLOCATE(ELAG)
      RETURN
      END

!----------------------------------------------------------------------!
!                                                                      !
!     O R B I T A L     O P T I M I Z A T I O N     B Y    S Q P       !
!                                                                      !
!   OrbOptSQP: Minimize w.r.t. MOs using a Sequential Quadratic Method !
!   CoefEnergy: Calculate Electronic Energy                            !
!   CoefConstr: Calculate Constraints w.r.t. COEF (CtSC=1)             !
!                                                                      !
!----------------------------------------------------------------------!

! OrbOptSQP
      SUBROUTINE OrbOptSQP(ITCALL,OVERLAP,AHCORE,IJKL,XIJKL,XIJKaux,QD, &
                           COEF,RO,CJ12,CK12,ELAG,DIPN,ILOOP,IPRINTOPT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL CONVGDELAG,SMCD,ERIACTIVATED
      COMMON/CONVERGENCE/DUMEL,PCONV,CONVGDELAG
      COMMON/ERITYPE/IERITYP,IRITYP,IGEN,ISTAR,MIXSTATE,SMCD
      COMMON/ERIACT/ERIACTIVATED,NIJKaux,NINTCRaux,NSTOREaux,IAUXDIM
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPFILE_NIJKL/NINTMX,NIJKL,NINTCR,NSTORE
      COMMON/PUNTEROSSQP/M1,M2,M3,M4,M5,M6,M7,M8,M9,M10,MUSE
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      COMMON/EHFEN/EHF,EN
      COMMON/INPNOF_THRESH/THRESHL,THRESHE,THRESHEC,THRESHEN
      COMMON/ENERGIAS/EELEC,EELEC_OLD,DIF_EELEC,EELEC_MIN
!
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
      INTEGER(8),ALLOCATABLE,DIMENSION(:) :: IUSER
      INTEGER,ALLOCATABLE,DIMENSION(:) :: IWORK,ISTATE
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)   :: USER,GRAD,WORK
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)   :: BL,BU,CONST,CLAMDA
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: R,AAA,CJAC(:,:)
!
!      EXTERNAL CoefConstr,CoefEnergy,E04UCF,E04UEF !nag
!-----------------------------------------------------------------------
!     IUSER
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(IUSER(NSTORE))
      CALL IXtoIX0(IJKL,IUSER,NSTORE)  ! IUSER = IJKL(NSTORE)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     USER
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      M1  = 1                    ! USER(M1) = AHCORE(NBF,NBF)
      M2  = M1  + NSQ            ! USER(M2) = XIJKL(NSTORE)
      M3  = M2  + NSTORE         ! USER(M3) = QD(NBF,NBF,NBF)
      M4  = M3  + NBF*NSQ        ! USER(M4) = RO(NBF5)
      M5  = M4  + NBF5           ! USER(M5) = CJ12(NBF5,NBF5)
      M6  = M5  + NSQ5           ! USER(M6) = CK12(NBF5,NBF5)
      M7  = M6  + NSQ5           ! USER(M7) = ELAG(NBF,NBF)
      M8  = M7  + NSQ            ! USER(M8) = DIPN(3)
      M9  = M8  +  3             ! USER(M9) = OVERLAP(NBF,NBF)
      if(IERITYP==1)then
       MUSE = M9 - M1 + NSQ
      else if(IERITYP==2 .or. IERITYP==3)then
       M10 = M9 +  NSQ           ! USER(M10)= XIJKaux(NSTOREaux)
       MUSE = M10 - M1 + NSTOREaux
      end if
      ALLOCATE (USER(MUSE))
!
      CALL XtoX0( AHCORE,USER(M1),NSQ)
      CALL XtoX0(  XIJKL,USER(M2),NSTORE)
      CALL XtoX0(     QD,USER(M3),NBF*NSQ)
      CALL XtoX0(     RO,USER(M4),NBF5)
      CALL XtoX0(   CJ12,USER(M5),NSQ5)
      CALL XtoX0(   CK12,USER(M6),NSQ5)
      CALL XtoX0(   ELAG,USER(M7),NSQ)
      CALL XtoX0(   DIPN,USER(M8),3)
      CALL XtoX0(OVERLAP,USER(M9),NSQ)
      if(IERITYP>1)CALL XtoX0(XIJKaux,USER(M10),NSTOREaux)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     First Call to OrbOptSQP
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(ITCALL==1)THEN
       IF(IPRINTOPT==1)WRITE(6,1)
       EELEC_OLD = EELEC
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Calculate Initial Electronic Energy (EELEC)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE (GRAD(NSQ))
      CALL CoefEnergy(0,NSQ,COEF,EELEC,GRAD,1,IUSER,USER)
!-.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.-
!                       Energy Minimization
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      NV = NSQ
      NRES = NBFT
!- - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Work Dimensions
!- - - - - - - - - - - - - - - - - - - - - - - - - - - -
      LIWORK = 3*NV + 2*NRES
      LWORK  = 2*NV*NV + 2*NV*NRES + 20*NV + 21*NRES
      IFAIL  = -1
!
      ALLOCATE (IWORK(3*NV+2*NRES),ISTATE(NV+NRES))
      ALLOCATE (R(NV,NV),AAA(1,1))
      ALLOCATE (BL(NV+NRES),BU(NV+NRES),CONST(NRES))
      ALLOCATE (CJAC(NRES,NV),CLAMDA(NV+NRES))
      ALLOCATE (WORK(2*NV*NV+2*NV*NRES+20*NV+21*NRES))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Lower and Upper bounds for COEF
!- - - - - - - - - - - - - - - - - - -
      K = 0
      DO IQ=1,NBF
       do i=1,nbf
        K = K+1
        BL(K) = COEF(i,IQ) - 1.0d0
        BU(K) = COEF(i,IQ) + 1.0d0
       enddo
      ENDDO
!- - - - - - - - - - - - - - - - - - -
!     Constraints as equalities
!- - - - - - - - - - - - - - - - - - -
      DO I=NV+1,NV+NRES
       BL(I) = 0.0d0
       BU(I) = 0.0d0
      ENDDO
      DO I=1,NBF
       II = I*(I+1)/2
       BL(NV+II) = 1.0d0
       BU(NV+II) = 1.0d0
      ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Send output of E04DGF to CGM file
!- - - - - - - - - - - - - - - - - - - - - - - - - - -
!      CALL X04ABF(1,2) !nag
!- - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Function Precision
!- - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(THRESHEC >= 1.0d-8)THEN
!       CALL E04UEF('Function Precision = 1.0d-8') !nag
      ELSEIF(THRESHEC == 1.0d-9)THEN
!       CALL E04UEF('Function Precision = 1.0d-9') !nag
      ELSEIF(THRESHEC <= 1.0d-10)THEN
!       CALL E04UEF('Function Precision = 1.0d-10') !nag
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Major Print Level (0 = No output)
!- - - - - - - - - - - - - - - - - - - - - - - - - - -
!      CALL E04UEF('Major Print Level = 0') !nag
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                  =-1 no verifica nada   , 0 Cheap Test (default)
!     Verify Level = 1 Objective Gradients, 2 Constraint Gradients
!                  = 3 Verify all Gradients
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      CALL E04UEF('Verify Level = -1') !nag
!-----------------------------------------------------------------------
!     MINIMIZATION PROCESS
!-----------------------------------------------------------------------
!      CALL E04UCF(NV,0,NRES,1,NRES,NV,AAA,BL,BU,CoefConstr,CoefEnergy,  & !nag
!                  ILOOP,ISTATE,CONST,CJAC,CLAMDA,EELEC,GRAD,R,COEF,     &
!                  IWORK,LIWORK,WORK,LWORK,IUSER,USER,IFAIL)
      DIF_EELEC = EELEC - EELEC_OLD
      EELEC_OLD = EELEC
      PCONV = ABS(DIF_EELEC)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Check for the Symmetry of Lagrangian (ELAG(IQJQ)-ELAG(JQIQ))
!     for ITCALL>1 (do not stop in the first call)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL XtoX0(USER(M7),ELAG,NSQ)
      CALL PCONVE(ELAG,DUMEL,MAXI,MAXJ,SUMDIF)
      IF(IPRINTOPT==1)WRITE(6,2)ITCALL,EELEC,EELEC+EN,DIF_EELEC,DUMEL
      IF(DUMEL<THRESHL .and. PCONV<THRESHE)CONVGDELAG=.TRUE.
!-----------------------------------------------------------------------
!     FORMAT STATEMENTS
!-----------------------------------------------------------------------
    1 FORMAT(//2X,'RESULTS OF OCCUPATION-COEFFICIENT SQP PROCEDURE'     &
            ,/1X,'=================================================',   &
        //2X,'Iter',5X,'Electronic Energy',6X,'Total Energy',           &
          3X,'Energy Convergency',4X,'Max Mul-Lag Diff',/)
    2 FORMAT(I5,'.',1X,F20.10,1X,F19.10,2X,F15.10,8X,F11.6)
!-----------------------------------------------------------------------
      DEALLOCATE (IUSER,IWORK,ISTATE)
      DEALLOCATE (GRAD,R,USER,AAA,BL,BU,CONST,CJAC,CLAMDA,WORK)
      RETURN
      END

! CoefEnergy
      SUBROUTINE CoefEnergy(MODE,NV,COEF,ENERGYEL,GG,NSTATE,IUSER,USER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL EFIELDL,SMCD,HighSpin
      COMMON/INPNOF_PNOF/IPNOF,NTWOPAR
      COMMON/INP_EFIELDL/EFX,EFY,EFZ,EFIELDL
      COMMON/ERITYPE/IERITYP,IRITYP,IGEN,ISTAR,MIXSTATE,SMCD
      COMMON/INPFILE_NIJKL/NINTMX,NIJKL,NINTCR,NSTORE
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/PUNTEROSSQP/M1,M2,M3,M4,M5,M6,M7,M8,M9,M10,MUSE
      COMMON/INPNOF_MOD/Imod
      INTEGER :: NSOC,NDNS,MSpin
      INTEGER(8),DIMENSION(NSTORE) :: IUSER
      DOUBLE PRECISION,DIMENSION(MUSE) :: USER
      DOUBLE PRECISION,DIMENSION(NV) :: GG
      DOUBLE PRECISION,DIMENSION(NBF,NBF) :: COEF
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: DEN,G
!-----------------------------------------------------------------------
      NSTATE = NSTATE
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Compute Gradient & Lagrangian Multipliers (ELAG)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(DEN(NBF,NBF),G(NBF,NBF5))
      CALL DENMATr(DEN,COEF,USER(M4),NBF,1,NBF5) ! Density Matrix
      IF( (IPNOF==5.or.IPNOF==7.or.(IPNOF==8.and.Imod==0)) .and.        &
              MSpin==0 .and. IERITYP==1 )THEN
       CALL ELAGr(USER(M1),IUSER,USER(M2),COEF,USER(M4),DEN,USER(M7),   &
                  G,IPNOF)
      ELSE
       CALL ELAG1r(USER(M1),IUSER,USER(M2),USER(M10),USER(M3),COEF,     &
                   USER(M4),USER(M5),USER(M6),USER(M7),G)
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Calculate Electronic Energy (ENERGYEL)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ENERGYEL = TRACE(DEN,USER(M1),NBF)
      DO iq=1,NBF5
       iqiq = M7-1 + iq+(iq-1)*nbf  ! ELAG(iq,iq)
       ENERGYEL = ENERGYEL + USER(iqiq)
      END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Include Nuclear Dipoles if electric field =/ 0
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(EFIELDL)THEN
       ENERGYEL = ENERGYEL- EFX*USER(M8)-EFY*USER(M8+1)-EFZ*USER(M8+2)
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Calculate GRADIENTS with respect to the Coefficients (COEF)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(MODE==1.or.MODE==2)THEN
       GG = 0.0d0
       DO IQ=1,NBF5
        do i=1,nbf
         iIQ = i+(IQ-1)*nbf
!        Use of real coefficients
         GG(iIQ)= 4.0 * G(i,IQ)
        enddo
       ENDDO
      ENDIF
!-----------------------------------------------------------------------
      DEALLOCATE(DEN,G)
      RETURN
      END

! CoefConstr
      SUBROUTINE CoefConstr(MODE,NRES,NV,LDCJ,NEEDC,COEF,CONST,CJAC,    &
                            NSTATE,IUSER,USER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/INPFILE_NIJKL/NINTMX,NIJKL,NINTCR,NSTORE
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/PUNTEROSSQP/M1,M2,M3,M4,M5,M6,M7,M8,M9,M10,MUSE
      INTEGER,DIMENSION(NRES) :: NEEDC
      INTEGER,DIMENSION(NSTORE) :: IUSER
      DOUBLE PRECISION,DIMENSION(MUSE) :: USER
      DOUBLE PRECISION,DIMENSION(NBF,NBF) :: COEF
      DOUBLE PRECISION,DIMENSION(NRES) :: CONST
      DOUBLE PRECISION,DIMENSION(LDCJ,NV) :: CJAC
!-----------------------------------------------------------------------
      IUSER(1)=IUSER(1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     First call to CONFUN.  Set all Jacobian elements to zero.
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(NSTATE==1)CJAC=0.0d0
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(MODE==0.or.MODE==2)THEN
       K = 0
       DO IQ=1,NBF
        DO JQ=1,IQ
         K = K+1
         CONST(K) = 0.0d0
         do i=1,nbf
          CONST(K) = CONST(K) + COEF(i,IQ)*FC(i,JQ,USER(M9),COEF)
         enddo
        ENDDO
       ENDDO
      ENDIF
!-----------------------------------------------------------------------
      IF(MODE==1.or.MODE==2)THEN
       DO I=1,NRES
        IF(NEEDC(I)>0)THEN
         DO J=1,NV
          CJAC(I,J)=0.0
         ENDDO
        ENDIF
       ENDDO
! - - - - - - - - - - - - - - - -
       K = 0
       DO IQ=1,NBF
        DO JQ=1,IQ
         K = K+1
         if(NEEDC(K)>0)then
          do i=1,nbf
           iIQ = i+(IQ-1)*nbf
           iJQ = i+(JQ-1)*nbf
           CJAC(K,iIQ) = CJAC(K,iIQ) + FC(i,JQ,USER(M9),COEF)
           IF(iJQ>=1)CJAC(K,iJQ) = CJAC(K,iJQ) + FC(i,IQ,USER(M9),COEF)
          enddo
         endif
        ENDDO
       ENDDO
      ENDIF
!-----------------------------------------------------------------------
      RETURN
      END
