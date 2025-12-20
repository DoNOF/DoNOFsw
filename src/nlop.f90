!-----------------------------------------------------------------------!
!                                                                       !
!  Polarizability (α), first (β) and second (γ) hyperpolarizabilities   !
!                                                                       !
!-----------------------------------------------------------------------!

! DIPNLOP
      SUBROUTINE DIPNLOP(NP,STEP,NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,    &
                         NSHELL,NPRIMI,ZAN,Cxyz,IAN,IMIN,IMAX,KSTART,   &
                         KATOM,KTYPE,KLOC,INTYP,KNG,KMIN,KMAX,ISH,ITYP, &
                         C1,C2,EX,CS,CP,CD,CF,CG,CH,CI,GRADS,IRUNTYP,   &
                         DIPS,SIZE_ENV,ENV,ATM,NBAS,BAS,IGTYP,NLOP,     &
                         ISOALPHA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!     ARGUMENTS
      INTEGER,INTENT(IN) :: NP,NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NLOP
      INTEGER,INTENT(IN) :: NSHELL,NPRIMI,IRUNTYP
      DOUBLE PRECISION,INTENT(IN) :: STEP
      DOUBLE PRECISION,DIMENSION(3,NAT) :: Cxyz
      DOUBLE PRECISION,DIMENSION(NAT),INTENT(IN) :: ZAN
      INTEGER,DIMENSION(NAT),INTENT(IN) :: IAN,IMIN,IMAX
      INTEGER,DIMENSION(NSHELL),INTENT(IN) :: KSTART,KATOM,KTYPE,KLOC
      INTEGER,DIMENSION(NSHELL),INTENT(IN) :: INTYP,KNG,KMIN,KMAX
      INTEGER,DIMENSION(NPRIMI),INTENT(IN) :: ISH,ITYP
      DOUBLE PRECISION,DIMENSION(NPRIMI),INTENT(IN) :: C1,C2,EX,CS,CP
      DOUBLE PRECISION,DIMENSION(NPRIMI),INTENT(IN) :: CD,CF,CG,CH,CI
      DOUBLE PRECISION,DIMENSION(3*NAT) :: GRADS
      DOUBLE PRECISION,DIMENSION(3) :: DIPS
      INTEGER :: SIZE_ENV,ATM(6,NAT),NBAS,BAS(8,NBAS),IGTYP
      DOUBLE PRECISION :: ENV(SIZE_ENV)
!-----------------------------------------------------------------------!
!                                                                       !
!     Compute tensor components of polarizability (α), first (β)        !
!     and second (γ) hyperpolarizabilities from dipole derivatives      !
!        w.r.t. the electric field Fz, using Romberg–Richardson         !
!                extrapolation of the dipole moment                     !
!                                                                       !
!     Flag NLOP controls which property to calculate:                   !
!                                                                       !
!        NLOP =  1 → ALPHA only                                         !
!        NLOP =  2 → BETA only                                          !
!        NLOP =  3 → GAMMA only                                         !
!        NLOP = -1 → All three                                          !
!                                                                       !
!     NPOINT: Number of steps used in the dyadic scaling of             !
!             the electric field. It represents the total number        !
!             of fields that will be considered in the calculations     !
!                                                                       !
!     STEP: Initial step size for the electric field. It defines        !
!           the base value of the field for the first step,             !
!           with subsequent values being scaled by powers of 2.         !
!                                                                       !
!     ISOALPHA: If=1, computes the diagonal components of the static    !
!               polarizability tensor (αxx, αyy, αzz) using the dyadic  !
!               Romberg–Richardson scheme looping the field direction   !
!               over x, y, z. Then reports the isotropic average and    !
!               the anisotropy (Raman convention):                      !
!                  α_iso = (αxx + αyy + αzz)/3                          !
!       α_aniso = sqrt( 0.5*((αxx-αyy)^2 + (αyy-αzz)^2 + (αzz-αxx)^2) ) !
!                                                                       !
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
      IF (NLOP==1) THEN    ! Polarizabilities
       IF (ISOALPHA==0) THEN
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
!                                  αzz                                  !
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
        CALL DIPNLOPz(NP,STEP,NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,       &
                       NSHELL,NPRIMI,ZAN,Cxyz,IAN,IMIN,IMAX,KSTART,     &
                       KATOM,KTYPE,KLOC,INTYP,KNG,KMIN,KMAX,ISH,ITYP,   &
                       C1,C2,EX,CS,CP,CD,CF,CG,CH,CI,GRADS,IRUNTYP,     &
                       DIPS,SIZE_ENV,ENV,ATM,NBAS,BAS,IGTYP,NLOP)
       ELSE IF (ISOALPHA==1) THEN
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
!                  α_iso = (αxx + αyy + αzz)/3                          !
!       α_aniso = sqrt( 0.5*((αxx-αyy)^2 + (αyy-αzz)^2 + (αzz-αxx)^2) ) !
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
        CALL DIPPOLARDIAG(NP,STEP,NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,   &
                          NSHELL,NPRIMI,ZAN,Cxyz,IAN,IMIN,IMAX,KSTART,  &
                          KATOM,KTYPE,KLOC,INTYP,KNG,KMIN,KMAX,ISH,ITYP,&
                          C1,C2,EX,CS,CP,CD,CF,CG,CH,CI,GRADS,IRUNTYP,  &
                          DIPS,SIZE_ENV,ENV,ATM,NBAS,BAS,IGTYP,NLOP)
       ENDIF
      END IF
!     Hyperpolarizabilities (βzzz,γzzzz) & Polarizability αzz (NLOP==-1)
      IF( NLOP==2 .or. NLOP==3 .or. NLOP==-1 )THEN
       CALL DIPNLOPz(NP,STEP,NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,        &
                     NSHELL,NPRIMI,ZAN,Cxyz,IAN,IMIN,IMAX,KSTART,       &
                     KATOM,KTYPE,KLOC,INTYP,KNG,KMIN,KMAX,ISH,ITYP,     &
                     C1,C2,EX,CS,CP,CD,CF,CG,CH,CI,GRADS,IRUNTYP,       &
                     DIPS,SIZE_ENV,ENV,ATM,NBAS,BAS,IGTYP,NLOP)
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
!             αxz,αyz,αzz, βxzz,βyzz,βzzz, γxzzz,γyzzz,γzzzz            !
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
!     CALL DIPNLOPxyz(NP,STEP,NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,       &
!                     NSHELL,NPRIMI,ZAN,Cxyz,IAN,IMIN,IMAX,KSTART,      &
!                     KATOM,KTYPE,KLOC,INTYP,KNG,KMIN,KMAX,ISH,ITYP,    &
!                     C1,C2,EX,CS,CP,CD,CF,CG,CH,CI,GRADS,IRUNTYP,      &
!                     DIPS,SIZE_ENV,ENV,ATM,NBAS,BAS,IGTYP,NLOP)
!-----------------------------------------------------------------------
      RETURN
      END

!-----------------------------------------------------------------------
!     Compute Z components of polarizability and hyperpolarizabilities
!                    Static Electric Field: (0,0,Fz)
!-----------------------------------------------------------------------

! DIPNLOPz
      SUBROUTINE DIPNLOPz(NP,STEP,NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,   &
                          NSHELL,NPRIMI,ZAN,Cxyz,IAN,IMIN,IMAX,KSTART,  &
                          KATOM,KTYPE,KLOC,INTYP,KNG,KMIN,KMAX,ISH,ITYP,&
                          C1,C2,EX,CS,CP,CD,CF,CG,CH,CI,GRADS,IRUNTYP,  &
                          DIPS,SIZE_ENV,ENV,ATM,NBAS,BAS,IGTYP,NLOP)
!-----------------------------------------------------------------------
!     Romberg–Richardson (full triangles) for derivatives of DIPz(F):
!     ALPHAzz = dDIPz/dF, BETAzzz = d2DIPz/dF2, GAMMAzzzz = d3DIPz/dF3
!
!       - Dyadic Scaling: h_k = STEP*2**(k-1), k=1..NP
!
!       - Centered schemes based O(h^2):
!           alpha(h) = [DIPz(+h)-DIPz(-h)]/(2h)
!           beta(h)  = [DIPz(+h)-2DIPz(0)+DIPz(-h)]/h^2
!           gamma(h) = [DIPz(+2h)-2DIPz(+h)+2DIPz(-h)-DIPz(-2h)]/(2h^3)
!
!       - Romberg triangle: p=2 for all three cases.
!
!       - Selecting the "best" value: minimum |R(k,j)-R(k,j-1)|;
!         If there is a tie, choose the one with the lowest |R(k,j)|.
!
!     Flag NLOP controls which property to calculate:
!
!        NLOP =  1 → ALPHAzz only
!        NLOP =  2 → BETAzzz only
!        NLOP =  3 → GAMMAzzzz only
!        NLOP = -1 → All three
!
!-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL EFIELDL,RESTART
      COMMON/INP_EFIELDL/EFX,EFY,EFZ,EFIELDL
      COMMON/INPNOF_RSTRT/RESTART,INPUTGAMMA,INPUTC,INPUTFMIUG,INPUTCXYZ
      COMMON/INPNOF_RHF/CONVRHFDM,IRHF,IRHFTYP,NCONVRHF,MAXITRHF
      COMMON/INPNOF_PRINT/NPRINT,IWRITEC,IMULPOP,IAIMPAC,IFCHK
      COMMON/INPNOF_MOLDEN/MOLDEN
      COMMON/EHFEN/EHF,EN
      COMMON/ENERGIAS/EELEC,EELEC_OLD,DIF_EELEC,EELEC_MIN
      COMMON/INPNOF_COEFOPT/MAXLOOP
!     ARGUMENTS
      INTEGER,INTENT(IN) :: NP,NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux
      INTEGER,INTENT(IN) :: NSHELL,NPRIMI,IRUNTYP,NLOP
      DOUBLE PRECISION,INTENT(IN) :: STEP
      DOUBLE PRECISION,DIMENSION(3,NAT) :: Cxyz
      DOUBLE PRECISION,DIMENSION(NAT),INTENT(IN) :: ZAN
      INTEGER,DIMENSION(NAT),INTENT(IN) :: IAN,IMIN,IMAX
      INTEGER,DIMENSION(NSHELL),INTENT(IN) :: KSTART,KATOM,KTYPE,KLOC
      INTEGER,DIMENSION(NSHELL),INTENT(IN) :: INTYP,KNG,KMIN,KMAX
      INTEGER,DIMENSION(NPRIMI),INTENT(IN) :: ISH,ITYP
      DOUBLE PRECISION,DIMENSION(NPRIMI),INTENT(IN) :: C1,C2,EX,CS,CP
      DOUBLE PRECISION,DIMENSION(NPRIMI),INTENT(IN) :: CD,CF,CG,CH,CI
      DOUBLE PRECISION,DIMENSION(3*NAT) :: GRADS
      DOUBLE PRECISION,DIMENSION(3) :: DIPS
      INTEGER :: SIZE_ENV,ATM(6,NAT),NBAS,BAS(8,NBAS),IGTYP
      DOUBLE PRECISION :: ENV(SIZE_ENV)
!
      INTEGER :: MAXM, M, I, I0, FINDX, IPrintTria
      INTEGER :: ibA, jbA, ibB, jbB, ibG, jbG
      DOUBLE PRECISION :: TOL, DMU(3)
      DOUBLE PRECISION :: bestA, bestB, bestG, errA, errB, errG
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: FV, DIPz, HPOS
!
      EXTERNAL FINDX
!-----------------------------------------------------------------------
      MAXM = 2*NP + 3
      ALLOCATE(FV(MAXM), DIPz(MAXM), HPOS(NP))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Build dyadic ladder and evaluate DIPz at each required field
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      TOL = 1.0D-16
      DMU  = 0.0D0
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Energy without Electric Field
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL OriginalOptions(0)   ! Save the original settings
      IF(EFX==0.0d0 .and. EFY==0.0d0 .and. EFZ==0.0d0)THEN
       CALL ENERGRAD(NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NSHELL,NPRIMI,  &
                     ZAN,Cxyz,IAN,IMIN,IMAX,KSTART,KATOM,KTYPE,KLOC,    &
                     INTYP,KNG,KMIN,KMAX,ISH,ITYP,C1,C2,EX,CS,CP,CD,CF, &
                     CG,CH,CI,GRADS,IRUNTYP,DIPS,SIZE_ENV,ENV,ATM,NBAS, &
                     BAS,IGTYP,0,1)
       WRITE(6,1)
       WRITE(6,2)
       WRITE(6,3)EFZ,DIPS(3),EELEC+EN
      ELSE
       WRITE(6,4)
       STOP
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Ensure F=0 (M=1) present
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      M = 1
      FV(M) = EFZ
      DIPz(M) = DIPS(3)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Energy and Dipole for I=1,NP+1 (EFZ>0)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL OriginalOptions(1)   ! Impose mandatory options
      DO I = 1, NP+1
       EFZ = STEP * 2.0D0**(I-1)
       M = M + 1
       MAXLOOP = 10
       CALL ENERGRAD(NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NSHELL,         &
                     NPRIMI,ZAN,Cxyz,IAN,IMIN,IMAX,KSTART,KATOM,        &
                     KTYPE,KLOC,INTYP,KNG,KMIN,KMAX,ISH,ITYP,C1,C2,EX,  &
                     CS,CP,CD,CF,CG,CH,CI,GRADS,IRUNTYP,DMU,SIZE_ENV,   &
                     ENV,ATM,NBAS,BAS,IGTYP,0,0)
       FV(M) = EFZ
       DIPz(M) = DMU(3)
       WRITE(6,3) FV(M), DIPz(M), EELEC+EN
      END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Recalculate the Energy without Electric Field
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL OriginalOptions(2)   ! Restore the original settings
      CALL ENERGRAD(NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NSHELL,NPRIMI,   &
                    ZAN,Cxyz,IAN,IMIN,IMAX,KSTART,KATOM,KTYPE,KLOC,     &
                    INTYP,KNG,KMIN,KMAX,ISH,ITYP,C1,C2,EX,CS,CP,CD,CF,  &
                    CG,CH,CI,GRADS,IRUNTYP,DMU,SIZE_ENV,ENV,ATM,NBAS,   &
                    BAS,IGTYP,0,0)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Energy and Dipole for I=1,NP+1 (EFZ<0)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL OriginalOptions(1)   ! Impose mandatory options
      DO I = 1, NP+1
       EFZ = - STEP * 2.0D0**(I-1)
       M = M + 1
       MAXLOOP = 10
       CALL ENERGRAD(NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NSHELL,         &
                     NPRIMI,ZAN,Cxyz,IAN,IMIN,IMAX,KSTART,KATOM,        &
                     KTYPE,KLOC,INTYP,KNG,KMIN,KMAX,ISH,ITYP,C1,C2,EX,  &
                     CS,CP,CD,CF,CG,CH,CI,GRADS,IRUNTYP,DMU,SIZE_ENV,   &
                     ENV,ATM,NBAS,BAS,IGTYP,0,0)
       FV(M) = EFZ
       DIPz(M) = DMU(3)
       WRITE(6,3) FV(M), DIPz(M), EELEC+EN
      END DO
      CALL OriginalOptions(2)   ! Restore the original settings
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Sort FV in ascending order
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL SORT2z(FV, DIPz, M)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Locate index of F=0 (after sorting)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      I0 = FINDX(FV, M, 0.0D0, TOL)
      IF (I0 <= 0) THEN
       WRITE(6,5)
       STOP
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Dyadic steps:  h, 2h, 4h, …  (ascending powers of two)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO I = 1, NP
       HPOS(I) = STEP * 2.0D0**(I-1)
      END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Build and process full Romberg–Richardson triangles for the first,
!     second and third derivatives (α, β, γ) of the dipole μ(F) w.r.t.
!     the electric field F.
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      WRITE(6,9)
      CALL ROMBERG_TRIANGLE(FV, DIPz, M, I0, NP, TOL, HPOS,             &
                            bestA, errA, ibA, jbA,                      &
                            bestB, errB, ibB, jbB,                      &
                            bestG, errG, ibG, jbG, NLOP, 1)
      IF (NLOP == -1) THEN
       WRITE(6,6) bestA, errA, ibA, jbA
       WRITE(6,7) bestB, errB, ibB, jbB
       WRITE(6,8) bestG, errG, ibG, jbG
      ELSE IF (NLOP == 1) THEN
       WRITE(6,6) bestA, errA, ibA, jbA
      ELSE IF (NLOP == 2) THEN
       WRITE(6,7) bestB, errB, ibB, jbB
      ELSE IF (NLOP == 3) THEN
       WRITE(6,8) bestG, errG, ibG, jbG
      END IF
      WRITE(6,9)
!-----------------------------------------------------------------------
!     Format definitions
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    1 FORMAT(//2X,'Non Linear Optical Properties (Romberg–Richardson)', &
            /1X,'====================================================')
    2 FORMAT(/,2X,'Field',16X,'Dipole_z',17X,'Energy',/)
    3 FORMAT(1X,F7.4,2X,F25.15,2X,F25.15)
    4 FORMAT(/,2X,'Stop: The Electric Field must be 0 for NLOP' )
    5 FORMAT(/,2X,'Error: F=0 not found after sorting')
    6 FORMAT(/,' Best ALPHAzz   =',1X,ES16.8,' ±',ES10.3,               &
               '  at (k=',I2,', j=',I2,')')
    7 FORMAT(/,' Best BETAzzz   =',1X,ES16.8,' ±',ES10.3,               &
               '  at (k=',I2,', j=',I2,')')
    8 FORMAT(/,' Best GAMMAzzzz =',1X,ES16.8,' ±',ES10.3,               &
               '  at (k=',I2,', j=',I2,')')
    9 FORMAT(/118('-'))
!-----------------------------------------------------------------------
      DEALLOCATE(FV,DIPz,HPOS)
      RETURN
      END

! OriginalOptions
      SUBROUTINE OriginalOptions(IOPTION)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER,INTENT(IN) :: IOPTION
      LOGICAL EFIELDL,EFIELDL0,RESTART,RESTART0
      COMMON/INP_EFIELDL/EFX,EFY,EFZ,EFIELDL
      COMMON/INPNOF_RSTRT/RESTART,INPUTGAMMA,INPUTC,INPUTFMIUG,INPUTCXYZ
      COMMON/INPNOF_RHF/CONVRHFDM,IRHF,IRHFTYP,NCONVRHF,MAXITRHF
      COMMON/INPNOF_HFCONVTECH/HFDAMP,HFEXTRAP,HFDIIS
      COMMON/INPNOF_PRINT/NPRINT,IWRITEC,IMULPOP,IAIMPAC,IFCHK
      COMMON/INPNOF_MOLDEN/MOLDEN
      COMMON/INPNOF_COEFOPT/MAXLOOP
      COMMON/OriOptions/RESTART0,INPUTGAMMA0,INPUTC0,INPUTFMIUG0,       &
                        INPUTCXYZ0,IRHF0,NPRINT0,IWRITEC0,IMULPOP0,     &
                        IAIMPAC0,IFCHK0,MOLDEN0,EFX0,EFY0,EFZ0,EFIELDL0
!-----------------------------------------------------------------------
      IF(IOPTION==0)THEN             ! Save the original settings
       EFX0 = EFX
       EFY0 = EFY
       EFZ0 = EFZ
       EFIELDL0 = EFIELDL
!
       RESTART0 = RESTART
       INPUTGAMMA0 = INPUTGAMMA
       INPUTC0 = INPUTC
       INPUTFMIUG0 = INPUTFMIUG
       INPUTCXYZ0 = INPUTCXYZ
!
       IRHF0 = IRHF
       CONVRHFDM0 = CONVRHFDM
       IRHFTYP0 = IRHFTYP
       NCONVRHF0 = NCONVRHF
       MAXITRHF0 = MAXITRHF
       HFDAMP0   = HFDAMP
       HFEXTRAP0 = HFEXTRAP
       HFDIIS0   = HFDIIS
!
       NPRINT0  = NPRINT
       IWRITEC0 = IWRITEC
       IMULPOP0 = IMULPOP
       IAIMPAC0 = IAIMPAC
       IFCHK0 = IFCHK
       MOLDEN0 = MOLDEN
      ELSE IF(IOPTION==1)THEN        ! Impose mandatory options
       EFX = EFX0
       EFY = EFY0
       EFIELDL = .TRUE.
!
       RESTART = .TRUE.
       INPUTGAMMA = 1
       INPUTC = 1
       INPUTFMIUG = 1
       INPUTCXYZ = 1
!
       IRHF  = 0
       NPRINT  = 0
       IWRITEC = 0
       IMULPOP = 0
       IAIMPAC = 0
       IFCHK = 0
       MOLDEN = 0
       MAXLOOP = 10                  ! ADAM Method
      ELSE IF(IOPTION==2)THEN        ! Restore the original settings
       EFX = EFX0
       EFY = EFY0
       EFZ = EFZ0
       EFIELDL = EFIELDL0
!
       RESTART = RESTART0
       INPUTGAMMA = INPUTGAMMA0
       INPUTC = INPUTC0
       INPUTFMIUG = INPUTFMIUG0
       INPUTCXYZ = INPUTCXYZ0
!
       IRHF = IRHF0
!      CONVRHFDM = CONVRHFDM0
!      IRHFTYP  = IRHFTYP0
!      NCONVRHF = NCONVRHF0
!      MAXITRHF = MAXITRHF0
!      HFDAMP   = HFDAMP0
!      HFEXTRAP = HFEXTRAP0
!      HFDIIS   = HFDIIS0
!
       NPRINT  = NPRINT0
       IWRITEC = IWRITEC0
       IMULPOP = IMULPOP0
       IAIMPAC = IAIMPAC0
       IFCHK = IFCHK0
       MOLDEN = MOLDEN0
       MAXLOOP = 10                  ! ADAM Method
      ENDIF
!-----------------------------------------------------------------------
      RETURN
      END

! ROMBERG_TRIANGLE
      SUBROUTINE ROMBERG_TRIANGLE(FV,DIPz,M,I0,NP,TOL,HPOS,             &
                                  bestA,errA,ibA,jbA,                   &
                                  bestB,errB,ibB,jbB,                   &
                                  bestG,errG,ibG,jbG,NLOP,IPrintTria)
!-----------------------------------------------------------------------
!     Uses centered finite differences (O(h^2)) for a dyadic ladder of
!     step sizes h, 2h, 4h, …  and constructs Romberg extrapolation
!     triangles (p=2) for α=dμ/dF, β=d²μ/dF², γ=d³μ/dF³.
!
!     Output:
!     bestA, bestB, bestG - Best α, β, γ selected from Romberg triangles
!     errA,  errB,  errG  - Corresponding errors (Δ between best pair)
!-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER M, I0, NP, I, K, J, ip, im, i2p, i2m, FINDX, IPrintTria
      INTEGER :: ibA, jbA, ibB, jbB, ibG, jbG
      DOUBLE PRECISION FV(M), DIPz(M), TOL
      DOUBLE PRECISION bestA, bestB, bestG, errA, errB, errG
!
      DOUBLE PRECISION HPOS(NP)
      DOUBLE PRECISION ABASE(NP), BBASE(NP), GBASE(NP)
      DOUBLE PRECISION RAL(NP,NP), RBE(NP,NP), RGA(NP,NP)
      DOUBLE PRECISION h, denom
!-----------------------------------------------------------------------
!     Centered finite-difference base derivatives O(h^2)
!-----------------------------------------------------------------------
      DO K=1,NP
       h = HPOS(K)
!      Find indices for ±h and ±2h
       ip  = FINDX(FV,M,+h ,TOL)
       im  = FINDX(FV,M,-h ,TOL)
       i2p = FINDX(FV,M,+2.0D0*h,TOL)
       i2m = FINDX(FV,M,-2.0D0*h,TOL)
!      First derivative: α ≈ [μ(+h) - μ(-h)] / (2h)
       IF (ip>0 .and. im>0) THEN
        ABASE(K) = (DIPz(ip)-DIPz(im))/(2.0D0*h)
       ELSE
        ABASE(K) = 0.0D0
       END IF
!      Second derivative: β ≈ [μ(+h) - 2μ(0) + μ(-h)] / h²
       IF (ip>0 .and. im>0) THEN
        BBASE(K) = (DIPz(ip)-2.0D0*DIPz(I0)+DIPz(im))/(h*h)
       ELSE
        BBASE(K) = 0.0D0
       END IF
!      Third derivative: γ ≈ [μ(2h)-2μ(h)+2μ(-h)-μ(-2h)] / (2h³)
       IF (i2p>0 .and. ip>0 .and. im>0 .and. i2m>0) THEN
        GBASE(K) = ( DIPz(i2p) - 2.0D0*DIPz(ip) + 2.0D0*DIPz(im)        &
                   - DIPz(i2m) ) / (2.0D0*h*h*h)
       ELSE
        GBASE(K) = 0.0D0
       END IF
      END DO
!-----------------------------------------------------------------------
!     Initialize Romberg triangles
!-----------------------------------------------------------------------
      RAL = 0.0D0
      RBE = 0.0D0
      RGA = 0.0D0
!
      DO K=1,NP
       RAL(K,1) = ABASE(K)
       RBE(K,1) = BBASE(K)
       RGA(K,1) = GBASE(K)
      END DO
!-----------------------------------------------------------------------
!     Romberg extrapolation (p=2) for α, β, γ
!-----------------------------------------------------------------------
      DO J=2,NP
       denom = 2.0D0**(2*(J-1)) - 1.0D0
       DO K=1,NP-J+1
        RAL(K,J) = RAL(K,J-1) + (RAL(K,J-1)-RAL(K+1,J-1))/denom
        RBE(K,J) = RBE(K,J-1) + (RBE(K,J-1)-RBE(K+1,J-1))/denom
        RGA(K,J) = RGA(K,J-1) + (RGA(K,J-1)-RGA(K+1,J-1))/denom
       END DO
      END DO
!-----------------------------------------------------------------------
!     Print full Romberg triangles for inspection
!-----------------------------------------------------------------------
      IF(IPrintTria==1)THEN
       IF (NLOP == -1) THEN
        WRITE(6,1) 'Alpha (dmu/dF) - Romberg triangle'
        CALL PRINT_TRIANGLE(RAL,NP)
        WRITE(6,1) 'Beta (d2mu/dF2) - Romberg triangle'
        CALL PRINT_TRIANGLE(RBE,NP)
        WRITE(6,1) 'Gamma (d3mu/dF3) - Romberg triangle'
        CALL PRINT_TRIANGLE(RGA,NP)
       ELSE IF (NLOP == 1) THEN
        WRITE(6,1) 'Alpha (dmu/dF) - Romberg triangle'
        CALL PRINT_TRIANGLE(RAL,NP)
       ELSE IF (NLOP == 2) THEN
        WRITE(6,1) 'Beta (d2mu/dF2) - Romberg triangle'
        CALL PRINT_TRIANGLE(RBE,NP)
       ELSE IF (NLOP == 3) THEN
        WRITE(6,1) 'Gamma (d3mu/dF3) - Romberg triangle'
        CALL PRINT_TRIANGLE(RGA,NP)
       END IF
      END IF
!-----------------------------------------------------------------------
!     Select best values by column difference criterion
!-----------------------------------------------------------------------
      CALL SELECT_BEST_COL(RAL,NP,bestA,errA,ibA,jbA)
      CALL SELECT_BEST_COL(RBE,NP,bestB,errB,ibB,jbB)
      CALL SELECT_BEST_COL(RGA,NP,bestG,errG,ibG,jbG)
!-----------------------------------------------------------------------
!     Format statements
!-----------------------------------------------------------------------
    1 FORMAT(/1X,A,/)
!-----------------------------------------------------------------------
      RETURN
      END

! PRINT_TRIANGLE
      SUBROUTINE PRINT_TRIANGLE(R,NP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION :: R(NP,NP)
      INTEGER :: NP, K, J
!-----------------------------------------------------------------------
!     Prints a full Romberg triangle row by row
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO K=1,NP
       WRITE(6,'(1X,*(1X,ES12.5))') (R(K,J), J=1,NP-K+1)
      END DO
!-----------------------------------------------------------------------
      RETURN
      END

! SELECT_BEST_COL
      SUBROUTINE SELECT_BEST_COL(R,NP,VAL,ERR,KB,JB)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Criterion:
!       - Scan columns j = 1,...,NP-1
!       - In each column, find vertical pair (k,k+1) minimizing
!           Δ = |R(k+1,j) - R(k,j)|
!       - Report VAL = R(k, j+1) (next column, upper element)
!       - Report ERR = Δ from column j
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER NP, KB, JB, K, J
      DOUBLE PRECISION R(NP,NP), VAL, ERR
      DOUBLE PRECISION delta, bestDelta, cand_next
!-----------------------------------------------------------------------
      bestDelta = 1.0D300
      VAL  = 0.0D0
      ERR  = 0.0D0
      KB   = 1
      JB   = 2
!
      DO j = 1, NP-1
       DO k = 1, NP-j
        delta = DABS( R(k+1,j) - R(k,j) )
        IF (delta < bestDelta) THEN
         bestDelta = delta
         cand_next = R(k, j+1)   ! next column, row k (superior)
         VAL  = cand_next
         ERR  = delta
         KB   = k
         JB   = j+1
        ELSE IF (delta == bestDelta) THEN
!        Tie-breaker: choose the smaller |R(k,j+1)|
         IF (DABS(R(k, j+1)) < DABS(VAL)) THEN
          VAL  = R(k, j+1)
          ERR  = delta
          KB   = k
          JB   = j+1
         END IF
        END IF
       END DO
      END DO
!-----------------------------------------------------------------------
      RETURN
      END

! SORT2z
      SUBROUTINE SORT2z(A,Bz,M)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER M, I, J, K
      DOUBLE PRECISION A(*), Bz(*), TA, TBz
      DO I=1,M-1
       K = I
       DO J=I+1,M
        IF (A(J) < A(K)) K = J
       END DO
       IF (K /= I) THEN
        TA = A(I)
        TBz = Bz(I)
        A(I) = A(K)
        Bz(I) = Bz(K)
        A(K) = TA
        Bz(K) = TBz
       END IF
      END DO
      END SUBROUTINE SORT2z

! FINDX
      INTEGER FUNCTION FINDX(A,N,X,TOL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION A(*), X, TOL
      INTEGER N, I
!-----------------------------------------------------------------------
      FINDX = -1
      DO I=1,N
       IF (DABS(A(I)-X) <= TOL) THEN
        FINDX = I
        RETURN
       END IF
      END DO
!-----------------------------------------------------------------------
      END FUNCTION FINDX

!-----------------------------------------------------------------------
!     Compute all components of polarizability and hyperpolarizabilities
!                    Static Electric Field: (0,0,Fz)
!
!     MUi = MU(0)i + Sum_j ALPHAij*Fj + 1/2 Sum_jk BETAijk*Fj*Fk
!                                     + 1/6 Sum_jkl GAMMAijkl*Fj*Fk*Fl
!
!     Afected Components: ALPHAiz, BETAizz, GAMMAizzz; i=x,y,z
!
!-----------------------------------------------------------------------

! DIPNLOPxyz
      SUBROUTINE DIPNLOPxyz(NP,STEP,NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux, &
                            NSHELL,NPRIMI,ZAN,Cxyz,IAN,IMIN,IMAX,KSTART,&
                            KATOM,KTYPE,KLOC,INTYP,KNG,KMIN,KMAX,ISH,   &
                            ITYP,C1,C2,EX,CS,CP,CD,CF,CG,CH,CI,GRADS,   &
                            IRUNTYP,DIPS,SIZE_ENV,ENV,ATM,NBAS,BAS,     &
                            IGTYP,NLOP)
!-----------------------------------------------------------------------
!     Compute components of polarizability and hyperpolarizabilities
!     from dipole derivatives induced by an electric field applied
!     exclusively along z (Fz), using Romberg–Richardson
!     extrapolation on a dyadic ladder of steps.
!
!     Calculates:
!      Alpha components: α_xz, α_yz, α_zz
!      Beta  components: β_xzz, β_yzz, β_zzz
!      Gamma components: γ_xzzz, γ_yzzz, γ_zzzz
!
!     Selecting the "best" value: minimum |R(k,j)-R(k,j-1)|;
!     If there is a tie, choose the one with the lowest |R(k,j)|.
!
!     Flag NLOP controls which property to calculate:
!
!        NLOP =  1 → ALPHA only
!        NLOP =  2 → BETA only
!        NLOP =  3 → GAMMA only
!        NLOP = -1 → All three
!
!-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL EFIELDL,RESTART
      COMMON/INP_EFIELDL/EFX,EFY,EFZ,EFIELDL
      COMMON/INPNOF_RSTRT/RESTART,INPUTGAMMA,INPUTC,INPUTFMIUG,INPUTCXYZ
      COMMON/INPNOF_RHF/CONVRHFDM,IRHF,IRHFTYP,NCONVRHF,MAXITRHF
      COMMON/INPNOF_PRINT/NPRINT,IWRITEC,IMULPOP,IAIMPAC,IFCHK
      COMMON/INPNOF_MOLDEN/MOLDEN
      COMMON/EHFEN/EHF,EN
      COMMON/ENERGIAS/EELEC,EELEC_OLD,DIF_EELEC,EELEC_MIN
      COMMON/INPNOF_COEFOPT/MAXLOOP
!
      INTEGER,INTENT(IN) :: NP,NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NLOP
      INTEGER,INTENT(IN) :: NSHELL,NPRIMI,IRUNTYP
      DOUBLE PRECISION,INTENT(IN) :: STEP
      DOUBLE PRECISION,DIMENSION(3,NAT) :: Cxyz
      DOUBLE PRECISION,DIMENSION(NAT),INTENT(IN) :: ZAN
      INTEGER,DIMENSION(NAT),INTENT(IN) :: IAN,IMIN,IMAX
      INTEGER,DIMENSION(NSHELL),INTENT(IN) :: KSTART,KATOM,KTYPE,KLOC
      INTEGER,DIMENSION(NSHELL),INTENT(IN) :: INTYP,KNG,KMIN,KMAX
      INTEGER,DIMENSION(NPRIMI),INTENT(IN) :: ISH,ITYP
      DOUBLE PRECISION,DIMENSION(NPRIMI),INTENT(IN) :: C1,C2,EX,CS,CP
      DOUBLE PRECISION,DIMENSION(NPRIMI),INTENT(IN) :: CD,CF,CG,CH,CI
      DOUBLE PRECISION,DIMENSION(3*NAT) :: GRADS
      DOUBLE PRECISION,DIMENSION(3) :: DIPS
      INTEGER :: SIZE_ENV,ATM(6,NAT),NBAS,BAS(8,NBAS),IGTYP
      DOUBLE PRECISION :: ENV(SIZE_ENV)
!
      INTEGER :: MAXM, M, I, I0
      DOUBLE PRECISION :: TOL
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: FV, HPOS
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: DIPx, DIPy, DIPz
      DOUBLE PRECISION :: DMU(3)
      DOUBLE PRECISION :: bestAx, bestAy, bestAz
      DOUBLE PRECISION :: bestBx, bestBy, bestBz
      DOUBLE PRECISION :: bestGx, bestGy, bestGz
      DOUBLE PRECISION :: errAx, errAy, errAz
      DOUBLE PRECISION :: errBx, errBy, errBz
      DOUBLE PRECISION :: errGx, errGy, errGz
!
      INTEGER FINDX
      EXTERNAL FINDX
!-----------------------------------------------------------------------
      MAXM = 2*NP + 3
      ALLOCATE(FV(MAXM), DIPx(MAXM), DIPy(MAXM), DIPz(MAXM), HPOS(NP))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Build dyadic ladder and evaluate DIP at each unique required field
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      TOL = 1D-16
      DMU  = 0.0D0
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Energy without Electric Field (F=0)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL OriginalOptions(0)   ! Save the original settings
      IF(EFX==0.0d0 .and. EFY==0.0d0 .and. EFZ==0.0d0)THEN
       CALL ENERGRAD(NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NSHELL,NPRIMI,  &
                     ZAN,Cxyz,IAN,IMIN,IMAX,KSTART,KATOM,KTYPE,KLOC,    &
                     INTYP,KNG,KMIN,KMAX,ISH,ITYP,C1,C2,EX,CS,CP,CD,CF, &
                     CG,CH,CI,GRADS,IRUNTYP,DIPS,SIZE_ENV,ENV,ATM,NBAS, &
                     BAS,IGTYP,0,1)
       WRITE(6,1)
       WRITE(6,2)
       WRITE(6,3)EFZ,DIPS(3),EELEC+EN
      ELSE
       WRITE(6,4)
       STOP
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Ensure F=0 (M=1) present
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      M = 1
      FV(M)  = EFZ
      DIPx(M) = DIPS(1)
      DIPy(M) = DIPS(2)
      DIPz(M) = DIPS(3)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Energy and Dipole for I=1,NP+1 (Positive Fields: EFZ>0)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL OriginalOptions(1)   ! Impose mandatory options
      DO I = 1, NP+1
       EFZ = STEP * 2.0D0**(I-1)
       M = M + 1
       MAXLOOP = 10
       CALL ENERGRAD(NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NSHELL,         &
                     NPRIMI,ZAN,Cxyz,IAN,IMIN,IMAX,KSTART,KATOM,        &
                     KTYPE,KLOC,INTYP,KNG,KMIN,KMAX,ISH,ITYP,C1,C2,EX,  &
                     CS,CP,CD,CF,CG,CH,CI,GRADS,IRUNTYP,DMU,SIZE_ENV,   &
                     ENV,ATM,NBAS,BAS,IGTYP,0,0)
       FV(M) = EFZ
       DIPx(M) = DMU(1)
       DIPy(M) = DMU(2)
       DIPz(M) = DMU(3)
       WRITE(6,3)FV(M), DIPz(M), EELEC+EN
      END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Recalculate the Energy without Electric Field
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL OriginalOptions(2)   ! Restore the original settings
      CALL ENERGRAD(NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NSHELL,NPRIMI,   &
                    ZAN,Cxyz,IAN,IMIN,IMAX,KSTART,KATOM,KTYPE,KLOC,     &
                    INTYP,KNG,KMIN,KMAX,ISH,ITYP,C1,C2,EX,CS,CP,CD,CF,  &
                    CG,CH,CI,GRADS,IRUNTYP,DMU,SIZE_ENV,ENV,ATM,NBAS,   &
                    BAS,IGTYP,0,0)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Energy and Dipole for I=1,NP+1 (EFZ<0)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL OriginalOptions(1)   ! Impose mandatory options
      DO I = 1, NP+1
       EFZ = - STEP * 2.0D0**(I-1)
       M = M + 1
       MAXLOOP = 10
       CALL ENERGRAD(NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NSHELL,         &
                     NPRIMI,ZAN,Cxyz,IAN,IMIN,IMAX,KSTART,KATOM,        &
                     KTYPE,KLOC,INTYP,KNG,KMIN,KMAX,ISH,ITYP,C1,C2,EX,  &
                     CS,CP,CD,CF,CG,CH,CI,GRADS,IRUNTYP,DMU,SIZE_ENV,   &
                     ENV,ATM,NBAS,BAS,IGTYP,0,0)
       FV(M) = EFZ
       DIPx(M) = DMU(1)
       DIPy(M) = DMU(2)
       DIPz(M) = DMU(3)
       WRITE(6,3)FV(M), DIPz(M), EELEC+EN
      END DO
      CALL OriginalOptions(2)   ! Restore the original settings
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Sort by Field in ascending order and reorder dipoles
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL SORT2xyz(FV, DIPx, DIPy, DIPz, M)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Locate index of F=0 (after sorting)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      I0 = FINDX(FV,M,0.0D0,TOL)
      IF (I0 <= 0) THEN
       WRITE(6,5)
       STOP
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Dyadic steps:  h, 2h, 4h, …  (ascending powers of two)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO I = 1, NP
       HPOS(I) = STEP * 2.0D0**(I-1)
      END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Process Romberg triangles for each dipole component
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IPrintTria = 1
!     X-component
      if(IPrintTria == 1)then
       write(6,16)
       write(6,'(/,1X,''X-component:'')')
      end if
      CALL ROMBERG_TRIANGLE(FV, DIPx, M, I0, NP, TOL, HPOS,             &
                            bestAx, errAx, ibAx, jbAx,                  &
                            bestBx, errBx, ibBx, jbBx,                  &
                            bestGx, errGx, ibGx, jbGx,                  &
                            NLOP, IPrintTria)
!     Y-component
      if(IPrintTria == 1)then
       write(6,16)
       write(6,'(/,1X,''Y-component:'')')
      end if
      CALL ROMBERG_TRIANGLE(FV, DIPy, M, I0, NP, TOL, HPOS,             &
                            bestAy, errAy, ibAy, jbAy,                  &
                            bestBy, errBy, ibBy, jbBy,                  &
                            bestGy, errGy, ibGy, jbGy,                  &
                            NLOP, IPrintTria)
!     Z-component
      if(IPrintTria == 1)then
       write(6,16)
       write(6,'(/,1X,''Z-component:'')')
      end if
      CALL ROMBERG_TRIANGLE(FV,DIPz,M,I0,NP,TOL,HPOS,                   &
                            bestAz, errAz, ibAz, jbAz,                  &
                            bestBz, errBz, ibBz, jbBz,                  &
                            bestGz, errGz, ibGz, jbGz,                  &
                            NLOP, IPrintTria)
!-----------------------------------------------------------------------
!     Print tensor summary
!-----------------------------------------------------------------------
      WRITE(6,6)
      IF (NLOP == -1) THEN
!      x
       WRITE(6,7)  bestAx, errAx, ibAx, jbAx
       WRITE(6,8)  bestBx, errBx, ibBx, jbBx
       WRITE(6,9)  bestGx, errGx, ibGx, jbGx
!      y
       WRITE(6,10) bestAy, errAy, ibAy, jbAy
       WRITE(6,11) bestBy, errBy, ibBy, jbBy
       WRITE(6,12) bestGy, errGy, ibGy, jbGy
!      z
       WRITE(6,13) bestAz, errAz, ibAz, jbAz
       WRITE(6,14) bestBz, errBz, ibBz, jbBz
       WRITE(6,15) bestGz, errGz, ibGz, jbGz
      ELSE IF (NLOP == 1) THEN
       WRITE(6,7)  bestAx, errAx, ibAx, jbAx
       WRITE(6,10) bestAy, errAy, ibAy, jbAy
       WRITE(6,13) bestAz, errAz, ibAz, jbAz
      ELSE IF (NLOP == 2) THEN
       WRITE(6,8)  bestBx, errBx, ibBx, jbBx
       WRITE(6,11) bestBy, errBy, ibBy, jbBy
       WRITE(6,14) bestBz, errBz, ibBz, jbBz
      ELSE IF (NLOP == 3) THEN
       WRITE(6,9)  bestGx, errGx, ibGx, jbGx
       WRITE(6,12) bestGy, errGy, ibGy, jbGy
       WRITE(6,15) bestGz, errGz, ibGz, jbGz
      END IF
      WRITE(6,16)
!-----------------------------------------------------------------------
!     Format definitions
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    1 FORMAT(//2X,'Non Linear Optical Properties (Romberg–Richardson)', &
             /1X,'====================================================')
    2 FORMAT(/,2X,'Field',16X,'Dipole_z',17X,'Energy',/)
    3 FORMAT(1X,F7.4,2X,F25.15,2X,F25.15)
    4 FORMAT(/,2X,'Stop: The Electric Field must be 0 for NLOP' )
    5 FORMAT(/,2X,'Error: F=0 not found after sorting')
    6 FORMAT(/,2X,'Summary of tensor components (field along z)',       &
             /1X,'==============================================')
    7 FORMAT(/,' Best ALPHAxz   =',1X,ES16.8,' ±',ES10.3,               &
               '  at (k=',I2,', j=',I2,')')
    8 FORMAT(/,' Best BETAxzz   =',1X,ES16.8,' ±',ES10.3,               &
               '  at (k=',I2,', j=',I2,')')
    9 FORMAT(/,' Best GAMMAxzzz =',1X,ES16.8,' ±',ES10.3,               &
               '  at (k=',I2,', j=',I2,')')
   10 FORMAT(/,' Best ALPHAyz   =',1X,ES16.8,' ±',ES10.3,               &
               '  at (k=',I2,', j=',I2,')')
   11 FORMAT(/,' Best BETAyzz   =',1X,ES16.8,' ±',ES10.3,               &
               '  at (k=',I2,', j=',I2,')')
   12 FORMAT(/,' Best GAMMAyzzz =',1X,ES16.8,' ±',ES10.3,               &
               '  at (k=',I2,', j=',I2,')')
   13 FORMAT(/,' Best ALPHAzz   =',1X,ES16.8,' ±',ES10.3,               &
               '  at (k=',I2,', j=',I2,')')
   14 FORMAT(/,' Best BETAzzz   =',1X,ES16.8,' ±',ES10.3,               &
               '  at (k=',I2,', j=',I2,')')
   15 FORMAT(/,' Best GAMMAzzzz =',1X,ES16.8,' ±',ES10.3,               &
               '  at (k=',I2,', j=',I2,')')
   16 FORMAT(/118('-'))
!-----------------------------------------------------------------------
      DEALLOCATE(FV,DIPx,DIPy,DIPz,HPOS)
      RETURN
      END

! SORT2xyz
      SUBROUTINE SORT2xyz(A,Bx,By,Bz,M)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER M, I, J, K
      DOUBLE PRECISION A(*),Bx(*),By(*),Bz(*),TA,TBx,TBy,TBz
      DO I=1,M-1
       K = I
       DO J=I+1,M
        IF (A(J) < A(K)) K = J
       END DO
       IF (K /= I) THEN
        TA = A(I)
        TBx = Bx(I)
        TBy = By(I)
        TBz = Bz(I)
        A(I) = A(K)
        Bx(I) = Bx(K)
        By(I) = By(K)
        Bz(I) = Bz(K)
        A(K) = TA
        Bx(K) = TBx
        By(K) = TBy
        Bz(K) = TBz
       END IF
      END DO
      END SUBROUTINE SORT2xyz

!-----------------------------------------------------------------------!
!                                                                       !
!     Computes the diagonal components of the static polarizability     !
!     tensor (αxx, αyy, αzz) using the dyadic Romberg–Richardson scheme !
!     looping the field direction over x, y, z. Then reports the        !
!     isotropic average and the anisotropy (Raman convention):          !
!                                                                       !
!                  α_iso = (αxx + αyy + αzz)/3                          !
!                                                                       !
!      α_aniso = sqrt( 0.5*((αxx-αyy)^2 + (αyy-αzz)^2 + (αzz-αxx)^2) )  !
!                                                                       !
!-----------------------------------------------------------------------!

! DIPPOLARDIAG
      SUBROUTINE DIPPOLARDIAG(NP,STEP,NINTEG,IDONTW,IEMOM,NAT,NBF,      &
                              NBFaux,NSHELL,NPRIMI,ZAN,Cxyz,IAN,IMIN,   &
                              IMAX,KSTART,KATOM,KTYPE,KLOC,INTYP,KNG,   &
                              KMIN,KMAX,ISH,ITYP,C1,C2,EX,CS,CP,CD,CF,  &
                              CG,CH,CI,GRADS,IRUNTYP,DIPS,SIZE_ENV,ENV, &
                              ATM,NBAS,BAS,IGTYP,NLOP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL EFIELDL,RESTART
      COMMON/INP_EFIELDL/EFX,EFY,EFZ,EFIELDL
      COMMON/INPNOF_RSTRT/RESTART,INPUTGAMMA,INPUTC,INPUTFMIUG,INPUTCXYZ
      COMMON/INPNOF_RHF/CONVRHFDM,IRHF,IRHFTYP,NCONVRHF,MAXITRHF
      COMMON/INPNOF_PRINT/NPRINT,IWRITEC,IMULPOP,IAIMPAC,IFCHK
      COMMON/INPNOF_MOLDEN/MOLDEN
      COMMON/EHFEN/EHF,EN
      COMMON/ENERGIAS/EELEC,EELEC_OLD,DIF_EELEC,EELEC_MIN
      COMMON/INPNOF_COEFOPT/MAXLOOP
!     ARGUMENTS
      INTEGER,INTENT(IN) :: NP,NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux
      INTEGER,INTENT(IN) :: NSHELL,NPRIMI,IRUNTYP,NLOP
      DOUBLE PRECISION,INTENT(IN) :: STEP
      DOUBLE PRECISION,DIMENSION(3,NAT) :: Cxyz
      DOUBLE PRECISION,DIMENSION(NAT),INTENT(IN) :: ZAN
      INTEGER,DIMENSION(NAT),INTENT(IN) :: IAN,IMIN,IMAX
      INTEGER,DIMENSION(NSHELL),INTENT(IN) :: KSTART,KATOM,KTYPE,KLOC
      INTEGER,DIMENSION(NSHELL),INTENT(IN) :: INTYP,KNG,KMIN,KMAX
      INTEGER,DIMENSION(NPRIMI),INTENT(IN) :: ISH,ITYP
      DOUBLE PRECISION,DIMENSION(NPRIMI),INTENT(IN) :: C1,C2,EX,CS,CP
      DOUBLE PRECISION,DIMENSION(NPRIMI),INTENT(IN) :: CD,CF,CG,CH,CI
      DOUBLE PRECISION,DIMENSION(3*NAT) :: GRADS
      DOUBLE PRECISION,DIMENSION(3) :: DIPS
      INTEGER :: SIZE_ENV,ATM(6,NAT),NBAS,BAS(8,NBAS),IGTYP
      DOUBLE PRECISION :: ENV(SIZE_ENV)
!     Locals
      INTEGER :: AX, MAXM, M, I, I0, FINDX, IPrintTria
      INTEGER :: ibA, jbA, ibB, jbB, ibG, jbG
      DOUBLE PRECISION :: TOL, h, DMU(3)
      DOUBLE PRECISION :: bestA, errA, dummyB, errB, dummyG, errG
      DOUBLE PRECISION :: axx, ayy, azz
      DOUBLE PRECISION :: ALPHAD(3), ALPHA_ISO, ALPHA_ANISO
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: FV, DIPc, HPOS
!
      EXTERNAL FINDX
!-----------------------------------------------------------------------
      MAXM = 2*NP + 3
      ALLOCATE(FV(MAXM), DIPc(MAXM), HPOS(NP))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Save the original settings
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      TOL = 1.0D-16
      DMU  = 0.0D0
      IF(EFX==0.0d0 .and. EFY==0.0d0 .and. EFZ==0.0d0)THEN
       CALL OriginalOptions(0)
      ELSE
       WRITE(6,4)
       STOP
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Dyadic steps:  h, 2h, 4h, …  (ascending powers of two)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO I = 1, NP
       HPOS(I) = STEP * 2.0D0**(I-1)
      END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Build dyadic ladder and evaluate DIPc at each required field
!     Loop over field axis: 1=x, 2=y, 3=z
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO AX = 1,3
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Energy without Electric Field
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       CALL OriginalOptions(2)   ! Restore the original settings
       if(AX==1)then
        CALL ENERGRAD(NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NSHELL,NPRIMI, &
                      ZAN,Cxyz,IAN,IMIN,IMAX,KSTART,KATOM,KTYPE,KLOC,   &
                      INTYP,KNG,KMIN,KMAX,ISH,ITYP,C1,C2,EX,CS,CP,CD,CF,&
                      CG,CH,CI,GRADS,IRUNTYP,DIPS,SIZE_ENV,ENV,ATM,NBAS,&
                      BAS,IGTYP,0,1)
        WRITE(6,1)
        WRITE(6,21)
       else if(AX==2)then
        CALL ENERGRAD(NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NSHELL,NPRIMI, &
                      ZAN,Cxyz,IAN,IMIN,IMAX,KSTART,KATOM,KTYPE,KLOC,   &
                      INTYP,KNG,KMIN,KMAX,ISH,ITYP,C1,C2,EX,CS,CP,CD,CF,&
                      CG,CH,CI,GRADS,IRUNTYP,DIPS,SIZE_ENV,ENV,ATM,NBAS,&
                      BAS,IGTYP,0,0)
        WRITE(6,22)
       else if(AX==3)then
        CALL ENERGRAD(NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NSHELL,NPRIMI, &
                      ZAN,Cxyz,IAN,IMIN,IMAX,KSTART,KATOM,KTYPE,KLOC,   &
                      INTYP,KNG,KMIN,KMAX,ISH,ITYP,C1,C2,EX,CS,CP,CD,CF,&
                      CG,CH,CI,GRADS,IRUNTYP,DIPS,SIZE_ENV,ENV,ATM,NBAS,&
                      BAS,IGTYP,0,0)
        WRITE(6,23)
       end if
!      Ensure F=0 (M=1) present
       M = 1
       FV(M) = 0.0d0
!      Take dipole component aligned with field axis
       DIPc(M) = DIPS(AX)
       WRITE(6,3)FV(M),DIPc(M),EELEC+EN
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Energy and Dipole for I=1,NP+1 [ EF(ax) > 0 ]
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       CALL OriginalOptions(1)   ! Impose mandatory options
       DO I = 1, NP+1
!       Electric Field
        EFX = 0.0D0
        EFY = 0.0D0
        EFZ = 0.0D0
        h = STEP * 2.0D0**(I-1)
        if(AX==1)then
         EFX = h
        else if(AX==2)then
         EFY = h
        else if(AX==3)then
         EFZ = h
        end if
!
        M = M + 1
        MAXLOOP = 10
        CALL ENERGRAD(NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NSHELL,        &
                      NPRIMI,ZAN,Cxyz,IAN,IMIN,IMAX,KSTART,KATOM,       &
                      KTYPE,KLOC,INTYP,KNG,KMIN,KMAX,ISH,ITYP,C1,C2,EX, &
                      CS,CP,CD,CF,CG,CH,CI,GRADS,IRUNTYP,DMU,SIZE_ENV,  &
                      ENV,ATM,NBAS,BAS,IGTYP,0,0)
        FV(M) = h
        DIPc(M) = DMU(AX)
        WRITE(6,3) FV(M), DIPc(M), EELEC+EN
       END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Recalculate the Energy without Electric Field
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       CALL OriginalOptions(2)   ! Restore the original settings
       CALL ENERGRAD(NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NSHELL,NPRIMI,  &
                     ZAN,Cxyz,IAN,IMIN,IMAX,KSTART,KATOM,KTYPE,KLOC,    &
                     INTYP,KNG,KMIN,KMAX,ISH,ITYP,C1,C2,EX,CS,CP,CD,CF, &
                     CG,CH,CI,GRADS,IRUNTYP,DMU,SIZE_ENV,ENV,ATM,NBAS,  &
                     BAS,IGTYP,0,0)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Energy and Dipole for I=1,NP+1 [ EF(ax) < 0 ]
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       CALL OriginalOptions(1)   ! Impose mandatory options
       DO I = 1, NP+1
!       Electric Field
        EFX = 0.0D0
        EFY = 0.0D0
        EFZ = 0.0D0
        h = STEP * 2.0D0**(I-1)
        if(AX==1)then
         EFX = -h
        else if(AX==2)then
         EFY = -h
        else if(AX==3)then
         EFZ = -h
        end if
!
        M = M + 1
        MAXLOOP = 10
        CALL ENERGRAD(NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NSHELL,        &
                      NPRIMI,ZAN,Cxyz,IAN,IMIN,IMAX,KSTART,KATOM,       &
                      KTYPE,KLOC,INTYP,KNG,KMIN,KMAX,ISH,ITYP,C1,C2,EX, &
                      CS,CP,CD,CF,CG,CH,CI,GRADS,IRUNTYP,DMU,SIZE_ENV,  &
                      ENV,ATM,NBAS,BAS,IGTYP,0,0)
        FV(M) = -h
        DIPc(M) = DMU(AX)
        WRITE(6,3) FV(M), DIPc(M), EELEC+EN
       END DO
       CALL OriginalOptions(2)   ! Restore the original settings
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Sort FV in ascending order
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       CALL SORT2z(FV, DIPc, M)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Locate index of F=0 (after sorting)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       I0 = FINDX(FV, M, 0.0D0, TOL)
       IF (I0 <= 0) THEN
        WRITE(6,5)
        STOP
       END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Build and process full Romberg–Richardson triangles for the
!      first derivatives (α) of the dipole μ(F) w.r.t. electric field F
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       WRITE(6,7)
       CALL ROMBERG_TRIANGLE(FV, DIPc, M, I0, NP, TOL, HPOS,            &
                             bestA, errA, ibA, jbA,                     &
                             dummyB, errB, ibB, jbB,                    &
                             dummyG, errG, ibG, jbG, NLOP, 1)
       if (AX==1) then
        WRITE(6,61) bestA, errA, ibA, jbA
       else if (AX==2) then
        WRITE(6,62) bestA, errA, ibA, jbA
       else if (AX==3) then
        WRITE(6,63) bestA, errA, ibA, jbA
       end if
       WRITE(6,7)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Store diagonal α for this axis
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       if (AX==1) then
        axx = bestA
       else if (AX==2) then
        ayy = bestA
       else if (AX==3) then
        azz = bestA
       end if
      END DO  ! AX = 1..3
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Isotropic average and the anisotropy (Raman convention)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALPHAD(1) = axx
      ALPHAD(2) = ayy
      ALPHAD(3) = azz
      ALPHA_ISO = (axx + ayy + azz) / 3.0D0
      ALPHA_ANISO = DSQRT( 0.5D0 * ( (axx-ayy)*(axx-ayy)                &
                     + (ayy-azz)*(ayy-azz) + (azz-axx)*(azz-axx) ) )

      WRITE(6,8) axx, ayy, azz
      WRITE(6,9) ALPHA_ISO
      WRITE(6,10)ALPHA_ANISO
      WRITE(6,7)
!-----------------------------------------------------------------------
!     Format definitions
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    1 FORMAT(//2X,'Non Linear Optical Properties (Romberg–Richardson)', &
            /1X,'====================================================')
   21 FORMAT(/,4X,'EFx',16X,'Dipole_x',17X,'Energy',/)
   22 FORMAT(/,4X,'EFy',16X,'Dipole_y',17X,'Energy',/)
   23 FORMAT(/,4X,'EFz',16X,'Dipole_z',17X,'Energy',/)
    3 FORMAT(1X,F7.4,2X,F25.15,2X,F25.15)
    4 FORMAT(/,2X,'Stop: The Electric Field must be 0 for NLOP' )
    5 FORMAT(/,2X,'Error: F=0 not found after sorting')
   61 FORMAT(/,' Best αxx   =',1X,ES16.8,' ±',ES10.3,                   &
               '  at (k=',I2,', j=',I2,')')
   62 FORMAT(/,' Best αyy   =',1X,ES16.8,' ±',ES10.3,                   &
               '  at (k=',I2,', j=',I2,')')
   63 FORMAT(/,' Best αzz   =',1X,ES16.8,' ±',ES10.3,                   &
               '  at (k=',I2,', j=',I2,')')
    7 FORMAT(/118('-'))
    8 FORMAT(/,' Alpha diagonal:  αxx=',ES16.8,'   αyy=',ES16.8,        &
                                               '   αzz=',ES16.8)
    9 FORMAT(/,' Alpha (isotropic average)  = ',ES16.8)
   10 FORMAT(/,' Alpha (anisotropy, Δα)     = ',ES16.8)
!-----------------------------------------------------------------------
      DEALLOCATE(FV,DIPc,HPOS)
      RETURN
      END
