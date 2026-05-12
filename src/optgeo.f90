!======================================================================!
!                                                                      !
!  G E O M E T R Y   O P T I M I Z A T I O N   S U B R O U T I N E S   !
!                                                                      !
!======================================================================!
!                                                                      !
!   OPTIMIZE: Routine to carry out geometry optimization (IRUNTYP=3)   !
!   OPTSUMSL: Optimize geometries using the CG SUMSL routine           !
!   POINTERSOPT : Define Pointers of the USER array                    !
!   CALCOPTE: Routine to compute energy needed by OPTSUMSL             !
!   CALCOPTG: Routine to compute gradients needed by OPTSUMSL          !
!                                                                      !
!   OPTCGNAG: Optimize geometries using the CG NAG routine E04DGF      !
!   ENERGYFUN: Routine to compute energy needed by OPTCGNAG            !
!                                                                      !
!======================================================================!

! OPTIMIZE
      SUBROUTINE OPTIMIZE(NINTEG,IDONTW,NAT,ZAN,Cxyz,IAN,IMIN,IMAX,     &
                          ZMASS,KSTART,KATOM,KTYPE,KLOC,INTYP,KNG,KMIN, &
                          KMAX,ISH,ITYP,C1,C2,EX1,CS,CP,CD,CF,CG,CH,CI, &
                          DIPS,IZCORE,SIZE_ENV,ENV,ATM,NBAS,BAS,IGTYP,  &
                          GRADS,IRUNTYP,IHSSCAL,IPROJECT,ISIGMA,KTYPEaux)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL RESTART,WRTHSSNOTE,WRTHSSCGO,WRXYZHDR
      INTEGER :: NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_RSTRT/RESTART,INPUTGAMMA,INPUTC,INPUTFMIUG,INPUTCXYZ
      COMMON/INPNOF_CGM/ICGMETHOD
      COMMON/INPNOF_PRINT/NPRINT,IWRITEC,IMULPOP,IAIMPAC,IFCHK
      COMMON/INPNOF_GENERALINF/ICOEF,ISOFTMAX,IORBOPT,MAXIT,MAXIT21
      COMMON/INPFILE_Naux/NBFaux,NSHELLaux
      COMMON/ENERGIAS/EELEC,EELEC_OLD,DIF_EELEC,EELEC_MIN
      COMMON/ELPROP/IEMOM
      COMMON/GEOSTAT/GEOMAXG,GEORMSG,IOPTIRC
!
      DOUBLE PRECISION,DIMENSION(NAT) :: ZAN,ZMASS
      INTEGER,DIMENSION(NAT) :: IAN,IMIN,IMAX,IZCORE
      INTEGER,DIMENSION(NSHELL) :: KSTART,KATOM,KTYPE,KLOC
      INTEGER,DIMENSION(NSHELL) :: INTYP,KNG,KMIN,KMAX
      INTEGER,DIMENSION(NSHELLaux) :: KTYPEaux
      INTEGER,DIMENSION(NPRIMI) :: ISH,ITYP
      DOUBLE PRECISION,DIMENSION(NPRIMI)::C1,C2,EX1,CS,CP,CD,CF,CG,CH,CI
      DOUBLE PRECISION,DIMENSION(3,NAT):: Cxyz
      DOUBLE PRECISION,DIMENSION(3,NAT):: CTRY
      DOUBLE PRECISION,DIMENSION(3*NAT):: GRADS
      DOUBLE PRECISION,DIMENSION(3):: DIPS

      INTEGER :: SIZE_ENV,NBAS,IGTYP
      DOUBLE PRECISION :: ENV(SIZE_ENV)
      INTEGER :: ATM(6,NAT), BAS(8,NBAS)
      DOUBLE PRECISION,PARAMETER :: GMAXTOL=1.0D-3
!-----------------------------------------------------------------------
      EELEC_MIN = 1.0d20
      IOPTIRC = 0
      GEOMAXG = 0.0D0
      GEORMSG = 0.0D0
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Generate an initial GCF if RESTART=F setting ICOEF=0
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(RESTART .eqv. .FALSE.)THEN
       ICOEFORI = ICOEF
       ICOEF = 0
!      Update coordinates of shells if use libcint library for ERIs
       CALL ENERGRAD(NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NSHELL,NPRIMI,  &
                     ZAN,Cxyz,IAN,IMIN,IMAX,KSTART,KATOM,KTYPE,KLOC,    &
                     INTYP,KNG,KMIN,KMAX,ISH,ITYP,C1,C2,EX1,CS,CP,CD,   &
                     CF,CG,CH,CI,GRADS,IRUNTYP,DIPS,IZCORE,SIZE_ENV,    &
                     ENV,ATM,NBAS,BAS,IGTYP,KTYPEaux,NSHELLaux,0,0)
       RESTART = .TRUE.
       ICOEF = ICOEFORI
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Select CG Method
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(ICGMETHOD==1)THEN
       CALL OPTSUMSL(NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NSHELL,NPRIMI,  &
                     ZAN,Cxyz,IAN,IMIN,IMAX,ZMASS,KSTART,KATOM,KTYPE,   &
                     KLOC,INTYP,KNG,KMIN,KMAX,ISH,ITYP,C1,C2,EX1,CS,CP, &
                     CD,CF,CG,CH,CI,IRUNTYP,SIZE_ENV,ENV,ATM,           &
                     NBAS,BAS,IGTYP,GRADS,DIPS,IZCORE,KTYPEaux,NSHELLaux)
      ELSE IF(ICGMETHOD==2)THEN                                         
       CALL OPTCGNAG(NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NSHELL,NPRIMI,  &
                     ZAN,Cxyz,IAN,IMIN,IMAX,ZMASS,KSTART,KATOM,KTYPE,   &
                     KLOC,INTYP,KNG,KMIN,KMAX,ISH,ITYP,C1,C2,EX1,CS,    &
                     CP,CD,CF,CG,CH,CI,IRUNTYP,GRADS,DIPS,IZCORE,       &
                     KTYPEaux,NSHELLaux)
      ELSE IF(ICGMETHOD==3)THEN
       WRITE(6,*)'Stop: To be implemented'
       Stop
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Compute Hessian from analytic gradients at stationary point
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(IHSSCAL==1)THEN
       IF(ICGMETHOD==1 .AND. GEOMAXG>GMAXTOL)THEN
        WRITE(6,'(/,1X,A,1PE12.4)')                                     &
        'Skipping HESSCAL: final max|Grad| is too large = ',GEOMAXG
        WRITE(6,'( 1X,A)')                                              &
        'Restart OPTGEO from the final GCF or run a GRAD check first.'
       ELSE
        CALL HESSCAL(NINTEG,IDONTW,NAT,ZAN,Cxyz,IAN,IMIN,IMAX,ZMASS,    &
                     KSTART,KATOM,KTYPE,KLOC,INTYP,KNG,KMIN,KMAX,ISH,   &
                     ITYP,C1,C2,EX1,CS,CP,CD,CF,CG,CH,CI,DIPS,          &
                     GRADS,IRUNTYP,IPROJECT,ISIGMA,SIZE_ENV,ENV,ATM,    &
                     NBAS,BAS,IGTYP,KTYPEaux,IZCORE)
       END IF
      END IF
!-----------------------------------------------------------------------
      RETURN
      END

! OPTIMIZETS
      SUBROUTINE OPTIMIZETS(NINTEG,IDONTW,NAT,ZAN,Cxyz,IAN,IMIN,IMAX,   &
                          ZMASS,KSTART,KATOM,KTYPE,KLOC,INTYP,KNG,KMIN, &
                          KMAX,ISH,ITYP,C1,C2,EX1,CS,CP,CD,CF,CG,CH,CI, &
                          DIPS,IZCORE,SIZE_ENV,ENV,ATM,NBAS,BAS,IGTYP,  &
                          GRADS,IRUNTYP,IHSSCAL,IPROJECT,ISIGMA,        &
                          KTYPEaux)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL RESTART,WRTHSSNOTE,WRTHSSCGO,WRXYZHDR
      INTEGER :: NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPFILE_Naux/NBFaux,NSHELLaux
      COMMON/INPNOF_RSTRT/RESTART,INPUTGAMMA,INPUTC,INPUTFMIUG,INPUTCXYZ
      COMMON/INPNOF_MOLDEN/MOLDEN
      COMMON/INPNOF_GENERALINF/ICOEF,ISOFTMAX,IORBOPT,MAXIT,MAXIT21
      COMMON/WRTGCF/IWRTGCF
      COMMON/EHFEN/EHF,EN
      COMMON/ENERGIAS/EELEC,EELEC_OLD,DIF_EELEC,EELEC_MIN
      COMMON/ELPROP/IEMOM
      COMMON/GEOSTAT/GEOMAXG,GEORMSG,IOPTIRC
      COMMON/HSSNOTE/WRTHSSNOTE
      COMMON/HSSCGO/WRTHSSCGO
      COMMON/XYZHDR/WRXYZHDR
!
      INTEGER,DIMENSION(NAT) :: IAN,IMIN,IMAX,IZCORE
      INTEGER,DIMENSION(NSHELL) :: KSTART,KATOM,KTYPE,KLOC
      INTEGER,DIMENSION(NSHELL) :: INTYP,KNG,KMIN,KMAX
      INTEGER,DIMENSION(NSHELLaux) :: KTYPEaux
      INTEGER,DIMENSION(NPRIMI) :: ISH,ITYP
      DOUBLE PRECISION,DIMENSION(NAT) :: ZAN,ZMASS
      DOUBLE PRECISION,DIMENSION(NPRIMI)::C1,C2,EX1,CS,CP,CD,CF,CG,CH,CI
      DOUBLE PRECISION,DIMENSION(3,NAT):: Cxyz
      DOUBLE PRECISION,DIMENSION(3,NAT):: CTRY
      DOUBLE PRECISION,DIMENSION(3*NAT):: GRADS
      DOUBLE PRECISION,DIMENSION(3):: DIPS
      INTEGER :: SIZE_ENV,NBAS,IGTYP
      DOUBLE PRECISION :: ENV(SIZE_ENV)
      INTEGER :: ATM(6,NAT), BAS(8,NBAS)
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: HESSIAN,VEC,DDM
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: EIG,WORK,STEPV,GPROJ
      DOUBLE PRECISION,PARAMETER :: BOHR = 0.52917724924D+00
      DOUBLE PRECISION,PARAMETER :: TSGTOL = 1.0D-3
      DOUBLE PRECISION,PARAMETER :: TSSTEPMAX = 0.15D0
      DOUBLE PRECISION,PARAMETER :: TSSTEPASC = 0.08D0
      DOUBLE PRECISION,PARAMETER :: TSSTEPDESC = 0.12D0
      DOUBLE PRECISION,PARAMETER :: TSFINEG = 2.0D-3
      DOUBLE PRECISION,PARAMETER :: TSSTEPMAXF = 0.05D0
      DOUBLE PRECISION,PARAMETER :: TSSTEPASCF = 0.03D0
      DOUBLE PRECISION,PARAMETER :: TSSTEPDESCF = 0.05D0
      DOUBLE PRECISION,PARAMETER :: TSCURVMIN = 1.0D-4
      DOUBLE PRECISION,PARAMETER :: TSTRYSHRINK = 0.5D0
      DOUBLE PRECISION,PARAMETER :: TSGRFACT = 1.10D0
      DOUBLE PRECISION,PARAMETER :: TSGRRELAX = 1.30D0
      DOUBLE PRECISION,PARAMETER :: TSLAMFACT = 1.50D0
      DOUBLE PRECISION,PARAMETER :: TSENERTOL = 5.0D-4
      DOUBLE PRECISION,PARAMETER :: TSRESCSTEP = 0.03D0
      DOUBLE PRECISION,PARAMETER :: TSRADSHRINK = 0.5D0
      DOUBLE PRECISION,PARAMETER :: TSMINSTEP = 5.0D-3
      INTEGER,PARAMETER :: TSTRYMAX = 6
      INTEGER,PARAMETER :: TSRADMAX = 4
      INTEGER,PARAMETER :: TSMAXIT = 40
      INTEGER :: I, J, NV, ITER, IMODE, NNEG, ICOEFORI, IZNUC
      INTEGER :: ITRY, IMTRY, NNEGTRY, IRADTRY
      DOUBLE PRECISION :: ENERGIA, STEP2, SCALE, CURV, COEF, STEPMAXI
      DOUBLE PRECISION :: DENOMASC, MUORTH, ORTH2, BASEI
      DOUBLE PRECISION :: MUH, MULOW, GMODE
      DOUBLE PRECISION :: SHRFAC, GEOMAXGTRY, ENERTRY, GTRYRAT
      DOUBLE PRECISION :: LAMREF, SHRRESC, STEPLEN, STEPASCI,STEPDESCI
      LOGICAL :: TSFAIL, ACCEPTSTEP
      INTEGER :: LUGEOCONV, MOLDEN, MOLDENSAVE
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
!-----------------------------------------------------------------------
      EELEC_MIN = 1.0D20
      IOPTIRC = 0
      GEOMAXG = 0.0D0
      GEORMSG = 0.0D0
      TSFAIL = .FALSE.
      NV = 3*NAT
      IF(RESTART .eqv. .FALSE.)THEN
       ICOEFORI = ICOEF
       ICOEF = 0
       CALL ENERGRAD(NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NSHELL,NPRIMI,  &
                     ZAN,Cxyz,IAN,IMIN,IMAX,KSTART,KATOM,KTYPE,KLOC,    &
                     INTYP,KNG,KMIN,KMAX,ISH,ITYP,C1,C2,EX1,CS,CP,CD,   &
                     CF,CG,CH,CI,GRADS,2,DIPS,IZCORE,SIZE_ENV,ENV,ATM,  &
                     NBAS,BAS,IGTYP,KTYPEaux,NSHELLaux,0,0)
       RESTART = .TRUE.
       ICOEF = ICOEFORI
      END IF
      ALLOCATE(HESSIAN(NV,NV),VEC(NV,NV),DDM(9,NAT))
      ALLOCATE(EIG(NV),WORK(NV),STEPV(NV),GPROJ(NV))
      OPEN(NEWUNIT=LUGEOCONV,STATUS='SCRATCH',FORM='FORMATTED')
      WRITE(6,'(/,1X,A29,/1X,29(1H=))')'TRANSITION-STATE OPTIMIZATION'
      WRITE(6,'(/,1X,A7,8X,A6,11X,A9,7X,A4,6X,A10,/)')                  &
              'TS Iter','Energy','max|Grad|','Nneg','lambda_min'
      WRXYZHDR = .FALSE.
      IF(MOLDEN==1)THEN
       WRITE(18,'(A18)')'[GEOMETRIES] (XYZ)'
       WRXYZHDR = .TRUE.
      END IF
      WRITE(11,'(/,72(1H*),/)')
      WRITE(11,*)'Initial Coordinates (Angs) for TS Optimization'
      WRITE(11,*)
      DO I=1,NAT
       WRITE(11,'(I5,3F15.4)')I,Cxyz(1,I)*BOHR,Cxyz(2,I)*BOHR,          &
                               Cxyz(3,I)*BOHR
      ENDDO
      WRITE(11,'(/,72(1H*),/)')
      DO ITER = 1,TSMAXIT
       DO I=1,NAT
        ENV(20+3*(I-1)+1) = Cxyz(1,I)
        ENV(20+3*(I-1)+2) = Cxyz(2,I)
        ENV(20+3*(I-1)+3) = Cxyz(3,I)
       ENDDO
       CALL ENERGRAD(NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NSHELL,NPRIMI,  &
                     ZAN,Cxyz,IAN,IMIN,IMAX,KSTART,KATOM,KTYPE,KLOC,    &
                     INTYP,KNG,KMIN,KMAX,ISH,ITYP,C1,C2,EX1,CS,CP,CD,   &
                     CF,CG,CH,CI,GRADS,2,DIPS,IZCORE,SIZE_ENV,ENV,ATM,  &
                     NBAS,BAS,IGTYP,KTYPEaux,NSHELLaux,0,0)
       ENERGIA = EELEC + EN
       GEOMAXG = 0.0D0
       GEORMSG = 0.0D0
       DO I=1,NV
        GEOMAXG = DMAX1(GEOMAXG,DABS(GRADS(I)))
        GEORMSG = GEORMSG + GRADS(I)*GRADS(I)
       ENDDO
       GEORMSG = DSQRT(GEORMSG/DBLE(NV))
       IWRTGCF = 0
       WRTHSSNOTE = .FALSE.
       WRTHSSCGO = .FALSE.
       MOLDENSAVE = MOLDEN
       MOLDEN = 0
       CALL HSSNUMd(HESSIAN,NV,Cxyz,GRADS,DIPS,DDM,NINTEG,IDONTW,ZAN,   &
                    IAN,IMIN,IMAX,KSTART,KATOM,KTYPE,KLOC,INTYP,KNG,    &
                    KMIN,KMAX,ISH,ITYP,C1,C2,EX1,CS,CP,CD,CF,CG,CH,CI,  &
                    2,SIZE_ENV,ENV,ATM,NBAS,BAS,IGTYP,KTYPEaux,IZCORE)
       MOLDEN = MOLDENSAVE
       CALL TSPROJMODES(NAT,NV,Cxyz,HESSIAN,GRADS,ZAN,ZMASS,IPROJECT,   &
                        VEC,EIG,IMODE,NNEG,IZCORE)
       WRITE(6,'(1X,I4,1X,F20.10,6X,1PE10.3,6X,I3,8X,1PE10.3)')         &
                   ITER,ENERGIA,GEOMAXG,NNEG,EIG(IMODE)
       WRITE(LUGEOCONV,'(F20.10)')ENERGIA
       IF(MOLDEN==1)THEN
        WRITE(18,'(I6)')NAT
        WRITE(18,'(I6,F20.10)')ITER,ENERGIA
        DO J=1,NAT
         IZNUC = INT(ZAN(J)) + IZCORE(J)
         WRITE(18,'(1X,A4,3F15.4)')ATMLAB(IZNUC),                       &
              Cxyz(1,J)*BOHR,Cxyz(2,J)*BOHR,Cxyz(3,J)*BOHR
        ENDDO
       END IF
       IF(GEOMAXG<TSGTOL .AND. NNEG==1) EXIT
       STEPMAXI = TSSTEPMAX
       STEPASCI = TSSTEPASC
       STEPDESCI = TSSTEPDESC
       IF(GEOMAXG<TSFINEG)THEN
        STEPMAXI = DMIN1(STEPMAXI,TSSTEPMAXF)
        STEPASCI = DMIN1(STEPASCI,TSSTEPASCF)
        STEPDESCI = DMIN1(STEPDESCI,TSSTEPDESCF)
       END IF
       DO I=1,NV
        GPROJ(I) = DDOT(NV,VEC(1,I),1,GRADS,1)
       ENDDO
       STEPV = 0.0D0
       GMODE = GPROJ(IMODE)
       DENOMASC = -DMAX1(DABS(EIG(IMODE)),DABS(GMODE)/STEPASCI,        &
                         TSCURVMIN)
       COEF = -GMODE/DENOMASC
       CALL DAXPY(NV,COEF,VEC(1,IMODE),1,STEPV,1)
       MUORTH = 0.0D0
       ORTH2 = 0.0D0
       DO I=1,NV
        IF(I==IMODE) CYCLE
        BASEI = DMAX1(EIG(I),TSCURVMIN)
        ORTH2 = ORTH2 + (GPROJ(I)/BASEI)**2
       ENDDO
       IF(DSQRT(ORTH2)>STEPDESCI)THEN
        MULOW = 0.0D0
        MUH = 1.0D0
  110   ORTH2 = 0.0D0
        DO I=1,NV
         IF(I==IMODE) CYCLE
         BASEI = DMAX1(EIG(I),TSCURVMIN)
         ORTH2 = ORTH2 + (GPROJ(I)/(BASEI+MUH))**2
        ENDDO
        IF(DSQRT(ORTH2)>STEPDESCI)THEN
         MUH = 2.0D0*MUH
         GO TO 110
        END IF
        DO I=1,40
         MUORTH = 0.5D0*(MULOW+MUH)
         ORTH2 = 0.0D0
         DO J=1,NV
          IF(J==IMODE) CYCLE
          BASEI = DMAX1(EIG(J),TSCURVMIN)
          ORTH2 = ORTH2 + (GPROJ(J)/(BASEI+MUORTH))**2
         ENDDO
         IF(DSQRT(ORTH2)>STEPDESCI)THEN
          MULOW = MUORTH
         ELSE
          MUH = MUORTH
         END IF
        ENDDO
       END IF
       DO I=1,NV
        IF(I==IMODE) CYCLE
        BASEI = DMAX1(EIG(I),TSCURVMIN)
        COEF = -GPROJ(I)/(BASEI+MUORTH)
        CALL DAXPY(NV,COEF,VEC(1,I),1,STEPV,1)
       ENDDO
       STEP2 = DDOT(NV,STEPV,1,STEPV,1)
       IF(STEP2>STEPMAXI*STEPMAXI)THEN
       SCALE = STEPMAXI/DSQRT(STEP2)
        CALL DSCAL(NV,SCALE,STEPV,1)
       END IF
  120  CONTINUE
       IRADTRY = 1
  130  CONTINUE
       SHRFAC = 1.0D0
       ACCEPTSTEP = .FALSE.
       DO ITRY = 1,TSTRYMAX
        DO I=1,NAT
         CTRY(1,I) = Cxyz(1,I) + SHRFAC*STEPV(3*(I-1)+1)
         CTRY(2,I) = Cxyz(2,I) + SHRFAC*STEPV(3*(I-1)+2)
         CTRY(3,I) = Cxyz(3,I) + SHRFAC*STEPV(3*(I-1)+3)
         ENV(20+3*(I-1)+1) = CTRY(1,I)
         ENV(20+3*(I-1)+2) = CTRY(2,I)
         ENV(20+3*(I-1)+3) = CTRY(3,I)
        ENDDO
        CALL ENERGRAD(NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NSHELL,NPRIMI, &
                      ZAN,CTRY,IAN,IMIN,IMAX,KSTART,KATOM,KTYPE,KLOC,   &
                      INTYP,KNG,KMIN,KMAX,ISH,ITYP,C1,C2,EX1,CS,CP,CD,  &
                      CF,CG,CH,CI,GRADS,2,DIPS,IZCORE,SIZE_ENV,ENV,ATM, &
                      NBAS,BAS,IGTYP,KTYPEaux,NSHELLaux,0,0)
        WRTHSSNOTE = .FALSE.
        WRTHSSCGO = .FALSE.
        MOLDENSAVE = MOLDEN
        MOLDEN = 0
        CALL HSSNUMd(HESSIAN,NV,CTRY,GRADS,DIPS,DDM,NINTEG,IDONTW,ZAN,  &
                     IAN,IMIN,IMAX,KSTART,KATOM,KTYPE,KLOC,INTYP,KNG,   &
                     KMIN,KMAX,ISH,ITYP,C1,C2,EX1,CS,CP,CD,CF,CG,CH,CI, &
                     2,SIZE_ENV,ENV,ATM,NBAS,BAS,IGTYP,KTYPEaux,IZCORE)
        MOLDEN = MOLDENSAVE
        CALL TSPROJMODES(NAT,NV,CTRY,HESSIAN,GRADS,ZAN,ZMASS,IPROJECT,  &
                         VEC,EIG,IMTRY,NNEGTRY,IZCORE)
        ENERTRY = EELEC + EN
        GEOMAXGTRY = 0.0D0
        DO J=1,NV
         GEOMAXGTRY = DMAX1(GEOMAXGTRY,DABS(GRADS(J)))
        ENDDO
        GTRYRAT = GEOMAXGTRY/DMAX1(GEOMAXG,1.0D-12)
        LAMREF = DABS(EIG(IMODE))
        IF(NNEGTRY==1 .AND.                                             &
           (GEOMAXGTRY<=TSGRFACT*GEOMAXG .OR.                           &
           (GTRYRAT<=TSGRRELAX .AND. DABS(ENERTRY-ENERGIA)<=TSENERTOL   &
            .AND. DABS(EIG(IMTRY))<=TSLAMFACT*DMAX1(LAMREF,TSCURVMIN))) &
          )THEN
         ACCEPTSTEP = .TRUE.
         IMODE = IMTRY
         NNEG = NNEGTRY
         GEOMAXG = GEOMAXGTRY
         Cxyz = CTRY
         EXIT
        END IF
        SHRFAC = SHRFAC*TSTRYSHRINK
      ENDDO
      IF(.NOT.ACCEPTSTEP)THEN
        SHRRESC = TSRESCSTEP/DMAX1(DABS(COEF),TSRESCSTEP)
        SHRRESC = DMIN1(1.0D0,SHRRESC)
        DO I=1,NAT
         CTRY(1,I) = Cxyz(1,I) + SHRRESC*COEF*VEC(3*(I-1)+1,IMODE)
         CTRY(2,I) = Cxyz(2,I) + SHRRESC*COEF*VEC(3*(I-1)+2,IMODE)
         CTRY(3,I) = Cxyz(3,I) + SHRRESC*COEF*VEC(3*(I-1)+3,IMODE)
         ENV(20+3*(I-1)+1) = CTRY(1,I)
         ENV(20+3*(I-1)+2) = CTRY(2,I)
         ENV(20+3*(I-1)+3) = CTRY(3,I)
        ENDDO
        CALL ENERGRAD(NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NSHELL,NPRIMI, &
                      ZAN,CTRY,IAN,IMIN,IMAX,KSTART,KATOM,KTYPE,KLOC,   &
                      INTYP,KNG,KMIN,KMAX,ISH,ITYP,C1,C2,EX1,CS,CP,CD,  &
                      CF,CG,CH,CI,GRADS,2,DIPS,IZCORE,SIZE_ENV,ENV,ATM, &
                      NBAS,BAS,IGTYP,KTYPEaux,NSHELLaux,0,0)
        WRTHSSNOTE = .FALSE.
        WRTHSSCGO = .FALSE.
        MOLDENSAVE = MOLDEN
        MOLDEN = 0
        CALL HSSNUMd(HESSIAN,NV,CTRY,GRADS,DIPS,DDM,NINTEG,IDONTW,ZAN,  &
                     IAN,IMIN,IMAX,KSTART,KATOM,KTYPE,KLOC,INTYP,KNG,   &
                     KMIN,KMAX,ISH,ITYP,C1,C2,EX1,CS,CP,CD,CF,CG,CH,CI, &
                     2,SIZE_ENV,ENV,ATM,NBAS,BAS,IGTYP,KTYPEaux,IZCORE)
        MOLDEN = MOLDENSAVE
        CALL TSPROJMODES(NAT,NV,CTRY,HESSIAN,GRADS,ZAN,ZMASS,IPROJECT,  &
                         VEC,EIG,IMTRY,NNEGTRY,IZCORE)
        GEOMAXGTRY = 0.0D0
        DO J=1,NV
         GEOMAXGTRY = DMAX1(GEOMAXGTRY,DABS(GRADS(J)))
        ENDDO
        IF(NNEGTRY==1 .AND. GEOMAXGTRY<=TSGRRELAX*GEOMAXG)THEN
         ACCEPTSTEP = .TRUE.
         IMODE = IMTRY
         NNEG = NNEGTRY
         GEOMAXG = GEOMAXGTRY
        Cxyz = CTRY
       ELSE
         STEP2 = DDOT(NV,STEPV,1,STEPV,1)
         STEPLEN = DSQRT(DMAX1(STEP2,0.0D0))
         IF(IRADTRY<TSRADMAX .AND. STEPLEN>TSMINSTEP)THEN
          IRADTRY = IRADTRY + 1
          CALL DSCAL(NV,TSRADSHRINK,STEPV,1)
          GO TO 130
         ELSE
          TSFAIL = .TRUE.
          EXIT
         END IF
        END IF
      END IF
      IF(ACCEPTSTEP)THEN
       DO I=1,NAT
        ENV(20+3*(I-1)+1) = Cxyz(1,I)
        ENV(20+3*(I-1)+2) = Cxyz(2,I)
        ENV(20+3*(I-1)+3) = Cxyz(3,I)
       ENDDO
!       Write the accepted TSOPT point to GCFe using the converged
!       electronic solution of the accepted geometry.
       IWRTGCF = 1
       CALL ENERGRAD(NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NSHELL,NPRIMI,  &
                     ZAN,Cxyz,IAN,IMIN,IMAX,KSTART,KATOM,KTYPE,KLOC,    &
                     INTYP,KNG,KMIN,KMAX,ISH,ITYP,C1,C2,EX1,CS,CP,CD,   &
                     CF,CG,CH,CI,GRADS,2,DIPS,IZCORE,SIZE_ENV,ENV,ATM,  &
                     NBAS,BAS,IGTYP,KTYPEaux,NSHELLaux,0,0)
       IWRTGCF = 0
      END IF
      ENDDO
      IOPTIRC = 0
      CALL COPYGCFE2GCFl
      WRITE(6,'(/1X,A24,2X,F20.10,/)')'Final objective value = ',ENERGIA
      WRITE(6,*)'New Coordinates after TS Optimization (Bohr)'
      WRITE(6,'(1X,51(1H-),/)')
      DO I=1,NAT
       WRITE(6,'(I5,3F15.4)')I,Cxyz(1,I),Cxyz(2,I),Cxyz(3,I)
      ENDDO
      IF(TSFAIL)THEN
       WRITE(6,'(/1X,A,A)')'TSOPT stopped: unable to keep one projected',&
                           ' imaginary mode after step backtracking'
      ELSE IF(GEOMAXG<TSGTOL .AND. NNEG==1)THEN
       WRITE(6,'(/1X,A,I4,2X,A,1PE10.3,2X,A,I3)')                       &
            'TSOPT converged in',ITER,'max|Grad| =',GEOMAXG,            &
            'Nneg =',NNEG
       WRITE(11,'(/,72(1H*),/)')
       WRITE(11,*)'New Coordinates after TS Optimization (Angs)'
       WRITE(11,*)
       DO I=1,NAT
        WRITE(11,'(I5,3F15.4)')                                         &
             I,Cxyz(1,I)*BOHR,Cxyz(2,I)*BOHR,Cxyz(3,I)*BOHR
       ENDDO
       CALL NUCDIST(NV,NAT,Cxyz)
       WRITE(11,'(/,72(1H*),/)')
       IF(IHSSCAL==1)THEN
        WRTHSSNOTE = .TRUE.
        WRTHSSCGO = .TRUE.
        CALL HESSCAL(NINTEG,IDONTW,NAT,ZAN,Cxyz,IAN,IMIN,IMAX,ZMASS,    &
                     KSTART,KATOM,KTYPE,KLOC,INTYP,KNG,KMIN,KMAX,ISH,   &
                     ITYP,C1,C2,EX1,CS,CP,CD,CF,CG,CH,CI,DIPS,GRADS,    &
                     IRUNTYP,IPROJECT,ISIGMA,SIZE_ENV,ENV,ATM,NBAS,BAS, &
                     IGTYP,KTYPEaux,IZCORE)
       END IF
      ELSE
       WRITE(6,'(/1X,A,I4,2X,A,1PE10.3,2X,A,I3)')                       &
            'TSOPT stopped after',TSMAXIT,'max|Grad| =',GEOMAXG,        &
            'Nneg =',NNEG
      END IF
      IF(.NOT.(GEOMAXG<TSGTOL .AND. NNEG==1))THEN
       WRITE(11,'(/,72(1H*),/)')
       WRITE(11,*)'New Coordinates after TS Optimization (Angs)'
       WRITE(11,*)
       DO I=1,NAT
        WRITE(11,'(I5,3F15.4)')                                         &
             I,Cxyz(1,I)*BOHR,Cxyz(2,I)*BOHR,Cxyz(3,I)*BOHR
       ENDDO
       CALL NUCDIST(NV,NAT,Cxyz)
       WRITE(11,'(/,72(1H*),/)')
      END IF
      IF(MOLDEN==1)WRITE(18,'(A9,/,A6)')'[GEOCONV]','energy'
   10 READ(LUGEOCONV,*,END=20) ENERTRY
      WRITE(18,'(F20.10)')ENERTRY
      GO TO 10
   20 CONTINUE
      CLOSE(LUGEOCONV)
      WRTHSSNOTE = .TRUE.
      DEALLOCATE(HESSIAN,VEC,DDM,EIG,WORK,STEPV,GPROJ)
      RETURN
      END

! TSPROJMODES
      SUBROUTINE TSPROJMODES(NAT,NC1,Cxyz,HESS,GRADS,ZAN,ZMASS,         &
                             IPROJECT,VEC,EIG,IMODE,NNEG,IZCORE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL EFIELDL,PRJGRD,PRJROT
      COMMON/INP_EFIELDL/EFX,EFY,EFZ,EFIELDL
      INTEGER,INTENT(IN) :: NAT,NC1,IPROJECT
      INTEGER,INTENT(OUT) :: IMODE,NNEG
      INTEGER,DIMENSION(NAT),INTENT(IN) :: IZCORE
      DOUBLE PRECISION,DIMENSION(3,NAT),INTENT(IN) :: Cxyz
      DOUBLE PRECISION,DIMENSION(NC1,NC1),INTENT(IN) :: HESS
      DOUBLE PRECISION,DIMENSION(NC1),INTENT(IN) :: GRADS
      DOUBLE PRECISION,DIMENSION(NAT),INTENT(IN) :: ZAN,ZMASS
      DOUBLE PRECISION,DIMENSION(NC1,NC1),INTENT(OUT) :: VEC
      DOUBLE PRECISION,DIMENSION(NC1),INTENT(OUT) :: EIG
      DOUBLE PRECISION,DIMENSION(NC1,NC1) :: HPROJ
      DOUBLE PRECISION,DIMENSION(NC1) :: RM,WORK
      DOUBLE PRECISION,DIMENSION(3,NAT) :: COM,GRAD3
      DOUBLE PRECISION,DIMENSION(3) :: CMASS
      DOUBLE PRECISION :: ZMASST
      INTEGER :: I,J,NNEGRAW,NSKIP,NLAST,NFIRST,MAYBE,NUCZ
!-----------------------------------------------------------------------
      HPROJ = HESS
      DO I=1,NAT
       DO J=1,3
        GRAD3(J,I) = GRADS(3*(I-1)+J)
       ENDDO
      ENDDO
      DO I=1,NAT
       IF(ZMASS(I)>0.0D0)THEN
        RM(3*(I-1)+1) = 1.0D0/DSQRT(ZMASS(I))
       ELSE
        RM(3*(I-1)+1) = 1.0D0
       END IF
       RM(3*(I-1)+2) = RM(3*(I-1)+1)
       RM(3*(I-1)+3) = RM(3*(I-1)+1)
      ENDDO
      CALL HESMASd(NC1,HPROJ,RM)
      CALL CENMASd(NAT,Cxyz,COM,ZMASST,CMASS,ZMASS)
      IF(IPROJECT==1)THEN
       PRJROT = .TRUE.
       IF(EFIELDL) PRJROT = .FALSE.
       PRJGRD = .FALSE.
       CALL PRJFCMd(PRJGRD,PRJROT,ZMASST,HPROJ,COM,GRAD3,RM,NAT,NC1)
      END IF
      CALL DIAG(NC1,HPROJ,VEC,EIG,WORK)
      CALL STFASEd(VEC,NC1,NC1,NC1)
      DO I = 1,NC1
       DO J = 1,NC1
        VEC(J,I) = VEC(J,I)*RM(J)
       END DO
      END DO
      NNEGRAW = 0
      DO I = 1,NC1
       IF(EIG(I)<0.0D0) NNEGRAW = NNEGRAW + 1
      ENDDO
      IF(IPROJECT==0)THEN
       NNEG = NNEGRAW
       IF(NNEG>=1)THEN
        IMODE = 1
       ELSE
        IMODE = 1
       END IF
       RETURN
      END IF
      NSKIP = 6
      IF(NAT==2) NSKIP = 5
      DO I = 1,NAT
       IF(INT(ZAN(I))+IZCORE(I)==0) NSKIP = NSKIP + 3
      ENDDO
      NLAST = NSKIP
      DO I = 1,NNEGRAW
       MAYBE = I + NSKIP
       IF(MAYBE>NC1) EXIT
       IF(DABS(EIG(I))>DABS(EIG(MAYBE))) NLAST = NLAST + 1
      ENDDO
      NFIRST = NLAST - NSKIP + 1
      NNEG = NFIRST - 1
      IF(NNEG>=1)THEN
       IMODE = 1
      ELSE
       IMODE = MIN(NC1,NSKIP+1)
      END IF
      RETURN
      END

! OPTSUMSL
      SUBROUTINE OPTSUMSL(NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NSHELL,    &
                          NPRIMI,ZAN,Cxyz,IAN,IMIN,IMAX,ZMASS,KSTART,   &
                          KATOM,KTYPE,KLOC,INTYP,KNG,KMIN,KMAX,ISH,     &
                          ITYP,C1,C2,EX,CS,CP,CD,CF,CG,CH,CI,IRUNTYP,   &
                          SIZE_ENV,ENV,ATM,NBAS,BAS,IGTYP,GRAD,DIPS,    &
                          IZCORE,KTYPEaux,NSHELLaux)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      PARAMETER (BOHR = 0.52917724924D+00) 
      INTEGER,DIMENSION(NAT)    :: IAN,IMIN,IMAX,IZCORE
      INTEGER,DIMENSION(NSHELL) :: KSTART,KATOM,KTYPE,KLOC
      INTEGER,DIMENSION(NSHELL) :: INTYP,KNG,KMIN,KMAX
      INTEGER,DIMENSION(NSHELLaux) :: KTYPEaux
      INTEGER,DIMENSION(NPRIMI) :: ISH,ITYP
      INTEGER,INTENT(IN)        :: IRUNTYP
      DOUBLE PRECISION,DIMENSION(NAT)   :: ZAN,ZMASS
      DOUBLE PRECISION,DIMENSION(NPRIMI):: C1,C2,EX,CS,CP,CD,CF,CG,CH,CI
      DOUBLE PRECISION,DIMENSION(3,NAT) :: Cxyz
      DOUBLE PRECISION,DIMENSION(3*NAT),INTENT(INOUT) :: GRAD
      DOUBLE PRECISION,DIMENSION(3) :: DIPS

      INTEGER::SIZE_ENV, IGTYP
      DOUBLE PRECISION,DIMENSION(SIZE_ENV) :: ENV
      INTEGER,DIMENSION(6*NAT) :: ATM
      INTEGER,DIMENSION(8*NBAS) :: BAS
!      
      COMMON/EHFEN/EHF,EN
      COMMON/ENERGIAS/EELEC_OLD,EELEC,DIF_EELEC,EELEC_MIN
      COMMON/INPNOF_THRESH/THRESHL,THRESHE,THRESHEC,THRESHEN
      COMMON/INPNOF_GENERALINF/ICOEF,ISOFTMAX,IORBOPT,MAXIT,MAXIT21
      COMMON/POINTERUSER/NIU1,NIU2,NIU3,NIU4,NIU5,NIU6,NIU7,NIU8,NIU9,  &
                         NIU10,NIU11,NIU12,NIU13,NIU14,NIU15,NIU16,     &
                         NIU17,NIU18,NIU19,NIU20,NIU21,NIU22,NIU23,     &
                         NIU24,NIU25,NIU26,NIU27,NIU28,NIU29,NIU30,     &
                         NIU31,NIU32,NIULAST,NU1,NU2,NU3,NU4,NU5,NU6,   &
                         NU7,NU8,NU9,NU10,NU11,NU12,NU13,NU14,NU15,    &
                         NU16,NULAST
      COMMON/INPNOF_MOLDEN/MOLDEN
      COMMON/GEOSTAT/GEOMAXG,GEORMSG,IOPTIRC
!
      INTEGER,ALLOCATABLE,DIMENSION(:) :: IUSER,IV
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: USER,D,V     
      EXTERNAL CALCOPTE,CALCOPTG
      DOUBLE PRECISION,PARAMETER :: GMAXWARN=1.0D-3
      INTEGER :: LUGEOCONV
!-----------------------------------------------------------------------
!     Define Pointers of the USER array
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL POINTERSOPT(NAT,NSHELL,NPRIMI,SIZE_ENV,NBAS,IGTYP,NSHELLaux)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Transfer working arrays to IUSER and USER
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(IUSER(NIULAST),USER(NULAST))
!
      IUSER(NIU1)          = NINTEG
      IUSER(NIU2)          = IDONTW  
      IUSER(NIU3)          = 0
      IUSER(NIU4)          = IEMOM  
      IUSER(NIU5)          = NAT
      IUSER(NIU6)          = NBF
      IUSER(NIU7)          = NSHELL
      IUSER(NIU8)          = NPRIMI
      IUSER(NIU9 :NIU10-1) = IAN
      IUSER(NIU10:NIU11-1) = IMIN
      IUSER(NIU11:NIU12-1) = IMAX
      IUSER(NIU12:NIU13-1) = KSTART 
      IUSER(NIU13:NIU14-1) = KATOM
      IUSER(NIU14:NIU15-1) = KTYPE
      IUSER(NIU15:NIU16-1) = KLOC
      IUSER(NIU16:NIU17-1) = INTYP
      IUSER(NIU17:NIU18-1) = KNG
      IUSER(NIU18:NIU19-1) = KMIN
      IUSER(NIU19:NIU20-1) = KMAX
      IUSER(NIU20:NIU21-1) = ISH
      IUSER(NIU21:NIU22-1) = ITYP
      IUSER(NIU22)         = IRUNTYP
      IUSER(NIU23)         = 0             ! ITCG
      IUSER(NIU24)         = NBFaux
      IUSER(NIU25)         = IGTYP
      IUSER(NIU26)         = NBAS
      IUSER(NIU27)         = SIZE_ENV
      IUSER(NIU28:NIU29-1) = ATM
      IUSER(NIU29:NIU30-1) = BAS
      IUSER(NIU30:NIU31-1) = KTYPEaux
      IUSER(NIU31)         = NSHELLaux
      IUSER(NIU32:NIULAST-1) = IZCORE
!
      USER(NU1 :NU2-1)     = ZAN
      USER(NU2 :NU3-1)     = ZMASS
      USER(NU3 :NU4-1)     = C1
      USER(NU4 :NU5-1)     = C2
      USER(NU5 :NU6-1)     = EX
      USER(NU6 :NU7-1)     = CS
      USER(NU7 :NU8-1)     = CP
      USER(NU8 :NU9-1)     = CD
      USER(NU9 :NU10-1)    = CF
      USER(NU10:NU11-1)    = CG
      USER(NU11:NU12-1)    = CH
      USER(NU12:NU13-1)    = CI
      USER(NU13:NU14-1)    = DIPS
      USER(NU14:NU15-1)    = ENV
      USER(NU15)           = 1.0D20
      USER(NU16:NULAST-1)  = 0.0D0
      OPEN(NEWUNIT=LUGEOCONV,STATUS='SCRATCH',FORM='FORMATTED')
      IUSER(NIU3)          = LUGEOCONV
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Write Initial Coordinates on File CGGRAD (Unit=11)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      WRITE(11,'(/,72(1H*),/)')
      WRITE(11,*)'Initial Coordinates (Angs) for Geometry Optimization'
      WRITE(11,*)
      DO I=1,NAT
       WRITE(11,'(I5,3F15.4)')  &
             I,Cxyz(1,I)*BOHR,Cxyz(2,I)*BOHR,Cxyz(3,I)*BOHR
      ENDDO
      WRITE(11,'(/,72(1H*),/)')      
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Minimization of the total energy with respect to Cxyz
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      WRITE(6,'(/,2X,A21,/1X,23(1H=))')'GEOMETRY OPTIMIZATION'
      WRITE(6,'(/,1X,A23,7X,A19,/)')'Call in CG Optimization',          &
                                    'Total Energy (a.u.)'
      IF(MOLDEN==1)WRITE(18,'(A18)')'[GEOMETRIES] (XYZ)'
!
      NV = 3*NAT
      LIV = 60
      LV = 71+NV*(NV+15)/2
!      
      ALLOCATE( D(NV),IV(LIV),V(LV) ) 
!      
      D(1:NV) = 1.0d0
      CALL DEFLT(2,IV,LIV,LV,V)
      IV(1) = 12
      ETOL = THRESHE
      ETOL = DMAX1(ETOL,THRESHEN)
      IF(IORBOPT==1) ETOL = DMAX1(ETOL,THRESHEC)
      V(31) = ETOL
      V(32) = ETOL
      EMINIMA = 1.0d6
!
      CALL SUMSL(NV,D,Cxyz,CALCOPTE,CALCOPTG,IV,LIV,LV,V,IUSER,USER)

!     iv(irc) is a return code having one of the following values
!
!        1 = switch models or try smaller step.
!        2 = switch models or accept step.
!        3 = accept step and determine v(radfac) by gradient
!             tests.
!        4 = accept step, v(radfac) has been determined.
!        5 = recompute step (using the same model).
!        6 = recompute step with radius = v(lmaxs) but do not
!             evaulate the objective function.
!        7 = x-convergence (see v(xctol)).
!        8 = relative function convergence (see v(rfctol)).
!        9 = both x- and relative function convergence.
!       10 = absolute function convergence (see v(afctol)).
!       11 = singular convergence (see v(lmaxs)).
!       12 = false convergence (see v(xftol)).
!       13 = iv(irc) was out of range on input.

      write(6,'(/,A10,I3)')' IV(irc) =',IV(29)
!      
      ENERGIA = EELEC + EN
      IF(EMINIMA>ENERGIA)EMINIMA = ENERGIA
      IF(USER(NU15)<ENERGIA-1.0D-12)THEN
       WRITE(6,'(/1X,A/)')                                              &
       'Promoting lowest-energy geometry reached during OPTGEO'
       CALL COPYGCFE2GCFl
       DO I=1,NAT
        Cxyz(1,I) = USER(NU16+3*(I-1))
        Cxyz(2,I) = USER(NU16+3*(I-1)+1)
        Cxyz(3,I) = USER(NU16+3*(I-1)+2)
        USER(NU14-1+20+3*(I-1)+1) = Cxyz(1,I)
        USER(NU14-1+20+3*(I-1)+2) = Cxyz(2,I)
        USER(NU14-1+20+3*(I-1)+3) = Cxyz(3,I)
       END DO
       CALL ENERGRAD(IUSER(NIU1),IUSER(NIU2),IUSER(NIU4),IUSER(NIU5),   &
                     IUSER(NIU6),IUSER(NIU24),IUSER(NIU7),IUSER(NIU8),  &
                     USER(NU1),Cxyz,IUSER(NIU9),IUSER(NIU10),           &
                     IUSER(NIU11),IUSER(NIU12),IUSER(NIU13),            &
                     IUSER(NIU14),IUSER(NIU15),IUSER(NIU16),            &
                     IUSER(NIU17),IUSER(NIU18),IUSER(NIU19),            &
                     IUSER(NIU20),IUSER(NIU21),USER(NU3),USER(NU4),     &
                     USER(NU5),USER(NU6),USER(NU7),USER(NU8),USER(NU9), &
                     USER(NU10),USER(NU11),USER(NU12),GRAD,             &
                     IUSER(NIU22),USER(NU13),IUSER(NIU32),IUSER(NIU27), &
                     USER(NU14),IUSER(NIU28),IUSER(NIU26),IUSER(NIU29), &
                     IUSER(NIU25),IUSER(NIU30),IUSER(NIU31),0,0)
       ENERGIA = EELEC + EN
       EMINIMA = ENERGIA
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Compute energy at solution to print data on the output file
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      WRITE(6,'(/1X,A24,2X,F20.10,/)')'Final objective value = ',EMINIMA      
      WRITE(6,*)'New Coordinates after Conjugate Gradient Opt (Bohr)'
      WRITE(6,'(1X,51(1H-),/)')
      DO I=1,NAT
       WRITE(6,'(I5,3F15.4)')I,Cxyz(1,I),Cxyz(2,I),Cxyz(3,I)
      ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Update coordinates of shells if use libcint library for ERIs
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO I=1,NAT
       USER(NU14-1+20+3*(I-1)+1) = Cxyz(1,I)
       USER(NU14-1+20+3*(I-1)+2) = Cxyz(2,I)
       USER(NU14-1+20+3*(I-1)+3) = Cxyz(3,I)
      END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL ENERGRAD(IUSER(NIU1),IUSER(NIU2),IUSER(NIU4),IUSER(NIU5),    &
                    IUSER(NIU6),IUSER(NIU24),IUSER(NIU7),IUSER(NIU8),   &
                    USER(NU1),Cxyz,IUSER(NIU9),IUSER(NIU10),            &
                    IUSER(NIU11),IUSER(NIU12),IUSER(NIU13),IUSER(NIU14),&
                    IUSER(NIU15),IUSER(NIU16),IUSER(NIU17),IUSER(NIU18),&
                    IUSER(NIU19),IUSER(NIU20),IUSER(NIU21),USER(NU3),   &
                    USER(NU4),USER(NU5),USER(NU6),USER(NU7),USER(NU8),  &
                    USER(NU9),USER(NU10),USER(NU11),USER(NU12),GRAD,    &
                    IUSER(NIU22),USER(NU13),IUSER(NIU32),IUSER(NIU27),  &
                    USER(NU14),IUSER(NIU28),IUSER(NIU26),IUSER(NIU29),  &
                    IUSER(NIU25),IUSER(NIU30),IUSER(NIU31),0,0)
      GEOMAXG = 0.0D0
      GEORMSG = 0.0D0
      DO I=1,NV
       GEOMAXG = DMAX1(GEOMAXG,DABS(GRAD(I)))
       GEORMSG = GEORMSG + GRAD(I)*GRAD(I)
      END DO
      GEORMSG = DSQRT(GEORMSG/DBLE(NV))
      IOPTIRC = IV(29)
      WRITE(6,'(/1X,A,1PE12.4,2X,A,1PE12.4)')                           &
      'Final max|Grad| = ',GEOMAXG,'RMS|Grad| = ',GEORMSG
      IF(IOPTIRC==12 .AND. GEOMAXG>GMAXWARN)THEN
       WRITE(6,'( 1X,A)')                                               &
       'Warning: SUMSL returned IV(irc)=12 with a non-small gradient.'
       WRITE(6,'( 1X,A)')                                               &
       'A restart from the final GCF is recommended before HESSCAL.'
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Write Final Coordinates on File CGGRAD (Unit=100)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      WRITE(11,'(/,72(1H*),/)')
      WRITE(11,*)'New Coordinates after Conjugate Gradient Opt (Angs)'
      WRITE(11,*)
      DO I=1,NAT
       WRITE(11,'(I5,3F15.4)')                                          &
            I,Cxyz(1,I)*BOHR,Cxyz(2,I)*BOHR,Cxyz(3,I)*BOHR
      ENDDO
      CALL NUCDIST(NV,NAT,Cxyz)
      WRITE(11,'(/,72(1H*),/)')  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Write geometry energies in [GEOCONV] section of the molden file
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -      
      IF(MOLDEN==1)WRITE(18,'(A9,/,A6)')'[GEOCONV]','energy'
      REWIND(IUSER(NIU3))
   10 READ(IUSER(NIU3),*,END=20) ENERGIA
      WRITE(18,'(F20.10)')ENERGIA
      GO TO 10
   20 CONTINUE
!-----------------------------------------------------------------------
      CLOSE(IUSER(NIU3))
      DEALLOCATE(IUSER,USER,D,IV,V)
      RETURN
      END

! POINTERSOPT
      SUBROUTINE POINTERSOPT(NAT,NSHELL,NPRIMI,SIZE_ENV,NBAS,IGTYP,     &
                             NSHELLaux)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      COMMON/POINTERUSER/NIU1,NIU2,NIU3,NIU4,NIU5,NIU6,NIU7,NIU8,NIU9,  &
                         NIU10,NIU11,NIU12,NIU13,NIU14,NIU15,NIU16,     &
                         NIU17,NIU18,NIU19,NIU20,NIU21,NIU22,NIU23,     &
                         NIU24,NIU25,NIU26,NIU27,NIU28,NIU29,NIU30,     &
                         NIU31,NIU32,NIULAST,NU1,NU2,NU3,NU4,NU5,NU6,   &
                         NU7,NU8,NU9,NU10,NU11,NU12,NU13,NU14,NU15,    &
                         NU16,NULAST
      INTEGER :: SIZE_ENV,NBAS
!-----------------------------------------------------------------------
!     Define Pointers of the USER array
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      NIU1  = 1                         ! NINTEG
      NIU2  = NIU1  + 1                 ! IDONTW 
      NIU3  = NIU2  + 1                 ! Scratch unit
      NIU4  = NIU3  + 1                 ! IEMOM  
      NIU5  = NIU4  + 1                 ! NAT
      NIU6  = NIU5  + 1                 ! NBF
      NIU7  = NIU6  + 1                 ! NSHELL
      NIU8  = NIU7  + 1                 ! NPRIMI
      NIU9  = NIU8  + 1                 ! IAN
      NIU10 = NIU9  + NAT               ! IMIN
      NIU11 = NIU10 + NAT               ! IMAX
      NIU12 = NIU11 + NAT               ! KSTART 
      NIU13 = NIU12 + NSHELL            ! KATOM
      NIU14 = NIU13 + NSHELL            ! KTYPE
      NIU15 = NIU14 + NSHELL            ! KLOC
      NIU16 = NIU15 + NSHELL            ! INTYP
      NIU17 = NIU16 + NSHELL            ! KNG
      NIU18 = NIU17 + NSHELL            ! KMIN
      NIU19 = NIU18 + NSHELL            ! KMAX 
      NIU20 = NIU19 + NSHELL            ! ISH
      NIU21 = NIU20 + NPRIMI            ! ITYP
      NIU22 = NIU21 + NPRIMI            ! IRUNTYP
      NIU23 = NIU22 + 1                 ! ITCG
      NIU24 = NIU23 + 1                 ! NBFaux
      NIU25 = NIU24 + 1                 ! IGTYP
      NIU26 = NIU25 + 1                 ! NBAS
      NIU27 = NIU26 + 1                 ! SIZE_ENV
      NIU28 = NIU27 + 1                 ! ATM
      NIU29 = NIU28 + 6*NAT             ! BAS
      NIU30 = NIU29 + 8*NBAS            ! KTYPEaux
      NIU31 = NIU30 + NSHELLaux         ! NSHELLaux
      NIU32 = NIU31 + 1                 ! IZCORE
      NIULAST = NIU32 + NAT
!
      NU1  = 1                          ! ZAN
      NU2  = NU1  + NAT                 ! ZMASS
      NU3  = NU2  + NAT                 ! C1
      NU4  = NU3  + NPRIMI              ! C2
      NU5  = NU4  + NPRIMI              ! EX   
      NU6  = NU5  + NPRIMI              ! CS
      NU7  = NU6  + NPRIMI              ! CP
      NU8  = NU7  + NPRIMI              ! CD
      NU9  = NU8  + NPRIMI              ! CF
      NU10 = NU9  + NPRIMI              ! CG
      NU11 = NU10 + NPRIMI              ! CH
      NU12 = NU11 + NPRIMI              ! CI
      NU13 = NU12 + NPRIMI              ! DIPS
      NU14 = NU13 + 3                   ! ENV
      NU15 = NU14 + SIZE_ENV            ! GEOBESTE
      NU16 = NU15 + 1                   ! GEOBESTXYZ
      NULAST = NU16 + 3*NAT
!-----------------------------------------------------------------------
      RETURN
      END      

! CALCOPTE      
      SUBROUTINE CALCOPTE(NV,Cxyz,NF,ENERGIA,IUSER,USER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      COMMON/EHFEN/EHF,EN
      COMMON/ENERGIAS/EELEC,EELEC_OLD,DIF_EELEC,EELEC_MIN
      COMMON/INPNOF_MOLDEN/MOLDEN
      COMMON/INPNOF_MOLDENGEO/MOLDENGEO
      COMMON/POINTERUSER/NIU1,NIU2,NIU3,NIU4,NIU5,NIU6,NIU7,NIU8,NIU9,  &
                         NIU10,NIU11,NIU12,NIU13,NIU14,NIU15,NIU16,     &
                         NIU17,NIU18,NIU19,NIU20,NIU21,NIU22,NIU23,     &
                         NIU24,NIU25,NIU26,NIU27,NIU28,NIU29,NIU30,     &
                         NIU31,NIU32,NIULAST,NU1,NU2,NU3,NU4,NU5,NU6,   &
                         NU7,NU8,NU9,NU10,NU11,NU12,NU13,NU14,NU15,    &
                         NU16,NULAST
!
      INTEGER,DIMENSION(NIULAST)         :: IUSER 
      DOUBLE PRECISION,DIMENSION(NV)     :: Cxyz
      DOUBLE PRECISION,DIMENSION(NULAST) :: USER
      PARAMETER (BOHR = 0.52917724924D+00) 
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
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: GRAD
!-----------------------------------------------------------------------
      IUSER(NIU23) = IUSER(NIU23) + 1                  ! ITCG = ITCG + 1
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -            
!     Update coordinates of shells if use libcint library for ERIs
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO I=1,NV
       USER(NU14+20+I-1) = Cxyz(I)
      END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(GRAD(NV))
      CALL ENERGRAD(IUSER(NIU1),IUSER(NIU2),IUSER(NIU4),IUSER(NIU5),    &
                    IUSER(NIU6),IUSER(NIU24),IUSER(NIU7),IUSER(NIU8),   &
                    USER(NU1),Cxyz,IUSER(NIU9),IUSER(NIU10),            &
                    IUSER(NIU11),IUSER(NIU12),IUSER(NIU13),IUSER(NIU14),&
                    IUSER(NIU15),IUSER(NIU16),IUSER(NIU17),IUSER(NIU18),&
                    IUSER(NIU19),IUSER(NIU20),IUSER(NIU21),USER(NU3),   &
                    USER(NU4),USER(NU5),USER(NU6),USER(NU7),USER(NU8),  &
                    USER(NU9),USER(NU10),USER(NU11),USER(NU12),GRAD,    &
                    IUSER(NIU22),USER(NU13),IUSER(NIU32),IUSER(NIU27),  &
                    USER(NU14),IUSER(NIU28),IUSER(NIU26),IUSER(NIU29),  &
                    IUSER(NIU25),IUSER(NIU30),IUSER(NIU31),0,0)
!
      ENERGIA = EELEC + EN     
      WRITE(6,'(8X,I3,16X,F20.10)')NF,ENERGIA
      WRITE(IUSER(NIU3),'(F20.10)')ENERGIA
      IF(ENERGIA<USER(NU15))THEN
       USER(NU15) = ENERGIA
       USER(NU16:NU16+NV-1) = Cxyz(1:NV)
      END IF
      DEALLOCATE(GRAD)
!
      if(MOLDEN==1)then 
       nat = IUSER(NIU5)
       write(18,'(I6)')nat
       write(18,'(I6,F20.10)')NF,ENERGIA
       do jat=1,nat
        j = (jat-1)*3
        IZNUC = INT(USER(NU1+jat-1))+IUSER(NIU32+jat-1)
        write(18,'(1X,A4,3F15.4)')ATMLAB(IZNUC),                        &
        Cxyz(1+j)*BOHR,Cxyz(2+j)*BOHR,Cxyz(3+j)*BOHR
       enddo
      endif
      IF(MOLDEN==1 .AND. MOLDENGEO==1) CALL WRITE_OPT_MLD_SNAPSHOT(NF)
!-----------------------------------------------------------------------
      RETURN
      END

! COPYGCFE2GCFl
      SUBROUTINE COPYGCFE2GCFl
      IMPLICIT NONE
      CHARACTER(LEN=256) :: LINE
!-----------------------------------------------------------------------
      REWIND(8)
      REWIND(3)
   10 CONTINUE
      READ(8,'(A)',END=20) LINE
      WRITE(3,'(A)') LINE
      GO TO 10
   20 CONTINUE
      ENDFILE(3)
      REWIND(3)
      REWIND(8)
!-----------------------------------------------------------------------
      RETURN
      END
      
! CALCOPTG      
      SUBROUTINE CALCOPTG(NV,Cxyz,NF,GRAD,IUSER,USER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      INTEGER :: NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/POINTERUSER/NIU1,NIU2,NIU3,NIU4,NIU5,NIU6,NIU7,NIU8,NIU9,  &
                         NIU10,NIU11,NIU12,NIU13,NIU14,NIU15,NIU16,     &
                         NIU17,NIU18,NIU19,NIU20,NIU21,NIU22,NIU23,     &
                         NIU24,NIU25,NIU26,NIU27,NIU28,NIU29,NIU30,     &
                         NIU31,NIU32,NIULAST,NU1,NU2,NU3,NU4,NU5,NU6,   &
                         NU7,NU8,NU9,NU10,NU11,NU12,NU13,NU14,NU15,    &
                         NU16,NULAST
      INTEGER,DIMENSION(NIULAST)         :: IUSER 
      DOUBLE PRECISION,DIMENSION(NV)     :: Cxyz
      DOUBLE PRECISION,DIMENSION(NV)     :: GRAD      
      DOUBLE PRECISION,DIMENSION(NULAST) :: USER
!-----------------------------------------------------------------------
!     Avoiding warnings
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      NF = NF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -      
!     Update coordinates of shells if use libcint library for ERIs
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO I=1,NV
       USER(NU14+20+I-1) = Cxyz(I)
      END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL ENERGRAD(IUSER(NIU1),IUSER(NIU2),IUSER(NIU4),IUSER(NIU5),    &
                    IUSER(NIU6),IUSER(NIU24),IUSER(NIU7),IUSER(NIU8),   &
                    USER(NU1),Cxyz,IUSER(NIU9),IUSER(NIU10),            &
                    IUSER(NIU11),IUSER(NIU12),IUSER(NIU13),IUSER(NIU14),&
                    IUSER(NIU15),IUSER(NIU16),IUSER(NIU17),IUSER(NIU18),&
                    IUSER(NIU19),IUSER(NIU20),IUSER(NIU21),USER(NU3),   &
                    USER(NU4),USER(NU5),USER(NU6),USER(NU7),USER(NU8),  &
                    USER(NU9),USER(NU10),USER(NU11),USER(NU12),GRAD,    &
                    IUSER(NIU22),USER(NU13),IUSER(NIU32),IUSER(NIU27),  &
                    USER(NU14),IUSER(NIU28),IUSER(NIU26),IUSER(NIU29),  &
                    IUSER(NIU25),IUSER(NIU30),IUSER(NIU31),0,0)
      ENERGIA = EELEC + EN
      IF(ENERGIA<USER(NU15))THEN
       USER(NU15) = ENERGIA
       USER(NU16:NU16+NV-1) = Cxyz(1:NV)
      END IF
!-----------------------------------------------------------------------
      RETURN
      END
      
! OPTCGNAG
      SUBROUTINE OPTCGNAG(NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NSHELL,    &
                          NPRIMI,ZAN,Cxyz,IAN,IMIN,IMAX,ZMASS,KSTART,   &
                          KATOM,KTYPE,KLOC,INTYP,KNG,KMIN,KMAX,ISH,     &
                          ITYP,C1,C2,EX,CS,CP,CD,CF,CG,CH,CI,IRUNTYP,   &
                          GRADIENT,DIPS,IZCORE,KTYPEaux,NSHELLaux)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      PARAMETER (BOHR = 0.52917724924D+00) 
      DOUBLE PRECISION,DIMENSION(NAT) :: ZAN,ZMASS
      INTEGER,DIMENSION(NAT):: IAN,IMIN,IMAX,IZCORE
      INTEGER,DIMENSION(NSHELL) :: KSTART,KATOM,KTYPE,KLOC
      INTEGER,DIMENSION(NSHELL) :: INTYP,KNG,KMIN,KMAX
      INTEGER,DIMENSION(NSHELLaux) :: KTYPEaux
      INTEGER,DIMENSION(NPRIMI) :: ISH,ITYP
      DOUBLE PRECISION,DIMENSION(NPRIMI):: C1,C2,EX,CS,CP,CD,CF,CG,CH,CI
      DOUBLE PRECISION,DIMENSION(3,NAT):: Cxyz
      DOUBLE PRECISION,DIMENSION(3*NAT),INTENT(OUT):: GRADIENT
      INTEGER,INTENT(IN):: IRUNTYP
      DOUBLE PRECISION,DIMENSION(3):: DIPS
      COMMON/POINTERUSER/NIU1,NIU2,NIU3,NIU4,NIU5,NIU6,NIU7,NIU8,NIU9,  &
                         NIU10,NIU11,NIU12,NIU13,NIU14,NIU15,NIU16,     &
                         NIU17,NIU18,NIU19,NIU20,NIU21,NIU22,NIU23,     &
                         NIU24,NIU25,NIU26,NIU27,NIU28,NIU29,NIU30,     &
                         NIU31,NIU32,NIULAST,NU1,NU2,NU3,NU4,NU5,NU6,   &
                         NU7,NU8,NU9,NU10,NU11,NU12,NU13,NU14,NU15,    &
                         NU16,NULAST
      COMMON/INPNOF_MOLDEN/MOLDEN
!                         
      INTEGER,ALLOCATABLE,DIMENSION(:) :: IUSER,IWORK
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: USER,WORK,GRADS
      EXTERNAL ENERGYFUN
      INTEGER :: LUGEOCONV
!-----------------------------------------------------------------------
!     Define Pointers of the user arrays
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     CALL POINTERSOPT(NAT,NSHELL,NPRIMI,SIZE_ENV,NBAS,IGTYP,NSHELLaux)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Transfer working arrays to IUSER and USER
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(IUSER(NIULAST),USER(NULAST))
!
      IUSER(NIU1)          = NINTEG
      IUSER(NIU2)          = IDONTW  
      IUSER(NIU3)          = 0
      IUSER(NIU4)          = IEMOM  
      IUSER(NIU5)          = NAT
      IUSER(NIU6)          = NBF
      IUSER(NIU7)          = NSHELL
      IUSER(NIU8)          = NPRIMI
      IUSER(NIU9 :NIU10-1) = IAN
      IUSER(NIU10:NIU11-1) = IMIN
      IUSER(NIU11:NIU12-1) = IMAX
      IUSER(NIU12:NIU13-1) = KSTART 
      IUSER(NIU13:NIU14-1) = KATOM
      IUSER(NIU14:NIU15-1) = KTYPE
      IUSER(NIU15:NIU16-1) = KLOC
      IUSER(NIU16:NIU17-1) = INTYP
      IUSER(NIU17:NIU18-1) = KNG
      IUSER(NIU18:NIU19-1) = KMIN
      IUSER(NIU19:NIU20-1) = KMAX
      IUSER(NIU20:NIU21-1) = ISH
      IUSER(NIU21:NIU22-1) = ITYP
      IUSER(NIU22)         = IRUNTYP
      IUSER(NIU23)         = 0             ! ITCG
      IUSER(NIU24)         = NBFaux
      IUSER(NIU25)         = IGTYP
      IUSER(NIU26)         = NBAS
      IUSER(NIU27)         = SIZE_ENV
      IUSER(NIU28:NIU29-1) = ATM
      IUSER(NIU29:NIU30-1) = BAS
      IUSER(NIU30:NIU31-1) = KTYPEaux
      IUSER(NIU31)         = NSHELLaux
      IUSER(NIU32:NIULAST-1) = IZCORE
!
      USER(NU1 :NU2-1)     = ZAN
      USER(NU2 :NU3-1)     = ZMASS 
      USER(NU3 :NU4-1)     = C1
      USER(NU4 :NU5-1)     = C2
      USER(NU5 :NU6-1)     = EX
      USER(NU6 :NU7-1)     = CS
      USER(NU7 :NU8-1)     = CP
      USER(NU8 :NU9-1)     = CD
      USER(NU9 :NU10-1)    = CF
      USER(NU10:NU11-1)    = CG
      USER(NU11:NU12-1)    = CH
      USER(NU12:NU13-1)    = CI
      USER(NU13:NU14-1)    = DIPS
      USER(NU14:NU15-1)    = ENV
      USER(NU15)           = 1.0D20
      USER(NU16:NULAST-1)  = 0.0D0
      OPEN(NEWUNIT=LUGEOCONV,STATUS='SCRATCH',FORM='FORMATTED')
      IUSER(NIU3)          = LUGEOCONV
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Write Initial Coordinates on File CGGRAD (Unit=11)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      WRITE(11,'(/,72(1H*),/)')
      WRITE(11,*)'Initial Coordinates (Angs) for Geometry Optimization'
      WRITE(11,*)
      DO I=1,NAT
       WRITE(11,'(I5,3F15.4)')  &
             I,Cxyz(1,I)*BOHR,Cxyz(2,I)*BOHR,Cxyz(3,I)*BOHR
      ENDDO
      WRITE(11,'(/,72(1H*),/)')      
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Minimization of the total energy with respect to Cxyz
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      WRITE(6,'(/,2X,A21,/1X,23(1H=))')'GEOMETRY OPTIMIZATION'
      WRITE(6,'(/,1X,A23,7X,A19,/)')'Call in CG Optimization',          &
                                    'Total Energy (a.u.)'
      IF(MOLDEN==1)WRITE(18,'(A18)')'[GEOMETRIES] (XYZ)'
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Send output of E04DGF to CGM file
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      CALL X04ABF(1,2)                                           !nag
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -      
!     Est. opt. function val.
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     CALL E04DKF ('Es = ....')                                   !nag
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Function Precision (machine precision**0.9)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      CALL E04DKF ('Function Precision = 1.0D-6')                !nag
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Linesearch Tolerance (0<r<1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      CALL E04DKF ('Linesearch Tolerance = 0.1')                 !nag
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Maximum Step Length (Default = 10**20))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      CALL E04DKF ('Maximum Step Length  = 0.01')                !nag
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Optimality Tolerance (deafult = relative precision**0.8)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      CALL E04DKF ('Optimality Tolerance = 1.0D-5')              !nag
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Verify Level (-1 = No checks, 0 = cheap test, 1 = 0 + gradients)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      CALL E04DKF ('Verify Level = -1')                          !nag
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Print Level (0 = No output, 1 = The final solution only)
!                 (5 = One line of summary output for each iteration)
!                 (10 = The final solution and one line for each iter.)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      CALL E04DKF ('Print Level = 0')                            !nag
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Calling to NAG Library for using the CG method
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      NV = 3*NAT
      IFAIL = 1
      EMINIMA = 1.0d6
      ALLOCATE(IWORK(NV+1),WORK(13*NV),GRADS(NV))
!      CALL E04DGF(NV,ENERGYFUN,ITER_E04DGF,ENERGIA,GRADS,Cxyz,          &  !nag
!                  IWORK,WORK,IUSER,USER,IFAIL)                             !nag
      IF(EMINIMA>ENERGIA)EMINIMA=ENERGIA
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Compute energy at solution to print data on the output file
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(IFAIL==0)then
       WRITE(6,'(/1X,A24,2X,F20.10,/)')'Optimal solution found !',EMINIMA
      else if(IFAIL==6)then       
       WRITE(6,'(/1X,A24,2X,F20.10,/)')'Final objective value = ',EMINIMA      
      end if
      WRITE(6,*)'New Coordinates after Conjugate Gradient Opt (Bohr)'
      WRITE(6,'(1X,51(1H-),/)')
      DO I=1,NAT
       WRITE(6,'(I5,3F15.4)')I,Cxyz(1,I),Cxyz(2,I),Cxyz(3,I)
      ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -      
!     Update coordinates of shells if use libcint library for ERIs
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -      
      CALL ENERGRAD(IUSER(NIU1),IUSER(NIU2),IUSER(NIU4),IUSER(NIU5),    &
                    IUSER(NIU6),IUSER(NIU24),IUSER(NIU7),IUSER(NIU8),   &
                    USER(NU1),Cxyz,IUSER(NIU9),IUSER(NIU10),            &
                    IUSER(NIU11),IUSER(NIU12),IUSER(NIU13),IUSER(NIU14),&
                    IUSER(NIU15),IUSER(NIU16),IUSER(NIU17),IUSER(NIU18),&
                    IUSER(NIU19),IUSER(NIU20),IUSER(NIU21),USER(NU3),   &
                    USER(NU4),USER(NU5),USER(NU6),USER(NU7),USER(NU8),  &
                    USER(NU9),USER(NU10),USER(NU11),USER(NU12),GRADS,   &
                    IUSER(NIU22),USER(NU13),IUSER(NIU32),IUSER(NIU27),  &
                    USER(NU14),IUSER(NIU28),IUSER(NIU26),IUSER(NIU29),  &
                    IUSER(NIU25),IUSER(NIU30),IUSER(NIU31),0,0)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Write Final Coordinates on File CGGRAD (Unit=100)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      WRITE(11,'(/,72(1H*),/)')
      WRITE(11,*)'New Coordinates after Conjugate Gradient Opt (Angs)'
      WRITE(11,*)
      DO I=1,NAT
       WRITE(11,'(I5,3F15.4)')                                          &
            I,Cxyz(1,I)*BOHR,Cxyz(2,I)*BOHR,Cxyz(3,I)*BOHR
      ENDDO
      CALL NUCDIST(NV,NAT,Cxyz)
      WRITE(11,'(/,72(1H*),/)')      
      GRADIENT = GRADS
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Write geometry energies in [GEOCONV] section of the molden file
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -      
      IF(MOLDEN==1)WRITE(18,'(A9,/,A6)')'[GEOCONV]','energy'
      REWIND(IUSER(NIU3))
   10 READ(IUSER(NIU3),*,END=20) ENERGIA
      WRITE(18,'(F20.10)')ENERGIA
      GO TO 10
   20 CONTINUE
!-----------------------------------------------------------------------
      CLOSE(IUSER(NIU3))
      DEALLOCATE(IUSER,USER,GRADS,IWORK,WORK)
      RETURN
      END
      
! ENERGYFUN
      SUBROUTINE ENERGYFUN(MODE,NV,Cxyz,ENERGIA,GRADS,NSTATE,IUSER,USER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/POINTERUSER/NIU1,NIU2,NIU3,NIU4,NIU5,NIU6,NIU7,NIU8,NIU9,  &
                         NIU10,NIU11,NIU12,NIU13,NIU14,NIU15,NIU16,     &
                         NIU17,NIU18,NIU19,NIU20,NIU21,NIU22,NIU23,     &
                         NIU24,NIU25,NIU26,NIU27,NIU28,NIU29,NIU30,     &
                         NIU31,NIU32,NIULAST,NU1,NU2,NU3,NU4,NU5,NU6,   &
                         NU7,NU8,NU9,NU10,NU11,NU12,NU13,NU14,NU15,    &
                         NU16,NULAST
      COMMON/EHFEN/EHF,EN
      COMMON/ENERGIAS/EELEC_OLD,EELEC,DIF_EELEC,EELEC_MIN
      COMMON/INPNOF_MOLDEN/MOLDEN
      COMMON/INPNOF_MOLDENGEO/MOLDENGEO
!
      INTEGER,DIMENSION(*) :: IUSER
      DOUBLE PRECISION,DIMENSION(*) :: USER
      DOUBLE PRECISION,DIMENSION(NV) :: Cxyz,GRADS
      PARAMETER (BOHR = 0.52917724924D+00)
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
!-----------------------------------------------------------------------
      NSTATE = 1
      MODE = 1
      IUSER(NIU23) = IUSER(NIU23) + 1                  ! ITCG = ITCG + 1
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -      
!     Update coordinates of shells if use libcint library for ERIs
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL ENERGRAD(IUSER(NIU1),IUSER(NIU2),IUSER(NIU4),IUSER(NIU5),    &
                    IUSER(NIU6),IUSER(NIU24),IUSER(NIU7),IUSER(NIU8),   &
                    USER(NU1),Cxyz,IUSER(NIU9),IUSER(NIU10),            &
                    IUSER(NIU11),IUSER(NIU12),IUSER(NIU13),IUSER(NIU14),&
                    IUSER(NIU15),IUSER(NIU16),IUSER(NIU17),IUSER(NIU18),&
                    IUSER(NIU19),IUSER(NIU20),IUSER(NIU21),USER(NU3),   &
                    USER(NU4),USER(NU5),USER(NU6),USER(NU7),USER(NU8),  &
                    USER(NU9),USER(NU10),USER(NU11),USER(NU12),GRADS,   &
                    IUSER(NIU22),USER(NU13),IUSER(NIU32),IUSER(NIU27),  &
                    USER(NU14),IUSER(NIU28),IUSER(NIU26),IUSER(NIU29),  &
                    IUSER(NIU25),IUSER(NIU30),IUSER(NIU31),0,0)
      ENERGIA = EELEC + EN
      WRITE(6,'(8X,I3,16X,F20.10)')IUSER(NIU23),ENERGIA
      WRITE(IUSER(NIU3),'(F20.10)')ENERGIA
!
      if(MOLDEN==1)then 
       nat = IUSER(NIU5)
       write(18,'(I6)')nat
       write(18,'(I6,F20.10)')IUSER(NIU23),ENERGIA
       do jat=1,nat
        j = (jat-1)*3
        IZNUC = INT(USER(NU1+jat-1))+IUSER(NIU32+jat-1)
        write(18,'(1X,A4,3F15.4)')ATMLAB(IZNUC),                        &
        Cxyz(1+j)*BOHR,Cxyz(2+j)*BOHR,Cxyz(3+j)*BOHR
       enddo
      endif
      IF(MOLDEN==1 .AND. MOLDENGEO==1)                                 &
     & CALL WRITE_OPT_MLD_SNAPSHOT(IUSER(NIU23))
!-----------------------------------------------------------------------
      RETURN
      END

      SUBROUTINE WRITE_OPT_MLD_SNAPSHOT(ISTEP)
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: ISTEP
      CHARACTER(LEN=32) :: SNAPFILE
      CHARACTER(LEN=512) :: LINE
      INTEGER :: IOS
!-----------------------------------------------------------------------
!     Copy the current full Molden file (unit 17, file MLD) into
!     snapshot-####.mld after each accepted OPTGEO/TSOPT evaluation.
      FLUSH(17)
      WRITE(SNAPFILE,'(A9,I4.4,A4)') 'snapshot-',ISTEP,'.mld'
      OPEN(97,FILE=SNAPFILE,STATUS='UNKNOWN',FORM='FORMATTED',          &
           ACCESS='SEQUENTIAL')
      REWIND(17)
   10 CONTINUE
      READ(17,'(A)',IOSTAT=IOS) LINE
      IF(IOS /= 0) GOTO 20
      WRITE(97,'(A)') TRIM(LINE)
      GOTO 10
   20 CONTINUE
      CLOSE(97)
      REWIND(17)
!-----------------------------------------------------------------------
      RETURN
      END
      
!======================================================================!
