!======================================================================!
!                                                                      !
!                    N O F     S U B R O U T I N E S                   !
!                                                                      !
!======================================================================!
!======================================================================!
      
! RunNOF
      SUBROUTINE RunNOF(NATOMSn,NBFn,NBFTn,NSHELLn,NPRIMIn,ZAN,Cxyz,    &
                        IAN,IMIN,IMAX,KSTART,KATOM,KTYPE,KLOC,INTYP,    &
                        KNG,KMIN,KMAX,ISH,ITYP,C1,C2,EX1,CS,CP,CD,CF,   &
                        CG,CH,CI,AHCORE,OVERLAP,CHF,EiHF,DIPN,QUADN,    &
                        OCTUN,NVAL,DQOInt,NINTMXn,NREC,IX2,BUFP2,       &
                        BUFP2aux,NINTEGt,NINTEGAUXtm,IDONTW,GRADS,      &
                        IRUNTYP,DIPS,XINTS,SIZE_ENV,ENV,ATM,NBAS,BAS,   &
                        IGTYP,IPRINTOPT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      LOGICAL CONVGDELAG,RESTART,ERIACTIVATED,HFID,HighSpin,RHF
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/CONVERGENCE/DUMEL,PCONV,CONVGDELAG
      COMMON/ERIACT/ERIACTIVATED,NIJKaux,NINTCRaux,NSTOREaux,IAUXDIM
      COMMON/INPNOF_RHF/IRHFTYP,NCONVRHF,CONVRHFDM,MAXITRHF,RHF
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
      COMMON/INPNOF_RSTRT/RESTART,INPUTGAMMA,INPUTC,INPUTFMIUG,INPUTCXYZ
      COMMON/INPNOF_HFID/HFID,NTHRESHEID,THRESHEID,MAXITID,KOOPMANS
      COMMON/INPNOF_PRINT/NPRINT,IWRITEC,IMULPOP,IAIMPAC,IFCHK
      COMMON/INPNOF_GENERALINF/ICOEF,ISOFTMAX,IORBOPT,MAXIT,MAXIT21
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      COMMON/INPFILE_NIJKL/NINTMX,NIJKL,NINTCR,NSTORE
      COMMON/EHFEN/EHF,EN
      COMMON/ENERGIAS/EELEC,EELEC_OLD,DIF_EELEC,EELEC_MIN
      LOGICAL SMCD
      COMMON/ERITYPE/IERITYP,IRITYP,IGEN,ISTAR,MIXSTATE,SMCD
      COMMON/ELPROP/IEMOM
      COMMON/ECP2/CLP(4004),ZLP(4004),NLP(4004),KFRST(1001,6),          &
                  KLAST(1001,6),LMAX(1001),LPSKIP(1001),IZCORE(1001)
      COMMON/PUNTEROSUSER/N1,N2,N3,N4,N5,N6,N7,N8,N9,N10,N11,N12,N13,   &
                          N14,N15,N16,N17,N18,N19,N20,N21,N22,N23,N24,  &
                          N25,N26,N27,N28,N29,N30,N31,N32,N33,N34,N35,  &
                          N36,N37,N38,N39,N40,N41,N42,N43,N44,N45,N46,  &
                          N47,N48,N49,N50,N51,NUSER
!
      INTEGER :: NATOMSn,NBFn,NBFTn,NSHELLn,NPRIMIn,NVAL,NINTMXn,NREC
      INTEGER :: NINTEGt,NINTEGAUXtm,IDONTW,IRUNTYP,IPRINTOPT,IGTYP
      DOUBLE PRECISION,DIMENSION(NATOMSn):: ZAN
      DOUBLE PRECISION,DIMENSION(3,NATOMSn):: Cxyz
      INTEGER,DIMENSION(NATOMSn):: IAN,IMIN,IMAX
      INTEGER,DIMENSION(NSHELLn):: KSTART,KATOM,KTYPE,KLOC
      INTEGER,DIMENSION(NSHELLn):: INTYP,KNG,KMIN,KMAX
      INTEGER,DIMENSION(NPRIMIn):: ISH,ITYP
      DOUBLE PRECISION,DIMENSION(NPRIMIn):: C1,C2,EX1
      DOUBLE PRECISION,DIMENSION(NPRIMIn):: CS,CP,CD,CF,CG,CH,CI
      DOUBLE PRECISION,DIMENSION(NBFn):: EiHF
      DOUBLE PRECISION,DIMENSION(NBFn,NBFn):: AHCORE,OVERLAP,CHF
      DOUBLE PRECISION,DIMENSION(3):: DIPN
      DOUBLE PRECISION,DIMENSION(6):: QUADN
      DOUBLE PRECISION,DIMENSION(10):: OCTUN
      DOUBLE PRECISION,DIMENSION(NVAL*NBFTn):: DQOInt
      INTEGER,DIMENSION(NINTEGt) :: IX2
      DOUBLE PRECISION,DIMENSION(NINTEGt) :: BUFP2
      DOUBLE PRECISION,DIMENSION(NINTEGAUXtm) :: BUFP2aux      
      DOUBLE PRECISION,DIMENSION(3*NATOMSn) :: GRADS
      DOUBLE PRECISION,DIMENSION((NSHELL*NSHELL+NSHELL)/2) :: XINTS
      DOUBLE PRECISION,DIMENSION(3):: DIPS
!
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:):: COEF
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:):: XIJKL,XIJKaux,USER
      CHARACTER*4,ALLOCATABLE,DIMENSION(:)::ATMNAME
      INTEGER,ALLOCATABLE,DIMENSION(:)::LIMLOW,LIMSUP,IJKL
      LOGICAL CONVG,COEF21
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:):: GAMMA,FMIUG0
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::   ELAGN,RON
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:):: ELAG,COEFN
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:):: ADIPx,ADIPy,ADIPz
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:):: AQUADxx,AQUADyy
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:):: AQUADzz,AQUADxy
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:):: AQUADxz,AQUADyz
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:):: AOCTxxx,AOCTyyy
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:):: AOCTzzz,AOCTxxz
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:):: AOCTxyy,AOCTyyz
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:):: AOCTxzz,AOCTxxy
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:):: AOCTyzz,AOCTxyz
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:):: XATOM,YATOM,ZATOM
      INTEGER :: CR, CM, TIMESTART, TIMEFINISH
      DOUBLE PRECISION :: RATE

      INTEGER :: SIZE_ENV,NBAS               !LIBCINT
      DOUBLE PRECISION :: ENV(SIZE_ENV)      !LIBCINT
      INTEGER :: ATM(6,NATOMSn), BAS(8,NBAS) !LIBCINT
!-----------------------------------------------------------------------
!     Initialization for system_clock
!-----------------------------------------------------------------------
      CALL SYSTEM_CLOCK(COUNT_RATE=CR)
      CALL SYSTEM_CLOCK(COUNT_MAX=CM)
      RATE = REAL(CR)
!-----------------------------------------------------------------------
      ALLOCATE(XATOM(NATOMS),YATOM(NATOMS),ZATOM(NATOMS))
      XATOM(1:NATOMS) = Cxyz(1,1:NATOMS)
      YATOM(1:NATOMS) = Cxyz(2,1:NATOMS)
      ZATOM(1:NATOMS) = Cxyz(3,1:NATOMS)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(COEF(NBF,NBF))
      COEF = CHF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Calculate Nuclear Energy (EN)
!     ZAN: Nuclear charge array (1,NATOMS)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL NUCLEARm(NATOMS,ZAN,Cxyz,EN)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IAN(1) = IAN(1)   ! Avoiding warnings
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Name of Atoms
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(ATMNAME(NATOMS))
      CALL ATOMNAMES(NATOMS,ZAN,IZCORE,ATMNAME,Cxyz,NPRINT,0)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Atomic Basis Set
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       IMIN: Index of minimal primitive for atom
!       IMAX: Index of maximal primitive for atom
!     KSTART: Tells the location of the first exponent and the first
!             contraction coefficient contained in a particular shell
!      KATOM: Tells which atom the shell is centered on
!      KTYPE: is 1,2,3,4,5,6,7 for S,P,D,F,G,H,I. For L shell is 2
!       KLOC: Gives the location of the shell in the total AO basis 
!      INTYP: Index for the type of the shell
!        KNG: Number of Gaussians in the shell
!       KMIN: Starting index of the shell
!       KMAX: Ending index of the shell
!        ISH: Shell for the primitive
!       ITYP: Index for the type of the primitive
!         EX: Gaussian exponents
!         C1: S,P,D,F,G,H,I CONTRACTION COEFFICIENTS.
!         C2: Extra L CONTRACTION COEFFICIENTS.
!      CS-CI: S,P,D,F,G,H,I CONTRACTION COEFFICIENTS.
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IPRINTATOMBASIS=0
      IF(NPRINT==1.and.IPRINTATOMBASIS==1.and.IPRINTOPT==1)             &
      CALL ATOMBASIS(NATOMS,ATMNAME,IMIN,IMAX,NPRIMI,ITYP,ISH,EX1,C1,C2)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Create the basis function symbol table
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(LIMLOW(NATOMS),LIMSUP(NATOMS))
      CALL SYMBOLTABLE(KATOM,ATMNAME,INTYP,KLOC,LIMLOW,LIMSUP)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     NINTMX:  Number of 2e-integrals (ERIs) per record (15000)
!     NINTCHK: Integral of NINTMX in each core
!     NPROCS:  Number of cores
!     NCHUNKS: Number of chunks
!     NINTCR:  Space needed to allocate 2e- integrals in Slaves
!     NSTORE:  Space needed to allocate 2e- integrals in Master
!              NIJKL if serial, NINTCR if parallel
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      NINTMX = NINTMXn
      NIJKL  = NINTEGt
      NIJKaux  = NINTEGAUXtm
      CALL DISTRIBUTION(IPRINTOPT,IDONTW)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Read two-electron Repulsion Integrals in AO basis (ERI)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(IJKL(NSTORE),XIJKL(NSTORE),XIJKaux(NSTOREaux))
      IF(IERITYP==1) THEN
       CALL READERIs(IJKL,XIJKL,IX2,BUFP2,NINTEGt,IDONTW,NREC)
      ELSE IF(IERITYP==2 .or. IERITYP==3) THEN
       CALL READERIsAUX(XIJKaux,BUFP2aux,NINTEGAUXtm)
      END IF   
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Allocate User array
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(USER(NUSER))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Pass atomic dipole, quadrupole and octupole matrices to USER array
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(IEMOM>=1)THEN
       ALLOCATE(ADIPx(NSQ),ADIPy(NSQ),ADIPz(NSQ))
       CALL CPYTSQ(DQOInt(1       ),ADIPx,NBF)
       CALL CPYTSQ(DQOInt(1+  NBFT),ADIPy,NBF)
       CALL CPYTSQ(DQOInt(1+2*NBFT),ADIPz,NBF)
       CALL PASSDIPUSER(DIPN,ADIPx,ADIPy,ADIPz,USER)
       DEALLOCATE(ADIPx,ADIPy,ADIPz)
      END IF
!
      IF(IEMOM>=2)THEN
       ALLOCATE(AQUADxx(NSQ),AQUADyy(NSQ),AQUADzz(NSQ),                 &
                AQUADxy(NSQ),AQUADxz(NSQ),AQUADyz(NSQ))
       CALL CPYTSQ(DQOInt(1+3*NBFT),AQUADxx,NBF)
       CALL CPYTSQ(DQOInt(1+4*NBFT),AQUADyy,NBF)
       CALL CPYTSQ(DQOInt(1+5*NBFT),AQUADzz,NBF)
       CALL CPYTSQ(DQOInt(1+6*NBFT),AQUADxy,NBF)
       CALL CPYTSQ(DQOInt(1+7*NBFT),AQUADxz,NBF)
       CALL CPYTSQ(DQOInt(1+8*NBFT),AQUADyz,NBF)
       CALL PASSQUADUSER(QUADN,AQUADxx,AQUADyy,AQUADzz,                 &
                         AQUADxy,AQUADxz,AQUADyz,USER)
       DEALLOCATE(AQUADxx,AQUADyy,AQUADzz,AQUADxy,AQUADxz,AQUADyz)
      END IF
!
      IF(IEMOM==3)THEN
       ALLOCATE(AOCTxxx(NSQ),AOCTyyy(NSQ),AOCTzzz(NSQ),AOCTxxy(NSQ),    &
                AOCTxxz(NSQ),AOCTxyy(NSQ),AOCTyyz(NSQ),AOCTxzz(NSQ),    &
                AOCTyzz(NSQ),AOCTxyz(NSQ))
       CALL CPYTSQ(DQOInt(1+ 9*NBFT),AOCTxxx,NBF)
       CALL CPYTSQ(DQOInt(1+10*NBFT),AOCTyyy,NBF)
       CALL CPYTSQ(DQOInt(1+11*NBFT),AOCTzzz,NBF)
       CALL CPYTSQ(DQOInt(1+12*NBFT),AOCTxxy,NBF)
       CALL CPYTSQ(DQOInt(1+13*NBFT),AOCTxxz,NBF)
       CALL CPYTSQ(DQOInt(1+14*NBFT),AOCTxyy,NBF)
       CALL CPYTSQ(DQOInt(1+15*NBFT),AOCTyyz,NBF)
       CALL CPYTSQ(DQOInt(1+16*NBFT),AOCTxzz,NBF)
       CALL CPYTSQ(DQOInt(1+17*NBFT),AOCTyzz,NBF)
       CALL CPYTSQ(DQOInt(1+18*NBFT),AOCTxyz,NBF)
       CALL PASSOCTUSER(OCTUN,AOCTxxx,AOCTyyy,AOCTzzz,AOCTxxy,AOCTxxz,  &
                        AOCTxyy,AOCTyyz,AOCTxzz,AOCTyzz,AOCTxyz,USER)
       DEALLOCATE(AOCTxxx,AOCTyyy,AOCTzzz,AOCTxxz,AOCTxyy)
       DEALLOCATE(AOCTyyz,AOCTxzz,AOCTxxy,AOCTyzz,AOCTxyz)
      END IF
!***********************************************************************
!     FIRSTCALL: Initialize variables    COEF,      GAMMA,     FMIUG0
!                         according to INPUTC, INPUTGAMMA, INPUTFMUIG
!***********************************************************************
      ALLOCATE(GAMMA(NBF5),FMIUG0(NBF))
      CALL INITr(COEF,OVERLAP,GAMMA,FMIUG0,IPRINTOPT)
!=======================================================================      
!     Restricted Hartree-Fock (RHF)                                    
!     Use the Iterative Diagonalization Method to generate the HF MOs  
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
      IF(HFID)THEN
       if(INPUTC==1)CHF = COEF
       CALL HFIDr(AHCORE,IJKL,XIJKL,XIJKaux,CHF,EiHF,USER,IPRINTOPT)
       if(INPUTC==0)COEF = CHF
       if(INPUTFMIUG==0)FMIUG0 = EiHF
      ELSE
       if(INPUTC==0.and.IPRINTOPT==1)then
        if(IRHFTYP==0)then
         WRITE(6,11)
        else
         WRITE(6,12)
        end if
       end if
      ENDIF
!=======================================================================      
!     INITIALIZE LOCAL VARIABLES                                       
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IT=0
      ILOOP=0
      ITTOTAL=0
      IFIRSTCALL=0
      CONVG=.FALSE.         
      COEF21=.FALSE.
!=======================================================================
!               OPTIMIZATION WITH RESPECT TO THE OCCUPATIONS           
!=======================================================================
      ALLOCATE(ELAG(NBF,NBF),ELAGN(NBF),COEFN(NBF,NBF),RON(NBF))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Note: CONVGDELAG is the fundamental criterion in the optimization, 
!           so it must be FALSE before minimizing respect to GAMMAs and 
!           being able to call the CG subroutine. Its value is 
!           determined in the orbital optimization subroutine.
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CONVGDELAG = .FALSE.
      DIF_EELEC = 0.0d0      
      IF(IRUNTYP/=3) EELEC_MIN = 1.0d20 ! GLOBAL FIRST CALL
      IF(ICOEF==21)THEN
       ICOEFori=21
       ICOEF=2
       COEF21=.TRUE.
      ENDIF
      CALL SYSTEM_CLOCK(TIMESTART)
      CALL OccOpt(IFIRSTCALL,CONVG,ATMNAME,ZAN,OVERLAP,LIMLOW,LIMSUP,   &
                  COEF,GAMMA,FMIUG0,AHCORE,IJKL,XIJKL,XIJKaux,ELAG,     &
                  USER,IZCORE,XATOM,YATOM,ZATOM,KSTART,KNG,KMIN,KMAX,   &
                  KATOM,KTYPE,KLOC,IMIN,IMAX,ISH,ITYP,EX1,C1,C2,CS,     &
                  CP,CD,CF,CG,CH,CI,ELAGN,COEFN,RON,IT,ITTOTAL,DIPS,    &
                  IPRINTOPT,IRUNTYP)
      CALL SYSTEM_CLOCK(TIMEFINISH)
      OCCTIME = (TIMEFINISH - TIMESTART)/RATE
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     END SINGLE-POINT CALCULATION (ICOEF=0)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(ICOEF==0)THEN
       if(IERITYP==3 .and. MIXSTATE==1)then  ! MIX: Change to FULL ERIs
        CALL READERIs(IJKL,XIJKL,IX2,BUFP2,NINTEGt,IDONTW,NREC)
        MIXSTATE = 2
       end if
       CALL PT2GRAD(ELAG,COEF,AHCORE,IJKL,XIJKL,XIJKaux,USER,IRUNTYP,   &
                    GRADS,ATMNAME,KATOM,KTYPE,KLOC,KMIN,KMAX,KSTART,KNG,&
                    XATOM,YATOM,ZATOM,ZAN,EX1,CS,CP,CD,CF,CG,XINTS,     &
                    SIZE_ENV,ENV,NATOMSn,ATM,NBAS,BAS,IGTYP,IPRINTOPT)
       IF(IERITYP==3 .and. MIXSTATE==2) MIXSTATE = 1              ! MIX
       IF(IPRINTOPT==1)WRITE(6,2)IT
       GOTO 10
      ENDIF
!======================================================================!
!     SINGLE-POINT CALCULATION: FULL OPTIMIZATION FOR A GIVEN GEOMETRY !
!     Optimization with respect to the Occupations and Orbitals (COEF) !
!     using the iterative diagonalization method                       !
!======================================================================!
      IFIRSTCALL=1
      ITLIM=1
      DO WHILE(IT<=MAXIT)
       IT=IT+1
       if(COEF21.and.IT>MAXIT21)then    ! ICOEF21=ICOEF2+ICOEF1
        ICOEF=1
        COEF21=.FALSE.
       end if        
!      Orbital Optimization
       IF(ICOEF==1.or.ICOEF==2)THEN
        CALL OrbOpt(IT,ITLIM,OVERLAP,AHCORE,IJKL,XIJKL,XIJKaux,         &
                    USER(N7),COEF,USER(N1),USER(N2),USER(N3),ELAG,      &
                    FMIUG0,USER(N11),USER(N12),USER(N13),USER(N14),     &
                    ILOOP,IORBOPT,OCCTIME,IPRINTOPT)
!      Core-Fragment Orbital Optimization
       ELSEIF(ICOEF==3)THEN
        if(MSpin==0)then          ! Singlet and Multiplet States
         CALL FragOrbOpt(IT,ITLIM,AHCORE,IJKL,XIJKL,USER(N7),COEF,      &
                         USER(N1),USER(N2),USER(N3),ELAG,FMIUG0,        &
                         USER(N11),USER(N12),USER(N13),USER(N14),       &
                         ILOOP,IRUNTYP)
        else if(MSpin>0)then      ! High-Spin States
         WRITE(6,13)
         STOP
        end if      
       ENDIF
       CALL FLUSH(6)
!      MIX: Change to FULL ERIs if RI is converged
       IF(IERITYP==3 .and. MIXSTATE==1 .and. CONVGDELAG) THEN
        IF(IPRINTOPT==1)WRITE(6,14)
        CONVGDELAG = .FALSE.
        CALL READERIs(IJKL,XIJKL,IX2,BUFP2,NINTEGt,IDONTW,NREC)
        MIXSTATE = 2
       END IF
       
!      Occupation Optimization
       ITTOTAL=ITTOTAL+ILOOP
       CALL SYSTEM_CLOCK(TIMESTART)
       CALL OccOpt(IFIRSTCALL,CONVG,ATMNAME,ZAN,OVERLAP,LIMLOW,LIMSUP,  &
                   COEF,GAMMA,FMIUG0,AHCORE,IJKL,XIJKL,XIJKaux,ELAG,    &
                   USER,IZCORE,XATOM,YATOM,ZATOM,KSTART,KNG,KMIN,KMAX,  &
                   KATOM,KTYPE,KLOC,IMIN,IMAX,ISH,ITYP,EX1,C1,C2,CS,    &
                   CP,CD,CF,CG,CH,CI,ELAGN,COEFN,RON,IT,ITTOTAL,DIPS,   &
                   IPRINTOPT,IRUNTYP)
       CALL SYSTEM_CLOCK(TIMEFINISH)
       OCCTIME = (TIMEFINISH - TIMESTART)/RATE
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      END SINGLE-POINT CALCULATION (CONVG=TRUE)=(CONVGDELAG & ICOEF/=0)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       IF(CONVG)THEN
        CALL PT2GRAD(ELAG,COEF,AHCORE,IJKL,XIJKL,XIJKaux,USER,IRUNTYP,  &
                     GRADS,ATMNAME,KATOM,KTYPE,KLOC,KMIN,KMAX,          &
                     KSTART,KNG,XATOM,YATOM,ZATOM,ZAN,EX1,CS,CP,        &
                     CD,CF,CG,XINTS,SIZE_ENV,ENV,NATOMSn,ATM,NBAS,BAS,  &
                     IGTYP,IPRINTOPT)
!       MIX: Return to RI ERIs if FULL is converged                    
        IF(IERITYP==3 .and. MIXSTATE==2 .and. CONVGDELAG)MIXSTATE = 1
!                    
        IF(IPRINTOPT==0)GOTO 10
        WRITE(6,4)
        WRITE(6,5)IT,ITTOTAL
        GOTO 10
       ENDIF
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
      ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     EXCESSIVE NUMBER OF ITERATIONS
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(IT>MAXIT)THEN
       WRITE(6,3)
       WRITE(6,5)IT,ITTOTAL
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     STOP PROGRAM, DEALLOCATE MEMORY, GIVES ELAPSED TIME if IT > MAXIT
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   10 CONTINUE 
      IF(IRUNTYP==5.and.ICOEFori==21)ICOEF=21
      DEALLOCATE(COEF,ATMNAME,LIMLOW,LIMSUP,XATOM,YATOM,ZATOM)
      DEALLOCATE(GAMMA,FMIUG0,ELAG,ELAGN,RON,USER,COEFN)
      DEALLOCATE(IJKL,XIJKL,XIJKaux)
!-----------------------------------------------------------------------
!   1 FORMAT(/,'  ELAPSED REAL TIME :',F10.2,'  (SECONDS)')
    2 FORMAT(//2X,'**************************************************', &
              /2X,'*                                                *', &
              /2x,'*       SINGLE-POINT DoNOF CALCULATION           *', &
              /2X,'*                                                *', &
              /2X,'*            No.ITER =',I6,'                     *', &
              /2X,'*         (Occupation Optimization)              *', &
              /2X,'*                                                *', &
              /2x,'*  FINAL RESULTS   FINAL RESULTS  FINAL RESULTS  *', &
              /2X,'*                                                *', &
              /2X,'**************************************************')
    3 FORMAT(/,10X,30(1H-),/,10X,'EXCESSIVE NUMBER OF ITERATIONS',      &
             /,10X,30(1H-))
    4 FORMAT(//2X,'**************************************************', &
              /2X,'*                                                *', &
              /2x,'*       SINGLE-POINT DoNOF CALCULATION           *', &
              /2X,'*                                                *', &
              /2X,'*             FULL OPTIMIZATION                  *', &
              /2X,'*                                                *', &
              /2x,'*  FINAL RESULTS   FINAL RESULTS  FINAL RESULTS  *', &
              /2X,'*                                                *', &
              /2X,'**************************************************')
    5 FORMAT(/2X,'**************************************************',  &
             /2X,'*         No. EXTERNAL ITER =',I6,'              *',  &
             /2X,'*         No. of TOTAL ITER =',I6,'              *',  &
             /2X,'**************************************************')
   11 FORMAT(/1X,'Input for Coefficients is HCORE')
   12 FORMAT(/1X,'Input for Coefficients is RHF')
   13 FORMAT(/1X,'Fragment Calculations are not possible with MSpin>0')
   14 FORMAT(/1X,'Change to FULL ERIs')
!-----------------------------------------------------------------------
      RETURN
      END

! PT2GRAD
      SUBROUTINE PT2GRAD(ELAG,COEF,AHCORE,IJKL,XIJKL,XIJKaux,USER,      &
                         IRUNTYP,GRADS,ATMNAME,KATOM,KTYPE,KLOC,KMIN,   &
                         KMAX,KSTART,KNG,XATOM,YATOM,ZATOM,ZAN,EX1,CS,  &
                         CP,CD,CF,CG,XINTS,SIZE_ENV,ENV,NAT,ATM,NBAS,   &
                         BAS,IGTYP,IPRINTOPT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      LOGICAL ERPA,MBPT,OIMP2,SC2MCPT,HighSpin
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
      COMMON/INPNOF_ERPA/ERPA
      COMMON/INPNOF_OIMP2/OIMP2,MBPT
      COMMON/INPNOF_SC2MCPT/SC2MCPT
      LOGICAL SMCD
      COMMON/ERITYPE/IERITYP,IRITYP,IGEN,ISTAR,MIXSTATE,SMCD
      COMMON/INPFILE_NIJKL/NINTMX,NIJKL,NINTCR,NSTORE
      COMMON/ERIACT/ERIACTIVATED,NIJKaux,NINTCRaux,NSTOREaux,IAUXDIM
      COMMON/PUNTEROSUSER/N1,N2,N3,N4,N5,N6,N7,N8,N9,N10,N11,N12,N13,   &
                          N14,N15,N16,N17,N18,N19,N20,N21,N22,N23,N24,  &
                          N25,N26,N27,N28,N29,N30,N31,N32,N33,N34,N35,  &
                          N36,N37,N38,N39,N40,N41,N42,N43,N44,N45,N46,  &
                          N47,N48,N49,N50,N51,NUSER
!
      INTEGER :: IRUNTYP,IPRINTOPT
      DOUBLE PRECISION,DIMENSION(NBF,NBF) :: ELAG,COEF,AHCORE
      INTEGER,DIMENSION(NSTORE) :: IJKL  
      DOUBLE PRECISION,DIMENSION(NSTORE):: XIJKL
      DOUBLE PRECISION,DIMENSION(NSTOREaux):: XIJKaux
      DOUBLE PRECISION,DIMENSION(NUSER) :: USER 
      DOUBLE PRECISION,DIMENSION(3*NATOMS) :: GRADS      
      CHARACTER*4 ATMNAME(NATOMS)
      INTEGER,DIMENSION(NSHELL):: KATOM,KTYPE,KLOC,KMIN,KMAX,KSTART,KNG
      DOUBLE PRECISION,DIMENSION(NATOMS) :: XATOM,YATOM,ZATOM,ZAN 
      DOUBLE PRECISION,DIMENSION(NPRIMI) :: EX1,CS,CP,CD,CF,CG
      DOUBLE PRECISION,DIMENSION((NSHELL*NSHELL+NSHELL)/2) :: XINTS

      INTEGER :: SIZE_ENV,NBAS,IGTYP
      DOUBLE PRECISION :: ENV(SIZE_ENV)
      INTEGER :: ATM(6,NAT), BAS(8,NBAS)  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Orbital-Invariant MP2 Perturbative Corrections
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(OIMP2 .OR. MBPT)THEN
#ifdef MPI
       WRITE(6,1)
       CALL ABRT       
!      avoiding warnings    
       AHCORE(1,1) = AHCORE(1,1)
       IJKL(1) = IJKL(1)
       XIJKL(1) = XIJKL(1)
#else
       if(IERITYP==1 .or. IERITYP==3)then
        IF(OIMP2) THEN
          CALL ORBINVMP2(ELAG,COEF,USER(N1),USER(N2),USER(N3),AHCORE,IJKL,&
                       XIJKL,USER(N12),USER(N13),USER(N14))
        ENDIF
        IF(MBPT)THEN
          CALL MBPT_RPA(ELAG,COEF,USER(N1),USER(N2),USER(N3),AHCORE,IJKL,&
                       XIJKL,USER(N12),USER(N13),USER(N14))
        ENDIF
       else if(IERITYP==2)then
        WRITE(6,4)
        CALL ABRT  
       end if       
#endif                                             
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     SC2-MCPT (Hartree-Fock Partition)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(SC2MCPT)THEN
       if(NSOC==0)then
#ifdef MPI
        WRITE(6,2)
        CALL ABRT       
!       avoiding warnings    
        AHCORE(1,1) = AHCORE(1,1)
        IJKL(1) = IJKL(1)
        XIJKL(1) = XIJKL(1)
#else
        if(IERITYP==1 .or. IERITYP==3)then
         CALL SC2MCPThf(USER(N1),COEF,AHCORE,IJKL,XIJKL,USER(N10))
        else if(IERITYP==2)then
         WRITE(6,5)
         CALL ABRT  
        end if
#endif 
       else
        WRITE(6,7)      
        CALL ABRT        
       end if
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Extended Random Phase Approximation
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(ERPA)THEN
#ifdef MPI
       WRITE(6,3)
       CALL ABRT
!      avoiding warnings    
       AHCORE(1,1) = AHCORE(1,1)
       IJKL(1) = IJKL(1)
       XIJKL(1) = XIJKL(1)
#else
       if(IERITYP==1)then
        WRITE(6,6)
        CALL ABRT
       else if(IERITYP==2 .or. IERITYP==3)then
        CALL PNOF_ERPA(USER(N1),COEF,AHCORE,XIJKaux)
       end if
#endif                                             
      END IF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Analytical Gradient Calculation for PNOF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(IRUNTYP>1)THEN
       IF(.not.HighSpin)THEN

        CALL PNOFGRAD(COEF,USER(N7),USER(N1),ELAG,GRADS,ATMNAME,KATOM, &
                      KTYPE,KLOC,KMIN,KMAX,KSTART,KNG,XATOM,YATOM,     &
                      ZATOM,ZAN,EX1,CS,CP,CD,CF,CG,USER(N2),USER(N3),  &
                      XINTS,SIZE_ENV,ENV,NAT,ATM,NBAS,BAS,IGTYP,       &
                      IPRINTOPT)
       ELSE
        WRITE(6,8)
        CALL ABRT       
       END IF
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -   
#ifdef MPI
    1 FORMAT(/1X,'Sorry: OIMP2/MBPT is not implemented with MPI.',      &
            //1X,'Use serial code for perturbative calculations.',/)           
    2 FORMAT(/1X,'Sorry: SC2MCPT is not implemented with MPI.',         &
            //1X,'Use serial code for perturbative calculations.',/)  
    3 FORMAT(/1X,'Sorry: ERPA is not implemented with MPI.',            &
            //1X,'Use serial code for perturbative calculations.',/)  
#else
    4 FORMAT(/1X,'Sorry: OIMP2/MBPT is not implemented with ERITYP=RI.',&
            //1X,'Use ERITYP=FULL for perturbative calculations.',/)
    5 FORMAT(/1X,'Sorry: SC2MCPT is not implemented with ERITYP=RI.',   &
            //1X,'Use ERITYP=FULL for perturbative calculations.',/) 
    6 FORMAT(/1X,'Sorry: ERPA is not implemented with ERITYP=FULL.',    &
            //1X,'Use ERITYP=RI for ERPA calclations.',/) 
#endif
    7 FORMAT(/1X,'Sorry: NSOC>0, SC2MCPT is not implemented for',       &
              1X,'Spin-uncompensated systems',/) 
    8 FORMAT(/1X,'Sorry: GRADIENT is not implemented for High Spin',/)
!-----------------------------------------------------------------------
      RETURN
      END

!----------------------------------------------------------------------!
!                                                                      !
!                    Nuclear related Subroutines                       !
!                                                                      !
!   NUCLEARm: Calculate the nuclear energy                             !
!   DQONuclear: Nuclear Dipole, Quadrupole, Octupole Elec. Moments     !
!   ATOMNAMES: Write atom name-coordinates                             !
!                                                                      !
!----------------------------------------------------------------------!

! NUCLEARm
      SUBROUTINE NUCLEARm(NATOMS,ZNUC,Cxyz,Enuc)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DOUBLE PRECISION,DIMENSION(NATOMS) :: ZNUC
      DOUBLE PRECISION,DIMENSION(3,NATOMS) :: Cxyz
!-----------------------------------------------------------------------
      Enuc = 0.0
      DO I=1,NATOMS-1
       DO J=I+1,NATOMS
        DISTNUC = (Cxyz(1,I)-Cxyz(1,J))**2                              &
                + (Cxyz(2,I)-Cxyz(2,J))**2 + (Cxyz(3,I)-Cxyz(3,J))**2
        Enuc = Enuc + ZNUC(I)*ZNUC(J)/SQRT(DISTNUC)
       ENDDO
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! DQONuclear
      SUBROUTINE DQONuclear(DIPN,QUADN,OCTUN,Cxyz,ZAN,NAT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CMCoord/Xcm,Ycm,Zcm
      DOUBLE PRECISION,DIMENSION(3) :: DIPN
      DOUBLE PRECISION,DIMENSION(6) :: QUADN
      DOUBLE PRECISION,DIMENSION(10):: OCTUN
      DOUBLE PRECISION,DIMENSION(3,NAT) :: Cxyz
      DOUBLE PRECISION,DIMENSION(NAT) :: ZAN
!-----------------------------------------------------------------------
      DIPN  = 0.0d0                                                     
      QUADN = 0.0d0                                                       
      OCTUN = 0.0d0
      DO I=1,NAT
       XN = Cxyz(1,I) - Xcm
       YN = Cxyz(2,I) - Ycm
       ZN = Cxyz(3,I) - Zcm
       DIPN(1)   = DIPN(1)   + ZAN(I)*XN                                        
       DIPN(2)   = DIPN(2)   + ZAN(I)*YN                                        
       DIPN(3)   = DIPN(3)   + ZAN(I)*ZN                                        
       QUADN(1)  = QUADN(1)  + ZAN(I)*XN*XN                               
       QUADN(2)  = QUADN(2)  + ZAN(I)*YN*YN                               
       QUADN(3)  = QUADN(3)  + ZAN(I)*ZN*ZN                               
       QUADN(4)  = QUADN(4)  + ZAN(I)*XN*YN                               
       QUADN(5)  = QUADN(5)  + ZAN(I)*XN*ZN                               
       QUADN(6)  = QUADN(6)  + ZAN(I)*YN*ZN                               
       OCTUN(1)  = OCTUN(1)  + ZAN(I)*XN*XN*XN                            
       OCTUN(2)  = OCTUN(2)  + ZAN(I)*YN*YN*YN                            
       OCTUN(3)  = OCTUN(3)  + ZAN(I)*ZN*ZN*ZN                            
       OCTUN(4)  = OCTUN(4)  + ZAN(I)*XN*XN*YN                            
       OCTUN(5)  = OCTUN(5)  + ZAN(I)*XN*XN*ZN                            
       OCTUN(6)  = OCTUN(6)  + ZAN(I)*XN*YN*YN                            
       OCTUN(7)  = OCTUN(7)  + ZAN(I)*YN*YN*ZN                            
       OCTUN(8)  = OCTUN(8)  + ZAN(I)*XN*ZN*ZN                            
       OCTUN(9)  = OCTUN(9)  + ZAN(I)*YN*ZN*ZN                            
       OCTUN(10) = OCTUN(10) + ZAN(I)*XN*YN*ZN                          
      END DO                                                            
!-----------------------------------------------------------------------
      RETURN
      END
      
! ATOMNAMES
      SUBROUTINE ATOMNAMES(NATOMS,ZNUC,IZCORE,ATMNAME,Cxyz,NPRINT,      &
                           IWRITECXYZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      INTEGER:: IWRITECXYZ      
      CHARACTER*4 ATMNAME(NATOMS)
      DIMENSION ZNUC(NATOMS),IZCORE(NATOMS)
      DIMENSION Cxyz(3,NATOMS)
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
      DO I=1,NATOMS
       IZNUC = INT(ZNUC(I))+IZCORE(I)
       ATMNAME(I)=ATMLAB(IZNUC)
      ENDDO
      IF(NPRINT==1.and.IWRITECXYZ==1)THEN
       WRITE(6,1)
       DO I=1,NATOMS
        WRITE(6,2)ATMNAME(I),ZNUC(I)+IZCORE(I),                         &
                  Cxyz(1,I),Cxyz(2,I),Cxyz(3,I)
       ENDDO
      END IF
!-----------------------------------------------------------------------
    1 FORMAT(/1X,'Atom',6X,'Charge',16X,'Coordinates (Bohr)'/           &
              27X,'x',13X,'y',13X,'z')                     
    2 FORMAT(1X,A4,2X,F10.5,3F14.4)
!-----------------------------------------------------------------------
      RETURN
      END

!----------------------------------------------------------------------!

! ATOMBASIS
      SUBROUTINE ATOMBASIS(NATOMS,ATMNAME,IMIN,IMAX,NPRIMI,ITYP,ISH,    &
                           EX1,C1,C2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      CHARACTER*4 ATMNAME(NATOMS)
      DIMENSION IMIN(NATOMS),IMAX(NATOMS),ITYP(NPRIMI),ISH(NPRIMI)
      DIMENSION EX1(NPRIMI),C1(NPRIMI),C2(NPRIMI)
      CHARACTER*2 LABEL(8)
      DATA LABEL/'S ','P ','D ','F ','G ','H ','I ','L'/
!-----------------------------------------------------------------------
      WRITE(6,1)
!-----------------------------------------------------------------------
!     Atomic Basis Set
!
!     NSHELL: THE TOTAL NUMBER OF SHELLS.  
!             P SHELL MEANS X,Y,Z,
!             D SHELL MEANS XX,YY,ZZ,XY,XY,YZ, AND SO ON FOR F,G,H,I.
!     NPRIMI: TOTAL NUMBER OF PRIMITIVE EXPONENTS
!       IMIN: INDEX OF MINIMAL PRIMITIVE for ATOM
!       IMAX: INDEX OF MAXIMAL PRIMITIVE for ATOM
!       ITYP: INDEX for the TYPE of the PRIMITIVE
!        ISH: SHELL for the PRIMITIVE
!        EX1: GAUSSIAN EXPONENTS, FOR EVERY SYMMETRY UNIQUE PRIMITIVE.
!         C1: S,P,D,F,G,H,I CONTRACTION COEFFICIENTS.
!         C2: NORMALLY ONLY C1 ARRAYS WILL BE NON-ZERO,
!             THE EXCEPTION IS "L" SHELLS,
!             WHERE BOTH C1 AND C2 WILL HAVE DIFFERENT VALUES.
!-----------------------------------------------------------------------
      DO I=1,NATOMS
       WRITE(6,2)ATMNAME(I)
       DO J=IMIN(I),IMAX(I)
        IF(ITYP(J)<8)THEN
         WRITE(6,3)ISH(J),LABEL(ITYP(J)),J,EX1(J),C1(J)
        ELSE
         WRITE(6,3)ISH(J),LABEL(ITYP(J)),J,EX1(J),C1(J),C2(J)
        ENDIF
       ENDDO
      ENDDO
!-----------------------------------------------------------------------
    1 FORMAT(//5X,'Atomic Basis Set'/5X,16(1H-),//,2X,                  &
             'THE CONTRACTED BASIS FUNCTIONS ARE NORMALIZED TO UNITY',  &
             //2X,'SHELL TYPE',2X,'PRIMITIVE',8X,'EXPONENT',6X,         &
             'CONTRACTION COEFFICIENT(S)')
    2 FORMAT(/1X,A4,/)
    3 FORMAT(1X,I6,3X,A2,I7,F22.7,2F18.12)
!-----------------------------------------------------------------------
      RETURN
      END

! SYMBOLTABLE
      SUBROUTINE SYMBOLTABLE(KATOM,ATMNAME,INTYP,KLOC,LIMLOW,LIMSUP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      CHARACTER*4 LABELAT,ATMNAME(NATOMS),BFNAM1(35)
      CHARACTER*6 BFNAM2(49)
      DIMENSION KATOM(NSHELL),INTYP(NSHELL),KLOC(NSHELL)
      DIMENSION LIMLOW(NATOMS),LIMSUP(NATOMS),KMIN(8),KMAX(8)
      DIMENSION INTYPS(8,NATOMS)
      DATA KMIN /1,2, 5,11,21,34,57,1/
      DATA KMAX /1,4,10,20,35,56,84,4/
      DATA BFNAM1/'  S ','  X ','  Y ','  Z ',                          &
                  ' XX ',' YY ',' ZZ ',' XY ',' XZ ',' YZ ',            &
                  ' XXX',' YYY',' ZZZ',' XXY',' XXZ',                   &
                  ' YYX',' YYZ',' ZZX',' ZZY',' XYZ',                   &
                  'XXXX','YYYY','ZZZZ','XXXY','XXXZ',                   &
                  'YYYX','YYYZ','ZZZX','ZZZY','XXYY',                   &
                  'XXZZ','YYZZ','XXYZ','YYXZ','ZZXY'/
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
!-----------------------------------------------------------------------
!         ----- BASIS FUNCTION SYMBOL TABLE -----
!         S,  X,Y,Z,  XX,YY,ZZ,XY,XZ,YZ,
!         1   2 3 4    5  6  7  8  9 10
!         XXX,YYY,ZZZ,XXY,XXZ,YYX,YYZ,ZZX,ZZY,XYZ, ... G,H,I
!         11  12  13,  14  15  16  17  18  19  20, ... G,H,I
!
!      KATOM: TELLS WHICH ATOM THE SHELL IS CENTERED ON,
!             NORMALLY MORE THAN ONE SHELL EXISTS ON EVERY ATOM.
!
!      KLOC: Gives the location of the shell in the total AO basis 
!
!      INTYP: INDEX for the TYPE of the SHELL
!
!     KMIN AND KMAX ARE THE STARTING AND ENDING INDICES OF THE SHELL.  
!     THESE ARE DEFINED AS
!                    S    P    D    F   G   H   I   L
!             KMIN   1    2    5   11  21  34  57   1
!             KMAX   1    4   10   20  35  56  84   4
!-----------------------------------------------------------------------
      WRITE(4,'(I5)')NBF
      DO J = 1,NSHELL
       IAT = KATOM(J)
       LABELAT = ATMNAME(IAT)
       MINI = KMIN(INTYP(J))
       MAXI = KMAX(INTYP(J))
       DO I = MINI,MAXI
        IF(I<=35)THEN
         WRITE(4,'(A2,I2,A4)')LABELAT,MOD(IAT,1000),BFNAM1(I)
        ELSE
         WRITE(4,'(A2,A6)')LABELAT,BFNAM2(I-35)
        END IF
       ENDDO
      ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     INTYPS: 8xNATOMS matrix containing the dimension of each shell
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      INTYPS = 0
      LAT = 1
      ISUMTYP = 0
      ITYP0 = INTYP(1)
      NTYP = 1
!
      DO J = 1,NSHELL
       IAT = KATOM(J)
       IF(IAT==LAT)THEN
        IF(INTYP(J)==ITYP0)THEN
         ISUMTYP = ISUMTYP + 1
        ELSE
         ISUMTYP = 1
         ITYP0 = INTYP(J)
        ENDIF
       ELSE
        LAT = IAT
        ISUMTYP = 1
        ITYP0 = INTYP(J)
       ENDIF       
       INTYPS(ITYP0,LAT) = ISUMTYP
      ENDDO
!
      DO J = 1,NATOMS
       WRITE(4,'(9I5)')J,(INTYPS(I,J),I=1,8)
      ENDDO
!-----------------------------------------------------------------------
!     NUMBER OF BASIS FUNCTIONS PER ATOM: LIMLOW(iat) ... LIMSUP(iat)
!-----------------------------------------------------------------------
      IF(NATOMS>1000)THEN
       WRITE(6,*)'Stop, NATOMS>1000, enlarge the dimensions of limlow'
       CALL ABRT       
      ENDIF
      LAT=1
      J=1
! - - - - - - - - - - - - - -  
      LIMLOW(1)=1  
      DO I=1,NSHELL
       IAT=KATOM(I)
       IF(LAT/=IAT)THEN
        LAT=IAT
        LIMSUP(J)=KLOC(I)-1
        J=J+1
        LIMLOW(J)=KLOC(I)
       ENDIF
      ENDDO
! - - - - - - - - - - - - - -      
      LIMSUP(J)=NBF
      IF(J<NATOMS)THEN
       JP1=J+1
       DO J=JP1,NATOMS
        LIMLOW(J)=NBF
        LIMSUP(J)=1
       ENDDO
      ENDIF
!-----------------------------------------------------------------------
      RETURN
      END

! DISTRIBUTION      
      SUBROUTINE DISTRIBUTION(IPRINTOPT,IDONTW)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      INTEGER:: IPRINTOPT,IDONTW     
      LOGICAL ERIACTIVATED       
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/ERIACT/ERIACTIVATED,NIJKaux,NINTCRaux,NSTOREaux,IAUXDIM
      LOGICAL SMCD
      COMMON/ERITYPE/IERITYP,IRITYP,IGEN,ISTAR,MIXSTATE,SMCD
      COMMON/INPFILE_Naux/NBFaux,NSHELLaux
      COMMON/INPFILE_NIJKL/NINTMX,NIJKL,NINTCR,NSTORE
!-----------------------------------------------------------------------      
#include "mpip.h"
!-.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.-
!     NINTMX:  Number of 2e-integrals (ERIs) per record (15000)
!     NINTCHK: Integral of NINTMX in each core
!     NPROCS:  Number of cores
!     NCHUNKS: Number of chunks
!     NINTCR:  Dimension of ERIs in Slaves
!     NSTORE:  Dimension of ERIs in Master
!              NIJKL if serial, NINTCR if parallel
!-.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.-
#ifdef MPI

      if(IERITYP==1 .or. IERITYP==3)then

       IF(IDONTW==0)THEN      
        NINTCHK = INT(NINTMX/(2*NPROCS))*2+2
        NCHUNKS = INT(NIJKL/NINTMX)+1+10
        NINTCR = NINTCHK*NCHUNKS
        NSTORE = NINTCR+NINTMX
        IF(IPRINTOPT==1)THEN
         Write(6,501)
         Write(6,502)NIJKL
         Write(6,503)NINTCHK
         Write(6,504)NCHUNKS
         Write(6,505)NINTCR
         Write(6,506)NSTORE
        ENDIF
       ELSE    ! NINTMX = NINTEGtm, only 1 record
        NINTCR = INT(NIJKL/(2*NPROCS))*2+2
        NSTORE = NINTCR+NIJKL
        IF(IPRINTOPT==1)THEN
         Write(6,501)
         Write(6,502)NIJKL
         Write(6,505)NINTCR
         Write(6,506)NSTORE
        ENDIF
       ENDIF
      end if

      if(IERITYP==2 .or. IERITYP==3)then       

       NSTOREaux = INT(NBFaux/NPROCS)*NBF*(NBF+1)/2
       IAUXDIM = INT(NBFaux/NPROCS)
       IEXTRAS = MOD(NBFaux,NPROCS)
       IF(IPRINTOPT==1)THEN
        Write(6,511)
        Write(6,512)NIJKaux
        IF(IEXTRAS==0) THEN
          Write(6,513) NSTOREaux
        ELSE
          Write(6,513) NSTOREaux+NBF*(NBF+1)/2
        END IF
        Write(6,514)NSTOREaux
       ENDIF      

      end if    
      
!-----------------------------------------------------------------------
!     Format Statements Full Center ERIs
!-----------------------------------------------------------------------
  501 FORMAT(/,1X,'Distribution Full ERIs:',/)
  502 FORMAT(1X,'Integrals:                ',I18)
  503 FORMAT(1X,'Integrals/chunk:          ',I18)
  504 FORMAT(1X,'Chunks:                   ',I18)
  505 FORMAT(1X,'Max size per slaves:      ',I18)
  506 FORMAT(1X,'Max size in master:       ',I18)
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     Format Statements RI b
!-----------------------------------------------------------------------
  511 FORMAT(/,1X,'Distribution b tensor:',/)
  512 FORMAT(1X,'Elements:                 ',I18)
  513 FORMAT(1X,'Max size per slaves:      ',I18)
  514 FORMAT(1X,'Max size in master:       ',I18)
!-----------------------------------------------------------------------
      
#else
      IPRINTOPT1 = IPRINTOPT
      if(IERITYP==1 .or. IERITYP==3)then
       NINTCR = NIJKL
       NSTORE = NINTCR
      end if

      if(IERITYP==2 .or. IERITYP==3)then              
       NINTCRaux = NIJKaux
       NSTOREaux = NINTCRaux
       IAUXDIM = NBFaux
      end if
#endif

!-----------------------------------------------------------------------
      RETURN
      END

! READERIs
      SUBROUTINE READERIs(IERI,ERI,IX2,BUFP2,NINTEGt,IDONTW,NREC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      LOGICAL ERIACTIVATED       
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/ERIACT/ERIACTIVATED,NIJKaux,NINTCRaux,NSTOREaux,IAUXDIM
      COMMON/INPFILE_NIJKL/NINTMX,NIJKL,NINTCR,NSTORE
      COMMON/USEHUBBARD/IHUB
#include "mpip.h"
      INTEGER :: NINTEGt,IDONTW,NREC
      INTEGER,DIMENSION(NINTEGt) :: IX2
      DOUBLE PRECISION,DIMENSION(NINTEGt) :: BUFP2
      INTEGER,DIMENSION(NSTORE) :: IERI
      DOUBLE PRECISION,DIMENSION(NSTORE) :: ERI
      INTEGER,ALLOCATABLE,DIMENSION(:) :: IX
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: XX
#ifdef MPI
      LOGICAL TROUBLE
#endif
!-----------------------------------------------------------------------
      INTTYPE = 1
!      
      IF(IDONTW==0)ALLOCATE(IX(NINTMX),XX(NINTMX))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     NINTMX: Total Number of 2e-integrals per record (15000)
!     NIJKL:  Total Number of 2e-integrals
!     IERI:   I,J,K,L Indeces of the integrals
!     ERI:    Value of the integrals
!     NINTCR: Dimension of ERIs in Slaves
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#ifdef MPI
      TROUBLE=.FALSE.
      IF(ERIACTIVATED) THEN
        DO I=1,NPROCS-1
          CALL MPI_SEND(0,1,MPI_INTEGER8,I,I,MPI_COMM_WORLD,IERR)
          CALL MPI_SEND(0,1,MPI_INTEGER8,I,I,MPI_COMM_WORLD,IERR)
        ENDDO
      END IF
      CALL MPI_BCAST(INTTYPE,1,MPI_INTEGER8,MASTER,MPI_COMM_WORLD,IERR) 
      DO I=1,NPROCS-1
        CALL MPI_SEND(NINTCR,1,MPI_INTEGER8,I,I,MPI_COMM_WORLD,IERR)
!HUB        CALL MPI_SEND(IHUB,1,MPI_INTEGER8,I,I,MPI_COMM_WORLD,IERR)        
      END DO
#endif
!-.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.-
!     III  : Last position of the integrals in the Master
!     JJ   : Pointer of the integrals
!     ITOT : Total integrals read so far
!-.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.-
      III  = 0
      JJ   = 0
      ITOT = 0
      IRECORD = 0
      IF(IDONTW==0)REWIND(1)
!-.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.-
   10 CONTINUE
      IF(IDONTW==1)THEN
       IRECORD = IRECORD + 1
       IJBUFi = (IRECORD-1)*NINTMX
       if(IRECORD<NREC)then
        NNINT = NINTMX
        NX = NNINT
       else
        NNINT = NINTEGt - (NREC-1)*NINTMX
        NX = -NNINT
       end if
       IF(NX==0)GOTO 20
       DO M=1,NNINT
        I = M + III
        IERI(I) = IX2(IJBUFi+M)
        ERI(I) = BUFP2(IJBUFi+M)
       ENDDO
      ELSE
       READ(1,END=100)NX,IX,XX
       IF(NX==0)GOTO 20
       NNINT = IABS(NX)
       DO M=1,NNINT
        I = M + III
        IERI(I) = IX(M)
        ERI(I) = XX(M)
       ENDDO
      END IF

#ifdef MPI
!-.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.-
!      Distributing integrals into chunks of size 2
!-.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.-
       IPARITY = MOD(NNINT,2)
       NPAIRS = (NNINT+IPARITY)/2
!-.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.-
!      If non-even number of integrals in a chunk other than last
!      (and only the last, complain!)
!-.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.-
       IF (IPARITY/=0.AND.TROUBLE) THEN
        Write(6,*) 'Non-even number of integrals, problem!'
        Call Abortx('Non-even number of integrals, problem!')
       ENDIF
       IF (IPARITY/=0) TROUBLE=.TRUE.
       NLEFT=MOD(NPAIRS,NPROCS)
       NCHSZ=(NPAIRS-NLEFT)/NPROCS
       IF (NLEFT>0) NCHSZ=NCHSZ+1
!-.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.-
!      Sends chunks of integrals to the slaves
!-.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.-
       III = III+NCHSZ*2
       JJ  = III+1
       NPAIRS=NPAIRS-NCHSZ

       DO I=1,NPROCS-1
        NLEFT=MOD(NPAIRS,NPROCS-I)
        NCHSZ=(NPAIRS-NLEFT)/(NPROCS-I)
        IF (NLEFT>0) NCHSZ=NCHSZ+1
        NPAIRS=NPAIRS-NCHSZ
        NSIZE=NCHSZ*2
        IF (I.EQ.(NPROCS-1).AND.IPARITY.NE.0) NSIZE=NSIZE-1
        CALL MPI_SEND(NSIZE,1,MPI_INTEGER8,I,I,MPI_COMM_WORLD,IERR)
        CALL MPI_SEND(ERI(JJ),NSIZE,MPI_REAL8,I,I,MPI_COMM_WORLD,IERR)
        CALL MPI_SEND(IERI(JJ),NSIZE,MPI_INTEGER8,I,I,MPI_COMM_WORLD,   &
                      IERR)
        JJ=JJ+NSIZE
        ITOT=ITOT+NSIZE
       ENDDO
#else
       III=III+NNINT
#endif
!-.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.-
      IF(NX>0)GOTO 10
   20 CONTINUE
!-.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.-
#ifdef MPI
      DO I=1,NPROCS-1
       CALL MPI_SEND(-I,1,MPI_INTEGER8,I,I,MPI_COMM_WORLD,IERR)
      ENDDO
      ITOT=ITOT+III
      NINTCR=III
      ERIACTIVATED = .TRUE.       
      IF(ITOT/=NIJKL)THEN
       Call Abortx("Problem with integrals. Not everything shared?!?")
      ENDIF
#endif
!-----------------------------------------------------------------------
      IF(IDONTW==0)DEALLOCATE(XX,IX)
      RETURN
  100 CONTINUE
      WRITE(6,200)
  200 FORMAT(1X,'ERROR - ENCOUNTERED END OF FILE ON UNIT 1',I2)
      STOP
      END

! READERIsAUX
      SUBROUTINE READERIsAUX(ERIaux,BUFP2aux,NINTEGAUXtm)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      LOGICAL ERIACTIVATED       
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/ERIACT/ERIACTIVATED,NIJKaux,NINTCRaux,NSTOREaux,IAUXDIM
      COMMON/INPFILE_Naux/NBFaux,NSHELLaux
      COMMON/INPFILE_NIJKL/NINTMX,NIJKL,NINTCR,NSTORE
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
#include "mpip.h"
      DOUBLE PRECISION,DIMENSION(NINTEGAUXtm) :: BUFP2aux
      DOUBLE PRECISION,DIMENSION(NSTOREaux) :: ERIaux
!-----------------------------------------------------------------------
!     NINTMX: Total Number of 2e-integrals per record (15000)
!     NIJKL:  Total Number of 2e-integrals
!     IERI:   I,J,K,L Indeces of the integrals
!     ERI:    Value of the integrals
!     NINTCR: Dimension of ERIs in Slaves
!-----------------------------------------------------------------------
      INTTYPE = 2
#ifdef MPI
      IF(ERIACTIVATED) THEN
        DO I=1,NPROCS-1
          CALL MPI_SEND(0,1,MPI_INTEGER8,I,I,MPI_COMM_WORLD,IERR)
          CALL MPI_SEND(0,1,MPI_INTEGER8,I,I,MPI_COMM_WORLD,IERR)
        ENDDO
      END IF
      CALL MPI_BCAST(INTTYPE,1,MPI_INTEGER8,MASTER,MPI_COMM_WORLD,IERR)
      L_PER_PROC = INT(NBFaux/NPROCS)
      IEXTRAS = MOD(NBFaux,NPROCS)
      DO I=1,NPROCS-1
        NSIZE = L_PER_PROC
        IF (IEXTRAS>0) THEN
         NSIZE = NSIZE + 1
         IEXTRAS = IEXTRAS -1
        END IF
        NINTCRaux = NSIZE*NBF*(NBF+1)/2
        CALL MPI_SEND(NINTCRaux,1,MPI_INTEGER8,I,I,MPI_COMM_WORLD,IERR)
      END DO

      CALL MPI_BCAST(NBF,1,MPI_INTEGER8,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NBF5,1,MPI_INTEGER8,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NBFaux,1,MPI_INTEGER8,MASTER,MPI_COMM_WORLD,IERR)

      IEXTRAS = MOD(NBFaux,NPROCS)
      JJ=1
      DO I=1,NPROCS-1
        NSIZE = L_PER_PROC
        IF (IEXTRAS>0) THEN
         NSIZE = NSIZE + 1
         IEXTRAS = IEXTRAS - 1
        END IF
        NINTCRaux = NSIZE*NBF*(NBF+1)/2
        CALL MPI_SEND(BUFP2aux(JJ),NINTCRaux,MPI_REAL8,I,I,             &
                      MPI_COMM_WORLD,IERR)
        JJ = JJ + NINTCRaux
      END DO
      ERIaux = BUFP2aux(JJ:)
      ERIACTIVATED = .TRUE.       
#else
      ERIaux = BUFP2aux
#endif
!-----------------------------------------------------------------------
      RETURN
      END      

! DEACTIVATEERIs
      SUBROUTINE DEACTIVATEERIs
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL ERIACTIVATED       
      COMMON/ERIACT/ERIACTIVATED,NIJKaux,NINTCRaux,NSTOREaux,IAUXDIM
#include "mpip.h"
#ifdef MPI
      ERIACTIVATED = .FALSE.
#endif
      END      
      
!----------------------------------------------------------------------!
!                                                                      !
!                     Initialize Variables in RunNOF                   !
!                                                                      !
!----------------------------------------------------------------------!

! INITr
      SUBROUTINE INITr(COEF,OVERLAP,GAMMA,FMIUG0,IPRINTOPT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      LOGICAL CHKORTHO,ORTHO,RESTART
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_ORTHOGONALITY/CHKORTHO,ORTHO
      COMMON/INPNOF_RSTRT/RESTART,INPUTGAMMA,INPUTC,INPUTFMIUG,INPUTCXYZ
      COMMON/INPNOF_ORBSPACE0/NO1,NDOC,NCO,NCWO,NVIR,NAC,NO0
      COMMON/INPNOF_GENERALINF/ICOEF,ISOFTMAX,IORBOPT,MAXIT,MAXIT21
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5      
!
      INTEGER :: IPRINTOPT
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::COEF,OVERLAP
      DOUBLE PRECISION,DIMENSION(NBF)::FMIUG0
      DOUBLE PRECISION,DIMENSION(NBF5)::GAMMA
!-----------------------------------------------------------------------
!     INPUTGAMMA=0: Initial Values for GAMMA close to Fermi-Dirac dist.
!     INPUTGAMMA=1: Read ONs on file 3 (GCF) and transform to GAMMA
! -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -
      IF(INPUTGAMMA==0)THEN
       IF(ISOFTMAX==0)THEN      
        do i=1,ndoc
         GAMMA(i)= DACOS(DSQRT(2.0d0*0.999d0-1.0d0))
         do iw=1,ncwo-1
          !ig = ndoc+(i-1)*(ncwo-1)+iw !old-sort
          ig = ndoc*(iw+1)-i+1         !new-sort
          GAMMA(ig) = dasin(dsqrt(1.0d0/dfloat(ncwo-iw+1)))
         enddo
        enddo
       ELSE IF(ISOFTMAX==1)THEN
        do i=1,ndoc
         GAMMA(i)= DLOG(0.999d0)
         do iw=1,ncwo
          !ig = ndoc+(i-1)*ncwo+iw !old-sort
          ig = ndoc*(iw+1)-i+1     !new-sort
          GAMMA(ig) = DLOG(0.001d0/dfloat(ncwo))
         enddo
        enddo
       END IF
      ELSEIF(INPUTGAMMA==1)THEN
       CALL READGAMMAr(GAMMA)
      ENDIF
! -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -
!     INPUTC=0: HF Coeffcient Matrix (COEF=CHF) from HFIDr or HCORE
!     INPUTC=1: Reading on file GCF (3)
! -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -
      IF(INPUTC==1)THEN
       CALL READCOEFMr(COEF,NSQ,NBF)
!      Orthonormalization of the input orbitals (ORTHO=T)
       if(ORTHO)then
!       Check the orthonormality of the input orbitals
        CALL CHECKORTHO(COEF,OVERLAP,IVIOORTHO,IPRINTOPT)
!       Orthonormalize input orbitals if necessary
        if(IVIOORTHO/=0)then
         IF(IPRINTOPT==1)WRITE(6,'(A27)')' Orthogonalize the orbitals'
         CALL ORTHONORMAL(NBF,NBF,NBFT,OVERLAP,COEF,1,IPRINTOPT)
        endif
!      Check the Orthonormality of the input orbitals (CHKORTHO=T)
       elseif(CHKORTHO)then
        CALL CHECKORTHO(COEF,OVERLAP,IVIOORTHO,IPRINTOPT)
       endif
      ENDIF
! -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -
!     INPUTFMIUG=0: Diagonal Elements of F (FMIUG0) = Eigenvalues (E)
!     INPUTFMIUG=1: Read Diagonal Elements of F on file GCF (3)
! -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -
      FMIUG0 = 0.0D0
      IF(INPUTFMIUG==1)CALL READFMIUG0(FMIUG0,NBF,NSQ)
!-----------------------------------------------------------------------
      RETURN
      END
      
! READGAMMAr
      SUBROUTINE READGAMMAr(GAMMA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ      
      COMMON/INPNOF_ORBSPACE0/NO1,NDOC,NCO,NCWO,NVIR,NAC,NO0
      COMMON/INPNOF_GENERALINF/ICOEF,ISOFTMAX,IORBOPT,MAXIT,MAXIT21      
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
!
      DOUBLE PRECISION,DIMENSION(NBF5)::GAMMA
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::RO,HR
      ALLOCATE (RO(NBF5),HR(ndoc*(ncwo-1)))
!-----------------------------------------------------------------------
!     Read Occupations from the GCF file
!-----------------------------------------------------------------------
      REWIND(3)
      DO I=1,NBF5
       READ(3,'(I6,F30.16)')II,RO(I)
      ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Obtain GAMMAs (i=1,ndoc*ncwo=nv) from the Occupation Numbers (RO)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(ISOFTMAX==0)THEN      
       do i=1,NDOC
        in = NO1+i                                     ! in=no1+1,nb
        GAMMA(i) = dacos(dsqrt(2.0d0*RO(in)-1.0d0))
        IF(NCWO/=1)THEN
         !ici = (ncwo-1)*(i-1)+1        !old-sort
         !icf = (ncwo-1)*i              !old-sort
         !HRin = 1.0d0 - RO(in)         !old-sort
         !HR(ici:icf) = HRin            !old-sort
         ici = ndoc-i+1                 !new-sort
         icf = (ncwo-2)*ndoc + ndoc-i+1 !new-sort
         HRin = 1.0d0 - RO(in)          !new-sort
         do ic=ici,icf,ndoc             !new-sort
           HR(ic)  = HRin               !new-sort
         end do                         !new-sort
         do iw=1,ncwo-1
          !ic = (ncwo-1)*(i-1)+iw       !old-sort       ! ic=1,ndoc*(ncwo-1)
          !ig = ndoc+ic                 !old-sort       ! ig=ndoc+1,ndoc*ncwo
          !in = na+ncwo*(ndoc-i)+iw     !old-sort       ! in=na+1,na+ncwo*ndoc-1        
          ic = (iw-1)*ndoc + ndoc-i+1   !new-sort
          ig = ndoc+ic                  !new-sort
          in = no1+(na-nb)+ig           !new-sort
          if(HR(ic)>0.0d0)then
           ARGUM=sqrt(RO(in)/HR(ic))
           if(ARGUM>1.0d0)ARGUM=1.0d0
           GAMMA(ig)=asin(ARGUM)
          else 
           GAMMA(ig) = 0.0d0
          endif
          if(iw<ncwo-1)then
           do ix=1,ncwo-1-iw
            !ic1 = ic+ix                 !old-sort  !ic < ic1 < i*(ncwo-1)
            ic1 = ic+ix*ndoc             !new-sort
            HR(ic1) = HR(ic1) - RO(in)
           enddo
          endif
         enddo
        ENDIF
       enddo
      ELSE IF(ISOFTMAX==1)THEN       
       do i=1,NDOC
        in = NO1+i                                     ! in=no1+1,nb
        GAMMA(i) = dlog(RO(in))
        do iw=1,ncwo
         !ig = ndoc+ncwo*(i-1)+iw   !old-sort     ! ig=ndoc+1,ndoc*ncwo                
         !in = na+ncwo*(ndoc-i)+iw  !old-sort     ! in=na+1,na+ncwo*ndoc (nbf5)
         ig = ndoc*(iw+1)-i+1       !new-sort
         in = no1+(na-nb)+ig        !new-sort
         GAMMA(ig) = dlog(RO(in) + 1.0D-8)
        enddo
       enddo      
      END IF
!-----------------------------------------------------------------------
      DEALLOCATE (RO,HR)
      RETURN
      END

! READCOEFMr
      SUBROUTINE READCOEFMr(C,NSQ,NBF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DOUBLE PRECISION,DIMENSION(NSQ)::C
!-----------------------------------------------------------------------
      REWIND(3)
      DO I=1,NBF
       READ(3,'(I6,F30.16)')II,ROI
      ENDDO
      READ(3,'(F30.16)')SS
      DO I = 1,NSQ
       READ(3,'(I6,F30.16)')II,C(I)
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! READFMIUG0
      SUBROUTINE READFMIUG0(F,NBF,NSQ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DOUBLE PRECISION,DIMENSION(NBF)::F
!-----------------------------------------------------------------------
!     Read diagonal elements of the Gen Fock Operator (FMIUG)
!-----------------------------------------------------------------------
      REWIND(3)
      DO I = 1,NBF
       READ(3,'(I6,F30.16)')II,ROI
      ENDDO
      READ(3,'(F30.16)')SS
      DO I = 1,NSQ+NBF
       READ(3,'(I6,F30.16)')II,AA
      ENDDO
      DO I = 1,NBF
       READ(3,'(I6,F30.16)')II,F(I)
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! READCXYZ
      SUBROUTINE READCXYZ(ZNUC,C,NAT,NBF,NSQ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      INTEGER,DIMENSION(NAT):: IZNUC
      DOUBLE PRECISION,DIMENSION(NAT):: ZNUC
      DOUBLE PRECISION,DIMENSION(3,NAT):: C
      DOUBLE PRECISION,PARAMETER:: BOHR = 0.52917724924D+00
!-----------------------------------------------------------------------
!     Read diagonal elements of the Gen Fock Operator (FMIUG)
!-----------------------------------------------------------------------
      REWIND(3)
      DO I = 1,NBF
       READ(3,'(I6,F30.16)')II,AA
      ENDDO
      READ(3,'(F30.16)')SS
      DO I = 1,NSQ+NBF
       READ(3,'(I6,F30.16)')II,AA
      ENDDO
      DO I = 1,NBF
       READ(3,'(I6,F30.16)')II,AA
      ENDDO
      READ(3,'(I6,F30.16)')IT,AA
      READ(3,'(I6,F30.16)')IT,AA
      DO I = 1,NAT
       READ(3,'(I6,3F30.16)')IZNUC(I),C(1,I),C(2,I),C(3,I)
      ENDDO
      ZNUC = REAL(IZNUC)
      C = C / BOHR
!-----------------------------------------------------------------------
      RETURN
      END

! DENMATr
      SUBROUTINE DENMATr(DM,C,RO,NBF,N1,N2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DOUBLE PRECISION,DIMENSION(N2)::RO
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::DM,C
!----------------------------------------------------------------------
      do i=1,NBF
       do j=i,NBF
        DM(i,j) = 0.0d0
        DO K=N1,N2
         DM(i,j) = DM(i,j) + RO(K)*C(i,K)*C(j,K)
        ENDDO
        DM(j,i) = DM(i,j)
       enddo
      enddo
      DM = 2.0d0*DM
!----------------------------------------------------------------------
      RETURN
      END

! FORM2JK
      SUBROUTINE FORM2JK(FM,PM,IERI,ERI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPFILE_NIJKL/NINTMX,NIJKL,NINTCR,NSTORE
      COMMON/USEHUBBARD/IHUB      
#include "mpip.h"
      INTEGER,DIMENSION(NSTORE)::IERI
      DOUBLE PRECISION,DIMENSION(NSTORE)::ERI
      DOUBLE PRECISION,DIMENSION(NBF,NBF):: FM,PM
      ALLOCATABLE::P(:),F(:)
#ifdef MPI
      ALLOCATABLE::FF(:)
#endif
      ALLOCATE (P(NBFT),F(NBFT))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Wake up the nodes for the task
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#ifdef MPI
      ALLOCATE (FF(NBFT))
      DO I=1,NPROCS-1
       NOPT=3
       CALL MPI_SEND(NOPT,1,MPI_INTEGER8,I,I,MPI_COMM_WORLD,IERR)
       CALL MPI_SEND(NBFT,1,MPI_INTEGER8,I,I,MPI_COMM_WORLD,IERR)
      ENDDO
#endif
!-----------------------------------------------------------------------
      CALL SQUARETRIAN(PM,P,NBF,NBFT)
!-----------------------------------------------------------------------
      F = 0.0d0
#ifdef MPI
      FF = 0.0d0
      CALL MPI_BCAST(NBF,1,MPI_INTEGER8,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(P,NBFT,MPI_REAL8,MASTER,MPI_COMM_WORLD,IERR)
#endif
      !$OMP PARALLEL DO PRIVATE(LABEL, XJ, XK, I, J, K, L, NIJ, NKL, NIK, NJL, NIL, NJK) REDUCTION(+:F)
      DO M=1,NINTCR
       LABEL = IERI(M)
       CALL LABELIJKL(LABEL,I,J,K,L)
!-----------------------------------------------------------------------
!      2*J
!-----------------------------------------------------------------------
       XJ = ERI(M)
       NIJ = I*(I-1)/2 + J
       NKL = K*(K-1)/2 + L

       IF(IHUB==0)CALL OTTOINTEGR(I,J,K,L,NIJ,NKL,XJ)

                      F(NIJ)=F(NIJ)+P(NKL)*XJ
       IF(NIJ/=NKL)   F(NKL)=F(NKL)+P(NIJ)*XJ
!-----------------------------------------------------------------------
!      -K
!-----------------------------------------------------------------------
       XJ = 0.25*XJ
       XK = XJ
       NIK = I*(I-1)/2 + K
       NJL = MAX0(J,L)*(MAX0(J,L)-1)/2 + MIN0(J,L)
       IF(I==K.OR.J==L) XK=XK+XK
                          F(NIK)=F(NIK)-P(NJL)*XK
       IF(NIK/=NJL)       F(NJL)=F(NJL)-P(NIK)*XK
       IF(I/=J.and.K/=L)THEN
        NIL = I*(I-1)/2 + L
        NJK = MAX0(J,K)*(MAX0(J,K)-1)/2 + MIN0(J,K)
        IF(I==L.OR.J==K) XJ=XJ+XJ
                           F(NIL)=F(NIL)-P(NJK)*XJ
        IF(NIL/=NJK)       F(NJK)=F(NJK)-P(NIL)*XJ
       ENDIF
      ENDDO
      !$OMP END PARALLEL DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Get the pieces from slaves
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#ifdef MPI
      CALL MPI_REDUCE(F,FF,NBFT,MPI_REAL8,MPI_SUM,MASTER,               &
                      MPI_COMM_WORLD,IERR)
      CALL TRIANSQUARE(FM,FF,NBF,NBFT)
      DEALLOCATE(P,F,FF)
#else
      CALL TRIANSQUARE(FM,F,NBF,NBFT)
      DEALLOCATE(P,F)
#endif
!----------------------------------------------------------------------
      RETURN
      END
           
!----------------------------------------------------------------------!
!                                                                      !
!            O C C U P A T I O N    O P T I M I Z A T I O N            !
!                                                                      !
!             Spin-compensated Systems (Restricted Shells)             !
!                                                                      !
!    'Restricted Closed' (rc) Case: MSpin=0 [Singlet and Multiplet]    !
!       Singlet States (S=0,Ms=0) and Multiplet States (S>0,Ms=0)      !
!                                                                      !
!        'Restricted Open' (ro) Case: MSpin>0 [High-Spin State]        !
!                                                                      !
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
!                                                                      !
!     NO1:  Number of inactive doubly occupied orbitals (OCC=1)        !
!     NDOC: Number of strongly doubly occupied MOs                     !
!     NSOC: Number of strongly singly occupied MOs                     !
!     NDNS: Number of strongly occupied MOs (NDNS=NDOC+NSOC)           !
!     NCWO: Number of coupled weakly occ. MOs per strongly doubly occ. !
!     NCWO*NDOC: Active orbitals in the virtual subspace               !
!     NO0:  Empty orbitals  (OCC=0)                                    !
!     NVIR: Number of weakly occupied MOs + empty MOs                  !
!                                                                      !
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
!                                                                      !
!         NO1 |  NDOC + NSOC  |   NCWO*NDOC + NO0  = NBF               !                                                                
!         NO1 |      NDNS     |          NVIR      = NBF               !
!             | -NAC- |       |  -   NAC  - |                          !
!                    NB      NA            NBF5                        !
!                                                                      !
!----------------------------------------------------------------------!
!                                                                      !
!    OccOpt: Minimize the energy with respect to the occupations (RO)  !
!                                                                      !
!----------------------------------------------------------------------!

! OccOpt
      SUBROUTINE OccOpt(IFIRSTCALL,CONVG,ATMNAME,ZNUC,OVERLAP,LIMLOW,   &
                        LIMSUP,COEF,GAMMA,FMIUG0,AHCORE,IJKL,XIJKL,     &
                        XIJKaux,ELAG,USER,IZCORE,CX0,CY0,CZ0,KSTART,    &
                        KNG,KKMIN,KKMAX,KATOM,KTYPE,KLOC,IMIN,IMAX,     &
                        ISH,ITYP,EX1,C1,C2,CS,CP,CD,CF,CG,CH,CI,ELAGN,  &
                        COEFN,RON,IT,ITTOTAL,DIPS,IPRINTOPT,IRUNTYP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      LOGICAL CONVG,CONVGDELAG,ERIACTIVATED,HighSpin
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/CONVERGENCE/DUMEL,PCONV,CONVGDELAG
      COMMON/ERIACT/ERIACTIVATED,NIJKaux,NINTCRaux,NSTOREaux,IAUXDIM
      COMMON/INPNOF_ORBSPACE0/NO1,NDOC,NCO,NCWO,NVIR,NAC,NO0
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
      COMMON/INPNOF_EINI/IEINI
      COMMON/INPNOF_CGM/ICGMETHOD
      COMMON/INPNOF_PRINT/NPRINT,IWRITEC,IMULPOP,IAIMPAC,IFCHK
      COMMON/INPNOF_GENERALINF/ICOEF,ISOFTMAX,IORBOPT,MAXIT,MAXIT21
      COMMON/SUMSZ/SUMS,SUMF
      COMMON/INPFILE_NIJKL/NINTMX,NIJKL,NINTCR,NSTORE
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      COMMON/WRTGCF/IWRTGCF
      COMMON/EHFEN/EHF,EN
      COMMON/ENERGIAS/EELEC,EELEC_OLD,DIF_EELEC,EELEC_MIN
      COMMON/PUNTEROSUSER/N1,N2,N3,N4,N5,N6,N7,N8,N9,N10,N11,N12,N13,   &
                          N14,N15,N16,N17,N18,N19,N20,N21,N22,N23,N24,  &
                          N25,N26,N27,N28,N29,N30,N31,N32,N33,N34,N35,  &
                          N36,N37,N38,N39,N40,N41,N42,N43,N44,N45,N46,  &
                          N47,N48,N49,N50,N51,NUSER
      COMMON/INPNOF_COEFOPT/MAXLOOP
!
      INTEGER :: IFIRSTCALL,IT,ITTOTAL,IPRINTOPT,IRUNTYP
      CHARACTER*4,DIMENSION(NATOMS)::ATMNAME
      INTEGER,DIMENSION(NATOMS)::LIMLOW,LIMSUP,IZCORE,IMIN,IMAX
      INTEGER,DIMENSION(NPRIMI)::ISH,ITYP
      INTEGER,DIMENSION(NSHELL)::KSTART,KNG,KKMIN,KKMAX,KATOM,KTYPE,KLOC
      INTEGER,DIMENSION(NSTORE)::IJKL
      DOUBLE PRECISION,DIMENSION(NATOMS)::ZNUC,CX0,CY0,CZ0
      DOUBLE PRECISION,DIMENSION(NPRIMI)::EX1,C1,C2,CS,CP,CD,CF,CG,CH,CI
      DOUBLE PRECISION,DIMENSION(NBF)::FMIUG0,ELAGN,RON
      DOUBLE PRECISION,DIMENSION(NBF5)::GAMMA
      DOUBLE PRECISION,DIMENSION(NSTORE)::XIJKL
      DOUBLE PRECISION,DIMENSION(NSTOREaux)::XIJKaux
      DOUBLE PRECISION,DIMENSION(NUSER)::USER
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::OVERLAP,AHCORE
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::COEF,ELAG,COEFN
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::GAMMA_OLD,EAHF,E
      DOUBLE PRECISION,DIMENSION(3)::DIPS
      INTEGER,ALLOCATABLE,DIMENSION(:) :: IUSER
!-----------------------------------------------------------------------
!     Define the number of variables in the occupation optimization
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(ISOFTMAX==0)THEN
       NV = NDOC*NCWO
      ELSE IF(ISOFTMAX==1)THEN
       NV = NDOC*(NCWO+1)
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     IFIRSTCALL
!          = 0  First Call to OccOpt (Occupation Optimization)
!          = 1  Iterative Procedure for Occupations and Coefficients
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(GAMMA_OLD(NBF5))
      IF(IFIRSTCALL==0)THEN
       GAMMA_OLD = 0.0d0
      ELSEIF(IFIRSTCALL==1)THEN
       CALL XtoX0(GAMMA,GAMMA_OLD,NBF5)      ! Keep GAMMA in GAMMA_OLD
      ENDIF
!----------------------------------------------------------------------!
!                     Pointers of the USER array                       !
!----------------------------------------------------------------------!
!                                                                      !
!     N1  = 1                  ! USER( N1) = RO(NBF5)                  !
!     N2  = N1  + NBF5         ! USER( N2) = CJ12(NBF5,NBF5)           !
!     N3  = N2  + NSQ5         ! USER( N3) = CK12(NBF5,NBF5)           !
!     N4  = N3  + NSQ5         ! USER( N4) = DR(NBF5,NBF5)             !
!     N5  = N4  + NSQ5         ! USER( N5) = DCJ12r(NBF5,NBF5,NBF5)    !
!     N6  = N5  + NSQ5*NBF5    ! USER( N6) = DCK12r(NBF5,NBF5,NBF5)    !
!     N7  = N6  + NSQ5         ! USER( N7) = QD(NBF,NBF,NBF)           !
!     N8  = N7  + NBF*NSQ      ! USER( N8) = HCORE(NBF5)               !
!     N9  = N8  + NBF5         ! USER( N9) = QJ(NBFT5)                 !
!     N10 = N9  + NBFT5        ! USER(N10) = QK(NBFT5)                 !
!     N11 = N10 + NBFT5        ! USER(N11) = DIPN                      !
!     N12 = N11 + 3            ! USER(N12) = ADIPx                     !
!     N13 = N12 + NSQ          ! USER(N13) = ADIPy                     !
!     N14 = N13 + NSQ          ! USER(N14) = ADIPz                     !
!     N15 = N14 + NSQ          ! USER(N15) = DIPx                      !
!     N16 = N15 + NSQ5         ! USER(N16) = DIPy                      !
!     N17 = N16 + NSQ5         ! USER(N17) = DIPz                      !
!     N18 = N17 + NSQ5         ! USER(N18) = QUADN(6)                  !
!     N19 = N18 + 6            ! USER(N19) = AQUADxx(NSQ)              !
!     N20 = N19 + NSQ          ! USER(N20) = AQUADyy(NSQ)              !
!     N21 = N20 + NSQ          ! USER(N21) = AQUADzz(NSQ)              !
!     N22 = N21 + NSQ          ! USER(N22) = AQUADxy(NSQ)              !
!     N23 = N22 + NSQ          ! USER(N23) = AQUADxz(NSQ)              !
!     N24 = N23 + NSQ          ! USER(N24) = AQUADyz(NSQ)              !
!     N25 = N24 + NSQ          ! USER(N25) = QUADxx(NSQ5)              !
!     N26 = N25 + NSQ5         ! USER(N26) = QUADyy(NSQ5)              !
!     N27 = N26 + NSQ5         ! USER(N27) = QUADzz(NSQ5)              !
!     N28 = N27 + NSQ5         ! USER(N28) = QUADxy(NSQ5)              !
!     N29 = N28 + NSQ5         ! USER(N29) = QUADxz(NSQ5)              !
!     N30 = N29 + NSQ5         ! USER(N30) = QUADyz(NSQ5)              !
!     N31 = N30 + NSQ5         ! USER(N31) = OCTUN(10)                 !
!     N32 = N31 + 10           ! USER(N32) = AOCTxxx(NSQ)              !
!     N33 = N32 + NSQ          ! USER(N33) = AOCTyyy(NSQ)              !
!     N34 = N33 + NSQ          ! USER(N34) = AOCTzzz(NSQ)              !
!     N35 = N34 + NSQ          ! USER(N35) = AOCTxxy(NSQ)              !
!     N36 = N35 + NSQ          ! USER(N36) = AOCTxxz(NSQ)              !
!     N37 = N36 + NSQ          ! USER(N37) = AOCTxyy(NSQ)              !
!     N38 = N37 + NSQ          ! USER(N38) = AOCTyyz(NSQ)              !
!     N39 = N38 + NSQ          ! USER(N39) = AOCTxzz(NSQ)              !
!     N40 = N39 + NSQ          ! USER(N40) = AOCTyzz(NSQ)              !
!     N41 = N40 + NSQ          ! USER(N41) = AOCTxyz(NSQ)              !
!     N42 = N41 + NSQ          ! USER(N42) = OCTXXX(NSQ5)              !
!     N43 = N42 + NSQ5         ! USER(N43) = OCTYYY(NSQ5)              !
!     N44 = N43 + NSQ5         ! USER(N44) = OCTZZZ(NSQ5)              !
!     N45 = N44 + NSQ5         ! USER(N45) = OCTXXY(NSQ5)              !
!     N46 = N45 + NSQ5         ! USER(N46) = OCTXXZ(NSQ5)              !
!     N47 = N46 + NSQ5         ! USER(N47) = OCTXYY(NSQ5)              !
!     N48 = N47 + NSQ5         ! USER(N48) = OCTYYZ(NSQ5)              !
!     N49 = N48 + NSQ5         ! USER(N49) = OCTXZZ(NSQ5)              !
!     N50 = N49 + NSQ5         ! USER(N50) = OCTYZZ(NSQ5)              !
!     N51 = N50 + NSQ5         ! USER(N51) = OCTXYZ(NSQ5)              !
!                                                                      !
!----------------------------------------------------------------------!
!     Transform atomic integrals using COEF to obtain HCORE, QJ and QK 
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL FORMQHJK(COEF,USER(N7),USER(N8),USER(N9),USER(N10),AHCORE,   &
                    IJKL,XIJKL,XIJKaux)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Output of the MO inter-electronic integrals
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      NPRINTQJQK = 0
      IF(NPRINTQJQK==1)CALL OUTQJQK(USER(N9),USER(N10))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Output Energy and Properties with the input, and Stop if IEINI=1
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(IEINI==1 .or. NV==0 .or. (ICOEF==2.and.IFIRSTCALL==0) )THEN
       CALL OUTINITIALSr(NV,GAMMA,COEF,                                 &
            USER(N1),USER(N2),USER(N3),USER(N4),USER(N5),USER(N6),      &
            USER(N7),USER(N8),USER(N9),USER(N10),USER(N11),USER(N12),   &
            USER(N13),USER(N14),USER(N15),USER(N16),USER(N17),          &
            USER(N18),USER(N19),USER(N20),USER(N21),USER(N22),          &
            USER(N23),USER(N24),USER(N25),USER(N26),USER(N27),          &
            USER(N28),USER(N29),USER(N30),USER(N31),USER(N32),          &
            USER(N33),USER(N34),USER(N35),USER(N36),USER(N37),          &
            USER(N38),USER(N39),USER(N40),USER(N41),USER(N42),          &
            USER(N43),USER(N44),USER(N45),USER(N46),USER(N47),          &
            USER(N48),USER(N49),USER(N50),USER(N51),                    &
            ATMNAME,ZNUC,LIMLOW,LIMSUP,OVERLAP,IT,ITTOTAL,              &
            IPRINTOPT,IRUNTYP)
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Minimization of the total energy with respect to GAMMAs
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      EELEC_OLD = EELEC
      IF( .NOT.CONVGDELAG .and. NV>0 .and. ICOEF/=2 )THEN
       IF(ICGMETHOD==1)THEN
        CALL CGOCUPSUMSL(NV,GAMMA,USER,EELEC)
       ELSE IF(ICGMETHOD==2)THEN
        CALL CGOCUPNAG(NV,GAMMA,USER,EELEC)
       ELSE IF(ICGMETHOD==3)THEN       
        CALL LBFGSOCUP(NV,GAMMA,USER,EELEC)
       ENDIF
      END IF
      DIF_EELEC = EELEC - EELEC_OLD                          
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Restore GAMMAs and regenerate variables if energy rised
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(DIF_EELEC>0.0d0.and.IFIRSTCALL==1) THEN
        GAMMA = GAMMA_OLD
        ALLOCATE(IUSER(1))
        CALL CALCOE(NV,GAMMA,0,ENERGY,IUSER,USER)
        DEALLOCATE(IUSER)
      END IF 
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Increase Orbital Iterations if OccOpt is not decreasing energy
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(DIF_EELEC>-1.0D-6.and.IFIRSTCALL==1) THEN
        MAXLOOP = MAXLOOP + 10
      END IF 
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     One-particle Energies
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(EAHF(NBF),E(NBF))
      IF(IFIRSTCALL==0)THEN
!      Determine the 1-energies if first call and pass to ELAG
       CALL ENENEWr(USER(N1),USER(N8),USER(N9),USER(N10),USER(N2),      &
                    USER(N3),USER(N15),USER(N16),USER(N17),EAHF,E)
       ELAG = 0.0d0
       DO I=1,NBF
        ELAG(I,I) = E(I)
       ENDDO
      ELSEIF(IFIRSTCALL==1)THEN
!      The 1-energies are taken as eigenvalues of the MO optimization
       DO I=1,NBF
        E(I) = ELAG(I,I)
       ENDDO
       if(ICOEF==3)then
!       Determine the 1-energies if fragment calculation
        CALL ENENEWr(USER(N1),USER(N8),USER(N9),USER(N10),USER(N2),     &
                     USER(N3),USER(N15),USER(N16),USER(N17),EAHF,E)
        DO I=1,NBF
         ELAG(I,I) = E(I)
        ENDDO
       endif
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Calculate Dipole Moment if IPRINTOPT=0
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(IPRINTOPT==0)THEN
       CALL DIPMOMr(USER(N11),USER(N12),USER(N13),USER(N14),USER(N15),  &
                    USER(N16),USER(N17),USER(N7),USER(N1),              &
                    DMXe,DMYe,DMZe,DIPS(1),DIPS(2),DIPS(3),DM)  
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Initial Output (IFIRSTCALL==0), Intermediate Output (NPRINT=2)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF((IFIRSTCALL==0.or.(IFIRSTCALL==1.and.NPRINT==2)).and.          &
          IPRINTOPT==1)THEN                                              
       CALL PRINTOCr(E,COEF,ATMNAME,ZNUC,LIMLOW,LIMSUP,                 &
                     OVERLAP,USER(N1),USER(N2),USER(N3),USER(N7),       &
                     USER(N10),USER(N11),USER(N12),                     &
                     USER(N13),USER(N14),USER(N15),USER(N16),USER(N17), &
                     USER(N18),USER(N19),USER(N20),USER(N21),USER(N22), &
                     USER(N23),USER(N24),USER(N25),USER(N26),USER(N27), &
                     USER(N28),USER(N29),USER(N30),USER(N31),USER(N32), &
                     USER(N33),USER(N34),USER(N35),USER(N36),USER(N37), &
                     USER(N38),USER(N39),USER(N40),USER(N41),USER(N42), &
                     USER(N43),USER(N44),USER(N45),USER(N46),USER(N47), &
                     USER(N48),USER(N49),USER(N50),USER(N51),IZCORE,    &
                     CX0,CY0,CZ0,KSTART,KNG,KKMIN,KKMAX,KATOM,KTYPE,    &
                     KLOC,IMIN,IMAX,ISH,ITYP,EX1,C1,C2,CS,CP,CD,CF,     &
                     CG,CH,CI,IFIRSTCALL,DIPS,IRUNTYP)
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Ordering RO,COEF,E,FMIUG0 below NB
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(ISOFTMAX==1)THEN
       CALL ORDERSM(USER(N1),E,ELAG,COEF,FMIUG0,USER(N7),USER(N2),      &
                   USER(N3),USER(N8),USER(N9),USER(N10),USER(N15),      &
                   USER(N16),USER(N17),GAMMA,NV)
      ENDIF
      IF(IFIRSTCALL==0 .AND. IORBOPT==1)THEN
       CALL ORDERm(USER(N1),E,ELAG,COEF,FMIUG0,USER(N7),USER(N2),       &
                   USER(N3),USER(N8),USER(N9),USER(N10),USER(N15),      &
                   USER(N16),USER(N17),2)
      ENDIF
      IF(CONVGDELAG.and.ICOEF/=0)THEN
       CALL ORDERm(USER(N1),E,ELAG,COEF,FMIUG0,USER(N7),USER(N2),       &
                   USER(N3),USER(N8),USER(N9),USER(N10),USER(N15),      &
                   USER(N16),USER(N17),1) 
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Write RO,COEF,E,FMIUG0 on GCF file
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(IWRTGCF==1)then
       IF(IRUNTYP<5)THEN
        CALL WRITEGCFr(3,USER(N1),SUMS,COEF,E,FMIUG0,NSQ,NBF,NBF5,IT,   &
                       EELEC,EN,NO1,NDOC,NSOC,NCWO,NAC,NO0,ZNUC,        &
                       CX0,CY0,CZ0,NATOMS,1) 
        IF(EELEC<EELEC_MIN.and.IPRINTOPT==1)THEN                          
         CALL WRITEGCFr(8,USER(N1),SUMS,COEF,E,FMIUG0,NSQ,NBF,NBF5,IT,  &
                        EELEC,EN,NO1,NDOC,NSOC,NCWO,NAC,NO0,ZNUC,       &
                        CX0,CY0,CZ0,NATOMS,1)
         EELEC_MIN = EELEC
        ENDIF
       ELSE IF(IRUNTYP==5)THEN
        CALL WRITEGCFr(3,USER(N1),SUMS,COEF,E,FMIUG0,NSQ,NBF,NBF5,IT,   &
                       EELEC,EN,NO1,NDOC,NSOC,NCWO,NAC,NO0,ZNUC,        &
                       CX0,CY0,CZ0,NATOMS,1)                        
!                      CX0,CY0,CZ0,NATOMS,0) !Don't write RO on GCF file
       END IF
      end if
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Exit if convergence achieved
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(CONVGDELAG.and.ICOEF/=0)THEN
       CONVG=.TRUE.
       CALL FINALOUTPUTr(E,USER(N1),USER(N2),USER(N3),                  &
                         USER(N7),USER(N8),USER(N9),USER(N10),          &
                         ELAG,COEF,ATMNAME,ZNUC,LIMLOW,LIMSUP,          &
                         OVERLAP,USER(N11),USER(N12),USER(N13),         &
                         USER(N14),USER(N15),USER(N16),USER(N17),       &
                         USER(N18),USER(N19),USER(N20),USER(N21),       &
                         USER(N22),USER(N23),USER(N24),USER(N25),       &
                         USER(N26),USER(N27),USER(N28),USER(N29),       &
                         USER(N30),USER(N31),USER(N32),USER(N33),       &
                         USER(N34),USER(N35),USER(N36),USER(N37),       &
                         USER(N38),USER(N39),USER(N40),USER(N41),       &
                         USER(N42),USER(N43),USER(N44),USER(N45),       &
                         USER(N46),USER(N47),USER(N48),USER(N49),       &
                         USER(N50),USER(N51),IZCORE,CX0,CY0,CZ0,        &
                         KSTART,KNG,KKMIN,KKMAX,KATOM,KTYPE,KLOC,       &
                         IMIN,IMAX,ISH,ITYP,EX1,C1,C2,CS,CP,CD,CF,      &
                         CG,CH,CI,ELAGN,COEFN,RON,AHCORE,IJKL,XIJKL,    &
                         IFIRSTCALL,DIPS,IPRINTOPT,IRUNTYP)
      ENDIF
!-----------------------------------------------------------------------
      DEALLOCATE(GAMMA_OLD,EAHF,E)
      RETURN
      END

! ORDERm
      SUBROUTINE ORDERm(RO,E,ELAG,COEF,FMIUG0,QD,CJ12,CK12,HCORE,QJ,QK, &
                        DIPx,DIPy,DIPz,MODE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_ORBSPACE0/NO1,NDOC,NCO,NCWO,NVIR,NAC,NO0      
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      COMMON/ELPROP/IEMOM      
      DOUBLE PRECISION,DIMENSION(NBF5)::RO,HCORE,DIPx,DIPy,DIPz
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5)::CJ12,CK12
      DOUBLE PRECISION,DIMENSION(NBFT5)::QJ,QK
      DOUBLE PRECISION,DIMENSION(NBF)::E,FMIUG0                          
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::ELAG,COEF
      DOUBLE PRECISION,DIMENSION(NBF,NBF,NBF)::QD
!----------------------------------------------------------------------!
!     Ordering by Occupation Numbers (NO1+1:NO1+NDOC)=(NO1+1:NB)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      do k=1,ndoc
       kmini = k
       im = NO1+k  ! im=no1+1,nb
       MINI = im
       if(MODE==1)then ! by RO
        DUM = RO(im)
        do i=im,nb
         IF(RO(i)>DUM)THEN
          DUM = RO(i)
          MINI = i  ! MINI -> MAXI
          kmini = MINI - NO1
         ENDIF
        end do
       else if(MODE==2)then ! by FMIUG0
        DUM = FMIUG0(im)
        do i=im,nb
         IF(FMIUG0(i)<DUM)THEN
          DUM = FMIUG0(i)
          MINI = i
          kmini = MINI - NO1
         ENDIF
        end do
       end if
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
       IF(MINI/=im)THEN
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       Occupancies
        DUM1 = RO(im)
        RO(im) = RO(MINI)
        RO(MINI) = DUM1
!       Energies
        DUM2 = E(im)
        E(im) = E(MINI)
        E(MINI) = DUM2
!       Diagonal F-elements
        DUM3 = FMIUG0(im)
        FMIUG0(im) = FMIUG0(MINI)
        FMIUG0(MINI) = DUM3
!       MO Coefficients 
        do j=1,NBF
         DUM4 = COEF(j,im)
         COEF(j,im) = COEF(j,MINI)
         COEF(j,MINI) = DUM4
        end do
!       QD
        do j=1,NBF
         do i=1,NBF
          DUM5 = QD(im,i,j)
          QD(im,i,j)=QD(MINI,i,j)
          QD(MINI,i,j)=DUM5
         end do
        end do
!       CJ12,CK12 columns
        do j=1,NBF5
         DUM6 = CJ12(j,im)
         DUM7 = CK12(j,im)
         CJ12(j,im) = CJ12(j,MINI)
         CK12(j,im) = CK12(j,MINI)
         CJ12(j,MINI) = DUM6
         CK12(j,MINI) = DUM7
        end do
!       CJ12,CK12 rows
        do j=1,NBF5
         DUM8 = CJ12(im,j)
         DUM9 = CK12(im,j)
         CJ12(im,j) = CJ12(MINI,j)
         CK12(im,j) = CK12(MINI,j)
         CJ12(MINI,j) = DUM8
         CK12(MINI,j) = DUM9
        end do
!       HCORE
        DUM10 = HCORE(im)
        HCORE(im) = HCORE(MINI)
        HCORE(MINI) = DUM10
!       DIP (No IEMOM>1: QUAD,OCT)
        DUM11 = DIPx(im)
        DUM12 = DIPy(im)
        DUM13 = DIPz(im)
        DIPx(im) = DIPx(MINI)
        DIPy(im) = DIPy(MINI)
        DIPz(im) = DIPz(MINI)
        DIPx(MINI) = DUM11
        DIPy(MINI) = DUM12
        DIPz(MINI) = DUM13
!       QJ,QK
        do j=1,NBF5
         if(j<=im)then
          jim = j+im*(im-1)/2
         else
          jim = im+j*(j-1)/2
         endif
         if(j<=mini)then
          jmini = j+mini*(mini-1)/2
         else
          jmini = mini+j*(j-1)/2
         endif
         DUM14 = QJ(jim)
         DUM15 = QK(jim)
         QJ(jim) = QJ(jmini)
         QK(jim) = QK(jmini)
         QJ(jmini) = DUM14
         QK(jmini) = DUM15
        end do
!       Lagrangian Multiplier columns
        do j=1,NBF
         DUM16 = ELAG(j,im)
         ELAG(j,im) = ELAG(j,MINI)
         ELAG(j,MINI) = DUM16
        end do
!       Lagrangian Multiplier rows
        do j=1,NBF
         DUM17 = ELAG(im,j)
         ELAG(im,j) = ELAG(MINI,j)
         ELAG(MINI,j) = DUM17
        end do
!- - - - - - - - - - - - - - - - - - - - - - - -
!       in = na+1,na+ncwo*ndoc         
!- - - - - - - - - - - - - - - - - - - - - - - -
        do iw=1,ncwo
         !in = na+ncwo*(ndoc-k)+iw                 !old-sort
         !inmini = na+ncwo*(ndoc-kmini)+iw         !old-sort         
         in = no1+(na-nb)+ndoc*(iw+1)-k+1          !new-sort
         inmini = no1+(na-nb)+ndoc*(iw+1)-kmini+1  !new-sort
!        Occupancies         
         DUM1 = RO(in)
         RO(in) = RO(inmini)
         RO(inmini) = DUM1
!        Energies         
         DUM2 = E(in)
         E(in) = E(inmini)
         E(inmini) = DUM2
!        Diagonal F-elements
         DUM3 = FMIUG0(in)
         FMIUG0(in) = FMIUG0(inmini)
         FMIUG0(inmini) = DUM3
!        MO Coefficients 
         do j=1,NBF
          DUM4 = COEF(j,in)
          COEF(j,in) = COEF(j,inmini)
          COEF(j,inmini) = DUM4
         end do
!        QD
         do j=1,NBF
          do i=1,NBF
           DUM5 = QD(in,i,j)
           QD(in,i,j)=QD(inmini,i,j)
           QD(inmini,i,j)=DUM5
          end do
         end do
!        CJ12,CK12 columns
         do j=1,NBF5
          DUM6 = CJ12(j,in)
          DUM7 = CK12(j,in)
          CJ12(j,in) = CJ12(j,inmini)
          CK12(j,in) = CK12(j,inmini)
          CJ12(j,inmini) = DUM6
          CK12(j,inmini) = DUM7
         end do
!        CJ12,CK12 rows
         do j=1,NBF5
          DUM8 = CJ12(in,j)
          DUM9 = CK12(in,j)
          CJ12(in,j) = CJ12(inmini,j)
          CK12(in,j) = CK12(inmini,j)
          CJ12(inmini,j) = DUM8
          CK12(inmini,j) = DUM9
         end do
!        HCORE
         DUM10 = HCORE(in)
         HCORE(in) = HCORE(inmini)
         HCORE(inmini) = DUM10
!        DIP (No IEMOM>1: QUAD,OCT)
         DUM11 = DIPx(in)
         DUM12 = DIPy(in)
         DUM13 = DIPz(in)
         DIPx(in) = DIPx(inmini)
         DIPy(in) = DIPy(inmini)
         DIPz(in) = DIPz(inmini)
         DIPx(inmini) = DUM11
         DIPy(inmini) = DUM12
         DIPz(inmini) = DUM13
!        QJ,QK
         do j=1,NBF5
          if(j<=in)then
           jin = j+in*(in-1)/2
          else
           jin = in+j*(j-1)/2
          endif
          if(j<=inmini)then
           jinmini = j+inmini*(inmini-1)/2
          else
           jinmini = inmini+j*(j-1)/2
          endif
          DUM14 = QJ(jin)
          DUM15 = QK(jin)
          QJ(jin) = QJ(jinmini)
          QK(jin) = QK(jinmini)
          QJ(jinmini) = DUM14
          QK(jinmini) = DUM15
         end do
!        Lagrangian Multiplier columns
         do j=1,NBF
          DUM16 = ELAG(j,in)
          ELAG(j,in) = ELAG(j,inmini)
          ELAG(j,inmini) = DUM16
         end do
!        Lagrangian Multiplier rows
         do j=1,NBF
          DUM17 = ELAG(in,j)
          ELAG(in,j) = ELAG(inmini,j)
          ELAG(inmini,j) = DUM17
         end do
        end do
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       END IF
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
      end do
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END

! ORDERSM
      SUBROUTINE ORDERSM(RO,E,ELAG,COEF,FMIUG0,QD,CJ12,CK12,HCORE,QJ,QK, &
                        DIPx,DIPy,DIPz,GAMMA,NV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_ORBSPACE0/NO1,NDOC,NCO,NCWO,NVIR,NAC,NO0      
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      COMMON/ELPROP/IEMOM      
      DOUBLE PRECISION,DIMENSION(NV)::GAMMA
      DOUBLE PRECISION,DIMENSION(NBF5)::RO,HCORE,DIPx,DIPy,DIPz
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5)::CJ12,CK12
      DOUBLE PRECISION,DIMENSION(NBFT5)::QJ,QK
      DOUBLE PRECISION,DIMENSION(NBF)::E,FMIUG0                          
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::ELAG,COEF
      DOUBLE PRECISION,DIMENSION(NBF,NBF,NBF)::QD
!----------------------------------------------------------------------!
!     Make sure the strongly occupied orbital is the one with the highest
!     occupancy in the subspace 
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      do k=1,ndoc
       im = NO1+k  ! im=no1+1,nb
       MINI = im
       DUM = RO(im)
       kw = 0
       do iw=1,ncwo
        !in = na+ncwo*(ndoc-k)+iw         !old-sort
        in = no1+(na-nb)+ndoc*(iw+1)-k+1  !new-sort
        IF(RO(in)>DUM)THEN
         DUM = RO(in)
         MINI = in  ! MINI -> MAXI
         kw = iw
        ENDIF
       end do
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
       IF(MINI/=im)THEN
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       Gamma
        !ig = ndoc+ncwo*(k-1)+kw  !old-sort
        ig = ndoc*(kw+1)-k+1      !new-sort
        DUM0 = GAMMA(k)
        GAMMA(k) = GAMMA(ig)
        GAMMA(ig) = DUM0
!       Occupancies
        DUM1 = RO(im)
        RO(im) = RO(MINI)
        RO(MINI) = DUM1
!       Energies
        DUM2 = E(im)
        E(im) = E(MINI)
        E(MINI) = DUM2
!       Diagonal F-elements
        DUM3 = FMIUG0(im)
        FMIUG0(im) = FMIUG0(MINI)
        FMIUG0(MINI) = DUM3
!       MO Coefficients 
        do j=1,NBF
         DUM4 = COEF(j,im)
         COEF(j,im) = COEF(j,MINI)
         COEF(j,MINI) = DUM4
        end do
!       QD
        do j=1,NBF
         do i=1,NBF
          DUM5 = QD(im,i,j)
          QD(im,i,j)=QD(MINI,i,j)
          QD(MINI,i,j)=DUM5
         end do
        end do
!       CJ12,CK12 columns
        do j=1,NBF5
         DUM6 = CJ12(j,im)
         DUM7 = CK12(j,im)
         CJ12(j,im) = CJ12(j,MINI)
         CK12(j,im) = CK12(j,MINI)
         CJ12(j,MINI) = DUM6
         CK12(j,MINI) = DUM7
        end do
!       CJ12,CK12 rows
        do j=1,NBF5
         DUM8 = CJ12(im,j)
         DUM9 = CK12(im,j)
         CJ12(im,j) = CJ12(MINI,j)
         CK12(im,j) = CK12(MINI,j)
         CJ12(MINI,j) = DUM8
         CK12(MINI,j) = DUM9
        end do
!       HCORE
        DUM10 = HCORE(im)
        HCORE(im) = HCORE(MINI)
        HCORE(MINI) = DUM10
!       DIP (No IEMOM>1: QUAD,OCT)
        DUM11 = DIPx(im)
        DUM12 = DIPy(im)
        DUM13 = DIPz(im)
        DIPx(im) = DIPx(MINI)
        DIPy(im) = DIPy(MINI)
        DIPz(im) = DIPz(MINI)
        DIPx(MINI) = DUM11
        DIPy(MINI) = DUM12
        DIPz(MINI) = DUM13
!       QJ,QK
        do j=1,NBF5
         if(j<=im)then
          jim = j+im*(im-1)/2
         else
          jim = im+j*(j-1)/2
         endif
         if(j<=mini)then
          jmini = j+mini*(mini-1)/2
         else
          jmini = mini+j*(j-1)/2
         endif
         DUM14 = QJ(jim)
         DUM15 = QK(jim)
         QJ(jim) = QJ(jmini)
         QK(jim) = QK(jmini)
         QJ(jmini) = DUM14
         QK(jmini) = DUM15
        end do
!       Lagrangian Multiplier columns
        do j=1,NBF
         DUM16 = ELAG(j,im)
         ELAG(j,im) = ELAG(j,MINI)
         ELAG(j,MINI) = DUM16
        end do
!       Lagrangian Multiplier rows
        do j=1,NBF
         DUM17 = ELAG(im,j)
         ELAG(im,j) = ELAG(MINI,j)
         ELAG(MINI,j) = DUM17
        end do
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       END IF
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
      end do
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END
 
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
!                                                                      !
!       Subroutines for transformations from atomic into molecular     !
!             integrals: Hjj, <ij|ij>, <ij|ji>, <ii|jj>                !
!                                                                      !
!                      HCORE, QJ and QK MATRICES                       !
!                                                                      !
!    FORMQHJK: Density for each j (QDj), HCORE, QJ and QK matrices     !
!    OUTQHJK: Print molecular integrals, that is, QJ and QK matrices   !
!    QHMATm: Calculate Dj matrix keeping in QD(J,:,:) and molec. Hcore !
!    DENMATj: Obtain Density matix for each j (called from QHMATm)     !
!    QJMATm: Coulomb Integrals QJ(i,j)                                 !
!    QKMATm: Exchange Integrals QK(i,j)                                !
!    QJKMATmRI: Coulomb and Exchange Integrals with the RI Approx.     !
!    HSTARJ: Determine the skeleton J from atomic integrals (AUX)      !
!    HSTARK: Determine the skeleton K from atomic integrals (AUX)      !
!    HSTARJKRI: Determine the skeleton J and K from atomic integrals   !
!               with the RI approximation                              !
!                                                                      !
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !

! FORMQHJK
      SUBROUTINE FORMQHJK(COEF,QD,HCORE,QJ,QK,AHCORE,IJKL,XIJKL,XIJKaux)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      LOGICAL ERIACTIVATED       
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/ERIACT/ERIACTIVATED,NIJKaux,NINTCRaux,NSTOREaux,IAUXDIM
      COMMON/INPFILE_NIJKL/NINTMX,NIJKL,NINTCR,NSTORE
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      LOGICAL SMCD
      COMMON/ERITYPE/IERITYP,IRITYP,IGEN,ISTAR,MIXSTATE,SMCD
!
      DIMENSION COEF(NBF,NBF)
      DIMENSION QD(NBF,NBF,NBF),HCORE(NBF5),QJ(NBFT5),QK(NBFT5)
      DIMENSION AHCORE(NBF,NBF),IJKL(NSTORE),XIJKL(NSTORE)
      DOUBLE PRECISION,DIMENSION(NSTOREaux)::XIJKaux
!-----------------------------------------------------------------------
      QJ=0.0d0
      QK=0.0d0
      CALL QHMATm(COEF,QD,HCORE,AHCORE)
      IF(IERITYP==1 .or. (IERITYP==3 .and. MIXSTATE==2) )THEN
       CALL QJMATm(QD,QJ,IJKL,XIJKL)
       CALL QKMATm(QD,QK,IJKL,XIJKL)
      ELSE IF(IERITYP==2 .or. (IERITYP==3 .and. MIXSTATE==1) )THEN 
       CALL QJKMATmRI(QJ,QK,XIJKaux,COEF)
      END IF
!-----------------------------------------------------------------------
      RETURN
      END

! OUTQJQK
      SUBROUTINE OUTQJQK(QJ,QK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      DIMENSION QJ(NBFT5),QK(NBFT5)
!-----------------------------------------------------------------------
!     Print molecular integrals (called from OccOpt)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - - - - -
      WRITE(6,1)
      DO J=1,NBF5
       DO I=1,J
        IJ=I+J*(J-1)/2
        WRITE(6,2)I,J,QJ(IJ),QK(IJ)
       ENDDO
      ENDDO
      DO J=1,NBF5
       JJ=J*(J+1)/2
       WRITE(6,3)J,QJ(JJ)/2.0d0
      ENDDO
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - - - - -
    1 FORMAT(/,                                                         &
      1X,46('-'),/,6X,'I',3X,'J',4X,'QJ(IJ)',2X,'QK(IJ)',/,1X,46('-'),/)
    2 FORMAT(3X,2I4,1X,2F8.3)
    3 FORMAT(3X,I4,1X,F8.3)
!-----------------------------------------------------------------------
      RETURN
      END

! QHMATm
      SUBROUTINE QHMATm(C,QD,HCORE,AHCORE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)  
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
!
      DIMENSION C(NBF,NBF),QD(NBF,NBF,NBF),HCORE(NBF5),AHCORE(NBF,NBF)
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: AUX
!-----------------------------------------------------------------------
!     Calculate D matrix for each value J and keep in QD(J,:,:)
!     Calculate molecular Hcore matrix (HCORE)
!-----------------------------------------------------------------------
      ALLOCATE(AUX(NBF,NBF))
      !$OMP PARALLEL DO PRIVATE(J, AUX)
      DO J=1,NBF
       CALL DENMATj(J,AUX,C,NBF)
       QD(J,1:NBF,1:NBF) = AUX(1:NBF,1:NBF)
      ENDDO
      !$OMP END PARALLEL DO
      !$OMP PARALLEL DO PRIVATE(J, AUX)
      DO J=1,NBF5
       AUX(1:NBF,1:NBF) = QD(J,1:NBF,1:NBF)
       CALL TRACEm(HCORE(J),AUX,AHCORE,NBF)
      ENDDO
      !$OMP END PARALLEL DO
!-----------------------------------------------------------------------
      DEALLOCATE(AUX)
      RETURN
      END

! DENMATj
      SUBROUTINE DENMATj(J,D,C,NBF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DIMENSION D(NBF,NBF),C(NBF,NBF)
      DO M=1,NBF
       DO N=M,NBF
        D(M,N)=2.0d0*C(M,J)*C(N,J)
        D(N,M)=D(M,N)
       ENDDO
      ENDDO
      RETURN
      END

! QJMATm
      SUBROUTINE QJMATm(QD,QJ,IJKL,XIJKL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPFILE_NIJKL/NINTMX,NIJKL,NINTCR,NSTORE
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
!
      DIMENSION QD(NBF,NBF,NBF),QJ(NBFT5),IJKL(NSTORE),XIJKL(NSTORE)
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: AUX1,AUX2
      ALLOCATE (AUX1(NBF,NBF),AUX2(NBF,NBF))       ! AUX1:Jj(i,l)
!-----------------------------------------------------------------------
!     Coulomb Integrals QJ(i,j)
!-----------------------------------------------------------------------
      DO J=1,NBF5
       AUX2(1:NBF,1:NBF) = QD(J,1:NBF,1:NBF)
       CALL HSTARJ(AUX1,AUX2,IJKL,XIJKL)
       !$OMP PARALLEL DO PRIVATE(I, IJ, AUX2)
       DO I=1,J
        IJ=I+J*(J-1)/2
        AUX2(1:NBF,1:NBF) = QD(I,1:NBF,1:NBF)
        CALL TRACEm(QJ(IJ),AUX2,AUX1,NBF)
       ENDDO
       !$OMP END PARALLEL DO
      ENDDO
!-----------------------------------------------------------------------
      DEALLOCATE (AUX1,AUX2)
      RETURN
      END

! QKMATm
      SUBROUTINE QKMATm(QD,QK,IJKL,XIJKL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPFILE_NIJKL/NINTMX,NIJKL,NINTCR,NSTORE
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
!
      DIMENSION QD(NBF,NBF,NBF),QK(NBFT5),IJKL(NSTORE),XIJKL(NSTORE)
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: AUX1,AUX2
      ALLOCATE (AUX1(NBF,NBF),AUX2(NBF,NBF))       ! AUX1:Kj(i,l)
!-----------------------------------------------------------------------
!     Exchange Integrals QK(i,j)
!-----------------------------------------------------------------------
      DO J=1,NBF5
       AUX2(1:NBF,1:NBF) = QD(J,1:NBF,1:NBF)
       CALL HSTARK(AUX1,AUX2,IJKL,XIJKL)
       !$OMP PARALLEL DO PRIVATE(I, IJ, AUX2)
       DO I=1,J
        IJ=I+J*(J-1)/2
        AUX2(1:NBF,1:NBF) = QD(I,1:NBF,1:NBF)
        CALL TRACEm(QK(IJ),AUX2,AUX1,NBF)
       ENDDO
       !$OMP END PARALLEL DO
      ENDDO
!-----------------------------------------------------------------------
      DEALLOCATE (AUX1,AUX2)
      RETURN
      END

! QJKMATmRI
      SUBROUTINE QJKMATmRI(QJ,QK,XIJKAUX,C)
!     Coulomb and Exchange Integrals QJ(i,j) and QK(i,j)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL ERIACTIVATED       
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/ERIACT/ERIACTIVATED,NIJKaux,NINTCRaux,NSTOREaux,IAUXDIM
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
#include "mpip.h"
      DOUBLE PRECISION::QJ(NBFT5),QK(NBFT5)
      DOUBLE PRECISION::XIJKAUX(NBF*(NBF+1)/2,IAUXDIM)
#ifdef MPI
      DOUBLE PRECISION,DIMENSION(NBFT5)::QQJ,QQK
#endif
!-----------------------------------------------------------------------
      DOUBLE PRECISION :: C(NBF,NBF)
      DOUBLE PRECISION,ALLOCATABLE :: B_IN(:,:), B_IJ(:,:)
      DOUBLE PRECISION,ALLOCATABLE :: QJJ(:), QKK(:)
      INTEGER :: I,J,K,M,N
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Wake up the nodes for the task
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#ifdef MPI
      DO I=1,NPROCS-1
       NOPT = 2
       QQJ = 0.0d0
       QQK = 0.0d0
       CALL MPI_SEND(NOPT,1,MPI_INTEGER8,I,I,MPI_COMM_WORLD,IERR)
      ENDDO
      CALL MPI_BCAST(C,NBF*NBF,MPI_REAL8,MASTER,MPI_COMM_WORLD,IERR)
#endif
      ALLOCATE(B_IN(NBF5,NBF), B_IJ(NBF5,NBF5))
      ALLOCATE(QJJ(NBFT5),QKK(NBFT5))
      QJJ = 0.0D0
      QKK = 0.0D0
      !$OMP PARALLEL DO PRIVATE(K, M, N, MN, I, J, IJ, B_IN, B_IJ)      &
      !$OMP REDUCTION(+:QJJ,QKK)
      DO K=1,IAUXDIM

        B_IN(1:NBF5,1:NBF) = 0.0d0
        DO I=1,NBF5
          DO N=1,NBF
            DO M=1,N
              MN = M + N*(N-1)/2
              B_IN(I,N) = B_IN(I,N) + C(M,I)*XIJKAUX(MN,K)
              IF(M.NE.N) B_IN(I,M) = B_IN(I,M)+C(N,I)*XIJKAUX(MN,K)
            END DO
          END DO
        END DO

        B_IJ(1:NBF5,1:NBF5) = 0.0d0
        DO J=1,NBF5
          DO N=1,NBF
            DO I=1,J
              B_IJ(I,J) = B_IJ(I,J) + C(N,J)*B_IN(I,N)
            END DO
          END DO
        END DO

        DO J=1,NBF5
          DO I=1,J
            IJ=I+J*(J-1)/2
            QJJ(IJ) = QJJ(IJ) + B_IJ(I,I)*B_IJ(J,J)
          END DO
        END DO

        DO J=1,NBF5
          DO I=1,J
            IJ=I+J*(J-1)/2
            QKK(IJ) = QKK(IJ) + B_IJ(I,J)*B_IJ(I,J)
          END DO
        END DO

      END DO
      !$OMP END PARALLEL DO
      DEALLOCATE(B_IN, B_IJ)
      QJ(1:NBFT5) = QJJ(1:NBFT5)
      QK(1:NBFT5) = QKK(1:NBFT5)
      DEALLOCATE(QJJ,QKK)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Get the pieces from slaves
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#ifdef MPI
      CALL MPI_REDUCE(QJ,QQJ,NBFT5,MPI_REAL8,MPI_SUM,MASTER,            &
                      MPI_COMM_WORLD,IERR)                               
      CALL MPI_REDUCE(QK,QQK,NBFT5,MPI_REAL8,MPI_SUM,MASTER,            &
                      MPI_COMM_WORLD,IERR)
      QJ = QQJ
      QK = QQK
#endif      
!-----------------------------------------------------------------------
      RETURN
      END

! HSTARJ
      SUBROUTINE HSTARJ(FM,PM,IERI,ERI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPFILE_NIJKL/NINTMX,NIJKL,NINTCR,NSTORE
      COMMON/USEHUBBARD/IHUB      
#include "mpip.h"
      INTEGER,DIMENSION(NSTORE)::IERI
      DOUBLE PRECISION,DIMENSION(NSTORE)::ERI
      DOUBLE PRECISION,DIMENSION(NBF,NBF):: FM,PM
      ALLOCATABLE::P(:),F(:)
#ifdef MPI
      ALLOCATABLE::FF(:)
#endif
      ALLOCATE (P(NBFT),F(NBFT))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Wake up the nodes for the task
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#ifdef MPI
      ALLOCATE (FF(NBFT))
      DO I=1,NPROCS-1
       NOPT=1
       CALL MPI_SEND(NOPT,1,MPI_INTEGER8,I,I,MPI_COMM_WORLD,IERR)
       CALL MPI_SEND(NBFT,1,MPI_INTEGER8,I,I,MPI_COMM_WORLD,IERR)
      ENDDO
#endif
      CALL SQUARETRIAN(PM,P,NBF,NBFT)
      F = 0.0d0
#ifdef MPI
      FF = 0.0d0
      CALL MPI_BCAST(NBF,1,MPI_INTEGER8,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(P,NBFT,MPI_REAL8,MASTER,MPI_COMM_WORLD,IERR)
#endif
      !$OMP PARALLEL DO PRIVATE(LABEL, I, J, K, L, XJ, NIJ, NKL) REDUCTION(+:F)
      DO M=1,NINTCR
       LABEL = IERI(M)
       CALL LABELIJKL(LABEL,I,J,K,L)
       XJ = ERI(M)
       NIJ = I*(I-1)/2 + J
       NKL = K*(K-1)/2 + L

       IF(IHUB==0)CALL OTTOINTEGR(I,J,K,L,NIJ,NKL,XJ)

                       F(NIJ)=F(NIJ)+0.5*P(NKL)*XJ
       IF(NIJ/=NKL)    F(NKL)=F(NKL)+0.5*P(NIJ)*XJ
      ENDDO
      !$OMP END PARALLEL DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Get the pieces from slaves
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#ifdef MPI
      CALL MPI_REDUCE(F,FF,NBFT,MPI_REAL8,MPI_SUM,MASTER,               &
                      MPI_COMM_WORLD,IERR)
      CALL TRIANSQUARE(FM,FF,NBF,NBFT)
      DEALLOCATE(P,F,FF)
#else
      CALL TRIANSQUARE(FM,F,NBF,NBFT)
      DEALLOCATE(P,F)
#endif
!----------------------------------------------------------------------
      RETURN
      END

! HSTARK
      SUBROUTINE HSTARK(FM,PM,IERI,ERI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPFILE_NIJKL/NINTMX,NIJKL,NINTCR,NSTORE
      COMMON/USEHUBBARD/IHUB      
#include "mpip.h"
      INTEGER,DIMENSION(NSTORE)::IERI
      DOUBLE PRECISION,DIMENSION(NSTORE)::ERI
      DOUBLE PRECISION,DIMENSION(NBF,NBF):: FM,PM
      ALLOCATABLE::P(:),F(:)
#ifdef MPI
      ALLOCATABLE::FF(:)
#endif
      ALLOCATE (P(NBFT),F(NBFT))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Wake up the nodes for the task
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#ifdef MPI
      ALLOCATE (FF(NBFT))
      DO I=1,NPROCS-1
       NOPT=2
       CALL MPI_SEND(NOPT,1,MPI_INTEGER8,I,I,MPI_COMM_WORLD,IERR)
       CALL MPI_SEND(NBFT,1,MPI_INTEGER8,I,I,MPI_COMM_WORLD,IERR)
      ENDDO
#endif
      CALL SQUARETRIAN(PM,P,NBF,NBFT)
      F = 0.0d0
#ifdef MPI
      FF = 0.0d0
      CALL MPI_BCAST(NBF,1,MPI_INTEGER8,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(P,NBFT,MPI_REAL8,MASTER,MPI_COMM_WORLD,IERR)
#endif
      !$OMP PARALLEL DO PRIVATE(LABEL, I, J, K, L, XJ, XK, NIJ, NKL, NIK, NJL, NIL, NJK) REDUCTION(+:F)
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
                          F(NIK)=F(NIK)+P(NJL)*XK
       IF(NIK/=NJL)       F(NJL)=F(NJL)+P(NIK)*XK
       IF(I/=J.and.K/=L)THEN
        NIL = I*(I-1)/2 + L
        NJK = MAX0(J,K)*(MAX0(J,K)-1)/2 + MIN0(J,K)
        IF(I==L.OR.J==K) XJ=XJ+XJ
                           F(NIL)=F(NIL)+P(NJK)*XJ
        IF(NIL/=NJK)       F(NJK)=F(NJK)+P(NIL)*XJ
       ENDIF
      ENDDO
      !$OMP END PARALLEL DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Get the pieces from slaves
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#ifdef MPI
      CALL MPI_REDUCE(F,FF,NBFT,MPI_REAL8,MPI_SUM,MASTER,               &
                      MPI_COMM_WORLD,IERR)
      CALL TRIANSQUARE(FM,FF,NBF,NBFT)
      DEALLOCATE (P,F,FF)
#else
      CALL TRIANSQUARE(FM,F,NBF,NBFT)
      DEALLOCATE (P,F)
#endif
!----------------------------------------------------------------------
      RETURN
      END

! HSTARJKRI
      SUBROUTINE HSTARJKRI(FMJ,FMK,ERIaux,C)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL ERIACTIVATED
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/ERIACT/ERIACTIVATED,NIJKaux,NINTCRaux,NSTOREaux,IAUXDIM
#include "mpip.h"
      DOUBLE PRECISION,DIMENSION(NBF*(NBF+1)/2,IAUXDIM)::ERIaux
      DOUBLE PRECISION,DIMENSION(NBF,NBF):: FMJ,FMK
      DOUBLE PRECISION,ALLOCATABLE :: FJ(:),FK(:)
#ifdef MPI
      DOUBLE PRECISION,DIMENSION(NBFT)::FFJ,FFK
#endif
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DOUBLE PRECISION,DIMENSION(NBF) :: C
      DOUBLE PRECISION, ALLOCATABLE :: B_IN(:)
      DOUBLE PRECISION :: B_II
      INTEGER :: M,N,K
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Wake up the nodes for the task
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#ifdef MPI
      DO I=1,NPROCS-1
       NOPT=1
       FFJ = 0.0d0
       FFK = 0.0d0
       CALL MPI_SEND(NOPT,1,MPI_INTEGER8,I,I,MPI_COMM_WORLD,IERR)
      ENDDO
      CALL MPI_BCAST(C,NBF,MPI_REAL8,MASTER,MPI_COMM_WORLD,IERR)
#endif

      ALLOCATE(FJ(NBFT), FK(NBFT), B_IN(NBF))

      FJ(1:NBFT) = 0.0d0
      FK(1:NBFT) = 0.0d0

      !$OMP PARALLEL DO PRIVATE(K, M, N, MN, B_II, B_IN)                &
      !$OMP REDUCTION(+:FJ, FK)
      DO K=1,IAUXDIM
        B_IN(1:NBF) = 0.0d0
        B_II = 0.0d0

        DO N=1,NBF
          DO M=1,N
            MN = M + N*(N-1)/2
            B_IN(N) = B_IN(N) + C(M)*ERIaux(MN,K)
            IF(M.NE.N) B_IN(M) = B_IN(M)+C(N)*ERIaux(MN,K)
          END DO
        END DO

        DO N=1,NBF
          B_II = B_II + C(N)*B_IN(N)
        END DO

        DO N=1,NBF
          DO M=1,N
            MN = M + N*(N-1)/2
            FJ(MN) = FJ(MN) + B_II*ERIaux(MN,K)
          END DO
        END DO

        DO N=1,NBF
          DO M=1,N
            MN = M + N*(N-1)/2
            FK(MN) = FK(MN) + B_IN(M)*B_IN(N)
          END DO
        END DO
      END DO
      !$OMP END PARALLEL DO

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Get the pieces from slaves
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#ifdef MPI
      CALL MPI_REDUCE(FJ,FFJ,NBFT,MPI_REAL8,MPI_SUM,MASTER,             &
                      MPI_COMM_WORLD,IERR)
      CALL MPI_REDUCE(FK,FFK,NBFT,MPI_REAL8,MPI_SUM,MASTER,             &
                      MPI_COMM_WORLD,IERR)
      CALL TRIANSQUARE(FMJ,FFJ,NBF,NBFT)
      CALL TRIANSQUARE(FMK,FFK,NBF,NBFT)
#else      
      CALL TRIANSQUARE(FMJ,FJ,NBF,NBFT)
      CALL TRIANSQUARE(FMK,FK,NBF,NBFT)
#endif      
      DEALLOCATE(FJ, FK, B_IN)
!----------------------------------------------------------------------
      RETURN
      END

! LABELIJKL
      SUBROUTINE LABELIJKL(LABEL,I,J,K,L)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
!-----------------------------------------------------------------------
!     Determine label (ijkl) (2**16-1=65535)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      I = ISHFT( LABEL, -48 )                                 
      J = IAND( ISHFT( LABEL, -32 ), 65535 )                  
      K = IAND( ISHFT( LABEL, -16 ), 65535 )                  
      L = IAND( LABEL, 65535 )                                
!-----------------------------------------------------------------------
      RETURN
      END

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
!                                                                      !
!    Minimization of Elec. Energy using CG or LBFGS Algorithms         !
!                                                                      !
!    OCUPACIONTFr: Calculate ONs and their Derivatives using Trigomet. !
!    OCUPACIONSMr: Calculate ONs and their Derivatives using Softmax   !
!    ENENEWr: Evaluate one-particle energies                           !
!                                                                      !
!    CGOCUPNAG: Prepare for calling the NAG subroutine E04DGF          !
!    ENERFUNr: External Energy subroutine that calls OCUPENERGYrc,ro   !
!    OCUPENERGYrc,ro: Calculate the electronic energy and gradient     !
!    CGOCUPSUMSL: Prepare for calling the subroutine SUMSL             !
!    CALCOE: Compute the Occupation Energy                             !
!    CALCOG: Compute the Occupation Energy Gradients                   !
!    LBFGSOCUP: Optimize occupations by LBFGS algorithm                !
!                                                                      !
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !

! OCUPACIONTFr
      SUBROUTINE OCUPACIONTFr(GAMMA,RO,CJ12,CK12,DR,DCJ12r,DCK12r,NV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      LOGICAL HighSpin
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_ORBSPACE0/NO1,NDOC,NCO,NCWO,NVIR,NAC,NO0
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
      COMMON/INPNOF_PNOF/IPNOF,NTWOPAR
      COMMON/INPNOF_EXSTA/NESt,OMEGA1            
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      COMMON/SUMSZ/SUMS,SUMF
!
      DOUBLE PRECISION,DIMENSION(NV)::GAMMA       ! NV = NCWO*NDOC
      DOUBLE PRECISION,DIMENSION(NBF5)::RO
      DOUBLE PRECISION,DIMENSION(NBF5,NV)::DR
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5)::CJ12,CK12
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5,NV)::DCJ12r,DCK12r
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::DRO,BETA,DBETA,HR
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::DB,DHR
!-----------------------------------------------------------------------
!                 Occupations and their Derivatives
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE (DRO(NBF5),BETA(NBF5),DBETA(NBF5),DB(NBF5,NV))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RO = 0.0d0
      BETA = 0.0d0
      DRO = 0.0d0
      DBETA = 0.0d0
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Occupations (1,NO1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(NO1>0)THEN
       DO in=1,NO1
        RO(in) = 1.0d0
        BETA(in) = 1.0d0
       ENDDO
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Occupations (NO1+1,NO1+NDOC)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO i=1,NDOC-NESt
       in = NO1+i
       RO(in)   = 0.5d0 + 0.5d0*DCOS(GAMMA(i))*DCOS(GAMMA(i))
       DRO(in)  = - 0.5d0*DSIN(2.0d0*GAMMA(i))
       BETA(in) = DSQRT(RO(in))
       DBETA(in)= 0.5d0*DRO(in)/BETA(in)
      ENDDO
!
      if(NESt==1)then   ! First Excited State
       i = NDOC
       in = NO1+i       ! 1/2 <= RO(in) <= (1+OMEGA1)/2
       RO(in)   = 0.5d0 + 0.5d0*OMEGA1*DCOS(GAMMA(i))*DCOS(GAMMA(i))
       DRO(in)  = - 0.5d0*OMEGA1*DSIN(2.0d0*GAMMA(i))
       BETA(in) = DSQRT(RO(in))
       DBETA(in)= 0.5d0*DRO(in)/BETA(in)
      end if
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Occupations (NO1+NDOC+1,NA=NO1+NDNS)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(NSOC>0)THEN
       if(.not.HighSpin)then
        DO i=NDOC+1,NDNS
         in = NO1+i
         MULT1 = NA - NB
         RO(in)  = 0.5d0*MULT1/NSOC
         DRO(in) = 0.0d0
         BETA(in) = DSQRT(RO(in))
         DBETA(in) = 0.0d0
        ENDDO
       else if(HighSpin)then
        DO i=NDOC+1,NDNS
         in = NO1+i
         RO(in)  = 1.0d0
         DRO(in) = 0.0d0
        ENDDO
       end if      
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Occupations (NA+1,NBF5)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(NCWO==1)THEN             ! PNOFi(1): Perfect Pairing (NCWO=1)

       DR = 0.0d0
       DB = 0.0d0
       DO i=1,NDOC
        in = NO1+i                               ! in=no1+1,nb 
        DR(in,i) = DRO(in)
        DB(in,i) = DBETA(in)
        icf = na+ndoc-i+1      ! icf=na+ncwo*(ndoc-i)+ncwo with ncwo=1
        RO(icf)   = 1.0d0 - RO(in)
        DRO(icf)  = - DRO(in)
        DR(icf,i) = DRO(icf)
        BETA(icf) = DSQRT(RO(icf))
        if(BETA(icf)>0.0d0)then
         DBETA(icf)= 0.5d0*DRO(icf)/BETA(icf)
        else
         DBETA(icf)= 0.0d0
        end if
        DB(icf,i) = DBETA(icf)
       ENDDO

      ELSE                        ! PNOFi(Nc): Extended PNOF (NCWO>1)

       ALLOCATE (HR(NV-NDOC),DHR(NV-NDOC,NV))
       DR = 0.0d0
       DB = 0.0d0
       HR = 0.0d0
       DHR = 0.0d0
       
       DO i=1,NDOC                                ! ig=i=1,ndoc
        in = NO1+i                                ! in=no1+1,nb
        DR(in,i) = DRO(in)
        DB(in,i) = DBETA(in)
        !ici = (ncwo-1)*(i-1)+1         !old-sort
        !icf = (ncwo-1)*i               !old-sort
        ici = ndoc-i+1                  !new-sort
        icf = (ncwo-2)*ndoc + ndoc-i+1  !new-sort
        ! HR(ici:icf)  = 1.0d0 - RO(in) !old-sort
        !DHR(ici:icf,i)= - DRO(in)      !old-sort
        do ic=ici,icf,ndoc              !new-sort
          HR(ic)  = 1.0d0 - RO(in)      !new-sort
         DHR(ic,i)= - DRO(in)           !new-sort
        end do                          !new-sort
!- - - -- - - - - - - - - - (i,iw) <-> ic,ig,im  - - - - - - - - - - - -
        do iw=1,ncwo-1
         !ic = (ncwo-1)*(i-1)+iw             ! ic=1,ndoc*(ncwo-1)      !old-sort
         !ig = ndoc+ic                       ! ig=ndoc+1,ndoc*ncwo     !old-sort
         !im = na+ncwo*(ndoc-i)+iw           ! im=na+1,na+ncwo*ndoc-1  !old-sort
         ic = (iw-1)*ndoc + ndoc-i+1                                   !new-sort
         ig = ndoc+ic                                                  !new-sort
         im = no1+(na-nb)+ig                                           !new-sort
       
         ROn = DSIN(GAMMA(ig))*DSIN(GAMMA(ig))
         DRO(im) = DSIN(2.0d0*GAMMA(ig))
         BETAn = DSQRT(ROn)
         if(BETAn>0.0d0)then
          DBETA(im) = 0.5d0*DRO(im)/BETAn
         else
          DBETA(im) = 0.0d0
         end if
         RO(im) =       HR(ic)*ROn
         RAIZic = DSQRT(HR(ic))
         BETA(im) = RAIZic*BETAn
         DR(im,i) = DHR(ic,i)*ROn
         if(RAIZic>0.0d0)then
          DB(im,i) = 0.5d0*DHR(ic,i)*BETAn/RAIZic
         else 
          DB(im,i) = 0.0d0
         endif
         !do ic1=ici,ic-1                          !old-sort
         do ic1=ici,icf-ndoc,ndoc                  !new-sort
          ig1 = ndoc+ic1                           !   i < ig1 < ig
          DR(im,ig1) =        DHR(ic,ig1)*ROn
          if(RAIZic>0.0d0)then
           DB(im,ig1) = 0.5d0*DHR(ic,ig1)*BETAn/RAIZic
          else 
           DB(im,ig1) = 0.0d0
          endif
         enddo
         DR(im,ig) = HR(ic)*DRO(im)
         DB(im,ig) = RAIZic*DBETA(im)

!- - - - HR(ic+1) - - - - - - - - - - - - - - - -
         if(iw<ncwo-1)then
          do ix=1,ncwo-1-iw
           !ic1 = ic+ix                      !old-sort
           ic1 = ic+ix*ndoc                  !new-sort
            HR(ic1)  =  HR(ic1) - RO(im)
           DHR(ic1,i)= DHR(ic1,i) - DR(im,i)
           !do icn=ici,ic-1                  !old-sort
           do icn=ici,icf-ndoc,ndoc          !new-sort
            ign = ndoc+icn
            DHR(ic1,ign)= DHR(ic1,ign) - DR(im,ign)
           enddo
           DHR(ic1,ig)= DHR(ic1,ig) - DR(im,ig)
          enddo
         endif
!- - - - HR(ic+1) - - - - - - - - - - - - - - - -
        enddo
!- - - -- - - - - - - - - - (i,iw) <-> ic,ig,im  - - - - - - - - - - - -

!- - - - ic = icf - last RO  - - - - - - - - - - - - - -
        ig = ndoc+icf               ! ig=ndoc+i*(ncwo-1)
        !im = na+ncwo*(ndoc-i)+ncwo !old-sort
        im = no1+(na-nb)+ig+ndoc    !new-sort
        Hn = DCOS(GAMMA(ig))*DCOS(GAMMA(ig))
        DRO(im) = -DSIN(2.0d0*GAMMA(ig))
        BETAn = DSQRT(Hn)
        if(BETAn>0.0d0)then
         DBETA(im) = 0.5d0*DRO(im)/BETAn
        else
         DBETA(im) = 0.0d0
        end if       
        RO(im)  =       HR(icf)*Hn
        RAIZicf = DSQRT(HR(icf))
        BETA(im)= RAIZicf*BETAn

        DR(im,i) =        DHR(icf,i)*Hn
        if(RAIZicf>0.0d0)then
         DB(im,i) = 0.5d0*DHR(icf,i)*BETAn/RAIZicf
        else 
         DB(im,i) = 0.0d0
        endif
 
        !do ic1=ici,icf-1           !old-sort ! ici < ic1 < icf-1
        do ic1=ici,icf-ndoc,ndoc    !new-sort
         ig1 = ndoc+ic1             !   i < ig1 < ig
         DR(im,ig1) =        DHR(icf,ig1)*Hn
         if(RAIZicf>0.0d0)then
          DB(im,ig1) = 0.5d0*DHR(icf,ig1)*BETAn/RAIZicf
         else 
          DB(im,ig1) = 0.0d0
         endif
        enddo
        DR(im,ig) = HR(icf) *DRO(im)
        DB(im,ig) = RAIZicf*DBETA(im)
!- - - - ic = icf - last RO  - - - - - - - - - - - - - -
       ENDDO
       DEALLOCATE(HR,DHR)       

      ENDIF
      DEALLOCATE(DRO,DBETA)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Sum of the Holes below the Fermi Level (SUMS)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      SUMS = DFLOAT(NB)
      do j=1,nb
       SUMS = SUMS - RO(j)
      enddo
!-----------------------------------------------------------------------
!                   CJ12, CK12, DCJ12r, DCK12r
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(IPNOF==3)THEN
       if(NSOC==0)CALL CJCKD3(NV,RO,DR,BETA,DB,CJ12,CK12,DCJ12r,DCK12r)
      ELSEIF(IPNOF==4)THEN
       if(NSOC==0)CALL CJCKD4(NV,RO,DR,BETA,DB,CJ12,CK12,DCJ12r,DCK12r)
      ELSEIF(IPNOF==5)THEN
       CALL CJCKD5(NV,RO,DR,BETA,DB,CJ12,CK12,DCJ12r,DCK12r)
      ELSEIF(IPNOF==6)THEN
       if(NSOC==0)CALL CJCKD6(NV,RO,DR,CJ12,CK12,DCJ12r,DCK12r)
      ELSEIF(IPNOF==7)THEN
       CALL CJCKD7(NV,RO,DR,BETA,DB,CJ12,CK12,DCJ12r,DCK12r)
      ELSEIF(IPNOF==8)THEN
       CALL CJCKD8(NV,RO,DR,BETA,DB,CJ12,CK12,DCJ12r,DCK12r)
      ENDIF
!-----------------------------------------------------------------------
      DEALLOCATE(BETA,DB)
      RETURN
      END

! OCUPACIONSMr
      SUBROUTINE OCUPACIONSMr(GAMMA,RO,CJ12,CK12,DR,DCJ12r,DCK12r,NV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      LOGICAL HighSpin
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_ORBSPACE0/NO1,NDOC,NCO,NCWO,NVIR,NAC,NO0
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
      COMMON/INPNOF_PNOF/IPNOF,NTWOPAR
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      COMMON/SUMSZ/SUMS,SUMF
!
      DOUBLE PRECISION,DIMENSION(NV)::GAMMA       ! NV = NDOC*(NCWO+1)
      DOUBLE PRECISION,DIMENSION(NBF5)::RO
      DOUBLE PRECISION,DIMENSION(NBF5,NV)::DR
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5)::CJ12,CK12
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5,NV)::DCJ12r,DCK12r
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::BETA,SUMRO
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::DB
!-----------------------------------------------------------------------
!                 Occupations and their Derivatives
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE (BETA(NBF5),DB(NBF5,NV),SUMRO(NDOC))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RO = 0.0d0
      BETA = 0.0d0
      DR = 0.0d0
      DB = 0.0d0
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Occupations (1,NO1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(NO1>0)THEN
       DO in=1,NO1
        RO(in) = 1.0d0
        BETA(in) = 1.0d0
       ENDDO
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Occupations (NO1+1,NO1+NDOC),(NA+1,NBF5)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO i=1,NDOC
       in = NO1+i                          ! in=no1+1,nb
       RO(in)   = DEXP(GAMMA(i))
       SUMRO(i) = RO(in)
       do iw=1,ncwo
        !ig = ndoc+ncwo*(i-1)+iw            !old-sort ! ig=ndoc+1,ndoc*(ncwo+1)
        !in = na+ncwo*(ndoc-i)+iw           !old-sort in=na+1,nbf5
        ig = ndoc*(iw+1)-i+1                !new-sort
        in = no1+(na-nb)+ig                 !new-sort        
        RO(in)   = DEXP(GAMMA(ig))
        SUMRO(i) = SUMRO(i) + RO(in)
       enddo
      ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Occupations (NO1+NDOC+1,NA=NO1+NDNS)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(NSOC>0)THEN
       if(.not.HighSpin)then
        DO i=NDOC+1,NDNS
         in = NO1+i
         MULT1 = NA - NB
         RO(in)  = 0.5d0*MULT1/NSOC
         BETA(in) = DSQRT(RO(in))
        ENDDO
       else if(HighSpin)then
        DO i=NDOC+1,NDNS
         in = NO1+i
         RO(in)  = 1.0d0
         BETA(in)= 1.0d0
        ENDDO
       end if      
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Enforce pairing conditions for each subspace: SUM(RO)=1; i=1,ndoc
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO i=1,NDOC
       in = NO1+i                          ! in=no1+1,nb
       RO(in) = RO(in)/SUMRO(i) 
       do iw=1,ncwo
        !ig = ndoc+ncwo*(i-1)+iw   !old-sort         ! ig=ndoc+1,ndoc*(ncwo+1)
        !im = na+ncwo*(ndoc-i)+iw  !old-sort         ! im=na+1,nbf5
        ig = ndoc*(iw+1)-i+1       !new-sort
        im = no1+(na-nb)+ig        !new-sort
        RO(im) = RO(im)/SUMRO(i)
       enddo 
      ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Derivatives of RO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO i=1,NDOC
       in = NO1+i                          
       DR(in,i) = RO(in)*(1.0d0-RO(in))
       do iw=1,ncwo
        !ig = ndoc+ncwo*(i-1)+iw  !old-sort           
        !im = na+ncwo*(ndoc-i)+iw !old-sort
        ig = ndoc*(iw+1)-i+1      !new-sort
        im = no1+(na-nb)+ig       !new-sort
        DR(in,ig) = - RO(in)*RO(im)                        
!        
        DR(im,ig) = RO(im)*(1.0d0-RO(im))                         
        DR(im,i)  = - RO(im)*RO(in)        
        do iw1=1,iw-1
         !ig1 = ndoc+ncwo*(i-1)+iw1  !old-sort
         !im1 = na+ncwo*(ndoc-i)+iw1 !old-sort
         ig1 = ndoc*(iw1+1)-i+1      !new-sort
         im1 = no1+(na-nb)+ig1       !new-sort
         DR(im,ig1)= - RO(im)*RO(im1)         
        enddo
        do iw1=iw+1,ncwo
         !ig1 = ndoc+ncwo*(i-1)+iw1  !old-sort
         !im1 = na+ncwo*(ndoc-i)+iw1 !old-sort
         ig1 = ndoc*(iw1+1)-i+1      !new-sort
         im1 = no1+(na-nb)+ig1       !new-sort
         DR(im,ig1)= - RO(im)*RO(im1)         
        enddo
       enddo 
      ENDDO        
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Beta
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO j=1,NV
       DO i=1,NDOC
        in = NO1+i                        ! in=no1+1,nb
        BETA(in) = DSQRT(RO(in))
        DB(in,j) = 0.5d0*DR(in,j)/BETA(in)
        do iw=1,ncwo
         !in = na+ncwo*(ndoc-i)+iw         !old-sort  ! in=na+1,nbf5
         in = no1+(na-nb)+ndoc*(iw+1)-i+1  !new-sort
         BETA(in) = DSQRT(RO(in))
         DB(in,j)= 0.5d0*DR(in,j)/BETA(in)
        enddo 
       ENDDO        
      ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Sum of the Holes below the Fermi Level (SUMS)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      SUMS = DFLOAT(NB)
      do j=1,nb
       SUMS = SUMS - RO(j)
      enddo
!-----------------------------------------------------------------------
!                   CJ12, CK12, DCJ12r, DCK12r
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(IPNOF==3)THEN
       if(NSOC==0)CALL CJCKD3(NV,RO,DR,BETA,DB,CJ12,CK12,DCJ12r,DCK12r)
      ELSEIF(IPNOF==4)THEN
       if(NSOC==0)CALL CJCKD4(NV,RO,DR,BETA,DB,CJ12,CK12,DCJ12r,DCK12r)
      ELSEIF(IPNOF==5)THEN
       CALL CJCKD5(NV,RO,DR,BETA,DB,CJ12,CK12,DCJ12r,DCK12r)
      ELSEIF(IPNOF==6)THEN
       if(NSOC==0)CALL CJCKD6(NV,RO,DR,CJ12,CK12,DCJ12r,DCK12r)
      ELSEIF(IPNOF==7)THEN
       CALL CJCKD7(NV,RO,DR,BETA,DB,CJ12,CK12,DCJ12r,DCK12r)
      ELSEIF(IPNOF==8)THEN
       CALL CJCKD8(NV,RO,DR,BETA,DB,CJ12,CK12,DCJ12r,DCK12r)
      ENDIF
!-----------------------------------------------------------------------
      DEALLOCATE(BETA,DB)
      RETURN
      END
      
! ENENEWr
      SUBROUTINE ENENEWr(RO,HCORE,QJ,QK,CJ12,CK12,DIPx,DIPy,DIPz,EAHF,E)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      LOGICAL EFIELDL,HighSpin
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      COMMON/INP_EFIELDL/EX,EY,EZ,EFIELDL
      COMMON/EHFEN/EHF,EN
!
      DOUBLE PRECISION,DIMENSION(NBF)::EAHF,E
      DOUBLE PRECISION,DIMENSION(NBF5)::RO,HCORE,DIPx,DIPy,DIPz
      DOUBLE PRECISION,DIMENSION(NBFT5)::QJ,QK
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5)::CJ12,CK12
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::AUX
!-----------------------------------------------------------------------
!     EAHF: New HF energies
!-----------------------------------------------------------------------
      ALLOCATE (AUX(NBF5,NBF5))
      DO j=1,NBF5
       DO i=1,NBF5
        AUX(j,i) = RO(j)*RO(i)
       ENDDO
      ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO j=1,na
       FI2j= RO(j)*(1.0-RO(j))
       FIj = SQRT(FI2j)
        jj = j*(j+1)/2
        ja = NBF5-j+1
       jja = j + ja*(ja-1)/2
       EAHF(j) =  RO(j)*HCORE(j) + 0.5d0*FI2j*QJ(jj)                    &
               - 0.5d0*FIj*QK(jja) - FI2j*(QJ(jja)-0.5d0*QK(jja))        
       do i=1,j                                                          
        ij = i + j*(j-1)/2                                               
        EAHF(j) = EAHF(j) + AUX(j,i)*(2.0d0*QJ(ij)-QK(ij))               
       enddo                                                             
       do i=j+1,NBF5                                                     
        ij = j + i*(i-1)/2                                               
        EAHF(j) = EAHF(j) + AUX(j,i)*(2.0d0*QJ(ij)-QK(ij))               
       enddo                                                             
      ENDDO                                                              
!                                                                        
      DO j=na+1,nbf5                                                     
       FI2j= RO(j)*(1.0-RO(j))                                           
       FIj = SQRT(FI2j)                                                  
        jj = j*(j+1)/2                                                   
        ja = NBF5-j+1                                                    
       jja = ja + j*(j-1)/2                                              
       EAHF(j) =  RO(j)*HCORE(j) + 0.5d0*FI2j*QJ(jj)                    &
               - 0.5d0*FIj*QK(jja) - FI2j*(QJ(jja)-0.5d0*QK(jja))
       do i=1,j
        ij = i + j*(j-1)/2
        EAHF(j) = EAHF(j) + AUX(j,i)*(2.0d0*QJ(ij)-QK(ij))
       enddo
       do i=j+1,NBF5
        ij = j + i*(i-1)/2
        EAHF(j) = EAHF(j) + AUX(j,i)*(2.0d0*QJ(ij)-QK(ij))
       enddo
      ENDDO
!      
      IF(NBF5<NBF)THEN
       DO J=NBF5+1,NBF
        EAHF(J) = 0.0d0
       ENDDO
      ENDIF
!-----------------------------------------------------------------------
!     E: Diagonal of the Lagrangian (ELAG)
!-----------------------------------------------------------------------
      if(MSpin==0)then
       DO J=1,NB
        JJ=J*(J+1)/2
        E(J) = RO(J) * ( HCORE(J) + QJ(JJ) )                            &
             + PRODCWQWj(J,CJ12,QJ) - PRODCWQWj(J,CK12,QK)               
       ENDDO                                                             
       DO J=NB+1,NA                                                      
        E(J) = RO(J) * HCORE(J)                                         &
             + PRODCWQWj(J,CJ12,QJ)-PRODCWQWj(J,CK12,QK)                 
       ENDDO                                                             
       DO J=NA+1,NBF5                                                    
        JJ=J*(J+1)/2                                                     
        E(J) = RO(J) * ( HCORE(J) + QJ(JJ) )                            &
             + PRODCWQWj(J,CJ12,QJ) - PRODCWQWj(J,CK12,QK)               
       ENDDO                                                             
      else if(MSpin>0)then                                               
       DO J=1,NB                                                         
        JJ=J*(J+1)/2                                                     
        E(J) = RO(J) * ( HCORE(J) + QJ(JJ) )                            &
             + PRODCWQWj1(J,CJ12,QJ) - PRODCWQWj1(J,CK12,QK)            &
             + 2.0d0*PRODROQWj1(J,RO,QJ)-PRODROQWj1(J,RO,QK)             
       ENDDO                                                             
       DO J=NB+1,NA                                                      
        E(J) = 0.5d0*( RO(J)*HCORE(J) + PRODROQWj0(J,RO,QJ)             &
                      - PRODROQWj0(J,RO,QK) )                            
       ENDDO                                                             
       DO J=NA+1,NBF5                                                    
        JJ=J*(J+1)/2                                                     
        E(J) = RO(J) * ( HCORE(J) + QJ(JJ) )                            &
             + PRODCWQWj2(J,CJ12,QJ) - PRODCWQWj2(J,CK12,QK)            &
             + 2.0d0*PRODROQWj2(J,RO,QJ)-PRODROQWj2(J,RO,QK)
       ENDDO
      end if 
!      
      IF(NBF5<NBF)THEN
       DO J=NBF5+1,NBF
        E(J) = 0.0d0                  ! /= ELAG(J,J) !!!
       ENDDO
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Including Electric Field
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(EFIELDL)THEN
       DO J=1,NBF5
        E(J) = E(J) + (EX*DIPx(J)+EY*DIPy(J)+EZ*DIPz(J))*RO(J)
       ENDDO
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Calculation of the Total Energy 
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      ETOTAL = 0.0
!      DO J=1,NB
!       ETOTAL = ETOTAL + E(J) + RO(J)*HCORE(J)
!      ENDDO
!      if(MSpin==0)then
!       DO J=NB+1,NA
!        ETOTAL = ETOTAL + E(J) + RO(J)*HCORE(J)       
!       ENDDO
!      else if(MSpin>0)then       
!       DO J=NB+1,NA
!        ETOTAL = ETOTAL + E(J) + 0.5d0*RO(J)*HCORE(J)       
!       ENDDO
!      end if
!      DO J=NA+1,NBF5
!       ETOTAL = ETOTAL + E(J) + RO(J)*HCORE(J)
!      ENDDO
!      IF(EFIELDL)THEN
!       DO J=1,NBF5
!        ETOTAL = ETOTAL + (EX*DIPx(J)+EY*DIPy(J)+EZ*DIPz(J))*RO(J)
!       ENDDO
!      ENDIF
!      write(6,*)'ETOTAL=',ETOTAL + EN
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DEALLOCATE (AUX)
      RETURN
      END

! CGOCUPSUMSL
      SUBROUTINE CGOCUPSUMSL(NV,GAMMA,USER,ENERGY)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      COMMON/ENERGIAS/EELEC,EELEC_OLD,DIF_EELEC,EELEC_MIN
      COMMON/PUNTEROSUSER/N1,N2,N3,N4,N5,N6,N7,N8,N9,N10,N11,N12,N13,   &
                          N14,N15,N16,N17,N18,N19,N20,N21,N22,N23,N24,  &
                          N25,N26,N27,N28,N29,N30,N31,N32,N33,N34,N35,  &
                          N36,N37,N38,N39,N40,N41,N42,N43,N44,N45,N46,  &
                          N47,N48,N49,N50,N51,NUSER
!
      DOUBLE PRECISION,DIMENSION(NV) :: GAMMA
      DOUBLE PRECISION,DIMENSION(NUSER) :: USER
      INTEGER,ALLOCATABLE,DIMENSION(:) :: IUSER,IV      
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: D,V     
      EXTERNAL CALCOE,CALCOG
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      LIV = 60
      LV = 71+NV*(NV+15)/2
      ALLOCATE( IUSER(1),IV(LIV),D(NV),V(LV) ) 
      IUSER(1) = 1
      IV = 0
      D(1:NV) = 1.0d-1       
      CALL SUMSL(NV,D,GAMMA,CALCOE,CALCOG,IV,LIV,LV,V,IUSER,USER)
      ENERGY = EELEC
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DEALLOCATE(IUSER,IV,D,V)
      RETURN
      END

! CALCOE      
      SUBROUTINE CALCOE(NV,GAMMA,NF,ENERGY,IUSER,USER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      LOGICAL EFIELDL,HighSpin
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_ORBSPACE0/NO1,NDOC,NCO,NCWO,NVIR,NAC,NO0
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
      COMMON/INPNOF_GENERALINF/ICOEF,ISOFTMAX,IORBOPT,MAXIT,MAXIT21      
      COMMON/INP_EFIELDL/EX,EY,EZ,EFIELDL
      COMMON/ENERGIAS/EELEC,EELEC_OLD,DIF_EELEC,EELEC_MIN
      COMMON/PUNTEROSUSER/N1,N2,N3,N4,N5,N6,N7,N8,N9,N10,N11,N12,N13,   &
                          N14,N15,N16,N17,N18,N19,N20,N21,N22,N23,N24,  &
                          N25,N26,N27,N28,N29,N30,N31,N32,N33,N34,N35,  &
                          N36,N37,N38,N39,N40,N41,N42,N43,N44,N45,N46,  &
                          N47,N48,N49,N50,N51,NUSER
!
      INTEGER,DIMENSION(1) :: IUSER 
      DOUBLE PRECISION,DIMENSION(NV)   :: GAMMA
      DOUBLE PRECISION,DIMENSION(NUSER):: USER
!-----------------------------------------------------------------------
!     Avoiding warnings
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!BORRAR NF = NF (Esta causando Segmentation Fault)
      IUSER(1) = IUSER(1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Occupations
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(ISOFTMAX==0)THEN
       CALL OCUPACIONTFr(GAMMA,USER(N1),USER(N2),USER(N3),              &
                               USER(N4),USER(N5),USER(N6),NV)
      ELSE IF(ISOFTMAX==1)THEN
       CALL OCUPACIONSMr(GAMMA,USER(N1),USER(N2),USER(N3),              &
                               USER(N4),USER(N5),USER(N6),NV)
      ENDIF
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!     Singlet State (S=0,Ms=0) and Multiplet States (S>0,Ms=0)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
      if(MSpin==0)then
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       IF(NCWO==1)THEN       ! PNOFi(1): Perfect Pairing (NCWO=1)
        ENERGY = 0.0d0
        do in=1,NO1
         ENERGY = ENERGY + USER(N1-1+in)                                &
                * ( 2.0*USER(N8-1+in) + USER(N9-1+in*(in+1)/2) )        &
                + PRODCWQWj(in,USER(N2),USER(N9))                       &
                - PRODCWQWj(in,USER(N3),USER(N10))                       
        enddo                                                            
        do i=1,NDOC                                                      
         in = NO1+i                                                      
         ENERGY = ENERGY + USER(N1-1+in)                                &
                * ( 2.0*USER(N8-1+in) + USER(N9-1+in*(in+1)/2) )        &
                + PRODCWQWj(in,USER(N2),USER(N9))                       &
                - PRODCWQWj(in,USER(N3),USER(N10))                       
         in = na+ndoc-i+1                                                
         ENERGY = ENERGY + USER(N1-1+in)                                &
                * ( 2.0*USER(N8-1+in) + USER(N9-1+in*(in+1)/2) )        &
                + PRODCWQWj(in,USER(N2),USER(N9))                       &
                - PRODCWQWj(in,USER(N3),USER(N10))                       
        enddo                                                            
       END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       IF(NCWO>1)THEN        ! PNOFi(Nc): Extended PNOF (NCWO>1)
        ENERGY = 0.0d0
        do in=1,NO1
         ENERGY = ENERGY + USER(N1-1+in)                                &
                * ( 2.0*USER(N8-1+in) + USER(N9-1+in*(in+1)/2) )        &
                + PRODCWQWj(in,USER(N2),USER(N9))                       &
                - PRODCWQWj(in,USER(N3),USER(N10))                       
        enddo                                                            
        do i=1,NDOC                                                      
         in = NO1+i                                                      
         ENERGY = ENERGY + USER(N1-1+in)                                &
                * ( 2.0*USER(N8-1+in) + USER(N9-1+in*(in+1)/2) )        &
                + PRODCWQWj(in,USER(N2),USER(N9))                       &
                - PRODCWQWj(in,USER(N3),USER(N10))                       
         do iw=1,ncwo-1                                                  
          !in = na+ncwo*(ndoc-i)+iw              !old-sort                                   
          in = no1+(na-nb)+ndoc*(iw+1)-i+1       !new-sort
          ENERGY = ENERGY + USER(N1-1+in)                               &
                 * ( 2.0*USER(N8-1+in) + USER(N9-1+in*(in+1)/2) )       &
                 + PRODCWQWj(in,USER(N2),USER(N9))                      &
                 - PRODCWQWj(in,USER(N3),USER(N10))                      
         enddo                                                           
         !in = na+ncwo*(ndoc-i)+ncwo              !old-sort                        
         in = no1+(na-nb)+ndoc*(ncwo+1)-i+1       !new-sort
         ENERGY = ENERGY + USER(N1-1+in)                                &
                * ( 2.0*USER(N8-1+in) + USER(N9-1+in*(in+1)/2) )        &
                + PRODCWQWj(in,USER(N2),USER(N9))                       &
                - PRODCWQWj(in,USER(N3),USER(N10))                       
        enddo                                                            
       END IF
       IF(NSOC>0)THEN                                                   
        do i=NDOC+1,NDNS                                                
         in = NO1+i                                                     
         ENERGY = ENERGY + 2.0d0*USER(N1-1+in)*USER(N8-1+in)           &
                + PRODCWQWj(in,USER(N2),USER(N9))                      &
                - PRODCWQWj(in,USER(N3),USER(N10))
        enddo
       END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       IF(EFIELDL)THEN       ! including electric field
        CALL DIPMOMr(USER(N11),USER(N12),USER(N13),USER(N14),          &
                     USER(N15),USER(N16),USER(N17),USER(N7),           &
                     USER(N1),DMXe,DMYe,DMZe,DMX,DMY,DMZ,DM)
        ENERGY = ENERGY - EX*DMX - EY*DMY - EZ*DMZ
       END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      end if
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --           
!     High-Spin Multiplet State (S>0,Ms=S)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
      if(MSpin>0)then
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       IF(NCWO==1)THEN       ! PNOFi(1): Perfect Pairing (NCWO=1)
        ENERGY = 0.0d0
        do in=1,NO1
         ENERGY = ENERGY + USER(N1-1+in)                                &
                * ( 2.0*USER(N8-1+in) + USER(N9-1+in*(in+1)/2) )        &
                + PRODCWQWj1(in,USER(N2),USER(N9))                      &
                - PRODCWQWj1(in,USER(N3),USER(N10))                     &
                + 2.0d0*PRODROQWj1(in,USER(N1),USER(N9))                &
                - PRODROQWj1(in,USER(N1),USER(N10))                      
        enddo                                                            
        do i=1,NDOC                                                      
         in = NO1+i                                                      
         ENERGY = ENERGY + USER(N1-1+in)                                &
                * ( 2.0*USER(N8-1+in) + USER(N9-1+in*(in+1)/2) )        &
                + PRODCWQWj1(in,USER(N2),USER(N9))                      &
                - PRODCWQWj1(in,USER(N3),USER(N10))                     &
                + 2.0d0*PRODROQWj1(in,USER(N1),USER(N9))                &
                - PRODROQWj1(in,USER(N1),USER(N10))                      
         in = na+ndoc-i+1                                                
         ENERGY = ENERGY + USER(N1-1+in)                                &
                * ( 2.0*USER(N8-1+in) + USER(N9-1+in*(in+1)/2) )        &
                + PRODCWQWj2(in,USER(N2),USER(N9))                      &
                - PRODCWQWj2(in,USER(N3),USER(N10))                     &
                + 2.0d0*PRODROQWj2(in,USER(N1),USER(N9))                &
                - PRODROQWj2(in,USER(N1),USER(N10))                     
        enddo                                                           
        do i=NDOC+1,NDNS                                                
         in = NO1+i                                                     
         ENERGY = ENERGY + USER(N1-1+in)*USER(N8-1+in)                  &
                + 0.5d0*(PRODROQWj0(in,USER(N1),USER(N9))               &
                - PRODROQWj0(in,USER(N1),USER(N10)))
        enddo
       END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       IF(NCWO>1)THEN        ! PNOFi(Nc): Extended PNOF (NCWO>1)
        ENERGY = 0.0d0
        do in=1,NO1
         ENERGY = ENERGY + USER(N1-1+in)                               &
                * ( 2.0*USER(N8-1+in) + USER(N9-1+in*(in+1)/2) )       &
                + PRODCWQWj1(in,USER(N2),USER(N9))                     &
                - PRODCWQWj1(in,USER(N3),USER(N10))                    &
                + 2.0d0*PRODROQWj1(in,USER(N1),USER(N9))               &
                - PRODROQWj1(in,USER(N1),USER(N10))                     
        enddo                                                           
        do i=1,NDOC                                                     
         in = NO1+i                                                     
         ENERGY = ENERGY + USER(N1-1+in)                               &
                * ( 2.0*USER(N8-1+in) + USER(N9-1+in*(in+1)/2) )       &
                + PRODCWQWj1(in,USER(N2),USER(N9))                     &
                - PRODCWQWj1(in,USER(N3),USER(N10))                    &
                + 2.0d0*PRODROQWj1(in,USER(N1),USER(N9))               &
                - PRODROQWj1(in,USER(N1),USER(N10))                     
         do iw=1,ncwo-1                                                 
          !in = na+ncwo*(ndoc-i)+iw              !old-sort                       
          in = no1+(na-nb)+ndoc*(iw+1)-i+1       !new-sort
          ENERGY = ENERGY + USER(N1-1+in)                              &
                 * ( 2.0*USER(N8-1+in) + USER(N9-1+in*(in+1)/2) )      &
                 + PRODCWQWj2(in,USER(N2),USER(N9))                    &
                 - PRODCWQWj2(in,USER(N3),USER(N10))                   &
                 + 2.0d0*PRODROQWj2(in,USER(N1),USER(N9))              &
                 - PRODROQWj2(in,USER(N1),USER(N10))                    
         enddo                                                          
         !in = na+ncwo*(ndoc-i)+ncwo              !old-sort                                    
         in = no1+(na-nb)+ndoc*(ncwo+1)-i+1       !new-sort
         ENERGY = ENERGY + USER(N1-1+in)                               &
                * ( 2.0*USER(N8-1+in) + USER(N9-1+in*(in+1)/2) )       &
                + PRODCWQWj2(in,USER(N2),USER(N9))                     &
                - PRODCWQWj2(in,USER(N3),USER(N10))                    &
                + 2.0d0*PRODROQWj2(in,USER(N1),USER(N9))               &
                - PRODROQWj2(in,USER(N1),USER(N10))          
        enddo
        do i=NDOC+1,NDNS
         in = NO1+i
         ENERGY = ENERGY + USER(N1-1+in)*USER(N8-1+in)                 &
                + 0.5d0*(PRODROQWj0(in,USER(N1),USER(N9))              &
                - PRODROQWj0(in,USER(N1),USER(N10)))
        enddo
       ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       IF(EFIELDL)THEN       ! including electric field
        CALL DIPMOMr(USER(N11),USER(N12),USER(N13),USER(N14),          &
                     USER(N15),USER(N16),USER(N17),USER(N7),           &
                     USER(N1),DMXe,DMYe,DMZe,DMX,DMY,DMZ,DM)
        ENERGY = ENERGY - EX*DMX - EY*DMY - EZ*DMZ
       END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      end if
!-----------------------------------------------------------------------
      EELEC = ENERGY
      RETURN
      END

! CALCOG      
      SUBROUTINE CALCOG(NV,GAMMA,NF,GRAD,IUSER,USER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      LOGICAL EFIELDL,HighSpin
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_ORBSPACE0/NO1,NDOC,NCO,NCWO,NVIR,NAC,NO0
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
      COMMON/INPNOF_GENERALINF/ICOEF,ISOFTMAX,IORBOPT,MAXIT,MAXIT21            
      COMMON/INP_EFIELDL/EX,EY,EZ,EFIELDL
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      COMMON/PUNTEROSUSER/N1,N2,N3,N4,N5,N6,N7,N8,N9,N10,N11,N12,N13,   &
                          N14,N15,N16,N17,N18,N19,N20,N21,N22,N23,N24,  &
                          N25,N26,N27,N28,N29,N30,N31,N32,N33,N34,N35,  &
                          N36,N37,N38,N39,N40,N41,N42,N43,N44,N45,N46,  &
                          N47,N48,N49,N50,N51,NUSER
!
      INTEGER,DIMENSION(1) :: IUSER 
      DOUBLE PRECISION,DIMENSION(NV) :: GAMMA,GRAD
      DOUBLE PRECISION,DIMENSION(NUSER) :: USER
!-----------------------------------------------------------------------
!     Avoiding warnings
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      NF = NF
      IUSER(1) = IUSER(1)
      GAMMA(1) = GAMMA(1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Occupations
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     IF(ISOFTMAX==0)THEN
!      CALL OCUPACIONTFr(GAMMA,USER(N1),USER(N2),USER(N3),              &
!                              USER(N4),USER(N5),USER(N6),NV)
!     ELSE IF(ISOFTMAX==1)THEN
!      CALL OCUPACIONSMr(GAMMA,USER(N1),USER(N2),USER(N3),              &
!                              USER(N4),USER(N5),USER(N6),NV)
!     ENDIF
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!     Singlet State (S=0,Ms=0) and Multiplet States (S>0,Ms=0)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
      if(MSpin==0)then
      
       IF(NCWO==1)THEN       ! PNOFi(1): Perfect Pairing (NCWO=1)
        GRAD = 0.0d0
        DO ig=1,NV
         do i=1,NDOC
          in = NO1+i
          GRAD(ig) = GRAD(ig) + USER(N4-1+in+(ig-1)*nbf5)               &
                   * ( 2.0d0*USER(N8-1+in) + USER(N9-1+in*(in+1)/2) )   &
                   + 2.0d0 * ( PRODCWQWjk(nv,in,ig,USER(N5),USER(N9))   &
                             - PRODCWQWjk(nv,in,ig,USER(N6),USER(N10)) ) 
          in = na+ndoc-i+1                                               
          GRAD(ig) = GRAD(ig) + USER(N4-1+in+(ig-1)*nbf5)               &
                   * ( 2.0d0*USER(N8-1+in) + USER(N9-1+in*(in+1)/2) )   &
                   + 2.0d0 * ( PRODCWQWjk(nv,in,ig,USER(N5),USER(N9))   &
                             - PRODCWQWjk(nv,in,ig,USER(N6),USER(N10)) )
         enddo
        ENDDO

        if(EFIELDL)then      ! including electric field
         DO ig=1,NV
          do i=1,NDOC
           in = NO1+i
           GRAD(ig) = GRAD(ig) + 2.0d0*USER(N4-1+in+(ig-1)*nbf5)       &
                    * ( EX*USER(N15-1+in) + EY*USER(N16-1+in)          &
                      + EZ*USER(N17-1+in) )                             
           in = na+ndoc-i+1                                             
           GRAD(ig) = GRAD(ig) + 2.0d0*USER(N4-1+in+(ig-1)*nbf5)       &
                    * ( EX*USER(N15-1+in) + EY*USER(N16-1+in)          &
                      + EZ*USER(N17-1+in) )
          enddo
         ENDDO
        endif
       END IF

       IF(NCWO>1)THEN        ! PNOFi(Nc): Extended PNOF (NCWO>1)

        GRAD = 0.0d0
        DO ig=1,NV
         do i=1,NDOC
          in = NO1+i
          GRAD(ig) = GRAD(ig) + USER(N4-1+in+(ig-1)*nbf5)               &
                   * ( 2.0d0*USER(N8-1+in) + USER(N9-1+in*(in+1)/2) )   &
                   + 2.0d0 * ( PRODCWQWjk(nv,in,ig,USER(N5),USER(N9))   &
                             - PRODCWQWjk(nv,in,ig,USER(N6),USER(N10)) ) 
          do iw=1,NCWO-1                                                 
           !in = na+ncwo*(ndoc-i)+iw              !old-sort                        
           in = no1+(na-nb)+ndoc*(iw+1)-i+1       !new-sort
           GRAD(ig) = GRAD(ig) + USER(N4-1+in+(ig-1)*nbf5)              &
                    * ( 2.0d0*USER(N8-1+in) + USER(N9-1+in*(in+1)/2) )  &
                    + 2.0d0 * (PRODCWQWjk(nv,in,ig,USER(N5),USER(N9))   &
                              -PRODCWQWjk(nv,in,ig,USER(N6),USER(N10)) ) 
          enddo                                                          
          !in = na+ncwo*(ndoc-i)+ncwo              !old-sort                       
          in = no1+(na-nb)+ndoc*(ncwo+1)-i+1       !new-sort
          GRAD(ig) = GRAD(ig) + USER(N4-1+in+(ig-1)*nbf5)               &
                   * ( 2.0d0*USER(N8-1+in) + USER(N9-1+in*(in+1)/2) )   &
                   + 2.0d0 * ( PRODCWQWjk(nv,in,ig,USER(N5),USER(N9))   &
                             - PRODCWQWjk(nv,in,ig,USER(N6),USER(N10)) )
         enddo
        ENDDO

        if(EFIELDL)then      ! including electric field
         DO ig=1,NV
          do i=1,NDOC
           in = NO1+i
           GRAD(ig) = GRAD(ig) + 2.0d0*USER(N4-1+in+(ig-1)*nbf5)       &
                    * ( EX*USER(N15-1+in) + EY*USER(N16-1+in)          &
                      + EZ*USER(N17-1+in) )                             
           do iw=1,NCWO-1                                               
            !in = na+ncwo*(ndoc-i)+iw              !old-sort                    
            in = no1+(na-nb)+ndoc*(iw+1)-i+1       !new-sort
            GRAD(ig) = GRAD(ig) + 2.0d0*USER(N4-1+in+(ig-1)*nbf5)      &
                     * ( EX*USER(N15-1+in) + EY*USER(N16-1+in)         &
                       + EZ*USER(N17-1+in) )                            
           enddo                                                        
           !in = na+ncwo*(ndoc-i)+ncwo              !old-sort                     
           in = no1+(na-nb)+ndoc*(ncwo+1)-i+1       !new-sort
           GRAD(ig) = GRAD(ig) + 2.0d0*USER(N4-1+in+(ig-1)*nbf5)       &
                    * ( EX*USER(N15-1+in) + EY*USER(N16-1+in)          &
                      + EZ*USER(N17-1+in) )
          enddo
         ENDDO
        endif

       END IF

      end if
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --           
!     High-Spin Multiplet State (S>0,Ms=S)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
      if(MSpin>0)then
      
       IF(NCWO==1)THEN       ! PNOFi(1): Perfect Pairing (NCWO=1)

        GRAD = 0.0d0
        DO ig=1,NV
         do i=1,NDOC
          in = NO1+i
          GRAD(ig) = GRAD(ig) + USER(N4-1+in+(ig-1)*nbf5)               &
                   * ( 2.0d0*USER(N8-1+in) + USER(N9-1+in*(in+1)/2) )   &
                   + 2.0d0 * ( PRODCWQWjk1(nv,in,ig,USER(N5),USER(N9))  &
                             - PRODCWQWjk1(nv,in,ig,USER(N6),USER(N10)))&
                   + 2.0d0 *   PRODDRQWjk1(nv,in,ig,USER(N4),USER(N9))  &
                   -           PRODDRQWjk1(nv,in,ig,USER(N4),USER(N10))  
          in = na+ndoc-i+1                                               
          GRAD(ig) = GRAD(ig) + USER(N4-1+in+(ig-1)*nbf5)               &
                   * ( 2.0d0*USER(N8-1+in) + USER(N9-1+in*(in+1)/2) )   &
                   + 2.0d0 * ( PRODCWQWjk2(nv,in,ig,USER(N5),USER(N9))  &
                             - PRODCWQWjk2(nv,in,ig,USER(N6),USER(N10)))&
                   + 2.0d0 *   PRODDRQWjk2(nv,in,ig,USER(N4),USER(N9))  &
                   -           PRODDRQWjk2(nv,in,ig,USER(N4),USER(N10))
         enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        ENDDO

        if(EFIELDL)then      ! including electric field   
         DO ig=1,NV
          do i=1,NDOC
           in = NO1+i
           GRAD(ig) = GRAD(ig) + 2.0d0*USER(N4-1+in+(ig-1)*nbf5)       &
                    * ( EX*USER(N15-1+in) + EY*USER(N16-1+in)          &
                      + EZ*USER(N17-1+in) )                             
           in = na+ndoc-i+1                                             
           GRAD(ig) = GRAD(ig) + 2.0d0*USER(N4-1+in+(ig-1)*nbf5)       &
                    * ( EX*USER(N15-1+in) + EY*USER(N16-1+in)          &
                      + EZ*USER(N17-1+in) )
          enddo
         ENDDO
        endif

       END IF
       
       IF(NCWO>1)THEN        ! PNOFi(Nc): Extended PNOF (NCWO>1)

        GRAD = 0.0d0
        DO ig=1,NV
         do i=1,NDOC
          in = NO1+i
          GRAD(ig) = GRAD(ig) + USER(N4-1+in+(ig-1)*nbf5)               &
                   * ( 2.0d0*USER(N8-1+in) + USER(N9-1+in*(in+1)/2) )   &
                   + 2.0d0 * ( PRODCWQWjk1(nv,in,ig,USER(N5),USER(N9))  &
                             - PRODCWQWjk1(nv,in,ig,USER(N6),USER(N10)))&
                   + 2.0d0 *   PRODDRQWjk1(nv,in,ig,USER(N4),USER(N9))  &
                   -           PRODDRQWjk1(nv,in,ig,USER(N4),USER(N10))  
          do iw=1,NCWO-1                                                 
           !in = na+ncwo*(ndoc-i)+iw              !old-sort                      
           in = no1+(na-nb)+ndoc*(iw+1)-i+1       !new-sort
           GRAD(ig) = GRAD(ig) + USER(N4-1+in+(ig-1)*nbf5)              &
                    * ( 2.0d0*USER(N8-1+in) + USER(N9-1+in*(in+1)/2) )  &
                    + 2.0d0 * (PRODCWQWjk2(nv,in,ig,USER(N5),USER(N9))  &
                              -PRODCWQWjk2(nv,in,ig,USER(N6),USER(N10)))&
                    + 2.0d0 *  PRODDRQWjk2(nv,in,ig,USER(N4),USER(N9))  &
                    -          PRODDRQWjk2(nv,in,ig,USER(N4),USER(N10))  
          enddo                                                          
          !in = na+ncwo*(ndoc-i)+ncwo              !old-sort                       
          in = no1+(na-nb)+ndoc*(ncwo+1)-i+1       !new-sort
          GRAD(ig) = GRAD(ig) + USER(N4-1+in+(ig-1)*nbf5)               &
                   * ( 2.0d0*USER(N8-1+in) + USER(N9-1+in*(in+1)/2) )   &
                   + 2.0d0 * ( PRODCWQWjk2(nv,in,ig,USER(N5),USER(N9))  &
                             - PRODCWQWjk2(nv,in,ig,USER(N6),USER(N10)))&
                   + 2.0d0 *   PRODDRQWjk2(nv,in,ig,USER(N4),USER(N9))  &
                   -           PRODDRQWjk2(nv,in,ig,USER(N4),USER(N10))
         enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -          
        ENDDO

        if(EFIELDL)then      ! including electric field   
         DO ig=1,NV
          do i=1,NDOC
           in = NO1+i
           GRAD(ig) = GRAD(ig) + 2.0d0*USER(N4-1+in+(ig-1)*nbf5)       &
                    * ( EX*USER(N15-1+in) + EY*USER(N16-1+in)          &
                      + EZ*USER(N17-1+in) )                             
           do iw=1,NCWO-1                                               
            !in = na+ncwo*(ndoc-i)+iw              !old-sort                                   
            in = no1+(na-nb)+ndoc*(iw+1)-i+1       !new-sort
            GRAD(ig) = GRAD(ig) + 2.0d0*USER(N4-1+in+(ig-1)*nbf5)      &
                     * ( EX*USER(N15-1+in) + EY*USER(N16-1+in)         &
                       + EZ*USER(N17-1+in) )                            
           enddo                                                        
           !in = na+ncwo*(ndoc-i)+ncwo              !old-sort                                   
           in = no1+(na-nb)+ndoc*(ncwo+1)-i+1       !new-sort
           GRAD(ig) = GRAD(ig) + 2.0d0*USER(N4-1+in+(ig-1)*nbf5)       &
                    * ( EX*USER(N15-1+in) + EY*USER(N16-1+in)          &
                      + EZ*USER(N17-1+in) )
          enddo
         ENDDO
        endif

       ENDIF

      end if
!-----------------------------------------------------------------------
      RETURN
      END
      
! CGOCUPNAG
      SUBROUTINE CGOCUPNAG(NV,GAMMA,USER,ENERGY)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_THRESH/THRESHL,THRESHE,THRESHEC,THRESHEN
      COMMON/PUNTEROSUSER/N1,N2,N3,N4,N5,N6,N7,N8,N9,N10,N11,N12,N13,   &
                          N14,N15,N16,N17,N18,N19,N20,N21,N22,N23,N24,  &
                          N25,N26,N27,N28,N29,N30,N31,N32,N33,N34,N35,  &
                          N36,N37,N38,N39,N40,N41,N42,N43,N44,N45,N46,  &
                          N47,N48,N49,N50,N51,NUSER
!
      DOUBLE PRECISION,DIMENSION(NV)::GAMMA
      DOUBLE PRECISION,DIMENSION(NUSER)::USER
      INTEGER,ALLOCATABLE,DIMENSION(:)::IUSER,IWORK
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::GRAD,WORK
      EXTERNAL ENERFUNr
      ALLOCATE (IUSER(1),IWORK(NV+1),GRAD(NV),WORK(13*NV))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Send output of E04DGF to CGM file
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      CALL X04ABF(1,2)                                        !nag
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Function Precision (machine precision**0.9)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(THRESHEN>=1.0d-10)THEN
!       CALL E04DKF ('Function Precision = 1.0D-10')           !nag
      ELSEIF(THRESHEN==1.0d-11)THEN                               
!       CALL E04DKF ('Function Precision = 1.0D-11')           !nag
      ELSEIF(THRESHEN==1.0d-12)THEN                                   
!       CALL E04DKF ('Function Precision = 1.0D-12')           !nag
      ELSEIF(THRESHEN==1.0d-13)THEN                                   
!       CALL E04DKF ('Function Precision = 1.0D-13')           !nag
      ELSEIF(THRESHEN==1.0d-14)THEN                                   
!       CALL E04DKF ('Function Precision = 1.0D-14')           !nag
      ELSEIF(THRESHEN==1.0d-15)THEN                                   
!       CALL E04DKF ('Function Precision = 1.0D-15')           !nag
      ELSEIF(THRESHEN<=1.0d-16)THEN                                   
!       CALL E04DKF ('Function Precision = 1.0D-16')           !nag
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Verify Level (-1 = No checks, 0 = cheap test, 1 = 0 + gradients)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      CALL E04DKF ('Verify Level = -1')                       !nag
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Print Level (0 = No output, 1 = The final solution only)
!                 (5 = One line of summary output for each iteration)
!                 (10 = The final solution and one line for each iter.)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      CALL E04DKF ('Print Level = 0')                         !nag
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Maximum Step Length (Default = 10**20))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      CALL E04DKF ('Maximum Step Length  = 0.1')              !nag
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Linesearch Tolerance (0<r<1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      CALL E04DKF ('Linesearch Tolerance = 0.01')             !nag
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Optimality Tolerance (default = relative precision**0.8)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      CALL E04DKF ('Optimality Tolerance = 1.0D-10')          !nag
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Calling to NAG Library for using the CG method
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IFAIL = 1
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -      
!     Avoiding warnings
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      GAMMA(1) = GAMMA(1)
!      CALL E04DGF(NV,ENERFUNr,ITER_E04DGF,ENERGY,GRAD,GAMMA,            &    !nag
!                  IWORK,WORK,IUSER,USER,IFAIL)                               !nag
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DEALLOCATE (IUSER,IWORK,GRAD,WORK)
      RETURN
      END

! ENERFUNr
      SUBROUTINE ENERFUNr(MODE,NV,X,ENERGY,GRAD,NSTATE,IUSER,USER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      LOGICAL HighSpin
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
      COMMON/PUNTEROSUSER/N1,N2,N3,N4,N5,N6,N7,N8,N9,N10,N11,N12,N13,   &
                          N14,N15,N16,N17,N18,N19,N20,N21,N22,N23,N24,  &
                          N25,N26,N27,N28,N29,N30,N31,N32,N33,N34,N35,  &
                          N36,N37,N38,N39,N40,N41,N42,N43,N44,N45,N46,  &
                          N47,N48,N49,N50,N51,NUSER
!
      INTEGER,DIMENSION(*)::IUSER
      DOUBLE PRECISION,DIMENSION(NV)::X,GRAD
      DOUBLE PRECISION,DIMENSION(NUSER)::USER
!-----------------------------------------------------------------------
      NSTATE=NSTATE
      IUSER(1)=IUSER(1)
!-----------------------------------------------------------------------
      if(MSpin==0)then
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -
!      Singlet State (S=0,Ms=0) and Multiplet States (S>0,Ms=0)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -
       CALL OCUPENERGYrc(MODE,X,USER(N1),USER(N2),USER(N3),USER(N4),    &
                         USER(N5),USER(N6),USER(N7),USER(N8),USER(N9),  &
                         USER(N10),USER(N11),USER(N12),USER(N13),       &
                         USER(N14),USER(N15),USER(N16),USER(N17),       &
                         ENERGY,GRAD,NV)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -   
      else if(MSpin>0)then
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -
!      High-Spin Multiplet State (S>0,Ms=S)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -
       CALL OCUPENERGYro(MODE,X,USER(N1),USER(N2),USER(N3),USER(N4),    &
                         USER(N5),USER(N6),USER(N7),USER(N8),USER(N9),  &
                         USER(N10),USER(N11),USER(N12),USER(N13),       &
                         USER(N14),USER(N15),USER(N16),USER(N17),       &
                         ENERGY,GRAD,NV)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -
      end if
!-----------------------------------------------------------------------
      RETURN
      END

! OCUPENERGYrc
      SUBROUTINE OCUPENERGYrc(MODE,GAMMA,RO,CJ12,CK12,DR,DCJ12r,DCK12r, &
                              QD,HCORE,QJ,QK,DIPN,ADIPx,ADIPy,ADIPz,    &
                              DIPx,DIPy,DIPz,ENERGY,GRAD,NV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      LOGICAL EFIELDL,HighSpin
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_ORBSPACE0/NO1,NDOC,NCO,NCWO,NVIR,NAC,NO0
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin 
      COMMON/INPNOF_GENERALINF/ICOEF,ISOFTMAX,IORBOPT,MAXIT,MAXIT21      
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      COMMON/INP_EFIELDL/EX,EY,EZ,EFIELDL
!
      DOUBLE PRECISION,DIMENSION(3)::DIPN
      DOUBLE PRECISION,DIMENSION(NV)::GAMMA,GRAD
      DOUBLE PRECISION,DIMENSION(NBF5)::RO,HCORE,DIPx,DIPy,DIPz
      DOUBLE PRECISION,DIMENSION(NBFT5)::QJ,QK
      DOUBLE PRECISION,DIMENSION(NBF5,NV)::DR
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::ADIPx,ADIPy,ADIPz
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5)::CJ12,CK12
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5,NV)::DCJ12r,DCK12r
      DOUBLE PRECISION,DIMENSION(NBF,NBF,NBF)::QD
!----------------------------------------------------------------------- 
!     Occupations
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(ISOFTMAX==0)THEN
       CALL OCUPACIONTFr(GAMMA,RO,CJ12,CK12,DR,DCJ12r,DCK12r,NV)
      ELSE IF(ISOFTMAX==1)THEN
       CALL OCUPACIONSMr(GAMMA,RO,CJ12,CK12,DR,DCJ12r,DCK12r,NV)
      ENDIF
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!                                ENERGY
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      IF(MODE==0.or.MODE==2)THEN
       IF(NCWO==1)THEN             ! PNOFi(1): Perfect Pairing (NCWO=1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ENERGY = 0.0d0
        do in=1,NO1
         ENERGY = ENERGY + RO(in)*(2.0*HCORE(in)+QJ(in*(in+1)/2))       &
                + PRODCWQWj(in,CJ12,QJ) - PRODCWQWj(in,CK12,QK)          
        enddo                                                            
        do i=1,NDOC                                                      
         in = NO1+i                                                      
         ENERGY = ENERGY + RO(in)*(2.0*HCORE(in)+QJ(in*(in+1)/2))       &
                + PRODCWQWj(in,CJ12,QJ) - PRODCWQWj(in,CK12,QK)          
         in = na+ndoc-i+1                                                
         ENERGY = ENERGY + RO(in)*(2.0*HCORE(in)+QJ(in*(in+1)/2))       &
                + PRODCWQWj(in,CJ12,QJ) - PRODCWQWj(in,CK12,QK)          
        enddo                                                            
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       ELSE                        ! PNOFi(Nc): Extended PNOF (NCWO>1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ENERGY = 0.0d0
        do in=1,NO1
         ENERGY = ENERGY + RO(in)*(2.0*HCORE(in)+QJ(in*(in+1)/2))       &
                + PRODCWQWj(in,CJ12,QJ) - PRODCWQWj(in,CK12,QK)          
        enddo                                                            
        do i=1,NDOC                                                      
         in = NO1+i                                                      
         ENERGY = ENERGY + RO(in)*(2.0*HCORE(in)+QJ(in*(in+1)/2))       &
                + PRODCWQWj(in,CJ12,QJ) - PRODCWQWj(in,CK12,QK)          
         do iw=1,ncwo-1                                                  
          !in = na+ncwo*(ndoc-i)+iw              !old-sort                         
          in = no1+(na-nb)+ndoc*(iw+1)-i+1       !new-sort
          ENERGY = ENERGY + RO(in)*(2.0*HCORE(in)+QJ(in*(in+1)/2))      &
                 + PRODCWQWj(in,CJ12,QJ) - PRODCWQWj(in,CK12,QK)         
         enddo                                                           
         !in = na+ncwo*(ndoc-i)+ncwo             !old-sort                         
         in = no1+(na-nb)+ndoc*(ncwo+1)-i+1      !new-sort
         ENERGY = ENERGY + RO(in)*(2.0*HCORE(in)+QJ(in*(in+1)/2))       &
                + PRODCWQWj(in,CJ12,QJ) - PRODCWQWj(in,CK12,QK)          
        enddo                                                            
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       ENDIF
       IF(NSOC>0)THEN                                                   
        do i=NDOC+1,NDNS                                                
         in = NO1+i       
         ENERGY = ENERGY + 2.0d0*RO(in)*HCORE(in)                      &
                + PRODCWQWj(in,CJ12,QJ) - PRODCWQWj(in,CK12,QK)
        enddo
       ENDIF
!- - - including Electric Field  - - - - - - - - - - - - - - - - -
       if(EFIELDL)then
        CALL DIPMOMr(DIPN,ADIPx,ADIPy,ADIPz,DIPx,DIPy,DIPz,             &
                     QD,RO,DMXe,DMYe,DMZe,DMX,DMY,DMZ,DM)
        ENERGY = ENERGY - EX*DMX - EY*DMY - EZ*DMZ
       endif
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ENDIF
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!                               GRADIENTS
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      IF(MODE==1.or.MODE==2)THEN
       GRAD = 0.0d0      
       IF(NCWO==1)THEN             ! PNOFi(1): Perfect Pairing (NCWO=1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        DO ig=1,NV
         do i=1,NDOC
          in = NO1+i
          GRAD(ig) = GRAD(ig)                                           &
                   + DR(in,ig) * ( 2.0d0*HCORE(in) + QJ(in*(in+1)/2) )  &
                   + 2.0d0 * ( PRODCWQWjk(nv,in,ig,DCJ12r,QJ)           &
                             - PRODCWQWjk(nv,in,ig,DCK12r,QK) )          
          in = na+ndoc-i+1                                               
          GRAD(ig) = GRAD(ig)                                           &
                   + DR(in,ig) * ( 2.0d0*HCORE(in) + QJ(in*(in+1)/2) )  &
                   + 2.0d0 * ( PRODCWQWjk(nv,in,ig,DCJ12r,QJ)           &
                             - PRODCWQWjk(nv,in,ig,DCK12r,QK) )
         enddo
        ENDDO
!- - -  including Electric Field  - - - - - - - - - - - - - - - -
        if(EFIELDL)then   
         DO ig=1,NV
          do i=1,NDOC
           in = NO1+i
           GRAD(ig) = GRAD(ig) + 2.0d0 * DR(in,ig)                      &
                    * ( EX*DIPx(in) + EY*DIPy(in) + EZ*DIPz(in) )        
           in = na+ndoc-i+1                                              
           GRAD(ig) = GRAD(ig) + 2.0d0 * DR(in,ig)                      &
                    * ( EX*DIPx(in) + EY*DIPy(in) + EZ*DIPz(in) )
          enddo
         ENDDO
        endif
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       ELSE                        ! PNOFi(Nc): Extended PNOF (NCWO>1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        DO ig=1,NV
         do i=1,NDOC
          in = NO1+i
          GRAD(ig) = GRAD(ig)                                           &
                   + DR(in,ig) * ( 2.0d0*HCORE(in) + QJ(in*(in+1)/2) )  &
                   + 2.0d0 * ( PRODCWQWjk(nv,in,ig,DCJ12r,QJ)           &
                             - PRODCWQWjk(nv,in,ig,DCK12r,QK) )          
          do iw=1,NCWO-1                                                 
           !in = na+ncwo*(ndoc-i)+iw             !old-sort                         
           in = no1+(na-nb)+ndoc*(iw+1)-i+1      !new-sort
           GRAD(ig) = GRAD(ig)                                          &
                    + DR(in,ig) * ( 2.0d0*HCORE(in) + QJ(in*(in+1)/2) ) &
                    + 2.0d0 * ( PRODCWQWjk(nv,in,ig,DCJ12r,QJ)          &
                              - PRODCWQWjk(nv,in,ig,DCK12r,QK) )         
          enddo                                                          
          !in = na+ncwo*(ndoc-i)+ncwo             !old-sort                        
          in = no1+(na-nb)+ndoc*(ncwo+1)-i+1      !new-sort
          GRAD(ig) = GRAD(ig)                                           &
                   + DR(in,ig) * ( 2.0d0*HCORE(in) + QJ(in*(in+1)/2) )  &
                   + 2.0d0 * ( PRODCWQWjk(nv,in,ig,DCJ12r,QJ)           &
                             - PRODCWQWjk(nv,in,ig,DCK12r,QK) )
         enddo
        ENDDO
!- - -  including Electric Field  - - - - - - - - - - - - - - - -
        if(EFIELDL)then   
         DO ig=1,NV
          do i=1,NDOC
           in = NO1+i
           GRAD(ig) = GRAD(ig) + 2.0d0 * DR(in,ig)                      &
                    * ( EX*DIPx(in) + EY*DIPy(in) + EZ*DIPz(in) )        
           do iw=1,NCWO-1                                                
            !in = na+ncwo*(ndoc-i)+iw             !old-sort                                     
            in = no1+(na-nb)+ndoc*(iw+1)-i+1      !new-sort
            GRAD(ig) = GRAD(ig) + 2.0d0 * DR(in,ig)                     &
                     * ( EX*DIPx(in) + EY*DIPy(in) + EZ*DIPz(in) )       
           enddo                                                         
           !in = na+ncwo*(ndoc-i)+ncwo             !old-sort                       
           in = no1+(na-nb)+ndoc*(ncwo+1)-i+1      !new-sort
           GRAD(ig) = GRAD(ig) + 2.0d0 * DR(in,ig)                      &
                    * ( EX*DIPx(in) + EY*DIPy(in) + EZ*DIPz(in) )
          enddo
         ENDDO
        endif
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       ENDIF
      ENDIF
!-----------------------------------------------------------------------
      RETURN
      END
      
! OCUPENERGYro
      SUBROUTINE OCUPENERGYro(MODE,GAMMA,RO,CJ12,CK12,DR,DCJ12r,DCK12r, &
                              QD,HCORE,QJ,QK,DIPN,ADIPx,ADIPy,ADIPz,    &
                              DIPx,DIPy,DIPz,ENERGY,GRAD,NV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      LOGICAL EFIELDL,HighSpin
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INP_EFIELDL/EX,EY,EZ,EFIELDL
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      COMMON/INPNOF_ORBSPACE0/NO1,NDOC,NCO,NCWO,NVIR,NAC,NO0
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
      COMMON/INPNOF_GENERALINF/ICOEF,ISOFTMAX,IORBOPT,MAXIT,MAXIT21            
!
      DOUBLE PRECISION,DIMENSION(3)::DIPN
      DOUBLE PRECISION,DIMENSION(NV)::GAMMA,GRAD
      DOUBLE PRECISION,DIMENSION(NBF5)::RO,HCORE,DIPx,DIPy,DIPz
      DOUBLE PRECISION,DIMENSION(NBFT5)::QJ,QK
      DOUBLE PRECISION,DIMENSION(NBF5,NV)::DR
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::ADIPx,ADIPy,ADIPz
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5)::CJ12,CK12
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5,NV)::DCJ12r,DCK12r
      DOUBLE PRECISION,DIMENSION(NBF,NBF,NBF)::QD
!----------------------------------------------------------------------- 
!     Occupations
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(ISOFTMAX==0)THEN
       CALL OCUPACIONTFr(GAMMA,RO,CJ12,CK12,DR,DCJ12r,DCK12r,NV)
      ELSE IF(ISOFTMAX==1)THEN
       CALL OCUPACIONSMr(GAMMA,RO,CJ12,CK12,DR,DCJ12r,DCK12r,NV)
      ENDIF
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!                                ENERGY
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      IF(MODE==0.or.MODE==2)THEN
       IF(NCWO==1)THEN             ! PNOFi(1): Perfect Pairing (NCWO=1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ENERGY = 0.0d0
        do in=1,NO1
         ENERGY = ENERGY + RO(in)*(2.0*HCORE(in)+QJ(in*(in+1)/2))       &
                + PRODCWQWj1(in,CJ12,QJ) - PRODCWQWj1(in,CK12,QK)       &
                + 2.0d0*PRODROQWj1(in,RO,QJ)-PRODROQWj1(in,RO,QK)        
        enddo                                                            
        do i=1,NDOC                                                      
         in = NO1+i                                                      
         ENERGY = ENERGY + RO(in)*(2.0*HCORE(in)+QJ(in*(in+1)/2))       &
                + PRODCWQWj1(in,CJ12,QJ) - PRODCWQWj1(in,CK12,QK)       &
                + 2.0d0*PRODROQWj1(in,RO,QJ)-PRODROQWj1(in,RO,QK)           
         in = na+ndoc-i+1                                                
         ENERGY = ENERGY + RO(in)*(2.0*HCORE(in)+QJ(in*(in+1)/2))       &
                + PRODCWQWj2(in,CJ12,QJ) - PRODCWQWj2(in,CK12,QK)       &
                + 2.0d0*PRODROQWj2(in,RO,QJ)-PRODROQWj2(in,RO,QK)     
        enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        do i=NDOC+1,NDNS
         in = NO1+i
         ENERGY = ENERGY + RO(in)*HCORE(in)                             &
                + 0.5d0*(PRODROQWj0(in,RO,QJ)-PRODROQWj0(in,RO,QK))
        enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       ELSE                        ! PNOFi(Nc): Extended PNOF (NCWO>1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ENERGY = 0.0d0
        do in=1,NO1
         ENERGY = ENERGY + RO(in)*(2.0*HCORE(in)+QJ(in*(in+1)/2))       &
                + PRODCWQWj1(in,CJ12,QJ) - PRODCWQWj1(in,CK12,QK)       &
                + 2.0d0*PRODROQWj1(in,RO,QJ)-PRODROQWj1(in,RO,QK)        
        enddo                                                            
        do i=1,NDOC                                                      
         in = NO1+i                                                      
         ENERGY = ENERGY + RO(in)*(2.0*HCORE(in)+QJ(in*(in+1)/2))       &
                + PRODCWQWj1(in,CJ12,QJ) - PRODCWQWj1(in,CK12,QK)       &
                + 2.0d0*PRODROQWj1(in,RO,QJ)-PRODROQWj1(in,RO,QK)        
         do iw=1,ncwo-1                                                  
          !in = na+ncwo*(ndoc-i)+iw             !old-sort                              
          in = no1+(na-nb)+ndoc*(iw+1)-i+1      !new-sort
          ENERGY = ENERGY + RO(in)*(2.0*HCORE(in)+QJ(in*(in+1)/2))      &
                 + PRODCWQWj2(in,CJ12,QJ) - PRODCWQWj2(in,CK12,QK)      &
                + 2.0d0*PRODROQWj2(in,RO,QJ)-PRODROQWj2(in,RO,QK)           
         enddo                                                           
         !in = na+ncwo*(ndoc-i)+ncwo             !old-sort                         
         in = no1+(na-nb)+ndoc*(ncwo+1)-i+1      !new-sort
         ENERGY = ENERGY + RO(in)*(2.0*HCORE(in)+QJ(in*(in+1)/2))       &
                + PRODCWQWj2(in,CJ12,QJ) - PRODCWQWj2(in,CK12,QK)       &
                + 2.0d0*PRODROQWj2(in,RO,QJ)-PRODROQWj2(in,RO,QK)          
        enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        do i=NDOC+1,NDNS
         in = NO1+i
         ENERGY = ENERGY + RO(in)*HCORE(in)                             &
                + 0.5d0*(PRODROQWj0(in,RO,QJ)-PRODROQWj0(in,RO,QK))
        enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       ENDIF
!- - - including Electric Field  - - - - - - - - - - - - - - - - -
       if(EFIELDL)then
        CALL DIPMOMr(DIPN,ADIPx,ADIPy,ADIPz,DIPx,DIPy,DIPz,             &
                     QD,RO,DMXe,DMYe,DMZe,DMX,DMY,DMZ,DM)
        ENERGY = ENERGY - EX*DMX - EY*DMY - EZ*DMZ
       endif
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ENDIF
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!                               GRADIENTS
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      IF(MODE==1.or.MODE==2)THEN
       IF(NCWO==1)THEN             ! PNOFi(1): Perfect Pairing (NCWO=1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        GRAD = 0.0d0
        DO ig=1,NV
         do i=1,NDOC
          in = NO1+i
          GRAD(ig) = GRAD(ig)                                           &
                   + DR(in,ig) * ( 2.0d0*HCORE(in) + QJ(in*(in+1)/2) )  &
                   + 2.0d0 * ( PRODCWQWjk1(nv,in,ig,DCJ12r,QJ)          &
                             - PRODCWQWjk1(nv,in,ig,DCK12r,QK) )        &
                   + 2.0d0 *   PRODDRQWjk1(nv,in,ig,DR,QJ)              &
                   -           PRODDRQWjk1(nv,in,ig,DR,QK)               
          in = na+ndoc-i+1                                               
          GRAD(ig) = GRAD(ig)                                           &
                   + DR(in,ig) * ( 2.0d0*HCORE(in) + QJ(in*(in+1)/2) )  &
                   + 2.0d0 * ( PRODCWQWjk2(nv,in,ig,DCJ12r,QJ)          &
                             - PRODCWQWjk2(nv,in,ig,DCK12r,QK) )        &
                   + 2.0d0 *   PRODDRQWjk2(nv,in,ig,DR,QJ)              &
                   -           PRODDRQWjk2(nv,in,ig,DR,QK)
         enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        ENDDO
! - - -  including Electric Field  - - - - - - - - - - - - - - - - - - - 
        if(EFIELDL)then   
         DO ig=1,NV
          do i=1,NDOC
           in = NO1+i
           GRAD(ig) = GRAD(ig) + 2.0d0 * DR(in,ig)                      &
                    * ( EX*DIPx(in) + EY*DIPy(in) + EZ*DIPz(in) )        
           in = na+ndoc-i+1                                              
           GRAD(ig) = GRAD(ig) + 2.0d0 * DR(in,ig)                      &
                    * ( EX*DIPx(in) + EY*DIPy(in) + EZ*DIPz(in) )
          enddo
         ENDDO
        endif
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       ELSE                        ! PNOFi(Nc): Extended PNOF (NCWO>1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        GRAD = 0.0d0
        DO ig=1,NV
         do i=1,NDOC
          in = NO1+i
          GRAD(ig) = GRAD(ig)                                           &
                   + DR(in,ig) * ( 2.0d0*HCORE(in) + QJ(in*(in+1)/2) )  &
                   + 2.0d0 * ( PRODCWQWjk1(nv,in,ig,DCJ12r,QJ)          &
                             - PRODCWQWjk1(nv,in,ig,DCK12r,QK) )        &
                   + 2.0d0 *   PRODDRQWjk1(nv,in,ig,DR,QJ)              &
                   -           PRODDRQWjk1(nv,in,ig,DR,QK)               
          do iw=1,NCWO-1                                                 
           !in = na+ncwo*(ndoc-i)+iw             !old-sort                                     
           in = no1+(na-nb)+ndoc*(iw+1)-i+1      !new-sort
           GRAD(ig) = GRAD(ig)                                          &
                    + DR(in,ig) * ( 2.0d0*HCORE(in) + QJ(in*(in+1)/2) ) &
                    + 2.0d0 * ( PRODCWQWjk2(nv,in,ig,DCJ12r,QJ)         &
                              - PRODCWQWjk2(nv,in,ig,DCK12r,QK) )       &
                    + 2.0d0 *   PRODDRQWjk2(nv,in,ig,DR,QJ)             &
                    -           PRODDRQWjk2(nv,in,ig,DR,QK)              
          enddo                                                          
          !in = na+ncwo*(ndoc-i)+ncwo             !old-sort                        
          in = no1+(na-nb)+ndoc*(ncwo+1)-i+1      !new-sort
          GRAD(ig) = GRAD(ig)                                           &
                   + DR(in,ig) * ( 2.0d0*HCORE(in) + QJ(in*(in+1)/2) )  &
                   + 2.0d0 * ( PRODCWQWjk2(nv,in,ig,DCJ12r,QJ)          &
                             - PRODCWQWjk2(nv,in,ig,DCK12r,QK) )        &
                   + 2.0d0 *   PRODDRQWjk2(nv,in,ig,DR,QJ)              &
                   -           PRODDRQWjk2(nv,in,ig,DR,QK)
         enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -          
        ENDDO
! - - -  including Electric Field  - - - - - - - - - - - - - - - - - - -
        if(EFIELDL)then   
         DO ig=1,NV
          do i=1,NDOC
           in = NO1+i
           GRAD(ig) = GRAD(ig) + 2.0d0 * DR(in,ig)                      &
                    * ( EX*DIPx(in) + EY*DIPy(in) + EZ*DIPz(in) )        
           do iw=1,NCWO-1                                                
            !in = na+ncwo*(ndoc-i)+iw             !old-sort                      
            in = no1+(na-nb)+ndoc*(iw+1)-i+1      !new-sort
            GRAD(ig) = GRAD(ig) + 2.0d0 * DR(in,ig)                     &
                     * ( EX*DIPx(in) + EY*DIPy(in) + EZ*DIPz(in) )       
           enddo                                                         
           !in = na+ncwo*(ndoc-i)+ncwo             !old-sort                   
           in = no1+(na-nb)+ndoc*(ncwo+1)-i+1      !new-sort
           GRAD(ig) = GRAD(ig) + 2.0d0 * DR(in,ig)                      &
                    * ( EX*DIPx(in) + EY*DIPy(in) + EZ*DIPz(in) )
          enddo
         ENDDO
        endif
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       ENDIF
      ENDIF
!-----------------------------------------------------------------------
      RETURN
      END

! LBFGSOCUP
      SUBROUTINE LBFGSOCUP(NV,GAMMA,USER,ENERGY)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL HighSpin
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
      COMMON /LB3/MP,LP,GTOL,STPMIN,STPMAX
      COMMON/PUNTEROSUSER/N1,N2,N3,N4,N5,N6,N7,N8,N9,N10,N11,N12,N13,   &
                          N14,N15,N16,N17,N18,N19,N20,N21,N22,N23,N24,  &
                          N25,N26,N27,N28,N29,N30,N31,N32,N33,N34,N35,  &
                          N36,N37,N38,N39,N40,N41,N42,N43,N44,N45,N46,  &
                          N47,N48,N49,N50,N51,NUSER
!
      LOGICAL::DIAGCO
      INTEGER,PARAMETER::MSAVE=7
      INTEGER::IFLAG,ICALL,N,M,MP,LP,NWORK
      INTEGER,DIMENSION(2)::IPRINT
      DOUBLE PRECISION :: F,EPS,XTOL,GTOL,STPMIN,STPMAX
      DOUBLE PRECISION :: X(NV),G(NV),DIAG(NV)
      DOUBLE PRECISION,DIMENSION(NV) :: GAMMA
      DOUBLE PRECISION,DIMENSION(NUSER) :: USER
      DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::W
!     The driver for LBFGS must always declare LB2 as EXTERNAL
      EXTERNAL LB2
!-----------------------------------------------------------------------
      NWORK=NV*(2*MSAVE +1)+2*MSAVE
      ALLOCATE(W(NWORK))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Calling to LBFGS SUBROUTINE
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     CHECK COMMON /LB3, MP SETS UNIT WHERE PRINTING OPTIMIZATION INFO, 
!     AND LP WHERE PRINTING INFO ABOUT ERRORS
      N=NV ! NUMBER OF VARIABLES
      M=5     ! 0 <= M <= 7
      IPRINT(1)= -1
!     IPRINT(1) < 0 : no output is generated,
!     IPRINT(1) = 0 : output only at first and last iteration,
!     IPRINT(1) > 0 : output every IPRINT(1) iterations.      
      IPRINT(2)= 0
!     IPRINT(2) = 0 : iteration count, number of function 
!                     evaluations, function value, norm of the
!                     gradient, and steplength,
!     IPRINT(2) = 1 : same as IPRINT(2)=0, plus vector of
!                     variables and  gradient vector at the
!                     initial point,
!     IPRINT(2) = 2 : same as IPRINT(2)=1, plus vector of
!                     variables,
!     IPRINT(2) = 3 : same as IPRINT(2)=2, plus gradient vector.      
!
!     We do not wish to provide the diagonal matrices Hk0, and 
!     therefore set DIAGCO to FALSE.
      DIAGCO= .FALSE.
      EPS= 1.0D-5
      XTOL= 1.0D-16
      ICALL=0
      IFLAG=0
      X = GAMMA ! INITIAL ESTIMATE OF THE SOLUTION VECTOR
      MODE = 2
      DO
       if(NSOC==0)then
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!       Singlet State (S=0,Ms=0) and Multiplet States (S>0,Ms=0)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
        CALL OCUPENERGYrc(MODE,X,USER(N1),USER(N2),USER(N3),            &
                          USER(N4),USER(N5),USER(N6),USER(N7),          &
                          USER(N8),USER(N9),USER(N10),USER(N11),        &
                          USER(N12),USER(N13),USER(N14),USER(N15),      &
                          USER(N16),USER(N17),ENERGY,G,NV)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --     
       else
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!       High-Spin Multiplet State (S>0,Ms=S)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
        CALL OCUPENERGYro(MODE,X,USER(N1),USER(N2),USER(N3),            &
                          USER(N4),USER(N5),USER(N6),USER(N7),          &
                          USER(N8),USER(N9),USER(N10),USER(N11),        &
                          USER(N12),USER(N13),USER(N14),USER(N15),      &
                          USER(N16),USER(N17),ENERGY,G,NV)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
       end if
!      F CONTAINS THE VALUE OF THE FUNCTION AT THE POINT X
!      G CONTAINS THE COMPONENTS OF GRADIENT AT X
       F = ENERGY
       CALL LBFGS(N,M,X,F,G,DIAGCO,DIAG,IPRINT,EPS,XTOL,W,IFLAG)
       GAMMA = X
       IF(IFLAG.LE.0) EXIT
       ICALL=ICALL + 1
!      We allow at most 1000 evaluations of F and G
       IF(ICALL.GT.1000) EXIT
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --       
      ENDDO
      DEALLOCATE(W)
!     FINAL CALL TO COMPUTE ENERGY AT EQUILIBRIUM POINT
      MODE = 1
      if(NSOC==0)then
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!      Singlet State (S=0,Ms=0) and Multiplet States (S>0,Ms=0)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
       CALL OCUPENERGYrc(MODE,X,USER(N1),USER(N2),USER(N3),             &
                         USER(N4),USER(N5),USER(N6),USER(N7),           &
                         USER(N8),USER(N9),USER(N10),USER(N11),         &
                         USER(N12),USER(N13),USER(N14),USER(N15),       &
                         USER(N16),USER(N17),ENERGY,G,NV)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --     
      else
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!      High-Spin Multiplet State (S>0,Ms=S)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
       CALL OCUPENERGYro(MODE,X,USER(N1),USER(N2),USER(N3),             &
                         USER(N4),USER(N5),USER(N6),USER(N7),           &
                         USER(N8),USER(N9),USER(N10),USER(N11),         &
                         USER(N12),USER(N13),USER(N14),USER(N15),       &
                         USER(N16),USER(N17),ENERGY,G,NV)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
      end if
!-----------------------------------------------------------------------
      RETURN
      END

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
!                                                                      !
!          Subroutines related to occupation optimization of           !
!          spin-compensated systems using different PNOFs              !
!                                                                      !
!  CJCKD3: Define the coefficientes in front J,K,L integrals for PNOF3 !
!          ( J. Chem. Phys. 132, 031103, 2010 )                        !
!  CJCKD4: Define the coefficientes in front J,K,L integrals for PNOF4 !
!          ( J. Chem. Phys. 133, 111101, 2010 )                        !
!  CJCKD5: Define the coefficientes in front J,K,L integrals for PNOF5 !
!          ( J. Chem. Phys. 134, 164102, 2011; 139, 234109, 2013 )     !
!  CJCKD6: Define the coefficientes in front J,K,L integrals for PNOF6 !
!          ( J. Chem. Phys. 141, 044107, 2014 )                        !
!  AACOMP: Add parallel spin components in PNOFi (i=4,5,6) if Daa=Dab  !
!  CJCKD7: Define the coefficientes in front J,K,L integrals for PNOF7 !
!          ( PRL 119, 063002, 2017; PRA 100, 032508, 2019 )            !
!  CJCKD8: Define the coefficientes in front J,K,L integrals for GNOF  !
!          ( Phys. Rev. Lett. 127, 233001, 2021 )                      !
!          Define the coefficientes in front J,K,L integrals for GNOFm !
!          ( Phys. Rev. Lett. 134, 206401, 2025 )                      !
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !

! CJCKD3 = PNOF3 + pairing conditions
      SUBROUTINE CJCKD3(NV,RO,DR,BETA,DB,CJ12,CK12,DCJ12r,DCK12r)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      COMMON/SUMSZ/SUMS,SUMF
!
      DOUBLE PRECISION,DIMENSION(NBF5)::RO,BETA
      DOUBLE PRECISION,DIMENSION(NBF5,NV)::DR,DB
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5)::CJ12,CK12
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5,NV)::DCJ12r,DCK12r
      ALLOCATABLE::DSUMS(:)
!-----------------------------------------------------------------------
!                  CJ12, CK12, DCJ12r and DCK12r
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!     CJpq = 2NpNq, CKpq = SQRT(NpNq)  (Note: in PNOF3 Daa = 0.0)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO j=1,NBF5
       DO i=1,NBF5
        CJ12(j,i) = 2.0d0*RO(j)*RO(i)
        CK12(j,i) = BETA(j)*BETA(i)
        do k=1,nv        
         DCJ12r(j,i,k) = 2.0d0*DR(j,k)*RO(i)
         DCK12r(j,i,k) = DB(j,k)*BETA(i)
        enddo
       ENDDO
      ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Derivatives of SUMS, FSUMS = (1-S)/S, and its derivatives
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(SUMS<=1.0d0)THEN
       ALLOCATE(DSUMS(NV))
       do k=1,nv
        DSUMS(k) = 0.0d0
        DO j=1,NB
         DSUMS(k) = DSUMS(k) - DR(j,k)
        ENDDO
       enddo
       if(SUMS>1.0d-20)then
        FSUMS = 1.0d0/SUMS - 1.0d0
        do k=1,nv
        DSUMS(k) = DSUMS(k)/(SUMS*SUMS)           ! DSUMS is changed
       enddo
       else
        FSUMS = 0.0d0
        DSUMS = 0.0d0
       endif
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Including interactions
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(SUMS<=1.0d0)THEN
       DO j=1,NB
        Hj = 1.0d0-RO(j)
        DO i=1,NB
         Hi = 1.0d0-RO(i)
         CJ12(j,i) = CJ12(j,i) - Hj*Hi
         do k=1,nv
          DCJ12r(j,i,k) = DCJ12r(j,i,k) + DR(j,k)*Hi
         enddo
        ENDDO
        DO i=NB+1,NBF5
         CJ12(j,i) = CJ12(j,i) - FSUMS*Hj*RO(i)
         CK12(j,i) = CK12(j,i) + DSQRT(Hj*RO(i))
         do k=1,nv        
          DCJ12r(j,i,k) = DCJ12r(j,i,k) + DSUMS(k)*Hj*RO(i)             &
                        + FSUMS*DR(j,k)*RO(i)                            
          if(Hj>0.0d0)DCK12r(j,i,k) = DCK12r(j,i,k)                     &
                      - 0.5d0*DSQRT(Hj*RO(i))*DR(j,k)/Hj                 
         enddo                                                           
        ENDDO                                                            
       ENDDO                                                             
       DO j=NB+1,NBF5                                                    
        DO i=1,NB                                                        
         Hi = 1.0d0-RO(i)                                                
         CJ12(j,i) = CJ12(j,i) - FSUMS*RO(j)*Hi                          
         CK12(j,i) = CK12(j,i) + DSQRT(RO(j)*Hi)                         
         do k=1,nv                                                       
          DCJ12r(j,i,k) = DCJ12r(j,i,k) - FSUMS*DR(j,k)*Hi               
          if(RO(j)>0.0d0)DCK12r(j,i,k) = DCK12r(j,i,k)                  &
                         + 0.5d0*DSQRT(RO(j)*Hi)*DR(j,k)/RO(j)
         enddo
        ENDDO
        DO i=NB+1,NBF5
         CJ12(j,i) = CJ12(j,i) - RO(j)*RO(i)
         CK12(j,i) = - CK12(j,i)
         do k=1,nv        
          DCJ12r(j,i,k) = DCJ12r(j,i,k) - DR(j,k)*RO(i)
          DCK12r(j,i,k) = - DCK12r(j,i,k)
         enddo
        ENDDO
       ENDDO
      ELSE
       WRITE(6,*)'STOP PNOF3: SUMS > 1'
       STOP
      ENDIF
!-----------------------------------------------------------------------
      IF(SUMS<=1.0d0)DEALLOCATE(DSUMS)
      RETURN
      END

! CJCKD4 = PNOF4 + pairing conditions + SUMS>1 implementation as PNOF6
      SUBROUTINE CJCKD4(NV,RO,DR,BETA,DB,CJ12,CK12,DCJ12r,DCK12r)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_PNOF/IPNOF,NTWOPAR
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      COMMON/SUMSZ/SUMS,SUMF
!
      DOUBLE PRECISION,DIMENSION(NBF5)::RO,BETA
      DOUBLE PRECISION,DIMENSION(NBF5,NV)::DR,DB
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5)::CJ12,CK12
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5,NV)::DCJ12r,DCK12r
      ALLOCATABLE::DSUMS(:),G(:),FI(:),DG(:,:),DF(:,:),DSUMF(:),DSUMG(:)
!-----------------------------------------------------------------------
!          Alpha-Beta Components of CJ12, CK12, DCJ12 and DCK12
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!                         CJpq = NpNq, CKpq = 0
!-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
      CK12 = 0.0d0
      DCK12r = 0.0d0
      DO j=1,NBF5
       DO i=1,NBF5
        CJ12(j,i) = RO(j)*RO(i)
        do k=1,nv        
         DCJ12r(j,i,k) = DR(j,k)*RO(i)
        enddo
       ENDDO
      ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Derivatives of SUMS
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(DSUMS(NV))
      do k=1,nv
       DSUMS(k) = 0.0d0
       DO j=1,NB
        DSUMS(k) = DSUMS(k) - DR(j,k)
       ENDDO
      enddo
!-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
!                 Including interactions in CJ and CK
!-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
      IF(SUMS<=1.0d0)THEN                                  ! SUMS <= 1
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!      FSUMS = (1-S)/S, and its derivatives
!- - - - - - - - - - - - - - - - - - - - - - 
       if(SUMS>1.0d-20)then
        FSUMS = 1.0d0/SUMS - 1.0d0
        do k=1,nv
!        DSUMS is changed
         DSUMS(k) = DSUMS(k)/(SUMS*SUMS)
        enddo
       else
        FSUMS = 0.0d0
        DSUMS = 0.0d0
       endif
!- - - - - - - - - - - - - - - - - - - - - - 
       DO j=1,NB
        Hj = 1.0d0-RO(j)
        DO i=1,NB
         Hi = 1.0d0-RO(i)
         CJ12(j,i) = CJ12(j,i) - Hj*Hi
         CK12(j,i) = DSQRT(Hj*Hi)
         do k=1,nv        
          DCJ12r(j,i,k) = DCJ12r(j,i,k) + DR(j,k)*Hi
          if(DABS(Hj)>0.0d0)DCK12r(j,i,k)=-0.5d0*CK12(j,i)*DR(j,k)/Hj
         enddo
        ENDDO
        DO i=NB+1,NBF5
         CJ12(j,i) = CJ12(j,i) - FSUMS*Hj*RO(i)
         CK12ji = Hj*RO(i)/SUMS
         CK12(j,i) = DSQRT( CK12ji*(RO(j)-RO(i)) + CK12ji*CK12ji )
         do k=1,nv        
          DCJ12r(j,i,k) = DCJ12r(j,i,k) + DSUMS(k)*Hj*RO(i)             &
                        + FSUMS*DR(j,k)*RO(i)                            
          DCK12jik = - DR(j,k)*RO(i)/SUMS - DSUMS(k)*Hj*RO(i)            
          if(DABS(CK12(j,i))>0.0d0)                                     &
           DCK12r(j,i,k) = ( DCK12jik*(RO(j)-RO(i)) + CK12ji*DR(j,k)    &
                         + 2.0d0*CK12ji*DCK12jik ) / (2.0d0*CK12(j,i))   
         enddo                                                           
        ENDDO                                                            
       ENDDO                                                             
       DEALLOCATE(DSUMS)                                                 
       DO j=NB+1,NBF5                                                    
        DO i=1,NB                                                        
         Hi = 1.0d0-RO(i)                                                
         CJ12(j,i) = CJ12(j,i) - FSUMS*RO(j)*Hi                          
         CK12ji = RO(j)*Hi/SUMS                                          
         CK12(j,i) = DSQRT( CK12ji*(RO(i)-RO(j)) + CK12ji*CK12ji )       
         do k=1,nv                                                       
          DCJ12r(j,i,k) = DCJ12r(j,i,k) - FSUMS*DR(j,k)*Hi               
          DCK12jik = DR(j,k)*Hi/SUMS                                     
          if(DABS(CK12(j,i))>0.0d0)                                     &
           DCK12r(j,i,k) = ( DCK12jik*(RO(i)-RO(j)) - CK12ji*DR(j,k)    &
                         + 2.0d0*CK12ji*DCK12jik ) / (2.0d0*CK12(j,i))
         enddo
        ENDDO
        DO i=NB+1,NBF5
         CJ12(j,i) = 0.0d0
         CK12(j,i) = - BETA(j)*BETA(i)
         do k=1,nv        
          DCJ12r(j,i,k) = 0.0d0
          DCK12r(j,i,k) = - DB(j,k)*BETA(i)
         enddo
        ENDDO
       ENDDO
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
      ELSE                                                 ! SUMS > 1
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!      Gq, FIq = Nq*Hq + Gq*Gq - Gq*S, SUMF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       ALLOCATE(G(NBF5),DG(NBF5,NV))
       RHc = 0.5d0
       SUMS1 = 1.0d0-SUMS                  ! SUMS1 is negative
       DO j=1,NB
        Hj = 1.0d0 - RO(j)
        G(j) = HEAV(SUMS1)*Hj + (1.0d0-HEAV(SUMS1))*Hj*DEXP(-Hj/RHc)
        do k=1,nv
         DG(j,k) = -HEAV(SUMS1)*DR(j,k) - (1.0d0-HEAV(SUMS1))           &
                 * (1.0d0-Hj/RHc)*DR(j,k)*DEXP(-Hj/RHc)                 &
                 + (DEXP(-Hj/RHc)-1.0d0)*Hj*DHEAV(SUMS1)*DSUMS(k)        
        enddo                                                            
       ENDDO                                                             
       DO j=NB+1,NBF5                                                    
        G(j) = HEAV(SUMS1)*RO(j)                                        &
             + (1.0d0-HEAV(SUMS1))*RO(j)*DEXP(-RO(j)/RHc)                
        do k=1,nv                                                        
         DG(j,k) = HEAV(SUMS1)*DR(j,k) + (1.0d0-HEAV(SUMS1))            &
                 * (1.0d0-RO(j)/RHc)*DR(j,k)*DEXP(-RO(j)/RHc)           &
                 +(DEXP(-RO(j)/RHc)-1.0d0)*RO(j)*DHEAV(SUMS1)*DSUMS(k)
        enddo
       ENDDO
       DEALLOCATE(DSUMS)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       SUMG = 0.0d0
       DO j=NB+1,NBF5
        SUMG = SUMG + G(j) 
       ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       ALLOCATE(DSUMG(NV))
       do k=1,nv
        DSUMG(k) = 0.0d0
        DO j=NB+1,NBF5
         DSUMG(k) = DSUMG(k) + DG(j,k)
        ENDDO
       enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       ALLOCATE(FI(NBF5),DF(NBF5,NV))
       DO j=1,NBF5
        FI(j) = RO(j)*(1.0d0-RO(j)) + G(j)*(G(j)-SUMG)
        do k=1,nv        
         DF(j,k) = (1.0d0-2.0d0*RO(j))*DR(j,k) + DG(j,k)*(G(j)-SUMG)    &
                 + G(j)*(DG(j,k)-DSUMG(k))
        enddo
       ENDDO
       DEALLOCATE(DSUMG)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       SUMF = 0.0d0
       DO j=NB+1,NBF5
        SUMF = SUMF + FI(j) 
       ENDDO
       IF(SUMF==0.0d0)SUMF=1.0d-6      ! to avoid dividing by zero
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       ALLOCATE(DSUMF(NV))
       do k=1,nv
        DSUMF(k) = 0.0d0
        DO j=NB+1,NBF5
         DSUMF(k) = DSUMF(k) + DF(j,k)
        ENDDO
        DSUMF(k) = DSUMF(k) / (SUMF*SUMF)
       enddo
!-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
       DO j=1,NB
        Hj = 1.0d0-RO(j)
        DO i=1,NB
         Hi = 1.0d0-RO(i)
         Deltaji = G(j)*G(i)
         CJ12(j,i) = CJ12(j,i) - Deltaji
         F1ji = RO(j)*Hi + Deltaji
         F2ji = Hj*RO(i) + Deltaji
         CK12(j,i) = DSQRT( F1ji*F2ji )
         do k=1,nv
          DDeltajik = DG(j,k)*G(i)
          DCJ12r(j,i,k) = DCJ12r(j,i,k) - DDeltajik
          DF1jik =  DR(j,k)*Hi    + DDeltajik
          DF2jik = -DR(j,k)*RO(i) + DDeltajik
          if(DABS(F1ji)>0.0d0)                                          &
           DCK12r(j,i,k) = DCK12r(j,i,k) + 0.5d0*CK12(j,i)*DF1jik/F1ji   
          if(DABS(F2ji)>0.0d0)                                          &
           DCK12r(j,i,k) = DCK12r(j,i,k) + 0.5d0*CK12(j,i)*DF2jik/F2ji   
         enddo                                                           
        ENDDO                                                            
        DO i=NB+1,NBF5                                                   
         Hi = 1.0d0-RO(i)                                                
         Deltaji = FI(j)*FI(i)/SUMF                                      
         CJ12(j,i) = CJ12(j,i) - Deltaji                                 
         F1ji = RO(j)*Hi + Deltaji                                       
         F2ji = Hj*RO(i) + Deltaji                                       
         CK12(j,i) = DSQRT( F1ji*F2ji )                                  
         do k=1,nv                                                       
          DDeltajik = DF(j,k)*FI(i)/SUMF                                 
          DCJ12r(j,i,k) = DCJ12r(j,i,k) - DDeltajik                      
          DF1jik =  DR(j,k)*Hi    + DDeltajik                            
          DF2jik = -DR(j,k)*RO(i) + DDeltajik                            
          if(DABS(F1ji)>0.0d0)                                          &
           DCK12r(j,i,k) = DCK12r(j,i,k) + 0.5d0*CK12(j,i)*DF1jik/F1ji   
          if(DABS(F2ji)>0.0d0)                                          &
           DCK12r(j,i,k) = DCK12r(j,i,k) + 0.5d0*CK12(j,i)*DF2jik/F2ji   
         enddo                                                           
        ENDDO                                                            
       ENDDO                                                             
       DO j=NB+1,NBF5                                                    
        Hj = 1.0d0-RO(j)                                                 
        DO i=1,NB                                                        
         Hi = 1.0d0-RO(i)                                                
         Deltaji = FI(j)*FI(i)/SUMF                                      
         CJ12(j,i) = CJ12(j,i) - Deltaji                                 
         F1ji = RO(j)*Hi + Deltaji                                       
         F2ji = Hj*RO(i) + Deltaji                                       
         CK12(j,i) = DSQRT( F1ji*F2ji )                                  
         do k=1,nv                                                       
          DDeltajik = DF(j,k)*FI(i)/SUMF - FI(j)*FI(i)*DSUMF(k)          
          DCJ12r(j,i,k) = DCJ12r(j,i,k) - DDeltajik                      
          DF1jik= DR(j,k)*Hi    + DDeltajik                              
          DF2jik=-DR(j,k)*RO(i) + DDeltajik                              
          if(DABS(F1ji)>0.0d0)                                          &
           DCK12r(j,i,k) = DCK12r(j,i,k) + 0.5d0*CK12(j,i)*DF1jik/F1ji   
          if(DABS(F2ji)>0.0d0)                                          &
           DCK12r(j,i,k) = DCK12r(j,i,k) + 0.5d0*CK12(j,i)*DF2jik/F2ji
         enddo
        ENDDO
        DO i=NB+1,NBF5
         Hi = 1.0d0-RO(i)
         Deltaji = G(j)*G(i)
         CJ12(j,i) = CJ12(j,i) - Deltaji
         F1ji = RO(j)*Hi + Deltaji
         F2ji = Hj*RO(i) + Deltaji
         CK12(j,i) = - DSQRT( F1ji*F2ji )
         do k=1,nv
          DDeltajik = DG(j,k)*G(i)
          DCJ12r(j,i,k) = DCJ12r(j,i,k) - DDeltajik
          DF1jik =  DR(j,k)*Hi    + DDeltajik
          DF2jik = -DR(j,k)*RO(i) + DDeltajik
          if(DABS(F1ji)>0.0d0)                                          &
           DCK12r(j,i,k) = DCK12r(j,i,k) + 0.5d0*CK12(j,i)*DF1jik/F1ji   
          if(DABS(F2ji)>0.0d0)                                          &
           DCK12r(j,i,k) = DCK12r(j,i,k) + 0.5d0*CK12(j,i)*DF2jik/F2ji
         enddo
        ENDDO
       ENDDO
!-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
       DEALLOCATE(G,DG,FI,DF,DSUMF)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
      ENDIF
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!                Alpha-Alpha Components (N > 2), Daa = Dab
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
      IF(NTWOPAR/=1)CALL AACOMP(NV,CJ12,CK12,DCJ12r,DCK12r)
!-----------------------------------------------------------------------
      RETURN
      END

! CJCKD6 = PNOF6(Nc)
      SUBROUTINE CJCKD6(NV,RO,DR,CJ12,CK12,DCJ12r,DCK12r)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_PNOF/IPNOF,NTWOPAR
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      COMMON/SUMSZ/SUMS,SUMF
!
      DOUBLE PRECISION,DIMENSION(NBF5)::RO
      DOUBLE PRECISION,DIMENSION(NBF5,NV)::DR
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5)::CJ12,CK12
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5,NV)::DCJ12r,DCK12r
      ALLOCATABLE::DSUMS(:),G(:),FI(:),DG(:,:),DF(:,:),DSUMF(:),DSUMG(:)
!-----------------------------------------------------------------------
!          Alpha-Beta Components of CJ12, CK12, DCJ12r and DCK12r
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!                         CJpq = NpNq, CKpq = 0
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CK12 = 0.0d0
      DCK12r = 0.0d0
      DO j=1,NBF5
       DO i=1,NBF5
        CJ12(j,i) = RO(j)*RO(i)
        do k=1,nv        
         DCJ12r(j,i,k) = DR(j,k)*RO(i)
        enddo
       ENDDO
      ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Derivatives of SUMS
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(DSUMS(NV))
      do k=1,nv
       DSUMS(k) = 0.0d0
       DO j=1,NB
        DSUMS(k) = DSUMS(k) - DR(j,k)
       ENDDO
      enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Gq, DGqk
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(G(NBF5),DG(NBF5,NV))
      DO j=1,NB
       Hj = 1.0d0 - RO(j)
       G(j) = Hj*DEXP(-SUMS)
       do k=1,nv
        DG(j,k) = - ( DR(j,k) + Hj*DSUMS(k) ) * DEXP(-SUMS)
       enddo
      ENDDO
      DO j=NB+1,NBF5
       G(j) = RO(j)*DEXP(-SUMS)
       do k=1,nv        
        DG(j,k) = ( DR(j,k) - RO(j)*DSUMS(k) ) * DEXP(-SUMS)
       enddo
      ENDDO
      DEALLOCATE(DSUMS)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     SUMG, DSUMG
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      SUMG = 0.0d0
      DO j=NB+1,NBF5
       SUMG = SUMG + G(j) 
      ENDDO
!- - - - - - - - - - - - - - - - - - - - 
      ALLOCATE(DSUMG(NV))
      do k=1,nv
       DSUMG(k) = 0.0d0
       DO j=NB+1,NBF5
        DSUMG(k) = DSUMG(k) + DG(j,k)
       ENDDO
      enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     FIq = Nq*Hq + Gq*Gq - Gq*S
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(FI(NBF5),DF(NBF5,NV))
      DO j=1,NBF5
       FI(j) = RO(j)*(1.0d0-RO(j)) + G(j)*(G(j)-SUMG)
       do k=1,nv        
        DF(j,k) = (1.0d0-2.0d0*RO(j))*DR(j,k) + DG(j,k)*(G(j)-SUMG)     &
                + G(j)*(DG(j,k)-DSUMG(k))
       enddo
      ENDDO
      DEALLOCATE(DSUMG)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     SUMF, DSUMF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      SUMF = 0.0d0
      DO j=NB+1,NBF5
       SUMF = SUMF + FI(j) 
      ENDDO
      IF(SUMF==0.0d0)SUMF=1.0d-6      ! to avoid dividing by zero
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(DSUMF(NV))
      do k=1,nv
       DSUMF(k) = 0.0d0
       DO j=NB+1,NBF5
        DSUMF(k) = DSUMF(k) + DF(j,k)
       ENDDO
       DSUMF(k) = DSUMF(k) / (SUMF*SUMF)
      enddo
!-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
!                 Including interactions in CJ and CK
!-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
      DO j=1,NB
       Hj = 1.0d0-RO(j)
       DO i=1,NB
        Hi = 1.0d0-RO(i)
        CJ12(j,i) = CJ12(j,i) - G(j)*G(i)
        CK12(j,i) = DSQRT(G(j)*G(i))
        do k=1,nv        
         DCJ12r(j,i,k) = DCJ12r(j,i,k) - DG(j,k)*G(i)
         if(G(j)>0.0d0)DCK12r(j,i,k) = 0.5d0*CK12(j,i)*DG(j,k)/G(j)
        enddo
       ENDDO
       DO i=NB+1,NBF5                    
        Hi = 1.0d0-RO(i)
        Deltaji = FI(j)*FI(i)/SUMF
        CJ12(j,i) = CJ12(j,i) - Deltaji
        F1ji = RO(j)*Hi + Deltaji
        F2ji = Hj*RO(i) + Deltaji
        CK12(j,i) = DSQRT( F1ji*F2ji )
        do k=1,nv        
         DCJ12r(j,i,k) = DCJ12r(j,i,k) - DF(j,k)*FI(i)/SUMF
         DF1jik =  DR(j,k)*Hi    + DF(j,k)*FI(i)/SUMF
         DF2jik = -DR(j,k)*RO(i) + DF(j,k)*FI(i)/SUMF
         if(DABS(F1ji)>0.0d0)                                           &
          DCK12r(j,i,k) = DCK12r(j,i,k) + 0.5d0*CK12(j,i)*DF1jik/F1ji    
         if(DABS(F2ji)>0.0d0)                                           &
          DCK12r(j,i,k) = DCK12r(j,i,k) + 0.5d0*CK12(j,i)*DF2jik/F2ji    
        enddo                                                            
       ENDDO                                                             
      ENDDO                                                              
      DO j=NB+1,NBF5                                                     
       Hj = 1.0d0-RO(j)                                                  
       DO i=1,NB                                                         
        Hi = 1.0d0-RO(i)                                                 
        Deltaji = FI(j)*FI(i)/SUMF                                       
        CJ12(j,i) = CJ12(j,i) - Deltaji                                  
        F1ji = RO(j)*Hi + Deltaji                                        
        F2ji = Hj*RO(i) + Deltaji                                        
        CK12(j,i) = DSQRT( F1ji*F2ji )                                   
        do k=1,nv                                                        
         DCJ12r(j,i,k) = DCJ12r(j,i,k) - DF(j,k)*FI(i)/SUMF             &
                       + FI(j)*FI(i)*DSUMF(k)                            
         DF1jik= DR(j,k)*Hi   +DF(j,k)*FI(i)/SUMF-FI(j)*FI(i)*DSUMF(k)   
         DF2jik=-DR(j,k)*RO(i)+DF(j,k)*FI(i)/SUMF-FI(j)*FI(i)*DSUMF(k)   
         if(DABS(F1ji)>0.0d0)                                           &
          DCK12r(j,i,k) = DCK12r(j,i,k) + 0.5d0*CK12(j,i)*DF1jik/F1ji    
         if(DABS(F2ji)>0.0d0)                                           &
          DCK12r(j,i,k) = DCK12r(j,i,k) + 0.5d0*CK12(j,i)*DF2jik/F2ji
        enddo
       ENDDO
       DO i=NB+1,NBF5
        Hi = 1.0d0-RO(i)
        CJ12(j,i) = CJ12(j,i) - G(j)*G(i)
        CK12(j,i) = - DSQRT(G(j)*G(i))
        do k=1,nv        
         DCJ12r(j,i,k) = DCJ12r(j,i,k) - DG(j,k)*G(i)
         if(G(j)>0.0d0)DCK12r(j,i,k) = 0.5d0*CK12(j,i)*DG(j,k)/G(j)
        enddo
       ENDDO
      ENDDO
!-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
      DEALLOCATE(G,DG,FI,DF,DSUMF)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!             Alpha-Alpha Components (N > 2), Daa = Dab
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
      IF(NTWOPAR/=1)CALL AACOMP(NV,CJ12,CK12,DCJ12r,DCK12r)
!-----------------------------------------------------------------------
      RETURN
      END

! AACOMP
      SUBROUTINE AACOMP(NV,CJ12,CK12,DCJ12r,DCK12r)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5)::CJ12,CK12
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5,NV)::DCJ12r,DCK12r
!-----------------------------------------------------------------------
      DO j=1,NBF5
       DO i=1,NBF5
        CK12(j,i) = CK12(j,i) +  CJ12(j,i)
        CJ12(j,i) = 2.0d0*CJ12(j,i)
        do k=1,nv
         DCK12r(j,i,k) = DCK12r(j,i,k) + DCJ12r(j,i,k)
         DCJ12r(j,i,k) = 2.0d0*DCJ12r(j,i,k) 
        enddo
       ENDDO
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! CJCKD5 = PNOF5(Nc)
      SUBROUTINE CJCKD5(NV,RO,DR,BETA,DB,CJ12,CK12,DCJ12r,DCK12r)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)  
      LOGICAL HighSpin
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_ORBSPACE0/NO1,NDOC,NCO,NCWO,NVIR,NAC,NO0
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Extended (Nc>1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DOUBLE PRECISION,DIMENSION(NBF5)::RO,BETA
      DOUBLE PRECISION,DIMENSION(NBF5,NV)::DR,DB
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5)::CJ12,CK12
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5,NV)::DCJ12r,DCK12r
!-----------------------------------------------------------------------
!                Inter-pair interactions for PNOF5(Nc)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!     CJpq = 2NpNq, CKpq = NpNq [ DCJpqk = 2DNpk*Nq, DCKpq = DNpk*Nq ]
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!           NO1 |     NDNS     |              NVIR         = NBF               
!           NO1 | NDOC + NSOC  |     NCWO*NDOC     + NO0   = NBF               
!                      | NSOC  | NBF5 - NDNS - NO1 | NO0   = NBF
!                      NB      NA                NBF5
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
      DO j=1,NBF5
       DO i=1,NBF5
        CJ12(j,i) = 2.0d0*RO(j)*RO(i)
        CK12(j,i) = RO(j)*RO(i)
        do k=1,nv
         DCJ12r(j,i,k) = 2.0d0*DR(j,k)*RO(i)
         DCK12r(j,i,k) = DR(j,k)*RO(i)
        enddo
       ENDDO
      ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - -              
      if(MSpin==0.and.NSOC>1)then                     
       DO j=NB+1,NA
        DO i=NB+1,NA
         CK12(j,i) = 2.0d0*RO(j)*RO(i)
         do k=1,nv
          DCK12r(j,i,k) = 2.0d0*DR(j,k)*RO(i)
         enddo
        ENDDO      
       ENDDO
      end if
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!                Intra-pair interactions for PNOF5(Nc)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
      DO l=1,NDOC
       ln = NO1+l
       !DO i=NDNS+NCWO*(NDOC-l)+1,NDNS+NCWO*(NDOC-l+1)      !old-sort
       DO i=NDNS+NDOC-l+1,NDNS+NDOC-l+1+(NCWO-1)*NDOC,NDOC  !new-sort
        in = NO1+i
        CJ12(ln,in) = 0.0d0
        CJ12(in,ln) = 0.0d0
        CK12(ln,in) = BETA(ln)*BETA(in)
        CK12(in,ln) = BETA(in)*BETA(ln)
        do k=1,nv
         DCJ12r(ln,in,k) = 0.0d0
         DCJ12r(in,ln,k) = 0.0d0
         DCK12r(ln,in,k) = DB(ln,k)*BETA(in)
         DCK12r(in,ln,k) = DB(in,k)*BETA(ln)
        enddo
        !DO j=NDNS+NCWO*(NDOC-l)+1,NDNS+NCWO*(NDOC-l+1)     !old-sort
        DO j=NDNS+NDOC-l+1,NDNS+NDOC-l+1+(NCWO-1)*NDOC,NDOC !new-sort
         jn = NO1+j
         CJ12(jn,in) = 0.0d0
         CK12(jn,in) = - BETA(jn)*BETA(in)
         do k=1,nv
          DCJ12r(jn,in,k) = 0.0d0
          DCK12r(jn,in,k) = - DB(jn,k)*BETA(in)
         enddo
        ENDDO
       ENDDO
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! CJCKD7 = PNOF7(Nc)(Ista=0) and PNOF7s(Nc)(Ista=1)
      SUBROUTINE CJCKD7(NV,RO,DR,BETA,DB,CJ12,CK12,DCJ12r,DCK12r)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      LOGICAL HighSpin
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_ORBSPACE0/NO1,NDOC,NCO,NCWO,NVIR,NAC,NO0
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
      COMMON/INPNOF_STATIC/Ista
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
!
      DOUBLE PRECISION,DIMENSION(NBF5)::RO,BETA
      DOUBLE PRECISION,DIMENSION(NBF5,NV)::DR,DB
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5)::CJ12,CK12
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5,NV)::DCJ12r,DCK12r
      ALLOCATABLE::FIs(:),DFIs(:,:)
!-----------------------------------------------------------------------
!     FIs
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(FIs(NBF5),DFIs(NBF5,NV))
      FIs = 0.0d0      
      DFIs = 0.0d0
      if(Ista==0)then
!- - - - - - - - - - - - - - - - - - - - - - - - - - -  
!      FIs = (Np*Hp)^1/2
!- - - - - - - - - - - - - - - - - - - - - - - - - - -
       DO j=NO1+1,NBF5
        FIs(j) = DSQRT( RO(j)*(1.0d0-RO(j)) )
        if(FIs(j)>1.0d-20)then
         do k=1,nv        
          DFIs(j,k) = (0.5d0-RO(j))*DR(j,k)/FIs(j)
         enddo
        endif
       ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - -              
      else if(Ista==1)then
!- - - - - - - - - - - - - - - - - - - - - - - - - - -  
!      FIs = (4*Np*Hp)^0.5*(Np*Hp)^0.5 = 2*Np*Hp
!- - - - - - - - - - - - - - - - - - - - - - - - - - -
       DO j=NO1+1,NBF5
        FIs(j) = 2.0d0*RO(j)*(1.0d0-RO(j))
        do k=1,nv        
         DFIs(j,k) = 2.0d0*(1.0d0-2.0d0*RO(j))*DR(j,k)
        enddo
       ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - -       
      end if
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                  - Interpair Electron Correlation -
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!           NO1 |     NDNS     |              NVIR         = NBF               
!           NO1 | NDOC + NSOC  |     NCWO*NDOC     + NO0   = NBF               
!                      | NSOC  | NBF5 - NDNS - NO1 | NO0   = NBF
!                      NB      NA                NBF5
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO j=1,NBF5
       DO i=1,NBF5
        CJ12(j,i) = 2.0d0*RO(j)*RO(i)
        CK12(j,i) = RO(j)*RO(i) + FIs(j)*FIs(i)        
        do k=1,nv
         DCJ12r(j,i,k) = 2.0d0*DR(j,k)*RO(i)
         DCK12r(j,i,k) = DR(j,k)*RO(i) + DFIs(j,k)*FIs(i)   
        enddo
       ENDDO
      ENDDO
      DEALLOCATE(FIs,DFIs)      
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(MSpin==0.and.NSOC>1)then
       DO j=NB+1,NA
        DO i=NB+1,NA
         CK12(j,i) = 2.0d0*RO(j)*RO(i)
         do k=1,nv
          DCK12r(j,i,k) = 2.0d0*DR(j,k)*RO(i)
         enddo
        ENDDO      
       ENDDO
      end if
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                  - Intrapair Electron Correlation -
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO l=1,NDOC
       ln = NO1+l
       !DO i=NDNS+NCWO*(NDOC-l)+1,NDNS+NCWO*(NDOC-l+1)      !old-sort
       DO i=NDNS+NDOC-l+1,NDNS+NDOC-l+1+(NCWO-1)*NDOC,NDOC  !new-sort
        in = NO1+i
        CJ12(ln,in) = 0.0d0
        CJ12(in,ln) = 0.0d0
        CK12(ln,in) = BETA(ln)*BETA(in)
        CK12(in,ln) = BETA(in)*BETA(ln)
        do k=1,nv
         DCJ12r(ln,in,k) = 0.0d0
         DCJ12r(in,ln,k) = 0.0d0
         DCK12r(ln,in,k) = DB(ln,k)*BETA(in)
         DCK12r(in,ln,k) = DB(in,k)*BETA(ln)
        enddo
        !DO j=NDNS+NCWO*(NDOC-l)+1,NDNS+NCWO*(NDOC-l+1)     !old-sort
        DO j=NDNS+NDOC-l+1,NDNS+NDOC-l+1+(NCWO-1)*NDOC,NDOC !new-sort
         jn = NO1+j
         CJ12(jn,in) = 0.0d0
         CK12(jn,in) = - BETA(jn)*BETA(in)
         do k=1,nv
          DCJ12r(jn,in,k) = 0.0d0
          DCK12r(jn,in,k) = - DB(jn,k)*BETA(in)
         enddo
        ENDDO
       ENDDO
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! CJCKD8 = GNOF(Imod=0) GNOFm(Imod=1)
      SUBROUTINE CJCKD8(NV,RO,DR,BETA,DB,CJ12,CK12,DCJ12r,DCK12r)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL HighSpin
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_ORBSPACE0/NO1,NDOC,NCO,NCWO,NVIR,NAC,NO0
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
      COMMON/INPNOF_MOD/Imod
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
!
      DOUBLE PRECISION,DIMENSION(NBF5)::RO,BETA
      DOUBLE PRECISION,DIMENSION(NBF5,NV)::DR,DB
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5)::CJ12,CK12
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5,NV)::DCJ12r,DCK12r
      ALLOCATABLE :: ROd(:),DROd(:,:),Rd(:),DRd(:,:),FIs(:),DFIs(:,:)
!-----------------------------------------------------------------------
      ALLOCATE (ROd(NBF5),DROd(NBF5,NV),Rd(NBF5),DRd(NBF5,NV))
      ALLOCATE (FIs(NBF5),DFIs(NBF5,NV))
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!          Dynamic and Static Occupation Numbers (ROd,Rd,FIs)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(Imod==0 .OR. Imod==1)then                   ! GNOFm
       Hcutd = 0.02d0*DSQRT(2.0d0)
      end if
      ROd  = 0.0d0
      DROd = 0.0d0
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     ROd (NO1+1:NBF5)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO i=1,NDOC
       in = NO1+i
       HOLEin = 1.0d0-RO(in)
       COC = HOLEin/Hcutd
       ARG  = - COC**2
       DARG = - 2*COC/Hcutd
       Fin = DEXP(ARG)                       ! Hd/Hole
!      ROd(NO1+1:NB)
       ROd(in) = RO(in) * Fin
       do k=1,nv
        DROd(in,k) = Fin * DR(in,k) * ( 1.0d0 - RO(in)*DARG )
       enddo
!      ROd(NA+1:NBF5)
       IF(NCWO>1)THEN                        ! extended PNOF
        do iw=1,ncwo
         ! above Fermi level
         !im = na+ncwo*(ndoc-i)+iw                       !old-sort
         im = no1+ndoc+(na-nb)+(ndoc-i+1)+ndoc*(iw-1)    !new-sort
         ROd(im) = RO(im) * Fin              ! ROd = RO*Hd/Hole
         do k=1,nv
          DROd(im,k) = Fin * ( DR(im,k) - RO(im)*DARG*DR(in,k) )
         enddo
        enddo
       ELSE                                  ! perfect-pairing
        icf = na+ndoc-i+1
        ROd(icf) = RO(icf) * Fin             ! ROd = RO*Hd/Hole
        do k=1,nv
         DROd(icf,k) = Fin * (DR(icf,k) - RO(icf)*DARG*DR(in,k) )
        enddo
       ENDIF
      ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Rd(1:NBF5) = SQRT(ROd)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Rd   = 0.0d0
      DRd  = 0.0d0
      DO j=1,NBF5
       Rd(j) = DSQRT(ROd(j))
       do k=1,nv
        if(Rd(j)>0.0d0)DRd(j,k) = 0.5d0 * DROd(j,k) / Rd(j)
       enddo
      ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     FIs
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      FIs = 0.0d0
      DFIs = 0.0d0
!- - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(Imod==0 .OR. Imod==1)then
!- - - - - - - - - - - - - - - - - - - - - - - - - - -
!      FIs = (Np*Hp)^1/2
!- - - - - - - - - - - - - - - - - - - - - - - - - - -
       DO j=NO1+1,NBF5
        FIs(j) = DSQRT( RO(j)*(1.0d0-RO(j)) )
        if(FIs(j)>0.0d0)then
         do k=1,nv
          DFIs(j,k) = (0.5d0-RO(j))*DR(j,k)/FIs(j)
         enddo
        endif
       ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - -
      end if
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!                  - Interpair Electron Correlation -
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(Imod==0) THEN
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       CJpq = 2NpNq , CKpq = NpNq
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        DO j=1,NBF5
         DO i=1,NBF5
          CJ12(j,i) = 2.0d0*RO(j)*RO(i)
          CK12(j,i) = RO(j)*RO(i)
          do k=1,nv
           DCJ12r(j,i,k) = 2.0d0*DR(j,k)*RO(i)
           DCK12r(j,i,k) = DR(j,k)*RO(i)
          enddo
         ENDDO
        ENDDO    
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       CKpq = NpNq - PIs_pq
!       Static component: PIs_pq = - FIs(p)*FIs(q)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        DO j=NO1+1,NB
         DO i=NA+1,NBF5     
          CK12(j,i) = CK12(j,i) + FIs(j)*FIs(i) 
          do k=1,nv
           DCK12r(j,i,k) = DCK12r(j,i,k) + DFIs(j,k)*FIs(i)
          enddo
         ENDDO
        ENDDO
        DO j=NA+1,NBF5           
         DO i=NO1+1,NB
          CK12(j,i) = CK12(j,i) + FIs(j)*FIs(i) 
          do k=1,nv
           DCK12r(j,i,k) = DCK12r(j,i,k) + DFIs(j,k)*FIs(i)
          enddo
         ENDDO
         DO i=NA+1,NBF5
          CK12(j,i) = CK12(j,i) + FIs(j)*FIs(i) 
          do k=1,nv
           DCK12r(j,i,k) = DCK12r(j,i,k) + DFIs(j,k)*FIs(i)
          enddo
         ENDDO
        ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if(MSpin==0.and.NSOC>0)then
         DO j=NO1+1,NB           
          DO i=NB+1,NA                                              
           CK12(j,i) = CK12(j,i) + 0.5d0*FIs(j)*0.5d0
           do k=1,nv                                                
            DCK12r(j,i,k) = DCK12r(j,i,k) + 0.5d0*DFIs(j,k)*0.5d0        
           enddo                                                    
          ENDDO                                                     
         ENDDO      
         DO j=NB+1,NA
          DO i=NO1+1,NB     
           CK12(j,i) = CK12(j,i) + 0.5d0*0.5d0*FIs(i) 
          ENDDO
         ENDDO
  !      
         DO j=NB+1,NA
          DO i=NA+1,NBF5     
           CK12(j,i) = CK12(j,i) + 0.5d0*FIs(i) 
          ENDDO
         ENDDO
         DO j=NA+1,NBF5           
          DO i=NB+1,NA
           CK12(j,i) = CK12(j,i) + FIs(j)*0.5d0
           do k=1,nv
            DCK12r(j,i,k) = DCK12r(j,i,k) + DFIs(j,k)*0.5d0
           enddo
          ENDDO
         ENDDO       
        end if
  
        if(MSpin==0.and.NSOC>1)then
         DO j=NB+1,NA
          DO i=NB+1,NA
           CK12(j,i) = 0.5d0              ! 2.0d0*RO(j)*RO(i)
           do k=1,nv
            DCK12r(j,i,k) = 0.0d0         ! 2.0d0*DR(j,k)*RO(i)
           enddo
          ENDDO      
         ENDDO
        end if            
      ELSE IF(Imod==1) THEN
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       CJpq = 2NpNq , CKpq = NpNq + FIs(p)*FIs(q)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        DO j=1,NBF5
         DO i=1,NBF5
          CJ12(j,i) = 2.0d0*RO(j)*RO(i)
          CK12(j,i) = RO(j)*RO(i) + FIs(j)*FIs(i)
          do k=1,nv
           DCJ12r(j,i,k) = 2.0d0*DR(j,k)*RO(i)
           DCK12r(j,i,k) = DR(j,k)*RO(i) + DFIs(j,k)*FIs(i)
          enddo
         ENDDO
        ENDDO
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     CKpq = CKpq - PId_pq
!     Dynamic component: PId_pq = - Rd(p)*Rd(q) + ROd(p)*ROd(q)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO j=NO1+1,NB
       DO i=NA+1,NBF5
        CK12(j,i) = CK12(j,i) + Rd(j)*Rd(i) - ROd(j)*ROd(i)
        do k=1,nv
         DCK12r(j,i,k) = DCK12r(j,i,k)+DRd(j,k)*Rd(i)-DROd(j,k)*ROd(i)
        enddo
       ENDDO
      ENDDO
      DO j=NA+1,NBF5
       DO i=NO1+1,NB
        CK12(j,i) = CK12(j,i) + Rd(j)*Rd(i) - ROd(j)*ROd(i)
        do k=1,nv
         DCK12r(j,i,k) = DCK12r(j,i,k)+DRd(j,k)*Rd(i)-DROd(j,k)*ROd(i)
        enddo
       ENDDO
       DO i=NA+1,NBF5
        CK12(j,i) = CK12(j,i) - Rd(j)*Rd(i) - ROd(j)*ROd(i)
        do k=1,nv
         DCK12r(j,i,k) = DCK12r(j,i,k)-DRd(j,k)*Rd(i)-DROd(j,k)*ROd(i)
        enddo
       ENDDO
      ENDDO
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!                  - Intrapair Electron Correlation -
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
      DO l=1,NDOC
       ln = NO1+l
       !DO i=NDNS+NCWO*(NDOC-l)+1,NDNS+NCWO*(NDOC-l+1)     !old-sort
       DO i=NDNS+NDOC-l+1,NDNS+NDOC-l+1+(NCWO-1)*NDOC,NDOC !new-sort
        in = NO1+i
        CJ12(ln,in) = 0.0d0
        CJ12(in,ln) = 0.0d0
        CK12(ln,in) = BETA(ln)*BETA(in)
        CK12(in,ln) = BETA(in)*BETA(ln)
        do k=1,nv
         DCJ12r(ln,in,k) = 0.0d0
         DCJ12r(in,ln,k) = 0.0d0
         DCK12r(ln,in,k) = DB(ln,k)*BETA(in)
         DCK12r(in,ln,k) = DB(in,k)*BETA(ln)
        enddo
        !DO j=NDNS+NCWO*(NDOC-l)+1,NDNS+NCWO*(NDOC-l+1)     !old-sort
        DO j=NDNS+NDOC-l+1,NDNS+NDOC-l+1+(NCWO-1)*NDOC,NDOC !new-sort
         jn = NO1+j
         CJ12(jn,in) = 0.0d0
         CK12(jn,in) = - BETA(jn)*BETA(in)
         do k=1,nv
          DCJ12r(jn,in,k) = 0.0d0
          DCK12r(jn,in,k) = - DB(jn,k)*BETA(in)
         enddo
        ENDDO
       ENDDO
      ENDDO
!-----------------------------------------------------------------------
      DEALLOCATE(ROd,DROd,Rd,DRd,FIs,DFIs)
      RETURN
      END

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
!                                                                      !
!                           Output Routines                            !
!                                                                      !
!  OUTINITIALSr: Output of properties with the input, stop if IEINI=1  !
!  GENOUTPUTr: Print Occupations, 1-Energies, Sum Occ                  !
!  PRINTOCr: Initial and Intermediate Outputs                          !
!  WRITEGCFr: Write ONs, NOs, and diag-elements of F on NFILE unix     !
!  FINALOUTPUTr: Final Output of the occupation optimization results   !
!  DIAGELAG: Diagonalization of Lagrange Multiplier Matrix if DIAGLAG  !
!  PRINTV: Print vectors V(N) without heading on NFILE unix            !
!  PRINTVERO: Print vectors V(N) with heading En-Occ on NFILE unix     !
!                                                                      !
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !

! OUTINITIALSr
      SUBROUTINE OUTINITIALSr(NV,GAMMA,COEF,RO,CJ12,CK12,               &
                             DR,DCJ12r,DCK12r,QD,HCORE,QJ,QK,           &
                             DIPN,ADIPx,ADIPy,ADIPz,DIPx,DIPy,DIPz,     &
                             QUADN,AQUADxx,AQUADyy,AQUADzz,AQUADxy,     &
                             AQUADxz,AQUADyz,QUADxx,QUADyy,QUADzz,      &
                             QUADxy,QUADxz,QUADyz,OCTUN,AOCTxxx,        &
                             AOCTyyy,AOCTzzz,AOCTxxy,AOCTxxz,AOCTxyy,   &
                             AOCTyyz,AOCTxzz,AOCTyzz,AOCTxyz,OCTxxx,    &
                             OCTyyy,OCTzzz,OCTxxy,OCTxxz,OCTxyy,OCTyyz, &
                             OCTxzz,OCTyzz,OCTxyz,ATMNAME,ZNUC,LIMLOW,  &
                             LIMSUP,OVERLAP,IT,ITTOTAL,IPRINTOPT,IRUNTYP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      LOGICAL EFIELDL,HighSpin,SCALING
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPSCALING/SCALING,NZEROS,NZEROSm,NZEROSr,ITZITER
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
      COMMON/INPNOF_PRINT/NPRINT,IWRITEC,IMULPOP,IAIMPAC,IFCHK
      COMMON/INPNOF_GENERALINF/ICOEF,ISOFTMAX,IORBOPT,MAXIT,MAXIT21
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      COMMON/INP_EFIELDL/EX,EY,EZ,EFIELDL
      COMMON/ELPROP/IEMOM
      COMMON/EHFEN/EHF,EN
      COMMON/ENERGIAS/EELEC,EELEC_OLD,DIF_EELEC,EELEC_MIN
!
      INTEGER :: NV,IT,ITTOTAL,IPRINTOPT,IRUNTYP
      CHARACTER*4,DIMENSION(NATOMS)::ATMNAME
      INTEGER,DIMENSION(NATOMS)::LIMLOW,LIMSUP
      DOUBLE PRECISION,DIMENSION(3)::DIPN
      DOUBLE PRECISION,DIMENSION(6)::QUADN
      DOUBLE PRECISION,DIMENSION(10)::OCTUN
      DOUBLE PRECISION,DIMENSION(NATOMS)::ZNUC
      DOUBLE PRECISION,DIMENSION(NBF5)::GAMMA,RO,HCORE,DIPx,DIPy,DIPz
      DOUBLE PRECISION,DIMENSION(NBF5)::QUADxx,QUADyy,QUADzz,QUADxy
      DOUBLE PRECISION,DIMENSION(NBF5)::QUADxz,QUADyz      
      DOUBLE PRECISION,DIMENSION(NBF5)::OCTxxx,OCTyyy,OCTzzz,OCTxxy
      DOUBLE PRECISION,DIMENSION(NBF5)::OCTxxz,OCTxyy,OCTyyz,OCTxzz
      DOUBLE PRECISION,DIMENSION(NBF5)::OCTyzz,OCTxyz
      DOUBLE PRECISION,DIMENSION(NBFT5)::QJ,QK
      DOUBLE PRECISION,DIMENSION(NBF5,NV)::DR
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::COEF,OVERLAP
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::ADIPx,ADIPy,ADIPz
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::AQUADxx,AQUADyy,AQUADzz
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::AQUADxy,AQUADxz,AQUADyz
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::AOCTxxx,AOCTyyy,AOCTzzz
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::AOCTxxy,AOCTxxz,AOCTxyy
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::AOCTyyz,AOCTxzz
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::AOCTyzz,AOCTxyz
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5)::CJ12,CK12
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5,NV)::DCJ12r,DCK12r
      DOUBLE PRECISION,DIMENSION(NBF,NBF,NBF)::QD
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::GRAD,EAHF,E
!
      DOUBLE PRECISION,PARAMETER::DFAC=2.54174D0    ! Debye
      DOUBLE PRECISION,PARAMETER::QFAC=1.345044D0   ! Buckinham
      DOUBLE PRECISION,PARAMETER::OFAC=7.117668D-01 ! X10**34 ESU-CM**3
!-----------------------------------------------------------------------
!     Evaluate the Electronic Energy
!-----------------------------------------------------------------------
      ALLOCATE(GRAD(NV))
!
      MODE = 0
      if(MSpin==0)then
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!      Singlet State (S=0,Ms=0) and Multiplet States (S>0,Ms=0)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
       CALL OCUPENERGYrc(MODE,GAMMA,RO,CJ12,CK12,DR,DCJ12r,DCK12r,      &
                         QD,HCORE,QJ,QK,DIPN,ADIPx,ADIPy,ADIPz,         &
                         DIPx,DIPy,DIPz,EELEC,GRAD,NV)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --    
      else if(MSpin>0)then
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!      High-Spin Multiplet State (S>0,Ms=S)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
       CALL OCUPENERGYro(MODE,GAMMA,RO,CJ12,CK12,DR,DCJ12r,DCK12r,      &
                         QD,HCORE,QJ,QK,DIPN,ADIPx,ADIPy,ADIPz,         &
                         DIPx,DIPy,DIPz,EELEC,GRAD,NV)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
      end if
!      
      DEALLOCATE(GRAD)
!-----------------------------------------------------------------------
      IF(NV>0)THEN
!-----------------------------------------------------------------------
!      Initial Values
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       IF(NZEROS<=NZEROSm .and. IRUNTYP/=3 )WRITE(6,1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      One-particle Energies
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       ALLOCATE(EAHF(NBF),E(NBF))
       CALL ENENEWr(RO,HCORE,QJ,QK,CJ12,CK12,DIPx,DIPy,DIPz,EAHF,E)
       DEALLOCATE(EAHF)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Write Coeficient Matrix (Natural Orbitals)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       IF(IWRITEC==1.AND.IPRINTOPT==1)THEN
        WRITE(6,11)
        CALL PRINTVERO(6,COEF,E,RO,NBF,NBF5)                
       ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Mulliken Population Analysis
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       IF(IMULPOP==1.AND.IPRINTOPT==1)THEN
        WRITE(6,12)
        if(MSpin==0)then
         CALL MULLIKENrc(ATMNAME,ZNUC,LIMLOW,LIMSUP,OVERLAP,RO,QD)
        else if(MSpin>0)then
!        CALL MULLIKENro        to do!
        end if
       ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      General Output (Occupations)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       IF(IPRINTOPT==1)CALL GENOUTPUTr(RO,E)
       DEALLOCATE(E)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Print Electric Field
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       IF(EFIELDL.AND.IPRINTOPT==1)WRITE(6,13)EX,EY,EZ
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Calculate Electrostatic Moments and print on the main output (6)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       IF(IEMOM>=1.AND.IPRINTOPT==1)THEN
        CALL DIPMOMr(DIPN,ADIPx,ADIPy,ADIPz,DIPx,DIPy,DIPz,             &
                     QD,RO,DMXe,DMYe,DMZe,DMX,DMY,DMZ,DM)                
        WRITE(6,14)DM*DFAC,DM,DMX,DMY,DMZ                                
       END IF                                                            
       IF(IEMOM>=2.AND.IPRINTOPT==1)THEN                                 
        CALL QUADMOMr(QUADN,AQUADxx,AQUADyy,AQUADzz,AQUADxy,AQUADxz,    &
                      AQUADyz,QUADxx,QUADyy,QUADzz,QUADxy,QUADxz,QUADyz,&
                      QD,RO,QMXXe,QMYYe,QMZZe,QMXYe,QMXZe,QMYZe,        &
                      QTxx,QTyy,QTzz,QTxy,QTxz,QTyz)                     
        WRITE(6,15)QTxx*QFAC,QTyy*QFAC,QTzz*QFAC,                       &
                   QTxy*QFAC,QTxz*QFAC,QTyz*QFAC,                       &
                   QTxx,QTyy,QTzz,QTxy,QTxz,QTyz                         
       ENDIF                                                             
       IF(IEMOM==3.AND.IPRINTOPT==1)THEN                                 
        CALL OCTMOMr(OCTUN,AOCTxxx,AOCTyyy,AOCTzzz,AOCTxxy,AOCTxxz,     &
                     AOCTxyy,AOCTyyz,AOCTxzz,AOCTyzz,AOCTxyz,OCTxxx,    &
                     OCTyyy,OCTzzz,OCTxxy,OCTxxz,OCTxyy,OCTyyz,OCTxzz,  &
                     OCTyzz,OCTxyz,QD,RO,OMXXXe,OMYYYe,OMZZZe,OMXXYe,   &
                     OMXXZe,OMXYYe,OMYYZe,OMXZZe,OMYZZe,OMXYZe,         &
                     OTXXX,OTYYY,OTZZZ,OTXXY,OTXXZ,                     &
                     OTXYY,OTYYZ,OTXZZ,OTYZZ,OTXYZ)                      
        WRITE(6,16)OTXXX*OFAC,OTYYY*OFAC,OTZZZ*OFAC,OTXXY*OFAC,         &
                   OTXXZ*OFAC,OTXYY*OFAC,OTYYZ*OFAC,OTXZZ*OFAC,         &
                   OTYZZ*OFAC,OTXYZ*OFAC,OTXXX,OTYYY,OTZZZ,             &
                   OTXXY,OTXXZ,OTXYY,OTYYZ,OTXZZ,OTYZZ,OTXYZ
       ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Energy Output
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       IF(ICOEF==2)RETURN
       IF(IPRINTOPT==0)STOP
       IF(NZEROS<=NZEROSm)THEN
        WRITE(6,100)
        WRITE(6,101)EHF+EN
        WRITE(6,102)EELEC+EN
        WRITE(6,103)EELEC-EHF
       ELSE
        WRITE(6,104)
        WRITE(6,101)EHF+EN
        WRITE(6,102)EELEC+EN
        WRITE(6,103)EELEC-EHF
        WRITE(6,105)IT,ITTOTAL
        WRITE(6,200)
       ENDIF
       STOP       
!-----------------------------------------------------------------------
      ENDIF
!-----------------------------------------------------------------------
    1 FORMAT(/3X,'INITIAL VALUES'                                       &
           ,/2X,'=================')
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   11 FORMAT(/,                                                         &
       18X,'------------------------------------------',/,              &
       18X,' NATURAL ORBITALS IN ATOMIC ORBITAL BASIS ',/,              &
       18X,'------------------------------------------')                 
   12 FORMAT(                                                           &
         2X,'------------------------------',/,                         &
         2X,' Mulliken Population Analysis ',/,                         &
         2X,'------------------------------')                            
   13 FORMAT(/,6X,'ELECTRIC FIELD (',D8.1,',',D8.1,',',D8.1,')')         
   14 FORMAT(/,2X,'---------------',                                    &
              /2X,' Dipole Moment',                                     &
              /2X,'---------------',                                    &
            //,3X,F9.4,' Debye',' [',F9.4,                              &
               2X,'(',F9.4,',',F9.4,',',F9.4,')',' ]')                   
   15 FORMAT(/,2X,'-------------------',                                &
              /2X,' Quadrupole Moment',                                 &
              /2X,'-------------------',                                &
             //6X,'QXX',6X,'QYY',6X,'QZZ',6X,'QXY',6X,'QXZ',6X,'QYZ',   &
             //2X,6F9.4,2X,'(Buckingham)',//1X,'[',6F9.4,1X,']')         
   16 FORMAT(/,2X,'-----------------',                                  &
              /2X,' Octupole Moment',                                   &
              /2X,'-----------------',                                  &
             //6X,'OXXX',5X,'OYYY',5X,'OZZZ',5X,'OXXY',5X,'OXXZ',       &
               5X,'OXYY',5X,'OYYZ',5X,'OXZZ',5X,'OYZZ',5X,'OXYZ',       &
             //2X,10F9.4,2X,'(X10**34 ESU-CM**3)',                      &
             //1X,'[',10F9.4,1X,']')
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  100 FORMAT(//2X,'INITIAL ENERGY (AU)',                                &
             /1X,'=====================')                                
  101 FORMAT(/,8X,' HF Total Energy =',F14.6)                            
  102 FORMAT(  8X,'NOF Total Energy =',F14.6)                            
  103 FORMAT(  6X,'Correlation Energy =',F14.6,/)                        
  104 FORMAT(//3X,'ENERGY (AU)',                                        &
             /2X,'=============')                                        
  105 FORMAT(/2X,'**************************************************',  &
             /2X,'*         No. EXTERNAL ITER =',I6,'              *',  &
             /2X,'*         No. of TOTAL ITER =',I6,'              *',  &
             /2X,'**************************************************')
!----------------------------------------------------------------------
  200 FORMAT(//2X,'**************************************************', &
              /2X,'*                                                *', &
              /2x,'*       SINGLE-POINT PNOF CALCULATION            *', &
              /2X,'*                                                *', &
              /2X,'*       BE CAREFUL NZEROS > NZEROSm!             *', &
              /2X,'*                                                *', &
              /2x,'*  FINAL RESULTS   FINAL RESULTS  FINAL RESULTS  *', &
              /2X,'*                                                *', &
              /2X,'**************************************************', &
              /)
!----------------------------------------------------------------------
      END

! GENOUTPUTr
      SUBROUTINE GENOUTPUTr(RO,E)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL HighSpin
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
!
      DOUBLE PRECISION,DIMENSION(NBF5)::RO,E
!-----------------------------------------------------------------------
!     Print the Occupations
!     SUMRO: Sum of the Occupations
!-----------------------------------------------------------------------
      WRITE(6,1)
      SUMRO = 0.0d0
      
      DO J=1,NB
       WRITE(6,2)J,2.0*RO(J),E(J)
       SUMRO = SUMRO + 2.0d0*RO(J)       
      ENDDO
      
      IF(NSOC>0)THEN      
       if(.not.HighSpin)then      
        DO J=NB+1,NA
         WRITE(6,2)J,2.0*RO(J),E(J)
         SUMRO = SUMRO + 2.0d0*RO(J)         
        ENDDO
       else if(HighSpin)then      
        DO J=NB+1,NA
         WRITE(6,2)J,RO(J),E(J)
         SUMRO = SUMRO + RO(J)         
        ENDDO
       end if
      END IF 
      
      DO J=NA+1,NBF5
       WRITE(6,2)J,2.0*RO(J),E(J)
       SUMRO = SUMRO + 2.0d0*RO(J)       
      ENDDO
      
      WRITE(6,3)NE,SUMRO
!-----------------------------------------------------------------------
    1 FORMAT(/,3X,'OM',5X,'Occupation',6X,'Elag Diag',/)
    2 FORMAT(2X,I3,2F15.7)
    3 FORMAT(/3X,'RO Sum (',I4,') =',F8.2)
!-----------------------------------------------------------------------     
      RETURN
      END

! PRINTOCr
      SUBROUTINE PRINTOCr(E,COEF,ATMNAME,ZNUC,LIMLOW,LIMSUP,            &
                          OVERLAP,RO,CJ12,CK12,QD,QK,                   &
                          DIPN,ADIPx,ADIPy,ADIPz,DIPx,DIPy,DIPz,        &
                          QUADN,AQUADxx,AQUADyy,AQUADzz,AQUADxy,        &
                          AQUADxz,AQUADyz,QUADxx,QUADyy,QUADzz,         &
                          QUADxy,QUADxz,QUADyz,OCTUN,                   &
                          AOCTxxx,AOCTyyy,AOCTzzz,AOCTxxy,AOCTxxz,      &
                          AOCTxyy,AOCTyyz,AOCTxzz,AOCTyzz,AOCTxyz,      &
                          OCTxxx,OCTyyy,OCTzzz,OCTxxy,OCTxxz,           &
                          OCTxyy,OCTyyz,OCTxzz,OCTyzz,OCTxyz,           &
                          IZCORE,CX0,CY0,CZ0,KSTART,KNG,                &
                          KKMIN,KKMAX,KATOM,KTYPE,KLOC,IMIN,IMAX,       &
                          ISH,ITYP,EX1,C1,C2,CS,CP,CD,CF,CG,CH,CI,      &
                          IFIRSTCALL,DIPS,IRUNTYP)                          
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)  
      LOGICAL EFIELDL,RESTART,APSG,OIMP2,MBPT,HighSpin
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INP_EFIELDL/EX,EY,EZ,EFIELDL
      COMMON/INPNOF_ORBSPACE0/NO1,NDOC,NCO,NCWO,NVIR,NAC,NO0
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
      COMMON/INPNOF_RSTRT/RESTART,INPUTGAMMA,INPUTC,INPUTFMIUG,INPUTCXYZ
      COMMON/INPNOF_APSG/APSG,NTHAPSG,THAPSG
      COMMON/INPNOF_OIMP2/OIMP2,MBPT
      COMMON/INPNOF_MOLDEN/MOLDEN
      COMMON/INPNOF_ARDM/THRESHDM,NOUTRDM,NSQT,NTHRESHDM      
      COMMON/INPNOF_CJK/NOUTCJK,NTHRESHCJK,THRESHCJK
      COMMON/INPNOF_PRINT/NPRINT,IWRITEC,IMULPOP,IAIMPAC,IFCHK
      COMMON/INPNOF_PNOF/IPNOF,NTWOPAR
      COMMON/INPNOF_GENERALINF/ICOEF,ISOFTMAX,IORBOPT,MAXIT,MAXIT21
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      COMMON/EHFEN/EHF,EN
      COMMON/ENERGIAS/EELEC,EELEC_OLD,DIF_EELEC,EELEC_MIN
      COMMON/CorrNonDynamic/ECnd,ECndl,ECndHF
      COMMON/ELPROP/IEMOM      
!
      CHARACTER*4,DIMENSION(NATOMS)::ATMNAME
      INTEGER,DIMENSION(NATOMS)::LIMLOW,LIMSUP,IZCORE,IMIN,IMAX
      INTEGER,DIMENSION(NSHELL)::KSTART,KNG,KKMIN,KKMAX,KATOM,KTYPE,KLOC
      INTEGER,DIMENSION(NPRIMI)::ISH,ITYP
      DOUBLE PRECISION,DIMENSION(3)::DIPN
      DOUBLE PRECISION,DIMENSION(6)::QUADN
      DOUBLE PRECISION,DIMENSION(10)::OCTUN
      DOUBLE PRECISION,DIMENSION(NATOMS)::ZNUC,CX0,CY0,CZ0
      DOUBLE PRECISION,DIMENSION(NPRIMI)::EX1,C1,C2,CS,CP,CD,CF,CG,CH,CI
      DOUBLE PRECISION,DIMENSION(NBF)::E
      DOUBLE PRECISION,DIMENSION(NBF5)::RO,DIPx,DIPy,DIPz
      DOUBLE PRECISION,DIMENSION(NBF5)::QUADxx,QUADyy,QUADzz
      DOUBLE PRECISION,DIMENSION(NBF5)::QUADxy,QUADxz,QUADyz
      DOUBLE PRECISION,DIMENSION(NBF5)::OCTxxx,OCTyyy,OCTzzz,OCTxxy
      DOUBLE PRECISION,DIMENSION(NBF5)::OCTxxz,OCTxyy,OCTyyz,OCTxzz
      DOUBLE PRECISION,DIMENSION(NBF5)::OCTyzz,OCTxyz
      DOUBLE PRECISION,DIMENSION(NBFT5)::QK
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5)::CJ12,CK12      
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::COEF,OVERLAP
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::ADIPx,ADIPy,ADIPz
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::AQUADxx,AQUADyy,AQUADzz
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::AQUADxy,AQUADxz,AQUADyz
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::AOCTxxx,AOCTyyy,AOCTzzz
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::AOCTxxy,AOCTxxz,AOCTxyy
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::AOCTyyz,AOCTxzz
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::AOCTyzz,AOCTxyz
      DOUBLE PRECISION,DIMENSION(NBF,NBF,NBF)::QD
      DOUBLE PRECISION,DIMENSION(3):: DIPS
!
      DOUBLE PRECISION,PARAMETER::DFAC=2.54174D0    ! Debye
      DOUBLE PRECISION,PARAMETER::QFAC=1.345044D0   ! Buckinham
      DOUBLE PRECISION,PARAMETER::OFAC=7.117668D-01 ! X10**34 ESU-CM**3
!-----------------------------------------------------------------------
      KLOC(1) = KLOC(1)
      IFIRSTCALL = IFIRSTCALL
!-----------------------------------------------------------------------
!     Skip Printing if NPRINT/=2
!-----------------------------------------------------------------------
      IF(ICOEF==0.or.NPRINT==2)THEN
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Write Header on output file
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       IF(ICOEF/=2)WRITE(6,1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      General Output (Occupations and N-representability conditions)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       CALL GENOUTPUTr(RO,E)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Electric Field and Electrostatic Moments
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       IF(EFIELDL)WRITE(6,7)EX,EY,EZ
       IF(IEMOM>=1)THEN
        CALL DIPMOMr(DIPN,ADIPx,ADIPy,ADIPz,DIPx,DIPy,DIPz,             &
                     QD,RO,DMXe,DMYe,DMZe,DMX,DMY,DMZ,DM)                
        DIPS(1) = DMX                                                    
        DIPS(2) = DMY                                                    
        DIPS(3) = DMZ                                                    
        WRITE(6,8)DM*DFAC,DM,DMX,DMY,DMZ                                 
       END IF                                                            
       IF(IEMOM>=2)THEN                                                  
        CALL QUADMOMr(QUADN,AQUADxx,AQUADyy,AQUADzz,AQUADxy,AQUADxz,    &
                      AQUADyz,QUADxx,QUADyy,QUADzz,QUADxy,QUADxz,QUADyz,&
                      QD,RO,QMXXe,QMYYe,QMZZe,QMXYe,QMXZe,QMYZe,        &
                      QTxx,QTyy,QTzz,QTxy,QTxz,QTyz)                     
        WRITE(6,9)QTxx*QFAC,QTyy*QFAC,QTzz*QFAC,                        &
                  QTxy*QFAC,QTxz*QFAC,QTyz*QFAC,                        &
                  QTxx,QTyy,QTzz,QTxy,QTxz,QTyz                          
       ENDIF                                                             
       IF(IEMOM==3)THEN                                                  
        CALL OCTMOMr(OCTUN,AOCTxxx,AOCTyyy,AOCTzzz,AOCTxxy,AOCTxxz,     &
                     AOCTxyy,AOCTyyz,AOCTxzz,AOCTyzz,AOCTxyz,OCTxxx,    &
                     OCTyyy,OCTzzz,OCTxxy,OCTxxz,OCTxyy,OCTyyz,OCTxzz,  &
                     OCTyzz,OCTxyz,QD,RO,OMXXXe,OMYYYe,OMZZZe,OMXXYe,   &
                     OMXXZe,OMXYYe,OMYYZe,OMXZZe,OMYZZe,OMXYZe,         &
                     OTXXX,OTYYY,OTZZZ,OTXXY,OTXXZ,                     &
                     OTXYY,OTYYZ,OTXZZ,OTYZZ,OTXYZ)                      
        WRITE(6,10)OTXXX*OFAC,OTYYY*OFAC,OTZZZ*OFAC,OTXXY*OFAC,         &
                   OTXXZ*OFAC,OTXYY*OFAC,OTYYZ*OFAC,OTXZZ*OFAC,         &
                   OTYZZ*OFAC,OTXYZ*OFAC,OTXXX,OTYYY,OTZZZ,             &
                   OTXXY,OTXXZ,OTXYY,OTYYZ,OTXZZ,OTYZZ,OTXYZ
       ENDIF       
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Coeficient matrix and one-particle occupations
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       IF(IWRITEC==1)THEN
        IF(ICOEF==0)THEN
         WRITE(6,3)
        ELSE
         WRITE(6,4)
        ENDIF
        CALL PRINTVERO(6,COEF,E,RO,NBF,NBF5)        
       ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Mulliken Population Analysis
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       IF(IMULPOP==1.and.MSpin==0)THEN
        WRITE(6,5)
        CALL MULLIKENrc(ATMNAME,ZNUC,LIMLOW,LIMSUP,OVERLAP,RO,QD)
       ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Energy Output
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       WRITE(6,100)
       IF(INPUTC==0.and.EHF<0.0)THEN
        WRITE(6,101)EHF+EN
        WRITE(6,102)EELEC+EN
        WRITE(6,103)EELEC-EHF
       ELSE
        WRITE(6,102)EELEC+EN
       ENDIF
!-----------------------------------------------------------------------
      END IF
!
      IF(ICOEF==0)THEN
!-----------------------------------------------------------------------
!      Bader's analysis: Write information into wavefunction file (WFN) 
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       IF(IAIMPAC==1.and.MSpin==0)THEN
        CALL AIMMEMrc(COEF,ZNUC,IZCORE,CX0,CY0,CZ0,KSTART,KNG,KKMIN,    &
                      KKMAX,KATOM,KTYPE,RO,E,EX1,CS,CP,CD,CF,CG,CH,CI)
       END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Write information into formatted checkpoint file (FCHK) [Unit=19]
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       IF(IFCHK==1.and.MSpin==0)THEN
        CALL FCHKrc(COEF,ZNUC,IZCORE,CX0,CY0,CZ0,KNG,KATOM,KTYPE,       &
                    RO,E,EX1,C1,C2,IMIN,IMAX,ISH,DIPS,IRUNTYP)
       END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Write information into a file in Molden Format (MLD) [Unit=17]
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       IF(MOLDEN==1.and.MSpin==0)THEN
        CALL MOLDENrc(ATMNAME,IZCORE,ZNUC,CX0,CY0,CZ0,IMIN,IMAX,        &
                      ISH,ITYP,EX1,C1,C2,RO,E,COEF)        
       END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      APSG File for the APSG generating wavefunction of PNOF5(Nc)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       IF(APSG.and.NSOC==0.and.IPNOF==5)THEN
        REWIND(9)
!       Coefficient matrix
        WRITE(9,'(A5)')' $VEC'
        CALL PUNCHVEC(COEF,NBF)
        WRITE(9,'(A5)')' $END'
!       APSG Expansion Coefficients of the Generating PNOF5-wavefunction
        WRITE(9,'(A6)')' $APSG'
        CALL PUNCHAPSG(NO1,NCWO,NA,NBF5,RO,SUMA,THAPSG)
        WRITE(9,'(A5)')' $END'
        WRITE(9,'(1X,A16,F15.10)')'EXP. COEF. SUM =',SUMA
       ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      PRINT CJ12 and CK12 (MOLECULAR RDMs) IN CJK FILE
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       IF(NOUTCJK==1.and.MSpin==0)CALL OUTPUTCJKrc(RO,CJ12,CK12)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      PRINT ATOMIC RDMs IN .1dm and .2dm FILES
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       IF(1<=NOUTRDM .and. NOUTRDM<=3 .and. MSpin==0)THEN
        CALL OUTPUTRDMrc(OVERLAP,RO,QD,CJ12,CK12)
       ENDIF
!-----------------------------------------------------------------------       
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Non-Dynamic Correction if OIMP2 or MBPT
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(OIMP2.or.MBPT)CALL ECorrNonDyn(RO,QK,ECndl)
!-----------------------------------------------------------------------
    1 FORMAT(//2X,'RESULTS OF THE OCCUPATION OPTIMIZATION'              &
            ,/1X,'========================================')             
!   2 FORMAT(/,3X,'Chem Pot =',F12.6,4X,F15.6)                           
    3 FORMAT(/,                                                         &
       18X,'------------------------------------------',/,              &
       18X,' NATURAL ORBITALS IN ATOMIC ORBITAL BASIS ',/,              &
       18X,'------------------------------------------')                 
    4 FORMAT(/,                                                         &
       2X,'--------------------',/,                                     &
       2X,' COEFFICIENT MATRIX',/,                                      &
       2X,'--------------------')                                        
    5 FORMAT(                                                           &
       /,2X,'------------------------------',/,                         &
         2X,' Mulliken Population Analysis ',/,                         &
         2X,'------------------------------')                            
!   6 FORMAT( 2X,'-------------',                                       &
!            /2X,' Occupations',                                        &
!            /2X,'-------------')                                        
    7 FORMAT(/,6X,'Electric Field (',D8.1,',',D8.1,',',D8.1,')')         
    8 FORMAT(/,2X,'---------------',                                    &
              /2X,' Dipole Moment',                                     &
              /2X,'---------------',                                    &
            //,3X,F9.4,' Debye',' [',F9.4,                              &
               2X,'(',F15.10,',',F15.10,',',F15.10,')',' ]')            
    9 FORMAT(/,2X,'-------------------',                                &
              /2X,' Quadrupole Moment',                                 &
              /2X,'-------------------',                                &
             //6X,'QXX',6X,'QYY',6X,'QZZ',6X,'QXY',6X,'QXZ',6X,'QYZ',   &
             //2X,6F9.4,2X,'(Buckingham)',//1X,'[',6F9.4,1X,']')         
   10 FORMAT(/,2X,'-----------------',                                  &
              /2X,' Octupole Moment',                                   &
              /2X,'-----------------',                                  &
             //6X,'OXXX',5X,'OYYY',5X,'OZZZ',5X,'OXXY',5X,'OXXZ',       &
               5X,'OXYY',5X,'OYYZ',5X,'OXZZ',5X,'OYZZ',5X,'OXYZ',       &
             //2X,10F9.4,2X,'(X10**34 ESU-CM**3)',                      &
             //1X,'[',10F9.4,1X,']')
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  100 FORMAT(/,                                                         &
       2X,'----------------------------------------',/,                 &
       2X,' Energies after Occupation Optimization ',/,                 &
       2X,'----------------------------------------')
  101 FORMAT(/,8X,' HF Total Energy =',F14.6)
  102 FORMAT(  8X,'NOF Total Energy =',F14.6)
  103 FORMAT(  6X,'Correlation Energy =',F14.6,/)
!-----------------------------------------------------------------------
      RETURN
      END

! WRITEGCFr
      SUBROUTINE WRITEGCFr(NFILE,RO,SUMS,C,E,F,NSQ,NBF,NBF5,IT,EELEC,EN,&
                           NO1,NDOC,NSOC,NCWO,NAC,NO0,ZNUC,CX0,CY0,CZ0, &
                           NAT,IWRTRO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DIMENSION RO(NBF5),C(NSQ),E(NBF),F(NBF)
      DOUBLE PRECISION,DIMENSION(NAT):: ZNUC,CX0,CY0,CZ0
      DOUBLE PRECISION,PARAMETER::BOHR = 0.52917724924D+00
!-----------------------------------------------------------------------
      REWIND(NFILE)
      if(IWRTRO==0)then  ! Don't write RO if RUNTYP=DYN (IRUNTYP==5)
       DO I=1,NBF5
        READ(NFILE,'(I6,F30.16)')II,RRII
       ENDDO
      else if(IWRTRO==1)then
       DO I=1,NBF5
        WRITE(NFILE,'(I6,F30.16)')I,RO(I)
       ENDDO
      endif
!      
      DO I=NBF5+1,NBF
       WRITE(NFILE,'(I6,F30.16)')I,0.0d0
      ENDDO
      WRITE(NFILE,'(F30.16)')SUMS
      DO J = 1,NBF
       DO I = 1,NBF
        WRITE(NFILE,'(I6,F30.16)')I,C(I+(J-1)*NBF)
       ENDDO
      ENDDO
      DO I = 1,NBF
       WRITE(NFILE,'(I6,F30.16)')I,E(I)
      ENDDO
      DO I = 1,NBF
       WRITE(NFILE,'(I6,F30.16)')I,F(I)
      ENDDO
      WRITE(NFILE,'(I6,F30.16)')IT,EELEC+EN
      WRITE(NFILE,'(6I6)')NO1,NDOC,NSOC,NCWO,NAC,NO0
      DO I = 1,NAT
       WRITE(NFILE,'(I6,3F30.16)')INT(ZNUC(I)),                         &
                                  BOHR*CX0(I),BOHR*CY0(I),BOHR*CZ0(I)
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! FINALOUTPUTr
      SUBROUTINE FINALOUTPUTr(E,RO,CJ12,CK12,QD,HCORE,QJ,QK,            &
                              ELAG,COEF,ATMNAME,ZNUC,LIMLOW,LIMSUP,     &
                              OVERLAP,DIPN,ADIPx,ADIPy,ADIPz,DIPx,DIPy, &
                              DIPz,QUADN,AQUADxx,AQUADyy,AQUADzz,       &
                              AQUADxy,AQUADxz,AQUADyz,QUADxx,QUADyy,    &
                              QUADzz,QUADxy,QUADxz,QUADyz,OCTUN,        &
                              AOCTxxx,AOCTyyy,AOCTzzz,AOCTxxy,AOCTxxz,  &
                              AOCTxyy,AOCTyyz,AOCTxzz,AOCTyzz,AOCTxyz,  &
                              OCTxxx,OCTyyy,OCTzzz,OCTxxy,OCTxxz,       &
                              OCTxyy,OCTyyz,OCTxzz,OCTyzz,OCTxyz,       &
                              IZCORE,CX0,CY0,CZ0,KSTART,KNG,KKMIN,      &
                              KKMAX,KATOM,KTYPE,KLOC,IMIN,IMAX,ISH,ITYP,&
                              EX1,C1,C2,CS,CP,CD,CF,CG,CH,CI,ELAGN,     &
                              COEFN,RON,AHCORE,IJKL,XIJKL,IFIRSTCALL,   &
                              DIPS,IPRINTOPT,IRUNTYP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      LOGICAL EFIELDL,CONVGDELAG,PRINTLAG,DIAGLAG,APSG,CHKORTHO,ORTHO
      LOGICAL OIMP2,MBPT,HighSpin
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INP_EFIELDL/EX,EY,EZ,EFIELDL
      COMMON/CONVERGENCE/DUMEL,PCONV,CONVGDELAG
      COMMON/INPNOF_ORBSPACE0/NO1,NDOC,NCO,NCWO,NVIR,NAC,NO0
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
      COMMON/INPNOF_LAGRANGE/PRINTLAG,DIAGLAG
      COMMON/INPNOF_APSG/APSG,NTHAPSG,THAPSG 
      COMMON/INPNOF_ORTHOGONALITY/CHKORTHO,ORTHO
      COMMON/INPNOF_OIMP2/OIMP2,MBPT
      COMMON/INPNOF_EKT/IEKT
      COMMON/INPNOF_MOLDEN/MOLDEN
      COMMON/INPNOF_ARDM/THRESHDM,NOUTRDM,NSQT,NTHRESHDM      
      COMMON/INPNOF_CJK/NOUTCJK,NTHRESHCJK,THRESHCJK
      COMMON/INPNOF_PRINT/NPRINT,IWRITEC,IMULPOP,IAIMPAC,IFCHK
      COMMON/INPNOF_PNOF/IPNOF,NTWOPAR
      COMMON/INPNOF_GENERALINF/ICOEF,ISOFTMAX,IORBOPT,MAXIT,MAXIT21
      COMMON/INPFILE_NIJKL/NINTMX,NIJKL,NINTCR,NSTORE
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      COMMON/ELPROP/IEMOM      
      COMMON/EHFEN/EHF,EN
      COMMON/ENERGIAS/EELEC,EELEC_OLD,DIF_EELEC,EELEC_MIN
      COMMON/CorrNonDynamic/ECnd,ECndl,ECndHF
!
      INTEGER :: IFIRSTCALL,IPRINTOPT,IRUNTYP
      CHARACTER*4,DIMENSION(NATOMS)::ATMNAME
      INTEGER,DIMENSION(NATOMS)::LIMLOW,LIMSUP,IZCORE,IMIN,IMAX
      INTEGER,DIMENSION(NSHELL)::KSTART,KNG,KKMIN,KKMAX,KATOM,KTYPE,KLOC
      INTEGER,DIMENSION(NPRIMI)::ISH,ITYP
      INTEGER,DIMENSION(NIJKL)::IJKL
      DOUBLE PRECISION,DIMENSION(3)::DIPN
      DOUBLE PRECISION,DIMENSION(6)::QUADN
      DOUBLE PRECISION,DIMENSION(10)::OCTUN
      DOUBLE PRECISION,DIMENSION(NATOMS)::ZNUC,CX0,CY0,CZ0
      DOUBLE PRECISION,DIMENSION(NPRIMI)::EX1,C1,C2,CS,CP,CD,CF,CG,CH,CI
      DOUBLE PRECISION,DIMENSION(NBF5)::HCORE,RO,DIPx,DIPy,DIPz
      DOUBLE PRECISION,DIMENSION(NBF5)::QUADxx,QUADyy,QUADzz
      DOUBLE PRECISION,DIMENSION(NBF5)::QUADxy,QUADxz,QUADyz
      DOUBLE PRECISION,DIMENSION(NBF5)::OCTxxx,OCTyyy,OCTzzz,OCTxxy
      DOUBLE PRECISION,DIMENSION(NBF5)::OCTxxz,OCTxyy,OCTyyz,OCTxzz
      DOUBLE PRECISION,DIMENSION(NBF5)::OCTyzz,OCTxyz
      DOUBLE PRECISION,DIMENSION(NBFT5)::QJ,QK
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5)::CJ12,CK12
      DOUBLE PRECISION,DIMENSION(NBF)::E,ELAGN,RON
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::ELAG,COEF,OVERLAP
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::ADIPx,ADIPy,ADIPz
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::AQUADxx,AQUADyy,AQUADzz
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::AQUADxy,AQUADxz,AQUADyz
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::AOCTxxx,AOCTyyy,AOCTzzz
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::AOCTxxy,AOCTxxz,AOCTxyy
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::AOCTyyz,AOCTxzz
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::AOCTyzz,AOCTxyz
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::COEFN,AHCORE
      DOUBLE PRECISION,DIMENSION(NBF,NBF,NBF)::QD
      DOUBLE PRECISION,DIMENSION(NIJKL)::XIJKL
      DOUBLE PRECISION,DIMENSION(3):: DIPS
!      
      DOUBLE PRECISION,PARAMETER::DFAC=2.54174D0    ! Debye
      DOUBLE PRECISION,PARAMETER::QFAC=1.345044D0   ! Buckinham
      DOUBLE PRECISION,PARAMETER::OFAC=7.117668D-01 ! X10**34 ESU-CM**3
!-----------------------------------------------------------------------
      KLOC(1) = KLOC(1)
      IFIRSTCALL = IFIRSTCALL
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Check the orthonormality
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(CHKORTHO)CALL CHECKORTHO(COEF,OVERLAP,IVIOORTHO,IPRINTOPT)
!-----------------------------------------------------------------------
      IF(IPRINTOPT==1)THEN
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Write Header on output file
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       IF(ICOEF/=2)WRITE(6,1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      General Output
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       CALL GENOUTPUTr(RO,E)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Electric Field and Electrostatic Moments
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       IF(EFIELDL)WRITE(6,7)EX,EY,EZ
       IF(IEMOM>=1)THEN
        CALL DIPMOMr(DIPN,ADIPx,ADIPy,ADIPz,DIPx,DIPy,DIPz,             &
                     QD,RO,DMXe,DMYe,DMZe,DMX,DMY,DMZ,DM)                
        DIPS(1) = DMX                                                    
        DIPS(2) = DMY                                                    
        DIPS(3) = DMZ                                                    
        WRITE(6,8)DM*DFAC,DM,DMX,DMY,DMZ                                 
       END IF                                                            
       IF(IEMOM>=2)THEN                                                  
        CALL QUADMOMr(QUADN,AQUADxx,AQUADyy,AQUADzz,AQUADxy,AQUADxz,    &
                      AQUADyz,QUADxx,QUADyy,QUADzz,QUADxy,QUADxz,QUADyz,&
                      QD,RO,QMXXe,QMYYe,QMZZe,QMXYe,QMXZe,QMYZe,        &
                      QTxx,QTyy,QTzz,QTxy,QTxz,QTyz)                     
        WRITE(6,9)QTxx*QFAC,QTyy*QFAC,QTzz*QFAC,                        &
                  QTxy*QFAC,QTxz*QFAC,QTyz*QFAC,                        &
                  QTxx,QTyy,QTzz,QTxy,QTxz,QTyz                          
       END IF                                                            
       IF(IEMOM==3)THEN                                                  
        CALL OCTMOMr(OCTUN,AOCTxxx,AOCTyyy,AOCTzzz,AOCTxxy,AOCTxxz,     &
                     AOCTxyy,AOCTyyz,AOCTxzz,AOCTyzz,AOCTxyz,OCTxxx,    &
                     OCTyyy,OCTzzz,OCTxxy,OCTxxz,OCTxyy,OCTyyz,OCTxzz,  &
                     OCTyzz,OCTxyz,QD,RO,OMXXXe,OMYYYe,OMZZZe,OMXXYe,   &
                     OMXXZe,OMXYYe,OMYYZe,OMXZZe,OMYZZe,OMXYZe,         &
                     OTXXX,OTYYY,OTZZZ,OTXXY,OTXXZ,                     &
                     OTXYY,OTYYZ,OTXZZ,OTYZZ,OTXYZ)                      
        WRITE(6,10)OTXXX*OFAC,OTYYY*OFAC,OTZZZ*OFAC,OTXXY*OFAC,         &
                   OTXXZ*OFAC,OTXYY*OFAC,OTYYZ*OFAC,OTXZZ*OFAC,         &
                   OTYZZ*OFAC,OTXYZ*OFAC,OTXXX,OTYYY,OTZZZ,             &
                   OTXXY,OTXXZ,OTXYY,OTYYZ,OTXZZ,OTYZZ,OTXYZ
       END IF
!-----------------------------------------------------------------------
!      Skip Printing if NPRINT=0, otherwise print EKT,ELAG,COEF,MulPop
!-----------------------------------------------------------------------
       IF(NPRINT>0)THEN
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       Coefficient matrix (NOs) and one-particle occupations
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        IF(IWRITEC==1)THEN
         WRITE(6,4)
         CALL PRINTVERO(6,COEF,E,RO,NBF,NBF5)        
        ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       Extended Koopmans Theorem (EKT)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        IF(IEKT==1.and.MSpin==0)THEN
         CALL EXTKOOPMANSrc(ELAG,COEF,OVERLAP,AHCORE,IJKL,XIJKL,RO)
        END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       Quasi-particle Energies
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        IQE=0
        IF(IQE==1.and.MSpin==0)CALL QUASIENEr(RO,HCORE,QJ,QK,CJ12,CK12)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       Print Lagrange Multipliers (ELAG Matrix) on main output (6)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        IF(PRINTLAG)THEN
         WRITE(6,3)
         CALL PRINTV(6,ELAG,NBF)
        ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       Diagonalization of Lagrange Multipliers (ELAG)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        IF(DIAGLAG)THEN
         CALL DIAGELAG(ELAG,COEF,RO,ELAGN,COEFN,RON)       
         RON = RON/2.0d0         
        ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       Mulliken Population Analysis 
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        IF(IMULPOP==1.and.MSpin==0)THEN
         WRITE(6,5)
         CALL MULLIKENrc(ATMNAME,ZNUC,LIMLOW,LIMSUP,OVERLAP,RO,QD)
        ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       Calculate the chemical potential (CHEMP) for PNOF5,7(Nc)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ICHEMPOT=0
        IF(ICHEMPOT==1.and.MSpin==0.and.(IPNOF==5.or.IPNOF==7))         &
         CALL CHEMPOTrc(HCORE,QJ,QK,RO,DIPx,DIPy,DIPz)       
       END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Energy Output
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       WRITE(6,100)
       IF(EHF<0.0)THEN
        WRITE(6,101)EHF+EN
        WRITE(6,102)EELEC+EN,DUMEL,ABS(EELEC_OLD-EELEC)
        WRITE(6,103)EELEC-EHF
       ELSE
        WRITE(6,102)EELEC+EN,DUMEL,ABS(EELEC_OLD-EELEC)
       ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Bader's analysis:Write information into a wavefunction file (WFN) 
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       IF(IAIMPAC==1.and.MSpin==0)THEN
        if(NPRINT>0.and.DIAGLAG)then   ! Canonical Orbitals
         CALL AIMMEMrc(COEFN,ZNUC,IZCORE,CX0,CY0,CZ0,KSTART,KNG,KKMIN,  &
                  KKMAX,KATOM,KTYPE,RON,ELAGN,EX1,CS,CP,CD,CF,CG,CH,CI)
        else
         CALL AIMMEMrc(COEF,ZNUC,IZCORE,CX0,CY0,CZ0,KSTART,KNG,KKMIN,   &
                       KKMAX,KATOM,KTYPE,RO,E,EX1,CS,CP,CD,CF,CG,CH,CI)
        endif                       
       END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Write information into formatted checkpoint file (FCHK) [Unit=19]
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       IF(IFCHK==1.and.MSpin==0)THEN
        if(NPRINT>0.and.DIAGLAG)then   ! Canonical Orbitals
         CALL FCHKrc(COEFN,ZNUC,IZCORE,CX0,CY0,CZ0,KNG,KATOM,KTYPE,     &
                     RON,ELAGN,EX1,C1,C2,IMIN,IMAX,ISH,DIPS,IRUNTYP)
        else
         CALL FCHKrc(COEF,ZNUC,IZCORE,CX0,CY0,CZ0,KNG,KATOM,KTYPE,       &
                     RO,E,EX1,C1,C2,IMIN,IMAX,ISH,DIPS,IRUNTYP)
        endif
       END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Write information into a file in Molden Format (MLD) [Unit=17]
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       IF(MOLDEN==1.and.MSpin==0)THEN
        if(NPRINT>0.and.DIAGLAG)then   ! Canonical Orbitals
         CALL MOLDENrc(ATMNAME,IZCORE,ZNUC,CX0,CY0,CZ0,IMIN,IMAX,       &
                       ISH,ITYP,EX1,C1,C2,RON,ELAGN,COEFN)               
        else                           ! Natural Orbitals                
         CALL MOLDENrc(ATMNAME,IZCORE,ZNUC,CX0,CY0,CZ0,IMIN,IMAX,       &
                       ISH,ITYP,EX1,C1,C2,RO,E,COEF)        
        endif
       END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      APSG File for the APSG generating wavefunction of PNOF5(Nc)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       IF(APSG .and. NSOC==0 .and. IPNOF==5)THEN
        REWIND(9)
!       Coefficient matrix
        WRITE(9,'(A5)')' $VEC'
        CALL PUNCHVEC(COEF,NBF)
        WRITE(9,'(A5)')' $END'
!       APSG Expansion Coefficients of the Generating PNOF5-wavefunction
        WRITE(9,'(A6)')' $APSG'
        CALL PUNCHAPSG(NO1,NCWO,NA,NBF5,RO,SUMA,THAPSG)
        WRITE(9,'(A5)')' $END'
        WRITE(9,'(1X,A16,F15.10)')'EXP. COEF. SUM =',SUMA
       ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      PRINT CJ12 and CK12 (MOLECULAR RDMs) IN CJK FILE
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       IF(NOUTCJK==1.and.MSpin==0)THEN
        CALL OUTPUTCJKrc(RO,CJ12,CK12)
       ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      PRINT ATOMIC RDMs IN .1dm and .2dm FILES
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       IF(1<=NOUTRDM .and. NOUTRDM<=3 .and. MSpin==0)THEN
        CALL OUTPUTRDMrc(OVERLAP,RO,QD,CJ12,CK12)
       ENDIF
!-----------------------------------------------------------------------       
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Non-Dynamic Correction if OIMP2 or MBPT
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(OIMP2.or.MBPT)CALL ECorrNonDyn(RO,QK,ECndl) 
!-----------------------------------------------------------------------
    1 FORMAT(//2X,'RESULTS OF THE OCCUPATION OPTIMIZATION'              &
            ,/1X,'========================================')
!   2 FORMAT(/,3X,'Chem Pot =',F12.6,4X,F15.6)
    3 FORMAT(/,                                                         &
       18X,'--------------------------------',/,                        &
       18X,' MATRIX OF LAGRANGE MULTIPLIERS ',/,                        &
       18X,'--------------------------------')                           
    4 FORMAT(/,                                                         &
       18X,'------------------------------------------',/,              &
       18X,' NATURAL ORBITALS IN ATOMIC ORBITAL BASIS ',/,              &
       18X,'------------------------------------------')                 
    5 FORMAT(/,                                                         &
       /,2X,'------------------------------',/,                         &
         2X,' Mulliken Population Analysis ',/,                         &
         2X,'------------------------------')
!   6 FORMAT(/2X,' Occupations ',                                       &
!             /2X,'-------------')
    7 FORMAT(/,6X,'Electric Field (',D8.1,',',D8.1,',',D8.1,')')
    8 FORMAT(/,2X,'---------------',                                    &
              /2X,' Dipole Moment',                                     &
              /2X,'---------------',                                    &
            //,3X,F9.4,' Debye',' [',F9.4,                              &
               2X,'(',F9.4,',',F9.4,',',F9.4,')',' ]')                   
    9 FORMAT(/,2X,'-------------------',                                &
              /2X,' Quadrupole Moment',                                 &
              /2X,'-------------------',                                &
             //6X,'QXX',6X,'QYY',6X,'QZZ',6X,'QXY',6X,'QXZ',6X,'QYZ',   &
             //2X,6F9.4,2X,'(Buckingham)',//1X,'[',6F9.4,1X,']')         
   10 FORMAT(/,2X,'-----------------',                                  &
              /2X,' Octupole Moment',                                   &
              /2X,'-----------------',                                  &
             //6X,'OXXX',5X,'OYYY',5X,'OZZZ',5X,'OXXY',5X,'OXXZ',       &
               5X,'OXYY',5X,'OYYZ',5X,'OXZZ',5X,'OYZZ',5X,'OXYZ',       &
             //2X,10F9.4,2X,'(X10**34 ESU-CM**3)',                      &
             //1X,'[',10F9.4,1X,']')
!-----------------------------------------------------------------------
  100 FORMAT(/,                                                         &
       1X,'----------------',/,                                         &
       1X,' Final Energies  ',/,                                        &
       1X,'----------------',/)                                          
  101 FORMAT(8X,' HF Total Energy =',F20.10)                             
  102 FORMAT(2X,'Final NOF Total Energy =',                             &
                   F20.10,' (',ES7.0,1X,',',ES7.0,' )')
  103 FORMAT(6X,'Correlation Energy =',F20.10,/)
!-----------------------------------------------------------------------
      RETURN
      END

! DIAGELAG 
      SUBROUTINE DIAGELAG(ELAG,COEF,RO,ELAGN,COEFN,RON)
!=======================================================================
!     DIAGONALIZATION OF LAGRANGE MULTIPLIERS (ELAG)
!=======================================================================
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
!
      DOUBLE PRECISION,DIMENSION(NBF5)::RO
      DOUBLE PRECISION,DIMENSION(NBF)::ELAGN,RON
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::ELAG,COEF,COEFN
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::TEMP
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::AUX,W,DENMAT
!-----------------------------------------------------------------------
!     INTERMEDIATE MATRICES
!-----------------------------------------------------------------------
      ALLOCATE (AUX(NBF,NBF),W(NBF,NBF),TEMP(NBF))
!-----------------------------------------------------------------------
!     DIAGONALIZATION OF THE LAGRANGE MULTIPLIERS (ELAG)
!-----------------------------------------------------------------------
!     ELAG -> SQUARE MATRIX (AUX)
!-----------------------------------------------------------------------
      DO I=1,NBF
       DO J=1,I
        AUX(I,J)=ELAG(I,J)
        AUX(J,I)=AUX(I,J)
       ENDDO
      ENDDO
!-----------------------------------------------------------------------
!     DIAGONALIZE SQUARE MATRIX (AUX) FOR REAL SYMMETRIC CASE
!     NOTE: ONLY LOWER TRIANGLE IS USED + THIS IS DESTROYED !!!
!     W - EIGENVECTORS, ELAGN - EIGENVALUES 
!-----------------------------------------------------------------------
      CALL DIAG(NBF,AUX,W,ELAGN,TEMP)
!-----------------------------------------------------------------------
!     New Density Matrix (D=Wt*RO*W)
!-----------------------------------------------------------------------
      ALLOCATE(DENMAT(NBF,NBF))
      DO IP=1,NBF
       DO IQ=1,NBF
        DENMAT(IP,IQ)=0.0d0
        do i=1,nbf5
         DENMAT(IP,IQ)=DENMAT(IP,IQ)+W(i,IP)*RO(i)*W(i,IQ)
        enddo
       ENDDO
      ENDDO
!-----------------------------------------------------------------------
!     WRITE ONE-PARTICLE ENERGIES and NEW AVERAGE OCCUPATIONS (OUTPUT)
!-----------------------------------------------------------------------
      WRITE(6,1)
      DO I=1,NBF
       RON(I)=2.0d0*DENMAT(I,I)
       IF(RON(I)>1.0d-6)THEN
        WRITE(6,2)I,ELAGN(I),ELAGN(I)*27.21138386,RON(I)
       ENDIF
      ENDDO
      CALL DMATMAX(DENMAT,NBF,MAXI,MAXJ,DUM)
      WRITE(6,3)DUM,MAXI,MAXJ
!-----------------------------------------------------------------------
!     Coefficients of Canonical Orbitals (COEFN=COEF*W)
!-----------------------------------------------------------------------
      COEFN = MATMUL(COEF,W)
      WRITE(6,4)
      CALL PRINTVERO(6,COEFN,ELAGN,RON,NBF,NBF5)
!-----------------------------------------------------------------------
    1 FORMAT(/2X,26('-'),/3X,'Canonical Representation',/2X,26('-'),/,  &
             /20X,'One-Particle Energies',11X,'1RDM Diag',              &
             //19X,'(aU)',14X,'(eV)')                                    
    2 FORMAT(2X,I4,4X,F15.6,4X,F15.6,9X,F8.6)                            
    3 FORMAT(/,15X,'Maximum 1RDM off-diagonal element:',F12.6,          &
               1X,'(',I3,',',I3,')')                                     
    4 FORMAT(/,                                                         &
       18X,'---------------------------------',/,                       &
       18X,' Eigenvectors of Lagrange Matrix ',/,                       &
       18X,'---------------------------------')
!-----------------------------------------------------------------------
      DEALLOCATE (AUX,W,TEMP,DENMAT)
      RETURN
      END

! PRINTV
      SUBROUTINE PRINTV(NFILE,V,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DOUBLE PRECISION,DIMENSION(N,N)::V
!-----------------------------------------------------------------------
      MAX=0
  10  MIN=MAX+1
      MAX=MAX+5
      IF(MAX>N)MAX=N
      WRITE(NFILE,1)(J,J=MIN,MAX)
      WRITE(NFILE,*)
      DO I=1,N
       WRITE(NFILE,2)I,(V(I,J),J=MIN,MAX)
      ENDDO
      IF(MAX<N) GOTO 10
      WRITE(NFILE,3)
!-----------------------------------------------------------------------
    1 FORMAT(/,7X,5(4X,I4,3X))
    2 FORMAT(I5,2X,5F11.6)
    3 FORMAT(/)
!-----------------------------------------------------------------------
      RETURN
      END

! PRINTVERO
      SUBROUTINE PRINTVERO(NFILE,V,E,RO,N,NL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      CHARACTER*2 LABELAT
      CHARACTER*4 BFNAM1
      CHARACTER*6 BFNAM2
      DOUBLE PRECISION,DIMENSION(N)::E,RO
      DOUBLE PRECISION,DIMENSION(N,N)::V
!-----------------------------------------------------------------------
      MAX=0
  10  MIN=MAX+1
      MAX=MAX+5
      IF(MAX>NL) MAX=NL
      WRITE(NFILE,1)(J,J=MIN,MAX)
      WRITE(NFILE,*)      
      WRITE(NFILE,2)(RO(J),J=MIN,MAX)
      WRITE(NFILE,2)(E(J),J=MIN,MAX)      
      REWIND(4)
      READ(4,'(I5)')N0
      DO I=1,N
       IF(I<=35)THEN
        READ(4,3)LABELAT,IAT,BFNAM1
        WRITE(NFILE,4) I,LABELAT,IAT,BFNAM1,(V(I,J),J=MIN,MAX)
       ELSE
        READ(4,5)LABELAT,BFNAM2
        WRITE(NFILE,6) I,LABELAT,BFNAM2,(V(I,J),J=MIN,MAX)
       ENDIF
      ENDDO
      IF(MAX<NL) GOTO 10
      WRITE(NFILE,7)
!-----------------------------------------------------------------------
    1 FORMAT(/,15X,5(4X,I4,3X))
    2 FORMAT(15X,5F11.4)
    3 FORMAT(A2,I2,A4)
    4 FORMAT(I5,2X,A2,I2,A4,5F11.6)
    5 FORMAT(A2,A6)
    6 FORMAT(I5,2X,A2,A6,5F11.6)
    7 FORMAT(/)
!-----------------------------------------------------------------------
      RETURN
      END

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
!                                                                      !
!             Multipurpose-Product Functions used in DoNOF             !
!                                                                      !
!   PRODCWQWj: PRODCWQW(j) = CW12j*QWj excluding i=j                   !
!   PRODCWQWjk: PRODCWQWjk = CW12jk*QWj excluding i=j                  !
!   PRODROQWj: PRODROQWj = RO*QWj                                      !
!   PRODWCWij: PRODWCWij = W(ij)*CW12                                  !
!   PRODWCWijq: PRODWCWijq = W(ij)*CW12(q) excluding q'=q              !
!                                                                      !
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !

! PRODCWQWj
      FUNCTION PRODCWQWj(J,CW12,QW)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5)::CW12
      DOUBLE PRECISION,DIMENSION(NBFT5)::QW
!-----------------------------------------------------------------------
!     PRODCWQWj = CW12j*QWj. Note: Term with I=J is not included
!-----------------------------------------------------------------------
      PRODCWQWj = 0.0d0
      DO I=1,J-1
       PRODCWQWj = PRODCWQWj + CW12(J,I)*QW(I+J*(J-1)/2)
      ENDDO
      DO I=J+1,NBF5
       PRODCWQWj = PRODCWQWj + CW12(J,I)*QW(J+I*(I-1)/2)
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! PRODCWQWj1
      FUNCTION PRODCWQWj1(J,CW12,QW)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5)::CW12
      DOUBLE PRECISION,DIMENSION(NBFT5)::QW
!-----------------------------------------------------------------------
!     PRODCWQWj1 = CW12j*QWj, j<NB (NSOC is excluded from the Sum)
!     Note: Term with I=J is not included
!-----------------------------------------------------------------------
      PRODCWQWj1 = 0.0d0
      DO I=1,J-1
       PRODCWQWj1 = PRODCWQWj1 + CW12(J,I)*QW(I+J*(J-1)/2)
      ENDDO
      DO I=J+1,NB
       PRODCWQWj1 = PRODCWQWj1 + CW12(J,I)*QW(J+I*(I-1)/2)
      ENDDO
      DO I=NA+1,NBF5
       PRODCWQWj1 = PRODCWQWj1 + CW12(J,I)*QW(J+I*(I-1)/2)
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END      

! PRODCWQWj2
      FUNCTION PRODCWQWj2(J,CW12,QW)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5)::CW12
      DOUBLE PRECISION,DIMENSION(NBFT5)::QW
!-----------------------------------------------------------------------
!     PRODCWQWj2 = CW12j*QWj, j>NA (NSOC is excluded from the Sum)
!     Note: Term with I=J is not included
!-----------------------------------------------------------------------
      PRODCWQWj2 = 0.0d0
      DO I=1,NB
       PRODCWQWj2 = PRODCWQWj2 + CW12(J,I)*QW(I+J*(J-1)/2)
      ENDDO
      DO I=NA+1,J-1
       PRODCWQWj2 = PRODCWQWj2 + CW12(J,I)*QW(I+J*(J-1)/2)
      ENDDO
      DO I=J+1,NBF5
       PRODCWQWj2 = PRODCWQWj2 + CW12(J,I)*QW(J+I*(I-1)/2)
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END
      
! PRODCWQWjk
      FUNCTION PRODCWQWjk(NV,J,K,CW12,QW)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5,NV)::CW12
      DOUBLE PRECISION,DIMENSION(NBFT5)::QW
!-----------------------------------------------------------------------
!     PRODCWQWjk = CW12jk*QWj. Note: Term with I=J is not included
!-----------------------------------------------------------------------
      PRODCWQWjk=0.0d0
      DO I=1,J-1
       PRODCWQWjk = PRODCWQWjk + CW12(J,I,K)*QW(I+J*(J-1)/2)
      ENDDO
      DO I=J+1,NBF5
       PRODCWQWjk = PRODCWQWjk + CW12(J,I,K)*QW(J+I*(I-1)/2)
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! PRODCWQWjk1
      FUNCTION PRODCWQWjk1(NV,J,K,CW12,QW)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5,NV)::CW12
      DOUBLE PRECISION,DIMENSION(NBFT5)::QW
!-----------------------------------------------------------------------
!     PRODCWQWjk = CW12jk*QWj, j<NB (NSOC is excluded from the Sum)
!     Note: Term with I=J is not included
!-----------------------------------------------------------------------
      PRODCWQWjk1 = 0.0d0
      DO I=1,J-1
       PRODCWQWjk1 = PRODCWQWjk1 + CW12(J,I,K)*QW(I+J*(J-1)/2)
      ENDDO
      DO I=J+1,NB
       PRODCWQWjk1 = PRODCWQWjk1 + CW12(J,I,K)*QW(J+I*(I-1)/2)
      ENDDO
      DO I=NA+1,NBF5
       PRODCWQWjk1 = PRODCWQWjk1 + CW12(J,I,K)*QW(J+I*(I-1)/2)
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! PRODCWQWjk2
      FUNCTION PRODCWQWjk2(NV,J,K,CW12,QW)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5,NV)::CW12
      DOUBLE PRECISION,DIMENSION(NBFT5)::QW
!-----------------------------------------------------------------------
!     PRODCWQWjk = CW12jk*QWj, j>NA (NSOC is excluded from the Sum)
!     Note: Term with I=J is not included
!-----------------------------------------------------------------------
      PRODCWQWjk2 = 0.0d0
      DO I=1,NB
       PRODCWQWjk2 = PRODCWQWjk2 + CW12(J,I,K)*QW(I+J*(J-1)/2)
      ENDDO
      DO I=NA+1,J-1
       PRODCWQWjk2 = PRODCWQWjk2 + CW12(J,I,K)*QW(I+J*(J-1)/2)
      ENDDO
      DO I=J+1,NBF5
       PRODCWQWjk2 = PRODCWQWjk2 + CW12(J,I,K)*QW(J+I*(I-1)/2)
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END
      
! PRODROQWj
      FUNCTION PRODROQWj(J,RO,QW)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      DOUBLE PRECISION,DIMENSION(NBF5)::RO
      DOUBLE PRECISION,DIMENSION(NBFT5)::QW
!-----------------------------------------------------------------------
!     PRODROQWj = RO*QWj    Note: Term with I=J is not included
!-----------------------------------------------------------------------
      PRODROQWj=0.0d0
      DO I=1,J-1
       PRODROQWj = PRODROQWj + RO(I)*QW(I+J*(J-1)/2)
      ENDDO
      DO I=J+1,NBF5
       PRODROQWj = PRODROQWj + RO(I)*QW(J+I*(I-1)/2)
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! PRODROQWj0
      FUNCTION PRODROQWj0(J,RO,QW)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      DOUBLE PRECISION,DIMENSION(NBF5)::RO
      DOUBLE PRECISION,DIMENSION(NBFT5)::QW
!-----------------------------------------------------------------------
!     PRODROQWj0 = RO*QWj  (Sum only for NSOC)
!     Note: Term with I=J is not included
!-----------------------------------------------------------------------
      PRODROQWj0 = 0.0d0
      DO I=NB+1,J-1
       PRODROQWj0 = PRODROQWj0 + RO(I)*QW(I+J*(J-1)/2)
      ENDDO
      DO I=J+1,NA
       PRODROQWj0 = PRODROQWj0 + RO(I)*QW(J+I*(I-1)/2)
      ENDDO
!--------------------------------------------------------------------      
      RETURN
      END
      
! PRODROQWj1
      FUNCTION PRODROQWj1(J,RO,QW)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      DOUBLE PRECISION,DIMENSION(NBF5)::RO
      DOUBLE PRECISION,DIMENSION(NBFT5)::QW
!-----------------------------------------------------------------------
!     PRODROQWj = RO*QWj, j<NB<i (Sum only for NSOC)
!-----------------------------------------------------------------------
      PRODROQWj1 = 0.0d0
      DO I=NB+1,NA
       PRODROQWj1 = PRODROQWj1 + RO(J)*QW(J+I*(I-1)/2)
      ENDDO
!-----------------------------------------------------------------------      
      RETURN
      END

! PRODROQWj2
      FUNCTION PRODROQWj2(J,RO,QW)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      DOUBLE PRECISION,DIMENSION(NBF5)::RO
      DOUBLE PRECISION,DIMENSION(NBFT5)::QW
!-----------------------------------------------------------------------
!     PRODROQWj = RO*QWj, j>NA>i (Sum only for NSOC)
!-----------------------------------------------------------------------
      PRODROQWj2 = 0.0d0
      DO I=NB+1,NA
       PRODROQWj2 = PRODROQWj2 + RO(J)*QW(I+J*(J-1)/2)
      ENDDO
!-----------------------------------------------------------------------      
      RETURN
      END

! PRODDRQWjk1
      FUNCTION PRODDRQWjk1(NV,J,K,DR,QW)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      DOUBLE PRECISION,DIMENSION(NBF5,NV)::DR
      DOUBLE PRECISION,DIMENSION(NBFT5)::QW
!-----------------------------------------------------------------------
!     PRODDRQWjk1 = DRjk*QWj, j<NB (NSOC is excluded from the Sum)
!     Note: Term with I=J is not included
!-----------------------------------------------------------------------
      PRODDRQWjk1=0.0d0
      DO I=NB+1,NA
       PRODDRQWjk1 = PRODDRQWjk1 + DR(J,K)*QW(J+I*(I-1)/2)
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! PRODDRQWjk2
      FUNCTION PRODDRQWjk2(NV,J,K,DR,QW)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      DOUBLE PRECISION,DIMENSION(NBF5,NV)::DR
      DOUBLE PRECISION,DIMENSION(NBFT5)::QW
!-----------------------------------------------------------------------
!     PRODDRQWjk2 = DRjk*QWj, j>NA (NSOC is excluded from the Sum)
!     Note: Term with I=J is not included
!-----------------------------------------------------------------------
      PRODDRQWjk2=0.0d0
      DO I=NB+1,NA
       PRODDRQWjk2 = PRODDRQWjk2 + DR(J,K)*QW(I+J*(J-1)/2)
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! PRODWCWij
      FUNCTION PRODWCWij(ij,W,CW12)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5)::CW12
      DOUBLE PRECISION,DIMENSION(NSQ,NBF5)::W
!-----------------------------------------------------------------------
!     Sum W*CW12(1,*) by IQP
!-----------------------------------------------------------------------
      PRODWCWij = 0.0d0
      DO IQP=1,NBF5
       PRODWCWij = PRODWCWij + W(ij,IQP)*CW12(1,IQP)
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! PRODWCWij1
      FUNCTION PRODWCWij1(ij,W,CW12)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5)::CW12
      DOUBLE PRECISION,DIMENSION(NSQ,NBF5)::W
!-----------------------------------------------------------------------
!     Sum W*CW12(1,*) by IQP. NSOC terms are excluded from the Sum.
!-----------------------------------------------------------------------
      PRODWCWij1 = 0.0d0
      DO IQP=1,NB
       PRODWCWij1 = PRODWCWij1 + W(ij,IQP)*CW12(1,IQP)
      ENDDO
      DO IQP=NA+1,NBF5
       PRODWCWij1 = PRODWCWij1 + W(ij,IQP)*CW12(1,IQP)
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END
      
! PRODWCWijq
      FUNCTION PRODWCWijq(ij,IQ,W,CW12)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5)::CW12
      DOUBLE PRECISION,DIMENSION(NSQ,NBF5)::W
!-----------------------------------------------------------------------
!     Sum W*CW12(IQ,*) by IQP, IQ is not considered
!-----------------------------------------------------------------------
      PRODWCWijq = 0.0d0
      DO IQP=1,IQ-1
       PRODWCWijq = PRODWCWijq + W(ij,IQP)*CW12(IQ,IQP)
      ENDDO
      DO IQP=IQ+1,NBF5
       PRODWCWijq = PRODWCWijq + W(ij,IQP)*CW12(IQ,IQP)
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! PRODWCWijq1
      FUNCTION PRODWCWijq1(ij,IQ,W,CW12)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5)::CW12
      DOUBLE PRECISION,DIMENSION(NSQ,NBF5)::W
!-----------------------------------------------------------------------
!     Sum W*CW12(IQ,*) by IQP, IQ is not considered (IQ<NB)
!     NSOC terms are excluded from the Sum.
!-----------------------------------------------------------------------
      PRODWCWijq1 = 0.0d0
      DO IQP=1,IQ-1
       PRODWCWijq1 = PRODWCWijq1 + W(ij,IQP)*CW12(IQ,IQP)
      ENDDO
      DO IQP=IQ+1,NB
       PRODWCWijq1 = PRODWCWijq1 + W(ij,IQP)*CW12(IQ,IQP)
      ENDDO
      DO IQP=NA+1,NBF5
       PRODWCWijq1 = PRODWCWijq1 + W(ij,IQP)*CW12(IQ,IQP)
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END
      
! PRODWCWijq2
      FUNCTION PRODWCWijq2(ij,IQ,W,CW12)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5)::CW12
      DOUBLE PRECISION,DIMENSION(NSQ,NBF5)::W
!-----------------------------------------------------------------------
!     Sum W*CW12(IQ,*) by IQP, IQ is not considered (IQ>NA)
!     NSOC terms are excluded from the Sum.
!-----------------------------------------------------------------------
      PRODWCWijq2 = 0.0d0
      DO IQP=1,NB
       PRODWCWijq2 = PRODWCWijq2 + W(ij,IQP)*CW12(IQ,IQP)
      ENDDO
      DO IQP=NA+1,IQ-1
       PRODWCWijq2 = PRODWCWijq2 + W(ij,IQP)*CW12(IQ,IQP)
      ENDDO
      DO IQP=IQ+1,NBF5
       PRODWCWijq2 = PRODWCWijq2 + W(ij,IQP)*CW12(IQ,IQP)
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! SUMWij
      FUNCTION SUMWij(ij,W)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      DOUBLE PRECISION,DIMENSION(NSQ,NBF5)::W      
!-----------------------------------------------------------------------
!     Sum Wij for NSOC terms
!-----------------------------------------------------------------------
      SUMWij = 0.0d0
      DO IQP=NB+1,NA
       SUMWij = SUMWij + W(ij,IQP)
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! SUMWijq
      FUNCTION SUMWijq(ij,IQ,W)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      DOUBLE PRECISION,DIMENSION(NSQ,NBF5)::W
!-----------------------------------------------------------------------
!     Sum Wij for NSOC terms, IQ is not considered
!-----------------------------------------------------------------------
      SUMWijq = 0.0d0
      DO IQP=NB+1,IQ-1
        SUMWijq = SUMWijq + W(ij,IQP)
      ENDDO
      DO IQP=IQ+1,NA
        SUMWijq = SUMWijq + W(ij,IQP)
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! PRODWROij
      FUNCTION PRODWROij(ij,W,RO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      DOUBLE PRECISION,DIMENSION(NBF5)::RO
      DOUBLE PRECISION,DIMENSION(NSQ,NBF5)::W
!-----------------------------------------------------------------------
!     Sum W*RO. NSOC terms are excluded from the Sum.
!-----------------------------------------------------------------------
      PRODWROij = 0.0d0
      DO IQP=1,NB
       PRODWROij = PRODWROij + W(ij,IQP)*RO(IQP)
      ENDDO
      DO IQP=NA+1,NBF5
       PRODWROij = PRODWROij + W(ij,IQP)*RO(IQP)
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END
