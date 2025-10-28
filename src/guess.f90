!======================================================================!
!                                                                      !
!             G U E S S   &  R H F    S U B R O U T I N E S            !
!                                                                      !
!======================================================================!
!                                                                      !
!    GuessHJKRHF: Calculate the RHF as guess                           !
!    GuessCore  : HCore as guess                                       !
!                                                                      !
!======================================================================!

! GuessHJKRHF
      SUBROUTINE GuessHJKRHF(KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX,     &
                             EX,CS,CP,CD,CF,CG,CH,CI,NPRIMI,Cxyz,H,S,   &
                             EIG,VEC,XInteg,IXInteg,XIntegaux,NSHELL,   &
                             NAT,NBF,NSQ,NBFT,NINTEGtm,NINTEGAUXtm,     &
                             NINTEGt,NREC,XINTS,NSH2,IDONTW,INPUTC,     &
                             IPRINTOPT,ZAN,SIZE_ENV,ENV,ATM,NBAS,BAS,   &
                             IGTYP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                
!
      COMMON/USELIBCINT/ILIBCINT
      COMMON/INPNOF_RHF/CONVRHFDM,IRHF,IRHFTYP,NCONVRHF,MAXITRHF
      COMMON/INFOA/NATOMS,ICH,MUL,NUM,NQMT,NE,NA,NB 
      COMMON/CONV/ACURCY,EN,Etot,EHF,EHF0,DIFF,ITER,ICALP,ICBET
      LOGICAL SMCD
      COMMON/ERITYPE/IERITYP,IRITYP,IGEN,ISTAR,MIXSTATE,SMCD
      COMMON/INPNOF_HFID/KOOPMANS
      COMMON/INTOPT/CUTOFF,ISCHWZ,IECP,NECP
!      
      INTEGER :: NPRIMI,NSHELL,NAT,NBF,NSQ,NBFT,NINTEGtm,NINTEGAUXtm
      INTEGER :: NINTEGt,NREC,IDONTW,INPUTC,IPRINTOPT
      INTEGER(8),DIMENSION(NINTEGtm) :: IXInteg      
      INTEGER,DIMENSION(NSHELL) :: KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX
      DOUBLE PRECISION,DIMENSION(NPRIMI) :: EX,CS,CP,CD,CF,CG,CH,CI
      DOUBLE PRECISION,DIMENSION(NAT) :: ZAN
      DOUBLE PRECISION,DIMENSION(3,NAT) :: Cxyz                                                         
      DOUBLE PRECISION,DIMENSION(NBF) :: EIG
      DOUBLE PRECISION,DIMENSION(NSQ) :: VEC
      DOUBLE PRECISION,DIMENSION(NBFT) :: H,S
      DOUBLE PRECISION,DIMENSION(NINTEGtm) :: XInteg
      DOUBLE PRECISION,DIMENSION(NINTEGAUXtm) :: XIntegaux
      DOUBLE PRECISION,DIMENSION(NSH2) :: XINTS
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: TKIN,DipoInt,Q
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: OCCa,OCCb,DENa,DENb
!     LIBCINT
      INTEGER :: SIZE_ENV,NBAS
      DOUBLE PRECISION :: ENV(SIZE_ENV)
      INTEGER :: ATM(6,NAT), BAS(8,NBAS)
      INTEGER :: IGTYP
!-----------------------------------------------------------------------
!     IRHFTYP=1: Restricted Closed HF
!     IRHFTYP=2: Restricted Open HF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(IRHF==1)THEN
       IF(NA==NB)THEN
        IRHFTYP=1
       ELSE IF(NA>NB)THEN
        IRHFTYP=2
       END IF
       if(IERITYP==2)then
        WRITE(6,*)
        WRITE(6,*)'Sorry: IRHF=1 is not possible yet with ERITYP=RI'
        WRITE(6,*)'       Changing to IRHF=2'
        IRHF=2
        IRHFTYP=0
       end if
      ELSE
       IRHFTYP=0      
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     1e Integrals
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(TKIN(NBFT),DipoInt(3*NBFT))
      CALL OneElecInt(Cxyz,H,S,TKIN,DipoInt,NBFT,IPRINTOPT,             &
                      EX,CS,CP,CD,CF,CG,CH,CI,NPRIMI,KSTART,            &
                      KATOM,KTYPE,KNG,KLOC,KMIN,KMAX,NSHELL,ZAN,NAT,    &
                      IECP,SIZE_ENV,ENV,ATM,NBAS,BAS,IGTYP)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Initial Orbitals
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(OCCa(NBF),OCCb(NBF),DENa(NBFT),DENb(NBFT),Q(NSQ))
      if(INPUTC==0.or.IRHFTYP>0)then
       CALL GuessCore(OCCa,OCCb,DENa,DENb,EIG,VEC,Q,H,S,NBF,NSQ,NBFT,   &
                      IPRINTOPT)
      endif
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     2e- Integrals
!     NINTEGtm & NINTEGAUXtm: Maximum numbers of distinct 2e- Integrals 
!     Note: NINTEGt<NINTEGtm due to the CUTOFF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(IERITYP==1 .or. IERITYP==3)then   ! FULL or MIX
       if(ILIBCINT==0)then
!HONDO  CALL JandK(XInteg,IXInteg,NINTEGtm,NINTEGt,NREC,XINTS,NSH2,     &
!HONDO             IDONTW,IPRINTOPT,EX,CS,CP,CD,CF,CG,CH,CI,NPRIMI,     &
!HONDO             KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX,NSHELL,        &
!HONDO             Cxyz,NAT,SIZE_ENV,ENV,ATM,NBAS,BAS,IGTYP,ISCHWZ)
       else if(ILIBCINT==1)then
        CALL JandKl(XInteg,IXInteg,NINTEGtm,NINTEGt,NREC,XINTS,NSH2,    &
                    IDONTW,IPRINTOPT,EX,CS,CP,CD,CF,CG,CH,CI,NPRIMI,    &
                    KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX,NSHELL,       &
                    Cxyz,NAT,SIZE_ENV,ENV,ATM,NBAS,BAS,IGTYP)
       end if
      end if
      if(IERITYP==2 .or. IERITYP==3)then   ! RI or MIX
       if(ILIBCINT==0)then
!HONDO  CALL JandKaux(XIntegaux,NINTEGAUXtm,IDONTW,IPRINTOPT,NBF,       &
!HONDO                EX,CS,CP,CD,CF,CG,CH,CI,NPRIMI,KSTART,KATOM,      &
!HONDO                KTYPE,KNG,KLOC,KMIN,KMAX,NSHELL,Cxyz,NAT,         &
!HONDO                SIZE_ENV,ENV,ATM,NBAS,BAS,IGTYP)
       else if(ILIBCINT==1)then
        CALL JandKauxl(XIntegaux,NINTEGAUXtm,IDONTW,IPRINTOPT,NBF,      &
                       EX,CS,CP,CD,CF,CG,CH,CI,NPRIMI,KSTART,KATOM,     &
                       KTYPE,KNG,KLOC,KMIN,KMAX,NSHELL,Cxyz,NAT,        &
                       SIZE_ENV,ENV,ATM,NBAS,BAS,IGTYP)
       end if
      end if
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Calculate the Restricted HF Energy of Molecule
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(IRHFTYP==1)THEN
       CALL RHFCL(OCCa,DENa,EIG,VEC,Q,H,S,NBF,NSQ,NBFT,XInteg,          &
                  IXInteg,NINTEGtm,NINTEGt,NREC,IDONTW,IPRINTOPT,       &
                  ZAN,Cxyz)        
       DENb = DENa                                                       
      END IF                                                             
      IF(IRHFTYP==2)THEN                                                 
       CALL RHFOP(OCCa,OCCb,DENa,DENb,EIG,VEC,Q,H,S,NBF,NSQ,NBFT,       &
                  XInteg,IXInteg,NINTEGtm,NINTEGt,NREC,IDONTW,          &
                  IPRINTOPT,ZAN,Cxyz)
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Energy Analysis
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(IRHFTYP>0)CALL EnergyComp(EN,Etot,KATOM,KLOC,DENa,DENb,H,TKIN, &
                                   S,NSHELL,NAT,NBF,NBFT,NE,IRHFTYP)
!-----------------------------------------------------------------------
      DEALLOCATE(TKIN,DipoInt,OCCa,OCCb,DENa,DENb,Q)
      RETURN
      END

! GuessCore
      SUBROUTINE GuessCore(OCCa,OCCb,DENa,DENb,EIG,VEC,                 &
                           Q,H,S,NBF,NSQ,NBFT,IPRINTOPT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      INTEGER :: NBF,NSQ,NBFT,IPRINTOPT
      DOUBLE PRECISION,DIMENSION(NBF) :: OCCa,OCCb,EIG
      DOUBLE PRECISION,DIMENSION(NSQ) :: VEC,Q
      DOUBLE PRECISION,DIMENSION(NBFT) :: DENa,DENb,H,S
      DOUBLE PRECISION, ALLOCATABLE,DIMENSION(:,:) :: OVERLAP
!-----------------------------------------------------------------------
!     Initial Molecular Orbitals
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(IPRINTOPT==1)WRITE(6,'(A22)')' Guess for MOs = HCORE'
      CALL HCORE(EIG,H,S,VEC,Q,NBF,NBFT,NSQ)
      CALL CLENMO(VEC,NBF,NBF,1.0D-08,1.0D-05,IPRINTOPT)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Orthonormalize initial MOs if necessary
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(OVERLAP(NBF,NBF))
      CALL CPYTSQ(S,OVERLAP,NBF)       
      CALL CHECKORTHO(VEC,OVERLAP,IVIOORTHO,IPRINTOPT)
      if(IVIOORTHO/=0)then
       IF(IPRINTOPT==1)WRITE(6,'(/,A29)')' Orthogonalizing the orbitals'
       CALL ORTHONORMAL(NBF,NBF,NBFT,OVERLAP,VEC,1,IPRINTOPT)
      endif
      DEALLOCATE(OVERLAP)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Initial Density Matrix
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL INIDEN(VEC,OCCa,OCCb,DENa,DENb,NBF,NBFT,NSQ)
!-----------------------------------------------------------------------
      RETURN
      END

! HCORE
      SUBROUTINE HCORE(EIG,H,S,VEC,Q,NBF,NBFT,NSQ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/INFOA/NAT,ICH,MUL,NUM,NQMT,NE,NA,NB
      DOUBLE PRECISION,DIMENSION(NBF) :: EIG
      DOUBLE PRECISION,DIMENSION(NBFT) :: H,S
      DOUBLE PRECISION,DIMENSION(NSQ) :: VEC,Q
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: H0,HH,W
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: SCR,Hsq
      ALLOCATE(H0(NBFT),HH(NBFT),SCR(NBF,NBF))
!-----------------------------------------------------------------------
!         Initial orbitals by diagonalization of 1e Hamiltonian
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Q-matrix and NQMT
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      H0 = S
      CALL QMTSYM(H0,VEC,Q,EIG,SCR,NBF,NQMT,NBFT,NSQ) 
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     HH = Q*H0*Q (H0=HCORE)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      H0 = H
      CALL TFTRI(HH,H0,Q,SCR,NQMT,NBF,NBF)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Diagonalize HH (HH->Hsq)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(Hsq(NBF,NBF),W(NBF))
      CALL CPYTSQ(HH,Hsq,NBF)      
      CALL DIAG(NBF,Hsq,VEC,EIG,W)
      DEALLOCATE(Hsq,W)      
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Back Transform
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL TFSQB(VEC,Q,SCR,NQMT,NBF,NBF)    
!-----------------------------------------------------------------------
      DEALLOCATE(H0,HH,SCR)
      RETURN
      END

! CLENMO                                           
      SUBROUTINE CLENMO(VEC,NAO,NMO,TOLZ,TOLE,IPRINTOPT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      INTEGER :: NAO,NMO,IPRINTOPT
      DOUBLE PRECISION :: TOLZ,TOLE     
      DOUBLE PRECISION,DIMENSION(NAO,NMO) :: VEC                                            
!-----------------------------------------------------------------------
      NZER = 0                                                            
      NEQU = 0                                                            
      DO J=1,NMO                                                    
       DO 1 I=1,NAO                                                 
        VAL1 = ABS(VEC(I,J))                                          
        IF(VAL1==0.0D0)GO TO 1                                  
!       Zero teeny MOs
        IF(VAL1<=TOLZ)THEN
         VEC(I,J)=0.0D0                                            
         NZER = NZER+1                                              
         GO TO 1                                                
        END IF
!       equal magnitudes for nonzero orbitals
        IF(I==NAO)GO TO 1                                      
        KMIN = I+1                                                    
        DO K=KMIN,NAO                                           
         VAL2 = ABS(VEC(K,J))                                       
         TEST = ABS(1.0d0-VAL2/VAL1)                                
         IF(TEST<=TOLE)THEN
          VEC(K,J) = SIGN(VAL1,VEC(K,J))                        
          NEQU = NEQU+1                                           
         END IF
        END DO                                                    
    1  CONTINUE                                                       
      END DO
      IF(IPRINTOPT==1)WRITE(6,10)NEQU,NZER
   10 FORMAT(I10,' MOs were matched and',I10,' equalized to zero')
!-----------------------------------------------------------------------       
      RETURN                                                            
      END                                                               
      
! INIDEN
      SUBROUTINE INIDEN(V,OCCa,OCCb,DENa,DENb,NBF,NBFT,NSQ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/INFOA/NAT,ICH,MUL,NUM,NQMT,NE,NA,NB
      CHARACTER*4,ALLOCATABLE,DIMENSION(:) :: RLABMO
      DOUBLE PRECISION,DIMENSION(NBF) :: OCCa,OCCb
      DOUBLE PRECISION,DIMENSION(NSQ) :: V
      DOUBLE PRECISION,DIMENSION(NBFT) :: DENa,DENb
      ALLOCATE(RLABMO(NBF))
!-----------------------------------------------------------------------
      IF(NA>NBF)THEN
       WRITE(6,*)'Stop: NA > NBF for initial density matrix'
       CALL ABRT
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF( NB==(NE+MUL-1)/2 .and. NA==(NE-MUL+1)/2 )THEN     ! Closed
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       CALL VCLR(OCCa,1,NBF)
       DO I=1,NA
        OCCa(I) = 2.0D+00      ! Occupation Numbers
       END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Symmetry (C1)
       CALL SYMMOS(RLABMO,V,NBF,NQMT,NBF)
       CALL DMTX(DENa,V,OCCa,NA,NBF,NBF)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ELSE IF( NB==(NE-MUL+1)/2 .and. NA==(NE+MUL-1)/2 )THEN   ! Open
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Alpha Orbitals
!- - - - - - - - - - - - - - - - - - - - - - - - - - - -
       CALL VCLR(OCCa,1,NBF)
       DO I=1,NA
        OCCa(I) = 1.0D+00
       END DO
!      Symmetry (C1)
       CALL SYMMOS(RLABMO,V,NBF,NQMT,NBF)
       CALL DMTX(DENa,V,OCCa,NA,NBF,NBF)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Beta Orbitals
!- - - - - - - - - - - - - - - - - - - - - - - - - - - -
       CALL VCLR(OCCb,1,NBF)
       DO I=1,NB
        OCCb(I) = 1.0D+00
       END DO
!      Symmetry (C1)
       CALL SYMMOS(RLABMO,V,NBF,NQMT,NBF)
       CALL DMTX(DENb,V,OCCb,NB,NBF,NBF)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - -       
      END IF
!-----------------------------------------------------------------------
      DEALLOCATE(RLABMO)
      RETURN
      END

! DMTX                                             
      SUBROUTINE DMTX(D,V,X,M,N,NDIM)                                   
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION D(*),V(NDIM,M),X(M)                                     
!-----------------------------------------------------------------------
      DO I = 1,N                                                    
       IJ = I*(I-1)/2                                                 
       DO J = 1,I                                                 
        IJ = IJ + 1                                                 
        DUM = 0.0D+00                                                  
        DO K = 1,M                                              
         DUM = DUM + X(K)*V(I,K)*V(J,K)                             
        END DO
        D(IJ) = DUM                                                 
       END DO
      END DO
!-----------------------------------------------------------------------
      RETURN                                                            
      END                                                               
     
! QMTSYM
      SUBROUTINE QMTSYM(S,WRK,Q,EE,SCR,NBF,NQMT,NBFT,NSQ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DOUBLE PRECISION,DIMENSION(NBF) :: EE
      DOUBLE PRECISION,DIMENSION(NBFT) :: S
      DOUBLE PRECISION,DIMENSION(NSQ) :: WRK
      DOUBLE PRECISION,DIMENSION(NBF,NBF) :: Q,SCR
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: WRKq
!-----------------------------------------------------------------------
!     Symmetry Adapted Orthonormal Orbitals
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Q = 0.0D+00
      DO I=1,NBF
       Q(I,I) = 1.0D+00
      END DO
      CALL TFTRI(WRK,S,Q,SCR,NBF,NBF,NBF)
!      
      ALLOCATE(WRKq(NBF,NBF))
      CALL CPYTSQ(WRK,WRKq,NBF)            
      CALL DIAG(NBF,WRKq,Q,EE,SCR)
      DEALLOCATE(WRKq)      
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Small Eigenvalues -> Eliminate eigenvectors
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      EWARN  = 1.0D-05
      QMTTOL = 1.0D-16      !  Linear Dependence Threshhold
      IRPZER=0
      DO ISALC=1,NBF
       IF(EE(ISALC)<QMTTOL)IRPZER=IRPZER+1
      END DO
!
      ESMALL= 1.0D+09
      NSMALL= 0
      NWARN = 0
      JSALC =0
      DO ISALC=1,NBF
       IF(EE(ISALC)<ESMALL)ESMALL=EE(ISALC)
       IF(EE(ISALC)<EWARN) NWARN = NWARN+1
       IF(EE(ISALC)>=QMTTOL)THEN
        JSALC=JSALC+1
        EE(JSALC)=1.0D+00/SQRT(EE(ISALC))
        DO I=1,NBF
         Q(I,JSALC)=Q(I,ISALC)
        END DO
       ELSE
        NSMALL=NSMALL+1
       END IF
      END DO
      NQMT = NBF-NSMALL
      NWARN = NWARN-NSMALL
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Canonical Orthonormal Orbitals
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO J=1,NBF
       DO I=1,NBF
        IF(J>NQMT)Q(I,J) = 0.0D+00
        Q(I,J) = Q(I,J)*EE(J)
       END DO
      END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Back Transform to AO Space
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      WRK = 0.0D+00
      DO I=1,NBF
       IDIAG = I + (I-1)*NBF
       WRK(IDIAG) = 1.0D+00
      END DO
      CALL TFSQB(Q,WRK,SCR,NBF,NBF,NBF)
!-----------------------------------------------------------------------
      RETURN
      END

! SYMMOS
      SUBROUTINE SYMMOS(RLABMO,V,NAO,NMO,LDQV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DIMENSION V(LDQV,NMO)
      CHARACTER*4,DIMENSION(NMO) :: RLABMO
!-----------------------------------------------------------------------
      CALL STFASEd(V,LDQV,NAO,NMO)
      DO I=1,NMO
       RLABMO(I) = 'A   '
      END DO
!-----------------------------------------------------------------------
      RETURN
      END
      
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ! 
!                         R H F : NA = NB                              !
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ! 

! RHFCL      
      SUBROUTINE RHFCL(OCC,DEN,EIG,VEC,Q,H,S,NBF,NSQ,NBFT,XInteg,       &
                       IXInteg,NINTEGtm,NINTEGt,NREC,IDONTW,IPRINTOPT,  &
                       ZAN,Cxyz)   
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                
      LOGICAL HFDAMP,HFEXTRAP,HFDIIS
      COMMON/INPNOF_HFCONVTECH/HFDAMP,HFEXTRAP,HFDIIS 
      COMMON/INPNOF_RHF/CONVRHFDM,IRHF,IRHFTYP,NCONVRHF,MAXITRHF
      COMMON/INPNOF_THRESH/THRESHL,THRESHE,THRESHEC,THRESHEN
      COMMON/INPNOF_NTHRESH/NTHRESHL,NTHRESHE,NTHRESHEC,NTHRESHEN
      COMMON/INFOA/NAT,ICH,MUL,NUM,NQMT,NE,NA,NB
      COMMON/DIISSO/ETHRSH,MAXDII,IRAF
      COMMON/ACONV/RRSHFT,EXTTOL,DMPTOL                
      COMMON/CONV/ACURCY,EN,Etot,EHF,EHF0,DIFF,ITER,ICALP,ICBET
!
      LOGICAL PRVEC,CVGED,CVDENS,CVENGY,CVDIIS,NOTOPN     
      INTEGER :: NBF,NSQ,NBFT,NINTEGtm,NINTEGt,NREC,IDONTW,IPRINTOPT
      INTEGER(8),DIMENSION(NINTEGtm) :: IXInteg
      DOUBLE PRECISION,DIMENSION(NAT) :: ZAN
      DOUBLE PRECISION,DIMENSION(3,NAT) :: Cxyz
      DOUBLE PRECISION,DIMENSION(NBF) :: OCC,EIG
      DOUBLE PRECISION,DIMENSION(NSQ) :: VEC,Q
      DOUBLE PRECISION,DIMENSION(NBFT) :: DEN,H,S
      DOUBLE PRECISION,DIMENSION(NINTEGtm) :: XInteg
!
      CHARACTER*4,ALLOCATABLE,DIMENSION(:) :: AWRK      
      INTEGER,ALLOCATABLE,DIMENSION(:) :: IPDIIS,IODII      
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: FAO,QFQ,WRK
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: Fsq,Ssq,ERR
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: ADIIS,XDIIS,BDIIS
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: WRK1,WRK2,WRK3
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: FAO0,FAO1,FAO2
!-----------------------------------------------------------------------
!     Nuclear Energy
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      EN = ENUC(NAT,ZAN,Cxyz)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -      
!     Initialize Variables
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      PRVEC  = .FALSE.
      CVGED  = .FALSE.                                                        
!
      CVDENS = .FALSE.
      ACURCY = CONVRHFDM                                                      
      DTOL   = 2.0d0 * ACURCY      
      DIFFP  = 0.0d0                                                       
      DIFF   = 0.0d0                                                       
!      
      CVENGY = .FALSE.
      ETOL   = THRESHE
      ENGTOL = 1.0D-10
      EHF    = 0.0d0                                                        
      Etot   = 0.0d0
      DELE   = 0.0d0                                                       
!      
      CVDIIS = .FALSE. 
      DIITOL = 1.0D-07      
      NOTOPN = .TRUE.                                                   
      ETHRSH = 0.5d0                                                         
      ITDIIS = 1 
      ERDIIS = 0.0d0      
      MAXDII = 10      
!      
      DAMP   = 0.0d0                                                       
      DAMP0  = 0.0d0                                                      
      ITERV  = 0                                                                
      RRSHFT = 0.0d0
      EXTTOL = 1.0D-03                                                  
      DMPTOL = 1.0D-04                                                  
      ICALP  = 0
      ICBET  = 0                                                               
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                 
!     Header
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      WRITE(6,10)'RHF SCF Calculation'
      WRITE(6,20)HFEXTRAP,HFDAMP,HFDIIS 
      WRITE(6,30)ACURCY,ETOL
      WRITE(6,35)ENGTOL,DTOL
      IF(HFDIIS)WRITE(6,36)DIITOL,DTOL                                                         
      IF(HFDAMP)DAMP=1.0d0
      IF(HFDIIS)THEN
       WRITE(6,40)
      ELSE
       WRITE(6,50)
      END IF                                                            
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                 
      ALLOCATE(FAO(NBFT),QFQ(NBFT),Fsq(NBF,NBF),WRK(NBF))
      IF(HFDIIS)THEN
       ALLOCATE(Ssq(NBF,NBF),ERR(NBF,NBF))
       MAXIT2 = (MAXITRHF*MAXITRHF+MAXITRHF)/2                                             
       ALLOCATE(ADIIS(MAXDII*MAXDII),XDIIS(MAXITRHF),IPDIIS(MAXITRHF))
       ALLOCATE(BDIIS(MAXIT2),IODII(4*MAXDII)) 
      END IF       
      IF(HFDAMP.or.HFEXTRAP)THEN      
       ALLOCATE(WRK1(NSQ),WRK2(NSQ),WRK3(NSQ))
       ALLOCATE(FAO0(NBFT),FAO1(NBFT),FAO2(NBFT))      
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                 
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                 
!                            RHF Iterations
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                 
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                                             
      DO ITER=1,MAXITRHF 
!      Skeleton 2e- Fock matrix                           
       CALL HSTAR(NBF,DEN,FAO,NBFT,XInteg,IXInteg,NINTEGt,NREC,IDONTW)
!      Add H to form Fock Matrix
       FAO = FAO + H
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
       IF(.not.CVGED)THEN
!- - - - - - - - - - - - - - - - - -
!       HF Energy
!- - - - - - - - - - - - - - - - - -
        EHF0 = EHF 
        EHF1 = TRACEP(DEN,FAO,NBF)                                   
        EHF2 = TRACEP(DEN,H,NBF)                                    
        EHF = (EHF1+EHF2)/2.0d0                                                        
        Etot0 = Etot                                                      
        Etot  = EHF + EN
        DELE0 = DELE                                                      
        DELE  = Etot - Etot0
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       DIIS Interpolation (ERR=F*D*S-S*D*F)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        IF(HFDIIS)THEN
         CALL CPYTSQ(FAO,Fsq,NBF)    ! FAO->Fsq
         CALL CPYTSQ(S,Ssq,NBF)      ! S->Ssq                 
         CALL DIISER(Fsq,DEN,Ssq,ERR,WRK,NBF,NBFT)                                      
         CALL DIISInter(ITDIIS,Q,FAO,ERR,Fsq,ADIIS,                     &
                        XDIIS,IPDIIS,BDIIS,IODII,WRK,                   &
                        NBF,NBFT,NSQ,MAXITRHF,MAXIT2,                   &
                        4*MAXDII,ERDIIS,NOTOPN)
         IF(ITDIIS>1)THEN                                           
          IF(HFDAMP)DAMP=1.0d0                                                   
          RRSHFT=0.0d0                                                 
         END IF                                                         
        END IF                                                            
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       Damping and Extrapolation of the Fock matrix
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        IF(ITER==2)DEAVG = ABS(DELE)                                   
        IF(ITDIIS==1)THEN
         IF(HFDAMP.and.ITER>2)THEN
          DEAVG = (ABS(DELE)+ABS(DELE0)+0.2d0*DEAVG)/2.2D0 
          DAMP0 = DAMP
          CALL DAMPD(DELE,DELE0,DEAVG,DAMP,ACURCY,DIFF,DIFFP,1.0D-02)
         END IF
         IF(DAMP<0.0d0)DAMP = 0.0d0
         IF(HFDAMP.or.HFEXTRAP)THEN
          CALL EXTRAPOL(DELE,DAMP,DAMP0,FAO,WRK1,WRK2,WRK3,FAO0,FAO1,   &
                        FAO2,NBF,NBFT,ITERV,1,1,ITER,ICALP,ICBET)
         END IF
        END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                              
!       Transform the Fock matrix to the Q-basis: QFQ = Q*FAO*Q
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
        CALL TFTRI(QFQ,FAO,Q,WRK,NQMT,NBF,NBF)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                                      
!       Diagonalize F' (QFQ->Fsq)                              
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                              
        CALL CPYTSQ(QFQ,Fsq,NBF)      
        CALL DIAG(NBF,Fsq,VEC,EIG,WRK)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                              
!       Back-transfrom the eigenvectors to the AO-basis
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                              
        CALL TFSQB(VEC,Q,WRK,NQMT,NBF,NBF) 
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                              
!       New Density Matrix (DEN), QFQ = Previous Density Matrix
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                              
        CALL DCOPY(NBFT,DEN,1,QFQ,1) 
        CALL DMTX(DEN,VEC,OCC,NA,NBF,NBF)   ! NOC=NA=NB
        DIFFP = DIFF                                                      
        CALL DDIFF(QFQ,DEN,NBFT,DIFF)                                
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                              
!       Printing iteration                                               
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                              
        IF(HFDIIS)THEN
         WRITE(6,60)ITER,Etot,DELE,DIFF,ERDIIS                              
        ELSE
         WRITE(6,70)ITER,Etot,DELE,DIFF
        END IF
        ICALP = ICALP+1                                                   
        ICBET = ICBET+1                                                   
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       Convergence: CVDENS, CVENGY, CVGED
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        CVDENS = DIFF<ACURCY  .and.  ABS(DELE)<ETOL                       
        CVENGY = ABS(DELE)<ENGTOL  .and.  DIFF<DTOL                
        IF(HFDIIS) CVDIIS = ERDIIS<DIITOL  .and.  DIFF<DTOL 
        CVGED  = CVDENS .or. CVENGY .or. CVDIIS
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                              
       END IF
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!      Check Convergence
!-- -- -- -- -- -- -- -- --
       IF(CVGED)THEN                                                    
        IF(CVDENS)THEN                                                
         WRITE(6,80)                                   
         GO TO 1                                                   
        END IF                                                         
        IF(CVDIIS) THEN                                                
         WRITE(6,90)                                   
         GO TO 1                                                   
        END IF                                                         
        IF(CVENGY) THEN                                                
         WRITE(6,100)                                   
         GO TO 1                                                   
        END IF                                                         
       END IF                                                            
      END DO
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!     SCF is unconverged
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
      WRITE(6,'(/1X,A24)')'Stop: SCF is unconverged'
      STOP                                                         
!-- -- -- -- -- -- -- -- -- --
!     Succesful Convergence
!-- -- -- -- -- -- -- -- -- --
    1 CONTINUE                                                          
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                 
!     Print Final Results
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(NTHRESHE<=10)then
       WRITE(6,110)Etot,ITER
      else
       WRITE(6,111)Etot,ITER
      end if
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                       
      ALLOCATE(AWRK(NQMT))
      CALL SYMMOS(AWRK,VEC,NBF,NQMT,NBF)      
      IF(PRVEC.and.IPRINTOPT==1)THEN
       WRITE(6,120)
       CALL PRINTEVECS(VEC,EIG,AWRK,NQMT,NBF)
      END IF      
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -       
      IF(.not.NOTOPN)CLOSE(UNIT=20,STATUS='DELETE')
      NOTOPN = .TRUE.     
!-----------------------------------------------------------------------
      DEALLOCATE(FAO,QFQ,Fsq,WRK,AWRK)
      IF(HFDIIS)DEALLOCATE(Ssq,ERR,ADIIS,XDIIS,IPDIIS,BDIIS,IODII)
      IF(HFDAMP.or.HFEXTRAP)THEN
       DEALLOCATE(WRK1,WRK2,WRK3,FAO0,FAO1,FAO2)
      ENDIF
      RETURN
!-----------------------------------------------------------------------     
   10 FORMAT(/1X,A19,/1X,19(1H-))  
   20 FORMAT(1X,'EXTRAP=',L1,'  DAMP=',L1,'  DIIS=',L1)  
   30 FORMAT(1X,'DMat Conv =',1P,E8.1,2X,                               &
                'Ener Tol =',1P,E8.1)                                    
   35 FORMAT(1X,'Ener Conv =',1P,E8.1,2X,                               &
                'DMat Tol =',1P,E8.1)                                    
   36 FORMAT(1X,'DIIS Conv =',1P,E8.1,2X,                               &
                'DMat Tol =',1P,E8.1)                                    
   40 FORMAT(/' ITER      TOTAL ENERGY        E CHANGE     ',           &
              'DEN. CHANGE    DIIS ERROR')                               
   50 FORMAT(/' ITER      TOTAL ENERGY        E CHANGE     ',           &
              'DEN. CHANGE')                                             
   60 FORMAT(1X,I3,F20.10,F17.10,2F14.9)                                 
   70 FORMAT(1X,I3,F20.10,F17.10, F14.9)                                 
   80 FORMAT(/10X,17(1H-)/10X,17HDENSITY CONVERGED/10X,17(1H-))          
   90 FORMAT(/10X,14(1H-)/10X,14HDIIS CONVERGED/10X,14(1H-))             
  100 FORMAT(/10X,16(1H-)/10X,16HENERGY CONVERGED/10X,16(1H-))           
  110 FORMAT(/1X,'FINAL RHF ENERGY IS',F20.10,' AFTER',I4,              &
                 ' ITERATIONS')
  111 FORMAT(/1X,'FINAL RHF ENERGY IS',F25.15,' AFTER',I4,              &
                 ' ITERATIONS')
  120 FORMAT(/10X,12HEigenvectors/10X,12(1H-))                   
!-----------------------------------------------------------------------           
      END

! HSTAR                                            
      SUBROUTINE HSTAR(NBF,Den,Skeleton,NBFT,XInteg,IXInteg,NINTEGt,    &
                       NREC,IDONTW)                  
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/INTFIL/NINTMX                 
      DOUBLE PRECISION,DIMENSION(NBFT) :: Den,Skeleton
      INTEGER(8),DIMENSION(NINTEGt) :: IXInteg
      DOUBLE PRECISION,DIMENSION(NINTEGt) :: XInteg    
      INTEGER(8),ALLOCATABLE,DIMENSION(:) :: IBUF
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: BUF
      INTEGER(8) :: LABEL, I, J, K, L
!-----------------------------------------------------------------------
!             Form the Skeleton of the Fock matrix for RHF
!-----------------------------------------------------------------------
      IF(IDONTW==0)THEN
       REWIND(1)
       ALLOCATE(IBUF(NINTMX),BUF(NINTMX))
      END IF
!      
      CALL VCLR(Skeleton,1,NBFT)                             
      I = 0                                                             
      J = 0                                                             
      K = 0                                                             
      L = 0                                                             
      NXInteg = 0                                                           
      IRECORD = 0 
!-.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.-
    1 CONTINUE
      IF(IDONTW==1)THEN        ! IDONTW=1  (IXInteg,XInteg)
       IRECORD = IRECORD + 1
       IJBUF = (IRECORD-1)*NINTMX
       if(IRECORD<NREC)then               
        NINTEG = NINTMX
        NXInteg = NINTEG
       else
        NINTEG = NINTEGt - (NREC-1)*NINTMX
        NXInteg = -NINTEG
       end if
       IF(NXInteg==0)GO TO 2 
       IF(NINTEG>NINTMX)CALL ABRT 
!       
       DO M = 1,NINTEG
        LABEL = IXInteg(IJBUF+M)
        CALL LABELIJKL(LABEL,I,J,K,L)
        VAL = XInteg(IJBUF+M)
        NIJ = (I*I-I)/2+J                                                  
        NKL = (K*K-K)/2+L                                                  
        NIK = (I*I-I)/2+K                                                  
        NIL = (I*I-I)/2+L                                                  
        IF(J<K)GO TO 32                                         
        NJK = (J*J-J)/2+K                                                  
        NJL = (J*J-J)/2+L                                                  
        GO TO 42                                                      
   32   NJK = (K*K-K)/2+J                                                  
        IF(J<L) GO TO 52                                         
        NJL = (J*J-J)/2+L                                                  
        GO TO 42                                                      
   52   NJL = (L*L-L)/2+J                                                  
   42   CONTINUE                                                       
        VAL4 = 4.0d0*VAL
        Skeleton(NIJ) = Skeleton(NIJ) + VAL4*Den(NKL)                                    
        Skeleton(NKL) = Skeleton(NKL) + VAL4*Den(NIJ)                                    
        Skeleton(NIK) = Skeleton(NIK) - VAL*Den(NJL)                                     
        Skeleton(NIL) = Skeleton(NIL) - VAL*Den(NJK)                                     
        Skeleton(NJK) = Skeleton(NJK) - VAL*Den(NIL)                                     
        Skeleton(NJL) = Skeleton(NJL) - VAL*Den(NIK)                                     
       END DO
      ELSE                    ! IDONTW=0 (IBUF,XBUF)
       READ(1)NXInteg,IBUF,BUF
       IF(NXInteg==0)GO TO 2                                        
       NINTEG = IABS(NXInteg)                                                
       IF(NINTEG>NINTMX)CALL ABRT
!       
       DO M = 1,NINTEG                                                 
        LABEL = IBUF(M) 
        CALL LABELIJKL(LABEL,I,J,K,L)                                          
        VAL = BUF(M)                                                
        NIJ = (I*I-I)/2+J                                                  
        NKL = (K*K-K)/2+L                                                  
        NIK = (I*I-I)/2+K                                                  
        NIL = (I*I-I)/2+L                                                  
        IF(J < K) GO TO 31                                         
        NJK = (J*J-J)/2+K                                                  
        NJL = (J*J-J)/2+L                                                  
        GO TO 41                                                      
   31   NJK = (K*K-K)/2+J                                                  
        IF(J<L) GO TO 51                                         
        NJL = (J*J-J)/2+L                                                  
        GO TO 41                                                      
   51   NJL = (L*L-L)/2+J                                                  
   41   CONTINUE                                                       
        VAL4 = (VAL+VAL)+(VAL+VAL)  
        Skeleton(NIJ) = Skeleton(NIJ) + VAL4*Den(NKL)                                    
        Skeleton(NKL) = Skeleton(NKL) + VAL4*Den(NIJ)                                    
        Skeleton(NIK) = Skeleton(NIK) - VAL*Den(NJL)                                     
        Skeleton(NIL) = Skeleton(NIL) - VAL*Den(NJK)                                     
        Skeleton(NJK) = Skeleton(NJK) - VAL*Den(NIL)                                     
        Skeleton(NJL) = Skeleton(NJL) - VAL*Den(NIK)                                     
       END DO
      END IF
!-.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.-
      IF(NXInteg>0)GO TO 1
    2 CONTINUE
!-----------------------------------------------------------------------                                                
      IF(IDONTW==0)THEN
       REWIND(1)
       DEALLOCATE(IBUF,BUF)
      END IF
!
      CALL DSCAL(NBFT,0.5d0,Skeleton,1)                                       
      II = 0                                                            
      DO I=1,NBF                                                       
       II = II + I                                                    
       Skeleton(II) = Skeleton(II) + Skeleton(II)                                            
      ENDDO 
!-----------------------------------------------------------------------  
      RETURN                                                            
      END                                                               

! DDIFF                                            
      SUBROUTINE DDIFF(D0,D1,L2,DIFF)                                   
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION D0(L2),D1(L2)                                             
!     Largest Absolute Change in Density Matrix
      DIFF = 0.0D+00                                                    
      DO I = 1,L2                                                   
       DIFIJ = ABS(D0(I)-D1(I))                                          
       IF (DIFIJ > DIFF) DIFF = DIFIJ
      END DO                                 
      RETURN                                                            
      END

! PRINTEVECS                                            
      SUBROUTINE PRINTEVECS(V,EE,RLABMO,NMO,NBF)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                      
      COMMON/RUNLAB/BFLAB(8192)            
      DIMENSION V(NBF,NMO),EE(NMO)
      CHARACTER*4 RLABMO(NMO)                                               
!-----------------------------------------------------------------------    
      MAX = 5                                                           
      IMAX = 0
    1 IMIN = IMAX+1                                                     
      IMAX = IMAX+MAX                                                   
      IF (IMAX > NMO) IMAX = NMO                                     
      WRITE (6,*)                                                   
      WRITE (6,'(15X,10(4X,I4,3X))')(I,I=IMIN,IMAX)                            
      WRITE (6,'(15X,10F11.4)')(EE(I),I=IMIN,IMAX)                            
      WRITE (6,'(16X,10(5X,A4,2X))')(RLABMO(I),I=IMIN,IMAX)                            
      DO J = 1,NBF                                                  
       WRITE (6,'(I5,2X,A8,10F11.6)')J,BFLAB(J),(V(J,I),I=IMIN,IMAX)              
      END DO
      IF(IMAX < NMO)GO TO 1                                      
!-----------------------------------------------------------------------
      RETURN                                                                                                                                   
      END                                                               

! EnergyComp
      SUBROUTINE EnergyComp(EN,Etot,KATOM,KLOC,DENa,DENb,H,TKIN,S,      &
                            NSHELL,NAT,NBF,NBFT,NE,ITYP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/VirialRatio/Virial                             ! RHF
      INTEGER,DIMENSION(NSHELL) :: KATOM,KLOC      
      DOUBLE PRECISION,DIMENSION(NBFT) :: DENa,DENb,H,TKIN,S
      INTEGER,ALLOCATABLE,DIMENSION(:) :: LIMLOW,LIMSUP            
!-----------------------------------------------------------------------
!     Basis Functions per Atom: Determine LIMLOW & LIMSUP
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(LIMLOW(NAT),LIMSUP(NAT))
      CALL AOLIM(NBF,NAT,NSHELL,LIMLOW,LIMSUP,KATOM,KLOC)

      Eone = TRACEP(DENa,H,NBF)
      IF(ITYP==2)Eone = Eone + TRACEP(DENb,H,NBF)

      Ecin = TRACEP(DENa,TKIN,NBF)
      IF(ITYP==2)Ecin = Ecin + TRACEP(DENb,TKIN,NBF)

      PSINRM = TRACEP(DENa,S,NBF)
      IF(ITYP==2)PSINRM = PSINRM + TRACEP(DENb,S,NBF)
      
      PSINRM = PSINRM/NE

      Vee = Etot - Eone - EN
      Vne = Eone - Ecin
      Vnn = EN
      Vtot = Vne + Vnn + Vee
      Virial = -Vtot/Ecin

      WRITE(6,10) PSINRM
      WRITE(6,20) Eone,Vee,Vnn,Etot
      WRITE(6,30) Vee,Vne,Vnn,Vtot,Ecin,Virial

      DEALLOCATE(LIMLOW,LIMSUP)
      RETURN
!-----------------------------------------------------------------------
   10 FORMAT(/1X,'        WAVEFUNCTION NORMALIZATION =',F19.10)
   20 FORMAT(/1X,'               ONE ELECTRON ENERGY =',F19.10/         &
              1X,'               TWO ELECTRON ENERGY =',F19.10/         &
              1X,'          NUCLEAR REPULSION ENERGY =',F19.10/         &
             38X,18(1H-)/                                               &
              1X,'                  RHF TOTAL ENERGY =',F19.10)          
   30 FORMAT(/1X,'ELECTRON-ELECTRON POTENTIAL ENERGY =',F19.10/         &
              1X,' NUCLEUS-ELECTRON POTENTIAL ENERGY =',F19.10/         &
              1X,'  NUCLEUS-NUCLEUS POTENTIAL ENERGY =',F19.10/         &
             38X,18(1H-)/                                               &
              1X,'            TOTAL POTENTIAL ENERGY =',F19.10/         &
              1X,'              TOTAL KINETIC ENERGY =',F19.10/         &
              1X,'                VIRIAL RATIO (V/T) =',F19.10)
!-----------------------------------------------------------------------              
      END

! AOLIM
      SUBROUTINE AOLIM(NBF,NAT,NSHELL,LIMLOW,LIMSUP,KATOM,KLOC)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER,DIMENSION(NAT) :: LIMLOW,LIMSUP
      INTEGER,DIMENSION(NSHELL) :: KATOM,KLOC
!-----------------------------------------------------------------------
      LIMLOW(1) = 1
      LAT = 1
      J = 1
      DO I = 1,NSHELL
       IAT = KATOM(I)
       IF (LAT == IAT) GO TO 10
       LAT = IAT
       LIMSUP(J) = KLOC(I)-1
       J = J+1
       LIMLOW(J) = KLOC(I)
   10  CONTINUE
      END DO
      LIMSUP(J) = NBF
      IF(J<NAT)THEN
       JP1=J+1
       DO J=JP1,NAT
        LIMLOW(J)=NBF
        LIMSUP(J)=1
       END DO
      END IF
!-----------------------------------------------------------------------
      RETURN
      END

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ! 
!                       R O H F : NA > NB                              !
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ! 

! RHFOP                                            
      SUBROUTINE RHFOP(OCCa,OCCb,DENa,DENb,EIG,VEC,Q,H,S,NBF,NSQ,NBFT,  &
                       XInteg,IXInteg,NINTEGtm,NINTEGt,NREC,            &
                       IDONTW,IPRINTOPT,ZAN,Cxyz)                                 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                 
      LOGICAL HFDAMP,HFEXTRAP,HFDIIS
      COMMON/INPNOF_HFCONVTECH/HFDAMP,HFEXTRAP,HFDIIS 
      COMMON/INPNOF_RHF/CONVRHFDM,IRHF,IRHFTYP,NCONVRHF,MAXITRHF
      COMMON/INFOA/NAT,ICH,MUL,NUM,NQMT,NE,NA,NB
      COMMON/DIISSO/ETHRSH,MAXDII,IRAF
      COMMON/ACONV/RRSHFT,EXTTOL,DMPTOL                
      COMMON/CONV  /ACURCY,EN,Etot,EHF,EHF0,DIFF,ITER,ICALP,ICBET
!
      LOGICAL PRVEC,CVGED,CVDENS,CVENGY,CVDIIS,NOTOPN          
      INTEGER(8),DIMENSION(NINTEGtm) :: IXInteg
      INTEGER :: NBF,NSQ,NBFT,NINTEGtm,NINTEGt,NREC,IDONTW,IPRINTOPT
      DOUBLE PRECISION,DIMENSION(NAT) :: ZAN
      DOUBLE PRECISION,DIMENSION(3,NAT) :: Cxyz
      DOUBLE PRECISION,DIMENSION(NBF) :: OCCa,OCCb,EIG
      DOUBLE PRECISION,DIMENSION(NSQ) :: VEC,Q
      DOUBLE PRECISION,DIMENSION(NBFT) :: DENa,DENb,H,S
      DOUBLE PRECISION,DIMENSION(NINTEGtm) :: XInteg
!
      CHARACTER*4,ALLOCATABLE,DIMENSION(:) :: AWRK      
      INTEGER,ALLOCATABLE,DIMENSION(:) :: IPDIIS,IODII            
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: FA,FB,FMO,QFQ,WRK
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: Fsq,Ssq,ERR
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: ADIIS,XDIIS
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: BDIIS
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: WRK1,WRK2,WRK3
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: FAO0,FAO1,FAO2
!-----------------------------------------------------------------------
!     Nuclear Energy
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      EN = ENUC(NAT,ZAN,Cxyz)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -      
!     Initialize Variables
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      PRVEC  = .FALSE.
      CVGED  = .FALSE.                                                        
!
      CVDENS = .FALSE.
      ACURCY = CONVRHFDM                                                      
      DTOL   = 2.0d0 * ACURCY      
      DIFFP  = 0.0d0                                                       
      DIFF   = 0.0d0
      DIFFA  = 0.0d0                                                       
      DIFFB  = 0.0d0                                                       
!      
      CVENGY = .FALSE.
      ETOL   = 1.0D-08      
      ENGTOL = 1.0D-10
      EHF    = 0.0d0                                                        
      Etot   = 0.0d0
      DELE   = 0.0d0                                                       
!      
      CVDIIS = .FALSE. 
      DIITOL = 1.0D-07      
      NOTOPN = .TRUE.                                                   
      ETHRSH = 0.5d0                                                         
      ITDIIS = 1 
      ERDIIS = 0.0d0      
      MAXDII = 10      
!      
      DAMP   = 0.0d0                                                       
      DAMP0  = 0.0d0                                                      
      ITERV  = 0                                                                
      RRSHFT = 0.0d0
      EXTTOL = 1.0D-03                                                  
      DMPTOL = 1.0D-04                                                  
      ICALP  = 0
      ICBET  = 0
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                 
!     Header
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      WRITE(6,10)'ROHF SCF Calculation'
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      WRITE(6,20)HFEXTRAP,HFDAMP,HFDIIS 
      WRITE(6,30)ACURCY,ETOL
      WRITE(6,35)ENGTOL,DTOL
      IF(HFDIIS)WRITE(6,36)DIITOL,DTOL                                                         
      IF(HFDAMP)DAMP=1.0d0
      IF(HFDIIS)THEN
       WRITE(6,40)
      ELSE
       WRITE(6,50)
      END IF                                                            
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                 
      ALLOCATE(FA(NBFT),FB(NBFT),FMO(NBFT),QFQ(NBFT))
      ALLOCATE(Fsq(NBF,NBF),WRK(NBF),WRK1(NSQ))
      IF(HFDIIS)THEN
       ALLOCATE(Ssq(NBF,NBF),ERR(NBF,NBF))
       MAXIT2 = (MAXITRHF*MAXITRHF+MAXITRHF)/2                                             
       ALLOCATE(ADIIS(MAXDII*MAXDII),XDIIS(MAXITRHF),IPDIIS(MAXITRHF))
       ALLOCATE(BDIIS(MAXIT2),IODII(4*MAXDII)) 
      END IF       
      IF(HFDAMP.or.HFEXTRAP)THEN      
       ALLOCATE(WRK2(NSQ),WRK3(NSQ))
       ALLOCATE(FAO0(NBFT),FAO1(NBFT),FAO2(NBFT))      
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                 
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                 
!                            ROHF Iterations
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                 
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO ITER=1,MAXITRHF
!      Skeletons 2e- Fock matrix  (FA,FB)                         
       CALL HSTARRO(NBF,DENa,FA,DENb,FB,NBFT,XInteg,IXInteg,            &
                    NINTEGt,NREC,IDONTW)    
!      Add H to form Alpha and Beta Fock Matrices in AO basis 
       FA = FA + H
       FB = FB + H       
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
       IF(.not.CVGED)THEN
!- - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       HF Energy
!- - - - - - - - - - - - - - - - - - - - - - - - - - - -
        EHF0 = EHF                                                        
        EHFA = TRACEP(DENa,H,NBF) + TRACEP(DENa,FA,NBF)        
        EHFB = TRACEP(DENb,H,NBF) + TRACEP(DENb,FB,NBF)        
        EHF  = (EHFA+EHFB)/2.0d0                                            
        Etot0= Etot                                                       
        Etot = EHF + EN                                                     
        DELE0 = DELE                                                       
        DELE = Etot-Etot0                                                 
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       High Spin Total MO Fock Matrix (FMO) -> Output: FA in AO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        CALL ROFOCK(FA,FB,FMO,VEC,WRK,S,WRK1,NA,NB,NQMT,NBF,NBFT,NSQ)                  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
!       DIIS Interpolation (ERR=F*D*S-S*D*F,D=DA+DB)   (F=FA)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
        IF(HFDIIS)THEN
         DENa = DENa + DENb         ! D = DENa
         CALL CPYTSQ(FA,Fsq,NBF)    ! FA->Fsq
         CALL CPYTSQ(S,Ssq,NBF)     ! S->Ssq                 
         CALL DIISER(Fsq,DENa,Ssq,ERR,WRK,NBF,NBFT)                                      
         DENa = DENa - DENb         
         CALL DIISInter(ITDIIS,Q,FA,ERR,Fsq,ADIIS,XDIIS,IPDIIS,BDIIS,   &
         IODII,WRK,NBF,NBFT,NSQ,MAXITRHF,MAXIT2,4*MAXDII,ERDIIS,NOTOPN)
         IF(ITDIIS>1)THEN                                           
          IF(HFDAMP)DAMP=1.0d0                                                   
          RRSHFT=0.0d0                                                 
         END IF                                                         
        END IF                                                            
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
!       Damping and Extrapolation of the Fock Matrix (FA)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        IF(ITER==2)DEAVG = ABS(DELE)                                   
        IF(ITDIIS==1)THEN
         IF(HFDAMP.and.ITER>2)THEN
          DEAVG = (ABS(DELE)+ABS(DELE0)+0.2d0*DEAVG)/2.2D0          
          DAMP0 = DAMP
          CALL DAMPD(DELE,DELE0,DEAVG,DAMP,ACURCY,DIFF,DIFFP,1.0D-02)
         END IF
         IF(DAMP<0.0d0)DAMP = 0.0d0
         IF(HFDAMP.or.HFEXTRAP)THEN
          CALL EXTRAPOL(DELE,DAMP,DAMP0,FA,WRK1,WRK2,WRK3,FAO0,FAO1,    &
                        FAO2,NBF,NBFT,ITERV,1,1,ITER,ICALP,ICBET)
         END IF
        END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                              
!       Transform the Fock matrix to the Q-basis: QFQ = Q*FA*Q
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
        CALL TFTRI(QFQ,FA,Q,WRK,NQMT,NBF,NBF)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                                      
!       Diagonalize F' (QFQ->Fsq)                              
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                              
        CALL CPYTSQ(QFQ,Fsq,NBF)      
        CALL DIAG(NBF,Fsq,VEC,EIG,WRK)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                              
!       Back-transfrom the eigenvectors to the AO-basis
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                              
        CALL TFSQB(VEC,Q,WRK,NQMT,NBF,NBF) 
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                              
!       New Density Matrix (DEN), QFQ = Previous Density Matrix
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                              
        CALL DCOPY(NBFT,DENa,1,QFQ,1) 
        CALL DMTX(DENa,VEC,OCCa,NA,NBF,NBF)
        CALL DDIFF(QFQ,DENa,NBFT,DIFFA) 
!        
        IF(NB>0)THEN                                                  
         CALL DCOPY(NBFT,DENb,1,WRK1,1)                             
         CALL DMTX(DENb,VEC,OCCb,NB,NBF,NBF)                      
         CALL DDIFF(WRK1,DENb,NBFT,DIFFB)                           
        ELSE                                                              
         WRK1  = 0.0d0
         DENb  = 0.0d0
         DIFFB = 0.0d0                                                     
        END IF
!        
        DIFFP = DIFF                                                      
        DIFF  = DIFFA + DIFFB
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                              
!       Printing iteration                                               
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                              
        IF(HFDIIS)THEN
         WRITE(6,60)ITER,Etot,DELE,DIFF,ERDIIS                              
        ELSE
         WRITE(6,70)ITER,Etot,DELE,DIFF
        END IF
        ICALP = ICALP+1                                                   
        ICBET = ICBET+1
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       Convergence: CVDENS, CVENGY, CVGED
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        CVDENS = DIFF<ACURCY  .and.  ABS(DELE)<ETOL                       
        CVENGY = ABS(DELE)<ENGTOL  .and.  DIFF<DTOL                
        IF(HFDIIS) CVDIIS = ERDIIS<DIITOL  .and.  DIFF<DTOL 
        CVGED  = CVDENS .or. CVENGY .or. CVDIIS
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                              
       END IF
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --       
!      Check Convergence
!-- -- -- -- -- -- -- -- --
       IF(CVGED)THEN                                                    
        IF(CVDENS)THEN                                                
         WRITE(6,80)                                   
         GO TO 1                                                   
        END IF                                                         
        IF(CVDIIS) THEN                                                
         WRITE(6,90)                                   
         GO TO 1                                                   
        END IF                                                         
        IF(CVENGY) THEN                                                
         WRITE(6,100)                                   
         GO TO 1                                                   
        END IF                                                         
       END IF                                                            
      END DO
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!     SCF is unconverged
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
      WRITE(6,'(/1X,A24)')'Stop: SCF is unconverged'
      STOP                                                         
!-- -- -- -- -- -- -- -- -- --
!     Succesful Convergence
!-- -- -- -- -- -- -- -- -- --
    1 CONTINUE                                                          
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                 
!     Print Final Results
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      WRITE(6,110)Etot,ITER  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                       
!     Spin values
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      SPINZ = 0.0d0                                                         
      S2 = 0.0d0
      CALL SPIN(SPINZ,S2,DENa,DENb,S,WRK1,WRK,NA,NB,NBF,NBFT)                                         
      WRITE(6,120)SPINZ,S2                                   
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      ALLOCATE(AWRK(NQMT))
      CALL SYMMOS(AWRK,VEC,NBF,NQMT,NBF)
      IF(PRVEC.and.IPRINTOPT==1)THEN
       WRITE(6,130)                                                 
       CALL PRINTEVECS(VEC,EIG,AWRK,NQMT,NBF)                  
      END IF                                                            
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(.not.NOTOPN)CLOSE(UNIT=20,STATUS='DELETE')           
      NOTOPN = .TRUE.                                                     
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DEALLOCATE(FA,FB,FMO,QFQ,Fsq,WRK,WRK1,AWRK)
      IF(HFDIIS)DEALLOCATE(Ssq,ERR,ADIIS,XDIIS,IPDIIS,BDIIS,IODII)
      IF(HFDAMP.or.HFEXTRAP)THEN
       DEALLOCATE(WRK2,WRK3,FAO0,FAO1,FAO2)
      ENDIF
      RETURN                                                            
!-----------------------------------------------------------------------
   10 FORMAT(/1X,A20,/1X,20(1H-)) 
   20 FORMAT(1X,'EXTRAP=',L1,'  DAMP=',L1,'  DIIS=',L1)  
   30 FORMAT(1X,'DMat Conv =',1P,E8.1,2X,                               &
                'Ener Tol =',1P,E8.1)                                    
   35 FORMAT(1X,'Ener Conv =',1P,E8.1,2X,                               &
                'DMat Tol =',1P,E8.1)                                    
   36 FORMAT(1X,'DIIS Conv =',1P,E8.1,2X,                               &
                'DMat Tol =',1P,E8.1)                                    
   40 FORMAT(/' ITER      TOTAL ENERGY        E CHANGE     ',           &
              'DEN. CHANGE    DIIS ERROR')                               
   50 FORMAT(/' ITER      TOTAL ENERGY        E CHANGE     ',           &
              'DEN. CHANGE')                                             
   60 FORMAT(1X,I3,F20.10,F17.10,2F14.9)                                 
   70 FORMAT(1X,I3,F20.10,F17.10, F14.9)                                 
   80 FORMAT(/10X,17(1H-)/10X,17HDENSITY CONVERGED/10X,17(1H-))          
   90 FORMAT(/10X,14(1H-)/10X,14HDIIS CONVERGED/10X,14(1H-))             
  100 FORMAT(/10X,16(1H-)/10X,16HENERGY CONVERGED/10X,16(1H-))           
  110 FORMAT(/1X,'FINAL ROHF ENERGY IS',F20.10,' AFTER',I4,             &
                 ' ITERATIONS')                                          
  120 FORMAT(/10X,20(1H-)/10X,12HSPIN SZ   = ,F8.3/                     &
              10X,12HS-SQUARED = ,F8.3/10X,20(1H-))                     
  130 FORMAT(/10X,12HEigenvectors/10X,12(1H-))              
!-----------------------------------------------------------------------
      END                                                               

! HSTARRO                                           
      SUBROUTINE HSTARRO(NBF,DAlpha,FA,DBeta,FB,NBFT,XInteg,IXInteg,    &
                         NINTEGt,NREC,IDONTW)     
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                     
      COMMON/INTFIL/NINTMX               
      DOUBLE PRECISION,DIMENSION(NBFT) :: DAlpha,FA,DBeta,FB
      INTEGER(8),DIMENSION(NINTEGt) :: IXInteg
      DOUBLE PRECISION,DIMENSION(NINTEGt) :: XInteg    
      INTEGER(8),ALLOCATABLE,DIMENSION(:) :: IBUF
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: BUF
      INTEGER(8) :: LABEL, I, J, K, L
!-----------------------------------------------------------------------
!     Form the Skeletons of the Fock matrices FA & FB for ROHF
!-----------------------------------------------------------------------
      IF(IDONTW==0)THEN
       REWIND(1)
       ALLOCATE(IBUF(NINTMX),BUF(NINTMX))
      END IF
!      
      CALL VCLR(FA,1,NBFT)                                            
      CALL VCLR(FB,1,NBFT)                                            
      I = 0                                                             
      J = 0                                                             
      K = 0                                                             
      L = 0                                                             
      NXInteg = 0                                                           
      IRECORD = 0 
!-.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.-      
    1 CONTINUE                                                          
      IF(IDONTW==1)THEN        ! IDONTW=1  (IXInteg,XInteg)      
       IRECORD = IRECORD + 1
       IJBUF = (IRECORD-1)*NINTMX
       if(IRECORD<NREC)then               
        NINTEG = NINTMX
        NXInteg = NINTEG
       else
        NINTEG = NINTEGt - (NREC-1)*NINTMX
        NXInteg = -NINTEG
       end if
       IF(NXInteg==0)GO TO 2                                        
       IF(NINTEG>NINTMX)CALL ABRT
!       
       DO M = 1,NINTEG
        LABEL = IXInteg(IJBUF+M)       
        CALL LABELIJKL(LABEL,I,J,K,L) 
        VAL = XInteg(IJBUF+M)        
        NIJ = (I*I-I)/2+J                                                  
        NKL = (K*K-K)/2+L                                                  
        NIK = (I*I-I)/2+K                                                  
        NIL = (I*I-I)/2+L                                                  
        IF(J<K) GO TO 32                                           
        NJK = (J*J-J)/2+K                                                  
        NJL = (J*J-J)/2+L                                                  
        GO TO 42                                                      
   32   NJK = (K*K-K)/2+J                                                  
        IF(J<L) GO TO 52                                           
        NJL = (J*J-J)/2+L                                                  
        GO TO 42                                                      
   52   NJL = (L*L-L)/2+J                                                  
   42   CONTINUE                                                       
        VAL2 = VAL+VAL                                                 
        VAL4 = VAL2+VAL2                                               
        DUM = VAL4*(DAlpha(NKL)+DBeta(NKL))                                   
        FA(NIJ) = FA(NIJ)+DUM                                          
        FB(NIJ) = FB(NIJ)+DUM                                          
        DUM = VAL4*(DAlpha(NIJ)+DBeta(NIJ))                                   
        FA(NKL) = FA(NKL)+DUM                                          
        FB(NKL) = FB(NKL)+DUM                                          
        FA(NIK) = FA(NIK)-VAL2*DAlpha(NJL)                                 
        FB(NIK) = FB(NIK)-VAL2*DBeta(NJL)                                 
        FA(NIL) = FA(NIL)-VAL2*DAlpha(NJK)                                 
        FB(NIL) = FB(NIL)-VAL2*DBeta(NJK)                                 
        FA(NJK) = FA(NJK)-VAL2*DAlpha(NIL)                                 
        FB(NJK) = FB(NJK)-VAL2*DBeta(NIL)                                 
        FA(NJL) = FA(NJL)-VAL2*DAlpha(NIK)                                 
        FB(NJL) = FB(NJL)-VAL2*DBeta(NIK)                                 
       END DO
      ELSE    
       READ(1)NXInteg,IBUF,BUF
       IF(NXInteg==0)GO TO 2                                               
       NINTEG = IABS(NXInteg)                                                
       IF(NINTEG>NINTMX)CALL ABRT 
!       
       DO M = 1,NINTEG                                                 
        LABEL = IBUF(M)
        CALL LABELIJKL(LABEL,I,J,K,L)                                                                                     
        VAL = BUF(M)                                                
        NIJ = (I*I-I)/2+J                                                  
        NKL = (K*K-K)/2+L                                                  
        NIK = (I*I-I)/2+K                                                  
        NIL = (I*I-I)/2+L                                                  
        IF(J<K) GO TO 31                                           
        NJK = (J*J-J)/2+K                                                  
        NJL = (J*J-J)/2+L                                                  
        GO TO 41                                                      
   31   NJK = (K*K-K)/2+J                                                  
        IF(J<L) GO TO 51                                           
        NJL = (J*J-J)/2+L                                                  
        GO TO 41                                                      
   51   NJL = (L*L-L)/2+J                                                  
   41   CONTINUE                                                       
        VAL2 = VAL+VAL                                                 
        VAL4 = VAL2+VAL2                                               
        DUM = VAL4*(DAlpha(NKL)+DBeta(NKL))                                   
        FA(NIJ) = FA(NIJ)+DUM                                          
        FB(NIJ) = FB(NIJ)+DUM                                          
        DUM = VAL4*(DAlpha(NIJ)+DBeta(NIJ))                                   
        FA(NKL) = FA(NKL)+DUM                                          
        FB(NKL) = FB(NKL)+DUM                                          
        FA(NIK) = FA(NIK)-VAL2*DAlpha(NJL)                                 
        FB(NIK) = FB(NIK)-VAL2*DBeta(NJL)                                 
        FA(NIL) = FA(NIL)-VAL2*DAlpha(NJK)                                 
        FB(NIL) = FB(NIL)-VAL2*DBeta(NJK)                                 
        FA(NJK) = FA(NJK)-VAL2*DAlpha(NIL)                                 
        FB(NJK) = FB(NJK)-VAL2*DBeta(NIL)                                 
        FA(NJL) = FA(NJL)-VAL2*DAlpha(NIK)                                 
        FB(NJL) = FB(NJL)-VAL2*DBeta(NIK)                                 
       END DO
      END IF       
!-.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.- -.-
      IF(NXInteg>0)GO TO 1
    2 CONTINUE                                                                
!----------------------------------------------------------------------- 
      IF(IDONTW==0)THEN
       REWIND(1)
       DEALLOCATE(IBUF,BUF)
      END IF
!
      CALL DSCAL(NBFT,0.5d0,FA,1)                                          
      CALL DSCAL(NBFT,0.5d0,FB,1)                                          
      II = 0                                                            
      DO I=1,NBF                                                 
       II = II + I                                                    
       FA(II) = FA(II)+FA(II)                                         
       FB(II) = FB(II)+FB(II)                                         
      END DO
!-----------------------------------------------------------------------
      RETURN                                                            
      END
                                                               
! ROFOCK                                           
      SUBROUTINE ROFOCK(FA,FB,FMO,V,WRK,S,WRK1,NA,NB,L0,L1,L2,L3)  
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)            
      DIMENSION FA(L2),FB(L2),FMO(L2),V(L3),WRK(L1),S(L2),WRK1(L3)         
      LOGICAL JA,JB,JC,JO,JV,IA,IB,IC,IO,IV
!-----------------------------------------------------------------------
!              FORM THE ROHF FOCK MATRIX IN MO BASIS
!-----------------------------------------------------------------------
!                                                                       
!     FA AND FB ARE ALPHA AND BETA FOCK MATRICES IN AO BASIS. 
!     V IS THE CURRENT MO VECTORS. S IS OVERLAP MATRIX.
!     FMO IS THE COMBINED FOCK OPERATOR IN THE MO BASIS.
!                                                                       
!     GUEST AND SAUNDERS FMO (MolPhys 28, 819, 1974) :
!
!                 CLOSED         OPEN         VIRTUAL                   
!                                                                       
!     CLOSED    (FA + FB)/2 |     FB      | (FA + FB)/2                 
!               ---------------------------------------                 
!     OPEN          FB      | (FA + FB)/2 |     FA                      
!               ---------------------------------------                 
!     VIRTUAL   (FA + FB)/2 |     FA      | (FA + FB)/2                  
!
!-----------------------------------------------------------------------
!     Transform Fock Matrices to MO basis 
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL TFTRI(FMO,FA,V,WRK,L0,L1,L1)   ! FMO=V*FA*V                              
      CALL DCOPY(L2,FMO,1,FA,1)           ! FMO -> FA                              
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(NB/=0)THEN
       CALL TFTRI(FMO,FB,V,WRK,L0,L1,L1)  ! FMO=V*FB*V
       CALL DCOPY(L2,FMO,1,FB,1)          ! FMO -> FB                              
       IJ = 0                                                            
       DO J = 1,L0                                                   
        JA = J<=NA                                                   
        JB = J<=NB                                                   
        JC = JA.and.JB                                                 
        JV = .not.JA.and..not.JB                                       
        JO = .not.JC.and..not.JV                                       
        DO I = 1,J                                                 
         IJ = IJ + 1                                                 
         IA = I<=NA                                                
         IB = I<=NB                                                
         IC = IA.and.IB                                              
         IV = .not.IA.and..not.IB                                    
         IO = .not.IC.and..not.IV                                    
         FMO(IJ) = 0.0d0                                                  
         IF(IC.and.JC) FMO(IJ) = (FA(IJ)+FB(IJ))/2.0d0               
         IF(IC.and.JO) FMO(IJ) = FB(IJ)                                  
         IF(IC.and.JV) FMO(IJ) = (FA(IJ)+FB(IJ))/2.0d0
         IF(IO.and.JO) FMO(IJ) = (FA(IJ)+FB(IJ))/2.0d0                 
         IF(IO.and.JV) FMO(IJ) = FA(IJ)                                  
         IF(IV.and.JV) FMO(IJ) = (FA(IJ)+FB(IJ))/2.0d0                 
        END DO
       END DO
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Back-transform to AO Basis FA = SV * FMO * (SV)+             
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -     
      CALL TFTRIB(FA,FMO,S,V,WRK1,WRK,L0,L1,L2,L3)   ! Output: FA                      
!-----------------------------------------------------------------------
      RETURN                                                            
      END                                                               

! SPIN                                             
      SUBROUTINE SPIN(SZ,S2,DAlpha,DBeta,S,D,T,NA,NB,L1,L2)                 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      DIMENSION DAlpha(L2),DBeta(L2),S(L2),D(L2),T(L1)                  
!-----------------------------------------------------------------------
!     Spin Expectation Values
!-----------------------------------------------------------------------
      IJ=0                                                              
      DO J=1,L1                                                     
       DO I=1,L1                                                  
        DUM = 0.0d0                                                  
        DO K = 1,L1                                             
         IK = (I*I-I)/2+K                                             
         IF(K>I) IK=(K*K-K)/2+I                                    
         JK = (J*J-J)/2+K                                             
         IF(K>J) JK = (K*K-K)/2+J                                  
         DUM = DUM+DAlpha(IK)*S(JK)                                   
        END DO
        T(I) = DUM                                                  
       END DO
       DO I = 1,J                                                 
        DUM = 0.0d0                                                  
        DO K = 1,L1                                             
         IK = (I*I-I)/2+K                                             
         IF(K>I) IK=(K*K-K)/2 + I                                  
         DUM = DUM+S(IK)*T(K)                                     
        END DO
        IJ = IJ+1                                                   
        D(IJ) = DUM                                                 
       END DO
      END DO
!                                                                       
      SZ = (NA-NB)/2.0d0
      S2 = SZ*SZ + (NA+NB)/2.0d0  -  TRACEP(DBeta,D,L1)                     
!-----------------------------------------------------------------------      
      RETURN                                                            
      END                                                               

! TFTRIB                                           
      SUBROUTINE TFTRIB(FAO,FMO,S,V,SV,WRK,L0,L1,L2,L3)                 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      DIMENSION FAO(L2),FMO(L2),S(L2),V(L3),SV(L3),WRK(L1)              
!-----------------------------------------------------------------------
!     Transform from MO basis to AO basis: FAO = SV * FMO * (SV)t                                
!-----------------------------------------------------------------------
      CALL MTARBR(S,L1,V,L0,SV,L1,1)        ! SV  = S*V                            
      CALL TRPOSQ(SV,L1)                    ! SV  = (SV)t                                              
      CALL TFTRI(FAO,FMO,SV,WRK,L1,L0,L1)   ! FAO = SV*FMO*(SV)t                           
!-----------------------------------------------------------------------
      RETURN                                                            
      END                                                               

! MTARBR                                           
      SUBROUTINE MTARBR(A,NA,B,MB,AB,NAB,INCA)                          
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      DIMENSION A(*),B(NA,MB),AB(NAB,MB)                                
!-----------------------------------------------------------------------
!     MULTIPLY SYMMETRIC MATRIX A TIMES RECTANGULAR MATRIX B                          
!     AND GET RECTANGULAR MATRIX AB   ( AB = A * B )
!-----------------------------------------------------------------------
      INC=INCA                                                          
!     PROCESS DIAGONAL ELEMENTS OF INPUT MATRIX A                    
      IJ=1-INC                                                          
      DO I=1,NA                                                     
       IJ=IJ+I*INC                                                    
       AIJ=A(IJ)                                                      
       DO K=1,MB                                                  
        AB(I,K)=AIJ*B(I,K)                                          
       END DO
      END DO
      IF(NA==1)RETURN                                                
!     PROCESS OFF-DIAGONAL ELEMENTS OF INPUT MATRIX A                
      IJ=1-INC                                                          
      DO I=2,NA                                                     
       IJ=IJ+INC                                                      
       IM1=I-1                                                        
       DO J=1,IM1                                                 
        IJ=IJ+INC                                                   
        AIJ=A(IJ)                                                   
        IF(AIJ/=0.0d0)THEN
         CALL DAXPY(MB,AIJ,B(I,1),NA,AB(J,1),NAB)                 
         CALL DAXPY(MB,AIJ,B(J,1),NA,AB(I,1),NAB)                 
        END IF
       END DO
      END DO
!-----------------------------------------------------------------------
      RETURN                                                            
      END                                                               

! TRPOSQ                                           
      SUBROUTINE TRPOSQ(A,N)                                            
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      DIMENSION A(N,N)
!-----------------------------------------------------------------------
!     TRANSPOSE SQUARE MATRIX IN PLACE                                  
!-----------------------------------------------------------------------
      DO J = 2,N                                                    
       JMO = J - 1                                                    
       DO I = 1,JMO                                               
        TMP = A(I,J)                                                
        A(I,J) = A(J,I)                                             
        A(J,I) = TMP
       END DO
      END DO                                                
!-----------------------------------------------------------------------
      RETURN                                                            
      END                                                               

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !                                                         
!                       D I I S Interpolation                          !
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !                                                 

! DIISER                                           
      SUBROUTINE DIISER(F,D,S,ERR,WRK,L1,L2)                
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION F(L1,L1),D(L2),S(L1,L1),ERR(L1,L1),WRK(L1)              
!-----------------------------------------------------------------------    
!     ERR = FDS - SDF
!-----------------------------------------------------------------------    
      DO I = 1,L1,5                                            
       IIMAX = MIN(L1,I+4)                                     
       DO II=I,IIMAX                                              
        KL = 0                                                      
        DO L = 1,L1                                             
         LM1 = L-1                                                
         DUM = 0.0d0                                               
         FDUM = F(L,II)                                           
         IF(LM1>0) THEN                                       
          DO K = 1,LM1                                      
           KL = KL+1                                          
           WRK(K) = WRK(K)+D(KL)*FDUM                         
           DUM = DUM+D(KL)*F(K,II)                            
          END DO                                                         
         END IF                                                   
         KL = KL+1                                                
         WRK(L) = DUM+D(KL)*FDUM                                  
        END DO                                                         
        DO J = 1,L1                                          
         DUM = DDOT(L1,WRK,1,S(1,J),1)                         
         IF(ABS(DUM)<1.0D-15)DUM=0.0d0                      
         ERR(II,J) = DUM                                       
        END DO                                                         
       END DO                                                         
      END DO                                                         
!     Substracting (F*D*S)t from itself
      DO J=1,L1                                                    
       DO I=1,J-1                                                 
        DIFF = ERR(I,J)-ERR(J,I)                                      
        ERR(I,J) =  DIFF                                               
        ERR(J,I) = -DIFF                                               
       END DO
       ERR(J,J) = 0.0d0                                                    
      END DO
!-----------------------------------------------------------------------    
      RETURN                                                            
      END 
      
! DIISInter                                             
      SUBROUTINE DIISInter(ITDIIS,Q,FCKA,ERR,WRK,A,X,IPVT,B,IODIIS,SCR, &
                           L1,L2,L3,MAXITRHF,MAXIT2,MAXIO,ERDIIS,NOTOPN)                       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      COMMON/DIISSO/ETHRSH,MAXDII,IRAF      
      DIMENSION Q(L3),FCKA(L2),ERR(L3),WRK(L3),IODIIS(MAXIO),           &
                X(MAXITRHF),IPVT(MAXITRHF),B(MAXIT2),SCR(L1)
      LOGICAL NOTOPN
      DIMENSION A(MAXDII,MAXDII)                        
!-----------------------------------------------------------------------
!             Direct Inversion in the Iterative Procedure
!     PULAY, JCC 3, 556 (1982); CPL 73, 393 (1980), JCP 84, 5728 (1986)     
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     DIIS DATA IS SAVED ON DIRECT FILE -IRAF- IN CYCLIC ORDER,         
!                 1,2,3,...,MAXDII = LAST MAXDII ERROR MATRICES         
!        MAXDII + 1,2,3,...,MAXDII = LAST MAXDII ERROR MATRICES (BETA)  
!      2*MAXDII + 1,2,3,...,MAXDII = LAST MAXDII ALPHA FOCK MATRICES    
!      3*MAXDII + 1,2,3,...,MAXDII = LAST MAXDII BETA  FOCK MATRICES    
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      LRAFEA = 0                                                         
      LRAFEB = MAXDII                                                    
      LRAFFA = MAXDII*2                                                  
      LRAFFB = MAXDII*3                                                  
!                                                                       
!     ----- IF WE ARE ON THE LAST SCF CYCLE, SHUT DIIS DOWN -----       
!                                                                       
      IF(ITDIIS>=MAXITRHF)THEN                                          
       IF(.not.NOTOPN)CLOSE(UNIT=IRAF,STATUS='DELETE')
       NOTOPN = .TRUE.                                                
       RETURN                                                         
      END IF                                                            
!                                                                       
      IF(NOTOPN) THEN                                                   
       NOTOPN = .FALSE.                                               
       IRAF   = 20                                                    
       NDAF20 = MAXIO                                                 
       LEN20  = L3                                                    
       CALL RAOPEN(IRAF,IODIIS,0,NDAF20,LEN20)                     
      END IF                                                            
!                                                                       
!     ----- PUT ERROR MATRIX INTO CONSISTENT O.N.B. -----               
!     PULAY USES S**-1/2, BUT HERE WE USE Q, Q OBEYS Q-DAGGER*S*Q=I     
!     E-ORTH = Q-DAGGER * E * Q, FCKA IS USED AS A SCRATCH -L1- VECTOR  
!                                                                       
      CALL DCOPY(L3,ERR,1,WRK,1)                                        
      CALL TFSQU(ERR,WRK,Q,SCR,L1,L1)                                   
!                                                                       
!     ----- START DIIS PROCEDURE IF ERDIIS < ETHRSH -----               
!                                                                       
      IMAX   = IDAMAX(L3,ERR,1)                                         
      ERDIIS = ABS(ERR(IMAX))                                           
      IF(ERDIIS>ETHRSH  .and.  ITDIIS==1) RETURN                   
!                                                                       
!     ----- SAVE THE CURRENT FOCK MATRIX -----                          
!                                                                       
      LFCKA = LRAFFA + MOD(ITDIIS-1,MAXDII) + 1                         
      LFCKB = LRAFFB + MOD(ITDIIS-1,MAXDII) + 1                         
      CALL RAWRIT(IRAF,IODIIS,FCKA,L2,LFCKA,0)                          
!                                                                       
!     ----- SAVE THE CURRENT ERROR MATRIX -----                         
!                                                                       
      LERR = LRAFEA + MOD(ITDIIS-1,MAXDII) + 1                          
      CALL RAWRIT(IRAF,IODIIS,ERR,L3,LERR,0)                            
      ITDIIS = ITDIIS+1                                                 
!                                                                       
!        FORM LAST ROW (ROW ITDIIS) OF TRIANGULAR MATRIX B,             
!        FIRST DOING THE DIAGONAL AND THE FIRST ELEMENT OF THE ROW.     
!        THE TRANSPOSE IN PULAY'S STEP (2) IS JUST A SIGN CHANGE        
!        BECAUSE ERROR MATRICES ARE ANTISYMMETRIC BY CONSTRUCTION.      
!        THE SIGN IS IGNORED BECAUSE IT DOESN'T AFFECT THE LINEAR       
!        EQUATION'S SOLUTION.                                           
!                                                                       
      BJJ=DDOT(L3,ERR,1,ERR,1)                                          
!                                                                       
      IF(ITDIIS==2) THEN                                              
       CALL VCLR(B,1,MAXIT2)                                          
       B(1)=0.0d0                                                      
       B(2)=-1.0d0                                                      
       B(3)= BJJ                                                      
       RETURN                                                         
      ELSE                                                              
       J1   = (ITDIIS*ITDIIS-ITDIIS)/2 + 1                            
       JJ   = (ITDIIS*ITDIIS+ITDIIS)/2                                
       B(J1)= -1.0d0                                                    
       B(JJ)= BJJ                                                     
      END IF                                                            
!                                                                       
!        THE REST OF THE BIJ'S.                                         
!        NOTE THAT WE ONLY COMPUTE THE NA-2 VALUES NEXT TO THE          
!        DIAGONAL, SO THAT THE OTHER ELEMENTS (EXCEPT THE FIRST         
!        COLUMN) ARE RANDOM VALUES.                                     
!                                                                       
      NA=MIN(ITDIIS,MAXDII)                                             
      NAM1=NA-1                                                         
      IJ=JJ-NA+1                                                        
      DO IX=2,NAM1                                                  
       LERR = LRAFEA + MOD(ITDIIS-NA+IX-2,MAXDII) + 1                 
       CALL RAREAD(IRAF,IODIIS,WRK,L3,LERR,0)                         
       BIJ=DDOT(L3,ERR,1,WRK,1)                                       
       IJ=IJ+1                                                        
       B(IJ)=BIJ                                                      
      END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     DIIS Linear Equations
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    1 CONTINUE                                                          
      CALL DIISEQ(A,NA,B,JJ,ITDIIS)                                     
      CALL VCLR(X,1,NA)                                                 
      X(1)=-1.0d0                                                         
      IERR=0                                                            
      CALL SLVLEQ(A,X,IPVT,NA,NA,0,IERR)                                
      IF(IERR/=0)THEN                                                
       NA = NA - 1                                                    
        WRITE(6,*)                                                      &
        'REDUCING DIIS EQUATION SIZE BY 1 FOR NUMERICAL STABILITY'  
       IF(NA/=1)GO TO 1                                       
       WRITE(6,*)'DIIS EQUATIONS ARE SINGULAR, BOMBING'  
       CALL ABRT                                                      
      END IF                                                            
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Interpolated Fock Matrix: Add C(I)*F(I)                          
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL VCLR(FCKA,1,L2)                                              
      DO IX=2,NA                                                    
       LFCKA = LRAFFA + MOD(ITDIIS-NA+IX-2,MAXDII) + 1                
       CALL RAREAD(IRAF,IODIIS,WRK,L2,LFCKA,0)                        
       CI = X(IX)                                                     
       CALL DAXPY(L2,CI,WRK,1,FCKA,1)                                 
      END DO
!-----------------------------------------------------------------------                                                                       
      RETURN                                                            
      END                                                               

! DIISEQ                                           
      SUBROUTINE DIISEQ(A,NA,B,NTT,NB)                                  
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      DIMENSION A(NA,NA),B(NTT)                                        
!-----------------------------------------------------------------------                                                                       
      NSKP1=NB-NA+1                                                     
      IB=0                                                              
      JA=0                                                              
      DO J=1,NB                                                     
       IA=0                                                           
       DO I=1,J                                                   
        IB=IB+1                                                     
        IF(I==1.or.I>NSKP1) THEN                             
         IF(J==1.or.J>NSKP1) THEN                          
          IA=IA+1                                               
          IF(IA==1) JA=JA+1                                   
          BIB=B(IB)                                             
          A(IA,JA)=BIB                                          
          A(JA,IA)=BIB                                          
         END IF                                                   
        END IF                                                      
       END DO
      END DO
!-----------------------------------------------------------------------
      RETURN                                                            
      END                                                               
            
! RAOPEN                                           
      SUBROUTINE RAOPEN(IRAFX,IORA,LPHYS,NUMREC,LENREC) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                      
      COMMON/DAIOLN/IRECLN
      COMMON/RAIOLN/JRECLN(10),JRECST(10)
      DIMENSION IORA(NUMREC)                                            
!-----------------------------------------------------------------------
!     Record Length
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                     
      IRECLN = 2048                                                 
      IF(LPHYS/=0.AND.LENREC<IRECLN)IRECLN=LENREC
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -             
!     Physical Records Needed
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      JRECST(IRAFX/10) = 1                                                 
      JRECLN(IRAFX/10) = IRECLN                                            
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Open File in UNIT=IRAFX
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      OPEN(UNIT=IRAFX,STATUS='UNKNOWN',ACCESS='DIRECT',                 &
           FORM='UNFORMATTED',RECL=8*IRECLN)                                              
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO I = 1,NUMREC                                               
       IORA(I) = -1                                                   
      END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN                                                            
      END                                                               
            
! RAREAD                                           
      SUBROUTINE RAREAD(IRAFX,IORA,V,LEN,NRECORD,NAVX)                       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                      
      COMMON/DAIOLN/IRECLN
      COMMON/RAIOLN/JRECLN(10),JRECST(10)      
      DIMENSION V(*),IORA(*)                                            
!     Read a logical record from IRAFX
      IF(IRAFX<0)WRITE(6,*) 'BOGUS RAREAD, NAV=',NAVX                 
      IRECLN = JRECLN(IRAFX/10)                                            
      N = IORA(NRECORD)                                                    
      IS = -IRECLN + 1                                                  
      NS = N                                                            
      LENT = LEN                                                        
    1 CONTINUE                                                          
      IS = IS + IRECLN                                               
      IF = IS + LENT - 1                                             
      IF ((IF-IS+1) > IRECLN) IF = IS + IRECLN - 1                
      NSP = NS                                                       
      LENW = IF - IS + 1                                             
      CALL RARD(V(IS),LENW,IRAFX,NSP)                                 
      LENT = LENT - IRECLN                                           
      NS = NS + 1                                                    
      N = NS                                                         
      IF(LENT>=1)GO TO 1                                        
      RETURN                                                            
      END                                                               

! RARD                                             
      SUBROUTINE RARD(V,LEN,IRAF,NS)                                    
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION V(LEN)                                                  
      READ(UNIT=IRAF,REC=NS) V                 
      RETURN                                                            
      END                                                               

! RAWRIT                                           
      SUBROUTINE RAWRIT(IRAFX,IORA,V,LEN,NRECORD,NAVM)                      
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                      
      COMMON/DAIOLN/IRECLN
      COMMON/RAIOLN/JRECLN(10),JRECST(10)      
      DIMENSION V(*),IORA(*)                                            
!     Write a logical record on IRAFX
      IF(IRAFX<0) WRITE(6,*) 'BOGUS RAWRIT, NAVM=',NAVM               
      IRECST=JRECST(IRAFX/10)                                            
      IRECLN=JRECLN(IRAFX/10)                                            
      N = IORA(NRECORD)                                                    
      IF(N<=0) THEN                                                   
       IORA(NRECORD) = IRECST                                            
       IRECST = IRECST + (LEN-1)/IRECLN + 1                           
       N = IORA(NRECORD)                                                 
      END IF                                                            
      JRECST(IRAFX/10)=IRECST                                            
      IST = -IRECLN + 1                                                 
      NS = N                                                            
      LENT = LEN                                                        
    1 CONTINUE                                                          
         IST = IST + IRECLN                                             
         IF = IST + LENT - 1                                            
         IF ((IF-IST+1) > IRECLN) IF = IST+IRECLN-1                  
         NSP = NS                                                       
         LENW = IF - IST + 1                                            
         CALL RAWRT(V(IST),LENW,IRAFX,NSP)                               
         LENT = LENT - IRECLN                                           
         NS = NS + 1                                                    
         N = NS                                                         
      IF (LENT >= 1) GO TO 1                                        
      RETURN                                                            
      END                                                               

! RAWRT                                            
      SUBROUTINE RAWRT(V,LEN,IRAF,NS)                                   
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION V(LEN)                                                  
      WRITE(UNIT=IRAF,REC=NS) V                                       
      RETURN                                                            
      END                                                               

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ! 
!                       Damping and Extrapolation                      !
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ! 
      
! DAMPD                                            
      SUBROUTINE DAMPD(DE,DEP,DEAVG,DAMP,ACURCY,DIFF,DIFFP,DMPTLC)      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      PARAMETER (FAC=16.0D0)              
!-----------------------------------------------------------------------
      DAMPO = DAMP                                                      
      ETEST = ACURCY*ACURCY                                             
      IF(ABS(DE) < ETEST .and. ABS(DEP) < ETEST) GO TO 300       
      IF(ABS(DE) < ETEST) GO TO 320                                 
      IF(ABS(DEP) < ETEST) GO TO 340                                
      IF((DIFFP-DIFF) < 0.0d0) GO TO 100                             
      IF( ABS(DE) >= ACURCY .or. DE > 0.0d0) GO TO 100             
!                                                                       
      DAMP = DAMP/FAC                                                   
      GO TO 280                                                         
!                                                                       
  100 CONTINUE                                                          
      IF ( DE > 0.0d0) GO TO 200                                      
      IF (DEP > 0.0d0) GO TO 180                                      
      IF ( DE > DEP) GO TO 140                                       
!                                                                       
!     ----- DE < 0. , DEP < 0. , DE < DEP -----                         
!                                                                       
      IF ( ABS(DE) < 2.0d0*DEAVG) GO TO 120                            
      DAMP = FAC*MAX(DAMP,DEAVG)                                        
      GO TO 280                                                         
!                                                                       
  120 IF ( ABS(DE) > 0.5d0*DEAVG) GO TO 280                            
      DAMP = DAMP/FAC                                                   
      GO TO 280                                                         
!                                                                       
  140 CONTINUE                                                          
!                                                                       
!     ----- DE < 0. , DEP < 0. , DE > DEP -----                         
!                                                                       
      IF (DE > 0.25d0*DEP) GO TO 160                                   
      DAMP = (DE/DEP)**2*MAX(DAMP,DEAVG)                                
      GO TO 280                                                         
!                                                                       
  160 DAMP = DAMP/FAC                                                   
      GO TO 280                                                         
!                                                                       
  180 CONTINUE                                                          
!                                                                       
!     ----- DE < 0. , DEP > 0. -----                                    
!                                                                       
      DAMP = 4.0d0*MAX(DAMP,DEAVG)                                       
      IF (-DE > DEAVG) DAMP = DAMP*FAC                               
      IF (-DE+DEP >= DEAVG) GO TO 280                                 
      DAMP = DAMP/FAC                                                   
      GO TO 280                                                         
!                                                                       
  200 CONTINUE                                                          
      IF (DEP > 0.0d0) GO TO 220                                      
!                                                                       
!     ----- DE > 0. , DEP < 0. -----                                    
!                                                                       
      DAMP = 4.0d0*MAX(DAMP,DEAVG)                                       
      IF (DE > 0.5d0*DEAVG) DAMP = DAMP*FAC                            
      IF (DE-DEP >= 0.2d0*DEAVG) GO TO 280                              
      DAMP = DAMP/FAC                                                   
      GO TO 280                                                         
!                                                                       
  220 CONTINUE                                                          
!                                                                       
!     ----- DE > 0. , DEP > 0. -----                                    
!                                                                       
      DAMP = 4.0d0*MAX(DAMP,DEAVG)                                       
      IF (DE < 4.0d0*DEP) GO TO 240                                   
      DAMP = FAC*MAX(DAMP,DEAVG)                                        
      GO TO 280                                                         
!                                                                       
  240 IF (DE > 0.25d0*DEP) GO TO 260                                   
      DAMP = DAMP/FAC                                                   
      GO TO 280                                                         
!                                                                       
  260 DAMP = (DE/DEP)**2*MAX(DAMP,DEAVG)                                
  280 CONTINUE                                                          
!                                                                       
!     ----- IF THE DENSITY CONVERGENCE WORSENED - MAKE SURE             
!           THAT THE DAMPING CAN'T DECREASE -----                       
!                                                                       
      IF ((DIFFP-DIFF) < 0.0d0) DAMP = MAX(DAMP,DAMPO)                
      GO TO 360                                                         
!                                                                       
  300 CONTINUE                                                          
!                                                                       
!        DE < ETEST AND DEP < ETEST                                     
      DAMP = DAMP/FAC                                                   
      GO TO 360                                                         
!                                                                       
  320 CONTINUE                                                          
!        DE < ETEST  DEP > ETEST                                        
      DAMP = DAMP/FAC                                                   
      GO TO 360                                                         
!                                                                       
  340 CONTINUE                                                          
!        DEP < ETEST  DE > ETEST                                        
      DAMP = DAMPO                                                      
      IF (DE > 0.0d0) DAMP = MAX(2.0d0*DAMP,DMPTLC)                     
!                                                                       
  360 CONTINUE 
!-----------------------------------------------------------------------                                                         
      RETURN                                                            
      END                                                               

! EXTRAPOL                                           
      SUBROUTINE EXTRAPOL(DE,DAMP,DAMP0,FA0,FA1,FA2,FA3,FAO0,FAO1,FAO2, &
                          NBF,NBFT,ITERV,NCALL,ITYP,ITER,ICALP,ICBET)                                     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      LOGICAL HFDAMP,HFEXTRAP,HFDIIS
      COMMON/INPNOF_HFCONVTECH/HFDAMP,HFEXTRAP,HFDIIS
      COMMON/ACONV/RRSHFT,EXTTOL,DMPTOL
      DIMENSION FA0(NBFT),FA1(NBFT),FA2(NBFT),FA3(NBFT)
      DIMENSION FAO0(NBFT),FAO1(NBFT),FAO2(NBFT)
      LOGICAL EXTPRE,DAMPRE                         
      SAVE DAMPRE,EXTPRE                                                
      DATA DAMPRE,EXTPRE/.FALSE.,.FALSE./                               
!-----------------------------------------------------------------------
      ICOUNT = ICALP                                                    
      IF(NCALL==2)ICOUNT = ICBET                                  
      IF(ITER>1)GO TO 140                                        
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     ITER = 1                                              
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO I = 1,NBFT                                                   
       FA1(I) = FA0(I)                                                     
       FA2(I) = FA0(I)                                                     
      END DO
      GO TO 660                                                         
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Previous Fock Matrices FAO0, FAO1, FAO2 -> FA1, FA2, FA3
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  140 CONTINUE
      FA1 = FAO0
      FA2 = FAO1
      FA3 = FAO2
      IF(ITER>2)GO TO 160                                        
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     ITER = 2
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF (HFDAMP) GO TO 320                                              
      GO TO 420                                      
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     ITER > 2
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  160 IF (DAMPRE) GO TO 200                                             
      IF ( .not. HFDAMP) GO TO 420                                       
      IF (ABS(DE) > EXTTOL) GO TO 180                                
      IF (DE > 0.0d0 .and. ICOUNT > 4) GO TO 180              
      GO TO 420                                                         
!                                                                       
  180 DMPTOL = DMPTOL/100.0d0                                            
      EXTTOL = EXTTOL/100.0d0                                            
      GO TO 300                                                         
!                                                                       
  200 IF (ABS(DE) > DMPTOL) GO TO 300                                
      IF (DE > 0.0d0) GO TO 300                                       
      IF (DAMP > 0.01d0) GO TO 300                                   
      GO TO 420                                                         
!                                                                       
      IF ( .not. EXTPRE) GO TO 260                                      
      IF (ABS(DE) > EXTTOL) GO TO 240                                
      IF (DE > 0.0d0 .and. ICOUNT > 4) GO TO 240              
      GO TO 420                                                         
!                                                                       
  240 EXTTOL = EXTTOL/100.0d0                                            
      DMPTOL = DMPTOL/100.0d0                                            
      RRSHFT = 0.8d0                                                    
      ITERV = 0                                                         
      IF (HFDAMP) GO TO 300                                              
      GO TO 620                                                         
!                                                                       
  260 IF (ABS(DE) > DMPTOL) GO TO 280                                
      IF (DE > 0.0d0) GO TO 280                                       
      IF (HFDAMP .and. DAMP > 0.01d0) GO TO 280                       
      IF ( .not. HFEXTRAP) GO TO 280                                      
      IF (RRSHFT >= 0.4d0) GO TO 280                                 
      IF (ITERV == 0) GO TO 280                                       
      GO TO 420                                                         
!                                                                       
  280 IF (HFDAMP) GO TO 300                                              
      IF (ITERV < 2) GO TO 620                                       
      IF ( .not. HFEXTRAP) GO TO 620                                      
      GO TO 440                                                         
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Davidson's Damping
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  300 IF (ICOUNT < 4) GO TO 320                                 
      GO TO 360                    
!                                                                       
  320 DO I = 1,NBFT                                                   
       FA0(I) = (FA0(I)+DAMP*FA1(I))/(1.0d0+DAMP)                             
      END DO
      GO TO 400                                                         
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     ITER > 2: Damping
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  360 CUTDAMP = 0.5d0-DAMP0
      DO I = 1,NBFT                                                   
       FA0(I) = (FA0(I)+DAMP*FA1(I))/(1.0d0+DAMP)                             
      END DO
      IF ( .not. HFEXTRAP .or. CUTDAMP < 0.0d0) GO TO 400
      DAMPRE = .TRUE.                                                   
      EXTPRE = .TRUE.                                                   
      GO TO 460                                                         
!                                                                       
  400 IF (NCALL /= 1) GO TO 660                                       
      DAMPRE = .TRUE.                                                   
      EXTPRE = .FALSE.                                                  
      GO TO 660                                                         
!                                                                       
  420 IF ( .not. HFEXTRAP) GO TO 620                                      
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Pople's Extrapolation
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF (NCALL /= 1) GO TO 460                                       
      DAMPRE = .FALSE.                                                  
      EXTPRE = .TRUE.                                                   
  440 DAMP = 0.0d0                                                       
  460 CONTINUE
      FAO0 = FA0
      FAO1 = FA1
      FAO2 = FA2
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Extrapolation
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO I = 1,NBFT                                                   
       FA3(I) = FA2(I) - FA3(I)                                               
       FA2(I) = FA1(I) - FA2(I)                                               
       FA1(I) = FA0(I) - FA1(I)                                               
      END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Return after Extrapolation
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(ICOUNT < 4 .or. ITER < 4) GO TO 680                
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Displacements
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(ITYP==2)GO TO 500                                        
      SP11 = TRACEP(FA1,FA1,NBF)                                           
      SP12 = TRACEP(FA2,FA1,NBF)                                           
      SP13 = TRACEP(FA3,FA1,NBF)                                           
      SP22 = TRACEP(FA2,FA2,NBF)                                           
      SP23 = TRACEP(FA3,FA2,NBF)                                           
      SP33 = TRACEP(FA3,FA3,NBF)                                           
      GO TO 520                                                         
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  500 CONTINUE                                                          
      SP11 = DDOT(NBFT,FA1,1,FA1,1)                                         
      SP12 = DDOT(NBFT,FA2,1,FA1,1)                                         
      SP13 = DDOT(NBFT,FA3,1,FA1,1)                                         
      SP22 = DDOT(NBFT,FA2,1,FA2,1)                                         
      SP23 = DDOT(NBFT,FA3,1,FA2,1)                                         
      SP33 = DDOT(NBFT,FA3,1,FA3,1)                                         
  520 CONTINUE                                                          
      DP1 = SQRT(SP11)                                                  
      DP2 = SQRT(SP22)                                                  
      DP3 = SQRT(SP33)                                                  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(DP1<1.0D-10) GO TO 680                                      
      IF(DP2<1.0D-10) GO TO 680                                      
      IF(DP3<1.0D-10) GO TO 680                                      
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Find angle PHI between successive displacements
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      COSPHI = SP12/(DP1*DP2)                                           
!
      Z = SP11*SP22-SP12*SP12                                           
      IF ( ABS(Z) < 1.0D-17 ) GO TO 680                                 
      X = (SP13*SP22-SP12*SP23)/Z                                       
      Y = (SP23*SP11-SP12*SP13)/Z                                       
      COSPSI = SQRT(X*X*SP11+Y*Y*SP22+2.0d0*X*Y*SP12)/DP3                 
!                                                                       
      IF(COSPSI<=1.0D-07)GO TO 680                                    
      IF(DAMP>0.01d0)GO TO 680                                   
!                                                                       
      Y = -Y/X                                                          
      X = 1.0d0/X                                                         
!                                                                       
      XY = Y*Y + 4.0d0*X                                                   
      IF (XY < 0.0d0) GO TO 680                                       
      XY = ABS(Y)+SQRT(XY)                                              
      IF (XY <= 1.9d0) GO TO 560                                       
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     If 4-point extrapolation is not possible, try 3-point extrap.
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF (ABS(COSPHI) <= 0.99d0) GO TO 680                              
      X = DP1/(DP2*COSPHI-DP1)                                          
      DO I = 1,NBFT                                                   
       FA0(I) = FA0(I) + X*FA1(I)                                             
      END DO
      GO TO 600                                                         
!                                                                       
  560 XXX = X/(1.0d0-X-Y)                                                 
      YYY = (X+Y)/(1.0d0-X-Y)                                             
      DO I = 1,NBFT                                                   
       FA0(I) = FA0(I) + XXX*FA2(I)+YYY*FA1(I)                                 
      END DO
!                                                                       
  600 CONTINUE                                                          
      IF(NCALL==1)ICALP=0                                            
      IF(NCALL==2)ICBET=0
      FAO0 = FA0
      GO TO 680                                                         
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     No Damping or Extrapolation
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  620 IF (NCALL /= 1) GO TO 640                                       
      DAMPRE = .FALSE.                                                  
      EXTPRE = .FALSE.                                                  
  640 DAMP = 0.0d0                                                       
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Save Modified Fockians in FAO-matrices
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                     
  660 CONTINUE
      FAO0 = FA0
      FAO1 = FA1
      FAO2 = FA2
!----------------------------------------------------------------------- 
  680 RETURN
      END 

!----------------------------------------------------------------------!
!                                                                      !
!   HFIDr Iterative Diagonalization method for the RHF using           !
!         the Iterative Diagonalization method (IRHF=3) proposed       !
!         in the J. Comp. Chem. 131, 021102, 2009.                     !
!                                                                      !
!----------------------------------------------------------------------!

! HFIDr
      SUBROUTINE HFIDr(AHCORE,IJKL,XIJKL,XIJKaux,CHF,EiHF,USER,IPRINTOPT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL HighSpin,HFDAMP,HFEXTRAP,HFDIIS
      COMMON/INPNOF_HFID/KOOPMANS
      COMMON/INPNOF_GENERALINF/ICOEF,ISOFTMAX,IORBOPT,MAXIT,MAXIT21
      COMMON/INPNOF_NTHRESH/NTHRESHL,NTHRESHE,NTHRESHEC,NTHRESHEN
      COMMON/INPNOF_THRESH/THRESHL,THRESHE,THRESHEC,THRESHEN
      COMMON/INPNOF_ORBSPACE0/NO1,NDOC,NCO,NCWO,NVIR,NAC,NO0
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
      COMMON/INPNOF_HFCONVTECH/HFDAMP,HFEXTRAP,HFDIIS
      LOGICAL EFIELDL,RESTART,CONVGDELAG,SMCD,ERIACTIVATED
      COMMON/INP_EFIELDL/EFX,EFY,EFZ,EFIELDL
      COMMON/INPNOF_RSTRT/RESTART,INPUTGAMMA,INPUTC,INPUTFMIUG,INPUTCXYZ
      COMMON/CONVERGENCE/DUMEL,PCONV,CONVGDELAG
      COMMON/ERITYPE/IERITYP,IRITYP,IGEN,ISTAR,MIXSTATE,SMCD
      COMMON/ERIACT/ERIACTIVATED,NIJKaux,NINTCRaux,NSTOREaux,IAUXDIM
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_COEFOPT/MAXLOOP
      COMMON/INPNOF_PRINT/NPRINT,IWRITEC,IMULPOP,IAIMPAC,IFCHK
      COMMON/INPFILE_NIJKL/NINTMX,NIJKL,NINTCR,NSTORE
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      COMMON/CONVERGESUM/SUMDIF,SUMDIF_OLD
      COMMON/ELPROP/IEMOM
      COMMON/EHFEN/EHF,EN
      COMMON/ENERGIAS/EELEC,EELEC_OLD,DIF_EELEC,EELEC_MIN
      COMMON/PUNTEROSUSER/N1,N2,N3,N4,N5,N6,N7,N8,N9,N10,N11,N12,N13,   &
                          N14,N15,N16,N17,N18,N19,N20,N21,N22,N23,N24,  &
                          N25,N26,N27,N28,N29,N30,N31,N32,N33,N34,N35,  &
                          N36,N37,N38,N39,N40,N41,N42,N43,N44,N45,N46,  &
                          N47,N48,N49,N50,N51,NUSER
      COMMON/USEHUBBARD/IHUB
!
      INTEGER :: IPRINTOPT
      INTEGER(8),DIMENSION(NSTORE)::IJKL
      DOUBLE PRECISION,DIMENSION(NSTORE)::XIJKL
      DOUBLE PRECISION,DIMENSION(NSTOREaux)::XIJKaux
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::AHCORE,CHF
      DOUBLE PRECISION,DIMENSION(NBF)::EiHF
      DOUBLE PRECISION,DIMENSION(NUSER)::USER
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::RO10,FMIUG0
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::EVA,TEMP,CFM
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::DMP,DM,WRK1,WRK2,WRK3
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::FAO,FAO0,FAO1,FAO2
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::CJ12HF,CK12HF,DEN
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::ELAG,G,FMIUG,W
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::CHFNEW,BFM
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:)::FK
!-----------------------------------------------------------------------
!     Initial Values
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      INPUTFMIUG_ORI = INPUTFMIUG
      INPUTFMIUG = 0
      NO1_ORI = NO1
      NO1 = NB
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(RO10(NBF5),CJ12HF(NBF5,NBF5),CK12HF(NBF5,NBF5))
      RO10 = 0.0d0
      DO i=1,NB
       RO10(i) = 1.0d0
      ENDDO
      IF(NSOC>0)THEN
       DO i=NB+1,NA
        RO10(i) = 0.5d0
       ENDDO
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Density Matrix
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(DEN(NBF,NBF))
      CALL DENMATr(DEN,CHF,RO10,NBF,1,NBF5)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(IPRINTOPT==1)WRITE(6,1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Calculate Initial HF Electronic Energy
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(ELAG(NBF,NBF),G(NBF,NBF5))
      DO j=1,NBF5
       DO i=1,NBF5
        CJ12HF(j,i) = 2.0d0*RO10(j)*RO10(i)
        CK12HF(j,i) = RO10(j)*RO10(i)
       ENDDO
      ENDDO
      if(MSpin==0.and.NSOC>1)then
       DO j=NB+1,NA
        DO i=NB+1,NA
         CK12HF(j,i) = 2.0d0*RO10(j)*RO10(i)
        ENDDO
       ENDDO
      end if
      CALL ELAG1r(AHCORE,IJKL,XIJKL,XIJKaux,USER(N7),CHF,RO10,CJ12HF,   &
                  CK12HF,ELAG,G)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Calculate HF Electronic Energy (EELEC=EHF)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      EELEC = TRACE(DEN,AHCORE,NBF)
      DO iq=1,NBF5
       EELEC = EELEC + ELAG(iq,iq)
      END DO
!     Include Nuclear Dipoles if electric field =/ 0
      IF(EFIELDL)THEN
       EELEC = EELEC - EFX*USER(N11) - EFY*USER(N11+1) - EFZ*USER(N11+2)
      END IF
      EHF = EELEC
      DELEHF = EHF
      CALL PCONVE(ELAG,DUMEL,MAXI,MAXJ,SUMDIF)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Itermediate Output of the external iteration
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(IPRINTOPT==1)THEN
       WRITE(6,5)1,EHF,EHF+EN,0.0,DUMEL
      ENDIF
!-----------------------------------------------------------------------
!                       START SCF-ITERATION CYCLE
!-----------------------------------------------------------------------
      ALLOCATE (FMIUG(NBF,NBF),W(NBF,NBF),EVA(NBF),TEMP(NBF))
      ALLOCATE (CHFNEW(NBF,NBF),FMIUG0(NBF))
      IF(HFDIIS)THEN
       ALLOCATE(CFM(MAXLOOP+1),BFM(MAXLOOP+1,MAXLOOP+1))
       ALLOCATE(FK(MAXLOOP,NBF,NBF))
      ENDIF
      IF(HFDAMP.or.HFEXTRAP)THEN
       ALLOCATE(DMP(NBFT),DM(NBFT),WRK1(NSQ),WRK2(NSQ),WRK3(NSQ))
       ALLOCATE(FAO(NBFT),FAO0(NBFT),FAO1(NBFT),FAO2(NBFT))
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO LOOPEXT=1,MAXIT
       IF(LOOPEXT==1)THEN
        MAXLP = 1
       ELSE
        MAXLP = MAXLOOP
       ENDIF
       EHF_OLD = EHF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Check for energy convergent solution and Output
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       IF( ABS(DELEHF) < THRESHE .and. DUMEL < THRESHL )THEN
        INPUTFMIUG = INPUTFMIUG_ORI
        NO1 = NO1_ORI
        IF(IPRINTOPT==0)RETURN
!       One-particle HF energies
        IF(IHUB==0)THEN
         DO i=1,NBF
          EiHF(i) = ELAG(i,i)
         ENDDO
         CALL PRINTEiHF(EiHF,NA,NBF)
        ENDIF
        IF(KOOPMANS==1 .and. MSpin==0)THEN
         CALL DIAGELAGHF(ELAG,CHF,RO10,EiHF,CHFNEW,TEMP)
        ENDIF
        if(NTHRESHE<=10)then
         WRITE(6,2)EHF
         WRITE(6,3)EHF+EN
        else
         WRITE(6,21)EHF
         WRITE(6,31)EHF+EN
        end if
        IF(EFIELDL)WRITE(6,4)EFX,EFY,EFZ
!       HF Electrostatic Moments
        IF(1<=IEMOM.and.IEMOM<=3)CALL DQOHF(IEMOM,USER,RO10)
        RETURN
       ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       ILOOP = 0
       DELE = 1.0d20
       IDIIS = 0
       DAMP   = 0.0d0
       DAMP0  = 0.0d0
       IF(HFDAMP)DAMP = 1.0d0
       DIFFP = 0.0d0
       DIFF  = 0.0d0
       ITERV = 0
       RRSHFT = 0.0d0
       EXTTOL = 1.0D-03
       DMPTOL = 1.0D-04
       ICALP  = 0
       ICBET  = 0
!-----------------------------------------------------------------------
       DO WHILE( ILOOP<MAXLP .and. DABS(DELE)>THRESHE )
        ILOOP=ILOOP+1
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       Fock Matrix (FMIUG)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        CALL FORMFMIUGHF(FMIUG,ELAG,FMIUG0,LOOPEXT)
        IF(HFDAMP.or.HFEXTRAP)CALL SQUARETRIAN(FMIUG,FAO,NBF,NBFT)
        IF(HFDIIS.and.DUMEL<1.0d-3)THEN
         CALL FFMIUG_DIIS(NBF,FMIUG,CFM,BFM,FK,IDIIS)
         IF(IDIIS>0)THEN
          IF(HFDAMP)DAMP = 1.0d0
          RRSHFT = 0.0d0
         END IF
        ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       Damping and Extrapolation
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        IF(ILOOP==2)DEAVG = ABS(DELE)
        IF(IDIIS==0)THEN
         IF(HFDAMP.and.ILOOP>2)THEN
          DEAVG = (ABS(DELE)+ABS(DELE0)+0.2d0*DEAVG)/2.2D0
          DAMP0 = DAMP
          CALL DAMPD(DELE,DELE0,DEAVG,DAMP,1.0D-05,DIFF,DIFFP,1.0D-02)
         END IF
         IF(DAMP<0.0d0)DAMP = 0.0d0
         IF(HFDAMP.or.HFEXTRAP)THEN
          CALL EXTRAPOL(DELE,DAMP,DAMP0,FAO,WRK1,WRK2,WRK3,FAO0,FAO1,   &
                        FAO2,NBF,NBFT,ITERV,1,1,ILOOP,ICALP,ICBET)
          CALL TRIANSQUARE(FMIUG,FAO,NBF,NBFT)
         ENDIF
        END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       DIAGONALIZE SQUARE MATRIX (FMIUG) FOR REAL SYMMETRIC CASE
!       W - EIGENVECTORS, EVA-EIGENVALUES IN ALGEBRAIC DESCENDING ORDER
!       HOUSEHOLDER METHOD
!       NOTE: ONLY LOWER TRIANGLE IS USED + THIS IS DESTROYED !!!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        CALL DIAG(NBF,FMIUG,W,EVA,TEMP)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       Move EVA -> FMIUG0
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        FMIUG0 = EVA
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       New Coefficients (CHFNEW=CHF*W), Move CHFNEW -> CHF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        CALL COEFW(NBF,NBF,CHFNEW,CHF,W)
        CHF = CHFNEW
        IF(HFDAMP.or.HFEXTRAP)CALL DCOPY(NBFT,DM,1,DMP,1)   ! DM -> DMP
        CALL DENMATr(DEN,CHF,RO10,NBF,1,NBF5)  ! Density Matrix
        IF(HFDAMP.or.HFEXTRAP)THEN
         DIFFP = DIFF
         CALL SQUARETRIAN(DEN,DM,NBF,NBFT)
         CALL DDIFF(DMP,DM,NBFT,DIFF)
         ICALP = ICALP+1
         ICBET = ICBET+1
        ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       Lagrangian Multipliers (ELAG) and one-energies (E)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        CALL ELAG1r(AHCORE,IJKL,XIJKL,XIJKaux,USER(N7),CHF,RO10,CJ12HF, &
                    CK12HF,ELAG,G)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       Calculate HF Electronic Energy (EELEC=EHF)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        EHF0 = EHF
        DELE0 = DELE
!       Eelec
        EELEC = TRACE(DEN,AHCORE,NBF)
        DO iq=1,NBF5
         EELEC = EELEC + ELAG(iq,iq)
        END DO
!       Include Nuclear Dipoles if electric field =/ 0
        IF(EFIELDL)THEN
         EELEC = EELEC - EFX*USER(N11)-EFY*USER(N11+1)-EFZ*USER(N11+2)
        END IF
        EHF = EELEC
!       Differences
        DELE = EHF - EHF0
        DELEHF = EHF - EHF_OLD
        CALL PCONVE(ELAG,DUMEL,MAXI,MAXJ,SUMDIF)
!       Intermediate Output at each internal interation (Nprint=2)
        IF(NPRINT==2.and.IPRINTOPT==1)THEN
         WRITE(6,5)ILOOP,EHF,EHF+EN,DELE,DUMEL
        ENDIF
!-----------------------------------------------------------------------
!                       LOOP-END OF SCF-ITERATION
!-----------------------------------------------------------------------
       ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Intermediate Output at each external iteration
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       IF(NPRINT<2.and.IPRINTOPT==1.and.LOOPEXT>1)THEN
        WRITE(6,5)LOOPEXT,EHF,EHF+EN,DELEHF,DUMEL
       ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ENDDO
!-----------------------------------------------------------------------
!     FORMAT STATEMENTS
!-----------------------------------------------------------------------
    1 FORMAT(//,1X,' HARTREE-FOCK ',/, 1X,'==============',//,          &
                2X,'ITER',5X,'ELECTRONIC ENERGY',5X,'TOTAL ENERGY',     &
                4X,'ENERGY CONVERGENCY',4X,'MAX MUL-LAG DIFF',/)
    2 FORMAT(/,4X,'ELECTRONIC HF ENERGY =',F20.10)
   21 FORMAT(/,4X,'ELECTRONIC HF ENERGY =',F25.15)
    3 FORMAT(/,8X,' HF TOTAL ENERGY =',F20.10)
   31 FORMAT(/,8X,' HF TOTAL ENERGY =',F25.15)
    4 FORMAT(/,6X,'ELECTRIC FIELD (',D8.1,',',D8.1,',',D8.1,')')
    5 FORMAT(2X,I3,'.',3X,F17.10,2X,F17.10,6X,F13.8,8X,F11.6)
!-----------------------------------------------------------------------
      DEALLOCATE(RO10,CJ12HF,CK12HF)
      DEALLOCATE(FMIUG,W,EVA,TEMP,CHFNEW,DEN,ELAG,G,FMIUG0)
      IF(HFDIIS)DEALLOCATE(CFM,BFM,FK)
      IF(HFDAMP.or.HFEXTRAP)THEN
       DEALLOCATE(DMP,DM,WRK1,WRK2,WRK3,FAO,FAO0,FAO1,FAO2)
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      INPUTFMIUG = INPUTFMIUG_ORI
      NO1 = NO1_ORI
!-----------------------------------------------------------------------
      RETURN
      END

! FORMFMIUGHF
      SUBROUTINE FORMFMIUGHF(FMIUG,ELAG,FMIUG0,ITCALL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      DOUBLE PRECISION,DIMENSION(NBF)::FMIUG0
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::FMIUG,ELAG
!-----------------------------------------------------------------------
!     Fock Matrix (FMIUG)
!-----------------------------------------------------------------------
      IF(ITCALL==1)THEN  ! only for itcall==1
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       DO IQ=1,NBF
        DO JQ=1,IQ-1
         FMIUG(IQ,JQ)=(ELAG(IQ,JQ)+ELAG(JQ,IQ))/2.0      ! Nondiagonal
         FMIUG(JQ,IQ)=FMIUG(IQ,JQ)                       ! Fji=Fij
        ENDDO
        FMIUG(IQ,IQ)=ELAG(IQ,IQ)                         ! Diagonal
       ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ELSE
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       DO IQ=1,NBF
        DO JQ=1,IQ-1
         FMIUG(IQ,JQ)=ELAG(IQ,JQ)-ELAG(JQ,IQ)            ! Nondiagonal
         FMIUG(JQ,IQ)=FMIUG(IQ,JQ)                       ! Fji=Fij
        ENDDO
        FMIUG(IQ,IQ)=FMIUG0(IQ)                          ! Diagonal
       ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ENDIF
!-----------------------------------------------------------------------
      RETURN
      END

! PRINTEiHF
      SUBROUTINE PRINTEiHF(EiHF,NA,NBF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION,DIMENSION(NBF)::EiHF
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::E1
!-----------------------------------------------------------------------
!     ORDERING ENERGIES
!-----------------------------------------------------------------------
      ALLOCATE(E1(NBF))
      DO I=1,NBF
       E1(I)=EiHF(I)
      ENDDO
      DO IM=1,NBF
       DUM1= E1(IM)
       MIN = IM
       DO I=IM,NBF
        IF(E1(I)<DUM1)THEN
         DUM1= E1(I)
         MIN = I
        ENDIF
       ENDDO
       IF(MIN/=IM)THEN
        DUM1   = E1(IM)
        E1(IM) = E1(MIN)
        E1(MIN)= DUM1
       ENDIF
      ENDDO
!-----------------------------------------------------------------------
!     WRITE ENERGIES ON THE OUTPUT FILE
!-----------------------------------------------------------------------
      WRITE(6,100)
      DO I=1,NA
       WRITE(6,101)I,E1(I),E1(I)*27.21138386
      ENDDO
!-----------------------------------------------------------------------
  100 FORMAT(/2X,'-------------',/2X,' HF Energies ',                   &
             /2X,'-------------',//19X,'(aU)',14X,'(eV)',20X)
  101 FORMAT(2X,I4,4X,F15.6,4X,F15.6)
!-----------------------------------------------------------------------
      DEALLOCATE(E1)
      RETURN
      END

! DIAGELAGHF
      SUBROUTINE DIAGELAGHF(ELAG,COEF,RO,ELAGN,COEFN,RON)
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
       IF(RON(I)>1.0d-1)THEN
        WRITE(6,2)I,-ELAGN(I),-ELAGN(I)*27.21138386,RON(I)
       ENDIF
      ENDDO
      CALL DMATMAX(DENMAT,NBF,MAXI,MAXJ,DUM)
      WRITE(6,3)DUM,MAXI,MAXJ
!-----------------------------------------------------------------------
!     Coefficients of Canonical Orbitals (COEFN=COEF*W)
!-----------------------------------------------------------------------
      COEFN = MATMUL(COEF,W)
      ICANHF=0
      IF(ICANHF==1)THEN
       WRITE(6,4)
       CALL PRINTVERO(6,COEFN,ELAGN,RON,NBF,NBF5)
      ENDIF
!-----------------------------------------------------------------------
    1 FORMAT(/2X,42('-'),/3X,'Canonical Representation: Koopmans Theo.',&
             /2X,42('-'),//20X,'Ionization Potentials',11X,'1RDM Diag', &
             //19X,'(aU)',14X,'(eV)')
    2 FORMAT(2X,I4,4X,F15.6,4X,F15.6,9X,F8.6)
    3 FORMAT(/,15X,'Maximum 1RDM off-diagonal element:',F12.6,          &
               1X,'(',I3,',',I3,')')
    4 FORMAT(/,                                                         &
       18X,'-----------------------',/,                                 &
       18X,' Canonical HF Orbitals ',/,                                 &
       18X,'-----------------------')
!-----------------------------------------------------------------------
      DEALLOCATE (AUX,W,TEMP,DENMAT)
      RETURN
      END

! DQOHF
      SUBROUTINE DQOHF(IEMOM,USER,RO10)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      COMMON/INPNOF_NTHRESH/NTHRESHL,NTHRESHE,NTHRESHEC,NTHRESHEN
      COMMON/PUNTEROSUSER/N1,N2,N3,N4,N5,N6,N7,N8,N9,N10,N11,N12,N13,   &
                          N14,N15,N16,N17,N18,N19,N20,N21,N22,N23,N24,  &
                          N25,N26,N27,N28,N29,N30,N31,N32,N33,N34,N35,  &
                          N36,N37,N38,N39,N40,N41,N42,N43,N44,N45,N46,  &
                          N47,N48,N49,N50,N51,NUSER
      DOUBLE PRECISION,DIMENSION(NBF5)::RO10
      DOUBLE PRECISION,DIMENSION(NUSER)::USER
      DOUBLE PRECISION,PARAMETER::DFAC=2.54174D0    ! Debye
      DOUBLE PRECISION,PARAMETER::QFAC=1.345044D0   ! Buckinham
      DOUBLE PRECISION,PARAMETER::OFAC=7.117668D-01 ! X10**34 ESU-CM**3
!-----------------------------------------------------------------------
      IF(IEMOM>=1)THEN  ! Dipole
       CALL DIPMOMr(USER(N11),USER(N12),USER(N13),USER(N14),USER(N15),  &
                    USER(N16),USER(N17),USER(N7),RO10,DMXe,DMYe,DMZe,   &
                    DTx,DTy,DTz,DT)
       IF(NTHRESHE<=10)THEN
        WRITE(6,1)DT*DFAC,DT,DTx,DTy,DTz
       ELSE
        WRITE(6,11)DT*DFAC,DT,DTx,DTy,DTz
       ENDIF
      ENDIF
      IF(IEMOM>=2)THEN  ! Quadrupole
       CALL QUADMOMr(USER(N18),USER(N19),USER(N20),USER(N21),           &
                     USER(N22),USER(N23),USER(N24),USER(N25),           &
                     USER(N26),USER(N27),USER(N28),USER(N29),           &
                     USER(N30),USER(N7),RO10,                           &
                     QMXXe,QMYYe,QMZZe,QMXYe,QMXZe,QMYZe,               &
                     QTxx,QTyy,QTzz,QTxy,QTxz,QTyz)
       WRITE(6,2)QTxx*QFAC,QTyy*QFAC,QTzz*QFAC,                         &
                 QTxy*QFAC,QTxz*QFAC,QTyz*QFAC,                         &
                 QTxx,QTyy,QTzz,QTxy,QTxz,QTyz

      ENDIF
      IF(IEMOM==3)THEN  ! Octupole
       CALL OCTMOMr(USER(N31),USER(N32),USER(N33),USER(N34),            &
                    USER(N35),USER(N36),USER(N37),USER(N38),            &
                    USER(N39),USER(N40),USER(N41),USER(N42),            &
                    USER(N43),USER(N44),USER(N45),USER(N46),            &
                    USER(N47),USER(N48),USER(N49),USER(N50),            &
                    USER(N51),USER(N7),RO10,                            &
                    OMXXXe,OMYYYe,OMZZZe,OMXXYe,OMXXZe,                 &
                    OMXYYe,OMYYZe,OMXZZe,OMYZZe,OMXYZe,                 &
                    OTXXX,OTYYY,OTZZZ,OTXXY,OTXXZ,                      &
                    OTXYY,OTYYZ,OTXZZ,OTYZZ,OTXYZ)
       WRITE(6,3)OTXXX*OFAC,OTYYY*OFAC,OTZZZ*OFAC,                      &
                 OTXXY*OFAC,OTXXZ*OFAC,OTXYY*OFAC,                      &
                 OTYYZ*OFAC,OTXZZ*OFAC,OTYZZ*OFAC,                      &
                 OTXYZ*OFAC,OTXXX,OTYYY,OTZZZ,OTXXY,                    &
                 OTXXZ,OTXYY,OTYYZ,OTXZZ,OTYZZ,OTXYZ
      ENDIF
!-----------------------------------------------------------------------
!     FORMAT STATEMENTS
!-----------------------------------------------------------------------
    1 FORMAT(/,2X,'------------------',                                 &
              /2X,' HF Dipole Moment',                                  &
              /2X,'------------------',                                 &
            //,3X,F9.4,' Debye',' [',F9.4,                              &
               2X,'(',F9.4,',',F9.4,',',F9.4,')',' ]')
   11 FORMAT(/,2X,'------------------',                                 &
              /2X,' HF Dipole Moment',                                  &
              /2X,'------------------',                                 &
            //,3X,F15.10,' Debye',' [',F20.15,                          &
               2X,'(',F20.15,',',F20.15,',',F20.15,')',' ]')
    2 FORMAT(/,2X,'----------------------',                             &
              /2X,' HF Quadrupole Moment',                              &
              /2X,'----------------------',                             &
             //6X,'QXX',6X,'QYY',6X,'QZZ',6X,'QXY',6X,'QXZ',6X,'QYZ',   &
             //2X,6F9.4,2X,'(Buckingham)',//1X,'[',6F9.4,1X,']')
    3 FORMAT(/,2X,'--------------------',                               &
              /2X,' HF Octupole Moment',                                &
              /2X,'--------------------',                               &
             //6X,'OXXX',5X,'OYYY',5X,'OZZZ',5X,'OXXY',5X,'OXXZ',       &
               5X,'OXYY',5X,'OYYZ',5X,'OXZZ',5X,'OYZZ',5X,'OXYZ',       &
             //2X,10F9.4,2X,'(X10**34 ESU-CM**3)',                      &
             //1X,'[',10F9.4,1X,']')
!-----------------------------------------------------------------------
      RETURN
      END

!----------------------------------------------------------------------!
!                                                                      !
!   Restricted Hartree-Fock Calculation by means of Orbital rotaions   !
!  through ADAM (IRHF=2) or Iterative Diagonalization (IRHF=3) Methods !
!                                                                      !
!----------------------------------------------------------------------!

! HFIDADAMr
      SUBROUTINE HFIDADAMr(AHCORE,IJKL,XIJKL,XIJKaux,CHF,EiHF,USER,     &
                           ELAG,OVERLAP,DIPN,IRHF,MAXIT,IPRINTOPT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL ERIACTIVATED,HighSpin,CONVGDELAG,EFIELDL
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPFILE_NIJKL/NINTMX,NIJKL,NINTCR,NSTORE
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
      COMMON/ERIACT/ERIACTIVATED,NIJKaux,NINTCRaux,NSTOREaux,IAUXDIM
      COMMON/INP_EFIELDL/EFX,EFY,EFZ,EFIELDL
      COMMON/ELPROP/IEMOM
      COMMON/PUNTEROSUSER/N1,N2,N3,N4,N5,N6,N7,N8,N9,N10,N11,N12,N13,   &
                          N14,N15,N16,N17,N18,N19,N20,N21,N22,N23,N24,  &
                          N25,N26,N27,N28,N29,N30,N31,N32,N33,N34,N35,  &
                          N36,N37,N38,N39,N40,N41,N42,N43,N44,N45,N46,  &
                          N47,N48,N49,N50,N51,NUSER
      COMMON/CONVERGENCE/DUMEL,PCONV,CONVGDELAG
      COMMON/INPNOF_NTHRESH/NTHRESHL,NTHRESHE,NTHRESHEC,NTHRESHEN
      COMMON/ENERGIAS/EELEC,EELEC_OLD,DIF_EELEC,EELEC_MIN
      COMMON/EHFEN/EHF,EN
      COMMON/SKIPRINT/ISKIP
      COMMON/INPNOF_COEFOPT/MAXLOOP
!
      INTEGER :: ITCALL,ITLIM,ILOOP,IPRINTOPT,IORBOPT
      INTEGER(8),DIMENSION(NSTORE)::IJKL
      DOUBLE PRECISION,DIMENSION(NSTORE)::XIJKL
      DOUBLE PRECISION,DIMENSION(NSTOREaux)::XIJKaux
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::AHCORE,CHF,ELAG,OVERLAP
      DOUBLE PRECISION,DIMENSION(NUSER)::USER
      DOUBLE PRECISION,DIMENSION(NBF)::EiHF
      DOUBLE PRECISION,DIMENSION(3)::DIPN
!
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::RO10
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::CJ12HF,CK12HF
!-----------------------------------------------------------------------
!     Initial Values
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(RO10(NBF5),CJ12HF(NBF5,NBF5),CK12HF(NBF5,NBF5))
      IT = 0
      ITLIM = 1
      ILOOP = 0
      OCCTIME = 0.0d0
      CONVGDELAG = .FALSE.
!- - - - - - - - - - - - - - - - - - - - - - -
      RO10 = 0.0d0
      DO i=1,NB
       RO10(i) = 1.0d0
      ENDDO
      IF(NSOC>0)THEN
       DO i=NB+1,NA
        RO10(i) = 0.5d0
       ENDDO
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - -
      DO j=1,NBF5
       DO i=1,NBF5
        CJ12HF(j,i) = 2.0d0*RO10(j)*RO10(i)
        CK12HF(j,i) = RO10(j)*RO10(i)
       ENDDO
      ENDDO
      if(MSpin==0.and.NSOC>1)then
       DO j=NB+1,NA
        DO i=NB+1,NA
         CK12HF(j,i) = 2.0d0*RO10(j)*RO10(i)
        ENDDO
       ENDDO
      end if
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Restricted Hartree-Fock Method (IRHF=2,3)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(IRHF==2)THEN
       ISKIP=1
       MAXLOOP_ORI = MAXLOOP
       MAXLOOP = 50
       DO WHILE(IT<=MAXIT)
        IT=IT+1
        CALL OrbOptADAM(IT,OVERLAP,AHCORE,IJKL,XIJKL,XIJKaux,           &
                        USER(N7),CHF,RO10,CJ12HF,CK12HF,ELAG,DIPN,      &
                        ILOOP,OCCTIME,IPRINTOPT)
        EHF = EELEC
!       DUMEL<THRESHL & PCONV<THRESHE
        IF(CONVGDELAG)THEN
         if(NTHRESHE<=10)then
          WRITE(6,2)EHF
          WRITE(6,3)EHF+EN
         else
          WRITE(6,4)EHF
          WRITE(6,5)EHF+EN
         end if
         if(EFIELDL)WRITE(6,6)EFX,EFY,EFZ
!        HF Electrostatic Moments
         if(1<=IEMOM.and.IEMOM<=3)CALL DQOHF(IEMOM,USER,RO10)
        RETURN
         ISKIP=0
         MAXLOOP = MAXLOOP_ORI
         RETURN
        END IF
       END DO
      ELSE IF(IRHF==3)THEN
       ISKIP=0
       CALL HFIDr(AHCORE,IJKL,XIJKL,XIJKaux,CHF,EiHF,USER,IPRINTOPT)
      END IF
!-----------------------------------------------------------------------
!     FORMAT STATEMENTS
!-----------------------------------------------------------------------
    2 FORMAT(/,4X,'ELECTRONIC HF ENERGY =',F20.10)
    3 FORMAT(/,8X,' HF TOTAL ENERGY =',F20.10)
    4 FORMAT(/,4X,'ELECTRONIC HF ENERGY =',F25.15)
    5 FORMAT(/,8X,' HF TOTAL ENERGY =',F25.15)
    6 FORMAT(/,6X,'ELECTRIC FIELD (',D8.1,',',D8.1,',',D8.1,')')
!-----------------------------------------------------------------------
      DEALLOCATE(RO10,CJ12HF,CK12HF)
      RETURN
      END

!=======================================================================
