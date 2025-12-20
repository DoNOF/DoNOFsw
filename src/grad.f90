!======================================================================!
!                                                                      !
!           N O F   G R A D I E N T   S U B R O U T I N E S            !
!                                                                      !
!             2017  Module implemented by Ion Mitxelena                !
!                                                                      !
!                ( J. Chem. Phys. 146, 014102, 2017 )                  !
!                                                                      !
!======================================================================!
!                                                                      !
!   PNOFGRAD:  Main driver to compute PNOF gradient at one geometry    !
!   STVDERNOF: Main driver to compute one-electron PNOF gradient       !
!   VNNDERNOF: Nuclear repulsion contribution                          !
!   AOLAGRAN: Calculate PNOF lagrangian in AO basis                    !
!                                                                      !
!   SDERNOFl: Density force contribution                               !
!   JKDERNOFl: Main driver to compute two-electron HF gradient         !
!   JKDATMNOF: Select centers for derivatives                          !
!   JKDSHLNOFl: Select indices for shell block                         !
!   JKDNDXNOFl: Select indices for shell block                         !
!   DABNOF2PRE: Contract CJ12 and CK12 with density matrix to later    !
!               use in DABNOF2                                         !
!   DABNOF2l: Obtain PNOF 2e density for one shell block, N**5         !
!   DABNOF5l: Obtain PNOF 2e density for one shell block, M*N**4,      !
!   for PNOF5(DABNOF5PRE) and PNOF7(DABNOF7PRE) by using separability  !
!                                                                      !
!   JKDSPDNOFl: Evaluate derivative integral                           !
!   JKDINVNOF: Process derivative gradient and add to total gradient   !
!   DSPDFSNOFl: Compute derivative integrals, in principle these       !
!          integrals are added to the gradient 'on the fly', but each  !
!          integral contribution can be stored just removing a 'return'!
!                                                                      !
!   HELFEYNOFl: Hellmann-Feynmann contribution                         !
!   TVDERNOFl: 1e contribution due to AO derivatives with respect to   !
!              nuclei motion                                           !
!                                                                      !
!======================================================================!

! PNOFGRAD
      SUBROUTINE PNOFGRAD(COEF,QD,RO,ELAG,GRADS,ATMNAME,KATOM,KTYPE,    &
                          KLOC,KKMIN,KKMAX,KSTART,KNG,CX0,CY0,CZ0,ZNUC, &
                          EX1,CS,CP,CD,CF,CG,CJ12,CK12,XINTS,SIZE_ENV,  &
                          ENV,NAT,ATM,NBAS,BAS,IGTYP,IPRINTOPT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      LOGICAL HighSpin
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
      COMMON/INPNOF_PNOF/IPNOF,NTWOPAR
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      COMMON/USELIBCINT/ILIBCINT      
!
      INTEGER,INTENT(IN)::IPRINTOPT
      CHARACTER*4 ATMNAME(NATOMS)
      INTEGER,DIMENSION(NSHELL),INTENT(IN)::KATOM,KTYPE,KLOC,KKMIN,KKMAX
      INTEGER,DIMENSION(NSHELL),INTENT(IN)::KSTART,KNG
      DOUBLE PRECISION,DIMENSION(NATOMS),INTENT(IN)::CX0,CY0,CZ0,ZNUC
      DOUBLE PRECISION,DIMENSION(NPRIMI),INTENT(IN)::EX1,CS,CP,CD,CF,CG
      DOUBLE PRECISION,DIMENSION(NBF,NBF),INTENT(IN)::COEF
      DOUBLE PRECISION,DIMENSION(NBF,NBF,NBF),INTENT(IN)::QD
      DOUBLE PRECISION,DIMENSION(NBF5),INTENT(IN)::RO
      DOUBLE PRECISION,DIMENSION(NBF,NBF),INTENT(IN)::ELAG
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5),INTENT(IN)::CJ12,CK12
      DOUBLE PRECISION,DIMENSION((NSHELL*NSHELL+NSHELL)/2),INTENT(IN)::XINTS
      DOUBLE PRECISION,DIMENSION(3,NATOMS),INTENT(INOUT):: GRADS

      INTEGER :: SIZE_ENV,NBAS,IGTYP
      DOUBLE PRECISION :: ENV(SIZE_ENV)
      INTEGER :: ATM(6,NAT), BAS(8,NBAS)      
!-----------------------------------------------------------------------
!     ONE-ELECTRON CONTRIBUTION TO THE GRADIENT
!-----------------------------------------------------------------------
      CALL STVDERNOF(COEF,QD,RO,ELAG,GRADS,KATOM,KTYPE,KLOC,KKMIN,KKMAX,&
                     KSTART,KNG,CX0,CY0,CZ0,ZNUC,EX1,CS,CP,CD,CF,CG,    &
                     SIZE_ENV,ENV,NAT,ATM,NBAS,BAS,IGTYP)
!-----------------------------------------------------------------------
!     TWO-ELECTRON CONTRIBUTION TO THE GRADIENT
!-----------------------------------------------------------------------
      if(ILIBCINT==0)then
!HONDO IF((IPNOF==5 .or. IPNOF==7) .and. NSOC==0) THEN
!HONDO  CALL JKDERNOF5(KATOM,KTYPE,KLOC,KKMIN,KKMAX,KSTART,KNG,QD,RO,   &
!HONDO                 GRADS,CX0,CY0,CZ0,EX1,CS,CP,CD,CF,CG,ATMNAME,    &
!HONDO                 XINTS,IPRINTOPT)
!HONDO ELSE
!HONDO  CALL JKDERNOF(KATOM,KTYPE,KLOC,KKMIN,KKMAX,KSTART,KNG,CJ12,     &
!HONDO                CK12,QD,RO,GRADS,CX0,CY0,CZ0,EX1,CS,CP,CD,CF,CG,  &
!HONDO                ATMNAME,XINTS,IPRINTOPT)
!HONDO ENDIF
      else if(ILIBCINT==1)then                    
       IF((IPNOF==5 .or. IPNOF==7) .and. NSOC==0) THEN
        CALL JKDERNOF5l(KATOM,KTYPE,KLOC,KKMIN,KKMAX,QD,RO,GRADS,       &
                          ATMNAME,XINTS,SIZE_ENV,ENV,NAT,ATM,NBAS,BAS,  &
                          IGTYP,IPRINTOPT)                  
       ELSE                                                            
        CALL JKDERNOFl(KATOM,KTYPE,KLOC,KKMIN,KKMAX,CJ12,CK12,          &
                         QD,RO,GRADS,ATMNAME,XINTS,SIZE_ENV,ENV,NAT,ATM,&
                         NBAS,BAS,IGTYP,IPRINTOPT)
       ENDIF
      end if
!-----------------------------------------------------------------------
      RETURN
      END      

! STVDERNOF
      SUBROUTINE STVDERNOF(COEF,QD,RO,ELAG,GRADS,KATOM,KTYPE,KLOC,      &
                           KKMIN,KKMAX,KSTART,KNG,CX0,CY0,CZ0,ZNUC,     &
                           EX1,CS,CP,CD,CF,CG,SIZE_ENV,ENV,NAT,ATM,NBAS,&
                           BAS,IGTYP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      LOGICAL FROZEN
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_FROZEN/FROZEN,IFROZEN(200) 
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      COMMON/USELIBCINT/ILIBCINT
!      
      INTEGER,DIMENSION(NSHELL),INTENT(IN)::KATOM,KTYPE,KLOC,KKMIN,KKMAX
      INTEGER,DIMENSION(NSHELL),INTENT(IN)::KSTART,KNG  
      DOUBLE PRECISION,DIMENSION(NATOMS),INTENT(IN)::CX0,CY0,CZ0,ZNUC
      DOUBLE PRECISION,DIMENSION(NPRIMI),INTENT(IN)::EX1,CS,CP,CD,CF,CG
      DOUBLE PRECISION,DIMENSION(NBF,NBF),INTENT(IN)::COEF
      DOUBLE PRECISION,DIMENSION(NBF,NBF,NBF),INTENT(IN)::QD
      DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::DM2
      DOUBLE PRECISION,DIMENSION(NBF5),INTENT(IN)::RO
      DOUBLE PRECISION,DIMENSION(NBF,NBF),INTENT(IN)::ELAG
      DOUBLE PRECISION,DIMENSION(3,NATOMS),INTENT(INOUT):: GRADS
      DOUBLE PRECISION,DIMENSION(NBFT)::LEPS
      DOUBLE PRECISION,PARAMETER::ZERO=0.0D+00

      INTEGER :: SIZE_ENV,NBAS,IGTYP
      DOUBLE PRECISION :: ENV(SIZE_ENV)
      INTEGER :: ATM(6,NAT), BAS(8,NBAS)     
!-----------------------------------------------------------------------
!     MAKE GRADIENT ZERO
      GRADS=ZERO
!-----------------------------------------------------------------------
!     Nuclear repulsion force
!-----------------------------------------------------------------------
      CALL VNNDERNOF(CX0,CY0,CZ0,ZNUC,GRADS)
!-----------------------------------------------------------------------
!     Density force contribution
!-----------------------------------------------------------------------      
!     GET LAGRANGIAN MATRIX (LEPS)    
      CALL AOLAGRAN(COEF,ELAG,LEPS,NBF,NBFT)
!     COMPUTE AND ADD OVERLAP DERIVATIVES
      if(ILIBCINT==0)then
!HONDO CALL SDERNOF(KATOM,KTYPE,KLOC,KKMIN,KKMAX,KSTART,KNG,            &
!HONDO              CX0,CY0,CZ0,EX1,CS,CP,CD,CF,CG,LEPS,GRADS)
      else if(ILIBCINT==1)then     
       CALL SDERNOFl(NSHELL,NBF,NBFT,LEPS,SIZE_ENV,ENV,                 &
                     NAT,ATM,NBAS,BAS,IGTYP,GRADS)
      end if      
!-----------------------------------------------------------------------
!     One-electron hamiltonian contribution
!-----------------------------------------------------------------------
!     CONTRACT OVER OCCUPATIONS TO OBTAIN DENSITY MATRIX
      ALLOCATE(DM2(NBF,NBF))
      DM2=ZERO
      DO J=1,NBF
       DO L=1,J
        DO I=1,NBF5
          DM2(J,L)=DM2(J,L)+QD(I,J,L)*RO(I)
        ENDDO
        DM2(L,J)=DM2(J,L)
       ENDDO
      ENDDO
      DM2=DM2+DM2
!      
!     HELLMANN-FEYNMAN FORCE      
!
      if(ILIBCINT==0)then
!HONDO CALL HELFEYNOF(KATOM,KTYPE,KLOC,KKMIN,KKMAX,KSTART,KNG,          &
!HONDO                CX0,CY0,CZ0,EX1,CS,CP,CD,CF,CG,ZNUC,DM2,GRADS)
      else if(ILIBCINT==1)then       
       CALL HELFEYNOFl(NBF,NBFT,NSHELL,CX0,CY0,CZ0,ZNUC,DM2,GRADS,      &
                       SIZE_ENV,ENV,NAT,ATM,NBAS,BAS,IGTYP)
      end if            
!                  
!     INTEGRAL FORCE (AO DERIVATIVE CONTRIBUTION)                  
!
      if(ILIBCINT==0)then
!HONDO CALL TVDERNOF(KATOM,KTYPE,KLOC,KKMIN,KKMAX,KSTART,KNG,           &
!HONDO               CX0,CY0,CZ0,EX1,CS,CP,CD,CF,CG,ZNUC,DM2,GRADS)
      else if(ILIBCINT==1)then                   
       CALL TVDERNOFl(NBF,NBFT,NSHELL,DM2,GRADS,SIZE_ENV,ENV,NAT,ATM,   &
                      NBAS,BAS,IGTYP)
      end if
!                  
      DEALLOCATE(DM2)
!     
!     MAKE ZERO GRADIENT CORRESPONDING TO FROZEN COORDINATES
!
      IF(FROZEN) THEN
        DO I=1,200,2
         IF(IFROZEN(I)==0) EXIT
         GRADS(IFROZEN(I),IFROZEN(I+1))=ZERO
        ENDDO
      ENDIF
!-----------------------------------------------------------------------
      RETURN
      END   

! VNNDERNOF
      SUBROUTINE VNNDERNOF(CX0,CY0,CZ0,ZNUC,DE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      DOUBLE PRECISION,DIMENSION(NATOMS),INTENT(IN)::CX0,CY0,CZ0,ZNUC
      DOUBLE PRECISION,DIMENSION(3,NATOMS),INTENT(INOUT)::DE
      DOUBLE PRECISION,DIMENSION(3,NATOMS)::C
      DOUBLE PRECISION,DIMENSION(NATOMS,NATOMS)::DRG
      DOUBLE PRECISION,PARAMETER::ZERO=0.0D+00,ONE=1.0D+00
!
      C(1,:)=CX0
      C(2,:)=CY0
      C(3,:)=CZ0
!
!     ----- FORM DISTANCE MATRIX -----
!
      DRG(1,1) = ZERO
      DO K = 2,NATOMS
         DRG(K,K) = ZERO
         K1 = K-1
         DO L = 1,K1
            RKL = ZERO
            DO I = 1,3
               RKL = RKL+(C(I,K)-C(I,L))**2
            ENDDO
            DRG(K,L) = -ONE/RKL
            DRG(L,K) = SQRT(RKL)
         ENDDO
      ENDDO
!
!     ----- NUCLEAR REPULSION CONTRIBUTION TO GRADIENT -----
!
      DO KK = 1,3
         DO K = 2,NATOMS
            ZAK = ZNUC(K)
            KM1 = K-1
            DO L = 1,KM1
               ZAL = ZNUC(L)
               PKL = (C(KK,K)-C(KK,L))/DRG(L,K)
               DE(KK,K) = DE(KK,K)+PKL*DRG(K,L)*ZAK*ZAL
            ENDDO
         ENDDO
!
         NAT1 = NATOMS-1
         DO K = 1,NAT1
            ZAK = ZNUC(K)
            KP1 = K+1
            DO L = KP1,NATOMS
               ZAL = ZNUC(L)
               PKL = (C(KK,K)-C(KK,L))/DRG(K,L)
               DE(KK,K) = DE(KK,K)+PKL*DRG(L,K)*ZAK*ZAL
            ENDDO
         ENDDO
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! AOLAGRAN
      SUBROUTINE AOLAGRAN(C,EE,LEPS,N,NL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      DOUBLE PRECISION,DIMENSION(N,N),INTENT(IN) :: C,EE
      DOUBLE PRECISION,DIMENSION(NL),INTENT(OUT) :: LEPS
      DOUBLE PRECISION,DIMENSION(N,N) :: AUXELG,AELG
      DOUBLE PRECISION :: AUX
      DOUBLE PRECISION,PARAMETER :: PT5=0.5D+00,ZERO=0.0D+00
!-----------------------------------------------------------------------
!     GET THE LAGRANGIAN MATRIX FOR PNOF (COEF*LAMBDA*COEF)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      LEPS=ZERO
!     SIMETRIZATION OF LAMBDA MATRIX
      DO I=1,NBF
       DO J=1,I
        AUXELG(I,J)=EE(I,J)
        AUXELG(J,I)=AUXELG(I,J)
       ENDDO
      ENDDO
!     SCALING AS N**3 CALCULATION, FIRST CONTRACTING OVER A SINGLE COEF
      AELG=ZERO
      DO L=1,N
       DO M=1,NBF
        DO NM=1,NBF
          AELG(M,L)=AELG(M,L)+AUXELG(M,NM)*C(L,NM)
        ENDDO
       ENDDO
      ENDDO
      KL=0
      DO K=1,N
       DO L=1,K
        KL=KL+1
        AUX=ZERO
        DO M=1,NBF
          AUX = AUX - AELG(M,L)*C(K,M)
        ENDDO
        LEPS(KL) = AUX + AUX
       ENDDO
!     DOUBLE COUNTING OF DIAGONAL IS AVOIDED IN TRACEs
      ENDDO
      LEPS = LEPS + LEPS
!-----------------------------------------------------------------------
      RETURN
      END

!----------------------------------------------------------------------!
!                                                                      !
!       2025 Use LIBCINT open source library for ERI calculation       !
!                                                                      !
!  implemented by Juan Felipe Huan Lew Yee and Jorge Martin del Campo  ! 
!                                                                      !
!----------------------------------------------------------------------!

! SDERNOFl
      SUBROUTINE SDERNOFl(NSHELL,NBF,NBFT,EPS,SIZE_ENV,ENV,             &
                          NAT,ATM,NBAS,BAS,IGTYP,DE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DOUBLE PRECISION, DIMENSION(NBFT), INTENT(IN) :: EPS
      DOUBLE PRECISION, DIMENSION(3, NAT), INTENT(INOUT) :: DE

      INTEGER :: II, JJ, I, J, LI, LJ, IAT, JAT, ISH, ERR, IGTYP
      INTEGER :: SIZE_ENV, NAT, NBAS
      DOUBLE PRECISION :: ENV(SIZE_ENV)
      INTEGER :: IA(NBF), LOC(NSHELL), Dcgto(NSHELL)

      INTEGER :: DI, DJ
      INTEGER(4) :: SHLS(2)
      INTEGER :: ATM(6, NAT), BAS(8, NBAS)

      DOUBLE PRECISION, ALLOCATABLE :: BLK(:,:,:)

      INTEGER, EXTERNAL :: CINTcgto_spheric, CINT1e_ipovlp_sph
      INTEGER, EXTERNAL :: CINTcgto_cart, CINT1e_ipovlp_cart
!-----------------------------------------------------------------------
! LOC indicates the starting location of each shell in the AO basis
!-----------------------------------------------------------------------
      LOC(1) = 1
      IF(IGTYP==1) DI = CINTcgto_cart(0, BAS)
      IF(IGTYP==2) DI = CINTcgto_spheric(0, BAS)
      Dcgto(1) = DI
      DO ISH = 2,NSHELL
        LOC(ISH) = LOC(ISH-1) + DI
        IF(IGTYP==1) DI = CINTcgto_cart(ISH-1, BAS)
        IF(IGTYP==2) DI = CINTcgto_spheric(ISH-1, BAS)
        Dcgto(ISH) = DI
      END DO

!
      DO I=1,NBF
        IA(I) = (I*I-I)/2
      ENDDO
!
!     ----- I SHELL
!
      DO II = 1,NSHELL
        IAT = BAS(1,II) + 1
        SHLS(1) = II - 1
        DI = Dcgto(II)
!
!     ----- J SHELL
!
        DO JJ = 1,II
          IF(II.EQ.JJ) CYCLE
          JAT = BAS(1, JJ) + 1
          SHLS(2) = JJ - 1
          DJ = Dcgto(JJ)

          ALLOCATE(BLK(DI, DJ, 3))
          BLK = 0.0D0

          IF(IGTYP==1) ERR = cint1e_ipovlp_cart(BLK, SHLS, ATM, NAT,   &
                  BAS, NBAS, ENV)
          IF(IGTYP==2) ERR = cint1e_ipovlp_sph(BLK, SHLS, ATM, NAT,    &
                  BAS, NBAS, ENV)
          DO I=1, DI
            LI = LOC(II) + I - 1
            DO J=1,DJ
              LJ = LOC(JJ) + J - 1
              IJ=IA(LI) + LJ
              DEN = EPS(IJ)
!
              DE(1, IAT) = DE(1, IAT) - BLK(I, J, 1)*DEN
              DE(2, IAT) = DE(2, IAT) - BLK(I, J, 2)*DEN
              DE(3, IAT) = DE(3, IAT) - BLK(I, J, 3)*DEN
              DE(1, JAT) = DE(1, JAT) + BLK(I, J, 1)*DEN
              DE(2, JAT) = DE(2, JAT) + BLK(I, J, 2)*DEN
              DE(3, JAT) = DE(3, JAT) + BLK(I, J, 3)*DEN
            END DO
          END DO
        DEALLOCATE(BLK)
        ENDDO
      ENDDO
!-----------------------------------------------------------------------      
      RETURN
      END

! HELFEYNOFl
      SUBROUTINE HELFEYNOFl(NBF,NBFT,NSHELL,CX0,CY0,CZ0,ZAN,PM,DE,      &
                            SIZE_ENV,ENV,NAT,ATM,NBAS,BAS,IGTYP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
      INTEGER, INTENT(IN) :: NAT,NSHELL,NBF,NBFT,SIZE_ENV,NBAS,IGTYP

      DOUBLE PRECISION, DIMENSION(NAT), INTENT(IN) :: CX0,CY0,CZ0,ZAN
      DOUBLE PRECISION, DIMENSION(NBF,NBF), INTENT(IN) :: PM
      DOUBLE PRECISION, DIMENSION(NBFT) :: P
      DOUBLE PRECISION, DIMENSION(3,NAT), INTENT(INOUT) :: DE

      INTEGER(4) :: DI, DJ
      INTEGER(4) :: SHLS(2)
      INTEGER :: II, JJ, I, J, LI, LJ, IC, ERR
      INTEGER,DIMENSION(NSHELL) :: LOC
      INTEGER, DIMENSION(NBF) :: IA
      INTEGER :: Dcgto(NSHELL)


      DOUBLE PRECISION :: ENV(SIZE_ENV)
      INTEGER :: ATM(6,NAT), BAS(8,NBAS)

      DOUBLE PRECISION, ALLOCATABLE :: BLK(:,:,:)

      DOUBLE PRECISION :: ZNUCC, CX, CY, CZ

      INTEGER, EXTERNAL :: CINTcgto_spheric, CINT1e_iprinv_sph
      INTEGER, EXTERNAL :: CINTcgto_cart, CINT1e_iprinv_cart

!-----------------------------------------------------------------------
! LOC indicates the starting location of each shell in the AO basis
!-----------------------------------------------------------------------
      LOC(1) = 1
      IF(IGTYP==1) DI = CINTcgto_cart(0, BAS)
      IF(IGTYP==2) DI = CINTcgto_spheric(0, BAS)
      Dcgto(1) = DI
      DO ISH = 2,NSHELL
        LOC(ISH) = LOC(ISH-1) + DI
        IF(IGTYP==1) DI = CINTcgto_cart(ISH-1, BAS)
        IF(IGTYP==2) DI = CINTcgto_spheric(ISH-1, BAS)
        Dcgto(ISH) = DI
      END DO

!     ----- HELMANN-FEYNMAN GRADIENT TERM -----
!     INTEGRAL TYPE IS <II/H'/JJ> = <II/V'/JJ>
!                     
      DO I=1, NBF
        IA(I) = (I*I-I)/2
      ENDDO

      CALL SQUARETRIAN2(PM, P, NBF, NBFT)
      P = P * 0.5D+00
!
!     ----- I SHELL
!
      DO II = 1, NSHELL
!
        SHLS(1) = II-1
        DI = Dcgto(II)
!
!     ----- J SHELL
!
        DO JJ = 1, NSHELL
!
          SHLS(2) = JJ-1
          DJ = Dcgto(JJ)

          DO IC = 1, NAT
            ZNUCC = ZAN(IC)
            CX = CX0(IC)
            CY = CY0(IC)
            CZ = CZ0(IC)

            ENV(5) = CX
            ENV(6) = CY
            ENV(7) = CZ
            ALLOCATE(BLK(DI,DJ,3))
            BLK = 0.0D0
            if(IGTYP==1) ERR = cint1e_iprinv_cart(BLK, SHLS, ATM, NAT, &
                    BAS, NBAS, ENV)
            if(IGTYP==2) ERR = cint1e_iprinv_sph(BLK, SHLS, ATM, NAT,  &
                    BAS, NBAS, ENV)

            DO I=1,DI
              LI = LOC(II) + I - 1
              DO J=1,DJ
                LJ = LOC(JJ) + J - 1
                DEN = PM(LI, LJ)
                DE(1, IC) = DE(1, IC) - BLK(I, J, 1) * DEN * ZNUCC
                DE(2, IC) = DE(2, IC) - BLK(I, J, 2) * DEN * ZNUCC
                DE(3, IC) = DE(3, IC) - BLK(I, J, 3) * DEN * ZNUCC
              ENDDO
            ENDDO
            DEALLOCATE(BLK)

          ENDDO
!
!     ----- END OF *ATOMS* LOOPS -----
!
        ENDDO
      ENDDO
      ENV(5:7) = 0.0D0
!
!     ----- END OF *SHELL* LOOPS -----
!             
!-----------------------------------------------------------------------
      RETURN
      END      

! TVDERNOFl
      SUBROUTINE TVDERNOFl(NBF,NBFT,NSHELL,PM,DE,SIZE_ENV,ENV,          &
                           NAT,ATM,NBAS,BAS, IGTYP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      
      INTEGER, INTENT(IN) :: NAT, NBAS, SIZE_ENV, IGTYP
      DOUBLE PRECISION, DIMENSION(NBF,NBF), INTENT(IN) :: PM
      DOUBLE PRECISION, DIMENSION(3,NAT), INTENT(INOUT) :: DE
      DOUBLE PRECISION, DIMENSION(SIZE_ENV), INTENT(IN) :: ENV
      INTEGER, DIMENSION(6,NAT), INTENT(IN) :: ATM
      INTEGER, DIMENSION(8,NBAS), INTENT(IN) :: BAS
    
      DOUBLE PRECISION, DIMENSION(NBFT) :: P
      INTEGER, DIMENSION(NBF) :: IA
      INTEGER, DIMENSION(NSHELL) :: LOC
      INTEGER :: Dcgto(NSHELL)
      DOUBLE PRECISION, ALLOCATABLE :: BLK(:,:,:)
      INTEGER(4) :: DI, DJ, SHLS(2)
      INTEGER :: II, JJ, I, J, LI, LJ, IAT, JAT, ERR

      INTEGER, EXTERNAL :: CINTcgto_spheric, CINT1e_ipkin_sph
      INTEGER, EXTERNAL :: CINT1e_nuc_sph
      INTEGER, EXTERNAL :: CINTcgto_cart, CINT1e_ipnuc_cart
      INTEGER, EXTERNAL :: CINT1e_nuc_cart
!-----------------------------------------------------------------------
! LOC indicates the starting location of each shell in the AO basis
!-----------------------------------------------------------------------
      LOC(1) = 1
      IF(IGTYP==1) DI = CINTcgto_cart(0, BAS)
      IF(IGTYP==2) DI = CINTcgto_spheric(0, BAS)
      Dcgto(1) = DI
      DO ISH = 2,NSHELL
        LOC(ISH) = LOC(ISH-1) + DI
        IF(IGTYP==1) DI = CINTcgto_cart(ISH-1, BAS)
        IF(IGTYP==2) DI = CINTcgto_spheric(ISH-1, BAS)
        Dcgto(ISH) = DI
      END DO
!-----------------------------------------------------------------------
!
!     ----- BASIS FUNCTION DERIVATIVE CONTRIBUTIONS TO GRADIENT -----
!     INTEGRALS ARE OF TYPE <II'/H/JJ> = <II'/T+V/JJ>
!
      DO I=1, NBF
        IA(I) = (I*I-I)/2
      ENDDO
      IAZ=0
      CALL SQUARETRIAN2(PM, P, NBF, NBFT)
      P = P * 0.5D+00  
!
!     ----- I SHELL
!
      DO II = 1, NSHELL
!      
        IAT = BAS(1,II) + 1
        SHLS(1) = II - 1
        DI = Dcgto(II)
!
!     ----- J SHELL
!
        DO JJ = 1,NSHELL
!
          JAT = BAS(1,JJ) + 1
          SHLS(2) = JJ - 1
          DJ = Dcgto(JJ)

          ALLOCATE(BLK(DI, DJ, 3))
          BLK = 0.0D0

          IF(IGTYP==1) ERR = cint1e_ipkin_cart(BLK, SHLS, ATM, NAT,    &
                  BAS, NBAS, ENV)
          IF(IGTYP==2) ERR = cint1e_ipkin_sph(BLK, SHLS, ATM, NAT,     &
                  BAS, NBAS, ENV)
          DO I=1,DI
            LI = LOC(II) + I - 1
            DO J=1,DJ
              LJ = LOC(JJ) + J - 1
              NN = IA(MAX0(LI, LJ))+MIN0(LI, LJ)
              DEN = P(NN)
              DE(1, IAT)=DE(1, IAT) - BLK(I, J, 1) * DEN
              DE(2, IAT)=DE(2, IAT) - BLK(I, J, 2) * DEN
              DE(3, IAT)=DE(3, IAT) - BLK(I, J, 3) * DEN
              IJ = IJ + 3
            ENDDO
          ENDDO
        DEALLOCATE(BLK)

!
!     ..... NUCLEAR ATTRACTION
!
        ALLOCATE(BLK(DI,DJ,3))
        BLK = 0.0D0
        IF(IGTYP==1) ERR = cint1e_ipnuc_cart(BLK, SHLS, ATM, NAT, BAS,&
                NBAS, ENV)
        IF(IGTYP==2) ERR = cint1e_ipnuc_sph(BLK, SHLS, ATM, NAT, BAS, &
                NBAS, ENV)
        DO I=1, DI
          LI = LOC(II) + I - 1
            DO J=1, DJ
              LJ = LOC(JJ) + J - 1
!             IF((IC.GT.NATOMS).AND.(IAT.EQ.IAZ)) CYCLE
              NN = IA(MAX0(LI, LJ)) + MIN0(LI, LJ)
              DEN = P(NN)
              DE(1, IAT)=DE(1, IAT) - BLK(I, J, 1) * DEN
              DE(2, IAT)=DE(2, IAT) - BLK(I, J, 2) * DEN
              DE(3, IAT)=DE(3, IAT) - BLK(I, J, 3) * DEN
              IJ = IJ + 3
            ENDDO
          ENDDO
          DEALLOCATE(BLK)
        END DO
      END DO
!
!     ----- END OF SHELL LOOPS -----
!
!-----------------------------------------------------------------------
      RETURN
      END

! JKDERNOFl
      SUBROUTINE JKDERNOFl(KATOM,KTYPE,KLOC,KKMIN,KKMAX,CJ12,CK12,      &
                           QD,RO,GRADS,ATMNAME,XINTS,SIZE_ENV,ENV,NAT,  &
                           ATM,NBAS,BAS,IGTYP,IPRINTOPT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      LOGICAL FROZEN
      COMMON/INPNOF_FROZEN/FROZEN,IFROZEN(200)      
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      COMMON/INTOPT/CUTOFF,ISCHWZ,IECP,NECP
!     INPUT-OUTPUT VARIABLES OR ARGUMENTS
      INTEGER,DIMENSION(NSHELL),INTENT(IN)::KATOM,KTYPE,KLOC,KKMIN,KKMAX
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5),INTENT(IN)::CJ12,CK12
      DOUBLE PRECISION,DIMENSION(NBF,NBF,NBF),INTENT(IN)::QD
      DOUBLE PRECISION,DIMENSION(NBF5),INTENT(IN)::RO
      DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::P,DAAUX,DAAUX2
      DOUBLE PRECISION,DIMENSION((NSHELL*NSHELL+NSHELL)/2),INTENT(IN)::XINTS
      CHARACTER*4 ATMNAME(NATOMS)
      DOUBLE PRECISION,DIMENSION(3,NATOMS),INTENT(INOUT)::GRADS
      INTEGER,INTENT(IN)::IPRINTOPT
!     INTERMEDIATE VARIABLES 
      INTEGER,DIMENSION(NBF)::IA
      INTEGER(4) :: DI,DJ,DK,DL
      INTEGER,DIMENSION(35)::IGXYZ,JGXYZ,KGXYZ,LGXYZ
      INTEGER,DIMENSION(:),ALLOCATABLE::IJKLG
      INTEGER::MAXNUM
      INTEGER::TOTCOUNT,INVTYP,KANG,I
      INTEGER::II,JJ,KK,LL,MAXLL,ISH,JSH,KSH,LSH,IIAT,JJAT,KKAT,LLAT
      INTEGER::MINJ,MAXJ,MINK,MAXK,MINL,MAXL,MINI,MAXI
      INTEGER::NUMI,NUMJ,NUMK,NUML
      LOGICAL::SKIPI,SKIPJ,SKIPK,SKIPL
      LOGICAL::IIEQJJ,KKEQLL,IJEQKL
      DOUBLE PRECISION::DABMAX,CUTOFF2,GMAX
      DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::DAB
      DOUBLE PRECISION,PARAMETER::ZERO=0.0D+00,ONE=1.0D+00
!     LIBCINT
      INTEGER :: SIZE_ENV,NBAS,IGTYP
      DOUBLE PRECISION :: ENV(SIZE_ENV)
      INTEGER :: ATM(6,NAT), BAS(8,NBAS)
!-----------------------------------------------------------------------
!     ROUTINE BASCHK DO:
      KANG=0
      DO N=1,NSHELL
       IF(KTYPE(N).GT.KANG) KANG = KTYPE(N)
      ENDDO
!
!     SET STARTING PARAMETERS
!
!     CUTOFF IS THE SCHWARZ SCREENING CUT OFF 
      CUTOFF2 = CUTOFF*0.5D+00
!     SET POINTERS FOR PARTITIONING MEMORY
      DO I=1,NBF
        IA(I) = (I*I-I)/2
      ENDDO   
!           
!     MAXNUM=NUMBER OF FUNCTIONS WITH ANGULAR MOMENTUM EQUAL TO MAXTYP
      MAXNUM=((KANG)*(KANG+1))/2
!     DO AT LEAST AN L SHELL      
      IF(MAXNUM.LT.4) MAXNUM=4
      MAXNUM=(MAXNUM**4)
!
      ALLOCATE(IJKLG(MAXNUM))
      ALLOCATE(DAB(MAXNUM))
      ALLOCATE(P(NBF,NBFT))
      CALL SQUARETRIAN3(QD,P,NBF,NBFT)
      P=P*0.5D+00
      ALLOCATE(DAAUX(NBF,NBFT))
      ALLOCATE(DAAUX2(NBF,NBFT))
      CALL DABNOF2PRE(CJ12,CK12,RO,P,DAAUX,DAAUX2)
!
      TOTCOUNT=1
!
!----I SHELL
!
      DO II=1,NSHELL
!      
!-----J SHELL
!
       DO JJ=1,II
!        
!-----K SHELL
!
        DO KK=1,II
!        
!-----L SHELL
!
        MAXLL=KK
        IF(KK.EQ.II) MAXLL=JJ
         DO LL=1,MAXLL
!
!         IMPLEMENT INTEGRAL SCREENING HERE USING EXCHANGE INTEGRALS
!                                                                       
          IJIJ=IA(MAX0(II,JJ))+MIN0(II,JJ)                          
          KLKL=IA(MAX0(KK,LL))+MIN0(KK,LL)                          
          GMAX=(XINTS(IJIJ)*XINTS(KLKL))
!
!         COARSE SCREENING, ON JUST THE INTEGRAL VALUE         
!         ONLY WORKS IF SCHWARZ SCREENING IS ON (NATOMS>5)
          IF (GMAX.LT.CUTOFF.AND.NATOMS>5) CYCLE
!
          ISH=II
          JSH=JJ
          KSH=KK
          LSH=LL
!          
!         SELECT CENTERS FOR DERIVATIVES
!
          CALL JKDATMNOF(ISH,JSH,KSH,LSH,SKIPI,SKIPJ,SKIPK,SKIPL,       &
                         INVTYP,KATOM,IIAT,JJAT,KKAT,LLAT)              
          IF(SKIPI.AND.SKIPJ.AND.SKIPK.AND.SKIPL) CYCLE                 
!                                                                       
!         SELECT INDICES FOR SHELL BLOCK                                
!                                                                       
          CALL JKDSHLNOFl(ISH,JSH,KSH,LSH,IIEQJJ,KKEQLL,IJEQKL,DI,DJ,   &
                          DK,DL,SIZE_ENV,ENV,NAT,ATM,NBAS,BAS,IGTYP)
          CALL JKDNDXNOFl(DI,DJ,DK,DL,IGXYZ,JGXYZ,KGXYZ,LGXYZ,IIEQJJ,   &
                          KKEQLL,IJEQKL,IJKLG,MAXNUM,ISH,JSH,KSH,LSH,   &
                          NBF,NSHELL,SIZE_ENV,ENV,NAT,ATM,NBAS,BAS,     &
                          IGTYP)
!                      
!         OBTAIN 2e DENSITY FOR THIS SHELL BLOCK                      
!
          CALL DABNOF2l(ISH,JSH,KSH,LSH,IGXYZ,JGXYZ,KGXYZ,LGXYZ,        &
                        IIEQJJ,KKEQLL,IJEQKL,IA,P,DAB,MAXNUM,DABMAX,    &
                        DAAUX,DAAUX2,NAT,ATM,NBAS,BAS,IGTYP)
!
!         FINE SCREENING, ON INTEGRAL VALUE TIMES DENSITY FACTOR
!         ONLY WORKS IF SCHWARZ SCREENING IS ON (NATOMS>5)
          IF(DABMAX*GMAX.LT.CUTOFF2.AND.NATOMS>5) CYCLE
!
!         EVALUATE DERIVATIVE INTEGRAL AND ADD TO THE GRADIENT
!
          CALL JKDSPDNOFl(TOTCOUNT,DAB,II,JJ,KK,LL,IIEQJJ,KKEQLL,       &
                          IJEQKL,SKIPI,SKIPJ,SKIPK,SKIPL,IJKLG,         &
                          MAXNUM,INVTYP,IIAT,JJAT,KKAT,LLAT,            &
                          SIZE_ENV,ENV,NAT,ATM,NBAS,BAS,IGTYP,GRADS)

!      
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      DEALLOCATE(P,IJKLG,DAB,DAAUX,DAAUX2)
!     
!     MAKE ZERO GRADIENT CORRESPONDING TO FROZEN COORDINATES
!
      IF(FROZEN) THEN
        DO I=1,200,2
         IF(IFROZEN(I).EQ.0) EXIT
         GRADS(IFROZEN(I),IFROZEN(I+1))=ZERO
        ENDDO
      ENDIF      
!
!     PRINT OUT TOTAL PNOF ENERGY GRADIENT
!
      IF(IPRINTOPT==1)THEN
       WRITE(6,1)
       DO I=1,NATOMS
        WRITE(6,2)I,ATMNAME(I),GRADS(1,I),GRADS(2,I),GRADS(3,I)
       ENDDO
      ENDIF
!-----------------------------------------------------------------------
    1 FORMAT( /1X,'----------------',                                   &
              /1X,' Total Gradient ',                                   &
              /1X,'----------------',                                   &
              //9X,'Atom',7X,'Ex',10X,'Ey',10X,'Ez' )
    2 FORMAT(/1X,I4,5X,A4,F10.4,2X,F10.4,2X,F10.4)
!-----------------------------------------------------------------------
      RETURN
      END

! JKDERNOF5l
      SUBROUTINE JKDERNOF5l(KATOM,KTYPE,KLOC,KKMIN,KKMAX,QD,RO,GRADS,   &
                            ATMNAME,XINTS,SIZE_ENV,ENV,NAT,ATM,NBAS,    &
                            BAS,IGTYP,IPRINTOPT)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      LOGICAL FROZEN
      COMMON/INPNOF_FROZEN/FROZEN,IFROZEN(200)
      COMMON/INPNOF_PNOF/IPNOF,NTWOPAR
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      COMMON/INTOPT/CUTOFF,ISCHWZ,IECP,NECP
!     INPUT-OUTPUT VARIABLES OR ARGUMENTS
      INTEGER,DIMENSION(NSHELL),INTENT(IN)::KATOM,KTYPE,KLOC,KKMIN,KKMAX
      DOUBLE PRECISION,DIMENSION(NBF,NBF,NBF),INTENT(IN)::QD
      DOUBLE PRECISION,DIMENSION(NBF5),INTENT(IN)::RO
      DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::P,CJAUX,CKAUX
      DOUBLE PRECISION,DIMENSION((NSHELL*NSHELL+NSHELL)/2),INTENT(IN)::XINTS
      CHARACTER*4 ATMNAME(NATOMS)
      DOUBLE PRECISION,DIMENSION(3,NATOMS),INTENT(INOUT)::GRADS
      INTEGER,INTENT(IN)::IPRINTOPT
!     INTERMEDIATE VARIABLES
      INTEGER,DIMENSION(NBF)::IA
      INTEGER,DIMENSION(35)::IGXYZ,JGXYZ,KGXYZ,LGXYZ
      INTEGER,DIMENSION(:),ALLOCATABLE::IJKLG
      INTEGER::TOTCOUNT,INVTYP,KANG,I
      INTEGER(4) :: DI,DJ,DK,DL
      INTEGER::II,JJ,KK,LL,MAXLL,ISH,JSH,KSH,LSH,IIAT,JJAT,KKAT,LLAT
      INTEGER::MINJ,MAXJ,MINK,MAXK,MINL,MAXL,MINI,MAXI
      INTEGER::NUMI,NUMJ,NUMK,NUML
      LOGICAL::SKIPI,SKIPJ,SKIPK,SKIPL
      LOGICAL::IIEQJJ,KKEQLL,IJEQKL
      DOUBLE PRECISION::DABMAX,CUTOFF2,GMAX
      DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::DAB
      DOUBLE PRECISION,PARAMETER::ZERO=0.0D+00,ONE=1.0D+00
!     LIBCINT
      INTEGER :: SIZE_ENV,NBAS,IGTYP
      DOUBLE PRECISION :: ENV(SIZE_ENV)
      INTEGER :: ATM(6,NAT), BAS(8,NBAS)
!-----------------------------------------------------------------------
!     ROUTINE BASCHK DO:
      KANG=0
      DO N=1,NSHELL
        IF(KTYPE(N).GT.KANG) KANG = KTYPE(N)
      ENDDO
!
!     SET STARTING PARAMETERS
!
!     CUTOFF IS THE SCHWARZ SCREENING CUT OFF
      CUTOFF2 = CUTOFF*0.5D+00
!     SET POINTERS FOR PARTITIONING MEMORY
      DO I=1,NBF
        IA(I) = (I*I-I)/2
      ENDDO
!
!     MAXNUM=NUMBER OF FUNCTIONS WITH ANGULAR MOMENTUM EQUAL TO MAXTYP
      MAXNUM=((KANG)*(KANG+1))/2
!     DO AT LEAST AN L SHELL
      IF(MAXNUM.LT.4) MAXNUM=4
      MAXNUM=(MAXNUM**4)
!
      ALLOCATE(IJKLG(MAXNUM))
      ALLOCATE(DAB(MAXNUM))
      ALLOCATE(P(NBF,NBFT))
      CALL SQUARETRIAN3(QD,P,NBF,NBFT)
      P=P*0.5D+00
      ALLOCATE(CJAUX(NBFT,NBFT))
      ALLOCATE(CKAUX(NBFT,NBFT))
      IF(IPNOF==5) CALL DABNOF5PRE(RO,P,CJAUX,CKAUX)
      IF(IPNOF==7) CALL DABNOF7PRE(RO,P,CJAUX,CKAUX)
      DEALLOCATE(P)
!
      TOTCOUNT=1
!
!----I SHELL
!
      DO II=1,NSHELL
!
!-----J SHELL
!
       DO JJ=1,II
!
!-----K SHELL
!
        DO KK=1,II
!
!-----L SHELL
!
        MAXLL=KK
        IF(KK.EQ.II) MAXLL=JJ
         DO LL=1,MAXLL
!
!         IMPLEMENT INTEGRAL SCREENING HERE USING EXCHANGE INTEGRALS
!
          IJIJ=IA(MAX0(II,JJ))+MIN0(II,JJ)
          KLKL=IA(MAX0(KK,LL))+MIN0(KK,LL)
          GMAX=(XINTS(IJIJ)*XINTS(KLKL))
!
!         COARSE SCREENING, ON JUST THE INTEGRAL VALUE
!         ONLY WORKS IF SCHWARZ SCREENING IS ON (NATOMS>5)
          IF (GMAX.LT.CUTOFF.AND.NATOMS>5) CYCLE
!
          ISH=II
          JSH=JJ
          KSH=KK
          LSH=LL
!
!         SELECT CENTERS FOR DERIVATIVES
!
          CALL JKDATMNOF(ISH,JSH,KSH,LSH,SKIPI,SKIPJ,SKIPK,SKIPL,       &
                         INVTYP,KATOM,IIAT,JJAT,KKAT,LLAT)
          IF(SKIPI.AND.SKIPJ.AND.SKIPK.AND.SKIPL) CYCLE
!
!         SELECT INDICES FOR SHELL BLOCK
!
          CALL JKDSHLNOFl(ISH,JSH,KSH,LSH,IIEQJJ,KKEQLL,IJEQKL,DI,DJ,   &
                          DK,DL,SIZE_ENV,ENV,NAT,ATM,NBAS,BAS,IGTYP)
          CALL JKDNDXNOFl(DI,DJ,DK,DL,IGXYZ,JGXYZ,KGXYZ,LGXYZ,IIEQJJ,   &
                          KKEQLL,IJEQKL,IJKLG,MAXNUM,ISH,JSH,KSH,LSH,   &
                          NBF,NSHELL,SIZE_ENV,ENV,NAT,ATM,NBAS,BAS,     &
                          IGTYP)
!
!         OBTAIN 2e DENSITY FOR THIS SHELL BLOCK
!
          CALL DABNOF5l(ISH,JSH,KSH,LSH,IGXYZ,JGXYZ,KGXYZ,LGXYZ,        &
                        IIEQJJ,KKEQLL,IJEQKL,IA,DAB,MAXNUM,DABMAX,      &
                        CJAUX,CKAUX,NAT,ATM,NBAS,BAS,IGTYP)
!
!         FINE SCREENING, ON INTEGRAL VALUE TIMES DENSITY FACTOR
!         ONLY WORKS IF SCHWARZ SCREENING IS ON (NATOMS>5)
          IF(DABMAX*GMAX.LT.CUTOFF2.AND.NATOMS>5) CYCLE
!
!         EVALUATE DERIVATIVE INTEGRAL AND ADD TO THE GRADIENT
!
          CALL JKDSPDNOFl(TOTCOUNT,DAB,II,JJ,KK,LL,IIEQJJ,KKEQLL,       &
                          IJEQKL,SKIPI,SKIPJ,SKIPK,SKIPL,IJKLG,         &
                          MAXNUM,INVTYP,IIAT,JJAT,KKAT,LLAT,            &
                          SIZE_ENV,ENV,NAT,ATM,NBAS,BAS,IGTYP,GRADS)
!
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      DEALLOCATE(IJKLG,DAB,CJAUX,CKAUX)
!
!     MAKE ZERO GRADIENT CORRESPONDING TO FROZEN COORDINATES
!
      IF(FROZEN) THEN
        DO I=1,200,2
         IF(IFROZEN(I).EQ.0) EXIT
         GRADS(IFROZEN(I),IFROZEN(I+1))=ZERO
        ENDDO
      ENDIF
!
!     PRINT OUT TOTAL PNOF ENERGY GRADIENT
!
      IF(IPRINTOPT==1)THEN
       WRITE(6,1)
       DO I=1,NATOMS
        WRITE(6,2)I,ATMNAME(I),GRADS(1,I),GRADS(2,I),GRADS(3,I)
       ENDDO
      ENDIF
!-----------------------------------------------------------------------
    1 FORMAT( /1X,'----------------',                                  &
              /1X,' Total Gradient ',                                  &
              /1X,'----------------',                                  &
              //9X,'Atom',7X,'Ex',10X,'Ey',10X,'Ez' )
    2 FORMAT(/1X,I4,5X,A4,F10.4,2X,F10.4,2X,F10.4)
!-----------------------------------------------------------------------
      RETURN
      END

! JKDATMNOF
      SUBROUTINE JKDATMNOF(II,JJ,KK,LL,SKIPI,SKIPJ,SKIPK,SKIPL,INVTYP,  &
                        KATOM,IIAT,JJAT,KKAT,LLAT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      INTEGER,DIMENSION(NSHELL),INTENT(IN)::KATOM
      INTEGER,INTENT(IN)::II,JJ,KK,LL
      INTEGER,INTENT(OUT)::INVTYP,IIAT,JJAT,KKAT,LLAT
      LOGICAL,INTENT(OUT)::SKIPI,SKIPJ,SKIPK,SKIPL
      LOGICAL::IANDJ,IANDK,IANDL,JANDK,JANDL,KANDL
!-----------------------------------------------------------------------
      SKIPI=.TRUE.
      SKIPJ=.TRUE.
      SKIPK=.TRUE.
      SKIPL=.TRUE.
      IIAT=KATOM(II)
      JJAT=KATOM(JJ)
      KKAT=KATOM(KK)
      LLAT=KATOM(LL)
      IANDJ=IIAT.EQ.JJAT
      IANDK=IIAT.EQ.KKAT
      IANDL=IIAT.EQ.LLAT
      JANDK=JJAT.EQ.KKAT
      JANDL=JJAT.EQ.LLAT
      KANDL=KKAT.EQ.LLAT
      IF(.NOT.IANDJ) GO TO 500
      IF(.NOT.IANDK) GO TO 200
      IF(.NOT.IANDL) GO TO 100
!     ----- IAT = JAT = KAT = LAT ----- (IAT,IAT/IAT,IAT) -----
      INVTYP=1
      GO TO 1500
  100 CONTINUE
!     ----- IAT = JAT = KAT ; LAT ----- (IAT,IAT/IAT,LAT) -----
      SKIPL=.FALSE.
      INVTYP=2
      GO TO 1500
  200 IF(.NOT.IANDL) GO TO 300
!     ----- IAT = JAT = LAT ; KAT ----- (IAT,IAT/KAT,IAT) -----
      SKIPK=.FALSE.
      INVTYP=3
      GO TO 1500
  300 IF(.NOT.KANDL) GO TO 400
!     ----- IAT = JAT ; KAT = LAT ----- (IAT,IAT/KAT,KAT) -----
      SKIPK=.FALSE.
      SKIPL=.FALSE.
      INVTYP=4
      GO TO 1500
  400 CONTINUE
!     ----- IAT = JAT ; KAT ; LAT ----- (IAT,IAT/KAT,LAT) -----
      SKIPK=.FALSE.
      SKIPL=.FALSE.
      INVTYP=5
      GO TO 1500
  500 IF(.NOT.IANDK) GO TO 800
      IF(.NOT.IANDL) GO TO 600
!     ----- IAT = KAT = LAT ; JAT ----- (IAT,JAT/IAT,IAT) -----
      SKIPJ=.FALSE.
      INVTYP=6
      GO TO 1500
  600 IF(.NOT.JANDL) GO TO 700
!     ----- IAT = KAT ; JAT = LAT ----- (IAT,JAT/IAT,JAT) -----
      SKIPJ=.FALSE.
      SKIPL=.FALSE.
      INVTYP=7
      GO TO 1500
  700 CONTINUE
!     ----- IAT = KAT ; JAT ; LAT ----- (IAT,JAT/IAT,LAT) -----
      SKIPJ=.FALSE.
      SKIPL=.FALSE.
      INVTYP=8
      GO TO 1500
  800 IF(.NOT.IANDL) GO TO 1000
      IF(.NOT.JANDK) GO TO 900
!     ----- IAT = LAT ; JAT = KAT ----- (IAT,JAT/JAT,IAT) -----
      SKIPJ=.FALSE.
      SKIPK=.FALSE.
      INVTYP=9
      GO TO 1500
  900 CONTINUE
!     ----- IAT = LAT ; JAT , KAT ----- (IAT,JAT/KAT,IAT) -----
      SKIPJ=.FALSE.
      SKIPK=.FALSE.
      INVTYP=10
      GO TO 1500
 1000 IF(.NOT.JANDK) GO TO 1200
      IF(.NOT.JANDL) GO TO 1100
!     ----- IAT ; JAT = JAT = JAT ----- (IAT,JAT/JAT,JAT) -----
      SKIPI=.FALSE.
      INVTYP=11
      GO TO 1500
 1100 CONTINUE
!     ----- IAT ; JAT = KAT ; LAT ----- (IAT,JAT/JAT,LAT) -----
      SKIPI=.FALSE.
      SKIPL=.FALSE.
      INVTYP=12
      GO TO 1500
 1200 IF(.NOT.JANDL) GO TO 1300
!     ----- JAT = LAT ; IAT ; KAT ----- (IAT,JAT/KAT,JAT) -----
      SKIPI=.FALSE.
      SKIPK=.FALSE.
      INVTYP=13
      GO TO 1500
 1300 IF(.NOT.KANDL) GO TO 1400
!     ----- KAT = LAT ; IAT ; JAT ----- (IAT,JAT/KAT,KAT) -----
      SKIPI=.FALSE.
      SKIPJ=.FALSE.
      INVTYP=14
      GO TO 1500
 1400 CONTINUE
!     ----- IAT ; JAT ; KAT ; LAT ----- (IAT,JAT/KAT,LAT) -----
      SKIPI=.FALSE.
      SKIPJ=.FALSE.
      SKIPK=.FALSE.
      INVTYP=15
 1500 CONTINUE
!-----------------------------------------------------------------------
      RETURN
      END

! JKDSHLNOFl
      SUBROUTINE JKDSHLNOFl(ISH,JSH,KSH,LSH,IIEQJJ,KKEQLL,IJEQKL,DI,DJ, &
                            DK,DL,SIZE_ENV,ENV,NAT,ATM,NBAS,BAS,IGTYP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      INTEGER,INTENT(IN)::ISH,JSH,KSH,LSH
      LOGICAL,INTENT(OUT)::IIEQJJ,KKEQLL,IJEQKL

      INTEGER(4) :: DI, DJ, DK, DL, I, J, K, L
      INTEGER(4) :: DJEFF, DKEFF, DLEFF

      INTEGER :: SIZE_ENV, NAT, NBAS, IGTYP
      DOUBLE PRECISION :: ENV(SIZE_ENV)
      INTEGER :: IA(NBF), LOC(NSHELL)

      INTEGER :: ATM(6, NAT), BAS(8, NBAS)

      DOUBLE PRECISION, ALLOCATABLE :: BLK(:,:,:)

      INTEGER, EXTERNAL :: CINTcgto_cart
      INTEGER, EXTERNAL :: CINTcgto_spheric
!-----------------------------------------------------------------------
      DABMAX=ZER         
!
!     GET 2e DENSITY FOR THIS SHELL BLOCK
!
!-----------------------------------------------------------------------
! LOC indicates the starting location of each shell in the AO basis
!-----------------------------------------------------------------------
      IF(IGTYP==1) THEN
        DI = CINTcgto_cart(ISH-1, BAS)
        DJ = CINTcgto_cart(JSH-1, BAS)
        DK = CINTcgto_cart(KSH-1, BAS)
        DL = CINTcgto_cart(LSH-1, BAS)
      ELSEIF(IGTYP==2) THEN
        DI = CINTcgto_spheric(ISH-1, BAS)
        DJ = CINTcgto_spheric(JSH-1, BAS)
        DK = CINTcgto_spheric(KSH-1, BAS)
        DL = CINTcgto_spheric(LSH-1, BAS)
      END IF
      
!-----------------------------------------------------------------------
      IIEQJJ=ISH.EQ.JSH
      KKEQLL=KSH.EQ.LSH
      IJEQKL=ISH.EQ.KSH.AND.JSH.EQ.LSH
!-----------------------------------------------------------------------
      RETURN
      END

! JKDNDXNOFl
      SUBROUTINE JKDNDXNOFl(DI,DJ,DK,DL,IGXYZ,JGXYZ,KGXYZ,              &
                            LGXYZ,IIEQJJ,KKEQLL,IJEQKL,IJKLG,MAXNUM,    &
                            II,JJ,KK,LL,NBF,NSHELL,SIZE_ENV,ENV,NAT,    &
                            ATM,NBAS,BAS,IGTYP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      INTEGER,DIMENSION(MAXNUM)::IJKLG
      INTEGER,DIMENSION(35),INTENT(OUT)::IGXYZ,JGXYZ,KGXYZ,LGXYZ
      LOGICAL,INTENT(IN)::IIEQJJ,KKEQLL,IJEQKL

      INTEGER(4) :: DI, DJ, DK, DL, I, J, K, L
      INTEGER(4) :: DJEFF, DKEFF, DLEFF

      INTEGER :: SIZE_ENV, NAT, NBAS
      DOUBLE PRECISION :: ENV(SIZE_ENV)
      INTEGER :: IA(NBF), LOC(NSHELL)

      INTEGER :: ATM(6, NAT), BAS(8, NBAS)

      DOUBLE PRECISION, ALLOCATABLE :: BLK(:,:,:)

      INTEGER, EXTERNAL :: CINTcgto_cart
      INTEGER, EXTERNAL :: CINTcgto_spheric
!-----------------------------------------------------------------------
      DABMAX=ZER         
!
!     GET 2e DENSITY FOR THIS SHELL BLOCK
!
!-----------------------------------------------------------------------
! LOC indicates the starting location of each shell in the AO basis
!-----------------------------------------------------------------------
!     ----- PREPARE INDICES FOR PAIRS OF (I,J) FUNCTIONS -----
      NI=DL*DK*DJ
      DO I=1,DI
        IGXYZ(I)=NI*(I-1)+1
      ENDDO

      NJ=DL*DK
      DO J=1,DJ
        JGXYZ(J)=NJ*(J-1)
      ENDDO
!     ----- PREPARE INDICES FOR PAIRS OF (K,L) FUNCTIONS -----
      NK=DL
      DO K=1,DK
        KGXYZ(K)=NK*(K-1)
      ENDDO
      NL=1
      DO L=1,DL
         LGXYZ(L)=NL*(L-1)
      ENDDO
!     ----- PREPARE INDICES FOR (IJ/KL) -----
      IJKL=0
      DO I=1,DI
        DJEFF = MERGE(I, DJ, IIEQJJ)
        DO J=1,DJEFF
          DKEFF = MERGE(I, DK, IJEQKL)
          DO K=1,DKEFF
            DLEFF = MERGE(K, DL, KKEQLL)
            DLEFF = MERGE(J, DLEFF, IJEQKL .AND. K.EQ.I)
            DO L=1,DLEFF
              IJKL=IJKL+1
              NN=((IGXYZ(I)+JGXYZ(J))+KGXYZ(K))+LGXYZ(L)
              IJKLG(IJKL)=   NN
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! DABNOF2PRE
      SUBROUTINE DABNOF2PRE(CJ12,CK12,RO,DA,DAAUX,DAAUX2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
!
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5),INTENT(IN)::CJ12,CK12
      DOUBLE PRECISION,DIMENSION(NBF5),INTENT(IN)::RO
      DOUBLE PRECISION,DIMENSION(NBF,NBFT),INTENT(IN)::DA
      DOUBLE PRECISION,DIMENSION(NBF,NBFT),INTENT(OUT)::DAAUX,DAAUX2
      DOUBLE PRECISION,PARAMETER::ZER=0.0D+00
!-----------------------------------------------------------------------
!
!     FIRST CONTRACT OVER ONE NATURAL ORBITAL COEFFICIENT
!
!-----------------------------------------------------------------------
      DAAUX=ZER
      DAAUX2=ZER
       DO IJ=1,NBFT
        DO LP=1,NBF5
         if(LP<=NB.or.LP>NA)DAAUX(LP,IJ)=DAAUX(LP,IJ)+RO(LP)*DA(LP,IJ)
         DO LQ=1,LP-1
           DAAUX(LQ,IJ) = DAAUX(LQ,IJ) + CJ12(LP,LQ)*DA(LP,IJ)
           DAAUX2(LQ,IJ) = DAAUX2(LQ,IJ) + CK12(LP,LQ)*DA(LP,IJ)
         ENDDO
         DO LQ=LP+1,NBF5
           DAAUX(LQ,IJ) = DAAUX(LQ,IJ) + CJ12(LP,LQ)*DA(LP,IJ)
           DAAUX2(LQ,IJ) = DAAUX2(LQ,IJ) + CK12(LP,LQ)*DA(LP,IJ)
         ENDDO
        ENDDO
       ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! DABNOF5PRE
      SUBROUTINE DABNOF5PRE(RO,DA,CJAUX,CKAUX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      COMMON/INPNOF_ORBSPACE0/NO1,NDOC,NCO,NCWO,NVIR,NAC,NO0
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Extended (Nc>1): NBF5 = NBF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DOUBLE PRECISION,DIMENSION(NBF5),INTENT(IN)::RO
      DOUBLE PRECISION,DIMENSION(NBF5)::BETA
      DOUBLE PRECISION,DIMENSION(NBF,NBFT),INTENT(IN)::DA
      DOUBLE PRECISION,DIMENSION(NBFT)::DENS5
      DOUBLE PRECISION,DIMENSION(NDOC,NBFT)::DAAUX2
      DOUBLE PRECISION,DIMENSION(NCO,NBFT)::DAAUX
      DOUBLE PRECISION,DIMENSION(NBFT,NBFT),INTENT(OUT)::CJAUX,CKAUX
      DOUBLE PRECISION,PARAMETER::ZER=0.0D+00
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO i=1,NBF5
        BETA(i)=DSQRT(RO(i))
      ENDDO
!     DENS5 AND DAAUX ARE RELATED WITH (NN-DELTA)
!     DAAUX2 IS RELATED TO PI MATRIX
      DENS5=ZER
      DAAUX=ZER
      DAAUX2=ZER
      CJAUX = ZER
      CKAUX = ZER
      DO IJ=1,NBFT
!     ADD HF-LIKE TERMS FROM FROZEN COEFFICIENTS
       DO LP=1,NO1
          DENS5(IJ) = DENS5(IJ) + RO(LP)*DA(LP,IJ)
          DAAUX(LP,IJ) = DAAUX(LP,IJ) + RO(LP)*DA(LP,IJ)
       ENDDO
       DO j=1,NDOC
        jn = NO1+j
        DENS5(IJ)=DENS5(IJ)+RO(jn)*DA(jn,IJ)
        DAAUX(jn,IJ)=DAAUX(jn,IJ)+RO(jn)*DA(jn,IJ)
        DAAUX2(j,IJ)=DAAUX2(j,IJ)-BETA(jn)*DA(jn,IJ)
        DO i=NDOC+NCWO*(NDOC-j)+1,NDOC+NCWO*(NDOC-j+1)
          in = NO1+i
          DENS5(IJ)=DENS5(IJ)+RO(in)*DA(in,IJ)
          DAAUX(jn,IJ)=DAAUX(jn,IJ)+RO(in)*DA(in,IJ)
          DAAUX2(j,IJ)=DAAUX2(j,IJ)+BETA(in)*DA(in,IJ)
        ENDDO
       ENDDO
!
       DO KL=1,IJ
        CJAUX(IJ,KL) = DENS5(IJ)*DENS5(KL)
        DO IG=1,NDOC
         IGG = IG + NO1
         CJAUX(IJ,KL) = CJAUX(IJ,KL) - DAAUX(IGG,IJ)*DAAUX(IGG,KL)
         CKAUX(IJ,KL) = CKAUX(IJ,KL) + DAAUX2(IG,IJ)*DAAUX2(IG,KL)
        ENDDO
!       AACOMP
        CKAUX(IJ,KL) = CJAUX(IJ,KL) - CKAUX(IJ,KL)
        CJAUX(IJ,KL) = CJAUX(IJ,KL) + CJAUX(IJ,KL)
!       SYMMETRY IN AO INDICES
        CKAUX(KL,IJ) = CKAUX(IJ,KL)
        CJAUX(KL,IJ) = CJAUX(IJ,KL)
       ENDDO
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! DABNOF7PRE
      SUBROUTINE DABNOF7PRE(RO,DA,CJAUX,CKAUX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_STATIC/Ista
      COMMON/INPNOF_ORBSPACE0/NO1,NDOC,NCO,NCWO,NVIR,NAC,NO0
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Extended (Nc>1): NBF5 = NBF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DOUBLE PRECISION,DIMENSION(NBF5),INTENT(IN)::RO
      DOUBLE PRECISION,DIMENSION(NBF5)::BETA
      DOUBLE PRECISION,DIMENSION(NBF,NBFT),INTENT(IN)::DA
      DOUBLE PRECISION,DIMENSION(NBFT)::DENS5,DENS7
      DOUBLE PRECISION,DIMENSION(NDOC,NBFT)::DAAUX2,DAAUX3
      DOUBLE PRECISION,DIMENSION(NCO,NBFT)::DAAUX
      DOUBLE PRECISION,DIMENSION(NBFT,NBFT),INTENT(OUT)::CJAUX,CKAUX
      DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::FIs
      DOUBLE PRECISION,PARAMETER::ZER=0.0D+00
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO i=1,NBF5
        BETA(i)=DSQRT(RO(i))
      ENDDO
      ALLOCATE(FIs(NBF5))
      FIs = 0.0d0
      if(Ista==0)then
       DO j=NO1+1,NBF5
        FIs(j) = DSQRT( RO(j)*(1.0d0-RO(j)) )
       ENDDO
      else if(Ista==1)then
       DO j=NO1+1,NBF5
        FIs(j) = 2.0d0*RO(j)*(1.0d0-RO(j))
       ENDDO
      end if
!
      DENS5=ZER
      DENS7=ZER
      DAAUX=ZER
      DAAUX2=ZER
      DAAUX3=ZER
      CJAUX = ZER
      CKAUX = ZER
      DO IJ=1,NBFT
       DO LP=1,NO1
          DENS5(IJ) = DENS5(IJ) + RO(LP)*DA(LP,IJ)
          DAAUX(LP,IJ) = DAAUX(LP,IJ) + RO(LP)*DA(LP,IJ)
       ENDDO
       DO j=1,NDOC
        jn = NO1+j
        DENS5(IJ)=DENS5(IJ)+RO(jn)*DA(jn,IJ)
        DENS7(IJ)=DENS7(IJ)+FIs(jn)*DA(jn,IJ)
        DAAUX(jn,IJ)=DAAUX(jn,IJ)+RO(jn)*DA(jn,IJ)
        DAAUX2(j,IJ)=DAAUX2(j,IJ)-BETA(jn)*DA(jn,IJ)
        DAAUX3(j,IJ)=DAAUX3(j,IJ)+FIs(jn)*DA(jn,IJ)
        DO i=NDOC+NCWO*(NDOC-j)+1,NDOC+NCWO*(NDOC-j+1)
          in = NO1+i
          DENS5(IJ)=DENS5(IJ)+RO(in)*DA(in,IJ)
          DENS7(IJ)=DENS7(IJ)+FIs(in)*DA(in,IJ)
          DAAUX(jn,IJ)=DAAUX(jn,IJ)+RO(in)*DA(in,IJ)
          DAAUX2(j,IJ)=DAAUX2(j,IJ)+BETA(in)*DA(in,IJ)
          DAAUX3(j,IJ)=DAAUX3(j,IJ)+FIs(in)*DA(in,IJ)
        ENDDO
       ENDDO
!
       DO KL=1,IJ
        CJAUX(IJ,KL) = DENS5(IJ)*DENS5(KL)
        CKAUX(IJ,KL) = - DENS7(IJ)*DENS7(KL)
        DO IG=1,NDOC
         IGG = IG + NO1
         CJAUX(IJ,KL) = CJAUX(IJ,KL) - DAAUX(IGG,IJ)*DAAUX(IGG,KL)
         CKAUX(IJ,KL) = CKAUX(IJ,KL) + DAAUX2(IG,IJ)*DAAUX2(IG,KL)
         CKAUX(IJ,KL) = CKAUX(IJ,KL) + DAAUX3(IG,IJ)*DAAUX3(IG,KL)
        ENDDO
!       AACOMP
        CKAUX(IJ,KL) = CJAUX(IJ,KL) - CKAUX(IJ,KL)
        CJAUX(IJ,KL) = CJAUX(IJ,KL) + CJAUX(IJ,KL)
!       SYMMETRY IN AO INDICES
        CKAUX(KL,IJ) = CKAUX(IJ,KL)
        CJAUX(KL,IJ) = CJAUX(IJ,KL)
       ENDDO
      ENDDO
      DEALLOCATE(FIs)
!-----------------------------------------------------------------------
      RETURN
      END

! DABNOF2l
      SUBROUTINE DABNOF2l(II,JJ,KK,LL,IGXYZ,JGXYZ,KGXYZ,LGXYZ,IIEQJJ,   &
                          KKEQLL,IJEQKL,IA,DA,DAB,MAXNUM,DABMAX,DAAUX,  &
                          DAAUX2,NAT,ATM,NBAS,BAS,IGTYP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
!
      INTEGER,INTENT(IN) :: II,JJ,KK,LL,MAXNUM,NAT,NBAS,IGTYP
      INTEGER,DIMENSION(35),INTENT(IN) :: IGXYZ,JGXYZ,KGXYZ,LGXYZ
      LOGICAL,INTENT(IN) :: IIEQJJ,KKEQLL,IJEQKL
      INTEGER,DIMENSION(NBF),INTENT(IN) :: IA
      DOUBLE PRECISION,DIMENSION(NBF,NBFT),INTENT(IN) :: DA,DAAUX,DAAUX2
      DOUBLE PRECISION,DIMENSION(MAXNUM),INTENT(OUT) :: DAB
      DOUBLE PRECISION,INTENT(OUT) :: DABMAX
      DOUBLE PRECISION,PARAMETER :: ZER=0.0D0, PT5=0.5D0


      INTEGER(4) :: DI, DJ, DK, DL, I ,J, K, L
      INTEGER(4) :: DJEFF, DKEFF, DLEFF

      INTEGER :: LOC(NSHELL)

      INTEGER :: ATM(6, NAT), BAS(8, NBAS)

      INTEGER, EXTERNAL :: CINTcgto_cart
      INTEGER, EXTERNAL :: CINTcgto_spheric
!-----------------------------------------------------------------------
      DABMAX=ZER         
!
!     GET 2e DENSITY FOR THIS SHELL BLOCK
!
!-----------------------------------------------------------------------
! LOC indicates the starting location of each shell in the AO basis
!-----------------------------------------------------------------------
      LOC(1) = 1
      IF(IGTYP==1) DI = CINTcgto_cart(0, BAS)
      IF(IGTYP==2) DI = CINTcgto_spheric(0, BAS)
      DO ISH = 2,NSHELL
        LOC(ISH) = LOC(ISH-1) + DI
        IF(IGTYP==1) DI = CINTcgto_cart(ISH-1, BAS)
        IF(IGTYP==2) DI = CINTcgto_spheric(ISH-1, BAS)
      END DO

      IF(IGTYP==1) THEN
        DI = CINTcgto_cart(II-1, BAS)
        DJ = CINTcgto_cart(JJ-1, BAS)
        DK = CINTcgto_cart(KK-1, BAS)
        DL = CINTcgto_cart(LL-1, BAS)
     ELSE IF(IGTYP==2) THEN
        DI = CINTcgto_spheric(II-1, BAS)
        DJ = CINTcgto_spheric(JJ-1, BAS)
        DK = CINTcgto_spheric(KK-1, BAS)
        DL = CINTcgto_spheric(LL-1, BAS)
      END IF

      LOCI=LOC(II) - 1
      LOCJ=LOC(JJ) - 1
      LOCK=LOC(KK) - 1
      LOCL=LOC(LL) - 1

      DO I=1,DI
        DJEFF = MERGE(I, DJ, IIEQJJ)
        DO J=1,DJEFF
          LI = LOCI + I
          LJ = LOCJ + J
          IAJ= MAX0(LI,LJ)
          IIJ= MIN0(LI,LJ)
          DKEFF = MERGE(I, DK, IJEQKL)
          DO K=1,DKEFF
            DLEFF = MERGE(K, DL, KKEQLL)
            DLEFF = MERGE(J, DLEFF, IJEQKL .AND. K.EQ.I)
            DO L=1,DLEFF
              LK = LOCK + K
              LLL = LOCL + L
              KAL= MAX0(LK,LLL)
              KIL= MIN0(LK,LLL)
              IN = IAJ
              JN = IIJ
              KN = KAL
              LN = KIL
              IF(IN.LT.KN .OR.(IN.EQ.KN .AND. JN.LT.LN)) THEN
                IN = KAL
                JN = KIL
                KN = IAJ
                LN = IIJ
              ENDIF
              IJ = IA(IN)+JN
              IK = IA(IN)+KN
              IL = IA(IN)+LN
              JK = IA(MAX0(JN,KN))+MIN0(JN,KN)
              JL = IA(JN)+LN
              IF(JN.LT.KN) JL = IA(MAX0(JN,LN))+MIN0(JN,LN)
              KL = IA(KN)+LN
!             CONTRACT OVER THE 2nd NO COEFFICIENT
              DF1=ZER
              DQ1=ZER
              DO LQ=1,NBF5
                DF1 = DF1 + DAAUX(LQ,KL)*DA(LQ,IJ)
                DQ1 = DQ1 + DAAUX2(LQ,JK)*DA(LQ,IL)
                DQ1 = DQ1 + DAAUX2(LQ,JL)*DA(LQ,IK)
              ENDDO
!             BUILD THE DENSITY TERM SUMMING COULOMB-LIKE
!             AND EXCHANGE-LIKE PARTS
              DF1 = DF1 + DF1 - DQ1
!             AVOID DOUBLE COUNTING OF DIAGONAL TERMS                     
              IF(JN.EQ.IN               ) DF1= DF1*PT5
              IF(LN.EQ.KN               ) DF1= DF1*PT5
              IF(KN.EQ.IN .AND. LN.EQ.JN) DF1= DF1*PT5
              IF(DABMAX.LT. ABS(DF1)) DABMAX= ABS(DF1)
!             IGXYZ AND J, K, AND L ARE SET UP IN JKDNDX
              IJKL=IGXYZ(I)+JGXYZ(J)+KGXYZ(K)+LGXYZ(L)
              DAB(IJKL)= DF1
              ENDDO
            ENDDO
          ENDDO
        ENDDO
!-----------------------------------------------------------------------
      RETURN
      END        

! DABNOF5l
      SUBROUTINE DABNOF5l(II,JJ,KK,LL,IGXYZ,JGXYZ,KGXYZ,LGXYZ,IIEQJJ,   &
                          KKEQLL,IJEQKL,IA,DAB,MAXNUM,DABMAX,CJAUX,     &
                          CKAUX,NAT,ATM,NBAS,BAS,IGTYP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      INTEGER,INTENT(IN)::II,JJ,KK,LL,IGTYP
      INTEGER(4) :: DI, DJ, DK, DL, I ,J, K, L
      INTEGER(4) :: DJEFF, DKEFF, DLEFF
      INTEGER,DIMENSION(35),INTENT(IN)::IGXYZ,JGXYZ,KGXYZ,LGXYZ
      INTEGER,DIMENSION(NBF),INTENT(IN)::IA
      INTEGER,INTENT(IN)::MAXNUM
      LOGICAL,INTENT(IN)::IIEQJJ,KKEQLL,IJEQKL
      DOUBLE PRECISION,DIMENSION(NBFT,NBFT),INTENT(IN)::CJAUX,CKAUX
      DOUBLE PRECISION,INTENT(OUT)::DABMAX
      DOUBLE PRECISION,DIMENSION(MAXNUM),INTENT(OUT)::DAB
      DOUBLE PRECISION,PARAMETER::ZER=0.0D+00,PT5=0.5D+00
      INTEGER :: NAT, NBAS
      INTEGER(4) :: LOC(NSHELL)

      INTEGER :: ATM(6, NAT), BAS(8, NBAS)

      INTEGER, EXTERNAL :: CINTcgto_cart
      INTEGER, EXTERNAL :: CINTcgto_spheric
!-----------------------------------------------------------------------
      DABMAX=ZER         
!
!     GET 2e DENSITY FOR THIS SHELL BLOCK
!
!-----------------------------------------------------------------------
! LOC indicates the starting location of each shell in the AO basis
!-----------------------------------------------------------------------
      LOC(1) = 1
      IF(IGTYP==1) DI = CINTcgto_cart(0, BAS)
      IF(IGTYP==2) DI = CINTcgto_spheric(0, BAS)
      DO ISH = 2,NSHELL
        LOC(ISH) = LOC(ISH-1) + DI
        IF(IGTYP==1) DI = CINTcgto_cart(ISH-1, BAS)
        IF(IGTYP==2) DI = CINTcgto_spheric(ISH-1, BAS)
      END DO

      IF(IGTYP==1) THEN
        DI = CINTcgto_cart(II-1, BAS)
        DJ = CINTcgto_cart(JJ-1, BAS)
        DK = CINTcgto_cart(KK-1, BAS)
        DL = CINTcgto_cart(LL-1, BAS)
      ELSE IF(IGTYP==2) THEN
        DI = CINTcgto_spheric(II-1, BAS)
        DJ = CINTcgto_spheric(JJ-1, BAS)
        DK = CINTcgto_spheric(KK-1, BAS)
        DL = CINTcgto_spheric(LL-1, BAS)
      END IF

      LOCI=LOC(II) - 1
      LOCJ=LOC(JJ) - 1
      LOCK=LOC(KK) - 1
      LOCL=LOC(LL) - 1
      DO I=1,DI
        DJEFF = MERGE(I, DJ, IIEQJJ)
        DO J=1,DJEFF
          LI = LOCI + I
          LJ = LOCJ + J
          IAJ= MAX0(LI,LJ)
          IIJ= MIN0(LI,LJ)
          DKEFF = MERGE(I, DK, IJEQKL)
          DO K=1,DKEFF
            DLEFF = MERGE(K, DL, KKEQLL)
            DLEFF = MERGE(J, DLEFF, IJEQKL .AND. K.EQ.I)
            DO L=1,DLEFF
              LK = LOCK + K
              LLL = LOCL + L
              KAL= MAX0(LK,LLL)
              KIL= MIN0(LK,LLL)
              IN = IAJ
              JN = IIJ
              KN = KAL
              LN = KIL
              IF(IN.LT.KN .OR.(IN.EQ.KN .AND. JN.LT.LN)) THEN
                IN = KAL
                JN = KIL
                KN = IAJ
                LN = IIJ
              ENDIF
              IJ = IA(IN)+JN
              IK = IA(IN)+KN
              IL = IA(IN)+LN
              JK = IA(MAX0(JN,KN))+MIN0(JN,KN)
              JL = IA(JN)+LN
              IF(JN.LT.KN) JL = IA(MAX0(JN,LN))+MIN0(JN,LN)
              KL = IA(KN)+LN

!             IGXYZ AND J, K, AND L ARE SET UP IN JKDNDX
              IJKL=IGXYZ(I)+JGXYZ(J)+KGXYZ(K)+LGXYZ(L)
              DAB(IJKL) = 2.0D0 * CJAUX(IJ,KL) - (CKAUX(JK,IL) + CKAUX(JL,IK))

!             AVOID DOUBLE COUNTING OF DIAGONAL TERMS                     
              IF (JN .EQ. IN) DAB(IJKL) = DAB(IJKL) * PT5
              IF (LN .EQ. KN) DAB(IJKL) = DAB(IJKL) * PT5
              IF (KN .EQ. IN .AND. LN .EQ. JN) DAB(IJKL) = DAB(IJKL) * PT5
              IF (ABS(DAB(IJKL)) > DABMAX) DABMAX = ABS(DAB(IJKL))

              ENDDO
            ENDDO
          ENDDO
        ENDDO
!-----------------------------------------------------------------------
      RETURN
      END      

! JKDSPDNOFl
      SUBROUTINE JKDSPDNOFl(TOTCOUNT,DAB,II,JJ,KK,LL,IIEQJJ,KKEQLL,     &
                            IJEQKL,SKIPI,SKIPJ,SKIPK,SKIPL,IJKLG,       &
                            MAXNUM,INVTYP,IIAT,JJAT,KKAT,LLAT,          &
                            SIZE_ENV,ENV,NAT,ATM,NBAS,BAS,IGTYP,GRADS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      INTEGER,INTENT(IN)::MAXNUM,IGTYP
      DOUBLE PRECISION,DIMENSION(MAXNUM),INTENT(IN)::DAB
      LOGICAL,INTENT(IN)::IIEQJJ,KKEQLL,IJEQKL
      LOGICAL,INTENT(IN)::SKIPI,SKIPJ,SKIPK,SKIPL
      INTEGER,INTENT(IN)::INVTYP,IIAT,JJAT,KKAT,LLAT
      INTEGER,INTENT(INOUT)::TOTCOUNT
      DOUBLE PRECISION,DIMENSION(3,NAT),INTENT(INOUT)::GRADS
      DOUBLE PRECISION,DIMENSION(12)::FD
      INTEGER,DIMENSION(MAXNUM),INTENT(IN)::IJKLG

      INTEGER :: SIZE_ENV,NBAS
      DOUBLE PRECISION :: ENV(SIZE_ENV)
      INTEGER :: ATM(6,NAT), BAS(8,NBAS)

      FD = 0.0d0
      CALL DSPDFSNOFl(TOTCOUNT,IJKLG,DAB,SKIPI,SKIPJ,SKIPK,SKIPL,       &
                      MAXNUM,IIEQJJ,KKEQLL,IJEQKL,FD,II,JJ,KK,LL,       &
                      SIZE_ENV,ENV,NAT,ATM,NBAS,BAS,IGTYP)
      CALL JKDINVNOF(INVTYP,FD,GRADS,IIAT,JJAT,KKAT,LLAT)
!-----------------------------------------------------------------------
      RETURN
      END

! JKDINVNOF
      SUBROUTINE JKDINVNOF(INVTYP,FD,GRADS,II,JJ,KK,LL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      INTEGER,INTENT(IN)::INVTYP,II,JJ,KK,LL
      DOUBLE PRECISION,DIMENSION(3,NATOMS),INTENT(INOUT)::GRADS
      DOUBLE PRECISION,DIMENSION(3,4)::FD
!     ----- TRANSLATIONAL INVARIANCE FOR GRADIENT ELEMENTS -----
      IF (INVTYP.EQ.2) THEN
        DO IXYZ=1,3
          FD(IXYZ,1)=- FD(IXYZ,4)
        ENDDO
      ELSE IF (INVTYP.EQ.3) THEN
        DO IXYZ=1,3
          FD(IXYZ,1)=- FD(IXYZ,3)
        ENDDO
      ELSE IF (INVTYP.EQ.4.OR.INVTYP.EQ.5) THEN
        DO IXYZ=1,3
          FD(IXYZ,1)=-(FD(IXYZ,3)+FD(IXYZ,4))
        ENDDO
      ELSE IF (INVTYP.EQ.6) THEN
        DO IXYZ=1,3
          FD(IXYZ,1)=- FD(IXYZ,2)
        ENDDO
      ELSE IF (INVTYP.EQ.7) THEN
        DO IXYZ=1,3
          FD(IXYZ,1)=-(FD(IXYZ,2)+FD(IXYZ,4))
        ENDDO
      ELSE IF (INVTYP.EQ.8) THEN
        DO IXYZ=1,3
          FD(IXYZ,1)=-(FD(IXYZ,2)+FD(IXYZ,4))
        ENDDO
      ELSE IF (INVTYP.EQ.9.OR.INVTYP.EQ.10) THEN
        DO IXYZ=1,3
          FD(IXYZ,1)=-(FD(IXYZ,2)+FD(IXYZ,3))
        ENDDO
      ELSE IF (INVTYP.EQ.11) THEN
        DO IXYZ=1,3
          FD(IXYZ,2)=- FD(IXYZ,1)
        ENDDO
      ELSE IF (INVTYP.EQ.12) THEN
        DO IXYZ=1,3
          FD(IXYZ,2)=-(FD(IXYZ,1)+FD(IXYZ,4))
        ENDDO
      ELSE IF (INVTYP.EQ.13) THEN
        DO IXYZ=1,3
          FD(IXYZ,2)=-(FD(IXYZ,1)+FD(IXYZ,3))
        ENDDO
      ELSE IF (INVTYP.EQ.14) THEN
        DO IXYZ=1,3
          FD(IXYZ,3)=-(FD(IXYZ,1)+FD(IXYZ,2))
        ENDDO
      ELSE IF (INVTYP.EQ.15) THEN
        DO IXYZ=1,3
          FD(IXYZ,4)=-(FD(IXYZ,1)+FD(IXYZ,2)+FD(IXYZ,3))
        ENDDO
      ENDIF
      DO IXYZ=1,3
        GRADS(IXYZ,II)=GRADS(IXYZ,II)+FD(IXYZ,1)
        GRADS(IXYZ,JJ)=GRADS(IXYZ,JJ)+FD(IXYZ,2)
        GRADS(IXYZ,KK)=GRADS(IXYZ,KK)+FD(IXYZ,3)
        GRADS(IXYZ,LL)=GRADS(IXYZ,LL)+FD(IXYZ,4)
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! DSPDFSNOFl
      SUBROUTINE DSPDFSNOFl(IIFINT,IJKLG,DAB,SKIPI,SKIPJ,SKIPK,SKIPL,   &
                            MAXNUM,IIEQJJ,KKEQLL,IJEQKL,FD,II,JJ,KK,LL, &
                            SIZE_ENV,ENV,NAT,ATM,NBAS,BAS,IGTYP)
                             
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INTEGER :: I, J, K, L, IGTYP
      INTEGER :: DJEFF, DKEFF, DLEFF
      LOGICAL, INTENT(IN) :: IIEQJJ, KKEQLL, IJEQKL
      LOGICAL, INTENT(IN) :: SKIPI, SKIPJ, SKIPK, SKIPL
      INTEGER, INTENT(INOUT) :: IIFINT
      INTEGER, DIMENSION(MAXNUM), INTENT(IN) :: IJKLG
      DOUBLE PRECISION, DIMENSION(12), INTENT(INOUT) :: FD
      DOUBLE PRECISION, DIMENSION(MAXNUM), INTENT(IN) :: DAB

      INTEGER :: SIZE_ENV, NBAS
      DOUBLE PRECISION :: ENV(SIZE_ENV)
      INTEGER :: ATM(6, NAT), BAS(8, NBAS)

      INTEGER(4) :: SHLS(4)
      INTEGER :: DI, DJ, DK, DL
      DOUBLE PRECISION, ALLOCATABLE :: BLK1(:,:,:,:,:), BLK2(:,:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: BLK3(:,:,:,:,:), BLK4(:,:,:,:,:)

      INTEGER, EXTERNAL :: CINTcgto_cart, cint2e_ip1_cart
      INTEGER, EXTERNAL :: CINTcgto_spheric, cint2e_ip1_sph
!-----------------------------------------------------------------------

      IF(IGTYP==1) THEN
        DI = CINTcgto_cart(II-1, BAS)
        DJ = CINTcgto_cart(JJ-1, BAS)
        DK = CINTcgto_cart(KK-1, BAS)
        DL = CINTcgto_cart(LL-1, BAS)
      ELSE IF(IGTYP==2) THEN
        DI = CINTcgto_spheric(II-1, BAS)
        DJ = CINTcgto_spheric(JJ-1, BAS)
        DK = CINTcgto_spheric(KK-1, BAS)
        DL = CINTcgto_spheric(LL-1, BAS)
      END IF

      SHLS = [II-1, JJ-1, KK-1, LL-1]
      ALLOCATE(BLK1(DI,DJ,DK,DL,3))
      IF(IGTYP==1) err = cint2e_ip1_cart(BLK1, SHLS, ATM, NAT, BAS, NBAS, ENV, 0_8)
      IF(IGTYP==2) err = cint2e_ip1_sph(BLK1, SHLS, ATM, NAT, BAS, NBAS, ENV, 0_8)

      SHLS = [JJ-1, II-1, KK-1, LL-1]
      ALLOCATE(BLK2(DJ,DI,DK,DL,3))
      IF(IGTYP==1) err = cint2e_ip1_cart(BLK2, SHLS, ATM, NAT, BAS, NBAS, ENV, 0_8)
      IF(IGTYP==2) err = cint2e_ip1_sph(BLK2, SHLS, ATM, NAT, BAS, NBAS, ENV, 0_8)

      SHLS = [KK-1, LL-1, II-1, JJ-1]
      ALLOCATE(BLK3(DK,DL,DI,DJ,3))
      IF(IGTYP==1) err = cint2e_ip1_cart(BLK3, SHLS, ATM, NAT, BAS, NBAS, ENV, 0_8)
      IF(IGTYP==2) err = cint2e_ip1_sph(BLK3, SHLS, ATM, NAT, BAS, NBAS, ENV, 0_8)

      SHLS = [LL-1, KK-1, II-1, JJ-1]
      ALLOCATE(BLK4(DL,DK,DI,DJ,3))
      IF(IGTYP==1) err = cint2e_ip1_cart(BLK4, SHLS, ATM, NAT, BAS, NBAS, ENV, 0_8)
      IF(IGTYP==2) err = cint2e_ip1_sph(BLK4, SHLS, ATM, NAT, BAS, NBAS, ENV, 0_8)

      IJKLN=0
      DO I=1, DI
        DJEFF = MERGE(I, DJ, IIEQJJ)
        DO J=1, DJEFF
          DKEFF = MERGE(I, DK, IJEQKL)
          DO K=1, DKEFF
            DLEFF = MERGE(K, DL, KKEQLL)
            DLEFF = MERGE(J, DLEFF, IJEQKL.AND.K.EQ.I)
            DO L=1, DLEFF
              IJKLN = IJKLN+1
              NN = IJKLG(IJKLN)
!
              IF (.NOT. SKIPI) THEN
                IIFINT = IIFINT + 3
                FD(1) = FD(1) - BLK1(I,J,K,L,1)*DAB(NN)
                FD(2) = FD(2) - BLK1(I,J,K,L,2)*DAB(NN)
                FD(3) = FD(3) - BLK1(I,J,K,L,3)*DAB(NN)
              END IF
!
              IF (.NOT. SKIPJ) THEN
!
              IIFINT=IIFINT+3 
              FD( 4)=FD( 4)-BLK2(J,I,K,L,1)*DAB(NN)
              FD( 5)=FD( 5)-BLK2(J,I,K,L,2)*DAB(NN)
              FD( 6)=FD( 6)-BLK2(J,I,K,L,3)*DAB(NN)
              END IF
              IF (.NOT. SKIPK) THEN
!
              IIFINT=IIFINT+3
              FD( 7)=FD( 7)-BLK3(K,L,I,J,1)*DAB(NN)
              FD( 8)=FD( 8)-BLK3(K,L,I,J,2)*DAB(NN)
              FD( 9)=FD( 9)-BLK3(K,L,I,J,3)*DAB(NN)
              END IF
              IF (.NOT. SKIPL) THEN
!
              IIFINT=IIFINT+3
              FD(10)=FD(10)-BLK4(L,K,I,J,1)*DAB(NN)
              FD(11)=FD(11)-BLK4(L,K,I,J,2)*DAB(NN)
              FD(12)=FD(12)-BLK4(L,K,I,J,3)*DAB(NN)
              END IF
!
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      DEALLOCATE(BLK1)
      DEALLOCATE(BLK2)
      DEALLOCATE(BLK3)
      DEALLOCATE(BLK4)

      END
