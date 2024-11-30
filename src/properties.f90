!======================================================================!
!                                                                      !
!              P R O P E R T Y       S U B R O U T I N E S             !
!                                                                      !
!======================================================================!
!                                                                      !
!   NUCDIST: Internuclear distances                                    !
!                                                                      !
!   MULLIKENrc: Calculate Mulliken populations on Atoms and every MO   !
!   WRITEMULLAPMOrc: Print Mulliken atomic populations in each MO      !
!   ATDENMATrc: Compute atomic density matrix                          !
!                                                                      !
!   Extended Koopmans' Theorem (EKT): J. Chem. Phys. 136, 174116, 2012 !
!   EXTKOOPMANSrc: Ionization Potentials by Ext. Koopmans' theorem     !
!   DYSONORB: Dyson molecular orbitals (related to EKT)                !
!   DYSVECOUTrc: Print Dyson eigenvectors (related to EKT)             !
!                                                                      !
!   PUNCHVEC,PUNCHAPSG,OneNCO,OddNCO,EvenNCO: JCP 139, 234109, 2013    !
!   Create an input APSG File (9) containing generating wfn of PNOF5   !
!                                                                      !
!   CHEMPOTrc: Chemical Potential ( IJQC 116, 805, 2016 )              !
!                                                                      !
!   OUTPUTCJKrc: Print CJ12 & CK12 (Molecular RDMs) in 'CJK' File (12) !
!   OUTPUTTijab_rc: Print non-dynamic CK12 & Tijab (MP2 amplitudes)    !
!                   to 'CND' File (13)                                 !
!                                                                      !
!   Print Atomic RDMs in 1DM and 2DM Files (NOUTRDM=1,2,3)             !
!   NSQT=0: 2DM (14), 1DM (15) [formatted files]                       !
!   NSQT=1: 2DM (14), 1DM (15), N2DM (16) [unformatted f iles]         !
!   OUTPUTRDMrc: Calculate, check norms and print atomic DMs           !
!   SUMDDL: Perform the trasformation from NOs to atomic MOs for 2DM   !
!   SUMDL: Perform the trasformation from NOs to atomic MOs for 1DM    !
!   RDM1NORM: Calculate and print the 1DM norm                         !
!   RDM2NORM: Calculate and print the 2DM norm                         !
!   DENSI: Calculate density for each 2DM(ijkl)                        !
!                                                                      !
!======================================================================!

! NUCDIST
      SUBROUTINE NUCDIST(NV,NAT,Cxyz)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NV,NAT
      DOUBLE PRECISION,DIMENSION(NV),INTENT(IN)::Cxyz
      DOUBLE PRECISION,DIMENSION(NAT)::DIST
      DOUBLE PRECISION::RR
      INTEGER::I,J,MAXX,MINN
      DOUBLE PRECISION,PARAMETER::BOHR = 0.52917724924D+00               
!-----------------------------------------------------------------------
      WRITE(11,'(1X)')
      WRITE(11,*)'Internuclear distances (Angs)'
      MAXX=0
      DO
      MINN = MAXX+1
      MAXX = MAXX+5
      IF(MAXX>NAT) MAXX=NAT
      WRITE(11,'(1X)')
      WRITE(11,'(10X,7(I4,8X))') (J,J=MINN,MAXX)
      WRITE(11,'(1X)')
!      
      DO I = 1,NAT
       DO J = MINN,MAXX
        RR = (Cxyz(1+(I-1)*3)-Cxyz(1+(J-1)*3))**2 +                     &
             (Cxyz(2+(I-1)*3)-Cxyz(2+(J-1)*3))**2 +                     &
             (Cxyz(3+(I-1)*3)-Cxyz(3+(J-1)*3))**2
        RR = DSQRT(RR)
        DIST(J) = RR*BOHR
       ENDDO
       WRITE(11,'(I4,7(F12.4))')I,(DIST(J),J=MINN,MAXX)
      ENDDO
      IF(.NOT.(MAXX<NAT)) EXIT
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! MULLIKENrc
      SUBROUTINE MULLIKENrc(ATMNAME,ZNUC,LIMLOW,LIMSUP,S,RO,QD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
!
      CHARACTER*4 ATMNAME(NATOMS)
      DIMENSION ZNUC(NATOMS),LIMLOW(NATOMS),LIMSUP(NATOMS)
      DIMENSION S(NBF,NBF),RO(NBF5),QD(NBF,NBF,NBF)
      ALLOCATABLE:: CMUL(:,:),QMOAT(:,:)
      ALLOCATABLE:: AOVLPOP(:,:),OVLPOP(:,:),POPAT(:)
!-----------------------------------------------------------------------
!     Mulliken Population:
!     CMUL:    Atomic Orbital populations in each Molecular Orbital
!     QMOAT:   Atomic populations in each Molecular Orbital
!     AOVLPOP: Atomic Overlap Populations
!     OVLPOP:  Overlap Populations
!     POPAT:   Populations on Atoms
!-----------------------------------------------------------------------
!     Compute Mulliken population of each AO in every MO
!-----------------------------------------------------------------------
      ALLOCATE(CMUL(NBF,NBF5))
      DO IQ=1,NBF5
       DO J=1,NBF
        SCSC=0.0d0
        DO K=1,NBF
         SCSC=SCSC+QD(IQ,K,J)*S(J,K)
        ENDDO
        CMUL(J,IQ)=SCSC*RO(IQ)
       ENDDO
      ENDDO
!-----------------------------------------------------------------------
!     Condense to atoms
!-----------------------------------------------------------------------
      ALLOCATE (QMOAT(NATOMS,NBF5))
      DO IQ=1,NBF5
       DO IAT=1,NATOMS
        SCMUL = 0.0d0
        IMIN=LIMLOW(IAT)
        IMAX=LIMSUP(IAT)
        DO J=IMIN,IMAX
         SCMUL = SCMUL + CMUL(J,IQ)
        ENDDO
        QMOAT(IAT,IQ) = SCMUL
       ENDDO
      ENDDO      
!-----------------------------------------------------------------------
!     Print Atomic populations in each Molecular Orbital
!-----------------------------------------------------------------------
      WRITE(6,1)
      CALL WRITEMULLAPMOrc(QMOAT,RO,NBF5,NATOMS)
!-----------------------------------------------------------------------
!     Calculate Overlap Atomic Populations
!-----------------------------------------------------------------------
      ALLOCATE(AOVLPOP(NBF,NBF))
      AOVLPOP=0.0d0
      DO K=1,NBF
       DO L=1,NBF
        AOVLPOP(K,L) = AOVLPOP(K,L)                                     &
                     + ATDENMATrc(K,L,RO,QD,NBF,NBF5)*S(K,L)
       ENDDO
      ENDDO
!-----------------------------------------------------------------------
!     Compress Overlap Atomic Populations to Atoms
!-----------------------------------------------------------------------
      ALLOCATE(OVLPOP(NATOMS,NATOMS))
      DO I=1,NATOMS
       IMIN=LIMLOW(I)
       IMAX=LIMSUP(I)
       DO J=1,NATOMS
        JMIN=LIMLOW(J)
        JMAX=LIMSUP(J)
        OVLPOP(I,J)=0.0d0
        IF(IMIN<=IMAX)THEN
         DO K=IMIN,IMAX
          IF(JMIN<=JMAX)THEN
           DO L=JMIN,JMAX
            OVLPOP(I,J) = OVLPOP(I,J) + AOVLPOP(K,L)
           ENDDO
          ENDIF
         ENDDO
        ENDIF
       ENDDO
      ENDDO
!-----------------------------------------------------------------------
!     Calculate total populations on Atoms
!-----------------------------------------------------------------------
      WRITE(6,2)
      WRITE(6,3)
      ALLOCATE(POPAT(NATOMS))
      DO I=1,NATOMS
       POPAT(I)=0.0d0
       DO J=1,NATOMS
        POPAT(I)=POPAT(I)+OVLPOP(I,J)
       ENDDO
!-----------------------------------------------------------------------
!      Print total Populations
!-----------------------------------------------------------------------
       WRITE(6,4)I,ATMNAME(I),POPAT(I),ZNUC(I)-POPAT(I)
      ENDDO
!-----------------------------------------------------------------------
    1 FORMAT(                                                           &
      /2X,' Atomic populations in each Molecular Orbital ',/,           &
       2X,'- - - - - - - - - - - - - - - - - - - - - - - -')             
    2 FORMAT(/,                                                         &
      /2X,' Total Populations on Atoms ',/,                             &
       2X,'- - - - - - - - - - - - - - -',/)
    3 FORMAT(9X,'Atom',8X,'Pop.',7X,'Charge',/)
    4 FORMAT(1X,I4,1X,A8,2F12.4)
!-----------------------------------------------------------------------
      DEALLOCATE(CMUL,QMOAT,AOVLPOP,OVLPOP)
      RETURN
      END

! WRITEMULLAPMOrc
      SUBROUTINE WRITEMULLAPMOrc(QMOAT,RO,NBF5,NATOMS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DIMENSION QMOAT(NATOMS,NBF5),RO(NBF5)
!-----------------------------------------------------------------------
!     Print Atomic populations in each Molecular Orbital
!-----------------------------------------------------------------------
      IMAX=0
  10  IMIN=IMAX+1
      IMAX=IMAX+10
      IF(IMAX>NBF5)IMAX=NBF5
      WRITE(6,1)
      WRITE(6,2)(I,I=IMIN,IMAX)
      WRITE(6,1)
      WRITE(6,3)(2*RO(I),I=IMIN,IMAX)
      WRITE(6,1)
      DO J=1,NATOMS
       WRITE(6,4)J,(QMOAT(J,I),I=IMIN,IMAX)
      ENDDO
      IF(IMAX<NBF5)GOTO 10
!-----------------------------------------------------------------------
    1 FORMAT(1X)
    2 FORMAT(9X,10(3X,I4,2X))
    3 FORMAT(9X,10F9.4)
    4 FORMAT(I5,4X,10F9.4)
!-----------------------------------------------------------------------
      RETURN
      END

! ATDENMATrc
      FUNCTION ATDENMATrc(K,L,RO,QD,NBF,NBF5)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DIMENSION RO(NBF5),QD(NBF,NBF,NBF)
!-----------------------------------------------------------------------
      ATDENMATrc = 0.0d0
      DO J =1,NBF5
       ATDENMATrc = ATDENMATrc + RO(J)*QD(J,K,L)
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! EXTKOOPMANSrc
      SUBROUTINE EXTKOOPMANSrc(ELAG,COEF,OVERLAP,AHCORE,IJKL,XIJKL,RO)
!-----------------------------------------------------------------------
!                Extended Koopmans' Theorem (EKT)
!                This subroutine is called when IEKT=1
!-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL HighSpin
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
      COMMON/INPFILE_NIJKL/NINTMX,NIJKL,NINTCR,NSTORE
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      COMMON/EHFEN/EHF,EN
      COMMON/ENERGIAS/EELEC,EELEC_OLD,DIF_EELEC,EELEC_MIN
!
      INTEGER,DIMENSION(NIJKL)::IJKL
      DOUBLE PRECISION,DIMENSION(NBF5)::RO
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::ELAG,COEF,OVERLAP,AHCORE
      DOUBLE PRECISION,DIMENSION(NIJKL)::XIJKL
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::EVA,TEMP,OCCD
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::AUX,W,DYSON
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::DM,GFOCK,COEFT
!-----------------------------------------------------------------------
!     Intermediate matrices
!-----------------------------------------------------------------------
      ALLOCATE (AUX(NBF,NBF),W(NBF,NBF),EVA(NBF),TEMP(NBF))
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - - - - -
!     Ionization Potentials:
!     Diagonalization of -ELAG/RAIZ[RO*RO], ELAG: Lagrangian
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - - - - -
!     ELAG/RAIZ[RO*RO] -> square symmetric matrix (AUX)
!-----------------------------------------------------------------------
      DO I=1,NBF5
       BETAi=DSQRT(RO(i))
       DO J=1,I
        BETAj=DSQRT(RO(j))
        IF( (DABS(BETAi)>1.0d-6).and.(DABS(BETAj)>1.0d-6) )THEN
         AUX(I,J)=ELAG(I,J)/(BETAi*BETAj)
        ELSE
         AUX(I,J)=0.0d0
        ENDIF
        AUX(J,I)=AUX(I,J)
       ENDDO
      ENDDO
      DO I=NBF5+1,NBF
       DO J=1,I
        AUX(I,J)=0.0d0
        AUX(J,I)=AUX(I,J)
       ENDDO
      ENDDO
!-----------------------------------------------------------------------
!     Diagonalize square matrix (AUX)
!     W - Eigenvectors, EVA - EIGENVALUES
!-----------------------------------------------------------------------
      CALL DIAG(NBF,AUX,W,EVA,TEMP)
!-----------------------------------------------------------------------
!     Write Ionization Energies (Output File)
!     Coefficients of Dyson Orbitals (DYSON=COEF*BETA*W)
!-----------------------------------------------------------------------
      ALLOCATE (DYSON(NBF,NBF5),OCCD(NBF5))
      CALL DYSONORB(OCCD,DYSON,COEF,RO,W,OVERLAP)
      WRITE(6,1)
      DO I=1,NBF5
       IF(EVA(I)<0.0d0.and.OCCD(I)>0.45d0)THEN
        WRITE(6,2)I,-EVA(I),-EVA(I)*27.21138386
       ENDIF
      ENDDO
!      
      IWRITEDYSON=0
      IF(IWRITEDYSON==1)THEN
       WRITE(6,4)
       CALL DYSVECOUTrc(DYSON,OCCD,NBF,NBF5,NBF5)
      ENDIF
      DEALLOCATE (DYSON)
!-----------------------------------------------------------------------
!     Write Cation Energy if ICATION=1 (Output File)
!     Ecation = Eelec + EN + MinVal IonPotential obtained with EKT
!-----------------------------------------------------------------------
      ICATION=0                                          ! ICATION=0
      IF(NSOC==0.and.ICATION==1)THEN
       TEMP = 1.0d06
       DO I=1,NBF5
        IF(-EVA(I)>0.0d0.and.OCCD(I)>0.45d0)TEMP(I)=-EVA(I)
       ENDDO
       PIMIN = MINVAL(TEMP,NBF5)
       WRITE(6,3)EELEC+EN+PIMIN
      ENDIF
      DEALLOCATE (OCCD)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - - - - -
!     Electron Affinities:
!     Diagonalization of (ELAG-GFOCK)/RAIZ[(1-RO)*(1-RO)]
!     GFOCK: Generalized Fock Matrix, ELAG: Lagrangian
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IEA=0                                                 ! IEA=0
      IF(IEA==1)THEN
!-----------------------------------------------------------------------
!      Form the Generalized Fock Matrix
!-----------------------------------------------------------------------
       ALLOCATE (DM(NBF,NBF),GFOCK(NBF,NBF))
       CALL DENMATr(DM,COEF,RO,NBF,1,NBF5)
       CALL FORM2JK(AUX,DM,IJKL,XIJKL)
       DEALLOCATE(DM)
       GFOCK = AHCORE + AUX
       AUX  = MATMUL(GFOCK,COEF)
       ALLOCATE (COEFT(NBF,NBF))
       COEFT = TRANSPOSE(COEF)     
       GFOCK = MATMUL(COEFT,AUX)
       DEALLOCATE (COEFT)
!-----------------------------------------------------------------------
!      (ELAG-GFOCK)/RAIZ[(1-RO)*(1-RO)] -> square symmetric matrix (AUX)
!-----------------------------------------------------------------------
       DO I=1,NB
        AUX(I,I) = -DABS(ELAG(I,I)-GFOCK(I,I)) / (1.0d0-RO(i))
       ENDDO

       DO I=NB+1,NBF5
        AUX(I,I) = (ELAG(I,I)-GFOCK(I,I)) / (1.0d0-RO(i))
       ENDDO

       DO I=1,NBF5
        BETAi=DSQRT(1.0d0-RO(i))
        DO J=1,I-1
         BETAj=DSQRT(1.0d0-RO(j))
         IF( (DABS(BETAi)>1.0d-3).and.(DABS(BETAj)>1.0d-3) )THEN
          AUX(I,J) = (ELAG(I,J)-GFOCK(I,J)) / (BETAi*BETAj)
         ELSE
          AUX(I,J)=0.0d0
         ENDIF
         AUX(J,I)=AUX(I,J)
        ENDDO
       ENDDO

       DO I=NBF5+1,NBF
        DO J=1,I
         AUX(I,J) = -GFOCK(I,J)
         AUX(J,I) = AUX(I,J)
        ENDDO
       ENDDO

       DEALLOCATE (GFOCK)
!-----------------------------------------------------------------------      
!      Diagonalize square matrix (AUX)
!      W - Eigenvectors, EVA - Eigenvalues
!-----------------------------------------------------------------------
       CALL DIAG(NBF,AUX,W,EVA,TEMP)
!-----------------------------------------------------------------------
!      Write electron affinities (Output File)
!-----------------------------------------------------------------------
       WRITE(6,5)
       DO I=1,NBF
        WRITE(6,2)I,EVA(I),EVA(I)*27.2107
       ENDDO
      ENDIF
!-----------------------------------------------------------------------
    1 FORMAT(/2X,'--------------------------------------------------',  &
             /2X,' Extended Koopmans'' Theorem (Ionization Energies) ', &
             /2X,'--------------------------------------------------',  &
            //4X,'OM',14X,'(aU)',14X,'(eV)',/)                           
    2 FORMAT(2X,I4,4X,F15.3,4X,F15.3,10X,F7.3)                           
    3 FORMAT(/,3X,'NOF Total Cation Energy (aU) =',F16.6)                
    4 FORMAT(/,                                                         &
       18X,'----------------------------------------',/,                &
       18X,' DYSON ORBITALS IN ATOMIC ORBITAL BASIS ',/,                &
       18X,'----------------------------------------')                   
    5 FORMAT(/2X,'-------------------------------------------------',   &
             /2X,' EXTENDED KOOPMANS THEOREM (ELECTRON AFFINITIES) ',   &
             /2X,'-------------------------------------------------',   &
            //20X,'(aU)',14X,'(eV)',/)
!-----------------------------------------------------------------------
      DEALLOCATE (AUX,W,EVA,TEMP)
      RETURN
      END

! QUASIENEr
      SUBROUTINE QUASIENEr(RO,HCORE,QJ,QK,CJ12,CK12)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
!
      DOUBLE PRECISION,DIMENSION(NBF5)::RO,HCORE
      DOUBLE PRECISION,DIMENSION(NBFT5)::QJ,QK
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5)::CJ12,CK12 
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::EHFa,EQP
!-----------------------------------------------------------------------
!     EAHF: New HF energies
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE (EHFa(NBF5))
      EHFa = 0.0d0
      DO j=na+1,nbf5 
       do i=1,j-1
        ij = i + j*(j-1)/2
        EHFa(j) = EHFa(j) + RO(i)*(2.0d0*QJ(ij)-QK(ij))
       enddo
       jj=j*(j+1)/2
       EHFa(j) = EHFa(j) + RO(j)*QJ(jj)
       do i=j+1,nbf5
        ij = j + i*(i-1)/2
        EHFa(j) = EHFa(j) + RO(i)*(2.0d0*QJ(ij)-QK(ij))
       enddo
       EHFa(j) = HCORE(j) + EHFa(j)/(1.0-RO(j))
      ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     EQP: Quasi-particle energies
!     MSpin=0 (Ms=0): Singlet States (S=0) & Multiplet States (S>0) 
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE (EQP(NBF5))
      DO j=1,nb
       jj=j*(j+1)/2
       EQP(j) = HCORE(j) + QJ(jj) + (1.0d0/RO(j)) *                     &
                ( PRODCWQWj(j,CJ12,QJ) - PRODCWQWj(j,CK12,QK) )
      ENDDO        
      DO j=nb+1,na  ! NSOC                                                    
       EQP(j) = HCORE(j) + (1.0d0/RO(j))                                &
              * ( PRODCWQWj(j,CJ12,QJ) - PRODCWQWj(j,CK12,QK) )           
      ENDDO                                                             
      DO j=na+1,nbf5                                                    
       jj=j*(j+1)/2                                                     
       EQP(j) = EHFa(j) - RO(j)*QJ(jj)/(1.0d0-RO(j))                    &
           - (PRODCWQWj(j,CJ12,QJ)-PRODCWQWj(j,CK12,QK))/(1.0d0-RO(j))
      ENDDO 
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL ORDEREQP(RO,EQP)
      WRITE(6,1)
      DO I=1,NBF5
       WRITE(6,2)I,EQP(I),EQP(I)*27.2107
      ENDDO
!-----------------------------------------------------------------------
    1 FORMAT(/2X,'--------------------------------------------------',  &
             /2X,' Quasi Particle Energies                           ', &
             /2X,'--------------------------------------------------',  &
            //4X,'OM',14X,'(aU)',14X,'(eV)',/) 
    2 FORMAT(2X,I4,4X,F15.3,4X,F15.3,10X,F7.3)            
!-----------------------------------------------------------------------
      DEALLOCATE (EHFa,EQP)
      RETURN
      END
      
! ORDEREQP
      SUBROUTINE ORDEREQP(RO,EQP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_ORBSPACE0/NO1,NDOC,NCO,NCWO,NVIR,NAC,NO0      
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      DOUBLE PRECISION,DIMENSION(NBF5)::RO
      DOUBLE PRECISION,DIMENSION(NBF)::EQP
!----------------------------------------------------------------------!
!     Ordering by Occupation Numbers (NO1+1:NO1+NDOC)=(NO1+1:NB)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      do k=1,ndoc
       kmini = k
       im = NO1+k  ! im=no1+1,nb
       MINI = im
       DUM = RO(im)
       do i=im,nb
        IF(RO(i)>DUM)THEN
         DUM = RO(i)
         MINI = i  ! MINI -> MAXI
         kmini = MINI - NO1
        ENDIF
       end do
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
       IF(MINI/=im)THEN
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       Occupancies
        DUM1 = RO(im)
        RO(im) = RO(MINI)
        RO(MINI) = DUM1
!       Energies
        DUM2 = EQP(im)
        EQP(im) = EQP(MINI)
        EQP(MINI) = DUM2
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
         DUM2 = EQP(in)
         EQP(in) = EQP(inmini)
         EQP(inmini) = DUM2
        end do
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       END IF
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
      end do
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END      
      
! DYSONORB
      SUBROUTINE DYSONORB(OCCD,DYSON,COEF,RO,W,OVERLAP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      DOUBLE PRECISION,DIMENSION(NBF5)::RO,OCCD
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::COEF,W,OVERLAP
      DOUBLE PRECISION,DIMENSION(NBF,NBF5)::DYSON
!-----------------------------------------------------------------------
!     DYSON=COEF*BETA*WC (unnormalized)
!-----------------------------------------------------------------------
      DO NIU=1,NBF
       DO IQ=1,NBF5
        DYSON(NIU,IQ)=0.0d0
        do i=1,nbf5
         BETAi = DSQRT(RO(i))
         DYSON(NIU,IQ)=DYSON(NIU,IQ)+COEF(NIU,i)*BETAi*W(i,IQ)
        enddo
       ENDDO
      ENDDO
!-----------------------------------------------------------------------
!     Occupation Number of Dyson Orbital (Normalization Factor)
!-----------------------------------------------------------------------
      DO IQ=1,NBF5
       OCCD(IQ)=0.0d0
       do niu=1,nbf
        OCCD(IQ) = OCCD(IQ) + DYSON(niu,IQ)*FC(niu,IQ,OVERLAP,DYSON)
       enddo
      ENDDO
!-----------------------------------------------------------------------
!     Normalized Dyson Orbitals
!-----------------------------------------------------------------------
      DO IQ=1,NBF5
       DO niu=1,nbf
        DYSON(niu,IQ)=DYSON(niu,IQ)/DSQRT(OCCD(IQ))
       ENDDO
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! DYSVECOUTrc
      SUBROUTINE DYSVECOUTrc(V,RO,NBF,NBF5,NL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
!-----------------------------------------------------------------------
!     Print Eigenvectors
!-----------------------------------------------------------------------
      CHARACTER*2 LABELAT
      CHARACTER*4 BFNAM1
      CHARACTER*6 BFNAM2
      DOUBLE PRECISION,DIMENSION(NBF5)::RO
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::V
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::DOSRO
      ALLOCATE(DOSRO(NBF))
!-----------------------------------------------------------------------
      DO J=1,NBF5
       DOSRO(J)=2.0*RO(J)
      ENDDO
      DO J=NBF5+1,NBF
       DOSRO(J)=0.0
      ENDDO
      MAX=0
  10  MIN=MAX+1
      MAX=MAX+5
      IF(MAX>NL) MAX=NL
      WRITE(6,1)(J,J=MIN,MAX)
      WRITE(6,2)(DOSRO(J),J=MIN,MAX)
      WRITE(6,*)
      REWIND(4)
      READ(4,'(I5)')NBF0
      DO I=1,NBF
       IF(I<=35)THEN
        READ(4,3)LABELAT,IAT,BFNAM1
        WRITE(6,4) I,LABELAT,IAT,BFNAM1,(V(I,J),J=MIN,MAX)
       ELSE
        READ(4,5)LABELAT,BFNAM2
        WRITE(6,6) I,LABELAT,BFNAM2,(V(I,J),J=MIN,MAX)
       ENDIF
      ENDDO
      IF(MAX<NL) GOTO 10
      WRITE(6,7)
!-----------------------------------------------------------------------
    1 FORMAT(/,15X,5(4X,I4,3X))
    2 FORMAT(/,15X,5F11.4)
    3 FORMAT(A2,I2,A4)
    4 FORMAT(I5,2X,A2,I2,A4,5F11.6)
    5 FORMAT(A2,A6)
    6 FORMAT(I5,2X,A2,A6,5F11.6)
    7 FORMAT(/)
!-----------------------------------------------------------------------
      DEALLOCATE(DOSRO)
      RETURN
      END

! PUNCHVEC
      SUBROUTINE PUNCHVEC(V,NBF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::V
!-----------------------------------------------------------------------
      DO J = 1,NBF
       K = 0
       MAX = 0
    1  MIN = MAX + 1
       MAX = MAX + 5
       K = K + 1
       IF(MAX>NBF)MAX = NBF
       MODJ = MOD(J,100)
       WRITE(9,10) MODJ,K,(V(I,J),I = MIN,MAX)
       IF (MAX<NBF)GOTO 1
      ENDDO
!-----------------------------------------------------------------------
   10 FORMAT(I2,I3,1P,5E15.8)
      RETURN
      END

! PUNCHAPSG
      SUBROUTINE PUNCHAPSG(NO1,NCWO,NCO,NBF5,RO,SUMA,THAPSG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DOUBLE PRECISION,DIMENSION(NBF5)::RO
      ALLOCATABLE::IOCU(:),BB(:,:)
!-----------------------------------------------------------------------
!     NCO:  Number of HF occupied MOs (OCC=1)
!     NVIR: Number of HF virtual  MOs (OCC=0)
!
!     NO1:  Number of inactive doubly occupied orbitals (OCC=1)
!     NDOC: Number of strongly occupied MOs
!     NCWO: Number of coupled weakly occupied MOs per strongly occupied
!     NCWO*NDOC: Active orbitals in the virtual subspace
!     NO0: Empty orbitals  (OCC=0)
!     NAC:  Dimension of the active natural orbital subspace
!
!           NCO     |       NVIR          = NBF
!       NO1 + NDOC  |  NCWO*NDOC + NO0    = NBF
!           |      NAC           |
!-----------------------------------------------------------------------
      NO1NAC=NCO+NCWO*(NCO-NO1)                             ! NO1+NAC
      ALLOCATE(IOCU(NO1NAC),BB(NCO,NO1NAC))
      IOCU = 0
      SUMA = 0.0d0

!     HF Fermi Vacuum
      D0 = 1.0d0
      do i=1,NCO
       IOCU(i) = i
       D0 = D0*DSQRT(RO(i))
      enddo
      WRITE(9,1)D0,(IOCU(i),-IOCU(i),i=1,NCO)
      SUMA = SUMA + D0*D0

      do i=NO1+1,NCO
       !do j=NCO+NCWO*(NCO-i)+1,NCO+NCWO*(NCO-i+1)  !old-sort
       do j=NCO+NCO-i+1,NCO+NCO-i+1+(NCO-1)*NCO,NCO !new-sort
        BB(i,j) = DSQRT(RO(j)/RO(i))
       enddo
      enddo

      CALL OneNCO(NO1,NCWO,NCO,NO1NAC,IOCU,BB,SUMA,D0,THAPSG)
      IF(NCO>1)THEN
       do ig=2,NCO
        if(mod(ig, 2)==0)then
         CALL EvenNCO(NO1,NCO-ig+2,NCWO,NCO,NO1NAC,                     &
                      IOCU,BB,SUMA,D0,THAPSG)                            
        else                                                             
         CALL  OddNCO(NO1,NCO-ig+3,NCWO,NCO,NO1NAC,                     &
                      IOCU,BB,SUMA,D0,THAPSG)
        endif
       enddo
      ENDIF
!-----------------------------------------------------------------------
    1 FORMAT(F20.16,2X,99I4)
      DEALLOCATE(IOCU,BB)
      RETURN
      END

! OneNCO
      SUBROUTINE OneNCO(I1,NCWO,NCO,NO1NAC,IOCU,BB,SUMA,DD,THAPSG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)         
      INTEGER,DIMENSION(NO1NAC)::IOCU
      DOUBLE PRECISION,DIMENSION(NCO,NO1NAC)::BB
!-----------------------------------------------------------------------
      do i=I1+1,NCO
       !do ip=NCO+NCWO*(NCO-i)+1,NCO+NCWO*(NCO-i+1)  !old-sort
       do ip=NCO+NCO-i+1,NCO+NCO-i+1+(NCO-1)*NCO,NCO !new-sort
        IOCU(i)=ip
        Di = -DD*BB(i,ip)
        IF(DABS(Di)>THAPSG)THEN
         WRITE(9,1)Di,(IOCU(ii),-IOCU(ii),ii=1,NCO)
         SUMA = SUMA + Di*Di
        ENDIF
        IOCU(i) = i
       enddo
      enddo
!-----------------------------------------------------------------------
    1 FORMAT(F20.16,2X,99I4)
      RETURN
      END

! OddNCO
      SUBROUTINE OddNCO(I1,I2,NCWO,NCO,NO1NAC,IOCU,BB,SUMA,DD,THAPSG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      INTEGER,DIMENSION(NO1NAC)::IOCU
      DOUBLE PRECISION,DIMENSION(NCO,NO1NAC)::BB
!-----------------------------------------------------------------------
      do i=I1+1,I2-2
       !do ip=NCO+NCWO*(NCO-i)+1,NCO+NCWO*(NCO-i+1)  !old-sort
       do ip=NCO+NCO-i+1,NCO+NCO-i+1+(NCO-1)*NCO,NCO !new-sort
        IOCU(i) = ip
        Di = -DD*BB(i,ip)
        CALL EvenNCO(i,I2,NCWO,NCO,NO1NAC,IOCU,BB,SUMA,Di,THAPSG)
        IOCU(i) = i
       enddo
      enddo
!-----------------------------------------------------------------------
      RETURN
      END

! EvenNCO
      SUBROUTINE EvenNCO(I1,I2,NCWO,NCO,NO1NAC,IOCU,BB,SUMA,DD,THAPSG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      INTEGER,DIMENSION(NO1NAC)::IOCU
      DOUBLE PRECISION,DIMENSION(NCO,NO1NAC)::BB
!-----------------------------------------------------------------------
      do i=I1+1,I2-1
       !do ip=NCO+NCWO*(NCO-i)+1,NCO+NCWO*(NCO-i+1)  !old-sort
       do ip=NCO+NCO-i+1,NCO+NCO-i+1+(NCO-1)*NCO,NCO !new-sort
        IOCU(i) = ip
        Di = -DD*BB(i,ip)
        if(I2==NCO)then
         CALL OneNCO(i,NCWO,NCO,NO1NAC,IOCU,BB,SUMA,Di,THAPSG)
        else
         CALL OddNCO(i,I2+2,NCWO,NCO,NO1NAC,IOCU,BB,SUMA,Di,THAPSG)
        endif
        IOCU(i) = i
       enddo
      enddo
!-----------------------------------------------------------------------
      RETURN
      END

! CHEMPOTrc
      SUBROUTINE CHEMPOTrc(HCORE,QJ,QK,RO,DIPx,DIPy,DIPz)
!-----------------------------------------------------------------------
!     This subroutine is called when ICHEMPOT = 1 for PNOF5 and PNOF7
!-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      LOGICAL EFIELDL,HighSpin
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_ORBSPACE0/NO1,NDOC,NCO,NCWO,NVIR,NAC,NO0
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
      COMMON/INPNOF_PNOF/IPNOF,NTWOPAR
      COMMON/INPNOF_STATIC/Ista
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      COMMON/INP_EFIELDL/EX,EY,EZ,EFIELDL
!
      DOUBLE PRECISION,DIMENSION(NBF5)::HCORE,RO,DIPx,DIPy,DIPz
      DOUBLE PRECISION,DIMENSION(NBFT5)::QJ,QK
      ALLOCATABLE::CMIU(:),DCJ12DRO(:,:),DCK12DRO(:,:)
!-----------------------------------------------------------------------
      ALLOCATE (DCJ12DRO(NBF5,NBF5),DCK12DRO(NBF5,NBF5))
      DCJ12DRO = 0.0d0
      DCK12DRO = 0.0d0
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!                      Inter-pair interactions
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
      IF(IPNOF==5)THEN                 ! PNOF5
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       DO j=1,NBF5
        DO i=1,NBF5
         DCJ12DRO(j,i) = 2.0d0*RO(i)
         DCK12DRO(j,i) = RO(i)
        ENDDO
       ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -       
      ELSE IF(IPNOF==7)THEN            ! PNOF7
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -      
       if(Ista==0)then
        DO j=1,NBF5
         DO i=1,NBF5
          DCJ12DRO(j,i) = 2.0d0*RO(i)
          FIj = DSQRT(RO(j)*(1.0d0-RO(j)))
          IF(FIj/=0.0d0)THEN          
           DCK12DRO(j,i) = RO(i)                                        &
                         + (0.5d0-RO(j))*DSQRT(RO(i)*(1.0d0-RO(i)))/FIj  
          ENDIF                                                          
         ENDDO                                                           
        ENDDO                                                            
       else if(Ista==1)then                                              
        DO j=1,NBF5                                                      
         DO i=1,NBF5                                                     
          DCJ12DRO(j,i) = 2.0d0*RO(i)                                    
          DCK12DRO(j,i) = RO(i)                                         &
                        + 4.0d0*(1.0d0-2.0d0*RO(j))*RO(i)*(1.0d0-RO(i))
         ENDDO
        ENDDO
       end if
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -      
      END IF       
!      
      if(MSpin==0.and.NSOC>1)then                     
       DO j=NB+1,NA
        DO i=NB+1,NA
         DCK12DRO(j,i) = 2.0d0*RO(i)        
        ENDDO      
       ENDDO
      endif
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!                      Intra-pair interactions
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!     below-above Fermi level interaction
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO j=1,NDOC
       jn = NO1+j
       !DO i=NDNS+NCWO*(NDOC-j)+1,NDNS+NCWO*(NDOC-j+1) !old-sort
       DO i=NCO+NCO-j+1,NCO+NCO-j+1+(NCO-1)*NCO,NCO   !new-sort
        in = NO1+i
        DCJ12DRO(jn,in) = 0.0d0
        DCJ12DRO(in,jn) = 0.0d0
        IF(RO(jn)/=0.0d0)THEN
         DCK12DRO(jn,in) = DSQRT(RO(in))/(2.0*DSQRT(RO(jn)))
        ELSE
         DCK12DRO(jn,in) = 0.0d0
        ENDIF
        IF(RO(in)/=0.0d0)THEN
         DCK12DRO(in,jn) = DSQRT(RO(jn))/(2.0*DSQRT(RO(in)))
        ELSE
         DCK12DRO(in,jn) = 0.0d0
        ENDIF
       ENDDO
      ENDDO
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     above-above Fermi level interaction for each geminal 'l'
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO l=1,NDOC
       !DO j=NDNS+NCWO*(NDOC-l)+1,NDNS+NCWO*(NDOC-l+1)!old-sort
       DO j=NCO+NCO-l+1,NCO+NCO-l+1+(NCO-1)*NCO,NCO   !new-sort
        jn = NO1+j
        !DO i=NDNS+NCWO*(NDOC-l)+1,NDNS+NCWO*(NDOC-l+1)!old-sort
        DO i=NCO+NCO-l+1,NCO+NCO-l+1+(NCO-1)*NCO,NCO   !new-sort
         in = NO1+i
         DCJ12DRO(jn,in) = 0.0d0
         IF(RO(jn)/=0.0d0)THEN
          DCK12DRO(jn,in) = - DSQRT(RO(in))/(2.0*DSQRT(RO(jn)))
         ELSE
          DCK12DRO(jn,in) = 0.0d0
         ENDIF
        ENDDO
       ENDDO
      ENDDO
!-----------------------------------------------------------------------      
!     Chemical Potential for each subspace      
!-----------------------------------------------------------------------
      ALLOCATE (CMIU(NBF5))
!
      DO j=1,NDOC
       jn = NO1+j      
       jj = jn*(jn+1)/2
       CMIU(jn) = HCORE(jn) + QJ(jj)/2.0 + PRODCWQWj(jn,DCJ12DRO,QJ)    &
                                         - PRODCWQWj(jn,DCK12DRO,QK)     
      ENDDO                                                              
      IF(MSpin==0.and.NSOC>1)THEN                                        
       DO j=NDOC+1,NDNS                                                  
        jn = NO1+j                                                       
        jj = jn*(jn+1)/2                                                 
        CMIU(jn) = HCORE(jn) + PRODCWQWj(jn,DCJ12DRO,QJ)                &
                             - PRODCWQWj(jn,DCK12DRO,QK)
       ENDDO
      ENDIF       
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Including Electric Field
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(EFIELDL)THEN
       DO j=1,NDNS
        jn = NO1+j      
        CMIU(jn) = CMIU(jn) - (EX*DIPx(jn)+EY*DIPy(jn)+EZ*DIPz(jn))
       ENDDO
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      write(6,'(/,11X,A19,/)')'Chemical Potentials'
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      DO j=1,NDNS
!       jn = NO1+j      
!       write(6,'(I6,2F19.6)')jn,CMIU(jn),CMIU(jn)*27.21138386
!      ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CHEPOT = -1.0d10
      DO j=1,NDNS
       jn = NO1+j
       if(CMIU(jn)>CHEPOT)CHEPOT=CMIU(jn)
      ENDDO
      write(6,1)CHEPOT,CHEPOT*27.21138386             
!-----------------------------------------------------------------------
    1 FORMAT(/,3X,'Chemical Potential =',F10.4,2X,'(',F10.4,1X,'eV )')
      DEALLOCATE (CMIU,DCJ12DRO,DCK12DRO)
      RETURN
      END

! OUTPUTCJKrc
      SUBROUTINE OUTPUTCJKrc(RO,CJ12,CK12)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      COMMON/INPNOF_ORBSPACE0/NO1,NDOC,NCO,NCWO,NVIR,NAC,NO0
      COMMON/INPNOF_CJK/NOUTCJK,NTHRESHCJK,THRESHCJK      
      COMMON/INPNOF_Tijab/NOUTTijab,NTHRESHTijab,THRESHTijab      
!
      DOUBLE PRECISION,DIMENSION(NBF5)::RO
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5)::CJ12,CK12
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::BETA,FIs
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::CK12nd
!----------------------------------------------------------------------- 
!     Print CJ12 and -CK12
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO i=1,NBF5
       if(RO(i)>THRESHCJK)WRITE(12)i,i,i,i,RO(i)    ! CLii = CJii
       DO j=1,i-1
        if(DABS(CJ12(i,j))>THRESHCJK)then
         WRITE(12)i,j,i,j,CJ12(i,j)
         WRITE(12)j,i,j,i,CJ12(i,j)
        end if
        if(DABS(CK12(i,j))>THRESHCJK)then            ! Print -CK12
         WRITE(12)i,j,j,i,-CK12(i,j)
         WRITE(12)j,i,i,j,-CK12(i,j)
        end if
       ENDDO
      ENDDO
      WRITE(12)0,0,0,0,0.0d0
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!                   CK12nd: Non-Dynamic (Static) CK12
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
      IF(NOUTTijab==1)THEN
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -      
       ALLOCATE(FIs(NBF5),BETA(NBF5),CK12nd(NBF5,NBF5))
       FIs = 0.0d0
       BETA = 0.0d0
       CK12nd = 0.0d0
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      FIs(j)  = 2*RO(j)*HOLE(j)
!      BETA(j) = DSQRT(2*HOLE(j)) * DSQRT(RO(j))
!      BETA(j) = DSQRT(2*RO(j))   * DSQRT(RO(j))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       DO j=NO1+1,NBF5
        FIs(j) = 2.0d0*RO(j)*(1.0d0-RO(j))
        Cj = 1.0d0 - DABS(1.0d0-2.0d0*RO(j))
        BETA(j) = DSQRT( Cj*RO(j) )
       ENDDO
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!      Inter-pair Non-Dynamic (Static) Electron Correlation
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       DO j=NO1+1,NBF5
        DO i=NO1+1,NBF5
         CK12nd(j,i) = FIs(j)*FIs(i)
        ENDDO
       ENDDO
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!      Intra-pair Non-Dynamic (Static) Electron Correlation
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!      below-above Fermi level interaction
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       DO j=1,NDOC
        jn = NO1+j
        !DO i=NDOC+NCWO*(NDOC-j)+1,NDOC+NCWO*(NDOC-j+1)!old-sort
        DO i=NCO+NCO-j+1,NCO+NCO-j+1+(NCO-1)*NCO,NCO   !new-sort
         in = NO1+i
         CK12nd(jn,in) = BETA(jn)*BETA(in)
         CK12nd(in,jn) = BETA(in)*BETA(jn)
        ENDDO
       ENDDO
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      above-above Fermi level interaction for each geminal 'l'
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       DO l=1,NDOC
        !DO j=NDOC+NCWO*(NDOC-l)+1,NDOC+NCWO*(NDOC-l+1)!old-sort
        DO j=NCO+NCO-l+1,NCO+NCO-l+1+(NCO-1)*NCO,NCO   !new-sort
         jn = NO1+j
         !DO i=NDOC+NCWO*(NDOC-l)+1,NDOC+NCWO*(NDOC-l+1)!old-sort
         DO i=NCO+NCO-l+1,NCO+NCO-l+1+(NCO-1)*NCO,NCO   !new-sort
          in = NO1+i
          CK12nd(jn,in) = - BETA(jn)*BETA(in)
         ENDDO
        ENDDO
       ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Print -CK12nd
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       DO i=1,NBF5
        DO j=1,i-1
         if(DABS(CK12nd(i,j))>THRESHCJK)then        
          WRITE(13)i,j,j,i,-CK12nd(i,j)
          WRITE(13)j,i,i,j,-CK12nd(i,j)
         end if
        ENDDO
       ENDDO
       WRITE(13)0,0,0,0,0.0d0       
       DEALLOCATE(FIs,BETA,CK12nd)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END IF
!-----------------------------------------------------------------------
      RETURN
      END

! OUTPUTTijab_rc
      SUBROUTINE OUTPUTTijab_rc(NOC,NVI,NN,Tijab)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/INPNOF_ORBSPACE0/NO1,NDOC,NCO,NCWO,NVIR,NAC,NO0
      COMMON/INPNOF_Tijab/NOUTTijab,NTHRESHTijab,THRESHTijab      
      DOUBLE PRECISION,DIMENSION(NN)::Tijab
!-----------------------------------------------------------------------  
      DO i=1,NOC
       DO j=1,NOC
         DO k=1,NVI
          ki = (k-1)*noc + i
          DO l=1,NVI
           ijkl =  i + (j-1)*noc + (k-1)*noc*noc + (l-1)*noc*noc*nvi
           ijlk =  i + (j-1)*noc + (l-1)*noc*noc + (k-1)*noc*noc*nvi
           TTijab = 2.0*Tijab(ijkl)-Tijab(ijlk)
           IF(DABS(TTijab)>THRESHTijab)THEN
            in = NO1+i
            jn = NO1+j
            kn = NO1+noc+k
            ln = NO1+noc+l
            WRITE(13)in,jn,kn,ln,TTijab
           ENDIF
          ENDDO
         ENDDO
       ENDDO
      ENDDO
      WRITE(13)0,0,0,0,0.0d0
!-----------------------------------------------------------------------
      RETURN
      END

! OUTPUTRDMrc
      SUBROUTINE OUTPUTRDMrc(OVERLAP,RO,QD,CJ12,CK12)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      COMMON/INPNOF_ARDM/THRESHDM,NOUTRDM,NSQT,NTHRESHDM      
!
      DOUBLE PRECISION,DIMENSION(NBF5)::RO
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5)::CJ12,CK12
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::OVERLAP
      DOUBLE PRECISION,DIMENSION(NBF,NBF,NBF)::QD
      ALLOCATABLE::IDXD(:),BUFD(:)
!-----------------------------------------------------------------------
!     2DM
!-----------------------------------------------------------------------
      IF( NOUTRDM==2 .or. NOUTRDM==3 )THEN
       LBUFDM = 10000       ! Length of buffer for Density Matrix values
       IF(NSQT==1)THEN
!       ifort -i8 in PNOF, MAX=2**64-1, MAX_NBF=2**16-1=65535 !!!
        NSHI1=2**16
!       NSH2=2**32
        NSHI2=NSHI1*NSHI1
!       NSH3=2**48
        NSHI3=NSHI2*NSHI1
        ALLOCATE (IDXD(LBUFDM),BUFD(LBUFDM))
        NBINT=0
        NREC=0
        DO IETA=1,NBF
         DO IMIU=1,NBF
          DO INIU=1,NBF
           DO ILAM=1,NBF
            RDM2 = SUMDDL(IETA,IMIU,INIU,ILAM,RO,QD,CJ12,CK12)
            IF(DABS(RDM2)>THRESHDM)THEN
             NBINT=NBINT+1
!            CHANGE INDEXES ( IMIU <-> INIU ) for Ugalde's Program
             IDXD(NBINT)=NSHI3*IETA+NSHI2*INIU+NSHI1*IMIU+ILAM
             BUFD(NBINT)=RDM2
!            Write: IDXD(IETA,INIU,IMIU,ILAM) ; BUFD(RDM2)
             IF(NBINT==LBUFDM)THEN
              WRITE(14)NBINT,(IDXD(M),M=1,NBINT),(BUFD(M),M=1,NBINT)
              NBINT=0
              NREC=NREC+1
             ENDIF
            ENDIF
           ENDDO
          ENDDO
         ENDDO
        ENDDO
!       Write the last record on file 14
        IF(NBINT/=0)THEN
         WRITE(14)NBINT,(IDXD(M),M=1,NBINT),(BUFD(M),M=1,NBINT)
         NREC=NREC+1
        ENDIF
!       Total number of written integrals
        NBINT=NBINT+(NREC-1)*LBUFDM
!       Write the number of records and 2DMS
        WRITE(16,3)NREC,NBINT
        DEALLOCATE (IDXD,BUFD)
       ELSEIF(NSQT==0)THEN
        DO IETA=1,NBF
         DO IMIU=1,NBF
          DO INIU=1,NBF
           DO ILAM=1,NBF
            RDM2 = SUMDDL(IETA,IMIU,INIU,ILAM,RO,QD,CJ12,CK12)
            IF(DABS(RDM2)>THRESHDM)THEN
!            Change indexes ( IMIU <-> INIU ) for Ugalde's Program
!            WRITE(14,1)IETA,INIU,IMIU,ILAM,RDM2
             WRITE(14,1)IETA,IMIU,INIU,ILAM,RDM2
            ENDIF
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDIF
      ENDIF
!-----------------------------------------------------------------------
!     1DM
!-----------------------------------------------------------------------
      IHAODE = 1
      HAODE = 0.0d0
      IF( NOUTRDM==1 .or. NOUTRDM==3 )THEN
       DO IETA=1,NBF
        DO IMIU=1,NBF
         RDM1 = SUMDL(IETA,IMIU,RO,QD)
         IF(IHAODE==1 .and. IETA/=IMIU)THEN     ! Off-diagonal elements
          HAODE = HAODE + RDM1*RDM1
         ENDIF
         IF(DABS(RDM1)>THRESHDM)THEN
          WRITE(15,2)IETA,IMIU,RDM1
         ENDIF
        ENDDO
       ENDDO
       NBF2 = NBF*(NBF-1)
       HAODE = DSQRT ( HAODE / DFLOAT(NBF2) )
      ENDIF
!     Check the normalization of the RDMs
      WRITE(6,4)
      IF(NOUTRDM==1.or.NOUTRDM==3)THEN
       CALL RDM1NORM(OVERLAP,RO,QD)
       WRITE(6,5)HAODE
      END IF
      IF(NOUTRDM==2.or.NOUTRDM==3)CALL RDM2NORM(OVERLAP,RO,QD,CJ12,CK12)
!-----------------------------------------------------------------------
1     FORMAT(4I4,D20.10)
2     FORMAT(2I4,D20.10)
3     FORMAT(I4,I20)
4     FORMAT(/' RDM Norms ',/,' ---------  ')
5     FORMAT(1X,'Harmonic average of atomic 1RDM =',F10.5,/)
      RETURN
      END
  
! SUMDDL (2RDM)
      FUNCTION SUMDDL(IETA,IMIU,INIU,ILAM,RO,QD,CJ12,CK12)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      COMMON/INPNOF_ARDM/THRESHDM,NOUTRDM,NSQT,NTHRESHDM      
!
      DOUBLE PRECISION,DIMENSION(NBF5)::RO
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5)::CJ12,CK12
      DOUBLE PRECISION,DIMENSION(NBF,NBF,NBF)::QD
      ALLOCATABLE::Dj(:,:),Di(:,:)
      ALLOCATE (Dj(NBF,NBF),Di(NBF,NBF))
!-----------------------------------------------------------------------
      SUMDDL=0.0d0
      DO j=1,NBF5
       Dj = QD(j,1:NBF,1:NBF)
       DO i=1,NBF5
        Di = QD(i,1:NBF,1:NBF)
        IF(i==j)THEN
         SUMDDL = SUMDDL + RO(j)*Dj(IETA,INIU)*Di(IMIU,ILAM)
        ELSE
         SUMDDL = SUMDDL + CJ12(j,i)*Dj(IETA,INIU)*Di(IMIU,ILAM)        &
                         - CK12(j,i)*Dj(IETA,ILAM)*Di(IMIU,INIU)        
!L                       - CK12(j,i)*Dj(IETA,IMIU)*Di(INIU,ILAM)
        ENDIF
       ENDDO
      ENDDO
      DO j=NB+1,NA
       Dj = QD(j,1:NBF,1:NBF)
       SUMDDL = SUMDDL - RO(j)*Dj(IETA,INIU)*Dj(IMIU,ILAM)
      ENDDO
      
      SUMDDL = SUMDDL/4.0d0
      IF(DABS(SUMDDL)<=THRESHDM)SUMDDL = 0.0d0
      DEALLOCATE (Dj,Di)
!-----------------------------------------------------------------------
      RETURN
      END

! SUMDL (1RDM)
      FUNCTION SUMDL(IETA,IMIU,RO,QD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL HighSpin
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
!
      DOUBLE PRECISION,DIMENSION(NBF5)::RO
      DOUBLE PRECISION,DIMENSION(NBF,NBF,NBF)::QD
!-----------------------------------------------------------------------
      SUMDL=0.0d0
      DO J=1,NB
       SUMDL = SUMDL + RO(J) * QD(J,IETA,IMIU)
      ENDDO
      IF(NSOC>0)THEN
       DO J=NB+1,NA
       if(HighSpin)THEN
        SUMDL = SUMDL + 0.5d0 * RO(J) * QD(J,IETA,IMIU)
       else
        SUMDL = SUMDL + RO(J) * QD(J,IETA,IMIU)
       endif
       ENDDO
      ENDIF
      DO J=NA+1,NBF5
       SUMDL = SUMDL + RO(J) * QD(J,IETA,IMIU)
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! RDM1NORM
      SUBROUTINE RDM1NORM(S,RO,QD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
!
      DOUBLE PRECISION,DIMENSION(NBF5)::RO
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::S
      DOUBLE PRECISION,DIMENSION(NBF,NBF,NBF)::QD
!-----------------------------------------------------------------------
      SUMNORM=0.0d0
      DO IETA=1,NBF
       DO IMIU=1,IETA
        SUMNORM=SUMNORM+2.0*S(IETA,IMIU)*SUMDL(IETA,IMIU,RO,QD)
       ENDDO
       SUMNORM=SUMNORM-S(IETA,IETA)*SUMDL(IETA,IETA,RO,QD)
      ENDDO
      WRITE(6,1)SUMNORM
!-----------------------------------------------------------------------     
1     FORMAT(/,1X,'1RDM Norm =',F10.3,/)
      RETURN
      END

! RDM2NORM
      SUBROUTINE RDM2NORM(S,RO,QD,CJ12,CK12)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
!
      DOUBLE PRECISION,DIMENSION(NBF5)::RO
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5)::CJ12,CK12
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::S
      DOUBLE PRECISION,DIMENSION(NBF,NBF,NBF)::QD
!-----------------------------------------------------------------------
      TWONORMA=0.D0
      DO I=1,NBF
       DO K=1,I
        DO J=1,NBF
         DO L=1,J
          DENSITY = DENSI(I,J,K,L,RO,QD,CJ12,CK12)
          TWONORMA=TWONORMA+DENSITY*S(I,K)*S(J,L)
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      WRITE(6,1)TWONORMA
1     FORMAT(1X,'2RDM Norm =',F7.3)
!-----------------------------------------------------------------------
      RETURN
      END

! DENSI
      DOUBLE PRECISION FUNCTION DENSI(II,JJ,KK,LL,RO,QD,CJ12,CK12)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
!
      DOUBLE PRECISION,DIMENSION(NBF5)::RO
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5)::CJ12,CK12
      DOUBLE PRECISION,DIMENSION(NBF,NBF,NBF)::QD
!-----------------------------------------------------------------------
      IF((II==KK).AND.(JJ==LL))THEN
       DENSI=SUMDDL(ii,jj,kk,ll,RO,QD,CJ12,CK12)
       RETURN
      ELSE IF((II==KK).AND.(JJ/=LL))THEN
       DENSI = SUMDDL(ii,jj,kk,ll,RO,QD,CJ12,CK12)                      &
             + SUMDDL(ii,ll,kk,jj,RO,QD,CJ12,CK12)                       
       RETURN                                                            
      ELSE IF((II/=KK).AND.(JJ==LL))THEN                                 
       DENSI = SUMDDL(ii,jj,kk,ll,RO,QD,CJ12,CK12)                      &
             + SUMDDL(kk,jj,ii,ll,RO,QD,CJ12,CK12)                       
       RETURN                                                            
      ELSE IF((II/=KK).AND.(JJ/=LL))THEN                                 
       DENSI = SUMDDL(ii,jj,kk,ll,RO,QD,CJ12,CK12)                      &
             + SUMDDL(ii,ll,kk,jj,RO,QD,CJ12,CK12)                      &
             + SUMDDL(kk,jj,ii,ll,RO,QD,CJ12,CK12)                      &
             + SUMDDL(kk,ll,ii,jj,RO,QD,CJ12,CK12)
       RETURN
      ELSE
       WRITE(6,*)' *** ERROR IN DENSI SUBROUTINE ***'
       STOP
      END IF
!----------------------------------------------------------------------- 
      END

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
!                                                                      !
!             E L E C T R O S T A T I C     M O M E N T S              !
!                                                                      !
!      2015 Electrostatic dipole, quadrupole and octupole moments      !
!                                                                      !
!                     implemented by Ion Mitxelena                     !
!                                                                      !
!                 ( J. Chem. Phys. 144, 204108, 2016 )                 !
!                                                                      !
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
!                                                                      !
!   DIPMOMr:  Electronic, nuclear and total dipole moments             !
!   QUADMOMr: Electronic, nuclear and total quadrupole moments         !
!   OCTMOMr:  Electronic, nuclear and total octupole moments           !
!                                                                      !
!   PASSDIPUSER:  Pass atomic dipole matrices to USER(13,14,15,16)     !
!   PASSUSERDIP:  Pass USER(13,14,15,16) to atomic dipole matrices     !
!   PASSQUADUSER: Pass atomic quadrupole matrices to USER(18-24)       !
!   PASSUSERQUAD: Pass USER(18-24) to atomic quadrupole matrices       !
!   PASSOCTUSER:  Pass atomic octupole matrices to USER(31-41)         !
!   PASSUSEROCT:  Pass USER(31-41) to atomic octupole matrices         !
!                                                                      !
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
                                                                       
! DIPMOMr
      SUBROUTINE DIPMOMr(DIPN,ADIPx,ADIPy,ADIPz,DIPx,DIPy,DIPz,QD,RO,   &
                         DM1e,DM2e,DM3e,DM1,DM2,DM3,DTOTAL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL HighSpin
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      !
      DOUBLE PRECISION,DIMENSION(3),INTENT(IN)::DIPN
      DOUBLE PRECISION,DIMENSION(NBF5),INTENT(IN)::RO,DIPx,DIPy,DIPz
      DOUBLE PRECISION,DIMENSION(NBF,NBF),INTENT(IN)::ADIPx,ADIPy,ADIPz
      DOUBLE PRECISION,DIMENSION(NBF,NBF,NBF),INTENT(IN)::QD
      DOUBLE PRECISION,INTENT(OUT)::DM1,DM2,DM3,DTOTAL
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::AUX
!-----------------------------------------------------------------------
!     Electron Contribution to Dipole Moment
!-----------------------------------------------------------------------
      DM1e = 0.0
      DM2e = 0.0
      DM3e = 0.0
!
      ALLOCATE (AUX(NBF,NBF))
      DO J=1,NB
       AUX(1:NBF,1:NBF) = QD(J,1:NBF,1:NBF)
       CALL TRACEm(DIPx(J),AUX,ADIPx,NBF)
       DM1e = DM1e - 2.0*RO(J)*DIPx(J)
       CALL TRACEm(DIPy(J),AUX,ADIPy,NBF)
       DM2e = DM2e - 2.0*RO(J)*DIPy(J)
       CALL TRACEm(DIPz(J),AUX,ADIPz,NBF)
       DM3e = DM3e - 2.0*RO(J)*DIPz(J)
      ENDDO
!      
      IF(NSOC>0)THEN
       if(.not.HighSpin)then
        DO J=NB+1,NA
         AUX(1:NBF,1:NBF) = QD(J,1:NBF,1:NBF)
         CALL TRACEm(DIPx(J),AUX,ADIPx,NBF)
         DM1e = DM1e - 2.0*RO(J)*DIPx(J)
         CALL TRACEm(DIPy(J),AUX,ADIPy,NBF)
         DM2e = DM2e - 2.0*RO(J)*DIPy(J)
         CALL TRACEm(DIPz(J),AUX,ADIPz,NBF)
         DM3e = DM3e - 2.0*RO(J)*DIPz(J)
        ENDDO
       else if(HighSpin)then
        DO J=NB+1,NA
         AUX(1:NBF,1:NBF) = QD(J,1:NBF,1:NBF)
         CALL TRACEm(DIPx(J),AUX,ADIPx,NBF)
         DM1e = DM1e - RO(J)*DIPx(J)
         CALL TRACEm(DIPy(J),AUX,ADIPy,NBF)
         DM2e = DM2e - RO(J)*DIPy(J)
         CALL TRACEm(DIPz(J),AUX,ADIPz,NBF)
         DM3e = DM3e - RO(J)*DIPz(J)
        ENDDO
       end if      
      ENDIF
!       
      DO J=NA+1,NBF5
       AUX(1:NBF,1:NBF) = QD(J,1:NBF,1:NBF)
       CALL TRACEm(DIPx(J),AUX,ADIPx,NBF)
       DM1e = DM1e - 2.0*RO(J)*DIPx(J)
       CALL TRACEm(DIPy(J),AUX,ADIPy,NBF)
       DM2e = DM2e - 2.0*RO(J)*DIPy(J)
       CALL TRACEm(DIPz(J),AUX,ADIPz,NBF)
       DM3e = DM3e - 2.0*RO(J)*DIPz(J)
      ENDDO
      DEALLOCATE (AUX)
!-----------------------------------------------------------------------
!     Nuclear Contribution to Dipole Moment
!-----------------------------------------------------------------------
      DM1 = DIPN(1) + DM1e
      DM2 = DIPN(2) + DM2e
      DM3 = DIPN(3) + DM3e
!-----------------------------------------------------------------------
!     Total Dipole Moment
!-----------------------------------------------------------------------
      DTOTAL = SQRT(DM1*DM1+DM2*DM2+DM3*DM3)
      RETURN
!-----------------------------------------------------------------------
      END SUBROUTINE DIPMOMr

! QUADMOMr
      SUBROUTINE QUADMOMr(QUADN,AQUADxx,AQUADyy,AQUADzz,AQUADxy,        &
                          AQUADxz,AQUADyz,QUADxx,QUADyy,QUADzz,QUADxy,  &
                          QUADxz,QUADyz,QD,RO,QM1e,QM2e,QM3e,QM4e,QM5e, &
                          QM6e,QM1,QM2,QM3,QM4,QM5,QM6)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      LOGICAL HighSpin
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
!
      INTEGER :: J
      DOUBLE PRECISION, DIMENSION(NBF,NBF,NBF), INTENT(IN) :: QD
      DOUBLE PRECISION, DIMENSION(NBF5), INTENT(IN) :: RO(NBF5)
      DOUBLE PRECISION, DIMENSION(6), INTENT(IN) :: QUADN
      DOUBLE PRECISION, DIMENSION(NBF,NBF), INTENT(IN) :: AQUADxx,AQUADyy,AQUADzz
      DOUBLE PRECISION, DIMENSION(NBF,NBF), INTENT(IN) :: AQUADxy,AQUADxz,AQUADyz
      DOUBLE PRECISION, DIMENSION(NBF5) :: QUADxx,QUADyy,QUADzz,QUADxy,QUADxz,QUADyz
      DOUBLE PRECISION, INTENT(OUT) :: QM1,QM2,QM3,QM4,QM5,QM6
      DOUBLE PRECISION :: QM1e,QM2e,QM3e,QM4e,QM5e,QM6e
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::AUX
!-----------------------------------------------------------------------
!     Electron Contribution to Quadrupole Moment
!-----------------------------------------------------------------------
      QM1e=0.0d0
      QM2e=0.0d0
      QM3e=0.0d0
      QM4e=0.0d0
      QM5e=0.0d0
      QM6e=0.0d0
!
      ALLOCATE (AUX(NBF,NBF))
      DO J=1,NB
       AUX(1:NBF,1:NBF) = QD(J,1:NBF,1:NBF)
       CALL TRACEm(QUADxx(J),AUX,AQUADxx,NBF)
       QM1e = QM1e - 2.d0*RO(J)*QUADxx(J)
       CALL TRACEm(QUADyy(J),AUX,AQUADyy,NBF)
       QM2e = QM2e - 2.d0*RO(J)*QUADyy(J)
       CALL TRACEm(QUADzz(J),AUX,AQUADzz,NBF)
       QM3e = QM3e - 2.d0*RO(J)*QUADzz(J)
       CALL TRACEm(QUADxy(J),AUX,AQUADxy,NBF)
       QM4e = QM4e - 2.d0*RO(J)*QUADxy(J)
       CALL TRACEm(QUADxz(J),AUX,AQUADxz,NBF)
       QM5e = QM5e - 2.d0*RO(J)*QUADxz(J)
       CALL TRACEm(QUADyz(J),AUX,AQUADyz,NBF)
       QM6e = QM6e - 2.d0*RO(J)*QUADyz(J)
      END DO
!      
      IF(NSOC>0)THEN
       if(.not.HighSpin)then
        DO J=NB+1,NA
         AUX(1:NBF,1:NBF) = QD(J,1:NBF,1:NBF)
         CALL TRACEm(QUADxx(J),AUX,AQUADxx,NBF)
         QM1e = QM1e - 2.d0*RO(J)*QUADxx(J)
         CALL TRACEm(QUADyy(J),AUX,AQUADyy,NBF)
         QM2e = QM2e - 2.d0*RO(J)*QUADyy(J)
         CALL TRACEm(QUADzz(J),AUX,AQUADzz,NBF)
         QM3e = QM3e - 2.d0*RO(J)*QUADzz(J)
         CALL TRACEm(QUADxy(J),AUX,AQUADxy,NBF)
         QM4e = QM4e - 2.d0*RO(J)*QUADxy(J)
         CALL TRACEm(QUADxz(J),AUX,AQUADxz,NBF)
         QM5e = QM5e - 2.d0*RO(J)*QUADxz(J)
         CALL TRACEm(QUADyz(J),AUX,AQUADyz,NBF)
         QM6e = QM6e - 2.d0*RO(J)*QUADyz(J)
        END DO
       else if(HighSpin)then        
        DO J=NB+1,NA
         AUX(1:NBF,1:NBF) = QD(J,1:NBF,1:NBF)
         CALL TRACEm(QUADxx(J),AUX,AQUADxx,NBF)
         QM1e = QM1e - RO(J)*QUADxx(J)
         CALL TRACEm(QUADyy(J),AUX,AQUADyy,NBF)
         QM2e = QM2e - RO(J)*QUADyy(J)
         CALL TRACEm(QUADzz(J),AUX,AQUADzz,NBF)
         QM3e = QM3e - RO(J)*QUADzz(J)
         CALL TRACEm(QUADxy(J),AUX,AQUADxy,NBF)
         QM4e = QM4e - RO(J)*QUADxy(J)
         CALL TRACEm(QUADxz(J),AUX,AQUADxz,NBF)
         QM5e = QM5e - RO(J)*QUADxz(J)
         CALL TRACEm(QUADyz(J),AUX,AQUADyz,NBF)
         QM6e = QM6e - RO(J)*QUADyz(J)
        END DO
       end if
      END IF
!
      DO J=NA+1,NBF5
       AUX(1:NBF,1:NBF) = QD(J,1:NBF,1:NBF)
       CALL TRACEm(QUADxx(J),AUX,AQUADxx,NBF)
       QM1e = QM1e - 2.d0*RO(J)*QUADxx(J)
       CALL TRACEm(QUADyy(J),AUX,AQUADyy,NBF)
       QM2e = QM2e - 2.d0*RO(J)*QUADyy(J)
       CALL TRACEm(QUADzz(J),AUX,AQUADzz,NBF)
       QM3e = QM3e - 2.d0*RO(J)*QUADzz(J)
       CALL TRACEm(QUADxy(J),AUX,AQUADxy,NBF)
       QM4e = QM4e - 2.d0*RO(J)*QUADxy(J)
       CALL TRACEm(QUADxz(J),AUX,AQUADxz,NBF)
       QM5e = QM5e - 2.d0*RO(J)*QUADxz(J)
       CALL TRACEm(QUADyz(J),AUX,AQUADyz,NBF)
       QM6e = QM6e - 2.d0*RO(J)*QUADyz(J)
      END DO
      DEALLOCATE (AUX)
!-----------------------------------------------------------------------
!     Add Nuclear Contribution to Quadrupole Moment
!-----------------------------------------------------------------------
      QM1e = QUADN(1) + QM1e
      QM2e = QUADN(2) + QM2e
      QM3e = QUADN(3) + QM3e
      QM4e = QUADN(4) + QM4e
      QM5e = QUADN(5) + QM5e
      QM6e = QUADN(6) + QM6e
!-----------------------------------------------------------------------
!     Form Quadrupole Tensor (BUCKINHAM)
!-----------------------------------------------------------------------
      QM1 = 0.5D0*(QM1e + QM1e - QM2e - QM3e)
      QM2 = 0.5D0*(QM2e + QM2e - QM1e - QM3e)
      QM3 = 0.5D0*(QM3e + QM3e - QM1e - QM2e)
      QM4 = 1.5D0*QM4e
      QM5 = 1.5D0*QM5e
      QM6 = 1.5D0*QM6e
!-----------------------------------------------------------------------    
      RETURN
      END SUBROUTINE QUADMOMr

! OCTMOMr
      SUBROUTINE OCTMOMr(OCTUN,AOCTxxx,AOCTyyy,AOCTzzz,                 &
                         AOCTxxy,AOCTxxz,AOCTxyy,AOCTyyz,AOCTxzz,       &
                         AOCTyzz,AOCTxyz,OCTxxx,OCTyyy,OCTzzz,          &
                         OCTxxy,OCTxxz,OCTxyy,OCTyyz,OCTxzz,            &
                         OCTyzz,OCTxyz,                                 &
                         QD,RO,OMXXXe,OMYYYe,OMZZZe,                    &
                         OMXXYe,OMXXZe,OMXYYe,OMYYZe,OMXZZe,            &
                         OMYZZe,OMXYZe,OXXX,OYYY,OZZZ,                  &
                         OXXY,OXXZ,OXYY,OYYZ,OXZZ,                      &
                         OYZZ,OXYZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL HighSpin
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
!
      INTEGER :: J
      DOUBLE PRECISION, DIMENSION(NBF,NBF,NBF), INTENT(IN) :: QD
      DOUBLE PRECISION, DIMENSION(NBF5), INTENT(IN) :: RO(NBF5)
      DOUBLE PRECISION, DIMENSION(10), INTENT(IN) :: OCTUN(10)
      DOUBLE PRECISION, DIMENSION(NBF,NBF), INTENT(IN) :: AOCTxxx,AOCTyyy,AOCTzzz
      DOUBLE PRECISION, DIMENSION(NBF,NBF), INTENT(IN) :: AOCTxyy,AOCTyyz,AOCTxzz
      DOUBLE PRECISION, DIMENSION(NBF,NBF), INTENT(IN) :: AOCTxxy,AOCTxxz
      DOUBLE PRECISION, DIMENSION(NBF,NBF), INTENT(IN) :: AOCTyzz,AOCTxyz
      DOUBLE PRECISION, DIMENSION(NBF5) :: OCTxxx,OCTyyy,OCTzzz,OCTxxy,OCTxxz
      DOUBLE PRECISION, DIMENSION(NBF5) :: OCTxyy,OCTyyz,OCTxzz,OCTyzz,OCTxyz
      DOUBLE PRECISION :: OMXXXe,OMYYYe,OMZZZe,OMXXYe,OMXXZe
      DOUBLE PRECISION :: OMXYYe,OMYYZe,OMXZZe,OMYZZe,OMXYZe
      DOUBLE PRECISION, INTENT(OUT) :: OXXX,OYYY,OZZZ,OXXY,OXXZ
      DOUBLE PRECISION, INTENT(OUT) :: OXYY,OYYZ,OXZZ,OYZZ,OXYZ
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::AUX
!-----------------------------------------------------------------------
!     Electron Contribution to Quadrupole Moment
!-----------------------------------------------------------------------
      OMXXXe=0.0d0
      OMYYYe=0.0d0
      OMZZZe=0.0d0
      OMXXYe=0.0d0
      OMXXZe=0.0d0
      OMXYYe=0.0d0
      OMYYZe=0.0d0
      OMXZZe=0.0d0
      OMYZZe=0.0d0
      OMXYZe=0.0d0
!
      ALLOCATE (AUX(NBF,NBF))
      DO J=1,NB
       AUX(1:NBF,1:NBF) = QD(J,1:NBF,1:NBF)
       CALL TRACEm(OCTxxx(J),AUX,AOCTxxx,NBF)
       OMXXXe = OMXXXe - 2.d0*RO(J)*OCTxxx(J)
       CALL TRACEm(OCTyyy(J),AUX,AOCTyyy,NBF)
       OMYYYe = OMYYYe - 2.d0*RO(J)*OCTyyy(J)
       CALL TRACEm(OCTzzz(J),AUX,AOCTzzz,NBF)
       OMZZZe = OMZZZe - 2.d0*RO(J)*OCTzzz(J)
       CALL TRACEm(OCTxxy(J),AUX,AOCTxxy,NBF)
       OMXXYe = OMXXYe - 2.d0*RO(J)*OCTxxy(J)
       CALL TRACEm(OCTxxz(J),AUX,AOCTxxz,NBF)
       OMXXZe = OMXXZe - 2.d0*RO(J)*OCTxxz(J)
       CALL TRACEm(OCTxyy(J),AUX,AOCTxyy,NBF)
       OMXYYe = OMXYYe - 2.d0*RO(J)*OCTxyy(J)
       CALL TRACEm(OCTyyz(J),AUX,AOCTyyz,NBF)
       OMYYZe = OMYYZe - 2.d0*RO(J)*OCTyyz(J)
       CALL TRACEm(OCTxzz(J),AUX,AOCTxzz,NBF)
       OMXZZe = OMXZZe - 2.d0*RO(J)*OCTxzz(J)
       CALL TRACEm(OCTyzz(J),AUX,AOCTyzz,NBF)
       OMYZZe = OMYZZe - 2.d0*RO(J)*OCTyzz(J)
       CALL TRACEm(OCTxyz(J),AUX,AOCTxyz,NBF)
       OMXYZe = OMXYZe - 2.d0*RO(J)*OCTxyz(J)
      END DO
!      
      IF(NSOC>0)THEN
       if(.not.HighSpin)then
        DO J=NB+1,NA
         AUX(1:NBF,1:NBF) = QD(J,1:NBF,1:NBF)
         CALL TRACEm(OCTxxx(J),AUX,AOCTxxx,NBF)
         OMXXXe = OMXXXe - 2.d0*RO(J)*OCTxxx(J)
         CALL TRACEm(OCTyyy(J),AUX,AOCTyyy,NBF)
         OMYYYe = OMYYYe - 2.d0*RO(J)*OCTyyy(J)
         CALL TRACEm(OCTzzz(J),AUX,AOCTzzz,NBF)
         OMZZZe = OMZZZe - 2.d0*RO(J)*OCTzzz(J)
         CALL TRACEm(OCTxxy(J),AUX,AOCTxxy,NBF)
         OMXXYe = OMXXYe - 2.d0*RO(J)*OCTxxy(J)
         CALL TRACEm(OCTxxz(J),AUX,AOCTxxz,NBF)
         OMXXZe = OMXXZe - 2.d0*RO(J)*OCTxxz(J)
         CALL TRACEm(OCTxyy(J),AUX,AOCTxyy,NBF)
         OMXYYe = OMXYYe - 2.d0*RO(J)*OCTxyy(J)
         CALL TRACEm(OCTyyz(J),AUX,AOCTyyz,NBF)
         OMYYZe = OMYYZe - 2.d0*RO(J)*OCTyyz(J)
         CALL TRACEm(OCTxzz(J),AUX,AOCTxzz,NBF)
         OMXZZe = OMXZZe - 2.d0*RO(J)*OCTxzz(J)
         CALL TRACEm(OCTyzz(J),AUX,AOCTyzz,NBF)
         OMYZZe = OMYZZe - 2.d0*RO(J)*OCTyzz(J)
         CALL TRACEm(OCTxyz(J),AUX,AOCTxyz,NBF)
         OMXYZe = OMXYZe - 2.d0*RO(J)*OCTxyz(J)
        END DO
       else if(HighSpin)then
        DO J=NB+1,NA
         AUX(1:NBF,1:NBF) = QD(J,1:NBF,1:NBF)
         CALL TRACEm(OCTxxx(J),AUX,AOCTxxx,NBF)
         OMXXXe = OMXXXe - RO(J)*OCTxxx(J)
         CALL TRACEm(OCTyyy(J),AUX,AOCTyyy,NBF)
         OMYYYe = OMYYYe - RO(J)*OCTyyy(J)
         CALL TRACEm(OCTzzz(J),AUX,AOCTzzz,NBF)
         OMZZZe = OMZZZe - RO(J)*OCTzzz(J)
         CALL TRACEm(OCTxxy(J),AUX,AOCTxxy,NBF)
         OMXXYe = OMXXYe - RO(J)*OCTxxy(J)
         CALL TRACEm(OCTxxz(J),AUX,AOCTxxz,NBF)
         OMXXZe = OMXXZe - RO(J)*OCTxxz(J)
         CALL TRACEm(OCTxyy(J),AUX,AOCTxyy,NBF)
         OMXYYe = OMXYYe - RO(J)*OCTxyy(J)
         CALL TRACEm(OCTyyz(J),AUX,AOCTyyz,NBF)
         OMYYZe = OMYYZe - RO(J)*OCTyyz(J)
         CALL TRACEm(OCTxzz(J),AUX,AOCTxzz,NBF)
         OMXZZe = OMXZZe - RO(J)*OCTxzz(J)
         CALL TRACEm(OCTyzz(J),AUX,AOCTyzz,NBF)
         OMYZZe = OMYZZe - RO(J)*OCTyzz(J)
         CALL TRACEm(OCTxyz(J),AUX,AOCTxyz,NBF)
         OMXYZe = OMXYZe - RO(J)*OCTxyz(J)
        END DO
       end if
      END IF
!       
      DO J=1,NBF5
       AUX(1:NBF,1:NBF) = QD(J,1:NBF,1:NBF)
       CALL TRACEm(OCTxxx(J),AUX,AOCTxxx,NBF)
       OMXXXe = OMXXXe - 2.d0*RO(J)*OCTxxx(J)
       CALL TRACEm(OCTyyy(J),AUX,AOCTyyy,NBF)
       OMYYYe = OMYYYe - 2.d0*RO(J)*OCTyyy(J)
       CALL TRACEm(OCTzzz(J),AUX,AOCTzzz,NBF)
       OMZZZe = OMZZZe - 2.d0*RO(J)*OCTzzz(J)
       CALL TRACEm(OCTxxy(J),AUX,AOCTxxy,NBF)
       OMXXYe = OMXXYe - 2.d0*RO(J)*OCTxxy(J)
       CALL TRACEm(OCTxxz(J),AUX,AOCTxxz,NBF)
       OMXXZe = OMXXZe - 2.d0*RO(J)*OCTxxz(J)
       CALL TRACEm(OCTxyy(J),AUX,AOCTxyy,NBF)
       OMXYYe = OMXYYe - 2.d0*RO(J)*OCTxyy(J)
       CALL TRACEm(OCTyyz(J),AUX,AOCTyyz,NBF)
       OMYYZe = OMYYZe - 2.d0*RO(J)*OCTyyz(J)
       CALL TRACEm(OCTxzz(J),AUX,AOCTxzz,NBF)
       OMXZZe = OMXZZe - 2.d0*RO(J)*OCTxzz(J)
       CALL TRACEm(OCTyzz(J),AUX,AOCTyzz,NBF)
       OMYZZe = OMYZZe - 2.d0*RO(J)*OCTyzz(J)
       CALL TRACEm(OCTxyz(J),AUX,AOCTxyz,NBF)
       OMXYZe = OMXYZe - 2.d0*RO(J)*OCTxyz(J)
      END DO
      DEALLOCATE (AUX)
!-----------------------------------------------------------------------
!     Nuclear Contribution to Quadrupole Moment
!-----------------------------------------------------------------------
      OMXXXe = OCTUN(1) + OMXXXe
      OMYYYe = OCTUN(2) + OMYYYe
      OMZZZe = OCTUN(3) + OMZZZe
      OMXXYe = OCTUN(4) + OMXXYe
      OMXXZe = OCTUN(5) + OMXXZe
      OMXYYe = OCTUN(6) + OMXYYe
      OMYYZe = OCTUN(7) + OMYYZe
      OMXZZe = OCTUN(8) + OMXZZe
      OMYZZe = OCTUN(9) + OMYZZe
      OMXYZe = OCTUN(10) + OMXYZe
!-----------------------------------------------------------------------
!     Form Octupole Tensor
!-----------------------------------------------------------------------
      OXXX = OMXXXe - 1.5D+00*OMXYYe - 1.5D+00*OMXZZe
      OYYY = OMYYYe - 1.5D+00*OMXXYe - 1.5D+00*OMYZZe
      OZZZ = OMZZZe - 1.5D+00*OMXXZe - 1.5D+00*OMYYZe
      OXXY = 2.0D+00*OMXXYe - 0.5D+00*OMYYYe - 0.5D+00*OMYZZe
      OXXZ = 2.0D+00*OMXXZe - 0.5D+00*OMYYZe - 0.5D+00*OMZZZe
      OXYY = 2.0D+00*OMXYYe - 0.5D+00*OMXXXe - 0.5D+00*OMXZZe
      OYYZ = 2.0D+00*OMYYZe - 0.5D+00*OMXXZe - 0.5D+00*OMZZZe
      OXZZ = 2.0D+00*OMXZZe - 0.5D+00*OMXXXe - 0.5D+00*OMXYYe
      OYZZ = 2.0D+00*OMYZZe - 0.5D+00*OMXXYe - 0.5D+00*OMYYYe
      OXYZ = 2.5D+00*OMXYZe
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE OCTMOMr
      
! PASSDIPUSER
      SUBROUTINE PASSDIPUSER(DIPN,ADIPx,ADIPy,ADIPz,USER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/PUNTEROSUSER/N1,N2,N3,N4,N5,N6,N7,N8,N9,N10,N11,N12,N13,   &
                          N14,N15,N16,N17,N18,N19,N20,N21,N22,N23,N24,  &
                          N25,N26,N27,N28,N29,N30,N31,N32,N33,N34,N35,  &
                          N36,N37,N38,N39,N40,N41,N42,N43,N44,N45,N46,  &
                          N47,N48,N49,N50,N51,NUSER
!
      DOUBLE PRECISION,DIMENSION(3)::DIPN
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::ADIPx,ADIPy,ADIPz
      DOUBLE PRECISION,DIMENSION(NUSER)::USER
!-----------------------------------------------------------------------
      CALL XtoX0(DIPN ,USER(N11),  3)
      CALL XtoX0(ADIPx,USER(N12),NSQ)
      CALL XtoX0(ADIPy,USER(N13),NSQ)
      CALL XtoX0(ADIPz,USER(N14),NSQ)
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE PASSDIPUSER
      
! PASSUSERDIP
      SUBROUTINE PASSUSERDIP(DIPN,ADIPx,ADIPy,ADIPz,USER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/PUNTEROSUSER/N1,N2,N3,N4,N5,N6,N7,N8,N9,N10,N11,N12,N13,   &
                          N14,N15,N16,N17,N18,N19,N20,N21,N22,N23,N24,  &
                          N25,N26,N27,N28,N29,N30,N31,N32,N33,N34,N35,  &
                          N36,N37,N38,N39,N40,N41,N42,N43,N44,N45,N46,  &
                          N47,N48,N49,N50,N51,NUSER
!
      DOUBLE PRECISION,DIMENSION(3)::DIPN
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::ADIPx,ADIPy,ADIPz
      DOUBLE PRECISION,DIMENSION(NUSER)::USER
!-----------------------------------------------------------------------
      CALL XtoX0(USER(N11),DIPN ,  3)
      CALL XtoX0(USER(N12),ADIPx,NSQ)
      CALL XtoX0(USER(N13),ADIPy,NSQ)
      CALL XtoX0(USER(N14),ADIPz,NSQ)
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE PASSUSERDIP
      
! PASSQUADUSER
      SUBROUTINE PASSQUADUSER(QUADN,AQUADxx,AQUADyy,AQUADzz,AQUADxy,    &
                              AQUADxz,AQUADyz,USER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/PUNTEROSUSER/N1,N2,N3,N4,N5,N6,N7,N8,N9,N10,N11,N12,N13,   &
                          N14,N15,N16,N17,N18,N19,N20,N21,N22,N23,N24,  &
                          N25,N26,N27,N28,N29,N30,N31,N32,N33,N34,N35,  &
                          N36,N37,N38,N39,N40,N41,N42,N43,N44,N45,N46,  &
                          N47,N48,N49,N50,N51,NUSER
!      
      DOUBLE PRECISION,DIMENSION(6)::QUADN
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::AQUADxx,AQUADyy,AQUADzz
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::AQUADxy,AQUADxz,AQUADyz
      DOUBLE PRECISION,DIMENSION(NUSER)::USER
!-----------------------------------------------------------------------
      CALL XtoX0(QUADN  ,USER(N18),  6)
      CALL XtoX0(AQUADxx,USER(N19),NSQ)
      CALL XtoX0(AQUADyy,USER(N20),NSQ)
      CALL XtoX0(AQUADzz,USER(N21),NSQ)
      CALL XtoX0(AQUADxy,USER(N22),NSQ)
      CALL XtoX0(AQUADxz,USER(N23),NSQ)
      CALL XtoX0(AQUADyz,USER(N24),NSQ)
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE PASSQUADUSER
      
! PASSUSERQUAD
      SUBROUTINE PASSUSERQUAD(QUADN,AQUADxx,AQUADyy,AQUADzz,            &
                              AQUADxy,AQUADxz,AQUADyz,USER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/PUNTEROSUSER/N1,N2,N3,N4,N5,N6,N7,N8,N9,N10,N11,N12,N13,   &
                          N14,N15,N16,N17,N18,N19,N20,N21,N22,N23,N24,  &
                          N25,N26,N27,N28,N29,N30,N31,N32,N33,N34,N35,  &
                          N36,N37,N38,N39,N40,N41,N42,N43,N44,N45,N46,  &
                          N47,N48,N49,N50,N51,NUSER
!
      DOUBLE PRECISION,DIMENSION(6)::QUADN
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::AQUADxx,AQUADyy,AQUADzz
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::AQUADxy,AQUADxz,AQUADyz
      DOUBLE PRECISION,DIMENSION(NUSER)::USER
!-----------------------------------------------------------------------
      CALL XtoX0(USER(N18),QUADN  ,  6)
      CALL XtoX0(USER(N19),AQUADxx,NSQ)
      CALL XtoX0(USER(N20),AQUADyy,NSQ)
      CALL XtoX0(USER(N21),AQUADzz,NSQ)
      CALL XtoX0(USER(N22),AQUADxy,NSQ)
      CALL XtoX0(USER(N23),AQUADxz,NSQ)
      CALL XtoX0(USER(N24),AQUADyz,NSQ)
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE PASSUSERQUAD
      
! PASSOCTUSER
      SUBROUTINE PASSOCTUSER(OCTUN,AOCTxxx,AOCTyyy,AOCTzzz,AOCTxxy,     &
                             AOCTxxz,AOCTxyy,AOCTyyz,AOCTxzz,AOCTyzz,   &
                             AOCTxyz,USER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/PUNTEROSUSER/N1,N2,N3,N4,N5,N6,N7,N8,N9,N10,N11,N12,N13,   &
                          N14,N15,N16,N17,N18,N19,N20,N21,N22,N23,N24,  &
                          N25,N26,N27,N28,N29,N30,N31,N32,N33,N34,N35,  &
                          N36,N37,N38,N39,N40,N41,N42,N43,N44,N45,N46,  &
                          N47,N48,N49,N50,N51,NUSER
!
      DOUBLE PRECISION,DIMENSION(10)::OCTUN
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::AOCTxxx,AOCTyyy,AOCTzzz,AOCTxxy,AOCTxxz
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::AOCTxyy,AOCTyyz,AOCTxzz,AOCTyzz,AOCTxyz
      DOUBLE PRECISION,DIMENSION(NUSER)::USER
!-----------------------------------------------------------------------
      CALL XtoX0(OCTUN  ,USER(N31), 10)
      CALL XtoX0(AOCTxxx,USER(N32),NSQ)
      CALL XtoX0(AOCTyyy,USER(N33),NSQ)
      CALL XtoX0(AOCTzzz,USER(N34),NSQ)
      CALL XtoX0(AOCTxxy,USER(N35),NSQ)
      CALL XtoX0(AOCTxxz,USER(N36),NSQ)
      CALL XtoX0(AOCTxyy,USER(N37),NSQ)
      CALL XtoX0(AOCTyyz,USER(N38),NSQ)
      CALL XtoX0(AOCTxzz,USER(N39),NSQ)
      CALL XtoX0(AOCTyzz,USER(N40),NSQ)
      CALL XtoX0(AOCTxyz,USER(N41),NSQ)
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE PASSOCTUSER
      
! PASSUSEROCT
      SUBROUTINE PASSUSEROCT(OCTUN,AOCTxxx,AOCTyyy,AOCTzzz,AOCTxxy,     &
                             AOCTxxz,AOCTxyy,AOCTyyz,AOCTxzz,AOCTyzz,   &
                             AOCTxyz,USER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/PUNTEROSUSER/N1,N2,N3,N4,N5,N6,N7,N8,N9,N10,N11,N12,N13,   &
                          N14,N15,N16,N17,N18,N19,N20,N21,N22,N23,N24,  &
                          N25,N26,N27,N28,N29,N30,N31,N32,N33,N34,N35,  &
                          N36,N37,N38,N39,N40,N41,N42,N43,N44,N45,N46,  &
                          N47,N48,N49,N50,N51,NUSER
!
      DOUBLE PRECISION,DIMENSION(10)::OCTUN
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::AOCTxxx,AOCTyyy,AOCTzzz,AOCTxxy,AOCTxxz
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::AOCTxyy,AOCTyyz,AOCTxzz,AOCTyzz,AOCTxyz
      DOUBLE PRECISION,DIMENSION(NUSER)::USER
!-----------------------------------------------------------------------
      CALL XtoX0(USER(N31),OCTUN  , 10)
      CALL XtoX0(USER(N32),AOCTxxx,NSQ)
      CALL XtoX0(USER(N33),AOCTyyy,NSQ)
      CALL XtoX0(USER(N34),AOCTzzz,NSQ)
      CALL XtoX0(USER(N35),AOCTxxy,NSQ)
      CALL XtoX0(USER(N36),AOCTxxz,NSQ)
      CALL XtoX0(USER(N37),AOCTxyy,NSQ)
      CALL XtoX0(USER(N38),AOCTyyz,NSQ)
      CALL XtoX0(USER(N39),AOCTxzz,NSQ)
      CALL XtoX0(USER(N40),AOCTyzz,NSQ)
      CALL XtoX0(USER(N41),AOCTxyz,NSQ)
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE PASSUSEROCT
      
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
!                                                                      !
!                     Thermodynamic Properties                         !
!                                                                      !
!      2018 Thermodynamic Properties implemented by Xabier Lopez       ! 
!                                                                      !
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
      
! THERMO    
      SUBROUTINE THERMO(NC1,NROTRA,FREQ,Cxyz,ZMASS,ISIGMA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/EHFEN/EHF,EN
      COMMON/ENERGIAS/EELEC,EELEC_OLD,DIF_EELEC,EELEC_MIN
!
!     THE FORMULAE USED HERE ARE BASED ON
!     OF 'MOLECULAR THERMODYNAMICS' BY McQuarrie
!
!     INPUT: NC1: number of coordinates 3*NATOMS 
!            FREQ(C1): Vibrational Frequencies
!            NROTRA: number of rotations and translations
!            Cxyz: Atom coordinates
!            ZMASS: Mass of the atoms
!            ISIGMA: Rotational symm number
!     OUTPUT: Thermochemistry is written in CGO file (11)
!
      INTEGER, INTENT(IN) :: NC1,ISIGMA
      INTEGER, INTENT(IN) :: NROTRA
      DOUBLE PRECISION,DIMENSION(NC1),INTENT(IN) :: FREQ, Cxyz
      DOUBLE PRECISION,DIMENSION(NATOMS),INTENT(IN) :: ZMASS
      DOUBLE PRECISION, DIMENSION(NC1) :: THETA,ZPVE,QVIB,SVIB,EVIB
      DOUBLE PRECISION, DIMENSION(3) :: VMOI,CMASS
      DOUBLE PRECISION, DIMENSION(NATOMS)   :: ZMS(NATOMS)
      DOUBLE PRECISION, DIMENSION(3,NATOMS) :: CC,COM
      
      LOGICAL LINEAR
      DOUBLE PRECISION KT,JPCAL
!
      PARAMETER (PLANCK=6.62606957D-34, BOLTZ=1.38064852D-23,           &
                 PASCAL=1.01325D+05,                                    &
                 AVOGAD=6.022140857D+23, TOANGS=0.52917724924D+00,      &
                 JPCAL=4.184D+00, CLIGHT=2.99792458D+08,                &
                 TOCM=2.194746D+05, TOKCAL=627.51D+00,                  &
                 TWO=2.0D+00, ONEPT5=1.5D+00, TWOPT5=2.5D+00,           &
                 ONE=1.0D+00, ZERO=0.0D+00,   TM3=1.0D-03,              &
                 EIGHT=8.0D+00, RMASSU=1.660539040D-27)
      
!     Definition of some constants and parameters 
      TEMP = 298.0 
      SCLFAC = ONE ! Posibility to scale the frequencies 
      PI = ACOS(-ONE)
      P = PASCAL
      R = BOLTZ * AVOGAD
      RKCAL = R/(JPCAL*1.0D03)
      KT = BOLTZ * TEMP
      RT = R * TEMP
      RTKCAL=RT/(JPCAL*1.0D03)
      BOHR = 5.29177249D-11
      CATOM = 12.011D+00
!
!     LOAD -CC- AND -ZMS- WITH COORDINATES AND MASSES ---
!     COMPUTING ALSO TOTAL MASS OF THE SYSTEM ZMASST
!
      DO I=1,NATOMS
      CC(1,I) = Cxyz(3*(I-1)+1)
      CC(2,I) = Cxyz(3*(I-1)+2)
      CC(3,I) = Cxyz(3*(I-1)+3)
      ZMS(I) = ZMASS(I)
      END DO
      CALL CENMASd(NATOMS,CC,COM,ZMASST,CMASS,ZMASS)
      
!     ----- ELECTRONIC -----
!     SINCE THE EXCITATION TO EXCITED STATES IS UNKNOWN, AND WE ALSO
!     KNOW NOTHING ABOUT SPIN-ORBIT SPLITTINGS, WE ASSUME THAT ONLY
!     ONE ELECTRONIC STATE CONTRIBUTES TO THE PARTITION FUNCTION.
!     IN ADDITION, SINCE WE DON'T KNOW WHAT THE SPATIAL DEGENERACY OF
!     THIS STATE IS, WE TAKE THE ELECTRONIC DEGENERATCY TO BE JUST
!     THE SPIN MULTIPLICITY.
!
      QELEC = MUL
      QELOG = LOG(QELEC)
      HELEC = EELEC ! EELEC is taken from the NOF calculation 
!
      SELEC = RKCAL*QELOG
      GELEC = HELEC - TEMP * SELEC /TOKCAL ! in a.u.
!
!     ----- TRANSLATION -----
!
!     qtrans = ((2*pi*M*kBT/h**2)**3/2)*V 
!            = (2*pi*M*kBT/h**2)**3/2*kbT/P
!            = ((2Pi)**3/2*(kBT)**5/2*M**3/2)/P
!
      FTRAN=(KT**2.5D+00/AVOGAD**1.5D+00)/PLANCK**3.0D+00 
      FTRAN=(1.0D+00/P)*FTRAN*(2.0D+00*PI)**1.5D+00 ! All in S.I. units
      QTRAN=FTRAN*(TM3*ZMASST)**1.5D+00
      
      WEIGHT = TM3 * ZMASST / AVOGAD
      TTEMP = TWO*PI*(WEIGHT/PLANCK)*(KT/PLANCK)
      QTRAN = TTEMP**ONEPT5 * KT/P
      
      QTLOG =DLOG(QTRAN)
      ETRAN = ONEPT5 * RTKCAL
      HTRAN = TWOPT5 * RTKCAL
      STRAN = RKCAL * (QTLOG + TWOPT5)
      GTRAN = HTRAN - TEMP * STRAN
!
!     ----- ROTATION -----
!           
!
!     ----- SET UP FOR POSSIBLE MONATOMIC -----
!
      IF(NATOMS.EQ.1) THEN
      QROT = ONE
      QRLOG = ZERO
      EROT = ZERO
      HROT = ZERO
      CVROT = ZERO
      CPROT = ZERO
      SROT = ZERO
      GROT = ZERO
      QVIB = ONE
      QVLOG = ZERO
      EVIB = ZERO
      HVIB = ZERO
      CVVIB = ZERO
      CPVIB = ZERO
      SVIB = ZERO
      GVIB = ZERO
      GO TO 700
      END IF
!
!     THE ROTATIONAL SYMMETRY NUMBER SHOULD BE SET CORRECTLY
!     C1    1
!     Cs    1  
!     C2    2 
!     C2v   2 
!     C3v   3 
!     C∞v   1 
!     D2h   4  
!     D3h   6 
!     D5h  10 
!     D∞h   2  
!     D3d   6 
!     Td   12 
!     Oh   24
!
!     At this time this is done through the input section since we 
!     do not consider symmetry in the program, there is the exception of
!     linear molecules, where we can easily detect whether there is a
!     center of symmetry (sigma=2) or not (sigma=1)
!
      SIGMA=DFLOAT(ISIGMA)
      
      CALL INRTIA(CC,COM,ZMS,VMOI,NATOMS)
      LINEAR = VMOI(1).LT.TM3
      IF(LINEAR) THEN
      SIGMA=ONE
      DO K=2,NATOMS
      IF(ABS(COM(1,1)+COM(1,K)).LT.TM3  .AND.                           &
         ABS(COM(2,1)+COM(2,K)).LT.TM3  .AND.                           &
         ABS(COM(3,1)+COM(3,K)).LT.TM3  .AND.                           &
         ZMS(1).EQ.ZMS(K)) SIGMA=TWO
      END DO
      END IF
      FACT1 = (CATOM*BOHR*BOHR)/(12.0D+03*AVOGAD)
      FACT2 = PLANCK/(8.0D+09*PI*PI)
      ACONST = ZERO
      BCONST = ZERO
      CCONST = ZERO
      IF(VMOI(1).GT.0.001D+00) ACONST = FACT2/(FACT1*VMOI(1))
      IF(VMOI(2).GT.0.001D+00) BCONST = FACT2/(FACT1*VMOI(2))
      IF(VMOI(3).GT.0.001D+00) CCONST = FACT2/(FACT1*VMOI(3))
      WRITE(11,9011) ACONST, BCONST, CCONST ! These are in G(HZ)
      
!     ROTATIONAL TEMPERATURES in K-1
      
      FTROT=4.79924466221135D-02       ! (h/kb)*10**9
      TROTA=FTROT*ACONST
      TROTB=FTROT*BCONST
      TROTC=FTROT*CCONST
      
      WRITE(11,9112) 'TROT(K): ', TROTA, TROTB, TROTC
      
      IF(LINEAR) THEN
      QROT = TEMP/(TROTC*SIGMA)
      QRLOG = LOG(QROT)
      EROT = RTKCAL
      SROT = RKCAL * (QRLOG + ONE)
      ELSE
      QROT =(SQRT(PI)/SIGMA)*((TEMP**1.5D0)/(DSQRT(TROTA*TROTB*TROTC)))
      QRLOG = LOG(QROT)
      EROT  = ONEPT5 * RTKCAL
      SROT = RKCAL * (QRLOG + ONEPT5)
      END IF
      HROT = EROT
      GROT = HROT - TEMP * SROT
!
!     ----- VIBRATION -----
!
      FACTHE=1.439773398663405D0 ! PLANCK*CLIGHT/BOLTZ*100 
                                 ! we apply to freq(cm-1)
                                 ! to get the vib. temp. (K)
      DO I=1,NC1
      THETA(i)=FACTHE*FREQ(i)*SCLFAC
      END DO
      QVIBT = ZERO
      EVIBT = ZERO
      SVIBT = ZERO
      ZPVET = ZERO
      DO I=NROTRA+1,NC1
      IF(THETA(I).EQ.ZERO) THEN
      WRITE(11,9200) I
      CALL ABRT
      END IF
!        Prefactor is 0.5hc*Nav/JPCAL         
      ZPVE(i) = 1.430561440340379D0*FREQ(i)*1.0D-03 
      QVIB(i)=DEXP(-0.5D0*THETA(i)/TEMP)/                               &
       (1.0D0 - DEXP(-1.0D0*THETA(i)/TEMP))                       
      THETAT=THETA(i)/TEMP                                            
      SVIB(i)= RKCAL*((THETAT)/(DEXP(THETAT)-1.0D0) -                   &
        DLOG(1-DEXP(-1.0D0*THETAT)))
      EVIB(i)=RKCAL*THETA(i)*(0.5D0+1.0D0/(DEXP(THETAT)-1.0D0))
      
      QVIBT = QVIBT + QVIB(I)
      SVIBT = SVIBT + SVIB(I)
      EVIBT = EVIBT + EVIB(I)
      ZPVET = ZPVET + ZPVE(I)
      END DO
      HVIBT=EVIBT
      GVIBT = HVIBT - TEMP * SVIBT
      
!
!     ----- TOTALS -----
!     Summing up Elec + Trans + Rot+ Vib
!
      700 CONTINUE
      QTOT = QELEC * QTRAN * QROT  * QVIBT
      WRITE(11,9071)
      WRITE(11,*) '####################################'
      WRITE(11,*) '#     THERMOCHEMISTRY ANALYSIS     #'
      WRITE(11,*) '####################################'
      write(11,9069) TEMP,SCLFAC
      WRITE(11,9071)
      WRITE(11,*)    '  Q   '
      WRITE(11,9080) 'ELEC. ',QELEC
      WRITE(11,9080) 'TRANS.',QTRAN
      WRITE(11,9081) 'ROT.  ',QROT,' with sigma ',SIGMA
      WRITE(11,9080) 'VIB.  ',QVIBT
      WRITE(11,9080) 'TOT.  ',QTOT
      
! 
!     thermal energies in kcal/mol
!
      WRITE(11,9090)
      WRITE(11,9100)
      WRITE(11,9110) 'TRANS.      ',ETRAN,HTRAN,GTRAN,STRAN*1.0D+03
      WRITE(11,9110) 'ROT.        ',EROT,HROT,GROT,SROT*1.0D+03
      WRITE(11,9110) 'VIB(Total)  ',EVIBT,HVIBT,GVIBT,SVIBT*1.0D+03
      WRITE(11,9114) 
      DO I=NROTRA+1,NC1
      WRITE(11,9113)I,FREQ(i),EVIB(i),SVIB(I)*1.0D+03 
      END DO   
!
!     Transforming to a.u. 
!
      ETRAN=ETRAN/TOKCAL
      HTRAN=HTRAN/TOKCAL
      GTRAN=GTRAN/TOKCAL
      
      EROT =EROT/TOKCAL
      HROT =HROT/TOKCAL
      GROT =GROT/TOKCAL
      
      ZPVET=ZPVET/TOKCAL
      EVIBT=EVIBT/TOKCAL
      HVIBT=HVIBT/TOKCAL
      GVIBT=GVIBT/TOKCAL
      
      CORRE = ETRAN + EROT + EVIBT
      CORRH = HTRAN + HROT + HVIBT
      CORRG = GTRAN + GROT + GVIBT
!
!  --- Final Results written to unit 100 ----
!
      WRITE(11,9071)
      WRITE(11,9101) 
      WRITE(11,9111) 'Zero-point correction                   =',ZPVET
      WRITE(11,9111) 'Thermal correction to Energy            =',CORRE
      WRITE(11,9111) 'Thermal correction to Enthalpy          =',CORRH
      WRITE(11,9111) 'Thermal correction to Gibbs Free Energy =',CORRG
      
      WRITE(11,9071)
      
      WRITE(11,9111) 'Ee + EN + ZPVE                   = ',             &
                EELEC+EN+ZPVET                                    
      WRITE(11,9111) 'Ee + EN + thermal Energies       = ',             &
                EELEC+EN+CORRE                                    
      WRITE(11,9111) 'Ee + EN + thermal Enthalpies     = ',             &
                EELEC+EN+CORRH                                    
      WRITE(11,9111) 'Ee + EN + thermal Free Energies  = ',             &
                EELEC+EN+CORRG
      WRITE(11,9101) 
      
      RETURN
!-----------------------------------------------------------------------
      9011 FORMAT(1X,'THE ROTATIONAL CONSTANTS ARE (IN GHZ)',/          &
      1X,3F12.5)                                                  
      9069 FORMAT(/1X,'TEMPERATURE ',F10.4, 4X,'SCALE FACTOR ',F10.4)
      9071 FORMAT(/1X)
      9080 FORMAT(1X,A6,1P,E15.5)
      9081 FORMAT(1X,A6,1P,E15.5,A15,E10.2)
      9090 FORMAT(/20X,'E',12X,'H',12X,'G',12X,'S')
      9100 FORMAT(17X,'KCAL/MOL',6X,'KCAL/MOL',6X,'KCAL/MOL',5X,'CAL/MOL-K')
      9101 FORMAT(1X,/,71(1H-),/)
      9110 FORMAT(1X,A10,3F14.6,3X,E14.6)
      9113 FORMAT(5X,I4,4X,F10.2,1X,F15.6,5X,E14.6)
      9114 FORMAT(/8X,'#',7X,'Freq(cm-1)   EVIB(KCAL/MOL)  SVIB(CAL/MOL-K)')
      9111 FORMAT(1X,A40,F14.6," (a.u.) ")
      9112 FORMAT(1X,A10,3F14.6)
      9200 FORMAT(1X,'FREQUENCY NUMBER',I4,' IS UNEXPECTEDLY ZERO!')
      END

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
