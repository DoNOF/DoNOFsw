!======================================================================!
!                                                                      !
!                    M C P T   S U B R O U T I N E S                   !
!                                                                      !
!                 ( J. Chem. Phys. 139, 064111, 2013 )                 !
!                                                                      !
!======================================================================!
!======================================================================!

! SC2MCPThf
      SUBROUTINE SC2MCPThf(RO,COEF,AHCORE,IERI,ERI,QK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/EHFEN/EHF,EN
      COMMON/INPNOF_PNOF/IPNOF,NTWOPAR
      COMMON/INPNOF_ORBSPACE0/NO1,NDOC,NCO,NCWO,NVIR,NAC,NO0
      COMMON/INPFILE_NO1PT2/NO1PT2,NEX
      COMMON/INPFILE_NIJKL/NINTMX,NIJKL,NINTCR,NSTORE
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      COMMON/NumLinIndOrb/NQMT
!
      INTEGER(8),DIMENSION(NIJKL)::IERI
      DOUBLE PRECISION,DIMENSION(NIJKL)::ERI
      DOUBLE PRECISION,DIMENSION(NBF5)::RO
      DOUBLE PRECISION,DIMENSION(NBFT5)::QK
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::COEF,AHCORE
!
      INTEGER,ALLOCATABLE,DIMENSION(:):: IPOS,IEX
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::  OCC,OCCn,QKv,EIG
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:):: VEC,DMhf,BBn
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:):: AUX,FOCK,TVEC,FOCKm
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:,:):: ERImol
!-----------------------------------------------------------------------
!     NA:  Number of HF occupied MOs (OCC=1)
!     NVI: Number of HF virtual  MOs (OCC=0)
!
!     NO1:  Number of inactive doubly occupied orbitals (OCC=1)
!     NDOC: Number of strongly occupied MOs
!     NCWO: Number of coupled weakly occupied MOs per strongly occupied
!     NCWO*NDOC: Active orbitals in the virtual subspace
!     NO0: Empty orbitals  (OCC=0)
!     NAC:  Dimension of the active natural orbital subspace
!
!           NA      |       NVI           = NBF
!       NO1 + NDOC  |  NCWO*NDOC + NO0    = NBF
!           |      NAC           |
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!                  EHFL: HF Energy with localized orbitals
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     DMhf: HF Density Matrix     ( DMhf -> FOCK = H + 2J-K )
!     FOCK: Fock Matrix in the atomic basis set
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE (DMhf(NBF,NBF),FOCK(NBF,NBF))
      CALL DENMATHFr(DMhf,COEF)
      CALL FORM2JK(FOCK,DMhf,IERI,ERI)
      FOCK = AHCORE + FOCK
      EHFL = TRACE(DMhf,AHCORE,NBF) + TRACE(DMhf,FOCK,NBF)
      DEALLOCATE(DMhf)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!     Substracting core orbitals (NO1)
!     If using ECPotentials (IECP/=0): NO1=0
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      NOC = NA - NO1                              ! NOC = NDOC
      NVI = NQMT - NA
      NORB = NQMT - NO1
      WRITE(6,1)NO1,NOC,NVI,NQMT
      NOCNC = NOC + NCWO*NOC
      NOCNCT = NOCNC*(NOCNC+1)/2
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Check the value of NCWO if NOCNC = NOC + NCWO*NOC > NORB
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(NOCNC>NORB)THEN
       NCWO1 = NORB/NOC - 1 
       IF(NCWO>NCWO1)THEN
        WRITE(6,2)NCWO1
        STOP
       ENDIF
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     OCC,VEC: Occupations and Molecular Orbitals without Core
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(OCC(NOCNC),VEC(NBF,NORB),QKv(NOCNCT))
      DO J=1,NOCNC
       OCC(J) = RO(J+NO1)
       DO I=1,J
        IJ = I + J*(J-1)/2
        IJold = I+NO1 + (J+NO1)*(J+NO1-1)/2
        QKv(IJ) = QK(IJold)
       ENDDO
      ENDDO
      DO J=1,NORB
       DO I=1,NBF
        VEC(I,J) = COEF(I,J+NO1)
       ENDDO
      ENDDO
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!     FOCKm: Fock Matrix in the molecular basis set
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
      ALLOCATE (FOCKm(NORB,NORB),AUX(NBF,NORB),TVEC(NORB,NBF))
      AUX   = MATMUL(FOCK,VEC)
      TVEC  = TRANSPOSE(VEC)     
      FOCKm = MATMUL(TVEC,AUX)
      DEALLOCATE(FOCK,AUX,TVEC)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     EIG = One-particle energies [E(i)=Ct(i)*FOCK*C(i)=FOCKm(i,i)]
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE (EIG(NORB))
      DO i=1,NORB
       EIG(i) = FOCKm(i,i)
      ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Print one-particle energies (EIG) if IPRINTEIG = 1
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IPRINTEIG=0
      IF(IPRINTEIG==1)THEN
       write(6,'(/,A22,/)')'One-particle energies:'
       do i=1,NORB
        write(6,'(4x,I5,F20.10)')NO1+i,EIG(i)
       enddo
      ENDIF
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!     Integral transformation ( Form <ab|rs> )
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!     Allocate <ab|rs>
      ALLOCATE(ERImol(NBF,NBF,NBF,NBF))
!     Form <al|mn> from <kl|mn>
      CALL ERIC1c(ERImol,IERI,ERI,VEC,NORB)
!     Form <ab|ms>  [ ERImol(ia,l,m,k) -> ERImol(ia,ib,is,k) ]
      CALL ERIC23c(ERImol,VEC,NORB)
!     Form <ab|rs>  [ ERImol(ia,ib,is,k) -> ERImol(ia,ib,is,ir) ]
      CALL ERIC4c(ERImol,VEC,NORB)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!                  PNOF5(Nc)-SC2-MCPT and PNOF5(Nc)-PT2
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      IF(IPNOF==5)THEN
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      BBn(i,iw) = BETA(i,iw)/SQRT(OCCi) = SQRT(OCCk)/SQRT(OCCi)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       ALLOCATE(BBn(NOC,NCWO))
       do i=1,NOC
        do iw=1,NCWO
         !k = noc+ncwo*(noc-i)+iw         !old-sort           ! k = (i,iw)
         k = noc+noc*(iw+1)-i+1           !new-sort
         BBn(i,iw) = DSQRT(OCC(k)/OCC(i))
        enddo
       enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      E0 = < 0'| H |0 > = EHFL + Ecorr   ( E1 = 0 )
!      Ecorr: Correlation Energy with respect to EHFL
!      E2x: Ec(2) correlation energy
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Ecorr = 0.0d0
       E2x = 0.0d0
       DO i=1,NOC
        do iw=1,NCWO
         !k = NOC+NCWO*(NOC-i)+iw         !old-sort
         k = noc+noc*(iw+1)-i+1           !new-sort
         ik = i + k*(k-1)/2
         BBQK = BBn(i,iw)*QKv(ik)
         Eki = 2.0d0*(EIG(k)-EIG(i))
         Ecorr = Ecorr - BBQK
         E2x = E2x - BBQK/Eki
        enddo
       ENDDO
       E0 = EHFL + Ecorr
       WRITE(6,4)EHFL+EN,Ecorr,E0+EN
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!      SC2-MCPT: E2 = Sum_k=1 [  < 0'| H |k >< k'| H |0 > / (E0-Ek) ]
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
       E2 = Ecorr*E2x                                                   &
          + E2HFs(NOC,NORB,EIG,FOCKm)                                   &
          + E2HFd(NOC,NORB,NBF,EIG,ERImol)                              &
          + E2F(NOC,NCWO,NORB,EIG,FOCKm,BBn)                            &
          + E2FERI(NOC,NCWO,NORB,NBF,EIG,FOCKm,ERImol,BBn)              &
          + E2ERIERI(NOC,NCWO,NORB,NBF,EIG,ERImol,BBn)                  &
          + E2dHF(NOC,NCWO,NORB,NBF,EIG,ERImol,BBn)                      
       IF(NOC/=1)THEN                                                    
        E2 = E2 - Ecorr*E2x/2.0d0                                       &
                + E2Ql(NOC,NCWO,NORB,NOCNCT,EIG,QKv,BBn)                &
                + E2Qd(NOC,NCWO,NORB,NBF,EIG,ERImol,BBn)
       ENDIF
       IF(NEX==0)WRITE(6,5)E2,E0+E2+EN
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!      PNOF5(Nc)-PT2: exclude excitations from the same spatial orbital
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
       E2ex = E2HFs(NOC,NORB,EIG,FOCKm)                                 &
            + E2HFdex(NOC,NORB,NBF,EIG,ERImol)                          &
            + E2F(NOC,NCWO,NORB,EIG,FOCKm,BBn)                          &
            + E2FERIex(NOC,NCWO,NORB,NBF,EIG,FOCKm,ERImol,BBn)          &
            + E2ERIERIex(NOC,NCWO,NORB,NBF,EIG,ERImol,BBn)              &
            + E2dHFex(NOC,NCWO,NORB,NBF,EIG,ERImol,BBn)
       IF(NOC/=1)THEN
        E2ex = E2ex + E2Qd(NOC,NCWO,NORB,NBF,EIG,ERImol,BBn) 
       ENDIF
       IF(NEX==0)WRITE(6,6)E2ex,E0+E2ex+EN
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!      PNOF5(Nc)-PT2-NEX: Remove subspaces associated to NEX->E2ex-E2exi
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
       IF(NEX>0)THEN
        IF(NEX>NOC)THEN
         WRITE(6,'(/1X,A32,I5)')'STOP: NEX must be lesser than ',NOC+1
         STOP
        ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       IEX: Position of the 'NEX' minimum strongly occupation numbers
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ALLOCATE(IEX(NEX),IPOS(1),OCCn(NOC))
        IEX = -1
        WRITE(6,3)NEX
        OCCn = OCC(1:NOC)
        do i=1,nex
         IPOS = MINLOC(OCCn)
         IEX(i) = IPOS(1)
         OCCn(IEX(i)) = 1.0d0
        enddo
        E2exi = E2HFs_nex(NOC,NCWO,NORB,EIG,FOCKm,IEX,NEX)              &
              + E2HFdex_nex(NOC,NCWO,NORB,NBF,EIG,ERImol,IEX,NEX)       &
              + E2F_nex(NOC,NCWO,NORB,EIG,FOCKm,BBn,IEX,NEX)            &
              + E2FERIex_nex(NOC,NCWO,NORB,NBF,EIG,FOCKm,ERImol,BBn,    &
                IEX,NEX) + E2ERIERIex_nex                               &
                (NOC,NCWO,NORB,NBF,EIG,ERImol,BBn,IEX,NEX)              &
              + E2dHFex_nex(NOC,NCWO,NORB,NBF,EIG,ERImol,BBn,IEX,NEX)    
        IF(NOC/=1)THEN                                                   
         E2exi = E2exi                                                  &
               + E2Qd_nex(NOC,NCWO,NORB,NBF,EIG,ERImol,BBn,IEX,NEX) 
        ENDIF
        WRITE(6,7)E2ex+E2exi,E0+E2ex+E2exi+EN
        DEALLOCATE(IEX,IPOS,OCCn)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       DEALLOCATE(BBn)
      ENDIF
!-----------------------------------------------------------------------
      DEALLOCATE(OCC,VEC,QKv,EIG,FOCKm,ERImol)
      RETURN
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    1 FORMAT(//2X,'PNOF5-MCPT',/1X,12('='),//,                          &
               1X,'NUMBER OF CORE ORBITALS    (NO1) =',I5,/,            &
               1X,'NUMBER OF OCCUPIED ORBS.   (NOC) =',I5,/,            &
               1X,'NUMBER OF VIRTUAL ORBS.    (NVI) =',I5,/,            &
               1X,'NUMBER OF LIN. IND. ORBS. (NQMT) =',I5)               
    2 FORMAT(/1X,'NCWO is too large, reduce the value at least to',I5)   
    3 FORMAT(/1X,'NUMBER OF EXCLUDED COUPLES (NEX) =',I5,/)              
    4 FORMAT(/11X,'Ehf(l) =',F20.10,/,                                  &
              11X,' Ecorr =',F20.10,/,                                  &
              11X,'  E(0) =',F20.10)                                     
    5 FORMAT(/6X,'       E(2) =',F20.10,/,                              &
              6X,'E(SC2-MCPT) =',F20.10)                                 
    6 FORMAT(/5X,'        E(2) =',F20.10,/,                             &
              5X,'E(PNOF5-PT2) =',F20.10)                                
    7 FORMAT(1X,'            E(2) =',F20.10,/,                          &
             1X,'E(PNOF5-PT2-NEX) =',F20.10)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
!                                                                      !
!                Integral transformation ( Form <ab|rs> )              !
!                                                                      !
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !

! ERIC1c
      SUBROUTINE ERIC1c(ERImol,IERI,ERI,VEC,NORB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPFILE_NIJKL/NINTMX,NIJKL,NINTCR,NSTORE
!
      INTEGER(8),DIMENSION(NSTORE)::IERI
      DOUBLE PRECISION,DIMENSION(NSTORE)::ERI
      DOUBLE PRECISION,DIMENSION(NBF,NBF,NBF,NBF)::ERImol
      DOUBLE PRECISION,DIMENSION(NBF,NORB)::VEC
!-----------------------------------------------------------------------
      ERImol = 0.0d0
      DO M=1,NINTCR
       XINT1 = ERI(M)
       XINT2 = XINT1
       LABEL = IERI(M)
       CALL LABELIJKL(LABEL,I,J,K,L)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       IF(I==J)XINT1 = XINT1 + XINT1
       IF(K==L)XINT2 = XINT2 + XINT2
       do ia=1,norb
        ERImol(ia,I,J,K) = ERImol(ia,I,J,K) + XINT1*VEC(L,ia)
        ERImol(ia,J,I,K) = ERImol(ia,I,J,K)
        ERImol(ia,I,J,L) = ERImol(ia,I,J,L) + XINT1*VEC(K,ia)
        ERImol(ia,J,I,L) = ERImol(ia,I,J,L)
        ERImol(ia,K,L,I) = ERImol(ia,K,L,I) + XINT2*VEC(J,ia)
        ERImol(ia,L,K,I) = ERImol(ia,K,L,I)
        ERImol(ia,K,L,J) = ERImol(ia,K,L,J) + XINT2*VEC(I,ia)
        ERImol(ia,L,K,J) = ERImol(ia,K,L,J)
       enddo
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! ERIC23c
      SUBROUTINE ERIC23c(ERImol,VEC,NORB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      DOUBLE PRECISION,DIMENSION(NBF,NBF,NBF,NBF)::ERImol
      DOUBLE PRECISION,DIMENSION(NBF,NORB)::VEC
      ALLOCATABLE::AUX1(:,:),AUX2(:,:)
      ALLOCATE(AUX1(NBF,NBF),AUX2(NBF,NORB))
!-----------------------------------------------------------------------
      do k=1,nbf
       do ia=1,norb
! - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       ERImol(ia,i,j,k) -> AUX1(i,j) for each (ia,k)
! - - - - - - - - - - - - - - - - - - - - - - - - - - -
        do i=1,nbf
         do j=1,nbf
          AUX1(i,j) = ERImol(ia,i,j,k)
         enddo
        enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       AUX2=AUX1*C
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        DO I=1,NBF
         do is=1,norb
          AUX2(I,is) = 0.0d0
          DO L=1,NBF
           AUX2(I,is) = AUX2(I,is) + AUX1(I,L)*VEC(L,is)
          ENDDO
         enddo
        ENDDO
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       AUX1=VECt*AUX2
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        CALL CeqAtB(AUX1,VEC,NBF,NORB,AUX2,NORB)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       AUX1 -> ERImol(ia,ib,is,k)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        do ib=1,norb
         do is=1,norb
          ERImol(ia,ib,is,k) = AUX1(ib,is)
         enddo
        enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       enddo
      enddo
!-----------------------------------------------------------------------
      DEALLOCATE(AUX1,AUX2)
      RETURN
      END

! ERIC4c
      SUBROUTINE ERIC4c(ERImol,VEC,NORB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      DOUBLE PRECISION,DIMENSION(NBF,NBF,NBF,NBF)::ERImol
      DOUBLE PRECISION,DIMENSION(NBF,NORB)::VEC
      ALLOCATABLE::AUX1(:,:),AUX2(:,:)
      ALLOCATE(AUX1(NBF,NORB),AUX2(NBF,NORB))
!-----------------------------------------------------------------------
      do ib=1,norb
       do is=1,norb
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       ERImol(ia,ib,is,k) -> AUX1(k,ia) for each (ib,is)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        do ia=1,norb
         do k=1,nbf
          AUX1(k,ia) = ERImol(ia,ib,is,k)
         enddo
        enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       AUX2 = VEC*AUX1
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        CALL CeqAtB(AUX2,VEC,NBF,NORB,AUX1,NORB)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       AUX2 -> ERImol(ia,ib,is,ir)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        do ir=1,norb
         do ia=1,norb
          ERImol(ia,ib,is,ir) = AUX2(ir,ia)
         enddo
        enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       enddo
      enddo
!-----------------------------------------------------------------------
      DEALLOCATE(AUX1,AUX2)
      RETURN
      END

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
!                                                                      !
!                     PNOF5(Nc)-SC2-MCPT -> E2                         !
!                                                                      !
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !

! E2HFs
      FUNCTION E2HFs(NOC,NORB,EIG,FOCKm)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DOUBLE PRECISION,DIMENSION(NORB)::EIG
      DOUBLE PRECISION,DIMENSION(NORB,NORB)::FOCKm
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      E2HFs = 0.0d0
      DO i=1,noc
       DO l=noc+1,norb
        Eil = EIG(l) - EIG(i)
        E2HFs = E2HFs - FOCKm(i,l)*FOCKm(l,i)/Eil
       ENDDO
      ENDDO
      E2HFs = 2.0d0*E2HFs
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END

! E2HFd
      FUNCTION E2HFd(NOC,NORB,NBF,EIG,ERImol)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DOUBLE PRECISION,DIMENSION(NORB)::EIG
      DOUBLE PRECISION,DIMENSION(NBF,NBF,NBF,NBF)::ERImol
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      E2HFd = 0.0d0
      DO i=1,noc
       DO j=1,noc
        DO k=noc+1,norb
         DO l=noc+1,norb
          Xijkl = ERImol(i,j,l,k)
          Xijlk = ERImol(i,j,k,l)
          Eklij = EIG(k) + EIG(l) - EIG(i) - EIG(j)
          E2HFd = E2HFd - Xijkl*(2.0*Xijkl-Xijlk)/Eklij
         ENDDO
        ENDDO
       ENDDO
      ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END

! E2F
      FUNCTION E2F(NOC,NCWO,NORB,EIG,FOCKm,BBn)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DOUBLE PRECISION,DIMENSION(NORB)::EIG
      DOUBLE PRECISION,DIMENSION(NOC,NCWO)::BBn
      DOUBLE PRECISION,DIMENSION(NORB,NORB)::FOCKm
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      E2F = 0.0d0
      DO i=1,noc
       do iw=1,ncwo
        !k = noc+ncwo*(noc-i)+iw         !old-sort     ! k = (i,iw)
        k = noc+noc*(iw+1)-i+1           !new-sort
        Eki = EIG(k) - EIG(i)
        E2F = E2F + BBn(i,iw)*FOCKm(i,k)*FOCKm(k,i)/Eki 
       enddo
      ENDDO
      E2F = 2.0d0*E2F
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END

! E2FERI
      FUNCTION E2FERI(NOC,NCWO,NORB,NBF,EIG,FOCKm,ERImol,BBn)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DOUBLE PRECISION,DIMENSION(NORB)::EIG
      DOUBLE PRECISION,DIMENSION(NOC,NCWO)::BBn
      DOUBLE PRECISION,DIMENSION(NORB,NORB)::FOCKm
      DOUBLE PRECISION,DIMENSION(NBF,NBF,NBF,NBF)::ERImol
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      E2FERI = 0.0d0

      DO i=1,noc
       DO l=noc+1,norb
        Eli = EIG(l) - EIG(i)
        do iw=1,ncwo
         !k = noc+ncwo*(noc-i)+iw         !old-sort     ! k = (i,iw)
         k = noc+noc*(iw+1)-i+1           !new-sort
         Eklii = EIG(k) - EIG(i) + Eli
         Xiikl = ERImol(i,i,l,k)
         Xilkk = ERImol(i,l,k,k)
         E2FERI = E2FERI + BBn(i,iw)*                                   &
                ( FOCKm(i,l)*Xilkk/Eli + FOCKm(k,l)*Xiikl/Eklii )
        enddo
       ENDDO
      ENDDO

      DO j=1,noc
       DO i=1,noc
        do iw=1,ncwo
         !k = noc+ncwo*(noc-i)+iw         !old-sort     ! k = (i,iw)
         k = noc+noc*(iw+1)-i+1           !new-sort
         Ekj = EIG(k) - EIG(j)
         Ekkij = 2.0d0*EIG(k) - EIG(i) - EIG(j)
         Xiijk = ERImol(i,i,k,j)
         Xijkk = ERImol(i,j,k,k)
         E2FERI = E2FERI - BBn(i,iw)*                                   &
                ( FOCKm(j,k)*Xiijk/Ekj + FOCKm(i,j)*Xijkk/Ekkij )
        enddo
       ENDDO
      ENDDO

      E2FERI = 2.0d0*E2FERI
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END

! E2ERIERI
      FUNCTION E2ERIERI(NOC,NCWO,NORB,NBF,EIG,ERImol,BBn)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DOUBLE PRECISION,DIMENSION(NORB)::EIG
      DOUBLE PRECISION,DIMENSION(NOC,NCWO)::BBn
      DOUBLE PRECISION,DIMENSION(NBF,NBF,NBF,NBF)::ERImol
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      E2ERIERI = 0.0d0
      DO i=1,noc
       do iw=1,ncwo
        !k = noc+ncwo*(noc-i)+iw         !old-sort     ! k = (i,iw)
        k = noc+noc*(iw+1)-i+1           !new-sort
        DO m=1,noc
         DO j=1,noc
          Ekkmj = 2.0d0*EIG(k) - EIG(m) - EIG(j)
          Xiimj = ERImol(i,i,j,m)
          Xmjkk = ERImol(m,j,k,k)
          E2ERIERI = E2ERIERI + BBn(i,iw)*Xiimj*Xmjkk/Ekkmj
         ENDDO
        ENDDO
        DO l=noc+1,norb
         DO n=noc+1,norb
          Enlii = EIG(n) + EIG(l) - 2.0d0*EIG(i)
          Xiiln = ERImol(i,i,l,n)
          Xlnkk = ERImol(l,n,k,k)
          E2ERIERI = E2ERIERI + BBn(i,iw)*Xiiln*Xlnkk/Enlii
         ENDDO
        ENDDO
       enddo
      ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END

! E2dHF
      FUNCTION E2dHF(NOC,NCWO,NORB,NBF,EIG,ERImol,BBn)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DOUBLE PRECISION,DIMENSION(NORB)::EIG
      DOUBLE PRECISION,DIMENSION(NOC,NCWO)::BBn
      DOUBLE PRECISION,DIMENSION(NBF,NBF,NBF,NBF)::ERImol
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      E2dHF = 0.0d0
      DO i=1,noc
       do iw=1,ncwo
        !k = noc+ncwo*(noc-i)+iw         !old-sort     ! k = (i,iw)
        k = noc+noc*(iw+1)-i+1           !new-sort
        DO j=1,noc
         DO l=noc+1,norb
          Elkij = EIG(l) + EIG(k) - EIG(i) - EIG(j)
          Xijlk = ERImol(i,j,k,l)
          Xiljk = ERImol(i,l,k,j)
          Xijkl = ERImol(i,j,l,k)
          Xilkj = ERImol(i,l,j,k)
          XXijlk = Xijlk - Xijkl
          XXiljk = Xiljk - Xilkj
          E2dHF = E2dHF - (BBn(i,iw)/Elkij) *                           &
                          ( XXijlk*Xilkj + Xijlk*Xiljk + Xijkl*XXiljk )
         ENDDO
        ENDDO
       enddo
      ENDDO
      E2dHF = 2.0d0*E2dHF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END

! E2Ql
      FUNCTION E2Ql(NOC,NCWO,NORB,NOCNCT,EIG,QKv,BBn)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DOUBLE PRECISION,DIMENSION(NORB)::EIG
      DOUBLE PRECISION,DIMENSION(NOCNCT)::QKv
      DOUBLE PRECISION,DIMENSION(NOC,NCWO)::BBn
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      E2Ql = 0.0d0
      DO i=1,noc
       DO j=1,noc
        if(j/=i)then
         do iw=1,ncwo
          !k = noc+ncwo*(noc-i)+iw         !old-sort     ! k = (i,iw)
          k = noc+noc*(iw+1)-i+1           !new-sort
          jk = j + k*(k-1)/2
          do jw=1,ncwo
           !l = noc+ncwo*(noc-j)+jw         !old-sort     ! k = (i,iw)
           l = noc+noc*(jw+1)-j+1           !new-sort
           il = i + l*(l-1)/2
           Eli = EIG(l)-EIG(i)
           BB2 = BBn(i,iw)*BBn(j,jw)
           E2Ql = E2Ql - BB2*QKv(il)*QKv(jk)/Eli
          enddo
         enddo
        endif
       ENDDO
      ENDDO
      E2Ql = E2Ql/4.0d0
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END

! E2Qd
      FUNCTION E2Qd(NOC,NCWO,NORB,NBF,EIG,ERImol,BBn)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DOUBLE PRECISION,DIMENSION(NORB)::EIG
      DOUBLE PRECISION,DIMENSION(NOC,NCWO)::BBn
      DOUBLE PRECISION,DIMENSION(NBF,NBF,NBF,NBF)::ERImol
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      E2Qd = 0.0d0
      DO i=1,noc
       DO j=1,noc
        do iw=1,ncwo
         !k = noc+ncwo*(noc-i)+iw         !old-sort     ! k = (i,iw)
         k = noc+noc*(iw+1)-i+1           !new-sort
         do jw=1,ncwo
          !l = noc+ncwo*(noc-j)+jw         !old-sort     ! k = (i,iw)
          l = noc+noc*(jw+1)-j+1           !new-sort
          BB2 = BBn(i,iw)*BBn(j,jw)
          if(j/=i)then
           Eklii = EIG(k) + EIG(l) - 2.0d0*EIG(i)
           Xiikl = ERImol(i,i,l,k)
           Xjjkl = ERImol(j,j,l,k)
           E2Qd = E2Qd + BB2*Xiikl*Xjjkl/Eklii
          endif
          if(l/=k)then
           Ekkij = 2.0d0*EIG(k) - EIG(i) - EIG(j)
           Xijkk = ERImol(i,j,k,k)
           Xijll = ERImol(i,j,l,l)
           E2Qd = E2Qd + BB2*Xijkk*Xijll/Ekkij
          endif
          Eklij = EIG(k) - EIG(i) + EIG(l) - EIG(j)
          Xijkl = ERImol(i,j,l,k)
          Xijlk = ERImol(i,j,k,l)
          XX = Xijkl*Xijkl - Xijkl*Xijlk + Xijlk*Xijlk
          E2Qd = E2Qd - BB2*XX/Eklij
         enddo
        enddo
       ENDDO
      ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
!                                                                      !
!                      PNOF5(Nc)-PT2 -> E2ex                           !
!                                                                      !
!         exclude excitations from the same spatial orbital            !
!                                                                      !
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !

! E2HFdex
      FUNCTION E2HFdex(NOC,NORB,NBF,EIG,ERImol)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DOUBLE PRECISION,DIMENSION(NORB)::EIG
      DOUBLE PRECISION,DIMENSION(NBF,NBF,NBF,NBF)::ERImol
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     i,j (i/=j) -> k,l
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      E2HFdex = 0.0d0
      DO i=1,noc
       DO j=1,noc
        if(i/=j)then
         DO k=noc+1,norb
          DO l=noc+1,norb
           Xijkl = ERImol(i,j,l,k)
           Xijlk = ERImol(i,j,k,l)
           Eklij = EIG(k) + EIG(l) - EIG(i) - EIG(j)
           E2HFdex = E2HFdex - Xijkl*(2.0*Xijkl-Xijlk)/Eklij
          ENDDO
         ENDDO
        endif
       ENDDO
      ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END

! E2FERIex
      FUNCTION E2FERIex(NOC,NCWO,NORB,NBF,EIG,FOCKm,ERImol,BBn)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DOUBLE PRECISION,DIMENSION(NORB)::EIG
      DOUBLE PRECISION,DIMENSION(NOC,NCWO)::BBn
      DOUBLE PRECISION,DIMENSION(NORB,NORB)::FOCKm
      DOUBLE PRECISION,DIMENSION(NBF,NBF,NBF,NBF)::ERImol
!-----------------------------------------------------------------------
      E2FERIex = 0.0d0
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     i -> l    ( BBn(i,iw)*FOCKm(k,l)*Xiikl/Eklii is excluded )
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO i=1,noc
       DO l=noc+1,norb
        Eli = EIG(l) - EIG(i)
        do iw=1,ncwo
         !k = noc+ncwo*(noc-i)+iw         !old-sort     ! k = (i,iw)
         k = noc+noc*(iw+1)-i+1           !new-sort
         Xilkk = ERImol(i,l,k,k)
         E2FERIex = E2FERIex + BBn(i,iw)*FOCKm(i,l)*Xilkk/Eli
        enddo
       ENDDO
      ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     j -> k and i/=j -> k,k
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO j=1,noc
       DO i=1,noc
        if(i/=j)then       
         do iw=1,ncwo
          !k = noc+ncwo*(noc-i)+iw         !old-sort     ! k = (i,iw)
          k = noc+noc*(iw+1)-i+1           !new-sort
          Ekj = EIG(k) - EIG(j)
          Ekkij = 2.0d0*EIG(k) - EIG(i) - EIG(j)
          Xiijk = ERImol(i,i,k,j)
          Xijkk = ERImol(i,j,k,k)
          E2FERIex = E2FERIex - BBn(i,iw)*                              &
                    ( FOCKm(j,k)*Xiijk/Ekj + FOCKm(i,j)*Xijkk/Ekkij )
         enddo
        endif
       ENDDO
      ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      E2FERIex = 2.0d0*E2FERIex
!-----------------------------------------------------------------------
      RETURN
      END

! E2ERIERIex
      FUNCTION E2ERIERIex(NOC,NCWO,NORB,NBF,EIG,ERImol,BBn)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DOUBLE PRECISION,DIMENSION(NORB)::EIG
      DOUBLE PRECISION,DIMENSION(NOC,NCWO)::BBn
      DOUBLE PRECISION,DIMENSION(NBF,NBF,NBF,NBF)::ERImol
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     j/=m -> k,k         ( BBn(i,iw)*Xiiln*Xlnkk/Enlii is excluded )
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      E2ERIERIex = 0.0d0
      DO i=1,noc
       do iw=1,ncwo
        !k = noc+ncwo*(noc-i)+iw         !old-sort     ! k = (i,iw)
        k = noc+noc*(iw+1)-i+1           !new-sort
        DO m=1,noc
         DO j=1,noc
          if(j/=m)then         
           Ekkmj = 2.0d0*EIG(k) - EIG(m) - EIG(j)
           Xiimj = ERImol(i,i,j,m)
           Xmjkk = ERImol(m,j,k,k)
           E2ERIERIex = E2ERIERIex + BBn(i,iw)*Xiimj*Xmjkk/Ekkmj
          endif
         ENDDO
        ENDDO
       enddo
      ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END
      
! E2dHFex
      FUNCTION E2dHFex(NOC,NCWO,NORB,NBF,EIG,ERImol,BBn)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DOUBLE PRECISION,DIMENSION(NORB)::EIG
      DOUBLE PRECISION,DIMENSION(NOC,NCWO)::BBn
      DOUBLE PRECISION,DIMENSION(NBF,NBF,NBF,NBF)::ERImol
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     i/=j -> k,l
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      E2dHFex = 0.0d0
      DO i=1,noc
       do iw=1,ncwo
        !k = noc+ncwo*(noc-i)+iw         !old-sort     ! k = (i,iw)
        k = noc+noc*(iw+1)-i+1           !new-sort
        DO j=1,noc
         if(j/=i)then
          DO l=noc+1,norb
           Elkij = EIG(l) + EIG(k) - EIG(i) - EIG(j)
           Xijlk = ERImol(i,j,k,l)
           Xiljk = ERImol(i,l,k,j)
           Xijkl = ERImol(i,j,l,k)
           Xilkj = ERImol(i,l,j,k)
           XXijlk = Xijlk - Xijkl
           XXiljk = Xiljk - Xilkj
           E2dHFex = E2dHFex - (BBn(i,iw)/Elkij) *                      &
                           ( XXijlk*Xilkj + Xijlk*Xiljk + Xijkl*XXiljk )
          ENDDO
         endif
        ENDDO
       enddo
      ENDDO
      E2dHFex = 2.0d0*E2dHFex
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
!                                                                      !
!                     PNOF5(Nc)-PT2-NEX - > E2exi                      !
!                                                                      !
!      exclude excitations associated only to the omitted subspaces    !
!                                                                      !
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !

! E2HFs_nex
      FUNCTION E2HFs_nex(NOC,NCWO,NORB,EIG,FOCKm,IEX,NEX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      INTEGER,DIMENSION(NEX)::IEX
      DOUBLE PRECISION,DIMENSION(NORB)::EIG
      DOUBLE PRECISION,DIMENSION(NORB,NORB)::FOCKm
!-----------------------------------------------------------------------
!     excluding ix -> lx < norb
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      E2HFs_nex = 0.0d0
      DO i=1,nex
       ix = iex(i)
       do j=1,nex
        jx = iex(j)
        do jxw=1,ncwo  
         !lx = noc+ncwo*(noc-jx)+jxw         !old-sort     ! k = (i,iw)
         lx = noc+noc*(jxw+1)-jx+1           !new-sort
         E2HFs_nex = E2HFs_nex                                          &
                   + FOCKm(ix,lx)*FOCKm(lx,ix)/(EIG(lx)-EIG(ix))
        enddo
       enddo
      ENDDO
      E2HFs_nex = 2.0d0*E2HFs_nex
!-----------------------------------------------------------------------
      RETURN
      END

! E2HFdex_nex
      FUNCTION E2HFdex_nex(NOC,NCWO,NORB,NBF,EIG,ERImol,IEX,NEX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      INTEGER,DIMENSION(NEX)::IEX
      DOUBLE PRECISION,DIMENSION(NORB)::EIG
      DOUBLE PRECISION,DIMENSION(NBF,NBF,NBF,NBF)::ERImol
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     excluding mx,nx (mx/=nx) -> kx,lx
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      E2HFdex_nex = 0.0d0
      do m=1,nex
       mx = iex(m)
       do n=1,nex
        nx = iex(n)
        if(mx/=nx)then
         do i=1,nex
          ix = iex(i)
          do ixw=1,ncwo  
           !kx = noc+ncwo*(noc-ix)+ixw         !old-sort     ! k = (i,iw)
           kx = noc+noc*(ixw+1)-ix+1           !new-sort
           do j=1,nex
            jx = iex(j)
            do jxw=1,ncwo  
             !lx = noc+ncwo*(noc-jx)+jxw         !old-sort     ! k = (i,iw)
             lx = noc+noc*(jxw+1)-jx+1           !new-sort
             Xmnkl = ERImol(mx,nx,lx,kx)
             Xmnlk = ERImol(mx,nx,kx,lx)
             Eklmn = EIG(kx) + EIG(lx) - EIG(mx) - EIG(nx)
             E2HFdex_nex = E2HFdex_nex + Xmnkl*(2.0*Xmnkl-Xmnlk)/Eklmn
            enddo
           enddo
          enddo
         enddo
        endif
       enddo
      enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END

! E2F_nex
      FUNCTION E2F_nex(NOC,NCWO,NORB,EIG,FOCKm,BBn,IEX,NEX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      INTEGER,DIMENSION(NEX)::IEX
      DOUBLE PRECISION,DIMENSION(NORB)::EIG
      DOUBLE PRECISION,DIMENSION(NOC,NCWO)::BBn
      DOUBLE PRECISION,DIMENSION(NORB,NORB)::FOCKm
!-----------------------------------------------------------------------
!     excluding ix -> kx                     
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      E2F_nex = 0.0d0
      do i=1,nex
       ix = iex(i)
       do ixw=1,ncwo
        !kx = noc+ncwo*(noc-ix)+ixw         !old-sort     ! k = (i,iw)
        kx = noc+noc*(ixw+1)-ix+1           !new-sort
        Eki = EIG(kx) - EIG(ix)
        E2F_nex = E2F_nex - BBn(ix,ixw)*FOCKm(ix,kx)*FOCKm(kx,ix)/Eki 
       enddo
      enddo
      E2F_nex = 2.0d0*E2F_nex
!-----------------------------------------------------------------------
      RETURN
      END

! E2FERIex_nex
      FUNCTION E2FERIex_nex(NOC,NCWO,NORB,NBF,EIG,FOCKm,ERImol,BBn,     &
                            IEX,NEX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      INTEGER,DIMENSION(NEX)::IEX
      DOUBLE PRECISION,DIMENSION(NORB)::EIG
      DOUBLE PRECISION,DIMENSION(NOC,NCWO)::BBn
      DOUBLE PRECISION,DIMENSION(NORB,NORB)::FOCKm
      DOUBLE PRECISION,DIMENSION(NBF,NBF,NBF,NBF)::ERImol
!-----------------------------------------------------------------------
      E2FERIex_nex = 0.0d0
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     excluding ix -> lx
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      do i=1,nex
       ix = iex(i)
       do j=1,nex
        jx = iex(j)
        do jxw=1,ncwo
         !lx = noc+ncwo*(noc-jx)+jxw         !old-sort     ! k = (i,iw)
         lx = noc+noc*(jxw+1)-jx+1           !new-sort
         Eli = EIG(lx) - EIG(ix)
         do ixw=1,ncwo
          !kx = noc+ncwo*(noc-ix)+ixw         !old-sort     ! k = (i,iw)
          kx = noc+noc*(ixw+1)-ix+1           !new-sort
          Xilkk = ERImol(ix,lx,kx,kx)
          E2FERIex_nex = E2FERIex_nex                                   &
                       - BBn(ix,ixw)*FOCKm(ix,lx)*Xilkk/Eli
         enddo
        enddo
       enddo
      enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     excluding ix,jx (ix/=jx) -> kx,kx
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO j=1,nex
       jx = iex(j)
       DO i=1,nex
        ix = iex(i)
        if(ix/=jx)then       
         do ixw=1,ncwo
          !kx = noc+ncwo*(noc-ix)+ixw         !old-sort     ! k = (i,iw)
          kx = noc+noc*(ixw+1)-ix+1           !new-sort
          Ekj = EIG(kx) - EIG(jx)          
          Ekkij = 2.0d0*EIG(kx) - EIG(ix) - EIG(jx) 
          Xiijk = ERImol(ix,ix,kx,jx)          
          Xijkk = ERImol(ix,jx,kx,kx)
          E2FERIex_nex = E2FERIex_nex + BBn(ix,ixw)*                    &
                  ( FOCKm(jx,kx)*Xiijk/Ekj + FOCKm(ix,jx)*Xijkk/Ekkij )
         enddo     
        endif
       ENDDO
      ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      E2FERIex_nex = 2.0d0*E2FERIex_nex
!-----------------------------------------------------------------------
      RETURN
      END

! E2ERIERIex_nex
      FUNCTION E2ERIERIex_nex(NOC,NCWO,NORB,NBF,EIG,ERImol,BBn,IEX,NEX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      INTEGER,DIMENSION(NEX)::IEX
      DOUBLE PRECISION,DIMENSION(NORB)::EIG
      DOUBLE PRECISION,DIMENSION(NOC,NCWO)::BBn
      DOUBLE PRECISION,DIMENSION(NBF,NBF,NBF,NBF)::ERImol
!-----------------------------------------------------------------------
!     excluding ix,jx -> lx,lx                        
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      E2ERIERIex_nex = 0.0d0
      do i=1,nex
       ix = iex(i)
       do j=1,nex
        jx = iex(j)
        if(jx/=ix)then
         do m=1,nex
          mx=iex(m)
          do mxw=1,ncwo  
           !lx = noc+ncwo*(noc-mx)+mxw         !old-sort     ! k = (i,iw)
           lx = noc+noc*(mxw+1)-mx+1           !new-sort
           Xmmij = ERImol(mx,mx,jx,ix)
           Xijll = ERImol(ix,jx,lx,lx)
           Ellij = 2.0d0*EIG(lx) - EIG(ix) - EIG(jx)
           E2ERIERIex_nex = E2ERIERIex_nex                              &
                          - BBn(mx,mxw)*Xmmij*Xijll/Ellij
          enddo
         ENDDO
        endif
       enddo
      enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END

! E2dHFex_nex
      FUNCTION E2dHFex_nex(NOC,NCWO,NORB,NBF,EIG,ERImol,BBn,IEX,NEX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      INTEGER,DIMENSION(NEX)::IEX
      DOUBLE PRECISION,DIMENSION(NORB)::EIG
      DOUBLE PRECISION,DIMENSION(NOC,NCWO)::BBn
      DOUBLE PRECISION,DIMENSION(NBF,NBF,NBF,NBF)::ERImol
!-----------------------------------------------------------------------
!     excluding ix,jx -> kx,lx
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      E2dHFex_nex = 0.0d0
      do i=1,nex
       ix = iex(i)
       do ixw=1,ncwo
        !kx = noc+ncwo*(noc-ix)+ixw         !old-sort     ! k = (i,iw)
        kx = noc+noc*(ixw+1)-ix+1           !new-sort
        do j=1,nex
         jx = iex(j)
         if(jx/=ix)then
          do m=1,nex
           mx = iex(m)
           do mxw=1,ncwo  
            !lx = noc+ncwo*(noc-mx)+mxw         !old-sort     ! k = (i,iw)
            lx = noc+noc*(mxw+1)-mx+1           !new-sort
            Elkij = EIG(lx) + EIG(kx) - EIG(ix) - EIG(jx)
            Xijlk = ERImol(ix,jx,kx,lx)
            Xiljk = ERImol(ix,lx,kx,jx)
            Xijkl = ERImol(ix,jx,lx,kx)
            Xilkj = ERImol(ix,lx,jx,kx)
            XXijlk = Xijlk - Xijkl
            XXiljk = Xiljk - Xilkj
            E2dHFex_nex = E2dHFex_nex + (BBn(ix,ixw)/Elkij) *           &
                          (XXijlk*Xilkj+Xijlk*Xiljk+Xijkl*XXiljk)
           enddo
          enddo
         endif
        enddo
       enddo
      enddo
      E2dHFex_nex = 2.0d0*E2dHFex_nex
!-----------------------------------------------------------------------
      RETURN
      END

! E2Qd_nex
      FUNCTION E2Qd_nex(NOC,NCWO,NORB,NBF,EIG,ERImol,BBn,IEX,NEX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      INTEGER,DIMENSION(NEX)::IEX
      DOUBLE PRECISION,DIMENSION(NORB)::EIG
      DOUBLE PRECISION,DIMENSION(NOC,NCWO)::BBn
      DOUBLE PRECISION,DIMENSION(NBF,NBF,NBF,NBF)::ERImol
!-----------------------------------------------------------------------
!     excluding ix,jx -> kx,lx
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      E2Qd_nex = 0.0d0
      do i=1,nex
       ix = iex(i)
       do j=1,nex
        jx = iex(j)
        if(jx/=ix)then
         do ixw=1,ncwo
          !kx = noc+ncwo*(noc-ix)+ixw         !old-sort     ! k = (i,iw)
          kx = noc+noc*(ixw+1)-ix+1           !new-sort
          do jxw=1,ncwo
           !lx = noc+ncwo*(noc-jx)+jxw         !old-sort     ! k = (i,iw)
           lx = noc+noc*(jxw+1)-jx+1           !new-sort
           BB2 = BBn(ix,ixw)*BBn(jx,jxw)
           Xijkl = ERImol(ix,jx,lx,kx)
           Xijlk = ERImol(ix,jx,kx,lx)
           XX = (2.0d0*Xijkl-Xijlk)*(2.0d0*Xijkl-Xijlk)
           Eklij = EIG(kx) + EIG(lx) - EIG(ix) - EIG(jx)
           E2Qd_nex = E2Qd_nex + BB2 * XX / Eklij
          enddo
         enddo
        endif
       enddo
      enddo
!-----------------------------------------------------------------------
      RETURN
      END
