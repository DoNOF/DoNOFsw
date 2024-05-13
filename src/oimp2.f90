!======================================================================!
!                                                                      !
!                 N O F - M P 2   S U B R O U T I N E S                !
!                                                                      !
!            ( PRL 119, 063002, 2017; PRA 98, 022504, 2018 )           !
!                                                                      !
!======================================================================!
!======================================================================!

! ECorrNonDyn
      SUBROUTINE ECorrNonDyn(RO,QK,ENonDEnergy)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL HighSpin
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_ORBSPACE0/NO1,NDOC,NCO,NCWO,NVIR,NAC,NO0
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      COMMON/CorrNonDynamic/ECnd,ECndl,ECndHF
      COMMON/INPNOF_PNOF/IPNOF,NTWOPAR      
!
      DOUBLE PRECISION,DIMENSION(NBF5)::RO
      DOUBLE PRECISION,DIMENSION(NBFT5)::QK
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::BETA,FIs
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::CK12nd
!-----------------------------------------------------------------------
      ALLOCATE(FIs(NBF5),CK12nd(NBF5,NBF5))
      CK12nd = 0.0d0      
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --      
!                       FIs(j)  = 2*RO(j)*HOLE(j)      
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
      FIs = 0.0d0
      DO j=NO1+1,NBF5
       FIs(j) = 2.0d0*RO(j)*(1.0d0-RO(j))
      ENDDO
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!          Inter-pair Non-Dynamic (Static) Electron Correlation
!                  CK12ji = +- FIs(j)*FIs(i) (-+PIs)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
      DO j=NO1+1,NBF5
       DO i=NO1+1,NBF5
        CK12nd(j,i) = FIs(j)*FIs(i)        
       ENDDO
      ENDDO
!- - - - - - - - - - - - - - - - - - - - - -
      IF(IPNOF==8)THEN
       DO j=NO1+1,NB
        DO i=NO1+1,NB
         CK12nd(j,i) = 0.0d0
        ENDDO
       ENDDO
       if(MSpin==0.and.NSOC>0)then       
        DO j=NO1+1,NB           
         DO i=NB+1,NA                                              
         CK12nd(j,i) = 0.5d0*CK12nd(j,i)
         ENDDO                                                     
        ENDDO      
        DO j=NB+1,NA
         DO i=NO1+1,NB     
          CK12nd(j,i) = 0.5d0*CK12nd(j,i)
         ENDDO
        ENDDO
       endif
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - -
      DEALLOCATE(FIs)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!           Intra-pair Non-Dynamic (Static) Electron Correlation
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!     BETA(j) = DSQRT(2*HOLE(j)) * DSQRT(RO(j))
!     BETA(j) = DSQRT(2*RO(j))   * DSQRT(RO(j))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(BETA(NBF5))
      DO j=1,NBF5
       Cj = 1.0d0 - DABS(1.0d0-2.0d0*RO(j))
       BETA(j) = DSQRT( Cj*RO(j) )
      ENDDO
!      
      DO l=1,NDOC
       ln = NO1+l
       DO i=NDNS+NCWO*(NDOC-l)+1,NDNS+NCWO*(NDOC-l+1)
        in = NO1+i
        CK12nd(ln,in) = BETA(ln)*BETA(in)
        CK12nd(in,ln) = BETA(in)*BETA(ln)
        DO j=NDNS+NCWO*(NDOC-l)+1,NDNS+NCWO*(NDOC-l+1)
         jn = NO1+j
         CK12nd(jn,in) = - BETA(jn)*BETA(in)         
        ENDDO
       ENDDO
      ENDDO
      DEALLOCATE(BETA)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!     Non-Dynamic Correlation Energy
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(MSpin==0)then
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -      
       ENonDEnergy = 0.0d0
       ECndHF = 0.0d0       
       do i=1,ndoc
        in = no1+i
        ENonDEnergy = ENonDEnergy - PRODCWQWj(in,CK12nd,QK)
        do iw=1,ncwo
         in = na+ncwo*(ndoc-i)+iw
         ENonDEnergy = ENonDEnergy -PRODCWQWj(in,CK12nd,QK)
        enddo
       enddo
       if(NSOC>0)then
        do i=ndoc+1,ndns
         in = no1+i
         ENonDEnergy = ENonDEnergy - PRODCWQWj(in,CK12nd,QK) 
         ECndHF = ECndHF - CK12nd(in,in)*QK(in*(in+1)/2)             
        end do
       end if
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      else if(MSpin>0)then
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -      
       ENonDEnergy = 0.0d0
       do i=1,ndoc
        in = no1+i
        ENonDEnergy = ENonDEnergy - PRODCWQWj1(in,CK12nd,QK)
        do iw=1,ncwo
         in = na+ncwo*(ndoc-i)+iw
         ENonDEnergy = ENonDEnergy -PRODCWQWj2(in,CK12nd,QK)
        enddo
       enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -      
      end if
!-----------------------------------------------------------------------
      DEALLOCATE(CK12nd)
      RETURN
      END

! ORBINVMP2 (oiMP2)
      SUBROUTINE ORBINVMP2(ELAG,COEF,RO,CJ12,CK12,AHCORE,IERI,ERI,      &
                           ADIPx,ADIPy,ADIPz)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL HighSpin
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_ORBSPACE0/NO1,NDOC,NCO,NCWO,NVIR,NAC,NO0
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
      COMMON/INPNOF_GENERALINF/ICOEF,ISOFTMAX,IORBOPT,MAXIT,MAXIT21
      COMMON/EHFEN/EHF,EN
      COMMON/INPFILE_NO1PT2/NO1PT2,NEX
      COMMON/INPFILE_NIJKL/NINTMX,NIJKL,NINTCR,NSTORE
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      COMMON/CorrNonDynamic/ECnd,ECndl,ECndHF
      COMMON/INPNOF_Tijab/NOUTTijab,NTHRESHTijab,THRESHTijab
      COMMON/NumLinIndOrb/NQMT
!
      INTEGER,DIMENSION(NIJKL) :: IERI
      DOUBLE PRECISION,DIMENSION(NIJKL) :: ERI
      DOUBLE PRECISION,DIMENSION(NBF5) :: RO
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5) :: CJ12,CK12
      DOUBLE PRECISION,DIMENSION(NBF,NBF) :: AHCORE,ADIPx,ADIPy,ADIPz
      DOUBLE PRECISION,DIMENSION(NBF,NBF) :: ELAG,COEF
!      
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: OCC,EIG,FI1,FI2,Tijab
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: VEC,FOCKm
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: ERImol
!-----------------------------------------------------------------------
!     NO1:  Number of inactive doubly occupied orbitals (OCC=1)         
!     NDOC: Number of strongly doubly occupied MOs                      
!     NSOC: Number of strongly singly occupied MOs                      
!     NDNS: Number of strongly occupied MOs (NDNS=NDOC+NSOC)                        
!     NCWO: Number of coupled weakly occ. MOs per strongly doubly occ.
!     NCWO*NDOC: Active orbitals in the virtual subspace                
!     NO0:  Empty orbitals  (OCC=0)                                      
!     NVIR: Number of weakly occupied MOs + empty MOs                   
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!           NO1 | NDOC  + NSOC  |   NCWO*NDOC + NO0  = NBF               
!           NO1 |      NDNS     |          NVIR      = NBF 
!               | -NAC- |       |  -   NAC  - |
!                      NB      NA            NBF5
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                               !  CLOSED (NB=NA=NCO,NSOC=0)
!     NDOC = NB - NO1           !  NDOC = NCO - NO1, NO1 <= NCO
!     NDNS = NDOC + NSOC        !  NDNS = NDOC      
!     NA   = NO1 + NDNS         !  NA = NB = NCO
!     NVIR = NBF - NA           !  NBF - NCO
!-----------------------------------------------------------------------
!     Calculation of the Matrix of Lagrange Multipliers if ICOEF=0
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(ICOEF==0)CALL ELAGCOEF0(ELAG,COEF,RO,CJ12,CK12,AHCORE,         &
                                 ADIPx,ADIPy,ADIPz,IERI,ERI)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Substracting core orbitals (NO1PT2)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      NOCB = NB - NO1PT2
      NOC = NA - NO1PT2
      NVI = NQMT - NA 
      NORB = NQMT - NO1PT2
      WRITE(6,1)NO1PT2,NOC,NVI,NORB
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     NOCNC = NOC + NCWO*NOCB <= NORB
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      NOCNC = NOC + NCWO*NOCB
      NOCNCT = NOCNC*(NOCNC+1)/2
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     If NOCNC > NORB: Check the value of NCWO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(NOCNC>NORB)THEN
       NCWO1 = (NORB-NSOC)/NOCB - 1 
       IF(NCWO>NCWO1)THEN
        WRITE(6,4)                                                      &
        'NCWO is too large, reduce the value at least to',NCWO1        
        STOP
       ENDIF
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     If NO1PT2<NO1: Check the value of NCWO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      NCWO1 = (NQMT-NO1-NSOC)/NOC - 1 
      IF(NO1PT2<NO1.and.NCWO>NCWO1)THEN
       WRITE(6,4)                                                       &
       'NCWO is too large, reduce the value at least to',NCWO1        
       STOP
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     OCC,VEC: Occupations and Molecular Orbitals without Core
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(OCC(NOCNC),VEC(NBF,NORB))
      DO i=1,NOC
       OCC(i) = RO(i+NO1PT2)
      ENDDO
      DO i=NOC+1,NOCNC
       if(NOCNC+NO1PT2<=NBF5)then      
        OCC(i) = RO(i+NO1PT2)
       else
        OCC(i) = 0.0d0       
       end if
      ENDDO
      DO J=1,NORB
       do i=1,NBF
        VEC(i,J) = COEF(i,J+NO1PT2)
       enddo
      ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     EIG: Energy Eigenvalues without Core
!     FOCKm: Fock Matrix in the Molecular Basis
!     EHFL = E(0) + E(1) : HF Energy with Non-HF Orbitals (COEF/=CHF)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(EIG(NORB),FOCKm(NORB,NORB))
      CALL FOCKMOL(NORB,COEF,VEC,ELAG,EIG,FOCKm,AHCORE,IERI,ERI,EHFL)      
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Integral transformation ( Form <ab|rs> )
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(ERImol(NOC,NBFT,NBF))
!     FORM (pq/rj)
      CALL ERIC1(ERImol,IERI,ERI,VEC,NOC,NORB)
!     FORM (ai/rj) for all a,i
      CALL ERIC23(ERImol,VEC,NVI,NOC,NORB)
!     FORM (ai/bj) for all b
      CALL ERIC4(ERImol,VEC(1,NOC+1),NOC,NVI,NORB)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!           Calculate second-order Dynamic Correlation E(2)
!         Calculate excitation amplitudes to determine Psi(1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      NN1 = NOC*NOC*NVI*NVI
      NN2 = NOC*NVI
      ALLOCATE (Tijab(NN1),FI1(NORB),FI2(NORB))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                              FI1 and FI2
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     like in PRL (NOF-MP2) (best result with PNOF7s)
!     Note: Ci=DSQRT[4*(1-OCCi)*OCCi],FI2i=1-Ci*Ci -> (1-2*OCCi)**2
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      FI1 = 1.0d0
      DO i=1,NOCNC             ! best for H2
       Ci = 1.0d0 - DABS(1.0d0-2.0d0*OCC(i))
       FI1(i) = 1.0d0 - Ci*Ci
      ENDDO
!- - - - - - - - - - - - - - - - - - - - - - -      
      FI2 = 1.0d0    
      DO i=NOC+1,NOCNC
       Ci = DABS(1.0d0-2.0d0*OCC(i))
       FI2(i) = Ci*Ci
      ENDDO
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!     Calculate Tijab amplitude solving Sparse Sym. Linear System
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
      CALL CALTijabIsym(NOCB,NOC,NVI,NORB,NN1,EIG,FOCKm,ERImol,Tijab,   &
                        FI1,FI2)                                         
      if(NOUTTijab==1)CALL OUTPUTTijab_rc(NOC,NVI,NN1,Tijab)             
      CALL ORBINVE2Totalsym(NOCB,NOC,NVI,NN1,NBF,NBFT,ERImol,           &
                            Tijab,ECd)
      WRITE(6,2)EHFL+EN+ECndHF
      WRITE(6,3)ECd,ECndl,ECd+ECndl,EHFL+ECd+ECndl+EN+ECndHF
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
      DEALLOCATE(OCC,VEC,FI1,FI2,EIG,FOCKm,ERImol,Tijab)
      RETURN
!-----------------------------------------------------------------------
    1 FORMAT(//,2X,'NOF-MP2',/1X,9('='),//,                             &
               1X,'NUMBER OF CORE ORBITALS  (NO1PT2) =',I5,/,           &
               1X,'NUMBER OF Doubly OCC. ORBS. (NOC) =',I5,/,           &
               1X,'NUMBER OF VIRTUAL ORBS.     (NVI) =',I5,/,           &
               1X,'NUMBER OF LIN. IND. ORBS.  (NORB) =',I5)              
    2 FORMAT(/3X,'          Ehfc =',F20.10)                              
    3 FORMAT(/3X,'           ECd =',F20.10,/,                           &
              3X,'          ECnd =',F20.10,/,                           &
              3X,'         ECorr =',F20.10,/,                           &
              3X,'     E(NOFMP2) =',F20.10)
    4 FORMAT(/1X,A47,I5)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END

! CALTijabIsym
      SUBROUTINE CALTijabIsym(NOCB,NOC,NVI,NORB,NN,EIG,FOCKm,           &
                              ERImol,Tijab,FI1,FI2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_CGM/ICGMETHOD
      COMMON/INPNOF_ORBSPACE0/NO1,NDOC,NCO,NCWO,NVIR,NAC,NO0
!
      DOUBLE PRECISION,DIMENSION(NN) :: Tijab
      DOUBLE PRECISION,DIMENSION(NORB) :: EIG,FI1,FI2
      DOUBLE PRECISION,DIMENSION(NORB,NORB) :: FOCKm
      DOUBLE PRECISION,DIMENSION(NOC,NBFT,NBF) :: ERImol
      INTEGER,ALLOCATABLE,DIMENSION(:) :: IROW,ICOL,NPAIR
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: A,B
!-----------------------------------------------------------------------
!     NPAIR: number of the pair to which the virtual orbital belongs
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(npair(nvi))
      do i=1,nocb
       do iw=1,ncwo
        k = ncwo*(nocb-i) + iw
        npair(k) = i
       end do
      end do
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     NN:  the order of the original sparse matrix
!     NNZ: the number of non-zero elements in the lower triangular of A
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      NNZM = 2*NOC*NOC*NVI*NVI*NORB
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     A: array (NNZ) that contains the non-zero elements of the lower 
!        triangular part of the original matrix, ordered by increasing
!        row index and by increasing column index within each row
!     IROW: the row indices of non-zero elements given in A
!     ICOL: the column indices of non-zero elements given in A
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(A(NNZM),IROW(NNZM),ICOL(NNZM))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Elements A(ijab,i'j'a'b')
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      nnz = 0
      nnz_end = 0

      do ib=1,nvi
       do ia=1,nvi
        do j=1,noc
         do i=1,noc

          jab =     (j-1)*noc + (ia-1)*noc*noc + (ib-1)*noc*noc*nvi
          iab = i             + (ia-1)*noc*noc + (ib-1)*noc*noc*nvi
          ijb = i + (j-1)*noc                  + (ib-1)*noc*noc*nvi
          ija = i + (j-1)*noc + (ia-1)*noc*noc
          ijab= i + (j-1)*noc + (ia-1)*noc*noc + (ib-1)*noc*noc*nvi

          nnz = nnz + 1
          nnz_ini = nnz_end + 1
          A(nnz) = FOCKm(ia+noc,ia+noc) + FOCKm(ib+noc,ib+noc)          &
                 - FOCKm(i,i) - FOCKm(j,j)    ! diagonal elem.

          IROW(nnz) = ijab
          ICOL(nnz) = i + jab

          do k=1,i-1
           IF(DABS(FOCKm(i,k))>1.0d-10)THEN
            nnz = nnz + 1
            Cki = FI2(k)*FI2(i)
            A(nnz) = - Cki*FOCKm(i,k)
            IROW(nnz) = ijab
            ICOL(nnz) = k + jab
           ENDIF
          end do

          do k=1,j-1
           IF(DABS(FOCKm(j,k))>1.0d-10)THEN
            nnz = nnz + 1
            Ckj = FI2(k)*FI2(j)
            A(nnz) = - Ckj*FOCKm(j,k)
            IROW(nnz) = ijab
            ICOL(nnz) = (k-1)*noc + iab
           ENDIF
          end do

          do k=1,ia-1
           IF(DABS(FOCKm(ia+noc,k+noc))>1.0d-10)THEN
            nnz = nnz + 1
            if(npair(k)==npair(ia))then
             Ckia = FI1(k+noc)*FI1(ia+noc)
            else
             Ckia = FI2(k+noc)*FI2(ia+noc)
            endif
            A(nnz) = Ckia*FOCKm(ia+noc,k+noc)
            IROW(nnz) = ijab
            ICOL(nnz) = (k-1)*noc*noc + ijb
           ENDIF
          end do

          do k=1,ib-1
           IF(DABS(FOCKm(ib+noc,k+noc))>1.0d-10)THEN
            nnz = nnz + 1
            if(npair(k)==npair(ib))then
             Ckib = FI1(k+noc)*FI1(ib+noc)
            else
             Ckib = FI2(k+noc)*FI2(ib+noc)
            endif
            A(nnz) = Ckib*FOCKm(ib+noc,k+noc)
            IROW(nnz) = ijab
            ICOL(nnz) = (k-1)*noc*noc*nvi + ija
           ENDIF
          end do

          nnz_end = nnz
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!         Sort Column Index for a given IROW = ijab
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          do k=nnz_ini,nnz_end-1
           do l=k+1,nnz_end
            if(ICOL(k)>ICOL(l))then
             ITEMPCOL = ICOL(l)
             TEMPVECT = A(l)
              ICOL(l) = ICOL(k)
                 A(l) = A(k)
              ICOL(k) = ITEMPCOL
                 A(k) = TEMPVECT
            endif
           end do
          end do
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         end do
        end do      
       end do
      end do
      DEALLOCATE(NPAIR)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Tijab: Initial approximation for the vector Tijab
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      iai = 0
      do ia=1,nvi
       do i=1,noc
        iai = iai+1
        do ib=1,nvi
         do j=1,noc
          ijab = i + (j-1)*noc + (ia-1)*noc*noc + (ib-1)*noc*noc*nvi
          Eijab = EIG(ib+noc) + EIG(ia+noc) - EIG(j) - EIG(i)
          Tijab(ijab) = - ERImol(j,iai,ib)/Eijab
         end do
        end do      
       end do
      end do
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     B: the right-hand side vector
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(B(NN))
      do i=1,noc
       lmin_i = noc+ncwo*(noc-i)+1
       lmax_i = noc+ncwo*(noc-i)+ncwo
       do j=1,noc
        if(j==i)then
!--------------------------------------------------------------------
!        excitations inside the subspace 'i'
!--------------------------------------------------------------------
         do k=1,nvi
          ik = i + (k-1)*noc
          kn = k + noc
          do l=1,nvi
           ln = l + noc
           if(      (lmin_i<=kn.and.kn<=lmax_i)                         &
              .and. (lmin_i<=ln.and.ln<=lmax_i) )then
            Ciikl = FI1(kn)*FI1(ln)*FI1(i)*FI1(i)
           else
            Ciikl = FI2(kn)*FI2(ln)*FI2(i)*FI2(i)
           endif
           iikl =  i + (i-1)*noc + (k-1)*noc*noc + (l-1)*noc*noc*nvi
           B(iikl) = - Ciikl*ERImol(i,ik,l)
          enddo
         enddo
!--------------------------------------------------------------------
        else
!--------------------------------------------------------------------
         do k=1,nvi
          ik = i + (k-1)*noc
          kn = k + noc
          do l=1,nvi
           ln = l + noc
           ijkl =  i + (j-1)*noc + (k-1)*noc*noc + (l-1)*noc*noc*nvi
           Cijkl = FI2(kn)*FI2(ln)*FI2(i)*FI2(j)
           B(ijkl) = - Cijkl*ERImol(j,ik,l)
          enddo
         enddo
!--------------------------------------------------------------------
        endif
       enddo
      enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Solving the Sparse Linear System AT=B using the CG method
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(ICGMETHOD==1)THEN
       CALL SparseSymLinearSystem_CG(NN,NNZ,A,IROW,ICOL,B,Tijab)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Prepare for calling the subroutine F11JEF of the NAG Library
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ELSE IF(ICGMETHOD==2)THEN
       CALL SparseSymLinearSystem_NAG(NN,NNZ,A,IROW,ICOL,B,Tijab)
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DEALLOCATE(A,IROW,ICOL,B)
      RETURN
!-----------------------------------------------------------------------
      END

! SparseSymLinearSystem_NAG
      SUBROUTINE SparseSymLinearSystem_NAG(NN,NNZ,A,IROW,ICOL,B,Tijab)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      CHARACTER PRECON
      CHARACTER*6 METHOD
      INTEGER,DIMENSION(NNZ) :: IROW,ICOL
      DOUBLE PRECISION,DIMENSION(NNZ) :: A
      DOUBLE PRECISION,DIMENSION(NN) :: B,Tijab
      INTEGER,ALLOCATABLE,DIMENSION(:) :: IWORK
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: WORK
!-----------------------------------------------------------------------
!     The iterative method to be used:
!     'CG'      ==> Conjugate gradient method
!     'SYMMLQ'  ==> Lanczos method (SYMMLQ).
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      METHOD = 'CG'
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     The type of preconditioning used:
!     'N'       ==> no preconditioning
!     'J'       ==> Jacobi preconditioning.
!     'S'       ==> symmetric SOR preconditioning.
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      PRECON = 'N'
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     OMEGA ( 0.0 < OMEGA < 2.0 )
!     if PRECON = 'S' then OMEGA is the relaxation parameter to be used
!     Otherwise OMEGA need not be initialised.
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      OMEGA = 1.0d0
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     TOL: the tolerance required ( TOL < 1.0 )
!          The iteration is judged to have converged at step k if:
!
!               || r || <= TOL*(|| b || + || A ||*|| x ||).
!                   k                                 k
!     If TOL = 0.0 the default value of SQRT(machine precision) is used.
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      TOL = 1.0d-10
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     MAXITN: the maximum number of iterations allowed.
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      MAXITN = 3000
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     LWORK: the dimension of the array WORK
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(METHOD=='CG')THEN
       IF(PRECON=='J'.or.PRECON=='S')THEN
        LWORK = 7*NN + 120
       ELSE
        LWORK = 6*NN + 120
       END IF
      ELSEIF(METHOD=='SYMMLQ')THEN
       IF(PRECON=='J'.or.PRECON=='S')THEN
        LWORK = 8*NN + 120
       ELSE
        LWORK = 7*NN + 120
       END IF
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                                 AX = B
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     NN:  the order of the original sparse matrix
!     NNZ: the number of non-zero elements in the lower triangular of A
!     A: array (NNZ) that contains the non-zero elements of the lower 
!        triangular part of the original matrix, ordered by increasing
!        row index and by increasing column index within each row
!     IROW: the row indices of non-zero elements given in A
!     ICOL: the column indices of non-zero elements given in A
!     B: array (NN), the right-hand side vector
!     RNORM: the final residual norm
!     ITN: the actual number of iterations taken
!     WORK: working REAL array of dimension LWORK
!     IWORK: working INTEGER array of dimension NN+1
!
!     IFAIL (On entry must be -1, 0, or 1)
!      On exit, the following values may occur:
!       IFAIL = 0 => no error detected.
!       IFAIL = 1 => METHOD invalid, or PRECON invalid, or
!                    NN < 1, or NNZ < 1, or NNZ > N*(N+1)/2, or
!                    OMEGA not in (0,2), or TOL >= 1.0, or
!                    MAXITN < 1, or LWORK too small.
!       IFAIL = 2 => the arrays IROW and ICOL fail to satisfy:
!                    1 <= IROW(i) <= N,
!                    1 <= ICOL(i) <= IROW(i), i = 1,2,...,NNZ.
!                    IROW(i-1) < IROW(i), or
!                    IROW(i-1) = IROW(i) and ICOL(i-1) < ICOL(i)
!                    for i = 2,3,...,NNZ.
!       IFAIL = 3 => A has a zero diagonal element.
!       IFAIL = 4 => A reasonable accuracy has been obtained and
!                    further iterations could not improve the result.
!       IFAIL = 5 => Required accuracy not obtained in MAXITN
!       IFAIL = 6 => The preconditioner is not positive definite.
!       IFAIL = 7 => The matrix A is not positive definite (CG only).
!       IFAIL = 8 => A serious error has occurred in an internal call.
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(WORK(LWORK),IWORK(NN+1))
      IFAIL = 1
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -      
!     Avoiding warnings
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      A(1) = A(1)
      IROW(1) = IROW(1)
      ICOL(1) = ICOL(1)
      B(1) = B(1)
      Tijab(1) = Tijab(1)
!      CALL F11JEF(METHOD,PRECON,NN,NNZ,A,IROW,ICOL,OMEGA,B,TOL,         &   !nag
!                  MAXITN,Tijab,RNORM,ITN,WORK,LWORK,IWORK,IFAIL)            !nag
      write(6,'(/,A8,E11.4,A7,I5)')' RNORM =',RNORM,', ITN =',ITN
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DEALLOCATE(WORK,IWORK)
      RETURN
!-----------------------------------------------------------------------
      END

! SparseSymLinearSystem_CG      
      SUBROUTINE SparseSymLinearSystem_CG(NN,NNZ,A,IROW,ICOL,B,Tijab)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      INTEGER,DIMENSION(NNZ) :: IROW,ICOL
      DOUBLE PRECISION,DIMENSION(NNZ) :: A
      DOUBLE PRECISION,DIMENSION(NN) :: B,Tijab
      INTEGER,ALLOCATABLE,DIMENSION(:) :: IIROW,IICOL
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: AA,AP,R,P
!-----------------------------------------------------------------------
!     Symmetric A -> Square AA
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      NNZ2 = NNZ*2
      ALLOCATE(IIROW(NNZ2),IICOL(NNZ2),AA(NNZ2))
      IIROW(1:NNZ) = IROW(1:NNZ) 
      IICOL(1:NNZ) = ICOL(1:NNZ)
      AA(1:NNZ) = A(1:NNZ)
      NZ = NNZ
      DO i=1,NNZ
       if(IIROW(i)>IICOL(i))then
        NZ = NZ + 1
        IIROW(NZ) = IICOL(i)
        IICOL(NZ) = IIROW(i)       
        AA(NZ)    = AA(i)
       end if 
      END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Initial Values for AP, R, P
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(AP(NN),R(NN),P(NN))
      CALL BeqProdAT(NN,NZ,IIROW,IICOL,AA,Tijab,AP)
      R(1:NN) = B(1:NN) - AP(1:NN)
      P(1:NN) = B(1:NN) - AP(1:NN)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Do Nsteps of CG method
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Nsteps = NN
      DO i = 1,Nsteps
       CALL BeqProdAT(NN,NZ,IIROW,IICOL,AA,P,AP)
       PR  = DOT_PRODUCT(P,R)                   
       PAP = DOT_PRODUCT(P,AP)                  
       IF(PAP==0.0d0)RETURN
       Tijab(1:NN) = Tijab(1:NN) + (PR/PAP)*P(1:NN)
       R(1:NN) = R(1:NN) - (PR/PAP)*AP(1:NN)       
       RAP = DOT_PRODUCT(R,AP)                
       P(1:NN) = R(1:NN) - (RAP/PAP)*P(1:NN)  
      END DO
!-----------------------------------------------------------------------
      DEALLOCATE(IIROW,IICOL,AA,AP,R,P)
      RETURN
      END

! BeqProdAT      
      SUBROUTINE BeqProdAT(N,NZ,IROW,ICOL,A,T,B )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      INTEGER,DIMENSION(NZ) :: ICOL,IROW
      DOUBLE PRECISION,DIMENSION(NZ) :: A
      DOUBLE PRECISION,DIMENSION(N) :: T,B
!-----------------------------------------------------------------------      
      B = 0.0d0
      do k = 1,NZ
       i = IROW(k)
       j = ICOL(k)
       B(i) = B(i) + A(k)*T(j)
      end do
!-----------------------------------------------------------------------
      RETURN
      END

! ORBINVE2Totalsym
      SUBROUTINE ORBINVE2Totalsym(NOCB,NOC,NVI,NN,NBF,NBFT,ERImol,      &
                                  Tijab,E2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DOUBLE PRECISION,DIMENSION(NN) :: Tijab
      DOUBLE PRECISION,DIMENSION(NOC,NBFT,NBF) :: ERImol
!-----------------------------------------------------------------------
      E2 = 0.0d0
      DO k=1,nvi
       DO l=1,nvi
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -       
        DO i=1,nocb
         ki = (k-1)*noc + i
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                  
         DO j=1,nocb
          Xijkl = ERImol(j,ki,l)
          ijkl = i+(j-1)*noc+(k-1)*noc*noc+(l-1)*noc*noc*nvi
          ijlk = i+(j-1)*noc+(l-1)*noc*noc+(k-1)*noc*noc*nvi
          E2 = E2 + Xijkl*(2.0d0*Tijab(ijkl)-Tijab(ijlk))
         ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
         DO j=nocb+1,noc
          Xijkl = ERImol(j,ki,l)
          ijkl = i+(j-1)*noc+(k-1)*noc*noc+(l-1)*noc*noc*nvi
          ijlk = i+(j-1)*noc+(l-1)*noc*noc+(k-1)*noc*noc*nvi
          E2 = E2 + Xijkl*(Tijab(ijkl)-0.5d0*Tijab(ijlk))
         ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                  
        ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        DO i=nocb+1,noc
         ki = (k-1)*noc + i
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                  
         DO j=1,nocb        
          Xijkl = ERImol(j,ki,l)
          ijkl = i+(j-1)*noc+(k-1)*noc*noc+(l-1)*noc*noc*nvi
          ijlk = i+(j-1)*noc+(l-1)*noc*noc+(k-1)*noc*noc*nvi
          E2 = E2 + Xijkl*(Tijab(ijkl)-0.5d0*Tijab(ijlk))
         ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
         DO j=nocb+1,noc
          Xijkl = ERImol(j,ki,l)
          ijkl = i+(j-1)*noc+(k-1)*noc*noc+(l-1)*noc*noc*nvi
          ijlk = i+(j-1)*noc+(l-1)*noc*noc+(k-1)*noc*noc*nvi
          if(j/=i)E2 = E2 + Xijkl*(Tijab(ijkl)-0.5d0*Tijab(ijlk))/2.0
         ENDDO
        ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -        
       ENDDO
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !        
!                                                                      !
!     ERIC1: FORM (PQ/RJ)                                              !
!    ERIC23: FORM (AI/RJ) for all A,I                                  !
!     ERIC4: FORM (AI/BJ) for all B                                    !
!   FOCKMOL: Form the Fock Matrix in the Molecular Basis               !
!   DENMATHFr: Calculate the density matrix (DM) using HF MOs (CHF)    !
!   FORM2JK: Calculate 2J-K in the MO basis to form the Fock matrix    ! 
!    FORMJK: Calculate J-K in the MO basis to form the Fock matrix     !
!   DENMATHF05ro: Calculate the open shell part of density matrix (DM) !
!   ELAGCOEF0: Calculate Matrix of Lagrange Multipliers if ICOEF=0     !
!                                                                      !
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !        

! ERIC1
      SUBROUTINE ERIC1(ERImol,IERI,ERI,VEC,NOC,NORB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)  
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPFILE_NIJKL/NINTMX,NIJKL,NINTCR,NSTORE
!
      INTEGER,DIMENSION(NSTORE)::IERI
      DOUBLE PRECISION,DIMENSION(NSTORE)::ERI
      DOUBLE PRECISION,DIMENSION(NOC,NBFT,NBF)::ERImol
      DOUBLE PRECISION,DIMENSION(NBF,NORB)::VEC
!-----------------------------------------------------------------------
      ERImol=0.0d0
      DO M=1,NINTCR
       XINT1 = ERI(M)
       XINT2 = XINT1
       LABEL = IERI(M)
       CALL LABELIJKL(LABEL,I,J,K,L)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       IJ = I*(I-1)/2 + J
       KL = K*(K-1)/2 + L
       IF(I==J)XINT1 = XINT1 + XINT1
       IF(K==L)XINT2 = XINT2 + XINT2
       DO N=1,NOC
        ERImol(N,IJ,K) = ERImol(N,IJ,K) + XINT1*VEC(L,N)
        ERImol(N,IJ,L) = ERImol(N,IJ,L) + XINT1*VEC(K,N)
        ERImol(N,KL,I) = ERImol(N,KL,I) + XINT2*VEC(J,N)
        ERImol(N,KL,J) = ERImol(N,KL,J) + XINT2*VEC(I,N)
       ENDDO
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! ERIC23
      SUBROUTINE ERIC23(ERImol,VEC,NVI,NOC,NORB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      DOUBLE PRECISION,DIMENSION(NOC,NBFT,NBF)::ERImol
      DOUBLE PRECISION,DIMENSION(NBF,NORB)::VEC
      ALLOCATABLE::AUX1(:,:),AUX2(:,:)
      ALLOCATE(AUX1(NBF,NBF),AUX2(NBF,NBF))
!-----------------------------------------------------------------------
      DO MR=1,NBF
       DO MJ=1,NOC
!       ERImol -> AUX1
        MPQ = 0
        DO MP=1,NBF
         DO MQ=1,MP
          MPQ = MPQ + 1
          AUX1(MP,MQ) = ERImol(MJ,MPQ,MR)
          AUX1(MQ,MP) = AUX1(MP,MQ)
         ENDDO
        ENDDO
!       AUX2=AUX1*C
        DO IQ=1,NBF
         DO JQ=1,NOC
          AUX2(IQ,JQ)=0.0d0
          do i=1,NBF
           AUX2(IQ,JQ)=AUX2(IQ,JQ)+AUX1(IQ,i)*VEC(i,JQ)
          enddo
         ENDDO
        ENDDO
!       AUX1=VECt*AUX2
        CALL CeqAtB(AUX1,VEC(1,NOC+1),NBF,NVI,AUX2,NOC)
!       AUX1 -> ERImol
        MAI=0
        DO MA=1,NVI
         DO MI=1,NOC
          MAI=MAI+1
          ERImol(MJ,MAI,MR) = AUX1(MA,MI)
         ENDDO
        ENDDO
       ENDDO
      ENDDO
!-----------------------------------------------------------------------
      DEALLOCATE(AUX1,AUX2)
      RETURN
      END

! ERIC4
      SUBROUTINE ERIC4(ERImol,VEC,NOC,NVI,NORB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      DOUBLE PRECISION,DIMENSION(NOC,NBFT,NBF)::ERImol
      DOUBLE PRECISION,DIMENSION(NBF,NORB)::VEC
      ALLOCATABLE::AUX1(:,:),AUX2(:,:)
      ALLOCATE(AUX1(NBF,NBF),AUX2(NBF,NBF))
!-----------------------------------------------------------------------
      MAI = 0
      DO MA=1,NVI
       DO MI=1,NOC
!       ERImol -> AUX1
        MAI = MAI+1
        DO MJ=1,NOC
         DO MR=1,NBF
          AUX1(MR,MJ) = ERImol(MJ,MAI,MR)
         ENDDO
        ENDDO
!       AUX2=VEC*AUX1
        CALL CeqAtB(AUX2,VEC,NBF,NVI,AUX1,NOC)
!       AUX2 -> ERImol  
        DO MB=1,NVI
         DO MJ=1,NOC
          ERImol(MJ,MAI,MB) = AUX2(MB,MJ)
         ENDDO
        ENDDO
       ENDDO
      ENDDO
!-----------------------------------------------------------------------
      DEALLOCATE(AUX1,AUX2)
      RETURN
      END

! FOCKMOL
      SUBROUTINE FOCKMOL(NORB,COEF,VEC,ELAG,EIG,FOCKm,AHCORE,IERI,ERI,  &
                         EHFL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      LOGICAL HighSpin
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
      COMMON/INPFILE_NO1PT2/NO1PT2,NEX
      COMMON/INPFILE_NIJKL/NINTMX,NIJKL,NINTCR,NSTORE
!
      INTEGER,DIMENSION(NIJKL) :: IERI
      DOUBLE PRECISION,DIMENSION(NIJKL) :: ERI
      DOUBLE PRECISION,DIMENSION(NORB) :: EIG
      DOUBLE PRECISION,DIMENSION(NBF,NBF) :: COEF,ELAG,AHCORE
      DOUBLE PRECISION,DIMENSION(NBF,NORB) :: VEC
      DOUBLE PRECISION,DIMENSION(NORB,NORB) :: FOCKm
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: DMhf,FOCK,DMa
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: FOCKa,AUX,TVEC
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     DMhf: HF Density Matrix
!     FOCK: Fock Matrix in the atomic basis set ( FOCK = H + 2J-K )
!     EHFL = E(0) + E(1) : HF Energy with Non-HF Orbitals (COEF/=CHF)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE (DMhf(NBF,NBF),FOCK(NBF,NBF))
      if(MSpin==0)then               ! Singlet and Multiplet States 
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       CALL DENMATHFr(DMhf,COEF)     ! 2*COEF(M,J)*COEF(N,J)
!- - - - - - - - - - - - - - - - - - - -      
       if(NSOC>0)then
        ALLOCATE (DMa(NBF,NBF))
        CALL DENMATHF05ro(DMa,COEF)                 ! NB+1,NA       
        DMhf = DMhf + DMa
       endif         
!- - - - - - - - - - - - - - - - - - - -
       CALL FORM2JK(FOCK,DMhf,IERI,ERI)
       FOCK = AHCORE + FOCK       
       EHFL = TRACE(DMhf,AHCORE,NBF) + TRACE(DMhf,FOCK,NBF) 
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -             
      else if(MSpin>0)then                    ! High-Spin State
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -            
       ALLOCATE (DMa(NBF,NBF))
       CALL DENMATHFr(DMhf,COEF)              ! 1,NB
       CALL DENMATHF05ro(DMa,COEF)            ! NB+1,NA       
!- - - - - - - - - - - - - - - - - - - - -       
       CALL FORM2JK(FOCK,DMhf,IERI,ERI)    
       EHFL=2.0*TRACE(DMhf+DMa,AHCORE,NBF)+TRACE(DMhf+2.0*DMa,FOCK,NBF)
       FOCK = AHCORE + FOCK       
       IF(NSOC>1)THEN
        ALLOCATE (FOCKa(NBF,NBF))       
        CALL FORMJK(FOCKa,DMa,IERI,ERI)
        EHFL = EHFL + 2.0*TRACE(DMa,FOCKa,NBF)
        FOCK = FOCK + FOCKa
        DEALLOCATE(FOCKa)                     
       ENDIF
       DEALLOCATE(DMa)             
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -             
      end if
      DEALLOCATE(DMhf)      
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     FOCKm: Fock Matrix in the molecular basis set
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE (AUX(NBF,NORB),TVEC(NORB,NBF))
      AUX   = MATMUL(FOCK,VEC)
      TVEC  = TRANSPOSE(VEC)     
      FOCKm = MATMUL(TVEC,AUX)
      DEALLOCATE(FOCK,AUX,TVEC)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     EIG: Energy Eigenvalues without Core
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IEIGeqELAG = 0
      IF(IEIGeqELAG==1)THEN
!- - - - - - - - - - - - - - - - - - - - 
       DO i=1,NORB
        EIG(i) = ELAG(i+NO1PT2,i+NO1PT2)
       ENDDO
!- - - - - - - - - - - - - - - - - - - - 
      ELSE
!- - - - - - - - - - - - - - - - - - - - 
!      EIG(i) -> FOCKm(i,i)
!- - - - - - - - - - - - - - - - - - - - 
       DO i=1,NORB
        EIG(i) = FOCKm(i,i)
       ENDDO
!- - - - - - - - - - - - - - - - - - - - 
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Print one-particle energies (EIG) if IPRINTEIG = 1
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IPRINTEIG = 0
      IF(IPRINTEIG==1)THEN
       write(6,'(/,A22,/)')'One-particle energies:'
       do i=1,NORB
        write(6,'(4x,I5,F20.10)')NO1PT2+i,EIG(i)
       enddo
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END

! DENMATHFr
      SUBROUTINE DENMATHFr(DM,CHF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::DM,CHF
!----------------------------------------------------------------------
      DO M=1,NBF
       DO N=M,NBF
        DM(M,N)=0.0d0
        DO J=1,NB                                   ! NB
         DM(M,N)=DM(M,N)+ CHF(M,J)*CHF(N,J)
        ENDDO
        DM(N,M)=DM(M,N)
       ENDDO
      ENDDO
      DM=2.0d0*DM
!----------------------------------------------------------------------
      RETURN
      END
      
! FORMJK
      SUBROUTINE FORMJK(FM,PM,IERI,ERI)
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
       NOPT=4
       CALL MPI_SEND(NOPT,1,MPI_INTEGER8,I,I,MPI_COMM_WORLD,IERR)
       CALL MPI_SEND(NBFT,1,MPI_INTEGER8,I,I,MPI_COMM_WORLD,IERR)
      ENDDO
#endif
!----------------------------------------------------------------------
      CALL SQUARETRIAN(PM,P,NBF,NBFT)
!----------------------------------------------------------------------
      F = 0.0d0
#ifdef MPI
      FF = 0.0d0
      CALL MPI_BCAST(NBF,1,MPI_INTEGER8,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(P,NBFT,MPI_REAL8,MASTER,MPI_COMM_WORLD,IERR)
#endif
      DO M=1,NINTCR
       LABEL = IERI(M)
       CALL LABELIJKL(LABEL,I,J,K,L)
!-----------------------------------------------------------------------
!      J
!-----------------------------------------------------------------------
       XJ = ERI(M)
       NIJ = I*(I-1)/2 + J
       NKL = K*(K-1)/2 + L

       IF(IHUB==0)CALL OTTOINTEGR(I,J,K,L,NIJ,NKL,XJ)

                      F(NIJ)=F(NIJ)+0.5*P(NKL)*XJ
       IF(NIJ/=NKL)   F(NKL)=F(NKL)+0.5*P(NIJ)*XJ
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

! DENMATHF05ro
      SUBROUTINE DENMATHF05ro(DMa,CHF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::DMa,CHF
!----------------------------------------------------------------------
      DO M=1,NBF
       DO N=M,NBF
       DMa(M,N)=0.0d0
       DO J=NB+1,NA
        DMa(M,N)=DMa(M,N)+ CHF(M,J)*CHF(N,J)
       ENDDO
       DMa(N,M)=DMa(M,N)
       ENDDO
      ENDDO
!----------------------------------------------------------------------
      RETURN
      END
                  
! ELAGCOEF0
      SUBROUTINE ELAGCOEF0(ELAG,COEF,RO,CJ12,CK12,AHCORE,               &
                           ADIPx,ADIPy,ADIPz,IERI,ERI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      LOGICAL HighSpin
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
      COMMON/INPFILE_NIJKL/NINTMX,NIJKL,NINTCR,NSTORE
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
!
      INTEGER,DIMENSION(NIJKL) :: IERI
      DOUBLE PRECISION,DIMENSION(NIJKL) :: ERI
      DOUBLE PRECISION,DIMENSION(NBF5) :: RO
      DOUBLE PRECISION,DIMENSION(NBF,NBF) :: ELAG,COEF,AHCORE
      DOUBLE PRECISION,DIMENSION(NBF,NBF) :: ADIPx,ADIPy,ADIPz
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5):: CJ12,CK12
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: WJj,WKj,WF,G,AUX1
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:):: AUX2
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Calculate AUX1: QD(j,miu,niu), Jj(miu,niu), Kj(miu,niu) (j=1,NBF5) 
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE (WJj(NSQ,NBF5),WKj(NSQ,NBF5),AUX1(NBF,NBF),AUX2(NSQ))
      DO j=1,NBF5
       CALL DENMATj(j,AUX1,COEF,NBF)
       CALL HSTARJ(AUX2,AUX1,IERI,ERI)
       WJj(1:NSQ,j) = AUX2(1:NSQ)
       CALL HSTARK(AUX2,AUX1,IERI,ERI)
       WKj(1:NSQ,j) = AUX2(1:NSQ)
      ENDDO
      DEALLOCATE (AUX1,AUX2)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Form F Matrix and keep it in WF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE (WF(NSQ,NBF5))
      if(MSpin==0)then               ! Singlet and Multiplet States
       CALL FFJMN1rc(RO,CJ12,CK12,AHCORE,WJj,WKj,WF,ADIPx,ADIPy,ADIPz)
      else if(MSpin>0)then           ! High-Spin States
       CALL FFJMN1ro(RO,CJ12,CK12,AHCORE,WJj,WKj,WF,ADIPx,ADIPy,ADIPz)      
      end if      
      DEALLOCATE (WJj,WKj)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Calculate G Matrix
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE (G(NBF,NBF5))
      DO IQ=1,NBF5
       do i=1,nbf
        G(i,IQ) = FC(i,IQ,WF(1,IQ),COEF)
       enddo
      ENDDO
      DEALLOCATE (WF)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Lagrangian Multipliers (ELAG)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL ELG(ELAG,COEF,G)
      DEALLOCATE (G)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END      
