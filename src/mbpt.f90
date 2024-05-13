!======================================================================!
!                                                                      !
!              N O F - M B P T  S U B R O U T I N E S                  !
!                                                                      !
!            ( J. Chem. Theory Comput. 17, 7562, 2021 )                !
!                                                                      !
!======================================================================!
!                                                                      !
!              Developed    Date: 10/01/2021                           !
!              by M. Rodríguez-Mayorga and I. Mitxelena                !
!              Adapted      Date: 30/10/2022                           !
!              by M. Rodríguez-Mayorga                                 !
!                                                                      !
!    mbpt.f90 contains the subroutines related to the RPA(GW)          ! 
!             calculations. RPA is trapezoidal rule of GW@GM           !
!                                                                      !
! ==================================================================== !

! MBPT (RPA) (NOF-c-X)
      SUBROUTINE MBPT_RPA(ELAG,COEF,RO,CJ12,CK12,AHCORE,IERI,ERI,      &
                           ADIPx,ADIPy,ADIPz)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL HighSpin
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_ORBSPACE0/NO1,NDOC,NCO,NCWO,NVIR,NAC,NO0
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
      COMMON/INPNOF_GENERALINF/ICOEF,MAXIT,MAXIT21
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
!-----------------------------------------------------------------------
!     Calculation RPA corr energy
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL ERPA_MBPT(ECd,EHFL,ELAG,COEF,RO,CJ12,CK12,AHCORE,ADIPx,ADIPy,ADIPz,IERI,ERI)
      WRITE(6,1)ECd,ECndl,ECd+ECndl,EHFL+ECd+ECndl+EN+ECndHF
    1 FORMAT(/3X,'           ECd =',F20.10,/,                           &
              3X,'          ECnd =',F20.10,/,                           &
              3X,'         ECorr =',F20.10,/,                           &
              3X,'     E(NOFRPA) =',F20.10)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END

! ERPA      
     SUBROUTINE ERPA_MBPT(ECd,EHFL,ELAG,COEF,RO,CJ12,CK12,AHCORE,ADIPx,ADIPy,ADIPz,IERI,ERI)
     COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
     COMMON/INPNOF_ORBSPACE0/NO1,NDOC,NCO,NCWO,NVIR,NAC,NO0
     COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
     COMMON/INPNOF_GENERALINF/ICOEF,MAXIT,MAXIT21
     COMMON/EHFEN/EHF,EN
     COMMON/INPFILE_NO1PT2/NO1PT2,NEX
     COMMON/INPFILE_NIJKL/NINTMX,NIJKL,NINTCR,NSTORE
     COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
     COMMON/CorrNonDynamic/ECnd,ECndl,ECndHF
     COMMON/INPNOF_Tijab/NOUTTijab,NTHRESHTijab,THRESHTijab
     COMMON/NumLinIndOrb/NQMT
     COMMON/ENERGIAS/EELEC,EELEC_OLD,DIF_EELEC,EELEC_MIN
     PARAMETER (ZERO=0.0d0)
     PARAMETER (FOURTH=0.25d0)
     PARAMETER (ONE=1.0d0)
     PARAMETER (TWO=2.0d0)
     PARAMETER (FOUR=4.0d0)
     LOGICAL::diagFOCK 
     INTEGER::i,j,k,l,a,b,info,last_coup
     INTEGER::order,Nab
     INTEGER,DIMENSION(NIJKL)::IERI
     DOUBLE PRECISION,DIMENSION(NIJKL)::ERI
     DOUBLE PRECISION::ECd,EPNOF,ESDc,EcRPA,EcMP2,EcGoWo
     DOUBLE PRECISION::mu,EHFL 
     DOUBLE PRECISION,DIMENSION(NBF5)::RO
     DOUBLE PRECISION,DIMENSION(NBF5,NBF5)::CJ12,CK12
     DOUBLE PRECISION,DIMENSION(NBF,NBF)::AHCORE,ADIPx,ADIPy,ADIPz
     DOUBLE PRECISION,DIMENSION(NBF,NBF)::ELAG,COEF
     DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::OCC,EIG,BIGOMEGA,TEMPV
     DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::CINTER,CINTRA 
     DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::TEMPM,FOCKm,COEF2
     DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::ApB,AmB
     DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::XmY,XpY 
     DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:,:)::ERImol,ERImol2
     DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:)::wmn
!----------------------------------------------------------------------
!    NCO:  Number of HF occupied MOs (OCC=1 in SD)
!    NVIR: Number of HF virtual  MOs (OCC=0 in SD)
!
!    NO1:  Number of inactive doubly occupied orbitals (OCC=1)
!    NDOC: Number of strongly occupied MOs
!    NCWO: Number of coupled weakly occupied MOs per strongly occupied
!    NCWO*NDOC: Active orbitals in the virtual subspace
!    NO0: Empty orbitals  (OCC=0)
!    NAC:  Dimension of the active natural orbital subspace
!
!    NB/NCO: Number of orbitals in the doubly-occupied space
!    NA: NB/NCO + Number of orbitals in the singly-occupied space (NSOC)
!
!          NCO + NSOC     |       NVIR          = NBF
!----------------------------------------------------------------------
     IF(ICOEF==0)CALL ELAGCOEF0(ELAG,COEF,RO,CJ12,CK12,AHCORE,ADIPx,ADIPy,ADIPz,IERI,ERI)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Initial allocation and def. of number of pairs (occ,virtual) i <-> a  
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     NVIR=NBF-NA
     Nab=NVIR*NA
     last_coup=NA+NCWO*(NCO-NO1PT2)
     WRITE(6,1)NBF,NO1PT2,NA,NVIR,last_coup,Nab
     ALLOCATE(OCC(NBF),TEMPM(NBF,NBF),COEF2(NBF,NBF))
     COEF2=COEF
     OCC=ZERO
     DO i=1,NBF5
      OCC(i) = RO(i)
     ENDDO
     DO j=1,NBF
      do i=1,NBF
       TEMPM(i,j) = COEF(i,j)
      enddo
     ENDDO
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! One orbital energies (EIG) and FOCK matrix (FOCKm) with COEF/=CHF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     ALLOCATE(EIG(NBF),FOCKm(NBF,NBF))
     CALL FOCKMOL(NBF,COEF,TEMPM,ELAG,EIG,FOCKm,AHCORE,IERI,ERI,ESD)
     EHFL=ESD
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Allocate 2e_integrals array(s) and form ERI in MO basis (see mp2.f)
! Prepare coef. factors CINTRA and CINTER.
! Tune Fock(i,a) and ERImol elements before computing W, wmn, etc. 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     ALLOCATE(ERImol(NBF,NBF,NBF,NBF))
     CALL ERIC1c(ERImol,IERI,ERI,TEMPM,NBF)
     CALL ERIC23c(ERImol,TEMPM,NBF)
     CALL ERIC4c(ERImol,TEMPM,NBF)
!     call print2eint(NBF,ERImol)
     ALLOCATE(CINTER(NBF),CINTRA(NBF),ERImol2(NBF,NBF,NBF,NBF))
     ERImol2=ERImol
     call tune_coefs(CINTER,CINTRA,OCC,NBF,NA)
     call tunefock(CINTRA,CINTER,FOCKm,NBF,NCO,NVIR,NCWO,NO1PT2,NA)
     call tuneerimol(CINTER,CINTRA,NBF,NCO,NVIR,NCWO,NO1PT2,ERImol,NA)
     DEALLOCATE(CINTER,CINTRA)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Transform FOCKm and ERIs from NO basis to "DIAG(FOCK)" basis (TEMPM)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     TEMPM=ZERO 
     diagFOCK=.true.
     call diafock(FOCKm,NBF,EIG,TEMPM,diagFOCK)
     DEALLOCATE(FOCKm)
     if(.not.diagFOCK) then
      COEF2=matmul(COEF,TEMPM)
      call transformERI(NBF,TEMPM,ERImol)
      call transformERI(NBF,TEMPM,ERImol2)
     endif
     DEALLOCATE(TEMPM)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Print orb. energies and occ numbers.
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     write(*,*) '          '
     write(*,*) ' LIST OF ORBITAL ENERGIES (A.U.) AND OCC NUMBERS IN USE'
     write(*,*) '          '
     mu=(EIG(NA+1)+EIG(NA))/TWO
     DO i=1,NBF
      EIG(i)=EIG(i)-mu
      write(6,2) EIG(i),OCC(i)
     ENDDO
     write(*,*) '          '
     write(*,'(a,f10.5)') ' COMPUTED CHEMICAL POTENTIAL',mu
     write(*,*) '          '
     DEALLOCATE(OCC) ! We do not need the OCCs anymore
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Build A+B and A-B matrices for CASIDA EQN and 1st contribution to EcRPA.
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     EcRPA=ZERO
     ALLOCATE(AmB(Nab,Nab),ApB(Nab,Nab))
     AmB=ZERO
     ApB=ZERO
     k=1      
     do i=1,NA
       do a=NA+1,nbf
         l=1
         do j=1,NA
          do b=NA+1,nbf
            ApB(k,l)=FOUR*ERImol(i,j,b,a) ! This is i-a and j-b, the 4 is from the spin sum 
            AmB(k,l)=ZERO 
           if(i.eq.j.and.a.eq.b) then      
            AmB(k,l)=EIG(a)-EIG(i)
            ApB(k,l)=ApB(k,l)+AmB(k,l)
           endif 
           if(k.eq.l) then ! 1st contrib. to RPA -> 0.25*diag_terms  
            EcRPA=EcRPA-FOURTH*(ApB(k,k)+AmB(k,k)) 
           endif  
           l=l+1
          enddo
         enddo
        k=k+1
       enddo
     enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Solve Casida Equation by Cholesky Decomposition (calling LAPACK):
! sqrt(BIGOMEGA)=Excitation energies
! the eigenvectors (F) of C*F=BIGOMEGA^2 *F are stored in columns in TEMPM
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     info=0
     CALL DPOTRF('L',Nab,ApB,Nab,info)
     if(info.ne.0) then
      write(*,*) ' ' 
      write(*,*) ' WARNING! DPOTRF FAILED TO CONVERGE ',info 
      write(*,*) ' '  
     else
      CALL DSYGST(3,'L',Nab,AmB,Nab,ApB,Nab,info)
      if(info.ne.0) then 
       write(*,*) ' ' 
       write(*,*) ' WARNING! DSYGST FAILED TO CONVERGE ',info
       write(*,*) ' '
      endif
     endif
     ALLOCATE(TEMPM(Nab,Nab),BIGOMEGA(Nab),TEMPV(Nab))
     ALLOCATE(XpY(Nab,Nab),XmY(Nab,Nab))
     do i=1,Nab
      do j=i+1,Nab
        AmB(i,j)=AmB(j,i)
      enddo
     enddo
     BIGOMEGA=ZERO
     TEMPV=ZERO 
     CALL DIAG(Nab,AmB,TEMPM,BIGOMEGA,TEMPV) 
     BIGOMEGA(:) = SQRT(BIGOMEGA(:))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Compute eigenvectors X+Y and X-Y as in MOLGW (calling LAPACK)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     XpY=TEMPM
     XmY=TEMPM
     if(info.eq.0) then
      call DTRMM('L','L','N','N',Nab,Nab,ONE,ApB,Nab,XmY,Nab)
      call DTRSM('L','L','T','N',Nab,Nab,ONE,ApB,Nab,XpY,Nab)
     endif
     forall(i=1:Nab)
      XpY(:,i)=XpY(:,i)*SQRT(BIGOMEGA(i))
      XmY(:,i)=XmY(:,i)/SQRT(BIGOMEGA(i))
     endforall
     DEALLOCATE(AmB,ApB,TEMPM,TEMPV,XmY) ! We do not need X-Y anymore 
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Compute omegas, oscillator strenghts and static polarizability (omega->0)       
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     call td_polarizability(NBF,NCO,Nab,NA,COEF2,XpY,BIGOMEGA,ADIPx,ADIPy,ADIPz,EcRPA)
     DEALLOCATE(COEF2)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Our RPA (trapezoidal rule) is 1/2 GoWo@Galitskii-Migdal EQN
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Compute w^s mn = sum_ia <im|an>* (X^s _ia + Y^s_ia) -> see MolGW paper (i is occ, a is virt) 
     ALLOCATE(wmn(NBF,NBF,Nab))
     write(*,*) ' '
     write(*,*) 'COMPUTING RPA CORRELATION ENERGY FOR NOF-c-RPA'
     wmn=ZERO; 
     call build_wmn(NBF,Nab,NA,wmn,ERImol,XpY) ! Building only w(i,a,s) and w(a,i,s) like terms
     EcGoWo=ZERO
     call gw_gm_eq(NA,NCO,NBF,Nab,wmn,EIG,EcGoWo,XpY,BIGOMEGA,ERImol2)
     ECd=EcGoWo
     DEALLOCATE(XpY,wmn,BIGOMEGA)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! MP2 Ec energy 
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     write(*,*) 'COMPUTING MP2 CORRELATION ENERGY FOR NOF-c-MP2'
     EcMP2=ZERO
     call mp2_eq(NA,NCO,NBF,EIG,ERImol,ERImol2,EcMP2)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Write Ec energies and final results
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     ESD=ESD+EN+ECndHF
     ESDc=ESD+ECndl
     EPNOF=EELEC+EN
     write(6,3)ESD,ESDc,EPNOF,ECndl,EcRPA,ECd,EcMP2,&
    & ESDc+EcRPA,ESDc+ECd,ESDc+EcMP2
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Deallocate and clean 
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     DEALLOCATE(ERImol2)
     DEALLOCATE(EIG,ERImol)
     RETURN
!---------------------------------------------------------------------
   1 FORMAT(2X,//,3X,'MBPT ENERGIES',/,                       &
           2X,15('='),//,                                    &
           1X,'NUMBER OF ORBITALS        (NBASIS) =',I5,/,   &
           1X,'NUMBER OF FROZEN PAIRS       (NFR) =',I5,/,   &
           1X,'NUMBER OF OCC ORBITALS       (NOC) =',I5,/,   &
           1X,'NUMBER OF VIR ORBITALS      (NVIR) =',I5,/,   &
           1X,'LAST COUPLED ORBITAL        (NLAS) =',I5,/,   &
           1X,'SIZE OF A+B AND A-B (NAB=NOCxNVIR) =',I5) 
   2 FORMAT(3X,F15.10,' ',F15.10,' ')
   3 FORMAT(/,1X,'E(HFL)                ',5X,F20.10,' a.u.',/, &
     1X,'E(HFL+ND)             ',5X,F20.10,' a.u.',/,         &
     1X,'E(PNOFi)              ',5X,F20.10,' a.u.',/,         &
     ' ',/,                                                   &
     1X,'Ec(ND)                ',5X,F20.10,' a.u.',/,         &
     1X,'Ec(RPA-FURCHE)        ',5X,F20.10,' a.u.',/,         &
     1X,'Ec(RPA)               ',5X,F20.10,' a.u.',/,         &
     1X,'Ec(MP2)               ',5X,F20.10,' a.u.',/,         &
     ' ',/,                                                   & 
     1X,'E(NOF-c-RPA-FURCHE)   ',5X,F20.10,' a.u.',/,         &
     1X,'E(NOF-c-RPA)          ',5X,F20.10,' a.u.',/,         &
     1X,'E(NOF-c-MP2)          ',5X,F20.10,' a.u.')
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     END SUBROUTINE ERPA_MBPT

! tune_coefs     
     subroutine tune_coefs(CINTER,CINTRA,OCC,NBF,NA)
     implicit none
     integer,intent(in)::NBF,NA
     double precision,dimension(nbf),intent(in)::OCC
     double precision,dimension(nbf),intent(inout)::CINTER,CINTRA
     integer::ii
     CINTER=1.0d0
     CINTRA=1.0d0
     do ii=1,NBF
      ! CINTRA 
      CINTRA(ii)=1.0d0-DABS(1.0d0-2.0d0*OCC(ii))
      CINTRA(ii)=1.0d0-CINTRA(ii)**2.0d0
      ! CINTER
      if(ii.le.NA) then ! eq. (34) PRA multiplets
        CINTER(ii)=1.0d0
      else
       CINTER(ii)=DABS(1.0d0-2.0d0*OCC(ii))
       CINTER(ii)=CINTER(ii)**2.0d0
      endif
     enddo
     write(*,*) ' '
     end subroutine tune_coefs

! tunefock     
     subroutine tunefock(CINTRA,CINTER,FOCKm,NBF,NCO,NVIR,NCWO,nfr,NA)
     implicit none
     integer,intent(in)::NBF,NCO,NVIR,NCWO,NFR,NA
     double precision,dimension(NBF),intent(inout)::CINTER,CINTRA
     double precision,dimension(NBF,NBF),intent(inout)::FOCKm
     integer,allocatable,dimension(:)::coup
     integer::i,k,in,kn,lmin_i,lmax_i,last_coup
     allocate(coup(NVIR))
     coup=0
     last_coup=NA+ncwo*(nco-nfr)
     do i=1,nco
      lmin_i = NA+ncwo*(nco-i)+1
      lmax_i = NA+ncwo*(nco-i)+ncwo
      do k=1,nvir
       kn=k+NA
       FOCKm(i,kn)=0.0d0 
       FOCKm(kn,i)=0.0d0 
       if((lmin_i<=kn.and.kn<=lmax_i).and.(lmax_i<=last_coup)) then
        coup(k)=i 
       endif
      enddo
     enddo
     do i=nco+1,NA
      do k=1,nvir
       kn=k+NA
       FOCKm(i,kn)=0.0d0
       FOCKm(kn,i)=0.0d0
      enddo
     enddo
     do i=1,NA
      do k=i+1,NA
       FOCKm(i,k)=FOCKm(i,k)*CINTER(i)*CINTER(k)
       FOCKm(k,i)=FOCKm(i,k)
      enddo
     enddo
     do i=1,nvir
      in=i+NA
      do k=i+1,nvir
       kn=k+NA
       if(coup(i)==coup(k)) then
        FOCKm(in,kn)=FOCKm(in,kn)*CINTRA(in)*CINTRA(kn)
       else
        FOCKm(in,kn)=FOCKm(in,kn)*CINTER(in)*CINTER(kn)
       endif
       FOCKm(kn,in)=FOCKm(in,kn)
      enddo
     enddo 
     deallocate(coup)
     end subroutine tunefock

! diafock     
     subroutine diafock(FOCKm,NBF,EIG,EIGENVE,diagF)
     implicit none
     logical,intent(inout)::diagF
     integer,intent(in)::NBF
     double precision,dimension(NBF),intent(inout)::EIG
     double precision,dimension(NBF,NBF),intent(inout)::FOCKm,EIGENVE
     double precision,allocatable,dimension(:)::TEMPV
     double precision::tol6
     integer::i,k
     tol6=1.0d-6
     do i=1,NBF-1
      do k=i+1,NBF
       if(abs(FOCKm(i,k))>tol6) then
        diagF=.false.
        goto 359
       endif
      enddo
     enddo
 359 continue
     if(diagF.eqv..false.) then
      allocate(TEMPV(NBF))
      CALL DIAG(NBF,FOCKm,EIGENVE,EIG,TEMPV) 
      deallocate(TEMPV)
     else
      EIGENVE=0.0d0
      do i=1,nbf
       EIGENVE(i,i)=1.0d0
      enddo
     endif
     end subroutine diafock

! tuneerimol     
     subroutine tuneerimol(CINTER,CINTRA,nbf,nco,nvir,ncwo,nfr,ERImol,NA)
     implicit none
     integer,intent(in)::NBF,NCO,NVIR,NCWO,NFR,NA
     double precision,dimension(NBF),intent(in)::CINTER,CINTRA
     double precision,dimension(nbf,nbf,nbf,nbf),intent(inout)::ERImol
     integer,dimension(:),allocatable::coup
     double precision,dimension(:,:,:,:),allocatable::ERImolTMP
     integer::i,j,k,l
     integer::a,b,c,d,an,bn,cn,dn,bmin_i,bmax_i,last_coup
     double precision::value1,Ciiab,Cijab,Cijkl,Cijka,Ciiia,Cabci,Cabcd
     allocate(ERImolTMP(NBF,NBF,NBF,NBF),coup(NA+1:NVIR+NA))
     ERImolTMP=0.0d0;coup(NA+1:NVIR+NA)=0
     last_coup=NA+ncwo*(nco-nfr)
     ! find the coupling
     do i=1,NA
      bmin_i=NA+ncwo*(nco-i)+1
      bmax_i=NA+ncwo*(nco-i)+ncwo
      do a=1,nvir
       an=a+NA
       if(      (bmin_i<=an.and.an<=bmax_i)  &
    &       .and. (bmax_i<=last_coup)        )then
        coup(an)=i
       endif
      enddo 
      ! <ov|o'v'> and <oo'|vv'> and equivalent terms
      do j=1,NA
       if(i==j .and. i<=nco) then ! Subspace 'i' occ, still could be inter or intra    
        do a=1,nvir
         an=a+NA
         do b=1,nvir
          bn=b+NA
          if((     (bmin_i<=an.and.an<=bmax_i)  &
    &        .and. (bmin_i<=bn.and.bn<=bmax_i)) & 
    &        .and. (bmax_i<=last_coup)        )then
           Ciiab=CINTRA(an)*CINTRA(bn)*CINTRA(i)*CINTRA(i)
          else
           Ciiab=CINTER(an)*CINTER(bn)*CINTER(i)*CINTER(i)
          endif
          value1=Ciiab*ERImol(i,i,bn,an)
          ERImolTMP(i,i,bn,an)=value1
          ERImolTMP(an,i,bn,i)=value1
          ERImolTMP(an,bn,i,i)=value1
          ERImolTMP(i,bn,i,an)=value1
          ERImolTMP(i,i,an,bn)=value1
          ERImolTMP(bn,i,an,i)=value1
          ERImolTMP(bn,an,i,i)=value1
          ERImolTMP(i,an,i,bn)=value1
          value1=Ciiab*ERImol(i,an,bn,i)
          ERImolTMP(i,an,bn,i)=value1
          ERImolTMP(i,bn,an,i)=value1
          ERImolTMP(bn,i,i,an)=value1
          ERImolTMP(an,i,i,bn)=value1
         enddo
        enddo
       else         ! Purely interspace  
        do a=1,nvir
         an=a+NA
         do b=1,nvir
          bn=b+NA
          Cijab=CINTER(an)*CINTER(bn)*CINTER(i)*CINTER(j)
          value1=Cijab*ERImol(i,j,bn,an)
          ERImolTMP(i,j,bn,an)=value1
          ERImolTMP(i,bn,j,an)=value1
          ERImolTMP(an,bn,j,i)=value1
          ERImolTMP(an,j,bn,i)=value1
          ERImolTMP(j,i,an,bn)=value1
          ERImolTMP(j,an,i,bn)=value1
          ERImolTMP(bn,an,i,j)=value1
          ERImolTMP(bn,i,an,j)=value1
          value1=Cijab*ERImol(i,an,bn,j)
          ERImolTMP(i,an,bn,j)=value1
          ERImolTMP(i,bn,an,j)=value1
          ERImolTMP(j,bn,an,i)=value1
          ERImolTMP(j,an,bn,i)=value1
          ERImolTMP(an,j,i,bn)=value1
          ERImolTMP(an,i,j,bn)=value1
          ERImolTMP(bn,i,j,an)=value1
          ERImolTMP(bn,j,i,an)=value1
         enddo
        enddo
       endif
      enddo
     enddo
     ERImol=ERImolTMP
     deallocate(ERImolTMP,coup)
     end subroutine tuneerimol      

! td_polarizability     
     subroutine td_polarizability(NBF,NCO,Nab,NA,COEF,XpY,BIGOMEGA,ADIPx,ADIPy,ADIPz,EcRPA)
     implicit none
     INTEGER,intent(in)::Nbf,Nab,NCO,NA
     DOUBLE PRECISION,intent(inout)::EcRPA
     DOUBLE PRECISION,DIMENSION(NBF,NBF),intent(in)::ADIPx,ADIPy,ADIPz
     DOUBLE PRECISION,DIMENSION(NBF,NBF),intent(in)::COEF
     DOUBLE PRECISION,DIMENSION(Nab),intent(in)::BIGOMEGA
     DOUBLE PRECISION,DIMENSION(Nab,Nab),intent(in)::XpY
     INTEGER::i,j,k,l,a
     DOUBLE PRECISION::tol10,tol6,AUtoEV
     DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) ::MDIPx,MDIPy,MDIPz
     DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) ::TEMPM
     DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) ::STATICPOL
     DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) ::OSCSTR,DIPSUM
     tol6=1.0d-6
     tol10=1.0d-10
     AUtoEV=27.211399
! Prepare dipole moment from AO to MO
     ALLOCATE(MDIPx(NBF,NBF),MDIPy(NBF,NBF),MDIPz(NBF,NBF))
     ALLOCATE(TEMPM(NBF,NBF))
!     x
     TEMPM=matmul(ADIPx,COEF) 
     MDIPx=matmul(transpose(COEF),TEMPM)
!     y
     TEMPM=matmul(ADIPy,COEF) 
     MDIPy=matmul(transpose(COEF),TEMPM)
!     z
     TEMPM=matmul(ADIPz,COEF) 
     MDIPz=matmul(transpose(COEF),TEMPM)
     DEALLOCATE(TEMPM)
! Compute OSCILLATOR STRENGHTS and 2nd contribution to EcRPA
     ALLOCATE(OSCSTR(Nab),DIPSUM(3),TEMPM(3,Nab))
     OSCSTR=0.0d0
     do k=1,Nab
      EcRPA=EcRPA+0.5d0*ABS(BIGOMEGA(k))  ! 2nd contribution to EcRPA
      DIPSUM=0.d0
      l=1
      do i=1,NA
       do a=NA+1,NBF
         DIPSUM(1)=DIPSUM(1)+MDIPx(i,a)*XpY(l,k)*DSQRT(2.0d0)
         DIPSUM(2)=DIPSUM(2)+MDIPy(i,a)*XpY(l,k)*DSQRT(2.0d0)
         DIPSUM(3)=DIPSUM(3)+MDIPz(i,a)*XpY(l,k)*DSQRT(2.0d0)
         l=l+1
       enddo
      enddo
      TEMPM(1,k)=DIPSUM(1)
      TEMPM(2,k)=DIPSUM(2)
      TEMPM(3,k)=DIPSUM(3)
      OSCSTR(k)=2.0d0*(DIPSUM(1)**2.0d0+DIPSUM(2)**2.0d0+DIPSUM(3)**2.0d0)*BIGOMEGA(k)/3.0d0
     enddo
! Print Omegas and Oscillator strenghts
     write(*,*) ' '
     write(*,*) 'TD-H (RPA) CASIDA EQ. SOLVED' 
     write(*,*) 'N. EXCITATION   A.U.         eV            nm      OSC. STRENGTH'
     do i=1,Nab
      if(OSCSTR(i)>tol6) then
       write(6,13) i,BIGOMEGA(i),BIGOMEGA(i)*AUtoEV,1239.84193/(BIGOMEGA(i)*AUtoEV),OSCSTR(i)
      endif
     enddo
     DEALLOCATE(OSCSTR,MDIPx,MDIPy,MDIPz,DIPSUM)
! Static polarizability alpha_xx' = <x|Xi(r,r',w=0)|x'> 
     ALLOCATE(STATICPOL(3,3))
     STATICPOL=0.0d0
     do i=1,Nab
      forall(j=1:3, k=1:3)
           STATICPOL(j,k)=STATICPOL(j,k)+2.0d0*TEMPM(j,i)*TEMPM(k,i)/(BIGOMEGA(i)+tol10)
      end forall
     enddo
     write(*,*) ' '
     write(6,14) STATICPOL(1,1),STATICPOL(1,2),STATICPOL(1,3),& 
    &   STATICPOL(2,1),STATICPOL(2,2),STATICPOL(2,3),         &
    &   STATICPOL(3,1),STATICPOL(3,2),STATICPOL(3,3)         
     write(6,15) ((STATICPOL(1,1)+STATICPOL(2,2)+STATICPOL(3,3))/3.0d0)
     DEALLOCATE(TEMPM,STATICPOL)

  13 FORMAT(3X,I5,3X,F12.5,1X,F12.5,1X,F12.4,1X,F12.6)
  14 FORMAT(1X,'STATIC POLARIZABILITY:',//,  &
           F15.10,3X,F15.10,3X,F15.10,/,     &
           F15.10,3X,F15.10,3X,F15.10,/,     &
           F15.10,3X,F15.10,3X,F15.10)
  15 FORMAT(/,1X,'TRACE OF THE STATIC POLARIZABILITY',F15.10)

     end subroutine td_polarizability

! build_wmn     
     subroutine build_wmn(NBF,Nab,NA,wmn,ERImol,XpY)
     implicit none
     integer,intent(in)::NBF,Nab,NA
     double precision,dimension(NBF,NBF,Nab),intent(inout)::wmn 
     double precision,dimension(Nab,Nab),intent(in)::XpY
     double precision,dimension(NBF,NBF,NBF,NBF),intent(in)::ERImol
     integer::s,ss,p,q,i,a
     do s=1,Nab
      ! Building only w(i,a,s) occ-vir 
      !do p=1,NBF 
      do p=NA,NBF
       !do q=1,p
       do q=1,p-1
        i=1
        a=NA+1
        do ss=1,Nab
         wmn(p,q,s)=wmn(p,q,s)+ERImol(i,p,q,a)*XpY(ss,s)
         a=a+1
         if(a>NBF) then
          a=NA+1
          i=i+1
         endif
        enddo
        wmn(p,q,s)=wmn(p,q,s)*DSQRT(2.0d0)
        if(p/=q) then
         wmn(q,p,s)=wmn(p,q,s)
        endif
       enddo
      enddo
     enddo
     end subroutine build_wmn

! gw_gm_eq     
!    Loop to compute the EcGoWo energy. Recall that Go is used in GM EQ.
!    Very symple Eq: Ec = 2 sum _i sum_a sum_s 
!                [ (wia)^s ]**2 /(e(i)-e(a)-Omega(s))
     subroutine gw_gm_eq(NA,NCO,NBF,Nab,wmn,EIG,EcGoWo,XpY,BIGOMEGA,ERImol)
     implicit none
     integer,intent(in)::NCO,NBF,Nab,NA
     double precision,intent(inout)::EcGoWo
     double precision,dimension(Nab),intent(in)::BIGOMEGA
     double precision,dimension(NBF),intent(in)::EIG
     double precision,dimension(Nab,Nab),intent(in)::XpY
     double precision,dimension(NBF,NBF,Nab),intent(in)::wmn
     double precision,dimension(NBF,NBF,NBF,NBF),intent(in)::ERImol
     integer::i,j,a,b,s,l,fst_virt
     double precision::tol10
     tol10=1.0d-10
     fst_virt=NA+1
     ! Doubly occ part (4 index version)
      do a=NA+1,NBF
       do b=NA+1,NBF
        do i=1,NCO
         do j=1,NCO
          l=(j-1)*(NBF-NA)+(b-fst_virt+1)
          do s=1,Nab 
           EcGoWo=EcGoWo+wmn(i,a,s)*ERImol(j,i,a,b)*XpY(l,s)*DSQRT(2.0d0)/(EIG(i)-EIG(a)-BIGOMEGA(s)+tol10)
          enddo
         enddo
        enddo
       enddo
      enddo
     ! Open-shell part
     if(NA/=NCO) then  
      do a=NA+1,NBF
       do b=NA+1,NBF
        do i=1,NCO
         do j=NCO+1,NA
          l=(j-1)*(NBF-NA)+(b-fst_virt+1)
          do s=1,Nab 
           EcGoWo=EcGoWo+0.5d0*wmn(i,a,s)*ERImol(j,i,a,b)*XpY(l,s)*DSQRT(2.0d0)/(EIG(i)-EIG(a)-BIGOMEGA(s)+tol10)
          enddo
         enddo
        enddo
        do i=NCO+1,NA
         do j=1,NCO
          l=(j-1)*(NBF-NA)+(b-fst_virt+1)
          do s=1,Nab 
           EcGoWo=EcGoWo+0.5d0*wmn(i,a,s)*ERImol(j,i,a,b)*XpY(l,s)*DSQRT(2.0d0)/(EIG(i)-EIG(a)-BIGOMEGA(s)+tol10)
          enddo
         enddo
         do j=NCO+1,NA
          if(i/=j) then
           l=(j-1)*(NBF-NA)+(b-fst_virt+1)
           do s=1,Nab 
            EcGoWo=EcGoWo+0.25d0*wmn(i,a,s)*ERImol(j,i,a,b)*XpY(l,s)*DSQRT(2.0d0)/(EIG(i)-EIG(a)-BIGOMEGA(s)+tol10)
           enddo
          endif
         enddo
        enddo
       enddo
      enddo
     endif
     end subroutine gw_gm_eq

! mp2_eq     
     subroutine mp2_eq(NA,NCO,NBF,EIG,ERImol,ERImol2,EcMP2)
     implicit none
     integer,intent(in)::NCO,NBF,NA
     double precision,intent(inout)::EcMP2
     double precision,dimension(NBF),intent(in)::EIG
     double precision,dimension(NBF,NBF,NBF,NBF),intent(in)::ERImol
     double precision,dimension(NBF,NBF,NBF,NBF),intent(in)::ERImol2
     integer::i,j,a,b
     double precision::tol10
     tol10=1.0d-10
     do a=NA+1,NBF
      do b=NA+1,NBF
       do i=1,NCO
        do j=1,NCO
         ! Use [2Tij,ab - Tij,ba] as 2-RDM element (the tuned one in Piris PRA).
         EcMP2=EcMP2+ERImol2(a,b,j,i)*(2.0d0*ERImol(a,b,j,i)-ERImol(a,b,i,j))/(EIG(i)+EIG(j)-EIG(a)-EIG(b)+tol10)
        enddo
        do j=NCO+1,NA
         ! Use [2Tij,ab - Tij,ba] as 2-RDM element (the tuned one in Piris PRA).
         EcMP2=EcMP2+ERImol2(a,b,j,i)*(ERImol(a,b,j,i)-0.5d0*ERImol(a,b,i,j))/(EIG(i)+EIG(j)-EIG(a)-EIG(b)+tol10)
        enddo
       enddo
       do i=NCO+1,NA
        do j=1,NCO
         ! Use [2Tij,ab - Tij,ba] as 2-RDM element (the tuned one in Piris PRA).
         EcMP2=EcMP2+ERImol2(a,b,j,i)*(ERImol(a,b,j,i)-0.5d0*ERImol(a,b,i,j))/(EIG(i)+EIG(j)-EIG(a)-EIG(b)+tol10)
        enddo
        do j=NCO+1,NA
         if(j/=i) then
          ! Use [2Tij,ab - Tij,ba] as 2-RDM element (the tuned one in Piris PRA).
          EcMP2=EcMP2+0.5d0*ERImol2(a,b,j,i)*(ERImol(a,b,j,i)-0.5d0*ERImol(a,b,i,j))/(EIG(i)+EIG(j)-EIG(a)-EIG(b)+tol10)
         endif
        enddo
       enddo
      enddo
     enddo 
     end subroutine mp2_eq

! transformERI     
     subroutine transformERI(NBF,TEMPM,ERImol)
     implicit none
     integer,intent(in)::NBF
     double precision,dimension(nbf,nbf,nbf,nbf),intent(inout)::ERImol
     double precision,dimension(nbf,nbf),intent(in)::TEMPM
     double precision,dimension(:,:,:,:),allocatable::TMPDM2
     integer::i,j,k,l,m
     write(*,*) ' '
     write(*,*) ' TRANSFORM <IJ|KL> -> <PQ|RS>'
     ! L -> S
     allocate(TMPDM2(NBF,NBF,NBF,NBF))
     do i=1,NBF
      do j=1,NBF
       do k=1,NBF
        do m=1,NBF
         TMPDM2(i,j,k,m)=0.0d0
         do l=1,NBF
          TMPDM2(i,j,k,m)=TMPDM2(i,j,k,m)+TEMPM(l,m)*ERImol(i,j,k,l)
         enddo
        enddo
       enddo
      enddo
     enddo
     write(*,*) ' L -> S DONE'
     ! K -> R
     do i=1,NBF
      do j=1,NBF
       do m=1,NBF
        do l=1,NBF
         ERImol(i,j,m,l)=0.0d0
         do k=1,NBF
          ERImol(i,j,m,l)=ERImol(i,j,m,l)+TEMPM(k,m)*TMPDM2(i,j,k,l)
         enddo
        enddo
       enddo
      enddo
     enddo
     write(*,*) ' K -> R DONE'
     ! J -> Q
     do i=1,NBF
      do m=1,NBF
       do k=1,NBF
        do l=1,NBF
         TMPDM2(i,m,k,l)=0.0d0
         do j=1,NBF
          TMPDM2(i,m,k,l)=TMPDM2(i,m,k,l)+TEMPM(j,m)*ERImol(i,j,k,l)
         enddo
        enddo
       enddo
      enddo
     enddo
     write(*,*) ' J -> Q DONE'
     ! I -> P
     do m=1,NBF
      do j=1,NBF
       do k=1,NBF
        do l=1,NBF
         ERImol(m,j,k,l)=0.0d0
         do i=1,NBF
          ERImol(m,j,k,l)=ERImol(m,j,k,l)+TEMPM(i,m)*TMPDM2(i,j,k,l)
         enddo
        enddo
       enddo
      enddo
     enddo
     write(*,*) ' I -> P DONE'
     write(*,*) ' TRANSFORMATION DONE'
     write(*,*) ' '
     deallocate(TMPDM2)
     end subroutine transformERI

! print2eint     
!     subroutine print2eint(nbf,erimol)
!     implicit none
!     integer,intent(in)::nbf
!     double precision,dimension(nbf,nbf,nbf,nbf),intent(in)::erimol
!     integer::i,j,k,l
!     do i=1,nbf
!      do j=1,nbf
!       do k=1,nbf
!        do l=1,nbf 
!          write(*,*) i,j,k,l,ERImol(i,j,k,l)
!        enddo
!       enddo
!      enddo
!     enddo
!     end subroutine print2eint

!======================================================================!
