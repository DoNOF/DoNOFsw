!======================================================================!
!                                                                      !
!           H U B B A R D   M O D E L  S U B R O U T I N E S           !
!                                                                      !
!======================================================================!
!======================================================================!

! HUBBARD
      SUBROUTINE HUBBARD(IRUNTYP,ICHARG,MULT,IDONTW)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL RESTART,HFID,CONVGDELAG,CONVG,COEF21
      COMMON/EFLDC_2/EVEC(3)
      COMMON/WRTGCF/IWRTGCF
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_RSTRT/RESTART,INPUTGAMMA,INPUTC,INPUTFMIUG,INPUTCXYZ
      COMMON/INPNOF_HFID/HFID,NTHRESHEID,THRESHEID,MAXITID,KOOPMANS
      COMMON/INPNOF_GENERALINF/ICOEF,ISOFTMAX,IORBOPT,MAXIT,MAXIT21
      COMMON/INPFILE_NIJKL/NINTMX,NIJKL,NINTCR,NSTORE
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5                                                        
      COMMON/EHFEN/EHF,EN
      COMMON/PUNTEROSUSER/N1,N2,N3,N4,N5,N6,N7,N8,N9,N10,N11,N12,N13,   &
                          N14,N15,N16,N17,N18,N19,N20,N21,N22,N23,N24,  &
                          N25,N26,N27,N28,N29,N30,N31,N32,N33,N34,N35,  &
                          N36,N37,N38,N39,N40,N41,N42,N43,N44,N45,N46,  &
                          N47,N48,N49,N50,N51,NUSER
      COMMON/ENERGIAS/EELEC,EELEC_OLD,DIF_EELEC,EELEC_MIN
#include "mpip.h"      
      INTEGER,ALLOCATABLE,DIMENSION(:)::IJKL
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::XIJKL,USER,EiHF
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::AHCORE,OVERLAP      
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::COEF,CHF
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:):: GAMMA,FMIUG0
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::   ELAGN,RON
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:):: ELAG,COEFN
!
      CHARACTER*4,ALLOCATABLE,DIMENSION(:)::CDUMMY
      INTEGER,ALLOCATABLE,DIMENSION(:)::IDUMMYA,IDUMMYP,IDUMMYS
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::RDUMMYA,RDUMMYP,DIPS
!----------------------------------------------------------------------!
!                   --- INPHUB NAMELIST VARIABLES ---                  !
!----------------------------------------------------------------------!
!
! NSITE            Number of sites in one dimension
!       = 1        (Default)
!
! NELEC            Number of electrons
!       = 1        (Default)
!
! NDIMH            Dimension considered in the Hubbard model
!       = 1        (Default)
!
! THOP             Near-neighbors hopping (t>0)
!       = 1.0d0    (Default)
!
! UONS             On-site energy = The site interaction parameter (U)
!       = 1.0d0    (Default)
!    
! IHTR             Householder transformation of the site 1RDM
!       = 0        (Default)
!-----------------------------------------------------------------------
      NAMELIST/INPHUB/NSITE,NELEC,NDIMH,THOP,UONS,IHTR
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Preset values to namelist variables
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      NSITE = 1
      NELEC = 1
      NDIMH = 1
      THOP  = 1.0d0
      UONS  = 1.0d0
      IHTR  = 0
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Read namelist variables
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#ifdef MPI
      CONTINUE
#else
      REWIND(5)
#endif
      READ(5,INPHUB,ERR=100,END=100)
!      
      IF(NDIMH==1)THEN
       WRITE(6,1)'1D'
       NBF = NSITE
      ELSE IF(NDIMH==2)THEN
       WRITE(6,1)'2D'
       NBF = NSITE*NSITE
      ELSE
       WRITE(6,2)NDIMH
       CALL ABRT                                                      
      END IF
!      
      IF(THOP<=0.0d0)THEN
       WRITE(6,3)THOP
       CALL ABRT                                                      
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Check Charge and Multiplicity
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      NA = (NELEC+MULT-1)/2            ! Number of Alpha Electrons
      NB = (NELEC-MULT+1)/2            ! Number of Beta Electrons
      IF(NA+NB/=NELEC)THEN
       WRITE(6,4)NELEC
       WRITE(6,5)ICHARG,MULT 
       CALL ABRT                                                   
      END IF
      IMULNE=MULT+NELEC                                                    
      IF(MULT>NELEC+1 .or. MULT<0 .or. 2*INT(IMULNE/2)==IMULNE)THEN  
       WRITE(6,6)    
       CALL ABRT                                                    
      END IF
      WRITE(6,7)NELEC,ICHARG,MULT,NSITE,NBF,THOP,UONS
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Print Input Run Options
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      WRITE(6,8)IRUNTYP,MULT,ICHARG
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Near-neighbors Hopping (t>0)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(AHCORE(NBF,NBF))
      AHCORE = 0.0d0
      IF(NDIMH==1)THEN
       do i=1,NSITE-1
        AHCORE(i,i+1) = -THOP
        AHCORE(i+1,i) = -THOP
       end do
!      Boundary conditions
       AHCORE(NSITE,1) = -THOP
       AHCORE(1,NSITE) = -THOP
      ELSE IF(NDIMH==2)THEN
       do j=1,NSITE-1
        do i=1,NSITE-1
         ij  = j + (i-1)*nsite           ! ij = (i,j)
         i1j = j + i*nsite               ! i1j = (i+1,j)
         AHCORE(ij,i1j) = -THOP          ! h[(i,j),(i+1,j)]
         AHCORE(i1j,ij) = -THOP          ! h[(i+1,j),(i,j)]
         ij1 = j+1 + (i-1)*nsite         ! ij1 = (i,j+1)
         AHCORE(ij,ij1) = -THOP          ! h[(i,j),(i,j+1)]
         AHCORE(ij1,ij) = -THOP          ! h[(i,j+1),(i,j)]
        end do
       end do               
!      Boundary conditions        
       do i=1,NSITE-1                    ! last column (j=nsite)
        in = i*nsite                     ! in = (i,nsite)       
        i1 = 1 + (i-1)*nsite             ! i1 = (i,1)       
        AHCORE(in,i1) = -THOP            ! h[(i,nsite),(i,1)]   
        AHCORE(i1,in) = -THOP            ! h[(i,1),(i,nsite)]
        i1n= (i+1)*nsite                 ! i1n = (i+1,nsite)        
        AHCORE(in,i1n) = -THOP           ! h[(i,nsite),(i+1,nsite)]
        AHCORE(i1n,in) = -THOP           ! h[(i+1,nsite),(i,nsite)]
       end do
       do j=1,NSITE-1                    ! last row (i=nsite)
        nj = j + (nsite-1)*nsite         ! nj = (nsite,j)
        AHCORE(nj,j) = -THOP             ! h[(nsite,j),(1,j)]
        AHCORE(j,nj) = -THOP             ! h[(1,j),(nsite,j)]
        nj1 = j+1 + (nsite-1)*nsite      ! nj1 = (nsite,j+1)
        AHCORE(nj,nj1) = -THOP           ! h[(nsite,j),(nsite,j+1)]
        AHCORE(nj1,nj) = -THOP           ! h[(j+1,nsite),(nsite,j)] 
       end do       
       nn = nsite*nsite                  ! last element (nsite,nsite)
       n1 = 1 + (nsite-1)*nsite          ! n1 = (nsite,1)
       AHCORE(nn,n1)    = -THOP          ! h[(nsite,nsite),(nsite,1)]
       AHCORE(n1,nn)    = -THOP          ! h[(nsite,1),(nsite,nsite)] 
       AHCORE(nn,nsite) = -THOP          ! h[(nsite,nsite),(1,nsite)]
       AHCORE(nsite,nn) = -THOP          ! h[(1,nsite),(nsite,nsite)]              
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     OVERLAP
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(OVERLAP(NBF,NBF))
      OVERLAP = 0.0d0
      DO I=1,NBF
       OVERLAP(I,I) = 1.0d0
      ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     On-site interaction parameter (U)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      NINTMX = NBF  ! Number of 2e-integrals (ERIs) per record
      NIJKL  = NBF  ! Total Number of 2e-integrals
      NINTCR = NBF  ! Space needed to allocate 2e- integrals in Slaves      
      NSTORE = NBF  ! Space needed to allocate 2e- integrals in Master
!
      ALLOCATE(IJKL(NIJKL),XIJKL(NIJKL))          
      DO I=1,NIJKL
       IJKL(I) = ISHFT(I,48) + ISHFT(I,32) + ISHFT(I,16) + I
       XIJKL(I) = UONS  
      END DO
      IF(IDONTW==0)WRITE(1)-NIJKL,IJKL,XIJKL
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Header on the output file
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(IDUMMYA(NSITE),RDUMMYA(3*NSITE))
      IDUMMYA = 0
      RDUMMYA = 0.0d0
      CALL RUNNOFHEADER(NSITE,ICHARG,MULT,NBF,0,NBF,NELEC,NA,NB,        &
                        EVEC(1),EVEC(2),EVEC(3),0,0,0,IDUMMYA,0,0,      &
                        IRUNTYP,RDUMMYA,RDUMMYA)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Write GCF file for Restart calculations
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IWRTGCF = 1
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Initial Orbitals
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(COEF(NBF,NBF))      
      COEF = 0.0d0      
      DO I=1,NBF
       COEF(I,I) = 1.0d0
      ENDDO
      EN = 0.0d0                               ! Nuclear Energy
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Create the basis function symbol table
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      WRITE(4,'(I5)')NBF
      DO I = 1,NBF
       WRITE(4,'(A2,I2,A4)')'    ',I,'  S '
      END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Read two-electron Repulsion Integrals in AO basis (ERI)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL READERIs(IJKL,XIJKL,IJKL,XIJKL,NBF,IDONTW,1)      
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Allocate User array
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(USER(NUSER))
      USER(N11:NUSER) = 0.0d0                  ! IEMOM = 0 
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     FIRSTCALL: Initialize variables    COEF,      GAMMA,     FMIUG0
!                         according to INPUTC, INPUTGAMMA, INPUTFMUIG
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(GAMMA(NBF5),FMIUG0(NBF))
      CALL INITr(COEF,OVERLAP,GAMMA,FMIUG0,1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Restricted Hartree-Fock (RHF)                                    
!     Use the Iterative Diagonalization Method to generate the HF MOs  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(HFID)THEN
       ALLOCATE(CHF(NBF,NBF),EiHF(NBF))
       CHF = COEF  
       CALL HFIDr(AHCORE,IJKL,XIJKL,XIJKL,CHF,EiHF,USER,1)
       if(INPUTC==0)COEF = CHF   ! Input for Coefficients is CHF
       if(INPUTFMIUG==0)FMIUG0 = EiHF
      ELSE
       if(INPUTC==0)WRITE(6,9)
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
      EELEC_MIN = 1.0d20        ! GLOBAL FIRST CALL
      IF(ICOEF==21)THEN
       ICOEF=2
       COEF21=.TRUE.
      ENDIF 
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Dummy Arrays
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(IDUMMYP(1),IDUMMYS(1))
      ALLOCATE(CDUMMY(NSITE),RDUMMYP(1),DIPS(3))
      IDUMMYA = 0
      IDUMMYP = 0
      IDUMMYS = 0      
      CDUMMY = '    '      
      RDUMMYP = 0.0d0      
      DIPS = 0.0d0
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -      
      CALL OccOpt(IFIRSTCALL,CONVG,CDUMMY,RDUMMYA,OVERLAP,IDUMMYA,      &
                  IDUMMYA,COEF,GAMMA,FMIUG0,AHCORE,IJKL,XIJKL,RDUMMYA,  &
                  ELAG,USER,IDUMMYA,RDUMMYA,RDUMMYA,RDUMMYA,IDUMMYS,    &
                  IDUMMYS,IDUMMYS,IDUMMYS,IDUMMYS,IDUMMYS,IDUMMYS,      &
                  IDUMMYA,IDUMMYA,IDUMMYP,IDUMMYP,RDUMMYP,RDUMMYP,      &
                  RDUMMYP,RDUMMYP,RDUMMYP,RDUMMYP,RDUMMYP,RDUMMYP,      &
                  RDUMMYP,RDUMMYP,ELAGN,COEFN,RON,IT,ITTOTAL,DIPS,1,1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     END SINGLE-POINT CALCULATION (ICOEF=0)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(ICOEF==0)THEN
       WRITE(6,10)IT
       GOTO 1000
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
        IF(IORBOPT==1)THEN       
         CALL OrbOptFMIUGr(IT,ITLIM,OVERLAP,AHCORE,IJKL,XIJKL,RDUMMYA,  &
                           USER(N7),COEF,USER(N1),USER(N2),USER(N3),    &
                           ELAG,FMIUG0,USER(N11),USER(N12),USER(N13),   &
                           USER(N14),ILOOP,1)
        ELSE IF(IORBOPT==2)THEN                           
         CALL OrbOptRot(IT,OVERLAP,AHCORE,IJKL,XIJKL,RDUMMYA,USER(N7),  &
                        COEF,USER(N1),USER(N2),USER(N3),ELAG,USER(N11), &
                        USER(N12),USER(N13),USER(N14),ILOOP,1)
        END IF                       
!      Core-Fragment Orbital Optimization
       ELSEIF(ICOEF==3)THEN
        WRITE(6,14)
        STOP       
       ENDIF 
       
!      Occupation Optimization
       ITTOTAL=ITTOTAL+ILOOP
       CALL OccOpt(IFIRSTCALL,CONVG,CDUMMY,RDUMMYA,OVERLAP,IDUMMYA,     &
                   IDUMMYA,COEF,GAMMA,FMIUG0,AHCORE,IJKL,XIJKL,RDUMMYA, &
                   ELAG,USER,IDUMMYA,RDUMMYA,RDUMMYA,RDUMMYA,IDUMMYS,   &
                   IDUMMYS,IDUMMYS,IDUMMYS,IDUMMYS,IDUMMYS,IDUMMYS,     &
                   IDUMMYA,IDUMMYA,IDUMMYP,IDUMMYP,RDUMMYP,RDUMMYP,     &
                   RDUMMYP,RDUMMYP,RDUMMYP,RDUMMYP,RDUMMYP,RDUMMYP,     &
                   RDUMMYP,RDUMMYP,ELAGN,COEFN,RON,IT,ITTOTAL,DIPS,1,1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      END SINGLE-POINT CALCULATION (CONVG=TRUE)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       IF(CONVG)THEN
!       Householder Transformation of 1RDM in the site basis
        IF(IHTR==1)THEN
         CALL HTR1RDM(OVERLAP,USER(N1),USER(N7),NBF5,NBF)
        ENDIF
!       Final Output
        WRITE(6,12)
        WRITE(6,13)IT,ITTOTAL
        GOTO 1000
       ENDIF
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
      ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     EXCESSIVE NUMBER OF ITERATIONS
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(IT>MAXIT)THEN
       WRITE(6,11)
       WRITE(6,13)IT,ITTOTAL
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     STOP PROGRAM, DEALLOCATE MEMORY, GIVES ELAPSED TIME if IT > MAXIT
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 1000 CONTINUE 
      DEALLOCATE(AHCORE,OVERLAP,IJKL,XIJKL,COEF,USER)
      DEALLOCATE(GAMMA,FMIUG0,ELAG,ELAGN,COEFN,RON)
      DEALLOCATE(IDUMMYA,IDUMMYP,IDUMMYS,CDUMMY,RDUMMYA,RDUMMYP,DIPS)
      if(HFID)DEALLOCATE(CHF,EiHF)
      NOPTCGMPI = 1      
#ifdef MPI
      DO I=1,NPROCS-1
       K=0
       CALL MPI_SEND(K,1,MPI_INTEGER8,I,I,MPI_COMM_WORLD,IERR)
       CALL MPI_SEND(NOPTCGMPI,1,MPI_INTEGER8,I,I,MPI_COMM_WORLD,IERR)
      ENDDO
#endif 
      RETURN      
!-----------------------------------------------------------------------
!     Format definitions
!-----------------------------------------------------------------------
    1 FORMAT(/,19X,'HUBBARD MODEL CALCULATION (',A2,')')
    2 FORMAT(/1X,'Stop: NDIMH =',I2,' and must be 1 or 2',/)    
    3 FORMAT(/1X,'Stop: THOP =',F6.2,' and must be greater than 0.0',/)
    4 FORMAT(/1X,'Number of Electrons =',I6)
    5 FORMAT(/1X,'Check your Charge and Multiplicity:',2I6)
    6 FORMAT(/1X,'Impossible Spin Multiplicity')
    7 FORMAT(/1X,'Number of Electrons   (NELEC)  =',I5/                 &
              1X,'Charge of the System  (ICHARG) =',I5/                 &
              1X,'Spin Multiplicity     (MULT)   =',I5/                 &
              1X,'Number of Sites in 1D (NSITE)  =',I5/                 &
             /1X,'Total number of sites             =',I5/              & 
              1X,'Near-neighbors hopping (t>0)      =',F7.2/            &
              1X,'On-site interaction parameter (U) =',F7.2)
    8 FORMAT(/1X,'Input Run Options',/,                                 &
              1X,'-----------------',/,                                 &
              1X,'IRUNTYP =',I2,2X,'MULT =',I2,2X,'ICHARG =',I3)
    9 FORMAT(/1X,'Input for Coefficients is Site Basis')
   10 FORMAT(//2X,'**************************************************', &
              /2X,'*                                                *', &
              /2x,'*       SINGLE-POINT DoNOF CALCULATION           *', &
              /2X,'*                                                *', &
              /2X,'*            No.ITER =',I6,'                     *', &
              /2X,'*         (Occupation Optimization)              *', &
              /2X,'*                                                *', &
              /2x,'*  FINAL RESULTS   FINAL RESULTS  FINAL RESULTS  *', &
              /2X,'*                                                *', &
              /2X,'**************************************************')
   11 FORMAT(/,10X,30(1H-),/,10X,'EXCESSIVE NUMBER OF ITERATIONS',      &
             /,10X,30(1H-))
   12 FORMAT(//2X,'**************************************************', &
              /2X,'*                                                *', &
              /2x,'*       SINGLE-POINT DoNOF CALCULATION           *', &
              /2X,'*                                                *', &
              /2X,'*             FULL OPTIMIZATION                  *', &
              /2X,'*                                                *', &
              /2x,'*  FINAL RESULTS   FINAL RESULTS  FINAL RESULTS  *', &
              /2X,'*                                                *', &
              /2X,'**************************************************')
   13 FORMAT(/2X,'**************************************************',  &
             /2X,'*         No. EXTERNAL ITER =',I6,'              *',  &
             /2X,'*         No. of TOTAL ITER =',I6,'              *',  &
             /2X,'**************************************************')
   14 FORMAT(/1X,'ICOEF cannot be equal 3')
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  100 WRITE(6,'(/2X,36A)')'Stop: Wrong INPHUB Namelist Variable'
      STOP
!-----------------------------------------------------------------------
      END

!======================================================================!
!                                                                      !
!         HUBBARD: HOUSEHOLDER TRANSFORMATION OF THE SITE 1RDM         ! 
!                                                                      !
!======================================================================!

! HTR1RDM
      SUBROUTINE HTR1RDM(OVERLAP,RO,QD,NBF5,NBF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION,DIMENSION(NBF5)::RO
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::OVERLAP
      DOUBLE PRECISION,DIMENSION(NBF,NBF,NBF)::QD
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::HLD
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::RDM1,RDM1BL,AUX
!-----------------------------------------------------------------------
      WRITE(6,1)
      
      ALLOCATE (RDM1(NBF,NBF),RDM1BL(NBF,NBF),HLD(NBF),AUX(NBF,NBF))
      DO ieta=1,NBF
       DO imiu=1,NBF
        RDM1(ieta,imiu) = SUMDL(ieta,imiu,RO,QD)
       ENDDO
      ENDDO

      write(6,2)'RDM1:  '
      do i = 1,nbf
       write(6,3)(RDM1(j,i),j=1,nbf)
      end do
      
!     Check the normalization of the 1RDM
      CALL RDM1NORM(OVERLAP,RO,QD)
!      
      CALL HOUSEHOLDER(NBF,RDM1,AUX,RDM1BL,HLD)
      
      write(6,2)'RDM1BL:'
      do i = 1,nbf
       write(6,3)(RDM1BL(j,i),j=1,nbf)
      end do
!-----------------------------------------------------------------------
    1 FORMAT(/' Householder Transformation of the 1RDM ',               &
             /' -------------------------------------- ')
    2 FORMAT(/,1X,A7,/)
    3 FORMAT(20F8.3)
!-----------------------------------------------------------------------    
      RETURN
      END

! HOUSEHOLDER
      SUBROUTINE HOUSEHOLDER(nsize,Gamma,P,Gamma_bl,V)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
!      Compute the Full Householder transformation of a given matrix   !
!      Proceed the Householder tridiagonal transformation of matrix    !
!      Gamma of nsize x nsize. Implemented by Sekaran Sajanthan (2022) !
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
      IMPLICIT NONE
      INTEGER :: nsize,i,j,k,l
      DOUBLE PRECISION :: sum_M, alpha_ns, alpha, r, x_sign      
      DOUBLE PRECISION,DIMENSION(nsize) :: V
      DOUBLE PRECISION,DIMENSION(nsize,nsize) :: Gamma,P,Gamma_bl,Work
!-----------------------------------------------------------------------
      DO k = 1, 1                                 ! nsize - 1

       sum_M = 0.d0
       do j = k + 1, nsize
        sum_M = sum_M + Gamma(j,k)*Gamma(j,k)
       end do

       alpha_ns = -SQRT(sum_M)
       x_sign = SIGN(1.d0, Gamma(k+1,k))
       alpha = x_sign*alpha_ns

       r = SQRT(5.d-1)*SQRT(alpha*alpha - alpha*Gamma(k+1,k))

!      Construct V_k
       do j = 1, k
        V(j) = 0.d0
       end do
       V(k+1) = (Gamma(k+1,k) - alpha) / (2.d0*r)
       do j = k+2,nsize
        V(j) = Gamma(j,k)/(2.d0*r)
       end do

!      Construct P_k
       P = 0.0d0
       do i = 1,nsize
        do j = 1,nsize
         if(i==j) P(i,j) = 1.0d0
         P(i,j) = P(i,j) - 2.0d0*V(i)*V(j)
        end do
       end do
      
!      Matrix product: Get Gamma_bl
       Work = 0.0d0
       do i = 1, nsize
        do j = 1, nsize
         do l = 1, nsize
          Work(i, j) = Work(i, j) + Gamma(i, l)*P(l, j)
         end do
        end do
       end do
       
       Gamma_bl = 0.0d0
       do i = 1, nsize
        do j = 1, nsize
         do l = 1, nsize
          Gamma_bl(i, j) = Gamma_bl(i, j) + P(i, l)*Work(l, j)
         end do
        end do
       end do
       
      END DO
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE HOUSEHOLDER

!----------------------------------------------------------------------!
