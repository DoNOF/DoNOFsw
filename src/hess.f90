!======================================================================!
!                                                                      !
!             N O F   H E S S I A N   S U B R O U T I N E S            !
!                                                                      !
!      Computation of Hessian numerically from analytic gradients      !
!                                                                      !
!              2018  Module implemented by Ion Mitxelena               !
!                                                                      !
!======================================================================!
!                                                                      !
!   HESSCAL: Calculate Hessian from analytic gradients (IRUNTYP=4)     !
!   HSSNUMd: Main routine, use the 6-grid numerical formula            !
!   GRADSPURIFY: Consider gradients only above the threshold = 10-6    !
!   SETFCMd: Form Hessian from the gradients                           !
!   SYMFCMd: Symmetrize the Hessian                                    !
!                                                                      !
!   Dipole derivative matrix calculation for obtaining IR intensities  !
!   SETDDMd: Numerical Dipole derivative matrix computation            !
!   SYMDDMd: Purify the dipole derivative tensor below a threshold     !
!                                                                      !
!   Harmonic Vibrational analysis                                      !
!   FGMTRXd: Main routine to do Harm Vib analysis                      !
!   HESMASd: Computes the mass-weighted Hessian                        !
!   CENMASd: Computes the center of mass with mass-weighting           !
!                                                                      !
!   PRJFCMd: Projection of the Force Constant Matrix (FCM)             !
!                                                                      !
!======================================================================!

! HESSCAL
      SUBROUTINE HESSCAL(NINTEG,IDONTW,NAT,ZAN,Cxyz,IAN,IMIN,IMAX,      &
                 ZMASS,KSTART,KATOM,KTYPE,KLOC,INTYP,KNG,KMIN,KMAX,     &
                 ISH,ITYP,C1,C2,EX1,CS,CP,CD,CF,CG,CH,CI,DIPS,          &
                 GRADS,IRUNTYP,IPROJECT,ISIGMA,SIZE_ENV,ENV,ATM,        &
                 NBAS,BAS,IGTYP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPFILE_Naux/NBFaux,NSHELLaux
      COMMON/WRTGCF/IWRTGCF
      COMMON/ELPROP/IEMOM
      INTEGER,DIMENSION(NAT) :: IAN,IMIN,IMAX
      INTEGER,DIMENSION(NSHELL) :: KSTART,KATOM,KTYPE,KLOC
      INTEGER,DIMENSION(NSHELL) :: INTYP,KNG,KMIN,KMAX
      INTEGER,DIMENSION(NPRIMI) :: ISH,ITYP
      DOUBLE PRECISION,DIMENSION(NAT) :: ZAN,ZMASS
      DOUBLE PRECISION,DIMENSION(NPRIMI)::C1,C2,EX1,CS,CP,CD,CF,CG,CH,CI
      DOUBLE PRECISION,DIMENSION(3,NAT):: Cxyz
      DOUBLE PRECISION,DIMENSION(3*NAT):: GRADS
      DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE:: HESSIANO,DDM
      DOUBLE PRECISION,DIMENSION(3):: DIPS

      INTEGER :: SIZE_ENV,NBAS,IGTYP
      DOUBLE PRECISION :: ENV(SIZE_ENV)
      INTEGER :: ATM(6,NAT), BAS(8,NBAS)
!-----------------------------------------------------------------------
      IF (IRUNTYP==4) THEN
!      Update coordinates of shells if use libint library for ERIs
       CALL ENERGRAD(NINTEG,IDONTW,IEMOM,NAT,NBF,NBFaux,NSHELL,NPRIMI,  &
                     ZAN,Cxyz,IAN,IMIN,IMAX,KSTART,KATOM,KTYPE,KLOC,    &
                     INTYP,KNG,KMIN,KMAX,ISH,ITYP,C1,C2,EX1,CS,CP,CD,   &
                     CF,CG,CH,CI,GRADS,IRUNTYP,DIPS,SIZE_ENV,ENV,ATM,   &
                     NBAS,BAS,IGTYP,0,1)
!      Write Coordinates on File CGGRAD (Unit=11)
       WRITE(11,1)
       DO I=1,NAT
        WRITE(11,2)I,Cxyz(1,I),Cxyz(2,I),Cxyz(3,I)
        ENV(20+3*(I-1)+1) = Cxyz(1,I)
        ENV(20+3*(I-1)+2) = Cxyz(2,I)
        ENV(20+3*(I-1)+3) = Cxyz(3,I)
       ENDDO
!      Internuclear distances       
       CALL NUCDIST(3*NAT,NAT,Cxyz)
      ENDIF
!     Compute Hessian from analytic gradients at stationary point
      IWRTGCF = 0
      ALLOCATE (HESSIANO(3*NAT,3*NAT),DDM(9,NAT))
      CALL HSSNUMd(HESSIANO,3*NAT,Cxyz,GRADS,DIPS,DDM,NINTEG,IDONTW,    &
                   ZAN,IAN,IMIN,IMAX,KSTART,KATOM,KTYPE,KLOC,INTYP,     &
                   KNG,KMIN,KMAX,ISH,ITYP,C1,C2,EX1,CS,CP,CD,CF,CG,     &
                   CH,CI,IRUNTYP,SIZE_ENV,ENV,ATM,NBAS,BAS,IGTYP)
!     Carry out normal mode vibrational analysis
      CALL FGMTRXd(Cxyz,HESSIANO,GRADS,3*NAT,ZAN,ZMASS,DDM,IPROJECT,    &
                   ISIGMA)
      DEALLOCATE (HESSIANO,DDM)
      RETURN
!-----------------------------------------------------------------------
1     FORMAT(/1X,'Atom',16X,'Coordinates (Bohr)',/16X,                  &
       'x',13X,'y',13X,'z')
2     FORMAT(1X,I4,3F14.4)
!-----------------------------------------------------------------------    
      END

! HSSNUMd      
      SUBROUTINE HSSNUMd(FCM,NC1,Cxyz,EG,DIP,DDM,NINTEG,IDONTW,         &
                         ZAN,IAN,IMIN,IMAX,KSTART,KATOM,KTYPE,          &
                         KLOC,INTYP,KNG,KMIN,KMAX,ISH,ITYP,C1,C2,       &
                         EX1,CS,CP,CD,CF,CG,CH,CI,IRUNTYP,              &
                         SIZE_ENV,ENV,ATM,NBAS,BAS,IGTYP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      LOGICAL CONVGDELAG,FROZEN
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/CONVERGENCE/DUMEL,PCONV,CONVGDELAG      
      COMMON/INPNOF_FROZEN/FROZEN,IFROZEN(200)
      COMMON/INPNOF_MOLDEN/MOLDEN
      COMMON/INPFILE_Naux/NBFaux,NSHELLaux
      COMMON/EHFEN/EHF,EN
      COMMON/ENERGIAS/EELEC,EELEC_OLD,DIF_EELEC,EELEC_MIN
      COMMON/ELPROP/IEMOM      
      COMMON/ECP2/CLP(4004),ZLP(4004),NLP(4004),KFRST(1001,6),          &
                  KLAST(1001,6),LMAX(1001),LPSKIP(1001),IZCORE(1001)
!     ARGUMENTS
      INTEGER,INTENT(IN) :: NC1,NINTEG,IDONTW
      DOUBLE PRECISION,DIMENSION(NC1,NC1),INTENT(OUT) :: FCM
      DOUBLE PRECISION,DIMENSION(NC1),INTENT(IN) :: Cxyz
      DOUBLE PRECISION,DIMENSION(NC1),INTENT(INOUT) :: EG
      DOUBLE PRECISION,DIMENSION(3),INTENT(INOUT) :: DIP
      DOUBLE PRECISION,DIMENSION(9,NC1/3),INTENT(OUT)::DDM
      DOUBLE PRECISION,DIMENSION(NC1/3),INTENT(IN) :: ZAN
      INTEGER,DIMENSION(NC1/3),INTENT(IN):: IAN,IMIN,IMAX
      INTEGER,DIMENSION(NSHELL),INTENT(IN) :: KSTART,KATOM,KTYPE,KLOC
      INTEGER,DIMENSION(NSHELL),INTENT(IN) :: INTYP,KNG,KMIN,KMAX
      INTEGER,DIMENSION(NPRIMI),INTENT(IN) :: ISH,ITYP
      DOUBLE PRECISION,DIMENSION(NPRIMI),INTENT(IN)::C1,C2,EX1,CS,CP
      DOUBLE PRECISION,DIMENSION(NPRIMI),INTENT(IN)::CD,CF,CG,CH,CI

      INTEGER :: SIZE_ENV,NBAS,IGTYP
      DOUBLE PRECISION :: ENV(SIZE_ENV)
      INTEGER :: ATM(6,NATOMS), BAS(8,NBAS)

!     VARIABLES      
      DOUBLE PRECISION,DIMENSION(NC1) :: CDISP,EGDISP
      DOUBLE PRECISION,DIMENSION(2) :: D
      DOUBLE PRECISION,DIMENSION(3) :: DEQ
      LOGICAL,DIMENSION(NATOMS)::SKIP
      DOUBLE PRECISION,PARAMETER :: UNIT=0.52917724924D+00
      DOUBLE PRECISION,PARAMETER :: ZERO=0.0D+00
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
!     NVIB   = The number of displacements in each Cartesian
!              direction for force field computation.C
!     VIBSIZ = Displacement size in Bohrs. Default=0.01     
!-----------------------------------------------------------------------
      WRITE(11,1)
      WRITE(6,'(/,72(1H-),/)')
      WRITE(6,*)'Note: Numerical Hessian and Frequencies in CGO file !!'
      WRITE(6,'(/,72(1H-),/)')
      VIBSIZ = 0.01D+00
      NVIB = 2
      D(1) =  VIBSIZ
      D(2) = -VIBSIZ
      DEL = VIBSIZ*NVIB*UNIT
      SKIP = .FALSE.
      CDISP = Cxyz
!     Purify Gradients by Householder method
      CALL GRADSPURIFY(EG,NC1)
      WRITE(11,3)
      EGDISP = EG
!     Print Energy & Convergence achieved by ELag
      WRITE(11,4)EELEC+EN,DUMEL
!     Print Gradients
      WRITE(11,5)
      DO I = 1,NATOMS
      WRITE(11,'(I4,3F15.4)')I,EG(1+(I-1)*3),EG(2+(I-1)*3),EG(3+(I-1)*3)
      ENDDO
      WRITE(11,*)      
      CALL SETDDMd(DDM,DIP,DEL,DEQ,0,NVIB,NC1,NVIB)
      CALL SETFCMd(FCM,NC1,NC1,EG,0)
!     IDENTIFY FROZEN ATOMS
      IF(FROZEN) THEN
       WRITE(11,40)
       DO I = 1,200,2
        IF(IFROZEN(I).EQ.0) EXIT
        SKIP(IFROZEN(I+1)) = .TRUE.
       ENDDO
      END IF
!
      NMAXSKIP = NATOMS
      DO I = NATOMS,1,-1
       IF(SKIP(I) .eqv. .FALSE.)THEN
        NMAXSKIP = I
        EXIT
       ENDIF
      ENDDO
!
      NDISPL=0
      NOPTCG=0
      III=0
      WRITE(11,6)
      IF(MOLDEN==1)THEN
       WRITE(18,7)
      ENDIF 
      DO IVIB = 1,NVIB
       DO IAT = 1,NATOMS
        IF(SKIP(IAT)) CYCLE
        NVA = 3*(IAT-1)
        DO ICOORD = 1,3
         III = III + 1        
         NV = NVA+ICOORD
         CDISP(NV) = CDISP(NV)+D(IVIB)
         DO I=1,NV
           ENV(20+I) = CDISP(I)
         END DO
!        ENERGY AND GRADIENT AT DISPLACED GEOM 
         IF(IVIB==NVIB.AND.IAT==NMAXSKIP.AND.ICOORD==3)NOPTCG=1
!        Update coordinates of shells if use libint library for ERIs
         CALL ENERGRAD(NINTEG,IDONTW,IEMOM,NATOMS,NBF,NBFaux,NSHELL,    &
                       NPRIMI,ZAN,CDISP,IAN,IMIN,IMAX,KSTART,KATOM,     &
                       KTYPE,KLOC,INTYP,KNG,KMIN,KMAX,ISH,ITYP,         &
                       C1,C2,EX1,CS,CP,CD,CF,CG,CH,CI,EGDISP,           &
                       IRUNTYP,DIP,SIZE_ENV,ENV,ATM,                    &
                       NBAS,BAS,IGTYP,NOPTCG,0)
!        PURIFY THE GRADIENT BY HOUSEHOLDER METHOD
         CALL GRADSPURIFY(EGDISP,NC1)
         ENERGY = EELEC + EN
!        Displaced Geometries for MOLDEN
         if(MOLDEN==1)then   
          write(18,8)NATOMS
          write(18,9)III,ENERGY
          do jat=1,natoms
           jo = (jat-1)*3
           IZNUC = INT(ZAN(jat))+IZCORE(jat)                 
           write(18,10)ATMLAB(IZNUC),                                   &
           CDISP(1+jo)*UNIT,CDISP(2+jo)*UNIT,CDISP(3+jo)*UNIT
          enddo
         endif
!        MOVE THIS DISPLACED ATOM BACK TO WHERE IT CAME FROM
         CDISP(NV) = CDISP(NV)-D(IVIB)
!        CHECK ENERGY
         WRITE(11,15)III,ENERGY
!        COMPUTE CORRESPONDING HESSIAN AND DIPOLE DERIVATIVES
         CALL SETFCMd(FCM,NC1,NV,EGDISP,IVIB)
         CALL SETDDMd(DDM,DIP,DEL,DEQ,NV,IVIB,NC1,NVIB)
         EGDISP = ZERO
         NDISPL=NDISPL+1
        ENDDO
       ENDDO
      ENDDO      
!     COMPLETE NUMERICAL DIFFERENTIATIONS
!     SYMMETRIZE THE FORCE CONSTANT AND DIPOLE DERIVATIVE MATRIX
      CALL SYMFCMd(FCM,NC1,VIBSIZ)
      CALL SYMDDMd(DDM,NATOMS)
!     ZERO OFF FCM ELEMENTS FOR PARTIAL HESSIAN ANALYSIS      
      IF (FROZEN) THEN
      DO IFF=1,200,2
       IF(IFROZEN(IFF).EQ.0) EXIT
       IFREEZZ = 3*(IFROZEN(IFF+1)-1)+IFROZEN(IFF)
!      SET DIAGONAL FCM ELEMENTS BE 1.0D-08 FOR PARTIAL HESSIAN       
       FCM(IFREEZZ,IFREEZZ)=1.0D-08
       DO I = 1, NC1
        FCM(I,IFREEZZ)=ZERO
        FCM(IFREEZZ,I)=ZERO
       ENDDO
      ENDDO
      ENDIF
!     PRINT OUT TOTAL PNOF HESSIAN
      WRITE(11,20)
      DO I=1,NATOMS
       DO J=1,I
        WRITE(11,30)'X',I,J,FCM(1+(I-1)*3,1+(J-1)*3),                   &
                    FCM(2+(I-1)*3,1+(J-1)*3),FCM(3+(I-1)*3,1+(J-1)*3)    
        WRITE(11,30)'Y',I,J,FCM(1+(I-1)*3,2+(J-1)*3),                   &
                    FCM(2+(I-1)*3,2+(J-1)*3),FCM(3+(I-1)*3,2+(J-1)*3)    
        WRITE(11,30)'Z',I,J,FCM(1+(I-1)*3,3+(J-1)*3),                   &
                    FCM(2+(I-1)*3,3+(J-1)*3),FCM(3+(I-1)*3,3+(J-1)*3)
       ENDDO
      ENDDO
!     END HESSIAN CALCULATION
      WRITE(11,2)
      RETURN
!-----------------------------------------------------------------------
    1 FORMAT(/4X,'- START NUMERICAL HESSIAN CALCULATION -')
    2 FORMAT(/4X,'- END OF NUMERICAL HESSIAN CALCULATION -')
    3 FORMAT(/1X,'1E-06 is used for Gradient Threshold')
    4 FORMAT(/1X,'Energy:',F16.6,' (',ES7.0,' )')
    5 FORMAT(/1X,'Energy Gradient',/)
    6 FORMAT(1X,'Energy at displaced Geometry',/)
    7 FORMAT('[GEOMETRIES] (XYZ)')
    8 FORMAT(I6)
    9 FORMAT('Displaced Geometry',I6,F20.10)
   10 FORMAT(1X,A4,3F15.4)
   15 FORMAT(I4,F20.10)
   20 FORMAT(//1X,'Hessian computed from analytic Gradients',           &         
              /1X,'----------------------------------------',           &
            //25X,'X',13X,'Y',13X,'Z',/)
   30 FORMAT(A4,2X,2I4,3F14.4)
   40 FORMAT(/1X,'Warning: Gradients related to a frozen Atom are 0')   
      END

! GRADSPURIFY
      SUBROUTINE GRADSPURIFY(GRADS,NV)
      IMPLICIT NONE
      INTEGER,INTENT(IN)::NV
      DOUBLE PRECISION,DIMENSION(NV),INTENT(INOUT)::GRADS
      DOUBLE PRECISION::GOLD
      INTEGER::I,J
      DOUBLE PRECISION,PARAMETER::THR=0.000001D+00,THR2=0.0000001D+00
      DOUBLE PRECISION,PARAMETER::ZERO=0.0D+00
!-----------------------------------------------------------------------
      DO I=1,NV
        IF(DABS(GRADS(I))<=THR) THEN
          GRADS(I) = ZERO
          CYCLE
        ENDIF
        DO J=1,I-1
          IF(DABS(DABS(GRADS(I))-DABS(GRADS(J)))<THR2) THEN
            GOLD = GRADS(J)
            GRADS(J) = DABS(GRADS(I))
            IF(GOLD<ZERO) GRADS(J) = - DABS(GRADS(I))
          ENDIF
        ENDDO
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! SETFCMd      
      SUBROUTINE SETFCMd(FCM,M,NV,EG,IVIB)
      IMPLICIT NONE
      INTEGER,INTENT(IN)::IVIB,M,NV
      DOUBLE PRECISION,DIMENSION(M,M),INTENT(OUT)::FCM
      DOUBLE PRECISION,DIMENSION(M),INTENT(IN)::EG
      INTEGER:: I
!-----------------------------------------------------------------------      
!     SET ELEMENTS OF THE FORCE CONSTANT MATRIX
!-----------------------------------------------------------------------
      IF(IVIB .EQ. 0) THEN
!       INITIALIZE FCM TO ZERO
        FCM = 0.0D+00
        RETURN
      ENDIF  
!     STORE COLUMN 'NV' IN THE FCM
      IF(IVIB.EQ.1) THEN
!     FIRST DIFFERENCING
       DO I = 1,M
         FCM(I,NV) = EG(I)
       ENDDO
       RETURN
      ENDIF
!     SECOND DIFFERENCING
      DO I = 1,M
         FCM(I,NV) = FCM(I,NV)-EG(I)
      ENDDO
!-----------------------------------------------------------------------      
      RETURN
      END

! SYMFCMd      
      SUBROUTINE SYMFCMd(FCM,NCOORD,VIBSIZ)
      IMPLICIT NONE
      INTEGER,INTENT(IN)::NCOORD
      DOUBLE PRECISION,INTENT(IN)::VIBSIZ
      DOUBLE PRECISION,DIMENSION(NCOORD,NCOORD),INTENT(INOUT)::FCM
      INTEGER:: I,J
      DOUBLE PRECISION:: DUM,AVE
      DOUBLE PRECISION,PARAMETER::HALF=0.5D+00
      DOUBLE PRECISION,PARAMETER::ONE=1.0D+00,TWO=2.0D+00
!-----------------------------------------------------------------------      
!      COMPLETE COMPUTATION OF THE FORCE CONSTANT MATRIX
!      MATRIX IS SYMMETRIZED.COMPLETE THE FINITE DIFFERENCING
!-----------------------------------------------------------------------
      DUM = ONE/(VIBSIZ*TWO)
      DO I = 1,NCOORD
         DO J = 1,NCOORD
            FCM(I,J) = DUM*FCM(I,J)
         ENDDO
      ENDDO
!     SYMMETRIZE THE FORCE CONSTANT MATRIX
      DO I = 2,NCOORD
         DO J = 1,I-1
            AVE = (FCM(I,J)+FCM(J,I))*HALF
            FCM(J,I) = AVE
            FCM(I,J) = AVE
         ENDDO
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! SETDDMd      
      SUBROUTINE SETDDMd(DDM,DIP,DEL,DEQ,NP,IVIB,NCOORD,NVIB)
      IMPLICIT NONE
!     ARGUMENTS      
      INTEGER,INTENT(IN):: NP,IVIB,NCOORD,NVIB
      DOUBLE PRECISION,DIMENSION(3,NCOORD),INTENT(INOUT):: DDM
      DOUBLE PRECISION,DIMENSION(3),INTENT(INOUT):: DEQ
      DOUBLE PRECISION,DIMENSION(3),INTENT(IN):: DIP
      DOUBLE PRECISION,INTENT(IN):: DEL
!     VARIABLES      
      DOUBLE PRECISION:: DELI
      INTEGER:: I
      DOUBLE PRECISION,PARAMETER ::ONE=1.0D+00,ZERO=0.0D+00
!      
      DELI = ONE/DEL
      IF(NP <= 0) THEN
!      INITIALIZE DIPOLE DERIVATIVE MATRIX
       DEQ(1) = DIP(1)
       DEQ(2) = DIP(2)
       DEQ(3) = DIP(3)
       DDM = ZERO
       RETURN
      ENDIF
!     UPDATE DIPOLE DERIVATIVE MATRIX
      IF(IVIB .EQ. 2) THEN
       DO I=1,3
        DDM(I,NP) = (DDM(I,NP) - DIP(I))*DELI
       ENDDO
       RETURN
      ENDIF
      IF(NVIB .EQ. 2) THEN
       DO I=1,3
        DDM(I,NP) = DIP(I)
       ENDDO
       RETURN
      ENDIF
      DO I=1,3
       DDM(I,NP) = (DIP(I) - DEQ(I))*DELI
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! SYMDDMd
      SUBROUTINE SYMDDMd(DDM,NAT)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NAT
      DOUBLE PRECISION,DIMENSION(9,NAT),INTENT(INOUT):: DDM
      INTEGER:: I,J
      DOUBLE PRECISION,PARAMETER ::ZERO=0.0D+00, TM7=1.0D-07
!     PURIFY DIPOLE DERIVATIVE TENSOR
      DO I=1,NAT
       DO J=1,9
        IF(DABS(DDM(J,I)).LT.TM7)DDM(J,I)=ZERO
       ENDDO
      ENDDO
      RETURN
      END

! FGMTRXd
      SUBROUTINE FGMTRXd(Cxyz,HESS,GRAD,NC1,ZNUC,ZMASS,DDM,IPROJECT,    &
                         ISIGMA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL EFIELDL,FROZEN
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INP_EFIELDL/EFX,EFY,EFZ,EFIELDL
      COMMON/INPNOF_MOLDEN/MOLDEN
      COMMON/INPNOF_INICON/INICOND
      COMMON/INPNOF_FROZEN/FROZEN,IFROZEN(200)
      COMMON/ECP2/CLP(4004),ZLP(4004),NLP(4004),KFRST(1001,6),          &
                  KLAST(1001,6),LMAX(1001),LPSKIP(1001),IZCORE(1001)
!
      CHARACTER*8 LETI,IBLANK
      CHARACTER*4 CLAB(3)
      DOUBLE PRECISION :: LAB(9)
      DOUBLE PRECISION,DIMENSION(NC1),INTENT(IN) :: Cxyz !NC1=3*NATOMS
      DOUBLE PRECISION,DIMENSION(NC1,NC1),INTENT(INOUT) :: HESS   
      DOUBLE PRECISION,DIMENSION(3,NATOMS),INTENT(IN) :: GRAD
      DOUBLE PRECISION,DIMENSION(3,NC1) :: SVTZR,SVTZT
      DOUBLE PRECISION,DIMENSION(NC1) :: SVTZTT,SVTZRT
      DOUBLE PRECISION,DIMENSION(NATOMS),INTENT(IN) :: ZNUC
      DOUBLE PRECISION,DIMENSION(NATOMS),INTENT(IN) :: ZMASS
      DOUBLE PRECISION,DIMENSION(3,NC1),INTENT(INOUT)::DDM
!     VARIABLES
      LOGICAL :: PRJGRD,PRJROT
      DOUBLE PRECISION :: ZMASST
      DOUBLE PRECISION,DIMENSION(3) :: CMASS,VMOI
      DOUBLE PRECISION,DIMENSION(NC1) :: FREQ,RM,E
      DOUBLE PRECISION,DIMENSION(3,NATOMS) :: COM
      DOUBLE PRECISION,DIMENSION(NC1,8) :: SCR
      INTEGER :: NIMAG,NLAST,NSKIP                                      
      DOUBLE PRECISION,PARAMETER :: ZERO=0.0D+00, ONE=1.0D+00 
      DOUBLE PRECISION,PARAMETER :: TWO=2.0D+00, THREE=3.0D+00
      DOUBLE PRECISION,PARAMETER :: FOUR=4.0D+00, SIX=6.0D+00
      DOUBLE PRECISION,PARAMETER :: SEVEN=7.0D+00
      DOUBLE PRECISION,PARAMETER :: TFACT=2.642461D+07
      DOUBLE PRECISION,PARAMETER :: AUtoCM1=219474.6314D+00
      DOUBLE PRECISION,PARAMETER :: PLANCK=6.62607015D-34
      DOUBLE PRECISION,PARAMETER :: BOHR = 5.2917720859D-11
      DOUBLE PRECISION,PARAMETER :: AVOGAD=6.02214076D+23 
      DOUBLE PRECISION,PARAMETER :: AUtoANG=0.52917720859d0
      DOUBLE PRECISION,PARAMETER :: AMUtoAU=1822.888485540950d0
      DOUBLE PRECISION,PARAMETER :: FStoAU=41.341373336561364d0
      DOUBLE PRECISION,PARAMETER :: AUtoKCAL=627.509391d0
      DOUBLE PRECISION,PARAMETER :: CATOM = 12.011D+00
      CHARACTER*4,DIMENSION(106) :: ATMLAB
      DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: VECt,FREQau
      DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE :: VEC
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
      DATA CLAB /'   X','   Y','   Z'/
      DATA LETI,IBLANK/' I','  '/
!
      DOUBLE PRECISION,PARAMETER::DFAC=2.54174D0    ! Debye
!-----------------------------------------------------------------------
!     WILSON -FG- MATRIX VIBRATIONAL ANALYSIS
!     W.D.GWINN   J.CHEM.PHYS.  55, 477-481 (1971)
!-----------------------------------------------------------------------
      ALLOCATE(VEC(NC1,NC1))
      NIMAG = 0
      WRITE(11,9000)
!     ZMASS=MASSES, RM=TRIPLES OF INVERSE SQUARE ROOTS OF MASSES.
      WRITE(11,9010)
      DO IAT = 1,NATOMS
       WRITE(11,9020) IAT,ZMASS(IAT)
      ENDDO
      I0=0
      DO I=1,NATOMS
       IF(ZMASS(I)>ZERO) THEN
        DMY = ONE/DSQRT(ZMASS(I))
       ELSE
        DMY = ONE
       END IF
       RM(I0+1) = DMY
       RM(I0+2) = DMY
       RM(I0+3) = DMY
       I0=I0+3
      ENDDO
!     ----- GENERATE MASS WEIGHTED GRADIENT AND HESSIAN
      I0 = 1
      DO I=1,NATOMS
       DO J=1,3
        SVTZT(J,I) = GRAD(J,I) 
        SVTZT(J,I) = RM(I0) * SVTZT(J,I)
       ENDDO
       I0 = I0+3
      ENDDO
      CALL HESMASd(NC1,HESS,RM)
!
!     ----- GET CENTER OF MASS IN MASS-WEIGHTED CARTESIAN COORDS
      CALL CENMASd(NATOMS,Cxyz,COM,ZMASST,CMASS,ZMASS)
!
!     PROJECT THE FORCE CONSTANT MATRIX (PROJECT=T)
      IF (IPROJECT==1) THEN
       WRITE(11,'(1X)')
       WRITE(11,9030)
       PRJROT = .TRUE.
       IF(EFIELDL) PRJROT = .FALSE.
       PRJGRD = .FALSE.
       CALL PRJFCMd(PRJGRD,PRJROT,ZMASST,HESS,COM,SVTZT,RM,NATOMS,NC1)
      END IF
!
!     ----- GET NORMAL MODES AND FREQUENCIES
      CALL DIAG(NC1,HESS,VEC,E,SCR(:,1))
      CALL STFASEd(VEC,NC1,NC1,NC1)
!
!     ----- TRANSLATIONAL AND ROTATIONAL SAYVETZ CONDITIONS
      SVTZR = ZERO
      SVTZT = ZERO
      DO I=1,NC1
       DO J=1,NATOMS
        JJ=MAX(3*(J-1),6*(J-1)-3*NATOMS)
        AMASS=ONE/RM((J-1)*3+1)
        DO K=1,3
         K1=MOD(K+1,4)+(K+1)/4
         K2=MOD(K+2,4)+(K+2)/4
         SVTZT(K,I)=SVTZT(K,I)+AMASS*VEC(JJ+K,I)
         SVTZR(K,I)=SVTZR(K,I)+COM(K1,J)*VEC(JJ+K2,I)                   &
                              -COM(K2,J)*VEC(JJ+K1,I)
        ENDDO
       ENDDO
       SVTZTT(I)=DSQRT(SVTZT(1,I)**2+SVTZT(2,I)**2+SVTZT(3,I)**2)
       SVTZRT(I)=DSQRT(SVTZR(1,I)**2+SVTZR(2,I)**2+SVTZR(3,I)**2)
      ENDDO
!
!     ----- CONVERT NORMAL MODE DISPLACEMENTS -----
      DO I = 1,NC1
       DO J = 1,NC1
        VEC(J,I) = VEC(J,I)*RM(J)
       END DO
      END DO
!
!     COMPUTE REDUCED MASS
      DO J = 1,NC1
       DD = DDOT(NC1,VEC(1,J),1,VEC(1,J),1)
       SCR(J,4) = ONE/DD
      END DO
!
!     COUNT NUMBER OF NEGATIVE EIGENVALUES
      DO I = 1,NC1
       IF(E(I)<ZERO)NIMAG=NIMAG+1
      ENDDO
!
!     DECIDE WHICH MODES ARE NOT TRUE VIBRATIONS
      NSKIP = 6
      IF(NATOMS==2)NSKIP = 5
      DO IAT=1,NATOMS
       NUCZ = INT(ZNUC(IAT))
       IF(NUCZ.EQ.0) NSKIP=NSKIP+3
      ENDDO
      IF(IPROJECT==1 .and. PRJGRD) NSKIP = NSKIP + 1
!
      NLAST = NSKIP
      DO I=1,NIMAG
       N2OF3 = 0
       MAYBE = I+NSKIP
       IF(MAYBE>NC1) EXIT
       IF(DABS(E(I))>DABS(E(MAYBE))) N2OF3=N2OF3+1
       IF(IPROJECT==1) THEN
        IF(N2OF3.EQ.1) NLAST = NLAST + 1
       ELSE
        IF(0.01D+00+SVTZRT(I)<SVTZRT(MAYBE)) N2OF3=N2OF3+1
        IF(0.01D+00+SVTZTT(I)<SVTZTT(MAYBE)) N2OF3=N2OF3+1
        IF(N2OF3.GE.2) NLAST = NLAST + 1
        IF(N2OF3.EQ.2) WRITE(11,9040) I,MAYBE
        IF(N2OF3.EQ.1) WRITE(11,9040) MAYBE,I
        IF(N2OF3.LE.1) EXIT
       END IF
      ENDDO
!
      NFIRST = NLAST - NSKIP + 1
      NIMAG = NFIRST-1
      WRITE(11,9050) NFIRST,NLAST
!
!     ----- PRINT MESSAGE FOR PARTIAL HESSIAN ANALYSIS -----
      IF(FROZEN) THEN
       NFRZ = 0
       DO I=1,200,2
        IF(IFROZEN(I).EQ.0) EXIT
        NFRZ = NFRZ + 1
       ENDDO      
       IF(NFRZ.GE.6)THEN
        IF(NFRZ.GE.9)THEN
         WRITE(11,9051)NLAST+1, NLAST+NFRZ-6, NLAST+NFRZ-5,             &
         NLAST+NFRZ-3, NLAST+1, NLAST+NFRZ-3                             
        ELSE                                                             
         WRITE(11,9052) NLAST+1, NLAST+NFRZ-3,                          &
         NLAST+1, NLAST+NFRZ-3
        END IF
       END IF
      ENDIF
!
!     ----- CONVERT FREQUENCIES TO WAVENUMBERS -----
      DO I = 1,NC1
       FREQ(I) = DSQRT(DABS(TFACT*E(I)))
      ENDDO
!
!     ----- PRINT WARNING FOR SKIPPING FROZEN FREQUENCIES -----
      IF(FROZEN) THEN
       IF(NFRZ.GE.6)THEN
        IF(FREQ(NLAST+NFRZ-3)>12.AND.FREQ(NLAST+NFRZ-2)<30)             &
        WRITE(11,9053)NLAST+NFRZ-3,NLAST+NFRZ-2
       END IF
      ENDIF
!
!     ----- COMPUTE IR INTENSITIES -----
!     PROJECT THE DIPOLE DERIVATIVE TENSOR ONTO EACH NORMAL MODE,
!     AND TAKE THE SQUARE OF THE NORM OF THIS 3 COMPONENT VECTOR
      DDM(:,:) = DDM(:,:)*DFAC
      DO J = 1,NC1
       DDX = DDOT(NC1,VEC(1,J),1,DDM(1,1),3)
       DDY = DDOT(NC1,VEC(1,J),1,DDM(2,1),3)
       DDZ = DDOT(NC1,VEC(1,J),1,DDM(3,1),3)
       SCR(J,1) = DDX*DDX + DDY*DDY + DDZ*DDZ
      ENDDO
!
!     PRINT OUT UNITS INFO
      WRITE(11,9055)
      MAXCOL = 0
      INCR = 5
!      
!     ----- PRINT THE FREQUENCY AND INTENSITY ----- 
!
      DO
       MINCOL = MAXCOL+1
       MAXCOL = MAXCOL+INCR
       IF (MAXCOL > NC1) MAXCOL = NC1
       WRITE (11,9060)
       WRITE (11,9100) (J,J = MINCOL,MAXCOL)
       DO J=MINCOL,MAXCOL
        JJ = J + 1 - MINCOL
        LAB(JJ) = transfer (LETI,LAB(JJ))            ! LAB(JJ)=LETI
        IF(J>NIMAG)LAB(JJ)=transfer(IBLANK,LAB(JJ))  ! LAB(JJ)=IBLANK
       ENDDO
       WRITE(11,9110)(FREQ(J),LAB(J+1-MINCOL),J = MINCOL,MAXCOL)
       WRITE(11,9120)(SCR(J,1),J=MINCOL,MAXCOL)
       WRITE(11,9115)(SCR(J,4),J=MINCOL,MAXCOL)
       WRITE(11,9060)
!      
!      PRINT AB INITIO NORMAL MODE COMPONENTS
       DO IAT = 1,NATOMS
        I0 = 3*(IAT-1)
        WRITE(11,9150)IAT,                                              &
                      CLAB(1),(VEC(I0+1,J),J=MINCOL,MAXCOL)
        WRITE(11,9160)CLAB(2),(VEC(I0+2,J),J=MINCOL,MAXCOL)
        WRITE(11,9160)CLAB(3),(VEC(I0+3,J),J=MINCOL,MAXCOL)
       ENDDO
!      
!      PRINT SAYVETZ CONDITIONS
       WRITE(11,9060)
       WRITE(11,9180)CLAB(1),(SVTZT(1,I),I=MINCOL,MAXCOL)
       WRITE(11,9160)CLAB(2),(SVTZT(2,I),I=MINCOL,MAXCOL)
       WRITE(11,9160)CLAB(3),(SVTZT(3,I),I=MINCOL,MAXCOL)
       WRITE(11,9200)(SVTZTT(I),I=MINCOL,MAXCOL)
       WRITE(11,9060)
       WRITE(11,9190)CLAB(1),(SVTZR(1,I),I=MINCOL,MAXCOL)
       WRITE(11,9160)CLAB(2),(SVTZR(2,I),I=MINCOL,MAXCOL)
       WRITE(11,9160)CLAB(3),(SVTZR(3,I),I=MINCOL,MAXCOL)
       WRITE(11,9200)(SVTZRT(I),I=MINCOL,MAXCOL)
       IF (MAXCOL>NC1 .OR. MAXCOL.EQ.NC1) EXIT
      ENDDO
      WRITE(11,9220)
      WRITE(11,9060)
      WRITE(11,9070)
!
!     Frequencies and corresponding normal coordinates for MOLDEN
!
       IF(MOLDEN==1)THEN
        WRITE(18,'(A6)')'[FREQ]'
        do j=1,nc1
         write(18,'(F10.2)')FREQ(j)
        enddo
        WRITE(18,'(A5)')'[INT]'
        do j=1,nc1
         write(18,'(F10.5)')SCR(j,1)
        enddo
        WRITE(18,'(A10)')'[FR-COORD]'
        do iat=1,NATOMS
         IZNUC = INT(ZNUC(iat))+IZCORE(iat)      
         io = 3*(iat-1)
         WRITE(18,'(1X,A4,3F15.4)')ATMLAB(IZNUC),                       &
                 Cxyz(1+io),Cxyz(2+io),Cxyz(3+io)         
        enddo
        WRITE(18,'(A15)')'[FR-NORM-COORD]'
        do j=1,nc1
         write(18,'(A9,1X,I5)')'vibration',j
         do iat = 1,NATOMS         
          io = 3*(iat-1)
          write(18,'(3F12.8)')VEC(1+io,j),VEC(2+io,j),VEC(3+io,j)
         enddo 
        enddo
       ENDIF
!
!      Initial Conditions according to normal modes
!      Equilibrium Geometry & Zero Point Energy Velocities (Ang/fs)
!
       IF(INICOND==1)THEN
        WRITE(11,9080)
        ALLOCATE(VECt(NC1),FREQau(NC1))       
        OPEN(33,FILE='ini.xyz',STATUS='UNKNOWN',ACTION='WRITE',         &
                FORM='FORMATTED',ACCESS='SEQUENTIAL')
!
!       v = sqrt [(h/2*pi) * (w/m) ]
!
        FREQau = DSQRT(FREQ/AUtoCM1)
        FREQau = FREQau * AUtoANG*FStoAU/DSQRT(AMUtoAU) !*0.512396293759
!
        VECt = 0.0d0

        IF(NATOMS>2)THEN
         do iat = 1,NATOMS
          io = 3*(iat-1)
          do j=7,nc1
           VECt(1+io) = VECt(1+io) + VEC(1+io,j)*FREQau(j)
           VECt(2+io) = VECt(2+io) + VEC(2+io,j)*FREQau(j)
           VECt(3+io) = VECt(3+io) + VEC(3+io,j)*FREQau(j)
          enddo
         enddo
        ELSE IF(NATOMS==2)THEN
         do iat = 1,NATOMS
          io = 3*(iat-1)
          VECt(1+io) = VECt(1+io) + VEC(1+io,6)*FREQau(6)
          VECt(2+io) = VECt(2+io) + VEC(2+io,6)*FREQau(6)
          VECt(3+io) = VECt(3+io) + VEC(3+io,6)*FREQau(6)
         enddo
        END IF
!
        ZPE = 0.0d0
        write(33,*) NATOMS
        write(33,'(f19.12)') 0.0
        do iat=1,NATOMS
         IZNUC = INT(ZNUC(iat))+IZCORE(iat)      
         io = 3*(iat-1)         
         write(33,'(a4, 6(1x,f14.8))')ATMLAB(IZNUC),                    &
              Cxyz(1+io)*AUtoANG,Cxyz(2+io)*AUtoANG,Cxyz(3+io)*AUtoANG, &
              VECt(1+io),VECt(2+io),VECt(3+io)
         write(11,'(a4, 6(1x,f14.8))')ATMLAB(IZNUC),                    &
              Cxyz(1+io)*AUtoANG,Cxyz(2+io)*AUtoANG,Cxyz(3+io)*AUtoANG, &
              VECt(1+io),VECt(2+io),VECt(3+io)
         ZPE = ZPE + ZMASS(iat) * ( VECt(1+io)*VECt(1+io)               &
                   + VECt(2+io)*VECt(2+io) + VECt(3+io)*VECt(3+io) )
        enddo       
        ZPE = 0.5d0*AMUtoAU*ZPE/(AUtoANG*FStoAU*AUtoANG*FStoAU)
        write(11,9090)CMASS(1)*AUtoANG,CMASS(2)*AUtoANG,CMASS(3)*AUtoANG
        write(11,9091)ZPE,ZPE*AUtoKCAL,ZPE*AUtoCM1        
!        
!       The excitation of a specific vibrational mode can be obtained as
!       Vj(k) = SUM_i {sqrt[(2k+1)*VEC^i_j(0)}
!
!jk     jk = 7  ! number of the excited mode (7 <= jk <= nc1)
!jk     VECt = 0.0d0
!jk     do iat = 1,NATOMS         
!jk      io = 3*(iat-1)
!jk      do j=7,nc1
!jk       if(j==jk)then 
!jk        VECt(1+io) = VECt(1+io)+VEC(1+io,j)*FREQau(j)*SQRT(5.0) !K=2
!jk        VECt(2+io) = VECt(2+io)+VEC(2+io,j)*FREQau(j)*SQRT(5.0) !K=2
!jk        VECt(3+io) = VECt(3+io)+VEC(3+io,j)*FREQau(j)*SQRT(5.0) !K=2
!jk       else
!jk        VECt(1+io) = VECt(1+io)+VEC(1+io,j)*FREQau(j) !*SQRT(1)  K=0
!jk        VECt(2+io) = VECt(2+io)+VEC(2+io,j)*FREQau(j) !*SQRT(1)  K=0
!jk        VECt(3+io) = VECt(3+io)+VEC(3+io,j)*FREQau(j) !*SQRT(1)  K=0
!jk       endif 
!jk      enddo
!jk     enddo
!jk
!jk     VibE = 0.0d0
!jk     write(33,*) NATOMS
!jk     write(33,'(f19.12)') 0.0
!jk     do iat=1,NATOMS
!jk      IZNUC = INT(ZNUC(iat))+IZCORE(iat)      
!jk      io = 3*(iat-1)         
!jk      write(33,'(a4, 6(1x,f14.8))')ATMLAB(IZNUC),                    &
!jk           Cxyz(1+io)*AUtoANG,Cxyz(2+io)*AUtoANG,Cxyz(3+io)*AUtoANG, &
!jk           VECt(1+io),VECt(2+io),VECt(3+io)
!jk      write(11,'(a4, 6(1x,f14.8))')ATMLAB(IZNUC),                    &
!jk           Cxyz(1+io)*AUtoANG,Cxyz(2+io)*AUtoANG,Cxyz(3+io)*AUtoANG, &
!jk           VECt(1+io),VECt(2+io),VECt(3+io)
!jk      VibE = VibE + ZMASS(iat) * ( VECt(1+io)*VECt(1+io)             &
!jk                  + VECt(2+io)*VECt(2+io) + VECt(3+io)*VECt(3+io) )
!jk     enddo       
!jk     VibE = 0.5d0*AMUtoAU*VibE/(AUtoANG*FStoAU*AUtoANG*FStoAU)
!jk     write(11,9092)VibE,VibE*AUtoKCAL,VibE*AUtoCM1        
!
        DEALLOCATE(VECt,FREQau)
        CLOSE(33)
       ENDIF       
!
!     ----- THERMOCHEMISTRY ANALYSIS ----- 
!
      CALL INRTIA(Cxyz,COM,ZMASS,VMOI,NATOMS)
      WRITE(11,9300) (VMOI(I),I=1,3),ONE
      FACT1 = (CATOM*BOHR*BOHR)/(12.0D+03*AVOGAD)
      PI = ACOS(-ONE)
      FACT2 = PLANCK/(8.0D+09*PI*PI)
      ACONST = ZERO
      BCONST = ZERO
      CCONST = ZERO
      IF(VMOI(1).GT.0.001D+00) ACONST = FACT2/(FACT1*VMOI(1))
      IF(VMOI(2).GT.0.001D+00) BCONST = FACT2/(FACT1*VMOI(2))
      IF(VMOI(3).GT.0.001D+00) CCONST = FACT2/(FACT1*VMOI(3))
      WRITE(11,9310) ACONST, BCONST, CCONST
!
      FACT = AUtoANG*AUtoANG/AVOGAD
      VMOI(1) = FACT * VMOI(1)
      VMOI(2) = FACT * VMOI(2)
      VMOI(3) = FACT * VMOI(3)
!
      NROTRA=NLAST
      CALL THERMO(NC1,NROTRA,FREQ,Cxyz,ZMASS,ISIGMA)
!
      DEALLOCATE(VEC)
      RETURN
!-----------------------------------------------------------------------      
 9000 FORMAT(/1X,                                                       &
       'Normal Coordinate Analysis in the Harmonic Approximation'/      &
       1X,56(1H-))                                                       
 9010 FORMAT(/10X,'ATOMIC WEIGHTS (AMU)'/)                               
 9020 FORMAT(I5,5X,F15.5)                                                
 9030 FORMAT(1X,'The Force Constant Matrix is projected to eliminate'   &
            /1X,'rotational and vibrational contaminants')               
 9040 FORMAT(/1X,'* * * WARNING, MODE',I2,' HAS BEEN CHOSEN AS A ',     &
                'VIBRATION'/10X,'WHILE MODE',I2,                        &
                ' IS ASSUMED TO BE A TRANSLATION/ROTATION.'/            &
             1X,'PLEASE VERIFY THE PROGRAM''S DECISION MANUALLY!'/)      
 9050 FORMAT(/1X,'MODES',I2,' TO',I2,' ARE TAKEN AS ROTATIONS',         &
             ' AND TRANSLATIONS.')                                       
 9051 FORMAT(/1X,'MODES',I3,' TO',I3,' ARE INTERNAL ',                  &
                'VIBRATIONS ',                                          &
                'OF FROZEN ATOMS.',                                     &
                /1X,'MODES',I3,' TO',I3,' ARE RELATIVE ',               &
                'VIBRATIONS ',                                          &
                'BETWEEN FROZEN AND UNFROZEN ATOMS.',                   &
                /1X,'MODES',I3,' TO',I3,' DO NOT CONTRIBUTE TO ',       &
                'VIBRATIONAL PARTITION FUNCTION',                       &
                ' AND ENERGY.',/)                                        
 9052 FORMAT(/1X,'MODES',I3,' TO',I3,' ARE RELATIVE ',                  &
                'VIBRATIONS ',                                          &
                'BETWEEN FROZEN AND UNFROZEN ATOMS.',                   &
                /1X,'MODES',I3,' TO',I3,' DO NOT CONTRIBUTE TO ',       &
                'VIBRATIONAL ',                                         &
                'ENERGIES AND PARTITION FUNCTION.',/)                    
 9053 FORMAT(//1X,'* * * WARNING ! * * *',//                            &
                1X,'MODE',I3,' IS TAKEN AS A RELATIVE VIBRATION ',      &
                'BETWEEN FROZEN AND UNFROZEN ATOMS ',                   &
                /1X,'WHILE MODE',I3,' IS TAKEN AS AN INTERNAL ',        &
                'VIBRATION OF UNFROZEN ATOMS.',//                       &
                1X,'PLEASE VERIFY THE PROGRAM''S DECISION MANUALLY !'/)  
 9055 FORMAT(/1X,'FREQUENCIES IN CM**-1',                               &
             /1X,'IR INTENSITIES IN DEBYE**2/AMU-ANGSTROM**2',          &
             /1X,'REDUCED MASSES IN AMU.')                               
 9060 FORMAT(1X)                                                         
 9070 FORMAT(5X,'END OF NORMAL MODES',/)
 9080 FORMAT(/1X,'INITIAL CONDITIONS ACCORDING TO NORMAL MODES',        &
             //17X,'Equilibrium Geometry (Ang)',                        &
               13X,'Zero Point Energy Velocities (Ang/fs)',/)
 9090 FORMAT(/,'CM',2X,3(1x,f14.8))                                
 9091 FORMAT(/1X,'Zero Point Energy =',F10.4,' a.u. (',                 &
             F7.2,' kcal/mol,',F9.1,' cm-1 )'/)
 9092 FORMAT(/1X,'Vibrational Energy =',F10.4,' a.u. (',                &
             F7.2,' kcal/mol,',F9.1,' cm-1 )'/)
 9100 FORMAT(20X,9(4X,I3,5X))                                            
 9110 FORMAT(1X,'      FREQUENCY:',3X,9(F10.2,A2))                       
 9115 FORMAT(1X,'   REDUCED MASS:',3X,9(F10.5,2X))                       
 9120 FORMAT(1X,'   IR INTENSITY:',3X,9(F10.5,2X))                       
 9150 FORMAT(I3,13X,A4,9F12.8)                                           
 9160 FORMAT(16X,A4,9F12.8)                                              
 9180 FORMAT(16H TRANS. SAYVETZ ,A4,9F12.8)                              
 9190 FORMAT(16H   ROT. SAYVETZ ,A4,9F12.8)                              
 9200 FORMAT(15X,5HTOTAL,9F12.8)                                         
 9220 FORMAT(/1X,'NOTE - THE MODES J,K ARE ORTHONORMALIZED',            &
                 ' ACCORDING TO'/                                       &
              1X,'SUM ON I   M(I) * (X(I,J)*X(I,K) + Y(I,J)*',          &
                 'Y(I,K) + Z(I,J)*Z(I,K)) = DELTA(J,K)')                 
 9300 FORMAT(1X,'THE MOMENTS OF INERTIA ARE (IN AMU*BOHR**2)'/          &
             1X,3F12.5/                                                 &
             1X,'THE ROTATIONAL SYMMETRY NUMBER IS',F5.1)                
 9310 FORMAT(1X,'THE ROTATIONAL CONSTANTS ARE (IN GHZ)',/               &
             1X,3F12.5)     
      END

! HESMASd      
      SUBROUTINE HESMASd(NCOORD,HESS,RTRMS)
      IMPLICIT NONE
      INTEGER,INTENT(IN)::NCOORD
      DOUBLE PRECISION,DIMENSION(NCOORD,NCOORD),INTENT(INOUT)::HESS
      DOUBLE PRECISION,DIMENSION(NCOORD),INTENT(IN)::RTRMS
      DOUBLE PRECISION,DIMENSION((NCOORD**2+NCOORD)/2):: A
      DOUBLE PRECISION:: RTRMSI
      INTEGER:: I,J,IJ
!-----------------------------------------------------------------------      
      IJ = 0
      DO I=1,NCOORD
       RTRMSI = RTRMS(I)
       DO J=1,I
        IJ = IJ + 1
        A(IJ) = RTRMSI * HESS(I,J) * RTRMS(J)
       END DO
      END DO
!     Symmetrize Hessian
      IJ=0
      DO I=1,NCOORD
       DO J=1,I
        IJ=IJ+1
        HESS(I,J)=A(IJ)
        HESS(J,I)=A(IJ)
       END DO
      END DO
!-----------------------------------------------------------------------      
      RETURN
      END

! CENMASd      
      SUBROUTINE CENMASd(NAT,C,COM,ZMASST,CMASS,ZMASS)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NAT
      DOUBLE PRECISION,INTENT(OUT):: ZMASST
      DOUBLE PRECISION,DIMENSION(3,NAT),INTENT(IN)::C
      DOUBLE PRECISION,DIMENSION(3,NAT),INTENT(OUT)::COM
      DOUBLE PRECISION,DIMENSION(3),INTENT(OUT)::CMASS
      DOUBLE PRECISION,DIMENSION(NAT),INTENT(IN)::ZMASS
      DOUBLE PRECISION:: AMASS
      INTEGER:: I,J
      DOUBLE PRECISION,PARAMETER::ZERO=0.0D+00,ONE=1.0D+00
!-----------------------------------------------------------------------      
!     CALCULATE TOTAL MASS AND CENTER OF MASS
!-----------------------------------------------------------------------
      ZMASST= ZERO
      CMASS = ZERO
!
      DO I=1,NAT
       AMASS=ZMASS(I)
       ZMASST=ZMASST+AMASS
       DO J=1,3
        CMASS(J)=CMASS(J)+AMASS*C(J,I)
       END DO
      END DO
!
      DO I=1,3
       CMASS(I)=CMASS(I)/ZMASST
      END DO
!
      DO I=1,NAT
       IF (ZMASS(I)>ZERO) THEN
        AMASS = DSQRT(ZMASS(I))
       ELSE
        AMASS = ONE
       END IF
       DO J=1,3
        COM(J,I) = AMASS * (C(J,I)-CMASS(J))
       END DO
      END DO
!-----------------------------------------------------------------------      
      RETURN
      END   

! PRJFCMd                                                              
      SUBROUTINE PRJFCMd(PRJGRD,PRJROT,TOTM,FCM,X,DX,RM,NATM,NCC)        
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!     ARGUMENTS
      LOGICAL,INTENT(IN):: PRJGRD,PRJROT
      INTEGER,INTENT(IN) :: NATM,NCC
      DOUBLE PRECISION,INTENT(IN) :: TOTM
      DOUBLE PRECISION,DIMENSION(NCC),INTENT(IN) :: X,RM
      DOUBLE PRECISION,DIMENSION(NCC),INTENT(INOUT) :: DX
      DOUBLE PRECISION,DIMENSION(NCC,NCC),INTENT(INOUT) :: FCM
!     VARIABLES 
      DOUBLE PRECISION,DIMENSION(NCC,NCC)::P,BUF
      DOUBLE PRECISION,DIMENSION(3,3,3)::TENS
      DOUBLE PRECISION,DIMENSION(3,3)::ROT,SCR
      INTEGER,DIMENSION(6)::ISCR
      DOUBLE PRECISION,DIMENSION(2)::DETERM
      DOUBLE PRECISION,PARAMETER::ZERO=0.0D+00,CUT8=1.0D-08
      DOUBLE PRECISION,PARAMETER::HALF=0.5D+00,ONE=1.0D+00
!         -TENS- IS "THE USUAL TOTALLY ASYMMETRIC CARTESIAN TENSOR"
      DATA TENS/ 0.0D+00,  0.0D+00,  0.0D+00, &  ! X
                 0.0D+00,  0.0D+00, -1.0D+00, &  ! X
                 0.0D+00,  1.0D+00,  0.0D+00, &  ! Y
                 0.0D+00,  0.0D+00,  1.0D+00, &  ! Y
                 0.0D+00,  0.0D+00,  0.0D+00, &  ! Y
                -1.0D+00,  0.0D+00,  0.0D+00, &  ! Z
                 0.0D+00, -1.0D+00,  0.0D+00, &  ! Z
                 1.0D+00,  0.0D+00,  0.0D+00, &  ! Z
                 0.0D+00,  0.0D+00,  0.0D+00  /
!
!     PROJECTION OF THE FORCE CONSTANT MATRIX AT AN EQUILIBRIUM
!     GEOMETRY IS DONE WITH AN IDENTICALLY ZERO GRADIENT VECTOR
!     SO THAT THE TRANSLATIONAL AND VIBRATIONAL CONTAMINANTS
!     ONLY ARE ELIMINATED.
!
!     PROJECTION OF THE F.C.M. AT A POINT WITH A NON-ZERO GRADIENT
!     ALSO ELIMINATES THESE CONTAMINANTS, BUT ALSO PROJECTS THE
!     MATRIX SO THAT ONE OF ITS NORMAL MODES LIES PARALLEL TO
!     THE MASS-WEIGHTED GRADIENT.  THIS MODE WILL HAVE A ZERO
!     FREQUENCY, AND THE OTHER 3N-7 MODES WILL BE ORTHOGONAL TO
!     THIS MODE.  THIS TYPE OF PROJECTION IS USEFUL IN DYNAMICS.
!
!     FOR A DESCRIPTION OF THE METHOD, SEE W.H.MILLER, N.C.HANDY,
!     J.E.ADAMS, IN J.CHEM.PHYS. 72, 99-112(1980).
!
!     PRJGRD    = IF .TRUE., THE PROJECTION WILL ZERO THE NORMAL
!                 MODE PARALLEL TO THE GRADIENT VECTOR.  USEFUL
!                 IN TRANSITION STATE THEORY.
!     PRJROT    = IF .TRUE., THE PROJECTION WILL ZERO THE NORMAL
!                 MODES CORRESPONDING TO ROTATIONS
!     TOTM      = TOTAL MASS OF SYSTEM (INCLUDING FRAGMENT)
!     FCM       = ON ENTRY, MASS-WEIGHTED FORCE CONSTANT MATRIX.
!                 ON EXIT, PROJECTED TO REMOVE T,R AND POSSIBLY
!                 THE GRADIENT'S DEGREE OF FREEDOM
!     X         = MASS-WEIGHTED COORDINATE
!     DX        = ON ENTRY, MASS-WEIGHTED GRADIENT VECTOR.
!                 ON EXIT, THIS WILL HAVE BEEN DESTROYED
!                 THIS IS NOT USED IF PRJGRD IS FALSE.
!     RM        = INVERSE OF THE SQUARE ROOT MASS.
!     P,BUF     = WORK BUFFERS (NCC*NCC)
!     NATM      = NUMBER OF ATOMS.
!     NCC       = NO OF CARTESIAN COORDINATES (3*NATM)
!
!     ----- NORMALIZE THE GRADIENT -----
      IF(PRJGRD) THEN
       GNORM = DDOT(NCC,DX,1,DX,1)
       GNORM = ONE/DSQRT(GNORM)
       CALL DSCAL(NCC,GNORM,DX,1)
      END IF
!      
!  ----- CALCULATE PROJECTED FORCE CONSTANT MATRIX -----
!  COPIED FROM POLYRATE BY S. KOSEKI, DEC. 2, 1985.
!
!       --- COMPUTE INERTIA TENSOR ---
!
      ROT = ZERO
      IF(.NOT.PRJROT) GO TO 200
!
      DO I=1,NATM
       L=3*(I-1)+1
       ROT(1,1)=ROT(1,1)+X(L+1)**2+X(L+2)**2
       ROT(1,2)=ROT(1,2)-X(L)*X(L+1)
       ROT(1,3)=ROT(1,3)-X(L)*X(L+2)
       ROT(2,2)=ROT(2,2)+X(L)**2+X(L+2)**2
       ROT(2,3)=ROT(2,3)-X(L+1)*X(L+2)
       ROT(3,3)=ROT(3,3)+X(L)**2+X(L+1)**2
      ENDDO
      ROT(2,1)=ROT(1,2)
      ROT(3,1)=ROT(1,3)
      ROT(3,2)=ROT(2,3)
!
      CHK=ROT(1,1)*ROT(2,2)*ROT(3,3)
      IF(DABS(CHK).GT.CUT8) GO TO 180
!
!       --- COMPUTE INVERSE IN CASE OF SOME ZERO MOMENTS ---
!
      IF(DABS(ROT(1,1)).GT.CUT8) GO TO 130
      IF(DABS(ROT(2,2)).GT.CUT8) GO TO 120
      IF(DABS(ROT(3,3)).GT.CUT8) GO TO 110
      WRITE(11,10)ROT(1,1),ROT(2,2),ROT(3,3)
      RETURN
!
!             X,Y=0 BUT Z.NE.0
  110 ROT(3,3)=ONE/ROT(3,3)
      GO TO 170
!
!             Y.NE.0
  120 IF(DABS(ROT(3,3)).GT.CUT8) GO TO 160
!             X,Z=0 BUT Y.NE.0
      ROT(2,2)=ONE/ROT(2,2)
      GO TO 170
!
!             X.NE.0
  130 IF(DABS(ROT(2,2)).GT.CUT8) GO TO 140
      IF(DABS(ROT(3,3)).GT.CUT8) GO TO 150
!
!             Y,Z=0 BUT X.NE.0
      ROT(1,1)=ONE/ROT(1,1)
      GO TO 170
!
!             X,Y.NE.0 BUT Z=0
  140 DET=ROT(1,1)*ROT(2,2)-ROT(1,2)*ROT(2,1)
      TRP=ROT(1,1)
      ROT(1,1)=ROT(2,2)/DET
      ROT(2,2)=TRP/DET
      ROT(1,2)=-ROT(1,2)/DET
      ROT(2,1)=-ROT(2,1)/DET
      GO TO 170
!
!             X,Z.NE.0 BUT Y=0
  150 DET=ROT(1,1)*ROT(3,3)-ROT(1,3)*ROT(3,1)
      TRP=ROT(1,1)
      ROT(1,1)=ROT(3,3)/DET
      ROT(3,3)=TRP/DET
      ROT(1,3)=-ROT(1,3)/DET
      ROT(3,1)=-ROT(3,1)/DET
      GO TO 170
!
!             Y,Z.NE.0 BUT X=0
  160 DET=ROT(3,3)*ROT(2,2)-ROT(3,2)*ROT(2,3)
      TRP=ROT(3,3)
      ROT(3,3)=ROT(2,2)/DET
      ROT(2,2)=TRP/DET
      ROT(3,2)=-ROT(3,2)/DET
      ROT(2,3)=-ROT(2,3)/DET
!
  170 CONTINUE
      GO TO 200
!
!       --- COMPUTE INVERSE FOR ALL MOMENTS NONZERO CASE ---
  180 CONTINUE
      INFO=0
      CALL DGEFA(ROT,3,3,ISCR,INFO)
      IF(INFO.NE.0)THEN
       WRITE(11,20)
       RETURN
      END IF
      CALL DGEDId(ROT,3,3,ISCR,DETERM,SCR,1)
!
  200 CONTINUE
!     COMPUTE P MATRIX by FORMULA (EQN 4.11) OF MILLER, HANDY, ADAMS.  
!     ATOM CENTER
      DO IP=1,NATM
       INDX=3*(IP-1)
       KNDX=MAX(3*(IP-1),6*(IP-1)-3*NATM)
       DO JP=1,IP
        JNDX=3*(JP-1)
        LNDX=MAX(3*(JP-1),6*(JP-1)-3*NATM)
        DO IC=1,3
         II=INDX+IC
         KK=KNDX+IC
         JEND=3
         IF(JP.EQ.IP) JEND=IC
         DO JC=1,JEND
          JJ=JNDX+JC
          LL=LNDX+JC
          XMSUM=ZERO
          IF(PRJROT) THEN
          DO IA=1,3
           DO IB=1,3
            IF(TENS(IA,IB,IC).EQ.ZERO) CYCLE
            DO JA=1,3
             DO JB=1,3
              IF(TENS(JA,JB,JC).EQ.ZERO) CYCLE
              XMSUM = XMSUM + TENS(IA,IB,IC)*TENS(JA,JB,JC)             &
                            * ROT(IA,JA)*X(INDX+IB)*X(JNDX+JB)
             ENDDO
            ENDDO
           ENDDO
          ENDDO
          END IF
          P(KK,LL)=XMSUM
          IF(PRJGRD)P(KK,LL)=P(KK,LL) + DX(KK)*DX(LL)
          IF(IC.EQ.JC)P(KK,LL)=P(KK,LL) + ONE/(RM(II)*RM(JJ)*TOTM)
         ENDDO
        ENDDO
       ENDDO
      ENDDO
!
!        --- COMPUTE (I-P) ---
      DO I=1,NCC
       DO J=1,I
        P(I,J)=-P(I,J)
        IF(I.EQ.J) P(I,J) = ONE + P(I,J)
        IF(DABS(P(I,J)).LT.CUT8) P(I,J)=ZERO
        P(J,I)=P(I,J)
       ENDDO
      ENDDO
!
!        --- PROJECT BY COMPUTING FCM = (I-P)*(FCM*(I-P)) ---
      CALL TFSQUd(BUF,FCM,P,DX,NCC,NCC)
      CALL DCOPY(NCC*NCC,BUF,1,FCM,1)
!
!     ----- ENFORCE SYMMETRY UPON F.C.M. -----
      DO I=1,NCC
       DO J=1,I
        AVRG   = HALF*(FCM(I,J)+FCM(J,I))
        FCM(J,I) = AVRG
        FCM(I,J) = AVRG
       ENDDO
      ENDDO
      RETURN
!-----------------------------------------------------------------------
   10 FORMAT(/1X,'PRJFCM: EVERY DIAGONAL ELEMENT IS ZERO?',3F20.10)
   20 FORMAT(/1X,                                                       &
       'INERTIA MATRIX IS SINGULAR HESSIAN WILL NOT BE PROJECTED!')
      END      

! TFSQUd      
      SUBROUTINE TFSQUd(H,FCM,T,WRK,N,M)                                                        
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: N,M
      DOUBLE PRECISION,DIMENSION(M,M),INTENT(OUT) :: H
      DOUBLE PRECISION,DIMENSION(N),INTENT(OUT) :: WRK
      DOUBLE PRECISION,DIMENSION(N,N),INTENT(IN) :: FCM
      DOUBLE PRECISION,DIMENSION(N,M),INTENT(IN) :: T
      REAL :: DDOT
      INTEGER :: I,II,L,J,IIMAX
      INTEGER,PARAMETER :: MXROWS=5                                                                  
!     ----- TRANSFORM THE SQUARE MATRIX FCM USING VECTORS T -----         
!                      H = T-DAGGER * FCM * T                             
!     THE ORDER OF THE SQUARE MATRICES H, FCM, AND T IS N.                                                                          
      DO I = 1,M,MXROWS                                             
       IIMAX = MIN(M,I+MXROWS-1)                                                                                             
       DO II=I,IIMAX                                              
        DO L=1,N                                                
         WRK(L) = DDOT(N,T(1,II),1,FCM(1,L),1)                      
        ENDDO
        DO J=1,M                                                
         H(II,J) = DDOT(N,WRK,1,T(1,J),1)                         
        ENDDO
       ENDDO
      ENDDO
!-----------------------------------------------------------------------                                                                   
      RETURN                                                            
      END
