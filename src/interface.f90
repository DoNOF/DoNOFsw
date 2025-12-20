!======================================================================!
!                                                                      !
!                   U T I L S   S U B R O U T I N E S                  !
!                                                                      !
!======================================================================!
!                                                                      !
!   AIMMEMrc: Create an input file for Bader's AIMPAC Program          !
!   AIMPACrc: Write information into a WFN file (7)                    !
!   FCHKrc: Create a formatted checkpoint (FCHK) input file (19)       !
!   MOLDENrc: Create an input file (MLD) in Molden Format (17)         !

!   SQUARETRIAN2: Put square matrix (HF density matrix) into diagonal  !
!              form without avoiding double counting of diagonal terms !                 
!   SQUARETRIAN3: Put square matrix (PNOF density matrix) to diagonal  !
!              form without avoiding double counting of diagonal terms !                                                    
!   TRACEs: Calculate the trace of lagrangian and one-electron         !
!           derivative integrals                                       !
!                                                                      !
!======================================================================!

! AIMMEMrc
      SUBROUTINE AIMMEMrc(COEF,ZNUC,IZCORE,CX0,CY0,CZ0,KSTART,KNG,KKMIN,&
                        KKMAX,KATOM,KTYPE,RO,E,EX1,CS,CP,CD,CF,CG,CH,CI)
!-----------------------------------------------------------------------
!     Create an input File for Bader's AIMPAC Program
!-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      INTEGER,DIMENSION(NATOMS)::IZCORE
      INTEGER,DIMENSION(NSHELL)::KSTART,KNG,KKMIN,KKMAX,KATOM,KTYPE
      DOUBLE PRECISION,DIMENSION(NBF)::E
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::COEF
      DOUBLE PRECISION,DIMENSION(NBF5)::RO
      DOUBLE PRECISION,DIMENSION(NATOMS)::ZNUC,CX0,CY0,CZ0
      DOUBLE PRECISION,DIMENSION(NPRIMI)::EX1,CS,CP,CD,CF,CG,CH,CI
!-----------------------------------------------------------------------
      NPRIMS = 0
      DO I = 1,NSHELL
       NPRIMS = NPRIMS + KNG(I) * (KKMAX(I)-KKMIN(I)+1)
      ENDDO
!-----------------------------------------------------------------------
      CALL AIMPACrc(NBF5,NPRIMS,COEF,ZNUC,IZCORE,CX0,CY0,CZ0,KSTART,    &
           KNG,KKMIN,KKMAX,KATOM,KTYPE,RO,E,EX1,CS,CP,CD,CF,CG,CH,CI)
!-----------------------------------------------------------------------
      RETURN
      END

! AIMPACrc
      SUBROUTINE AIMPACrc(NMO,NPRIMS,VEC,ZNUC,IZCORE,CX0,CY0,CZ0,       &
                          KSTART,KNG,KKMIN,KKMAX,KATOM,KTYPE,RO,E,      &
                          EX1,CS,CP,CD,CF,CG,CH,CI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      CHARACTER(80) :: TITLE
      COMMON/TIT/TITLE
      COMMON/EHFEN/EHF,EN
      COMMON/ENERGIAS/EELEC,EELEC_OLD,DIF_EELEC,EELEC_MIN
      COMMON/VirialRatio/Virial                             ! RHF
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
!
      CHARACTER*4,DIMENSION(106)::ATMLAB
      INTEGER,DIMENSION(NATOMS)::IZCORE
      INTEGER,DIMENSION(NSHELL)::KSTART,KNG,KKMIN,KKMAX,KATOM,KTYPE
      INTEGER,DIMENSION(NPRIMS)::ICENT,ITYPE
      DOUBLE PRECISION,DIMENSION(NBF)::E
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::VEC
      DOUBLE PRECISION,DIMENSION(NBF5)::RO
      DOUBLE PRECISION,DIMENSION(NATOMS)::ZNUC,CX0,CY0,CZ0
      DOUBLE PRECISION,DIMENSION(NPRIMI)::EX1,CS,CP,CD,CF,CG,CH,CI
      DOUBLE PRECISION,DIMENSION(NPRIMS)::EXPON,AIC,ACO
!-----------------------------------------------------------------------
      DATA ATMLAB/'  H ','  HE','  LI','  BE','  B ','  C ',            &
                  '  N ','  O ','  F ','  NE','  NA','  MG',            &
                  '  AL','  SI','  P ','  S ','  CL','  AR',            &
                  '  K ','  CA','  SC','  TI','  V ','  CR',            &
                  '  MN','  FE','  CO','  NI','  CU','  ZN',            &
                  '  GA','  GE','  AS','  SE','  BR','  KR',            &
                  '  RB','  SR','  Y ','  ZR','  NB','  MO',            &
                  '  TC','  RU','  RH','  PD','  AG','  CD',            &
                  '  IN','  SN','  SB','  TE','  I ','  XE',            &
                  '  CS','  BA','  LA','  CE','  PR','  ND',            &
                  '  PM','  SM','  EU','  GD','  TB','  DY',            &
                  '  HO','  ER','  TM','  YB','  LU','  HF',            &
                  '  TA','  W ','  RE','  OS','  IR','  PT',            &
                  '  AU','  HG','  TL','  PB','  BI','  PO',            &
                  '  AT','  RN','  FR','  RA','  AC','  TH',            &
                  '  PA','  U ','  NP','  PU','  AM','  CM',            &
                  '  BK','  CF','  ES','  FM','  MD','  NO',            &
                  '  LR','    ','    ','    '/
!-----------------------------------------------------------------------
      REWIND(7)
      WRITE(7,1)TITLE
      WRITE(7,2)NMO,NPRIMS,NATOMS
!-----------------------------------------------------------------------
!     CX0,CY0,CZ0: Cartesian nuclear coordinate arrays
!     ZNUC: Nuclear charge array (1,NATOMS)
!     IZCORE: Number of electrons removed from each atom (ECP)
!-----------------------------------------------------------------------
      DO I=1,NATOMS
       IZNUC = INT(ZNUC(I))+IZCORE(I)
       WRITE(7,3)ATMLAB(IZNUC),I,I,CX0(I),CY0(I),CZ0(I),                &
                 ZNUC(I)+IZCORE(I)
      ENDDO
!-----------------------------------------------------------------------
!     Fill AIC (CONTRACTION COEFFS)
!-----------------------------------------------------------------------
      L = 0
      DO I = 1,NSHELL
       JST = KSTART(I)
       JEN = JST + KNG(I) - 1
       MINK = KKMIN(I)
       MAXK = KKMAX(I)
       DO K = MINK,MAXK
        DO J = JST,JEN
         L = L + 1
         ICENT(L) = KATOM(I)
         ITYPE(L) = K
         EXPON(L) = EX1(J)
         IF( K==1          )AIC(L) = CS(J)
         IF( 2<=K.and.K<= 4)AIC(L) = CP(J)
         IF( 5<=K.and.K<=10)AIC(L) = CD(J)
         IF(11<=K.and.K<=20)AIC(L) = CF(J)
         IF(21<=K.and.K<=35)AIC(L) = CG(J)
         IF(36<=K.and.K<=56)AIC(L) = CH(J)
         IF(57<=K.and.K<=84)AIC(L) = CI(J)
        ENDDO
       ENDDO
      ENDDO
!-----------------------------------------------------------------------
!     Write primitive information: ICENT, ITYPE, EXPON
!-----------------------------------------------------------------------
      WRITE(7,4) (ICENT(I),I=1,NPRIMS)
      WRITE(7,5) (ITYPE(I),I=1,NPRIMS)
      WRITE(7,6) (EXPON(I),I=1,NPRIMS)
!-----------------------------------------------------------------------
!     Expand the Natural Orbitals in the Atomic Orbital Basis
!     into the primitive representation. Write NOs and their OCC and E.
!-----------------------------------------------------------------------
      SQRT3 = SQRT(3.0D0)
      SQRT5 = SQRT(5.0D0)
      SQRT7 = SQRT(7.0D0)
      DO II=1,NMO
       L = 0
       M = 0
       DO I=1,NSHELL
        KST = KSTART(I)
        KEN = KST + KNG(I) - 1
        JMIN= KKMIN(I)
        JMAX= KKMAX(I)
        ITYP= KTYPE(I)
        DO J=JMIN,JMAX
         M = M+1
         DO K=KST,KEN
          L = L+1
          ACO(L) = VEC(M,II)*AIC(L)
          IF(ITYP==3)THEN
           IF(J>=8)ACO(L) = ACO(L)*SQRT3
          END IF
          IF(ITYP==4)THEN
           IF(J>=14)ACO(L)=ACO(L)*SQRT5
           IF(J==17)THEN
            DUMMY    = ACO(L-1)
            ACO(L-1) = ACO(L)
            ACO(L)   = DUMMY
           ENDIF
           IF(J==20)ACO(L) = ACO(L)*SQRT3
          ENDIF
          IF(ITYP==5)THEN
           IF(J>=24)ACO(L) = ACO(L)*SQRT7
           IF(J>=30)ACO(L) = ACO(L)*SQRT5/SQRT3
           IF(J>=33)ACO(L) = ACO(L)*SQRT3
          ENDIF
         ENDDO
        ENDDO
       ENDDO
       WRITE(7,7) II,2.0*RO(II),E(II)
       WRITE(7,8) (ACO(L),L=1,NPRIMS)
      ENDDO
!-----------------------------------------------------------------------
!     Terminate data, add SCF Energy and Virial
!-----------------------------------------------------------------------
      WRITE(7,9)
      WRITE(7,10)'THE RHF ',EELEC+EN,Virial
!-----------------------------------------------------------------------
    1 FORMAT(A80)
    2 FORMAT ('GAUSSIAN',10X,I5,' MOL ORBITALS',1X,I6,' PRIMITIVES',    &
              4X,I5,' NUCLEI')                                           
    3 FORMAT(A4,I4,4X,'(CENTRE',I3,')',1X,3F12.8,'  CHARGE =',F5.1)      
    4 FORMAT('CENTRE ASSIGNMENTS',2X,20I3)                               
    5 FORMAT('TYPE ASSIGNMENTS',4X,20I3)                                 
    6 FORMAT('EXPONENTS',1X,1P,5E14.7)                                   
    7 Format('MO',I5,5X,'MO 0.0',8X,'OCC NO = ',F12.7,                  &
             '  ORB. ENERGY =', F12.6)
    8 FORMAT(1P,5E16.8)
    9 FORMAT('END DATA')
   10 FORMAT(A8,' ENERGY =',F20.12,' THE VIRIAL(-V/T)=',F13.8)
!-----------------------------------------------------------------------
      RETURN
      END

! FCHKrc      
      SUBROUTINE FCHKrc(COEF,ZNUC,IZCORE,CX0,CY0,CZ0,KNG,KATOM,KTYPE,   &
                        RO,E,EX1,C1,C2,IMIN,IMAX,ISH,DIPS,IRUNTYP,IGTYP)
!-----------------------------------------------------------------------
!     Create a Formatted Checkpoint input File for Gaussview
!-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      CHARACTER(80) :: TITLE
      COMMON/TIT/TITLE
      CHARACTER(80) :: BASIS_FILE
      COMMON/BASIS_FILE/BASIS_FILE
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_PNOF/IPNOF,NTWOPAR
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      COMMON/CONV/ACURCY,EN,Etot,EHF,EHF0,DIFF,ITER,ICALP,ICBET
      COMMON/ENERGIAS/EELEC,EELEC_OLD,DIF_EELEC,EELEC_MIN
      COMMON/VirialRatio/Virial                             ! RHF
      INTEGER,DIMENSION(NATOMS)::IZCORE,IMIN,IMAX
      INTEGER,DIMENSION(NPRIMI)::ISH
      INTEGER,DIMENSION(NSHELL)::KNG,KATOM,KTYPE
      DOUBLE PRECISION,DIMENSION(NBF)::E
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::COEF
      DOUBLE PRECISION,DIMENSION(NBF5)::RO
      DOUBLE PRECISION,DIMENSION(NATOMS)::ZNUC,CX0,CY0,CZ0
      DOUBLE PRECISION,DIMENSION(NPRIMI)::EX1,C1,C2
      DOUBLE PRECISION,DIMENSION(3):: DIPS

      INTEGER :: IGTYP
      INTEGER,ALLOCATABLE,DIMENSION(:) :: NUMPRIMi,IZNUC,ISHELLTYP
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: DMs
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: DM
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: TMP !Fix: JFHLewYee
      CHARACTER(6),DIMENSION(4) :: RUNTYP
      DATA RUNTYP/'ENERGY','GRAD  ','OPTGEO','HESS  '/
      INTEGER,DIMENSION(10) :: FORD
      DATA FORD / 1,2,3,6,4,5,8,9,7,10/                 
      INTEGER,DIMENSION(15) :: GORD
      DATA GORD / 3,9,12,7,2,8,15,14,6,11,13,10,5,4,1/
      INTEGER,DIMENSION(21) :: HORD
      DATA HORD / 3,9,15,13,7,2,8,18,21,17,6,14,20,19,12,11,16,10,5,4,1/                 
      INTEGER,DIMENSION(28) :: IORD
      DATA IORD / 3,9,15,21,13,7,2,8,18,27,25,17,6,14,26,28,24,12,20,   &
                  23,22,19,11,16,10,5,4,1/
!----------------------------------------------------------------------- 
      REWIND(19)
      WRITE(19,1)TITLE
!     Run Type, Functional, Basis
      WRITE(19,2)RUNTYP(IRUNTYP),IPNOF,BASIS_FILE
!     
      WRITE(19,3)NATOMS
      WRITE(19,4)ICH
      WRITE(19,5)MUL
      WRITE(19,6)NE
      WRITE(19,7)NA
      WRITE(19,8)NB
      WRITE(19,9)NBF
      WRITE(19,10)NBF5
      WRITE(19,11)NSHELL
!     Highest angular momentum 
      CALL BASCHK(MAXANG,KTYPE,NSHELL)
      WRITE(19,12)MAXANG  
!     Largest degree of contraction 
      ALLOCATE(NUMPRIMi(NSHELL))
      LDCONT = 1
      DO I=1,NATOMS
       ISH1 = ISH(IMIN(I))
       NSHELLi = 1
       NUMPRIMi(NSHELLi) = 0
!      Number of primitives per shell of Atom 'i'
       DO J=IMIN(I),IMAX(I)
        if(ISH(J)==ISH1)then
         NUMPRIMi(NSHELLi) = NUMPRIMi(NSHELLi) + 1
        else
         NSHELLi = NSHELLi + 1
         NUMPRIMi(NSHELLi) = 1
         ISH1 = ISH(J)
        endif
       END DO
       LDCONTi = MAXVAL(NUMPRIMi(1:NSHELLi))
       if( LDCONTi > LDCONT ) LDCONT = LDCONTi
      END DO
      DEALLOCATE(NUMPRIMi)
      WRITE(19,13)LDCONT
!
      WRITE(19,14)NPRIMI
      WRITE(19,15)Virial
      ALLOCATE(IZNUC(NATOMS))
      DO i=1,NATOMS
       IZNUC(i) = INT(ZNUC(i))+IZCORE(i)
      END DO
      WRITE(19,16)NATOMS
      WRITE(19,17)(IZNUC(i),i=1,NATOMS)
      DEALLOCATE(IZNUC)
      WRITE(19,18)NATOMS
      WRITE(19,19)(ZNUC(i)+IZCORE(i),i=1,NATOMS)
      WRITE(19,20)3*NATOMS
      WRITE(19,19)(CX0(i),CY0(i),CZ0(i),i=1,NATOMS)
      WRITE(19,21)EELEC + EN
      WRITE(19,22)NSHELL
      ALLOCATE(ISHELLTYP(NSHELL))
!     Shell Types (KTYPE-1) 
      LSHELL = 0
      DO I=1,NSHELL
       if(KTYPE(I)<8)then
        ISHELLTYP(I) = KTYPE(I) - 1
       else
        ISHELLTYP(I) = -1
        LSHELL = 1
       end if
      END DO
      DO i=1,NSHELL
        IF(IGTYP .EQ. 2 .AND. ISHELLTYP(i) > 1) THEN
          WRITE(19,17) -ISHELLTYP(i)
        ELSE
          WRITE(19,17) ISHELLTYP(i)
        END IF
      END DO
      WRITE(19,17)(ISHELLTYP(i),i=1,NSHELL)

      IF (IGTYP.EQ.1) THEN
        !Start - block fix JFHLewYee
        ALLOCATE(TMP(NBF,NBF))
        TMP = COEF
        LOC = 0
        DO I=1,NSHELL
         IF(ISHELLTYP(I).EQ.0) THEN
           COEF(LOC+1, :) = TMP(LOC+1, :) 
         ELSE IF(ISHELLTYP(I).EQ.1) THEN
           COEF(LOC+1, :) = TMP(LOC+1, :) 
           COEF(LOC+2, :) = TMP(LOC+2, :) 
           COEF(LOC+3, :) = TMP(LOC+3, :) 
         ELSE IF(ISHELLTYP(I).EQ.2) THEN
           COEF(LOC+1, :) = TMP(LOC+1, :) 
           COEF(LOC+2, :) = TMP(LOC+4, :) 
           COEF(LOC+3, :) = TMP(LOC+6, :) 
           COEF(LOC+4, :) = TMP(LOC+2, :)  / sqrt(3.0D0)
           COEF(LOC+5, :) = TMP(LOC+3, :)  / sqrt(3.0D0)
           COEF(LOC+6, :) = TMP(LOC+5, :)  / sqrt(3.0D0)
         ELSE IF(ISHELLTYP(I).EQ.3) THEN
           COEF(LOC+1, :) = TMP(LOC+1, :) 
           COEF(LOC+2, :) = TMP(LOC+7, :) 
           COEF(LOC+3, :) = TMP(LOC+10, :) 
           COEF(LOC+4, :) = TMP(LOC+4, :)  / sqrt(5.0D0)
           COEF(LOC+5, :) = TMP(LOC+2, :)  / sqrt(5.0D0)
           COEF(LOC+6, :) = TMP(LOC+3, :)  / sqrt(5.0D0)
           COEF(LOC+7, :) = TMP(LOC+6, :)  / sqrt(5.0D0)
           COEF(LOC+8, :) = TMP(LOC+9, :)  / sqrt(5.0D0)
           COEF(LOC+9, :) = TMP(LOC+8, :)  / sqrt(5.0D0)
           COEF(LOC+10,:) = TMP(LOC+5, :)  / sqrt(15.0D0)
         ELSE IF(ISHELLTYP(I).EQ.4) THEN
           COEF(LOC+1, :) = TMP(LOC+1, :)
           COEF(LOC+2, :) = TMP(LOC+11, :)
           COEF(LOC+3, :) = TMP(LOC+15, :)
           COEF(LOC+4, :) = TMP(LOC+2, :)  / sqrt(7.0D0)
           COEF(LOC+5, :) = TMP(LOC+3, :)  / sqrt(7.0D0)
           COEF(LOC+6, :) = TMP(LOC+7, :)  / sqrt(7.0D0)
           COEF(LOC+7, :) = TMP(LOC+12, :) / sqrt(7.0D0)
           COEF(LOC+8, :) = TMP(LOC+10, :) / sqrt(7.0D0)
           COEF(LOC+9, :) = TMP(LOC+14, :) / sqrt(7.0D0)
           COEF(LOC+10,:) = TMP(LOC+4, :)  * sqrt(3.0D0) / sqrt(35.0D0)
           COEF(LOC+11,:) = TMP(LOC+6, :)  * sqrt(3.0D0) / sqrt(35.0D0)
           COEF(LOC+12,:) = TMP(LOC+13, :) * sqrt(3.0D0) / sqrt(35.0D0)
           COEF(LOC+13,:) = TMP(LOC+5, :)  / sqrt(35.0D0)
           COEF(LOC+14,:) = TMP(LOC+8, :)  / sqrt(35.0D0)
           COEF(LOC+15,:) = TMP(LOC+9, :)  / sqrt(35.0D0)
         ELSE IF(ISHELLTYP(I).GT.4) THEN
           !TODO
         END IF
         LOC = LOC + (ISHELLTYP(i)+1)*(ISHELLTYP(i)+2)/2 
        END DO
        DEALLOCATE(TMP)
        !End - block fix JFHLewYee
      END IF

      DEALLOCATE(ISHELLTYP)
!
      WRITE(19,23)NSHELL
      WRITE(19,17)(KNG(i),i=1,NSHELL)
      WRITE(19,24)NSHELL
      WRITE(19,17)(KATOM(i),i=1,NSHELL)
      WRITE(19,25)3*NSHELL
      WRITE(19,19)(CX0(KATOM(i)),CY0(KATOM(i)),CZ0(KATOM(i)),i=1,NSHELL)
      WRITE(19,26)NPRIMI
      WRITE(19,19)(EX1(i),i=1,NPRIMI)
      WRITE(19,27)NPRIMI
      WRITE(19,19)(C1(i),i=1,NPRIMI)
      if(LSHELL==1)then
       WRITE(19,28)NPRIMI
       WRITE(19,19)(C2(i),i=1,NPRIMI)
      end if
      
      WRITE(19,29)NBF5
      WRITE(19,19)(E(i),i=1,NBF5)
      WRITE(19,30)NBF*NBF5
      WRITE(19,19)((COEF(i,j),i=1,NBF),j=1,NBF5)

!     Total SCF Density  
      WRITE(19,31)NBFT
      ALLOCATE (DM(NBF,NBF),DMs(NBFT))
      CALL DENMATr(DM,COEF,RO,NBF,1,NBF5)
      NZ=0
      DO I=1,NBF
       DO J=1,I
        NZ=NZ+1
        DMs(NZ)=DM(I,J)
       ENDDO
      ENDDO
      WRITE(19,19)(DMs(i),i=1,NBFT)
      DEALLOCATE (DM,DMs)
      WRITE(19,32)NBF5
      WRITE(19,19)(RO(i),i=1,NBF5)
      WRITE(19,33)
      WRITE(19,19)(DIPS(i),i=1,3)
!-----------------------------------------------------------------------
    1 FORMAT(A80)
    2 FORMAT(A6,4X,'PNOF',I1,4X,A80)
    3 FORMAT('Number of atoms                ',12X,'I',5X,I12)
    4 FORMAT('Charge                         ',12X,'I',5X,I12)
    5 FORMAT('Multiplicity                   ',12X,'I',5X,I12)
    6 FORMAT('Number of electrons            ',12X,'I',5X,I12)
    7 FORMAT('Number of alpha electrons      ',12X,'I',5X,I12)
    8 FORMAT('Number of beta electrons       ',12X,'I',5X,I12)
    9 FORMAT('Number of basis functions      ',12X,'I',5X,I12)
   10 FORMAT('Number of independent functions',12X,'I',5X,I12)
   11 FORMAT('Number of contracted shells    ',12X,'I',5X,I12)
   12 FORMAT('Highest angular momentum       ',12X,'I',5X,I12)
   13 FORMAT('Largest degree of contraction  ',12X,'I',5X,I12)
   14 FORMAT('Number of primitive shells     ',12X,'I',5X,I12)
   15 FORMAT('Virial Ratio                   ',12X,'R',5X,ES22.15)
   16 FORMAT('Atomic numbers                 ',12X,'I   N=',I12)
   17 FORMAT(6I12) 
   18 FORMAT('Nuclear charges                ',12X,'R   N=',I12)
   19 FORMAT(5ES16.8)   
   20 FORMAT('Current cartesian coordinates  ',12X,'R   N=',I12)
   21 FORMAT('Total Energy                   ',12X,'R',5X,ES22.15)
   22 FORMAT('Shell types                    ',12X,'I   N=',I12)
   23 FORMAT('Number of primitives per shell ',12X,'I   N=',I12)
   24 FORMAT('Shell to atom map              ',12X,'I   N=',I12)
   25 FORMAT('Coordinates of each shell      ',12X,'R   N=',I12)
   26 FORMAT('Primitive exponents            ',12X,'R   N=',I12)            
   27 FORMAT('Contraction coefficients       ',12X,'R   N=',I12)
   28 FORMAT('P(S=P) Contraction coefficients',12X,'R   N=',I12)
   29 FORMAT('Alpha Orbital Energies         ',12X,'R   N=',I12)
   30 FORMAT('Alpha MO coefficients          ',12X,'R   N=',I12)
   31 FORMAT('Total SCF Density              ',12X,'R   N=',I12)
   32 FORMAT('Natural orbital occupancies    ',12X,'R   N=',I12)
   33 FORMAT('Dipole Moment                  ',12X,'R   N=',11X,'3')
!-----------------------------------------------------------------------              
      RETURN
      END
      
! MOLDENrc
      SUBROUTINE MOLDENrc(ATMNAME,IZCORE,ZNUC,CX0,CY0,CZ0,IMIN,IMAX,ISH,&
                          ITYP,EX1,C1,C2,RO,E,COEF)
!-----------------------------------------------------------------------
!     Create an input File to interface to the Molden program
!-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      CHARACTER(80) :: TITLE
      COMMON/TIT/TITLE
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      CHARACTER*4,DIMENSION(NATOMS)::ATMNAME
      INTEGER,DIMENSION(NATOMS)::IZCORE,IMIN,IMAX
      INTEGER,DIMENSION(NPRIMI)::ISH,ITYP
      DOUBLE PRECISION,DIMENSION(NATOMS)::ZNUC,CX0,CY0,CZ0
      DOUBLE PRECISION,DIMENSION(NPRIMI)::EX1,C1,C2
      DOUBLE PRECISION,DIMENSION(NBF5)::RO
      DOUBLE PRECISION,DIMENSION(NBF)::E
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::COEF
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
      CHARACTER*2 LABEL(8)
      DATA LABEL/'s ','p ','d ','f ','g ','h ','i ','sp'/
      INTEGER,ALLOCATABLE,DIMENSION(:)::NUMPRIM
!-----------------------------------------------------------------------
      REWIND(17)
      WRITE(17,1)
      WRITE(17,2)TITLE
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Coordinates for the Electron Density / Molecular orbitals
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      WRITE(17,3)      
      DO I=1,NATOMS
       IZNUC = INT(ZNUC(I))+IZCORE(I)
       ATMNAME(I)=ATMLAB(IZNUC)
      ENDDO
      DO I=1,NATOMS
       IZNUC = INT(ZNUC(I))+IZCORE(I)
       WRITE(17,4)ATMNAME(I),I,IZNUC,CX0(I),CY0(I),CZ0(I)
      ENDDO      
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Specification of the basis-set consisting of Contracted GTOs
!     Recognized shell_labels by Molden: s,p,d,f,g,sp
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      WRITE(17,5)
      ALLOCATE(NUMPRIM(NSHELL))
      DO I=1,NATOMS
       WRITE(17,6)I,0
!      Number of primitives per shell of Atom 'i'
       ISH1 = ISH(IMIN(I))
       NSHELLi = 1
       NUMPRIM(NSHELLi) = 0
       DO J=IMIN(I),IMAX(I)
        if(ISH(J)==ISH1)then
         NUMPRIM(NSHELLi) = NUMPRIM(NSHELLi) + 1
        else
         NSHELLi = NSHELLi + 1
         NUMPRIM(NSHELLi) = 1
         ISH1 = ISH(J)
        endif
       ENDDO
!      shell_label, number_of_primitives
       I1 = IMIN(I)
       DO K=1,NSHELLi
        WRITE(17,7)LABEL(ITYP(I1)),NUMPRIM(K),1.00
!       exponent_primitive, contraction_coefficient
        do j=I1,I1+NUMPRIM(K)-1
         if(ITYP(J)<8)then
          WRITE(17,8)EX1(j),C1(j)
         else
          WRITE(17,8)EX1(j),C1(j),C2(j)
         endif
        enddo        
        I1 = I1 + NUMPRIM(K) 
       ENDDO
       write(17,*)  ! empty line
      ENDDO
      DEALLOCATE(NUMPRIM)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Specification of the molecular orbitals    
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      WRITE(17,9)
      DO I=1,NBF5
       WRITE(17,10)I
       WRITE(17,11)E(I)
       WRITE(17,12)
       WRITE(17,13)2.0d0*RO(I)
       DO J=1,NBF
        WRITE(17,14)J,COEF(J,I)
       ENDDO
      ENDDO
      DO I=NBF5+1,NBF
       WRITE(17,10)I
       WRITE(17,11)E(I)
       WRITE(17,12)
       WRITE(17,13)0.0d0
       DO J=1,NBF
        WRITE(17,14)J,COEF(J,I)
       ENDDO
      ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
!-----------------------------------------------------------------------
    1 FORMAT('[Molden Format]')
    2 FORMAT(A80)
    3 FORMAT('[Atoms] AU')
    4 FORMAT(1X,A4,2I6,3F14.4)                              
    5 FORMAT('[GTO]')
    6 FORMAT(I5,I2)
    7 FORMAT(A2,I4,F5.2)
    8 FORMAT(1X,F20.12,2E20.12)
    9 FORMAT('[MO]')
   10 FORMAT('Sym=',I6,'a')
   11 FORMAT('Ene=',E16.6)
   12 FORMAT('Spin= Alpha')
   13 FORMAT('Occup=',F9.6)
   14 FORMAT(I6,1X,F12.6)
!-----------------------------------------------------------------------
      END
