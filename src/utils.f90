!======================================================================!
!                                                                      !
!         M A T H   U T I L I T I E S    S U B R O U T I N E S         !
!                                                                      !
!======================================================================!
!                                                                      !
!   ORTHONORMAL: Orthogonalize and Orthonormalize MOs (COEF)           !
!   CANORTORB: Form Canonical Orthonormal Orbitals that diagonalize S  !
!   ORTHOGONALIZE: Project COEF into W vectors                         !
!   ORTHONORMALIZE: Orthonormalize COEF                                !
!   CHECKORTHO: Check the Orthonormality of the MOs                    !
!   PRODWCWFik: PRODWCWFik = W(ik)*CW12(q)                             !
!   PRODWCWFikq: PRODWCWFikq = W(ik)*CW12(q) excluding q'=q            !
!   SORT_IDX_DECREASING: Indexes for decreasing sorting an array       !
!                                                                      !
!======================================================================!


! ORTHONORMAL
      SUBROUTINE ORTHONORMAL(NV,NBF,NBFT,OVERLAP,COEF,ICHECK,IPRINTOPT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      INTEGER::NV,NBF,NBFT,ICHECK,IPRINTOPT      
      DOUBLE PRECISION,DIMENSION(NBF,NV)::COEF
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::OVERLAP
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::EVA,AUXE,S
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::AUX,W,COEFN
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Form Canonical Orthonormal Orbitals: Diagonalize the Overlap (AUX)
!     W - eigenvectors, EVA - eigenvalues
!     LIV - number of linearly independent eigenvectors 
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(AUX(NBF,NBF),W(NBF,NBF),EVA(NBF),AUXE(NBF))
      AUX = OVERLAP
      CALL CANORTORB(AUX,W,EVA,AUXE,LIV,NBF)       ! LIV <= NBF
      DEALLOCATE(AUX)
!     Square OVERLAP -> Symmetric S
      ALLOCATE(S(NBFT))
      NZ=0
      DO I=1,NBF
       DO J=1,I
        NZ=NZ+1
        S(NZ)=OVERLAP(I,J)
       ENDDO
      ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Orthogonalize orbitals
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Project COEF into W vectors
      CALL ORTHOGONALIZE(W,S,COEF,AUXE,NBF,LIV,NBFT)
!     Orthonormalize COEF
      CALL ORTHONORMALIZE(COEF,NBF,LIV)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     New Coefficients: Back-Transform to first basis (COEFN=COEF*W)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(COEFN(NBF,NV))
      CALL COEFW(NBF,LIV,COEFN,W,COEF)
      COEF=COEFN
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Checking the orthonormality of the new orbitals, if they are not
!     orthonormal then the coefficients are orthonormalized once again
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(ICHECK==1)THEN
       CALL CHECKORTHO(COEF,OVERLAP,IVIOORTHO,IPRINTOPT)
       IF(IVIOORTHO/=0)THEN
        CALL ORTHOGONALIZE(W,S,COEF,AUXE,NBF,LIV,NBFT)
        CALL ORTHONORMALIZE(COEF,NBF,LIV)
        CALL COEFW(NBF,LIV,COEFN,W,COEF)
        COEF=COEFN
        CALL CHECKORTHO(COEF,OVERLAP,IVIOORTHO,IPRINTOPT)
       ENDIF
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DEALLOCATE(COEFN,W,EVA,AUXE,S)
      RETURN
      END

! CANORTORB
      SUBROUTINE CANORTORB(S,Q,EIG,AUX,LIV,N)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION S(N,N),Q(N,N),EIG(N),AUX(N)
!-----------------------------------------------------------------------
      CALL DIAG(N,S,Q,EIG,AUX)
      DUM = EIG(1)
      J = 0
      K = 0
      DO I = 1,N
       IF(EIG(I)<1.0D-16)THEN
        K = K+1
       ELSE
        J = J+1
        EIG(J) = 1.0D0/SQRT(EIG(I))
        Q(1,J)=Q(1,I)
       ENDIF
      ENDDO
      LIV = N-K                                ! LIV <= N

      DO J=1,LIV
       DO i=1,N
        QiJ=Q(i,J)*EIG(J)
        IF(ABS(QiJ)<1.0D-16)QiJ=0.0d0
        Q(i,J)=QiJ
       ENDDO
      ENDDO

!      IF(LIV/=N)THEN
!       WRITE(6,'(/,50A,I5,/)')
!     & 'The number of linearly independent eigenvectors is',LIV
!       STOP
!      ENDIF
!-----------------------------------------------------------------------
      RETURN
      END

! ORTHOGONALIZE
      SUBROUTINE ORTHOGONALIZE(Q,S,V,T,N,LIV,NT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      DIMENSION Q(N,LIV),S(NT),V(N,N),T(N)
!-----------------------------------------------------------------------
      DO J=1,N
       IK=0
       DO I=1,N
        DUM = 0.0d0
        VDUM = V(I,J)
        IF(I>1) THEN
         DO K=1,I-1
          IK = IK+1
          DUM = DUM+S(IK)*V(K,J)
          T(K) = T(K)+S(IK)*VDUM
         ENDDO
        ENDIF
        IK = IK+1
        T(I) = DUM+S(IK)*VDUM
       ENDDO
       DO I=1,N
        V(I,J) = T(I)
       ENDDO
      ENDDO
!-----------------------------------------------------------------------
      DO J=1,N
       DO I=1,LIV
        T(I) = DDOT(N,Q(1,I),1,V(1,J),1)
       ENDDO
       DO I=1,LIV
        V(I,J) = T(I)
       ENDDO
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! ORTHONORMALIZE
      SUBROUTINE ORTHONORMALIZE(V,N,LIV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DIMENSION V(N,N)
!-----------------------------------------------------------------------
!     Orthonormalize first LIV Orbitals
!-----------------------------------------------------------------------
      IF(LIV>0)THEN
       DO I=1,LIV
        DUMI = 0.0d0
        DO K=1,N
         DUMI = DUMI + V(K,I)*V(K,I)
        ENDDO
        DUMI = 1.0d0/ SQRT(DUMI)
        DO K=1,N
         V(K,I) = V(K,I)*DUMI
        ENDDO
        IF(I<LIV)THEN
         DO J=I+1,LIV
          DUM = -DDOT(N,V(1,J),1,V(1,I),1)
          CALL DAXPY(N,DUM,V(1,I),1,V(1,J),1)
         ENDDO
        ENDIF
       ENDDO
      ENDIF
      IF(LIV==N)RETURN
!-----------------------------------------------------------------------
!     Get Orthogonal Space
!-----------------------------------------------------------------------
      I = LIV
      J = 0
    1 I0 = I
      I = I+1
      IF(I>N)RETURN
    2 J = J+1
      IF(J>N)THEN
       WRITE(6,10)I0,N
       STOP
      ENDIF
      DO K=1,N
       V(K,I) = 0.0d0
      ENDDO
      V(J,I) = 1.0d0
      DO L=1,I0
       DUM = -DDOT(N,V(1,L),1,V(1,I),1)
       CALL DAXPY(N,DUM,V(1,L),1,V(1,I),1)
      ENDDO
      DUMI = 0.0d0
      DO K=1,N
       DUMI = DUMI+V(K,I)*V(K,I)
      ENDDO
      IF(ABS(DUMI)<1.0D-16)GOTO 2
      DUMI = 1.0d0/ SQRT(DUMI)
      DO K=1,N
       V(K,I) = V(K,I)*DUMI
      ENDDO
      GOTO 1
   10 FORMAT(1X,'Orthogonalizing found a redundant set of vectors'/     &
             1X,'only',I5,' independent vectors found but there'/       &
             1X,'should be',I5,' MOs')
!-----------------------------------------------------------------------
      END
      
! CHECKORTHO
      SUBROUTINE CHECKORTHO(COEF,OVERLAP,IVIOORTHO,IPRINTOPT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      INTEGER :: IVIOORTHO,IPRINTOPT      
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_NTHRESH/NTHRESHL,NTHRESHE,NTHRESHEC,NTHRESHEN
      DIMENSION COEF(NBF,NBF),OVERLAP(NBF,NBF)
!-----------------------------------------------------------------------
      IF(NTHRESHE<=8)THEN
       THRESHORTHO=1.0d-08
      ELSE
       THRESHORTHO=10.0**(-NTHRESHE)
      ENDIF
      IVIOORTHO=0
      CSCMAX=0.0d0
      !$OMP PARALLEL DO PRIVATE(IQ,JQ,CSC,i) REDUCTION(+:IVIOORTHO) REDUCTION(MAX: CSCMAX)
      DO IQ=1,NBF
       DO JQ=1,IQ
        CSC=0.0d0
        do i=1,nbf
         CSC = CSC + COEF(i,IQ)*FC(i,JQ,OVERLAP,COEF)
        enddo
        IF(IQ==JQ.and.ABS(1.0-CSC)>THRESHORTHO)THEN
         IVIOORTHO=IVIOORTHO+1
         if(CSCMAX<ABS(1.0-CSC))CSCMAX=ABS(1.0-CSC)
!         write(6,*)IQ,JQ,CSC         
        ELSEIF(IQ/=JQ.and.ABS(CSC)>THRESHORTHO)THEN
         IVIOORTHO=IVIOORTHO+1
         if(CSCMAX<ABS(CSC))CSCMAX=ABS(CSC)
!         write(6,*)IQ,JQ,CSC         
        ENDIF
       ENDDO
      ENDDO
      !$OMP END PARALLEL DO
      IF(IPRINTOPT==0)RETURN
      IF(IVIOORTHO==0)THEN
       WRITE(6,1)
      ELSE
       WRITE(6,2)IVIOORTHO,CSCMAX
      ENDIF
      RETURN
!-----------------------------------------------------------------------
!     Format definitions
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    1 FORMAT(/,' No violations of the orthonormality')
    2 FORMAT(/,' Orthonormality violations:',I8,                        &
               ', Maximum Violation:',E8.1)
!-----------------------------------------------------------------------
      END

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
!                                                                      !
!                  U T I L       S U B R O U T I N E S                 !
!                                                                      !
!   ABRT: Generic Stop Subroutine                                      !
!   XtoX0: It keeps the linear arrangement X(N) in X0(N)               !
!   TRACE: Calculate the trace of a symmetric matrix A(N,N)            !
!   TRACEm: Calculate the trace (E) of a symmetric matrix A(N,N)       !
!   SQUARETRIAN: Put square matrix FM(N,N) in triangular form F(NT)    !
!   TRIANSQUARE: Put triangular F(NT) in square matrix form FM(N,N)    !
!   CeqAtB: C = At*B where C(NA,NB), A(NA,MA) and B((NA,NB))           !
!   COEFW: COEFN = COEF*W  where COEFN(N,L), COEF(N,N) and W(N,L)      !
!   COEFW1: COEFN = COEF*W where COEFN(L,N), COEF(L,N) and W(N,N)      !
!   FC: FC(i,IQ) = Sum_j[F(i,j)*C(j,IQ)], Sym. F(nbf,nbf), C(nbf,nbf)  !
!   DMATMAX: Determine the maximum off-diagonal element of symmetric A !
!   HEAV and DHEAV: Heaviside function and its derivative              !
!                                                                      !
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !

! ABRT
      SUBROUTINE ABRT  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      CHARACTER*24 STAMP                                                
      DOUBLE PRECISION TIMSTR(3)  
!-----------------------------------------------------------------------
      CALL FDATE(STAMP)                                                 
      READ(UNIT=STAMP,FMT='(3A8)')TIMSTR                               
      WRITE(6,1)TIMSTR
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    
!     GENERIC STOP, IN CASE YOU GET THIS FAR
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    
      STOP 'IN ABRT'                                                    
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    
    1 FORMAT(1X,'Execution terminated abnormally at ',3A8)                                           
      END                                                               

! XtoX0
      SUBROUTINE XtoX0(X,X0,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DOUBLE PRECISION,DIMENSION(N)::X,X0
!-----------------------------------------------------------------------
      DO I=1,N
       X0(I) = X(I)
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! IXtoIX0
      SUBROUTINE IXtoIX0(IX,IX0,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      INTEGER(8),DIMENSION(N)::IX,IX0
!-----------------------------------------------------------------------
      DO I=1,N
       IX0(I) = IX(I)
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END
      
! TRACE
      FUNCTION TRACE(D,A,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DOUBLE PRECISION,DIMENSION(N,N)::D,A
!-----------------------------------------------------------------------
!     Calculate the trace of a symmetric matrix A(N,N)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      TRACE = 0.0
      DO K=1,N
       DO L=1,K
        TRACE = TRACE + A(K,L)*D(K,L)
       ENDDO
       TRACE = TRACE - A(K,K)*D(K,K)*0.5
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! TRACEm
      SUBROUTINE TRACEm(E,D,AA,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DIMENSION D(N,N),AA(N,N)
!-----------------------------------------------------------------------
!     Calculate the trace (E) of a symmetric matrix A(N,N)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      E = 0.0d0
      DO K=1,N
       DO L=1,K
        E = E + AA(K,L)*D(K,L)
       ENDDO
       E = E - AA(K,K)*D(K,K)*0.5d0
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! TRACEP                                           
      DOUBLE PRECISION FUNCTION TRACEP(A,B,N)                           
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      DIMENSION A(*),B(*)                                               
      N2 = (N*N+N)/2                                                    
      TRACE = DDOT(N2,A,1,B,1)                                          
      TRACE = TRACE+TRACE                                               
      K = 0                                                             
      DO I = 1,N                                                    
       K = K+I                                                        
       TRACE = TRACE - A(K)*B(K)                                      
      END DO
      TRACEP = TRACE                                                    
      RETURN                                                            
      END                                                               
      
! SQUARETRIAN
      SUBROUTINE SQUARETRIAN(FM,F,N,NT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DOUBLE PRECISION,DIMENSION(N,N)::FM
      DOUBLE PRECISION,DIMENSION(NT)::F
!-----------------------------------------------------------------------
!     PUT SQUARE FM IN TRIANGULAR FORM F
!-----------------------------------------------------------------------
      NZ=0
      DO I=1,N
       DO J=1,I
        NZ=NZ+1
        F(NZ)=FM(I,J)+FM(I,J)
       ENDDO
       F(NZ)=F(NZ)-FM(I,I)
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! TRIANSQUARE
      SUBROUTINE TRIANSQUARE(FM,F,N,NT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DOUBLE PRECISION,DIMENSION(N,N)::FM
      DOUBLE PRECISION,DIMENSION(NT)::F
!-----------------------------------------------------------------------
!     PUT TRIANGULAR F IN SQUARE FORM FM
!-----------------------------------------------------------------------
      NZ=0
      DO I=1,N
       DO J=1,I
        NZ=NZ+1
        FM(I,J)=F(NZ)
        FM(J,I)=F(NZ)
       ENDDO
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! CeqAtB
      SUBROUTINE CeqAtB(C,A,NA,MA,B,NB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DOUBLE PRECISION,DIMENSION(NA,MA)::A
      DOUBLE PRECISION,DIMENSION(NA,NB)::B
      DOUBLE PRECISION,DIMENSION(NA,NB)::C
!-----------------------------------------------------------------------
!     C = At * B
!-----------------------------------------------------------------------
      DO i=1,MA
       DO j=1,NB
        C(i,j) = 0.0d0
        do k=1,NA
         C(i,j) = C(i,j) + A(k,i)*B(k,j)
        enddo
       ENDDO
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! COEFW
      SUBROUTINE COEFW(N,L,COEFN,COEF,W)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DOUBLE PRECISION,DIMENSION(N,L)::COEFN,W
      DOUBLE PRECISION,DIMENSION(N,N)::COEF
!-----------------------------------------------------------------------
!     COEFN = COEF*W
!-----------------------------------------------------------------------
      DO i=1,N
       DO j=1,L
        COEFN(i,j)=0.0d0
        do k=1,N
         COEFN(i,j)=COEFN(i,j)+COEF(i,k)*W(k,j)
        enddo
       ENDDO
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! COEFW1
      SUBROUTINE COEFW1(L,N,COEFN,COEF,W)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DOUBLE PRECISION,DIMENSION(L,N)::COEFN,COEF
      DOUBLE PRECISION,DIMENSION(N,N)::W
!-----------------------------------------------------------------------
!     COEFN = COEF*W
!-----------------------------------------------------------------------
      DO i=1,L
       DO j=1,N
        COEFN(i,j)=0.0d0
        do k=1,N
         COEFN(i,j) = COEFN(i,j) + COEF(i,k)*W(k,j)
        enddo
       ENDDO
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! FC
      FUNCTION FC(i,IQ,F,C)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::F,C
!-----------------------------------------------------------------------
!     FC(i,IQ) = Sum_j[F(i,j)*C(j,IQ)], F:Sym, C(nbf,nbf)
!-----------------------------------------------------------------------
      FC = 0.0d0
      do j=1,NBF
       FC = FC + F(i,j) * C(j,IQ)
      enddo
!-----------------------------------------------------------------------
      RETURN
      END

! DMATMAX
      SUBROUTINE DMATMAX(A,N,MAXI,MAXJ,DUM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DOUBLE PRECISION,DIMENSION(N,N)::A
!-----------------------------------------------------------------------
!     Determine the maximum off-diagonal element of symmetric matrix A 
!-----------------------------------------------------------------------
      DUM=0.0d0
      DO IQ=1,N
       DO JQ=1,IQ-1
        DMij=ABS(A(IQ,JQ))
        IF(DMij>DUM)THEN
         DUM=DMij
         MAXI=IQ
         MAXJ=JQ
        ENDIF
       ENDDO
      ENDDO
!      WRITE(6,1)MAXI,MAXJ,DUM
!-----------------------------------------------------------------------
!    1 FORMAT(4X,'(',I3,',',I3,')',4X,F15.6)
      RETURN
      END

! Heaviside
      FUNCTION HEAV(X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)         
      IF(X>=0.0d0)THEN
       HEAV = 1.0d0
      ELSE
       HEAV = 0.0d0
      ENDIF
      RETURN
      END

! Derivative of Heaviside
      FUNCTION DHEAV(X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)               
      IF(X>=0.0d0)THEN
       DHEAV = 0.0d0
      ELSE
       DHEAV = 0.0d0
      ENDIF
      RETURN
      END

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
!                                                                      !
!   DIAG: Diagonalization of real symmetric matrix using Houselholder  !
!         method. NAN is the maximum possible size = 5000.             !
!                                                                      !
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !

! DIAG
      SUBROUTINE DIAG(NN,A,VEC,EIG,W)                                  
!**********************************************************************!
!     MATRIX DIAGNOLIZATION ROUTINE FOR REAL SYMMETRIC CASE            !
!     HOUSEHOLDER METHOD                                               !
!     RHOSQ = UPPER LIMIT FOR OFF-DIAGONAL ELEMENT                     !
!     NN= SIZE OF MATRIX                                               !
!     A = MATRIX (ONLY LOWER TRIANGLE IS USED + THIS IS DESTROYED)     !
!     EIG = RETURNED EIGENVALUES IN ALGEBRAIC DESCENDING ORDER         !
!     VEC = RETURNED EIGENVECTORS IN COLUMNS                           !
!**********************************************************************!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      PARAMETER (NAN = 5000)
      INTEGER,DIMENSION(NAN)::IPOSV,IORD,IVPOS
      DOUBLE PRECISION,DIMENSION(NAN)::GAMMA,BETA,P,Q,BETASQ
      DOUBLE PRECISION,DIMENSION(NN)::EIG,W
      DOUBLE PRECISION,DIMENSION(NN,NN)::A,VEC
      EQUIVALENCE (P,Q,IVPOS,BETA),(IPOSV,GAMMA),(IORD,BETASQ)
      DATA ZERO/0.0d0/,PT5/0.5d0/,ONE/1.0d0/,TWO/2.0d0/,RHOSQ/1.0D-30/
!-----------------------------------------------------------------------
      N=NN
      IF (N == 0) GO TO 550
      N1=N-1
      N2=N-2
      GAMMA(1)=A(1,1)
      IF (N2<0)GO TO 180
      IF (N2==0)GO TO 170
      DO NR=1,N2
       B=A(NR+1,NR)
       S=ZERO
       DO I=NR,N2
        S=S+A(I+2,NR)**2
       END DO
!      PREPARE FOR POSSIBLE BYPASS OF TRANSFORMATION
       A(NR+1,NR)=ZERO
       IF (S <= ZERO) GO TO 150
       S=S+B*B
       SGN=+ONE
       IF (B < ZERO) SGN=-ONE
       SQRTS= SQRT(S)
       D=SGN/(SQRTS+SQRTS)
       TEMP= SQRT(PT5+B*D)
       W(NR)=TEMP
       A(NR+1,NR)=TEMP
       D=D/TEMP
       B=-SGN*SQRTS
!      D IS FACTOR OF PROPORTIONALITY. NOW COMPUTE AND SAVE W VECTOR.
!      EXTRA SINGLY SUBSCRIPTED W VECTOR USED FOR SPEED.
       DO I=NR,N2
        TEMP=D*A(I+2,NR)
        W(I+1)=TEMP
        A(I+2,NR)=TEMP
       END DO
!      PREMULTIPLY VECTOR W BY MATRIX A TO OBTAIN P VECTOR.
!      SIMULTANEOUSLY ACCUMULATE DOT PRODUCT WP,(THE SCALAR K)
       WTAW=ZERO
       DO I=NR,N1
        SUM=ZERO
        DO J=NR,I
         SUM=SUM+A(I+1,J+1)*W(J)
        END DO
        I1=I+1
        IF( (N1-I1) >= 0 ) THEN
         DO J=I1,N1
          SUM=SUM+A(J+1,I+1)*W(J)
         END DO
        END IF
        P(I)=SUM
        WTAW=WTAW+SUM*W(I)
       END DO
!      P VECTOR AND SCALAR K  NOW STORED. NEXT COMPUTE Q VECTOR
       DO I=NR,N1
        Q(I)=P(I)-WTAW*W(I)
       END DO
!      NOW FORM PAP MATRIX, REQUIRED PART
       DO J=NR,N1
        QJ=Q(J)
        WJ=W(J)
        DO I=J,N1
         A(I+1,J+1)=A(I+1,J+1)-TWO*(W(I)*QJ+WJ*Q(I))
        END DO
       END DO
  150  BETA(NR)=B
       BETASQ(NR)=B*B
       GAMMA(NR+1)=A(NR+1,NR+1)
      END DO
  170 B=A(N,N-1)
      BETA(N-1)=B
      BETASQ(N-1)=B*B
      GAMMA(N)=A(N,N)
  180 BETASQ(N)=ZERO
!     ADJOIN AN IDENTIFY MATRIX TO BE POSTMULTIPLIED BY ROTATIONS.
      DO I=1,N
       DO J=1,N
        VEC(I,J)=ZERO
       END DO
       VEC(I,I)=ONE
      END DO
      M=N
      SUM=ZERO
      NPAS=1
      GO TO 330
  210 SUM=SUM+SHIFT
      COSA=ONE
      G=GAMMA(1)-SHIFT
      PP=G
      PPBS=PP*PP+BETASQ(1)
      PPBR= SQRT(PPBS)
      DO J=1,M
       COSAP=COSA
       IF(PPBS == ZERO)THEN
        SINA=ZERO
        SINA2=ZERO
        COSA=ONE
        GO TO 270
       END IF
       SINA=BETA(J)/PPBR
       SINA2=BETASQ(J)/PPBS
       COSA=PP/PPBR
!      POSTMULTIPLY IDENTITY BY P-TRANSPOSE MATRIX
       NT=J+NPAS
       IF( (NT-N) >= 0 ) NT=N
       DO I=1,NT
        TEMP=COSA*VEC(I,J)+SINA*VEC(I,J+1)
        VEC(I,J+1)=-SINA*VEC(I,J)+COSA*VEC(I,J+1)
        VEC(I,J)=TEMP
       END DO
  270  DIA=GAMMA(J+1)-SHIFT
       U=SINA2*(G+DIA)
       GAMMA(J)=G+U
       G=DIA-U
       PP=DIA*COSA-SINA*COSAP*BETA(J)
       IF( (J-M) == 0 )THEN
        BETA(J)=SINA*PP
        BETASQ(J)=SINA2*PP*PP
        GO TO 310
       ENDIF
       PPBS=PP*PP+BETASQ(J+1)
       PPBR= SQRT(PPBS)
       BETA(J)=SINA*PPBR
       BETASQ(J)=SINA2*PPBS
      END DO
  310 GAMMA(M+1)=G
!     TEST FOR CONVERGENCE OF LAST DIAGONAL ELEMENT
      NPAS=NPAS+1
      IF( (BETASQ(M)-RHOSQ) > ZERO) GO TO 350
  320 EIG(M+1)=GAMMA(M+1)+SUM
  330 BETA(M)=ZERO
      BETASQ(M)=ZERO
      M=M-1
      IF(M == 0) GO TO 380
      IF( (BETASQ(M)-RHOSQ) <= ZERO )GO TO 320
!     TAKE ROOT OF CORNER 2 BY 2 NEAREST TO LOWER DIAGONAL IN VALUE
!     AS ESTIMATE OF EIGENVALUE TO USE FOR SHIFT
  350 A2=GAMMA(M+1)
      R2=PT5*A2
      R1=PT5*GAMMA(M)
      R12=R1+R2
      DIF=R1-R2
      TEMP= SQRT(DIF*DIF+BETASQ(M))
      R1=R12+TEMP
      R2=R12-TEMP
      DIF= ABS(A2-R1)- ABS(A2-R2)
      IF(DIF >= ZERO)THEN
       SHIFT=R2
       GO TO 210
      END IF
      SHIFT=R1
      GO TO 210
  380 EIG(1)=GAMMA(1)+SUM
!     INITIALIZE AUXILIARY TABLES REQUIRED FOR REARRANGING THE VECTORS
      DO J=1,N
       IPOSV(J)=J
       IVPOS(J)=J
       IORD(J)=J
      END DO
!     USE A TRANSPOSITION SORT TO ORDER THE EIGENVALUES
      M=N
      GO TO 430
  400 DO J=1,M
       IF( (EIG(J)-EIG(J+1)) > ZERO )THEN
        TEMP=EIG(J)
        EIG(J)=EIG(J+1)
        EIG(J+1)=TEMP
        ITEMP=IORD(J)
        IORD(J)=IORD(J+1)
        IORD(J+1)=ITEMP
       END IF
      END DO
  430 M=M-1
      IF(M /= 0)GO TO 400
      IF(N1 == 0)GO TO 490
      DO L=1,N1
       NV=IORD(L)
       NP=IPOSV(NV)
       IF( (NP-L) == 0 )GO TO 480
       LV=IVPOS(L)
       IVPOS(NP)=LV
       IPOSV(LV)=NP
       DO I=1,N
        TEMP=VEC(I,L)
        VEC(I,L)=VEC(I,NP)
        VEC(I,NP)=TEMP
       END DO
  480  CONTINUE
      END DO
!     BACK TRANSFORM THE VECTORS OF THE TRIPLE DIAGONAL MATRIX
  490 DO NRR=1,N
       K=N1
  500  K=K-1
       IF( K <= 0)GO TO 540
       SUM=ZERO
       DO I=K,N1
        SUM=SUM+VEC(I+1,NRR)*A(I+1,K)
       END DO
       SUM=SUM+SUM
       DO I=K,N1
        VEC(I+1,NRR)=VEC(I+1,NRR)-SUM*A(I+1,K)
       END DO
       GO TO 500
  540  CONTINUE
      END DO
  550 RETURN
      END

! PRODWCWFik
      FUNCTION PRODWCWFik(ik,INDf,W,CW12)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/INPFILE_FRAG/NO1f,NBFf,NBF5f
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
!
      INTEGER,DIMENSION(NBFf)::INDf
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5)::CW12
      DOUBLE PRECISION,DIMENSION(NSQ,NBF5)::W
!-----------------------------------------------------------------------
!     Sum W*CW12(1,*) por IQP
!-----------------------------------------------------------------------
      IQ = INDf(1)
      PRODWCWFik = 0.0d0
      DO j=1,NBF5f
       IQP = INDf(j)
       PRODWCWFik = PRODWCWFik + W(ik,IQP)*CW12(IQ,IQP)
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! PRODWCWFikq
      FUNCTION PRODWCWFikq(ik,jf,INDf,W,CW12)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/INPFILE_FRAG/NO1f,NBFf,NBF5f
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
!
      INTEGER,DIMENSION(NBFf)::INDf
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5)::CW12
      DOUBLE PRECISION,DIMENSION(NSQ,NBF5)::W
!-----------------------------------------------------------------------
!     Sum W*CW12(IQ,*) por IQP, IQ is not considered
!-----------------------------------------------------------------------
      PRODWCWFikq = 0.0d0
      IQ = INDf(jf)
      DO j=1,jf-1
       IQP = INDf(j)
       PRODWCWFikq = PRODWCWFikq + W(ik,IQP)*CW12(IQ,IQP)
      ENDDO
      DO j=jf+1,NBF5f
       IQP = INDf(j)
       PRODWCWFikq = PRODWCWFikq + W(ik,IQP)*CW12(IQ,IQP)
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! SQUARETRIAN2
      SUBROUTINE SQUARETRIAN2(FM,F,N,NT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION,DIMENSION(N,N)::FM
      DOUBLE PRECISION,DIMENSION(NT)::F
!-----------------------------------------------------------------------
!     PUT SQUARE FM IN TRIANGULAR FORM F
!-----------------------------------------------------------------------
      NZ=0 
      DO I=1,N
       DO J=1,I
        NZ=NZ+1
        F(NZ)=FM(I,J)+FM(I,J)
       ENDDO
!      DOUBLE COUNTING OF DIAGONAL TERMS IS SKIPPED IN DABCLU
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! SQUARETRIAN3
      SUBROUTINE SQUARETRIAN3(FM,F,N,NT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION,DIMENSION(N,N,N)::FM
      DOUBLE PRECISION,DIMENSION(N,NT)::F
      NZ=0
      DO I=1,N 
       DO J=1,I
        NZ=NZ+1
        F(:,NZ)=FM(:,I,J)+FM(:,J,I)
       ENDDO
!      AVOIDING DOUBLE DIAGONAL COUNTING IN DABNOF       
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! TRACEs
      SUBROUTINE TRACEs(E,AA,EEPS,N,NL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION,DIMENSION(N,N),INTENT(IN)::AA
      DOUBLE PRECISION,DIMENSION(NL),INTENT(IN)::EEPS
      DOUBLE PRECISION,INTENT(OUT)::E
      E = 0.0D+0
      KL = 0
      DO K=1,N
       DO L=1,K
        KL = KL + 1
        E = E + AA(K,L)*EEPS(KL)
       ENDDO
       E = E - AA(K,K)*EEPS(KL)*0.5d0
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! TFTRI                                            
      SUBROUTINE TFTRI(H,FMAT,TMAT,WRK,M,N,LDT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION H(*),FMAT(*),TMAT(LDT,M),WRK(N)
!-----------------------------------------------------------------------
!     Transform triangular matrix FMAT: H = TMAT+ * FMAT * TMAT
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IJ = 0
      DO J = 1,M,5
       JJMAX = MIN(M,J+4)
       DO JJ=J,JJMAX
        IK = 0
        DO I = 1,N
         IM1 = I-1
         DUM = 0.0D+00
         TDUM = TMAT(I,JJ)
         IF(IM1>0)THEN
          DO K = 1,IM1
           IK = IK+1
           WRK(K) = WRK(K)+FMAT(IK)*TDUM
           DUM = DUM+FMAT(IK)*TMAT(K,JJ)
          END DO
         END IF
         IK = IK+1
         WRK(I) = DUM+FMAT(IK)*TDUM
        END DO
        DO I = 1,JJ
         IJ = IJ+1
         HIJ = DDOT(N,TMAT(1,I),1,WRK,1)
         IF(ABS(HIJ)<1.0D-11) HIJ = 0.0D+00
         H(IJ) = HIJ
        END DO
       END DO
      END DO
!-----------------------------------------------------------------------
      RETURN
      END

! TFSQB                                            
      SUBROUTINE TFSQB(V,Q,WRK,M,N,LDQV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION V(LDQV,M),Q(LDQV,M),WRK(N)
!-----------------------------------------------------------------------
!     Back Transform: V = Q*V                            
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO J = 1,M,5
       JJMAX = MIN(M,J+4)
       DO JJ=J,JJMAX
        DO I = 1,N
         W = DDOT(M,Q(I,1),LDQV,V(1,JJ),1)
         IF(ABS(W)<1.0D-11)W = 0.0D+00
         WRK(I)=W
        END DO
        DO I = 1,N
         V(I,JJ) = WRK(I)
        END DO
       END DO
      END DO
!-----------------------------------------------------------------------
      RETURN
      END

! STFASEd        
      SUBROUTINE STFASEd(A,LDA,N,M)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER,INTENT(IN) :: LDA,M,N
      DOUBLE PRECISION,DIMENSION(LDA,M),INTENT(INOUT):: A
      INTEGER:: I,J
!-----------------------------------------------------------------------      
!     SET THE PHASE OF EACH COLUMN OF A MATRIX SO THE LARGEST
!     ELEMENT IS POSITIVE
!-----------------------------------------------------------------------
      DO I = 1,M
       LARGE = IDAMAX(N,A(1,I),1)
       IF(LARGE<=0)LARGE=1
       IF(LARGE>N)LARGE=1
       IF(A(LARGE,I)<0.0D+00)THEN
        DO J = 1,N
         A(J,I) = -A(J,I)
        END DO
       END IF
      END DO
!-----------------------------------------------------------------------            
      RETURN
      END

! CPYTSQ                                           
      SUBROUTINE CPYTSQ(A,B,NA)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(*),B(NA,NA)
!-----------------------------------------------------------------------    
!     ---- COPY TRIANGULAR A TO SQUARE B (NA BY NA) ----                
!-----------------------------------------------------------------------    
      IJ=1
      DO I=1,NA
       DO J=1,I
        B(I,J) = A(IJ)
        B(J,I) = A(IJ)
        IJ = IJ + 1
       END DO
      END DO
!-----------------------------------------------------------------------    
      RETURN
      END

! VCLR                                             
      SUBROUTINE VCLR(A,INCA,N)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(*)
!-----------------------------------------------------------------------
!     ZERO OUT VECTOR -A-, USING INCREMENT -INCA- 
!-----------------------------------------------------------------------
      IF(INCA==1)THEN
       DO L=1,N
        A(L) = 0.0d0
       ENDDO
       RETURN
      ENDIF
      LA = 1-INCA
      DO L=1,N
       LA = LA+INCA
       A(LA) = 0.0d0
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END

! TFSQU                                            
      SUBROUTINE TFSQU(H,F,T,WRK,N,M)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION H(M,M),F(N,N),T(N,M),WRK(N)
!-----------------------------------------------------------------------               
!     H = T+ * F * T
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO I = 1,M,5
       IIMAX = MIN(M,I+4)
       DO II=I,IIMAX
        DO L=1,N
         WRK(L) = DDOT(N,T(1,II),1,F(1,L),1)
        END DO
        DO J=1,M
         H(II,J) = DDOT(N,WRK,1,T(1,J),1)
        END DO
       END DO
      END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -      
      RETURN
      END


! SLVLEQ                                           
      SUBROUTINE SLVLEQ(A,X,IPVT,N,LDA,JOB,IERR)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(LDA,N),X(N),IPVT(N)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Solve Linear equations: A*X=B (JOB=0) or At*X=B (JOB/=0)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IERR=0
      CALL DGEFA(A,LDA,N,IPVT,IERR)
      IF(IERR/=0) RETURN
      CALL DGESL(A,LDA,N,IPVT,X,JOB)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -      
      RETURN
      END
! DGEFA                                            
      SUBROUTINE DGEFA(A,LDA,N,IPVT,INFO)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(LDA,*),IPVT(*)
!                                                                       
!     DGEFA FACTORS A DOUBLE PRECISION MATRIX BY GAUSSIAN ELIMINATION.  
!                                                                       
!     DGEFA IS USUALLY CALLED BY DGECO, BUT IT CAN BE CALLED            
!     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.          
!     (TIME FOR DGECO) = (1 + 9/N)*(TIME FOR DGEFA) .                   
!                                                                       
!     ON ENTRY                                                          
!                                                                       
!        A       DOUBLE PRECISION(LDA, N)                               
!                THE MATRIX TO BE FACTORED.                             
!                                                                       
!        LDA     INTEGER                                                
!                THE LEADING DIMENSION OF THE ARRAY  A .                
!                                                                       
!        N       INTEGER                                                
!                THE ORDER OF THE MATRIX  A .                           
!                                                                       
!     ON RETURN                                                         
!                                                                       
!        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS         
!                WHICH WERE USED TO OBTAIN IT.                          
!                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE       
!                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER          
!                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.       
!                                                                       
!        IPVT    INTEGER(N)                                             
!                AN INTEGER VECTOR OF PIVOT INDICES.                    
!                                                                       
!        INFO    INTEGER                                                
!                = 0  NORMAL VALUE.                                     
!                = K  IF  U(K,K) == 0.0 .  THIS IS NOT AN ERROR       
!                     CONDITION FOR THIS ROUTINE, BUT IT DOES           
!                     INDICATE THAT DGESL OR DGEDI WILL DIVIDE BY ZERO  
!                     IF CALLED.  USE  RCOND  IN DGECO FOR A RELIABLE   
!                     INDICATION OF SINGULARITY.                        
!                                                                       
!     LINPACK. THIS VERSION DATED 08/14/78 .                            
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      
!                                                                       
!     ROUTINES AND FUNCTIONS                                            
!                                                                       
!     BLAS DAXPY,DSCAL,IDAMAX                                           
!                                                                       
!     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING                        
!                                                                       
      INFO = 0
      NM1 = N - 1
      IF (NM1 < 1) GO TO 70
      DO 60 K = 1, NM1
         KP1 = K + 1
!                                                                       
!        FIND L = PIVOT INDEX                                           
!                                                                       
         L = IDAMAX(N-K+1,A(K,K),1) + K - 1
         IPVT(K) = L
!                                                                       
!        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED          
!                                                                       
         IF (A(L,K) == 0.0D+00) GO TO 40
!                                                                       
!           INTERCHANGE IF NECESSARY                                    
!                                                                       
            IF (L == K) GO TO 10
               T = A(L,K)
               A(L,K) = A(K,K)
               A(K,K) = T
   10       CONTINUE
!                                                                       
!           COMPUTE MULTIPLIERS                                         
!                                                                       
            T = -1.0D+00/A(K,K)
            CALL DSCAL(N-K,T,A(K+1,K),1)
!                                                                       
!           ROW ELIMINATION WITH COLUMN INDEXING                        
!                                                                       
            DO 30 J = KP1, N
               T = A(L,J)
               IF (L == K) GO TO 20
                  A(L,J) = A(K,J)
                  A(K,J) = T
   20          CONTINUE
               CALL DAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)
   30       CONTINUE
         GO TO 50
   40    CONTINUE
            INFO = K
   50    CONTINUE
   60 CONTINUE
   70 CONTINUE
      IPVT(N) = N
      IF (A(N,N) == 0.0D+00) INFO = N
      RETURN
      END

! DGESL                                            
      SUBROUTINE DGESL(A,LDA,N,IPVT,B,JOB)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(LDA,N),B(N),IPVT(N)
!                                                                       
!     DGESL SOLVES THE DOUBLE PRECISION SYSTEM                          
!     A * X = B  OR  TRANS(A) * X = B                                   
!     USING THE FACTORS COMPUTED BY DGECO OR DGEFA.                     
!                                                                       
!     ON ENTRY                                                          
!                                                                       
!        A       DOUBLE PRECISION(LDA, N)                               
!                THE OUTPUT FROM DGECO OR DGEFA.                        
!                                                                       
!        LDA     INTEGER                                                
!                THE LEADING DIMENSION OF THE ARRAY  A .                
!                                                                       
!        N       INTEGER                                                
!                THE ORDER OF THE MATRIX  A .                           
!                                                                       
!        IPVT    INTEGER(N)                                             
!                THE PIVOT VECTOR FROM DGECO OR DGEFA.                  
!                                                                       
!        B       DOUBLE PRECISION(N)                                    
!                THE RIGHT HAND SIDE VECTOR.                            
!                                                                       
!        JOB     INTEGER                                                
!                = 0         TO SOLVE  A*X = B ,                        
!                = NONZERO   TO SOLVE  TRANS(A)*X = B  WHERE            
!                            TRANS(A)  IS THE TRANSPOSE.                
!                                                                       
!     ON RETURN                                                         
!                                                                       
!        B       THE SOLUTION VECTOR  X .                               
!                                                                       
!     ERROR CONDITION                                                   
!                                                                       
!        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A   
!        ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES SINGULARITY  
!        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER       
!        SETTING OF LDA .  IT WILL NOT OCCUR IF THE ROUTINES ARE        
!        CALLED CORRECTLY AND IF DGECO HAS SET RCOND > 0.0           
!        OR DGEFA HAS SET INFO == 0 .                                 
!                                                                       
!     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX                 
!     WITH  P  COLUMNS                                                  
!           CALL DGECO(A,LDA,N,IPVT,RCOND,Z)                            
!           IF (RCOND IS TOO SMALL) GO TO ...                           
!           DO 10 J = 1, P                                              
!              CALL DGESL(A,LDA,N,IPVT,C(1,J),0)                        
!        10 CONTINUE                                                    
!                                                                       
!     LINPACK. THIS VERSION DATED 08/14/78 .                            
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      
!                                                                       
!     ROUTINES AND FUNCTIONS                                            
!                                                                       
!     BLAS DAXPY,DDOT                                                   
!                                                                       
      NM1 = N - 1
      IF (JOB /= 0) GO TO 50
!                                                                       
!        JOB = 0 , SOLVE  A * X = B                                     
!        FIRST SOLVE  L*Y = B                                           
!                                                                       
         IF (NM1 < 1) GO TO 30
         DO 20 K = 1, NM1
            L = IPVT(K)
            T = B(L)
            IF (L == K) GO TO 10
               B(L) = B(K)
               B(K) = T
   10       CONTINUE
            CALL DAXPY(N-K,T,A(K+1,K),1,B(K+1),1)
   20    CONTINUE
   30    CONTINUE
!                                                                       
!        NOW SOLVE  U*X = Y                                             
!                                                                       
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K)/A(K,K)
            T = -B(K)
            CALL DAXPY(K-1,T,A(1,K),1,B(1),1)
   40    CONTINUE
      GO TO 100
   50 CONTINUE
!                                                                       
!        JOB = NONZERO, SOLVE  TRANS(A) * X = B                         
!        FIRST SOLVE  TRANS(U)*Y = B                                    
!                                                                       
         DO 60 K = 1, N
            T = DDOT(K-1,A(1,K),1,B(1),1)
            B(K) = (B(K) - T)/A(K,K)
   60    CONTINUE
!                                                                       
!        NOW SOLVE TRANS(L)*X = Y                                       
!                                                                       
         IF (NM1 < 1) GO TO 90
         DO 80 KB = 1, NM1
            K = N - KB
            B(K) = B(K) + DDOT(N-K,A(K+1,K),1,B(K+1),1)
            L = IPVT(K)
            IF (L == K) GO TO 70
               T = B(L)
               B(L) = B(K)
               B(K) = T
   70       CONTINUE
   80    CONTINUE
   90    CONTINUE
  100 CONTINUE
      RETURN
      END

! DGEDId
      SUBROUTINE DGEDId(A,LDA,N,IPVT,DET,WORK,JOB)                       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      DIMENSION A(LDA,*),DET(2),WORK(*),IPVT(*)                         
!                                                                       
!     DGEDI COMPUTES THE DETERMINANT AND INVERSE OF A MATRIX            
!     USING THE FACTORS COMPUTED BY DGECO OR DGEFA.                     
!                                                                       
!     ON ENTRY                                                          
!                                                                       
!        A       DOUBLE PRECISION(LDA, N)                               
!                THE OUTPUT FROM DGECO OR DGEFA.                        
!                                                                       
!        LDA     INTEGER                                                
!                THE LEADING DIMENSION OF THE ARRAY  A .                
!                                                                       
!        N       INTEGER                                                
!                THE ORDER OF THE MATRIX  A .                           
!                                                                       
!        IPVT    INTEGER(N)                                             
!                THE PIVOT VECTOR FROM DGECO OR DGEFA.                  
!                                                                       
!        WORK    DOUBLE PRECISION(N)                                    
!                WORK VECTOR.  CONTENTS DESTROYED.                      
!                                                                       
!        JOB     INTEGER                                                
!                = 11   BOTH DETERMINANT AND INVERSE.                   
!                = 01   INVERSE ONLY.                                   
!                = 10   DETERMINANT ONLY.                               
!                                                                       
!     ON RETURN                                                         
!                                                                       
!        A       INVERSE OF ORIGINAL MATRIX IF REQUESTED.               
!                OTHERWISE UNCHANGED.                                   
!                                                                       
!        DET     DOUBLE PRECISION(2)                                    
!                DETERMINANT OF ORIGINAL MATRIX IF REQUESTED.           
!                OTHERWISE NOT REFERENCED.                              
!                DETERMINANT = DET(1) * 10.0**DET(2)                    
!                WITH  1.0 .LE. ABS(DET(1)) .LT. 10.0                   
!                OR  DET(1) .EQ. 0.0 .                                  
!                                                                       
!     ERROR CONDITION                                                   
!                                                                       
!        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS     
!        A ZERO ON THE DIAGONAL AND THE INVERSE IS REQUESTED.           
!        IT WILL NOT OCCUR IF THE ROUTINES ARE CALLED CORRECTLY         
!        AND IF DGECO HAS SET RCOND .GT. 0.0 OR DGEFA HAS SET           
!        INFO .EQ. 0 .                                                  
!                                                                       
!     LINPACK. THIS VERSION DATED 08/14/78 .                            
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      
!                                                                       
!     ROUTINES AND FUNCTIONS                                            
!                                                                       
!     BLAS DAXPY,DSCAL,DSWAP                                            
!     FORTRAN ABS,MOD                                                   
!                                                                       
!     COMPUTE DETERMINANT                                               
!                                                                       
      IF (JOB/10 .EQ. 0) GO TO 70                                       
         DET(1) = 1.0D+00                                               
         DET(2) = 0.0D+00                                               
         TEN = 10.0D+00                                                 
         DO 50 I = 1, N                                                 
            IF (IPVT(I) .NE. I) DET(1) = -DET(1)                        
            DET(1) = A(I,I)*DET(1)                                      
!        ...EXIT                                                        
            IF (DET(1) .EQ. 0.0D+00) GO TO 60                           
   10       IF (ABS(DET(1)) .GE. 1.0D+00) GO TO 20                      
               DET(1) = TEN*DET(1)                                      
               DET(2) = DET(2) - 1.0D+00                                
            GO TO 10                                                    
   20       CONTINUE                                                    
   30       IF (ABS(DET(1)) .LT. TEN) GO TO 40                          
               DET(1) = DET(1)/TEN                                      
               DET(2) = DET(2) + 1.0D+00                                
            GO TO 30                                                    
   40       CONTINUE                                                    
   50    CONTINUE                                                       
   60    CONTINUE                                                       
   70 CONTINUE                                                          
!                                                                       
!     COMPUTE INVERSE(U)                                                
!                                                                       
      IF (MOD(JOB,10) .EQ. 0) GO TO 150                                 
         DO 100 K = 1, N                                                
            A(K,K) = 1.0D+00/A(K,K)                                     
            T = -A(K,K)                                                 
            CALL DSCAL(K-1,T,A(1,K),1)                                  
            KP1 = K + 1                                                 
            IF (N .LT. KP1) GO TO 90                                    
            DO 80 J = KP1, N                                            
               T = A(K,J)                                               
               A(K,J) = 0.0D+00                                         
               CALL DAXPY(K,T,A(1,K),1,A(1,J),1)                        
   80       CONTINUE                                                    
   90       CONTINUE                                                    
  100    CONTINUE                                                       
!                                                                       
!        FORM INVERSE(U)*INVERSE(L)                                     
!                                                                       
         NM1 = N - 1                                                    
         IF (NM1 .LT. 1) GO TO 140                                      
         DO 130 KB = 1, NM1                                             
            K = N - KB                                                  
            KP1 = K + 1                                                 
            DO 110 I = KP1, N                                           
               WORK(I) = A(I,K)                                         
               A(I,K) = 0.0D+00                                         
  110       CONTINUE                                                    
            DO 120 J = KP1, N                                           
               T = WORK(J)                                              
               CALL DAXPY(N,T,A(1,J),1,A(1,K),1)                        
  120       CONTINUE                                                    
            L = IPVT(K)                                                 
            IF (L .NE. K) CALL DSWAP(N,A(1,K),1,A(1,L),1)               
  130    CONTINUE                                                       
  140    CONTINUE                                                       
  150 CONTINUE                                                          
      RETURN                                                            
      END

! M_inv
      SUBROUTINE M_inv(GMAT,N,IPRINTOPT)
      INTEGER :: N,I,J,IDX
      DOUBLE PRECISION :: GMAT(N,N),TOL
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: VEC,AUX
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: EIG,W

      ALLOCATE(VEC(N,N),AUX(N,N))
      ALLOCATE(EIG(N),W(N))
!
      TOL = 1.0d-5
      AUX = GMAT

      CALL DIAG(N,AUX,VEC,EIG,W)

!
      IDX = 0
      AUX = 0.0d0
      DO I=1,N
       IF (ABS(EIG(I))<=TOL) THEN
        EIG(I) = SIGN(TOL,EIG(I))
        IF(EIG(I)==0) EIG(I) = TOL
        IDX = IDX + 1
       END IF
       AUX(I,I) = 1.0/EIG(I)
      END DO

      IF(IPRINTOPT==1) write(6,*)                                       &
       "Number of corrected eigenvalues at inverse",IDX

      DEALLOCATE(EIG,W)
      GMAT = MATMUL(MATMUL(VEC,AUX),TRANSPOSE(VEC))
      DEALLOCATE(VEC,AUX)

      END

      SUBROUTINE SORT_IDX_DECREASING(VEC0, SORT_IDX, N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER :: N
      INTEGER,DIMENSION(N)::SORT_IDX
      DOUBLE PRECISION,DIMENSION(N)::VEC0

      INTEGER :: I, J, MAX_IDX, ITEMP
      DOUBLE PRECISION :: DTEMP
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::VEC

      ALLOCATE(VEC(N))
      VEC=VEC0

      DO I=1,N
        SORT_IDX(I) = I
      END DO
  
      DO I=1, N-1
        MAX_IDX = I
        DO J=I+1,N
          IF (VEC(J) > VEC(MAX_IDX)) THEN
            MAX_IDX = J
          END IF
        END DO

        DTEMP = VEC(I)
        VEC(I) = VEC(MAX_IDX)
        VEC(MAX_IDX) = DTEMP
        ITEMP = SORT_IDX(I)
        SORT_IDX(I) = SORT_IDX(MAX_IDX)
        SORT_IDX(MAX_IDX) = ITEMP
      END DO

      DEALLOCATE(VEC)

      END

