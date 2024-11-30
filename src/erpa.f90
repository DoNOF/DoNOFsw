!======================================================================!
!                                                                      !
!                    E R P A   S U B R O U T I N E S                   !
!                                                                      !
!             2024  Module implemented by J.F.H. Lew-Yee,              !
!          I.A. Bonfil-Rivera, M. Piris, J. Martin del Campo           !
!                                                                      !
!            J. Chem. Theory Comput. 20(5), 2140–2151 (2024)           !
!                                                                      !
!======================================================================!
!                                                                      !
!   PNOF_ERPA     : Compute ERPA0/ERPA1/ERPA2                          !
!   Compute_2RDM  : Daa and Dab of PNOF5/PNOF7/GNOF                    !
!   DIAG_ERPA     : Diagonalize EPRA matrices and print exc. energies  ! 
!                                                                      !
!======================================================================!

! PNOF_ERPA
      SUBROUTINE PNOF_ERPA(RO,COEF,AHCORE,XIJKAUX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/ERIACT/ERIACTIVATED,NIJKaux,NINTCRaux,NSTOREaux,IAUXDIM
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
      DOUBLE PRECISION,DIMENSION(NBF*(NBF+1)/2,IAUXDIM)::XIJKAUX
      DOUBLE PRECISION,DIMENSION(NBF5) :: RO
      DOUBLE PRECISION,DIMENSION(NBF,NBF) :: COEF,AHCORE

      INTEGER I, J, K, L, P, Q, R, S, T, U, W
      INTEGER MAX_IDX, N_S, N_NS, ITMP
      DOUBLE PRECISION :: DTMP, DTMP2
      INTEGER,ALLOCATABLE,DIMENSION(:) :: dN_SORT_IDX
      INTEGER,ALLOCATABLE,DIMENSION(:) :: SIGNIF, NO_SIGNIF
      INTEGER,ALLOCATABLE,DIMENSION(:) :: EMPTY
      INTEGER,ALLOCATABLE,DIMENSION(:) :: M_SORT_IDX

      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: dN, dNm1
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: VV
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: ROsqrt
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: Hij, C_NBF5
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: M_ERPA, M_ERPA0
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: M_ERPA1, M_ERPA2
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: AA, BB, ApB, AmB
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: C, E, F, CFm1E
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: B_IN, B_IJ
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:,:) :: Daa, Dab
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:,:) :: Xijkl
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:,:) :: A

!-----------------------------------------------------------------------
!     This code solve the ERPA equations using the PNOFi (i=5,7,8) RDMs
!     requiring to solve the eigenproblem:
!
!     ( A_rspq    B_rspq    A_rspp) (Xpq)     ( dN_rs    0    0)(Xpq)
!     ( B_rspq    A_rspq    B_rspp) (Ypq) = w (   0   -dN_rs  0)(Ypq)
!     ( A_rrpq    A_rrqp    A_rrpp) (Zp )     (   0      0    0)(Zp )
!
!     In ERPA0, only Xpq and Ypq are used.
!     In ERPA1, Zp is used.
!     In ERPA2, Zp is replaced with Vp, and the matrices are modified.
!
!     The following variables are used in the code
!
!     Hij:        One electron core-hamiltonian integrals in NO Basis
!     A:          Array for A_rspq
!     B:          In literature, B is defined as B_rspq = A_rsqp
!     M_ERPA:     The ERPA matrix at the left side of the
!                 ERPA Generalized eigenvalue problem
!     dN:         n_s - n_r
!     Daa:        D^(alpha alpha) of the PNOF
!     Dab:        D^(alpha beta) of the PNOF
!     B_IN, B_IJ: Arrays for AO->NO ERIs transformation
!     Xijkl:      NO ERIs (currently created from RI, but reusing the
!                 original array is work in progress)
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Transform Integrals
!-----------------------------------------------------------------------
      ALLOCATE(Hij(NBF5,NBF5), C_NBF5(NBF,NBF5))
      C_NBF5 = COEF(:,1:NBF5)
      Hij = MATMUL(TRANSPOSE(C_NBF5),MATMUL(AHCORE,C_NBF5))
      DEALLOCATE(C_NBF5)

      ALLOCATE(Xijkl(NBF5,NBF5,NBF5,NBF5))
      ALLOCATE(B_IN(NBF5,NBF), B_IJ(NBF5,NBF5))
      Xijkl(1:NBF5,1:NBF5,1:NBF5,1:NBF5) = 0.0D0

      DO K=1,IAUXDIM
        B_IN(1:NBF5,1:NBF) = 0.0D0
        DO I=1,NBF5
          DO N=1,NBF
            DO M=1,N
              MN = M + N*(N-1)/2
              B_IN(I,N) = B_IN(I,N) + COEF(M,I)*XIJKAUX(MN,K)
              IF(M.NE.N) B_IN(I,M) = B_IN(I,M)+COEF(N,I)*XIJKAUX(MN,K)
            END DO
          END DO
        END DO

        B_IJ(1:NBF5,1:NBF5) = 0.0D0
        DO J=1,NBF5
          DO I=1,J
            DO N=1,NBF
              B_IJ(I,J) = B_IJ(I,J) + COEF(N,J)*B_IN(I,N)
            END DO
            B_IJ(J,I) = B_IJ(I,J)
          END DO
        END DO

        DO I=1,NBF5
          DO J=1,NBF5
            DO KK=1,NBF5
              DO L=1,NBF5
                Xijkl(I,KK,J,L) = Xijkl(I,KK,J,L) + B_IJ(I,J)*B_IJ(KK,L)
              END DO
            END DO
          END DO
        END DO

      END DO
      DEALLOCATE(B_IN, B_IJ)

      ALLOCATE(A(NBF5,NBF5,NBF5,NBF5))
      A(1:NBF5,1:NBF5,1:NBF5,1:NBF5) = 0.0D0
      DO P=1,NBF5
        DO Q=1,NBF5
          DO S=1,NBF5
            A(P,S,P,Q) = A(P,S,P,Q) + Hij(S,Q)*RO(P)
            A(P,S,P,Q) = A(P,S,P,Q) - Hij(S,Q)*RO(S)
          END DO
          DO R=1,NBF5
            A(R,Q,P,Q) = A(R,Q,P,Q) + Hij(P,R)*RO(Q)
            A(R,Q,P,Q) = A(R,Q,P,Q) - Hij(P,R)*RO(R)
          END DO
        END DO
      END DO
      DEALLOCATE(Hij)

!-----------------------------------------------------------------------
!  Compute PNOF 2RDM (Daa and Dab subblocks)
!-----------------------------------------------------------------------
      ALLOCATE(Daa(NBF5,NBF5,NBF5,NBF5),Dab(NBF5,NBF5,NBF5,NBF5))
      CALL Compute_2RDM(RO,Daa,Dab)

!-----------------------------------------------------------------------
!  Build A_rspq (Eq. 10) J. Chem. Theory Comput. 20(5), 2140–2151 (2024)
!-----------------------------------------------------------------------
      DO R=1,NBF5
        DO S=1,NBF5
          DO P=1,NBF5
            DO Q=1,NBF5
              DO T=1,NBF5
                DO U=1,NBF5
                  A(R,S,P,Q) = A(R,S,P,Q) + Xijkl(S,T,Q,U) * Daa(P,U,R,T)
                  A(R,S,P,Q) = A(R,S,P,Q) - Xijkl(S,T,U,Q) * Daa(P,U,R,T)
                  A(R,S,P,Q) = A(R,S,P,Q) + Xijkl(S,T,Q,U) * Dab(P,U,R,T)
                  A(R,S,P,Q) = A(R,S,P,Q) + Xijkl(S,T,U,Q) * Dab(U,P,R,T)
                  !----
                  A(R,S,P,Q) = A(R,S,P,Q) + Xijkl(U,P,T,R) * Daa(S,T,Q,U)
                  A(R,S,P,Q) = A(R,S,P,Q) - Xijkl(U,P,R,T) * Daa(S,T,Q,U)
                  A(R,S,P,Q) = A(R,S,P,Q) + Xijkl(U,P,T,R) * Dab(S,T,Q,U)
                  A(R,S,P,Q) = A(R,S,P,Q) + Xijkl(U,P,R,T) * Dab(S,T,U,Q)
                  !----
                  A(R,S,P,Q) = A(R,S,P,Q) + Xijkl(P,S,T,U) * Daa(T,U,R,Q)
                  A(R,S,P,Q) = A(R,S,P,Q) - Xijkl(P,S,T,U) * Dab(U,T,R,Q)
                  !----
                  A(R,S,P,Q) = A(R,S,P,Q) + Xijkl(T,U,Q,R) * Daa(S,P,T,U)
                  A(R,S,P,Q) = A(R,S,P,Q) - Xijkl(T,U,Q,R) * Dab(P,S,T,U)
                END DO
              END DO
            END DO
          END DO
        END DO 
      END DO

      DO R=1,NBF5
        DO P=1,NBF5
          DTMP  = 0.0D0
          DTMP2 = 0.0D0
          DO T=1,NBF5
            DO U=1,NBF5
              DO W=1,NBF5
                DTMP  = DTMP  + Xijkl(T,P,W,U) * Daa(W,U,R,T)
                DTMP2 = DTMP2 + Xijkl(T,P,W,U) * Dab(U,W,R,T)
              END DO
            END DO
          END DO
          DO S=1,NBF5
            A(R,S,P,S) = A(R,S,P,S) + DTMP 
            A(R,S,P,S) = A(R,S,P,S) - DTMP2
          END DO
        END DO
      END DO

      DO S=1,NBF5
        DO Q=1,NBF5
          DTMP  = 0.0D0
          DTMP2 = 0.0D0
          DO T=1,NBF5
            DO U=1,NBF5
              DO W=1,NBF5
                DTMP  = DTMP  + Xijkl(T,U,W,Q) * Daa(S,W,T,U)
                DTMP2 = DTMP2 + Xijkl(T,U,W,Q) * Dab(W,S,T,U)
              END DO
            END DO
          END DO
          DO R=1,NBF5
            A(R,S,R,Q) = A(R,S,R,Q) + DTMP 
            A(R,S,R,Q) = A(R,S,R,Q) - DTMP2
          END DO
        END DO
      END DO

      DEALLOCATE(Xijkl)
      DEALLOCATE(Daa,Dab)

!-----------------------------------------------------------------------
!  Build M_ERPA
!-----------------------------------------------------------------------
      ALLOCATE(M_ERPA(NBF5**2,NBF5**2))
      M_ERPA = 0.0D0
      I = 0
      ! LHS Eq.(21) J. Chem. Theory Comput. 20(5), 2140–2151 (2024).
      DO S=1,NBF5
        DO R=S+1,NBF5
          I = I + 1
          J = 0
          DO Q=1,NBF5
            DO P=Q+1,NBF5
              J = J + 1
              M_ERPA(I,J) = A(R,S,P,Q)
            END DO
          END DO
          DO Q=1,NBF5
            DO P=Q+1,NBF5
              J = J + 1
              M_ERPA(I,J) = A(R,S,Q,P)
            END DO
          END DO
          DO P=1,NBF5
            J = J + 1
            M_ERPA(I,J) = A(R,S,P,P)
          END DO
        END DO
      END DO
      ! LHS Eq.(22) J. Chem. Theory Comput. 20(5), 2140–2151 (2024).
      DO S=1,NBF5
        DO R=S+1,NBF5
          I = I + 1
          J = 0
          DO Q=1,NBF5
            DO P=Q+1,NBF5
              J = J + 1
              M_ERPA(I,J) = A(R,S,Q,P)
            END DO
          END DO
          DO Q=1,NBF5
            DO P=Q+1,NBF5
              J = J + 1
              M_ERPA(I,J) = A(R,S,P,Q)
            END DO
          END DO
          DO P=1,NBF5
            J = J + 1
            M_ERPA(I,J) = A(R,S,P,P)
          END DO
        END DO
      END DO
      ! LHS Eq.(23) J. Chem. Theory Comput. 20(5), 2140–2151 (2024).
      DO R=1,NBF5
        I = I + 1
        J = 0
        DO Q=1,NBF5
          DO P=Q+1,NBF5
            J = J + 1
            M_ERPA(I,J) = A(R,R,P,Q)
          END DO
        END DO
        DO Q=1,NBF5
          DO P=Q+1,NBF5
            J = J + 1
            M_ERPA(I,J) = A(R,R,Q,P)
          END DO
        END DO
        DO P=1,NBF5
          J = J + 1
          M_ERPA(I,J) = A(R,R,P,P)
        END DO
      END DO

!-----------------------------------------------------------------------
!  Build dN
!-----------------------------------------------------------------------
      ALLOCATE(dN(NBF5*(NBF5-1)/2))
      I = 0
      DO S=1,NBF5
        DO R=S+1,NBF5
          I = I + 1
          dN(I) = +(RO(S) - RO(R))
        END DO
      END DO

!-----------------------------------------------------------------------
!  Sort dN in decreasing order of absolute values
!-----------------------------------------------------------------------
      ALLOCATE(dN_SORT_IDX(NBF5*(NBF5-1)/2))
      CALL SORT_IDX_DECREASING(ABS(dN), dN_SORT_IDX, NBF5*(NBF5-1)/2)
      dN = dN(dN_SORT_IDX)

      ! Keep elements of dN > TOL (cut small values of dN) 
      N_S = NBF5*(NBF5-1)/2
      TOL = 1E-15
      DO I=1, NBF5*(NBF5-1)/2
        IF(ABS(dN(I)) > TOL) N_S = I !Index of the last allowed dN
      END DO
      N_NS = NBF5*(NBF5-1)/2 - N_S

!-----------------------------------------------------------------------
!  Create a new axis array with elements given by:
!  M_SORT_IDX: | Signif | Signif | Empty | No_Signif | No_Signif |
!  Signif    : Values of dN >  TOL
!  Empty     : Elments of V_ERPA with zeros in the diagonal
!  No_Signif : Values of dN <= TOL
!-----------------------------------------------------------------------
      ALLOCATE(SIGNIF(N_S), NO_SIGNIF(N_NS))
      ALLOCATE(EMPTY(NBF5))
      SIGNIF = dN_SORT_IDX(1:N_S)
      NO_SIGNIF = dN_SORT_IDX(N_S+1:NBF5*(NBF5-1)/2)
      DO I=1,NBF5
        EMPTY(I) = NBF5*(NBF5-1)+I
      END DO
      DEALLOCATE(dN_SORT_IDX)

      ALLOCATE(M_SORT_IDX(NBF5**2))
      M_SORT_IDX(1:N_S) = SIGNIF
      M_SORT_IDX(N_S+1:2*N_S) = NBF5*(NBF5-1)/2 + SIGNIF
      M_SORT_IDX(2*N_S+1:2*N_S+NBF5) = EMPTY
      M_SORT_IDX(2*N_S+NBF5+1:2*N_S+NBF5+N_NS) = NO_SIGNIF
      M_SORT_IDX(2*N_S+NBF5+N_NS+1:NBF5**2) = NBF5*(NBF5-1)/2 + NO_SIGNIF
      DEALLOCATE(SIGNIF, NO_SIGNIF, EMPTY)
!-----------------------------------------------------------------------
!  Permute M_ERPA matrix as given by M_SORT_IDX axis
!-----------------------------------------------------------------------
      M_ERPA = M_ERPA(M_SORT_IDX,M_SORT_IDX)

!-----------------------------------------------------------------------
!                               ERPA0
!-----------------------------------------------------------------------
! Take the operator as:
!
! O^t =   \sum_{p>q} X_{pq} (a_{p_a}^t a_{q_a} + a_{p_b}^t a_{q_b})
!       + \sum_{p>q} Y_{pq} (a_{q_a}^t a_{p_a} + a_{q_b}^t a_{p_b})
!
! Here we solve the Generalized eigenproblem:
!
!     ( A_rspq    B_rspq) (Xpq)     ( (n_r - n_s)     0     )(Xpq)
!     ( B_rspq    A_rspq) (Ypq) = w (      0    -(n_r - n_s))(Ypq)
!
! Lets solve this as:
!
!      (dN)^-1 (A+B) (dN)^-1 (A-B) (X-Y) = w^2 dN (X-Y)
!
! Note. we are removing columns and rows asociated with small dN values
!-----------------------------------------------------------------------

      ALLOCATE(AA(N_S, N_S), BB(N_S, N_S))
      ALLOCATE(ApB(N_S, N_S), AmB(N_S, N_S))

      AA = M_ERPA(1:N_S, 1:N_S)
      BB = M_ERPA(1:N_S, N_S+1:2*N_S)

      ApB = AA + BB
      AmB = AA - BB
      DEALLOCATE(AA, BB)

      ALLOCATE(dNm1(N_S))

      dNm1 = 1.0D0 / dN(1:N_S)
      DEALLOCATE(dN)

      ALLOCATE(M_ERPA0(N_S,N_S))
      M_ERPA0 = 0.0D0
      DO I=1, N_S
        DO J=1, N_S
          DTMP = dNm1(I)*ApB(I,J)*dNm1(J)
          DO K=1, N_S
            M_ERPA0(I,K) = M_ERPA0(I,K) + DTMP*AmB(J,K)
          END DO
        END DO
      END DO

      WRITE(6,1)
      CALL DIAG_ERPA(M_ERPA0, N_S, 0)
      DEALLOCATE(M_ERPA0)

!-----------------------------------------------------------------------
!                               ERPA1
!-----------------------------------------------------------------------
! Take the operator as:
!
! O^t =   \sum_{p>q} X_{pq} (a_{p_a}^t a_{q_a} + a_{p_b}^t a_{q_b})
!       + \sum_{p>q} Y_{pq} (a_{q_a}^t a_{p_a} + a_{q_b}^t a_{p_b})
!       + \sum_{p}   Z_{p}  (a_{p_a}^t a_{p_a} + a_{p_b}^t a_{p_b})
!
! We have to solve the Generalized eigenvalue problem:
!
!     ( A_rspq    B_rspq    A_rspp) (Xpq)     (dN_rs   0    0)(Xpq)
!     ( B_rspq    A_rspq    B_rspp) (Ypq) = w (  0  -dN_rs  0)(Ypq)
!     ( A_rrpq    B_rrpq    A_rrpp) (Zp )     (  0     0    0)(Zp )
!
! Lets formulate this as
!
!     ( A   B   C ) (X)     ( dN   0   0)(X)
!     ( B   A   C ) (Y) = w ( 0  -dN   0)(Y)
!     ( E   E   F ) (Z)     ( 0    0   0)(Z)
!
! Solve for
!
!      (dN)^-1 (A+B - C F^-1 E) (dN)^-1 (A-B) (X-Y) = w^2 dN (X-Y)
!
!-----------------------------------------------------------------------

      ALLOCATE(C(N_S,NBF5), E(NBF5,N_S), F(NBF5,NBF5))
      C = M_ERPA(1:N_S,2*N_S+1:2*N_S+NBF5)
      E = M_ERPA(2*N_S+1:2*N_S+NBF5,1:N_S)
      F = M_ERPA(2*N_S+1:2*N_S+NBF5,2*N_S+1:2*N_S+NBF5)
      DEALLOCATE(M_ERPA)

      CALL M_inv(F,NBF5,0)

      ALLOCATE(CFm1E(N_S, N_S))
      CFm1E = 2*MATMUL(MATMUL(C,F),E)
      DEALLOCATE(C, E, F)

      ALLOCATE(M_ERPA1(N_S, N_S))
      M_ERPA1 = 0.0D0
      DO I=1, N_S
        DO J=1, N_S
          DTMP = dNm1(I)*(ApB(I,J)-CFm1E(I,J))*dNm1(J)
          DO K=1, N_S
            M_ERPA1(I,K) = M_ERPA1(I,K) + DTMP*AmB(J,K)
          END DO
        END DO
      END DO
      DEALLOCATE(ApB, AmB)
      DEALLOCATE(CFm1E, dNm1)

      WRITE(6,2)
      CALL DIAG_ERPA(M_ERPA1, N_S, 1)
      DEALLOCATE(M_ERPA1)

!-----------------------------------------------------------------------
!                               ERPA2
!-----------------------------------------------------------------------
! Take the operator as:
!
! O^t =   \sum_{p>q} X_{pq} (a_{p_a}^t a_{q_a} + a_{p_b}^t a_{q_b})
!       + \sum_{p>q} Y_{pq} (a_{q_a}^t a_{p_a} + a_{q_b}^t a_{p_b})
!       + \sum_{pq}  D_{pq} (a_{p_b}^t a_{q_b} a_{p_a}^t a_{p_a})
!
! We have to solve the Generalized eigenvalue problem:
!
!(A_rspq        B_rspq        c_s/2c_p(c_r+c_s)*C_rsp)(Xpq)  (dN_rs  0   0)(Xpq)
!(B_rspq        A_rspq        c_r/2c_p(c_r+c_s)*C_rsp)(Ypq)=w(  0 -dN_rs 0)(Ypq)
!((1/c_r)*E_rpq (1/c_r)*E_rpq (1/(4c_p*c_r))F_pq     )(Vp )  (  0    0   1)(Vp )
!
!-----------------------------------------------------------------------

      ALLOCATE(ROsqrt(NBF5))
      ROsqrt(1:NDNS) = SQRT(RO(1:NDNS))
      ROsqrt(NDNS+1:NBF5) = -SQRT(RO(NDNS+1:NBF5))
      
      ALLOCATE(M_ERPA(NBF5**2,NBF5**2))
      M_ERPA = 0.0D0
      I = 0
      DO S=1,NBF5
        DO R=S+1,NBF5
          I = I + 1
          J = 0
          DO Q=1,NBF5
            DO P=Q+1,NBF5
              J = J + 1
              M_ERPA(I,J) = A(R,S,P,Q)
            END DO
          END DO
          DO Q=1,NBF5
            DO P=Q+1,NBF5
              J = J + 1
              M_ERPA(I,J) = A(R,S,Q,P)
            END DO
          END DO
          DO P=1,NBF5
            J = J + 1
            M_ERPA(I,J) = ROsqrt(S)/(2*ROsqrt(P)*(ROsqrt(R)+ROsqrt(S)))*A(R,S,P,P)
          END DO
        END DO
      END DO
      !
      DO S=1,NBF5
        DO R=S+1,NBF5
          I = I + 1
          J = 0
          DO Q=1,NBF5
            DO P=Q+1,NBF5
              J = J + 1
              M_ERPA(I,J) = A(R,S,Q,P)
            END DO
          END DO
          DO Q=1,NBF5
            DO P=Q+1,NBF5
              J = J + 1
              M_ERPA(I,J) = A(R,S,P,Q)
            END DO
          END DO
          DO P=1,NBF5
            J = J + 1
            M_ERPA(I,J) = ROsqrt(R)/(2*ROsqrt(P)*(ROsqrt(R)+ROsqrt(S)))*A(R,S,P,P)
          END DO
        END DO
      END DO
      !
      DO R=1,NBF5
        I = I + 1
        J = 0
        DO Q=1,NBF5
          DO P=Q+1,NBF5
            J = J + 1
            M_ERPA(I,J) = 1/ROsqrt(R)*A(R,R,P,Q)
          END DO
        END DO
        DO Q=1,NBF5
          DO P=Q+1,NBF5
            J = J + 1
            M_ERPA(I,J) = 1/ROsqrt(R)*A(R,R,P,Q)
          END DO
        END DO
        DO P=1,NBF5
          J = J + 1
          M_ERPA(I,J) = 1/(4*ROsqrt(P)*ROsqrt(R))*A(R,R,P,P)
        END DO
      END DO
      DEALLOCATE(A, ROsqrt)

!-----------------------------------------------------------------------
!  Build dN
!-----------------------------------------------------------------------
      ALLOCATE(dN(NBF5*(NBF5-1)/2))
      I = 0
      DO S=1,NBF5
        DO R=S+1,NBF5
          I = I + 1
          dN(I) = +(RO(S) - RO(R))
        END DO
      END DO

      ALLOCATE(VV(NBF5**2))
      VV(1:NBF5*(NBF5-1)/2)= dN
      VV(NBF5*(NBF5-1)/2+1:NBF5*(NBF5-1)) = -dN
      VV(NBF5*(NBF5-1)+1:NBF5**2) = 1
      VV = VV(M_SORT_IDX)
      DEALLOCATE(dN)

!-----------------------------------------------------------------------
!  Sort ERPA2 Matrix and take the slice of significant values
!-----------------------------------------------------------------------
      M_ERPA = M_ERPA(M_SORT_IDX,M_SORT_IDX)
      DEALLOCATE(M_SORT_IDX)

      ALLOCATE(M_ERPA2(2*N_S+NBF5, 2*N_S+NBF5))
      M_ERPA2 = M_ERPA(1:2*N_S+NBF5,1:2*N_S+NBF5)
      DEALLOCATE(M_ERPA)
     
      DO I=1,2*N_S+NBF5
        M_ERPA2(I,:) = M_ERPA2(I,:)/VV(I)
      END DO
      DEALLOCATE(VV)

      WRITE(6,3)
      CALL DIAG_ERPA(M_ERPA2, 2*N_S+NBF5, 2)
      DEALLOCATE(M_ERPA2)
      
      RETURN
!-----------------------------------------------------------------------
    1 FORMAT(/,                                                       &
       1X,'---------------------------',/,                            &
       1X,' ERPA0 Excitation Energies' ,/,                            &
       1X,'---------------------------',/,                            &
       3X,' ID ',1X,'       w(eV)     ',/)
    2 FORMAT(/,                                                       &
       1X,'---------------------------',/,                            &
       1X,' ERPA1  Excitation Energies ',/,                            &
       1X,'---------------------------',/,                            &
       3X,' ID ',1X,'       w(eV)     ',/)
    3 FORMAT(/,                                                       &
       1X,'---------------------------',/,                            &
       1X,' ERPA2 Excitation Energies' ,/,                            &
       1X,'---------------------------',/,                            &
       3X,' ID ',1X,'       w(eV)     ',/)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END

! 2RDM
      SUBROUTINE Compute_2RDM(RO,Daa,Dab)
      COMMON/MAIN/NATOMS,ICH,MUL,NE,NA,NB,NSHELL,NPRIMI,NBF,NBFT,NSQ
      COMMON/INPNOF_ORBSPACE0/NO1,NDOC,NCO,NCWO,NVIR,NAC,NO0
      COMMON/INPNOF_ORBSPACE1/NSOC,NDNS,MSpin,HighSpin
      COMMON/INPFILE_NBF5/NBF5,NBFT5,NSQ5
      COMMON/INPNOF_PNOF/IPNOF,NTWOPAR
      DOUBLE PRECISION,DIMENSION(NBF5) :: RO
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5,NBF5,NBF5) :: Daa,Dab

      INTEGER :: I,J,P,Q,R,T
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: ROd,Rd
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: FI
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: Pi_s,Pi_d,INTRA,INTER,INTER2

      ALLOCATE (ROd(NBF5),Rd(NBF5))
      ALLOCATE(FI(NBF5))
      ALLOCATE(Pi_s(NBF5,NBF5),Pi_d(NBF5,NBF5),INTRA(NBF5,NBF5),INTER(NBF5,NBF5),INTER2(NBF5,NBF5))

      IF(IPNOF.NE.5 .AND. IPNOF.NE.7 .AND. IPNOF.NE.8) THEN
        WRITE(*,*) "ERROR in ERPA: Compute_2RDM Subroutine has been coded only for PNOF5, PNOF7 and GNOF"
      END IF

      Daa(1:NBF5,1:NBF5,1:NBF5,1:NBF5) = 0.0D0
      Dab(1:NBF5,1:NBF5,1:NBF5,1:NBF5) = 0.0D0

      ! Form D(alpha,alpha) and D(alpha,beta) for PNOF5

      INTRA(1:NBF5,1:NBF5) = 0.0D0
      INTER(1:NBF5,1:NBF5) = 0.0D0

      !PNOF5
      DO I=1,NBF5
        DO J=1,NBF5
          INTER=RO(I)*RO(J)
        END DO
      END DO

      ! PNOF5 eq. (18)  J. Chem. Phys. 139, 234109 (2013). 
      DO l=1,NDOC
       ln = NO1+l
       INTER(ln,ln) = 0.0D0
       INTRA(ln,ln) = SQRT(RO(ln))*SQRT(RO(ln))
       !DO i=NDNS+NCWO*(NDOC-l)+1,NDNS+NCWO*(NDOC-l+1)     !old-sort
       DO i=NDNS+NDOC-l+1,NDNS+NDOC-l+1+(NCWO-1)*NDOC,NDOC !new-sort
        in = NO1+i
        INTER(ln,in) = 0.0D0
        INTER(in,ln) = 0.0D0
        INTRA(ln,in) = -SQRT(RO(ln))*SQRT(RO(in))
        INTRA(in,ln) = -SQRT(RO(in))*SQRT(RO(ln))
        !DO j=NDNS+NCWO*(NDOC-l)+1,NDNS+NCWO*(NDOC-l+1)     !old-sort
        DO j=NDNS+NDOC-l+1,NDNS+NDOC-l+1+(NCWO-1)*NDOC,NDOC !new-sort
         jn = NO1+j
         INTER(jn,in) = 0.0D0
         INTRA(jn,in) = SQRT(RO(jn))*SQRT(RO(in))
        ENDDO
       ENDDO
      ENDDO

      DO I=NDOC+1,NDNS
        INTER(I,I) = 0.0D0
      END DO 

      ! PNOF5 eqs. (23) and (29) J. Chem. Phys. 139, 234109 (2013).
      DO P=1,NBF5
        DO Q=1,NBF5
          Daa(P,Q,P,Q) = Daa(P,Q,P,Q) + INTER(P,Q)
          Daa(P,Q,Q,P) = Daa(P,Q,Q,P) - INTER(P,Q)
          Dab(P,Q,P,Q) = Dab(P,Q,P,Q) + INTER(P,Q)
          Dab(P,P,Q,Q) = Dab(P,P,Q,Q) + INTRA(P,Q)
        END DO
      END DO

      ! Add PNOF7 contributions to D(alpha,beta)
      ! eq. (21) Phys. Rev. A 100, 032508 (2019).
      IF(IPNOF.EQ.7 .OR. IPNOF.EQ.8) THEN
        FI = RO*(1-RO)
        FI = SQRT(MAX(FI,0.0D0))
        DO I=1,NBF5
          DO J=1,NBF5
            Pi_s(I,J)=FI(I)*FI(J)
          END DO
        END DO

        DO l=1,NDOC
         ln = NO1+l
         Pi_s(ln,ln) = 0.0D0
         !DO j=NDNS+NCWO*(NDOC-l)+1,NDNS+NCWO*(NDOC-l+1)     !old-sort
         DO i=NDNS+NDOC-l+1,NDNS+NDOC-l+1+(NCWO-1)*NDOC,NDOC !new-sort
          in = NO1+i
          Pi_s(ln,in) = 0.0D0
          Pi_s(in,ln) = 0.0D0
          !DO j=NDNS+NCWO*(NDOC-l)+1,NDNS+NCWO*(NDOC-l+1)     !old-sort
          DO j=NDNS+NDOC-l+1,NDNS+NDOC-l+1+(NCWO-1)*NDOC,NDOC !new-sort
           jn = NO1+j
           Pi_s(jn,in) = 0.0D0
          ENDDO
         ENDDO
        ENDDO
        DO I=NDOC+1,NDNS
          Pi_s(I,I) = 0.0D0
        END DO 

        INTER(1:NBF5,1:NBF5) = 0.0
        INTER(NDOC+1:NDNS,NDOC+1:NDNS) = Pi_s(NDOC+1:NDNS,NDOC+1:NDNS)

        ! In case of GNOF, modify Pi_s for the D(alpha,beta) contributions
        ! as in eq. (4)  Phys. Rev. Lett. 127, 233001 (2021).
        IF(IPNOF.EQ.8) THEN
          Pi_s(1:NDOC,1:NDOC) = 0.0D0
          Pi_s(1:NDOC,NDOC+1:NDNS) = Pi_s(1:NDOC,NDOC+1:NDNS)*0.5
          Pi_s(NDOC+1:NDNS,1:NDOC) = Pi_s(NDOC+1:NDNS,1:NDOC)*0.5
        END IF

        DO P=1,NBF5
          DO Q=1,NBF5
            Dab(P,Q,Q,P) = Dab(P,Q,Q,P) - INTER(P,Q)
          END DO
        END DO

        ! GNOF Pi_d for the D(alpha,beta) contributions
        ! as in eq. (4)  Phys. Rev. Lett. 127, 233001 (2021).
        IF(IPNOF.EQ.8) THEN
          Hcut = 0.02d0*DSQRT(2.0d0)
          ROd  = 0.0D0
          DROd = 0.0D0
          Rd   = 0.0D0
          DRd  = 0.0D0
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !     ROd (NO1+1:NBF5)
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          DO i=1,NDOC
           in = NO1+i
           HOLEin = 1.0d0-RO(in)
           COC = HOLEin/Hcut
           ARG  = - COC**2
           DARG = - 2*COC/Hcut
           Fin = DEXP(ARG)                       ! Hd/Hole
    !      ROd(NO1+1:NB)
           ROd(in) = RO(in) * Fin
    !      ROd(NA+1:NBF5)
           IF(NCWO>1)THEN                        ! extended PNOF
            do iw=1,ncwo
             ! above Fermi level
             !im = na+ncwo*(ndoc-i)+iw                       !old-sort
             im = no1+ndoc+(na-nb)+(ndoc-i+1)+ndoc*(iw-1)    !new-sort
             ROd(im) = RO(im) * Fin              ! ROd = RO*Hd/Hole
            enddo
           ELSE                                  ! perfect-pairing
            icf = na+ndoc-i+1
            ROd(icf) = RO(icf) * Fin             ! ROd = RO*Hd/Hole
           ENDIF
          ENDDO
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !     Rd(1:NBF5) = SQRT(ROd)
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          DO j=1,NBF5
           Rd(j) = DSQRT(ROd(j))
          ENDDO

          DO I=1,NBF5
            DO J=1,NBF5
              INTER(I,J) = Rd(I)*Rd(J) - ROd(I)*ROd(J)
              INTER2(I,J) = Rd(I)*Rd(J) + ROd(I)*ROd(J)
            END DO
          END DO

          DO l=1,NDOC
           ln = NO1+l
           INTER(ln,ln) = 0.0D0
           INTER2(ln,ln) = 0.0D0
           !DO i=NDNS+NCWO*(NDOC-l)+1,NDNS+NCWO*(NDOC-l+1)     !old-sort
           DO i=NDNS+NDOC-l+1,NDNS+NDOC-l+1+(NCWO-1)*NDOC,NDOC !new-sort
            in = NO1+i
            INTER(ln,in) = 0.0D0
            INTER(in,ln) = 0.0D0
            INTER2(ln,in) = 0.0D0
            INTER2(in,ln) = 0.0D0
            !DO j=NDNS+NCWO*(NDOC-l)+1,NDNS+NCWO*(NDOC-l+1)     !old-sort
            DO j=NDNS+NDOC-l+1,NDNS+NDOC-l+1+(NCWO-1)*NDOC,NDOC !new-sort
             jn = NO1+j
             INTER(jn,in) = 0.0D0
             INTER2(jn,in) = 0.0D0
            ENDDO
           ENDDO
          ENDDO

          INTER(NDOC+1:,NDOC+1:) = 0.0D0
          INTER(1:NDNS,1:NDNS) = 0.0D0
          INTER2(1:NDNS,:) = 0.0D0
          INTER2(:,1:NDNS) = 0.0D0

          DO I=1,NBF5
            DO J=1,NBF5
              Pi_d(I,J) = INTER(I,J) + INTER2(I,J)
            END DO
          END DO
          DO P=1,NBF5
            DO Q=1,NBF5
              Dab(P,P,Q,Q) = Dab(P,P,Q,Q) - Pi_d(P,Q)
            END DO
          END DO

        END IF

      END IF

      DEALLOCATE (ROd,Rd)
      DEALLOCATE(FI)
      DEALLOCATE(Pi_s,Pi_d,INTRA,INTER,INTER2)

      END

! DIAG_ERPA
      SUBROUTINE DIAG_ERPA(M,N,ERPA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER :: N, ERPA
      DOUBLE PRECISION,DIMENSION(N,N) :: M

      INTEGER :: I,J, LWORK, INFO, NVALS
      DOUBLE PRECISION :: DTEMP, TOL
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: ALPHAI,ALPHAR,BETA
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: VALS,WORK
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: ID,VL,VR

      LWORK = 8*N
      ALLOCATE(WORK(LWORK))
      ALLOCATE(ALPHAR(N),ALPHAI(N),BETA(N))

      ! This IDentity matrix is to use DGGEV, we can remove it 
      ! by changing DGGEV to DGEEV
      ALLOCATE(ID(N,N))
      ID = 0.0D0
      DO I=1,N
        ID(I,I) = 1.0D0
      END DO

      CALL DGGEV('N','N',N,M,N,ID,N,ALPHAR,ALPHAI,BETA,VL,N,VR,N,WORK,  &
                 LWORK,INFO)
      DEALLOCATE(ID, WORK)

      IF (INFO.NE.0) WRITE(*,*) "WARNING!!!! DGGEV INFO = ", INFO

      ! Remove small values due to numeric instabilities (TOL 0.5eV)
      IF(ERPA.EQ.0 .OR. ERPA.EQ. 1) TOL = 3.37E-4
      IF(ERPA.EQ.2) TOL = 0.01837

      ! Check number of real and complexes eigenvalues
      NVALS = 0
      DO J=1,N
        IF((ALPHAR(J) .GT. 1E-7) .AND. (ABS(ALPHAI(J)) .LT. 1E-7)) THEN
          IF (ALPHAR(J)/BETA(J) > TOL) THEN
            NVALS = NVALS + 1
          END IF
        END IF
      END DO
      ALLOCATE(VALS(NVALS))

      I = 0
      DO J=1,N
        IF((ALPHAR(J) .GT. 1E-7) .AND. (ABS(ALPHAI(J)) .LT. 1E-7)) THEN
          IF (ALPHAR(J)/BETA(J) > TOL) THEN
            I = I + 1
            VALS(I) = ALPHAR(J)/BETA(J)
          END IF
        END IF
      END DO
      DEALLOCATE(ALPHAR, ALPHAI, BETA)

      ! Sort
      DO I=1, NVALS-1
        DO J=1, NVALS-I
          IF (VALS(J) > VALS(J+1)) THEN
            DTEMP = VALS(J)
            VALS(J) = VALS(J+1)
            VALS(J+1) = DTEMP
           END IF
        END DO
      END DO

      DO I=1,MIN(NVALS,10)
        IF(ERPA.EQ.0) WRITE(6,1) I,SQRT(VALS(I))*27.2114
        IF(ERPA.EQ.1) WRITE(6,1) I,SQRT(VALS(I))*27.2114
        IF(ERPA.EQ.2) WRITE(6,1) I,VALS(I)*27.2114
      END DO
      DEALLOCATE(VALS)

!----------------------------------------------------------------------!
    1 FORMAT(1X,I5,F14.3)
!----------------------------------------------------------------------!

      END

