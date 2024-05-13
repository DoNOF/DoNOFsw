!======================================================================!
!                                                                      !
!                   LAPACK SUBROUTINES USED IN DoNOF                   !
!                                                                      !
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !

! DLAMCH
      DOUBLE PRECISION FUNCTION dlamch( CMACH )
*
*  -- LAPACK auxiliary routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     April 2012
*
*     .. Scalar Arguments ..
      CHARACTER          cmach
*     ..
*     .. Parameters ..
      DOUBLE PRECISION   one, zero
      parameter( one = 1.0d+0, zero = 0.0d+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            first, lrnd
      INTEGER            beta, imax, imin, it
      DOUBLE PRECISION   base, emax, emin, eps, prec, rmach, rmax, rmin,
     $                   rnd, sfmin, small, t
*     ..
*     .. External Functions ..
      LOGICAL            lsame
      EXTERNAL           lsame
*     ..
*     .. External Subroutines ..
      EXTERNAL           dlamc2
*     ..
*     .. Save statement ..
      SAVE               first, eps, sfmin, base, t, rnd, emin, rmin,
     $                   emax, rmax, prec
*     ..
*     .. Data statements ..
      DATA               first / .true. /
*     ..
*     .. Executable Statements ..
*
      IF( first ) THEN
         CALL dlamc2( beta, it, lrnd, eps, imin, rmin, imax, rmax )
         base = beta
         t = it
         IF( lrnd ) THEN
            rnd = one
            eps = ( base**( 1-it ) ) / 2
         ELSE
            rnd = zero
            eps = base**( 1-it )
         END IF
         prec = eps*base
         emin = imin
         emax = imax
         sfmin = rmin
         small = one / rmax
         IF( small.GE.sfmin ) THEN
*
*           Use SMALL plus a bit, to avoid the possibility of rounding
*           causing overflow when computing  1/sfmin.
*
            sfmin = small*( one+eps )
         END IF
      END IF
*
      IF( lsame( cmach, 'E' ) ) THEN
         rmach = eps
      ELSE IF( lsame( cmach, 'S' ) ) THEN
         rmach = sfmin
      ELSE IF( lsame( cmach, 'B' ) ) THEN
         rmach = base
      ELSE IF( lsame( cmach, 'P' ) ) THEN
         rmach = prec
      ELSE IF( lsame( cmach, 'N' ) ) THEN
         rmach = t
      ELSE IF( lsame( cmach, 'R' ) ) THEN
         rmach = rnd
      ELSE IF( lsame( cmach, 'M' ) ) THEN
         rmach = emin
      ELSE IF( lsame( cmach, 'U' ) ) THEN
         rmach = rmin
      ELSE IF( lsame( cmach, 'L' ) ) THEN
         rmach = emax
      ELSE IF( lsame( cmach, 'O' ) ) THEN
         rmach = rmax
      END IF
*
      dlamch = rmach
      first  = .false.
      RETURN
*
*     End of DLAMCH
*
      END

! DLAMC1
      SUBROUTINE DLAMC1( BETA, T, RND, IEEE1 )
*
*  -- LAPACK auxiliary routine (version 3.7.0) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2010
*
*     .. Scalar Arguments ..
      LOGICAL            IEEE1, RND
      INTEGER            BETA, T
*     ..
* =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            FIRST, LIEEE1, LRND
      INTEGER            LBETA, LT
      DOUBLE PRECISION   A, B, C, F, ONE, QTR, SAVEC, T1, T2
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           dlamc3
*     ..
*     .. Save statement ..
      SAVE               first, lieee1, lbeta, lrnd, lt
*     ..
*     .. Data statements ..
      DATA               first / .true. /
*     ..
*     .. Executable Statements ..
*
      IF( first ) THEN
         one = 1
*
*        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA,
*        IEEE1, T and RND.
*
*        Throughout this routine  we use the function  DLAMC3  to ensure
*        that relevant values are  stored and not held in registers,  or
*        are not affected by optimizers.
*
*        Compute  a = 2.0**m  with the  smallest positive integer m such
*        that
*
*           fl( a + 1.0 ) = a.
*
         a = 1
         c = 1
*
*+       WHILE( C.EQ.ONE )LOOP
   10    CONTINUE
         IF( c.EQ.one ) THEN
            a = 2*a
            c = dlamc3( a, one )
            c = dlamc3( c, -a )
            GO TO 10
         END IF
*+       END WHILE
*
*        Now compute  b = 2.0**m  with the smallest positive integer m
*        such that
*
*           fl( a + b ) .gt. a.
*
         b = 1
         c = dlamc3( a, b )
*
*+       WHILE( C.EQ.A )LOOP
   20    CONTINUE
         IF( c.EQ.a ) THEN
            b = 2*b
            c = dlamc3( a, b )
            GO TO 20
         END IF
*+       END WHILE
*
*        Now compute the base.  a and c  are neighbouring floating point
*        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so
*        their difference is beta. Adding 0.25 to c is to ensure that it
*        is truncated to beta and not ( beta - 1 ).
*
         qtr = one / 4
         savec = c
         c = dlamc3( c, -a )
         lbeta = c + qtr
*
*        Now determine whether rounding or chopping occurs,  by adding a
*        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a.
*
         b = lbeta
         f = dlamc3( b / 2, -b / 100 )
         c = dlamc3( f, a )
         IF( c.EQ.a ) THEN
            lrnd = .true.
         ELSE
            lrnd = .false.
         END IF
         f = dlamc3( b / 2, b / 100 )
         c = dlamc3( f, a )
         IF( ( lrnd ) .AND. ( c.EQ.a ) )
     $      lrnd = .false.
*
*        Try and decide whether rounding is done in the  IEEE  'round to
*        nearest' style. B/2 is half a unit in the last place of the two
*        numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit
*        zero, and SAVEC is odd. Thus adding B/2 to A should not  change
*        A, but adding B/2 to SAVEC should change SAVEC.
*
         t1 = dlamc3( b / 2, a )
         t2 = dlamc3( b / 2, savec )
         lieee1 = ( t1.EQ.a ) .AND. ( t2.GT.savec ) .AND. lrnd
*
*        Now find  the  mantissa, t.  It should  be the  integer part of
*        log to the base beta of a,  however it is safer to determine  t
*        by powering.  So we find t as the smallest positive integer for
*        which
*
*           fl( beta**t + 1.0 ) = 1.0.
*
         lt = 0
         a = 1
         c = 1
*
*+       WHILE( C.EQ.ONE )LOOP
   30    CONTINUE
         IF( c.EQ.one ) THEN
            lt = lt + 1
            a = a*lbeta
            c = dlamc3( a, one )
            c = dlamc3( c, -a )
            GO TO 30
         END IF
*+       END WHILE
*
      END IF
*
      beta = lbeta
      t = lt
      rnd = lrnd
      ieee1 = lieee1
      first = .false.
      RETURN
*
*     End of DLAMC1
*
      END

! DLAMC2      
      SUBROUTINE DLAMC2( BETA, T, RND, EPS, EMIN, RMIN, EMAX, RMAX )
*
*  -- LAPACK auxiliary routine (version 3.7.0) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2010
*
*     .. Scalar Arguments ..
      LOGICAL            RND
      INTEGER            BETA, EMAX, EMIN, T
      DOUBLE PRECISION   EPS, RMAX, RMIN
*     ..
* =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            FIRST, IEEE, IWARN, LIEEE1, LRND
      INTEGER            GNMIN, GPMIN, I, LBETA, LEMAX, LEMIN, LT,
     $                   NGNMIN, NGPMIN
      DOUBLE PRECISION   A, B, C, HALF, LEPS, LRMAX, LRMIN, ONE, RBASE,
     $                   SIXTH, SMALL, THIRD, TWO, ZERO
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           dlamc3
*     ..
*     .. External Subroutines ..
      EXTERNAL           dlamc1, dlamc4, dlamc5
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, max, min
*     ..
*     .. Save statement ..
      SAVE               first, iwarn, lbeta, lemax, lemin, leps, lrmax,
     $                   lrmin, lt
*     ..
*     .. Data statements ..
      DATA               first / .true. / , iwarn / .false. /
*     ..
*     .. Executable Statements ..
*
      IF( first ) THEN
         zero = 0
         one = 1
         two = 2
*
*        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of
*        BETA, T, RND, EPS, EMIN and RMIN.
*
*        Throughout this routine  we use the function  DLAMC3  to ensure
*        that relevant values are stored  and not held in registers,  or
*        are not affected by optimizers.
*
*        DLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1.
*
         CALL dlamc1( lbeta, lt, lrnd, lieee1 )
*
*        Start to find EPS.
*
         b = lbeta
         a = b**( -lt )
         leps = a
*
*        Try some tricks to see whether or not this is the correct  EPS.
*
         b = two / 3
         half = one / 2
         sixth = dlamc3( b, -half )
         third = dlamc3( sixth, sixth )
         b = dlamc3( third, -half )
         b = dlamc3( b, sixth )
         b = abs( b )
         IF( b.LT.leps )
     $      b = leps
*
         leps = 1
*
*+       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP
   10    CONTINUE
         IF( ( leps.GT.b ) .AND. ( b.GT.zero ) ) THEN
            leps = b
            c = dlamc3( half*leps, ( two**5 )*( leps**2 ) )
            c = dlamc3( half, -c )
            b = dlamc3( half, c )
            c = dlamc3( half, -b )
            b = dlamc3( half, c )
            GO TO 10
         END IF
*+       END WHILE
*
         IF( a.LT.leps )
     $      leps = a
*
*        Computation of EPS complete.
*
*        Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3)).
*        Keep dividing  A by BETA until (gradual) underflow occurs. This
*        is detected when we cannot recover the previous A.
*
         rbase = one / lbeta
         small = one
         DO 20 i = 1, 3
            small = dlamc3( small*rbase, zero )
   20    CONTINUE
         a = dlamc3( one, small )
         CALL dlamc4( ngpmin, one, lbeta )
         CALL dlamc4( ngnmin, -one, lbeta )
         CALL dlamc4( gpmin, a, lbeta )
         CALL dlamc4( gnmin, -a, lbeta )
         ieee = .false.
*
         IF( ( ngpmin.EQ.ngnmin ) .AND. ( gpmin.EQ.gnmin ) ) THEN
            IF( ngpmin.EQ.gpmin ) THEN
               lemin = ngpmin
*            ( Non twos-complement machines, no gradual underflow;
*              e.g.,  VAX )
            ELSE IF( ( gpmin-ngpmin ).EQ.3 ) THEN
               lemin = ngpmin - 1 + lt
               ieee = .true.
*            ( Non twos-complement machines, with gradual underflow;
*              e.g., IEEE standard followers )
            ELSE
               lemin = min( ngpmin, gpmin )
*            ( A guess; no known machine )
               iwarn = .true.
            END IF
*
         ELSE IF( ( ngpmin.EQ.gpmin ) .AND. ( ngnmin.EQ.gnmin ) ) THEN
            IF( abs( ngpmin-ngnmin ).EQ.1 ) THEN
               lemin = max( ngpmin, ngnmin )
*            ( Twos-complement machines, no gradual underflow;
*              e.g., CYBER 205 )
            ELSE
               lemin = min( ngpmin, ngnmin )
*            ( A guess; no known machine )
               iwarn = .true.
            END IF
*
         ELSE IF( ( abs( ngpmin-ngnmin ).EQ.1 ) .AND.
     $            ( gpmin.EQ.gnmin ) ) THEN
            IF( ( gpmin-min( ngpmin, ngnmin ) ).EQ.3 ) THEN
               lemin = max( ngpmin, ngnmin ) - 1 + lt
*            ( Twos-complement machines with gradual underflow;
*              no known machine )
            ELSE
               lemin = min( ngpmin, ngnmin )
*            ( A guess; no known machine )
               iwarn = .true.
            END IF
*
         ELSE
            lemin = min( ngpmin, ngnmin, gpmin, gnmin )
*         ( A guess; no known machine )
            iwarn = .true.
         END IF
         first = .false.
***
* Comment out this if block if EMIN is ok
         IF( iwarn ) THEN
            first = .true.
            WRITE( 6, fmt = 9999 )lemin
         END IF
***
*
*        Assume IEEE arithmetic if we found denormalised  numbers above,
*        or if arithmetic seems to round in the  IEEE style,  determined
*        in routine DLAMC1. A true IEEE machine should have both  things
*        true; however, faulty machines may have one or the other.
*
         ieee = ieee .OR. lieee1
*
*        Compute  RMIN by successive division by  BETA. We could compute
*        RMIN as BASE**( EMIN - 1 ),  but some machines underflow during
*        this computation.
*
         lrmin = 1
         DO 30 i = 1, 1 - lemin
            lrmin = dlamc3( lrmin*rbase, zero )
   30    CONTINUE
*
*        Finally, call DLAMC5 to compute EMAX and RMAX.
*
         CALL dlamc5( lbeta, lt, lemin, ieee, lemax, lrmax )
      END IF
*
      beta = lbeta
      t = lt
      rnd = lrnd
      eps = leps
      emin = lemin
      rmin = lrmin
      emax = lemax
      rmax = lrmax
*
      RETURN
*
 9999 FORMAT( / / ' WARNING. The value EMIN may be incorrect:-',
     $      '  EMIN = ', i8, /
     $      ' If, after inspection, the value EMIN looks',
     $      ' acceptable please comment out ',
     $      / ' the IF block as marked within the code of routine',
     $      ' DLAMC2,', / ' otherwise supply EMIN explicitly.', / )
*
*     End of DLAMC2
*
      END

! DLAMC3
      DOUBLE PRECISION FUNCTION dlamc3( A, B )
*
*  -- LAPACK auxiliary routine (version 3.7.0) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2010
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   a, b
*     ..
* =====================================================================
*
*     .. Executable Statements ..
*
      dlamc3 = a + b
*
      RETURN
*
*     End of DLAMC3
*
      END

! DLAMC4
      SUBROUTINE DLAMC4( EMIN, START, BASE )
*
*  -- LAPACK auxiliary routine (version 3.7.0) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2010
*
*     .. Scalar Arguments ..
      INTEGER            BASE, EMIN
      DOUBLE PRECISION   START
*     ..
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION   A, B1, B2, C1, C2, D1, D2, ONE, RBASE, ZERO
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           dlamc3
*     ..
*     .. Executable Statements ..
*
      a = start
      one = 1
      rbase = one / base
      zero = 0
      emin = 1
      b1 = dlamc3( a*rbase, zero )
      c1 = a
      c2 = a
      d1 = a
      d2 = a
*+    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND.
*    $       ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP
   10 CONTINUE
      IF( ( c1.EQ.a ) .AND. ( c2.EQ.a ) .AND. ( d1.EQ.a ) .AND.
     $    ( d2.EQ.a ) ) THEN
         emin = emin - 1
         a = b1
         b1 = dlamc3( a / base, zero )
         c1 = dlamc3( b1*base, zero )
         d1 = zero
         DO 20 i = 1, base
            d1 = d1 + b1
   20    CONTINUE
         b2 = dlamc3( a*rbase, zero )
         c2 = dlamc3( b2 / rbase, zero )
         d2 = zero
         DO 30 i = 1, base
            d2 = d2 + b2
   30    CONTINUE
         GO TO 10
      END IF
*+    END WHILE
*
      RETURN
*
*     End of DLAMC4
*
      END

! DLAMC5
      SUBROUTINE DLAMC5( BETA, P, EMIN, IEEE, EMAX, RMAX )
*
*  -- LAPACK auxiliary routine (version 3.7.0) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2010
*
*     .. Scalar Arguments ..
      LOGICAL            IEEE
      INTEGER            BETA, EMAX, EMIN, P
      DOUBLE PRECISION   RMAX
*     ..
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      parameter( zero = 0.0d0, one = 1.0d0 )
*     ..
*     .. Local Scalars ..
      INTEGER            EXBITS, EXPSUM, I, LEXP, NBITS, TRY, UEXP
      DOUBLE PRECISION   OLDY, RECBAS, Y, Z
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           dlamc3
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          mod
*     ..
*     .. Executable Statements ..
*
*     First compute LEXP and UEXP, two powers of 2 that bound
*     abs(EMIN). We then assume that EMAX + abs(EMIN) will sum
*     approximately to the bound that is closest to abs(EMIN).
*     (EMAX is the exponent of the required number RMAX).
*
      lexp = 1
      exbits = 1
   10 CONTINUE
      try = lexp*2
      IF( try.LE.( -emin ) ) THEN
         lexp = try
         exbits = exbits + 1
         GO TO 10
      END IF
      IF( lexp.EQ.-emin ) THEN
         uexp = lexp
      ELSE
         uexp = try
         exbits = exbits + 1
      END IF
*
*     Now -LEXP is less than or equal to EMIN, and -UEXP is greater
*     than or equal to EMIN. EXBITS is the number of bits needed to
*     store the exponent.
*
      IF( ( uexp+emin ).GT.( -lexp-emin ) ) THEN
         expsum = 2*lexp
      ELSE
         expsum = 2*uexp
      END IF
*
*     EXPSUM is the exponent range, approximately equal to
*     EMAX - EMIN + 1 .
*
      emax = expsum + emin - 1
      nbits = 1 + exbits + p
*
*     NBITS is the total number of bits needed to store a
*     floating-point number.
*
      IF( ( mod( nbits, 2 ).EQ.1 ) .AND. ( beta.EQ.2 ) ) THEN
*
*        Either there are an odd number of bits used to store a
*        floating-point number, which is unlikely, or some bits are
*        not used in the representation of numbers, which is possible,
*        (e.g. Cray machines) or the mantissa has an implicit bit,
*        (e.g. IEEE machines, Dec Vax machines), which is perhaps the
*        most likely. We have to assume the last alternative.
*        If this is true, then we need to reduce EMAX by one because
*        there must be some way of representing zero in an implicit-bit
*        system. On machines like Cray, we are reducing EMAX by one
*        unnecessarily.
*
         emax = emax - 1
      END IF
*
      IF( ieee ) THEN
*
*        Assume we are on an IEEE machine which reserves one exponent
*        for infinity and NaN.
*
         emax = emax - 1
      END IF
*
*     Now create RMAX, the largest machine number, which should
*     be equal to (1.0 - BETA**(-P)) * BETA**EMAX .
*
*     First compute 1.0 - BETA**(-P), being careful that the
*     result is less than 1.0 .
*
      recbas = one / beta
      z = beta - one
      y = zero
      DO 20 i = 1, p
         z = z*recbas
         IF( y.LT.one )oldy = y
         y = dlamc3( y, z )
   20 CONTINUE
      IF( y.GE.one )y = oldy
*
*     Now multiply by BETA**EMAX to get RMAX.
*
      DO 30 i = 1, emax
         y = dlamc3( y*beta, zero )
   30 CONTINUE
*
      rmax = y
      RETURN
*
*     End of DLAMC5
*
      END

! DLAPY2
      DOUBLE PRECISION FUNCTION DLAPY2( X, Y )
C
C  -- LAPACK AUXILIARY ROUTINE (VERSION 1.1) --
!     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
!     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
!     OCTOBER 31, 1992
C
!     .. SCALAR ARGUMENTS ..
      DOUBLE PRECISION   X, Y
!     ..
C
C  PURPOSE
C  =======
C
C  DLAPY2 RETURNS SQRT(X**2+Y**2), TAKING CARE NOT TO CAUSE UNNECESSARY
C  OVERFLOW.
C
C  ARGUMENTS
C  =========
C
C  X       (INPUT) DOUBLE PRECISION
C  Y       (INPUT) DOUBLE PRECISION
!          X AND Y SPECIFY THE VALUES X AND Y.
C
C  =====================================================================
C
!     .. PARAMETERS ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. LOCAL SCALARS ..
      DOUBLE PRECISION   W, XABS, YABS, Z
!     ..
!     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          ABS, MAX, MIN, SQRT
!     ..
!     .. EXECUTABLE STATEMENTS ..
C
      XABS = ABS( X )
      YABS = ABS( Y )
      W = MAX( XABS, YABS )
      Z = MIN( XABS, YABS )
      IF( Z==ZERO ) THEN
         DLAPY2 = W
      ELSE
         DLAPY2 = W*SQRT( ONE+( Z / W )**2 )
      END IF
      RETURN
C
!     END OF DLAPY2
C
      END

! DLAPY3
      DOUBLE PRECISION FUNCTION DLAPY3( X, Y, Z )
C
C  -- LAPACK AUXILIARY ROUTINE (VERSION 1.1) --
!     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
!     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
!     OCTOBER 31, 1992
C
!     .. SCALAR ARGUMENTS ..
      DOUBLE PRECISION   X, Y, Z
!     ..
C
C  PURPOSE
C  =======
C
C  DLAPY3 RETURNS SQRT(X**2+Y**2+Z**2), TAKING CARE NOT TO CAUSE
C  UNNECESSARY OVERFLOW.
C
C  ARGUMENTS
C  =========
C
C  X       (INPUT) DOUBLE PRECISION
C  Y       (INPUT) DOUBLE PRECISION
C  Z       (INPUT) DOUBLE PRECISION
!          X, Y AND Z SPECIFY THE VALUES X, Y AND Z.
C
C  =====================================================================
C
!     .. PARAMETERS ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
!     ..
!     .. LOCAL SCALARS ..
      DOUBLE PRECISION   W, XABS, YABS, ZABS
!     ..
!     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          ABS, MAX, SQRT
!     ..
!     .. EXECUTABLE STATEMENTS ..
C
      XABS = ABS( X )
      YABS = ABS( Y )
      ZABS = ABS( Z )
      W = MAX( XABS, YABS, ZABS )
      IF( W==ZERO ) THEN
         DLAPY3 = ZERO
      ELSE
         DLAPY3 = W*SQRT( ( XABS / W )**2+( YABS / W )**2+
     $            ( ZABS / W )**2 )
      END IF
      RETURN
C
!     END OF DLAPY3
C
      END

! DLADIV
      SUBROUTINE DLADIV( A, B, C, D, P, Q )
!
!  -- LAPACK AUXILIARY ROUTINE (VERSION 1.1) --
!     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
!     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
!     OCTOBER 31, 1992
!
!     .. SCALAR ARGUMENTS ..
      DOUBLE PRECISION   A, B, C, D, P, Q
!     ..
!
!  PURPOSE
!  =======
!
!  DLADIV PERFORMS COMPLEX DIVISION IN  REAL ARITHMETIC
!
!                        A + I*B
!             P + I*Q = ---------
!                        C + I*D
!
!  THE ALGORITHM IS DUE TO ROBERT L. SMITH AND CAN BE FOUND
!  IN D. KNUTH, THE ART OF COMPUTER PROGRAMMING, VOL.2, P.195
!
!  ARGUMENTS
!  =========
!
!  A       (INPUT) DOUBLE PRECISION
!  B       (INPUT) DOUBLE PRECISION
!  C       (INPUT) DOUBLE PRECISION
!  D       (INPUT) DOUBLE PRECISION
!          THE SCALARS A, B, C, AND D IN THE ABOVE EXPRESSION.
!
!  P       (OUTPUT) DOUBLE PRECISION
!  Q       (OUTPUT) DOUBLE PRECISION
!          THE SCALARS P AND Q IN THE ABOVE EXPRESSION.
!
!  =====================================================================
!
!     .. LOCAL SCALARS ..
      DOUBLE PRECISION   E, F
!     ..
!     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          ABS
!     ..
!     .. EXECUTABLE STATEMENTS ..
!
      IF( ABS( D )<ABS( C ) ) THEN
         E = D / C
         F = C + D*E
         P = ( A+B*E ) / F
         Q = ( B-A*E ) / F
      ELSE
         E = C / D
         F = D + C*E
         P = ( B+A*E ) / F
         Q = ( -A+B*E ) / F
      END IF
!
      RETURN
!
!     END OF DLADIV
!
      END

! DLARTG
      SUBROUTINE DLARTG( F, G, CS, SN, R )
!
!  -- LAPACK AUXILIARY ROUTINE (VERSION 1.1) --
!     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
!     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
!     OCTOBER 31, 1992
!
!     .. SCALAR ARGUMENTS ..
      DOUBLE PRECISION   CS, F, G, R, SN
!     ..
!
!  PURPOSE
!  =======
!
!  DLARTG GENERATE A PLANE ROTATION SO THAT
!
!     [  CS  SN  ]  .  [ F ]  =  [ R ]   WHERE CS**2 + SN**2 = 1.
!     [ -SN  CS  ]     [ G ]     [ 0 ]
!
!  THIS IS A FASTER VERSION OF THE BLAS ROUTINE DROTG, EXCEPT FOR
!  THE FOLLOWING DIFFERENCES:
!     F AND G ARE UNCHANGED ON RETURN.
!     IF G=0, THEN CS=1 AND SN=0.
!     IF F=0 AND (G /= 0), THEN CS=0 AND SN=1 WITHOUT DOING ANY
!        FLOATING POINT OPERATIONS (SAVES WORK IN DBDSQR WHEN
!        THERE ARE ZEROS ON THE DIAGONAL).
!
!  ARGUMENTS
!  =========
!
!  F       (INPUT) DOUBLE PRECISION
!          THE FIRST COMPONENT OF VECTOR TO BE ROTATED.
!
!  G       (INPUT) DOUBLE PRECISION
!          THE SECOND COMPONENT OF VECTOR TO BE ROTATED.
!
!  CS      (OUTPUT) DOUBLE PRECISION
!          THE COSINE OF THE ROTATION.
!
!  SN      (OUTPUT) DOUBLE PRECISION
!          THE SINE OF THE ROTATION.
!
!  R       (OUTPUT) DOUBLE PRECISION
!          THE NONZERO COMPONENT OF THE ROTATED VECTOR.
!
!  =====================================================================
!
!     .. PARAMETERS ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. LOCAL SCALARS ..
      DOUBLE PRECISION   T, TT
!     ..
!     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          ABS, SQRT
!     ..
!     .. EXECUTABLE STATEMENTS ..
!
      IF( G==ZERO ) THEN
         CS = ONE
         SN = ZERO
         R = F
      ELSE IF( F==ZERO ) THEN
         CS = ZERO
         SN = ONE
         R = G
      ELSE
         IF( ABS( F )>ABS( G ) ) THEN
            T = G / F
            TT = SQRT( ONE+T*T )
            CS = ONE / TT
            SN = T*CS
            R = F*TT
         ELSE
            T = F / G
            TT = SQRT( ONE+T*T )
            SN = ONE / TT
            CS = T*SN
            R = G*TT
         END IF
      END IF
      RETURN
C
!     END OF DLARTG
C
      END

! DGESV
      SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
*      
* DGESV computes the solution to system of linear equations A * X = B
*
*  Definition:
*  ===========
*
*       SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
* 
*       .. Scalar Arguments ..
*       INTEGER            INFO, LDA, LDB, N, NRHS
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*       ..
*  
*
* Purpose:
* =============
*
* DGESV computes the solution to a real system of linear equations
*    A * X = B,
* where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
*
* The LU decomposition with partial pivoting and row interchanges is
* used to factor A as
*    A = P * L * U,
* where P is a permutation matrix, L is unit lower triangular, and U is
* upper triangular.  The factored form of A is then used to solve the
* system of equations A * X = B.
*
* Arguments:
* ==========
*
* \param[in] N
*          N is INTEGER
*          The number of linear equations, i.e., the order of the
*          matrix A.  N >= 0.
*
* \param[in] NRHS
*          NRHS is INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
* \param[in,out] A
*          A is DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the N-by-N coefficient matrix A.
*          On exit, the factors L and U from the factorization
*          A = P*L*U; the unit diagonal elements of L are not stored.
*
* \param[in] LDA
*          LDA is INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
* \param[out] IPIV
*          IPIV is INTEGER array, dimension (N)
*          The pivot indices that define the permutation matrix P;
*          row i of the matrix was interchanged with row IPIV(i).
*
* \param[in,out] B
*          B is DOUBLE PRECISION array, dimension (LDB,NRHS)
*          On entry, the N-by-NRHS matrix of right hand side matrix B.
*          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
*
* \param[in] LDB
*          LDB is INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
* \param[out] INFO
*          INFO is INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
*                has been completed, but the factor U is exactly
*                singular, so the solution could not be computed.
*
* Authors:
* ========
*
* \author Univ. of Tennessee 
* \author Univ. of California Berkeley 
* \author Univ. of Colorado Denver 
* \author NAG Ltd. 
*
* \date November 2011
*
* \ingroup doubleGEsolve
*
*  =====================================================================
*
*  -- LAPACK driver routine (version 3.4.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
*     November 2011
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  =====================================================================
*
*     .. External Subroutines ..
      EXTERNAL           DGETRF, DGETRS, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGESV ', -INFO )
         RETURN
      END IF
*
*     Compute the LU factorization of A.
*
      CALL DGETRF( N, N, A, LDA, IPIV, INFO )
      IF( INFO.EQ.0 ) THEN
*
*        Solve the system A*X = B, overwriting B with X.
*
         CALL DGETRS( 'No transpose', N, NRHS, A, LDA, IPIV, B, LDB,
     $                INFO )
      END IF
      RETURN
*
*     End of DGESV
*
      END

! DGETRF
      SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
*
*  Definition:
*  ===========
*
*       SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
* 
*       .. Scalar Arguments ..
*       INTEGER            INFO, LDA, M, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       DOUBLE PRECISION   A( LDA, * )
*       ..
*  
*  Purpose:
* ==========
*
* DGETRF computes an LU factorization of a general M-by-N matrix A
* using partial pivoting with row interchanges.
*
* The factorization has the form
*    A = P * L * U
* where P is a permutation matrix, L is lower triangular with unit
* diagonal elements (lower trapezoidal if m > n), and U is upper
* triangular (upper trapezoidal if m < n).
*
* This is the right-looking Level 3 BLAS version of the algorithm.
*
* Arguments:
* ==========
*
* \param[in] M
* 
*          M is INTEGER
*          The number of rows of the matrix A.  M >= 0.
* 
*
* \param[in] N
* 
*          N is INTEGER
*          The number of columns of the matrix A.  N >= 0.
* 
*
* \param[in,out] A
* 
*          A is DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the M-by-N matrix to be factored.
*          On exit, the factors L and U from the factorization
*          A = P*L*U; the unit diagonal elements of L are not stored.
* 
*
* \param[in] LDA
* 
*          LDA is INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
* 
*
* \param[out] IPIV
* 
*          IPIV is INTEGER array, dimension (min(M,N))
*          The pivot indices; for 1 <= i <= min(M,N), row i of the
*          matrix was interchanged with row IPIV(i).
* 
*
* \param[out] INFO
* 
*          INFO is INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
*                has been completed, but the factor U is exactly
*                singular, and division by zero will occur if it is used
*                to solve a system of equations.
* 
*
* -- LAPACK computational routine (version 3.4.0) --
* -- LAPACK is a software package provided by Univ. of Tennessee,    --
* -- Univ. of California Berkeley, Univ. of Colorado Denver
*    November 2011
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IINFO, J, JB, NB
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMM, DGETF2, DLASWP, DTRSM, XERBLA
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRF', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     Determine the block size for this environment.
*
      NB = ILAENV( 1, 'DGETRF', ' ', M, N, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN
*
*        Use unblocked code.
*
         CALL DGETF2( M, N, A, LDA, IPIV, INFO )
      ELSE
*
*        Use blocked code.
*
         DO 20 J = 1, MIN( M, N ), NB
            JB = MIN( MIN( M, N )-J+1, NB )
*
*           Factor diagonal and subdiagonal blocks and test for exact
*           singularity.
*
            CALL DGETF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
*
*           Adjust INFO and the pivot indices.
*
            IF( INFO.EQ.0 .AND. IINFO.GT.0 )
     $         INFO = IINFO + J - 1
            DO 10 I = J, MIN( M, J+JB-1 )
               IPIV( I ) = J - 1 + IPIV( I )
   10       CONTINUE
*
*           Apply interchanges to columns 1:J-1.
*
            CALL DLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
*
            IF( J+JB.LE.N ) THEN
*
*              Apply interchanges to columns J+JB:N.
*
               CALL DLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1,
     $                      IPIV, 1 )
*
*              Compute block row of U.
*
               CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB,
     $                     N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ),
     $                     LDA )
               IF( J+JB.LE.M ) THEN
*
*                 Update trailing submatrix.
*
                  CALL DGEMM( 'No transpose', 'No transpose', M-J-JB+1,
     $                        N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA,
     $                        A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ),
     $                        LDA )
               END IF
            END IF
   20    CONTINUE
      END IF
      RETURN
*
*     End of DGETRF
*
      END

! DGETRS
      SUBROUTINE DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
*      
*  Definition:
*  ===========
*
*       SUBROUTINE DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
* 
*       .. Scalar Arguments ..
*       CHARACTER          TRANS
*       INTEGER            INFO, LDA, LDB, N, NRHS
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*       ..
*  
*
*  Purpose:
*  ========
*
* DGETRS solves a system of linear equations
*    A * X = B  or  A**T * X = B
* with a general N-by-N matrix A using the LU factorization computed
* by DGETRF.
*
* Arguments:
* ==========
*
* \param[in] TRANS
* 
*          TRANS is CHARACTER*1
*          Specifies the form of the system of equations:
*          = 'N':  A * X = B  (No transpose)
*          = 'T':  A**T* X = B  (Transpose)
*          = 'C':  A**T* X = B  (Conjugate transpose = Transpose)
*
* \param[in] N
* 
*          N is INTEGER
*          The order of the matrix A.  N >= 0.
*
* \param[in] NRHS
* 
*          NRHS is INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
* \param[in] A
* 
*          A is DOUBLE PRECISION array, dimension (LDA,N)
*          The factors L and U from the factorization A = P*L*U
*          as computed by DGETRF.
*
* \param[in] LDA
* 
*          LDA is INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
* 
* \param[in] IPIV
* 
*          IPIV is INTEGER array, dimension (N)
*          The pivot indices from DGETRF; for 1<=i<=N, row i of the
*          matrix was interchanged with row IPIV(i).
*
* \param[in,out] B
* 
*          B is DOUBLE PRECISION array, dimension (LDB,NRHS)
*          On entry, the right hand side matrix B.
*          On exit, the solution matrix X.
*
* \param[in] LDB
* 
*          LDB is INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
* 
* \param[out] INFO
* 
*          INFO is INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*  -- LAPACK computational routine (version 3.4.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
*     November 2011
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOTRAN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLASWP, DTRSM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.
     $    LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 )
     $   RETURN
*
      IF( NOTRAN ) THEN
*
*        Solve A * X = B.
*
*        Apply row interchanges to the right hand sides.
*
         CALL DLASWP( NRHS, B, LDB, 1, N, IPIV, 1 )
*
*        Solve L*X = B, overwriting B with X.
*
         CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', N, NRHS,
     $               ONE, A, LDA, B, LDB )
*
*        Solve U*X = B, overwriting B with X.
*
         CALL DTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N,
     $               NRHS, ONE, A, LDA, B, LDB )
      ELSE
*
*        Solve A**T * X = B.
*
*        Solve U**T *X = B, overwriting B with X.
*
         CALL DTRSM( 'Left', 'Upper', 'Transpose', 'Non-unit', N, NRHS,
     $               ONE, A, LDA, B, LDB )
*
*        Solve L**T *X = B, overwriting B with X.
*
         CALL DTRSM( 'Left', 'Lower', 'Transpose', 'Unit', N, NRHS, ONE,
     $               A, LDA, B, LDB )
*
*        Apply row interchanges to the solution vectors.
*
         CALL DLASWP( NRHS, B, LDB, 1, N, IPIV, -1 )
      END IF
*
      RETURN
*
*     End of DGETRS
*
      END

! DGETF2
      SUBROUTINE DGETF2( M, N, A, LDA, IPIV, INFO )
*
*  Definition:
*  ===========
*
*       SUBROUTINE DGETF2( M, N, A, LDA, IPIV, INFO )
* 
*       .. Scalar Arguments ..
*       INTEGER            INFO, LDA, M, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       DOUBLE PRECISION   A( LDA, * )
*       ..
*
* Purpose:
* ========
*
* DGETF2 computes an LU factorization of a general m-by-n matrix A
* using partial pivoting with row interchanges.
*
* The factorization has the form
*    A = P * L * U
* where P is a permutation matrix, L is lower triangular with unit
* diagonal elements (lower trapezoidal if m > n), and U is upper
* triangular (upper trapezoidal if m < n).
*
* This is the right-looking Level 2 BLAS version of the algorithm.
*
* Arguments:
* ==========
*
* \param[in] M
* 
*          M is INTEGER
*          The number of rows of the matrix A.  M >= 0.
* 
*
* \param[in] N
* 
*          N is INTEGER
*          The number of columns of the matrix A.  N >= 0.
* 
*
* \param[in,out] A
* 
*          A is DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the m by n matrix to be factored.
*          On exit, the factors L and U from the factorization
*          A = P*L*U; the unit diagonal elements of L are not stored.
* 
*
* \param[in] LDA
* 
*          LDA is INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
* 
*
* \param[out] IPIV
* 
*          IPIV is INTEGER array, dimension (min(M,N))
*          The pivot indices; for 1 <= i <= min(M,N), row i of the
*          matrix was interchanged with row IPIV(i).
* 
*
* \param[out] INFO
* 
*          INFO is INTEGER
*          = 0: successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
*               has been completed, but the factor U is exactly
*               singular, and division by zero will occur if it is used
*               to solve a system of equations.
*
*  =====================================================================
*
*  -- LAPACK computational routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   SFMIN 
      INTEGER            I, J, JP
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH      
      INTEGER            IDAMAX
      EXTERNAL           DLAMCH, IDAMAX
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGER, DSCAL, DSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETF2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     Compute machine safe minimum 
* 
      SFMIN = DLAMCH('S')  
*
      DO 10 J = 1, MIN( M, N )
*
*        Find pivot and test for singularity.
*
         JP = J - 1 + IDAMAX( M-J+1, A( J, J ), 1 )
         IPIV( J ) = JP
         IF( A( JP, J ).NE.ZERO ) THEN
*
*           Apply the interchange to columns 1:N.
*
            IF( JP.NE.J )
     $         CALL DSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )
*
*           Compute elements J+1:M of J-th column.
*
            IF( J.LT.M ) THEN 
               IF( ABS(A( J, J )) .GE. SFMIN ) THEN 
                  CALL DSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 ) 
               ELSE 
                 DO 20 I = 1, M-J 
                    A( J+I, J ) = A( J+I, J ) / A( J, J ) 
   20            CONTINUE 
               END IF 
            END IF 
*
         ELSE IF( INFO.EQ.0 ) THEN
*
            INFO = J
         END IF
*
         IF( J.LT.MIN( M, N ) ) THEN
*
*           Update trailing submatrix.
*
            CALL DGER( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ), LDA,
     $                 A( J+1, J+1 ), LDA )
         END IF
   10 CONTINUE
      RETURN
*
*     End of DGETF2
*
      END

! DLASWP
      SUBROUTINE DLASWP( N, A, LDA, K1, K2, IPIV, INCX )
*
*  Definition:
*  ===========
*
*       SUBROUTINE DLASWP( N, A, LDA, K1, K2, IPIV, INCX )
* 
*       .. Scalar Arguments ..
*       INTEGER            INCX, K1, K2, LDA, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       DOUBLE PRECISION   A( LDA, * )
*       ..
*
* \par Purpose:
* =============
*
* DLASWP performs a series of row interchanges on the matrix A.
* One row interchange is initiated for each of rows K1 through K2 of A.
* 
*
* Arguments:
* ==========
*
* \param[in] N
* 
*          N is INTEGER
*          The number of columns of the matrix A.
* 
*
* \param[in,out] A
* 
*          A is DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the matrix of column dimension N to which the row
*          interchanges will be applied.
*          On exit, the permuted matrix.
* 
*
* \param[in] LDA
* 
*          LDA is INTEGER
*          The leading dimension of the array A.
* 
*
* \param[in] K1
* 
*          K1 is INTEGER
*          The first element of IPIV for which a row interchange will
*          be done.
* 
*
* \param[in] K2
* 
*          K2 is INTEGER
*          The last element of IPIV for which a row interchange will
*          be done.
* 
*
* \param[in] IPIV
* 
*          IPIV is INTEGER array, dimension (K2*abs(INCX))
*          The vector of pivot indices.  Only the elements in positions
*          K1 through K2 of IPIV are accessed.
*          IPIV(K) = L implies rows K and L are to be interchanged.
* 
*
* \param[in] INCX
* 
*          INCX is INTEGER
*          The increment between successive values of IPIV.  If IPIV
*          is negative, the pivots are applied in reverse order.
* 
*
* Authors:
* ========
*
* \author Univ. of Tennessee 
* \author Univ. of California Berkeley 
* \author Univ. of Colorado Denver 
* \author NAG Ltd. 
*
* \date September 2012
*
* \ingroup doubleOTHERauxiliary
*
* \par Further Details:
* =====================
*
* 
*
*  Modified by
*   R. C. Whaley, Computer Science Dept., Univ. of Tenn., Knoxville, USA
* 
*
*  =====================================================================
*
*  -- LAPACK auxiliary routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      INTEGER            INCX, K1, K2, LDA, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, I1, I2, INC, IP, IX, IX0, J, K, N32
      DOUBLE PRECISION   TEMP
*     ..
*     .. Executable Statements ..
*
*     Interchange row I with row IPIV(I) for each of rows K1 through K2.
*
      IF( INCX.GT.0 ) THEN
         IX0 = K1
         I1 = K1
         I2 = K2
         INC = 1
      ELSE IF( INCX.LT.0 ) THEN
         IX0 = 1 + ( 1-K2 )*INCX
         I1 = K2
         I2 = K1
         INC = -1
      ELSE
         RETURN
      END IF
*
      N32 = ( N / 32 )*32
      IF( N32.NE.0 ) THEN
         DO 30 J = 1, N32, 32
            IX = IX0
            DO 20 I = I1, I2, INC
               IP = IPIV( IX )
               IF( IP.NE.I ) THEN
                  DO 10 K = J, J + 31
                     TEMP = A( I, K )
                     A( I, K ) = A( IP, K )
                     A( IP, K ) = TEMP
   10             CONTINUE
               END IF
               IX = IX + INCX
   20       CONTINUE
   30    CONTINUE
      END IF
      IF( N32.NE.N ) THEN
         N32 = N32 + 1
         IX = IX0
         DO 50 I = I1, I2, INC
            IP = IPIV( IX )
            IF( IP.NE.I ) THEN
               DO 40 K = N32, N
                  TEMP = A( I, K )
                  A( I, K ) = A( IP, K )
                  A( IP, K ) = TEMP
   40          CONTINUE
            END IF
            IX = IX + INCX
   50    CONTINUE
      END IF
*
      RETURN
*
*     End of DLASWP
*
      END

! DTRSM
      SUBROUTINE DTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
*
*  Definition:
*  ===========
*
*       SUBROUTINE DTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
* 
*       .. Scalar Arguments ..
*       DOUBLE PRECISION ALPHA
*       INTEGER LDA,LDB,M,N
*       CHARACTER DIAG,SIDE,TRANSA,UPLO
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION A(LDA,*),B(LDB,*)
*       ..
*
* Purpose:
* ========
*
* DTRSM  solves one of the matrix equations
*
*    op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
*
* where alpha is a scalar, X and B are m by n matrices, A is a unit, or
* non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
*
*    op( A ) = A   or   op( A ) = A**T.
*
* The matrix X is overwritten on B.
*
* Arguments:
* ==========
*
* \param[in] SIDE
* 
*          SIDE is CHARACTER*1
*           On entry, SIDE specifies whether op( A ) appears on the left
*           or right of X as follows:
*
*              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
*
*              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
* 
*
* \param[in] UPLO
* 
*          UPLO is CHARACTER*1
*           On entry, UPLO specifies whether the matrix A is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
* 
*
* \param[in] TRANSA
* 
*          TRANSA is CHARACTER*1
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n'   op( A ) = A.
*
*              TRANSA = 'T' or 't'   op( A ) = A**T.
*
*              TRANSA = 'C' or 'c'   op( A ) = A**T.
* 
*
* \param[in] DIAG
* 
*          DIAG is CHARACTER*1
*           On entry, DIAG specifies whether or not A is unit triangular
*           as follows:
*
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*
* \param[in] M
* 
*          M is INTEGER
*           On entry, M specifies the number of rows of B. M must be at
*           least zero.
* 
*
* \param[in] N
* 
*          N is INTEGER
*           On entry, N specifies the number of columns of B.  N must be
*           at least zero.
* 
*
* \param[in] ALPHA
* 
*          ALPHA is DOUBLE PRECISION.
*           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
*           zero then  A is not referenced and  B need not be set before
*           entry.
* 
*
* \param[in] A
* 
*          A is DOUBLE PRECISION array of DIMENSION ( LDA, k ),
*           where k is m when SIDE = 'L' or 'l'  
*             and k is n when SIDE = 'R' or 'r'.
*           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
*           upper triangular part of the array  A must contain the upper
*           triangular matrix  and the strictly lower triangular part of
*           A is not referenced.
*           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
*           lower triangular part of the array  A must contain the lower
*           triangular matrix  and the strictly upper triangular part of
*           A is not referenced.
*           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
*           A  are not referenced either,  but are assumed to be  unity.
* 
*
* \param[in] LDA
* 
*          LDA is INTEGER
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
*           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
*           then LDA must be at least max( 1, n ).
* 
*
* \param[in,out] B
* 
*          B is DOUBLE PRECISION array of DIMENSION ( LDB, n ).
*           Before entry,  the leading  m by n part of the array  B must
*           contain  the  right-hand  side  matrix  B,  and  on exit  is
*           overwritten by the solution matrix  X.
* 
*
* \param[in] LDB
* 
*          LDB is INTEGER
*           On entry, LDB specifies the first dimension of B as declared
*           in  the  calling  (sub)  program.   LDB  must  be  at  least
*           max( 1, m ).
*
*  Level 3 Blas routine.
*
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Iain Duff, Jeremy Du Croz, Sven Hammarling
*
*  =====================================================================
*
*  -- Reference BLAS level3 routine (version 3.4.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA
      INTEGER LDA,LDB,M,N
      CHARACTER DIAG,SIDE,TRANSA,UPLO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),B(LDB,*)
*     ..
*
*  =====================================================================
*
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MAX
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,J,K,NROWA
      LOGICAL LSIDE,NOUNIT,UPPER
*     ..
*     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
*     ..
*
*     Test the input parameters.
*
      LSIDE = LSAME(SIDE,'L')
      IF (LSIDE) THEN
          NROWA = M
      ELSE
          NROWA = N
      END IF
      NOUNIT = LSAME(DIAG,'N')
      UPPER = LSAME(UPLO,'U')
*
      INFO = 0
      IF ((.NOT.LSIDE) .AND. (.NOT.LSAME(SIDE,'R'))) THEN
          INFO = 1
      ELSE IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
          INFO = 2
      ELSE IF ((.NOT.LSAME(TRANSA,'N')) .AND.
     +         (.NOT.LSAME(TRANSA,'T')) .AND.
     +         (.NOT.LSAME(TRANSA,'C'))) THEN
          INFO = 3
      ELSE IF ((.NOT.LSAME(DIAG,'U')) .AND. (.NOT.LSAME(DIAG,'N'))) THEN
          INFO = 4
      ELSE IF (M.LT.0) THEN
          INFO = 5
      ELSE IF (N.LT.0) THEN
          INFO = 6
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
          INFO = 9
      ELSE IF (LDB.LT.MAX(1,M)) THEN
          INFO = 11
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DTRSM ',INFO)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF (M.EQ.0 .OR. N.EQ.0) RETURN
*
*     And when  alpha.eq.zero.
*
      IF (ALPHA.EQ.ZERO) THEN
          DO 20 J = 1,N
              DO 10 I = 1,M
                  B(I,J) = ZERO
   10         CONTINUE
   20     CONTINUE
          RETURN
      END IF
*
*     Start the operations.
*
      IF (LSIDE) THEN
          IF (LSAME(TRANSA,'N')) THEN
*
*           Form  B := alpha*inv( A )*B.
*
              IF (UPPER) THEN
                  DO 60 J = 1,N
                      IF (ALPHA.NE.ONE) THEN
                          DO 30 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
   30                     CONTINUE
                      END IF
                      DO 50 K = M,1,-1
                          IF (B(K,J).NE.ZERO) THEN
                              IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                              DO 40 I = 1,K - 1
                                  B(I,J) = B(I,J) - B(K,J)*A(I,K)
   40                         CONTINUE
                          END IF
   50                 CONTINUE
   60             CONTINUE
              ELSE
                  DO 100 J = 1,N
                      IF (ALPHA.NE.ONE) THEN
                          DO 70 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
   70                     CONTINUE
                      END IF
                      DO 90 K = 1,M
                          IF (B(K,J).NE.ZERO) THEN
                              IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                              DO 80 I = K + 1,M
                                  B(I,J) = B(I,J) - B(K,J)*A(I,K)
   80                         CONTINUE
                          END IF
   90                 CONTINUE
  100             CONTINUE
              END IF
          ELSE
*
*           Form  B := alpha*inv( A**T )*B.
*
              IF (UPPER) THEN
                  DO 130 J = 1,N
                      DO 120 I = 1,M
                          TEMP = ALPHA*B(I,J)
                          DO 110 K = 1,I - 1
                              TEMP = TEMP - A(K,I)*B(K,J)
  110                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/A(I,I)
                          B(I,J) = TEMP
  120                 CONTINUE
  130             CONTINUE
              ELSE
                  DO 160 J = 1,N
                      DO 150 I = M,1,-1
                          TEMP = ALPHA*B(I,J)
                          DO 140 K = I + 1,M
                              TEMP = TEMP - A(K,I)*B(K,J)
  140                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/A(I,I)
                          B(I,J) = TEMP
  150                 CONTINUE
  160             CONTINUE
              END IF
          END IF
      ELSE
          IF (LSAME(TRANSA,'N')) THEN
*
*           Form  B := alpha*B*inv( A ).
*
              IF (UPPER) THEN
                  DO 210 J = 1,N
                      IF (ALPHA.NE.ONE) THEN
                          DO 170 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
  170                     CONTINUE
                      END IF
                      DO 190 K = 1,J - 1
                          IF (A(K,J).NE.ZERO) THEN
                              DO 180 I = 1,M
                                  B(I,J) = B(I,J) - A(K,J)*B(I,K)
  180                         CONTINUE
                          END IF
  190                 CONTINUE
                      IF (NOUNIT) THEN
                          TEMP = ONE/A(J,J)
                          DO 200 I = 1,M
                              B(I,J) = TEMP*B(I,J)
  200                     CONTINUE
                      END IF
  210             CONTINUE
              ELSE
                  DO 260 J = N,1,-1
                      IF (ALPHA.NE.ONE) THEN
                          DO 220 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
  220                     CONTINUE
                      END IF
                      DO 240 K = J + 1,N
                          IF (A(K,J).NE.ZERO) THEN
                              DO 230 I = 1,M
                                  B(I,J) = B(I,J) - A(K,J)*B(I,K)
  230                         CONTINUE
                          END IF
  240                 CONTINUE
                      IF (NOUNIT) THEN
                          TEMP = ONE/A(J,J)
                          DO 250 I = 1,M
                              B(I,J) = TEMP*B(I,J)
  250                     CONTINUE
                      END IF
  260             CONTINUE
              END IF
          ELSE
*
*           Form  B := alpha*B*inv( A**T ).
*
              IF (UPPER) THEN
                  DO 310 K = N,1,-1
                      IF (NOUNIT) THEN
                          TEMP = ONE/A(K,K)
                          DO 270 I = 1,M
                              B(I,K) = TEMP*B(I,K)
  270                     CONTINUE
                      END IF
                      DO 290 J = 1,K - 1
                          IF (A(J,K).NE.ZERO) THEN
                              TEMP = A(J,K)
                              DO 280 I = 1,M
                                  B(I,J) = B(I,J) - TEMP*B(I,K)
  280                         CONTINUE
                          END IF
  290                 CONTINUE
                      IF (ALPHA.NE.ONE) THEN
                          DO 300 I = 1,M
                              B(I,K) = ALPHA*B(I,K)
  300                     CONTINUE
                      END IF
  310             CONTINUE
              ELSE
                  DO 360 K = 1,N
                      IF (NOUNIT) THEN
                          TEMP = ONE/A(K,K)
                          DO 320 I = 1,M
                              B(I,K) = TEMP*B(I,K)
  320                     CONTINUE
                      END IF
                      DO 340 J = K + 1,N
                          IF (A(J,K).NE.ZERO) THEN
                              TEMP = A(J,K)
                              DO 330 I = 1,M
                                  B(I,J) = B(I,J) - TEMP*B(I,K)
  330                         CONTINUE
                          END IF
  340                 CONTINUE
                      IF (ALPHA.NE.ONE) THEN
                          DO 350 I = 1,M
                              B(I,K) = ALPHA*B(I,K)
  350                     CONTINUE
                      END IF
  360             CONTINUE
              END IF
          END IF
      END IF
*
      RETURN
*
*     End of DTRSM .
*
      END

! LSAME
      LOGICAL FUNCTION LSAME( CA, CB )
C
C  -- LAPACK AUXILIARY ROUTINE (VERSION 1.1) --
!     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
!     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
!     FEBRUARY 29, 1992 
C
!     .. SCALAR ARGUMENTS ..
      CHARACTER          CA, CB
!     ..
!
!  PURPOSE
!  =======
!
!  LSAME RETURNS .TRUE. IF CA IS THE SAME LETTER AS CB REGARDLESS OF
!  CASE.
!
!  ARGUMENTS
!  =========
!
!  CA      (INPUT) CHARACTER*1
!  CB      (INPUT) CHARACTER*1
!          CA AND CB SPECIFY THE SINGLE CHARACTERS TO BE COMPARED.
!
!     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          ICHAR
!     ..
!     .. LOCAL SCALARS ..
      INTEGER            INTA, INTB, ZCODE
!     ..
!     .. EXECUTABLE STATEMENTS ..
!
!     TEST IF THE CHARACTERS ARE EQUAL
!
      LSAME = CA==CB
      IF( LSAME )
     $   RETURN
!
!     NOW TEST FOR EQUIVALENCE IF BOTH CHARACTERS ARE ALPHABETIC.
!
      ZCODE = ICHAR( 'Z' )
!
!     USE 'Z' RATHER THAN 'A' SO THAT ASCII CAN BE DETECTED ON PRIME
!     MACHINES, ON WHICH ICHAR RETURNS A VALUE WITH BIT 8 SET.
!     ICHAR('A') ON PRIME MACHINES RETURNS 193 WHICH IS THE SAME AS
!     ICHAR('A') ON AN EBCDIC MACHINE.
!
      INTA = ICHAR( CA )
      INTB = ICHAR( CB )
!
      IF( ZCODE==90 .OR. ZCODE==122 ) THEN
!
!        ASCII IS ASSUMED - ZCODE IS THE ASCII CODE OF EITHER LOWER OR
!        UPPER CASE 'Z'.
!
         IF( INTA.GE.97 .AND. INTA<=122 ) INTA = INTA - 32
         IF( INTB.GE.97 .AND. INTB<=122 ) INTB = INTB - 32
!
      ELSE IF( ZCODE==233 .OR. ZCODE==169 ) THEN
!
!        EBCDIC IS ASSUMED - ZCODE IS THE EBCDIC CODE OF EITHER LOWER OR
!        UPPER CASE 'Z'.
!
         IF( INTA.GE.129 .AND. INTA<=137 .OR.
     $       INTA.GE.145 .AND. INTA<=153 .OR.
     $       INTA.GE.162 .AND. INTA<=169 ) INTA = INTA + 64
         IF( INTB.GE.129 .AND. INTB<=137 .OR.
     $       INTB.GE.145 .AND. INTB<=153 .OR.
     $       INTB.GE.162 .AND. INTB<=169 ) INTB = INTB + 64
!
      ELSE IF( ZCODE==218 .OR. ZCODE==250 ) THEN
!
!        ASCII IS ASSUMED, ON PRIME MACHINES - ZCODE IS THE ASCII CODE
!        PLUS 128 OF EITHER LOWER OR UPPER CASE 'Z'.
!
         IF( INTA.GE.225 .AND. INTA<=250 ) INTA = INTA - 32
         IF( INTB.GE.225 .AND. INTB<=250 ) INTB = INTB - 32
      END IF
      LSAME = INTA==INTB
!
!     RETURN
!
!     END OF LSAME
!
      END

! XERBLA
      SUBROUTINE XERBLA( SRNAME, INFO )
C
C  -- LAPACK AUXILIARY ROUTINE (VERSION 1.1) --
!     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
!     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
!     FEBRUARY 29, 1992
C
!     .. SCALAR ARGUMENTS ..
      CHARACTER*6        SRNAME
      INTEGER            INFO
!     ..
C
C  PURPOSE
C  =======
C
C  XERBLA  IS AN ERROR HANDLER FOR THE LAPACK ROUTINES.
C  IT IS CALLED BY AN LAPACK ROUTINE IF AN INPUT PARAMETER HAS AN
C  INVALID VALUE.  A MESSAGE IS PRINTED AND EXECUTION STOPS.
C
C  INSTALLERS MAY CONSIDER MODIFYING THE STOP STATEMENT IN ORDER TO
C  CALL SYSTEM-SPECIFIC EXCEPTION-HANDLING FACILITIES.
C
C  ARGUMENTS
C  =========
C
C  SRNAME  (INPUT) CHARACTER*6
!          THE NAME OF THE ROUTINE WHICH CALLED XERBLA.
C
C  INFO    (INPUT) INTEGER
!          THE POSITION OF THE INVALID PARAMETER IN THE PARAMETER LIST
!          OF THE CALLING ROUTINE.
C
!     .. EXECUTABLE STATEMENTS ..
C
      WRITE( *, FMT = 9999) SRNAME,INFO
      CALL ABRT
      STOP
C
 9999 FORMAT( ' ** ON ENTRY TO ', A6, ' PARAMETER NUMBER ', I2, ' HAD ',
     $      'AN ILLEGAL VALUE' )
C
!     END OF XERBLA
C
      END

! ILAENV
      INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
C
C  -- LAPACK AUXILIARY ROUTINE (PRELIMINARY VERSION) --
!     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
!     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
!     FEBRUARY 20, 1992
C
!     .. SCALAR ARGUMENTS ..
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ISPEC, N1, N2, N3, N4
!     ..
C
C  PURPOSE
C  =======
C
C  ILAENV IS CALLED FROM THE LAPACK ROUTINES TO CHOOSE PROBLEM-DEPENDENT
C  PARAMETERS FOR THE LOCAL ENVIRONMENT.  SEE ISPEC FOR A DESCRIPTION OF
C  THE PARAMETERS.
C
C  THIS VERSION PROVIDES A SET OF PARAMETERS WHICH SHOULD GIVE GOOD,
C  BUT NOT OPTIMAL, PERFORMANCE ON MANY OF THE CURRENTLY AVAILABLE
C  COMPUTERS.  USERS ARE ENCOURAGED TO MODIFY THIS SUBROUTINE TO SET
C  THE TUNING PARAMETERS FOR THEIR PARTICULAR MACHINE USING THE OPTION
C  AND PROBLEM SIZE INFORMATION IN THE ARGUMENTS.
C
C  THIS ROUTINE WILL NOT FUNCTION CORRECTLY IF IT IS CONVERTED TO ALL
C  LOWER CASE.  CONVERTING IT TO ALL UPPER CASE IS ALLOWED.
C
C  ARGUMENTS
C  =========
C
C  ISPEC   (INPUT) INTEGER
!          SPECIFIES THE PARAMETER TO BE RETURNED AS THE VALUE OF
!          ILAENV.
!          = 1: THE OPTIMAL BLOCKSIZE; IF THIS VALUE IS 1, AN UNBLOCKED
!               ALGORITHM WILL GIVE THE BEST PERFORMANCE.
!          = 2: THE MINIMUM BLOCK SIZE FOR WHICH THE BLOCK ROUTINE
!               SHOULD BE USED; IF THE USABLE BLOCK SIZE IS LESS THAN
!               THIS VALUE, AN UNBLOCKED ROUTINE SHOULD BE USED.
!          = 3: THE CROSSOVER POINT (IN A BLOCK ROUTINE, FOR N LESS
!               THAN THIS VALUE, AN UNBLOCKED ROUTINE SHOULD BE USED)
!          = 4: THE NUMBER OF SHIFTS, USED IN THE NONSYMMETRIC
!               EIGENVALUE ROUTINES
!          = 5: THE MINIMUM COLUMN DIMENSION FOR BLOCKING TO BE USED;
!               RECTANGULAR BLOCKS MUST HAVE DIMENSION AT LEAST K BY M,
!               WHERE K IS GIVEN BY ILAENV(2,...) AND M BY ILAENV(5,...)
!          = 6: THE CROSSOVER POINT FOR THE SVD (WHEN REDUCING AN M BY N
!               MATRIX TO BIDIAGONAL FORM, IF MAX(M,N)/MIN(M,N) EXCEEDS
!               THIS VALUE, A QR FACTORIZATION IS USED FIRST TO REDUCE
!               THE MATRIX TO A TRIANGULAR FORM.)
!          = 7: THE NUMBER OF PROCESSORS
!          = 8: THE CROSSOVER POINT FOR THE MULTISHIFT QR AND QZ METHODS
!               FOR NONSYMMETRIC EIGENVALUE PROBLEMS.
C
C  NAME    (INPUT) CHARACTER*(*)
!          THE NAME OF THE CALLING SUBROUTINE, IN EITHER UPPER CASE OR
!          LOWER CASE.
C
C  OPTS    (INPUT) CHARACTER*(*)
!          THE CHARACTER OPTIONS TO THE SUBROUTINE NAME, CONCATENATED
!          INTO A SINGLE CHARACTER STRING.  FOR EXAMPLE, UPLO = 'U',
!          TRANS = 'T', AND DIAG = 'N' FOR A TRIANGULAR ROUTINE WOULD
!          BE SPECIFIED AS OPTS = 'UTN'.
C
C  N1      (INPUT) INTEGER
C  N2      (INPUT) INTEGER
C  N3      (INPUT) INTEGER
C  N4      (INPUT) INTEGER
!          PROBLEM DIMENSIONS FOR THE SUBROUTINE NAME; THESE MAY NOT ALL
!          BE REQUIRED.
C
C (ILAENV) (OUTPUT) INTEGER
!          >= 0: THE VALUE OF THE PARAMETER SPECIFIED BY ISPEC
!          < 0:  IF ILAENV = -K, THE K-TH ARGUMENT HAD AN ILLEGAL VALUE.
C
C  FURTHER DETAILS
C  ===============
C
C  THE FOLLOWING CONVENTIONS HAVE BEEN USED WHEN CALLING ILAENV FROM THE
C  LAPACK ROUTINES:
C  1)  OPTS IS A CONCATENATION OF ALL OF THE CHARACTER OPTIONS TO
!      SUBROUTINE NAME, IN THE SAME ORDER THAT THEY APPEAR IN THE
!      ARGUMENT LIST FOR NAME, EVEN IF THEY ARE NOT USED IN DETERMINING
!      THE VALUE OF THE PARAMETER SPECIFIED BY ISPEC.
C  2)  THE PROBLEM DIMENSIONS N1, N2, N3, N4 ARE SPECIFIED IN THE ORDER
!      THAT THEY APPEAR IN THE ARGUMENT LIST FOR NAME.  N1 IS USED
!      FIRST, N2 SECOND, AND SO ON, AND UNUSED PROBLEM DIMENSIONS ARE
!      PASSED A VALUE OF -1.
C  3)  THE PARAMETER VALUE RETURNED BY ILAENV IS CHECKED FOR VALIDITY IN
!      THE CALLING SUBROUTINE.  FOR EXAMPLE, ILAENV IS USED TO RETRIEVE
!      THE OPTIMAL BLOCKSIZE FOR STRTRI AS FOLLOWS:
C
!      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
!      IF( NB<=1 ) NB = MAX( 1, N )
C
C  =====================================================================
C
!     .. LOCAL SCALARS ..
      LOGICAL            CNAME, SNAME
      CHARACTER*1        C1
      CHARACTER*2        C2, C4
      CHARACTER*3        C3
      CHARACTER*6        SUBNAM
      INTEGER            I, IC, IZ, NB, NBMIN, NX
!     ..
!     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
!     ..
!     .. EXECUTABLE STATEMENTS ..
C
!         NEXT BIT ADDED BY MWS TO SUPPRESS FTNCHEK WARNINGS
      IF(ISPEC>8) THEN
!change WRITE(6,*) 'DUMMY MESSAGE TO USE ARGS OPTS AND N3',N3,OPTS
      END IF
C
      GO TO ( 100, 100, 100, 400, 500, 600, 700, 800 ) ISPEC
C
!     INVALID VALUE FOR ISPEC
C
      ILAENV = -1
      RETURN
C
  100 CONTINUE
C
!     CONVERT NAME TO UPPER CASE IF THE FIRST CHARACTER IS LOWER CASE.
C
      ILAENV = 1
      SUBNAM = NAME
      IC = ICHAR( SUBNAM( 1:1 ) )
      IZ = ICHAR( 'Z' )
      IF( IZ==90 .OR. IZ==122 ) THEN
C
!        ASCII CHARACTER SET
C
         IF( IC.GE.97 .AND. IC<=122 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 10 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.97 .AND. IC<=122 )
     $            SUBNAM( I:I ) = CHAR( IC-32 )
   10       CONTINUE
         END IF
C
      ELSE IF( IZ==233 .OR. IZ==169 ) THEN
C
!        EBCDIC CHARACTER SET
C
         IF( ( IC.GE.129 .AND. IC<=137 ) .OR.
     $       ( IC.GE.145 .AND. IC<=153 ) .OR.
     $       ( IC.GE.162 .AND. IC<=169 ) ) THEN
            SUBNAM( 1:1 ) = CHAR( IC+64 )
            DO 20 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( ( IC.GE.129 .AND. IC<=137 ) .OR.
     $             ( IC.GE.145 .AND. IC<=153 ) .OR.
     $             ( IC.GE.162 .AND. IC<=169 ) )
     $            SUBNAM( I:I ) = CHAR( IC+64 )
   20       CONTINUE
         END IF
C
      ELSE IF( IZ==218 .OR. IZ==250 ) THEN
C
!        PRIME MACHINES:  ASCII+128
C
         IF( IC.GE.225 .AND. IC<=250 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 30 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.225 .AND. IC<=250 )
     $            SUBNAM( I:I ) = CHAR( IC-32 )
   30       CONTINUE
         END IF
      END IF
C
      C1 = SUBNAM( 1:1 )
      SNAME = C1=='S' .OR. C1=='D'
      CNAME = C1=='C' .OR. C1=='Z'
      IF( .NOT.( CNAME .OR. SNAME ) )
     $   RETURN
      C2 = SUBNAM( 2:3 )
      C3 = SUBNAM( 4:6 )
      C4 = C3( 2:3 )
C
      GO TO ( 110, 200, 300 ) ISPEC
C
  110 CONTINUE
C
!     ISPEC = 1:  BLOCK SIZE
C
!     IN THESE EXAMPLES, SEPARATE CODE IS PROVIDED FOR SETTING NB FOR
!     REAL AND COMPLEX.  WE ASSUME THAT NB WILL TAKE THE SAME VALUE IN
!     SINGLE OR DOUBLE PRECISION.
C
      NB = 1
C
      IF( C2=='GE' ) THEN
         IF( C3=='TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( C3=='QRF' .OR. C3=='RQF' .OR. C3=='LQF' .OR.
     $            C3=='QLF' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3=='HRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3=='BRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3=='TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2=='PO' ) THEN
         IF( C3=='TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2=='SY' ) THEN
         IF( C3=='TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( SNAME .AND. C3=='TRD' ) THEN
            NB = 1
         ELSE IF( SNAME .AND. C3=='GST' ) THEN
            NB = 64
         END IF
      ELSE IF( CNAME .AND. C2=='HE' ) THEN
         IF( C3=='TRF' ) THEN
            NB = 64
         ELSE IF( C3=='TRD' ) THEN
            NB = 1
         ELSE IF( C3=='GST' ) THEN
            NB = 64
         END IF
      ELSE IF( SNAME .AND. C2=='OR' ) THEN
         IF( C3( 1:1 )=='G' ) THEN
            IF( C4=='QR' .OR. C4=='RQ' .OR. C4=='LQ' .OR.
     $          C4=='QL' .OR. C4=='HR' .OR. C4=='TR' .OR.
     $          C4=='BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 )=='M' ) THEN
            IF( C4=='QR' .OR. C4=='RQ' .OR. C4=='LQ' .OR.
     $          C4=='QL' .OR. C4=='HR' .OR. C4=='TR' .OR.
     $          C4=='BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( CNAME .AND. C2=='UN' ) THEN
         IF( C3( 1:1 )=='G' ) THEN
            IF( C4=='QR' .OR. C4=='RQ' .OR. C4=='LQ' .OR.
     $          C4=='QL' .OR. C4=='HR' .OR. C4=='TR' .OR.
     $          C4=='BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 )=='M' ) THEN
            IF( C4=='QR' .OR. C4=='RQ' .OR. C4=='LQ' .OR.
     $          C4=='QL' .OR. C4=='HR' .OR. C4=='TR' .OR.
     $          C4=='BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( C2=='GB' ) THEN
         IF( C3=='TRF' ) THEN
            IF( SNAME ) THEN
               IF( N4<=64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N4<=64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2=='PB' ) THEN
         IF( C3=='TRF' ) THEN
            IF( SNAME ) THEN
               IF( N2<=64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N2<=64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2=='TR' ) THEN
         IF( C3=='TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2=='LA' ) THEN
         IF( C3=='UUM' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( SNAME .AND. C2=='ST' ) THEN
         IF( C3=='EBZ' ) THEN
            NB = 1
         END IF
      END IF
      ILAENV = NB
      RETURN
C
  200 CONTINUE
C
!     ISPEC = 2:  MINIMUM BLOCK SIZE
C
      NBMIN = 2
      IF( C2=='GE' ) THEN
         IF( C3=='QRF' .OR. C3=='RQF' .OR. C3=='LQF' .OR.
     $       C3=='QLF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3=='HRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3=='BRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3=='TRI' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         END IF
      ELSE IF( C2=='SY' ) THEN
         IF( C3=='TRF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( SNAME .AND. C3=='TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( CNAME .AND. C2=='HE' ) THEN
         IF( C3=='TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( SNAME .AND. C2=='OR' ) THEN
         IF( C3( 1:1 )=='G' ) THEN
            IF( C4=='QR' .OR. C4=='RQ' .OR. C4=='LQ' .OR.
     $          C4=='QL' .OR. C4=='HR' .OR. C4=='TR' .OR.
     $          C4=='BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 )=='M' ) THEN
            IF( C4=='QR' .OR. C4=='RQ' .OR. C4=='LQ' .OR.
     $          C4=='QL' .OR. C4=='HR' .OR. C4=='TR' .OR.
     $          C4=='BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      ELSE IF( CNAME .AND. C2=='UN' ) THEN
         IF( C3( 1:1 )=='G' ) THEN
            IF( C4=='QR' .OR. C4=='RQ' .OR. C4=='LQ' .OR.
     $          C4=='QL' .OR. C4=='HR' .OR. C4=='TR' .OR.
     $          C4=='BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 )=='M' ) THEN
            IF( C4=='QR' .OR. C4=='RQ' .OR. C4=='LQ' .OR.
     $          C4=='QL' .OR. C4=='HR' .OR. C4=='TR' .OR.
     $          C4=='BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      END IF
      ILAENV = NBMIN
      RETURN
C
  300 CONTINUE
C
!     ISPEC = 3:  CROSSOVER POINT
C
      NX = 0
      IF( C2=='GE' ) THEN
         IF( C3=='QRF' .OR. C3=='RQF' .OR. C3=='LQF' .OR.
     $       C3=='QLF' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3=='HRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3=='BRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         END IF
      ELSE IF( C2=='SY' ) THEN
         IF( SNAME .AND. C3=='TRD' ) THEN
            NX = 1
         END IF
      ELSE IF( CNAME .AND. C2=='HE' ) THEN
         IF( C3=='TRD' ) THEN
            NX = 1
         END IF
      ELSE IF( SNAME .AND. C2=='OR' ) THEN
         IF( C3( 1:1 )=='G' ) THEN
            IF( C4=='QR' .OR. C4=='RQ' .OR. C4=='LQ' .OR.
     $          C4=='QL' .OR. C4=='HR' .OR. C4=='TR' .OR.
     $          C4=='BR' ) THEN
               NX = 128
            END IF
         END IF
      ELSE IF( CNAME .AND. C2=='UN' ) THEN
         IF( C3( 1:1 )=='G' ) THEN
            IF( C4=='QR' .OR. C4=='RQ' .OR. C4=='LQ' .OR.
     $          C4=='QL' .OR. C4=='HR' .OR. C4=='TR' .OR.
     $          C4=='BR' ) THEN
               NX = 128
            END IF
         END IF
      END IF
      ILAENV = NX
      RETURN
C
  400 CONTINUE
C
!     ISPEC = 4:  NUMBER OF SHIFTS (USED BY XHSEQR)
C
      ILAENV = 6
      RETURN
C
  500 CONTINUE
C
!     ISPEC = 5:  MINIMUM COLUMN DIMENSION (NOT USED)
C
      ILAENV = 2
      RETURN
C
  600 CONTINUE 
C
!     ISPEC = 6:  CROSSOVER POINT FOR SVD (USED BY XGELSS AND XGESVD)
C
      ILAENV = INT( MIN(N1,N2)*1.6D+00 )
      RETURN
C
  700 CONTINUE
C
!     ISPEC = 7:  NUMBER OF PROCESSORS (NOT USED)
C
      ILAENV = 1
      RETURN
C
  800 CONTINUE
C
!     ISPEC = 8:  CROSSOVER POINT FOR MULTISHIFT (USED BY XHSEQR)
C
      ILAENV = 50
      RETURN
C
!     END OF ILAENV
C
      END

! DASUM
      DOUBLE PRECISION FUNCTION DASUM(N,DX,INCX)
*
*  -- Reference BLAS level1 routine (version 3.8.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2017
*
*     .. Scalar Arguments ..
      INTEGER incx,n
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION dx(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION dtemp
      INTEGER i,m,mp1,nincx
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC dabs,mod
*     ..
      dasum = 0.0d0
      dtemp = 0.0d0
      IF (n.LE.0 .OR. incx.LE.0) RETURN
      IF (incx.EQ.1) THEN
*        code for increment equal to 1
*
*
*        clean-up loop
*
         m = mod(n,6)
         IF (m.NE.0) THEN
            DO i = 1,m
               dtemp = dtemp + dabs(dx(i))
            END DO
            IF (n.LT.6) THEN
               dasum = dtemp
               RETURN
            END IF
         END IF
         mp1 = m + 1
         DO i = mp1,n,6
            dtemp = dtemp + dabs(dx(i)) + dabs(dx(i+1)) +
     $              dabs(dx(i+2)) + dabs(dx(i+3)) +
     $              dabs(dx(i+4)) + dabs(dx(i+5))
         END DO
      ELSE
*
*        code for increment not equal to 1
*
         nincx = n*incx
         DO i = 1,nincx,incx
            dtemp = dtemp + dabs(dx(i))
         END DO
      END IF
      dasum = dtemp
      RETURN
      END

! DROT
      SUBROUTINE DROT(N,DX,INCX,DY,INCY,C,S)
*
*  -- Reference BLAS level1 routine (version 3.8.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2017
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION C,S
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION DTEMP
      INTEGER I,IX,IY
*     ..
      IF (n.LE.0) RETURN
      IF (incx.EQ.1 .AND. incy.EQ.1) THEN
*
*       code for both increments equal to 1
*
         DO i = 1,n
            dtemp = c*dx(i) + s*dy(i)
            dy(i) = c*dy(i) - s*dx(i)
            dx(i) = dtemp
         END DO
      ELSE
*
*       code for unequal increments or equal increments not equal
*         to 1
*
         ix = 1
         iy = 1
         IF (incx.LT.0) ix = (-n+1)*incx + 1
         IF (incy.LT.0) iy = (-n+1)*incy + 1
         DO i = 1,n
            dtemp = c*dx(ix) + s*dy(iy)
            dy(iy) = c*dy(iy) - s*dx(ix)
            dx(ix) = dtemp
            ix = ix + incx
            iy = iy + incy
         END DO
      END IF
      RETURN
      END

! DROTG
      SUBROUTINE DROTG(DA,DB,C,S)
C
!     CONSTRUCT GIVENS PLANE ROTATION.
!     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DA,DB,C,S,ROE,SCALE,R,Z
      DOUBLE PRECISION ZERO, ONE
C
      PARAMETER (ZERO=0.0D+00, ONE=1.0D+00)
C
!-----------------------------------------------------------------------
C
C
      ROE = DB
      IF( ABS(DA) > ABS(DB) ) ROE = DA
      SCALE = ABS(DA) + ABS(DB)
      IF( SCALE /= ZERO ) GO TO 10
         C = ONE
         S = ZERO
         R = ZERO
         GO TO 20
C
   10 R = SCALE*SQRT((DA/SCALE)**2 + (DB/SCALE)**2)
      R = SIGN(ONE,ROE)*R
      C = DA/R
      S = DB/R
   20 Z = ONE
      IF( ABS(DA) > ABS(DB) ) Z = S
      IF( ABS(DB) .GE. ABS(DA) .AND. C /= ZERO ) Z = ONE/C
      DA = R
      DB = Z
      RETURN
      END

! DTRMM
      SUBROUTINE DTRMM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,
     $                   B, LDB )
!     .. SCALAR ARGUMENTS ..
      CHARACTER*1      SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, LDA, LDB
      DOUBLE PRECISION   ALPHA
!     .. ARRAY ARGUMENTS ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
!     ..
!
!  PURPOSE
!  =======
!
!  DTRMM  PERFORMS ONE OF THE MATRIX-MATRIX OPERATIONS
!
!     B := ALPHA*OP( A )*B,   OR   B := ALPHA*B*OP( A ),
!
!  WHERE  ALPHA  IS A SCALAR,  B  IS AN M BY N MATRIX,  A  IS A UNIT, OR
!  NON-UNIT,  UPPER OR LOWER TRIANGULAR MATRIX  AND  OP( A )  IS ONE  OF
!
!     OP( A ) = A   OR   OP( A ) = A'.
!
!  PARAMETERS
!  ==========
!
!  SIDE   - CHARACTER*1.
!           ON ENTRY,  SIDE SPECIFIES WHETHER  OP( A ) MULTIPLIES B FROM
!           THE LEFT OR RIGHT AS FOLLOWS:
!
!              SIDE = 'L' OR 'L'   B := ALPHA*OP( A )*B.
!
!              SIDE = 'R' OR 'R'   B := ALPHA*B*OP( A ).
!
!           UNCHANGED ON EXIT.
!
!  UPLO   - CHARACTER*1.
!           ON ENTRY, UPLO SPECIFIES WHETHER THE MATRIX A IS AN UPPER OR
!           LOWER TRIANGULAR MATRIX AS FOLLOWS:
!
!              UPLO = 'U' OR 'U'   A IS AN UPPER TRIANGULAR MATRIX.
!
!              UPLO = 'L' OR 'L'   A IS A LOWER TRIANGULAR MATRIX.
!
!           UNCHANGED ON EXIT.
!
!  TRANSA - CHARACTER*1.
!           ON ENTRY, TRANSA SPECIFIES THE FORM OF OP( A ) TO BE USED IN
!           THE MATRIX MULTIPLICATION AS FOLLOWS:
!
!              TRANSA = 'N' OR 'N'   OP( A ) = A.
!
!              TRANSA = 'T' OR 'T'   OP( A ) = A'.
!
!              TRANSA = 'C' OR 'C'   OP( A ) = A'.
!
!           UNCHANGED ON EXIT.
!
!  DIAG   - CHARACTER*1.
!           ON ENTRY, DIAG SPECIFIES WHETHER OR NOT A IS UNIT TRIANGULAR
!           AS FOLLOWS:
!
!              DIAG = 'U' OR 'U'   A IS ASSUMED TO BE UNIT TRIANGULAR.
!
!              DIAG = 'N' OR 'N'   A IS NOT ASSUMED TO BE UNIT
!                                  TRIANGULAR.
!
!           UNCHANGED ON EXIT.
!
!  M      - INTEGER.
!           ON ENTRY, M SPECIFIES THE NUMBER OF ROWS OF B. M MUST BE AT
!           LEAST ZERO.
!           UNCHANGED ON EXIT.
!
!  N      - INTEGER.
!           ON ENTRY, N SPECIFIES THE NUMBER OF COLUMNS OF B.  N MUST BE
!           AT LEAST ZERO.
!           UNCHANGED ON EXIT.
!
!  ALPHA  - DOUBLE PRECISION.
!           ON ENTRY,  ALPHA SPECIFIES THE SCALAR  ALPHA. WHEN  ALPHA IS
!           ZERO THEN  A IS NOT REFERENCED AND  B NEED NOT BE SET BEFORE
!           ENTRY.
!           UNCHANGED ON EXIT.
!
!  A      - DOUBLE PRECISION ARRAY OF DIMENSION ( LDA, K ), WHERE K IS M
!           WHEN  SIDE = 'L' OR 'L'  AND IS  N  WHEN  SIDE = 'R' OR 'R'.
!           BEFORE ENTRY  WITH  UPLO = 'U' OR 'U',  THE  LEADING  K BY K
!           UPPER TRIANGULAR PART OF THE ARRAY  A MUST CONTAIN THE UPPER
!           TRIANGULAR MATRIX  AND THE STRICTLY LOWER TRIANGULAR PART OF
!           A IS NOT REFERENCED.
!           BEFORE ENTRY  WITH  UPLO = 'L' OR 'L',  THE  LEADING  K BY K
!           LOWER TRIANGULAR PART OF THE ARRAY  A MUST CONTAIN THE LOWER
!           TRIANGULAR MATRIX  AND THE STRICTLY UPPER TRIANGULAR PART OF
!           A IS NOT REFERENCED.
!           NOTE THAT WHEN  DIAG = 'U' OR 'U',  THE DIAGONAL ELEMENTS OF
!           A  ARE NOT REFERENCED EITHER,  BUT ARE ASSUMED TO BE  UNITY.
!           UNCHANGED ON EXIT.
!
!  LDA    - INTEGER.
!           ON ENTRY, LDA SPECIFIES THE FIRST DIMENSION OF A AS DECLARED
!           IN THE CALLING (SUB) PROGRAM.  WHEN  SIDE = 'L' OR 'L'  THEN
!           LDA  MUST BE AT LEAST  MAX( 1, M ),  WHEN  SIDE = 'R' OR 'R'
!           THEN LDA MUST BE AT LEAST MAX( 1, N ).
!           UNCHANGED ON EXIT.
!
!  B      - DOUBLE PRECISION ARRAY OF DIMENSION ( LDB, N ).
!           BEFORE ENTRY,  THE LEADING  M BY N PART OF THE ARRAY  B MUST
!           CONTAIN THE MATRIX  B,  AND  ON EXIT  IS OVERWRITTEN  BY THE
!           TRANSFORMED MATRIX.
!
!  LDB    - INTEGER.
!           ON ENTRY, LDB SPECIFIES THE FIRST DIMENSION OF B AS DECLARED
!           IN  THE  CALLING  (SUB)  PROGRAM.   LDB  MUST  BE  AT  LEAST
!           MAX( 1, M ).
!           UNCHANGED ON EXIT.
!
!
!  LEVEL 3 BLAS ROUTINE.
!
!  -- WRITTEN ON 8-FEBRUARY-1989.
!     JACK DONGARRA, ARGONNE NATIONAL LABORATORY.
!     IAIN DUFF, AERE HARWELL.
!     JEREMY DU CROZ, NUMERICAL ALGORITHMS GROUP LTD.
!     SVEN HAMMARLING, NUMERICAL ALGORITHMS GROUP LTD.
!
!
!     .. EXTERNAL FUNCTIONS ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     .. EXTERNAL ROUTINES ..
      EXTERNAL           XERBLA
!     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          MAX
!     .. LOCAL SCALARS ..
      LOGICAL            LSIDE, NOUNIT, UPPER
      INTEGER            I, INFO, J, K, NROWA
      DOUBLE PRECISION   TEMP
!     .. PARAMETERS ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. EXECUTABLE STATEMENTS ..
!
!     TEST THE INPUT PARAMETERS.
!
      LSIDE  = LSAME( SIDE  , 'L' )
      IF( LSIDE )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      NOUNIT = LSAME( DIAG  , 'N' )
      UPPER  = LSAME( UPLO  , 'U' )
!
      INFO   = 0
      IF(      ( .NOT.LSIDE                ).AND.
     $         ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER                ).AND.
     $         ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
         INFO = 2
      ELSE IF( ( .NOT.LSAME( TRANSA, 'N' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
         INFO = 3
      ELSE IF( ( .NOT.LSAME( DIAG  , 'U' ) ).AND.
     $         ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
         INFO = 4
      ELSE IF( M  <0               )THEN
         INFO = 5
      ELSE IF( N  <0               )THEN
         INFO = 6
      ELSE IF( LDA<MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDB<MAX( 1, M     ) )THEN
         INFO = 11
      END IF
      IF( INFO/=0 )THEN
         CALL XERBLA( 'DTRMM ', INFO )
         RETURN
      END IF
!
!     QUICK RETURN IF POSSIBLE.
!
      IF( N==0 )
     $   RETURN
!
!     AND WHEN  ALPHA==ZERO.
!
      IF( ALPHA==ZERO )THEN
         DO 20, J = 1, N
            DO 10, I = 1, M
               B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
         RETURN
      END IF
!
!     START THE OPERATIONS.
!
      IF( LSIDE )THEN
         IF( LSAME( TRANSA, 'N' ) )THEN
!
!           FORM  B := ALPHA*A*B.
!
            IF( UPPER )THEN
               DO 50, J = 1, N
                  DO 40, K = 1, M
                     IF( B( K, J )/=ZERO )THEN
                        TEMP = ALPHA*B( K, J )
                        DO 30, I = 1, K - 1
                           B( I, J ) = B( I, J ) + TEMP*A( I, K )
   30                   CONTINUE
                        IF( NOUNIT )
     $                     TEMP = TEMP*A( K, K )
                        B( K, J ) = TEMP
                     END IF
   40             CONTINUE
   50          CONTINUE
            ELSE
               DO 80, J = 1, N
                  DO 70 K = M, 1, -1
                     IF( B( K, J )/=ZERO )THEN
                        TEMP      = ALPHA*B( K, J )
                        B( K, J ) = TEMP
                        IF( NOUNIT )
     $                     B( K, J ) = B( K, J )*A( K, K )
                        DO 60, I = K + 1, M
                           B( I, J ) = B( I, J ) + TEMP*A( I, K )
   60                   CONTINUE
                     END IF
   70             CONTINUE
   80          CONTINUE
            END IF
         ELSE
!
!           FORM  B := ALPHA*B*A'.
!
            IF( UPPER )THEN
               DO 110, J = 1, N
                  DO 100, I = M, 1, -1
                     TEMP = B( I, J )
                     IF( NOUNIT )
     $                  TEMP = TEMP*A( I, I )
                     DO 90, K = 1, I - 1
                        TEMP = TEMP + A( K, I )*B( K, J )
   90                CONTINUE
                     B( I, J ) = ALPHA*TEMP
  100             CONTINUE
  110          CONTINUE
            ELSE
               DO 140, J = 1, N
                  DO 130, I = 1, M
                     TEMP = B( I, J )
                     IF( NOUNIT )
     $                  TEMP = TEMP*A( I, I )
                     DO 120, K = I + 1, M
                        TEMP = TEMP + A( K, I )*B( K, J )
  120                CONTINUE
                     B( I, J ) = ALPHA*TEMP
  130             CONTINUE
  140          CONTINUE
            END IF
         END IF
      ELSE
         IF( LSAME( TRANSA, 'N' ) )THEN
!
!           FORM  B := ALPHA*B*A.
!
            IF( UPPER )THEN
               DO 180, J = N, 1, -1
                  TEMP = ALPHA
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 150, I = 1, M
                     B( I, J ) = TEMP*B( I, J )
  150             CONTINUE
                  DO 170, K = 1, J - 1
                     IF( A( K, J )/=ZERO )THEN
                        TEMP = ALPHA*A( K, J )
                        DO 160, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  160                   CONTINUE
                     END IF
  170             CONTINUE
  180          CONTINUE
            ELSE
               DO 220, J = 1, N
                  TEMP = ALPHA
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 190, I = 1, M
                     B( I, J ) = TEMP*B( I, J )
  190             CONTINUE
                  DO 210, K = J + 1, N
                     IF( A( K, J )/=ZERO )THEN
                        TEMP = ALPHA*A( K, J )
                        DO 200, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  200                   CONTINUE
                     END IF
  210             CONTINUE
  220          CONTINUE
            END IF
         ELSE
!
!           FORM  B := ALPHA*B*A'.
!
            IF( UPPER )THEN
               DO 260, K = 1, N
                  DO 240, J = 1, K - 1
                     IF( A( J, K )/=ZERO )THEN
                        TEMP = ALPHA*A( J, K )
                        DO 230, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  230                   CONTINUE
                     END IF
  240             CONTINUE
                  TEMP = ALPHA
                  IF( NOUNIT )
     $               TEMP = TEMP*A( K, K )
                  IF( TEMP/=ONE )THEN
                     DO 250, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  250                CONTINUE
                  END IF
  260          CONTINUE
            ELSE
               DO 300, K = N, 1, -1
                  DO 280, J = K + 1, N
                     IF( A( J, K )/=ZERO )THEN
                        TEMP = ALPHA*A( J, K )
                        DO 270, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  270                   CONTINUE
                     END IF
  280             CONTINUE
                  TEMP = ALPHA
                  IF( NOUNIT )
     $               TEMP = TEMP*A( K, K )
                  IF( TEMP/=ONE )THEN
                     DO 290, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  290                CONTINUE
                  END IF
  300          CONTINUE
            END IF
         END IF
      END IF
!
      RETURN
!
!     END OF DTRMM .
!
      END

! IDAMAX
      INTEGER FUNCTION IDAMAX(N,DX,INCX)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION DX(*)
!
!     FINDS THE INDEX OF ELEMENT HAVING MAX. ABSOLUTE VALUE.
!     JACK DONGARRA, LINPACK, 3/11/78.
!
      IDAMAX = 0
      IF( N < 1 ) RETURN
      IDAMAX = 1
      IF(N==1) RETURN
      IF(INCX==1)GO TO 20
!
!        CODE FOR INCREMENT NOT EQUAL TO 1
!
      IX = 1
      RMAX = ABS(DX(1))
      IX = IX + INCX
      DO 10 I = 2,N
         IF(ABS(DX(IX))<=RMAX) GO TO 5
         IDAMAX = I
         RMAX = ABS(DX(IX))
    5    IX = IX + INCX
   10 CONTINUE
      RETURN
!
!        CODE FOR INCREMENT EQUAL TO 1
!
   20 RMAX = ABS(DX(1))
      DO 30 I = 2,N
         IF(ABS(DX(I))<=RMAX) GO TO 30
         IDAMAX = I
         RMAX = ABS(DX(I))
   30 CONTINUE
      RETURN
      END

! DAXPY
      SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
!
!  -- Reference BLAS level1 routine (version 3.8.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2017
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION DA
      INTEGER INCX,INCY,N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER I,IX,IY,M,MP1
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC mod
!     ..
      IF (n.LE.0) RETURN
      IF (da.EQ.0.0d0) RETURN
      IF (incx.EQ.1 .AND. incy.EQ.1) THEN
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
         m = mod(n,4)
         IF (m.NE.0) THEN
            DO i = 1,m
               dy(i) = dy(i) + da*dx(i)
            END DO
         END IF
         IF (n.LT.4) RETURN
         mp1 = m + 1
         DO i = mp1,n,4
            dy(i) = dy(i) + da*dx(i)
            dy(i+1) = dy(i+1) + da*dx(i+1)
            dy(i+2) = dy(i+2) + da*dx(i+2)
            dy(i+3) = dy(i+3) + da*dx(i+3)
         END DO
      ELSE
!
!        code for unequal increments or equal increments
!          not equal to 1
!
         ix = 1
         iy = 1
         IF (incx.LT.0) ix = (-n+1)*incx + 1
         IF (incy.LT.0) iy = (-n+1)*incy + 1
         DO i = 1,n
          dy(iy) = dy(iy) + da*dx(ix)
          ix = ix + incx
          iy = iy + incy
         END DO
      END IF
      RETURN
      END
       
! DDOT
       DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
!
!  -- Reference BLAS level1 routine (version 3.8.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2017
!
!     .. Scalar Arguments ..
      INTEGER incx,incy,n
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION dx(*),dy(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      DOUBLE PRECISION dtemp
      INTEGER i,ix,iy,m,mp1
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC mod
!     ..
      ddot = 0.0d0
      dtemp = 0.0d0
      IF (n.LE.0) RETURN
      IF (incx.EQ.1 .AND. incy.EQ.1) THEN
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
         m = mod(n,5)
         IF (m.NE.0) THEN
            DO i = 1,m
               dtemp = dtemp + dx(i)*dy(i)
            END DO
            IF (n.LT.5) THEN
               ddot=dtemp
            RETURN
            END IF
         END IF
         mp1 = m + 1
         DO i = mp1,n,5
          dtemp = dtemp + dx(i)*dy(i) + dx(i+1)*dy(i+1) +
     $            dx(i+2)*dy(i+2) + dx(i+3)*dy(i+3) + dx(i+4)*dy(i+4)
         END DO
      ELSE
!
!        code for unequal increments or equal increments
!          not equal to 1
!
         ix = 1
         iy = 1
         IF (incx.LT.0) ix = (-n+1)*incx + 1
         IF (incy.LT.0) iy = (-n+1)*incy + 1
         DO i = 1,n
            dtemp = dtemp + dx(ix)*dy(iy)
            ix = ix + incx
            iy = iy + incy
         END DO
      END IF
      ddot = dtemp
      RETURN
      END
       
! DNRM2
       DOUBLE PRECISION FUNCTION dnrm2(N,X,INCX)
!
!  -- Reference BLAS level1 routine (version 3.8.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2017
!
!     .. Scalar Arguments ..
      INTEGER incx,n
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION x(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION one,zero
      parameter(one=1.0d+0,zero=0.0d+0)
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION absxi,norm,scale,ssq
      INTEGER ix
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC abs,sqrt
!     ..
      IF (n.LT.1 .OR. incx.LT.1) THEN
          norm = zero
      ELSE IF (n.EQ.1) THEN
          norm = abs(x(1))
      ELSE
          scale = zero
          ssq = one
!        The following loop is equivalent to this call to the LAPACK
!        auxiliary routine:
!        CALL DLASSQ( N, X, INCX, SCALE, SSQ )
!
          DO 10 ix = 1,1 + (n-1)*incx,incx
              IF (x(ix).NE.zero) THEN
                  absxi = abs(x(ix))
                  IF (scale.LT.absxi) THEN
                      ssq = one + ssq* (scale/absxi)**2
                      scale = absxi
                  ELSE
                      ssq = ssq + (absxi/scale)**2
                  END IF
              END IF
   10     CONTINUE
          norm = scale*sqrt(ssq)
      END IF
!
      dnrm2 = norm
      RETURN
!
!     End of DNRM2.
!
      END

! DGEMV
      SUBROUTINE DGEMV(FORMA,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*1 FORMA
      DIMENSION A(LDA,*),X(*),Y(*)
      PARAMETER (ZERO=0.0D+00, ONE=1.0D+00)
!
!        CLONE OF -DGEMV- WRITTEN BY MIKE SCHMIDT
!
      LOCY = 1
!
!                  Y = ALPHA * A * X + BETA * Y
!
      IF(FORMA=='N') THEN
         IF(ALPHA==ONE  .AND.  BETA==ZERO) THEN
            DO 110 I=1,M
               Y(LOCY) =       DDOT(N,A(I,1),LDA,X,INCX)
               LOCY = LOCY+INCY
  110       CONTINUE
         ELSE
            DO 120 I=1,M
               Y(LOCY) = ALPHA*DDOT(N,A(I,1),LDA,X,INCX) + BETA*Y(LOCY)
               LOCY = LOCY+INCY
  120       CONTINUE
         END IF
         RETURN
      END IF
!
!                  Y = ALPHA * A-TRANSPOSE * X + BETA * Y
!
      IF(FORMA=='T') THEN
         IF(ALPHA==ONE  .AND.  BETA==ZERO) THEN
            DO 210 I=1,N
               Y(LOCY) =       DDOT(M,A(1,I),1,X,INCX)
               LOCY = LOCY+INCY
  210       CONTINUE
         ELSE
            DO 220 I=1,N
               Y(LOCY) = ALPHA*DDOT(M,A(1,I),1,X,INCX) + BETA*Y(LOCY)
               LOCY = LOCY+INCY
  220       CONTINUE
         END IF
         RETURN
      END IF
!
      WRITE(6,900) FORMA
      CALL ABRT
      RETURN
  900 FORMAT(1X,'ERROR IN -DGEMV- ... UNRECOGNIZED FORMA=',A1)
      END

! DTRMV
      SUBROUTINE DTRMV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
!     .. SCALAR ARGUMENTS ..
      INTEGER            INCX, LDA, N
      CHARACTER*1      DIAG, TRANS, UPLO
!     .. ARRAY ARGUMENTS ..
      DOUBLE PRECISION   A( LDA, * ), X( * )
!     ..
!
!  PURPOSE
!  =======
!
!  DTRMV  PERFORMS ONE OF THE MATRIX-VECTOR OPERATIONS
!
!     X := A*X,   OR   X := A'*X,
!
!  WHERE X IS AN N ELEMENT VECTOR AND  A IS AN N BY N UNIT, OR NON-UNIT,
!  UPPER OR LOWER TRIANGULAR MATRIX.
!
!  PARAMETERS
!  ==========
!
!  UPLO   - CHARACTER*1.
!           ON ENTRY, UPLO SPECIFIES WHETHER THE MATRIX IS AN UPPER OR
!           LOWER TRIANGULAR MATRIX AS FOLLOWS:
!
!              UPLO = 'U' OR 'U'   A IS AN UPPER TRIANGULAR MATRIX.
!
!              UPLO = 'L' OR 'L'   A IS A LOWER TRIANGULAR MATRIX.
!
!           UNCHANGED ON EXIT.
!
!  TRANS  - CHARACTER*1.
!           ON ENTRY, TRANS SPECIFIES THE OPERATION TO BE PERFORMED AS
!           FOLLOWS:
!
!              TRANS = 'N' OR 'N'   X := A*X.
!
!              TRANS = 'T' OR 'T'   X := A'*X.
!
!              TRANS = 'C' OR 'C'   X := A'*X.
!
!           UNCHANGED ON EXIT.
!
!  DIAG   - CHARACTER*1.
!           ON ENTRY, DIAG SPECIFIES WHETHER OR NOT A IS UNIT
!           TRIANGULAR AS FOLLOWS:
!
!              DIAG = 'U' OR 'U'   A IS ASSUMED TO BE UNIT TRIANGULAR.
!
!              DIAG = 'N' OR 'N'   A IS NOT ASSUMED TO BE UNIT
!                                  TRIANGULAR.
!
!           UNCHANGED ON EXIT.
!
!  N      - INTEGER.
!           ON ENTRY, N SPECIFIES THE ORDER OF THE MATRIX A.
!           N MUST BE AT LEAST ZERO.
!           UNCHANGED ON EXIT.
!
!  A      - DOUBLE PRECISION ARRAY OF DIMENSION ( LDA, N ).
!           BEFORE ENTRY WITH  UPLO = 'U' OR 'U', THE LEADING N BY N
!           UPPER TRIANGULAR PART OF THE ARRAY A MUST CONTAIN THE UPPER
!           TRIANGULAR MATRIX AND THE STRICTLY LOWER TRIANGULAR PART OF
!           A IS NOT REFERENCED.
!           BEFORE ENTRY WITH UPLO = 'L' OR 'L', THE LEADING N BY N
!           LOWER TRIANGULAR PART OF THE ARRAY A MUST CONTAIN THE LOWER
!           TRIANGULAR MATRIX AND THE STRICTLY UPPER TRIANGULAR PART OF
!           A IS NOT REFERENCED.
!           NOTE THAT WHEN  DIAG = 'U' OR 'U', THE DIAGONAL ELEMENTS OF
!           A ARE NOT REFERENCED EITHER, BUT ARE ASSUMED TO BE UNITY.
!           UNCHANGED ON EXIT.
!
!  LDA    - INTEGER.
!           ON ENTRY, LDA SPECIFIES THE FIRST DIMENSION OF A AS DECLARED
!           IN THE CALLING (SUB) PROGRAM. LDA MUST BE AT LEAST
!           MAX( 1, N ).
!           UNCHANGED ON EXIT.
!
!  X      - DOUBLE PRECISION ARRAY OF DIMENSION AT LEAST
!           ( 1 + ( N - 1 )*ABS( INCX ) ).
!           BEFORE ENTRY, THE INCREMENTED ARRAY X MUST CONTAIN THE N
!           ELEMENT VECTOR X. ON EXIT, X IS OVERWRITTEN WITH THE
!           TRANFORMED VECTOR X.
!
!  INCX   - INTEGER.
!           ON ENTRY, INCX SPECIFIES THE INCREMENT FOR THE ELEMENTS OF
!           X. INCX MUST NOT BE ZERO.
!           UNCHANGED ON EXIT.
!
!
!  LEVEL 2 BLAS ROUTINE.
!
!  -- WRITTEN ON 22-OCTOBER-1986.
!     JACK DONGARRA, JEREMY DU CROZ, SVEN HAMMARLING, RICHARD HANSON
!
!     .. PARAMETERS ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
!     .. LOCAL SCALARS ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, J, JX, KX
      LOGICAL            NOUNIT
!     .. EXTERNAL FUNCTIONS ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     .. EXTERNAL ROUTINES ..
      EXTERNAL           XERBLA
!     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          MAX
!     ..
!     .. EXECUTABLE STATEMENTS ..
!
!     TEST THE INPUT PARAMETERS.
!
      INFO = 0
      IF     ( .NOT.LSAME( UPLO , 'U' ).AND.
     $         .NOT.LSAME( UPLO , 'L' )      )THEN
         INFO = 1
      ELSE IF( .NOT.LSAME( TRANS, 'N' ).AND.
     $         .NOT.LSAME( TRANS, 'T' ).AND.
     $         .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 2
      ELSE IF( .NOT.LSAME( DIAG , 'U' ).AND.
     $         .NOT.LSAME( DIAG , 'N' )      )THEN
         INFO = 3
      ELSE IF( N<0 )THEN
         INFO = 4
      ELSE IF( LDA<MAX( 1, N ) )THEN
         INFO = 6
      ELSE IF( INCX==0 )THEN
         INFO = 8
      END IF
      IF( INFO/=0 )THEN
         CALL XERBLA( 'DTRMV ', INFO )
         RETURN
      END IF
!
!     QUICK RETURN IF POSSIBLE.
!
      IF( N==0 )
     $   RETURN
!
      NOUNIT = LSAME( DIAG, 'N' )
!
!     SET UP THE START POINT IN X IF THE INCREMENT IS NOT UNITY. THIS
!     WILL BE  ( N - 1 )*INCX  TOO SMALL FOR DESCENDING LOOPS.
!
      IF( INCX<=0 )THEN
         KX = 1 - ( N - 1 )*INCX
      ELSE IF( INCX/=1 )THEN
         KX = 1
      END IF
!
!     START THE OPERATIONS. IN THIS VERSION THE ELEMENTS OF A ARE
!     ACCESSED SEQUENTIALLY WITH ONE PASS THROUGH A.
!
      IF( LSAME( TRANS, 'N' ) )THEN
!
!        FORM  X := A*X.
!
         IF( LSAME( UPLO, 'U' ) )THEN
            IF( INCX==1 )THEN
               DO 20, J = 1, N
                  IF( X( J )/=ZERO )THEN
                     TEMP = X( J )
                     DO 10, I = 1, J - 1
                        X( I ) = X( I ) + TEMP*A( I, J )
   10                CONTINUE
                     IF( NOUNIT )
     $                  X( J ) = X( J )*A( J, J )
                  END IF
   20          CONTINUE
            ELSE
               JX = KX
               DO 40, J = 1, N
                  IF( X( JX )/=ZERO )THEN
                     TEMP = X( JX )
                     IX   = KX
                     DO 30, I = 1, J - 1
                        X( IX ) = X( IX ) + TEMP*A( I, J )
                        IX      = IX      + INCX
   30                CONTINUE
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )*A( J, J )
                  END IF
                  JX = JX + INCX
   40          CONTINUE
            END IF
         ELSE
            IF( INCX==1 )THEN
               DO 60, J = N, 1, -1
                  IF( X( J )/=ZERO )THEN
                     TEMP = X( J )
                     DO 50, I = N, J + 1, -1
                        X( I ) = X( I ) + TEMP*A( I, J )
   50                CONTINUE
                     IF( NOUNIT )
     $                  X( J ) = X( J )*A( J, J )
                  END IF
   60          CONTINUE
            ELSE
               KX = KX + ( N - 1 )*INCX
               JX = KX
               DO 80, J = N, 1, -1
                  IF( X( JX )/=ZERO )THEN
                     TEMP = X( JX )
                     IX   = KX
                     DO 70, I = N, J + 1, -1
                        X( IX ) = X( IX ) + TEMP*A( I, J )
                        IX      = IX      - INCX
   70                CONTINUE
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )*A( J, J )
                  END IF
                  JX = JX - INCX
   80          CONTINUE
            END IF
         END IF
      ELSE
!
!        FORM  X := A'*X.
!
         IF( LSAME( UPLO, 'U' ) )THEN
            IF( INCX==1 )THEN
               DO 100, J = N, 1, -1
                  TEMP = X( J )
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 90, I = J - 1, 1, -1
                     TEMP = TEMP + A( I, J )*X( I )
   90             CONTINUE
                  X( J ) = TEMP
  100          CONTINUE
            ELSE
               JX = KX + ( N - 1 )*INCX
               DO 120, J = N, 1, -1
                  TEMP = X( JX )
                  IX   = JX
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 110, I = J - 1, 1, -1
                     IX   = IX   - INCX
                     TEMP = TEMP + A( I, J )*X( IX )
  110             CONTINUE
                  X( JX ) = TEMP
                  JX      = JX   - INCX
  120          CONTINUE
            END IF
         ELSE
            IF( INCX==1 )THEN
               DO 140, J = 1, N
                  TEMP = X( J )
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 130, I = J + 1, N
                     TEMP = TEMP + A( I, J )*X( I )
  130             CONTINUE
                  X( J ) = TEMP
  140          CONTINUE
            ELSE
               JX = KX
               DO 160, J = 1, N
                  TEMP = X( JX )
                  IX   = JX
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 150, I = J + 1, N
                     IX   = IX   + INCX
                     TEMP = TEMP + A( I, J )*X( IX )
  150             CONTINUE
                  X( JX ) = TEMP
                  JX      = JX   + INCX
  160          CONTINUE
            END IF
         END IF
      END IF
!
      RETURN
!
!     END OF DTRMV .
!
      END

! DGEMM
      SUBROUTINE DGEMM(FORMA,FORMB,L,N,M,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*1 FORMA,FORMB
      DOUBLE PRECISION A(LDA,*), B(LDB,*), C(LDC,*)
      PARAMETER (ZERO=0.0D+00, ONE=1.0D+00)
!
!     THIS IS A PLAIN VANILLA FORTRAN CLONE OF DGEMM
!
      IF (FORMA=='N' .AND. FORMB=='N') THEN
         IF(ALPHA==ONE  .AND.  BETA==ZERO) THEN
            DO 30 IL = 1, L
               DO 20 IN = 1, N
                  T = ZERO
                  DO 10 IM = 1, M
                     T = T + A(IL,IM)*B(IM,IN)
   10             CONTINUE
                  C(IL,IN) = T
   20          CONTINUE
   30       CONTINUE
         ELSE
            DO 80 IL = 1, L
               DO 70 IN = 1, N
                  T = ZERO
                  DO 60 IM = 1, M
                     T = T + A(IL,IM)*B(IM,IN)
   60             CONTINUE
                  C(IL,IN) = ALPHA*T + BETA*C(IL,IN)
   70          CONTINUE
   80       CONTINUE
         END IF
         RETURN
      END IF
!
      IF (FORMA=='T' .AND. FORMB=='N') THEN
         IF(ALPHA==ONE  .AND.  BETA==ZERO) THEN
            DO 130 IL = 1, L
               DO 120 IN = 1, N
                  T = ZERO
                  DO 110 IM = 1, M
                     T = T + A(IM,IL)*B(IM,IN)
  110             CONTINUE
                  C(IL,IN) = T
  120          CONTINUE
  130       CONTINUE
         ELSE
            DO 180 IL = 1, L
               DO 170 IN = 1, N
                  T = ZERO
                  DO 160 IM = 1, M
                     T = T + A(IM,IL)*B(IM,IN)
  160             CONTINUE
                  C(IL,IN) = ALPHA*T + BETA*C(IL,IN)
  170          CONTINUE
  180       CONTINUE
         END IF
         RETURN
      END IF
!
      IF (FORMA=='N' .AND. FORMB=='T') THEN
         IF(ALPHA==ONE  .AND.  BETA==ZERO) THEN
            DO 230 IL = 1, L
               DO 220 IN = 1, N
                  T = ZERO
                  DO 210 IM = 1, M
                     T = T + A(IL,IM)*B(IN,IM)
  210             CONTINUE
                  C(IL,IN) = T
  220          CONTINUE
  230       CONTINUE
         ELSE
            DO 280 IL = 1, L
               DO 270 IN = 1, N
                  T = ZERO
                  DO 260 IM = 1, M
                     T = T + A(IL,IM)*B(IN,IM)
  260             CONTINUE
                  C(IL,IN) = ALPHA*T + BETA*C(IL,IN)
  270          CONTINUE
  280       CONTINUE
         END IF
         RETURN
      END IF
!
      IF (FORMA=='T' .AND. FORMB=='T') THEN
         IF(ALPHA==ONE  .AND.  BETA==ZERO) THEN
            DO 330 IL = 1, L
               DO 320 IN = 1, N
                  T = ZERO
                  DO 310 IM = 1, M
                     T = T + A(IM,IL)*B(IN,IM)
  310             CONTINUE
                  C(IL,IN) = T
  320          CONTINUE
  330       CONTINUE
         ELSE
            DO 380 IL = 1, L
               DO 370 IN = 1, N
                  T = ZERO
                  DO 360 IM = 1, M
                     T = T + A(IM,IL)*B(IN,IM)
  360             CONTINUE
                  C(IL,IN) = ALPHA*T + BETA*C(IL,IN)
  370          CONTINUE
  380       CONTINUE
         END IF
         RETURN
      END IF
!
      WRITE(6,900) FORMA,FORMB
      CALL ABRT
      RETURN
  900 FORMAT(1X,'ERROR IN -DGEMM- ... ILLEGAL FORMA/FORMB=',A1,1X,A1)
      END

! DCOPY
      SUBROUTINE DCOPY(N,DX,INCX,DY,INCY)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION DX(*),DY(*)
!
!     COPIES A VECTOR.
!           DY(I) <== DX(I)
!     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
!     JACK DONGARRA, LINPACK, 3/11/78.
!
      IF(N<=0) RETURN
      IF(INCX==1.AND.INCY==1)GO TO 20
!
!        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
!          NOT EQUAL TO 1
!
      IX = 1
      IY = 1
      IF(INCX<0)IX = (-N+1)*INCX + 1
      IF(INCY<0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
!
!        CODE FOR BOTH INCREMENTS EQUAL TO 1
!
!
!        CLEAN-UP LOOP
!
   20 M = MOD(N,7)
      IF( M == 0 ) GO TO 40
      DO 30 I = 1,M
        DY(I) = DX(I)
   30 CONTINUE
      IF( N < 7 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,7
        DY(I) = DX(I)
        DY(I + 1) = DX(I + 1)
        DY(I + 2) = DX(I + 2)
        DY(I + 3) = DX(I + 3)
        DY(I + 4) = DX(I + 4)
        DY(I + 5) = DX(I + 5)
        DY(I + 6) = DX(I + 6)
   50 CONTINUE
      RETURN
      END

! DGER
      SUBROUTINE DGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
!     .. SCALAR ARGUMENTS ..
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, INCY, LDA, M, N
!     .. ARRAY ARGUMENTS ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
!     ..
!
!  PURPOSE
!  =======
!
!  DGER   PERFORMS THE RANK 1 OPERATION
!
!     A := ALPHA*X*Y' + A,
!
!  WHERE ALPHA IS A SCALAR, X IS AN M ELEMENT VECTOR, Y IS AN N ELEMENT
!  VECTOR AND A IS AN M BY N MATRIX.
!
!  PARAMETERS
!  ==========
!
!  M      - INTEGER.
!           ON ENTRY, M SPECIFIES THE NUMBER OF ROWS OF THE MATRIX A.
!           M MUST BE AT LEAST ZERO.
!           UNCHANGED ON EXIT.
!
!  N      - INTEGER.
!           ON ENTRY, N SPECIFIES THE NUMBER OF COLUMNS OF THE MATRIX A.
!           N MUST BE AT LEAST ZERO.
!           UNCHANGED ON EXIT.
!
!  ALPHA  - DOUBLE PRECISION.
!           ON ENTRY, ALPHA SPECIFIES THE SCALAR ALPHA.
!           UNCHANGED ON EXIT.
!
!  X      - DOUBLE PRECISION ARRAY OF DIMENSION AT LEAST
!           ( 1 + ( M - 1 )*ABS( INCX ) ).
!           BEFORE ENTRY, THE INCREMENTED ARRAY X MUST CONTAIN THE M
!           ELEMENT VECTOR X.
!           UNCHANGED ON EXIT.
!
!  INCX   - INTEGER.
!           ON ENTRY, INCX SPECIFIES THE INCREMENT FOR THE ELEMENTS OF
!           X. INCX MUST NOT BE ZERO.
!           UNCHANGED ON EXIT.
!
!  Y      - DOUBLE PRECISION ARRAY OF DIMENSION AT LEAST
!           ( 1 + ( N - 1 )*ABS( INCY ) ).
!           BEFORE ENTRY, THE INCREMENTED ARRAY Y MUST CONTAIN THE N
!           ELEMENT VECTOR Y.
!           UNCHANGED ON EXIT.
!
!  INCY   - INTEGER.
!           ON ENTRY, INCY SPECIFIES THE INCREMENT FOR THE ELEMENTS OF
!           Y. INCY MUST NOT BE ZERO.
!           UNCHANGED ON EXIT.
!
!  A      - DOUBLE PRECISION ARRAY OF DIMENSION ( LDA, N ).
!           BEFORE ENTRY, THE LEADING M BY N PART OF THE ARRAY A MUST
!           CONTAIN THE MATRIX OF COEFFICIENTS. ON EXIT, A IS
!           OVERWRITTEN BY THE UPDATED MATRIX.
!
!  LDA    - INTEGER.
!           ON ENTRY, LDA SPECIFIES THE FIRST DIMENSION OF A AS DECLARED
!           IN THE CALLING (SUB) PROGRAM. LDA MUST BE AT LEAST
!           MAX( 1, M ).
!           UNCHANGED ON EXIT.
!
!
!  LEVEL 2 BLAS ROUTINE.
!
!  -- WRITTEN ON 22-OCTOBER-1986.
!     JACK DONGARRA, JEREMY DU CROZ, SVEN HAMMARLING, RICHARD HANSON
!
!
!     .. PARAMETERS ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
!     .. LOCAL SCALARS ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, J, JY, KX
!     .. EXTERNAL ROUTINES ..
      EXTERNAL           XERBLA
!     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          MAX
!     ..
!     .. EXECUTABLE STATEMENTS ..
!
!     TEST THE INPUT PARAMETERS.
!
      INFO = 0
      IF     ( M<0 )THEN
         INFO = 1
      ELSE IF( N<0 )THEN
         INFO = 2
      ELSE IF( INCX==0 )THEN
         INFO = 5
      ELSE IF( INCY==0 )THEN
         INFO = 7
      ELSE IF( LDA<MAX( 1, M ) )THEN
         INFO = 9
      END IF
      IF( INFO/=0 )THEN
         CALL XERBLA( 'DGER  ', INFO )
         RETURN
      END IF
!
!     QUICK RETURN IF POSSIBLE.
!
      IF( ( M==0 ).OR.( N==0 ).OR.( ALPHA==ZERO ) )
     $   RETURN
!
!     START THE OPERATIONS. IN THIS VERSION THE ELEMENTS OF A ARE
!     ACCESSED SEQUENTIALLY WITH ONE PASS THROUGH A.
!
      IF( INCY>0 )THEN
         JY = 1
      ELSE
         JY = 1 - ( N - 1 )*INCY
      END IF
      IF( INCX==1 )THEN
         DO 20, J = 1, N
            IF( Y( JY )/=ZERO )THEN
               TEMP = ALPHA*Y( JY )
               DO 10, I = 1, M
                  A( I, J ) = A( I, J ) + X( I )*TEMP
   10          CONTINUE
            END IF
            JY = JY + INCY
   20    CONTINUE
      ELSE
         IF( INCX>0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( M - 1 )*INCX
         END IF
         DO 40, J = 1, N
            IF( Y( JY )/=ZERO )THEN
               TEMP = ALPHA*Y( JY )
               IX   = KX
               DO 30, I = 1, M
                  A( I, J ) = A( I, J ) + X( IX )*TEMP
                  IX        = IX        + INCX
   30          CONTINUE
            END IF
            JY = JY + INCY
   40    CONTINUE
      END IF
!
      RETURN
!
!     END OF DGER  .
!
      END

! DSCAL
      SUBROUTINE DSCAL(N,DA,DX,INCX)
!
!  -- Reference BLAS level1 routine (version 3.8.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2017
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION DA
      INTEGER INCX,N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION DX(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER I,M,MP1,NINCX
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC mod
!     ..
      IF (n.LE.0 .OR. incx.LE.0) RETURN
      IF (incx.EQ.1) THEN
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
         m = mod(n,5)
         IF (m.NE.0) THEN
            DO i = 1,m
               dx(i) = da*dx(i)
            END DO
            IF (n.LT.5) RETURN
         END IF
         mp1 = m + 1
         DO i = mp1,n,5
            dx(i) = da*dx(i)
            dx(i+1) = da*dx(i+1)
            dx(i+2) = da*dx(i+2)
            dx(i+3) = da*dx(i+3)
            dx(i+4) = da*dx(i+4)
         END DO
      ELSE
!
!        code for increment not equal to 1
!
         nincx = n*incx
         DO i = 1,nincx,incx
            dx(i) = da*dx(i)
         END DO
      END IF
      RETURN
      END

! DSWAP
      SUBROUTINE DSWAP(N,DX,INCX,DY,INCY)
!
!  -- Reference BLAS level1 routine (version 3.8.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2017
!
!     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      DOUBLE PRECISION DTEMP
      INTEGER I,IX,IY,M,MP1
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC mod
!     ..
      IF (n.LE.0) RETURN
      IF (incx.EQ.1 .AND. incy.EQ.1) THEN
!
!       code for both increments equal to 1
!
!
!       clean-up loop
!
         m = mod(n,3)
         IF (m.NE.0) THEN
            DO i = 1,m
               dtemp = dx(i)
               dx(i) = dy(i)
               dy(i) = dtemp
            END DO
            IF (n.LT.3) RETURN
         END IF
         mp1 = m + 1
         DO i = mp1,n,3
            dtemp = dx(i)
            dx(i) = dy(i)
            dy(i) = dtemp
            dtemp = dx(i+1)
            dx(i+1) = dy(i+1)
            dy(i+1) = dtemp
            dtemp = dx(i+2)
            dx(i+2) = dy(i+2)
            dy(i+2) = dtemp
         END DO
      ELSE
!
!       code for unequal increments or equal increments not equal
!         to 1
!
         ix = 1
         iy = 1
         IF (incx.LT.0) ix = (-n+1)*incx + 1
         IF (incy.LT.0) iy = (-n+1)*incy + 1
         DO i = 1,n
            dtemp = dx(ix)
            dx(ix) = dy(iy)
            dy(iy) = dtemp
            ix = ix + incx
            iy = iy + incy
         END DO
      END IF
      RETURN
      END

! DSYEV
      SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
!
!  -- LAPACK driver routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, LDA, LWORK, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), W( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DSYEV computes all eigenvalues and, optionally, eigenvectors of a
!  real symmetric matrix A.
!
!  Arguments
!  =========
!
!  JOBZ    (input) CHARACTER*1
!          = 'N':  Compute eigenvalues only;
!          = 'V':  Compute eigenvalues and eigenvectors.
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
!          On entry, the symmetric matrix A.  If UPLO = 'U', the
!          leading N-by-N upper triangular part of A contains the
!          upper triangular part of the matrix A.  If UPLO = 'L',
!          the leading N-by-N lower triangular part of A contains
!          the lower triangular part of the matrix A.
!          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
!          orthonormal eigenvectors of the matrix A.
!          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
!          or the upper triangle (if UPLO='U') of A, including the
!          diagonal, is destroyed.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  W       (output) DOUBLE PRECISION array, dimension (N)
!          If INFO = 0, the eigenvalues in ascending order.
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The length of the array WORK.  LWORK >= max(1,3*N-1).
!          For optimal efficiency, LWORK >= (NB+2)*N,
!          where NB is the blocksize for DSYTRD returned by ILAENV.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, the algorithm failed to converge; i
!                off-diagonal elements of an intermediate tridiagonal
!                form did not converge to zero.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LOWER, LQUERY, WANTZ
      INTEGER            IINFO, IMAX, INDE, INDTAU, INDWRK, ISCALE,
     $                   LLWORK, LWKOPT, NB
      DOUBLE PRECISION   ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA,
     $                   SMLNUM
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH, DLANSY
      EXTERNAL           LSAME, ILAENV, DLAMCH, DLANSY
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLASCL, DORGTR, DSCAL, DSTEQR, DSTERF, DSYTRD,
     $                   XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      WANTZ = LSAME( JOBZ, 'V' )
      LOWER = LSAME( UPLO, 'L' )
      LQUERY = ( LWORK.EQ.-1 )
!
      INFO = 0
      IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( LOWER .OR. LSAME( UPLO, 'U' ) ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
!
      IF( INFO.EQ.0 ) THEN
         NB = ILAENV( 1, 'DSYTRD', UPLO, N, -1, -1, -1 )
         LWKOPT = MAX( 1, ( NB+2 )*N )
         WORK( 1 ) = LWKOPT
!
         IF( LWORK.LT.MAX( 1, 3*N-1 ) .AND. .NOT.LQUERY )
     $      INFO = -8
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYEV ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) THEN
         RETURN
      END IF
!
      IF( N.EQ.1 ) THEN
         W( 1 ) = A( 1, 1 )
         WORK( 1 ) = 2
         IF( WANTZ )
     $      A( 1, 1 ) = ONE
         RETURN
      END IF
!
!     Get machine constants.
!
      SAFMIN = DLAMCH( 'Safe minimum' )
      EPS = DLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN = SQRT( SMLNUM )
      RMAX = SQRT( BIGNUM )
!
!     Scale matrix to allowable range, if necessary.
!
      ANRM = DLANSY( 'M', UPLO, N, A, LDA, WORK )
      ISCALE = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) THEN
         ISCALE = 1
         SIGMA = RMIN / ANRM
      ELSE IF( ANRM.GT.RMAX ) THEN
         ISCALE = 1
         SIGMA = RMAX / ANRM
      END IF
      IF( ISCALE.EQ.1 )
     $   CALL DLASCL( UPLO, 0, 0, ONE, SIGMA, N, N, A, LDA, INFO )
!
!     Call DSYTRD to reduce symmetric matrix to tridiagonal form.
!
      INDE = 1
      INDTAU = INDE + N
      INDWRK = INDTAU + N
      LLWORK = LWORK - INDWRK + 1
      CALL DSYTRD( UPLO, N, A, LDA, W, WORK( INDE ), WORK( INDTAU ),
     $             WORK( INDWRK ), LLWORK, IINFO )
!
!     For eigenvalues only, call DSTERF.  For eigenvectors, first call
!     DORGTR to generate the orthogonal matrix, then call DSTEQR.
!
      IF( .NOT.WANTZ ) THEN
         CALL DSTERF( N, W, WORK( INDE ), INFO )
      ELSE
         CALL DORGTR( UPLO, N, A, LDA, WORK( INDTAU ), WORK( INDWRK ),
     $                LLWORK, IINFO )
         CALL DSTEQR( JOBZ, N, W, WORK( INDE ), A, LDA, WORK( INDTAU ),
     $                INFO )
      END IF
!
!     If matrix was scaled, then rescale eigenvalues appropriately.
!
      IF( ISCALE.EQ.1 ) THEN
         IF( INFO.EQ.0 ) THEN
            IMAX = N
         ELSE
            IMAX = INFO - 1
         END IF
         CALL DSCAL( IMAX, ONE / SIGMA, W, 1 )
      END IF
!
!     Set WORK(1) to optimal workspace size.
!
      WORK( 1 ) = LWKOPT
!
      RETURN
!
!     End of DSYEV
!
      END

! DLANSY
      DOUBLE PRECISION FUNCTION DLANSY( NORM, UPLO, N, A, LDA, WORK )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          NORM, UPLO
      INTEGER            LDA, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DLANSY  returns the value of the one norm,  or the Frobenius norm, or
*  the  infinity norm,  or the  element of  largest absolute value  of a
*  real symmetric matrix A.
*
*  Description
*  ===========
*
*  DLANSY returns the value
*
*     DLANSY = ( max(abs(A(i,j))), NORM = 'M' or 'm'
*              (
*              ( norm1(A),         NORM = '1', 'O' or 'o'
*              (
*              ( normI(A),         NORM = 'I' or 'i'
*              (
*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
*
*  where  norm1  denotes the  one norm of a matrix (maximum column sum),
*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
*  normF  denotes the  Frobenius norm of a matrix (square root of sum of
*  squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.
*
*  Arguments
*  =========
*
*  NORM    (input) CHARACTER*1
*          Specifies the value to be returned in DLANSY as described
*          above.
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the upper or lower triangular part of the
*          symmetric matrix A is to be referenced.
*          = 'U':  Upper triangular part of A is referenced
*          = 'L':  Lower triangular part of A is referenced
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.  When N = 0, DLANSY is
*          set to zero.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The symmetric matrix A.  If UPLO = 'U', the leading n by n
*          upper triangular part of A contains the upper triangular part
*          of the matrix A, and the strictly lower triangular part of A
*          is not referenced.  If UPLO = 'L', the leading n by n lower
*          triangular part of A contains the lower triangular part of
*          the matrix A, and the strictly upper triangular part of A is
*          not referenced.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(N,1).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (MAX(1,LWORK)),
*          where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,
*          WORK is not referenced.
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   ABSA, SCALE, SUM, VALUE
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLASSQ
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
*     ..
*     .. Executable Statements ..
*
      IF( N.EQ.0 ) THEN
         VALUE = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
*
*        Find max(abs(A(i,j))).
*
         VALUE = ZERO
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 20 J = 1, N
               DO 10 I = 1, J
                  VALUE = MAX( VALUE, ABS( A( I, J ) ) )
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40 J = 1, N
               DO 30 I = J, N
                  VALUE = MAX( VALUE, ABS( A( I, J ) ) )
   30          CONTINUE
   40       CONTINUE
         END IF
      ELSE IF( ( LSAME( NORM, 'I' ) ) .OR. ( LSAME( NORM, 'O' ) ) .OR.
     $         ( NORM.EQ.'1' ) ) THEN
*
*        Find normI(A) ( = norm1(A), since A is symmetric).
*
         VALUE = ZERO
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 60 J = 1, N
               SUM = ZERO
               DO 50 I = 1, J - 1
                  ABSA = ABS( A( I, J ) )
                  SUM = SUM + ABSA
                  WORK( I ) = WORK( I ) + ABSA
   50          CONTINUE
               WORK( J ) = SUM + ABS( A( J, J ) )
   60       CONTINUE
            DO 70 I = 1, N
               VALUE = MAX( VALUE, WORK( I ) )
   70       CONTINUE
         ELSE
            DO 80 I = 1, N
               WORK( I ) = ZERO
   80       CONTINUE
            DO 100 J = 1, N
               SUM = WORK( J ) + ABS( A( J, J ) )
               DO 90 I = J + 1, N
                  ABSA = ABS( A( I, J ) )
                  SUM = SUM + ABSA
                  WORK( I ) = WORK( I ) + ABSA
   90          CONTINUE
               VALUE = MAX( VALUE, SUM )
  100       CONTINUE
         END IF
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
*
*        Find normF(A).
*
         SCALE = ZERO
         SUM = ONE
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 110 J = 2, N
               CALL DLASSQ( J-1, A( 1, J ), 1, SCALE, SUM )
  110       CONTINUE
         ELSE
            DO 120 J = 1, N - 1
               CALL DLASSQ( N-J, A( J+1, J ), 1, SCALE, SUM )
  120       CONTINUE
         END IF
         SUM = 2*SUM
         CALL DLASSQ( N, A, LDA+1, SCALE, SUM )
         VALUE = SCALE*SQRT( SUM )
      END IF
*
      DLANSY = VALUE
      RETURN
*
*     End of DLANSY
*
      END

! DLASCL
      SUBROUTINE DLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          TYPE
      INTEGER            INFO, KL, KU, LDA, M, N
      DOUBLE PRECISION   CFROM, CTO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DLASCL multiplies the M by N real matrix A by the real scalar
*  CTO/CFROM.  This is done without over/underflow as long as the final
*  result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that
*  A may be full, upper triangular, lower triangular, upper Hessenberg,
*  or banded.
*
*  Arguments
*  =========
*
*  TYPE    (input) CHARACTER*1
*          TYPE indices the storage type of the input matrix.
*          = 'G':  A is a full matrix.
*          = 'L':  A is a lower triangular matrix.
*          = 'U':  A is an upper triangular matrix.
*          = 'H':  A is an upper Hessenberg matrix.
*          = 'B':  A is a symmetric band matrix with lower bandwidth KL
*                  and upper bandwidth KU and with the only the lower
*                  half stored.
*          = 'Q':  A is a symmetric band matrix with lower bandwidth KL
*                  and upper bandwidth KU and with the only the upper
*                  half stored.
*          = 'Z':  A is a band matrix with lower bandwidth KL and upper
*                  bandwidth KU.
*
*  KL      (input) INTEGER
*          The lower bandwidth of A.  Referenced only if TYPE = 'B',
*          'Q' or 'Z'.
*
*  KU      (input) INTEGER
*          The upper bandwidth of A.  Referenced only if TYPE = 'B',
*          'Q' or 'Z'.
*
*  CFROM   (input) DOUBLE PRECISION
*  CTO     (input) DOUBLE PRECISION
*          The matrix A is multiplied by CTO/CFROM. A(I,J) is computed
*          without over/underflow if the final result CTO*A(I,J)/CFROM
*          can be represented without over/underflow.  CFROM must be
*          nonzero.
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          The matrix to be multiplied by CTO/CFROM.  See TYPE for the
*          storage type.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  INFO    (output) INTEGER
*          0  - successful exit
*          <0 - if INFO = -i, the i-th argument had an illegal value.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            DONE
      INTEGER            I, ITYPE, J, K1, K2, K3, K4
      DOUBLE PRECISION   BIGNUM, CFROM1, CFROMC, CTO1, CTOC, MUL, SMLNUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, DLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
*
      IF( LSAME( TYPE, 'G' ) ) THEN
         ITYPE = 0
      ELSE IF( LSAME( TYPE, 'L' ) ) THEN
         ITYPE = 1
      ELSE IF( LSAME( TYPE, 'U' ) ) THEN
         ITYPE = 2
      ELSE IF( LSAME( TYPE, 'H' ) ) THEN
         ITYPE = 3
      ELSE IF( LSAME( TYPE, 'B' ) ) THEN
         ITYPE = 4
      ELSE IF( LSAME( TYPE, 'Q' ) ) THEN
         ITYPE = 5
      ELSE IF( LSAME( TYPE, 'Z' ) ) THEN
         ITYPE = 6
      ELSE
         ITYPE = -1
      END IF
*
      IF( ITYPE.EQ.-1 ) THEN
         INFO = -1
      ELSE IF( CFROM.EQ.ZERO ) THEN
         INFO = -4
      ELSE IF( M.LT.0 ) THEN
         INFO = -6
      ELSE IF( N.LT.0 .OR. ( ITYPE.EQ.4 .AND. N.NE.M ) .OR.
     $         ( ITYPE.EQ.5 .AND. N.NE.M ) ) THEN
         INFO = -7
      ELSE IF( ITYPE.LE.3 .AND. LDA.LT.MAX( 1, M ) ) THEN
         INFO = -9
      ELSE IF( ITYPE.GE.4 ) THEN
         IF( KL.LT.0 .OR. KL.GT.MAX( M-1, 0 ) ) THEN
            INFO = -2
         ELSE IF( KU.LT.0 .OR. KU.GT.MAX( N-1, 0 ) .OR.
     $            ( ( ITYPE.EQ.4 .OR. ITYPE.EQ.5 ) .AND. KL.NE.KU ) )
     $             THEN
            INFO = -3
         ELSE IF( ( ITYPE.EQ.4 .AND. LDA.LT.KL+1 ) .OR.
     $            ( ITYPE.EQ.5 .AND. LDA.LT.KU+1 ) .OR.
     $            ( ITYPE.EQ.6 .AND. LDA.LT.2*KL+KU+1 ) ) THEN
            INFO = -9
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLASCL', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. M.EQ.0 )
     $   RETURN
*
*     Get machine parameters
*
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
*
      CFROMC = CFROM
      CTOC = CTO
*
   10 CONTINUE
      CFROM1 = CFROMC*SMLNUM
      CTO1 = CTOC / BIGNUM
      IF( ABS( CFROM1 ).GT.ABS( CTOC ) .AND. CTOC.NE.ZERO ) THEN
         MUL = SMLNUM
         DONE = .FALSE.
         CFROMC = CFROM1
      ELSE IF( ABS( CTO1 ).GT.ABS( CFROMC ) ) THEN
         MUL = BIGNUM
         DONE = .FALSE.
         CTOC = CTO1
      ELSE
         MUL = CTOC / CFROMC
         DONE = .TRUE.
      END IF
*
      IF( ITYPE.EQ.0 ) THEN
*
*        Full matrix
*
         DO 30 J = 1, N
            DO 20 I = 1, M
               A( I, J ) = A( I, J )*MUL
   20       CONTINUE
   30    CONTINUE
*
      ELSE IF( ITYPE.EQ.1 ) THEN
*
*        Lower triangular matrix
*
         DO 50 J = 1, N
            DO 40 I = J, M
               A( I, J ) = A( I, J )*MUL
   40       CONTINUE
   50    CONTINUE
*
      ELSE IF( ITYPE.EQ.2 ) THEN
*
*        Upper triangular matrix
*
         DO 70 J = 1, N
            DO 60 I = 1, MIN( J, M )
               A( I, J ) = A( I, J )*MUL
   60       CONTINUE
   70    CONTINUE
*
      ELSE IF( ITYPE.EQ.3 ) THEN
*
*        Upper Hessenberg matrix
*
         DO 90 J = 1, N
            DO 80 I = 1, MIN( J+1, M )
               A( I, J ) = A( I, J )*MUL
   80       CONTINUE
   90    CONTINUE
*
      ELSE IF( ITYPE.EQ.4 ) THEN
*
*        Lower half of a symmetric band matrix
*
         K3 = KL + 1
         K4 = N + 1
         DO 110 J = 1, N
            DO 100 I = 1, MIN( K3, K4-J )
               A( I, J ) = A( I, J )*MUL
  100       CONTINUE
  110    CONTINUE
*
      ELSE IF( ITYPE.EQ.5 ) THEN
*
*        Upper half of a symmetric band matrix
*
         K1 = KU + 2
         K3 = KU + 1
         DO 130 J = 1, N
            DO 120 I = MAX( K1-J, 1 ), K3
               A( I, J ) = A( I, J )*MUL
  120       CONTINUE
  130    CONTINUE
*
      ELSE IF( ITYPE.EQ.6 ) THEN
*
*        Band matrix
*
         K1 = KL + KU + 2
         K2 = KL + 1
         K3 = 2*KL + KU + 1
         K4 = KL + KU + 1 + M
         DO 150 J = 1, N
            DO 140 I = MAX( K1-J, K2 ), MIN( K3, K4-J )
               A( I, J ) = A( I, J )*MUL
  140       CONTINUE
  150    CONTINUE
*
      END IF
*
      IF( .NOT.DONE )
     $   GO TO 10
*
      RETURN
*
*     End of DLASCL
*
      END

! DSYTRD
      SUBROUTINE DSYTRD( UPLO, N, A, LDA, D, E, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), D( * ), E( * ), TAU( * ),
     $                   WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DSYTRD reduces a real symmetric matrix A to real symmetric
*  tridiagonal form T by an orthogonal similarity transformation:
*  Q**T * A * Q = T.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
*          N-by-N upper triangular part of A contains the upper
*          triangular part of the matrix A, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading N-by-N lower triangular part of A contains the lower
*          triangular part of the matrix A, and the strictly upper
*          triangular part of A is not referenced.
*          On exit, if UPLO = 'U', the diagonal and first superdiagonal
*          of A are overwritten by the corresponding elements of the
*          tridiagonal matrix T, and the elements above the first
*          superdiagonal, with the array TAU, represent the orthogonal
*          matrix Q as a product of elementary reflectors; if UPLO
*          = 'L', the diagonal and first subdiagonal of A are over-
*          written by the corresponding elements of the tridiagonal
*          matrix T, and the elements below the first subdiagonal, with
*          the array TAU, represent the orthogonal matrix Q as a product
*          of elementary reflectors. See Further Details.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  D       (output) DOUBLE PRECISION array, dimension (N)
*          The diagonal elements of the tridiagonal matrix T:
*          D(i) = A(i,i).
*
*  E       (output) DOUBLE PRECISION array, dimension (N-1)
*          The off-diagonal elements of the tridiagonal matrix T:
*          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.
*
*  TAU     (output) DOUBLE PRECISION array, dimension (N-1)
*          The scalar factors of the elementary reflectors (see Further
*          Details).
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.  LWORK >= 1.
*          For optimum performance LWORK >= N*NB, where NB is the
*          optimal blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  Further Details
*  ===============
*
*  If UPLO = 'U', the matrix Q is represented as a product of elementary
*  reflectors
*
*     Q = H(n-1) . . . H(2) H(1).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a real scalar, and v is a real vector with
*  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
*  A(1:i-1,i+1), and tau in TAU(i).
*
*  If UPLO = 'L', the matrix Q is represented as a product of elementary
*  reflectors
*
*     Q = H(1) H(2) . . . H(n-1).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a real scalar, and v is a real vector with
*  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
*  and tau in TAU(i).
*
*  The contents of A on exit are illustrated by the following examples
*  with n = 5:
*
*  if UPLO = 'U':                       if UPLO = 'L':
*
*    (  d   e   v2  v3  v4 )              (  d                  )
*    (      d   e   v3  v4 )              (  e   d              )
*    (          d   e   v4 )              (  v1  e   d          )
*    (              d   e  )              (  v1  v2  e   d      )
*    (                  d  )              (  v1  v2  v3  e   d  )
*
*  where d and e denote diagonal and off-diagonal elements of T, and vi
*  denotes an element of the vector defining H(i).
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, UPPER
      INTEGER            I, IINFO, IWS, J, KK, LDWORK, LWKOPT, NB,
     $                   NBMIN, NX
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLATRD, DSYR2K, DSYTD2, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK.EQ.-1 )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.1 .AND. .NOT.LQUERY ) THEN
         INFO = -9
      END IF
*
      IF( INFO.EQ.0 ) THEN
*
*        Determine the block size.
*
         NB = ILAENV( 1, 'DSYTRD', UPLO, N, -1, -1, -1 )
         LWKOPT = N*NB
         WORK( 1 ) = LWKOPT
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYTRD', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      NX = N
      IWS = 1
      IF( NB.GT.1 .AND. NB.LT.N ) THEN
*
*        Determine when to cross over from blocked to unblocked code
*        (last block is always handled by unblocked code).
*
         NX = MAX( NB, ILAENV( 3, 'DSYTRD', UPLO, N, -1, -1, -1 ) )
         IF( NX.LT.N ) THEN
*
*           Determine if workspace is large enough for blocked code.
*
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
*
*              Not enough workspace to use optimal NB:  determine the
*              minimum value of NB, and reduce NB or force use of
*              unblocked code by setting NX = N.
*
               NB = MAX( LWORK / LDWORK, 1 )
               NBMIN = ILAENV( 2, 'DSYTRD', UPLO, N, -1, -1, -1 )
               IF( NB.LT.NBMIN )
     $            NX = N
            END IF
         ELSE
            NX = N
         END IF
      ELSE
         NB = 1
      END IF
*
      IF( UPPER ) THEN
*
*        Reduce the upper triangle of A.
*        Columns 1:kk are handled by the unblocked method.
*
         KK = N - ( ( N-NX+NB-1 ) / NB )*NB
         DO 20 I = N - NB + 1, KK + 1, -NB
*
*           Reduce columns i:i+nb-1 to tridiagonal form and form the
*           matrix W which is needed to update the unreduced part of
*           the matrix
*
            CALL DLATRD( UPLO, I+NB-1, NB, A, LDA, E, TAU, WORK,
     $                   LDWORK )
*
*           Update the unreduced submatrix A(1:i-1,1:i-1), using an
*           update of the form:  A := A - V*W' - W*V'
*
            CALL DSYR2K( UPLO, 'No transpose', I-1, NB, -ONE, A( 1, I ),
     $                   LDA, WORK, LDWORK, ONE, A, LDA )
*
*           Copy superdiagonal elements back into A, and diagonal
*           elements into D
*
            DO 10 J = I, I + NB - 1
               A( J-1, J ) = E( J-1 )
               D( J ) = A( J, J )
   10       CONTINUE
   20    CONTINUE
*
*        Use unblocked code to reduce the last or only block
*
         CALL DSYTD2( UPLO, KK, A, LDA, D, E, TAU, IINFO )
      ELSE
*
*        Reduce the lower triangle of A
*
         DO 40 I = 1, N - NX, NB
*
*           Reduce columns i:i+nb-1 to tridiagonal form and form the
*           matrix W which is needed to update the unreduced part of
*           the matrix
*
            CALL DLATRD( UPLO, N-I+1, NB, A( I, I ), LDA, E( I ),
     $                   TAU( I ), WORK, LDWORK )
*
*           Update the unreduced submatrix A(i+ib:n,i+ib:n), using
*           an update of the form:  A := A - V*W' - W*V'
*
            CALL DSYR2K( UPLO, 'No transpose', N-I-NB+1, NB, -ONE,
     $                   A( I+NB, I ), LDA, WORK( NB+1 ), LDWORK, ONE,
     $                   A( I+NB, I+NB ), LDA )
*
*           Copy subdiagonal elements back into A, and diagonal
*           elements into D
*
            DO 30 J = I, I + NB - 1
               A( J+1, J ) = E( J )
               D( J ) = A( J, J )
   30       CONTINUE
   40    CONTINUE
*
*        Use unblocked code to reduce the last or only block
*
         CALL DSYTD2( UPLO, N-I+1, A( I, I ), LDA, D( I ), E( I ),
     $                TAU( I ), IINFO )
      END IF
*
      WORK( 1 ) = LWKOPT
      RETURN
*
*     End of DSYTRD
*
      END

! DSTERF
      SUBROUTINE DSTERF( N, D, E, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * )
*     ..
*
*  Purpose
*  =======
*
*  DSTERF computes all eigenvalues of a symmetric tridiagonal matrix
*  using the Pal-Walker-Kahan variant of the QL or QR algorithm.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix.  N >= 0.
*
*  D       (input/output) DOUBLE PRECISION array, dimension (N)
*          On entry, the n diagonal elements of the tridiagonal matrix.
*          On exit, if INFO = 0, the eigenvalues in ascending order.
*
*  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
*          On entry, the (n-1) subdiagonal elements of the tridiagonal
*          matrix.
*          On exit, E has been destroyed.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  the algorithm failed to find all of the eigenvalues in
*                a total of 30*N iterations; if INFO = i, then i
*                elements of E have not converged to zero.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO, THREE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,
     $                   THREE = 3.0D0 )
      INTEGER            MAXIT
      PARAMETER          ( MAXIT = 30 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, ISCALE, JTOT, L, L1, LEND, LENDSV, LSV, M,
     $                   NMAXIT
      DOUBLE PRECISION   ALPHA, ANORM, BB, C, EPS, EPS2, GAMMA, OLDC,
     $                   OLDGAM, P, R, RT1, RT2, RTE, S, SAFMAX, SAFMIN,
     $                   SIGMA, SSFMAX, SSFMIN
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLANST, DLAPY2
      EXTERNAL           DLAMCH, DLANST, DLAPY2
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAE2, DLASCL, DLASRT, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*
*     Quick return if possible
*
      IF( N.LT.0 ) THEN
         INFO = -1
         CALL XERBLA( 'DSTERF', -INFO )
         RETURN
      END IF
      IF( N.LE.1 )
     $   RETURN
*
*     Determine the unit roundoff for this environment.
*
      EPS = DLAMCH( 'E' )
      EPS2 = EPS**2
      SAFMIN = DLAMCH( 'S' )
      SAFMAX = ONE / SAFMIN
      SSFMAX = SQRT( SAFMAX ) / THREE
      SSFMIN = SQRT( SAFMIN ) / EPS2
*
*     Compute the eigenvalues of the tridiagonal matrix.
*
      NMAXIT = N*MAXIT
      SIGMA = ZERO
      JTOT = 0
*
*     Determine where the matrix splits and choose QL or QR iteration
*     for each block, according to whether top or bottom diagonal
*     element is smaller.
*
      L1 = 1
*
   10 CONTINUE
      IF( L1.GT.N )
     $   GO TO 170
      IF( L1.GT.1 )
     $   E( L1-1 ) = ZERO
      DO 20 M = L1, N - 1
         IF( ABS( E( M ) ).LE.( SQRT( ABS( D( M ) ) )*SQRT( ABS( D( M+
     $       1 ) ) ) )*EPS ) THEN
            E( M ) = ZERO
            GO TO 30
         END IF
   20 CONTINUE
      M = N
*
   30 CONTINUE
      L = L1
      LSV = L
      LEND = M
      LENDSV = LEND
      L1 = M + 1
      IF( LEND.EQ.L )
     $   GO TO 10
*
*     Scale submatrix in rows and columns L to LEND
*
      ANORM = DLANST( 'I', LEND-L+1, D( L ), E( L ) )
      ISCALE = 0
      IF( ANORM.GT.SSFMAX ) THEN
         ISCALE = 1
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L+1, 1, D( L ), N,
     $                INFO )
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L, 1, E( L ), N,
     $                INFO )
      ELSE IF( ANORM.LT.SSFMIN ) THEN
         ISCALE = 2
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L+1, 1, D( L ), N,
     $                INFO )
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L, 1, E( L ), N,
     $                INFO )
      END IF
*
      DO 40 I = L, LEND - 1
         E( I ) = E( I )**2
   40 CONTINUE
*
*     Choose between QL and QR iteration
*
      IF( ABS( D( LEND ) ).LT.ABS( D( L ) ) ) THEN
         LEND = LSV
         L = LENDSV
      END IF
*
      IF( LEND.GE.L ) THEN
*
*        QL Iteration
*
*        Look for small subdiagonal element.
*
   50    CONTINUE
         IF( L.NE.LEND ) THEN
            DO 60 M = L, LEND - 1
               IF( ABS( E( M ) ).LE.EPS2*ABS( D( M )*D( M+1 ) ) )
     $            GO TO 70
   60       CONTINUE
         END IF
         M = LEND
*
   70    CONTINUE
         IF( M.LT.LEND )
     $      E( M ) = ZERO
         P = D( L )
         IF( M.EQ.L )
     $      GO TO 90
*
*        If remaining matrix is 2 by 2, use DLAE2 to compute its
*        eigenvalues.
*
         IF( M.EQ.L+1 ) THEN
            RTE = SQRT( E( L ) )
            CALL DLAE2( D( L ), RTE, D( L+1 ), RT1, RT2 )
            D( L ) = RT1
            D( L+1 ) = RT2
            E( L ) = ZERO
            L = L + 2
            IF( L.LE.LEND )
     $         GO TO 50
            GO TO 150
         END IF
*
         IF( JTOT.EQ.NMAXIT )
     $      GO TO 150
         JTOT = JTOT + 1
*
*        Form shift.
*
         RTE = SQRT( E( L ) )
         SIGMA = ( D( L+1 )-P ) / ( TWO*RTE )
         R = DLAPY2( SIGMA, ONE )
         SIGMA = P - ( RTE / ( SIGMA+SIGN( R, SIGMA ) ) )
*
         C = ONE
         S = ZERO
         GAMMA = D( M ) - SIGMA
         P = GAMMA*GAMMA
*
*        Inner loop
*
         DO 80 I = M - 1, L, -1
            BB = E( I )
            R = P + BB
            IF( I.NE.M-1 )
     $         E( I+1 ) = S*R
            OLDC = C
            C = P / R
            S = BB / R
            OLDGAM = GAMMA
            ALPHA = D( I )
            GAMMA = C*( ALPHA-SIGMA ) - S*OLDGAM
            D( I+1 ) = OLDGAM + ( ALPHA-GAMMA )
            IF( C.NE.ZERO ) THEN
               P = ( GAMMA*GAMMA ) / C
            ELSE
               P = OLDC*BB
            END IF
   80    CONTINUE
*
         E( L ) = S*P
         D( L ) = SIGMA + GAMMA
         GO TO 50
*
*        Eigenvalue found.
*
   90    CONTINUE
         D( L ) = P
*
         L = L + 1
         IF( L.LE.LEND )
     $      GO TO 50
         GO TO 150
*
      ELSE
*
*        QR Iteration
*
*        Look for small superdiagonal element.
*
  100    CONTINUE
         DO 110 M = L, LEND + 1, -1
            IF( ABS( E( M-1 ) ).LE.EPS2*ABS( D( M )*D( M-1 ) ) )
     $         GO TO 120
  110    CONTINUE
         M = LEND
*
  120    CONTINUE
         IF( M.GT.LEND )
     $      E( M-1 ) = ZERO
         P = D( L )
         IF( M.EQ.L )
     $      GO TO 140
*
*        If remaining matrix is 2 by 2, use DLAE2 to compute its
*        eigenvalues.
*
         IF( M.EQ.L-1 ) THEN
            RTE = SQRT( E( L-1 ) )
            CALL DLAE2( D( L ), RTE, D( L-1 ), RT1, RT2 )
            D( L ) = RT1
            D( L-1 ) = RT2
            E( L-1 ) = ZERO
            L = L - 2
            IF( L.GE.LEND )
     $         GO TO 100
            GO TO 150
         END IF
*
         IF( JTOT.EQ.NMAXIT )
     $      GO TO 150
         JTOT = JTOT + 1
*
*        Form shift.
*
         RTE = SQRT( E( L-1 ) )
         SIGMA = ( D( L-1 )-P ) / ( TWO*RTE )
         R = DLAPY2( SIGMA, ONE )
         SIGMA = P - ( RTE / ( SIGMA+SIGN( R, SIGMA ) ) )
*
         C = ONE
         S = ZERO
         GAMMA = D( M ) - SIGMA
         P = GAMMA*GAMMA
*
*        Inner loop
*
         DO 130 I = M, L - 1
            BB = E( I )
            R = P + BB
            IF( I.NE.M )
     $         E( I-1 ) = S*R
            OLDC = C
            C = P / R
            S = BB / R
            OLDGAM = GAMMA
            ALPHA = D( I+1 )
            GAMMA = C*( ALPHA-SIGMA ) - S*OLDGAM
            D( I ) = OLDGAM + ( ALPHA-GAMMA )
            IF( C.NE.ZERO ) THEN
               P = ( GAMMA*GAMMA ) / C
            ELSE
               P = OLDC*BB
            END IF
  130    CONTINUE
*
         E( L-1 ) = S*P
         D( L ) = SIGMA + GAMMA
         GO TO 100
*
*        Eigenvalue found.
*
  140    CONTINUE
         D( L ) = P
*
         L = L - 1
         IF( L.GE.LEND )
     $      GO TO 100
         GO TO 150
*
      END IF
*
*     Undo scaling if necessary
*
  150 CONTINUE
      IF( ISCALE.EQ.1 )
     $   CALL DLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV+1, 1,
     $                D( LSV ), N, INFO )
      IF( ISCALE.EQ.2 )
     $   CALL DLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV+1, 1,
     $                D( LSV ), N, INFO )
*
*     Check for no convergence to an eigenvalue after a total
*     of N*MAXIT iterations.
*
      IF( JTOT.LT.NMAXIT )
     $   GO TO 10
      DO 160 I = 1, N - 1
         IF( E( I ).NE.ZERO )
     $      INFO = INFO + 1
  160 CONTINUE
      GO TO 180
*
*     Sort eigenvalues in increasing order.
*
  170 CONTINUE
      CALL DLASRT( 'I', N, D, INFO )
*
  180 CONTINUE
      RETURN
*
*     End of DSTERF
*
      END

! DORGTR
      SUBROUTINE DORGTR( UPLO, N, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DORGTR generates a real orthogonal matrix Q which is defined as the
*  product of n-1 elementary reflectors of order N, as returned by
*  DSYTRD:
*
*  if UPLO = 'U', Q = H(n-1) . . . H(2) H(1),
*
*  if UPLO = 'L', Q = H(1) H(2) . . . H(n-1).
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U': Upper triangle of A contains elementary reflectors
*                 from DSYTRD;
*          = 'L': Lower triangle of A contains elementary reflectors
*                 from DSYTRD.
*
*  N       (input) INTEGER
*          The order of the matrix Q. N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the vectors which define the elementary reflectors,
*          as returned by DSYTRD.
*          On exit, the N-by-N orthogonal matrix Q.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= max(1,N).
*
*  TAU     (input) DOUBLE PRECISION array, dimension (N-1)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by DSYTRD.
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK. LWORK >= max(1,N-1).
*          For optimum performance LWORK >= (N-1)*NB, where NB is
*          the optimal blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, UPPER
      INTEGER            I, IINFO, J, LWKOPT, NB
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           DORGQL, DORGQR, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.MAX( 1, N-1 ) .AND. .NOT.LQUERY ) THEN
         INFO = -7
      END IF
*
      IF( INFO.EQ.0 ) THEN
         IF( UPPER ) THEN
            NB = ILAENV( 1, 'DORGQL', ' ', N-1, N-1, N-1, -1 )
         ELSE
            NB = ILAENV( 1, 'DORGQR', ' ', N-1, N-1, N-1, -1 )
         END IF
         LWKOPT = MAX( 1, N-1 )*NB
         WORK( 1 ) = LWKOPT
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORGTR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      IF( UPPER ) THEN
*
*        Q was determined by a call to DSYTRD with UPLO = 'U'
*
*        Shift the vectors which define the elementary reflectors one
*        column to the left, and set the last row and column of Q to
*        those of the unit matrix
*
         DO 20 J = 1, N - 1
            DO 10 I = 1, J - 1
               A( I, J ) = A( I, J+1 )
   10       CONTINUE
            A( N, J ) = ZERO
   20    CONTINUE
         DO 30 I = 1, N - 1
            A( I, N ) = ZERO
   30    CONTINUE
         A( N, N ) = ONE
*
*        Generate Q(1:n-1,1:n-1)
*
         CALL DORGQL( N-1, N-1, N-1, A, LDA, TAU, WORK, LWORK, IINFO )
*
      ELSE
*
*        Q was determined by a call to DSYTRD with UPLO = 'L'.
*
*        Shift the vectors which define the elementary reflectors one
*        column to the right, and set the first row and column of Q to
*        those of the unit matrix
*
         DO 50 J = N, 2, -1
            A( 1, J ) = ZERO
            DO 40 I = J + 1, N
               A( I, J ) = A( I, J-1 )
   40       CONTINUE
   50    CONTINUE
         A( 1, 1 ) = ONE
         DO 60 I = 2, N
            A( I, 1 ) = ZERO
   60    CONTINUE
         IF( N.GT.1 ) THEN
*
*           Generate Q(2:n,2:n)
*
            CALL DORGQR( N-1, N-1, N-1, A( 2, 2 ), LDA, TAU, WORK,
     $                   LWORK, IINFO )
         END IF
      END IF
      WORK( 1 ) = LWKOPT
      RETURN
*
*     End of DORGTR
*
      END

! DSTEQR
      SUBROUTINE DSTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          COMPZ
      INTEGER            INFO, LDZ, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * ), WORK( * ), Z( LDZ, * )
*     ..
*
*  Purpose
*  =======
*
*  DSTEQR computes all eigenvalues and, optionally, eigenvectors of a
*  symmetric tridiagonal matrix using the implicit QL or QR method.
*  The eigenvectors of a full or band symmetric matrix can also be found
*  if DSYTRD or DSPTRD or DSBTRD has been used to reduce this matrix to
*  tridiagonal form.
*
*  Arguments
*  =========
*
*  COMPZ   (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only.
*          = 'V':  Compute eigenvalues and eigenvectors of the original
*                  symmetric matrix.  On entry, Z must contain the
*                  orthogonal matrix used to reduce the original matrix
*                  to tridiagonal form.
*          = 'I':  Compute eigenvalues and eigenvectors of the
*                  tridiagonal matrix.  Z is initialized to the identity
*                  matrix.
*
*  N       (input) INTEGER
*          The order of the matrix.  N >= 0.
*
*  D       (input/output) DOUBLE PRECISION array, dimension (N)
*          On entry, the diagonal elements of the tridiagonal matrix.
*          On exit, if INFO = 0, the eigenvalues in ascending order.
*
*  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
*          On entry, the (n-1) subdiagonal elements of the tridiagonal
*          matrix.
*          On exit, E has been destroyed.
*
*  Z       (input/output) DOUBLE PRECISION array, dimension (LDZ, N)
*          On entry, if  COMPZ = 'V', then Z contains the orthogonal
*          matrix used in the reduction to tridiagonal form.
*          On exit, if INFO = 0, then if  COMPZ = 'V', Z contains the
*          orthonormal eigenvectors of the original symmetric matrix,
*          and if COMPZ = 'I', Z contains the orthonormal eigenvectors
*          of the symmetric tridiagonal matrix.
*          If COMPZ = 'N', then Z is not referenced.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ >= 1, and if
*          eigenvectors are desired, then  LDZ >= max(1,N).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (max(1,2*N-2))
*          If COMPZ = 'N', then WORK is not referenced.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  the algorithm has failed to find all the eigenvalues in
*                a total of 30*N iterations; if INFO = i, then i
*                elements of E have not converged to zero; on exit, D
*                and E contain the elements of a symmetric tridiagonal
*                matrix which is orthogonally similar to the original
*                matrix.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO, THREE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,
     $                   THREE = 3.0D0 )
      INTEGER            MAXIT
      PARAMETER          ( MAXIT = 30 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, ICOMPZ, II, ISCALE, J, JTOT, K, L, L1, LEND,
     $                   LENDM1, LENDP1, LENDSV, LM1, LSV, M, MM, MM1,
     $                   NM1, NMAXIT
      DOUBLE PRECISION   ANORM, B, C, EPS, EPS2, F, G, P, R, RT1, RT2,
     $                   S, SAFMAX, SAFMIN, SSFMAX, SSFMIN, TST
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANST, DLAPY2
      EXTERNAL           LSAME, DLAMCH, DLANST, DLAPY2
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAE2, DLAEV2, DLARTG, DLASCL, DLASET, DLASR,
     $                   DLASRT, DSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SIGN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*
      IF( LSAME( COMPZ, 'N' ) ) THEN
         ICOMPZ = 0
      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
         ICOMPZ = 1
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ICOMPZ = 2
      ELSE
         ICOMPZ = -1
      END IF
      IF( ICOMPZ.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( ( LDZ.LT.1 ) .OR. ( ICOMPZ.GT.0 .AND. LDZ.LT.MAX( 1,
     $         N ) ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSTEQR', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( N.EQ.1 ) THEN
         IF( ICOMPZ.EQ.2 )
     $      Z( 1, 1 ) = ONE
         RETURN
      END IF
*
*     Determine the unit roundoff and over/underflow thresholds.
*
      EPS = DLAMCH( 'E' )
      EPS2 = EPS**2
      SAFMIN = DLAMCH( 'S' )
      SAFMAX = ONE / SAFMIN
      SSFMAX = SQRT( SAFMAX ) / THREE
      SSFMIN = SQRT( SAFMIN ) / EPS2
*
*     Compute the eigenvalues and eigenvectors of the tridiagonal
*     matrix.
*
      IF( ICOMPZ.EQ.2 )
     $   CALL DLASET( 'Full', N, N, ZERO, ONE, Z, LDZ )
*
      NMAXIT = N*MAXIT
      JTOT = 0
*
*     Determine where the matrix splits and choose QL or QR iteration
*     for each block, according to whether top or bottom diagonal
*     element is smaller.
*
      L1 = 1
      NM1 = N - 1
*
   10 CONTINUE
      IF( L1.GT.N )
     $   GO TO 160
      IF( L1.GT.1 )
     $   E( L1-1 ) = ZERO
      IF( L1.LE.NM1 ) THEN
         DO 20 M = L1, NM1
            TST = ABS( E( M ) )
            IF( TST.EQ.ZERO )
     $         GO TO 30
            IF( TST.LE.( SQRT( ABS( D( M ) ) )*SQRT( ABS( D( M+
     $          1 ) ) ) )*EPS ) THEN
               E( M ) = ZERO
               GO TO 30
            END IF
   20    CONTINUE
      END IF
      M = N
*
   30 CONTINUE
      L = L1
      LSV = L
      LEND = M
      LENDSV = LEND
      L1 = M + 1
      IF( LEND.EQ.L )
     $   GO TO 10
*
*     Scale submatrix in rows and columns L to LEND
*
      ANORM = DLANST( 'I', LEND-L+1, D( L ), E( L ) )
      ISCALE = 0
      IF( ANORM.EQ.ZERO )
     $   GO TO 10
      IF( ANORM.GT.SSFMAX ) THEN
         ISCALE = 1
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L+1, 1, D( L ), N,
     $                INFO )
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L, 1, E( L ), N,
     $                INFO )
      ELSE IF( ANORM.LT.SSFMIN ) THEN
         ISCALE = 2
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L+1, 1, D( L ), N,
     $                INFO )
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L, 1, E( L ), N,
     $                INFO )
      END IF
*
*     Choose between QL and QR iteration
*
      IF( ABS( D( LEND ) ).LT.ABS( D( L ) ) ) THEN
         LEND = LSV
         L = LENDSV
      END IF
*
      IF( LEND.GT.L ) THEN
*
*        QL Iteration
*
*        Look for small subdiagonal element.
*
   40    CONTINUE
         IF( L.NE.LEND ) THEN
            LENDM1 = LEND - 1
            DO 50 M = L, LENDM1
               TST = ABS( E( M ) )**2
               IF( TST.LE.( EPS2*ABS( D( M ) ) )*ABS( D( M+1 ) )+
     $             SAFMIN )GO TO 60
   50       CONTINUE
         END IF
*
         M = LEND
*
   60    CONTINUE
         IF( M.LT.LEND )
     $      E( M ) = ZERO
         P = D( L )
         IF( M.EQ.L )
     $      GO TO 80
*
*        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
*        to compute its eigensystem.
*
         IF( M.EQ.L+1 ) THEN
            IF( ICOMPZ.GT.0 ) THEN
               CALL DLAEV2( D( L ), E( L ), D( L+1 ), RT1, RT2, C, S )
               WORK( L ) = C
               WORK( N-1+L ) = S
               CALL DLASR( 'R', 'V', 'B', N, 2, WORK( L ),
     $                     WORK( N-1+L ), Z( 1, L ), LDZ )
            ELSE
               CALL DLAE2( D( L ), E( L ), D( L+1 ), RT1, RT2 )
            END IF
            D( L ) = RT1
            D( L+1 ) = RT2
            E( L ) = ZERO
            L = L + 2
            IF( L.LE.LEND )
     $         GO TO 40
            GO TO 140
         END IF
*
         IF( JTOT.EQ.NMAXIT )
     $      GO TO 140
         JTOT = JTOT + 1
*
*        Form shift.
*
         G = ( D( L+1 )-P ) / ( TWO*E( L ) )
         R = DLAPY2( G, ONE )
         G = D( M ) - P + ( E( L ) / ( G+SIGN( R, G ) ) )
*
         S = ONE
         C = ONE
         P = ZERO
*
*        Inner loop
*
         MM1 = M - 1
         DO 70 I = MM1, L, -1
            F = S*E( I )
            B = C*E( I )
            CALL DLARTG( G, F, C, S, R )
            IF( I.NE.M-1 )
     $         E( I+1 ) = R
            G = D( I+1 ) - P
            R = ( D( I )-G )*S + TWO*C*B
            P = S*R
            D( I+1 ) = G + P
            G = C*R - B
*
*           If eigenvectors are desired, then save rotations.
*
            IF( ICOMPZ.GT.0 ) THEN
               WORK( I ) = C
               WORK( N-1+I ) = -S
            END IF
*
   70    CONTINUE
*
*        If eigenvectors are desired, then apply saved rotations.
*
         IF( ICOMPZ.GT.0 ) THEN
            MM = M - L + 1
            CALL DLASR( 'R', 'V', 'B', N, MM, WORK( L ), WORK( N-1+L ),
     $                  Z( 1, L ), LDZ )
         END IF
*
         D( L ) = D( L ) - P
         E( L ) = G
         GO TO 40
*
*        Eigenvalue found.
*
   80    CONTINUE
         D( L ) = P
*
         L = L + 1
         IF( L.LE.LEND )
     $      GO TO 40
         GO TO 140
*
      ELSE
*
*        QR Iteration
*
*        Look for small superdiagonal element.
*
   90    CONTINUE
         IF( L.NE.LEND ) THEN
            LENDP1 = LEND + 1
            DO 100 M = L, LENDP1, -1
               TST = ABS( E( M-1 ) )**2
               IF( TST.LE.( EPS2*ABS( D( M ) ) )*ABS( D( M-1 ) )+
     $             SAFMIN )GO TO 110
  100       CONTINUE
         END IF
*
         M = LEND
*
  110    CONTINUE
         IF( M.GT.LEND )
     $      E( M-1 ) = ZERO
         P = D( L )
         IF( M.EQ.L )
     $      GO TO 130
*
*        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
*        to compute its eigensystem.
*
         IF( M.EQ.L-1 ) THEN
            IF( ICOMPZ.GT.0 ) THEN
               CALL DLAEV2( D( L-1 ), E( L-1 ), D( L ), RT1, RT2, C, S )
               WORK( M ) = C
               WORK( N-1+M ) = S
               CALL DLASR( 'R', 'V', 'F', N, 2, WORK( M ),
     $                     WORK( N-1+M ), Z( 1, L-1 ), LDZ )
            ELSE
               CALL DLAE2( D( L-1 ), E( L-1 ), D( L ), RT1, RT2 )
            END IF
            D( L-1 ) = RT1
            D( L ) = RT2
            E( L-1 ) = ZERO
            L = L - 2
            IF( L.GE.LEND )
     $         GO TO 90
            GO TO 140
         END IF
*
         IF( JTOT.EQ.NMAXIT )
     $      GO TO 140
         JTOT = JTOT + 1
*
*        Form shift.
*
         G = ( D( L-1 )-P ) / ( TWO*E( L-1 ) )
         R = DLAPY2( G, ONE )
         G = D( M ) - P + ( E( L-1 ) / ( G+SIGN( R, G ) ) )
*
         S = ONE
         C = ONE
         P = ZERO
*
*        Inner loop
*
         LM1 = L - 1
         DO 120 I = M, LM1
            F = S*E( I )
            B = C*E( I )
            CALL DLARTG( G, F, C, S, R )
            IF( I.NE.M )
     $         E( I-1 ) = R
            G = D( I ) - P
            R = ( D( I+1 )-G )*S + TWO*C*B
            P = S*R
            D( I ) = G + P
            G = C*R - B
*
*           If eigenvectors are desired, then save rotations.
*
            IF( ICOMPZ.GT.0 ) THEN
               WORK( I ) = C
               WORK( N-1+I ) = S
            END IF
*
  120    CONTINUE
*
*        If eigenvectors are desired, then apply saved rotations.
*
         IF( ICOMPZ.GT.0 ) THEN
            MM = L - M + 1
            CALL DLASR( 'R', 'V', 'F', N, MM, WORK( M ), WORK( N-1+M ),
     $                  Z( 1, M ), LDZ )
         END IF
*
         D( L ) = D( L ) - P
         E( LM1 ) = G
         GO TO 90
*
*        Eigenvalue found.
*
  130    CONTINUE
         D( L ) = P
*
         L = L - 1
         IF( L.GE.LEND )
     $      GO TO 90
         GO TO 140
*
      END IF
*
*     Undo scaling if necessary
*
  140 CONTINUE
      IF( ISCALE.EQ.1 ) THEN
         CALL DLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV+1, 1,
     $                D( LSV ), N, INFO )
         CALL DLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV, 1, E( LSV ),
     $                N, INFO )
      ELSE IF( ISCALE.EQ.2 ) THEN
         CALL DLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV+1, 1,
     $                D( LSV ), N, INFO )
         CALL DLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV, 1, E( LSV ),
     $                N, INFO )
      END IF
*
*     Check for no convergence to an eigenvalue after a total
*     of N*MAXIT iterations.
*
      IF( JTOT.LT.NMAXIT )
     $   GO TO 10
      DO 150 I = 1, N - 1
         IF( E( I ).NE.ZERO )
     $      INFO = INFO + 1
  150 CONTINUE
      GO TO 190
*
*     Order eigenvalues and eigenvectors.
*
  160 CONTINUE
      IF( ICOMPZ.EQ.0 ) THEN
*
*        Use Quick Sort
*
         CALL DLASRT( 'I', N, D, INFO )
*
      ELSE
*
*        Use Selection Sort to minimize swaps of eigenvectors
*
         DO 180 II = 2, N
            I = II - 1
            K = I
            P = D( I )
            DO 170 J = II, N
               IF( D( J ).LT.P ) THEN
                  K = J
                  P = D( J )
               END IF
  170       CONTINUE
            IF( K.NE.I ) THEN
               D( K ) = D( I )
               D( I ) = P
               CALL DSWAP( N, Z( 1, I ), 1, Z( 1, K ), 1 )
            END IF
  180    CONTINUE
      END IF
*
  190 CONTINUE
      RETURN
*
*     End of DSTEQR
*
      END

! DSYTRF_RK
      SUBROUTINE DSYTRF_RK( UPLO, N, A, LDA, E, IPIV, WORK, LWORK,
     $                      INFO )

*  -- LAPACK computational routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*> \par Purpose:
*  =============
*>
*> \verbatim
*> DSYTRF_RK computes the factorization of a real symmetric matrix A
*> using the bounded Bunch-Kaufman (rook) diagonal pivoting method:
*>
*>    A = P*U*D*(U**T)*(P**T) or A = P*L*D*(L**T)*(P**T),
*>
*> where U (or L) is unit upper (or lower) triangular matrix,
*> U**T (or L**T) is the transpose of U (or L), P is a permutation
*> matrix, P**T is the transpose of P, and D is symmetric and block
*> diagonal with 1-by-1 and 2-by-2 diagonal blocks.
*>
*> This is the blocked version of the algorithm, calling Level 3 BLAS.
*> For more information see Further Details section.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies whether the upper or lower triangular part of the
*>          symmetric matrix A is stored:
*>          = 'U':  Upper triangular
*>          = 'L':  Lower triangular
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          On entry, the symmetric matrix A.
*>            If UPLO = 'U': the leading N-by-N upper triangular part
*>            of A contains the upper triangular part of the matrix A,
*>            and the strictly lower triangular part of A is not
*>            referenced.
*>
*>            If UPLO = 'L': the leading N-by-N lower triangular part
*>            of A contains the lower triangular part of the matrix A,
*>            and the strictly upper triangular part of A is not
*>            referenced.
*>
*>          On exit, contains:
*>            a) ONLY diagonal elements of the symmetric block diagonal
*>               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k);
*>               (superdiagonal (or subdiagonal) elements of D
*>                are stored on exit in array E), and
*>            b) If UPLO = 'U': factor U in the superdiagonal part of A.
*>               If UPLO = 'L': factor L in the subdiagonal part of A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[out] E
*> \verbatim
*>          E is DOUBLE PRECISION array, dimension (N)
*>          On exit, contains the superdiagonal (or subdiagonal)
*>          elements of the symmetric block diagonal matrix D
*>          with 1-by-1 or 2-by-2 diagonal blocks, where
*>          If UPLO = 'U': E(i) = D(i-1,i), i=2:N, E(1) is set to 0;
*>          If UPLO = 'L': E(i) = D(i+1,i), i=1:N-1, E(N) is set to 0.
*>
*>          NOTE: For 1-by-1 diagonal block D(k), where
*>          1 <= k <= N, the element E(k) is set to 0 in both
*>          UPLO = 'U' or UPLO = 'L' cases.
*> \endverbatim
*>
*> \param[out] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (N)
*>          IPIV describes the permutation matrix P in the factorization
*>          of matrix A as follows. The absolute value of IPIV(k)
*>          represents the index of row and column that were
*>          interchanged with the k-th row and column. The value of UPLO
*>          describes the order in which the interchanges were applied.
*>          Also, the sign of IPIV represents the block structure of
*>          the symmetric block diagonal matrix D with 1-by-1 or 2-by-2
*>          diagonal blocks which correspond to 1 or 2 interchanges
*>          at each factorization step. For more info see Further
*>          Details section.
*>
*>          If UPLO = 'U',
*>          ( in factorization order, k decreases from N to 1 ):
*>            a) A single positive entry IPIV(k) > 0 means:
*>               D(k,k) is a 1-by-1 diagonal block.
*>               If IPIV(k) != k, rows and columns k and IPIV(k) were
*>               interchanged in the matrix A(1:N,1:N);
*>               If IPIV(k) = k, no interchange occurred.
*>
*>            b) A pair of consecutive negative entries
*>               IPIV(k) < 0 and IPIV(k-1) < 0 means:
*>               D(k-1:k,k-1:k) is a 2-by-2 diagonal block.
*>               (NOTE: negative entries in IPIV appear ONLY in pairs).
*>               1) If -IPIV(k) != k, rows and columns
*>                  k and -IPIV(k) were interchanged
*>                  in the matrix A(1:N,1:N).
*>                  If -IPIV(k) = k, no interchange occurred.
*>               2) If -IPIV(k-1) != k-1, rows and columns
*>                  k-1 and -IPIV(k-1) were interchanged
*>                  in the matrix A(1:N,1:N).
*>                  If -IPIV(k-1) = k-1, no interchange occurred.
*>
*>            c) In both cases a) and b), always ABS( IPIV(k) ) <= k.
*>
*>            d) NOTE: Any entry IPIV(k) is always NONZERO on output.
*>
*>          If UPLO = 'L',
*>          ( in factorization order, k increases from 1 to N ):
*>            a) A single positive entry IPIV(k) > 0 means:
*>               D(k,k) is a 1-by-1 diagonal block.
*>               If IPIV(k) != k, rows and columns k and IPIV(k) were
*>               interchanged in the matrix A(1:N,1:N).
*>               If IPIV(k) = k, no interchange occurred.
*>
*>            b) A pair of consecutive negative entries
*>               IPIV(k) < 0 and IPIV(k+1) < 0 means:
*>               D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
*>               (NOTE: negative entries in IPIV appear ONLY in pairs).
*>               1) If -IPIV(k) != k, rows and columns
*>                  k and -IPIV(k) were interchanged
*>                  in the matrix A(1:N,1:N).
*>                  If -IPIV(k) = k, no interchange occurred.
*>               2) If -IPIV(k+1) != k+1, rows and columns
*>                  k-1 and -IPIV(k-1) were interchanged
*>                  in the matrix A(1:N,1:N).
*>                  If -IPIV(k+1) = k+1, no interchange occurred.
*>
*>            c) In both cases a) and b), always ABS( IPIV(k) ) >= k.
*>
*>            d) NOTE: Any entry IPIV(k) is always NONZERO on output.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension ( MAX(1,LWORK) ).
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The length of WORK.  LWORK >=1.  For best performance
*>          LWORK >= N*NB, where NB is the block size returned
*>          by ILAENV.
*>
*>          If LWORK = -1, then a workspace query is assumed;
*>          the routine only calculates the optimal size of the WORK
*>          array, returns this value as the first entry of the WORK
*>          array, and no error message related to LWORK is issued
*>          by XERBLA.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0: successful exit
*>
*>          < 0: If INFO = -k, the k-th argument had an illegal value
*>
*>          > 0: If INFO = k, the matrix A is singular, because:
*>                 If UPLO = 'U': column k in the upper
*>                 triangular part of A contains all zeros.
*>                 If UPLO = 'L': column k in the lower
*>                 triangular part of A contains all zeros.
*>
*>               Therefore D(k,k) is exactly zero, and superdiagonal
*>               elements of column k of U (or subdiagonal elements of
*>               column k of L ) are all zeros. The factorization has
*>               been completed, but the block diagonal matrix D is
*>               exactly singular, and division by zero will occur if
*>               it is used to solve a system of equations.
*>
*>               NOTE: INFO only stores the first occurrence of
*>               a singularity, any subsequent occurrence of singularity
*>               is not stored in INFO even though the factorization
*>               always completes.
*> \endverbatim
*
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), E( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            LQUERY, UPPER
      INTEGER            I, IINFO, IP, IWS, K, KB, LDWORK, LWKOPT,
     $                   NB, NBMIN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLASYF_RK, DSYTF2_RK, DSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK.EQ.-1 )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.1 .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
*
      IF( INFO.EQ.0 ) THEN
*
*        Determine the block size
*
         NB = ILAENV( 1, 'DSYTRF_RK', UPLO, N, -1, -1, -1 )
         LWKOPT = N*NB
         WORK( 1 ) = LWKOPT
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYTRF_RK', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
      NBMIN = 2
      LDWORK = N
      IF( NB.GT.1 .AND. NB.LT.N ) THEN
         IWS = LDWORK*NB
         IF( LWORK.LT.IWS ) THEN
            NB = MAX( LWORK / LDWORK, 1 )
            NBMIN = MAX( 2, ILAENV( 2, 'DSYTRF_RK',
     $                              UPLO, N, -1, -1, -1 ) )
         END IF
      ELSE
         IWS = 1
      END IF
      IF( NB.LT.NBMIN )
     $   NB = N
*
      IF( UPPER ) THEN
*
*        Factorize A as U*D*U**T using the upper triangle of A
*
*        K is the main loop index, decreasing from N to 1 in steps of
*        KB, where KB is the number of columns factorized by DLASYF_RK;
*        KB is either NB or NB-1, or K for the last block
*
         K = N
   10    CONTINUE
*
*        If K < 1, exit from loop
*
         IF( K.LT.1 )
     $      GO TO 15
*
         IF( K.GT.NB ) THEN
*
*           Factorize columns k-kb+1:k of A and use blocked code to
*           update columns 1:k-kb
*
            CALL DLASYF_RK( UPLO, K, NB, KB, A, LDA, E,
     $                      IPIV, WORK, LDWORK, IINFO )
         ELSE
*
*           Use unblocked code to factorize columns 1:k of A
*
            CALL DSYTF2_RK( UPLO, K, A, LDA, E, IPIV, IINFO )
            KB = K
         END IF
*
*        Set INFO on the first occurrence of a zero pivot
*
         IF( INFO.EQ.0 .AND. IINFO.GT.0 )
     $      INFO = IINFO
*
*        No need to adjust IPIV
*
*
*        Apply permutations to the leading panel 1:k-1
*
*        Read IPIV from the last block factored, i.e.
*        indices  k-kb+1:k and apply row permutations to the
*        last k+1 colunms k+1:N after that block
*        (We can do the simple loop over IPIV with decrement -1,
*        since the ABS value of IPIV( I ) represents the row index
*        of the interchange with row i in both 1x1 and 2x2 pivot cases)
*
         IF( K.LT.N ) THEN
            DO I = K, ( K - KB + 1 ), -1
               IP = ABS( IPIV( I ) )
               IF( IP.NE.I ) THEN
                  CALL DSWAP( N-K, A( I, K+1 ), LDA,
     $                        A( IP, K+1 ), LDA )
               END IF
            END DO
         END IF
*
*        Decrease K and return to the start of the main loop
*
         K = K - KB
         GO TO 10
*
*        This label is the exit from main loop over K decreasing
*        from N to 1 in steps of KB
*
   15    CONTINUE
*
      ELSE
*
*        Factorize A as L*D*L**T using the lower triangle of A
*
*        K is the main loop index, increasing from 1 to N in steps of
*        KB, where KB is the number of columns factorized by DLASYF_RK;
*        KB is either NB or NB-1, or N-K+1 for the last block
*
         K = 1
   20    CONTINUE
*
*        If K > N, exit from loop
*
         IF( K.GT.N )
     $      GO TO 35
*
         IF( K.LE.N-NB ) THEN
*
*           Factorize columns k:k+kb-1 of A and use blocked code to
*           update columns k+kb:n
*
            CALL DLASYF_RK( UPLO, N-K+1, NB, KB, A( K, K ), LDA, E( K ),
     $                        IPIV( K ), WORK, LDWORK, IINFO )


         ELSE
*
*           Use unblocked code to factorize columns k:n of A
*
            CALL DSYTF2_RK( UPLO, N-K+1, A( K, K ), LDA, E( K ),
     $                      IPIV( K ), IINFO )
            KB = N - K + 1
*
         END IF
*
*        Set INFO on the first occurrence of a zero pivot
*
         IF( INFO.EQ.0 .AND. IINFO.GT.0 )
     $      INFO = IINFO + K - 1
*
*        Adjust IPIV
*
         DO I = K, K + KB - 1
            IF( IPIV( I ).GT.0 ) THEN
               IPIV( I ) = IPIV( I ) + K - 1
            ELSE
               IPIV( I ) = IPIV( I ) - K + 1
            END IF
         END DO
*
*        Apply permutations to the leading panel 1:k-1
*
*        Read IPIV from the last block factored, i.e.
*        indices  k:k+kb-1 and apply row permutations to the
*        first k-1 colunms 1:k-1 before that block
*        (We can do the simple loop over IPIV with increment 1,
*        since the ABS value of IPIV( I ) represents the row index
*        of the interchange with row i in both 1x1 and 2x2 pivot cases)
*
         IF( K.GT.1 ) THEN
            DO I = K, ( K + KB - 1 ), 1
               IP = ABS( IPIV( I ) )
               IF( IP.NE.I ) THEN
                  CALL DSWAP( K-1, A( I, 1 ), LDA,
     $                        A( IP, 1 ), LDA )
               END IF
            END DO
         END IF
*
*        Increase K and return to the start of the main loop
*
         K = K + KB
         GO TO 20
*
*        This label is the exit from main loop over K increasing
*        from 1 to N in steps of KB
*
   35    CONTINUE
*
*     End Lower
*
      END IF
*
      WORK( 1 ) = LWKOPT
      RETURN
*
*     End of DSYTRF_RK
*
      END

! DORGQL
      SUBROUTINE DORGQL( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, J, KK, L, LDWORK, LWKOPT,
     $                   NB, NBMIN, NX
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARFB, DLARFT, DORG2L, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      END IF
*
      IF( INFO.EQ.0 ) THEN
         IF( N.EQ.0 ) THEN
            LWKOPT = 1
         ELSE
            NB = ILAENV( 1, 'DORGQL', ' ', M, N, K, -1 )
            LWKOPT = N*NB
         END IF
         WORK( 1 ) = LWKOPT
*
         IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
            INFO = -8
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORGQL', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.0 ) THEN
         RETURN
      END IF
*
      NBMIN = 2
      NX = 0
      IWS = N
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
*
*        Determine when to cross over from blocked to unblocked code.
*
         NX = MAX( 0, ILAENV( 3, 'DORGQL', ' ', M, N, K, -1 ) )
         IF( NX.LT.K ) THEN
*
*           Determine if workspace is large enough for blocked code.
*
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
*
*              Not enough workspace to use optimal NB:  reduce NB and
*              determine the minimum value of NB.
*
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'DORGQL', ' ', M, N, K, -1 ) )
            END IF
         END IF
      END IF
*
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
*
*        Use blocked code after the first block.
*        The last kk columns are handled by the block method.
*
         KK = MIN( K, ( ( K-NX+NB-1 ) / NB )*NB )
*
*        Set A(m-kk+1:m,1:n-kk) to zero.
*
         DO 20 J = 1, N - KK
            DO 10 I = M - KK + 1, M
               A( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
      ELSE
         KK = 0
      END IF
*
*     Use unblocked code for the first or only block.
*
      CALL DORG2L( M-KK, N-KK, K-KK, A, LDA, TAU, WORK, IINFO )
*
      IF( KK.GT.0 ) THEN
*
*        Use blocked code
*
         DO 50 I = K - KK + 1, K, NB
            IB = MIN( NB, K-I+1 )
            IF( N-K+I.GT.1 ) THEN
*
*              Form the triangular factor of the block reflector
*              H = H(i+ib-1) . . . H(i+1) H(i)
*
               CALL DLARFT( 'Backward', 'Columnwise', M-K+I+IB-1, IB,
     $                      A( 1, N-K+I ), LDA, TAU( I ), WORK, LDWORK )
*
*              Apply H to A(1:m-k+i+ib-1,1:n-k+i-1) from the left
*
               CALL DLARFB( 'Left', 'No transpose', 'Backward',
     $                      'Columnwise', M-K+I+IB-1, N-K+I-1, IB,
     $                      A( 1, N-K+I ), LDA, WORK, LDWORK, A, LDA,
     $                      WORK( IB+1 ), LDWORK )
            END IF
*
*           Apply H to rows 1:m-k+i+ib-1 of current block
*
            CALL DORG2L( M-K+I+IB-1, IB, IB, A( 1, N-K+I ), LDA,
     $                   TAU( I ), WORK, IINFO )
*
*           Set rows m-k+i+ib:m of current block to zero
*
            DO 40 J = N - K + I, N - K + I + IB - 1
               DO 30 L = M - K + I + IB, M
                  A( L, J ) = ZERO
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
      END IF
*
      WORK( 1 ) = IWS
      RETURN
*
*     End of DORGQL
*
      END

! DORGQR
      SUBROUTINE DORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, J, KI, KK, L, LDWORK,
     $                   LWKOPT, NB, NBMIN, NX
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARFB, DLARFT, DORG2R, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      NB = ILAENV( 1, 'DORGQR', ' ', M, N, K, -1 )
      LWKOPT = MAX( 1, N )*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORGQR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      NBMIN = 2
      NX = 0
      IWS = N
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
*
*        Determine when to cross over from blocked to unblocked code.
*
         NX = MAX( 0, ILAENV( 3, 'DORGQR', ' ', M, N, K, -1 ) )
         IF( NX.LT.K ) THEN
*
*           Determine if workspace is large enough for blocked code.
*
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
*
*              Not enough workspace to use optimal NB:  reduce NB and
*              determine the minimum value of NB.
*
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'DORGQR', ' ', M, N, K, -1 ) )
            END IF
         END IF
      END IF
*
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
*
*        Use blocked code after the last block.
*        The first kk columns are handled by the block method.
*
         KI = ( ( K-NX-1 ) / NB )*NB
         KK = MIN( K, KI+NB )
*
*        Set A(1:kk,kk+1:n) to zero.
*
         DO 20 J = KK + 1, N
            DO 10 I = 1, KK
               A( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
      ELSE
         KK = 0
      END IF
*
*     Use unblocked code for the last or only block.
*
      IF( KK.LT.N )
     $   CALL DORG2R( M-KK, N-KK, K-KK, A( KK+1, KK+1 ), LDA,
     $                TAU( KK+1 ), WORK, IINFO )
*
      IF( KK.GT.0 ) THEN
*
*        Use blocked code
*
         DO 50 I = KI + 1, 1, -NB
            IB = MIN( NB, K-I+1 )
            IF( I+IB.LE.N ) THEN
*
*              Form the triangular factor of the block reflector
*              H = H(i) H(i+1) . . . H(i+ib-1)
*
               CALL DLARFT( 'Forward', 'Columnwise', M-I+1, IB,
     $                      A( I, I ), LDA, TAU( I ), WORK, LDWORK )
*
*              Apply H to A(i:m,i+ib:n) from the left
*
               CALL DLARFB( 'Left', 'No transpose', 'Forward',
     $                      'Columnwise', M-I+1, N-I-IB+1, IB,
     $                      A( I, I ), LDA, WORK, LDWORK, A( I, I+IB ),
     $                      LDA, WORK( IB+1 ), LDWORK )
            END IF
*
*           Apply H to rows i:m of current block
*
            CALL DORG2R( M-I+1, IB, IB, A( I, I ), LDA, TAU( I ), WORK,
     $                   IINFO )
*
*           Set rows 1:i-1 of current block to zero
*
            DO 40 J = I, I + IB - 1
               DO 30 L = 1, I - 1
                  A( L, J ) = ZERO
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
      END IF
*
      WORK( 1 ) = IWS
      RETURN
*
*     End of DORGQR
*
      END

! DSYR2K
      SUBROUTINE DSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
*
*  -- Reference BLAS level3 routine (version 3.7.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA,BETA
      INTEGER K,LDA,LDB,LDC,N
      CHARACTER TRANS,UPLO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
*     ..
*
*  =====================================================================
*
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MAX
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION TEMP1,TEMP2
      INTEGER I,INFO,J,L,NROWA
      LOGICAL UPPER
*     ..
*     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
*     ..
*
*     Test the input parameters.
*
      IF (LSAME(TRANS,'N')) THEN
          NROWA = N
      ELSE
          NROWA = K
      END IF
      UPPER = LSAME(UPLO,'U')
*
      INFO = 0
      IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
          INFO = 1
      ELSE IF ((.NOT.LSAME(TRANS,'N')) .AND.
     +         (.NOT.LSAME(TRANS,'T')) .AND.
     +         (.NOT.LSAME(TRANS,'C'))) THEN
          INFO = 2
      ELSE IF (N.LT.0) THEN
          INFO = 3
      ELSE IF (K.LT.0) THEN
          INFO = 4
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
          INFO = 7
      ELSE IF (LDB.LT.MAX(1,NROWA)) THEN
          INFO = 9
      ELSE IF (LDC.LT.MAX(1,N)) THEN
          INFO = 12
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DSYR2K',INFO)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((N.EQ.0) .OR. (((ALPHA.EQ.ZERO).OR.
     +    (K.EQ.0)).AND. (BETA.EQ.ONE))) RETURN
*
*     And when  alpha.eq.zero.
*
      IF (ALPHA.EQ.ZERO) THEN
          IF (UPPER) THEN
              IF (BETA.EQ.ZERO) THEN
                  DO 20 J = 1,N
                      DO 10 I = 1,J
                          C(I,J) = ZERO
   10                 CONTINUE
   20             CONTINUE
              ELSE
                  DO 40 J = 1,N
                      DO 30 I = 1,J
                          C(I,J) = BETA*C(I,J)
   30                 CONTINUE
   40             CONTINUE
              END IF
          ELSE
              IF (BETA.EQ.ZERO) THEN
                  DO 60 J = 1,N
                      DO 50 I = J,N
                          C(I,J) = ZERO
   50                 CONTINUE
   60             CONTINUE
              ELSE
                  DO 80 J = 1,N
                      DO 70 I = J,N
                          C(I,J) = BETA*C(I,J)
   70                 CONTINUE
   80             CONTINUE
              END IF
          END IF
          RETURN
      END IF
*
*     Start the operations.
*
      IF (LSAME(TRANS,'N')) THEN
*
*        Form  C := alpha*A*B**T + alpha*B*A**T + C.
*
          IF (UPPER) THEN
              DO 130 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                      DO 90 I = 1,J
                          C(I,J) = ZERO
   90                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 100 I = 1,J
                          C(I,J) = BETA*C(I,J)
  100                 CONTINUE
                  END IF
                  DO 120 L = 1,K
                      IF ((A(J,L).NE.ZERO) .OR. (B(J,L).NE.ZERO)) THEN
                          TEMP1 = ALPHA*B(J,L)
                          TEMP2 = ALPHA*A(J,L)
                          DO 110 I = 1,J
                              C(I,J) = C(I,J) + A(I,L)*TEMP1 +
     +                                 B(I,L)*TEMP2
  110                     CONTINUE
                      END IF
  120             CONTINUE
  130         CONTINUE
          ELSE
              DO 180 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                      DO 140 I = J,N
                          C(I,J) = ZERO
  140                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 150 I = J,N
                          C(I,J) = BETA*C(I,J)
  150                 CONTINUE
                  END IF
                  DO 170 L = 1,K
                      IF ((A(J,L).NE.ZERO) .OR. (B(J,L).NE.ZERO)) THEN
                          TEMP1 = ALPHA*B(J,L)
                          TEMP2 = ALPHA*A(J,L)
                          DO 160 I = J,N
                              C(I,J) = C(I,J) + A(I,L)*TEMP1 +
     +                                 B(I,L)*TEMP2
  160                     CONTINUE
                      END IF
  170             CONTINUE
  180         CONTINUE
          END IF
      ELSE
*
*        Form  C := alpha*A**T*B + alpha*B**T*A + C.
*
          IF (UPPER) THEN
              DO 210 J = 1,N
                  DO 200 I = 1,J
                      TEMP1 = ZERO
                      TEMP2 = ZERO
                      DO 190 L = 1,K
                          TEMP1 = TEMP1 + A(L,I)*B(L,J)
                          TEMP2 = TEMP2 + B(L,I)*A(L,J)
  190                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP1 + ALPHA*TEMP2
                      ELSE
                          C(I,J) = BETA*C(I,J) + ALPHA*TEMP1 +
     +                             ALPHA*TEMP2
                      END IF
  200             CONTINUE
  210         CONTINUE
          ELSE
              DO 240 J = 1,N
                  DO 230 I = J,N
                      TEMP1 = ZERO
                      TEMP2 = ZERO
                      DO 220 L = 1,K
                          TEMP1 = TEMP1 + A(L,I)*B(L,J)
                          TEMP2 = TEMP2 + B(L,I)*A(L,J)
  220                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP1 + ALPHA*TEMP2
                      ELSE
                          C(I,J) = BETA*C(I,J) + ALPHA*TEMP1 +
     +                             ALPHA*TEMP2
                      END IF
  230             CONTINUE
  240         CONTINUE
          END IF
      END IF
*
      RETURN
*
*     End of DSYR2K.
*
      END

! DORG2R
      SUBROUTINE DORG2R( M, N, K, A, LDA, TAU, WORK, INFO )
*
*  -- LAPACK computational routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, L
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARF, DSCAL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORG2R', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.0 )
     $   RETURN
*
*     Initialise columns k+1:n to columns of the unit matrix
*
      DO 20 J = K + 1, N
         DO 10 L = 1, M
            A( L, J ) = ZERO
   10    CONTINUE
         A( J, J ) = ONE
   20 CONTINUE
*
      DO 40 I = K, 1, -1
*
*        Apply H(i) to A(i:m,i:n) from the left
*
         IF( I.LT.N ) THEN
            A( I, I ) = ONE
            CALL DLARF( 'Left', M-I+1, N-I, A( I, I ), 1, TAU( I ),
     $                  A( I, I+1 ), LDA, WORK )
         END IF
         IF( I.LT.M )
     $      CALL DSCAL( M-I, -TAU( I ), A( I+1, I ), 1 )
         A( I, I ) = ONE - TAU( I )
*
*        Set A(1:i-1,i) to zero
*
         DO 30 L = 1, I - 1
            A( L, I ) = ZERO
   30    CONTINUE
   40 CONTINUE
      RETURN
*
*     End of DORG2R
*
      END

! DLARFT
      SUBROUTINE DLARFT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )
*
*  -- LAPACK auxiliary routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      CHARACTER          DIRECT, STOREV
      INTEGER            K, LDT, LDV, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   T( LDT, * ), TAU( * ), V( LDV, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, PREVLASTV, LASTV
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMV, DTRMV
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( LSAME( DIRECT, 'F' ) ) THEN
         PREVLASTV = N
         DO I = 1, K
            PREVLASTV = MAX( I, PREVLASTV )
            IF( TAU( I ).EQ.ZERO ) THEN
*
*              H(i)  =  I
*
               DO J = 1, I
                  T( J, I ) = ZERO
               END DO
            ELSE
*
*              general case
*
               IF( LSAME( STOREV, 'C' ) ) THEN
*                 Skip any trailing zeros.
                  DO LASTV = N, I+1, -1
                     IF( V( LASTV, I ).NE.ZERO ) EXIT
                  END DO
                  DO J = 1, I-1
                     T( J, I ) = -TAU( I ) * V( I , J )
                  END DO
                  J = MIN( LASTV, PREVLASTV )
*
*                 T(1:i-1,i) := - tau(i) * V(i:j,1:i-1)**T * V(i:j,i)
*
                  CALL DGEMV( 'Transpose', J-I, I-1, -TAU( I ),
     $                        V( I+1, 1 ), LDV, V( I+1, I ), 1, ONE,
     $                        T( 1, I ), 1 )
               ELSE
*                 Skip any trailing zeros.
                  DO LASTV = N, I+1, -1
                     IF( V( I, LASTV ).NE.ZERO ) EXIT
                  END DO
                  DO J = 1, I-1
                     T( J, I ) = -TAU( I ) * V( J , I )
                  END DO
                  J = MIN( LASTV, PREVLASTV )
*
*                 T(1:i-1,i) := - tau(i) * V(1:i-1,i:j) * V(i,i:j)**T
*
                  CALL DGEMV( 'No transpose', I-1, J-I, -TAU( I ),
     $                        V( 1, I+1 ), LDV, V( I, I+1 ), LDV, ONE,
     $                        T( 1, I ), 1 )
               END IF
*
*              T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i)
*
               CALL DTRMV( 'Upper', 'No transpose', 'Non-unit', I-1, T,
     $                     LDT, T( 1, I ), 1 )
               T( I, I ) = TAU( I )
               IF( I.GT.1 ) THEN
                  PREVLASTV = MAX( PREVLASTV, LASTV )
               ELSE
                  PREVLASTV = LASTV
               END IF
            END IF
         END DO
      ELSE
         PREVLASTV = 1
         DO I = K, 1, -1
            IF( TAU( I ).EQ.ZERO ) THEN
*
*              H(i)  =  I
*
               DO J = I, K
                  T( J, I ) = ZERO
               END DO
            ELSE
*
*              general case
*
               IF( I.LT.K ) THEN
                  IF( LSAME( STOREV, 'C' ) ) THEN
*                    Skip any leading zeros.
                     DO LASTV = 1, I-1
                        IF( V( LASTV, I ).NE.ZERO ) EXIT
                     END DO
                     DO J = I+1, K
                        T( J, I ) = -TAU( I ) * V( N-K+I , J )
                     END DO
                     J = MAX( LASTV, PREVLASTV )
*
*                    T(i+1:k,i) = -tau(i) * V(j:n-k+i,i+1:k)**T * V(j:n-k+i,i)
*
                     CALL DGEMV( 'Transpose', N-K+I-J, K-I, -TAU( I ),
     $                           V( J, I+1 ), LDV, V( J, I ), 1, ONE,
     $                           T( I+1, I ), 1 )
                  ELSE
*                    Skip any leading zeros.
                     DO LASTV = 1, I-1
                        IF( V( I, LASTV ).NE.ZERO ) EXIT
                     END DO
                     DO J = I+1, K
                        T( J, I ) = -TAU( I ) * V( J, N-K+I )
                     END DO
                     J = MAX( LASTV, PREVLASTV )
*
*                    T(i+1:k,i) = -tau(i) * V(i+1:k,j:n-k+i) * V(i,j:n-k+i)**T
*
                     CALL DGEMV( 'No transpose', K-I, N-K+I-J,
     $                    -TAU( I ), V( I+1, J ), LDV, V( I, J ), LDV,
     $                    ONE, T( I+1, I ), 1 )
                  END IF
*
*                 T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i)
*
                  CALL DTRMV( 'Lower', 'No transpose', 'Non-unit', K-I,
     $                        T( I+1, I+1 ), LDT, T( I+1, I ), 1 )
                  IF( I.GT.1 ) THEN
                     PREVLASTV = MIN( PREVLASTV, LASTV )
                  ELSE
                     PREVLASTV = LASTV
                  END IF
               END IF
               T( I, I ) = TAU( I )
            END IF
         END DO
      END IF
      RETURN
*
*     End of DLARFT
*
      END

! DLARFB
      SUBROUTINE DLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV,
     $                   T, LDT, C, LDC, WORK, LDWORK )
*
*  -- LAPACK auxiliary routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     June 2013
*
*     .. Scalar Arguments ..
      CHARACTER          DIRECT, SIDE, STOREV, TRANS
      INTEGER            K, LDC, LDT, LDV, LDWORK, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   C( LDC, * ), T( LDT, * ), V( LDV, * ),
     $                   WORK( LDWORK, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      CHARACTER          TRANST
      INTEGER            I, J
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEMM, DTRMM
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( M.LE.0 .OR. N.LE.0 )
     $   RETURN
*
      IF( LSAME( TRANS, 'N' ) ) THEN
         TRANST = 'T'
      ELSE
         TRANST = 'N'
      END IF
*
      IF( LSAME( STOREV, 'C' ) ) THEN
*
         IF( LSAME( DIRECT, 'F' ) ) THEN
*
*           Let  V =  ( V1 )    (first K rows)
*                     ( V2 )
*           where  V1  is unit lower triangular.
*
            IF( LSAME( SIDE, 'L' ) ) THEN
*
*              Form  H * C  or  H**T * C  where  C = ( C1 )
*                                                    ( C2 )
*
*              W := C**T * V  =  (C1**T * V1 + C2**T * V2)  (stored in WORK)
*
*              W := C1**T
*
               DO 10 J = 1, K
                  CALL DCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
   10          CONTINUE
*
*              W := W * V1
*
               CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', N,
     $                     K, ONE, V, LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
*
*                 W := W + C2**T * V2
*
                  CALL DGEMM( 'Transpose', 'No transpose', N, K, M-K,
     $                        ONE, C( K+1, 1 ), LDC, V( K+1, 1 ), LDV,
     $                        ONE, WORK, LDWORK )
               END IF
*
*              W := W * T**T  or  W * T
*
               CALL DTRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - V * W**T
*
               IF( M.GT.K ) THEN
*
*                 C2 := C2 - V2 * W**T
*
                  CALL DGEMM( 'No transpose', 'Transpose', M-K, N, K,
     $                        -ONE, V( K+1, 1 ), LDV, WORK, LDWORK, ONE,
     $                        C( K+1, 1 ), LDC )
               END IF
*
*              W := W * V1**T
*
               CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', N, K,
     $                     ONE, V, LDV, WORK, LDWORK )
*
*              C1 := C1 - W**T
*
               DO 30 J = 1, K
                  DO 20 I = 1, N
                     C( J, I ) = C( J, I ) - WORK( I, J )
   20             CONTINUE
   30          CONTINUE
*
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*              Form  C * H  or  C * H**T  where  C = ( C1  C2 )
*
*              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
*
*              W := C1
*
               DO 40 J = 1, K
                  CALL DCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
   40          CONTINUE
*
*              W := W * V1
*
               CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', M,
     $                     K, ONE, V, LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
*
*                 W := W + C2 * V2
*
                  CALL DGEMM( 'No transpose', 'No transpose', M, K, N-K,
     $                        ONE, C( 1, K+1 ), LDC, V( K+1, 1 ), LDV,
     $                        ONE, WORK, LDWORK )
               END IF
*
*              W := W * T  or  W * T**T
*
               CALL DTRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - W * V**T
*
               IF( N.GT.K ) THEN
*
*                 C2 := C2 - W * V2**T
*
                  CALL DGEMM( 'No transpose', 'Transpose', M, N-K, K,
     $                        -ONE, WORK, LDWORK, V( K+1, 1 ), LDV, ONE,
     $                        C( 1, K+1 ), LDC )
               END IF
*
*              W := W * V1**T
*
               CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', M, K,
     $                     ONE, V, LDV, WORK, LDWORK )
*
*              C1 := C1 - W
*
               DO 60 J = 1, K
                  DO 50 I = 1, M
                     C( I, J ) = C( I, J ) - WORK( I, J )
   50             CONTINUE
   60          CONTINUE
            END IF
*
         ELSE
*
*           Let  V =  ( V1 )
*                     ( V2 )    (last K rows)
*           where  V2  is unit upper triangular.
*
            IF( LSAME( SIDE, 'L' ) ) THEN
*
*              Form  H * C  or  H**T * C  where  C = ( C1 )
*                                                    ( C2 )
*
*              W := C**T * V  =  (C1**T * V1 + C2**T * V2)  (stored in WORK)
*
*              W := C2**T
*
               DO 70 J = 1, K
                  CALL DCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
   70          CONTINUE
*
*              W := W * V2
*
               CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', N,
     $                     K, ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
*
*                 W := W + C1**T * V1
*
                  CALL DGEMM( 'Transpose', 'No transpose', N, K, M-K,
     $                        ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
*
*              W := W * T**T  or  W * T
*
               CALL DTRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - V * W**T
*
               IF( M.GT.K ) THEN
*
*                 C1 := C1 - V1 * W**T
*
                  CALL DGEMM( 'No transpose', 'Transpose', M-K, N, K,
     $                        -ONE, V, LDV, WORK, LDWORK, ONE, C, LDC )
               END IF
*
*              W := W * V2**T
*
               CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', N, K,
     $                     ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK )
*
*              C2 := C2 - W**T
*
               DO 90 J = 1, K
                  DO 80 I = 1, N
                     C( M-K+J, I ) = C( M-K+J, I ) - WORK( I, J )
   80             CONTINUE
   90          CONTINUE
*
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*              Form  C * H  or  C * H**T  where  C = ( C1  C2 )
*
*              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
*
*              W := C2
*
               DO 100 J = 1, K
                  CALL DCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
  100          CONTINUE
*
*              W := W * V2
*
               CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', M,
     $                     K, ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
*
*                 W := W + C1 * V1
*
                  CALL DGEMM( 'No transpose', 'No transpose', M, K, N-K,
     $                        ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
*
*              W := W * T  or  W * T**T
*
               CALL DTRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - W * V**T
*
               IF( N.GT.K ) THEN
*
*                 C1 := C1 - W * V1**T
*
                  CALL DGEMM( 'No transpose', 'Transpose', M, N-K, K,
     $                        -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC )
               END IF
*
*              W := W * V2**T
*
               CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', M, K,
     $                     ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK )
*
*              C2 := C2 - W
*
               DO 120 J = 1, K
                  DO 110 I = 1, M
                     C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
  110             CONTINUE
  120          CONTINUE
            END IF
         END IF
*
      ELSE IF( LSAME( STOREV, 'R' ) ) THEN
*
         IF( LSAME( DIRECT, 'F' ) ) THEN
*
*           Let  V =  ( V1  V2 )    (V1: first K columns)
*           where  V1  is unit upper triangular.
*
            IF( LSAME( SIDE, 'L' ) ) THEN
*
*              Form  H * C  or  H**T * C  where  C = ( C1 )
*                                                    ( C2 )
*
*              W := C**T * V**T  =  (C1**T * V1**T + C2**T * V2**T) (stored in WORK)
*
*              W := C1**T
*
               DO 130 J = 1, K
                  CALL DCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
  130          CONTINUE
*
*              W := W * V1**T
*
               CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', N, K,
     $                     ONE, V, LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
*
*                 W := W + C2**T * V2**T
*
                  CALL DGEMM( 'Transpose', 'Transpose', N, K, M-K, ONE,
     $                        C( K+1, 1 ), LDC, V( 1, K+1 ), LDV, ONE,
     $                        WORK, LDWORK )
               END IF
*
*              W := W * T**T  or  W * T
*
               CALL DTRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - V**T * W**T
*
               IF( M.GT.K ) THEN
*
*                 C2 := C2 - V2**T * W**T
*
                  CALL DGEMM( 'Transpose', 'Transpose', M-K, N, K, -ONE,
     $                        V( 1, K+1 ), LDV, WORK, LDWORK, ONE,
     $                        C( K+1, 1 ), LDC )
               END IF
*
*              W := W * V1
*
               CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', N,
     $                     K, ONE, V, LDV, WORK, LDWORK )
*
*              C1 := C1 - W**T
*
               DO 150 J = 1, K
                  DO 140 I = 1, N
                     C( J, I ) = C( J, I ) - WORK( I, J )
  140             CONTINUE
  150          CONTINUE
*
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*              Form  C * H  or  C * H**T  where  C = ( C1  C2 )
*
*              W := C * V**T  =  (C1*V1**T + C2*V2**T)  (stored in WORK)
*
*              W := C1
*
               DO 160 J = 1, K
                  CALL DCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
  160          CONTINUE
*
*              W := W * V1**T
*
               CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', M, K,
     $                     ONE, V, LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
*
*                 W := W + C2 * V2**T
*
                  CALL DGEMM( 'No transpose', 'Transpose', M, K, N-K,
     $                        ONE, C( 1, K+1 ), LDC, V( 1, K+1 ), LDV,
     $                        ONE, WORK, LDWORK )
               END IF
*
*              W := W * T  or  W * T**T
*
               CALL DTRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - W * V
*
               IF( N.GT.K ) THEN
*
*                 C2 := C2 - W * V2
*
                  CALL DGEMM( 'No transpose', 'No transpose', M, N-K, K,
     $                        -ONE, WORK, LDWORK, V( 1, K+1 ), LDV, ONE,
     $                        C( 1, K+1 ), LDC )
               END IF
*
*              W := W * V1
*
               CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', M,
     $                     K, ONE, V, LDV, WORK, LDWORK )
*
*              C1 := C1 - W
*
               DO 180 J = 1, K
                  DO 170 I = 1, M
                     C( I, J ) = C( I, J ) - WORK( I, J )
  170             CONTINUE
  180          CONTINUE
*
            END IF
*
         ELSE
*
*           Let  V =  ( V1  V2 )    (V2: last K columns)
*           where  V2  is unit lower triangular.
*
            IF( LSAME( SIDE, 'L' ) ) THEN
*
*              Form  H * C  or  H**T * C  where  C = ( C1 )
*                                                    ( C2 )
*
*              W := C**T * V**T  =  (C1**T * V1**T + C2**T * V2**T) (stored in WORK)
*
*              W := C2**T
*
               DO 190 J = 1, K
                  CALL DCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
  190          CONTINUE
*
*              W := W * V2**T
*
               CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', N, K,
     $                     ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
*
*                 W := W + C1**T * V1**T
*
                  CALL DGEMM( 'Transpose', 'Transpose', N, K, M-K, ONE,
     $                        C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
*
*              W := W * T**T  or  W * T
*
               CALL DTRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - V**T * W**T
*
               IF( M.GT.K ) THEN
*
*                 C1 := C1 - V1**T * W**T
*
                  CALL DGEMM( 'Transpose', 'Transpose', M-K, N, K, -ONE,
     $                        V, LDV, WORK, LDWORK, ONE, C, LDC )
               END IF
*
*              W := W * V2
*
               CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', N,
     $                     K, ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK )
*
*              C2 := C2 - W**T
*
               DO 210 J = 1, K
                  DO 200 I = 1, N
                     C( M-K+J, I ) = C( M-K+J, I ) - WORK( I, J )
  200             CONTINUE
  210          CONTINUE
*
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*              Form  C * H  or  C * H'  where  C = ( C1  C2 )
*
*              W := C * V**T  =  (C1*V1**T + C2*V2**T)  (stored in WORK)
*
*              W := C2
*
               DO 220 J = 1, K
                  CALL DCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
  220          CONTINUE
*
*              W := W * V2**T
*
               CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', M, K,
     $                     ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
*
*                 W := W + C1 * V1**T
*
                  CALL DGEMM( 'No transpose', 'Transpose', M, K, N-K,
     $                        ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
*
*              W := W * T  or  W * T**T
*
               CALL DTRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - W * V
*
               IF( N.GT.K ) THEN
*
*                 C1 := C1 - W * V1
*
                  CALL DGEMM( 'No transpose', 'No transpose', M, N-K, K,
     $                        -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC )
               END IF
*
*              W := W * V2
*
               CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', M,
     $                     K, ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK )
*
*              C1 := C1 - W
*
               DO 240 J = 1, K
                  DO 230 I = 1, M
                     C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
  230             CONTINUE
  240          CONTINUE
*
            END IF
*
         END IF
      END IF
*
      RETURN
*
*     End of DLARFB
*
      END

! DORG2L
      SUBROUTINE DORG2L( M, N, K, A, LDA, TAU, WORK, INFO )
*
*  -- LAPACK computational routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, II, J, L
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARF, DSCAL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORG2L', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.0 )
     $   RETURN
*
*     Initialise columns 1:n-k to columns of the unit matrix
*
      DO 20 J = 1, N - K
         DO 10 L = 1, M
            A( L, J ) = ZERO
   10    CONTINUE
         A( M-N+J, J ) = ONE
   20 CONTINUE
*
      DO 40 I = 1, K
         II = N - K + I
*
*        Apply H(i) to A(1:m-k+i,1:n-k+i) from the left
*
         A( M-N+II, II ) = ONE
         CALL DLARF( 'Left', M-N+II, II-1, A( 1, II ), 1, TAU( I ), A,
     $               LDA, WORK )
         CALL DSCAL( M-N+II-1, -TAU( I ), A( 1, II ), 1 )
         A( M-N+II, II ) = ONE - TAU( I )
*
*        Set A(m-k+i+1:m,n-k+i) to zero
*
         DO 30 L = M - N + II + 1, M
            A( L, II ) = ZERO
   30    CONTINUE
   40 CONTINUE
      RETURN
*
*     End of DORG2L
*
      END

! DLATRD
      SUBROUTINE DLATRD( UPLO, N, NB, A, LDA, E, TAU, W, LDW )
*
*  -- LAPACK auxiliary routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, LDW, N, NB
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), E( * ), TAU( * ), W( LDW, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, HALF
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, HALF = 0.5D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IW
      DOUBLE PRECISION   ALPHA
*     ..
*     .. External Subroutines ..
      EXTERNAL           DAXPY, DGEMV, DLARFG, DSCAL, DSYMV
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DDOT
      EXTERNAL           LSAME, DDOT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( N.LE.0 )
     $   RETURN
*
      IF( LSAME( UPLO, 'U' ) ) THEN
*
*        Reduce last NB columns of upper triangle
*
         DO 10 I = N, N - NB + 1, -1
            IW = I - N + NB
            IF( I.LT.N ) THEN
*
*              Update A(1:i,i)
*
               CALL DGEMV( 'No transpose', I, N-I, -ONE, A( 1, I+1 ),
     $                     LDA, W( I, IW+1 ), LDW, ONE, A( 1, I ), 1 )
               CALL DGEMV( 'No transpose', I, N-I, -ONE, W( 1, IW+1 ),
     $                     LDW, A( I, I+1 ), LDA, ONE, A( 1, I ), 1 )
            END IF
            IF( I.GT.1 ) THEN
*
*              Generate elementary reflector H(i) to annihilate
*              A(1:i-2,i)
*
               CALL DLARFG( I-1, A( I-1, I ), A( 1, I ), 1, TAU( I-1 ) )
               E( I-1 ) = A( I-1, I )
               A( I-1, I ) = ONE
*
*              Compute W(1:i-1,i)
*
               CALL DSYMV( 'Upper', I-1, ONE, A, LDA, A( 1, I ), 1,
     $                     ZERO, W( 1, IW ), 1 )
               IF( I.LT.N ) THEN
                  CALL DGEMV( 'Transpose', I-1, N-I, ONE, W( 1, IW+1 ),
     $                        LDW, A( 1, I ), 1, ZERO, W( I+1, IW ), 1 )
                  CALL DGEMV( 'No transpose', I-1, N-I, -ONE,
     $                        A( 1, I+1 ), LDA, W( I+1, IW ), 1, ONE,
     $                        W( 1, IW ), 1 )
                  CALL DGEMV( 'Transpose', I-1, N-I, ONE, A( 1, I+1 ),
     $                        LDA, A( 1, I ), 1, ZERO, W( I+1, IW ), 1 )
                  CALL DGEMV( 'No transpose', I-1, N-I, -ONE,
     $                        W( 1, IW+1 ), LDW, W( I+1, IW ), 1, ONE,
     $                        W( 1, IW ), 1 )
               END IF
               CALL DSCAL( I-1, TAU( I-1 ), W( 1, IW ), 1 )
               ALPHA = -HALF*TAU( I-1 )*DDOT( I-1, W( 1, IW ), 1,
     $                 A( 1, I ), 1 )
               CALL DAXPY( I-1, ALPHA, A( 1, I ), 1, W( 1, IW ), 1 )
            END IF
*
   10    CONTINUE
      ELSE
*
*        Reduce first NB columns of lower triangle
*
         DO 20 I = 1, NB
*
*           Update A(i:n,i)
*
            CALL DGEMV( 'No transpose', N-I+1, I-1, -ONE, A( I, 1 ),
     $                  LDA, W( I, 1 ), LDW, ONE, A( I, I ), 1 )
            CALL DGEMV( 'No transpose', N-I+1, I-1, -ONE, W( I, 1 ),
     $                  LDW, A( I, 1 ), LDA, ONE, A( I, I ), 1 )
            IF( I.LT.N ) THEN
*
*              Generate elementary reflector H(i) to annihilate
*              A(i+2:n,i)
*
               CALL DLARFG( N-I, A( I+1, I ), A( MIN( I+2, N ), I ), 1,
     $                      TAU( I ) )
               E( I ) = A( I+1, I )
               A( I+1, I ) = ONE
*
*              Compute W(i+1:n,i)
*
               CALL DSYMV( 'Lower', N-I, ONE, A( I+1, I+1 ), LDA,
     $                     A( I+1, I ), 1, ZERO, W( I+1, I ), 1 )
               CALL DGEMV( 'Transpose', N-I, I-1, ONE, W( I+1, 1 ), LDW,
     $                     A( I+1, I ), 1, ZERO, W( 1, I ), 1 )
               CALL DGEMV( 'No transpose', N-I, I-1, -ONE, A( I+1, 1 ),
     $                     LDA, W( 1, I ), 1, ONE, W( I+1, I ), 1 )
               CALL DGEMV( 'Transpose', N-I, I-1, ONE, A( I+1, 1 ), LDA,
     $                     A( I+1, I ), 1, ZERO, W( 1, I ), 1 )
               CALL DGEMV( 'No transpose', N-I, I-1, -ONE, W( I+1, 1 ),
     $                     LDW, W( 1, I ), 1, ONE, W( I+1, I ), 1 )
               CALL DSCAL( N-I, TAU( I ), W( I+1, I ), 1 )
               ALPHA = -HALF*TAU( I )*DDOT( N-I, W( I+1, I ), 1,
     $                 A( I+1, I ), 1 )
               CALL DAXPY( N-I, ALPHA, A( I+1, I ), 1, W( I+1, I ), 1 )
            END IF
*
   20    CONTINUE
      END IF
*
      RETURN
*
*     End of DLATRD
*
      END

! DSYTD2
      SUBROUTINE DSYTD2( UPLO, N, A, LDA, D, E, TAU, INFO )
*
*  -- LAPACK computational routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), D( * ), E( * ), TAU( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO, HALF
      PARAMETER          ( ONE = 1.0D0, ZERO = 0.0D0,
     $                   HALF = 1.0D0 / 2.0D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I
      DOUBLE PRECISION   ALPHA, TAUI
*     ..
*     .. External Subroutines ..
      EXTERNAL           DAXPY, DLARFG, DSYMV, DSYR2, XERBLA
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DDOT
      EXTERNAL           LSAME, DDOT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYTD2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.0 )
     $   RETURN
*
      IF( UPPER ) THEN
*
*        Reduce the upper triangle of A
*
         DO 10 I = N - 1, 1, -1
*
*           Generate elementary reflector H(i) = I - tau * v * v**T
*           to annihilate A(1:i-1,i+1)
*
            CALL DLARFG( I, A( I, I+1 ), A( 1, I+1 ), 1, TAUI )
            E( I ) = A( I, I+1 )
*
            IF( TAUI.NE.ZERO ) THEN
*
*              Apply H(i) from both sides to A(1:i,1:i)
*
               A( I, I+1 ) = ONE
*
*              Compute  x := tau * A * v  storing x in TAU(1:i)
*
               CALL DSYMV( UPLO, I, TAUI, A, LDA, A( 1, I+1 ), 1, ZERO,
     $                     TAU, 1 )
*
*              Compute  w := x - 1/2 * tau * (x**T * v) * v
*
               ALPHA = -HALF*TAUI*DDOT( I, TAU, 1, A( 1, I+1 ), 1 )
               CALL DAXPY( I, ALPHA, A( 1, I+1 ), 1, TAU, 1 )
*
*              Apply the transformation as a rank-2 update:
*                 A := A - v * w**T - w * v**T
*
               CALL DSYR2( UPLO, I, -ONE, A( 1, I+1 ), 1, TAU, 1, A,
     $                     LDA )
*
               A( I, I+1 ) = E( I )
            END IF
            D( I+1 ) = A( I+1, I+1 )
            TAU( I ) = TAUI
   10    CONTINUE
         D( 1 ) = A( 1, 1 )
      ELSE
*
*        Reduce the lower triangle of A
*
         DO 20 I = 1, N - 1
*
*           Generate elementary reflector H(i) = I - tau * v * v**T
*           to annihilate A(i+2:n,i)
*
            CALL DLARFG( N-I, A( I+1, I ), A( MIN( I+2, N ), I ), 1,
     $                   TAUI )
            E( I ) = A( I+1, I )
*
            IF( TAUI.NE.ZERO ) THEN
*
*              Apply H(i) from both sides to A(i+1:n,i+1:n)
*
               A( I+1, I ) = ONE
*
*              Compute  x := tau * A * v  storing y in TAU(i:n-1)
*
               CALL DSYMV( UPLO, N-I, TAUI, A( I+1, I+1 ), LDA,
     $                     A( I+1, I ), 1, ZERO, TAU( I ), 1 )
*
*              Compute  w := x - 1/2 * tau * (x**T * v) * v
*
               ALPHA = -HALF*TAUI*DDOT( N-I, TAU( I ), 1, A( I+1, I ),
     $                 1 )
               CALL DAXPY( N-I, ALPHA, A( I+1, I ), 1, TAU( I ), 1 )
*
*              Apply the transformation as a rank-2 update:
*                 A := A - v * w**T - w * v**T
*
               CALL DSYR2( UPLO, N-I, -ONE, A( I+1, I ), 1, TAU( I ), 1,
     $                     A( I+1, I+1 ), LDA )
*
               A( I+1, I ) = E( I )
            END IF
            D( I ) = A( I, I )
            TAU( I ) = TAUI
   20    CONTINUE
         D( N ) = A( N, N )
      END IF
*
      RETURN
*
*     End of DSYTD2
*
      END

! DLASSQ
      SUBROUTINE DLASSQ( N, X, INCX, SCALE, SUMSQ )
*
*  -- LAPACK auxiliary routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
      DOUBLE PRECISION   SCALE, SUMSQ
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
*     ..
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            IX
      DOUBLE PRECISION   ABSXI
*     ..
*     .. External Functions ..
      LOGICAL            DISNAN
      EXTERNAL           DISNAN
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
*
      IF( N.GT.0 ) THEN
         DO 10 IX = 1, 1 + ( N-1 )*INCX, INCX
            ABSXI = ABS( X( IX ) )
            IF( ABSXI.GT.ZERO.OR.DISNAN( ABSXI ) ) THEN
               IF( SCALE.LT.ABSXI ) THEN
                  SUMSQ = 1 + SUMSQ*( SCALE / ABSXI )**2
                  SCALE = ABSXI
               ELSE
                  SUMSQ = SUMSQ + ( ABSXI / SCALE )**2
               END IF
            END IF
   10    CONTINUE
      END IF
      RETURN
*
*     End of DLASSQ
*
      END

! DLARF
      SUBROUTINE DLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
*
*  -- LAPACK auxiliary routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE
      INTEGER            INCV, LDC, M, N
      DOUBLE PRECISION   TAU
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   C( LDC, * ), V( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            APPLYLEFT
      INTEGER            I, LASTV, LASTC
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMV, DGER
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILADLR, ILADLC
      EXTERNAL           LSAME, ILADLR, ILADLC
*     ..
*     .. Executable Statements ..
*
      APPLYLEFT = LSAME( SIDE, 'L' )
      LASTV = 0
      LASTC = 0
      IF( TAU.NE.ZERO ) THEN
!     Set up variables for scanning V.  LASTV begins pointing to the end
!     of V.
         IF( APPLYLEFT ) THEN
            LASTV = M
         ELSE
            LASTV = N
         END IF
         IF( INCV.GT.0 ) THEN
            I = 1 + (LASTV-1) * INCV
         ELSE
            I = 1
         END IF
!     Look for the last non-zero row in V.
         DO WHILE( LASTV.GT.0 .AND. V( I ).EQ.ZERO )
            LASTV = LASTV - 1
            I = I - INCV
         END DO
         IF( APPLYLEFT ) THEN
!     Scan for the last non-zero column in C(1:lastv,:).
            LASTC = ILADLC(LASTV, N, C, LDC)
         ELSE
!     Scan for the last non-zero row in C(:,1:lastv).
            LASTC = ILADLR(M, LASTV, C, LDC)
         END IF
      END IF
!     Note that lastc.eq.0 renders the BLAS operations null; no special
!     case is needed at this level.
      IF( APPLYLEFT ) THEN
*
*        Form  H * C
*
         IF( LASTV.GT.0 ) THEN
*
*           w(1:lastc,1) := C(1:lastv,1:lastc)**T * v(1:lastv,1)
*
            CALL DGEMV( 'Transpose', LASTV, LASTC, ONE, C, LDC, V, INCV,
     $           ZERO, WORK, 1 )
*
*           C(1:lastv,1:lastc) := C(...) - v(1:lastv,1) * w(1:lastc,1)**T
*
            CALL DGER( LASTV, LASTC, -TAU, V, INCV, WORK, 1, C, LDC )
         END IF
      ELSE
*
*        Form  C * H
*
         IF( LASTV.GT.0 ) THEN
*
*           w(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1)
*
            CALL DGEMV( 'No transpose', LASTC, LASTV, ONE, C, LDC,
     $           V, INCV, ZERO, WORK, 1 )
*
*           C(1:lastc,1:lastv) := C(...) - w(1:lastc,1) * v(1:lastv,1)**T
*
            CALL DGER( LASTC, LASTV, -TAU, WORK, 1, V, INCV, C, LDC )
         END IF
      END IF
      RETURN
*
*     End of DLARF
*
      END

! ILADLC
      INTEGER FUNCTION ILADLC( M, N, A, LDA )
*
*  -- LAPACK auxiliary routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      INTEGER            M, N, LDA
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER I
*     ..
*     .. Executable Statements ..
*
*     Quick test for the common case where one corner is non-zero.
      IF( N.EQ.0 ) THEN
         ILADLC = N
      ELSE IF( A(1, N).NE.ZERO .OR. A(M, N).NE.ZERO ) THEN
         ILADLC = N
      ELSE
*     Now scan each column from the end, returning with the first non-zero.
         DO ILADLC = N, 1, -1
            DO I = 1, M
               IF( A(I, ILADLC).NE.ZERO ) RETURN
            END DO
         END DO
      END IF
      RETURN
      END

! ILADLR
      INTEGER FUNCTION ILADLR( M, N, A, LDA )
*
*  -- LAPACK auxiliary routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      INTEGER            M, N, LDA
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER I, J
*     ..
*     .. Executable Statements ..
*
*     Quick test for the common case where one corner is non-zero.
      IF( M.EQ.0 ) THEN
         ILADLR = M
      ELSE IF( A(M, 1).NE.ZERO .OR. A(M, N).NE.ZERO ) THEN
         ILADLR = M
      ELSE
*     Scan up each column tracking the last zero row seen.
         ILADLR = 0
         DO J = 1, N
            I=M
            DO WHILE((A(MAX(I,1),J).EQ.ZERO).AND.(I.GE.1))
               I=I-1
            ENDDO
            ILADLR = MAX( ILADLR, I )
         END DO
      END IF
      RETURN
      END

! DLARFG
      SUBROUTINE DLARFG( N, ALPHA, X, INCX, TAU )
*
*  -- LAPACK auxiliary routine (version 3.8.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2017
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
      DOUBLE PRECISION   ALPHA, TAU
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            J, KNT
      DOUBLE PRECISION   BETA, RSAFMN, SAFMIN, XNORM
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLAPY2, DNRM2
      EXTERNAL           DLAMCH, DLAPY2, DNRM2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN
*     ..
*     .. External Subroutines ..
      EXTERNAL           DSCAL
*     ..
*     .. Executable Statements ..
*
      IF( N.LE.1 ) THEN
         TAU = ZERO
         RETURN
      END IF
*
      XNORM = DNRM2( N-1, X, INCX )
*
      IF( XNORM.EQ.ZERO ) THEN
*
*        H  =  I
*
         TAU = ZERO
      ELSE
*
*        general case
*
         BETA = -SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
         SAFMIN = DLAMCH( 'S' ) / DLAMCH( 'E' )
         KNT = 0
         IF( ABS( BETA ).LT.SAFMIN ) THEN
*
*           XNORM, BETA may be inaccurate; scale X and recompute them
*
            RSAFMN = ONE / SAFMIN
   10       CONTINUE
            KNT = KNT + 1
            CALL DSCAL( N-1, RSAFMN, X, INCX )
            BETA = BETA*RSAFMN
            ALPHA = ALPHA*RSAFMN
            IF( (ABS( BETA ).LT.SAFMIN) .AND. (KNT .LT. 20) )
     $         GO TO 10
*
*           New BETA is at most 1, at least SAFMIN
*
            XNORM = DNRM2( N-1, X, INCX )
            BETA = -SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
         END IF
         TAU = ( BETA-ALPHA ) / BETA
         CALL DSCAL( N-1, ONE / ( ALPHA-BETA ), X, INCX )
*
*        If ALPHA is subnormal, it may lose relative accuracy
*
         DO 20 J = 1, KNT
            BETA = BETA*SAFMIN
 20      CONTINUE
         ALPHA = BETA
      END IF
*
      RETURN
*
*     End of DLARFG
*
      END

! DSYMV
      SUBROUTINE DSYMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
*
*  -- Reference BLAS level2 routine (version 3.7.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA,BETA
      INTEGER INCX,INCY,LDA,N
      CHARACTER UPLO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),X(*),Y(*)
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION TEMP1,TEMP2
      INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY
*     ..
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MAX
*     ..
*
*     Test the input parameters.
*
      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
          INFO = 1
      ELSE IF (N.LT.0) THEN
          INFO = 2
      ELSE IF (LDA.LT.MAX(1,N)) THEN
          INFO = 5
      ELSE IF (INCX.EQ.0) THEN
          INFO = 7
      ELSE IF (INCY.EQ.0) THEN
          INFO = 10
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DSYMV ',INFO)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((N.EQ.0) .OR. ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN
*
*     Set up the start points in  X  and  Y.
*
      IF (INCX.GT.0) THEN
          KX = 1
      ELSE
          KX = 1 - (N-1)*INCX
      END IF
      IF (INCY.GT.0) THEN
          KY = 1
      ELSE
          KY = 1 - (N-1)*INCY
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through the triangular part
*     of A.
*
*     First form  y := beta*y.
*
      IF (BETA.NE.ONE) THEN
          IF (INCY.EQ.1) THEN
              IF (BETA.EQ.ZERO) THEN
                  DO 10 I = 1,N
                      Y(I) = ZERO
   10             CONTINUE
              ELSE
                  DO 20 I = 1,N
                      Y(I) = BETA*Y(I)
   20             CONTINUE
              END IF
          ELSE
              IY = KY
              IF (BETA.EQ.ZERO) THEN
                  DO 30 I = 1,N
                      Y(IY) = ZERO
                      IY = IY + INCY
   30             CONTINUE
              ELSE
                  DO 40 I = 1,N
                      Y(IY) = BETA*Y(IY)
                      IY = IY + INCY
   40             CONTINUE
              END IF
          END IF
      END IF
      IF (ALPHA.EQ.ZERO) RETURN
      IF (LSAME(UPLO,'U')) THEN
*
*        Form  y  when A is stored in upper triangle.
*
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
              DO 60 J = 1,N
                  TEMP1 = ALPHA*X(J)
                  TEMP2 = ZERO
                  DO 50 I = 1,J - 1
                      Y(I) = Y(I) + TEMP1*A(I,J)
                      TEMP2 = TEMP2 + A(I,J)*X(I)
   50             CONTINUE
                  Y(J) = Y(J) + TEMP1*A(J,J) + ALPHA*TEMP2
   60         CONTINUE
          ELSE
              JX = KX
              JY = KY
              DO 80 J = 1,N
                  TEMP1 = ALPHA*X(JX)
                  TEMP2 = ZERO
                  IX = KX
                  IY = KY
                  DO 70 I = 1,J - 1
                      Y(IY) = Y(IY) + TEMP1*A(I,J)
                      TEMP2 = TEMP2 + A(I,J)*X(IX)
                      IX = IX + INCX
                      IY = IY + INCY
   70             CONTINUE
                  Y(JY) = Y(JY) + TEMP1*A(J,J) + ALPHA*TEMP2
                  JX = JX + INCX
                  JY = JY + INCY
   80         CONTINUE
          END IF
      ELSE
*
*        Form  y  when A is stored in lower triangle.
*
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
              DO 100 J = 1,N
                  TEMP1 = ALPHA*X(J)
                  TEMP2 = ZERO
                  Y(J) = Y(J) + TEMP1*A(J,J)
                  DO 90 I = J + 1,N
                      Y(I) = Y(I) + TEMP1*A(I,J)
                      TEMP2 = TEMP2 + A(I,J)*X(I)
   90             CONTINUE
                  Y(J) = Y(J) + ALPHA*TEMP2
  100         CONTINUE
          ELSE
              JX = KX
              JY = KY
              DO 120 J = 1,N
                  TEMP1 = ALPHA*X(JX)
                  TEMP2 = ZERO
                  Y(JY) = Y(JY) + TEMP1*A(J,J)
                  IX = JX
                  IY = JY
                  DO 110 I = J + 1,N
                      IX = IX + INCX
                      IY = IY + INCY
                      Y(IY) = Y(IY) + TEMP1*A(I,J)
                      TEMP2 = TEMP2 + A(I,J)*X(IX)
  110             CONTINUE
                  Y(JY) = Y(JY) + ALPHA*TEMP2
                  JX = JX + INCX
                  JY = JY + INCY
  120         CONTINUE
          END IF
      END IF
*
      RETURN
*
*     End of DSYMV .
*
      END

! DSYR2
      SUBROUTINE DSYR2(UPLO,N,ALPHA,X,INCX,Y,INCY,A,LDA)
*
*  -- Reference BLAS level2 routine (version 3.7.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA
      INTEGER INCX,INCY,LDA,N
      CHARACTER UPLO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),X(*),Y(*)
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D+0)
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION TEMP1,TEMP2
      INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY
*     ..
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MAX
*     ..
*
*     Test the input parameters.
*
      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
          INFO = 1
      ELSE IF (N.LT.0) THEN
          INFO = 2
      ELSE IF (INCX.EQ.0) THEN
          INFO = 5
      ELSE IF (INCY.EQ.0) THEN
          INFO = 7
      ELSE IF (LDA.LT.MAX(1,N)) THEN
          INFO = 9
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DSYR2 ',INFO)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((N.EQ.0) .OR. (ALPHA.EQ.ZERO)) RETURN
*
*     Set up the start points in X and Y if the increments are not both
*     unity.
*
      IF ((INCX.NE.1) .OR. (INCY.NE.1)) THEN
          IF (INCX.GT.0) THEN
              KX = 1
          ELSE
              KX = 1 - (N-1)*INCX
          END IF
          IF (INCY.GT.0) THEN
              KY = 1
          ELSE
              KY = 1 - (N-1)*INCY
          END IF
          JX = KX
          JY = KY
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through the triangular part
*     of A.
*
      IF (LSAME(UPLO,'U')) THEN
*
*        Form  A  when A is stored in the upper triangle.
*
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
              DO 20 J = 1,N
                  IF ((X(J).NE.ZERO) .OR. (Y(J).NE.ZERO)) THEN
                      TEMP1 = ALPHA*Y(J)
                      TEMP2 = ALPHA*X(J)
                      DO 10 I = 1,J
                          A(I,J) = A(I,J) + X(I)*TEMP1 + Y(I)*TEMP2
   10                 CONTINUE
                  END IF
   20         CONTINUE
          ELSE
              DO 40 J = 1,N
                  IF ((X(JX).NE.ZERO) .OR. (Y(JY).NE.ZERO)) THEN
                      TEMP1 = ALPHA*Y(JY)
                      TEMP2 = ALPHA*X(JX)
                      IX = KX
                      IY = KY
                      DO 30 I = 1,J
                          A(I,J) = A(I,J) + X(IX)*TEMP1 + Y(IY)*TEMP2
                          IX = IX + INCX
                          IY = IY + INCY
   30                 CONTINUE
                  END IF
                  JX = JX + INCX
                  JY = JY + INCY
   40         CONTINUE
          END IF
      ELSE
*
*        Form  A  when A is stored in the lower triangle.
*
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
              DO 60 J = 1,N
                  IF ((X(J).NE.ZERO) .OR. (Y(J).NE.ZERO)) THEN
                      TEMP1 = ALPHA*Y(J)
                      TEMP2 = ALPHA*X(J)
                      DO 50 I = J,N
                          A(I,J) = A(I,J) + X(I)*TEMP1 + Y(I)*TEMP2
   50                 CONTINUE
                  END IF
   60         CONTINUE
          ELSE
              DO 80 J = 1,N
                  IF ((X(JX).NE.ZERO) .OR. (Y(JY).NE.ZERO)) THEN
                      TEMP1 = ALPHA*Y(JY)
                      TEMP2 = ALPHA*X(JX)
                      IX = JX
                      IY = JY
                      DO 70 I = J,N
                          A(I,J) = A(I,J) + X(IX)*TEMP1 + Y(IY)*TEMP2
                          IX = IX + INCX
                          IY = IY + INCY
   70                 CONTINUE
                  END IF
                  JX = JX + INCX
                  JY = JY + INCY
   80         CONTINUE
          END IF
      END IF
*
      RETURN
*
*     End of DSYR2 .
*
      END

! DISNAN
      LOGICAL FUNCTION DISNAN( DIN )
*
*  -- LAPACK auxiliary routine (version 3.7.1) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     June 2017
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION, INTENT(IN) :: DIN
*     ..
*
*  =====================================================================
*
*  .. External Functions ..
      LOGICAL DLAISNAN
      EXTERNAL DLAISNAN
*  ..
*  .. Executable Statements ..
      DISNAN = DLAISNAN(DIN,DIN)
      RETURN
      END

! DLANST
      DOUBLE PRECISION FUNCTION DLANST( NORM, N, D, E )
*
*  -- LAPACK auxiliary routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER            N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION   ANORM, SCALE, SUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME, DISNAN
      EXTERNAL           LSAME, DISNAN
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLASSQ
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SQRT
*     ..
*     .. Executable Statements ..
*
      IF( N.LE.0 ) THEN
         ANORM = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
*
*        Find max(abs(A(i,j))).
*
         ANORM = ABS( D( N ) )
         DO 10 I = 1, N - 1
            SUM = ABS( D( I ) )
            IF( ANORM .LT. SUM .OR. DISNAN( SUM ) ) ANORM = SUM
            SUM = ABS( E( I ) )
            IF( ANORM .LT. SUM .OR. DISNAN( SUM ) ) ANORM = SUM
   10    CONTINUE
      ELSE IF( LSAME( NORM, 'O' ) .OR. NORM.EQ.'1' .OR.
     $         LSAME( NORM, 'I' ) ) THEN
*
*        Find norm1(A).
*
         IF( N.EQ.1 ) THEN
            ANORM = ABS( D( 1 ) )
         ELSE
            ANORM = ABS( D( 1 ) )+ABS( E( 1 ) )
            SUM = ABS( E( N-1 ) )+ABS( D( N ) )
            IF( ANORM .LT. SUM .OR. DISNAN( SUM ) ) ANORM = SUM
            DO 20 I = 2, N - 1
               SUM = ABS( D( I ) )+ABS( E( I ) )+ABS( E( I-1 ) )
               IF( ANORM .LT. SUM .OR. DISNAN( SUM ) ) ANORM = SUM
   20       CONTINUE
         END IF
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
*
*        Find normF(A).
*
         SCALE = ZERO
         SUM = ONE
         IF( N.GT.1 ) THEN
            CALL DLASSQ( N-1, E, 1, SCALE, SUM )
            SUM = 2*SUM
         END IF
         CALL DLASSQ( N, D, 1, SCALE, SUM )
         ANORM = SCALE*SQRT( SUM )
      END IF
*
      DLANST = ANORM
      RETURN
*
*     End of DLANST
*
      END

! DLAEV2
      SUBROUTINE DLAEV2( A, B, C, RT1, RT2, CS1, SN1 )
*
*  -- LAPACK auxiliary routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, CS1, RT1, RT2, SN1
*     ..
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   HALF
      PARAMETER          ( HALF = 0.5D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            SGN1, SGN2
      DOUBLE PRECISION   AB, ACMN, ACMX, ACS, ADF, CS, CT, DF, RT, SM,
     $                   TB, TN
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SQRT
*     ..
*     .. Executable Statements ..
*
*     Compute the eigenvalues
*
      SM = A + C
      DF = A - C
      ADF = ABS( DF )
      TB = B + B
      AB = ABS( TB )
      IF( ABS( A ).GT.ABS( C ) ) THEN
         ACMX = A
         ACMN = C
      ELSE
         ACMX = C
         ACMN = A
      END IF
      IF( ADF.GT.AB ) THEN
         RT = ADF*SQRT( ONE+( AB / ADF )**2 )
      ELSE IF( ADF.LT.AB ) THEN
         RT = AB*SQRT( ONE+( ADF / AB )**2 )
      ELSE
*
*        Includes case AB=ADF=0
*
         RT = AB*SQRT( TWO )
      END IF
      IF( SM.LT.ZERO ) THEN
         RT1 = HALF*( SM-RT )
         SGN1 = -1
*
*        Order of execution important.
*        To get fully accurate smaller eigenvalue,
*        next line needs to be executed in higher precision.
*
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE IF( SM.GT.ZERO ) THEN
         RT1 = HALF*( SM+RT )
         SGN1 = 1
*
*        Order of execution important.
*        To get fully accurate smaller eigenvalue,
*        next line needs to be executed in higher precision.
*
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE
*
*        Includes case RT1 = RT2 = 0
*
         RT1 = HALF*RT
         RT2 = -HALF*RT
         SGN1 = 1
      END IF
*
*     Compute the eigenvector
*
      IF( DF.GE.ZERO ) THEN
         CS = DF + RT
         SGN2 = 1
      ELSE
         CS = DF - RT
         SGN2 = -1
      END IF
      ACS = ABS( CS )
      IF( ACS.GT.AB ) THEN
         CT = -TB / CS
         SN1 = ONE / SQRT( ONE+CT*CT )
         CS1 = CT*SN1
      ELSE
         IF( AB.EQ.ZERO ) THEN
            CS1 = ONE
            SN1 = ZERO
         ELSE
            TN = -CS / TB
            CS1 = ONE / SQRT( ONE+TN*TN )
            SN1 = TN*CS1
         END IF
      END IF
      IF( SGN1.EQ.SGN2 ) THEN
         TN = CS1
         CS1 = -SN1
         SN1 = TN
      END IF
      RETURN
*
*     End of DLAEV2
*
      END

! DLASR
      SUBROUTINE DLASR( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA )
*
*  -- LAPACK auxiliary routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      CHARACTER          DIRECT, PIVOT, SIDE
      INTEGER            LDA, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( * ), S( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, INFO, J
      DOUBLE PRECISION   CTEMP, STEMP, TEMP
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      IF( .NOT.( LSAME( SIDE, 'L' ) .OR. LSAME( SIDE, 'R' ) ) ) THEN
         INFO = 1
      ELSE IF( .NOT.( LSAME( PIVOT, 'V' ) .OR. LSAME( PIVOT,
     $         'T' ) .OR. LSAME( PIVOT, 'B' ) ) ) THEN
         INFO = 2
      ELSE IF( .NOT.( LSAME( DIRECT, 'F' ) .OR. LSAME( DIRECT, 'B' ) ) )
     $          THEN
         INFO = 3
      ELSE IF( M.LT.0 ) THEN
         INFO = 4
      ELSE IF( N.LT.0 ) THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLASR ', INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( ( M.EQ.0 ) .OR. ( N.EQ.0 ) )
     $   RETURN
      IF( LSAME( SIDE, 'L' ) ) THEN
*
*        Form  P * A
*
         IF( LSAME( PIVOT, 'V' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 20 J = 1, M - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 10 I = 1, N
                        TEMP = A( J+1, I )
                        A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                        A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
   10                CONTINUE
                  END IF
   20          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 40 J = M - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 30 I = 1, N
                        TEMP = A( J+1, I )
                        A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                        A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
   30                CONTINUE
                  END IF
   40          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'T' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 60 J = 2, M
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 50 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                        A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
   50                CONTINUE
                  END IF
   60          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 80 J = M, 2, -1
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 70 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                        A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
   70                CONTINUE
                  END IF
   80          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'B' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 100 J = 1, M - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 90 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                        A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
   90                CONTINUE
                  END IF
  100          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 120 J = M - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 110 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                        A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
  110                CONTINUE
                  END IF
  120          CONTINUE
            END IF
         END IF
      ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*        Form A * P**T
*
         IF( LSAME( PIVOT, 'V' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 140 J = 1, N - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 130 I = 1, M
                        TEMP = A( I, J+1 )
                        A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  130                CONTINUE
                  END IF
  140          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 160 J = N - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 150 I = 1, M
                        TEMP = A( I, J+1 )
                        A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  150                CONTINUE
                  END IF
  160          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'T' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 180 J = 2, N
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 170 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                        A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
  170                CONTINUE
                  END IF
  180          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 200 J = N, 2, -1
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 190 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                        A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
  190                CONTINUE
                  END IF
  200          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'B' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 220 J = 1, N - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 210 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                        A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
  210                CONTINUE
                  END IF
  220          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 240 J = N - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 230 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                        A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
  230                CONTINUE
                  END IF
  240          CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of DLASR
*
      END

! DLAE2
      SUBROUTINE DLAE2( A, B, C, RT1, RT2 )
*
*  -- LAPACK auxiliary routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, RT1, RT2
*     ..
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   HALF
      PARAMETER          ( HALF = 0.5D0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   AB, ACMN, ACMX, ADF, DF, RT, SM, TB
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SQRT
*     ..
*     .. Executable Statements ..
*
*     Compute the eigenvalues
*
      SM = A + C
      DF = A - C
      ADF = ABS( DF )
      TB = B + B
      AB = ABS( TB )
      IF( ABS( A ).GT.ABS( C ) ) THEN
         ACMX = A
         ACMN = C
      ELSE
         ACMX = C
         ACMN = A
      END IF
      IF( ADF.GT.AB ) THEN
         RT = ADF*SQRT( ONE+( AB / ADF )**2 )
      ELSE IF( ADF.LT.AB ) THEN
         RT = AB*SQRT( ONE+( ADF / AB )**2 )
      ELSE
*
*        Includes case AB=ADF=0
*
         RT = AB*SQRT( TWO )
      END IF
      IF( SM.LT.ZERO ) THEN
         RT1 = HALF*( SM-RT )
*
*        Order of execution important.
*        To get fully accurate smaller eigenvalue,
*        next line needs to be executed in higher precision.
*
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE IF( SM.GT.ZERO ) THEN
         RT1 = HALF*( SM+RT )
*
*        Order of execution important.
*        To get fully accurate smaller eigenvalue,
*        next line needs to be executed in higher precision.
*
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE
*
*        Includes case RT1 = RT2 = 0
*
         RT1 = HALF*RT
         RT2 = -HALF*RT
      END IF
      RETURN
*
*     End of DLAE2
*
      END

! DLASET
      SUBROUTINE DLASET( UPLO, M, N, ALPHA, BETA, A, LDA )
*
*  -- LAPACK auxiliary routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, M, N
      DOUBLE PRECISION   ALPHA, BETA
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, J
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
      IF( LSAME( UPLO, 'U' ) ) THEN
*
*        Set the strictly upper triangular or trapezoidal part of the
*        array to ALPHA.
*
         DO 20 J = 2, N
            DO 10 I = 1, MIN( J-1, M )
               A( I, J ) = ALPHA
   10       CONTINUE
   20    CONTINUE
*
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
*
*        Set the strictly lower triangular or trapezoidal part of the
*        array to ALPHA.
*
         DO 40 J = 1, MIN( M, N )
            DO 30 I = J + 1, M
               A( I, J ) = ALPHA
   30       CONTINUE
   40    CONTINUE
*
      ELSE
*
*        Set the leading m-by-n submatrix to ALPHA.
*
         DO 60 J = 1, N
            DO 50 I = 1, M
               A( I, J ) = ALPHA
   50       CONTINUE
   60    CONTINUE
      END IF
*
*     Set the first min(M,N) diagonal elements to BETA.
*
      DO 70 I = 1, MIN( M, N )
         A( I, I ) = BETA
   70 CONTINUE
*
      RETURN
*
*     End of DLASET
*
      END

! DSYTF2_RK
      SUBROUTINE DSYTF2_RK( UPLO, N, A, LDA, E, IPIV, INFO )
*
*  -- LAPACK computational routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), E( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION   EIGHT, SEVTEN
      PARAMETER          ( EIGHT = 8.0D+0, SEVTEN = 17.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER, DONE
      INTEGER            I, IMAX, J, JMAX, ITEMP, K, KK, KP, KSTEP,
     $                   P, II
      DOUBLE PRECISION   ABSAKK, ALPHA, COLMAX, D11, D12, D21, D22,
     $                   ROWMAX, DTEMP, T, WK, WKM1, WKP1, SFMIN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IDAMAX
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, IDAMAX, DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           DSCAL, DSWAP, DSYR, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYTF2_RK', -INFO )
         RETURN
      END IF
*
*     Initialize ALPHA for use in choosing pivot block size.
*
      ALPHA = ( ONE+SQRT( SEVTEN ) ) / EIGHT
*
*     Compute machine safe minimum
*
      SFMIN = DLAMCH( 'S' )
*
      IF( UPPER ) THEN
*
*        Factorize A as U*D*U**T using the upper triangle of A
*
*        Initialize the first entry of array E, where superdiagonal
*        elements of D are stored
*
         E( 1 ) = ZERO
*
*        K is the main loop index, decreasing from N to 1 in steps of
*        1 or 2
*
         K = N
   10    CONTINUE
*
*        If K < 1, exit from loop
*
         IF( K.LT.1 )
     $      GO TO 34
         KSTEP = 1
         P = K
*
*        Determine rows and columns to be interchanged and whether
*        a 1-by-1 or 2-by-2 pivot block will be used
*
         ABSAKK = ABS( A( K, K ) )
*
*        IMAX is the row-index of the largest off-diagonal element in
*        column K, and COLMAX is its absolute value.
*        Determine both COLMAX and IMAX.
*
         IF( K.GT.1 ) THEN
            IMAX = IDAMAX( K-1, A( 1, K ), 1 )
            COLMAX = ABS( A( IMAX, K ) )
         ELSE
            COLMAX = ZERO
         END IF
*
         IF( (MAX( ABSAKK, COLMAX ).EQ.ZERO) ) THEN
*
*           Column K is zero or underflow: set INFO and continue
*
            IF( INFO.EQ.0 )
     $         INFO = K
            KP = K
*
*           Set E( K ) to zero
*
            IF( K.GT.1 )
     $         E( K ) = ZERO
*
         ELSE
*
*           Test for interchange
*
*           Equivalent to testing for (used to handle NaN and Inf)
*           ABSAKK.GE.ALPHA*COLMAX
*
            IF( .NOT.( ABSAKK.LT.ALPHA*COLMAX ) ) THEN
*
*              no interchange,
*              use 1-by-1 pivot block
*
               KP = K
            ELSE
*
               DONE = .FALSE.
*
*              Loop until pivot found
*
   12          CONTINUE
*
*                 Begin pivot search loop body
*
*                 JMAX is the column-index of the largest off-diagonal
*                 element in row IMAX, and ROWMAX is its absolute value.
*                 Determine both ROWMAX and JMAX.
*
                  IF( IMAX.NE.K ) THEN
                     JMAX = IMAX + IDAMAX( K-IMAX, A( IMAX, IMAX+1 ),
     $                                    LDA )
                     ROWMAX = ABS( A( IMAX, JMAX ) )
                  ELSE
                     ROWMAX = ZERO
                  END IF
*
                  IF( IMAX.GT.1 ) THEN
                     ITEMP = IDAMAX( IMAX-1, A( 1, IMAX ), 1 )
                     DTEMP = ABS( A( ITEMP, IMAX ) )
                     IF( DTEMP.GT.ROWMAX ) THEN
                        ROWMAX = DTEMP
                        JMAX = ITEMP
                     END IF
                  END IF
*
*                 Equivalent to testing for (used to handle NaN and Inf)
*                 ABS( A( IMAX, IMAX ) ).GE.ALPHA*ROWMAX
*
                  IF( .NOT.( ABS( A( IMAX, IMAX ) ).LT.ALPHA*ROWMAX ) )
     $            THEN
*
*                    interchange rows and columns K and IMAX,
*                    use 1-by-1 pivot block
*
                     KP = IMAX
                     DONE = .TRUE.
*
*                 Equivalent to testing for ROWMAX .EQ. COLMAX,
*                 used to handle NaN and Inf
*
                  ELSE IF( ( P.EQ.JMAX ).OR.( ROWMAX.LE.COLMAX ) ) THEN
*
*                    interchange rows and columns K+1 and IMAX,
*                    use 2-by-2 pivot block
*
                     KP = IMAX
                     KSTEP = 2
                     DONE = .TRUE.
                  ELSE
*
*                    Pivot NOT found, set variables and repeat
*
                     P = IMAX
                     COLMAX = ROWMAX
                     IMAX = JMAX
                  END IF
*
*                 End pivot search loop body
*
               IF( .NOT. DONE ) GOTO 12
*
            END IF
*
*           Swap TWO rows and TWO columns
*
*           First swap
*
            IF( ( KSTEP.EQ.2 ) .AND. ( P.NE.K ) ) THEN
*
*              Interchange rows and column K and P in the leading
*              submatrix A(1:k,1:k) if we have a 2-by-2 pivot
*
               IF( P.GT.1 )
     $            CALL DSWAP( P-1, A( 1, K ), 1, A( 1, P ), 1 )
               IF( P.LT.(K-1) )
     $            CALL DSWAP( K-P-1, A( P+1, K ), 1, A( P, P+1 ),
     $                     LDA )
               T = A( K, K )
               A( K, K ) = A( P, P )
               A( P, P ) = T
*
*              Convert upper triangle of A into U form by applying
*              the interchanges in columns k+1:N.
*
               IF( K.LT.N )
     $            CALL DSWAP( N-K, A( K, K+1 ), LDA, A( P, K+1 ), LDA )
*
            END IF
*
*           Second swap
*
            KK = K - KSTEP + 1
            IF( KP.NE.KK ) THEN
*
*              Interchange rows and columns KK and KP in the leading
*              submatrix A(1:k,1:k)
*
               IF( KP.GT.1 )
     $            CALL DSWAP( KP-1, A( 1, KK ), 1, A( 1, KP ), 1 )
               IF( ( KK.GT.1 ) .AND. ( KP.LT.(KK-1) ) )
     $            CALL DSWAP( KK-KP-1, A( KP+1, KK ), 1, A( KP, KP+1 ),
     $                     LDA )
               T = A( KK, KK )
               A( KK, KK ) = A( KP, KP )
               A( KP, KP ) = T
               IF( KSTEP.EQ.2 ) THEN
                  T = A( K-1, K )
                  A( K-1, K ) = A( KP, K )
                  A( KP, K ) = T
               END IF
*
*              Convert upper triangle of A into U form by applying
*              the interchanges in columns k+1:N.
*
               IF( K.LT.N )
     $            CALL DSWAP( N-K, A( KK, K+1 ), LDA, A( KP, K+1 ),
     $                        LDA )
*
            END IF
*
*           Update the leading submatrix
*
            IF( KSTEP.EQ.1 ) THEN
*
*              1-by-1 pivot block D(k): column k now holds
*
*              W(k) = U(k)*D(k)
*
*              where U(k) is the k-th column of U
*
               IF( K.GT.1 ) THEN
*
*                 Perform a rank-1 update of A(1:k-1,1:k-1) and
*                 store U(k) in column k
*
                  IF( ABS( A( K, K ) ).GE.SFMIN ) THEN
*
*                    Perform a rank-1 update of A(1:k-1,1:k-1) as
*                    A := A - U(k)*D(k)*U(k)**T
*                       = A - W(k)*1/D(k)*W(k)**T
*
                     D11 = ONE / A( K, K )
                     CALL DSYR( UPLO, K-1, -D11, A( 1, K ), 1, A, LDA )
*
*                    Store U(k) in column k
*
                     CALL DSCAL( K-1, D11, A( 1, K ), 1 )
                  ELSE
*
*                    Store L(k) in column K
*
                     D11 = A( K, K )
                     DO 16 II = 1, K - 1
                        A( II, K ) = A( II, K ) / D11
   16                CONTINUE
*
*                    Perform a rank-1 update of A(k+1:n,k+1:n) as
*                    A := A - U(k)*D(k)*U(k)**T
*                       = A - W(k)*(1/D(k))*W(k)**T
*                       = A - (W(k)/D(k))*(D(k))*(W(k)/D(K))**T
*
                     CALL DSYR( UPLO, K-1, -D11, A( 1, K ), 1, A, LDA )
                  END IF
*
*                 Store the superdiagonal element of D in array E
*
                  E( K ) = ZERO
*
               END IF
*
            ELSE
*
*              2-by-2 pivot block D(k): columns k and k-1 now hold
*
*              ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)
*
*              where U(k) and U(k-1) are the k-th and (k-1)-th columns
*              of U
*
*              Perform a rank-2 update of A(1:k-2,1:k-2) as
*
*              A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )**T
*                 = A - ( ( A(k-1)A(k) )*inv(D(k)) ) * ( A(k-1)A(k) )**T
*
*              and store L(k) and L(k+1) in columns k and k+1
*
               IF( K.GT.2 ) THEN
*
                  D12 = A( K-1, K )
                  D22 = A( K-1, K-1 ) / D12
                  D11 = A( K, K ) / D12
                  T = ONE / ( D11*D22-ONE )
*
                  DO 30 J = K - 2, 1, -1
*
                     WKM1 = T*( D11*A( J, K-1 )-A( J, K ) )
                     WK = T*( D22*A( J, K )-A( J, K-1 ) )
*
                     DO 20 I = J, 1, -1
                        A( I, J ) = A( I, J ) - (A( I, K ) / D12 )*WK -
     $                              ( A( I, K-1 ) / D12 )*WKM1
   20                CONTINUE
*
*                    Store U(k) and U(k-1) in cols k and k-1 for row J
*
                     A( J, K ) = WK / D12
                     A( J, K-1 ) = WKM1 / D12
*
   30             CONTINUE
*
               END IF
*
*              Copy superdiagonal elements of D(K) to E(K) and
*              ZERO out superdiagonal entry of A
*
               E( K ) = A( K-1, K )
               E( K-1 ) = ZERO
               A( K-1, K ) = ZERO
*
            END IF
*
*           End column K is nonsingular
*
         END IF
*
*        Store details of the interchanges in IPIV
*
         IF( KSTEP.EQ.1 ) THEN
            IPIV( K ) = KP
         ELSE
            IPIV( K ) = -P
            IPIV( K-1 ) = -KP
         END IF
*
*        Decrease K and return to the start of the main loop
*
         K = K - KSTEP
         GO TO 10
*
   34    CONTINUE
*
      ELSE
*
*        Factorize A as L*D*L**T using the lower triangle of A
*
*        Initialize the unused last entry of the subdiagonal array E.
*
         E( N ) = ZERO
*
*        K is the main loop index, increasing from 1 to N in steps of
*        1 or 2
*
         K = 1
   40    CONTINUE
*
*        If K > N, exit from loop
*
         IF( K.GT.N )
     $      GO TO 64
         KSTEP = 1
         P = K
*
*        Determine rows and columns to be interchanged and whether
*        a 1-by-1 or 2-by-2 pivot block will be used
*
         ABSAKK = ABS( A( K, K ) )
*
*        IMAX is the row-index of the largest off-diagonal element in
*        column K, and COLMAX is its absolute value.
*        Determine both COLMAX and IMAX.
*
         IF( K.LT.N ) THEN
            IMAX = K + IDAMAX( N-K, A( K+1, K ), 1 )
            COLMAX = ABS( A( IMAX, K ) )
         ELSE
            COLMAX = ZERO
         END IF
*
         IF( ( MAX( ABSAKK, COLMAX ).EQ.ZERO ) ) THEN
*
*           Column K is zero or underflow: set INFO and continue
*
            IF( INFO.EQ.0 )
     $         INFO = K
            KP = K
*
*           Set E( K ) to zero
*
            IF( K.LT.N )
     $         E( K ) = ZERO
*
         ELSE
*
*           Test for interchange
*
*           Equivalent to testing for (used to handle NaN and Inf)
*           ABSAKK.GE.ALPHA*COLMAX
*
            IF( .NOT.( ABSAKK.LT.ALPHA*COLMAX ) ) THEN
*
*              no interchange, use 1-by-1 pivot block
*
               KP = K
*
            ELSE
*
               DONE = .FALSE.
*
*              Loop until pivot found
*
   42          CONTINUE
*
*                 Begin pivot search loop body
*
*                 JMAX is the column-index of the largest off-diagonal
*                 element in row IMAX, and ROWMAX is its absolute value.
*                 Determine both ROWMAX and JMAX.
*
                  IF( IMAX.NE.K ) THEN
                     JMAX = K - 1 + IDAMAX( IMAX-K, A( IMAX, K ), LDA )
                     ROWMAX = ABS( A( IMAX, JMAX ) )
                  ELSE
                     ROWMAX = ZERO
                  END IF
*
                  IF( IMAX.LT.N ) THEN
                     ITEMP = IMAX + IDAMAX( N-IMAX, A( IMAX+1, IMAX ),
     $                                     1 )
                     DTEMP = ABS( A( ITEMP, IMAX ) )
                     IF( DTEMP.GT.ROWMAX ) THEN
                        ROWMAX = DTEMP
                        JMAX = ITEMP
                     END IF
                  END IF
*
*                 Equivalent to testing for (used to handle NaN and Inf)
*                 ABS( A( IMAX, IMAX ) ).GE.ALPHA*ROWMAX
*
                  IF( .NOT.( ABS( A( IMAX, IMAX ) ).LT.ALPHA*ROWMAX ) )
     $            THEN
*
*                    interchange rows and columns K and IMAX,
*                    use 1-by-1 pivot block
*
                     KP = IMAX
                     DONE = .TRUE.
*
*                 Equivalent to testing for ROWMAX .EQ. COLMAX,
*                 used to handle NaN and Inf
*
                  ELSE IF( ( P.EQ.JMAX ).OR.( ROWMAX.LE.COLMAX ) ) THEN
*
*                    interchange rows and columns K+1 and IMAX,
*                    use 2-by-2 pivot block
*
                     KP = IMAX
                     KSTEP = 2
                     DONE = .TRUE.
                  ELSE
*
*                    Pivot NOT found, set variables and repeat
*
                     P = IMAX
                     COLMAX = ROWMAX
                     IMAX = JMAX
                  END IF
*
*                 End pivot search loop body
*
               IF( .NOT. DONE ) GOTO 42
*
            END IF
*
*           Swap TWO rows and TWO columns
*
*           First swap
*
            IF( ( KSTEP.EQ.2 ) .AND. ( P.NE.K ) ) THEN
*
*              Interchange rows and column K and P in the trailing
*              submatrix A(k:n,k:n) if we have a 2-by-2 pivot
*
               IF( P.LT.N )
     $            CALL DSWAP( N-P, A( P+1, K ), 1, A( P+1, P ), 1 )
               IF( P.GT.(K+1) )
     $            CALL DSWAP( P-K-1, A( K+1, K ), 1, A( P, K+1 ), LDA )
               T = A( K, K )
               A( K, K ) = A( P, P )
               A( P, P ) = T
*
*              Convert lower triangle of A into L form by applying
*              the interchanges in columns 1:k-1.
*
               IF ( K.GT.1 )
     $            CALL DSWAP( K-1, A( K, 1 ), LDA, A( P, 1 ), LDA )
*
            END IF
*
*           Second swap
*
            KK = K + KSTEP - 1
            IF( KP.NE.KK ) THEN
*
*              Interchange rows and columns KK and KP in the trailing
*              submatrix A(k:n,k:n)
*
               IF( KP.LT.N )
     $            CALL DSWAP( N-KP, A( KP+1, KK ), 1, A( KP+1, KP ), 1 )
               IF( ( KK.LT.N ) .AND. ( KP.GT.(KK+1) ) )
     $            CALL DSWAP( KP-KK-1, A( KK+1, KK ), 1, A( KP, KK+1 ),
     $                     LDA )
               T = A( KK, KK )
               A( KK, KK ) = A( KP, KP )
               A( KP, KP ) = T
               IF( KSTEP.EQ.2 ) THEN
                  T = A( K+1, K )
                  A( K+1, K ) = A( KP, K )
                  A( KP, K ) = T
               END IF
*
*              Convert lower triangle of A into L form by applying
*              the interchanges in columns 1:k-1.
*
               IF ( K.GT.1 )
     $            CALL DSWAP( K-1, A( KK, 1 ), LDA, A( KP, 1 ), LDA )
*
            END IF
*
*           Update the trailing submatrix
*
            IF( KSTEP.EQ.1 ) THEN
*
*              1-by-1 pivot block D(k): column k now holds
*
*              W(k) = L(k)*D(k)
*
*              where L(k) is the k-th column of L
*
               IF( K.LT.N ) THEN
*
*              Perform a rank-1 update of A(k+1:n,k+1:n) and
*              store L(k) in column k
*
                  IF( ABS( A( K, K ) ).GE.SFMIN ) THEN
*
*                    Perform a rank-1 update of A(k+1:n,k+1:n) as
*                    A := A - L(k)*D(k)*L(k)**T
*                       = A - W(k)*(1/D(k))*W(k)**T
*
                     D11 = ONE / A( K, K )
                     CALL DSYR( UPLO, N-K, -D11, A( K+1, K ), 1,
     $                          A( K+1, K+1 ), LDA )
*
*                    Store L(k) in column k
*
                     CALL DSCAL( N-K, D11, A( K+1, K ), 1 )
                  ELSE
*
*                    Store L(k) in column k
*
                     D11 = A( K, K )
                     DO 46 II = K + 1, N
                        A( II, K ) = A( II, K ) / D11
   46                CONTINUE
*
*                    Perform a rank-1 update of A(k+1:n,k+1:n) as
*                    A := A - L(k)*D(k)*L(k)**T
*                       = A - W(k)*(1/D(k))*W(k)**T
*                       = A - (W(k)/D(k))*(D(k))*(W(k)/D(K))**T
*
                     CALL DSYR( UPLO, N-K, -D11, A( K+1, K ), 1,
     $                          A( K+1, K+1 ), LDA )
                  END IF
*
*                 Store the subdiagonal element of D in array E
*
                  E( K ) = ZERO
*
               END IF
*
            ELSE
*
*              2-by-2 pivot block D(k): columns k and k+1 now hold
*
*              ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k)
*
*              where L(k) and L(k+1) are the k-th and (k+1)-th columns
*              of L
*
*
*              Perform a rank-2 update of A(k+2:n,k+2:n) as
*
*              A := A - ( L(k) L(k+1) ) * D(k) * ( L(k) L(k+1) )**T
*                 = A - ( ( A(k)A(k+1) )*inv(D(k) ) * ( A(k)A(k+1) )**T
*
*              and store L(k) and L(k+1) in columns k and k+1
*
               IF( K.LT.N-1 ) THEN
*
                  D21 = A( K+1, K )
                  D11 = A( K+1, K+1 ) / D21
                  D22 = A( K, K ) / D21
                  T = ONE / ( D11*D22-ONE )
*
                  DO 60 J = K + 2, N
*
*                    Compute  D21 * ( W(k)W(k+1) ) * inv(D(k)) for row J
*
                     WK = T*( D11*A( J, K )-A( J, K+1 ) )
                     WKP1 = T*( D22*A( J, K+1 )-A( J, K ) )
*
*                    Perform a rank-2 update of A(k+2:n,k+2:n)
*
                     DO 50 I = J, N
                        A( I, J ) = A( I, J ) - ( A( I, K ) / D21 )*WK -
     $                              ( A( I, K+1 ) / D21 )*WKP1
   50                CONTINUE
*
*                    Store L(k) and L(k+1) in cols k and k+1 for row J
*
                     A( J, K ) = WK / D21
                     A( J, K+1 ) = WKP1 / D21
*
   60             CONTINUE
*
               END IF
*
*              Copy subdiagonal elements of D(K) to E(K) and
*              ZERO out subdiagonal entry of A
*
               E( K ) = A( K+1, K )
               E( K+1 ) = ZERO
               A( K+1, K ) = ZERO
*
            END IF
*
*           End column K is nonsingular
*
         END IF
*
*        Store details of the interchanges in IPIV
*
         IF( KSTEP.EQ.1 ) THEN
            IPIV( K ) = KP
         ELSE
            IPIV( K ) = -P
            IPIV( K+1 ) = -KP
         END IF
*
*        Increase K and return to the start of the main loop
*
         K = K + KSTEP
         GO TO 40
*
   64    CONTINUE
*
      END IF
*
      RETURN
*
*     End of DSYTF2_RK
*
      END

! DLASYF_RK
      SUBROUTINE DLASYF_RK( UPLO, N, NB, KB, A, LDA, E, IPIV, W, LDW,
     $                      INFO )
*
*  -- LAPACK computational routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, KB, LDA, LDW, N, NB
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), E( * ), W( LDW, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION   EIGHT, SEVTEN
      PARAMETER          ( EIGHT = 8.0D+0, SEVTEN = 17.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            DONE
      INTEGER            IMAX, ITEMP, J, JB, JJ, JMAX, K, KK, KW, KKW,
     $                   KP, KSTEP, P, II
      DOUBLE PRECISION   ABSAKK, ALPHA, COLMAX, D11, D12, D21, D22,
     $                   DTEMP, R1, ROWMAX, T, SFMIN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IDAMAX
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, IDAMAX, DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEMM, DGEMV, DSCAL, DSWAP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
      INFO = 0
*
*     Initialize ALPHA for use in choosing pivot block size.
*
      ALPHA = ( ONE+SQRT( SEVTEN ) ) / EIGHT
*
*     Compute machine safe minimum
*
      SFMIN = DLAMCH( 'S' )
*
      IF( LSAME( UPLO, 'U' ) ) THEN
*
*        Factorize the trailing columns of A using the upper triangle
*        of A and working backwards, and compute the matrix W = U12*D
*        for use in updating A11
*
*        Initialize the first entry of array E, where superdiagonal
*        elements of D are stored
*
         E( 1 ) = ZERO
*
*        K is the main loop index, decreasing from N in steps of 1 or 2
*
         K = N
   10    CONTINUE
*
*        KW is the column of W which corresponds to column K of A
*
         KW = NB + K - N
*
*        Exit from loop
*
         IF( ( K.LE.N-NB+1 .AND. NB.LT.N ) .OR. K.LT.1 )
     $      GO TO 30
*
         KSTEP = 1
         P = K
*
*        Copy column K of A to column KW of W and update it
*
         CALL DCOPY( K, A( 1, K ), 1, W( 1, KW ), 1 )
         IF( K.LT.N )
     $      CALL DGEMV( 'No transpose', K, N-K, -ONE, A( 1, K+1 ),
     $                  LDA, W( K, KW+1 ), LDW, ONE, W( 1, KW ), 1 )
*
*        Determine rows and columns to be interchanged and whether
*        a 1-by-1 or 2-by-2 pivot block will be used
*
         ABSAKK = ABS( W( K, KW ) )
*
*        IMAX is the row-index of the largest off-diagonal element in
*        column K, and COLMAX is its absolute value.
*        Determine both COLMAX and IMAX.
*
         IF( K.GT.1 ) THEN
            IMAX = IDAMAX( K-1, W( 1, KW ), 1 )
            COLMAX = ABS( W( IMAX, KW ) )
         ELSE
            COLMAX = ZERO
         END IF
*
         IF( MAX( ABSAKK, COLMAX ).EQ.ZERO ) THEN
*
*           Column K is zero or underflow: set INFO and continue
*
            IF( INFO.EQ.0 )
     $         INFO = K
            KP = K
            CALL DCOPY( K, W( 1, KW ), 1, A( 1, K ), 1 )
*
*           Set E( K ) to zero
*
            IF( K.GT.1 )
     $         E( K ) = ZERO
*
         ELSE
*
*           ============================================================
*
*           Test for interchange
*
*           Equivalent to testing for ABSAKK.GE.ALPHA*COLMAX
*           (used to handle NaN and Inf)
*
            IF( .NOT.( ABSAKK.LT.ALPHA*COLMAX ) ) THEN
*
*              no interchange, use 1-by-1 pivot block
*
               KP = K
*
            ELSE
*
               DONE = .FALSE.
*
*              Loop until pivot found
*
   12          CONTINUE
*
*                 Begin pivot search loop body
*
*
*                 Copy column IMAX to column KW-1 of W and update it
*
                  CALL DCOPY( IMAX, A( 1, IMAX ), 1, W( 1, KW-1 ), 1 )
                  CALL DCOPY( K-IMAX, A( IMAX, IMAX+1 ), LDA,
     $                        W( IMAX+1, KW-1 ), 1 )
*
                  IF( K.LT.N )
     $               CALL DGEMV( 'No transpose', K, N-K, -ONE,
     $                           A( 1, K+1 ), LDA, W( IMAX, KW+1 ), LDW,
     $                           ONE, W( 1, KW-1 ), 1 )
*
*                 JMAX is the column-index of the largest off-diagonal
*                 element in row IMAX, and ROWMAX is its absolute value.
*                 Determine both ROWMAX and JMAX.
*
                  IF( IMAX.NE.K ) THEN
                     JMAX = IMAX + IDAMAX( K-IMAX, W( IMAX+1, KW-1 ),
     $                                     1 )
                     ROWMAX = ABS( W( JMAX, KW-1 ) )
                  ELSE
                     ROWMAX = ZERO
                  END IF
*
                  IF( IMAX.GT.1 ) THEN
                     ITEMP = IDAMAX( IMAX-1, W( 1, KW-1 ), 1 )
                     DTEMP = ABS( W( ITEMP, KW-1 ) )
                     IF( DTEMP.GT.ROWMAX ) THEN
                        ROWMAX = DTEMP
                        JMAX = ITEMP
                     END IF
                  END IF
*
*                 Equivalent to testing for
*                 ABS( W( IMAX, KW-1 ) ).GE.ALPHA*ROWMAX
*                 (used to handle NaN and Inf)
*
                  IF( .NOT.(ABS( W( IMAX, KW-1 ) ).LT.ALPHA*ROWMAX ) )
     $            THEN
*
*                    interchange rows and columns K and IMAX,
*                    use 1-by-1 pivot block
*
                     KP = IMAX
*
*                    copy column KW-1 of W to column KW of W
*
                     CALL DCOPY( K, W( 1, KW-1 ), 1, W( 1, KW ), 1 )
*
                     DONE = .TRUE.
*
*                 Equivalent to testing for ROWMAX.EQ.COLMAX,
*                 (used to handle NaN and Inf)
*
                  ELSE IF( ( P.EQ.JMAX ) .OR. ( ROWMAX.LE.COLMAX ) )
     $            THEN
*
*                    interchange rows and columns K-1 and IMAX,
*                    use 2-by-2 pivot block
*
                     KP = IMAX
                     KSTEP = 2
                     DONE = .TRUE.
                  ELSE
*
*                    Pivot not found: set params and repeat
*
                     P = IMAX
                     COLMAX = ROWMAX
                     IMAX = JMAX
*
*                    Copy updated JMAXth (next IMAXth) column to Kth of W
*
                     CALL DCOPY( K, W( 1, KW-1 ), 1, W( 1, KW ), 1 )
*
                  END IF
*
*                 End pivot search loop body
*
               IF( .NOT. DONE ) GOTO 12
*
            END IF
*
*           ============================================================
*
            KK = K - KSTEP + 1
*
*           KKW is the column of W which corresponds to column KK of A
*
            KKW = NB + KK - N
*
            IF( ( KSTEP.EQ.2 ) .AND. ( P.NE.K ) ) THEN
*
*              Copy non-updated column K to column P
*
               CALL DCOPY( K-P, A( P+1, K ), 1, A( P, P+1 ), LDA )
               CALL DCOPY( P, A( 1, K ), 1, A( 1, P ), 1 )
*
*              Interchange rows K and P in last N-K+1 columns of A
*              and last N-K+2 columns of W
*
               CALL DSWAP( N-K+1, A( K, K ), LDA, A( P, K ), LDA )
               CALL DSWAP( N-KK+1, W( K, KKW ), LDW, W( P, KKW ), LDW )
            END IF
*
*           Updated column KP is already stored in column KKW of W
*
            IF( KP.NE.KK ) THEN
*
*              Copy non-updated column KK to column KP
*
               A( KP, K ) = A( KK, K )
               CALL DCOPY( K-1-KP, A( KP+1, KK ), 1, A( KP, KP+1 ),
     $                     LDA )
               CALL DCOPY( KP, A( 1, KK ), 1, A( 1, KP ), 1 )
*
*              Interchange rows KK and KP in last N-KK+1 columns
*              of A and W
*
               CALL DSWAP( N-KK+1, A( KK, KK ), LDA, A( KP, KK ), LDA )
               CALL DSWAP( N-KK+1, W( KK, KKW ), LDW, W( KP, KKW ),
     $                     LDW )
            END IF
*
            IF( KSTEP.EQ.1 ) THEN
*
*              1-by-1 pivot block D(k): column KW of W now holds
*
*              W(k) = U(k)*D(k)
*
*              where U(k) is the k-th column of U
*
*              Store U(k) in column k of A
*
               CALL DCOPY( K, W( 1, KW ), 1, A( 1, K ), 1 )
               IF( K.GT.1 ) THEN
                  IF( ABS( A( K, K ) ).GE.SFMIN ) THEN
                     R1 = ONE / A( K, K )
                     CALL DSCAL( K-1, R1, A( 1, K ), 1 )
                  ELSE IF( A( K, K ).NE.ZERO ) THEN
                     DO 14 II = 1, K - 1
                        A( II, K ) = A( II, K ) / A( K, K )
   14                CONTINUE
                  END IF
*
*                 Store the superdiagonal element of D in array E
*
                  E( K ) = ZERO
*
               END IF
*
            ELSE
*
*              2-by-2 pivot block D(k): columns KW and KW-1 of W now
*              hold
*
*              ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)
*
*              where U(k) and U(k-1) are the k-th and (k-1)-th columns
*              of U
*
               IF( K.GT.2 ) THEN
*
*                 Store U(k) and U(k-1) in columns k and k-1 of A
*
                  D12 = W( K-1, KW )
                  D11 = W( K, KW ) / D12
                  D22 = W( K-1, KW-1 ) / D12
                  T = ONE / ( D11*D22-ONE )
                  DO 20 J = 1, K - 2
                     A( J, K-1 ) = T*( (D11*W( J, KW-1 )-W( J, KW ) ) /
     $                             D12 )
                     A( J, K ) = T*( ( D22*W( J, KW )-W( J, KW-1 ) ) /
     $                           D12 )
   20             CONTINUE
               END IF
*
*              Copy diagonal elements of D(K) to A,
*              copy superdiagonal element of D(K) to E(K) and
*              ZERO out superdiagonal entry of A
*
               A( K-1, K-1 ) = W( K-1, KW-1 )
               A( K-1, K ) = ZERO
               A( K, K ) = W( K, KW )
               E( K ) = W( K-1, KW )
               E( K-1 ) = ZERO
*
            END IF
*
*           End column K is nonsingular
*
         END IF
*
*        Store details of the interchanges in IPIV
*
         IF( KSTEP.EQ.1 ) THEN
            IPIV( K ) = KP
         ELSE
            IPIV( K ) = -P
            IPIV( K-1 ) = -KP
         END IF
*
*        Decrease K and return to the start of the main loop
*
         K = K - KSTEP
         GO TO 10
*
   30    CONTINUE
*
*        Update the upper triangle of A11 (= A(1:k,1:k)) as
*
*        A11 := A11 - U12*D*U12**T = A11 - U12*W**T
*
*        computing blocks of NB columns at a time
*
         DO 50 J = ( ( K-1 ) / NB )*NB + 1, 1, -NB
            JB = MIN( NB, K-J+1 )
*
*           Update the upper triangle of the diagonal block
*
            DO 40 JJ = J, J + JB - 1
               CALL DGEMV( 'No transpose', JJ-J+1, N-K, -ONE,
     $                     A( J, K+1 ), LDA, W( JJ, KW+1 ), LDW, ONE,
     $                     A( J, JJ ), 1 )
   40       CONTINUE
*
*           Update the rectangular superdiagonal block
*
            IF( J.GE.2 )
     $         CALL DGEMM( 'No transpose', 'Transpose', J-1, JB,
     $                  N-K, -ONE, A( 1, K+1 ), LDA, W( J, KW+1 ),
     $                  LDW, ONE, A( 1, J ), LDA )
   50    CONTINUE
*
*        Set KB to the number of columns factorized
*
         KB = N - K
*
      ELSE
*
*        Factorize the leading columns of A using the lower triangle
*        of A and working forwards, and compute the matrix W = L21*D
*        for use in updating A22
*
*        Initialize the unused last entry of the subdiagonal array E.
*
         E( N ) = ZERO
*
*        K is the main loop index, increasing from 1 in steps of 1 or 2
*
         K = 1
   70   CONTINUE
*
*        Exit from loop
*
         IF( ( K.GE.NB .AND. NB.LT.N ) .OR. K.GT.N )
     $      GO TO 90
*
         KSTEP = 1
         P = K
*
*        Copy column K of A to column K of W and update it
*
         CALL DCOPY( N-K+1, A( K, K ), 1, W( K, K ), 1 )
         IF( K.GT.1 )
     $      CALL DGEMV( 'No transpose', N-K+1, K-1, -ONE, A( K, 1 ),
     $                  LDA, W( K, 1 ), LDW, ONE, W( K, K ), 1 )
*
*        Determine rows and columns to be interchanged and whether
*        a 1-by-1 or 2-by-2 pivot block will be used
*
         ABSAKK = ABS( W( K, K ) )
*
*        IMAX is the row-index of the largest off-diagonal element in
*        column K, and COLMAX is its absolute value.
*        Determine both COLMAX and IMAX.
*
         IF( K.LT.N ) THEN
            IMAX = K + IDAMAX( N-K, W( K+1, K ), 1 )
            COLMAX = ABS( W( IMAX, K ) )
         ELSE
            COLMAX = ZERO
         END IF
*
         IF( MAX( ABSAKK, COLMAX ).EQ.ZERO ) THEN
*
*           Column K is zero or underflow: set INFO and continue
*
            IF( INFO.EQ.0 )
     $         INFO = K
            KP = K
            CALL DCOPY( N-K+1, W( K, K ), 1, A( K, K ), 1 )
*
*           Set E( K ) to zero
*
            IF( K.LT.N )
     $         E( K ) = ZERO
*
         ELSE
*
*           ============================================================
*
*           Test for interchange
*
*           Equivalent to testing for ABSAKK.GE.ALPHA*COLMAX
*           (used to handle NaN and Inf)
*
            IF( .NOT.( ABSAKK.LT.ALPHA*COLMAX ) ) THEN
*
*              no interchange, use 1-by-1 pivot block
*
               KP = K
*
            ELSE
*
               DONE = .FALSE.
*
*              Loop until pivot found
*
   72          CONTINUE
*
*                 Begin pivot search loop body
*
*
*                 Copy column IMAX to column K+1 of W and update it
*
                  CALL DCOPY( IMAX-K, A( IMAX, K ), LDA, W( K, K+1 ), 1)
                  CALL DCOPY( N-IMAX+1, A( IMAX, IMAX ), 1,
     $                        W( IMAX, K+1 ), 1 )
                  IF( K.GT.1 )
     $               CALL DGEMV( 'No transpose', N-K+1, K-1, -ONE,
     $                           A( K, 1 ), LDA, W( IMAX, 1 ), LDW,
     $                           ONE, W( K, K+1 ), 1 )
*
*                 JMAX is the column-index of the largest off-diagonal
*                 element in row IMAX, and ROWMAX is its absolute value.
*                 Determine both ROWMAX and JMAX.
*
                  IF( IMAX.NE.K ) THEN
                     JMAX = K - 1 + IDAMAX( IMAX-K, W( K, K+1 ), 1 )
                     ROWMAX = ABS( W( JMAX, K+1 ) )
                  ELSE
                     ROWMAX = ZERO
                  END IF
*
                  IF( IMAX.LT.N ) THEN
                     ITEMP = IMAX + IDAMAX( N-IMAX, W( IMAX+1, K+1 ), 1)
                     DTEMP = ABS( W( ITEMP, K+1 ) )
                     IF( DTEMP.GT.ROWMAX ) THEN
                        ROWMAX = DTEMP
                        JMAX = ITEMP
                     END IF
                  END IF
*
*                 Equivalent to testing for
*                 ABS( W( IMAX, K+1 ) ).GE.ALPHA*ROWMAX
*                 (used to handle NaN and Inf)
*
                  IF( .NOT.( ABS( W( IMAX, K+1 ) ).LT.ALPHA*ROWMAX ) )
     $            THEN
*
*                    interchange rows and columns K and IMAX,
*                    use 1-by-1 pivot block
*
                     KP = IMAX
*
*                    copy column K+1 of W to column K of W
*
                     CALL DCOPY( N-K+1, W( K, K+1 ), 1, W( K, K ), 1 )
*
                     DONE = .TRUE.
*
*                 Equivalent to testing for ROWMAX.EQ.COLMAX,
*                 (used to handle NaN and Inf)
*
                  ELSE IF( ( P.EQ.JMAX ) .OR. ( ROWMAX.LE.COLMAX ) )
     $            THEN
*
*                    interchange rows and columns K+1 and IMAX,
*                    use 2-by-2 pivot block
*
                     KP = IMAX
                     KSTEP = 2
                     DONE = .TRUE.
                  ELSE
*
*                    Pivot not found: set params and repeat
*
                     P = IMAX
                     COLMAX = ROWMAX
                     IMAX = JMAX
*
*                    Copy updated JMAXth (next IMAXth) column to Kth of W
*
                     CALL DCOPY( N-K+1, W( K, K+1 ), 1, W( K, K ), 1 )
*
                  END IF
*
*                 End pivot search loop body
*
               IF( .NOT. DONE ) GOTO 72
*
            END IF
*
*           ============================================================
*
            KK = K + KSTEP - 1
*
            IF( ( KSTEP.EQ.2 ) .AND. ( P.NE.K ) ) THEN
*
*              Copy non-updated column K to column P
*
               CALL DCOPY( P-K, A( K, K ), 1, A( P, K ), LDA )
               CALL DCOPY( N-P+1, A( P, K ), 1, A( P, P ), 1 )
*
*              Interchange rows K and P in first K columns of A
*              and first K+1 columns of W
*
               CALL DSWAP( K, A( K, 1 ), LDA, A( P, 1 ), LDA )
               CALL DSWAP( KK, W( K, 1 ), LDW, W( P, 1 ), LDW )
            END IF
*
*           Updated column KP is already stored in column KK of W
*
            IF( KP.NE.KK ) THEN
*
*              Copy non-updated column KK to column KP
*
               A( KP, K ) = A( KK, K )
               CALL DCOPY( KP-K-1, A( K+1, KK ), 1, A( KP, K+1 ), LDA )
               CALL DCOPY( N-KP+1, A( KP, KK ), 1, A( KP, KP ), 1 )
*
*              Interchange rows KK and KP in first KK columns of A and W
*
               CALL DSWAP( KK, A( KK, 1 ), LDA, A( KP, 1 ), LDA )
               CALL DSWAP( KK, W( KK, 1 ), LDW, W( KP, 1 ), LDW )
            END IF
*
            IF( KSTEP.EQ.1 ) THEN
*
*              1-by-1 pivot block D(k): column k of W now holds
*
*              W(k) = L(k)*D(k)
*
*              where L(k) is the k-th column of L
*
*              Store L(k) in column k of A
*
               CALL DCOPY( N-K+1, W( K, K ), 1, A( K, K ), 1 )
               IF( K.LT.N ) THEN
                  IF( ABS( A( K, K ) ).GE.SFMIN ) THEN
                     R1 = ONE / A( K, K )
                     CALL DSCAL( N-K, R1, A( K+1, K ), 1 )
                  ELSE IF( A( K, K ).NE.ZERO ) THEN
                     DO 74 II = K + 1, N
                        A( II, K ) = A( II, K ) / A( K, K )
   74                CONTINUE
                  END IF
*
*                 Store the subdiagonal element of D in array E
*
                  E( K ) = ZERO
*
               END IF
*
            ELSE
*
*              2-by-2 pivot block D(k): columns k and k+1 of W now hold
*
*              ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k)
*
*              where L(k) and L(k+1) are the k-th and (k+1)-th columns
*              of L
*
               IF( K.LT.N-1 ) THEN
*
*                 Store L(k) and L(k+1) in columns k and k+1 of A
*
                  D21 = W( K+1, K )
                  D11 = W( K+1, K+1 ) / D21
                  D22 = W( K, K ) / D21
                  T = ONE / ( D11*D22-ONE )
                  DO 80 J = K + 2, N
                     A( J, K ) = T*( ( D11*W( J, K )-W( J, K+1 ) ) /
     $                           D21 )
                     A( J, K+1 ) = T*( ( D22*W( J, K+1 )-W( J, K ) ) /
     $                             D21 )
   80             CONTINUE
               END IF
*
*              Copy diagonal elements of D(K) to A,
*              copy subdiagonal element of D(K) to E(K) and
*              ZERO out subdiagonal entry of A
*
               A( K, K ) = W( K, K )
               A( K+1, K ) = ZERO
               A( K+1, K+1 ) = W( K+1, K+1 )
               E( K ) = W( K+1, K )
               E( K+1 ) = ZERO
*
            END IF
*
*           End column K is nonsingular
*
         END IF
*
*        Store details of the interchanges in IPIV
*
         IF( KSTEP.EQ.1 ) THEN
            IPIV( K ) = KP
         ELSE
            IPIV( K ) = -P
            IPIV( K+1 ) = -KP
         END IF
*
*        Increase K and return to the start of the main loop
*
         K = K + KSTEP
         GO TO 70
*
   90    CONTINUE
*
*        Update the lower triangle of A22 (= A(k:n,k:n)) as
*
*        A22 := A22 - L21*D*L21**T = A22 - L21*W**T
*
*        computing blocks of NB columns at a time
*
         DO 110 J = K, N, NB
            JB = MIN( NB, N-J+1 )
*
*           Update the lower triangle of the diagonal block
*
            DO 100 JJ = J, J + JB - 1
               CALL DGEMV( 'No transpose', J+JB-JJ, K-1, -ONE,
     $                     A( JJ, 1 ), LDA, W( JJ, 1 ), LDW, ONE,
     $                     A( JJ, JJ ), 1 )
  100       CONTINUE
*
*           Update the rectangular subdiagonal block
*
            IF( J+JB.LE.N )
     $         CALL DGEMM( 'No transpose', 'Transpose', N-J-JB+1, JB,
     $                     K-1, -ONE, A( J+JB, 1 ), LDA, W( J, 1 ),
     $                     LDW, ONE, A( J+JB, J ), LDA )
  110    CONTINUE
*
*        Set KB to the number of columns factorized
*
         KB = K - 1
*
      END IF
*
      RETURN
*
*     End of DLASYF_RK
*
      END

! DLAISNAN
      LOGICAL FUNCTION DLAISNAN( DIN1, DIN2 )
*
*  -- LAPACK auxiliary routine (version 3.7.1) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     June 2017
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION, INTENT(IN) :: DIN1, DIN2
*     ..
*
*  =====================================================================
*
*  .. Executable Statements ..
      DLAISNAN = (DIN1.NE.DIN2)
      RETURN
      END

! DSYR
      SUBROUTINE DSYR(UPLO,N,ALPHA,X,INCX,A,LDA)
*
*  -- Reference BLAS level2 routine (version 3.7.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA
      INTEGER INCX,LDA,N
      CHARACTER UPLO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),X(*)
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D+0)
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,IX,J,JX,KX
*     ..
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MAX
*     ..
*
*     Test the input parameters.
*
      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
          INFO = 1
      ELSE IF (N.LT.0) THEN
          INFO = 2
      ELSE IF (INCX.EQ.0) THEN
          INFO = 5
      ELSE IF (LDA.LT.MAX(1,N)) THEN
          INFO = 7
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DSYR  ',INFO)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((N.EQ.0) .OR. (ALPHA.EQ.ZERO)) RETURN
*
*     Set the start point in X if the increment is not unity.
*
      IF (INCX.LE.0) THEN
          KX = 1 - (N-1)*INCX
      ELSE IF (INCX.NE.1) THEN
          KX = 1
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through the triangular part
*     of A.
*
      IF (LSAME(UPLO,'U')) THEN
*
*        Form  A  when A is stored in upper triangle.
*
          IF (INCX.EQ.1) THEN
              DO 20 J = 1,N
                  IF (X(J).NE.ZERO) THEN
                      TEMP = ALPHA*X(J)
                      DO 10 I = 1,J
                          A(I,J) = A(I,J) + X(I)*TEMP
   10                 CONTINUE
                  END IF
   20         CONTINUE
          ELSE
              JX = KX
              DO 40 J = 1,N
                  IF (X(JX).NE.ZERO) THEN
                      TEMP = ALPHA*X(JX)
                      IX = KX
                      DO 30 I = 1,J
                          A(I,J) = A(I,J) + X(IX)*TEMP
                          IX = IX + INCX
   30                 CONTINUE
                  END IF
                  JX = JX + INCX
   40         CONTINUE
          END IF
      ELSE
*
*        Form  A  when A is stored in lower triangle.
*
          IF (INCX.EQ.1) THEN
              DO 60 J = 1,N
                  IF (X(J).NE.ZERO) THEN
                      TEMP = ALPHA*X(J)
                      DO 50 I = J,N
                          A(I,J) = A(I,J) + X(I)*TEMP
   50                 CONTINUE
                  END IF
   60         CONTINUE
          ELSE
              JX = KX
              DO 80 J = 1,N
                  IF (X(JX).NE.ZERO) THEN
                      TEMP = ALPHA*X(JX)
                      IX = JX
                      DO 70 I = J,N
                          A(I,J) = A(I,J) + X(IX)*TEMP
                          IX = IX + INCX
   70                 CONTINUE
                  END IF
                  JX = JX + INCX
   80         CONTINUE
          END IF
      END IF
*
      RETURN
*
*     End of DSYR  .
*
      END

! DLASRT
      SUBROUTINE DLASRT( ID, N, D, INFO )
*
*  -- LAPACK computational routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     June 2016
*
*     .. Scalar Arguments ..
      CHARACTER          ID
      INTEGER            INFO, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            SELECT
      PARAMETER          ( SELECT = 20 )
*     ..
*     .. Local Scalars ..
      INTEGER            DIR, ENDD, I, J, START, STKPNT
      DOUBLE PRECISION   D1, D2, D3, DMNMX, TMP
*     ..
*     .. Local Arrays ..
      INTEGER            STACK( 2, 32 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      DIR = -1
      IF( LSAME( ID, 'D' ) ) THEN
         DIR = 0
      ELSE IF( LSAME( ID, 'I' ) ) THEN
         DIR = 1
      END IF
      IF( DIR.EQ.-1 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLASRT', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.1 )
     $   RETURN
*
      STKPNT = 1
      STACK( 1, 1 ) = 1
      STACK( 2, 1 ) = N
   10 CONTINUE
      START = STACK( 1, STKPNT )
      ENDD = STACK( 2, STKPNT )
      STKPNT = STKPNT - 1
      IF( ENDD-START.LE.SELECT .AND. ENDD-START.GT.0 ) THEN
*
*        Do Insertion sort on D( START:ENDD )
*
         IF( DIR.EQ.0 ) THEN
*
*           Sort into decreasing order
*
            DO 30 I = START + 1, ENDD
               DO 20 J = I, START + 1, -1
                  IF( D( J ).GT.D( J-1 ) ) THEN
                     DMNMX = D( J )
                     D( J ) = D( J-1 )
                     D( J-1 ) = DMNMX
                  ELSE
                     GO TO 30
                  END IF
   20          CONTINUE
   30       CONTINUE
*
         ELSE
*
*           Sort into increasing order
*
            DO 50 I = START + 1, ENDD
               DO 40 J = I, START + 1, -1
                  IF( D( J ).LT.D( J-1 ) ) THEN
                     DMNMX = D( J )
                     D( J ) = D( J-1 )
                     D( J-1 ) = DMNMX
                  ELSE
                     GO TO 50
                  END IF
   40          CONTINUE
   50       CONTINUE
*
         END IF
*
      ELSE IF( ENDD-START.GT.SELECT ) THEN
*
*        Partition D( START:ENDD ) and stack parts, largest one first
*
*        Choose partition entry as median of 3
*
         D1 = D( START )
         D2 = D( ENDD )
         I = ( START+ENDD ) / 2
         D3 = D( I )
         IF( D1.LT.D2 ) THEN
            IF( D3.LT.D1 ) THEN
               DMNMX = D1
            ELSE IF( D3.LT.D2 ) THEN
               DMNMX = D3
            ELSE
               DMNMX = D2
            END IF
         ELSE
            IF( D3.LT.D2 ) THEN
               DMNMX = D2
            ELSE IF( D3.LT.D1 ) THEN
               DMNMX = D3
            ELSE
               DMNMX = D1
            END IF
         END IF
*
         IF( DIR.EQ.0 ) THEN
*
*           Sort into decreasing order
*
            I = START - 1
            J = ENDD + 1
   60       CONTINUE
   70       CONTINUE
            J = J - 1
            IF( D( J ).LT.DMNMX )
     $         GO TO 70
   80       CONTINUE
            I = I + 1
            IF( D( I ).GT.DMNMX )
     $         GO TO 80
            IF( I.LT.J ) THEN
               TMP = D( I )
               D( I ) = D( J )
               D( J ) = TMP
               GO TO 60
            END IF
            IF( J-START.GT.ENDD-J-1 ) THEN
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
            ELSE
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
            END IF
         ELSE
*
*           Sort into increasing order
*
            I = START - 1
            J = ENDD + 1
   90       CONTINUE
  100       CONTINUE
            J = J - 1
            IF( D( J ).GT.DMNMX )
     $         GO TO 100
  110       CONTINUE
            I = I + 1
            IF( D( I ).LT.DMNMX )
     $         GO TO 110
            IF( I.LT.J ) THEN
               TMP = D( I )
               D( I ) = D( J )
               D( J ) = TMP
               GO TO 90
            END IF
            IF( J-START.GT.ENDD-J-1 ) THEN
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
            ELSE
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
            END IF
         END IF
      END IF
      IF( STKPNT.GT.0 )
     $   GO TO 10
      RETURN
*
*     End of DLASRT
*
      END

! DTRTRS
      SUBROUTINE DTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB,
     $                   INFO )
*
*  -- LAPACK computational routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, TRANS, UPLO
      INTEGER            INFO, LDA, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOUNIT
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DTRSM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      NOUNIT = LSAME( DIAG, 'N' )
      IF( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LSAME( TRANS, 'N' ) .AND. .NOT.
     $         LSAME( TRANS, 'T' ) .AND. .NOT.LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DTRTRS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Check for singularity.
*
      IF( NOUNIT ) THEN
         DO 10 INFO = 1, N
            IF( A( INFO, INFO ).EQ.ZERO )
     $         RETURN
   10    CONTINUE
      END IF
      INFO = 0
*
*     Solve A * x = b  or  A**T * x = b.
*
      CALL DTRSM( 'Left', UPLO, TRANS, DIAG, N, NRHS, ONE, A, LDA, B,
     $            LDB )
*
      RETURN
*
*     End of DTRTRS
*
      END      

! DPOTRF
      SUBROUTINE dpotrf( UPLO, N, A, LDA, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      parameter( one = 1.0d+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J, JB, NB
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           lsame, ilaenv
*     ..
*     .. External Subroutines ..
      EXTERNAL           dgemm, dpotrf2, dsyrk, dtrsm, xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          max, min
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      info = 0
      upper = lsame( uplo, 'U' )
      IF( .NOT.upper .AND. .NOT.lsame( uplo, 'L' ) ) THEN
         info = -1
      ELSE IF( n.LT.0 ) THEN
         info = -2
      ELSE IF( lda.LT.max( 1, n ) ) THEN
         info = -4
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'DPOTRF', -info )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( n.EQ.0 )
     $   RETURN
*
*     Determine the block size for this environment.
*
      nb = ilaenv( 1, 'DPOTRF', uplo, n, -1, -1, -1 )
      IF( nb.LE.1 .OR. nb.GE.n ) THEN
*
*        Use unblocked code.
*
         CALL dpotrf2( uplo, n, a, lda, info )
      ELSE
*
*        Use blocked code.
*
         IF( upper ) THEN
*
*           Compute the Cholesky factorization A = U**T*U.
*
            DO 10 j = 1, n, nb
*
*              Update and factorize the current diagonal block and test
*              for non-positive-definiteness.
*
               jb = min( nb, n-j+1 )
               CALL dsyrk( 'Upper', 'Transpose', jb, j-1, -one,
     $                     a( 1, j ), lda, one, a( j, j ), lda )
               CALL dpotrf2( 'Upper', jb, a( j, j ), lda, info )
               IF( info.NE.0 )
     $            GO TO 30
               IF( j+jb.LE.n ) THEN
*
*                 Compute the current block row.
*
                  CALL dgemm( 'Transpose', 'No transpose', jb, n-j-jb+1,
     $                        j-1, -one, a( 1, j ), lda, a( 1, j+jb ),
     $                        lda, one, a( j, j+jb ), lda )
                  CALL dtrsm( 'Left', 'Upper', 'Transpose', 'Non-unit',
     $                        jb, n-j-jb+1, one, a( j, j ), lda,
     $                        a( j, j+jb ), lda )
               END IF
   10       CONTINUE
*
         ELSE
*
*           Compute the Cholesky factorization A = L*L**T.
*
            DO 20 j = 1, n, nb
*
*              Update and factorize the current diagonal block and test
*              for non-positive-definiteness.
*
               jb = min( nb, n-j+1 )
               CALL dsyrk( 'Lower', 'No transpose', jb, j-1, -one,
     $                     a( j, 1 ), lda, one, a( j, j ), lda )
               CALL dpotrf2( 'Lower', jb, a( j, j ), lda, info )
               IF( info.NE.0 )
     $            GO TO 30
               IF( j+jb.LE.n ) THEN
*
*                 Compute the current block column.
*
                  CALL dgemm( 'No transpose', 'Transpose', n-j-jb+1, jb,
     $                        j-1, -one, a( j+jb, 1 ), lda, a( j, 1 ),
     $                        lda, one, a( j+jb, j ), lda )
                  CALL dtrsm( 'Right', 'Lower', 'Transpose', 'Non-unit',
     $                        n-j-jb+1, jb, one, a( j, j ), lda,
     $                        a( j+jb, j ), lda )
               END IF
   20       CONTINUE
         END IF
      END IF
      GO TO 40
*
   30 CONTINUE
      info = info + j - 1
*
   40 CONTINUE
      RETURN
*
*     End of DPOTRF
*
      END
      
! DSYRK
      SUBROUTINE dsyrk(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
*
*  -- Reference BLAS level3 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA,BETA
      INTEGER K,LDA,LDC,N
      CHARACTER TRANS,UPLO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),C(LDC,*)
*     ..
*
*  =====================================================================
*
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL lsame
*     ..
*     .. External Subroutines ..
      EXTERNAL xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC max
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,J,L,NROWA
      LOGICAL UPPER
*     ..
*     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      parameter(one=1.0d+0,zero=0.0d+0)
*     ..
*
*     Test the input parameters.
*
      IF (lsame(trans,'N')) THEN
          nrowa = n
      ELSE
          nrowa = k
      END IF
      upper = lsame(uplo,'U')
*
      info = 0
      IF ((.NOT.upper) .AND. (.NOT.lsame(uplo,'L'))) THEN
          info = 1
      ELSE IF ((.NOT.lsame(trans,'N')) .AND.
     +         (.NOT.lsame(trans,'T')) .AND.
     +         (.NOT.lsame(trans,'C'))) THEN
          info = 2
      ELSE IF (n.LT.0) THEN
          info = 3
      ELSE IF (k.LT.0) THEN
          info = 4
      ELSE IF (lda.LT.max(1,nrowa)) THEN
          info = 7
      ELSE IF (ldc.LT.max(1,n)) THEN
          info = 10
      END IF
      IF (info.NE.0) THEN
          CALL xerbla('DSYRK ',info)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((n.EQ.0) .OR. (((alpha.EQ.zero).OR.
     +    (k.EQ.0)).AND. (beta.EQ.one))) RETURN
*
*     And when  alpha.eq.zero.
*
      IF (alpha.EQ.zero) THEN
          IF (upper) THEN
              IF (beta.EQ.zero) THEN
                  DO 20 j = 1,n
                      DO 10 i = 1,j
                          c(i,j) = zero
   10                 CONTINUE
   20             CONTINUE
              ELSE
                  DO 40 j = 1,n
                      DO 30 i = 1,j
                          c(i,j) = beta*c(i,j)
   30                 CONTINUE
   40             CONTINUE
              END IF
          ELSE
              IF (beta.EQ.zero) THEN
                  DO 60 j = 1,n
                      DO 50 i = j,n
                          c(i,j) = zero
   50                 CONTINUE
   60             CONTINUE
              ELSE
                  DO 80 j = 1,n
                      DO 70 i = j,n
                          c(i,j) = beta*c(i,j)
   70                 CONTINUE
   80             CONTINUE
              END IF
          END IF
          RETURN
      END IF
*
*     Start the operations.
*
      IF (lsame(trans,'N')) THEN
*
*        Form  C := alpha*A*A**T + beta*C.
*
          IF (upper) THEN
              DO 130 j = 1,n
                  IF (beta.EQ.zero) THEN
                      DO 90 i = 1,j
                          c(i,j) = zero
   90                 CONTINUE
                  ELSE IF (beta.NE.one) THEN
                      DO 100 i = 1,j
                          c(i,j) = beta*c(i,j)
  100                 CONTINUE
                  END IF
                  DO 120 l = 1,k
                      IF (a(j,l).NE.zero) THEN
                          temp = alpha*a(j,l)
                          DO 110 i = 1,j
                              c(i,j) = c(i,j) + temp*a(i,l)
  110                     CONTINUE
                      END IF
  120             CONTINUE
  130         CONTINUE
          ELSE
              DO 180 j = 1,n
                  IF (beta.EQ.zero) THEN
                      DO 140 i = j,n
                          c(i,j) = zero
  140                 CONTINUE
                  ELSE IF (beta.NE.one) THEN
                      DO 150 i = j,n
                          c(i,j) = beta*c(i,j)
  150                 CONTINUE
                  END IF
                  DO 170 l = 1,k
                      IF (a(j,l).NE.zero) THEN
                          temp = alpha*a(j,l)
                          DO 160 i = j,n
                              c(i,j) = c(i,j) + temp*a(i,l)
  160                     CONTINUE
                      END IF
  170             CONTINUE
  180         CONTINUE
          END IF
      ELSE
*
*        Form  C := alpha*A**T*A + beta*C.
*
          IF (upper) THEN
              DO 210 j = 1,n
                  DO 200 i = 1,j
                      temp = zero
                      DO 190 l = 1,k
                          temp = temp + a(l,i)*a(l,j)
  190                 CONTINUE
                      IF (beta.EQ.zero) THEN
                          c(i,j) = alpha*temp
                      ELSE
                          c(i,j) = alpha*temp + beta*c(i,j)
                      END IF
  200             CONTINUE
  210         CONTINUE
          ELSE
              DO 240 j = 1,n
                  DO 230 i = j,n
                      temp = zero
                      DO 220 l = 1,k
                          temp = temp + a(l,i)*a(l,j)
  220                 CONTINUE
                      IF (beta.EQ.zero) THEN
                          c(i,j) = alpha*temp
                      ELSE
                          c(i,j) = alpha*temp + beta*c(i,j)
                      END IF
  230             CONTINUE
  240         CONTINUE
          END IF
      END IF
*
      RETURN
*
*     End of DSYRK
*
      END

! DPOTRF2
      RECURSIVE SUBROUTINE dpotrf2( UPLO, N, A, LDA, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          uplo
      INTEGER            info, lda, n
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   a( lda, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   one, zero
      parameter( one = 1.0d+0, zero = 0.0d+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            upper
      INTEGER            n1, n2, iinfo
*     ..
*     .. External Functions ..
      LOGICAL            lsame, disnan
      EXTERNAL           lsame, disnan
*     ..
*     .. External Subroutines ..
      EXTERNAL           dsyrk, dtrsm, xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          max, sqrt
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      info = 0
      upper = lsame( uplo, 'U' )
      IF( .NOT.upper .AND. .NOT.lsame( uplo, 'L' ) ) THEN
         info = -1
      ELSE IF( n.LT.0 ) THEN
         info = -2
      ELSE IF( lda.LT.max( 1, n ) ) THEN
         info = -4
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'DPOTRF2', -info )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( n.EQ.0 )
     $   RETURN
*
*     N=1 case
*
      IF( n.EQ.1 ) THEN
*
*        Test for non-positive-definiteness
*
         IF( a( 1, 1 ).LE.zero.OR.disnan( a( 1, 1 ) ) ) THEN
            info = 1
            RETURN
         END IF
*
*        Factor
*
         a( 1, 1 ) = sqrt( a( 1, 1 ) )
*
*     Use recursive code
*
      ELSE
         n1 = n/2
         n2 = n-n1
*
*        Factor A11
*
         CALL dpotrf2( uplo, n1, a( 1, 1 ), lda, iinfo )
         IF ( iinfo.NE.0 ) THEN
            info = iinfo
            RETURN
         END IF
*
*        Compute the Cholesky factorization A = U**T*U
*
         IF( upper ) THEN
*
*           Update and scale A12
*
            CALL dtrsm( 'L', 'U', 'T', 'N', n1, n2, one,
     $                  a( 1, 1 ), lda, a( 1, n1+1 ), lda )
*
*           Update and factor A22
*
            CALL dsyrk( uplo, 'T', n2, n1, -one, a( 1, n1+1 ), lda,
     $                  one, a( n1+1, n1+1 ), lda )
            CALL dpotrf2( uplo, n2, a( n1+1, n1+1 ), lda, iinfo )
            IF ( iinfo.NE.0 ) THEN
               info = iinfo + n1
               RETURN
            END IF
*
*        Compute the Cholesky factorization A = L*L**T
*
         ELSE
*
*           Update and scale A21
*
            CALL dtrsm( 'R', 'L', 'T', 'N', n2, n1, one,
     $                  a( 1, 1 ), lda, a( n1+1, 1 ), lda )
*
*           Update and factor A22
*
            CALL dsyrk( uplo, 'N', n2, n1, -one, a( n1+1, 1 ), lda,
     $                  one, a( n1+1, n1+1 ), lda )
            CALL dpotrf2( uplo, n2, a( n1+1, n1+1 ), lda, iinfo )
            IF ( iinfo.NE.0 ) THEN
               info = iinfo + n1
               RETURN
            END IF
         END IF
      END IF
      RETURN
*
*     End of DPOTRF2
*
      END
      
! DSYGST
      SUBROUTINE dsygst( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, ITYPE, LDA, LDB, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, HALF
      parameter( one = 1.0d0, half = 0.5d0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            K, KB, NB
*     ..
*     .. External Subroutines ..
      EXTERNAL           dsygs2, dsymm, dsyr2k, dtrmm, dtrsm, xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          max, min
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           lsame, ilaenv
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      info = 0
      upper = lsame( uplo, 'U' )
      IF( itype.LT.1 .OR. itype.GT.3 ) THEN
         info = -1
      ELSE IF( .NOT.upper .AND. .NOT.lsame( uplo, 'L' ) ) THEN
         info = -2
      ELSE IF( n.LT.0 ) THEN
         info = -3
      ELSE IF( lda.LT.max( 1, n ) ) THEN
         info = -5
      ELSE IF( ldb.LT.max( 1, n ) ) THEN
         info = -7
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'DSYGST', -info )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( n.EQ.0 )
     $   RETURN
*
*     Determine the block size for this environment.
*
      nb = ilaenv( 1, 'DSYGST', uplo, n, -1, -1, -1 )
*
      IF( nb.LE.1 .OR. nb.GE.n ) THEN
*
*        Use unblocked code
*
         CALL dsygs2( itype, uplo, n, a, lda, b, ldb, info )
      ELSE
*
*        Use blocked code
*
         IF( itype.EQ.1 ) THEN
            IF( upper ) THEN
*
*              Compute inv(U**T)*A*inv(U)
*
               DO 10 k = 1, n, nb
                  kb = min( n-k+1, nb )
*
*                 Update the upper triangle of A(k:n,k:n)
*
                  CALL dsygs2( itype, uplo, kb, a( k, k ), lda,
     $                         b( k, k ), ldb, info )
                  IF( k+kb.LE.n ) THEN
                     CALL dtrsm( 'Left', uplo, 'Transpose', 'Non-unit',
     $                           kb, n-k-kb+1, one, b( k, k ), ldb,
     $                           a( k, k+kb ), lda )
                     CALL dsymm( 'Left', uplo, kb, n-k-kb+1, -half,
     $                           a( k, k ), lda, b( k, k+kb ), ldb, one,
     $                           a( k, k+kb ), lda )
                     CALL dsyr2k( uplo, 'Transpose', n-k-kb+1, kb, -one,
     $                            a( k, k+kb ), lda, b( k, k+kb ), ldb,
     $                            one, a( k+kb, k+kb ), lda )
                     CALL dsymm( 'Left', uplo, kb, n-k-kb+1, -half,
     $                           a( k, k ), lda, b( k, k+kb ), ldb, one,
     $                           a( k, k+kb ), lda )
                     CALL dtrsm( 'Right', uplo, 'No transpose',
     $                           'Non-unit', kb, n-k-kb+1, one,
     $                           b( k+kb, k+kb ), ldb, a( k, k+kb ),
     $                           lda )
                  END IF
   10          CONTINUE
            ELSE
*
*              Compute inv(L)*A*inv(L**T)
*
               DO 20 k = 1, n, nb
                  kb = min( n-k+1, nb )
*
*                 Update the lower triangle of A(k:n,k:n)
*
                  CALL dsygs2( itype, uplo, kb, a( k, k ), lda,
     $                         b( k, k ), ldb, info )
                  IF( k+kb.LE.n ) THEN
                     CALL dtrsm( 'Right', uplo, 'Transpose', 'Non-unit',
     $                           n-k-kb+1, kb, one, b( k, k ), ldb,
     $                           a( k+kb, k ), lda )
                     CALL dsymm( 'Right', uplo, n-k-kb+1, kb, -half,
     $                           a( k, k ), lda, b( k+kb, k ), ldb, one,
     $                           a( k+kb, k ), lda )
                     CALL dsyr2k( uplo, 'No transpose', n-k-kb+1, kb,
     $                            -one, a( k+kb, k ), lda, b( k+kb, k ),
     $                            ldb, one, a( k+kb, k+kb ), lda )
                     CALL dsymm( 'Right', uplo, n-k-kb+1, kb, -half,
     $                           a( k, k ), lda, b( k+kb, k ), ldb, one,
     $                           a( k+kb, k ), lda )
                     CALL dtrsm( 'Left', uplo, 'No transpose',
     $                           'Non-unit', n-k-kb+1, kb, one,
     $                           b( k+kb, k+kb ), ldb, a( k+kb, k ),
     $                           lda )
                  END IF
   20          CONTINUE
            END IF
         ELSE
            IF( upper ) THEN
*
*              Compute U*A*U**T
*
               DO 30 k = 1, n, nb
                  kb = min( n-k+1, nb )
*
*                 Update the upper triangle of A(1:k+kb-1,1:k+kb-1)
*
                  CALL dtrmm( 'Left', uplo, 'No transpose', 'Non-unit',
     $                        k-1, kb, one, b, ldb, a( 1, k ), lda )
                  CALL dsymm( 'Right', uplo, k-1, kb, half, a( k, k ),
     $                        lda, b( 1, k ), ldb, one, a( 1, k ), lda )
                  CALL dsyr2k( uplo, 'No transpose', k-1, kb, one,
     $                         a( 1, k ), lda, b( 1, k ), ldb, one, a,
     $                         lda )
                  CALL dsymm( 'Right', uplo, k-1, kb, half, a( k, k ),
     $                        lda, b( 1, k ), ldb, one, a( 1, k ), lda )
                  CALL dtrmm( 'Right', uplo, 'Transpose', 'Non-unit',
     $                        k-1, kb, one, b( k, k ), ldb, a( 1, k ),
     $                        lda )
                  CALL dsygs2( itype, uplo, kb, a( k, k ), lda,
     $                         b( k, k ), ldb, info )
   30          CONTINUE
            ELSE
*
*              Compute L**T*A*L
*
               DO 40 k = 1, n, nb
                  kb = min( n-k+1, nb )
*
*                 Update the lower triangle of A(1:k+kb-1,1:k+kb-1)
*
                  CALL dtrmm( 'Right', uplo, 'No transpose', 'Non-unit',
     $                        kb, k-1, one, b, ldb, a( k, 1 ), lda )
                  CALL dsymm( 'Left', uplo, kb, k-1, half, a( k, k ),
     $                        lda, b( k, 1 ), ldb, one, a( k, 1 ), lda )
                  CALL dsyr2k( uplo, 'Transpose', k-1, kb, one,
     $                         a( k, 1 ), lda, b( k, 1 ), ldb, one, a,
     $                         lda )
                  CALL dsymm( 'Left', uplo, kb, k-1, half, a( k, k ),
     $                        lda, b( k, 1 ), ldb, one, a( k, 1 ), lda )
                  CALL dtrmm( 'Left', uplo, 'Transpose', 'Non-unit', kb,
     $                        k-1, one, b( k, k ), ldb, a( k, 1 ), lda )
                  CALL dsygs2( itype, uplo, kb, a( k, k ), lda,
     $                         b( k, k ), ldb, info )
   40          CONTINUE
            END IF
         END IF
      END IF
      RETURN
*
*     End of DSYGST
*
      END

! DSYMM
      SUBROUTINE dsymm(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
*
*  -- Reference BLAS level3 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA,BETA
      INTEGER LDA,LDB,LDC,M,N
      CHARACTER SIDE,UPLO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
*     ..
*
*  =====================================================================
*
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL lsame
*     ..
*     .. External Subroutines ..
      EXTERNAL xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC max
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION TEMP1,TEMP2
      INTEGER I,INFO,J,K,NROWA
      LOGICAL UPPER
*     ..
*     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      parameter(one=1.0d+0,zero=0.0d+0)
*     ..
*
*     Set NROWA as the number of rows of A.
*
      IF (lsame(side,'L')) THEN
          nrowa = m
      ELSE
          nrowa = n
      END IF
      upper = lsame(uplo,'U')
*
*     Test the input parameters.
*
      info = 0
      IF ((.NOT.lsame(side,'L')) .AND. (.NOT.lsame(side,'R'))) THEN
          info = 1
      ELSE IF ((.NOT.upper) .AND. (.NOT.lsame(uplo,'L'))) THEN
          info = 2
      ELSE IF (m.LT.0) THEN
          info = 3
      ELSE IF (n.LT.0) THEN
          info = 4
      ELSE IF (lda.LT.max(1,nrowa)) THEN
          info = 7
      ELSE IF (ldb.LT.max(1,m)) THEN
          info = 9
      ELSE IF (ldc.LT.max(1,m)) THEN
          info = 12
      END IF
      IF (info.NE.0) THEN
          CALL xerbla('DSYMM ',info)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((m.EQ.0) .OR. (n.EQ.0) .OR.
     +    ((alpha.EQ.zero).AND. (beta.EQ.one))) RETURN
*
*     And when  alpha.eq.zero.
*
      IF (alpha.EQ.zero) THEN
          IF (beta.EQ.zero) THEN
              DO 20 j = 1,n
                  DO 10 i = 1,m
                      c(i,j) = zero
   10             CONTINUE
   20         CONTINUE
          ELSE
              DO 40 j = 1,n
                  DO 30 i = 1,m
                      c(i,j) = beta*c(i,j)
   30             CONTINUE
   40         CONTINUE
          END IF
          RETURN
      END IF
*
*     Start the operations.
*
      IF (lsame(side,'L')) THEN
*
*        Form  C := alpha*A*B + beta*C.
*
          IF (upper) THEN
              DO 70 j = 1,n
                  DO 60 i = 1,m
                      temp1 = alpha*b(i,j)
                      temp2 = zero
                      DO 50 k = 1,i - 1
                          c(k,j) = c(k,j) + temp1*a(k,i)
                          temp2 = temp2 + b(k,j)*a(k,i)
   50                 CONTINUE
                      IF (beta.EQ.zero) THEN
                          c(i,j) = temp1*a(i,i) + alpha*temp2
                      ELSE
                          c(i,j) = beta*c(i,j) + temp1*a(i,i) +
     +                             alpha*temp2
                      END IF
   60             CONTINUE
   70         CONTINUE
          ELSE
              DO 100 j = 1,n
                  DO 90 i = m,1,-1
                      temp1 = alpha*b(i,j)
                      temp2 = zero
                      DO 80 k = i + 1,m
                          c(k,j) = c(k,j) + temp1*a(k,i)
                          temp2 = temp2 + b(k,j)*a(k,i)
   80                 CONTINUE
                      IF (beta.EQ.zero) THEN
                          c(i,j) = temp1*a(i,i) + alpha*temp2
                      ELSE
                          c(i,j) = beta*c(i,j) + temp1*a(i,i) +
     +                             alpha*temp2
                      END IF
   90             CONTINUE
  100         CONTINUE
          END IF
      ELSE
*
*        Form  C := alpha*B*A + beta*C.
*
          DO 170 j = 1,n
              temp1 = alpha*a(j,j)
              IF (beta.EQ.zero) THEN
                  DO 110 i = 1,m
                      c(i,j) = temp1*b(i,j)
  110             CONTINUE
              ELSE
                  DO 120 i = 1,m
                      c(i,j) = beta*c(i,j) + temp1*b(i,j)
  120             CONTINUE
              END IF
              DO 140 k = 1,j - 1
                  IF (upper) THEN
                      temp1 = alpha*a(k,j)
                  ELSE
                      temp1 = alpha*a(j,k)
                  END IF
                  DO 130 i = 1,m
                      c(i,j) = c(i,j) + temp1*b(i,k)
  130             CONTINUE
  140         CONTINUE
              DO 160 k = j + 1,n
                  IF (upper) THEN
                      temp1 = alpha*a(j,k)
                  ELSE
                      temp1 = alpha*a(k,j)
                  END IF
                  DO 150 i = 1,m
                      c(i,j) = c(i,j) + temp1*b(i,k)
  150             CONTINUE
  160         CONTINUE
  170     CONTINUE
      END IF
*
      RETURN
*
*     End of DSYMM
*
      END

! DSYGS2
      SUBROUTINE dsygs2( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, ITYPE, LDA, LDB, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, HALF
      parameter( one = 1.0d0, half = 0.5d0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            K
      DOUBLE PRECISION   AKK, BKK, CT
*     ..
*     .. External Subroutines ..
      EXTERNAL           daxpy, dscal, dsyr2, dtrmv, dtrsv, xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          max
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           lsame
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      info = 0
      upper = lsame( uplo, 'U' )
      IF( itype.LT.1 .OR. itype.GT.3 ) THEN
         info = -1
      ELSE IF( .NOT.upper .AND. .NOT.lsame( uplo, 'L' ) ) THEN
         info = -2
      ELSE IF( n.LT.0 ) THEN
         info = -3
      ELSE IF( lda.LT.max( 1, n ) ) THEN
         info = -5
      ELSE IF( ldb.LT.max( 1, n ) ) THEN
         info = -7
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'DSYGS2', -info )
         RETURN
      END IF
*
      IF( itype.EQ.1 ) THEN
         IF( upper ) THEN
*
*           Compute inv(U**T)*A*inv(U)
*
            DO 10 k = 1, n
*
*              Update the upper triangle of A(k:n,k:n)
*
               akk = a( k, k )
               bkk = b( k, k )
               akk = akk / bkk**2
               a( k, k ) = akk
               IF( k.LT.n ) THEN
                  CALL dscal( n-k, one / bkk, a( k, k+1 ), lda )
                  ct = -half*akk
                  CALL daxpy( n-k, ct, b( k, k+1 ), ldb, a( k, k+1 ),
     $                        lda )
                  CALL dsyr2( uplo, n-k, -one, a( k, k+1 ), lda,
     $                        b( k, k+1 ), ldb, a( k+1, k+1 ), lda )
                  CALL daxpy( n-k, ct, b( k, k+1 ), ldb, a( k, k+1 ),
     $                        lda )
                  CALL dtrsv( uplo, 'Transpose', 'Non-unit', n-k,
     $                        b( k+1, k+1 ), ldb, a( k, k+1 ), lda )
               END IF
   10       CONTINUE
         ELSE
*
*           Compute inv(L)*A*inv(L**T)
*
            DO 20 k = 1, n
*
*              Update the lower triangle of A(k:n,k:n)
*
               akk = a( k, k )
               bkk = b( k, k )
               akk = akk / bkk**2
               a( k, k ) = akk
               IF( k.LT.n ) THEN
                  CALL dscal( n-k, one / bkk, a( k+1, k ), 1 )
                  ct = -half*akk
                  CALL daxpy( n-k, ct, b( k+1, k ), 1, a( k+1, k ), 1 )
                  CALL dsyr2( uplo, n-k, -one, a( k+1, k ), 1,
     $                        b( k+1, k ), 1, a( k+1, k+1 ), lda )
                  CALL daxpy( n-k, ct, b( k+1, k ), 1, a( k+1, k ), 1 )
                  CALL dtrsv( uplo, 'No transpose', 'Non-unit', n-k,
     $                        b( k+1, k+1 ), ldb, a( k+1, k ), 1 )
               END IF
   20       CONTINUE
         END IF
      ELSE
         IF( upper ) THEN
*
*           Compute U*A*U**T
*
            DO 30 k = 1, n
*
*              Update the upper triangle of A(1:k,1:k)
*
               akk = a( k, k )
               bkk = b( k, k )
               CALL dtrmv( uplo, 'No transpose', 'Non-unit', k-1, b,
     $                     ldb, a( 1, k ), 1 )
               ct = half*akk
               CALL daxpy( k-1, ct, b( 1, k ), 1, a( 1, k ), 1 )
               CALL dsyr2( uplo, k-1, one, a( 1, k ), 1, b( 1, k ), 1,
     $                     a, lda )
               CALL daxpy( k-1, ct, b( 1, k ), 1, a( 1, k ), 1 )
               CALL dscal( k-1, bkk, a( 1, k ), 1 )
               a( k, k ) = akk*bkk**2
   30       CONTINUE
         ELSE
*
*           Compute L**T *A*L
*
            DO 40 k = 1, n
*
*              Update the lower triangle of A(1:k,1:k)
*
               akk = a( k, k )
               bkk = b( k, k )
               CALL dtrmv( uplo, 'Transpose', 'Non-unit', k-1, b, ldb,
     $                     a( k, 1 ), lda )
               ct = half*akk
               CALL daxpy( k-1, ct, b( k, 1 ), ldb, a( k, 1 ), lda )
               CALL dsyr2( uplo, k-1, one, a( k, 1 ), lda, b( k, 1 ),
     $                     ldb, a, lda )
               CALL daxpy( k-1, ct, b( k, 1 ), ldb, a( k, 1 ), lda )
               CALL dscal( k-1, bkk, a( k, 1 ), lda )
               a( k, k ) = akk*bkk**2
   40       CONTINUE
         END IF
      END IF
      RETURN
*
*     End of DSYGS2
*
      END

! DTRSV
      SUBROUTINE dtrsv(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
*
*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER INCX,LDA,N
      CHARACTER DIAG,TRANS,UPLO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),X(*)
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION ZERO
      parameter(zero=0.0d+0)
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,IX,J,JX,KX
      LOGICAL NOUNIT
*     ..
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL lsame
*     ..
*     .. External Subroutines ..
      EXTERNAL xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC max
*     ..
*
*     Test the input parameters.
*
      info = 0
      IF (.NOT.lsame(uplo,'U') .AND. .NOT.lsame(uplo,'L')) THEN
          info = 1
      ELSE IF (.NOT.lsame(trans,'N') .AND. .NOT.lsame(trans,'T') .AND.
     +         .NOT.lsame(trans,'C')) THEN
          info = 2
      ELSE IF (.NOT.lsame(diag,'U') .AND. .NOT.lsame(diag,'N')) THEN
          info = 3
      ELSE IF (n.LT.0) THEN
          info = 4
      ELSE IF (lda.LT.max(1,n)) THEN
          info = 6
      ELSE IF (incx.EQ.0) THEN
          info = 8
      END IF
      IF (info.NE.0) THEN
          CALL xerbla('DTRSV ',info)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF (n.EQ.0) RETURN
*
      nounit = lsame(diag,'N')
*
*     Set up the start point in X if the increment is not unity. This
*     will be  ( N - 1 )*INCX  too small for descending loops.
*
      IF (incx.LE.0) THEN
          kx = 1 - (n-1)*incx
      ELSE IF (incx.NE.1) THEN
          kx = 1
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF (lsame(trans,'N')) THEN
*
*        Form  x := inv( A )*x.
*
          IF (lsame(uplo,'U')) THEN
              IF (incx.EQ.1) THEN
                  DO 20 j = n,1,-1
                      IF (x(j).NE.zero) THEN
                          IF (nounit) x(j) = x(j)/a(j,j)
                          temp = x(j)
                          DO 10 i = j - 1,1,-1
                              x(i) = x(i) - temp*a(i,j)
   10                     CONTINUE
                      END IF
   20             CONTINUE
              ELSE
                  jx = kx + (n-1)*incx
                  DO 40 j = n,1,-1
                      IF (x(jx).NE.zero) THEN
                          IF (nounit) x(jx) = x(jx)/a(j,j)
                          temp = x(jx)
                          ix = jx
                          DO 30 i = j - 1,1,-1
                              ix = ix - incx
                              x(ix) = x(ix) - temp*a(i,j)
   30                     CONTINUE
                      END IF
                      jx = jx - incx
   40             CONTINUE
              END IF
          ELSE
              IF (incx.EQ.1) THEN
                  DO 60 j = 1,n
                      IF (x(j).NE.zero) THEN
                          IF (nounit) x(j) = x(j)/a(j,j)
                          temp = x(j)
                          DO 50 i = j + 1,n
                              x(i) = x(i) - temp*a(i,j)
   50                     CONTINUE
                      END IF
   60             CONTINUE
              ELSE
                  jx = kx
                  DO 80 j = 1,n
                      IF (x(jx).NE.zero) THEN
                          IF (nounit) x(jx) = x(jx)/a(j,j)
                          temp = x(jx)
                          ix = jx
                          DO 70 i = j + 1,n
                              ix = ix + incx
                              x(ix) = x(ix) - temp*a(i,j)
   70                     CONTINUE
                      END IF
                      jx = jx + incx
   80             CONTINUE
              END IF
          END IF
      ELSE
*
*        Form  x := inv( A**T )*x.
*
          IF (lsame(uplo,'U')) THEN
              IF (incx.EQ.1) THEN
                  DO 100 j = 1,n
                      temp = x(j)
                      DO 90 i = 1,j - 1
                          temp = temp - a(i,j)*x(i)
   90                 CONTINUE
                      IF (nounit) temp = temp/a(j,j)
                      x(j) = temp
  100             CONTINUE
              ELSE
                  jx = kx
                  DO 120 j = 1,n
                      temp = x(jx)
                      ix = kx
                      DO 110 i = 1,j - 1
                          temp = temp - a(i,j)*x(ix)
                          ix = ix + incx
  110                 CONTINUE
                      IF (nounit) temp = temp/a(j,j)
                      x(jx) = temp
                      jx = jx + incx
  120             CONTINUE
              END IF
          ELSE
              IF (incx.EQ.1) THEN
                  DO 140 j = n,1,-1
                      temp = x(j)
                      DO 130 i = n,j + 1,-1
                          temp = temp - a(i,j)*x(i)
  130                 CONTINUE
                      IF (nounit) temp = temp/a(j,j)
                      x(j) = temp
  140             CONTINUE
              ELSE
                  kx = kx + (n-1)*incx
                  jx = kx
                  DO 160 j = n,1,-1
                      temp = x(jx)
                      ix = kx
                      DO 150 i = n,j + 1,-1
                          temp = temp - a(i,j)*x(ix)
                          ix = ix - incx
  150                 CONTINUE
                      IF (nounit) temp = temp/a(j,j)
                      x(jx) = temp
                      jx = jx - incx
  160             CONTINUE
              END IF
          END IF
      END IF
*
      RETURN
*
*     End of DTRSV
*
      END

! DGGEV
      SUBROUTINE dggev( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR, ALPHAI,
     $                  BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          JOBVL, JOBVR
      INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), ALPHAI( * ), ALPHAR( * ),
     $                   b( ldb, * ), beta( * ), vl( ldvl, * ),
     $                   vr( ldvr, * ), work( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      parameter( zero = 0.0d+0, one = 1.0d+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ILASCL, ILBSCL, ILV, ILVL, ILVR, LQUERY
      CHARACTER          CHTEMP
      INTEGER            ICOLS, IERR, IHI, IJOBVL, IJOBVR, ILEFT, ILO,
     $                   in, iright, irows, itau, iwrk, jc, jr, maxwrk,
     $                   minwrk
      DOUBLE PRECISION   ANRM, ANRMTO, BIGNUM, BNRM, BNRMTO, EPS,
     $                   smlnum, temp
*     ..
*     .. Local Arrays ..
      LOGICAL            LDUMMA( 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           dgeqrf, dggbak, dggbal, dgghrd, dhgeqz, dlabad,
     $                   dlacpy,dlascl, dlaset, dorgqr, dormqr, dtgevc,
     $                   xerbla
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH, DLANGE
      EXTERNAL           lsame, ilaenv, dlamch, dlange
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, max, sqrt
*     ..
*     .. Executable Statements ..
*
*     Decode the input arguments
*
      IF( lsame( jobvl, 'N' ) ) THEN
         ijobvl = 1
         ilvl = .false.
      ELSE IF( lsame( jobvl, 'V' ) ) THEN
         ijobvl = 2
         ilvl = .true.
      ELSE
         ijobvl = -1
         ilvl = .false.
      END IF
*
      IF( lsame( jobvr, 'N' ) ) THEN
         ijobvr = 1
         ilvr = .false.
      ELSE IF( lsame( jobvr, 'V' ) ) THEN
         ijobvr = 2
         ilvr = .true.
      ELSE
         ijobvr = -1
         ilvr = .false.
      END IF
      ilv = ilvl .OR. ilvr
*
*     Test the input arguments
*
      info = 0
      lquery = ( lwork.EQ.-1 )
      IF( ijobvl.LE.0 ) THEN
         info = -1
      ELSE IF( ijobvr.LE.0 ) THEN
         info = -2
      ELSE IF( n.LT.0 ) THEN
         info = -3
      ELSE IF( lda.LT.max( 1, n ) ) THEN
         info = -5
      ELSE IF( ldb.LT.max( 1, n ) ) THEN
         info = -7
      ELSE IF( ldvl.LT.1 .OR. ( ilvl .AND. ldvl.LT.n ) ) THEN
         info = -12
      ELSE IF( ldvr.LT.1 .OR. ( ilvr .AND. ldvr.LT.n ) ) THEN
         info = -14
      END IF
*
*     Compute workspace
*      (Note: Comments in the code beginning "Workspace:" describe the
*       minimal amount of workspace needed at that point in the code,
*       as well as the preferred amount for good performance.
*       NB refers to the optimal block size for the immediately
*       following subroutine, as returned by ILAENV. The workspace is
*       computed assuming ILO = 1 and IHI = N, the worst case.)
*
      IF( info.EQ.0 ) THEN
         minwrk = max( 1, 8*n )
         maxwrk = max( 1, n*( 7 +
     $                 ilaenv( 1, 'DGEQRF', ' ', n, 1, n, 0 ) ) )
         maxwrk = max( maxwrk, n*( 7 +
     $                 ilaenv( 1, 'DORMQR', ' ', n, 1, n, 0 ) ) )
         IF( ilvl ) THEN
            maxwrk = max( maxwrk, n*( 7 +
     $                 ilaenv( 1, 'DORGQR', ' ', n, 1, n, -1 ) ) )
         END IF
         work( 1 ) = maxwrk
*
         IF( lwork.LT.minwrk .AND. .NOT.lquery )
     $      info = -16
      END IF
*
      IF( info.NE.0 ) THEN
         CALL xerbla( 'DGGEV ', -info )
         RETURN
      ELSE IF( lquery ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( n.EQ.0 )
     $   RETURN
*
*     Get machine constants
*
      eps = dlamch( 'P' )
      smlnum = dlamch( 'S' )
      bignum = one / smlnum
      CALL dlabad( smlnum, bignum )
      smlnum = sqrt( smlnum ) / eps
      bignum = one / smlnum
*
*     Scale A if max element outside range [SMLNUM,BIGNUM]
*
      anrm = dlange( 'M', n, n, a, lda, work )
      ilascl = .false.
      IF( anrm.GT.zero .AND. anrm.LT.smlnum ) THEN
         anrmto = smlnum
         ilascl = .true.
      ELSE IF( anrm.GT.bignum ) THEN
         anrmto = bignum
         ilascl = .true.
      END IF
      IF( ilascl )
     $   CALL dlascl( 'G', 0, 0, anrm, anrmto, n, n, a, lda, ierr )
*
*     Scale B if max element outside range [SMLNUM,BIGNUM]
*
      bnrm = dlange( 'M', n, n, b, ldb, work )
      ilbscl = .false.
      IF( bnrm.GT.zero .AND. bnrm.LT.smlnum ) THEN
         bnrmto = smlnum
         ilbscl = .true.
      ELSE IF( bnrm.GT.bignum ) THEN
         bnrmto = bignum
         ilbscl = .true.
      END IF
      IF( ilbscl )
     $   CALL dlascl( 'G', 0, 0, bnrm, bnrmto, n, n, b, ldb, ierr )
*
*     Permute the matrices A, B to isolate eigenvalues if possible
*     (Workspace: need 6*N)
*
      ileft = 1
      iright = n + 1
      iwrk = iright + n
      CALL dggbal( 'P', n, a, lda, b, ldb, ilo, ihi, work( ileft ),
     $             work( iright ), work( iwrk ), ierr )
*
*     Reduce B to triangular form (QR decomposition of B)
*     (Workspace: need N, prefer N*NB)
*
      irows = ihi + 1 - ilo
      IF( ilv ) THEN
         icols = n + 1 - ilo
      ELSE
         icols = irows
      END IF
      itau = iwrk
      iwrk = itau + irows
      CALL dgeqrf( irows, icols, b( ilo, ilo ), ldb, work( itau ),
     $             work( iwrk ), lwork+1-iwrk, ierr )
*
*     Apply the orthogonal transformation to matrix A
*     (Workspace: need N, prefer N*NB)
*
      CALL dormqr( 'L', 'T', irows, icols, irows, b( ilo, ilo ), ldb,
     $             work( itau ), a( ilo, ilo ), lda, work( iwrk ),
     $             lwork+1-iwrk, ierr )
*
*     Initialize VL
*     (Workspace: need N, prefer N*NB)
*
      IF( ilvl ) THEN
         CALL dlaset( 'Full', n, n, zero, one, vl, ldvl )
         IF( irows.GT.1 ) THEN
            CALL dlacpy( 'L', irows-1, irows-1, b( ilo+1, ilo ), ldb,
     $                   vl( ilo+1, ilo ), ldvl )
         END IF
         CALL dorgqr( irows, irows, irows, vl( ilo, ilo ), ldvl,
     $                work( itau ), work( iwrk ), lwork+1-iwrk, ierr )
      END IF
*
*     Initialize VR
*
      IF( ilvr )
     $   CALL dlaset( 'Full', n, n, zero, one, vr, ldvr )
*
*     Reduce to generalized Hessenberg form
*     (Workspace: none needed)
*
      IF( ilv ) THEN
*
*        Eigenvectors requested -- work on whole matrix.
*
         CALL dgghrd( jobvl, jobvr, n, ilo, ihi, a, lda, b, ldb, vl,
     $                ldvl, vr, ldvr, ierr )
      ELSE
         CALL dgghrd( 'N', 'N', irows, 1, irows, a( ilo, ilo ), lda,
     $                b( ilo, ilo ), ldb, vl, ldvl, vr, ldvr, ierr )
      END IF
*
*     Perform QZ algorithm (Compute eigenvalues, and optionally, the
*     Schur forms and Schur vectors)
*     (Workspace: need N)
*
      iwrk = itau
      IF( ilv ) THEN
         chtemp = 'S'
      ELSE
         chtemp = 'E'
      END IF
      CALL dhgeqz( chtemp, jobvl, jobvr, n, ilo, ihi, a, lda, b, ldb,
     $             alphar, alphai, beta, vl, ldvl, vr, ldvr,
     $             work( iwrk ), lwork+1-iwrk, ierr )
      IF( ierr.NE.0 ) THEN
         IF( ierr.GT.0 .AND. ierr.LE.n ) THEN
            info = ierr
         ELSE IF( ierr.GT.n .AND. ierr.LE.2*n ) THEN
            info = ierr - n
         ELSE
            info = n + 1
         END IF
         GO TO 110
      END IF
*
*     Compute Eigenvectors
*     (Workspace: need 6*N)
*
      IF( ilv ) THEN
         IF( ilvl ) THEN
            IF( ilvr ) THEN
               chtemp = 'B'
            ELSE
               chtemp = 'L'
            END IF
         ELSE
            chtemp = 'R'
         END IF
         CALL dtgevc( chtemp, 'B', ldumma, n, a, lda, b, ldb, vl, ldvl,
     $                vr, ldvr, n, in, work( iwrk ), ierr )
         IF( ierr.NE.0 ) THEN
            info = n + 2
            GO TO 110
         END IF
*
*        Undo balancing on VL and VR and normalization
*        (Workspace: none needed)
*
         IF( ilvl ) THEN
            CALL dggbak( 'P', 'L', n, ilo, ihi, work( ileft ),
     $                   work( iright ), n, vl, ldvl, ierr )
            DO 50 jc = 1, n
               IF( alphai( jc ).LT.zero )
     $            GO TO 50
               temp = zero
               IF( alphai( jc ).EQ.zero ) THEN
                  DO 10 jr = 1, n
                     temp = max( temp, abs( vl( jr, jc ) ) )
   10             CONTINUE
               ELSE
                  DO 20 jr = 1, n
                     temp = max( temp, abs( vl( jr, jc ) )+
     $                      abs( vl( jr, jc+1 ) ) )
   20             CONTINUE
               END IF
               IF( temp.LT.smlnum )
     $            GO TO 50
               temp = one / temp
               IF( alphai( jc ).EQ.zero ) THEN
                  DO 30 jr = 1, n
                     vl( jr, jc ) = vl( jr, jc )*temp
   30             CONTINUE
               ELSE
                  DO 40 jr = 1, n
                     vl( jr, jc ) = vl( jr, jc )*temp
                     vl( jr, jc+1 ) = vl( jr, jc+1 )*temp
   40             CONTINUE
               END IF
   50       CONTINUE
         END IF
         IF( ilvr ) THEN
            CALL dggbak( 'P', 'R', n, ilo, ihi, work( ileft ),
     $                   work( iright ), n, vr, ldvr, ierr )
            DO 100 jc = 1, n
               IF( alphai( jc ).LT.zero )
     $            GO TO 100
               temp = zero
               IF( alphai( jc ).EQ.zero ) THEN
                  DO 60 jr = 1, n
                     temp = max( temp, abs( vr( jr, jc ) ) )
   60             CONTINUE
               ELSE
                  DO 70 jr = 1, n
                     temp = max( temp, abs( vr( jr, jc ) )+
     $                      abs( vr( jr, jc+1 ) ) )
   70             CONTINUE
               END IF
               IF( temp.LT.smlnum )
     $            GO TO 100
               temp = one / temp
               IF( alphai( jc ).EQ.zero ) THEN
                  DO 80 jr = 1, n
                     vr( jr, jc ) = vr( jr, jc )*temp
   80             CONTINUE
               ELSE
                  DO 90 jr = 1, n
                     vr( jr, jc ) = vr( jr, jc )*temp
                     vr( jr, jc+1 ) = vr( jr, jc+1 )*temp
   90             CONTINUE
               END IF
  100       CONTINUE
         END IF
*
*        End of eigenvector calculation
*
      END IF
*
*     Undo scaling if necessary
*
  110 CONTINUE
*
      IF( ilascl ) THEN
         CALL dlascl( 'G', 0, 0, anrmto, anrm, n, 1, alphar, n, ierr )
         CALL dlascl( 'G', 0, 0, anrmto, anrm, n, 1, alphai, n, ierr )
      END IF
*
      IF( ilbscl ) THEN
         CALL dlascl( 'G', 0, 0, bnrmto, bnrm, n, 1, beta, n, ierr )
      END IF
*
      work( 1 ) = maxwrk
      RETURN
*
*     End of DGGEV
*
      END

! DLABAD
      SUBROUTINE dlabad( SMALL, LARGE )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   LARGE, SMALL
*     ..
*
*  =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC          log10, sqrt
*     ..
*     .. Executable Statements ..
*
*     If it looks like we're on a Cray, take the square root of
*     SMALL and LARGE to avoid overflow and underflow problems.
*
      IF( log10( large ).GT.2000.d0 ) THEN
         small = sqrt( small )
         large = sqrt( large )
      END IF
*
      RETURN
*
*     End of DLABAD
*
      END

! DLANGE
      DOUBLE PRECISION FUNCTION dlange( NORM, M, N, A, LDA, WORK )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          norm
      INTEGER            lda, m, n
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   a( lda, * ), work( * )
*     ..
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   one, zero
      parameter( one = 1.0d+0, zero = 0.0d+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            i, j
      DOUBLE PRECISION   scale, sum, VALUE, temp
*     ..
*     .. External Subroutines ..
      EXTERNAL           dlassq
*     ..
*     .. External Functions ..
      LOGICAL            lsame, disnan
      EXTERNAL           lsame, disnan
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, min, sqrt
*     ..
*     .. Executable Statements ..
*
      IF( min( m, n ).EQ.0 ) THEN
         VALUE = zero
      ELSE IF( lsame( norm, 'M' ) ) THEN
*
*        Find max(abs(A(i,j))).
*
         VALUE = zero
         DO 20 j = 1, n
            DO 10 i = 1, m
               temp = abs( a( i, j ) )
               IF( VALUE.LT.temp .OR. disnan( temp ) ) VALUE = temp
   10       CONTINUE
   20    CONTINUE
      ELSE IF( ( lsame( norm, 'O' ) ) .OR. ( norm.EQ.'1' ) ) THEN
*
*        Find norm1(A).
*
         VALUE = zero
         DO 40 j = 1, n
            sum = zero
            DO 30 i = 1, m
               sum = sum + abs( a( i, j ) )
   30       CONTINUE
            IF( VALUE.LT.sum .OR. disnan( sum ) ) VALUE = sum
   40    CONTINUE
      ELSE IF( lsame( norm, 'I' ) ) THEN
*
*        Find normI(A).
*
         DO 50 i = 1, m
            work( i ) = zero
   50    CONTINUE
         DO 70 j = 1, n
            DO 60 i = 1, m
               work( i ) = work( i ) + abs( a( i, j ) )
   60       CONTINUE
   70    CONTINUE
         VALUE = zero
         DO 80 i = 1, m
            temp = work( i )
            IF( VALUE.LT.temp .OR. disnan( temp ) ) VALUE = temp
   80    CONTINUE
      ELSE IF( ( lsame( norm, 'F' ) ) .OR. ( lsame( norm, 'E' ) ) ) THEN
*
*        Find normF(A).
*
         scale = zero
         sum = one
         DO 90 j = 1, n
            CALL dlassq( m, a( 1, j ), 1, scale, sum )
   90    CONTINUE
         VALUE = scale*sqrt( sum )
      END IF
*
      dlange = VALUE
      RETURN
*
*     End of DLANGE
*
      END

! DGGBAL
      SUBROUTINE dggbal( JOB, N, A, LDA, B, LDB, ILO, IHI, LSCALE,
     $                   RSCALE, WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          JOB
      INTEGER            IHI, ILO, INFO, LDA, LDB, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), LSCALE( * ),
     $                   rscale( * ), work( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, HALF, ONE
      parameter( zero = 0.0d+0, half = 0.5d+0, one = 1.0d+0 )
      DOUBLE PRECISION   THREE, SCLFAC
      parameter( three = 3.0d+0, sclfac = 1.0d+1 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, ICAB, IFLOW, IP1, IR, IRAB, IT, J, JC, JP1,
     $                   k, kount, l, lcab, lm1, lrab, lsfmax, lsfmin,
     $                   m, nr, nrp2
      DOUBLE PRECISION   ALPHA, BASL, BETA, CAB, CMAX, COEF, COEF2,
     $                   coef5, cor, ew, ewc, gamma, pgamma, rab, sfmax,
     $                   sfmin, sum, t, ta, tb, tc
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IDAMAX
      DOUBLE PRECISION   DDOT, DLAMCH
      EXTERNAL           lsame, idamax, ddot, dlamch
*     ..
*     .. External Subroutines ..
      EXTERNAL           daxpy, dscal, dswap, xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, dble, int, log10, max, min, sign
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      info = 0
      IF( .NOT.lsame( job, 'N' ) .AND. .NOT.lsame( job, 'P' ) .AND.
     $    .NOT.lsame( job, 'S' ) .AND. .NOT.lsame( job, 'B' ) ) THEN
         info = -1
      ELSE IF( n.LT.0 ) THEN
         info = -2
      ELSE IF( lda.LT.max( 1, n ) ) THEN
         info = -4
      ELSE IF( ldb.LT.max( 1, n ) ) THEN
         info = -6
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'DGGBAL', -info )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( n.EQ.0 ) THEN
         ilo = 1
         ihi = n
         RETURN
      END IF
*
      IF( n.EQ.1 ) THEN
         ilo = 1
         ihi = n
         lscale( 1 ) = one
         rscale( 1 ) = one
         RETURN
      END IF
*
      IF( lsame( job, 'N' ) ) THEN
         ilo = 1
         ihi = n
         DO 10 i = 1, n
            lscale( i ) = one
            rscale( i ) = one
   10    CONTINUE
         RETURN
      END IF
*
      k = 1
      l = n
      IF( lsame( job, 'S' ) )
     $   GO TO 190
*
      GO TO 30
*
*     Permute the matrices A and B to isolate the eigenvalues.
*
*     Find row with one nonzero in columns 1 through L
*
   20 CONTINUE
      l = lm1
      IF( l.NE.1 )
     $   GO TO 30
*
      rscale( 1 ) = one
      lscale( 1 ) = one
      GO TO 190
*
   30 CONTINUE
      lm1 = l - 1
      DO 80 i = l, 1, -1
         DO 40 j = 1, lm1
            jp1 = j + 1
            IF( a( i, j ).NE.zero .OR. b( i, j ).NE.zero )
     $         GO TO 50
   40    CONTINUE
         j = l
         GO TO 70
*
   50    CONTINUE
         DO 60 j = jp1, l
            IF( a( i, j ).NE.zero .OR. b( i, j ).NE.zero )
     $         GO TO 80
   60    CONTINUE
         j = jp1 - 1
*
   70    CONTINUE
         m = l
         iflow = 1
         GO TO 160
   80 CONTINUE
      GO TO 100
*
*     Find column with one nonzero in rows K through N
*
   90 CONTINUE
      k = k + 1
*
  100 CONTINUE
      DO 150 j = k, l
         DO 110 i = k, lm1
            ip1 = i + 1
            IF( a( i, j ).NE.zero .OR. b( i, j ).NE.zero )
     $         GO TO 120
  110    CONTINUE
         i = l
         GO TO 140
  120    CONTINUE
         DO 130 i = ip1, l
            IF( a( i, j ).NE.zero .OR. b( i, j ).NE.zero )
     $         GO TO 150
  130    CONTINUE
         i = ip1 - 1
  140    CONTINUE
         m = k
         iflow = 2
         GO TO 160
  150 CONTINUE
      GO TO 190
*
*     Permute rows M and I
*
  160 CONTINUE
      lscale( m ) = i
      IF( i.EQ.m )
     $   GO TO 170
      CALL dswap( n-k+1, a( i, k ), lda, a( m, k ), lda )
      CALL dswap( n-k+1, b( i, k ), ldb, b( m, k ), ldb )
*
*     Permute columns M and J
*
  170 CONTINUE
      rscale( m ) = j
      IF( j.EQ.m )
     $   GO TO 180
      CALL dswap( l, a( 1, j ), 1, a( 1, m ), 1 )
      CALL dswap( l, b( 1, j ), 1, b( 1, m ), 1 )
*
  180 CONTINUE
      GO TO ( 20, 90 )iflow
*
  190 CONTINUE
      ilo = k
      ihi = l
*
      IF( lsame( job, 'P' ) ) THEN
         DO 195 i = ilo, ihi
            lscale( i ) = one
            rscale( i ) = one
  195    CONTINUE
         RETURN
      END IF
*
      IF( ilo.EQ.ihi )
     $   RETURN
*
*     Balance the submatrix in rows ILO to IHI.
*
      nr = ihi - ilo + 1
      DO 200 i = ilo, ihi
         rscale( i ) = zero
         lscale( i ) = zero
*
         work( i ) = zero
         work( i+n ) = zero
         work( i+2*n ) = zero
         work( i+3*n ) = zero
         work( i+4*n ) = zero
         work( i+5*n ) = zero
  200 CONTINUE
*
*     Compute right side vector in resulting linear equations
*
      basl = log10( sclfac )
      DO 240 i = ilo, ihi
         DO 230 j = ilo, ihi
            tb = b( i, j )
            ta = a( i, j )
            IF( ta.EQ.zero )
     $         GO TO 210
            ta = log10( abs( ta ) ) / basl
  210       CONTINUE
            IF( tb.EQ.zero )
     $         GO TO 220
            tb = log10( abs( tb ) ) / basl
  220       CONTINUE
            work( i+4*n ) = work( i+4*n ) - ta - tb
            work( j+5*n ) = work( j+5*n ) - ta - tb
  230    CONTINUE
  240 CONTINUE
*
      coef = one / dble( 2*nr )
      coef2 = coef*coef
      coef5 = half*coef2
      nrp2 = nr + 2
      beta = zero
      it = 1
*
*     Start generalized conjugate gradient iteration
*
  250 CONTINUE
*
      gamma = ddot( nr, work( ilo+4*n ), 1, work( ilo+4*n ), 1 ) +
     $        ddot( nr, work( ilo+5*n ), 1, work( ilo+5*n ), 1 )
*
      ew = zero
      ewc = zero
      DO 260 i = ilo, ihi
         ew = ew + work( i+4*n )
         ewc = ewc + work( i+5*n )
  260 CONTINUE
*
      gamma = coef*gamma - coef2*( ew**2+ewc**2 ) - coef5*( ew-ewc )**2
      IF( gamma.EQ.zero )
     $   GO TO 350
      IF( it.NE.1 )
     $   beta = gamma / pgamma
      t = coef5*( ewc-three*ew )
      tc = coef5*( ew-three*ewc )
*
      CALL dscal( nr, beta, work( ilo ), 1 )
      CALL dscal( nr, beta, work( ilo+n ), 1 )
*
      CALL daxpy( nr, coef, work( ilo+4*n ), 1, work( ilo+n ), 1 )
      CALL daxpy( nr, coef, work( ilo+5*n ), 1, work( ilo ), 1 )
*
      DO 270 i = ilo, ihi
         work( i ) = work( i ) + tc
         work( i+n ) = work( i+n ) + t
  270 CONTINUE
*
*     Apply matrix to vector
*
      DO 300 i = ilo, ihi
         kount = 0
         sum = zero
         DO 290 j = ilo, ihi
            IF( a( i, j ).EQ.zero )
     $         GO TO 280
            kount = kount + 1
            sum = sum + work( j )
  280       CONTINUE
            IF( b( i, j ).EQ.zero )
     $         GO TO 290
            kount = kount + 1
            sum = sum + work( j )
  290    CONTINUE
         work( i+2*n ) = dble( kount )*work( i+n ) + sum
  300 CONTINUE
*
      DO 330 j = ilo, ihi
         kount = 0
         sum = zero
         DO 320 i = ilo, ihi
            IF( a( i, j ).EQ.zero )
     $         GO TO 310
            kount = kount + 1
            sum = sum + work( i+n )
  310       CONTINUE
            IF( b( i, j ).EQ.zero )
     $         GO TO 320
            kount = kount + 1
            sum = sum + work( i+n )
  320    CONTINUE
         work( j+3*n ) = dble( kount )*work( j ) + sum
  330 CONTINUE
*
      sum = ddot( nr, work( ilo+n ), 1, work( ilo+2*n ), 1 ) +
     $      ddot( nr, work( ilo ), 1, work( ilo+3*n ), 1 )
      alpha = gamma / sum
*
*     Determine correction to current iteration
*
      cmax = zero
      DO 340 i = ilo, ihi
         cor = alpha*work( i+n )
         IF( abs( cor ).GT.cmax )
     $      cmax = abs( cor )
         lscale( i ) = lscale( i ) + cor
         cor = alpha*work( i )
         IF( abs( cor ).GT.cmax )
     $      cmax = abs( cor )
         rscale( i ) = rscale( i ) + cor
  340 CONTINUE
      IF( cmax.LT.half )
     $   GO TO 350
*
      CALL daxpy( nr, -alpha, work( ilo+2*n ), 1, work( ilo+4*n ), 1 )
      CALL daxpy( nr, -alpha, work( ilo+3*n ), 1, work( ilo+5*n ), 1 )
*
      pgamma = gamma
      it = it + 1
      IF( it.LE.nrp2 )
     $   GO TO 250
*
*     End generalized conjugate gradient iteration
*
  350 CONTINUE
      sfmin = dlamch( 'S' )
      sfmax = one / sfmin
      lsfmin = int( log10( sfmin ) / basl+one )
      lsfmax = int( log10( sfmax ) / basl )
      DO 360 i = ilo, ihi
         irab = idamax( n-ilo+1, a( i, ilo ), lda )
         rab = abs( a( i, irab+ilo-1 ) )
         irab = idamax( n-ilo+1, b( i, ilo ), ldb )
         rab = max( rab, abs( b( i, irab+ilo-1 ) ) )
         lrab = int( log10( rab+sfmin ) / basl+one )
         ir = int(lscale( i ) + sign( half, lscale( i ) ))
         ir = min( max( ir, lsfmin ), lsfmax, lsfmax-lrab )
         lscale( i ) = sclfac**ir
         icab = idamax( ihi, a( 1, i ), 1 )
         cab = abs( a( icab, i ) )
         icab = idamax( ihi, b( 1, i ), 1 )
         cab = max( cab, abs( b( icab, i ) ) )
         lcab = int( log10( cab+sfmin ) / basl+one )
         jc = int(rscale( i ) + sign( half, rscale( i ) ))
         jc = min( max( jc, lsfmin ), lsfmax, lsfmax-lcab )
         rscale( i ) = sclfac**jc
  360 CONTINUE
*
*     Row scaling of matrices A and B
*
      DO 370 i = ilo, ihi
         CALL dscal( n-ilo+1, lscale( i ), a( i, ilo ), lda )
         CALL dscal( n-ilo+1, lscale( i ), b( i, ilo ), ldb )
  370 CONTINUE
*
*     Column scaling of matrices A and B
*
      DO 380 j = ilo, ihi
         CALL dscal( ihi, rscale( j ), a( 1, j ), 1 )
         CALL dscal( ihi, rscale( j ), b( 1, j ), 1 )
  380 CONTINUE
*
      RETURN
*
*     End of DGGBAL
*
      END

! DGEQRF
      SUBROUTINE dgeqrf( M, N, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, K, LDWORK, LWKOPT, NB,
     $                   NBMIN, NX
*     ..
*     .. External Subroutines ..
      EXTERNAL           dgeqr2, dlarfb, dlarft, xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          max, min
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ilaenv
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      k = min( m, n )
      info = 0
      nb = ilaenv( 1, 'DGEQRF', ' ', m, n, -1, -1 )
      lquery = ( lwork.EQ.-1 )
      IF( m.LT.0 ) THEN
         info = -1
      ELSE IF( n.LT.0 ) THEN
         info = -2
      ELSE IF( lda.LT.max( 1, m ) ) THEN
         info = -4
      ELSE IF( .NOT.lquery ) THEN
         IF( lwork.LE.0 .OR. ( m.GT.0 .AND. lwork.LT.max( 1, n ) ) )
     $      info = -7
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'DGEQRF', -info )
         RETURN
      ELSE IF( lquery ) THEN
         IF( k.EQ.0 ) THEN
            lwkopt = 1
         ELSE
            lwkopt = n*nb
         END IF
         work( 1 ) = lwkopt
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( k.EQ.0 ) THEN
         work( 1 ) = 1
         RETURN
      END IF
*
      nbmin = 2
      nx = 0
      iws = n
      IF( nb.GT.1 .AND. nb.LT.k ) THEN
*
*        Determine when to cross over from blocked to unblocked code.
*
         nx = max( 0, ilaenv( 3, 'DGEQRF', ' ', m, n, -1, -1 ) )
         IF( nx.LT.k ) THEN
*
*           Determine if workspace is large enough for blocked code.
*
            ldwork = n
            iws = ldwork*nb
            IF( lwork.LT.iws ) THEN
*
*              Not enough workspace to use optimal NB:  reduce NB and
*              determine the minimum value of NB.
*
               nb = lwork / ldwork
               nbmin = max( 2, ilaenv( 2, 'DGEQRF', ' ', m, n, -1,
     $                 -1 ) )
            END IF
         END IF
      END IF
*
      IF( nb.GE.nbmin .AND. nb.LT.k .AND. nx.LT.k ) THEN
*
*        Use blocked code initially
*
         DO 10 i = 1, k - nx, nb
            ib = min( k-i+1, nb )
*
*           Compute the QR factorization of the current block
*           A(i:m,i:i+ib-1)
*
            CALL dgeqr2( m-i+1, ib, a( i, i ), lda, tau( i ), work,
     $                   iinfo )
            IF( i+ib.LE.n ) THEN
*
*              Form the triangular factor of the block reflector
*              H = H(i) H(i+1) . . . H(i+ib-1)
*
               CALL dlarft( 'Forward', 'Columnwise', m-i+1, ib,
     $                      a( i, i ), lda, tau( i ), work, ldwork )
*
*              Apply H**T to A(i:m,i+ib:n) from the left
*
               CALL dlarfb( 'Left', 'Transpose', 'Forward',
     $                      'Columnwise', m-i+1, n-i-ib+1, ib,
     $                      a( i, i ), lda, work, ldwork, a( i, i+ib ),
     $                      lda, work( ib+1 ), ldwork )
            END IF
   10    CONTINUE
      ELSE
         i = 1
      END IF
*
*     Use unblocked code to factor the last or only block.
*
      IF( i.LE.k )
     $   CALL dgeqr2( m-i+1, n-i+1, a( i, i ), lda, tau( i ), work,
     $                iinfo )
*
      work( 1 ) = iws
      RETURN
*
*     End of DGEQRF
*
      END

! DORMQR
      SUBROUTINE dormqr( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
     $                   WORK, LWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, LWORK, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NBMAX, LDT, TSIZE
      parameter( nbmax = 64, ldt = nbmax+1,
     $                     tsize = ldt*nbmax )
*     ..
*     .. Local Scalars ..
      LOGICAL            LEFT, LQUERY, NOTRAN
      INTEGER            I, I1, I2, I3, IB, IC, IINFO, IWT, JC, LDWORK,
     $                   lwkopt, mi, nb, nbmin, ni, nq, nw
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           lsame, ilaenv
*     ..
*     .. External Subroutines ..
      EXTERNAL           dlarfb, dlarft, dorm2r, xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          max, min
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      info = 0
      left = lsame( side, 'L' )
      notran = lsame( trans, 'N' )
      lquery = ( lwork.EQ.-1 )
*
*     NQ is the order of Q and NW is the minimum dimension of WORK
*
      IF( left ) THEN
         nq = m
         nw = max( 1, n )
      ELSE
         nq = n
         nw = max( 1, m )
      END IF
      IF( .NOT.left .AND. .NOT.lsame( side, 'R' ) ) THEN
         info = -1
      ELSE IF( .NOT.notran .AND. .NOT.lsame( trans, 'T' ) ) THEN
         info = -2
      ELSE IF( m.LT.0 ) THEN
         info = -3
      ELSE IF( n.LT.0 ) THEN
         info = -4
      ELSE IF( k.LT.0 .OR. k.GT.nq ) THEN
         info = -5
      ELSE IF( lda.LT.max( 1, nq ) ) THEN
         info = -7
      ELSE IF( ldc.LT.max( 1, m ) ) THEN
         info = -10
      ELSE IF( lwork.LT.nw .AND. .NOT.lquery ) THEN
         info = -12
      END IF
*
      IF( info.EQ.0 ) THEN
*
*        Compute the workspace requirements
*
         nb = min( nbmax, ilaenv( 1, 'DORMQR', side // trans, m, n, k,
     $        -1 ) )
         lwkopt = nw*nb + tsize
         work( 1 ) = lwkopt
      END IF
*
      IF( info.NE.0 ) THEN
         CALL xerbla( 'DORMQR', -info )
         RETURN
      ELSE IF( lquery ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( m.EQ.0 .OR. n.EQ.0 .OR. k.EQ.0 ) THEN
         work( 1 ) = 1
         RETURN
      END IF
*
      nbmin = 2
      ldwork = nw
      IF( nb.GT.1 .AND. nb.LT.k ) THEN
         IF( lwork.LT.lwkopt ) THEN
            nb = (lwork-tsize) / ldwork
            nbmin = max( 2, ilaenv( 2, 'DORMQR', side // trans, m, n, k,
     $              -1 ) )
         END IF
      END IF
*
      IF( nb.LT.nbmin .OR. nb.GE.k ) THEN
*
*        Use unblocked code
*
         CALL dorm2r( side, trans, m, n, k, a, lda, tau, c, ldc, work,
     $                iinfo )
      ELSE
*
*        Use blocked code
*
         iwt = 1 + nw*nb
         IF( ( left .AND. .NOT.notran ) .OR.
     $       ( .NOT.left .AND. notran ) ) THEN
            i1 = 1
            i2 = k
            i3 = nb
         ELSE
            i1 = ( ( k-1 ) / nb )*nb + 1
            i2 = 1
            i3 = -nb
         END IF
*
         IF( left ) THEN
            ni = n
            jc = 1
         ELSE
            mi = m
            ic = 1
         END IF
*
         DO 10 i = i1, i2, i3
            ib = min( nb, k-i+1 )
*
*           Form the triangular factor of the block reflector
*           H = H(i) H(i+1) . . . H(i+ib-1)
*
            CALL dlarft( 'Forward', 'Columnwise', nq-i+1, ib, a( i, i ),
     $                   lda, tau( i ), work( iwt ), ldt )
            IF( left ) THEN
*
*              H or H**T is applied to C(i:m,1:n)
*
               mi = m - i + 1
               ic = i
            ELSE
*
*              H or H**T is applied to C(1:m,i:n)
*
               ni = n - i + 1
               jc = i
            END IF
*
*           Apply H or H**T
*
            CALL dlarfb( side, trans, 'Forward', 'Columnwise', mi, ni,
     $                   ib, a( i, i ), lda, work( iwt ), ldt,
     $                   c( ic, jc ), ldc, work, ldwork )
   10    CONTINUE
      END IF
      work( 1 ) = lwkopt
      RETURN
*
*     End of DORMQR
*
      END

! DGGHRD
      SUBROUTINE dgghrd( COMPQ, COMPZ, N, ILO, IHI, A, LDA, B, LDB, Q,
     $                   LDQ, Z, LDZ, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          COMPQ, COMPZ
      INTEGER            IHI, ILO, INFO, LDA, LDB, LDQ, LDZ, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), Q( LDQ, * ),
     $                   z( ldz, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      parameter( one = 1.0d+0, zero = 0.0d+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ILQ, ILZ
      INTEGER            ICOMPQ, ICOMPZ, JCOL, JROW
      DOUBLE PRECISION   C, S, TEMP
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           lsame
*     ..
*     .. External Subroutines ..
      EXTERNAL           dlartg, dlaset, drot, xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          max
*     ..
*     .. Executable Statements ..
*
*     Decode COMPQ
*
      IF( lsame( compq, 'N' ) ) THEN
         ilq = .false.
         icompq = 1
      ELSE IF( lsame( compq, 'V' ) ) THEN
         ilq = .true.
         icompq = 2
      ELSE IF( lsame( compq, 'I' ) ) THEN
         ilq = .true.
         icompq = 3
      ELSE
         icompq = 0
      END IF
*
*     Decode COMPZ
*
      IF( lsame( compz, 'N' ) ) THEN
         ilz = .false.
         icompz = 1
      ELSE IF( lsame( compz, 'V' ) ) THEN
         ilz = .true.
         icompz = 2
      ELSE IF( lsame( compz, 'I' ) ) THEN
         ilz = .true.
         icompz = 3
      ELSE
         icompz = 0
      END IF
*
*     Test the input parameters.
*
      info = 0
      IF( icompq.LE.0 ) THEN
         info = -1
      ELSE IF( icompz.LE.0 ) THEN
         info = -2
      ELSE IF( n.LT.0 ) THEN
         info = -3
      ELSE IF( ilo.LT.1 ) THEN
         info = -4
      ELSE IF( ihi.GT.n .OR. ihi.LT.ilo-1 ) THEN
         info = -5
      ELSE IF( lda.LT.max( 1, n ) ) THEN
         info = -7
      ELSE IF( ldb.LT.max( 1, n ) ) THEN
         info = -9
      ELSE IF( ( ilq .AND. ldq.LT.n ) .OR. ldq.LT.1 ) THEN
         info = -11
      ELSE IF( ( ilz .AND. ldz.LT.n ) .OR. ldz.LT.1 ) THEN
         info = -13
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'DGGHRD', -info )
         RETURN
      END IF
*
*     Initialize Q and Z if desired.
*
      IF( icompq.EQ.3 )
     $   CALL dlaset( 'Full', n, n, zero, one, q, ldq )
      IF( icompz.EQ.3 )
     $   CALL dlaset( 'Full', n, n, zero, one, z, ldz )
*
*     Quick return if possible
*
      IF( n.LE.1 )
     $   RETURN
*
*     Zero out lower triangle of B
*
      DO 20 jcol = 1, n - 1
         DO 10 jrow = jcol + 1, n
            b( jrow, jcol ) = zero
   10    CONTINUE
   20 CONTINUE
*
*     Reduce A and B
*
      DO 40 jcol = ilo, ihi - 2
*
         DO 30 jrow = ihi, jcol + 2, -1
*
*           Step 1: rotate rows JROW-1, JROW to kill A(JROW,JCOL)
*
            temp = a( jrow-1, jcol )
            CALL dlartg( temp, a( jrow, jcol ), c, s,
     $                   a( jrow-1, jcol ) )
            a( jrow, jcol ) = zero
            CALL drot( n-jcol, a( jrow-1, jcol+1 ), lda,
     $                 a( jrow, jcol+1 ), lda, c, s )
            CALL drot( n+2-jrow, b( jrow-1, jrow-1 ), ldb,
     $                 b( jrow, jrow-1 ), ldb, c, s )
            IF( ilq )
     $         CALL drot( n, q( 1, jrow-1 ), 1, q( 1, jrow ), 1, c, s )
*
*           Step 2: rotate columns JROW, JROW-1 to kill B(JROW,JROW-1)
*
            temp = b( jrow, jrow )
            CALL dlartg( temp, b( jrow, jrow-1 ), c, s,
     $                   b( jrow, jrow ) )
            b( jrow, jrow-1 ) = zero
            CALL drot( ihi, a( 1, jrow ), 1, a( 1, jrow-1 ), 1, c, s )
            CALL drot( jrow-1, b( 1, jrow ), 1, b( 1, jrow-1 ), 1, c,
     $                 s )
            IF( ilz )
     $         CALL drot( n, z( 1, jrow ), 1, z( 1, jrow-1 ), 1, c, s )
   30    CONTINUE
   40 CONTINUE
*
      RETURN
*
*     End of DGGHRD
*
      END

! DHGEQZ
      SUBROUTINE dhgeqz( JOB, COMPQ, COMPZ, N, ILO, IHI, H, LDH, T, LDT,
     $                   ALPHAR, ALPHAI, BETA, Q, LDQ, Z, LDZ, WORK,
     $                   LWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          COMPQ, COMPZ, JOB
      INTEGER            IHI, ILO, INFO, LDH, LDQ, LDT, LDZ, LWORK, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   ALPHAI( * ), ALPHAR( * ), BETA( * ),
     $                   H( LDH, * ), Q( LDQ, * ), T( LDT, * ),
     $                   work( * ), z( ldz, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
*    $                     SAFETY = 1.0E+0 )
      DOUBLE PRECISION   HALF, ZERO, ONE, SAFETY
      PARAMETER          ( HALF = 0.5d+0, zero = 0.0d+0, one = 1.0d+0,
     $                   safety = 1.0d+2 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ILAZR2, ILAZRO, ILPIVT, ILQ, ILSCHR, ILZ,
     $                   LQUERY
      INTEGER            ICOMPQ, ICOMPZ, IFIRST, IFRSTM, IITER, ILAST,
     $                   ILASTM, IN, ISCHUR, ISTART, J, JC, JCH, JITER,
     $                   jr, maxit
      DOUBLE PRECISION   A11, A12, A1I, A1R, A21, A22, A2I, A2R, AD11,
     $                   AD11L, AD12, AD12L, AD21, AD21L, AD22, AD22L,
     $                   ad32l, an, anorm, ascale, atol, b11, b1a, b1i,
     $                   b1r, b22, b2a, b2i, b2r, bn, bnorm, bscale,
     $                   btol, c, c11i, c11r, c12, c21, c22i, c22r, cl,
     $                   cq, cr, cz, eshift, s, s1, s1inv, s2, safmax,
     $                   safmin, scale, sl, sqi, sqr, sr, szi, szr, t1,
     $                   tau, temp, temp2, tempi, tempr, u1, u12, u12l,
     $                   u2, ulp, vs, w11, w12, w21, w22, wabs, wi, wr,
     $                   wr2
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   V( 3 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANHS, DLAPY2, DLAPY3
      EXTERNAL           lsame, dlamch, dlanhs, dlapy2, dlapy3
*     ..
*     .. External Subroutines ..
      EXTERNAL           dlag2, dlarfg, dlartg, dlaset, dlasv2, drot,
     $                   xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, dble, max, min, sqrt
*     ..
*     .. Executable Statements ..
*
*     Decode JOB, COMPQ, COMPZ
*
      IF( lsame( job, 'E' ) ) THEN
         ilschr = .false.
         ischur = 1
      ELSE IF( lsame( job, 'S' ) ) THEN
         ilschr = .true.
         ischur = 2
      ELSE
         ischur = 0
      END IF
*
      IF( lsame( compq, 'N' ) ) THEN
         ilq = .false.
         icompq = 1
      ELSE IF( lsame( compq, 'V' ) ) THEN
         ilq = .true.
         icompq = 2
      ELSE IF( lsame( compq, 'I' ) ) THEN
         ilq = .true.
         icompq = 3
      ELSE
         icompq = 0
      END IF
*
      IF( lsame( compz, 'N' ) ) THEN
         ilz = .false.
         icompz = 1
      ELSE IF( lsame( compz, 'V' ) ) THEN
         ilz = .true.
         icompz = 2
      ELSE IF( lsame( compz, 'I' ) ) THEN
         ilz = .true.
         icompz = 3
      ELSE
         icompz = 0
      END IF
*
*     Check Argument Values
*
      info = 0
      work( 1 ) = max( 1, n )
      lquery = ( lwork.EQ.-1 )
      IF( ischur.EQ.0 ) THEN
         info = -1
      ELSE IF( icompq.EQ.0 ) THEN
         info = -2
      ELSE IF( icompz.EQ.0 ) THEN
         info = -3
      ELSE IF( n.LT.0 ) THEN
         info = -4
      ELSE IF( ilo.LT.1 ) THEN
         info = -5
      ELSE IF( ihi.GT.n .OR. ihi.LT.ilo-1 ) THEN
         info = -6
      ELSE IF( ldh.LT.n ) THEN
         info = -8
      ELSE IF( ldt.LT.n ) THEN
         info = -10
      ELSE IF( ldq.LT.1 .OR. ( ilq .AND. ldq.LT.n ) ) THEN
         info = -15
      ELSE IF( ldz.LT.1 .OR. ( ilz .AND. ldz.LT.n ) ) THEN
         info = -17
      ELSE IF( lwork.LT.max( 1, n ) .AND. .NOT.lquery ) THEN
         info = -19
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'DHGEQZ', -info )
         RETURN
      ELSE IF( lquery ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( n.LE.0 ) THEN
         work( 1 ) = dble( 1 )
         RETURN
      END IF
*
*     Initialize Q and Z
*
      IF( icompq.EQ.3 )
     $   CALL dlaset( 'Full', n, n, zero, one, q, ldq )
      IF( icompz.EQ.3 )
     $   CALL dlaset( 'Full', n, n, zero, one, z, ldz )
*
*     Machine Constants
*
      in = ihi + 1 - ilo
      safmin = dlamch( 'S' )
      safmax = one / safmin
      ulp = dlamch( 'E' )*dlamch( 'B' )
      anorm = dlanhs( 'F', in, h( ilo, ilo ), ldh, work )
      bnorm = dlanhs( 'F', in, t( ilo, ilo ), ldt, work )
      atol = max( safmin, ulp*anorm )
      btol = max( safmin, ulp*bnorm )
      ascale = one / max( safmin, anorm )
      bscale = one / max( safmin, bnorm )
*
*     Set Eigenvalues IHI+1:N
*
      DO 30 j = ihi + 1, n
         IF( t( j, j ).LT.zero ) THEN
            IF( ilschr ) THEN
               DO 10 jr = 1, j
                  h( jr, j ) = -h( jr, j )
                  t( jr, j ) = -t( jr, j )
   10          CONTINUE
            ELSE
               h( j, j ) = -h( j, j )
               t( j, j ) = -t( j, j )
            END IF
            IF( ilz ) THEN
               DO 20 jr = 1, n
                  z( jr, j ) = -z( jr, j )
   20          CONTINUE
            END IF
         END IF
         alphar( j ) = h( j, j )
         alphai( j ) = zero
         beta( j ) = t( j, j )
   30 CONTINUE
*
*     If IHI < ILO, skip QZ steps
*
      IF( ihi.LT.ilo )
     $   GO TO 380
*
*     MAIN QZ ITERATION LOOP
*
*     Initialize dynamic indices
*
*     Eigenvalues ILAST+1:N have been found.
*        Column operations modify rows IFRSTM:whatever.
*        Row operations modify columns whatever:ILASTM.
*
*     If only eigenvalues are being computed, then
*        IFRSTM is the row of the last splitting row above row ILAST;
*        this is always at least ILO.
*     IITER counts iterations since the last eigenvalue was found,
*        to tell when to use an extraordinary shift.
*     MAXIT is the maximum number of QZ sweeps allowed.
*
      ilast = ihi
      IF( ilschr ) THEN
         ifrstm = 1
         ilastm = n
      ELSE
         ifrstm = ilo
         ilastm = ihi
      END IF
      iiter = 0
      eshift = zero
      maxit = 30*( ihi-ilo+1 )
*
      DO 360 jiter = 1, maxit
*
*        Split the matrix if possible.
*
*        Two tests:
*           1: H(j,j-1)=0  or  j=ILO
*           2: T(j,j)=0
*
         IF( ilast.EQ.ilo ) THEN
*
*           Special case: j=ILAST
*
            GO TO 80
         ELSE
            IF( abs( h( ilast, ilast-1 ) ).LE.max( safmin, ulp*( 
     $         abs( h( ilast, ilast ) ) + abs( h( ilast-1, ilast-1 ) ) 
     $         ) ) ) THEN
               h( ilast, ilast-1 ) = zero
               GO TO 80
            END IF
         END IF
*
         IF( abs( t( ilast, ilast ) ).LE.btol ) THEN
            t( ilast, ilast ) = zero
            GO TO 70
         END IF
*
*        General case: j<ILAST
*
         DO 60 j = ilast - 1, ilo, -1
*
*           Test 1: for H(j,j-1)=0 or j=ILO
*
            IF( j.EQ.ilo ) THEN
               ilazro = .true.
            ELSE
               IF( abs( h( j, j-1 ) ).LE.max( safmin, ulp*( 
     $         abs( h( j, j ) ) + abs( h( j-1, j-1 ) ) 
     $         ) ) ) THEN
                  h( j, j-1 ) = zero
                  ilazro = .true.
               ELSE
                  ilazro = .false.
               END IF
            END IF
*
*           Test 2: for T(j,j)=0
*
            IF( abs( t( j, j ) ).LT.btol ) THEN
               t( j, j ) = zero
*
*              Test 1a: Check for 2 consecutive small subdiagonals in A
*
               ilazr2 = .false.
               IF( .NOT.ilazro ) THEN
                  temp = abs( h( j, j-1 ) )
                  temp2 = abs( h( j, j ) )
                  tempr = max( temp, temp2 )
                  IF( tempr.LT.one .AND. tempr.NE.zero ) THEN
                     temp = temp / tempr
                     temp2 = temp2 / tempr
                  END IF
                  IF( temp*( ascale*abs( h( j+1, j ) ) ).LE.temp2*
     $                ( ascale*atol ) )ilazr2 = .true.
               END IF
*
*              If both tests pass (1 & 2), i.e., the leading diagonal
*              element of B in the block is zero, split a 1x1 block off
*              at the top. (I.e., at the J-th row/column) The leading
*              diagonal element of the remainder can also be zero, so
*              this may have to be done repeatedly.
*
               IF( ilazro .OR. ilazr2 ) THEN
                  DO 40 jch = j, ilast - 1
                     temp = h( jch, jch )
                     CALL dlartg( temp, h( jch+1, jch ), c, s,
     $                            h( jch, jch ) )
                     h( jch+1, jch ) = zero
                     CALL drot( ilastm-jch, h( jch, jch+1 ), ldh,
     $                          h( jch+1, jch+1 ), ldh, c, s )
                     CALL drot( ilastm-jch, t( jch, jch+1 ), ldt,
     $                          t( jch+1, jch+1 ), ldt, c, s )
                     IF( ilq )
     $                  CALL drot( n, q( 1, jch ), 1, q( 1, jch+1 ), 1,
     $                             c, s )
                     IF( ilazr2 )
     $                  h( jch, jch-1 ) = h( jch, jch-1 )*c
                     ilazr2 = .false.
                     IF( abs( t( jch+1, jch+1 ) ).GE.btol ) THEN
                        IF( jch+1.GE.ilast ) THEN
                           GO TO 80
                        ELSE
                           ifirst = jch + 1
                           GO TO 110
                        END IF
                     END IF
                     t( jch+1, jch+1 ) = zero
   40             CONTINUE
                  GO TO 70
               ELSE
*
*                 Only test 2 passed -- chase the zero to T(ILAST,ILAST)
*                 Then process as in the case T(ILAST,ILAST)=0
*
                  DO 50 jch = j, ilast - 1
                     temp = t( jch, jch+1 )
                     CALL dlartg( temp, t( jch+1, jch+1 ), c, s,
     $                            t( jch, jch+1 ) )
                     t( jch+1, jch+1 ) = zero
                     IF( jch.LT.ilastm-1 )
     $                  CALL drot( ilastm-jch-1, t( jch, jch+2 ), ldt,
     $                             t( jch+1, jch+2 ), ldt, c, s )
                     CALL drot( ilastm-jch+2, h( jch, jch-1 ), ldh,
     $                          h( jch+1, jch-1 ), ldh, c, s )
                     IF( ilq )
     $                  CALL drot( n, q( 1, jch ), 1, q( 1, jch+1 ), 1,
     $                             c, s )
                     temp = h( jch+1, jch )
                     CALL dlartg( temp, h( jch+1, jch-1 ), c, s,
     $                            h( jch+1, jch ) )
                     h( jch+1, jch-1 ) = zero
                     CALL drot( jch+1-ifrstm, h( ifrstm, jch ), 1,
     $                          h( ifrstm, jch-1 ), 1, c, s )
                     CALL drot( jch-ifrstm, t( ifrstm, jch ), 1,
     $                          t( ifrstm, jch-1 ), 1, c, s )
                     IF( ilz )
     $                  CALL drot( n, z( 1, jch ), 1, z( 1, jch-1 ), 1,
     $                             c, s )
   50             CONTINUE
                  GO TO 70
               END IF
            ELSE IF( ilazro ) THEN
*
*              Only test 1 passed -- work on J:ILAST
*
               ifirst = j
               GO TO 110
            END IF
*
*           Neither test passed -- try next J
*
   60    CONTINUE
*
*        (Drop-through is "impossible")
*
         info = n + 1
         GO TO 420
*
*        T(ILAST,ILAST)=0 -- clear H(ILAST,ILAST-1) to split off a
*        1x1 block.
*
   70    CONTINUE
         temp = h( ilast, ilast )
         CALL dlartg( temp, h( ilast, ilast-1 ), c, s,
     $                h( ilast, ilast ) )
         h( ilast, ilast-1 ) = zero
         CALL drot( ilast-ifrstm, h( ifrstm, ilast ), 1,
     $              h( ifrstm, ilast-1 ), 1, c, s )
         CALL drot( ilast-ifrstm, t( ifrstm, ilast ), 1,
     $              t( ifrstm, ilast-1 ), 1, c, s )
         IF( ilz )
     $      CALL drot( n, z( 1, ilast ), 1, z( 1, ilast-1 ), 1, c, s )
*
*        H(ILAST,ILAST-1)=0 -- Standardize B, set ALPHAR, ALPHAI,
*                              and BETA
*
   80    CONTINUE
         IF( t( ilast, ilast ).LT.zero ) THEN
            IF( ilschr ) THEN
               DO 90 j = ifrstm, ilast
                  h( j, ilast ) = -h( j, ilast )
                  t( j, ilast ) = -t( j, ilast )
   90          CONTINUE
            ELSE
               h( ilast, ilast ) = -h( ilast, ilast )
               t( ilast, ilast ) = -t( ilast, ilast )
            END IF
            IF( ilz ) THEN
               DO 100 j = 1, n
                  z( j, ilast ) = -z( j, ilast )
  100          CONTINUE
            END IF
         END IF
         alphar( ilast ) = h( ilast, ilast )
         alphai( ilast ) = zero
         beta( ilast ) = t( ilast, ilast )
*
*        Go to next block -- exit if finished.
*
         ilast = ilast - 1
         IF( ilast.LT.ilo )
     $      GO TO 380
*
*        Reset counters
*
         iiter = 0
         eshift = zero
         IF( .NOT.ilschr ) THEN
            ilastm = ilast
            IF( ifrstm.GT.ilast )
     $         ifrstm = ilo
         END IF
         GO TO 350
*
*        QZ step
*
*        This iteration only involves rows/columns IFIRST:ILAST. We
*        assume IFIRST < ILAST, and that the diagonal of B is non-zero.
*
  110    CONTINUE
         iiter = iiter + 1
         IF( .NOT.ilschr ) THEN
            ifrstm = ifirst
         END IF
*
*        Compute single shifts.
*
*        At this point, IFIRST < ILAST, and the diagonal elements of
*        T(IFIRST:ILAST,IFIRST,ILAST) are larger than BTOL (in
*        magnitude)
*
         IF( ( iiter / 10 )*10.EQ.iiter ) THEN
*
*           Exceptional shift.  Chosen for no particularly good reason.
*           (Single shift only.)
*
            IF( ( dble( maxit )*safmin )*abs( h( ilast, ilast-1 ) ).LT.
     $          abs( t( ilast-1, ilast-1 ) ) ) THEN
               eshift = h( ilast, ilast-1 ) /
     $                  t( ilast-1, ilast-1 )
            ELSE
               eshift = eshift + one / ( safmin*dble( maxit ) )
            END IF
            s1 = one
            wr = eshift
*
         ELSE
*
*           Shifts based on the generalized eigenvalues of the
*           bottom-right 2x2 block of A and B. The first eigenvalue
*           returned by DLAG2 is the Wilkinson shift (AEP p.512),
*
            CALL dlag2( h( ilast-1, ilast-1 ), ldh,
     $                  t( ilast-1, ilast-1 ), ldt, safmin*safety, s1,
     $                  s2, wr, wr2, wi )
*
            IF ( abs( (wr/s1)*t( ilast, ilast ) - h( ilast, ilast ) )
     $         .GT. abs( (wr2/s2)*t( ilast, ilast )
     $         - h( ilast, ilast ) ) ) THEN
               temp = wr
               wr = wr2
               wr2 = temp
               temp = s1
               s1 = s2
               s2 = temp
            END IF
            temp = max( s1, safmin*max( one, abs( wr ), abs( wi ) ) )
            IF( wi.NE.zero )
     $         GO TO 200
         END IF
*
*        Fiddle with shift to avoid overflow
*
         temp = min( ascale, one )*( half*safmax )
         IF( s1.GT.temp ) THEN
            scale = temp / s1
         ELSE
            scale = one
         END IF
*
         temp = min( bscale, one )*( half*safmax )
         IF( abs( wr ).GT.temp )
     $      scale = min( scale, temp / abs( wr ) )
         s1 = scale*s1
         wr = scale*wr
*
*        Now check for two consecutive small subdiagonals.
*
         DO 120 j = ilast - 1, ifirst + 1, -1
            istart = j
            temp = abs( s1*h( j, j-1 ) )
            temp2 = abs( s1*h( j, j )-wr*t( j, j ) )
            tempr = max( temp, temp2 )
            IF( tempr.LT.one .AND. tempr.NE.zero ) THEN
               temp = temp / tempr
               temp2 = temp2 / tempr
            END IF
            IF( abs( ( ascale*h( j+1, j ) )*temp ).LE.( ascale*atol )*
     $          temp2 )GO TO 130
  120    CONTINUE
*
         istart = ifirst
  130    CONTINUE
*
*        Do an implicit single-shift QZ sweep.
*
*        Initial Q
*
         temp = s1*h( istart, istart ) - wr*t( istart, istart )
         temp2 = s1*h( istart+1, istart )
         CALL dlartg( temp, temp2, c, s, tempr )
*
*        Sweep
*
         DO 190 j = istart, ilast - 1
            IF( j.GT.istart ) THEN
               temp = h( j, j-1 )
               CALL dlartg( temp, h( j+1, j-1 ), c, s, h( j, j-1 ) )
               h( j+1, j-1 ) = zero
            END IF
*
            DO 140 jc = j, ilastm
               temp = c*h( j, jc ) + s*h( j+1, jc )
               h( j+1, jc ) = -s*h( j, jc ) + c*h( j+1, jc )
               h( j, jc ) = temp
               temp2 = c*t( j, jc ) + s*t( j+1, jc )
               t( j+1, jc ) = -s*t( j, jc ) + c*t( j+1, jc )
               t( j, jc ) = temp2
  140       CONTINUE
            IF( ilq ) THEN
               DO 150 jr = 1, n
                  temp = c*q( jr, j ) + s*q( jr, j+1 )
                  q( jr, j+1 ) = -s*q( jr, j ) + c*q( jr, j+1 )
                  q( jr, j ) = temp
  150          CONTINUE
            END IF
*
            temp = t( j+1, j+1 )
            CALL dlartg( temp, t( j+1, j ), c, s, t( j+1, j+1 ) )
            t( j+1, j ) = zero
*
            DO 160 jr = ifrstm, min( j+2, ilast )
               temp = c*h( jr, j+1 ) + s*h( jr, j )
               h( jr, j ) = -s*h( jr, j+1 ) + c*h( jr, j )
               h( jr, j+1 ) = temp
  160       CONTINUE
            DO 170 jr = ifrstm, j
               temp = c*t( jr, j+1 ) + s*t( jr, j )
               t( jr, j ) = -s*t( jr, j+1 ) + c*t( jr, j )
               t( jr, j+1 ) = temp
  170       CONTINUE
            IF( ilz ) THEN
               DO 180 jr = 1, n
                  temp = c*z( jr, j+1 ) + s*z( jr, j )
                  z( jr, j ) = -s*z( jr, j+1 ) + c*z( jr, j )
                  z( jr, j+1 ) = temp
  180          CONTINUE
            END IF
  190    CONTINUE
*
         GO TO 350
*
*        Use Francis double-shift
*
*        Note: the Francis double-shift should work with real shifts,
*              but only if the block is at least 3x3.
*              This code may break if this point is reached with
*              a 2x2 block with real eigenvalues.
*
  200    CONTINUE
         IF( ifirst+1.EQ.ilast ) THEN
*
*           Special case -- 2x2 block with complex eigenvectors
*
*           Step 1: Standardize, that is, rotate so that
*
*                       ( B11  0  )
*                   B = (         )  with B11 non-negative.
*                       (  0  B22 )
*
            CALL dlasv2( t( ilast-1, ilast-1 ), t( ilast-1, ilast ),
     $                   t( ilast, ilast ), b22, b11, sr, cr, sl, cl )
*
            IF( b11.LT.zero ) THEN
               cr = -cr
               sr = -sr
               b11 = -b11
               b22 = -b22
            END IF
*
            CALL drot( ilastm+1-ifirst, h( ilast-1, ilast-1 ), ldh,
     $                 h( ilast, ilast-1 ), ldh, cl, sl )
            CALL drot( ilast+1-ifrstm, h( ifrstm, ilast-1 ), 1,
     $                 h( ifrstm, ilast ), 1, cr, sr )
*
            IF( ilast.LT.ilastm )
     $         CALL drot( ilastm-ilast, t( ilast-1, ilast+1 ), ldt,
     $                    t( ilast, ilast+1 ), ldt, cl, sl )
            IF( ifrstm.LT.ilast-1 )
     $         CALL drot( ifirst-ifrstm, t( ifrstm, ilast-1 ), 1,
     $                    t( ifrstm, ilast ), 1, cr, sr )
*
            IF( ilq )
     $         CALL drot( n, q( 1, ilast-1 ), 1, q( 1, ilast ), 1, cl,
     $                    sl )
            IF( ilz )
     $         CALL drot( n, z( 1, ilast-1 ), 1, z( 1, ilast ), 1, cr,
     $                    sr )
*
            t( ilast-1, ilast-1 ) = b11
            t( ilast-1, ilast ) = zero
            t( ilast, ilast-1 ) = zero
            t( ilast, ilast ) = b22
*
*           If B22 is negative, negate column ILAST
*
            IF( b22.LT.zero ) THEN
               DO 210 j = ifrstm, ilast
                  h( j, ilast ) = -h( j, ilast )
                  t( j, ilast ) = -t( j, ilast )
  210          CONTINUE
*
               IF( ilz ) THEN
                  DO 220 j = 1, n
                     z( j, ilast ) = -z( j, ilast )
  220             CONTINUE
               END IF
               b22 = -b22
            END IF
*
*           Step 2: Compute ALPHAR, ALPHAI, and BETA (see refs.)
*
*           Recompute shift
*
            CALL dlag2( h( ilast-1, ilast-1 ), ldh,
     $                  t( ilast-1, ilast-1 ), ldt, safmin*safety, s1,
     $                  temp, wr, temp2, wi )
*
*           If standardization has perturbed the shift onto real line,
*           do another (real single-shift) QR step.
*
            IF( wi.EQ.zero )
     $         GO TO 350
            s1inv = one / s1
*
*           Do EISPACK (QZVAL) computation of alpha and beta
*
            a11 = h( ilast-1, ilast-1 )
            a21 = h( ilast, ilast-1 )
            a12 = h( ilast-1, ilast )
            a22 = h( ilast, ilast )
*
*           Compute complex Givens rotation on right
*           (Assume some element of C = (sA - wB) > unfl )
*                            __
*           (sA - wB) ( CZ   -SZ )
*                     ( SZ    CZ )
*
            c11r = s1*a11 - wr*b11
            c11i = -wi*b11
            c12 = s1*a12
            c21 = s1*a21
            c22r = s1*a22 - wr*b22
            c22i = -wi*b22
*
            IF( abs( c11r )+abs( c11i )+abs( c12 ).GT.abs( c21 )+
     $          abs( c22r )+abs( c22i ) ) THEN
               t1 = dlapy3( c12, c11r, c11i )
               cz = c12 / t1
               szr = -c11r / t1
               szi = -c11i / t1
            ELSE
               cz = dlapy2( c22r, c22i )
               IF( cz.LE.safmin ) THEN
                  cz = zero
                  szr = one
                  szi = zero
               ELSE
                  tempr = c22r / cz
                  tempi = c22i / cz
                  t1 = dlapy2( cz, c21 )
                  cz = cz / t1
                  szr = -c21*tempr / t1
                  szi = c21*tempi / t1
               END IF
            END IF
*
*           Compute Givens rotation on left
*
*           (  CQ   SQ )
*           (  __      )  A or B
*           ( -SQ   CQ )
*
            an = abs( a11 ) + abs( a12 ) + abs( a21 ) + abs( a22 )
            bn = abs( b11 ) + abs( b22 )
            wabs = abs( wr ) + abs( wi )
            IF( s1*an.GT.wabs*bn ) THEN
               cq = cz*b11
               sqr = szr*b22
               sqi = -szi*b22
            ELSE
               a1r = cz*a11 + szr*a12
               a1i = szi*a12
               a2r = cz*a21 + szr*a22
               a2i = szi*a22
               cq = dlapy2( a1r, a1i )
               IF( cq.LE.safmin ) THEN
                  cq = zero
                  sqr = one
                  sqi = zero
               ELSE
                  tempr = a1r / cq
                  tempi = a1i / cq
                  sqr = tempr*a2r + tempi*a2i
                  sqi = tempi*a2r - tempr*a2i
               END IF
            END IF
            t1 = dlapy3( cq, sqr, sqi )
            cq = cq / t1
            sqr = sqr / t1
            sqi = sqi / t1
*
*           Compute diagonal elements of QBZ
*
            tempr = sqr*szr - sqi*szi
            tempi = sqr*szi + sqi*szr
            b1r = cq*cz*b11 + tempr*b22
            b1i = tempi*b22
            b1a = dlapy2( b1r, b1i )
            b2r = cq*cz*b22 + tempr*b11
            b2i = -tempi*b11
            b2a = dlapy2( b2r, b2i )
*
*           Normalize so beta > 0, and Im( alpha1 ) > 0
*
            beta( ilast-1 ) = b1a
            beta( ilast ) = b2a
            alphar( ilast-1 ) = ( wr*b1a )*s1inv
            alphai( ilast-1 ) = ( wi*b1a )*s1inv
            alphar( ilast ) = ( wr*b2a )*s1inv
            alphai( ilast ) = -( wi*b2a )*s1inv
*
*           Step 3: Go to next block -- exit if finished.
*
            ilast = ifirst - 1
            IF( ilast.LT.ilo )
     $         GO TO 380
*
*           Reset counters
*
            iiter = 0
            eshift = zero
            IF( .NOT.ilschr ) THEN
               ilastm = ilast
               IF( ifrstm.GT.ilast )
     $            ifrstm = ilo
            END IF
            GO TO 350
         ELSE
*
*           Usual case: 3x3 or larger block, using Francis implicit
*                       double-shift
*
*                                    2
*           Eigenvalue equation is  w  - c w + d = 0,
*
*                                         -1 2        -1
*           so compute 1st column of  (A B  )  - c A B   + d
*           using the formula in QZIT (from EISPACK)
*
*           We assume that the block is at least 3x3
*
            ad11 = ( ascale*h( ilast-1, ilast-1 ) ) /
     $             ( bscale*t( ilast-1, ilast-1 ) )
            ad21 = ( ascale*h( ilast, ilast-1 ) ) /
     $             ( bscale*t( ilast-1, ilast-1 ) )
            ad12 = ( ascale*h( ilast-1, ilast ) ) /
     $             ( bscale*t( ilast, ilast ) )
            ad22 = ( ascale*h( ilast, ilast ) ) /
     $             ( bscale*t( ilast, ilast ) )
            u12 = t( ilast-1, ilast ) / t( ilast, ilast )
            ad11l = ( ascale*h( ifirst, ifirst ) ) /
     $              ( bscale*t( ifirst, ifirst ) )
            ad21l = ( ascale*h( ifirst+1, ifirst ) ) /
     $              ( bscale*t( ifirst, ifirst ) )
            ad12l = ( ascale*h( ifirst, ifirst+1 ) ) /
     $              ( bscale*t( ifirst+1, ifirst+1 ) )
            ad22l = ( ascale*h( ifirst+1, ifirst+1 ) ) /
     $              ( bscale*t( ifirst+1, ifirst+1 ) )
            ad32l = ( ascale*h( ifirst+2, ifirst+1 ) ) /
     $              ( bscale*t( ifirst+1, ifirst+1 ) )
            u12l = t( ifirst, ifirst+1 ) / t( ifirst+1, ifirst+1 )
*
            v( 1 ) = ( ad11-ad11l )*( ad22-ad11l ) - ad12*ad21 +
     $               ad21*u12*ad11l + ( ad12l-ad11l*u12l )*ad21l
            v( 2 ) = ( ( ad22l-ad11l )-ad21l*u12l-( ad11-ad11l )-
     $               ( ad22-ad11l )+ad21*u12 )*ad21l
            v( 3 ) = ad32l*ad21l
*
            istart = ifirst
*
            CALL dlarfg( 3, v( 1 ), v( 2 ), 1, tau )
            v( 1 ) = one
*
*           Sweep
*
            DO 290 j = istart, ilast - 2
*
*              All but last elements: use 3x3 Householder transforms.
*
*              Zero (j-1)st column of A
*
               IF( j.GT.istart ) THEN
                  v( 1 ) = h( j, j-1 )
                  v( 2 ) = h( j+1, j-1 )
                  v( 3 ) = h( j+2, j-1 )
*
                  CALL dlarfg( 3, h( j, j-1 ), v( 2 ), 1, tau )
                  v( 1 ) = one
                  h( j+1, j-1 ) = zero
                  h( j+2, j-1 ) = zero
               END IF
*
               DO 230 jc = j, ilastm
                  temp = tau*( h( j, jc )+v( 2 )*h( j+1, jc )+v( 3 )*
     $                   h( j+2, jc ) )
                  h( j, jc ) = h( j, jc ) - temp
                  h( j+1, jc ) = h( j+1, jc ) - temp*v( 2 )
                  h( j+2, jc ) = h( j+2, jc ) - temp*v( 3 )
                  temp2 = tau*( t( j, jc )+v( 2 )*t( j+1, jc )+v( 3 )*
     $                    t( j+2, jc ) )
                  t( j, jc ) = t( j, jc ) - temp2
                  t( j+1, jc ) = t( j+1, jc ) - temp2*v( 2 )
                  t( j+2, jc ) = t( j+2, jc ) - temp2*v( 3 )
  230          CONTINUE
               IF( ilq ) THEN
                  DO 240 jr = 1, n
                     temp = tau*( q( jr, j )+v( 2 )*q( jr, j+1 )+v( 3 )*
     $                      q( jr, j+2 ) )
                     q( jr, j ) = q( jr, j ) - temp
                     q( jr, j+1 ) = q( jr, j+1 ) - temp*v( 2 )
                     q( jr, j+2 ) = q( jr, j+2 ) - temp*v( 3 )
  240             CONTINUE
               END IF
*
*              Zero j-th column of B (see DLAGBC for details)
*
*              Swap rows to pivot
*
               ilpivt = .false.
               temp = max( abs( t( j+1, j+1 ) ), abs( t( j+1, j+2 ) ) )
               temp2 = max( abs( t( j+2, j+1 ) ), abs( t( j+2, j+2 ) ) )
               IF( max( temp, temp2 ).LT.safmin ) THEN
                  scale = zero
                  u1 = one
                  u2 = zero
                  GO TO 250
               ELSE IF( temp.GE.temp2 ) THEN
                  w11 = t( j+1, j+1 )
                  w21 = t( j+2, j+1 )
                  w12 = t( j+1, j+2 )
                  w22 = t( j+2, j+2 )
                  u1 = t( j+1, j )
                  u2 = t( j+2, j )
               ELSE
                  w21 = t( j+1, j+1 )
                  w11 = t( j+2, j+1 )
                  w22 = t( j+1, j+2 )
                  w12 = t( j+2, j+2 )
                  u2 = t( j+1, j )
                  u1 = t( j+2, j )
               END IF
*
*              Swap columns if nec.
*
               IF( abs( w12 ).GT.abs( w11 ) ) THEN
                  ilpivt = .true.
                  temp = w12
                  temp2 = w22
                  w12 = w11
                  w22 = w21
                  w11 = temp
                  w21 = temp2
               END IF
*
*              LU-factor
*
               temp = w21 / w11
               u2 = u2 - temp*u1
               w22 = w22 - temp*w12
               w21 = zero
*
*              Compute SCALE
*
               scale = one
               IF( abs( w22 ).LT.safmin ) THEN
                  scale = zero
                  u2 = one
                  u1 = -w12 / w11
                  GO TO 250
               END IF
               IF( abs( w22 ).LT.abs( u2 ) )
     $            scale = abs( w22 / u2 )
               IF( abs( w11 ).LT.abs( u1 ) )
     $            scale = min( scale, abs( w11 / u1 ) )
*
*              Solve
*
               u2 = ( scale*u2 ) / w22
               u1 = ( scale*u1-w12*u2 ) / w11
*
  250          CONTINUE
               IF( ilpivt ) THEN
                  temp = u2
                  u2 = u1
                  u1 = temp
               END IF
*
*              Compute Householder Vector
*
               t1 = sqrt( scale**2+u1**2+u2**2 )
               tau = one + scale / t1
               vs = -one / ( scale+t1 )
               v( 1 ) = one
               v( 2 ) = vs*u1
               v( 3 ) = vs*u2
*
*              Apply transformations from the right.
*
               DO 260 jr = ifrstm, min( j+3, ilast )
                  temp = tau*( h( jr, j )+v( 2 )*h( jr, j+1 )+v( 3 )*
     $                   h( jr, j+2 ) )
                  h( jr, j ) = h( jr, j ) - temp
                  h( jr, j+1 ) = h( jr, j+1 ) - temp*v( 2 )
                  h( jr, j+2 ) = h( jr, j+2 ) - temp*v( 3 )
  260          CONTINUE
               DO 270 jr = ifrstm, j + 2
                  temp = tau*( t( jr, j )+v( 2 )*t( jr, j+1 )+v( 3 )*
     $                   t( jr, j+2 ) )
                  t( jr, j ) = t( jr, j ) - temp
                  t( jr, j+1 ) = t( jr, j+1 ) - temp*v( 2 )
                  t( jr, j+2 ) = t( jr, j+2 ) - temp*v( 3 )
  270          CONTINUE
               IF( ilz ) THEN
                  DO 280 jr = 1, n
                     temp = tau*( z( jr, j )+v( 2 )*z( jr, j+1 )+v( 3 )*
     $                      z( jr, j+2 ) )
                     z( jr, j ) = z( jr, j ) - temp
                     z( jr, j+1 ) = z( jr, j+1 ) - temp*v( 2 )
                     z( jr, j+2 ) = z( jr, j+2 ) - temp*v( 3 )
  280             CONTINUE
               END IF
               t( j+1, j ) = zero
               t( j+2, j ) = zero
  290       CONTINUE
*
*           Last elements: Use Givens rotations
*
*           Rotations from the left
*
            j = ilast - 1
            temp = h( j, j-1 )
            CALL dlartg( temp, h( j+1, j-1 ), c, s, h( j, j-1 ) )
            h( j+1, j-1 ) = zero
*
            DO 300 jc = j, ilastm
               temp = c*h( j, jc ) + s*h( j+1, jc )
               h( j+1, jc ) = -s*h( j, jc ) + c*h( j+1, jc )
               h( j, jc ) = temp
               temp2 = c*t( j, jc ) + s*t( j+1, jc )
               t( j+1, jc ) = -s*t( j, jc ) + c*t( j+1, jc )
               t( j, jc ) = temp2
  300       CONTINUE
            IF( ilq ) THEN
               DO 310 jr = 1, n
                  temp = c*q( jr, j ) + s*q( jr, j+1 )
                  q( jr, j+1 ) = -s*q( jr, j ) + c*q( jr, j+1 )
                  q( jr, j ) = temp
  310          CONTINUE
            END IF
*
*           Rotations from the right.
*
            temp = t( j+1, j+1 )
            CALL dlartg( temp, t( j+1, j ), c, s, t( j+1, j+1 ) )
            t( j+1, j ) = zero
*
            DO 320 jr = ifrstm, ilast
               temp = c*h( jr, j+1 ) + s*h( jr, j )
               h( jr, j ) = -s*h( jr, j+1 ) + c*h( jr, j )
               h( jr, j+1 ) = temp
  320       CONTINUE
            DO 330 jr = ifrstm, ilast - 1
               temp = c*t( jr, j+1 ) + s*t( jr, j )
               t( jr, j ) = -s*t( jr, j+1 ) + c*t( jr, j )
               t( jr, j+1 ) = temp
  330       CONTINUE
            IF( ilz ) THEN
               DO 340 jr = 1, n
                  temp = c*z( jr, j+1 ) + s*z( jr, j )
                  z( jr, j ) = -s*z( jr, j+1 ) + c*z( jr, j )
                  z( jr, j+1 ) = temp
  340          CONTINUE
            END IF
*
*           End of Double-Shift code
*
         END IF
*
         GO TO 350
*
*        End of iteration loop
*
  350    CONTINUE
  360 CONTINUE
*
*     Drop-through = non-convergence
*
      info = ilast
      GO TO 420
*
*     Successful completion of all QZ steps
*
  380 CONTINUE
*
*     Set Eigenvalues 1:ILO-1
*
      DO 410 j = 1, ilo - 1
         IF( t( j, j ).LT.zero ) THEN
            IF( ilschr ) THEN
               DO 390 jr = 1, j
                  h( jr, j ) = -h( jr, j )
                  t( jr, j ) = -t( jr, j )
  390          CONTINUE
            ELSE
               h( j, j ) = -h( j, j )
               t( j, j ) = -t( j, j )
            END IF
            IF( ilz ) THEN
               DO 400 jr = 1, n
                  z( jr, j ) = -z( jr, j )
  400          CONTINUE
            END IF
         END IF
         alphar( j ) = h( j, j )
         alphai( j ) = zero
         beta( j ) = t( j, j )
  410 CONTINUE
*
*     Normal Termination
*
      info = 0
*
*     Exit (other than argument error) -- return optimal workspace size
*
  420 CONTINUE
      work( 1 ) = dble( n )
      RETURN
*
*     End of DHGEQZ
*
      END

! DTGEVC
      SUBROUTINE dtgevc( SIDE, HOWMNY, SELECT, N, S, LDS, P, LDP, VL,
     $                   LDVL, VR, LDVR, MM, M, WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          HOWMNY, SIDE
      INTEGER            INFO, LDP, LDS, LDVL, LDVR, M, MM, N
*     ..
*     .. Array Arguments ..
      LOGICAL            SELECT( * )
      DOUBLE PRECISION   P( LDP, * ), S( LDS, * ), VL( LDVL, * ),
     $                   vr( ldvr, * ), work( * )
*     ..
*
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, SAFETY
      parameter( zero = 0.0d+0, one = 1.0d+0,
     $                   safety = 1.0d+2 )
*     ..
*     .. Local Scalars ..
      LOGICAL            COMPL, COMPR, IL2BY2, ILABAD, ILALL, ILBACK,
     $                   ilbbad, ilcomp, ilcplx, lsa, lsb
      INTEGER            I, IBEG, IEIG, IEND, IHWMNY, IINFO, IM, ISIDE,
     $                   j, ja, jc, je, jr, jw, na, nw
      DOUBLE PRECISION   ACOEF, ACOEFA, ANORM, ASCALE, BCOEFA, BCOEFI,
     $                   bcoefr, big, bignum, bnorm, bscale, cim2a,
     $                   cim2b, cimaga, cimagb, cre2a, cre2b, creala,
     $                   crealb, dmin, safmin, salfar, sbeta, scale,
     $                   small, temp, temp2, temp2i, temp2r, ulp, xmax,
     $                   xscale
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   BDIAG( 2 ), SUM( 2, 2 ), SUMS( 2, 2 ),
     $                   sump( 2, 2 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           lsame, dlamch
*     ..
*     .. External Subroutines ..
      EXTERNAL           dgemv, dlabad, dlacpy, dlag2, dlaln2, xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, max, min
*     ..
*     .. Executable Statements ..
*
*     Decode and Test the input parameters
*
      IF( lsame( howmny, 'A' ) ) THEN
         ihwmny = 1
         ilall = .true.
         ilback = .false.
      ELSE IF( lsame( howmny, 'S' ) ) THEN
         ihwmny = 2
         ilall = .false.
         ilback = .false.
      ELSE IF( lsame( howmny, 'B' ) ) THEN
         ihwmny = 3
         ilall = .true.
         ilback = .true.
      ELSE
         ihwmny = -1
         ilall = .true.
      END IF
*
      IF( lsame( side, 'R' ) ) THEN
         iside = 1
         compl = .false.
         compr = .true.
      ELSE IF( lsame( side, 'L' ) ) THEN
         iside = 2
         compl = .true.
         compr = .false.
      ELSE IF( lsame( side, 'B' ) ) THEN
         iside = 3
         compl = .true.
         compr = .true.
      ELSE
         iside = -1
      END IF
*
      info = 0
      IF( iside.LT.0 ) THEN
         info = -1
      ELSE IF( ihwmny.LT.0 ) THEN
         info = -2
      ELSE IF( n.LT.0 ) THEN
         info = -4
      ELSE IF( lds.LT.max( 1, n ) ) THEN
         info = -6
      ELSE IF( ldp.LT.max( 1, n ) ) THEN
         info = -8
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'DTGEVC', -info )
         RETURN
      END IF
*
*     Count the number of eigenvectors to be computed
*
      IF( .NOT.ilall ) THEN
         im = 0
         ilcplx = .false.
         DO 10 j = 1, n
            IF( ilcplx ) THEN
               ilcplx = .false.
               GO TO 10
            END IF
            IF( j.LT.n ) THEN
               IF( s( j+1, j ).NE.zero )
     $            ilcplx = .true.
            END IF
            IF( ilcplx ) THEN
               IF( SELECT( j ) .OR. SELECT( j+1 ) )
     $            im = im + 2
            ELSE
               IF( SELECT( j ) )
     $            im = im + 1
            END IF
   10    CONTINUE
      ELSE
         im = n
      END IF
*
*     Check 2-by-2 diagonal blocks of A, B
*
      ilabad = .false.
      ilbbad = .false.
      DO 20 j = 1, n - 1
         IF( s( j+1, j ).NE.zero ) THEN
            IF( p( j, j ).EQ.zero .OR. p( j+1, j+1 ).EQ.zero .OR.
     $          p( j, j+1 ).NE.zero )ilbbad = .true.
            IF( j.LT.n-1 ) THEN
               IF( s( j+2, j+1 ).NE.zero )
     $            ilabad = .true.
            END IF
         END IF
   20 CONTINUE
*
      IF( ilabad ) THEN
         info = -5
      ELSE IF( ilbbad ) THEN
         info = -7
      ELSE IF( compl .AND. ldvl.LT.n .OR. ldvl.LT.1 ) THEN
         info = -10
      ELSE IF( compr .AND. ldvr.LT.n .OR. ldvr.LT.1 ) THEN
         info = -12
      ELSE IF( mm.LT.im ) THEN
         info = -13
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'DTGEVC', -info )
         RETURN
      END IF
*
*     Quick return if possible
*
      m = im
      IF( n.EQ.0 )
     $   RETURN
*
*     Machine Constants
*
      safmin = dlamch( 'Safe minimum' )
      big = one / safmin
      CALL dlabad( safmin, big )
      ulp = dlamch( 'Epsilon' )*dlamch( 'Base' )
      small = safmin*n / ulp
      big = one / small
      bignum = one / ( safmin*n )
*
*     Compute the 1-norm of each column of the strictly upper triangular
*     part (i.e., excluding all elements belonging to the diagonal
*     blocks) of A and B to check for possible overflow in the
*     triangular solver.
*
      anorm = abs( s( 1, 1 ) )
      IF( n.GT.1 )
     $   anorm = anorm + abs( s( 2, 1 ) )
      bnorm = abs( p( 1, 1 ) )
      work( 1 ) = zero
      work( n+1 ) = zero
*
      DO 50 j = 2, n
         temp = zero
         temp2 = zero
         IF( s( j, j-1 ).EQ.zero ) THEN
            iend = j - 1
         ELSE
            iend = j - 2
         END IF
         DO 30 i = 1, iend
            temp = temp + abs( s( i, j ) )
            temp2 = temp2 + abs( p( i, j ) )
   30    CONTINUE
         work( j ) = temp
         work( n+j ) = temp2
         DO 40 i = iend + 1, min( j+1, n )
            temp = temp + abs( s( i, j ) )
            temp2 = temp2 + abs( p( i, j ) )
   40    CONTINUE
         anorm = max( anorm, temp )
         bnorm = max( bnorm, temp2 )
   50 CONTINUE
*
      ascale = one / max( anorm, safmin )
      bscale = one / max( bnorm, safmin )
*
*     Left eigenvectors
*
      IF( compl ) THEN
         ieig = 0
*
*        Main loop over eigenvalues
*
         ilcplx = .false.
         DO 220 je = 1, n
*
*           Skip this iteration if (a) HOWMNY='S' and SELECT=.FALSE., or
*           (b) this would be the second of a complex pair.
*           Check for complex eigenvalue, so as to be sure of which
*           entry(-ies) of SELECT to look at.
*
            IF( ilcplx ) THEN
               ilcplx = .false.
               GO TO 220
            END IF
            nw = 1
            IF( je.LT.n ) THEN
               IF( s( je+1, je ).NE.zero ) THEN
                  ilcplx = .true.
                  nw = 2
               END IF
            END IF
            IF( ilall ) THEN
               ilcomp = .true.
            ELSE IF( ilcplx ) THEN
               ilcomp = SELECT( je ) .OR. SELECT( je+1 )
            ELSE
               ilcomp = SELECT( je )
            END IF
            IF( .NOT.ilcomp )
     $         GO TO 220
*
*           Decide if (a) singular pencil, (b) real eigenvalue, or
*           (c) complex eigenvalue.
*
            IF( .NOT.ilcplx ) THEN
               IF( abs( s( je, je ) ).LE.safmin .AND.
     $             abs( p( je, je ) ).LE.safmin ) THEN
*
*                 Singular matrix pencil -- return unit eigenvector
*
                  ieig = ieig + 1
                  DO 60 jr = 1, n
                     vl( jr, ieig ) = zero
   60             CONTINUE
                  vl( ieig, ieig ) = one
                  GO TO 220
               END IF
            END IF
*
*           Clear vector
*
            DO 70 jr = 1, nw*n
               work( 2*n+jr ) = zero
   70       CONTINUE
*                                                 T
*           Compute coefficients in  ( a A - b B )  y = 0
*              a  is  ACOEF
*              b  is  BCOEFR + i*BCOEFI
*
            IF( .NOT.ilcplx ) THEN
*
*              Real eigenvalue
*
               temp = one / max( abs( s( je, je ) )*ascale,
     $                abs( p( je, je ) )*bscale, safmin )
               salfar = ( temp*s( je, je ) )*ascale
               sbeta = ( temp*p( je, je ) )*bscale
               acoef = sbeta*ascale
               bcoefr = salfar*bscale
               bcoefi = zero
*
*              Scale to avoid underflow
*
               scale = one
               lsa = abs( sbeta ).GE.safmin .AND. abs( acoef ).LT.small
               lsb = abs( salfar ).GE.safmin .AND. abs( bcoefr ).LT.
     $               small
               IF( lsa )
     $            scale = ( small / abs( sbeta ) )*min( anorm, big )
               IF( lsb )
     $            scale = max( scale, ( small / abs( salfar ) )*
     $                    min( bnorm, big ) )
               IF( lsa .OR. lsb ) THEN
                  scale = min( scale, one /
     $                    ( safmin*max( one, abs( acoef ),
     $                    abs( bcoefr ) ) ) )
                  IF( lsa ) THEN
                     acoef = ascale*( scale*sbeta )
                  ELSE
                     acoef = scale*acoef
                  END IF
                  IF( lsb ) THEN
                     bcoefr = bscale*( scale*salfar )
                  ELSE
                     bcoefr = scale*bcoefr
                  END IF
               END IF
               acoefa = abs( acoef )
               bcoefa = abs( bcoefr )
*
*              First component is 1
*
               work( 2*n+je ) = one
               xmax = one
            ELSE
*
*              Complex eigenvalue
*
               CALL dlag2( s( je, je ), lds, p( je, je ), ldp,
     $                     safmin*safety, acoef, temp, bcoefr, temp2,
     $                     bcoefi )
               bcoefi = -bcoefi
               IF( bcoefi.EQ.zero ) THEN
                  info = je
                  RETURN
               END IF
*
*              Scale to avoid over/underflow
*
               acoefa = abs( acoef )
               bcoefa = abs( bcoefr ) + abs( bcoefi )
               scale = one
               IF( acoefa*ulp.LT.safmin .AND. acoefa.GE.safmin )
     $            scale = ( safmin / ulp ) / acoefa
               IF( bcoefa*ulp.LT.safmin .AND. bcoefa.GE.safmin )
     $            scale = max( scale, ( safmin / ulp ) / bcoefa )
               IF( safmin*acoefa.GT.ascale )
     $            scale = ascale / ( safmin*acoefa )
               IF( safmin*bcoefa.GT.bscale )
     $            scale = min( scale, bscale / ( safmin*bcoefa ) )
               IF( scale.NE.one ) THEN
                  acoef = scale*acoef
                  acoefa = abs( acoef )
                  bcoefr = scale*bcoefr
                  bcoefi = scale*bcoefi
                  bcoefa = abs( bcoefr ) + abs( bcoefi )
               END IF
*
*              Compute first two components of eigenvector
*
               temp = acoef*s( je+1, je )
               temp2r = acoef*s( je, je ) - bcoefr*p( je, je )
               temp2i = -bcoefi*p( je, je )
               IF( abs( temp ).GT.abs( temp2r )+abs( temp2i ) ) THEN
                  work( 2*n+je ) = one
                  work( 3*n+je ) = zero
                  work( 2*n+je+1 ) = -temp2r / temp
                  work( 3*n+je+1 ) = -temp2i / temp
               ELSE
                  work( 2*n+je+1 ) = one
                  work( 3*n+je+1 ) = zero
                  temp = acoef*s( je, je+1 )
                  work( 2*n+je ) = ( bcoefr*p( je+1, je+1 )-acoef*
     $                             s( je+1, je+1 ) ) / temp
                  work( 3*n+je ) = bcoefi*p( je+1, je+1 ) / temp
               END IF
               xmax = max( abs( work( 2*n+je ) )+abs( work( 3*n+je ) ),
     $                abs( work( 2*n+je+1 ) )+abs( work( 3*n+je+1 ) ) )
            END IF
*
            dmin = max( ulp*acoefa*anorm, ulp*bcoefa*bnorm, safmin )
*
*                                           T
*           Triangular solve of  (a A - b B)  y = 0
*
*                                   T
*           (rowwise in  (a A - b B) , or columnwise in (a A - b B) )
*
            il2by2 = .false.
*
            DO 160 j = je + nw, n
               IF( il2by2 ) THEN
                  il2by2 = .false.
                  GO TO 160
               END IF
*
               na = 1
               bdiag( 1 ) = p( j, j )
               IF( j.LT.n ) THEN
                  IF( s( j+1, j ).NE.zero ) THEN
                     il2by2 = .true.
                     bdiag( 2 ) = p( j+1, j+1 )
                     na = 2
                  END IF
               END IF
*
*              Check whether scaling is necessary for dot products
*
               xscale = one / max( one, xmax )
               temp = max( work( j ), work( n+j ),
     $                acoefa*work( j )+bcoefa*work( n+j ) )
               IF( il2by2 )
     $            temp = max( temp, work( j+1 ), work( n+j+1 ),
     $                   acoefa*work( j+1 )+bcoefa*work( n+j+1 ) )
               IF( temp.GT.bignum*xscale ) THEN
                  DO 90 jw = 0, nw - 1
                     DO 80 jr = je, j - 1
                        work( ( jw+2 )*n+jr ) = xscale*
     $                     work( ( jw+2 )*n+jr )
   80                CONTINUE
   90             CONTINUE
                  xmax = xmax*xscale
               END IF
*
*              Compute dot products
*
*                    j-1
*              SUM = sum  conjg( a*S(k,j) - b*P(k,j) )*x(k)
*                    k=je
*
*              To reduce the op count, this is done as
*
*              _        j-1                  _        j-1
*              a*conjg( sum  S(k,j)*x(k) ) - b*conjg( sum  P(k,j)*x(k) )
*                       k=je                          k=je
*
*              which may cause underflow problems if A or B are close
*              to underflow.  (E.g., less than SMALL.)
*
*
               DO 120 jw = 1, nw
                  DO 110 ja = 1, na
                     sums( ja, jw ) = zero
                     sump( ja, jw ) = zero
*
                     DO 100 jr = je, j - 1
                        sums( ja, jw ) = sums( ja, jw ) +
     $                                   s( jr, j+ja-1 )*
     $                                   work( ( jw+1 )*n+jr )
                        sump( ja, jw ) = sump( ja, jw ) +
     $                                   p( jr, j+ja-1 )*
     $                                   work( ( jw+1 )*n+jr )
  100                CONTINUE
  110             CONTINUE
  120          CONTINUE
*
               DO 130 ja = 1, na
                  IF( ilcplx ) THEN
                     sum( ja, 1 ) = -acoef*sums( ja, 1 ) +
     $                              bcoefr*sump( ja, 1 ) -
     $                              bcoefi*sump( ja, 2 )
                     sum( ja, 2 ) = -acoef*sums( ja, 2 ) +
     $                              bcoefr*sump( ja, 2 ) +
     $                              bcoefi*sump( ja, 1 )
                  ELSE
                     sum( ja, 1 ) = -acoef*sums( ja, 1 ) +
     $                              bcoefr*sump( ja, 1 )
                  END IF
  130          CONTINUE
*
*                                  T
*              Solve  ( a A - b B )  y = SUM(,)
*              with scaling and perturbation of the denominator
*
               CALL dlaln2( .true., na, nw, dmin, acoef, s( j, j ), lds,
     $                      bdiag( 1 ), bdiag( 2 ), sum, 2, bcoefr,
     $                      bcoefi, work( 2*n+j ), n, scale, temp,
     $                      iinfo )
               IF( scale.LT.one ) THEN
                  DO 150 jw = 0, nw - 1
                     DO 140 jr = je, j - 1
                        work( ( jw+2 )*n+jr ) = scale*
     $                     work( ( jw+2 )*n+jr )
  140                CONTINUE
  150             CONTINUE
                  xmax = scale*xmax
               END IF
               xmax = max( xmax, temp )
  160       CONTINUE
*
*           Copy eigenvector to VL, back transforming if
*           HOWMNY='B'.
*
            ieig = ieig + 1
            IF( ilback ) THEN
               DO 170 jw = 0, nw - 1
                  CALL dgemv( 'N', n, n+1-je, one, vl( 1, je ), ldvl,
     $                        work( ( jw+2 )*n+je ), 1, zero,
     $                        work( ( jw+4 )*n+1 ), 1 )
  170          CONTINUE
               CALL dlacpy( ' ', n, nw, work( 4*n+1 ), n, vl( 1, je ),
     $                      ldvl )
               ibeg = 1
            ELSE
               CALL dlacpy( ' ', n, nw, work( 2*n+1 ), n, vl( 1, ieig ),
     $                      ldvl )
               ibeg = je
            END IF
*
*           Scale eigenvector
*
            xmax = zero
            IF( ilcplx ) THEN
               DO 180 j = ibeg, n
                  xmax = max( xmax, abs( vl( j, ieig ) )+
     $                   abs( vl( j, ieig+1 ) ) )
  180          CONTINUE
            ELSE
               DO 190 j = ibeg, n
                  xmax = max( xmax, abs( vl( j, ieig ) ) )
  190          CONTINUE
            END IF
*
            IF( xmax.GT.safmin ) THEN
               xscale = one / xmax
*
               DO 210 jw = 0, nw - 1
                  DO 200 jr = ibeg, n
                     vl( jr, ieig+jw ) = xscale*vl( jr, ieig+jw )
  200             CONTINUE
  210          CONTINUE
            END IF
            ieig = ieig + nw - 1
*
  220    CONTINUE
      END IF
*
*     Right eigenvectors
*
      IF( compr ) THEN
         ieig = im + 1
*
*        Main loop over eigenvalues
*
         ilcplx = .false.
         DO 500 je = n, 1, -1
*
*           Skip this iteration if (a) HOWMNY='S' and SELECT=.FALSE., or
*           (b) this would be the second of a complex pair.
*           Check for complex eigenvalue, so as to be sure of which
*           entry(-ies) of SELECT to look at -- if complex, SELECT(JE)
*           or SELECT(JE-1).
*           If this is a complex pair, the 2-by-2 diagonal block
*           corresponding to the eigenvalue is in rows/columns JE-1:JE
*
            IF( ilcplx ) THEN
               ilcplx = .false.
               GO TO 500
            END IF
            nw = 1
            IF( je.GT.1 ) THEN
               IF( s( je, je-1 ).NE.zero ) THEN
                  ilcplx = .true.
                  nw = 2
               END IF
            END IF
            IF( ilall ) THEN
               ilcomp = .true.
            ELSE IF( ilcplx ) THEN
               ilcomp = SELECT( je ) .OR. SELECT( je-1 )
            ELSE
               ilcomp = SELECT( je )
            END IF
            IF( .NOT.ilcomp )
     $         GO TO 500
*
*           Decide if (a) singular pencil, (b) real eigenvalue, or
*           (c) complex eigenvalue.
*
            IF( .NOT.ilcplx ) THEN
               IF( abs( s( je, je ) ).LE.safmin .AND.
     $             abs( p( je, je ) ).LE.safmin ) THEN
*
*                 Singular matrix pencil -- unit eigenvector
*
                  ieig = ieig - 1
                  DO 230 jr = 1, n
                     vr( jr, ieig ) = zero
  230             CONTINUE
                  vr( ieig, ieig ) = one
                  GO TO 500
               END IF
            END IF
*
*           Clear vector
*
            DO 250 jw = 0, nw - 1
               DO 240 jr = 1, n
                  work( ( jw+2 )*n+jr ) = zero
  240          CONTINUE
  250       CONTINUE
*
*           Compute coefficients in  ( a A - b B ) x = 0
*              a  is  ACOEF
*              b  is  BCOEFR + i*BCOEFI
*
            IF( .NOT.ilcplx ) THEN
*
*              Real eigenvalue
*
               temp = one / max( abs( s( je, je ) )*ascale,
     $                abs( p( je, je ) )*bscale, safmin )
               salfar = ( temp*s( je, je ) )*ascale
               sbeta = ( temp*p( je, je ) )*bscale
               acoef = sbeta*ascale
               bcoefr = salfar*bscale
               bcoefi = zero
*
*              Scale to avoid underflow
*
               scale = one
               lsa = abs( sbeta ).GE.safmin .AND. abs( acoef ).LT.small
               lsb = abs( salfar ).GE.safmin .AND. abs( bcoefr ).LT.
     $               small
               IF( lsa )
     $            scale = ( small / abs( sbeta ) )*min( anorm, big )
               IF( lsb )
     $            scale = max( scale, ( small / abs( salfar ) )*
     $                    min( bnorm, big ) )
               IF( lsa .OR. lsb ) THEN
                  scale = min( scale, one /
     $                    ( safmin*max( one, abs( acoef ),
     $                    abs( bcoefr ) ) ) )
                  IF( lsa ) THEN
                     acoef = ascale*( scale*sbeta )
                  ELSE
                     acoef = scale*acoef
                  END IF
                  IF( lsb ) THEN
                     bcoefr = bscale*( scale*salfar )
                  ELSE
                     bcoefr = scale*bcoefr
                  END IF
               END IF
               acoefa = abs( acoef )
               bcoefa = abs( bcoefr )
*
*              First component is 1
*
               work( 2*n+je ) = one
               xmax = one
*
*              Compute contribution from column JE of A and B to sum
*              (See "Further Details", above.)
*
               DO 260 jr = 1, je - 1
                  work( 2*n+jr ) = bcoefr*p( jr, je ) -
     $                             acoef*s( jr, je )
  260          CONTINUE
            ELSE
*
*              Complex eigenvalue
*
               CALL dlag2( s( je-1, je-1 ), lds, p( je-1, je-1 ), ldp,
     $                     safmin*safety, acoef, temp, bcoefr, temp2,
     $                     bcoefi )
               IF( bcoefi.EQ.zero ) THEN
                  info = je - 1
                  RETURN
               END IF
*
*              Scale to avoid over/underflow
*
               acoefa = abs( acoef )
               bcoefa = abs( bcoefr ) + abs( bcoefi )
               scale = one
               IF( acoefa*ulp.LT.safmin .AND. acoefa.GE.safmin )
     $            scale = ( safmin / ulp ) / acoefa
               IF( bcoefa*ulp.LT.safmin .AND. bcoefa.GE.safmin )
     $            scale = max( scale, ( safmin / ulp ) / bcoefa )
               IF( safmin*acoefa.GT.ascale )
     $            scale = ascale / ( safmin*acoefa )
               IF( safmin*bcoefa.GT.bscale )
     $            scale = min( scale, bscale / ( safmin*bcoefa ) )
               IF( scale.NE.one ) THEN
                  acoef = scale*acoef
                  acoefa = abs( acoef )
                  bcoefr = scale*bcoefr
                  bcoefi = scale*bcoefi
                  bcoefa = abs( bcoefr ) + abs( bcoefi )
               END IF
*
*              Compute first two components of eigenvector
*              and contribution to sums
*
               temp = acoef*s( je, je-1 )
               temp2r = acoef*s( je, je ) - bcoefr*p( je, je )
               temp2i = -bcoefi*p( je, je )
               IF( abs( temp ).GE.abs( temp2r )+abs( temp2i ) ) THEN
                  work( 2*n+je ) = one
                  work( 3*n+je ) = zero
                  work( 2*n+je-1 ) = -temp2r / temp
                  work( 3*n+je-1 ) = -temp2i / temp
               ELSE
                  work( 2*n+je-1 ) = one
                  work( 3*n+je-1 ) = zero
                  temp = acoef*s( je-1, je )
                  work( 2*n+je ) = ( bcoefr*p( je-1, je-1 )-acoef*
     $                             s( je-1, je-1 ) ) / temp
                  work( 3*n+je ) = bcoefi*p( je-1, je-1 ) / temp
               END IF
*
               xmax = max( abs( work( 2*n+je ) )+abs( work( 3*n+je ) ),
     $                abs( work( 2*n+je-1 ) )+abs( work( 3*n+je-1 ) ) )
*
*              Compute contribution from columns JE and JE-1
*              of A and B to the sums.
*
               creala = acoef*work( 2*n+je-1 )
               cimaga = acoef*work( 3*n+je-1 )
               crealb = bcoefr*work( 2*n+je-1 ) -
     $                  bcoefi*work( 3*n+je-1 )
               cimagb = bcoefi*work( 2*n+je-1 ) +
     $                  bcoefr*work( 3*n+je-1 )
               cre2a = acoef*work( 2*n+je )
               cim2a = acoef*work( 3*n+je )
               cre2b = bcoefr*work( 2*n+je ) - bcoefi*work( 3*n+je )
               cim2b = bcoefi*work( 2*n+je ) + bcoefr*work( 3*n+je )
               DO 270 jr = 1, je - 2
                  work( 2*n+jr ) = -creala*s( jr, je-1 ) +
     $                             crealb*p( jr, je-1 ) -
     $                             cre2a*s( jr, je ) + cre2b*p( jr, je )
                  work( 3*n+jr ) = -cimaga*s( jr, je-1 ) +
     $                             cimagb*p( jr, je-1 ) -
     $                             cim2a*s( jr, je ) + cim2b*p( jr, je )
  270          CONTINUE
            END IF
*
            dmin = max( ulp*acoefa*anorm, ulp*bcoefa*bnorm, safmin )
*
*           Columnwise triangular solve of  (a A - b B)  x = 0
*
            il2by2 = .false.
            DO 370 j = je - nw, 1, -1
*
*              If a 2-by-2 block, is in position j-1:j, wait until
*              next iteration to process it (when it will be j:j+1)
*
               IF( .NOT.il2by2 .AND. j.GT.1 ) THEN
                  IF( s( j, j-1 ).NE.zero ) THEN
                     il2by2 = .true.
                     GO TO 370
                  END IF
               END IF
               bdiag( 1 ) = p( j, j )
               IF( il2by2 ) THEN
                  na = 2
                  bdiag( 2 ) = p( j+1, j+1 )
               ELSE
                  na = 1
               END IF
*
*              Compute x(j) (and x(j+1), if 2-by-2 block)
*
               CALL dlaln2( .false., na, nw, dmin, acoef, s( j, j ),
     $                      lds, bdiag( 1 ), bdiag( 2 ), work( 2*n+j ),
     $                      n, bcoefr, bcoefi, sum, 2, scale, temp,
     $                      iinfo )
               IF( scale.LT.one ) THEN
*
                  DO 290 jw = 0, nw - 1
                     DO 280 jr = 1, je
                        work( ( jw+2 )*n+jr ) = scale*
     $                     work( ( jw+2 )*n+jr )
  280                CONTINUE
  290             CONTINUE
               END IF
               xmax = max( scale*xmax, temp )
*
               DO 310 jw = 1, nw
                  DO 300 ja = 1, na
                     work( ( jw+1 )*n+j+ja-1 ) = sum( ja, jw )
  300             CONTINUE
  310          CONTINUE
*
*              w = w + x(j)*(a S(*,j) - b P(*,j) ) with scaling
*
               IF( j.GT.1 ) THEN
*
*                 Check whether scaling is necessary for sum.
*
                  xscale = one / max( one, xmax )
                  temp = acoefa*work( j ) + bcoefa*work( n+j )
                  IF( il2by2 )
     $               temp = max( temp, acoefa*work( j+1 )+bcoefa*
     $                      work( n+j+1 ) )
                  temp = max( temp, acoefa, bcoefa )
                  IF( temp.GT.bignum*xscale ) THEN
*
                     DO 330 jw = 0, nw - 1
                        DO 320 jr = 1, je
                           work( ( jw+2 )*n+jr ) = xscale*
     $                        work( ( jw+2 )*n+jr )
  320                   CONTINUE
  330                CONTINUE
                     xmax = xmax*xscale
                  END IF
*
*                 Compute the contributions of the off-diagonals of
*                 column j (and j+1, if 2-by-2 block) of A and B to the
*                 sums.
*
*
                  DO 360 ja = 1, na
                     IF( ilcplx ) THEN
                        creala = acoef*work( 2*n+j+ja-1 )
                        cimaga = acoef*work( 3*n+j+ja-1 )
                        crealb = bcoefr*work( 2*n+j+ja-1 ) -
     $                           bcoefi*work( 3*n+j+ja-1 )
                        cimagb = bcoefi*work( 2*n+j+ja-1 ) +
     $                           bcoefr*work( 3*n+j+ja-1 )
                        DO 340 jr = 1, j - 1
                           work( 2*n+jr ) = work( 2*n+jr ) -
     $                                      creala*s( jr, j+ja-1 ) +
     $                                      crealb*p( jr, j+ja-1 )
                           work( 3*n+jr ) = work( 3*n+jr ) -
     $                                      cimaga*s( jr, j+ja-1 ) +
     $                                      cimagb*p( jr, j+ja-1 )
  340                   CONTINUE
                     ELSE
                        creala = acoef*work( 2*n+j+ja-1 )
                        crealb = bcoefr*work( 2*n+j+ja-1 )
                        DO 350 jr = 1, j - 1
                           work( 2*n+jr ) = work( 2*n+jr ) -
     $                                      creala*s( jr, j+ja-1 ) +
     $                                      crealb*p( jr, j+ja-1 )
  350                   CONTINUE
                     END IF
  360             CONTINUE
               END IF
*
               il2by2 = .false.
  370       CONTINUE
*
*           Copy eigenvector to VR, back transforming if
*           HOWMNY='B'.
*
            ieig = ieig - nw
            IF( ilback ) THEN
*
               DO 410 jw = 0, nw - 1
                  DO 380 jr = 1, n
                     work( ( jw+4 )*n+jr ) = work( ( jw+2 )*n+1 )*
     $                                       vr( jr, 1 )
  380             CONTINUE
*
*                 A series of compiler directives to defeat
*                 vectorization for the next loop
*
*
                  DO 400 jc = 2, je
                     DO 390 jr = 1, n
                        work( ( jw+4 )*n+jr ) = work( ( jw+4 )*n+jr ) +
     $                     work( ( jw+2 )*n+jc )*vr( jr, jc )
  390                CONTINUE
  400             CONTINUE
  410          CONTINUE
*
               DO 430 jw = 0, nw - 1
                  DO 420 jr = 1, n
                     vr( jr, ieig+jw ) = work( ( jw+4 )*n+jr )
  420             CONTINUE
  430          CONTINUE
*
               iend = n
            ELSE
               DO 450 jw = 0, nw - 1
                  DO 440 jr = 1, n
                     vr( jr, ieig+jw ) = work( ( jw+2 )*n+jr )
  440             CONTINUE
  450          CONTINUE
*
               iend = je
            END IF
*
*           Scale eigenvector
*
            xmax = zero
            IF( ilcplx ) THEN
               DO 460 j = 1, iend
                  xmax = max( xmax, abs( vr( j, ieig ) )+
     $                   abs( vr( j, ieig+1 ) ) )
  460          CONTINUE
            ELSE
               DO 470 j = 1, iend
                  xmax = max( xmax, abs( vr( j, ieig ) ) )
  470          CONTINUE
            END IF
*
            IF( xmax.GT.safmin ) THEN
               xscale = one / xmax
               DO 490 jw = 0, nw - 1
                  DO 480 jr = 1, iend
                     vr( jr, ieig+jw ) = xscale*vr( jr, ieig+jw )
  480             CONTINUE
  490          CONTINUE
            END IF
  500    CONTINUE
      END IF
*
      RETURN
*
*     End of DTGEVC
*
      END

! DGGBAK
      SUBROUTINE dggbak( JOB, SIDE, N, ILO, IHI, LSCALE, RSCALE, M, V,
     $                   LDV, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          JOB, SIDE
      INTEGER            IHI, ILO, INFO, LDV, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   LSCALE( * ), RSCALE( * ), V( LDV, * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            LEFTV, RIGHTV
      INTEGER            I, K
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           lsame
*     ..
*     .. External Subroutines ..
      EXTERNAL           dscal, dswap, xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          max, int
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      rightv = lsame( side, 'R' )
      leftv = lsame( side, 'L' )
*
      info = 0
      IF( .NOT.lsame( job, 'N' ) .AND. .NOT.lsame( job, 'P' ) .AND.
     $    .NOT.lsame( job, 'S' ) .AND. .NOT.lsame( job, 'B' ) ) THEN
         info = -1
      ELSE IF( .NOT.rightv .AND. .NOT.leftv ) THEN
         info = -2
      ELSE IF( n.LT.0 ) THEN
         info = -3
      ELSE IF( ilo.LT.1 ) THEN
         info = -4
      ELSE IF( n.EQ.0 .AND. ihi.EQ.0 .AND. ilo.NE.1 ) THEN
         info = -4
      ELSE IF( n.GT.0 .AND. ( ihi.LT.ilo .OR. ihi.GT.max( 1, n ) ) )
     $   THEN
         info = -5
      ELSE IF( n.EQ.0 .AND. ilo.EQ.1 .AND. ihi.NE.0 ) THEN
         info = -5
      ELSE IF( m.LT.0 ) THEN
         info = -8
      ELSE IF( ldv.LT.max( 1, n ) ) THEN
         info = -10
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'DGGBAK', -info )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( n.EQ.0 )
     $   RETURN
      IF( m.EQ.0 )
     $   RETURN
      IF( lsame( job, 'N' ) )
     $   RETURN
*
      IF( ilo.EQ.ihi )
     $   GO TO 30
*
*     Backward balance
*
      IF( lsame( job, 'S' ) .OR. lsame( job, 'B' ) ) THEN
*
*        Backward transformation on right eigenvectors
*
         IF( rightv ) THEN
            DO 10 i = ilo, ihi
               CALL dscal( m, rscale( i ), v( i, 1 ), ldv )
   10       CONTINUE
         END IF
*
*        Backward transformation on left eigenvectors
*
         IF( leftv ) THEN
            DO 20 i = ilo, ihi
               CALL dscal( m, lscale( i ), v( i, 1 ), ldv )
   20       CONTINUE
         END IF
      END IF
*
*     Backward permutation
*
   30 CONTINUE
      IF( lsame( job, 'P' ) .OR. lsame( job, 'B' ) ) THEN
*
*        Backward permutation on right eigenvectors
*
         IF( rightv ) THEN
            IF( ilo.EQ.1 )
     $         GO TO 50
*
            DO 40 i = ilo - 1, 1, -1
               k = int(rscale( i ))
               IF( k.EQ.i )
     $            GO TO 40
               CALL dswap( m, v( i, 1 ), ldv, v( k, 1 ), ldv )
   40       CONTINUE
*
   50       CONTINUE
            IF( ihi.EQ.n )
     $         GO TO 70
            DO 60 i = ihi + 1, n
               k = int(rscale( i ))
               IF( k.EQ.i )
     $            GO TO 60
               CALL dswap( m, v( i, 1 ), ldv, v( k, 1 ), ldv )
   60       CONTINUE
         END IF
*
*        Backward permutation on left eigenvectors
*
   70    CONTINUE
         IF( leftv ) THEN
            IF( ilo.EQ.1 )
     $         GO TO 90
            DO 80 i = ilo - 1, 1, -1
               k = int(lscale( i ))
               IF( k.EQ.i )
     $            GO TO 80
               CALL dswap( m, v( i, 1 ), ldv, v( k, 1 ), ldv )
   80       CONTINUE
*
   90       CONTINUE
            IF( ihi.EQ.n )
     $         GO TO 110
            DO 100 i = ihi + 1, n
               k = int(lscale( i ))
               IF( k.EQ.i )
     $            GO TO 100
               CALL dswap( m, v( i, 1 ), ldv, v( k, 1 ), ldv )
  100       CONTINUE
         END IF
      END IF
*
  110 CONTINUE
*
      RETURN
*
*     End of DGGBAK
*
      END

! DLACPY
      SUBROUTINE dlacpy( UPLO, M, N, A, LDA, B, LDB )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, LDB, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, J
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           lsame
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          min
*     ..
*     .. Executable Statements ..
*
      IF( lsame( uplo, 'U' ) ) THEN
         DO 20 j = 1, n
            DO 10 i = 1, min( j, m )
               b( i, j ) = a( i, j )
   10       CONTINUE
   20    CONTINUE
      ELSE IF( lsame( uplo, 'L' ) ) THEN
         DO 40 j = 1, n
            DO 30 i = j, m
               b( i, j ) = a( i, j )
   30       CONTINUE
   40    CONTINUE
      ELSE
         DO 60 j = 1, n
            DO 50 i = 1, m
               b( i, j ) = a( i, j )
   50       CONTINUE
   60    CONTINUE
      END IF
      RETURN
*
*     End of DLACPY
*
      END

! DGEQR2
      SUBROUTINE dgeqr2( M, N, A, LDA, TAU, WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      parameter( one = 1.0d+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, K
      DOUBLE PRECISION   AII
*     ..
*     .. External Subroutines ..
      EXTERNAL           dlarf, dlarfg, xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          max, min
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      info = 0
      IF( m.LT.0 ) THEN
         info = -1
      ELSE IF( n.LT.0 ) THEN
         info = -2
      ELSE IF( lda.LT.max( 1, m ) ) THEN
         info = -4
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'DGEQR2', -info )
         RETURN
      END IF
*
      k = min( m, n )
*
      DO 10 i = 1, k
*
*        Generate elementary reflector H(i) to annihilate A(i+1:m,i)
*
         CALL dlarfg( m-i+1, a( i, i ), a( min( i+1, m ), i ), 1,
     $                tau( i ) )
         IF( i.LT.n ) THEN
*
*           Apply H(i) to A(i:m,i+1:n) from the left
*
            aii = a( i, i )
            a( i, i ) = one
            CALL dlarf( 'Left', m-i+1, n-i, a( i, i ), 1, tau( i ),
     $                  a( i, i+1 ), lda, work )
            a( i, i ) = aii
         END IF
   10 CONTINUE
      RETURN
*
*     End of DGEQR2
*
      END

! DORM2R
      SUBROUTINE dorm2r( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
     $                   WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      parameter( one = 1.0d+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LEFT, NOTRAN
      INTEGER            I, I1, I2, I3, IC, JC, MI, NI, NQ
      DOUBLE PRECISION   AII
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           lsame
*     ..
*     .. External Subroutines ..
      EXTERNAL           dlarf, xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          max
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      info = 0
      left = lsame( side, 'L' )
      notran = lsame( trans, 'N' )
*
*     NQ is the order of Q
*
      IF( left ) THEN
         nq = m
      ELSE
         nq = n
      END IF
      IF( .NOT.left .AND. .NOT.lsame( side, 'R' ) ) THEN
         info = -1
      ELSE IF( .NOT.notran .AND. .NOT.lsame( trans, 'T' ) ) THEN
         info = -2
      ELSE IF( m.LT.0 ) THEN
         info = -3
      ELSE IF( n.LT.0 ) THEN
         info = -4
      ELSE IF( k.LT.0 .OR. k.GT.nq ) THEN
         info = -5
      ELSE IF( lda.LT.max( 1, nq ) ) THEN
         info = -7
      ELSE IF( ldc.LT.max( 1, m ) ) THEN
         info = -10
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'DORM2R', -info )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( m.EQ.0 .OR. n.EQ.0 .OR. k.EQ.0 )
     $   RETURN
*
      IF( ( left .AND. .NOT.notran ) .OR. ( .NOT.left .AND. notran ) )
     $     THEN
         i1 = 1
         i2 = k
         i3 = 1
      ELSE
         i1 = k
         i2 = 1
         i3 = -1
      END IF
*
      IF( left ) THEN
         ni = n
         jc = 1
      ELSE
         mi = m
         ic = 1
      END IF
*
      DO 10 i = i1, i2, i3
         IF( left ) THEN
*
*           H(i) is applied to C(i:m,1:n)
*
            mi = m - i + 1
            ic = i
         ELSE
*
*           H(i) is applied to C(1:m,i:n)
*
            ni = n - i + 1
            jc = i
         END IF
*
*        Apply H(i)
*
         aii = a( i, i )
         a( i, i ) = one
         CALL dlarf( side, mi, ni, a( i, i ), 1, tau( i ), c( ic, jc ),
     $               ldc, work )
         a( i, i ) = aii
   10 CONTINUE
      RETURN
*
*     End of DORM2R
*
      END

! DLANHS
      DOUBLE PRECISION FUNCTION dlanhs( NORM, N, A, LDA, WORK )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          norm
      INTEGER            lda, n
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   a( lda, * ), work( * )
*     ..
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   one, zero
      parameter( one = 1.0d+0, zero = 0.0d+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            i, j
      DOUBLE PRECISION   scale, sum, value
*     ..
*     .. External Subroutines ..
      EXTERNAL           dlassq
*     ..
*     .. External Functions ..
      LOGICAL            lsame, disnan
      EXTERNAL           lsame, disnan
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, min, sqrt
*     ..
*     .. Executable Statements ..
*
      IF( n.EQ.0 ) THEN
         VALUE = zero
      ELSE IF( lsame( norm, 'M' ) ) THEN
*
*        Find max(abs(A(i,j))).
*
         VALUE = zero
         DO 20 j = 1, n
            DO 10 i = 1, min( n, j+1 )
               sum = abs( a( i, j ) )
               IF( VALUE .LT. sum .OR. disnan( sum ) ) VALUE = sum
   10       CONTINUE
   20    CONTINUE
      ELSE IF( ( lsame( norm, 'O' ) ) .OR. ( norm.EQ.'1' ) ) THEN
*
*        Find norm1(A).
*
         VALUE = zero
         DO 40 j = 1, n
            sum = zero
            DO 30 i = 1, min( n, j+1 )
               sum = sum + abs( a( i, j ) )
   30       CONTINUE
            IF( VALUE .LT. sum .OR. disnan( sum ) ) VALUE = sum
   40    CONTINUE
      ELSE IF( lsame( norm, 'I' ) ) THEN
*
*        Find normI(A).
*
         DO 50 i = 1, n
            work( i ) = zero
   50    CONTINUE
         DO 70 j = 1, n
            DO 60 i = 1, min( n, j+1 )
               work( i ) = work( i ) + abs( a( i, j ) )
   60       CONTINUE
   70    CONTINUE
         VALUE = zero
         DO 80 i = 1, n
            sum = work( i )
            IF( VALUE .LT. sum .OR. disnan( sum ) ) VALUE = sum
   80    CONTINUE
      ELSE IF( ( lsame( norm, 'F' ) ) .OR. ( lsame( norm, 'E' ) ) ) THEN
*
*        Find normF(A).
*
         scale = zero
         sum = one
         DO 90 j = 1, n
            CALL dlassq( min( n, j+1 ), a( 1, j ), 1, scale, sum )
   90    CONTINUE
         VALUE = scale*sqrt( sum )
      END IF
*
      dlanhs = VALUE
      RETURN
*
*     End of DLANHS
*
      END

! DLAG2
      SUBROUTINE dlag2( A, LDA, B, LDB, SAFMIN, SCALE1, SCALE2, WR1,
     $                  WR2, WI )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            LDA, LDB
      DOUBLE PRECISION   SAFMIN, SCALE1, SCALE2, WI, WR1, WR2
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO
      parameter( zero = 0.0d+0, one = 1.0d+0, two = 2.0d+0 )
      DOUBLE PRECISION   HALF
      parameter( half = one / two )
      DOUBLE PRECISION   FUZZY1
      parameter( fuzzy1 = one+1.0d-5 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   A11, A12, A21, A22, ABI22, ANORM, AS11, AS12,
     $                   as22, ascale, b11, b12, b22, binv11, binv22,
     $                   bmin, bnorm, bscale, bsize, c1, c2, c3, c4, c5,
     $                   diff, discr, pp, qq, r, rtmax, rtmin, s1, s2,
     $                   safmax, shift, ss, sum, wabs, wbig, wdet,
     $                   wscale, wsize, wsmall
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, max, min, sign, sqrt
*     ..
*     .. Executable Statements ..
*
      rtmin = sqrt( safmin )
      rtmax = one / rtmin
      safmax = one / safmin
*
*     Scale A
*
      anorm = max( abs( a( 1, 1 ) )+abs( a( 2, 1 ) ),
     $        abs( a( 1, 2 ) )+abs( a( 2, 2 ) ), safmin )
      ascale = one / anorm
      a11 = ascale*a( 1, 1 )
      a21 = ascale*a( 2, 1 )
      a12 = ascale*a( 1, 2 )
      a22 = ascale*a( 2, 2 )
*
*     Perturb B if necessary to insure non-singularity
*
      b11 = b( 1, 1 )
      b12 = b( 1, 2 )
      b22 = b( 2, 2 )
      bmin = rtmin*max( abs( b11 ), abs( b12 ), abs( b22 ), rtmin )
      IF( abs( b11 ).LT.bmin )
     $   b11 = sign( bmin, b11 )
      IF( abs( b22 ).LT.bmin )
     $   b22 = sign( bmin, b22 )
*
*     Scale B
*
      bnorm = max( abs( b11 ), abs( b12 )+abs( b22 ), safmin )
      bsize = max( abs( b11 ), abs( b22 ) )
      bscale = one / bsize
      b11 = b11*bscale
      b12 = b12*bscale
      b22 = b22*bscale
*
*     Compute larger eigenvalue by method described by C. van Loan
*
*     ( AS is A shifted by -SHIFT*B )
*
      binv11 = one / b11
      binv22 = one / b22
      s1 = a11*binv11
      s2 = a22*binv22
      IF( abs( s1 ).LE.abs( s2 ) ) THEN
         as12 = a12 - s1*b12
         as22 = a22 - s1*b22
         ss = a21*( binv11*binv22 )
         abi22 = as22*binv22 - ss*b12
         pp = half*abi22
         shift = s1
      ELSE
         as12 = a12 - s2*b12
         as11 = a11 - s2*b11
         ss = a21*( binv11*binv22 )
         abi22 = -ss*b12
         pp = half*( as11*binv11+abi22 )
         shift = s2
      END IF
      qq = ss*as12
      IF( abs( pp*rtmin ).GE.one ) THEN
         discr = ( rtmin*pp )**2 + qq*safmin
         r = sqrt( abs( discr ) )*rtmax
      ELSE
         IF( pp**2+abs( qq ).LE.safmin ) THEN
            discr = ( rtmax*pp )**2 + qq*safmax
            r = sqrt( abs( discr ) )*rtmin
         ELSE
            discr = pp**2 + qq
            r = sqrt( abs( discr ) )
         END IF
      END IF
*
*     Note: the test of R in the following IF is to cover the case when
*           DISCR is small and negative and is flushed to zero during
*           the calculation of R.  On machines which have a consistent
*           flush-to-zero threshold and handle numbers above that
*           threshold correctly, it would not be necessary.
*
      IF( discr.GE.zero .OR. r.EQ.zero ) THEN
         sum = pp + sign( r, pp )
         diff = pp - sign( r, pp )
         wbig = shift + sum
*
*        Compute smaller eigenvalue
*
         wsmall = shift + diff
         IF( half*abs( wbig ).GT.max( abs( wsmall ), safmin ) ) THEN
            wdet = ( a11*a22-a12*a21 )*( binv11*binv22 )
            wsmall = wdet / wbig
         END IF
*
*        Choose (real) eigenvalue closest to 2,2 element of A*B**(-1)
*        for WR1.
*
         IF( pp.GT.abi22 ) THEN
            wr1 = min( wbig, wsmall )
            wr2 = max( wbig, wsmall )
         ELSE
            wr1 = max( wbig, wsmall )
            wr2 = min( wbig, wsmall )
         END IF
         wi = zero
      ELSE
*
*        Complex eigenvalues
*
         wr1 = shift + pp
         wr2 = wr1
         wi = r
      END IF
*
*     Further scaling to avoid underflow and overflow in computing
*     SCALE1 and overflow in computing w*B.
*
*     This scale factor (WSCALE) is bounded from above using C1 and C2,
*     and from below using C3 and C4.
*        C1 implements the condition  s A  must never overflow.
*        C2 implements the condition  w B  must never overflow.
*        C3, with C2,
*           implement the condition that s A - w B must never overflow.
*        C4 implements the condition  s    should not underflow.
*        C5 implements the condition  max(s,|w|) should be at least 2.
*
      c1 = bsize*( safmin*max( one, ascale ) )
      c2 = safmin*max( one, bnorm )
      c3 = bsize*safmin
      IF( ascale.LE.one .AND. bsize.LE.one ) THEN
         c4 = min( one, ( ascale / safmin )*bsize )
      ELSE
         c4 = one
      END IF
      IF( ascale.LE.one .OR. bsize.LE.one ) THEN
         c5 = min( one, ascale*bsize )
      ELSE
         c5 = one
      END IF
*
*     Scale first eigenvalue
*
      wabs = abs( wr1 ) + abs( wi )
      wsize = max( safmin, c1, fuzzy1*( wabs*c2+c3 ),
     $        min( c4, half*max( wabs, c5 ) ) )
      IF( wsize.NE.one ) THEN
         wscale = one / wsize
         IF( wsize.GT.one ) THEN
            scale1 = ( max( ascale, bsize )*wscale )*
     $               min( ascale, bsize )
         ELSE
            scale1 = ( min( ascale, bsize )*wscale )*
     $               max( ascale, bsize )
         END IF
         wr1 = wr1*wscale
         IF( wi.NE.zero ) THEN
            wi = wi*wscale
            wr2 = wr1
            scale2 = scale1
         END IF
      ELSE
         scale1 = ascale*bsize
         scale2 = scale1
      END IF
*
*     Scale second eigenvalue (if real)
*
      IF( wi.EQ.zero ) THEN
         wsize = max( safmin, c1, fuzzy1*( abs( wr2 )*c2+c3 ),
     $           min( c4, half*max( abs( wr2 ), c5 ) ) )
         IF( wsize.NE.one ) THEN
            wscale = one / wsize
            IF( wsize.GT.one ) THEN
               scale2 = ( max( ascale, bsize )*wscale )*
     $                  min( ascale, bsize )
            ELSE
               scale2 = ( min( ascale, bsize )*wscale )*
     $                  max( ascale, bsize )
            END IF
            wr2 = wr2*wscale
         ELSE
            scale2 = ascale*bsize
         END IF
      END IF
*
*     End of DLAG2
*
      RETURN
      END

! DLASV2
      SUBROUTINE dlasv2( F, G, H, SSMIN, SSMAX, SNR, CSR, SNL, CSL )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   CSL, CSR, F, G, H, SNL, SNR, SSMAX, SSMIN
*     ..
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      parameter( zero = 0.0d0 )
      DOUBLE PRECISION   HALF
      parameter( half = 0.5d0 )
      DOUBLE PRECISION   ONE
      parameter( one = 1.0d0 )
      DOUBLE PRECISION   TWO
      parameter( two = 2.0d0 )
      DOUBLE PRECISION   FOUR
      parameter( four = 4.0d0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            GASMAL, SWAP
      INTEGER            PMAX
      DOUBLE PRECISION   A, CLT, CRT, D, FA, FT, GA, GT, HA, HT, L, M,
     $                   MM, R, S, SLT, SRT, T, TEMP, TSIGN, TT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, sign, sqrt
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           dlamch
*     ..
*     .. Executable Statements ..
*
      ft = f
      fa = abs( ft )
      ht = h
      ha = abs( h )
*
*     PMAX points to the maximum absolute element of matrix
*       PMAX = 1 if F largest in absolute values
*       PMAX = 2 if G largest in absolute values
*       PMAX = 3 if H largest in absolute values
*
      pmax = 1
      swap = ( ha.GT.fa )
      IF( swap ) THEN
         pmax = 3
         temp = ft
         ft = ht
         ht = temp
         temp = fa
         fa = ha
         ha = temp
*
*        Now FA .ge. HA
*
      END IF
      gt = g
      ga = abs( gt )
      IF( ga.EQ.zero ) THEN
*
*        Diagonal matrix
*
         ssmin = ha
         ssmax = fa
         clt = one
         crt = one
         slt = zero
         srt = zero
      ELSE
         gasmal = .true.
         IF( ga.GT.fa ) THEN
            pmax = 2
            IF( ( fa / ga ).LT.dlamch( 'EPS' ) ) THEN
*
*              Case of very large GA
*
               gasmal = .false.
               ssmax = ga
               IF( ha.GT.one ) THEN
                  ssmin = fa / ( ga / ha )
               ELSE
                  ssmin = ( fa / ga )*ha
               END IF
               clt = one
               slt = ht / gt
               srt = one
               crt = ft / gt
            END IF
         END IF
         IF( gasmal ) THEN
*
*           Normal case
*
            d = fa - ha
            IF( d.EQ.fa ) THEN
*
*              Copes with infinite F or H
*
               l = one
            ELSE
               l = d / fa
            END IF
*
*           Note that 0 .le. L .le. 1
*
            m = gt / ft
*
*           Note that abs(M) .le. 1/macheps
*
            t = two - l
*
*           Note that T .ge. 1
*
            mm = m*m
            tt = t*t
            s = sqrt( tt+mm )
*
*           Note that 1 .le. S .le. 1 + 1/macheps
*
            IF( l.EQ.zero ) THEN
               r = abs( m )
            ELSE
               r = sqrt( l*l+mm )
            END IF
*
*           Note that 0 .le. R .le. 1 + 1/macheps
*
            a = half*( s+r )
*
*           Note that 1 .le. A .le. 1 + abs(M)
*
            ssmin = ha / a
            ssmax = fa*a
            IF( mm.EQ.zero ) THEN
*
*              Note that M is very tiny
*
               IF( l.EQ.zero ) THEN
                  t = sign( two, ft )*sign( one, gt )
               ELSE
                  t = gt / sign( d, ft ) + m / t
               END IF
            ELSE
               t = ( m / ( s+t )+m / ( r+l ) )*( one+a )
            END IF
            l = sqrt( t*t+four )
            crt = two / l
            srt = t / l
            clt = ( crt+srt*m ) / a
            slt = ( ht / ft )*srt / a
         END IF
      END IF
      IF( swap ) THEN
         csl = srt
         snl = crt
         csr = slt
         snr = clt
      ELSE
         csl = clt
         snl = slt
         csr = crt
         snr = srt
      END IF
*
*     Correct signs of SSMAX and SSMIN
*
      IF( pmax.EQ.1 )
     $   tsign = sign( one, csr )*sign( one, csl )*sign( one, f )
      IF( pmax.EQ.2 )
     $   tsign = sign( one, snr )*sign( one, csl )*sign( one, g )
      IF( pmax.EQ.3 )
     $   tsign = sign( one, snr )*sign( one, snl )*sign( one, h )
      ssmax = sign( ssmax, tsign )
      ssmin = sign( ssmin, tsign*sign( one, f )*sign( one, h ) )
      RETURN
*
*     End of DLASV2
*
      END

! DLALN2
      SUBROUTINE dlaln2( LTRANS, NA, NW, SMIN, CA, A, LDA, D1, D2, B,
     $                   LDB, WR, WI, X, LDX, SCALE, XNORM, INFO )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      LOGICAL            LTRANS
      INTEGER            INFO, LDA, LDB, LDX, NA, NW
      DOUBLE PRECISION   CA, D1, D2, SCALE, SMIN, WI, WR, XNORM
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), X( LDX, * )
*     ..
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      parameter( zero = 0.0d0, one = 1.0d0 )
      DOUBLE PRECISION   TWO
      parameter( two = 2.0d0 )
*     ..
*     .. Local Scalars ..
      INTEGER            ICMAX, J
      DOUBLE PRECISION   BBND, BI1, BI2, BIGNUM, BNORM, BR1, BR2, CI21,
     $                   ci22, cmax, cnorm, cr21, cr22, csi, csr, li21,
     $                   lr21, smini, smlnum, temp, u22abs, ui11, ui11r,
     $                   ui12, ui12s, ui22, ur11, ur11r, ur12, ur12s,
     $                   ur22, xi1, xi2, xr1, xr2
*     ..
*     .. Local Arrays ..
      LOGICAL            RSWAP( 4 ), ZSWAP( 4 )
      INTEGER            IPIVOT( 4, 4 )
      DOUBLE PRECISION   CI( 2, 2 ), CIV( 4 ), CR( 2, 2 ), CRV( 4 )
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           dlamch
*     ..
*     .. External Subroutines ..
      EXTERNAL           dladiv
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, max
*     ..
*     .. Equivalences ..
      equivalence( ci( 1, 1 ), civ( 1 ) ),
     $                   ( cr( 1, 1 ), crv( 1 ) )
*     ..
*     .. Data statements ..
      DATA               zswap / .false., .false., .true., .true. /
      DATA               rswap / .false., .true., .false., .true. /
      DATA               ipivot / 1, 2, 3, 4, 2, 1, 4, 3, 3, 4, 1, 2, 4,
     $                   3, 2, 1 /
*     ..
*     .. Executable Statements ..
*
*     Compute BIGNUM
*
      smlnum = two*dlamch( 'Safe minimum' )
      bignum = one / smlnum
      smini = max( smin, smlnum )
*
*     Don't check for input errors
*
      info = 0
*
*     Standard Initializations
*
      scale = one
*
      IF( na.EQ.1 ) THEN
*
*        1 x 1  (i.e., scalar) system   C X = B
*
         IF( nw.EQ.1 ) THEN
*
*           Real 1x1 system.
*
*           C = ca A - w D
*
            csr = ca*a( 1, 1 ) - wr*d1
            cnorm = abs( csr )
*
*           If | C | < SMINI, use C = SMINI
*
            IF( cnorm.LT.smini ) THEN
               csr = smini
               cnorm = smini
               info = 1
            END IF
*
*           Check scaling for  X = B / C
*
            bnorm = abs( b( 1, 1 ) )
            IF( cnorm.LT.one .AND. bnorm.GT.one ) THEN
               IF( bnorm.GT.bignum*cnorm )
     $            scale = one / bnorm
            END IF
*
*           Compute X
*
            x( 1, 1 ) = ( b( 1, 1 )*scale ) / csr
            xnorm = abs( x( 1, 1 ) )
         ELSE
*
*           Complex 1x1 system (w is complex)
*
*           C = ca A - w D
*
            csr = ca*a( 1, 1 ) - wr*d1
            csi = -wi*d1
            cnorm = abs( csr ) + abs( csi )
*
*           If | C | < SMINI, use C = SMINI
*
            IF( cnorm.LT.smini ) THEN
               csr = smini
               csi = zero
               cnorm = smini
               info = 1
            END IF
*
*           Check scaling for  X = B / C
*
            bnorm = abs( b( 1, 1 ) ) + abs( b( 1, 2 ) )
            IF( cnorm.LT.one .AND. bnorm.GT.one ) THEN
               IF( bnorm.GT.bignum*cnorm )
     $            scale = one / bnorm
            END IF
*
*           Compute X
*
            CALL dladiv( scale*b( 1, 1 ), scale*b( 1, 2 ), csr, csi,
     $                   x( 1, 1 ), x( 1, 2 ) )
            xnorm = abs( x( 1, 1 ) ) + abs( x( 1, 2 ) )
         END IF
*
      ELSE
*
*        2x2 System
*
*        Compute the real part of  C = ca A - w D  (or  ca A**T - w D )
*
         cr( 1, 1 ) = ca*a( 1, 1 ) - wr*d1
         cr( 2, 2 ) = ca*a( 2, 2 ) - wr*d2
         IF( ltrans ) THEN
            cr( 1, 2 ) = ca*a( 2, 1 )
            cr( 2, 1 ) = ca*a( 1, 2 )
         ELSE
            cr( 2, 1 ) = ca*a( 2, 1 )
            cr( 1, 2 ) = ca*a( 1, 2 )
         END IF
*
         IF( nw.EQ.1 ) THEN
*
*           Real 2x2 system  (w is real)
*
*           Find the largest element in C
*
            cmax = zero
            icmax = 0
*
            DO 10 j = 1, 4
               IF( abs( crv( j ) ).GT.cmax ) THEN
                  cmax = abs( crv( j ) )
                  icmax = j
               END IF
   10       CONTINUE
*
*           If norm(C) < SMINI, use SMINI*identity.
*
            IF( cmax.LT.smini ) THEN
               bnorm = max( abs( b( 1, 1 ) ), abs( b( 2, 1 ) ) )
               IF( smini.LT.one .AND. bnorm.GT.one ) THEN
                  IF( bnorm.GT.bignum*smini )
     $               scale = one / bnorm
               END IF
               temp = scale / smini
               x( 1, 1 ) = temp*b( 1, 1 )
               x( 2, 1 ) = temp*b( 2, 1 )
               xnorm = temp*bnorm
               info = 1
               RETURN
            END IF
*
*           Gaussian elimination with complete pivoting.
*
            ur11 = crv( icmax )
            cr21 = crv( ipivot( 2, icmax ) )
            ur12 = crv( ipivot( 3, icmax ) )
            cr22 = crv( ipivot( 4, icmax ) )
            ur11r = one / ur11
            lr21 = ur11r*cr21
            ur22 = cr22 - ur12*lr21
*
*           If smaller pivot < SMINI, use SMINI
*
            IF( abs( ur22 ).LT.smini ) THEN
               ur22 = smini
               info = 1
            END IF
            IF( rswap( icmax ) ) THEN
               br1 = b( 2, 1 )
               br2 = b( 1, 1 )
            ELSE
               br1 = b( 1, 1 )
               br2 = b( 2, 1 )
            END IF
            br2 = br2 - lr21*br1
            bbnd = max( abs( br1*( ur22*ur11r ) ), abs( br2 ) )
            IF( bbnd.GT.one .AND. abs( ur22 ).LT.one ) THEN
               IF( bbnd.GE.bignum*abs( ur22 ) )
     $            scale = one / bbnd
            END IF
*
            xr2 = ( br2*scale ) / ur22
            xr1 = ( scale*br1 )*ur11r - xr2*( ur11r*ur12 )
            IF( zswap( icmax ) ) THEN
               x( 1, 1 ) = xr2
               x( 2, 1 ) = xr1
            ELSE
               x( 1, 1 ) = xr1
               x( 2, 1 ) = xr2
            END IF
            xnorm = max( abs( xr1 ), abs( xr2 ) )
*
*           Further scaling if  norm(A) norm(X) > overflow
*
            IF( xnorm.GT.one .AND. cmax.GT.one ) THEN
               IF( xnorm.GT.bignum / cmax ) THEN
                  temp = cmax / bignum
                  x( 1, 1 ) = temp*x( 1, 1 )
                  x( 2, 1 ) = temp*x( 2, 1 )
                  xnorm = temp*xnorm
                  scale = temp*scale
               END IF
            END IF
         ELSE
*
*           Complex 2x2 system  (w is complex)
*
*           Find the largest element in C
*
            ci( 1, 1 ) = -wi*d1
            ci( 2, 1 ) = zero
            ci( 1, 2 ) = zero
            ci( 2, 2 ) = -wi*d2
            cmax = zero
            icmax = 0
*
            DO 20 j = 1, 4
               IF( abs( crv( j ) )+abs( civ( j ) ).GT.cmax ) THEN
                  cmax = abs( crv( j ) ) + abs( civ( j ) )
                  icmax = j
               END IF
   20       CONTINUE
*
*           If norm(C) < SMINI, use SMINI*identity.
*
            IF( cmax.LT.smini ) THEN
               bnorm = max( abs( b( 1, 1 ) )+abs( b( 1, 2 ) ),
     $                 abs( b( 2, 1 ) )+abs( b( 2, 2 ) ) )
               IF( smini.LT.one .AND. bnorm.GT.one ) THEN
                  IF( bnorm.GT.bignum*smini )
     $               scale = one / bnorm
               END IF
               temp = scale / smini
               x( 1, 1 ) = temp*b( 1, 1 )
               x( 2, 1 ) = temp*b( 2, 1 )
               x( 1, 2 ) = temp*b( 1, 2 )
               x( 2, 2 ) = temp*b( 2, 2 )
               xnorm = temp*bnorm
               info = 1
               RETURN
            END IF
*
*           Gaussian elimination with complete pivoting.
*
            ur11 = crv( icmax )
            ui11 = civ( icmax )
            cr21 = crv( ipivot( 2, icmax ) )
            ci21 = civ( ipivot( 2, icmax ) )
            ur12 = crv( ipivot( 3, icmax ) )
            ui12 = civ( ipivot( 3, icmax ) )
            cr22 = crv( ipivot( 4, icmax ) )
            ci22 = civ( ipivot( 4, icmax ) )
            IF( icmax.EQ.1 .OR. icmax.EQ.4 ) THEN
*
*              Code when off-diagonals of pivoted C are real
*
               IF( abs( ur11 ).GT.abs( ui11 ) ) THEN
                  temp = ui11 / ur11
                  ur11r = one / ( ur11*( one+temp**2 ) )
                  ui11r = -temp*ur11r
               ELSE
                  temp = ur11 / ui11
                  ui11r = -one / ( ui11*( one+temp**2 ) )
                  ur11r = -temp*ui11r
               END IF
               lr21 = cr21*ur11r
               li21 = cr21*ui11r
               ur12s = ur12*ur11r
               ui12s = ur12*ui11r
               ur22 = cr22 - ur12*lr21
               ui22 = ci22 - ur12*li21
            ELSE
*
*              Code when diagonals of pivoted C are real
*
               ur11r = one / ur11
               ui11r = zero
               lr21 = cr21*ur11r
               li21 = ci21*ur11r
               ur12s = ur12*ur11r
               ui12s = ui12*ur11r
               ur22 = cr22 - ur12*lr21 + ui12*li21
               ui22 = -ur12*li21 - ui12*lr21
            END IF
            u22abs = abs( ur22 ) + abs( ui22 )
*
*           If smaller pivot < SMINI, use SMINI
*
            IF( u22abs.LT.smini ) THEN
               ur22 = smini
               ui22 = zero
               info = 1
            END IF
            IF( rswap( icmax ) ) THEN
               br2 = b( 1, 1 )
               br1 = b( 2, 1 )
               bi2 = b( 1, 2 )
               bi1 = b( 2, 2 )
            ELSE
               br1 = b( 1, 1 )
               br2 = b( 2, 1 )
               bi1 = b( 1, 2 )
               bi2 = b( 2, 2 )
            END IF
            br2 = br2 - lr21*br1 + li21*bi1
            bi2 = bi2 - li21*br1 - lr21*bi1
            bbnd = max( ( abs( br1 )+abs( bi1 ) )*
     $             ( u22abs*( abs( ur11r )+abs( ui11r ) ) ),
     $             abs( br2 )+abs( bi2 ) )
            IF( bbnd.GT.one .AND. u22abs.LT.one ) THEN
               IF( bbnd.GE.bignum*u22abs ) THEN
                  scale = one / bbnd
                  br1 = scale*br1
                  bi1 = scale*bi1
                  br2 = scale*br2
                  bi2 = scale*bi2
               END IF
            END IF
*
            CALL dladiv( br2, bi2, ur22, ui22, xr2, xi2 )
            xr1 = ur11r*br1 - ui11r*bi1 - ur12s*xr2 + ui12s*xi2
            xi1 = ui11r*br1 + ur11r*bi1 - ui12s*xr2 - ur12s*xi2
            IF( zswap( icmax ) ) THEN
               x( 1, 1 ) = xr2
               x( 2, 1 ) = xr1
               x( 1, 2 ) = xi2
               x( 2, 2 ) = xi1
            ELSE
               x( 1, 1 ) = xr1
               x( 2, 1 ) = xr2
               x( 1, 2 ) = xi1
               x( 2, 2 ) = xi2
            END IF
            xnorm = max( abs( xr1 )+abs( xi1 ), abs( xr2 )+abs( xi2 ) )
*
*           Further scaling if  norm(A) norm(X) > overflow
*
            IF( xnorm.GT.one .AND. cmax.GT.one ) THEN
               IF( xnorm.GT.bignum / cmax ) THEN
                  temp = cmax / bignum
                  x( 1, 1 ) = temp*x( 1, 1 )
                  x( 2, 1 ) = temp*x( 2, 1 )
                  x( 1, 2 ) = temp*x( 1, 2 )
                  x( 2, 2 ) = temp*x( 2, 2 )
                  xnorm = temp*xnorm
                  scale = temp*scale
               END IF
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of DLALN2
*
      END

! ZGEEV      
      SUBROUTINE ZGEEV( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR,
     $                  WORK, LWORK, RWORK, INFO )
      implicit none
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          JOBVL, JOBVR
      INTEGER            INFO, LDA, LDVL, LDVR, LWORK, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ),
     $                   W( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, SCALEA, WANTVL, WANTVR
      CHARACTER          SIDE
      INTEGER            HSWORK, I, IBAL, IERR, IHI, ILO, IRWORK, ITAU,
     $                   IWRK, K, LWORK_TREVC, MAXWRK, MINWRK, NOUT
      DOUBLE PRECISION   ANRM, BIGNUM, CSCALE, EPS, SCL, SMLNUM
      COMPLEX*16         TMP
*     ..
*     .. Local Arrays ..
      LOGICAL            SELECT( 1 )
      DOUBLE PRECISION   DUM( 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLABAD, XERBLA, ZDSCAL, ZGEBAK, ZGEBAL, ZGEHRD,
     $                   ZHSEQR, ZLACPY, ZLASCL, ZSCAL, ZTREVC3, ZUNGHR
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IDAMAX, ILAENV
      DOUBLE PRECISION   DLAMCH, DZNRM2, ZLANGE
      EXTERNAL           LSAME, IDAMAX, ILAENV, DLAMCH, DZNRM2, ZLANGE
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCMPLX, CONJG, AIMAG, MAX, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      WANTVL = LSAME( JOBVL, 'V' )
      WANTVR = LSAME( JOBVR, 'V' )
      IF( ( .NOT.WANTVL ) .AND. ( .NOT.LSAME( JOBVL, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( ( .NOT.WANTVR ) .AND. ( .NOT.LSAME( JOBVR, 'N' ) ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDVL.LT.1 .OR. ( WANTVL .AND. LDVL.LT.N ) ) THEN
         INFO = -8
      ELSE IF( LDVR.LT.1 .OR. ( WANTVR .AND. LDVR.LT.N ) ) THEN
         INFO = -10
      END IF
*
*     Compute workspace
*      (Note: Comments in the code beginning "Workspace:" describe the
*       minimal amount of workspace needed at that point in the code,
*       as well as the preferred amount for good performance.
*       CWorkspace refers to complex workspace, and RWorkspace to real
*       workspace. NB refers to the optimal block size for the
*       immediately following subroutine, as returned by ILAENV.
*       HSWORK refers to the workspace preferred by ZHSEQR, as
*       calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
*       the worst case.)
*
      IF( INFO.EQ.0 ) THEN
         IF( N.EQ.0 ) THEN
            MINWRK = 1
            MAXWRK = 1
         ELSE
            MAXWRK = N + N*ILAENV( 1, 'ZGEHRD', ' ', N, 1, N, 0 )
            MINWRK = 2*N
            IF( WANTVL ) THEN
               MAXWRK = MAX( MAXWRK, N + ( N - 1 )*ILAENV( 1, 'ZUNGHR',
     $                       ' ', N, 1, N, -1 ) )
               CALL ZTREVC3( 'L', 'B', SELECT, N, A, LDA,
     $                       VL, LDVL, VR, LDVR,
     $                       N, NOUT, WORK, -1, RWORK, -1, IERR )
               LWORK_TREVC = INT( WORK(1) )
               MAXWRK = MAX( MAXWRK, N + LWORK_TREVC )
               CALL ZHSEQR( 'S', 'V', N, 1, N, A, LDA, W, VL, LDVL,
     $                      WORK, -1, INFO )
            ELSE IF( WANTVR ) THEN
               MAXWRK = MAX( MAXWRK, N + ( N - 1 )*ILAENV( 1, 'ZUNGHR',
     $                       ' ', N, 1, N, -1 ) )
               CALL ZTREVC3( 'R', 'B', SELECT, N, A, LDA,
     $                       VL, LDVL, VR, LDVR,
     $                       N, NOUT, WORK, -1, RWORK, -1, IERR )
               LWORK_TREVC = INT( WORK(1) )
               MAXWRK = MAX( MAXWRK, N + LWORK_TREVC )
               CALL ZHSEQR( 'S', 'V', N, 1, N, A, LDA, W, VR, LDVR,
     $                      WORK, -1, INFO )
            ELSE
               CALL ZHSEQR( 'E', 'N', N, 1, N, A, LDA, W, VR, LDVR,
     $                      WORK, -1, INFO )
            END IF
            HSWORK = INT( WORK(1) )
            MAXWRK = MAX( MAXWRK, HSWORK, MINWRK )
         END IF
         WORK( 1 ) = MAXWRK
*
         IF( LWORK.LT.MINWRK .AND. .NOT.LQUERY ) THEN
            INFO = -12
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGEEV ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Get machine constants
*
      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
      SMLNUM = SQRT( SMLNUM ) / EPS
      BIGNUM = ONE / SMLNUM
*
*     Scale A if max element outside range [SMLNUM,BIGNUM]
*
      ANRM = ZLANGE( 'M', N, N, A, LDA, DUM )
      SCALEA = .FALSE.
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
         SCALEA = .TRUE.
         CSCALE = SMLNUM
      ELSE IF( ANRM.GT.BIGNUM ) THEN
         SCALEA = .TRUE.
         CSCALE = BIGNUM
      END IF
      IF( SCALEA )
     $   CALL ZLASCL( 'G', 0, 0, ANRM, CSCALE, N, N, A, LDA, IERR )
*
*     Balance the matrix
*     (CWorkspace: none)
*     (RWorkspace: need N)
*
      IBAL = 1
      CALL ZGEBAL( 'B', N, A, LDA, ILO, IHI, RWORK( IBAL ), IERR )
*
*     Reduce to upper Hessenberg form
*     (CWorkspace: need 2*N, prefer N+N*NB)
*     (RWorkspace: none)
*
      ITAU = 1
      IWRK = ITAU + N
      CALL ZGEHRD( N, ILO, IHI, A, LDA, WORK( ITAU ), WORK( IWRK ),
     $             LWORK-IWRK+1, IERR )
*
      IF( WANTVL ) THEN
*
*        Want left eigenvectors
*        Copy Householder vectors to VL
*
         SIDE = 'L'
         CALL ZLACPY( 'L', N, N, A, LDA, VL, LDVL )
*
*        Generate unitary matrix in VL
*        (CWorkspace: need 2*N-1, prefer N+(N-1)*NB)
*        (RWorkspace: none)
*
         CALL ZUNGHR( N, ILO, IHI, VL, LDVL, WORK( ITAU ), WORK( IWRK ),
     $                LWORK-IWRK+1, IERR )
*
*        Perform QR iteration, accumulating Schur vectors in VL
*        (CWorkspace: need 1, prefer HSWORK (see comments) )
*        (RWorkspace: none)
*
         IWRK = ITAU
         CALL ZHSEQR( 'S', 'V', N, ILO, IHI, A, LDA, W, VL, LDVL,
     $                WORK( IWRK ), LWORK-IWRK+1, INFO )
*
         IF( WANTVR ) THEN
*
*           Want left and right eigenvectors
*           Copy Schur vectors to VR
*
            SIDE = 'B'
            CALL ZLACPY( 'F', N, N, VL, LDVL, VR, LDVR )
         END IF
*
      ELSE IF( WANTVR ) THEN
*
*        Want right eigenvectors
*        Copy Householder vectors to VR
*
         SIDE = 'R'
         CALL ZLACPY( 'L', N, N, A, LDA, VR, LDVR )
*
*        Generate unitary matrix in VR
*        (CWorkspace: need 2*N-1, prefer N+(N-1)*NB)
*        (RWorkspace: none)
*
         CALL ZUNGHR( N, ILO, IHI, VR, LDVR, WORK( ITAU ), WORK( IWRK ),
     $                LWORK-IWRK+1, IERR )
*
*        Perform QR iteration, accumulating Schur vectors in VR
*        (CWorkspace: need 1, prefer HSWORK (see comments) )
*        (RWorkspace: none)
*
         IWRK = ITAU
         CALL ZHSEQR( 'S', 'V', N, ILO, IHI, A, LDA, W, VR, LDVR,
     $                WORK( IWRK ), LWORK-IWRK+1, INFO )
*
      ELSE
*
*        Compute eigenvalues only
*        (CWorkspace: need 1, prefer HSWORK (see comments) )
*        (RWorkspace: none)
*
         IWRK = ITAU
         CALL ZHSEQR( 'E', 'N', N, ILO, IHI, A, LDA, W, VR, LDVR,
     $                WORK( IWRK ), LWORK-IWRK+1, INFO )
      END IF
*
*     If INFO .NE. 0 from ZHSEQR, then quit
*
      IF( INFO.NE.0 )
     $   GO TO 50
*
      IF( WANTVL .OR. WANTVR ) THEN
*
*        Compute left and/or right eigenvectors
*        (CWorkspace: need 2*N, prefer N + 2*N*NB)
*        (RWorkspace: need 2*N)
*
         IRWORK = IBAL + N
         CALL ZTREVC3( SIDE, 'B', SELECT, N, A, LDA, VL, LDVL, VR, LDVR,
     $                 N, NOUT, WORK( IWRK ), LWORK-IWRK+1,
     $                 RWORK( IRWORK ), N, IERR )
      END IF
*
      IF( WANTVL ) THEN
*
*        Undo balancing of left eigenvectors
*        (CWorkspace: none)
*        (RWorkspace: need N)
*
         CALL ZGEBAK( 'B', 'L', N, ILO, IHI, RWORK( IBAL ), N, VL, LDVL,
     $                IERR )
*
*        Normalize left eigenvectors and make largest component real
*
         DO 20 I = 1, N
            SCL = ONE / DZNRM2( N, VL( 1, I ), 1 )
            CALL ZDSCAL( N, SCL, VL( 1, I ), 1 )
            DO 10 K = 1, N
               RWORK( IRWORK+K-1 ) = DBLE( VL( K, I ) )**2 +
     $                               AIMAG( VL( K, I ) )**2
   10       CONTINUE
            K = IDAMAX( N, RWORK( IRWORK ), 1 )
            TMP = CONJG( VL( K, I ) ) / SQRT( RWORK( IRWORK+K-1 ) )
            CALL ZSCAL( N, TMP, VL( 1, I ), 1 )
            VL( K, I ) = DCMPLX( DBLE( VL( K, I ) ), ZERO )
   20    CONTINUE
      END IF
*
      IF( WANTVR ) THEN
*
*        Undo balancing of right eigenvectors
*        (CWorkspace: none)
*        (RWorkspace: need N)
*
         CALL ZGEBAK( 'B', 'R', N, ILO, IHI, RWORK( IBAL ), N, VR, LDVR,
     $                IERR )
*
*        Normalize right eigenvectors and make largest component real
*
         DO 40 I = 1, N
            SCL = ONE / DZNRM2( N, VR( 1, I ), 1 )
            CALL ZDSCAL( N, SCL, VR( 1, I ), 1 )
            DO 30 K = 1, N
               RWORK( IRWORK+K-1 ) = DBLE( VR( K, I ) )**2 +
     $                               AIMAG( VR( K, I ) )**2
   30       CONTINUE
            K = IDAMAX( N, RWORK( IRWORK ), 1 )
            TMP = CONJG( VR( K, I ) ) / SQRT( RWORK( IRWORK+K-1 ) )
            CALL ZSCAL( N, TMP, VR( 1, I ), 1 )
            VR( K, I ) = DCMPLX( DBLE( VR( K, I ) ), ZERO )
   40    CONTINUE
      END IF
*
*     Undo scaling if necessary
*
   50 CONTINUE
      IF( SCALEA ) THEN
         CALL ZLASCL( 'G', 0, 0, CSCALE, ANRM, N-INFO, 1, W( INFO+1 ),
     $                MAX( N-INFO, 1 ), IERR )
         IF( INFO.GT.0 ) THEN
            CALL ZLASCL( 'G', 0, 0, CSCALE, ANRM, ILO-1, 1, W, N, IERR )
         END IF
      END IF
*
      WORK( 1 ) = MAXWRK
      RETURN
*
*     End of ZGEEV
*
      END

! ZGETRF      
      SUBROUTINE ZGETRF( M, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IINFO, J, JB, NB
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZGEMM, ZGETRF2, ZLASWP, ZTRSM
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGETRF', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     Determine the block size for this environment.
*
      NB = ILAENV( 1, 'ZGETRF', ' ', M, N, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN
*
*        Use unblocked code.
*
         CALL ZGETRF2( M, N, A, LDA, IPIV, INFO )
      ELSE
*
*        Use blocked code.
*
         DO 20 J = 1, MIN( M, N ), NB
            JB = MIN( MIN( M, N )-J+1, NB )
*
*           Factor diagonal and subdiagonal blocks and test for exact
*           singularity.
*
            CALL ZGETRF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
*
*           Adjust INFO and the pivot indices.
*
            IF( INFO.EQ.0 .AND. IINFO.GT.0 )
     $         INFO = IINFO + J - 1
            DO 10 I = J, MIN( M, J+JB-1 )
               IPIV( I ) = J - 1 + IPIV( I )
   10       CONTINUE
*
*           Apply interchanges to columns 1:J-1.
*
            CALL ZLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
*
            IF( J+JB.LE.N ) THEN
*
*              Apply interchanges to columns J+JB:N.
*
               CALL ZLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1,
     $                      IPIV, 1 )
*
*              Compute block row of U.
*
               CALL ZTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB,
     $                     N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ),
     $                     LDA )
               IF( J+JB.LE.M ) THEN
*
*                 Update trailing submatrix.
*
                  CALL ZGEMM( 'No transpose', 'No transpose', M-J-JB+1,
     $                        N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA,
     $                        A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ),
     $                        LDA )
               END IF
            END IF
   20    CONTINUE
      END IF
      RETURN
*
*     End of ZGETRF
*
      END
      
! ZGETRI      
      SUBROUTINE ZGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ),
     $                   ONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IWS, J, JB, JJ, JP, LDWORK, LWKOPT, NB,
     $                   NBMIN, NN
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZGEMM, ZGEMV, ZSWAP, ZTRSM, ZTRTRI
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      NB = ILAENV( 1, 'ZGETRI', ' ', N, -1, -1, -1 )
      LWKOPT = N*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -3
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGETRI', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Form inv(U).  If INFO > 0 from ZTRTRI, then U is singular,
*     and the inverse is not computed.
*
      CALL ZTRTRI( 'Upper', 'Non-unit', N, A, LDA, INFO )
      IF( INFO.GT.0 )
     $   RETURN
*
      NBMIN = 2
      LDWORK = N
      IF( NB.GT.1 .AND. NB.LT.N ) THEN
         IWS = MAX( LDWORK*NB, 1 )
         IF( LWORK.LT.IWS ) THEN
            NB = LWORK / LDWORK
            NBMIN = MAX( 2, ILAENV( 2, 'ZGETRI', ' ', N, -1, -1, -1 ) )
         END IF
      ELSE
         IWS = N
      END IF
*
*     Solve the equation inv(A)*L = inv(U) for inv(A).
*
      IF( NB.LT.NBMIN .OR. NB.GE.N ) THEN
*
*        Use unblocked code.
*
         DO 20 J = N, 1, -1
*
*           Copy current column of L to WORK and replace with zeros.
*
            DO 10 I = J + 1, N
               WORK( I ) = A( I, J )
               A( I, J ) = ZERO
   10       CONTINUE
*
*           Compute current column of inv(A).
*
            IF( J.LT.N )
     $         CALL ZGEMV( 'No transpose', N, N-J, -ONE, A( 1, J+1 ),
     $                     LDA, WORK( J+1 ), 1, ONE, A( 1, J ), 1 )
   20    CONTINUE
      ELSE
*
*        Use blocked code.
*
         NN = ( ( N-1 ) / NB )*NB + 1
         DO 50 J = NN, 1, -NB
            JB = MIN( NB, N-J+1 )
*
*           Copy current block column of L to WORK and replace with
*           zeros.
*
            DO 40 JJ = J, J + JB - 1
               DO 30 I = JJ + 1, N
                  WORK( I+( JJ-J )*LDWORK ) = A( I, JJ )
                  A( I, JJ ) = ZERO
   30          CONTINUE
   40       CONTINUE
*
*           Compute current block column of inv(A).
*
            IF( J+JB.LE.N )
     $         CALL ZGEMM( 'No transpose', 'No transpose', N, JB,
     $                     N-J-JB+1, -ONE, A( 1, J+JB ), LDA,
     $                     WORK( J+JB ), LDWORK, ONE, A( 1, J ), LDA )
            CALL ZTRSM( 'Right', 'Lower', 'No transpose', 'Unit', N, JB,
     $                  ONE, WORK( J ), LDWORK, A( 1, J ), LDA )
   50    CONTINUE
      END IF
*
*     Apply column interchanges.
*
      DO 60 J = N - 1, 1, -1
         JP = IPIV( J )
         IF( JP.NE.J )
     $      CALL ZSWAP( N, A( 1, J ), 1, A( 1, JP ), 1 )
   60 CONTINUE
*
      WORK( 1 ) = IWS
      RETURN
*
*     End of ZGETRI
*
      END

! ZTREVC3
      SUBROUTINE ZTREVC3( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR,
     $                    LDVR, MM, M, WORK, LWORK, RWORK, LRWORK, INFO)
      IMPLICIT NONE
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          HOWMNY, SIDE
      INTEGER            INFO, LDT, LDVL, LDVR, LWORK, LRWORK, M, MM, N
*     ..
*     .. Array Arguments ..
      LOGICAL            SELECT( * )
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ),
     $                   WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ),
     $                     CONE  = ( 1.0D+0, 0.0D+0 ) )
      INTEGER            NBMIN, NBMAX
      PARAMETER          ( NBMIN = 8, NBMAX = 128 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ALLV, BOTHV, LEFTV, LQUERY, OVER, RIGHTV, SOMEV
      INTEGER            I, II, IS, J, K, KI, IV, MAXWRK, NB
      DOUBLE PRECISION   OVFL, REMAX, SCALE, SMIN, SMLNUM, ULP, UNFL
      COMPLEX*16         CDUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV, IZAMAX
      DOUBLE PRECISION   DLAMCH, DZASUM
      EXTERNAL           LSAME, ILAENV, IZAMAX, DLAMCH, DZASUM
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZCOPY, ZDSCAL, ZGEMV, ZLATRS,
     $                   ZGEMM, DLABAD, ZLASET, ZLACPY
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCMPLX, CONJG, AIMAG, MAX
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
*     ..
*     .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( DBLE( CDUM ) ) + ABS( AIMAG( CDUM ) )
*     ..
*     .. Executable Statements ..
*
*     Decode and test the input parameters
*
      BOTHV  = LSAME( SIDE, 'B' )
      RIGHTV = LSAME( SIDE, 'R' ) .OR. BOTHV
      LEFTV  = LSAME( SIDE, 'L' ) .OR. BOTHV
*
      ALLV  = LSAME( HOWMNY, 'A' )
      OVER  = LSAME( HOWMNY, 'B' )
      SOMEV = LSAME( HOWMNY, 'S' )
*
*     Set M to the number of columns required to store the selected
*     eigenvectors.
*
      IF( SOMEV ) THEN
         M = 0
         DO 10 J = 1, N
            IF( SELECT( J ) )
     $         M = M + 1
   10    CONTINUE
      ELSE
         M = N
      END IF
*
      INFO = 0
      NB = ILAENV( 1, 'ZTREVC', SIDE // HOWMNY, N, -1, -1, -1 )
      MAXWRK = N + 2*N*NB
      WORK(1) = MAXWRK
      RWORK(1) = N
      LQUERY = ( LWORK.EQ.-1 .OR. LRWORK.EQ.-1 )
      IF( .NOT.RIGHTV .AND. .NOT.LEFTV ) THEN
         INFO = -1
      ELSE IF( .NOT.ALLV .AND. .NOT.OVER .AND. .NOT.SOMEV ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDT.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDVL.LT.1 .OR. ( LEFTV .AND. LDVL.LT.N ) ) THEN
         INFO = -8
      ELSE IF( LDVR.LT.1 .OR. ( RIGHTV .AND. LDVR.LT.N ) ) THEN
         INFO = -10
      ELSE IF( MM.LT.M ) THEN
         INFO = -11
      ELSE IF( LWORK.LT.MAX( 1, 2*N ) .AND. .NOT.LQUERY ) THEN
         INFO = -14
      ELSE IF ( LRWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -16
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZTREVC3', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Use blocked version of back-transformation if sufficient workspace.
*     Zero-out the workspace to avoid potential NaN propagation.
*
      IF( OVER .AND. LWORK .GE. N + 2*N*NBMIN ) THEN
         NB = (LWORK - N) / (2*N)
         NB = MIN( NB, NBMAX )
         CALL ZLASET( 'F', N, 1+2*NB, CZERO, CZERO, WORK, N )
      ELSE
         NB = 1
      END IF
*
*     Set the constants to control overflow.
*
      UNFL = DLAMCH( 'Safe minimum' )
      OVFL = ONE / UNFL
      CALL DLABAD( UNFL, OVFL )
      ULP = DLAMCH( 'Precision' )
      SMLNUM = UNFL*( N / ULP )
*
*     Store the diagonal elements of T in working array WORK.
*
      DO 20 I = 1, N
         WORK( I ) = T( I, I )
   20 CONTINUE
*
*     Compute 1-norm of each column of strictly upper triangular
*     part of T to control overflow in triangular solver.
*
      RWORK( 1 ) = ZERO
      DO 30 J = 2, N
         RWORK( J ) = DZASUM( J-1, T( 1, J ), 1 )
   30 CONTINUE
*
      IF( RIGHTV ) THEN
*
*        ============================================================
*        Compute right eigenvectors.
*
*        IV is index of column in current block.
*        Non-blocked version always uses IV=NB=1;
*        blocked     version starts with IV=NB, goes down to 1.
*        (Note the "0-th" column is used to store the original diagonal.)
         IV = NB
         IS = M
         DO 80 KI = N, 1, -1
            IF( SOMEV ) THEN
               IF( .NOT.SELECT( KI ) )
     $            GO TO 80
            END IF
            SMIN = MAX( ULP*( CABS1( T( KI, KI ) ) ), SMLNUM )
*
*           --------------------------------------------------------
*           Complex right eigenvector
*
            WORK( KI + IV*N ) = CONE
*
*           Form right-hand side.
*
            DO 40 K = 1, KI - 1
               WORK( K + IV*N ) = -T( K, KI )
   40       CONTINUE
*
*           Solve upper triangular system:
*           [ T(1:KI-1,1:KI-1) - T(KI,KI) ]*X = SCALE*WORK.
*
            DO 50 K = 1, KI - 1
               T( K, K ) = T( K, K ) - T( KI, KI )
               IF( CABS1( T( K, K ) ).LT.SMIN )
     $            T( K, K ) = SMIN
   50       CONTINUE
*
            IF( KI.GT.1 ) THEN
               CALL ZLATRS( 'Upper', 'No transpose', 'Non-unit', 'Y',
     $                      KI-1, T, LDT, WORK( 1 + IV*N ), SCALE,
     $                      RWORK, INFO )
               WORK( KI + IV*N ) = SCALE
            END IF
*
*           Copy the vector x or Q*x to VR and normalize.
*
            IF( .NOT.OVER ) THEN
*              ------------------------------
*              no back-transform: copy x to VR and normalize.
               CALL ZCOPY( KI, WORK( 1 + IV*N ), 1, VR( 1, IS ), 1 )
*
               II = IZAMAX( KI, VR( 1, IS ), 1 )
               REMAX = ONE / CABS1( VR( II, IS ) )
               CALL ZDSCAL( KI, REMAX, VR( 1, IS ), 1 )
*
               DO 60 K = KI + 1, N
                  VR( K, IS ) = CZERO
   60          CONTINUE
*
            ELSE IF( NB.EQ.1 ) THEN
*              ------------------------------
*              version 1: back-transform each vector with GEMV, Q*x.
               IF( KI.GT.1 )
     $            CALL ZGEMV( 'N', N, KI-1, CONE, VR, LDVR,
     $                        WORK( 1 + IV*N ), 1, DCMPLX( SCALE ),
     $                        VR( 1, KI ), 1 )
*
               II = IZAMAX( N, VR( 1, KI ), 1 )
               REMAX = ONE / CABS1( VR( II, KI ) )
               CALL ZDSCAL( N, REMAX, VR( 1, KI ), 1 )
*
            ELSE
*              ------------------------------
*              version 2: back-transform block of vectors with GEMM
*              zero out below vector
               DO K = KI + 1, N
                  WORK( K + IV*N ) = CZERO
               END DO
*
*              Columns IV:NB of work are valid vectors.
*              When the number of vectors stored reaches NB,
*              or if this was last vector, do the GEMM
               IF( (IV.EQ.1) .OR. (KI.EQ.1) ) THEN
                  CALL ZGEMM( 'N', 'N', N, NB-IV+1, KI+NB-IV, CONE,
     $                        VR, LDVR,
     $                        WORK( 1 + (IV)*N    ), N,
     $                        CZERO,
     $                        WORK( 1 + (NB+IV)*N ), N )
*                 normalize vectors
                  DO K = IV, NB
                     II = IZAMAX( N, WORK( 1 + (NB+K)*N ), 1 )
                     REMAX = ONE / CABS1( WORK( II + (NB+K)*N ) )
                     CALL ZDSCAL( N, REMAX, WORK( 1 + (NB+K)*N ), 1 )
                  END DO
                  CALL ZLACPY( 'F', N, NB-IV+1,
     $                         WORK( 1 + (NB+IV)*N ), N,
     $                         VR( 1, KI ), LDVR )
                  IV = NB
               ELSE
                  IV = IV - 1
               END IF
            END IF
*
*           Restore the original diagonal elements of T.
*
            DO 70 K = 1, KI - 1
               T( K, K ) = WORK( K )
   70       CONTINUE
*
            IS = IS - 1
   80    CONTINUE
      END IF
*
      IF( LEFTV ) THEN
*
*        ============================================================
*        Compute left eigenvectors.
*
*        IV is index of column in current block.
*        Non-blocked version always uses IV=1;
*        blocked     version starts with IV=1, goes up to NB.
*        (Note the "0-th" column is used to store the original diagonal.)
         IV = 1
         IS = 1
         DO 130 KI = 1, N
*
            IF( SOMEV ) THEN
               IF( .NOT.SELECT( KI ) )
     $            GO TO 130
            END IF
            SMIN = MAX( ULP*( CABS1( T( KI, KI ) ) ), SMLNUM )
*
*           --------------------------------------------------------
*           Complex left eigenvector
*
            WORK( KI + IV*N ) = CONE
*
*           Form right-hand side.
*
            DO 90 K = KI + 1, N
               WORK( K + IV*N ) = -CONJG( T( KI, K ) )
   90       CONTINUE
*
*           Solve conjugate-transposed triangular system:
*           [ T(KI+1:N,KI+1:N) - T(KI,KI) ]**H * X = SCALE*WORK.
*
            DO 100 K = KI + 1, N
               T( K, K ) = T( K, K ) - T( KI, KI )
               IF( CABS1( T( K, K ) ).LT.SMIN )
     $            T( K, K ) = SMIN
  100       CONTINUE
*
            IF( KI.LT.N ) THEN
               CALL ZLATRS( 'Upper', 'Conjugate transpose', 'Non-unit',
     $                      'Y', N-KI, T( KI+1, KI+1 ), LDT,
     $                      WORK( KI+1 + IV*N ), SCALE, RWORK, INFO )
               WORK( KI + IV*N ) = SCALE
            END IF
*
*           Copy the vector x or Q*x to VL and normalize.
*
            IF( .NOT.OVER ) THEN
*              ------------------------------
*              no back-transform: copy x to VL and normalize.
               CALL ZCOPY( N-KI+1, WORK( KI + IV*N ), 1, VL(KI,IS), 1 )
*
               II = IZAMAX( N-KI+1, VL( KI, IS ), 1 ) + KI - 1
               REMAX = ONE / CABS1( VL( II, IS ) )
               CALL ZDSCAL( N-KI+1, REMAX, VL( KI, IS ), 1 )
*
               DO 110 K = 1, KI - 1
                  VL( K, IS ) = CZERO
  110          CONTINUE
*
            ELSE IF( NB.EQ.1 ) THEN
*              ------------------------------
*              version 1: back-transform each vector with GEMV, Q*x.
               IF( KI.LT.N )
     $            CALL ZGEMV( 'N', N, N-KI, CONE, VL( 1, KI+1 ), LDVL,
     $                        WORK( KI+1 + IV*N ), 1, DCMPLX( SCALE ),
     $                        VL( 1, KI ), 1 )
*
               II = IZAMAX( N, VL( 1, KI ), 1 )
               REMAX = ONE / CABS1( VL( II, KI ) )
               CALL ZDSCAL( N, REMAX, VL( 1, KI ), 1 )
*
            ELSE
*              ------------------------------
*              version 2: back-transform block of vectors with GEMM
*              zero out above vector
*              could go from KI-NV+1 to KI-1
               DO K = 1, KI - 1
                  WORK( K + IV*N ) = CZERO
               END DO
*
*              Columns 1:IV of work are valid vectors.
*              When the number of vectors stored reaches NB,
*              or if this was last vector, do the GEMM
               IF( (IV.EQ.NB) .OR. (KI.EQ.N) ) THEN
                  CALL ZGEMM( 'N', 'N', N, IV, N-KI+IV, CONE,
     $                        VL( 1, KI-IV+1 ), LDVL,
     $                        WORK( KI-IV+1 + (1)*N ), N,
     $                        CZERO,
     $                        WORK( 1 + (NB+1)*N ), N )
*                 normalize vectors
                  DO K = 1, IV
                     II = IZAMAX( N, WORK( 1 + (NB+K)*N ), 1 )
                     REMAX = ONE / CABS1( WORK( II + (NB+K)*N ) )
                     CALL ZDSCAL( N, REMAX, WORK( 1 + (NB+K)*N ), 1 )
                  END DO
                  CALL ZLACPY( 'F', N, IV,
     $                         WORK( 1 + (NB+1)*N ), N,
     $                         VL( 1, KI-IV+1 ), LDVL )
                  IV = 1
               ELSE
                  IV = IV + 1
               END IF
            END IF
*
*           Restore the original diagonal elements of T.
*
            DO 120 K = KI + 1, N
               T( K, K ) = WORK( K )
  120       CONTINUE
*
            IS = IS + 1
  130    CONTINUE
      END IF
*
      RETURN
*
*     End of ZTREVC3
*
      END

! ZHSEQR
      SUBROUTINE zhseqr( JOB, COMPZ, N, ILO, IHI, H, LDH, W, Z, LDZ,
     $                   WORK, LWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            IHI, ILO, INFO, LDH, LDZ, LWORK, N
      CHARACTER          COMPZ, JOB
*     ..
*     .. Array Arguments ..
      COMPLEX*16         H( LDH, * ), W( * ), WORK( * ), Z( LDZ, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
*
*     ==== Matrices of order NTINY or smaller must be processed by
*     .    ZLAHQR because of insufficient subdiagonal scratch space.
*     .    (This is a hard limit.) ====
      INTEGER            NTINY
      parameter( ntiny = 15 )
*
*     ==== NL allocates some local workspace to help small matrices
*     .    through a rare ZLAHQR failure.  NL > NTINY = 15 is
*     .    required and NL <= NMIN = ILAENV(ISPEC=12,...) is recom-
*     .    mended.  (The default value of NMIN is 75.)  Using NL = 49
*     .    allows up to six simultaneous shifts and a 16-by-16
*     .    deflation window.  ====
      INTEGER            NL
      parameter( nl = 49 )
      COMPLEX*16         ZERO, ONE
      parameter( zero = ( 0.0d0, 0.0d0 ),
     $                   one = ( 1.0d0, 0.0d0 ) )
      DOUBLE PRECISION   RZERO
      parameter( rzero = 0.0d0 )
*     ..
*     .. Local Arrays ..
      COMPLEX*16         HL( NL, NL ), WORKL( NL )
*     ..
*     .. Local Scalars ..
      INTEGER            KBOT, NMIN
      LOGICAL            INITZ, LQUERY, WANTT, WANTZ
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      LOGICAL            LSAME
      EXTERNAL           ilaenv, lsame
*     ..
*     .. External Subroutines ..
      EXTERNAL           xerbla, zcopy, zlacpy, zlahqr, zlaqr0, zlaset
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          dble, dcmplx, max, min
*     ..
*     .. Executable Statements ..
*
*     ==== Decode and check the input parameters. ====
*
      wantt = lsame( job, 'S' )
      initz = lsame( compz, 'I' )
      wantz = initz .OR. lsame( compz, 'V' )
      work( 1 ) = dcmplx( dble( max( 1, n ) ), rzero )
      lquery = lwork.EQ.-1
*
      info = 0
      IF( .NOT.lsame( job, 'E' ) .AND. .NOT.wantt ) THEN
         info = -1
      ELSE IF( .NOT.lsame( compz, 'N' ) .AND. .NOT.wantz ) THEN
         info = -2
      ELSE IF( n.LT.0 ) THEN
         info = -3
      ELSE IF( ilo.LT.1 .OR. ilo.GT.max( 1, n ) ) THEN
         info = -4
      ELSE IF( ihi.LT.min( ilo, n ) .OR. ihi.GT.n ) THEN
         info = -5
      ELSE IF( ldh.LT.max( 1, n ) ) THEN
         info = -7
      ELSE IF( ldz.LT.1 .OR. ( wantz .AND. ldz.LT.max( 1, n ) ) ) THEN
         info = -10
      ELSE IF( lwork.LT.max( 1, n ) .AND. .NOT.lquery ) THEN
         info = -12
      END IF
*
      IF( info.NE.0 ) THEN
*
*        ==== Quick return in case of invalid argument. ====
*
         CALL xerbla( 'ZHSEQR', -info )
         RETURN
*
      ELSE IF( n.EQ.0 ) THEN
*
*        ==== Quick return in case N = 0; nothing to do. ====
*
         RETURN
*
      ELSE IF( lquery ) THEN
*
*        ==== Quick return in case of a workspace query ====
*
         CALL zlaqr0( wantt, wantz, n, ilo, ihi, h, ldh, w, ilo, ihi, z,
     $                ldz, work, lwork, info )
*        ==== Ensure reported workspace size is backward-compatible with
*        .    previous LAPACK versions. ====
         work( 1 ) = dcmplx( max( dble( work( 1 ) ), dble( max( 1,
     $               n ) ) ), rzero )
         RETURN
*
      ELSE
*
*        ==== copy eigenvalues isolated by ZGEBAL ====
*
         IF( ilo.GT.1 )
     $      CALL zcopy( ilo-1, h, ldh+1, w, 1 )
         IF( ihi.LT.n )
     $      CALL zcopy( n-ihi, h( ihi+1, ihi+1 ), ldh+1, w( ihi+1 ), 1 )
*
*        ==== Initialize Z, if requested ====
*
         IF( initz )
     $      CALL zlaset( 'A', n, n, zero, one, z, ldz )
*
*        ==== Quick return if possible ====
*
         IF( ilo.EQ.ihi ) THEN
            w( ilo ) = h( ilo, ilo )
            RETURN
         END IF
*
*        ==== ZLAHQR/ZLAQR0 crossover point ====
*
         nmin = ilaenv( 12, 'ZHSEQR', job( : 1 ) // compz( : 1 ), n,
     $          ilo, ihi, lwork )
         nmin = max( ntiny, nmin )
*
*        ==== ZLAQR0 for big matrices; ZLAHQR for small ones ====
*
         IF( n.GT.nmin ) THEN
            CALL zlaqr0( wantt, wantz, n, ilo, ihi, h, ldh, w, ilo, ihi,
     $                   z, ldz, work, lwork, info )
         ELSE
*
*           ==== Small matrix ====
*
            CALL zlahqr( wantt, wantz, n, ilo, ihi, h, ldh, w, ilo, ihi,
     $                   z, ldz, info )
*
            IF( info.GT.0 ) THEN
*
*              ==== A rare ZLAHQR failure!  ZLAQR0 sometimes succeeds
*              .    when ZLAHQR fails. ====
*
               kbot = info
*
               IF( n.GE.nl ) THEN
*
*                 ==== Larger matrices have enough subdiagonal scratch
*                 .    space to call ZLAQR0 directly. ====
*
                  CALL zlaqr0( wantt, wantz, n, ilo, kbot, h, ldh, w,
     $                         ilo, ihi, z, ldz, work, lwork, info )
*
               ELSE
*
*                 ==== Tiny matrices don't have enough subdiagonal
*                 .    scratch space to benefit from ZLAQR0.  Hence,
*                 .    tiny matrices must be copied into a larger
*                 .    array before calling ZLAQR0. ====
*
                  CALL zlacpy( 'A', n, n, h, ldh, hl, nl )
                  hl( n+1, n ) = zero
                  CALL zlaset( 'A', nl, nl-n, zero, zero, hl( 1, n+1 ),
     $                         nl )
                  CALL zlaqr0( wantt, wantz, nl, ilo, kbot, hl, nl, w,
     $                         ilo, ihi, z, ldz, workl, nl, info )
                  IF( wantt .OR. info.NE.0 )
     $               CALL zlacpy( 'A', n, n, hl, nl, h, ldh )
               END IF
            END IF
         END IF
*
*        ==== Clear out the trash, if necessary. ====
*
         IF( ( wantt .OR. info.NE.0 ) .AND. n.GT.2 )
     $      CALL zlaset( 'L', n-2, n-2, zero, zero, h( 3, 1 ), ldh )
*
*        ==== Ensure reported workspace size is backward-compatible with
*        .    previous LAPACK versions. ====
*
         work( 1 ) = dcmplx( max( dble( max( 1, n ) ),
     $               dble( work( 1 ) ) ), rzero )
      END IF
*
*     ==== End of ZHSEQR ====
*
      END

! ZLANGE      
      DOUBLE PRECISION FUNCTION zlange( NORM, M, N, A, LDA, WORK )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          norm
      INTEGER            lda, m, n
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   work( * )
      COMPLEX*16         a( lda, * )
*     ..
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   one, zero
      parameter( one = 1.0d+0, zero = 0.0d+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            i, j
      DOUBLE PRECISION   scale, sum, VALUE, temp
*     ..
*     .. External Functions ..
      LOGICAL            lsame, disnan
      EXTERNAL           lsame, disnan
*     ..
*     .. External Subroutines ..
      EXTERNAL           zlassq
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, min, sqrt
*     ..
*     .. Executable Statements ..
*
      IF( min( m, n ).EQ.0 ) THEN
         VALUE = zero
      ELSE IF( lsame( norm, 'M' ) ) THEN
*
*        Find max(abs(A(i,j))).
*
         VALUE = zero
         DO 20 j = 1, n
            DO 10 i = 1, m
               temp = abs( a( i, j ) )
               IF( VALUE.LT.temp .OR. disnan( temp ) ) VALUE = temp
   10       CONTINUE
   20    CONTINUE
      ELSE IF( ( lsame( norm, 'O' ) ) .OR. ( norm.EQ.'1' ) ) THEN
*
*        Find norm1(A).
*
         VALUE = zero
         DO 40 j = 1, n
            sum = zero
            DO 30 i = 1, m
               sum = sum + abs( a( i, j ) )
   30       CONTINUE
            IF( VALUE.LT.sum .OR. disnan( sum ) ) VALUE = sum
   40    CONTINUE
      ELSE IF( lsame( norm, 'I' ) ) THEN
*
*        Find normI(A).
*
         DO 50 i = 1, m
            work( i ) = zero
   50    CONTINUE
         DO 70 j = 1, n
            DO 60 i = 1, m
               work( i ) = work( i ) + abs( a( i, j ) )
   60       CONTINUE
   70    CONTINUE
         VALUE = zero
         DO 80 i = 1, m
            temp = work( i )
            IF( VALUE.LT.temp .OR. disnan( temp ) ) VALUE = temp
   80    CONTINUE
      ELSE IF( ( lsame( norm, 'F' ) ) .OR. ( lsame( norm, 'E' ) ) ) THEN
*
*        Find normF(A).
*
         scale = zero
         sum = one
         DO 90 j = 1, n
            CALL zlassq( m, a( 1, j ), 1, scale, sum )
   90    CONTINUE
         VALUE = scale*sqrt( sum )
      END IF
*
      zlange = VALUE
      RETURN
*
*     End of ZLANGE
*
      END

! ZLASCL
      SUBROUTINE zlascl( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          TYPE
      INTEGER            INFO, KL, KU, LDA, M, N
      DOUBLE PRECISION   CFROM, CTO
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      parameter( zero = 0.0d0, one = 1.0d0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            DONE
      INTEGER            I, ITYPE, J, K1, K2, K3, K4
      DOUBLE PRECISION   BIGNUM, CFROM1, CFROMC, CTO1, CTOC, MUL, SMLNUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME, DISNAN
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           lsame, dlamch, disnan
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, max, min
*     ..
*     .. External Subroutines ..
      EXTERNAL           xerbla
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      info = 0
*
      IF( lsame( TYPE, 'G' ) ) then
         itype = 0
      ELSE IF( lsame( TYPE, 'L' ) ) then
         itype = 1
      ELSE IF( lsame( TYPE, 'U' ) ) then
         itype = 2
      ELSE IF( lsame( TYPE, 'H' ) ) then
         itype = 3
      ELSE IF( lsame( TYPE, 'B' ) ) then
         itype = 4
      ELSE IF( lsame( TYPE, 'Q' ) ) then
         itype = 5
      ELSE IF( lsame( TYPE, 'Z' ) ) then
         itype = 6
      ELSE
         itype = -1
      END IF
*
      IF( itype.EQ.-1 ) THEN
         info = -1
      ELSE IF( cfrom.EQ.zero .OR. disnan(cfrom) ) THEN
         info = -4
      ELSE IF( disnan(cto) ) THEN
         info = -5
      ELSE IF( m.LT.0 ) THEN
         info = -6
      ELSE IF( n.LT.0 .OR. ( itype.EQ.4 .AND. n.NE.m ) .OR.
     $         ( itype.EQ.5 .AND. n.NE.m ) ) THEN
         info = -7
      ELSE IF( itype.LE.3 .AND. lda.LT.max( 1, m ) ) THEN
         info = -9
      ELSE IF( itype.GE.4 ) THEN
         IF( kl.LT.0 .OR. kl.GT.max( m-1, 0 ) ) THEN
            info = -2
         ELSE IF( ku.LT.0 .OR. ku.GT.max( n-1, 0 ) .OR.
     $            ( ( itype.EQ.4 .OR. itype.EQ.5 ) .AND. kl.NE.ku ) )
     $             THEN
            info = -3
         ELSE IF( ( itype.EQ.4 .AND. lda.LT.kl+1 ) .OR.
     $            ( itype.EQ.5 .AND. lda.LT.ku+1 ) .OR.
     $            ( itype.EQ.6 .AND. lda.LT.2*kl+ku+1 ) ) THEN
            info = -9
         END IF
      END IF
*
      IF( info.NE.0 ) THEN
         CALL xerbla( 'ZLASCL', -info )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( n.EQ.0 .OR. m.EQ.0 )
     $   RETURN
*
*     Get machine parameters
*
      smlnum = dlamch( 'S' )
      bignum = one / smlnum
*
      cfromc = cfrom
      ctoc = cto
*
   10 CONTINUE
      cfrom1 = cfromc*smlnum
      IF( cfrom1.EQ.cfromc ) THEN
!        CFROMC is an inf.  Multiply by a correctly signed zero for
!        finite CTOC, or a NaN if CTOC is infinite.
         mul = ctoc / cfromc
         done = .true.
         cto1 = ctoc
      ELSE
         cto1 = ctoc / bignum
         IF( cto1.EQ.ctoc ) THEN
!           CTOC is either 0 or an inf.  In both cases, CTOC itself
!           serves as the correct multiplication factor.
            mul = ctoc
            done = .true.
            cfromc = one
         ELSE IF( abs( cfrom1 ).GT.abs( ctoc ) .AND. ctoc.NE.zero ) THEN
            mul = smlnum
            done = .false.
            cfromc = cfrom1
         ELSE IF( abs( cto1 ).GT.abs( cfromc ) ) THEN
            mul = bignum
            done = .false.
            ctoc = cto1
         ELSE
            mul = ctoc / cfromc
            done = .true.
            IF (mul .EQ. one)
     $         RETURN
         END IF
      END IF
*
      IF( itype.EQ.0 ) THEN
*
*        Full matrix
*
         DO 30 j = 1, n
            DO 20 i = 1, m
               a( i, j ) = a( i, j )*mul
   20       CONTINUE
   30    CONTINUE
*
      ELSE IF( itype.EQ.1 ) THEN
*
*        Lower triangular matrix
*
         DO 50 j = 1, n
            DO 40 i = j, m
               a( i, j ) = a( i, j )*mul
   40       CONTINUE
   50    CONTINUE
*
      ELSE IF( itype.EQ.2 ) THEN
*
*        Upper triangular matrix
*
         DO 70 j = 1, n
            DO 60 i = 1, min( j, m )
               a( i, j ) = a( i, j )*mul
   60       CONTINUE
   70    CONTINUE
*
      ELSE IF( itype.EQ.3 ) THEN
*
*        Upper Hessenberg matrix
*
         DO 90 j = 1, n
            DO 80 i = 1, min( j+1, m )
               a( i, j ) = a( i, j )*mul
   80       CONTINUE
   90    CONTINUE
*
      ELSE IF( itype.EQ.4 ) THEN
*
*        Lower half of a symmetric band matrix
*
         k3 = kl + 1
         k4 = n + 1
         DO 110 j = 1, n
            DO 100 i = 1, min( k3, k4-j )
               a( i, j ) = a( i, j )*mul
  100       CONTINUE
  110    CONTINUE
*
      ELSE IF( itype.EQ.5 ) THEN
*
*        Upper half of a symmetric band matrix
*
         k1 = ku + 2
         k3 = ku + 1
         DO 130 j = 1, n
            DO 120 i = max( k1-j, 1 ), k3
               a( i, j ) = a( i, j )*mul
  120       CONTINUE
  130    CONTINUE
*
      ELSE IF( itype.EQ.6 ) THEN
*
*        Band matrix
*
         k1 = kl + ku + 2
         k2 = kl + 1
         k3 = 2*kl + ku + 1
         k4 = kl + ku + 1 + m
         DO 150 j = 1, n
            DO 140 i = max( k1-j, k2 ), min( k3, k4-j )
               a( i, j ) = a( i, j )*mul
  140       CONTINUE
  150    CONTINUE
*
      END IF
*
      IF( .NOT.done )
     $   GO TO 10
*
      RETURN
*
*     End of ZLASCL
*
      END

! ZGEBAL
      SUBROUTINE zgebal( JOB, N, A, LDA, ILO, IHI, SCALE, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          JOB
      INTEGER            IHI, ILO, INFO, LDA, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   SCALE( * )
      COMPLEX*16         A( LDA, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      parameter( zero = 0.0d+0, one = 1.0d+0 )
      DOUBLE PRECISION   SCLFAC
      parameter( sclfac = 2.0d+0 )
      DOUBLE PRECISION   FACTOR
      parameter( factor = 0.95d+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOCONV
      INTEGER            I, ICA, IEXC, IRA, J, K, L, M
      DOUBLE PRECISION   C, CA, F, G, R, RA, S, SFMAX1, SFMAX2, SFMIN1,
     $                   SFMIN2
*     ..
*     .. External Functions ..
      LOGICAL            DISNAN, LSAME
      INTEGER            IZAMAX
      DOUBLE PRECISION   DLAMCH, DZNRM2
      EXTERNAL           disnan, lsame, izamax, dlamch, dznrm2
*     ..
*     .. External Subroutines ..
      EXTERNAL           xerbla, zdscal, zswap
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, dble, dimag, max, min
*
*     Test the input parameters
*
      info = 0
      IF( .NOT.lsame( job, 'N' ) .AND. .NOT.lsame( job, 'P' ) .AND.
     $    .NOT.lsame( job, 'S' ) .AND. .NOT.lsame( job, 'B' ) ) THEN
         info = -1
      ELSE IF( n.LT.0 ) THEN
         info = -2
      ELSE IF( lda.LT.max( 1, n ) ) THEN
         info = -4
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'ZGEBAL', -info )
         RETURN
      END IF
*
      k = 1
      l = n
*
      IF( n.EQ.0 )
     $   GO TO 210
*
      IF( lsame( job, 'N' ) ) THEN
         DO 10 i = 1, n
            scale( i ) = one
   10    CONTINUE
         GO TO 210
      END IF
*
      IF( lsame( job, 'S' ) )
     $   GO TO 120
*
*     Permutation to isolate eigenvalues if possible
*
      GO TO 50
*
*     Row and column exchange.
*
   20 CONTINUE
      scale( m ) = j
      IF( j.EQ.m )
     $   GO TO 30
*
      CALL zswap( l, a( 1, j ), 1, a( 1, m ), 1 )
      CALL zswap( n-k+1, a( j, k ), lda, a( m, k ), lda )
*
   30 CONTINUE
      GO TO ( 40, 80 )iexc
*
*     Search for rows isolating an eigenvalue and push them down.
*
   40 CONTINUE
      IF( l.EQ.1 )
     $   GO TO 210
      l = l - 1
*
   50 CONTINUE
      DO 70 j = l, 1, -1
*
         DO 60 i = 1, l
            IF( i.EQ.j )
     $         GO TO 60
            IF( dble( a( j, i ) ).NE.zero .OR. dimag( a( j, i ) ).NE.
     $          zero )GO TO 70
   60    CONTINUE
*
         m = l
         iexc = 1
         GO TO 20
   70 CONTINUE
*
      GO TO 90
*
*     Search for columns isolating an eigenvalue and push them left.
*
   80 CONTINUE
      k = k + 1
*
   90 CONTINUE
      DO 110 j = k, l
*
         DO 100 i = k, l
            IF( i.EQ.j )
     $         GO TO 100
            IF( dble( a( i, j ) ).NE.zero .OR. dimag( a( i, j ) ).NE.
     $          zero )GO TO 110
  100    CONTINUE
*
         m = k
         iexc = 2
         GO TO 20
  110 CONTINUE
*
  120 CONTINUE
      DO 130 i = k, l
         scale( i ) = one
  130 CONTINUE
*
      IF( lsame( job, 'P' ) )
     $   GO TO 210
*
*     Balance the submatrix in rows K to L.
*
*     Iterative loop for norm reduction
*
      sfmin1 = dlamch( 'S' ) / dlamch( 'P' )
      sfmax1 = one / sfmin1
      sfmin2 = sfmin1*sclfac
      sfmax2 = one / sfmin2
  140 CONTINUE
      noconv = .false.
*
      DO 200 i = k, l
*
         c = dznrm2( l-k+1, a( k, i ), 1 )
         r = dznrm2( l-k+1, a( i, k ), lda )
         ica = izamax( l, a( 1, i ), 1 )
         ca = abs( a( ica, i ) )
         ira = izamax( n-k+1, a( i, k ), lda )
         ra = abs( a( i, ira+k-1 ) )
*
*        Guard against zero C or R due to underflow.
*
         IF( c.EQ.zero .OR. r.EQ.zero )
     $      GO TO 200
         g = r / sclfac
         f = one
         s = c + r
  160    CONTINUE
         IF( c.GE.g .OR. max( f, c, ca ).GE.sfmax2 .OR.
     $       min( r, g, ra ).LE.sfmin2 )GO TO 170
            IF( disnan( c+f+ca+r+g+ra ) ) THEN
*
*           Exit if NaN to avoid infinite loop
*
            info = -3
            CALL xerbla( 'ZGEBAL', -info )
            RETURN
         END IF
         f = f*sclfac
         c = c*sclfac
         ca = ca*sclfac
         r = r / sclfac
         g = g / sclfac
         ra = ra / sclfac
         GO TO 160
*
  170    CONTINUE
         g = c / sclfac
  180    CONTINUE
         IF( g.LT.r .OR. max( r, ra ).GE.sfmax2 .OR.
     $       min( f, c, g, ca ).LE.sfmin2 )GO TO 190
         f = f / sclfac
         c = c / sclfac
         g = g / sclfac
         ca = ca / sclfac
         r = r*sclfac
         ra = ra*sclfac
         GO TO 180
*
*        Now balance.
*
  190    CONTINUE
         IF( ( c+r ).GE.factor*s )
     $      GO TO 200
         IF( f.LT.one .AND. scale( i ).LT.one ) THEN
            IF( f*scale( i ).LE.sfmin1 )
     $         GO TO 200
         END IF
         IF( f.GT.one .AND. scale( i ).GT.one ) THEN
            IF( scale( i ).GE.sfmax1 / f )
     $         GO TO 200
         END IF
         g = one / f
         scale( i ) = scale( i )*f
         noconv = .true.
*
         CALL zdscal( n-k+1, g, a( i, k ), lda )
         CALL zdscal( l, f, a( 1, i ), 1 )
*
  200 CONTINUE
*
      IF( noconv )
     $   GO TO 140
*
  210 CONTINUE
      ilo = k
      ihi = l
*
      RETURN
*
*     End of ZGEBAL
*
      END

! ZGEHRD
      SUBROUTINE zgehrd( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            IHI, ILO, INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16        A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NBMAX, LDT, TSIZE
      parameter( nbmax = 64, ldt = nbmax+1,
     $                     tsize = ldt*nbmax )
      COMPLEX*16        ZERO, ONE
      parameter( zero = ( 0.0d+0, 0.0d+0 ),
     $                     one = ( 1.0d+0, 0.0d+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWT, J, LDWORK, LWKOPT, NB,
     $                   NBMIN, NH, NX
      COMPLEX*16        EI
*     ..
*     .. External Subroutines ..
      EXTERNAL           zaxpy, zgehd2, zgemm, zlahr2, zlarfb, ztrmm,
     $                   xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          max, min
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ilaenv
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      info = 0
      lquery = ( lwork.EQ.-1 )
      IF( n.LT.0 ) THEN
         info = -1
      ELSE IF( ilo.LT.1 .OR. ilo.GT.max( 1, n ) ) THEN
         info = -2
      ELSE IF( ihi.LT.min( ilo, n ) .OR. ihi.GT.n ) THEN
         info = -3
      ELSE IF( lda.LT.max( 1, n ) ) THEN
         info = -5
      ELSE IF( lwork.LT.max( 1, n ) .AND. .NOT.lquery ) THEN
         info = -8
      END IF
*
      IF( info.EQ.0 ) THEN
*
*        Compute the workspace requirements
*
         nb = min( nbmax, ilaenv( 1, 'ZGEHRD', ' ', n, ilo, ihi, -1 ) )
         lwkopt = n*nb + tsize
         work( 1 ) = lwkopt
      ENDIF
*
      IF( info.NE.0 ) THEN
         CALL xerbla( 'ZGEHRD', -info )
         RETURN
      ELSE IF( lquery ) THEN
         RETURN
      END IF
*
*     Set elements 1:ILO-1 and IHI:N-1 of TAU to zero
*
      DO 10 i = 1, ilo - 1
         tau( i ) = zero
   10 CONTINUE
      DO 20 i = max( 1, ihi ), n - 1
         tau( i ) = zero
   20 CONTINUE
*
*     Quick return if possible
*
      nh = ihi - ilo + 1
      IF( nh.LE.1 ) THEN
         work( 1 ) = 1
         RETURN
      END IF
*
*     Determine the block size
*
      nb = min( nbmax, ilaenv( 1, 'ZGEHRD', ' ', n, ilo, ihi, -1 ) )
      nbmin = 2
      IF( nb.GT.1 .AND. nb.LT.nh ) THEN
*
*        Determine when to cross over from blocked to unblocked code
*        (last block is always handled by unblocked code)
*
         nx = max( nb, ilaenv( 3, 'ZGEHRD', ' ', n, ilo, ihi, -1 ) )
         IF( nx.LT.nh ) THEN
*
*           Determine if workspace is large enough for blocked code
*
            IF( lwork.LT.n*nb+tsize ) THEN
*
*              Not enough workspace to use optimal NB:  determine the
*              minimum value of NB, and reduce NB or force use of
*              unblocked code
*
               nbmin = max( 2, ilaenv( 2, 'ZGEHRD', ' ', n, ilo, ihi,
     $                 -1 ) )
               IF( lwork.GE.(n*nbmin + tsize) ) THEN
                  nb = (lwork-tsize) / n
               ELSE
                  nb = 1
               END IF
            END IF
         END IF
      END IF
      ldwork = n
*
      IF( nb.LT.nbmin .OR. nb.GE.nh ) THEN
*
*        Use unblocked code below
*
         i = ilo
*
      ELSE
*
*        Use blocked code
*
         iwt = 1 + n*nb
         DO 40 i = ilo, ihi - 1 - nx, nb
            ib = min( nb, ihi-i )
*
*           Reduce columns i:i+ib-1 to Hessenberg form, returning the
*           matrices V and T of the block reflector H = I - V*T*V**H
*           which performs the reduction, and also the matrix Y = A*V*T
*
            CALL zlahr2( ihi, i, ib, a( 1, i ), lda, tau( i ),
     $                   work( iwt ), ldt, work, ldwork )
*
*           Apply the block reflector H to A(1:ihi,i+ib:ihi) from the
*           right, computing  A := A - Y * V**H. V(i+ib,ib-1) must be set
*           to 1
*
            ei = a( i+ib, i+ib-1 )
            a( i+ib, i+ib-1 ) = one
            CALL zgemm( 'No transpose', 'Conjugate transpose',
     $                  ihi, ihi-i-ib+1,
     $                  ib, -one, work, ldwork, a( i+ib, i ), lda, one,
     $                  a( 1, i+ib ), lda )
            a( i+ib, i+ib-1 ) = ei
*
*           Apply the block reflector H to A(1:i,i+1:i+ib-1) from the
*           right
*
            CALL ztrmm( 'Right', 'Lower', 'Conjugate transpose',
     $                  'Unit', i, ib-1,
     $                  one, a( i+1, i ), lda, work, ldwork )
            DO 30 j = 0, ib-2
               CALL zaxpy( i, -one, work( ldwork*j+1 ), 1,
     $                     a( 1, i+j+1 ), 1 )
   30       CONTINUE
*
*           Apply the block reflector H to A(i+1:ihi,i+ib:n) from the
*           left
*
            CALL zlarfb( 'Left', 'Conjugate transpose', 'Forward',
     $                   'Columnwise',
     $                   ihi-i, n-i-ib+1, ib, a( i+1, i ), lda,
     $                   work( iwt ), ldt, a( i+1, i+ib ), lda,
     $                   work, ldwork )
   40    CONTINUE
      END IF
*
*     Use unblocked code to reduce the rest of the matrix
*
      CALL zgehd2( n, i, ihi, a, lda, tau, work, iinfo )
      work( 1 ) = lwkopt
*
      RETURN
*
*     End of ZGEHRD
*
      END
      
! ZLACPY
      SUBROUTINE zlacpy( UPLO, M, N, A, LDA, B, LDB )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, LDB, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, J
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           lsame
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          min
*     ..
*     .. Executable Statements ..
*
      IF( lsame( uplo, 'U' ) ) THEN
         DO 20 j = 1, n
            DO 10 i = 1, min( j, m )
               b( i, j ) = a( i, j )
   10       CONTINUE
   20    CONTINUE
*
      ELSE IF( lsame( uplo, 'L' ) ) THEN
         DO 40 j = 1, n
            DO 30 i = j, m
               b( i, j ) = a( i, j )
   30       CONTINUE
   40    CONTINUE
*
      ELSE
         DO 60 j = 1, n
            DO 50 i = 1, m
               b( i, j ) = a( i, j )
   50       CONTINUE
   60    CONTINUE
      END IF
*
      RETURN
*
*     End of ZLACPY
*
      END

! ZUNGHR
      SUBROUTINE zunghr( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            IHI, ILO, INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ZERO, ONE
      parameter( zero = ( 0.0d+0, 0.0d+0 ),
     $                   one = ( 1.0d+0, 0.0d+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IINFO, J, LWKOPT, NB, NH
*     ..
*     .. External Subroutines ..
      EXTERNAL           xerbla, zungqr
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ilaenv
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          max, min
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      info = 0
      nh = ihi - ilo
      lquery = ( lwork.EQ.-1 )
      IF( n.LT.0 ) THEN
         info = -1
      ELSE IF( ilo.LT.1 .OR. ilo.GT.max( 1, n ) ) THEN
         info = -2
      ELSE IF( ihi.LT.min( ilo, n ) .OR. ihi.GT.n ) THEN
         info = -3
      ELSE IF( lda.LT.max( 1, n ) ) THEN
         info = -5
      ELSE IF( lwork.LT.max( 1, nh ) .AND. .NOT.lquery ) THEN
         info = -8
      END IF
*
      IF( info.EQ.0 ) THEN
         nb = ilaenv( 1, 'ZUNGQR', ' ', nh, nh, nh, -1 )
         lwkopt = max( 1, nh )*nb
         work( 1 ) = lwkopt
      END IF
*
      IF( info.NE.0 ) THEN
         CALL xerbla( 'ZUNGHR', -info )
         RETURN
      ELSE IF( lquery ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( n.EQ.0 ) THEN
         work( 1 ) = 1
         RETURN
      END IF
*
*     Shift the vectors which define the elementary reflectors one
*     column to the right, and set the first ilo and the last n-ihi
*     rows and columns to those of the unit matrix
*
      DO 40 j = ihi, ilo + 1, -1
         DO 10 i = 1, j - 1
            a( i, j ) = zero
   10    CONTINUE
         DO 20 i = j + 1, ihi
            a( i, j ) = a( i, j-1 )
   20    CONTINUE
         DO 30 i = ihi + 1, n
            a( i, j ) = zero
   30    CONTINUE
   40 CONTINUE
      DO 60 j = 1, ilo
         DO 50 i = 1, n
            a( i, j ) = zero
   50    CONTINUE
         a( j, j ) = one
   60 CONTINUE
      DO 80 j = ihi + 1, n
         DO 70 i = 1, n
            a( i, j ) = zero
   70    CONTINUE
         a( j, j ) = one
   80 CONTINUE
*
      IF( nh.GT.0 ) THEN
*
*        Generate Q(ilo+1:ihi,ilo+1:ihi)
*
         CALL zungqr( nh, nh, nh, a( ilo+1, ilo+1 ), lda, tau( ilo ),
     $                work, lwork, iinfo )
      END IF
      work( 1 ) = lwkopt
      RETURN
*
*     End of ZUNGHR
*
      END
      
! ZGEBAK
      SUBROUTINE zgebak( JOB, SIDE, N, ILO, IHI, SCALE, M, V, LDV,
     $                   INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          JOB, SIDE
      INTEGER            IHI, ILO, INFO, LDV, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   SCALE( * )
      COMPLEX*16         V( LDV, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      parameter( one = 1.0d+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LEFTV, RIGHTV
      INTEGER            I, II, K
      DOUBLE PRECISION   S
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           lsame
*     ..
*     .. External Subroutines ..
      EXTERNAL           xerbla, zdscal, zswap
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          max, min
*     ..
*     .. Executable Statements ..
*
*     Decode and Test the input parameters
*
      rightv = lsame( side, 'R' )
      leftv = lsame( side, 'L' )
*
      info = 0
      IF( .NOT.lsame( job, 'N' ) .AND. .NOT.lsame( job, 'P' ) .AND.
     $    .NOT.lsame( job, 'S' ) .AND. .NOT.lsame( job, 'B' ) ) THEN
         info = -1
      ELSE IF( .NOT.rightv .AND. .NOT.leftv ) THEN
         info = -2
      ELSE IF( n.LT.0 ) THEN
         info = -3
      ELSE IF( ilo.LT.1 .OR. ilo.GT.max( 1, n ) ) THEN
         info = -4
      ELSE IF( ihi.LT.min( ilo, n ) .OR. ihi.GT.n ) THEN
         info = -5
      ELSE IF( m.LT.0 ) THEN
         info = -7
      ELSE IF( ldv.LT.max( 1, n ) ) THEN
         info = -9
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'ZGEBAK', -info )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( n.EQ.0 )
     $   RETURN
      IF( m.EQ.0 )
     $   RETURN
      IF( lsame( job, 'N' ) )
     $   RETURN
*
      IF( ilo.EQ.ihi )
     $   GO TO 30
*
*     Backward balance
*
      IF( lsame( job, 'S' ) .OR. lsame( job, 'B' ) ) THEN
*
         IF( rightv ) THEN
            DO 10 i = ilo, ihi
               s = scale( i )
               CALL zdscal( m, s, v( i, 1 ), ldv )
   10       CONTINUE
         END IF
*
         IF( leftv ) THEN
            DO 20 i = ilo, ihi
               s = one / scale( i )
               CALL zdscal( m, s, v( i, 1 ), ldv )
   20       CONTINUE
         END IF
*
      END IF
*
*     Backward permutation
*
*     For  I = ILO-1 step -1 until 1,
*              IHI+1 step 1 until N do --
*
   30 CONTINUE
      IF( lsame( job, 'P' ) .OR. lsame( job, 'B' ) ) THEN
         IF( rightv ) THEN
            DO 40 ii = 1, n
               i = ii
               IF( i.GE.ilo .AND. i.LE.ihi )
     $            GO TO 40
               IF( i.LT.ilo )
     $            i = ilo - ii
               k = int( scale( i ) )
               IF( k.EQ.i )
     $            GO TO 40
               CALL zswap( m, v( i, 1 ), ldv, v( k, 1 ), ldv )
   40       CONTINUE
         END IF
*
         IF( leftv ) THEN
            DO 50 ii = 1, n
               i = ii
               IF( i.GE.ilo .AND. i.LE.ihi )
     $            GO TO 50
               IF( i.LT.ilo )
     $            i = ilo - ii
               k = int( scale( i ) )
               IF( k.EQ.i )
     $            GO TO 50
               CALL zswap( m, v( i, 1 ), ldv, v( k, 1 ), ldv )
   50       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of ZGEBAK
*
      END

! DZNRM2      
      DOUBLE PRECISION FUNCTION DZNRM2(N,X,INCX)
*     .. Scalar Arguments ..
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      DOUBLE COMPLEX X(*)
*     ..
*
*  Purpose
*  =======
*
*  DZNRM2 returns the euclidean norm of a vector via the function
*  name, so that
*
*     DZNRM2 := sqrt( conjg( x' )*x )
*
*
*  -- This version written on 25-October-1982.
*     Modified on 14-October-1993 to inline the call to ZLASSQ.
*     Sven Hammarling, Nag Ltd.
*
*
*     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION NORM,SCALE,SSQ,TEMP
      INTEGER IX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC ABS,DBLE,DIMAG,SQRT
*     ..
      IF (N.LT.1 .OR. INCX.LT.1) THEN
          NORM = ZERO
      ELSE
          SCALE = ZERO
          SSQ = ONE
*        The following loop is equivalent to this call to the LAPACK
*        auxiliary routine:
*        CALL ZLASSQ( N, X, INCX, SCALE, SSQ )
*
          DO 10 IX = 1,1 + (N-1)*INCX,INCX
              IF (DBLE(X(IX)).NE.ZERO) THEN
                  TEMP = ABS(DBLE(X(IX)))
                  IF (SCALE.LT.TEMP) THEN
                      SSQ = ONE + SSQ* (SCALE/TEMP)**2
                      SCALE = TEMP
                  ELSE
                      SSQ = SSQ + (TEMP/SCALE)**2
                  END IF
              END IF
              IF (DIMAG(X(IX)).NE.ZERO) THEN
                  TEMP = ABS(DIMAG(X(IX)))
                  IF (SCALE.LT.TEMP) THEN
                      SSQ = ONE + SSQ* (SCALE/TEMP)**2
                      SCALE = TEMP
                  ELSE
                      SSQ = SSQ + (TEMP/SCALE)**2
                  END IF
              END IF
   10     CONTINUE
          NORM = SCALE*SQRT(SSQ)
      END IF
*
      DZNRM2 = NORM
      RETURN
*
*     End of DZNRM2.
*
      END

! ZDSCAL
      SUBROUTINE zdscal(N,DA,ZX,INCX)
*
*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION DA
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      COMPLEX*16 ZX(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER I,NINCX
*     .. Parameters ..
      DOUBLE PRECISION ONE
      parameter(one=1.0d+0)
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC dble, dcmplx, dimag
*     ..
      IF (n.LE.0 .OR. incx.LE.0 .OR. da.EQ.one) RETURN
      IF (incx.EQ.1) THEN
*
*        code for increment equal to 1
*
         DO i = 1,n
            zx(i) = dcmplx(da*dble(zx(i)),da*dimag(zx(i)))
         END DO
      ELSE
*
*        code for increment not equal to 1
*
         nincx = n*incx
         DO i = 1,nincx,incx
            zx(i) = dcmplx(da*dble(zx(i)),da*dimag(zx(i)))
         END DO
      END IF
      RETURN
*
*     End of ZDSCAL
*
      END
      
! ZSCAL
      SUBROUTINE zscal(N,ZA,ZX,INCX)
*
*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      COMPLEX*16 ZA
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      COMPLEX*16 ZX(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER I,NINCX
*     ..
*     .. Parameters ..
      COMPLEX*16 ONE
      parameter(one= (1.0d+0,0.0d+0))
*     ..
      IF (n.LE.0 .OR. incx.LE.0 .OR. za.EQ.one) RETURN
      IF (incx.EQ.1) THEN
*
*        code for increment equal to 1
*
         DO i = 1,n
            zx(i) = za*zx(i)
         END DO
      ELSE
*
*        code for increment not equal to 1
*
         nincx = n*incx
         DO i = 1,nincx,incx
            zx(i) = za*zx(i)
         END DO
      END IF
      RETURN
*
*     End of ZSCAL
*
      END

! ZGETRF2
      RECURSIVE SUBROUTINE zgetrf2( M, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            info, lda, m, n
*     ..
*     .. Array Arguments ..
      INTEGER            ipiv( * )
      COMPLEX*16         a( lda, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         one, zero
      parameter( one = ( 1.0d+0, 0.0d+0 ),
     $                     zero = ( 0.0d+0, 0.0d+0 ) )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   sfmin
      COMPLEX*16         temp
      INTEGER            i, iinfo, n1, n2
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   dlamch
      INTEGER            izamax
      EXTERNAL           dlamch, izamax
*     ..
*     .. External Subroutines ..
      EXTERNAL           zgemm, zscal, zlaswp, ztrsm, xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          max, min
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      info = 0
      IF( m.LT.0 ) THEN
         info = -1
      ELSE IF( n.LT.0 ) THEN
         info = -2
      ELSE IF( lda.LT.max( 1, m ) ) THEN
         info = -4
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'ZGETRF2', -info )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( m.EQ.0 .OR. n.EQ.0 )
     $   RETURN
 
      IF ( m.EQ.1 ) THEN
*
*        Use unblocked code for one row case
*        Just need to handle IPIV and INFO
*
         ipiv( 1 ) = 1
         IF ( a(1,1).EQ.zero )
     $      info = 1
*
      ELSE IF( n.EQ.1 ) THEN
*
*        Use unblocked code for one column case
*
*
*        Compute machine safe minimum
*
         sfmin = dlamch('S')
*
*        Find pivot and test for singularity
*
         i = izamax( m, a( 1, 1 ), 1 )
         ipiv( 1 ) = i
         IF( a( i, 1 ).NE.zero ) THEN
*
*           Apply the interchange
*
            IF( i.NE.1 ) THEN
               temp = a( 1, 1 )
               a( 1, 1 ) = a( i, 1 )
               a( i, 1 ) = temp
            END IF
*
*           Compute elements 2:M of the column
*
            IF( abs(a( 1, 1 )) .GE. sfmin ) THEN
               CALL zscal( m-1, one / a( 1, 1 ), a( 2, 1 ), 1 )
            ELSE
               DO 10 i = 1, m-1
                  a( 1+i, 1 ) = a( 1+i, 1 ) / a( 1, 1 )
   10          CONTINUE
            END IF
*
         ELSE
            info = 1
         END IF
 
      ELSE
*
*        Use recursive code
*
         n1 = min( m, n ) / 2
         n2 = n-n1
*
*               [ A11 ]
*        Factor [ --- ]
*               [ A21 ]
*
         CALL zgetrf2( m, n1, a, lda, ipiv, iinfo )
 
         IF ( info.EQ.0 .AND. iinfo.GT.0 )
     $      info = iinfo
*
*                              [ A12 ]
*        Apply interchanges to [ --- ]
*                              [ A22 ]
*
         CALL zlaswp( n2, a( 1, n1+1 ), lda, 1, n1, ipiv, 1 )
*
*        Solve A12
*
         CALL ztrsm( 'L', 'L', 'N', 'U', n1, n2, one, a, lda,
     $               a( 1, n1+1 ), lda )
*
*        Update A22
*
         CALL zgemm( 'N', 'N', m-n1, n2, n1, -one, a( n1+1, 1 ), lda,
     $               a( 1, n1+1 ), lda, one, a( n1+1, n1+1 ), lda )
*
*        Factor A22
*
         CALL zgetrf2( m-n1, n2, a( n1+1, n1+1 ), lda, ipiv( n1+1 ),
     $                 iinfo )
*
*        Adjust INFO and the pivot indices
*
         IF ( info.EQ.0 .AND. iinfo.GT.0 )
     $      info = iinfo + n1
         DO 20 i = n1+1, min( m, n )
            ipiv( i ) = ipiv( i ) + n1
   20    CONTINUE
*
*        Apply interchanges to A21
*
         CALL zlaswp( n1, a( 1, 1 ), lda, n1+1, min( m, n), ipiv, 1 )
*
      END IF
      RETURN
*
*     End of ZGETRF2
*
      END
      
! ZLASWP
      SUBROUTINE zlaswp( N, A, LDA, K1, K2, IPIV, INCX )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INCX, K1, K2, LDA, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * )
*     ..
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, I1, I2, INC, IP, IX, IX0, J, K, N32
      COMPLEX*16         TEMP
*     ..
*     .. Executable Statements ..
*
*     Interchange row I with row IPIV(K1+(I-K1)*abs(INCX)) for each of rows
*     K1 through K2.
*
      IF( incx.GT.0 ) THEN
         ix0 = k1
         i1 = k1
         i2 = k2
         inc = 1
      ELSE IF( incx.LT.0 ) THEN
         ix0 = k1 + ( k1-k2 )*incx
         i1 = k2
         i2 = k1
         inc = -1
      ELSE
         RETURN
      END IF
*
      n32 = ( n / 32 )*32
      IF( n32.NE.0 ) THEN
         DO 30 j = 1, n32, 32
            ix = ix0
            DO 20 i = i1, i2, inc
               ip = ipiv( ix )
               IF( ip.NE.i ) THEN
                  DO 10 k = j, j + 31
                     temp = a( i, k )
                     a( i, k ) = a( ip, k )
                     a( ip, k ) = temp
   10             CONTINUE
               END IF
               ix = ix + incx
   20       CONTINUE
   30    CONTINUE
      END IF
      IF( n32.NE.n ) THEN
         n32 = n32 + 1
         ix = ix0
         DO 50 i = i1, i2, inc
            ip = ipiv( ix )
            IF( ip.NE.i ) THEN
               DO 40 k = n32, n
                  temp = a( i, k )
                  a( i, k ) = a( ip, k )
                  a( ip, k ) = temp
   40          CONTINUE
            END IF
            ix = ix + incx
   50    CONTINUE
      END IF
*
      RETURN
*
*     End of ZLASWP
*
      END

! ZTRSM
      SUBROUTINE ztrsm(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
*
*  -- Reference BLAS level3 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      COMPLEX*16 ALPHA
      INTEGER LDA,LDB,M,N
      CHARACTER DIAG,SIDE,TRANSA,UPLO
*     ..
*     .. Array Arguments ..
      COMPLEX*16 A(LDA,*),B(LDB,*)
*     ..
*
*  =====================================================================
*
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL lsame
*     ..
*     .. External Subroutines ..
      EXTERNAL xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC dconjg,max
*     ..
*     .. Local Scalars ..
      COMPLEX*16 TEMP
      INTEGER I,INFO,J,K,NROWA
      LOGICAL LSIDE,NOCONJ,NOUNIT,UPPER
*     ..
*     .. Parameters ..
      COMPLEX*16 ONE
      parameter(one= (1.0d+0,0.0d+0))
      COMPLEX*16 ZERO
      parameter(zero= (0.0d+0,0.0d+0))
*     ..
*
*     Test the input parameters.
*
      lside = lsame(side,'L')
      IF (lside) THEN
          nrowa = m
      ELSE
          nrowa = n
      END IF
      noconj = lsame(transa,'T')
      nounit = lsame(diag,'N')
      upper = lsame(uplo,'U')
*
      info = 0
      IF ((.NOT.lside) .AND. (.NOT.lsame(side,'R'))) THEN
          info = 1
      ELSE IF ((.NOT.upper) .AND. (.NOT.lsame(uplo,'L'))) THEN
          info = 2
      ELSE IF ((.NOT.lsame(transa,'N')) .AND.
     +         (.NOT.lsame(transa,'T')) .AND.
     +         (.NOT.lsame(transa,'C'))) THEN
          info = 3
      ELSE IF ((.NOT.lsame(diag,'U')) .AND. (.NOT.lsame(diag,'N'))) THEN
          info = 4
      ELSE IF (m.LT.0) THEN
          info = 5
      ELSE IF (n.LT.0) THEN
          info = 6
      ELSE IF (lda.LT.max(1,nrowa)) THEN
          info = 9
      ELSE IF (ldb.LT.max(1,m)) THEN
          info = 11
      END IF
      IF (info.NE.0) THEN
          CALL xerbla('ZTRSM ',info)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF (m.EQ.0 .OR. n.EQ.0) RETURN
*
*     And when  alpha.eq.zero.
*
      IF (alpha.EQ.zero) THEN
          DO 20 j = 1,n
              DO 10 i = 1,m
                  b(i,j) = zero
   10         CONTINUE
   20     CONTINUE
          RETURN
      END IF
*
*     Start the operations.
*
      IF (lside) THEN
          IF (lsame(transa,'N')) THEN
*
*           Form  B := alpha*inv( A )*B.
*
              IF (upper) THEN
                  DO 60 j = 1,n
                      IF (alpha.NE.one) THEN
                          DO 30 i = 1,m
                              b(i,j) = alpha*b(i,j)
   30                     CONTINUE
                      END IF
                      DO 50 k = m,1,-1
                          IF (b(k,j).NE.zero) THEN
                              IF (nounit) b(k,j) = b(k,j)/a(k,k)
                              DO 40 i = 1,k - 1
                                  b(i,j) = b(i,j) - b(k,j)*a(i,k)
   40                         CONTINUE
                          END IF
   50                 CONTINUE
   60             CONTINUE
              ELSE
                  DO 100 j = 1,n
                      IF (alpha.NE.one) THEN
                          DO 70 i = 1,m
                              b(i,j) = alpha*b(i,j)
   70                     CONTINUE
                      END IF
                      DO 90 k = 1,m
                          IF (b(k,j).NE.zero) THEN
                              IF (nounit) b(k,j) = b(k,j)/a(k,k)
                              DO 80 i = k + 1,m
                                  b(i,j) = b(i,j) - b(k,j)*a(i,k)
   80                         CONTINUE
                          END IF
   90                 CONTINUE
  100             CONTINUE
              END IF
          ELSE
*
*           Form  B := alpha*inv( A**T )*B
*           or    B := alpha*inv( A**H )*B.
*
              IF (upper) THEN
                  DO 140 j = 1,n
                      DO 130 i = 1,m
                          temp = alpha*b(i,j)
                          IF (noconj) THEN
                              DO 110 k = 1,i - 1
                                  temp = temp - a(k,i)*b(k,j)
  110                         CONTINUE
                              IF (nounit) temp = temp/a(i,i)
                          ELSE
                              DO 120 k = 1,i - 1
                                  temp = temp - dconjg(a(k,i))*b(k,j)
  120                         CONTINUE
                              IF (nounit) temp = temp/dconjg(a(i,i))
                          END IF
                          b(i,j) = temp
  130                 CONTINUE
  140             CONTINUE
              ELSE
                  DO 180 j = 1,n
                      DO 170 i = m,1,-1
                          temp = alpha*b(i,j)
                          IF (noconj) THEN
                              DO 150 k = i + 1,m
                                  temp = temp - a(k,i)*b(k,j)
  150                         CONTINUE
                              IF (nounit) temp = temp/a(i,i)
                          ELSE
                              DO 160 k = i + 1,m
                                  temp = temp - dconjg(a(k,i))*b(k,j)
  160                         CONTINUE
                              IF (nounit) temp = temp/dconjg(a(i,i))
                          END IF
                          b(i,j) = temp
  170                 CONTINUE
  180             CONTINUE
              END IF
          END IF
      ELSE
          IF (lsame(transa,'N')) THEN
*
*           Form  B := alpha*B*inv( A ).
*
              IF (upper) THEN
                  DO 230 j = 1,n
                      IF (alpha.NE.one) THEN
                          DO 190 i = 1,m
                              b(i,j) = alpha*b(i,j)
  190                     CONTINUE
                      END IF
                      DO 210 k = 1,j - 1
                          IF (a(k,j).NE.zero) THEN
                              DO 200 i = 1,m
                                  b(i,j) = b(i,j) - a(k,j)*b(i,k)
  200                         CONTINUE
                          END IF
  210                 CONTINUE
                      IF (nounit) THEN
                          temp = one/a(j,j)
                          DO 220 i = 1,m
                              b(i,j) = temp*b(i,j)
  220                     CONTINUE
                      END IF
  230             CONTINUE
              ELSE
                  DO 280 j = n,1,-1
                      IF (alpha.NE.one) THEN
                          DO 240 i = 1,m
                              b(i,j) = alpha*b(i,j)
  240                     CONTINUE
                      END IF
                      DO 260 k = j + 1,n
                          IF (a(k,j).NE.zero) THEN
                              DO 250 i = 1,m
                                  b(i,j) = b(i,j) - a(k,j)*b(i,k)
  250                         CONTINUE
                          END IF
  260                 CONTINUE
                      IF (nounit) THEN
                          temp = one/a(j,j)
                          DO 270 i = 1,m
                              b(i,j) = temp*b(i,j)
  270                     CONTINUE
                      END IF
  280             CONTINUE
              END IF
          ELSE
*
*           Form  B := alpha*B*inv( A**T )
*           or    B := alpha*B*inv( A**H ).
*
              IF (upper) THEN
                  DO 330 k = n,1,-1
                      IF (nounit) THEN
                          IF (noconj) THEN
                              temp = one/a(k,k)
                          ELSE
                              temp = one/dconjg(a(k,k))
                          END IF
                          DO 290 i = 1,m
                              b(i,k) = temp*b(i,k)
  290                     CONTINUE
                      END IF
                      DO 310 j = 1,k - 1
                          IF (a(j,k).NE.zero) THEN
                              IF (noconj) THEN
                                  temp = a(j,k)
                              ELSE
                                  temp = dconjg(a(j,k))
                              END IF
                              DO 300 i = 1,m
                                  b(i,j) = b(i,j) - temp*b(i,k)
  300                         CONTINUE
                          END IF
  310                 CONTINUE
                      IF (alpha.NE.one) THEN
                          DO 320 i = 1,m
                              b(i,k) = alpha*b(i,k)
  320                     CONTINUE
                      END IF
  330             CONTINUE
              ELSE
                  DO 380 k = 1,n
                      IF (nounit) THEN
                          IF (noconj) THEN
                              temp = one/a(k,k)
                          ELSE
                              temp = one/dconjg(a(k,k))
                          END IF
                          DO 340 i = 1,m
                              b(i,k) = temp*b(i,k)
  340                     CONTINUE
                      END IF
                      DO 360 j = k + 1,n
                          IF (a(j,k).NE.zero) THEN
                              IF (noconj) THEN
                                  temp = a(j,k)
                              ELSE
                                  temp = dconjg(a(j,k))
                              END IF
                              DO 350 i = 1,m
                                  b(i,j) = b(i,j) - temp*b(i,k)
  350                         CONTINUE
                          END IF
  360                 CONTINUE
                      IF (alpha.NE.one) THEN
                          DO 370 i = 1,m
                              b(i,k) = alpha*b(i,k)
  370                     CONTINUE
                      END IF
  380             CONTINUE
              END IF
          END IF
      END IF
*
      RETURN
*
*     End of ZTRSM
*
      END
      
! ZGEMM
      SUBROUTINE zgemm(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
*
*  -- Reference BLAS level3 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      COMPLEX*16 ALPHA,BETA
      INTEGER K,LDA,LDB,LDC,M,N
      CHARACTER TRANSA,TRANSB
*     ..
*     .. Array Arguments ..
      COMPLEX*16 A(LDA,*),B(LDB,*),C(LDC,*)
*     ..
*
*  =====================================================================
*
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL lsame
*     ..
*     .. External Subroutines ..
      EXTERNAL xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC dconjg,max
*     ..
*     .. Local Scalars ..
      COMPLEX*16 TEMP
      INTEGER I,INFO,J,L,NROWA,NROWB
      LOGICAL CONJA,CONJB,NOTA,NOTB
*     ..
*     .. Parameters ..
      COMPLEX*16 ONE
      parameter(one= (1.0d+0,0.0d+0))
      COMPLEX*16 ZERO
      parameter(zero= (0.0d+0,0.0d+0))
*     ..
*
*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
*     conjugated or transposed, set  CONJA and CONJB  as true if  A  and
*     B  respectively are to be  transposed but  not conjugated  and set
*     NROWA and NROWB  as the number of rows  of  A  and  B  respectively.
*
      nota = lsame(transa,'N')
      notb = lsame(transb,'N')
      conja = lsame(transa,'C')
      conjb = lsame(transb,'C')
      IF (nota) THEN
          nrowa = m
      ELSE
          nrowa = k
      END IF
      IF (notb) THEN
          nrowb = k
      ELSE
          nrowb = n
      END IF
*
*     Test the input parameters.
*
      info = 0
      IF ((.NOT.nota) .AND. (.NOT.conja) .AND.
     +    (.NOT.lsame(transa,'T'))) THEN
          info = 1
      ELSE IF ((.NOT.notb) .AND. (.NOT.conjb) .AND.
     +         (.NOT.lsame(transb,'T'))) THEN
          info = 2
      ELSE IF (m.LT.0) THEN
          info = 3
      ELSE IF (n.LT.0) THEN
          info = 4
      ELSE IF (k.LT.0) THEN
          info = 5
      ELSE IF (lda.LT.max(1,nrowa)) THEN
          info = 8
      ELSE IF (ldb.LT.max(1,nrowb)) THEN
          info = 10
      ELSE IF (ldc.LT.max(1,m)) THEN
          info = 13
      END IF
      IF (info.NE.0) THEN
          CALL xerbla('ZGEMM ',info)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((m.EQ.0) .OR. (n.EQ.0) .OR.
     +    (((alpha.EQ.zero).OR. (k.EQ.0)).AND. (beta.EQ.one))) RETURN
*
*     And when  alpha.eq.zero.
*
      IF (alpha.EQ.zero) THEN
          IF (beta.EQ.zero) THEN
              DO 20 j = 1,n
                  DO 10 i = 1,m
                      c(i,j) = zero
   10             CONTINUE
   20         CONTINUE
          ELSE
              DO 40 j = 1,n
                  DO 30 i = 1,m
                      c(i,j) = beta*c(i,j)
   30             CONTINUE
   40         CONTINUE
          END IF
          RETURN
      END IF
*
*     Start the operations.
*
      IF (notb) THEN
          IF (nota) THEN
*
*           Form  C := alpha*A*B + beta*C.
*
              DO 90 j = 1,n
                  IF (beta.EQ.zero) THEN
                      DO 50 i = 1,m
                          c(i,j) = zero
   50                 CONTINUE
                  ELSE IF (beta.NE.one) THEN
                      DO 60 i = 1,m
                          c(i,j) = beta*c(i,j)
   60                 CONTINUE
                  END IF
                  DO 80 l = 1,k
                      temp = alpha*b(l,j)
                      DO 70 i = 1,m
                          c(i,j) = c(i,j) + temp*a(i,l)
   70                 CONTINUE
   80             CONTINUE
   90         CONTINUE
          ELSE IF (conja) THEN
*
*           Form  C := alpha*A**H*B + beta*C.
*
              DO 120 j = 1,n
                  DO 110 i = 1,m
                      temp = zero
                      DO 100 l = 1,k
                          temp = temp + dconjg(a(l,i))*b(l,j)
  100                 CONTINUE
                      IF (beta.EQ.zero) THEN
                          c(i,j) = alpha*temp
                      ELSE
                          c(i,j) = alpha*temp + beta*c(i,j)
                      END IF
  110             CONTINUE
  120         CONTINUE
          ELSE
*
*           Form  C := alpha*A**T*B + beta*C
*
              DO 150 j = 1,n
                  DO 140 i = 1,m
                      temp = zero
                      DO 130 l = 1,k
                          temp = temp + a(l,i)*b(l,j)
  130                 CONTINUE
                      IF (beta.EQ.zero) THEN
                          c(i,j) = alpha*temp
                      ELSE
                          c(i,j) = alpha*temp + beta*c(i,j)
                      END IF
  140             CONTINUE
  150         CONTINUE
          END IF
      ELSE IF (nota) THEN
          IF (conjb) THEN
*
*           Form  C := alpha*A*B**H + beta*C.
*
              DO 200 j = 1,n
                  IF (beta.EQ.zero) THEN
                      DO 160 i = 1,m
                          c(i,j) = zero
  160                 CONTINUE
                  ELSE IF (beta.NE.one) THEN
                      DO 170 i = 1,m
                          c(i,j) = beta*c(i,j)
  170                 CONTINUE
                  END IF
                  DO 190 l = 1,k
                      temp = alpha*dconjg(b(j,l))
                      DO 180 i = 1,m
                          c(i,j) = c(i,j) + temp*a(i,l)
  180                 CONTINUE
  190             CONTINUE
  200         CONTINUE
          ELSE
*
*           Form  C := alpha*A*B**T + beta*C
*
              DO 250 j = 1,n
                  IF (beta.EQ.zero) THEN
                      DO 210 i = 1,m
                          c(i,j) = zero
  210                 CONTINUE
                  ELSE IF (beta.NE.one) THEN
                      DO 220 i = 1,m
                          c(i,j) = beta*c(i,j)
  220                 CONTINUE
                  END IF
                  DO 240 l = 1,k
                      temp = alpha*b(j,l)
                      DO 230 i = 1,m
                          c(i,j) = c(i,j) + temp*a(i,l)
  230                 CONTINUE
  240             CONTINUE
  250         CONTINUE
          END IF
      ELSE IF (conja) THEN
          IF (conjb) THEN
*
*           Form  C := alpha*A**H*B**H + beta*C.
*
              DO 280 j = 1,n
                  DO 270 i = 1,m
                      temp = zero
                      DO 260 l = 1,k
                          temp = temp + dconjg(a(l,i))*dconjg(b(j,l))
  260                 CONTINUE
                      IF (beta.EQ.zero) THEN
                          c(i,j) = alpha*temp
                      ELSE
                          c(i,j) = alpha*temp + beta*c(i,j)
                      END IF
  270             CONTINUE
  280         CONTINUE
          ELSE
*
*           Form  C := alpha*A**H*B**T + beta*C
*
              DO 310 j = 1,n
                  DO 300 i = 1,m
                      temp = zero
                      DO 290 l = 1,k
                          temp = temp + dconjg(a(l,i))*b(j,l)
  290                 CONTINUE
                      IF (beta.EQ.zero) THEN
                          c(i,j) = alpha*temp
                      ELSE
                          c(i,j) = alpha*temp + beta*c(i,j)
                      END IF
  300             CONTINUE
  310         CONTINUE
          END IF
      ELSE
          IF (conjb) THEN
*
*           Form  C := alpha*A**T*B**H + beta*C
*
              DO 340 j = 1,n
                  DO 330 i = 1,m
                      temp = zero
                      DO 320 l = 1,k
                          temp = temp + a(l,i)*dconjg(b(j,l))
  320                 CONTINUE
                      IF (beta.EQ.zero) THEN
                          c(i,j) = alpha*temp
                      ELSE
                          c(i,j) = alpha*temp + beta*c(i,j)
                      END IF
  330             CONTINUE
  340         CONTINUE
          ELSE
*
*           Form  C := alpha*A**T*B**T + beta*C
*
              DO 370 j = 1,n
                  DO 360 i = 1,m
                      temp = zero
                      DO 350 l = 1,k
                          temp = temp + a(l,i)*b(j,l)
  350                 CONTINUE
                      IF (beta.EQ.zero) THEN
                          c(i,j) = alpha*temp
                      ELSE
                          c(i,j) = alpha*temp + beta*c(i,j)
                      END IF
  360             CONTINUE
  370         CONTINUE
          END IF
      END IF
*
      RETURN
*
*     End of ZGEMM
*
      END

! ZTRTRI
      SUBROUTINE ztrtri( UPLO, DIAG, N, A, LDA, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      parameter( one = ( 1.0d+0, 0.0d+0 ),
     $                   zero = ( 0.0d+0, 0.0d+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOUNIT, UPPER
      INTEGER            J, JB, NB, NN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           lsame, ilaenv
*     ..
*     .. External Subroutines ..
      EXTERNAL           xerbla, ztrmm, ztrsm, ztrti2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          max, min
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      info = 0
      upper = lsame( uplo, 'U' )
      nounit = lsame( diag, 'N' )
      IF( .NOT.upper .AND. .NOT.lsame( uplo, 'L' ) ) THEN
         info = -1
      ELSE IF( .NOT.nounit .AND. .NOT.lsame( diag, 'U' ) ) THEN
         info = -2
      ELSE IF( n.LT.0 ) THEN
         info = -3
      ELSE IF( lda.LT.max( 1, n ) ) THEN
         info = -5
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'ZTRTRI', -info )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( n.EQ.0 )
     $   RETURN
*
*     Check for singularity if non-unit.
*
      IF( nounit ) THEN
         DO 10 info = 1, n
            IF( a( info, info ).EQ.zero )
     $         RETURN
   10    CONTINUE
         info = 0
      END IF
*
*     Determine the block size for this environment.
*
      nb = ilaenv( 1, 'ZTRTRI', uplo // diag, n, -1, -1, -1 )
      IF( nb.LE.1 .OR. nb.GE.n ) THEN
*
*        Use unblocked code
*
         CALL ztrti2( uplo, diag, n, a, lda, info )
      ELSE
*
*        Use blocked code
*
         IF( upper ) THEN
*
*           Compute inverse of upper triangular matrix
*
            DO 20 j = 1, n, nb
               jb = min( nb, n-j+1 )
*
*              Compute rows 1:j-1 of current block column
*
               CALL ztrmm( 'Left', 'Upper', 'No transpose', diag, j-1,
     $                     jb, one, a, lda, a( 1, j ), lda )
               CALL ztrsm( 'Right', 'Upper', 'No transpose', diag, j-1,
     $                     jb, -one, a( j, j ), lda, a( 1, j ), lda )
*
*              Compute inverse of current diagonal block
*
               CALL ztrti2( 'Upper', diag, jb, a( j, j ), lda, info )
   20       CONTINUE
         ELSE
*
*           Compute inverse of lower triangular matrix
*
            nn = ( ( n-1 ) / nb )*nb + 1
            DO 30 j = nn, 1, -nb
               jb = min( nb, n-j+1 )
               IF( j+jb.LE.n ) THEN
*
*                 Compute rows j+jb:n of current block column
*
                  CALL ztrmm( 'Left', 'Lower', 'No transpose', diag,
     $                        n-j-jb+1, jb, one, a( j+jb, j+jb ), lda,
     $                        a( j+jb, j ), lda )
                  CALL ztrsm( 'Right', 'Lower', 'No transpose', diag,
     $                        n-j-jb+1, jb, -one, a( j, j ), lda,
     $                        a( j+jb, j ), lda )
               END IF
*
*              Compute inverse of current diagonal block
*
               CALL ztrti2( 'Lower', diag, jb, a( j, j ), lda, info )
   30       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of ZTRTRI
*
      END

! ZGEMV
      SUBROUTINE zgemv(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
*
*  -- Reference BLAS level2 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      COMPLEX*16 ALPHA,BETA
      INTEGER INCX,INCY,LDA,M,N
      CHARACTER TRANS
*     ..
*     .. Array Arguments ..
      COMPLEX*16 A(LDA,*),X(*),Y(*)
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16 ONE
      parameter(one= (1.0d+0,0.0d+0))
      COMPLEX*16 ZERO
      parameter(zero= (0.0d+0,0.0d+0))
*     ..
*     .. Local Scalars ..
      COMPLEX*16 TEMP
      INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY,LENX,LENY
      LOGICAL NOCONJ
*     ..
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL lsame
*     ..
*     .. External Subroutines ..
      EXTERNAL xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC dconjg,max
*     ..
*
*     Test the input parameters.
*
      info = 0
      IF (.NOT.lsame(trans,'N') .AND. .NOT.lsame(trans,'T') .AND.
     +    .NOT.lsame(trans,'C')) THEN
          info = 1
      ELSE IF (m.LT.0) THEN
          info = 2
      ELSE IF (n.LT.0) THEN
          info = 3
      ELSE IF (lda.LT.max(1,m)) THEN
          info = 6
      ELSE IF (incx.EQ.0) THEN
          info = 8
      ELSE IF (incy.EQ.0) THEN
          info = 11
      END IF
      IF (info.NE.0) THEN
          CALL xerbla('ZGEMV ',info)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((m.EQ.0) .OR. (n.EQ.0) .OR.
     +    ((alpha.EQ.zero).AND. (beta.EQ.one))) RETURN
*
      noconj = lsame(trans,'T')
*
*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
*     up the start points in  X  and  Y.
*
      IF (lsame(trans,'N')) THEN
          lenx = n
          leny = m
      ELSE
          lenx = m
          leny = n
      END IF
      IF (incx.GT.0) THEN
          kx = 1
      ELSE
          kx = 1 - (lenx-1)*incx
      END IF
      IF (incy.GT.0) THEN
          ky = 1
      ELSE
          ky = 1 - (leny-1)*incy
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
*     First form  y := beta*y.
*
      IF (beta.NE.one) THEN
          IF (incy.EQ.1) THEN
              IF (beta.EQ.zero) THEN
                  DO 10 i = 1,leny
                      y(i) = zero
   10             CONTINUE
              ELSE
                  DO 20 i = 1,leny
                      y(i) = beta*y(i)
   20             CONTINUE
              END IF
          ELSE
              iy = ky
              IF (beta.EQ.zero) THEN
                  DO 30 i = 1,leny
                      y(iy) = zero
                      iy = iy + incy
   30             CONTINUE
              ELSE
                  DO 40 i = 1,leny
                      y(iy) = beta*y(iy)
                      iy = iy + incy
   40             CONTINUE
              END IF
          END IF
      END IF
      IF (alpha.EQ.zero) RETURN
      IF (lsame(trans,'N')) THEN
*
*        Form  y := alpha*A*x + y.
*
          jx = kx
          IF (incy.EQ.1) THEN
              DO 60 j = 1,n
                  temp = alpha*x(jx)
                  DO 50 i = 1,m
                      y(i) = y(i) + temp*a(i,j)
   50             CONTINUE
                  jx = jx + incx
   60         CONTINUE
          ELSE
              DO 80 j = 1,n
                  temp = alpha*x(jx)
                  iy = ky
                  DO 70 i = 1,m
                      y(iy) = y(iy) + temp*a(i,j)
                      iy = iy + incy
   70             CONTINUE
                  jx = jx + incx
   80         CONTINUE
          END IF
      ELSE
*
*        Form  y := alpha*A**T*x + y  or  y := alpha*A**H*x + y.
*
          jy = ky
          IF (incx.EQ.1) THEN
              DO 110 j = 1,n
                  temp = zero
                  IF (noconj) THEN
                      DO 90 i = 1,m
                          temp = temp + a(i,j)*x(i)
   90                 CONTINUE
                  ELSE
                      DO 100 i = 1,m
                          temp = temp + dconjg(a(i,j))*x(i)
  100                 CONTINUE
                  END IF
                  y(jy) = y(jy) + alpha*temp
                  jy = jy + incy
  110         CONTINUE
          ELSE
              DO 140 j = 1,n
                  temp = zero
                  ix = kx
                  IF (noconj) THEN
                      DO 120 i = 1,m
                          temp = temp + a(i,j)*x(ix)
                          ix = ix + incx
  120                 CONTINUE
                  ELSE
                      DO 130 i = 1,m
                          temp = temp + dconjg(a(i,j))*x(ix)
                          ix = ix + incx
  130                 CONTINUE
                  END IF
                  y(jy) = y(jy) + alpha*temp
                  jy = jy + incy
  140         CONTINUE
          END IF
      END IF
*
      RETURN
*
*     End of ZGEMV
*
      END
      
! ZSWAP
      SUBROUTINE zswap(N,ZX,INCX,ZY,INCY)
*
*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      COMPLEX*16 ZX(*),ZY(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      COMPLEX*16 ZTEMP
      INTEGER I,IX,IY
*     ..
      IF (n.LE.0) RETURN
      IF (incx.EQ.1 .AND. incy.EQ.1) THEN
*
*       code for both increments equal to 1
         DO i = 1,n
            ztemp = zx(i)
            zx(i) = zy(i)
            zy(i) = ztemp
         END DO
      ELSE
*
*       code for unequal increments or equal increments not equal
*         to 1
*
         ix = 1
         iy = 1
         IF (incx.LT.0) ix = (-n+1)*incx + 1
         IF (incy.LT.0) iy = (-n+1)*incy + 1
         DO i = 1,n
            ztemp = zx(ix)
            zx(ix) = zy(iy)
            zy(iy) = ztemp
            ix = ix + incx
            iy = iy + incy
         END DO
      END IF
      RETURN
*
*     End of ZSWAP
*
      END

! ZLASET
      SUBROUTINE zlaset( UPLO, M, N, ALPHA, BETA, A, LDA )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, M, N
      COMPLEX*16         ALPHA, BETA
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, J
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           lsame
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          min
*     ..
*     .. Executable Statements ..
*
      IF( lsame( uplo, 'U' ) ) THEN
*
*        Set the diagonal to BETA and the strictly upper triangular
*        part of the array to ALPHA.
*
         DO 20 j = 2, n
            DO 10 i = 1, min( j-1, m )
               a( i, j ) = alpha
   10       CONTINUE
   20    CONTINUE
         DO 30 i = 1, min( n, m )
            a( i, i ) = beta
   30    CONTINUE
*
      ELSE IF( lsame( uplo, 'L' ) ) THEN
*
*        Set the diagonal to BETA and the strictly lower triangular
*        part of the array to ALPHA.
*
         DO 50 j = 1, min( m, n )
            DO 40 i = j + 1, m
               a( i, j ) = alpha
   40       CONTINUE
   50    CONTINUE
         DO 60 i = 1, min( n, m )
            a( i, i ) = beta
   60    CONTINUE
*
      ELSE
*
*        Set the array to BETA on the diagonal and ALPHA on the
*        offdiagonal.
*
         DO 80 j = 1, n
            DO 70 i = 1, m
               a( i, j ) = alpha
   70       CONTINUE
   80    CONTINUE
         DO 90 i = 1, min( m, n )
            a( i, i ) = beta
   90    CONTINUE
      END IF
*
      RETURN
*
*     End of ZLASET
*
      END
      
! DZASUM
      DOUBLE PRECISION FUNCTION dzasum(N,ZX,INCX)
*
*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER incx,n
*     ..
*     .. Array Arguments ..
      COMPLEX*16 zx(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION stemp
      INTEGER i,nincx
*     ..
*     .. External Functions ..
      DOUBLE PRECISION dcabs1
      EXTERNAL dcabs1
*     ..
      dzasum = 0.0d0
      stemp = 0.0d0
      IF (n.LE.0 .OR. incx.LE.0) RETURN
      IF (incx.EQ.1) THEN
*
*        code for increment equal to 1
*
         DO i = 1,n
            stemp = stemp + dcabs1(zx(i))
         END DO
      ELSE
*
*        code for increment not equal to 1
*
         nincx = n*incx
         DO i = 1,nincx,incx
            stemp = stemp + dcabs1(zx(i))
         END DO
      END IF
      dzasum = stemp
      RETURN
*
*     End of DZASUM
*
      END
      
! ZLATRS
      SUBROUTINE zlatrs( UPLO, TRANS, DIAG, NORMIN, N, A, LDA, X, SCALE,
     $                   CNORM, INFO )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, NORMIN, TRANS, UPLO
      INTEGER            INFO, LDA, N
      DOUBLE PRECISION   SCALE
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   CNORM( * )
      COMPLEX*16         A( LDA, * ), X( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, HALF, ONE, TWO
      parameter( zero = 0.0d+0, half = 0.5d+0, one = 1.0d+0,
     $                   two = 2.0d+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOTRAN, NOUNIT, UPPER
      INTEGER            I, IMAX, J, JFIRST, JINC, JLAST
      DOUBLE PRECISION   BIGNUM, GROW, REC, SMLNUM, TJJ, TMAX, TSCAL,
     $                   xbnd, xj, xmax
      COMPLEX*16         CSUMJ, TJJS, USCAL, ZDUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IDAMAX, IZAMAX
      DOUBLE PRECISION   DLAMCH, DZASUM
      COMPLEX*16         ZDOTC, ZDOTU, ZLADIV
      EXTERNAL           lsame, idamax, izamax, dlamch, dzasum, zdotc,
     $                   zdotu, zladiv
*     ..
*     .. External Subroutines ..
      EXTERNAL           dscal, xerbla, zaxpy, zdscal, ztrsv
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, dble, dcmplx, dconjg, dimag, max, min
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   CABS1, CABS2
*     ..
*     .. Statement Function definitions ..
      cabs1( zdum ) = abs( dble( zdum ) ) + abs( dimag( zdum ) )
      cabs2( zdum ) = abs( dble( zdum ) / 2.d0 ) +
     $                abs( dimag( zdum ) / 2.d0 )
*     ..
*     .. Executable Statements ..
*
      info = 0
      upper = lsame( uplo, 'U' )
      notran = lsame( trans, 'N' )
      nounit = lsame( diag, 'N' )
*
*     Test the input parameters.
*
      IF( .NOT.upper .AND. .NOT.lsame( uplo, 'L' ) ) THEN
         info = -1
      ELSE IF( .NOT.notran .AND. .NOT.lsame( trans, 'T' ) .AND. .NOT.
     $         lsame( trans, 'C' ) ) THEN
         info = -2
      ELSE IF( .NOT.nounit .AND. .NOT.lsame( diag, 'U' ) ) THEN
         info = -3
      ELSE IF( .NOT.lsame( normin, 'Y' ) .AND. .NOT.
     $         lsame( normin, 'N' ) ) THEN
         info = -4
      ELSE IF( n.LT.0 ) THEN
         info = -5
      ELSE IF( lda.LT.max( 1, n ) ) THEN
         info = -7
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'ZLATRS', -info )
         RETURN
      END IF
*
*     Quick return if possible
*
      scale = one
      IF( n.EQ.0 )
     $   RETURN
*
*     Determine machine dependent parameters to control overflow.
*
      smlnum = dlamch( 'Safe minimum' ) / dlamch( 'Precision' )
      bignum = one / smlnum
*
      IF( lsame( normin, 'N' ) ) THEN
*
*        Compute the 1-norm of each column, not including the diagonal.
*
         IF( upper ) THEN
*
*           A is upper triangular.
*
            DO 10 j = 1, n
               cnorm( j ) = dzasum( j-1, a( 1, j ), 1 )
   10       CONTINUE
         ELSE
*
*           A is lower triangular.
*
            DO 20 j = 1, n - 1
               cnorm( j ) = dzasum( n-j, a( j+1, j ), 1 )
   20       CONTINUE
            cnorm( n ) = zero
         END IF
      END IF
*
*     Scale the column norms by TSCAL if the maximum element in CNORM is
*     greater than BIGNUM/2.
*
      imax = idamax( n, cnorm, 1 )
      tmax = cnorm( imax )
      IF( tmax.LE.bignum*half ) THEN
         tscal = one
      ELSE
*
*        Avoid NaN generation if entries in CNORM exceed the
*        overflow threshold
*
         IF ( tmax.LE.dlamch('Overflow') ) THEN
*           Case 1: All entries in CNORM are valid floating-point numbers
            tscal = half / ( smlnum*tmax )
            CALL dscal( n, tscal, cnorm, 1 )
         ELSE
*           Case 2: At least one column norm of A cannot be
*           represented as a floating-point number. Find the
*           maximum offdiagonal absolute value
*           max( |Re(A(I,J))|, |Im(A(I,J)| ). If this entry is
*           not +/- Infinity, use this value as TSCAL.
            tmax = zero
            IF( upper ) THEN
*
*              A is upper triangular.
*
               DO j = 2, n
                  DO i = 1, j - 1
                     tmax = max( tmax, abs( dble( a( i, j ) ) ),
     $                           abs( dimag(a( i, j ) ) ) )
                  END DO
               END DO
            ELSE
*
*              A is lower triangular.
*
               DO j = 1, n - 1
                  DO i = j + 1, n
                     tmax = max( tmax, abs( dble( a( i, j ) ) ),
     $                           abs( dimag(a( i, j ) ) ) )
                  END DO
               END DO
            END IF
*
            IF( tmax.LE.dlamch('Overflow') ) THEN
               tscal = one / ( smlnum*tmax )
               DO j = 1, n
                  IF( cnorm( j ).LE.dlamch('Overflow') ) THEN
                     cnorm( j ) = cnorm( j )*tscal
                  ELSE
*                    Recompute the 1-norm of each column without
*                    introducing Infinity in the summation.
                     tscal = two * tscal
                     cnorm( j ) = zero
                     IF( upper ) THEN
                        DO i = 1, j - 1
                           cnorm( j ) = cnorm( j ) +
     $                                  tscal * cabs2( a( i, j ) )
                        END DO
                     ELSE
                        DO i = j + 1, n
                           cnorm( j ) = cnorm( j ) +
     $                                  tscal * cabs2( a( i, j ) )
                        END DO
                     END IF
                     tscal = tscal * half
                  END IF
               END DO
            ELSE
*              At least one entry of A is not a valid floating-point
*              entry. Rely on TRSV to propagate Inf and NaN.
               CALL ztrsv( uplo, trans, diag, n, a, lda, x, 1 )
               RETURN
            END IF
         END IF
      END IF
*
*     Compute a bound on the computed solution vector to see if the
*     Level 2 BLAS routine ZTRSV can be used.
*
      xmax = zero
      DO 30 j = 1, n
         xmax = max( xmax, cabs2( x( j ) ) )
   30 CONTINUE
      xbnd = xmax
*
      IF( notran ) THEN
*
*        Compute the growth in A * x = b.
*
         IF( upper ) THEN
            jfirst = n
            jlast = 1
            jinc = -1
         ELSE
            jfirst = 1
            jlast = n
            jinc = 1
         END IF
*
         IF( tscal.NE.one ) THEN
            grow = zero
            GO TO 60
         END IF
*
         IF( nounit ) THEN
*
*           A is non-unit triangular.
*
*           Compute GROW = 1/G(j) and XBND = 1/M(j).
*           Initially, G(0) = max{x(i), i=1,...,n}.
*
            grow = half / max( xbnd, smlnum )
            xbnd = grow
            DO 40 j = jfirst, jlast, jinc
*
*              Exit the loop if the growth factor is too small.
*
               IF( grow.LE.smlnum )
     $            GO TO 60
*
               tjjs = a( j, j )
               tjj = cabs1( tjjs )
*
               IF( tjj.GE.smlnum ) THEN
*
*                 M(j) = G(j-1) / abs(A(j,j))
*
                  xbnd = min( xbnd, min( one, tjj )*grow )
               ELSE
*
*                 M(j) could overflow, set XBND to 0.
*
                  xbnd = zero
               END IF
*
               IF( tjj+cnorm( j ).GE.smlnum ) THEN
*
*                 G(j) = G(j-1)*( 1 + CNORM(j) / abs(A(j,j)) )
*
                  grow = grow*( tjj / ( tjj+cnorm( j ) ) )
               ELSE
*
*                 G(j) could overflow, set GROW to 0.
*
                  grow = zero
               END IF
   40       CONTINUE
            grow = xbnd
         ELSE
*
*           A is unit triangular.
*
*           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
*
            grow = min( one, half / max( xbnd, smlnum ) )
            DO 50 j = jfirst, jlast, jinc
*
*              Exit the loop if the growth factor is too small.
*
               IF( grow.LE.smlnum )
     $            GO TO 60
*
*              G(j) = G(j-1)*( 1 + CNORM(j) )
*
               grow = grow*( one / ( one+cnorm( j ) ) )
   50       CONTINUE
         END IF
   60    CONTINUE
*
      ELSE
*
*        Compute the growth in A**T * x = b  or  A**H * x = b.
*
         IF( upper ) THEN
            jfirst = 1
            jlast = n
            jinc = 1
         ELSE
            jfirst = n
            jlast = 1
            jinc = -1
         END IF
*
         IF( tscal.NE.one ) THEN
            grow = zero
            GO TO 90
         END IF
*
         IF( nounit ) THEN
*
*           A is non-unit triangular.
*
*           Compute GROW = 1/G(j) and XBND = 1/M(j).
*           Initially, M(0) = max{x(i), i=1,...,n}.
*
            grow = half / max( xbnd, smlnum )
            xbnd = grow
            DO 70 j = jfirst, jlast, jinc
*
*              Exit the loop if the growth factor is too small.
*
               IF( grow.LE.smlnum )
     $            GO TO 90
*
*              G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) )
*
               xj = one + cnorm( j )
               grow = min( grow, xbnd / xj )
*
               tjjs = a( j, j )
               tjj = cabs1( tjjs )
*
               IF( tjj.GE.smlnum ) THEN
*
*                 M(j) = M(j-1)*( 1 + CNORM(j) ) / abs(A(j,j))
*
                  IF( xj.GT.tjj )
     $               xbnd = xbnd*( tjj / xj )
               ELSE
*
*                 M(j) could overflow, set XBND to 0.
*
                  xbnd = zero
               END IF
   70       CONTINUE
            grow = min( grow, xbnd )
         ELSE
*
*           A is unit triangular.
*
*           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
*
            grow = min( one, half / max( xbnd, smlnum ) )
            DO 80 j = jfirst, jlast, jinc
*
*              Exit the loop if the growth factor is too small.
*
               IF( grow.LE.smlnum )
     $            GO TO 90
*
*              G(j) = ( 1 + CNORM(j) )*G(j-1)
*
               xj = one + cnorm( j )
               grow = grow / xj
   80       CONTINUE
         END IF
   90    CONTINUE
      END IF
*
      IF( ( grow*tscal ).GT.smlnum ) THEN
*
*        Use the Level 2 BLAS solve if the reciprocal of the bound on
*        elements of X is not too small.
*
         CALL ztrsv( uplo, trans, diag, n, a, lda, x, 1 )
      ELSE
*
*        Use a Level 1 BLAS solve, scaling intermediate results.
*
         IF( xmax.GT.bignum*half ) THEN
*
*           Scale X so that its components are less than or equal to
*           BIGNUM in absolute value.
*
            scale = ( bignum*half ) / xmax
            CALL zdscal( n, scale, x, 1 )
            xmax = bignum
         ELSE
            xmax = xmax*two
         END IF
*
         IF( notran ) THEN
*
*           Solve A * x = b
*
            DO 120 j = jfirst, jlast, jinc
*
*              Compute x(j) = b(j) / A(j,j), scaling x if necessary.
*
               xj = cabs1( x( j ) )
               IF( nounit ) THEN
                  tjjs = a( j, j )*tscal
               ELSE
                  tjjs = tscal
                  IF( tscal.EQ.one )
     $               GO TO 110
               END IF
               tjj = cabs1( tjjs )
               IF( tjj.GT.smlnum ) THEN
*
*                    abs(A(j,j)) > SMLNUM:
*
                  IF( tjj.LT.one ) THEN
                     IF( xj.GT.tjj*bignum ) THEN
*
*                          Scale x by 1/b(j).
*
                        rec = one / xj
                        CALL zdscal( n, rec, x, 1 )
                        scale = scale*rec
                        xmax = xmax*rec
                     END IF
                  END IF
                  x( j ) = zladiv( x( j ), tjjs )
                  xj = cabs1( x( j ) )
               ELSE IF( tjj.GT.zero ) THEN
*
*                    0 < abs(A(j,j)) <= SMLNUM:
*
                  IF( xj.GT.tjj*bignum ) THEN
*
*                       Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM
*                       to avoid overflow when dividing by A(j,j).
*
                     rec = ( tjj*bignum ) / xj
                     IF( cnorm( j ).GT.one ) THEN
*
*                          Scale by 1/CNORM(j) to avoid overflow when
*                          multiplying x(j) times column j.
*
                        rec = rec / cnorm( j )
                     END IF
                     CALL zdscal( n, rec, x, 1 )
                     scale = scale*rec
                     xmax = xmax*rec
                  END IF
                  x( j ) = zladiv( x( j ), tjjs )
                  xj = cabs1( x( j ) )
               ELSE
*
*                    A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
*                    scale = 0, and compute a solution to A*x = 0.
*
                  DO 100 i = 1, n
                     x( i ) = zero
  100             CONTINUE
                  x( j ) = one
                  xj = one
                  scale = zero
                  xmax = zero
               END IF
  110          CONTINUE
*
*              Scale x if necessary to avoid overflow when adding a
*              multiple of column j of A.
*
               IF( xj.GT.one ) THEN
                  rec = one / xj
                  IF( cnorm( j ).GT.( bignum-xmax )*rec ) THEN
*
*                    Scale x by 1/(2*abs(x(j))).
*
                     rec = rec*half
                     CALL zdscal( n, rec, x, 1 )
                     scale = scale*rec
                  END IF
               ELSE IF( xj*cnorm( j ).GT.( bignum-xmax ) ) THEN
*
*                 Scale x by 1/2.
*
                  CALL zdscal( n, half, x, 1 )
                  scale = scale*half
               END IF
*
               IF( upper ) THEN
                  IF( j.GT.1 ) THEN
*
*                    Compute the update
*                       x(1:j-1) := x(1:j-1) - x(j) * A(1:j-1,j)
*
                     CALL zaxpy( j-1, -x( j )*tscal, a( 1, j ), 1, x,
     $                           1 )
                     i = izamax( j-1, x, 1 )
                     xmax = cabs1( x( i ) )
                  END IF
               ELSE
                  IF( j.LT.n ) THEN
*
*                    Compute the update
*                       x(j+1:n) := x(j+1:n) - x(j) * A(j+1:n,j)
*
                     CALL zaxpy( n-j, -x( j )*tscal, a( j+1, j ), 1,
     $                           x( j+1 ), 1 )
                     i = j + izamax( n-j, x( j+1 ), 1 )
                     xmax = cabs1( x( i ) )
                  END IF
               END IF
  120       CONTINUE
*
         ELSE IF( lsame( trans, 'T' ) ) THEN
*
*           Solve A**T * x = b
*
            DO 170 j = jfirst, jlast, jinc
*
*              Compute x(j) = b(j) - sum A(k,j)*x(k).
*                                    k<>j
*
               xj = cabs1( x( j ) )
               uscal = tscal
               rec = one / max( xmax, one )
               IF( cnorm( j ).GT.( bignum-xj )*rec ) THEN
*
*                 If x(j) could overflow, scale x by 1/(2*XMAX).
*
                  rec = rec*half
                  IF( nounit ) THEN
                     tjjs = a( j, j )*tscal
                  ELSE
                     tjjs = tscal
                  END IF
                  tjj = cabs1( tjjs )
                  IF( tjj.GT.one ) THEN
*
*                       Divide by A(j,j) when scaling x if A(j,j) > 1.
*
                     rec = min( one, rec*tjj )
                     uscal = zladiv( uscal, tjjs )
                  END IF
                  IF( rec.LT.one ) THEN
                     CALL zdscal( n, rec, x, 1 )
                     scale = scale*rec
                     xmax = xmax*rec
                  END IF
               END IF
*
               csumj = zero
               IF( uscal.EQ.dcmplx( one ) ) THEN
*
*                 If the scaling needed for A in the dot product is 1,
*                 call ZDOTU to perform the dot product.
*
                  IF( upper ) THEN
                     csumj = zdotu( j-1, a( 1, j ), 1, x, 1 )
                  ELSE IF( j.LT.n ) THEN
                     csumj = zdotu( n-j, a( j+1, j ), 1, x( j+1 ), 1 )
                  END IF
               ELSE
*
*                 Otherwise, use in-line code for the dot product.
*
                  IF( upper ) THEN
                     DO 130 i = 1, j - 1
                        csumj = csumj + ( a( i, j )*uscal )*x( i )
  130                CONTINUE
                  ELSE IF( j.LT.n ) THEN
                     DO 140 i = j + 1, n
                        csumj = csumj + ( a( i, j )*uscal )*x( i )
  140                CONTINUE
                  END IF
               END IF
*
               IF( uscal.EQ.dcmplx( tscal ) ) THEN
*
*                 Compute x(j) := ( x(j) - CSUMJ ) / A(j,j) if 1/A(j,j)
*                 was not used to scale the dotproduct.
*
                  x( j ) = x( j ) - csumj
                  xj = cabs1( x( j ) )
                  IF( nounit ) THEN
                     tjjs = a( j, j )*tscal
                  ELSE
                     tjjs = tscal
                     IF( tscal.EQ.one )
     $                  GO TO 160
                  END IF
*
*                    Compute x(j) = x(j) / A(j,j), scaling if necessary.
*
                  tjj = cabs1( tjjs )
                  IF( tjj.GT.smlnum ) THEN
*
*                       abs(A(j,j)) > SMLNUM:
*
                     IF( tjj.LT.one ) THEN
                        IF( xj.GT.tjj*bignum ) THEN
*
*                             Scale X by 1/abs(x(j)).
*
                           rec = one / xj
                           CALL zdscal( n, rec, x, 1 )
                           scale = scale*rec
                           xmax = xmax*rec
                        END IF
                     END IF
                     x( j ) = zladiv( x( j ), tjjs )
                  ELSE IF( tjj.GT.zero ) THEN
*
*                       0 < abs(A(j,j)) <= SMLNUM:
*
                     IF( xj.GT.tjj*bignum ) THEN
*
*                          Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.
*
                        rec = ( tjj*bignum ) / xj
                        CALL zdscal( n, rec, x, 1 )
                        scale = scale*rec
                        xmax = xmax*rec
                     END IF
                     x( j ) = zladiv( x( j ), tjjs )
                  ELSE
*
*                       A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
*                       scale = 0 and compute a solution to A**T *x = 0.
*
                     DO 150 i = 1, n
                        x( i ) = zero
  150                CONTINUE
                     x( j ) = one
                     scale = zero
                     xmax = zero
                  END IF
  160             CONTINUE
               ELSE
*
*                 Compute x(j) := x(j) / A(j,j) - CSUMJ if the dot
*                 product has already been divided by 1/A(j,j).
*
                  x( j ) = zladiv( x( j ), tjjs ) - csumj
               END IF
               xmax = max( xmax, cabs1( x( j ) ) )
  170       CONTINUE
*
         ELSE
*
*           Solve A**H * x = b
*
            DO 220 j = jfirst, jlast, jinc
*
*              Compute x(j) = b(j) - sum A(k,j)*x(k).
*                                    k<>j
*
               xj = cabs1( x( j ) )
               uscal = tscal
               rec = one / max( xmax, one )
               IF( cnorm( j ).GT.( bignum-xj )*rec ) THEN
*
*                 If x(j) could overflow, scale x by 1/(2*XMAX).
*
                  rec = rec*half
                  IF( nounit ) THEN
                     tjjs = dconjg( a( j, j ) )*tscal
                  ELSE
                     tjjs = tscal
                  END IF
                  tjj = cabs1( tjjs )
                  IF( tjj.GT.one ) THEN
*
*                       Divide by A(j,j) when scaling x if A(j,j) > 1.
*
                     rec = min( one, rec*tjj )
                     uscal = zladiv( uscal, tjjs )
                  END IF
                  IF( rec.LT.one ) THEN
                     CALL zdscal( n, rec, x, 1 )
                     scale = scale*rec
                     xmax = xmax*rec
                  END IF
               END IF
*
               csumj = zero
               IF( uscal.EQ.dcmplx( one ) ) THEN
*
*                 If the scaling needed for A in the dot product is 1,
*                 call ZDOTC to perform the dot product.
*
                  IF( upper ) THEN
                     csumj = zdotc( j-1, a( 1, j ), 1, x, 1 )
                  ELSE IF( j.LT.n ) THEN
                     csumj = zdotc( n-j, a( j+1, j ), 1, x( j+1 ), 1 )
                  END IF
               ELSE
*
*                 Otherwise, use in-line code for the dot product.
*
                  IF( upper ) THEN
                     DO 180 i = 1, j - 1
                        csumj = csumj + ( dconjg( a( i, j ) )*uscal )*
     $                          x( i )
  180                CONTINUE
                  ELSE IF( j.LT.n ) THEN
                     DO 190 i = j + 1, n
                        csumj = csumj + ( dconjg( a( i, j ) )*uscal )*
     $                          x( i )
  190                CONTINUE
                  END IF
               END IF
*
               IF( uscal.EQ.dcmplx( tscal ) ) THEN
*
*                 Compute x(j) := ( x(j) - CSUMJ ) / A(j,j) if 1/A(j,j)
*                 was not used to scale the dotproduct.
*
                  x( j ) = x( j ) - csumj
                  xj = cabs1( x( j ) )
                  IF( nounit ) THEN
                     tjjs = dconjg( a( j, j ) )*tscal
                  ELSE
                     tjjs = tscal
                     IF( tscal.EQ.one )
     $                  GO TO 210
                  END IF
*
*                    Compute x(j) = x(j) / A(j,j), scaling if necessary.
*
                  tjj = cabs1( tjjs )
                  IF( tjj.GT.smlnum ) THEN
*
*                       abs(A(j,j)) > SMLNUM:
*
                     IF( tjj.LT.one ) THEN
                        IF( xj.GT.tjj*bignum ) THEN
*
*                             Scale X by 1/abs(x(j)).
*
                           rec = one / xj
                           CALL zdscal( n, rec, x, 1 )
                           scale = scale*rec
                           xmax = xmax*rec
                        END IF
                     END IF
                     x( j ) = zladiv( x( j ), tjjs )
                  ELSE IF( tjj.GT.zero ) THEN
*
*                       0 < abs(A(j,j)) <= SMLNUM:
*
                     IF( xj.GT.tjj*bignum ) THEN
*
*                          Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.
*
                        rec = ( tjj*bignum ) / xj
                        CALL zdscal( n, rec, x, 1 )
                        scale = scale*rec
                        xmax = xmax*rec
                     END IF
                     x( j ) = zladiv( x( j ), tjjs )
                  ELSE
*
*                       A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
*                       scale = 0 and compute a solution to A**H *x = 0.
*
                     DO 200 i = 1, n
                        x( i ) = zero
  200                CONTINUE
                     x( j ) = one
                     scale = zero
                     xmax = zero
                  END IF
  210             CONTINUE
               ELSE
*
*                 Compute x(j) := x(j) / A(j,j) - CSUMJ if the dot
*                 product has already been divided by 1/A(j,j).
*
                  x( j ) = zladiv( x( j ), tjjs ) - csumj
               END IF
               xmax = max( xmax, cabs1( x( j ) ) )
  220       CONTINUE
         END IF
         scale = scale / tscal
      END IF
*
*     Scale the column norms by 1/TSCAL for return.
*
      IF( tscal.NE.one ) THEN
         CALL dscal( n, one / tscal, cnorm, 1 )
      END IF
*
      RETURN
*
*     End of ZLATRS
*
      END      
      
! ZCOPY
      SUBROUTINE zcopy(N,ZX,INCX,ZY,INCY)
*
*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      COMPLEX*16 ZX(*),ZY(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER I,IX,IY
*     ..
      IF (n.LE.0) RETURN
      IF (incx.EQ.1 .AND. incy.EQ.1) THEN
*
*        code for both increments equal to 1
*
         DO i = 1,n
          zy(i) = zx(i)
         END DO
      ELSE
*
*        code for unequal increments or equal increments
*          not equal to 1
*
         ix = 1
         iy = 1
         IF (incx.LT.0) ix = (-n+1)*incx + 1
         IF (incy.LT.0) iy = (-n+1)*incy + 1
         DO i = 1,n
            zy(iy) = zx(ix)
            ix = ix + incx
            iy = iy + incy
         END DO
      END IF
      RETURN
*
*     End of ZCOPY
*
      END      
      
! IZAMAX
      INTEGER FUNCTION izamax(N,ZX,INCX)
*
*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER incx,n
*     ..
*     .. Array Arguments ..
      COMPLEX*16 zx(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION dmax
      INTEGER i,ix
*     ..
*     .. External Functions ..
      DOUBLE PRECISION dcabs1
      EXTERNAL dcabs1
*     ..
      izamax = 0
      IF (n.LT.1 .OR. incx.LE.0) RETURN
      izamax = 1
      IF (n.EQ.1) RETURN
      IF (incx.EQ.1) THEN
*
*        code for increment equal to 1
*
         dmax = dcabs1(zx(1))
         DO i = 2,n
            IF (dcabs1(zx(i)).GT.dmax) THEN
               izamax = i
               dmax = dcabs1(zx(i))
            END IF
         END DO
      ELSE
*
*        code for increment not equal to 1
*
         ix = 1
         dmax = dcabs1(zx(1))
         ix = ix + incx
         DO i = 2,n
            IF (dcabs1(zx(ix)).GT.dmax) THEN
               izamax = i
               dmax = dcabs1(zx(ix))
            END IF
            ix = ix + incx
         END DO
      END IF
      RETURN
*
*     End of IZAMAX
*
      END
      
! ZLAQR0
      SUBROUTINE zlaqr0( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ,
     $                   IHIZ, Z, LDZ, WORK, LWORK, INFO )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, LWORK, N
      LOGICAL            WANTT, WANTZ
*     ..
*     .. Array Arguments ..
      COMPLEX*16         H( LDH, * ), W( * ), WORK( * ), Z( LDZ, * )
*     ..
*
*  ================================================================
*
*     .. Parameters ..
*
*     ==== Matrices of order NTINY or smaller must be processed by
*     .    ZLAHQR because of insufficient subdiagonal scratch space.
*     .    (This is a hard limit.) ====
      INTEGER            NTINY
      parameter( ntiny = 15 )
*
*     ==== Exceptional deflation windows:  try to cure rare
*     .    slow convergence by varying the size of the
*     .    deflation window after KEXNW iterations. ====
      INTEGER            KEXNW
      parameter( kexnw = 5 )
*
*     ==== Exceptional shifts: try to cure rare slow convergence
*     .    with ad-hoc exceptional shifts every KEXSH iterations.
*     .    ====
      INTEGER            KEXSH
      parameter( kexsh = 6 )
*
*     ==== The constant WILK1 is used to form the exceptional
*     .    shifts. ====
      DOUBLE PRECISION   WILK1
      parameter( wilk1 = 0.75d0 )
      COMPLEX*16         ZERO, ONE
      parameter( zero = ( 0.0d0, 0.0d0 ),
     $                   one = ( 1.0d0, 0.0d0 ) )
      DOUBLE PRECISION   TWO
      parameter( two = 2.0d0 )
*     ..
*     .. Local Scalars ..
      COMPLEX*16         AA, BB, CC, CDUM, DD, DET, RTDISC, SWAP, TR2
      DOUBLE PRECISION   S
      INTEGER            I, INF, IT, ITMAX, K, KACC22, KBOT, KDU, KS,
     $                   kt, ktop, ku, kv, kwh, kwtop, kwv, ld, ls,
     $                   lwkopt, ndec, ndfl, nh, nho, nibble, nmin, ns,
     $                   nsmax, nsr, nve, nw, nwmax, nwr, nwupbd
      LOGICAL            SORTED
      CHARACTER          JBCMPZ*2
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ilaenv
*     ..
*     .. Local Arrays ..
      COMPLEX*16         ZDUM( 1, 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           zlacpy, zlahqr, zlaqr3, zlaqr4, zlaqr5
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, dble, dcmplx, dimag, int, max, min, mod,
     $                   sqrt
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
*     ..
*     .. Statement Function definitions ..
      cabs1( cdum ) = abs( dble( cdum ) ) + abs( dimag( cdum ) )
*     ..
*     .. Executable Statements ..
      info = 0
*
*     ==== Quick return for N = 0: nothing to do. ====
*
      IF( n.EQ.0 ) THEN
         work( 1 ) = one
         RETURN
      END IF
*
      IF( n.LE.ntiny ) THEN
*
*        ==== Tiny matrices must use ZLAHQR. ====
*
         lwkopt = 1
         IF( lwork.NE.-1 )
     $      CALL zlahqr( wantt, wantz, n, ilo, ihi, h, ldh, w, iloz,
     $                   ihiz, z, ldz, info )
      ELSE
*
*        ==== Use small bulge multi-shift QR with aggressive early
*        .    deflation on larger-than-tiny matrices. ====
*
*        ==== Hope for the best. ====
*
         info = 0
*
*        ==== Set up job flags for ILAENV. ====
*
         IF( wantt ) THEN
            jbcmpz( 1: 1 ) = 'S'
         ELSE
            jbcmpz( 1: 1 ) = 'E'
         END IF
         IF( wantz ) THEN
            jbcmpz( 2: 2 ) = 'V'
         ELSE
            jbcmpz( 2: 2 ) = 'N'
         END IF
*
*        ==== NWR = recommended deflation window size.  At this
*        .    point,  N .GT. NTINY = 15, so there is enough
*        .    subdiagonal workspace for NWR.GE.2 as required.
*        .    (In fact, there is enough subdiagonal space for
*        .    NWR.GE.4.) ====
*
         nwr = ilaenv( 13, 'ZLAQR0', jbcmpz, n, ilo, ihi, lwork )
         nwr = max( 2, nwr )
         nwr = min( ihi-ilo+1, ( n-1 ) / 3, nwr )
*
*        ==== NSR = recommended number of simultaneous shifts.
*        .    At this point N .GT. NTINY = 15, so there is at
*        .    enough subdiagonal workspace for NSR to be even
*        .    and greater than or equal to two as required. ====
*
         nsr = ilaenv( 15, 'ZLAQR0', jbcmpz, n, ilo, ihi, lwork )
         nsr = min( nsr, ( n-3 ) / 6, ihi-ilo )
         nsr = max( 2, nsr-mod( nsr, 2 ) )
*
*        ==== Estimate optimal workspace ====
*
*        ==== Workspace query call to ZLAQR3 ====
*
         CALL zlaqr3( wantt, wantz, n, ilo, ihi, nwr+1, h, ldh, iloz,
     $                ihiz, z, ldz, ls, ld, w, h, ldh, n, h, ldh, n, h,
     $                ldh, work, -1 )
*
*        ==== Optimal workspace = MAX(ZLAQR5, ZLAQR3) ====
*
         lwkopt = max( 3*nsr / 2, int( work( 1 ) ) )
*
*        ==== Quick return in case of workspace query. ====
*
         IF( lwork.EQ.-1 ) THEN
            work( 1 ) = dcmplx( lwkopt, 0 )
            RETURN
         END IF
*
*        ==== ZLAHQR/ZLAQR0 crossover point ====
*
         nmin = ilaenv( 12, 'ZLAQR0', jbcmpz, n, ilo, ihi, lwork )
         nmin = max( ntiny, nmin )
*
*        ==== Nibble crossover point ====
*
         nibble = ilaenv( 14, 'ZLAQR0', jbcmpz, n, ilo, ihi, lwork )
         nibble = max( 0, nibble )
*
*        ==== Accumulate reflections during ttswp?  Use block
*        .    2-by-2 structure during matrix-matrix multiply? ====
*
         kacc22 = ilaenv( 16, 'ZLAQR0', jbcmpz, n, ilo, ihi, lwork )
         kacc22 = max( 0, kacc22 )
         kacc22 = min( 2, kacc22 )
*
*        ==== NWMAX = the largest possible deflation window for
*        .    which there is sufficient workspace. ====
*
         nwmax = min( ( n-1 ) / 3, lwork / 2 )
         nw = nwmax
*
*        ==== NSMAX = the Largest number of simultaneous shifts
*        .    for which there is sufficient workspace. ====
*
         nsmax = min( ( n-3 ) / 6, 2*lwork / 3 )
         nsmax = nsmax - mod( nsmax, 2 )
*
*        ==== NDFL: an iteration count restarted at deflation. ====
*
         ndfl = 1
*
*        ==== ITMAX = iteration limit ====
*
         itmax = max( 30, 2*kexsh )*max( 10, ( ihi-ilo+1 ) )
*
*        ==== Last row and column in the active block ====
*
         kbot = ihi
*
*        ==== Main Loop ====
*
         DO 70 it = 1, itmax
*
*           ==== Done when KBOT falls below ILO ====
*
            IF( kbot.LT.ilo )
     $         GO TO 80
*
*           ==== Locate active block ====
*
            DO 10 k = kbot, ilo + 1, -1
               IF( h( k, k-1 ).EQ.zero )
     $            GO TO 20
   10       CONTINUE
            k = ilo
   20       CONTINUE
            ktop = k
*
*           ==== Select deflation window size:
*           .    Typical Case:
*           .      If possible and advisable, nibble the entire
*           .      active block.  If not, use size MIN(NWR,NWMAX)
*           .      or MIN(NWR+1,NWMAX) depending upon which has
*           .      the smaller corresponding subdiagonal entry
*           .      (a heuristic).
*           .
*           .    Exceptional Case:
*           .      If there have been no deflations in KEXNW or
*           .      more iterations, then vary the deflation window
*           .      size.   At first, because, larger windows are,
*           .      in general, more powerful than smaller ones,
*           .      rapidly increase the window to the maximum possible.
*           .      Then, gradually reduce the window size. ====
*
            nh = kbot - ktop + 1
            nwupbd = min( nh, nwmax )
            IF( ndfl.LT.kexnw ) THEN
               nw = min( nwupbd, nwr )
            ELSE
               nw = min( nwupbd, 2*nw )
            END IF
            IF( nw.LT.nwmax ) THEN
               IF( nw.GE.nh-1 ) THEN
                  nw = nh
               ELSE
                  kwtop = kbot - nw + 1
                  IF( cabs1( h( kwtop, kwtop-1 ) ).GT.
     $                cabs1( h( kwtop-1, kwtop-2 ) ) )nw = nw + 1
               END IF
            END IF
            IF( ndfl.LT.kexnw ) THEN
               ndec = -1
            ELSE IF( ndec.GE.0 .OR. nw.GE.nwupbd ) THEN
               ndec = ndec + 1
               IF( nw-ndec.LT.2 )
     $            ndec = 0
               nw = nw - ndec
            END IF
*
*           ==== Aggressive early deflation:
*           .    split workspace under the subdiagonal into
*           .      - an nw-by-nw work array V in the lower
*           .        left-hand-corner,
*           .      - an NW-by-at-least-NW-but-more-is-better
*           .        (NW-by-NHO) horizontal work array along
*           .        the bottom edge,
*           .      - an at-least-NW-but-more-is-better (NHV-by-NW)
*           .        vertical work array along the left-hand-edge.
*           .        ====
*
            kv = n - nw + 1
            kt = nw + 1
            nho = ( n-nw-1 ) - kt + 1
            kwv = nw + 2
            nve = ( n-nw ) - kwv + 1
*
*           ==== Aggressive early deflation ====
*
            CALL zlaqr3( wantt, wantz, n, ktop, kbot, nw, h, ldh, iloz,
     $                   ihiz, z, ldz, ls, ld, w, h( kv, 1 ), ldh, nho,
     $                   h( kv, kt ), ldh, nve, h( kwv, 1 ), ldh, work,
     $                   lwork )
*
*           ==== Adjust KBOT accounting for new deflations. ====
*
            kbot = kbot - ld
*
*           ==== KS points to the shifts. ====
*
            ks = kbot - ls + 1
*
*           ==== Skip an expensive QR sweep if there is a (partly
*           .    heuristic) reason to expect that many eigenvalues
*           .    will deflate without it.  Here, the QR sweep is
*           .    skipped if many eigenvalues have just been deflated
*           .    or if the remaining active block is small.
*
            IF( ( ld.EQ.0 ) .OR. ( ( 100*ld.LE.nw*nibble ) .AND. ( kbot-
     $          ktop+1.GT.min( nmin, nwmax ) ) ) ) THEN
*
*              ==== NS = nominal number of simultaneous shifts.
*              .    This may be lowered (slightly) if ZLAQR3
*              .    did not provide that many shifts. ====
*
               ns = min( nsmax, nsr, max( 2, kbot-ktop ) )
               ns = ns - mod( ns, 2 )
*
*              ==== If there have been no deflations
*              .    in a multiple of KEXSH iterations,
*              .    then try exceptional shifts.
*              .    Otherwise use shifts provided by
*              .    ZLAQR3 above or from the eigenvalues
*              .    of a trailing principal submatrix. ====
*
               IF( mod( ndfl, kexsh ).EQ.0 ) THEN
                  ks = kbot - ns + 1
                  DO 30 i = kbot, ks + 1, -2
                     w( i ) = h( i, i ) + wilk1*cabs1( h( i, i-1 ) )
                     w( i-1 ) = w( i )
   30             CONTINUE
               ELSE
*
*                 ==== Got NS/2 or fewer shifts? Use ZLAQR4 or
*                 .    ZLAHQR on a trailing principal submatrix to
*                 .    get more. (Since NS.LE.NSMAX.LE.(N-3)/6,
*                 .    there is enough space below the subdiagonal
*                 .    to fit an NS-by-NS scratch array.) ====
*
                  IF( kbot-ks+1.LE.ns / 2 ) THEN
                     ks = kbot - ns + 1
                     kt = n - ns + 1
                     CALL zlacpy( 'A', ns, ns, h( ks, ks ), ldh,
     $                            h( kt, 1 ), ldh )
                     IF( ns.GT.nmin ) THEN
                        CALL zlaqr4( .false., .false., ns, 1, ns,
     $                               h( kt, 1 ), ldh, w( ks ), 1, 1,
     $                               zdum, 1, work, lwork, inf )
                     ELSE
                        CALL zlahqr( .false., .false., ns, 1, ns,
     $                               h( kt, 1 ), ldh, w( ks ), 1, 1,
     $                               zdum, 1, inf )
                     END IF
                     ks = ks + inf
*
*                    ==== In case of a rare QR failure use
*                    .    eigenvalues of the trailing 2-by-2
*                    .    principal submatrix.  Scale to avoid
*                    .    overflows, underflows and subnormals.
*                    .    (The scale factor S can not be zero,
*                    .    because H(KBOT,KBOT-1) is nonzero.) ====
*
                     IF( ks.GE.kbot ) THEN
                        s = cabs1( h( kbot-1, kbot-1 ) ) +
     $                      cabs1( h( kbot, kbot-1 ) ) +
     $                      cabs1( h( kbot-1, kbot ) ) +
     $                      cabs1( h( kbot, kbot ) )
                        aa = h( kbot-1, kbot-1 ) / s
                        cc = h( kbot, kbot-1 ) / s
                        bb = h( kbot-1, kbot ) / s
                        dd = h( kbot, kbot ) / s
                        tr2 = ( aa+dd ) / two
                        det = ( aa-tr2 )*( dd-tr2 ) - bb*cc
                        rtdisc = sqrt( -det )
                        w( kbot-1 ) = ( tr2+rtdisc )*s
                        w( kbot ) = ( tr2-rtdisc )*s
*
                        ks = kbot - 1
                     END IF
                  END IF
*
                  IF( kbot-ks+1.GT.ns ) THEN
*
*                    ==== Sort the shifts (Helps a little) ====
*
                     sorted = .false.
                     DO 50 k = kbot, ks + 1, -1
                        IF( sorted )
     $                     GO TO 60
                        sorted = .true.
                        DO 40 i = ks, k - 1
                           IF( cabs1( w( i ) ).LT.cabs1( w( i+1 ) ) )
     $                          THEN
                              sorted = .false.
                              swap = w( i )
                              w( i ) = w( i+1 )
                              w( i+1 ) = swap
                           END IF
   40                   CONTINUE
   50                CONTINUE
   60                CONTINUE
                  END IF
               END IF
*
*              ==== If there are only two shifts, then use
*              .    only one.  ====
*
               IF( kbot-ks+1.EQ.2 ) THEN
                  IF( cabs1( w( kbot )-h( kbot, kbot ) ).LT.
     $                cabs1( w( kbot-1 )-h( kbot, kbot ) ) ) THEN
                     w( kbot-1 ) = w( kbot )
                  ELSE
                     w( kbot ) = w( kbot-1 )
                  END IF
               END IF
*
*              ==== Use up to NS of the the smallest magnitude
*              .    shifts.  If there aren't NS shifts available,
*              .    then use them all, possibly dropping one to
*              .    make the number of shifts even. ====
*
               ns = min( ns, kbot-ks+1 )
               ns = ns - mod( ns, 2 )
               ks = kbot - ns + 1
*
*              ==== Small-bulge multi-shift QR sweep:
*              .    split workspace under the subdiagonal into
*              .    - a KDU-by-KDU work array U in the lower
*              .      left-hand-corner,
*              .    - a KDU-by-at-least-KDU-but-more-is-better
*              .      (KDU-by-NHo) horizontal work array WH along
*              .      the bottom edge,
*              .    - and an at-least-KDU-but-more-is-better-by-KDU
*              .      (NVE-by-KDU) vertical work WV arrow along
*              .      the left-hand-edge. ====
*
               kdu = 2*ns
               ku = n - kdu + 1
               kwh = kdu + 1
               nho = ( n-kdu+1-4 ) - ( kdu+1 ) + 1
               kwv = kdu + 4
               nve = n - kdu - kwv + 1
*
*              ==== Small-bulge multi-shift QR sweep ====
*
               CALL zlaqr5( wantt, wantz, kacc22, n, ktop, kbot, ns,
     $                      w( ks ), h, ldh, iloz, ihiz, z, ldz, work,
     $                      3, h( ku, 1 ), ldh, nve, h( kwv, 1 ), ldh,
     $                      nho, h( ku, kwh ), ldh )
            END IF
*
*           ==== Note progress (or the lack of it). ====
*
            IF( ld.GT.0 ) THEN
               ndfl = 1
            ELSE
               ndfl = ndfl + 1
            END IF
*
*           ==== End of main loop ====
   70    CONTINUE
*
*        ==== Iteration limit exceeded.  Set INFO to show where
*        .    the problem occurred and exit. ====
*
         info = kbot
   80    CONTINUE
      END IF
*
*     ==== Return the optimal value of LWORK. ====
*
      work( 1 ) = dcmplx( lwkopt, 0 )
*
*     ==== End of ZLAQR0 ====
*
      END
      
! ZLAHQR
      SUBROUTINE zlahqr( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ,
     $                   IHIZ, Z, LDZ, INFO )
      IMPLICIT NONE
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, N
      LOGICAL            WANTT, WANTZ
*     ..
*     .. Array Arguments ..
      COMPLEX*16         H( LDH, * ), W( * ), Z( LDZ, * )
*     ..
*
*  =========================================================
*
*     .. Parameters ..
      COMPLEX*16         ZERO, ONE
      parameter( zero = ( 0.0d0, 0.0d0 ),
     $                   one = ( 1.0d0, 0.0d0 ) )
      DOUBLE PRECISION   RZERO, RONE, HALF
      parameter( rzero = 0.0d0, rone = 1.0d0, half = 0.5d0 )
      DOUBLE PRECISION   DAT1
      parameter( dat1 = 3.0d0 / 4.0d0 )
      INTEGER            KEXSH
      parameter( kexsh = 10 )
*     ..
*     .. Local Scalars ..
      COMPLEX*16         CDUM, H11, H11S, H22, SC, SUM, T, T1, TEMP, U,
     $                   v2, x, y
      DOUBLE PRECISION   AA, AB, BA, BB, H10, H21, RTEMP, S, SAFMAX,
     $                   safmin, smlnum, sx, t2, tst, ulp
      INTEGER            I, I1, I2, ITS, ITMAX, J, JHI, JLO, K, L, M,
     $                   nh, nz, kdefl
*     ..
*     .. Local Arrays ..
      COMPLEX*16         V( 2 )
*     ..
*     .. External Functions ..
      COMPLEX*16         ZLADIV
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           zladiv, dlamch
*     ..
*     .. External Subroutines ..
      EXTERNAL           dlabad, zcopy, zlarfg, zscal
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, dble, dconjg, dimag, max, min, sqrt
*     ..
*     .. Statement Function definitions ..
      cabs1( cdum ) = abs( dble( cdum ) ) + abs( dimag( cdum ) )
*     ..
*     .. Executable Statements ..
*
      info = 0
*
*     Quick return if possible
*
      IF( n.EQ.0 )
     $   RETURN
      IF( ilo.EQ.ihi ) THEN
         w( ilo ) = h( ilo, ilo )
         RETURN
      END IF
*
*     ==== clear out the trash ====
      DO 10 j = ilo, ihi - 3
         h( j+2, j ) = zero
         h( j+3, j ) = zero
   10 CONTINUE
      IF( ilo.LE.ihi-2 )
     $   h( ihi, ihi-2 ) = zero
*     ==== ensure that subdiagonal entries are real ====
      IF( wantt ) THEN
         jlo = 1
         jhi = n
      ELSE
         jlo = ilo
         jhi = ihi
      END IF
      DO 20 i = ilo + 1, ihi
         IF( dimag( h( i, i-1 ) ).NE.rzero ) THEN
*           ==== The following redundant normalization
*           .    avoids problems with both gradual and
*           .    sudden underflow in ABS(H(I,I-1)) ====
            sc = h( i, i-1 ) / cabs1( h( i, i-1 ) )
            sc = dconjg( sc ) / abs( sc )
            h( i, i-1 ) = abs( h( i, i-1 ) )
            CALL zscal( jhi-i+1, sc, h( i, i ), ldh )
            CALL zscal( min( jhi, i+1 )-jlo+1, dconjg( sc ),
     $                  h( jlo, i ), 1 )
            IF( wantz )
     $         CALL zscal( ihiz-iloz+1, dconjg( sc ), z( iloz, i ), 1 )
         END IF
   20 CONTINUE
*
      nh = ihi - ilo + 1
      nz = ihiz - iloz + 1
*
*     Set machine-dependent constants for the stopping criterion.
*
      safmin = dlamch( 'SAFE MINIMUM' )
      safmax = rone / safmin
      CALL dlabad( safmin, safmax )
      ulp = dlamch( 'PRECISION' )
      smlnum = safmin*( dble( nh ) / ulp )
*
*     I1 and I2 are the indices of the first row and last column of H
*     to which transformations must be applied. If eigenvalues only are
*     being computed, I1 and I2 are set inside the main loop.
*
      IF( wantt ) THEN
         i1 = 1
         i2 = n
      END IF
*
*     ITMAX is the total number of QR iterations allowed.
*
      itmax = 30 * max( 10, nh )
*
*     KDEFL counts the number of iterations since a deflation
*
      kdefl = 0
*
*     The main loop begins here. I is the loop index and decreases from
*     IHI to ILO in steps of 1. Each iteration of the loop works
*     with the active submatrix in rows and columns L to I.
*     Eigenvalues I+1 to IHI have already converged. Either L = ILO, or
*     H(L,L-1) is negligible so that the matrix splits.
*
      i = ihi
   30 CONTINUE
      IF( i.LT.ilo )
     $   GO TO 150
*
*     Perform QR iterations on rows and columns ILO to I until a
*     submatrix of order 1 splits off at the bottom because a
*     subdiagonal element has become negligible.
*
      l = ilo
      DO 130 its = 0, itmax
*
*        Look for a single small subdiagonal element.
*
         DO 40 k = i, l + 1, -1
            IF( cabs1( h( k, k-1 ) ).LE.smlnum )
     $         GO TO 50
            tst = cabs1( h( k-1, k-1 ) ) + cabs1( h( k, k ) )
            IF( tst.EQ.zero ) THEN
               IF( k-2.GE.ilo )
     $            tst = tst + abs( dble( h( k-1, k-2 ) ) )
               IF( k+1.LE.ihi )
     $            tst = tst + abs( dble( h( k+1, k ) ) )
            END IF
*           ==== The following is a conservative small subdiagonal
*           .    deflation criterion due to Ahues & Tisseur (LAWN 122,
*           .    1997). It has better mathematical foundation and
*           .    improves accuracy in some examples.  ====
            IF( abs( dble( h( k, k-1 ) ) ).LE.ulp*tst ) THEN
               ab = max( cabs1( h( k, k-1 ) ), cabs1( h( k-1, k ) ) )
               ba = min( cabs1( h( k, k-1 ) ), cabs1( h( k-1, k ) ) )
               aa = max( cabs1( h( k, k ) ),
     $              cabs1( h( k-1, k-1 )-h( k, k ) ) )
               bb = min( cabs1( h( k, k ) ),
     $              cabs1( h( k-1, k-1 )-h( k, k ) ) )
               s = aa + ab
               IF( ba*( ab / s ).LE.max( smlnum,
     $             ulp*( bb*( aa / s ) ) ) )GO TO 50
            END IF
   40    CONTINUE
   50    CONTINUE
         l = k
         IF( l.GT.ilo ) THEN
*
*           H(L,L-1) is negligible
*
            h( l, l-1 ) = zero
         END IF
*
*        Exit from loop if a submatrix of order 1 has split off.
*
         IF( l.GE.i )
     $      GO TO 140
         kdefl = kdefl + 1
*
*        Now the active submatrix is in rows and columns L to I. If
*        eigenvalues only are being computed, only the active submatrix
*        need be transformed.
*
         IF( .NOT.wantt ) THEN
            i1 = l
            i2 = i
         END IF
*
         IF( mod(kdefl,2*kexsh).EQ.0 ) THEN
*
*           Exceptional shift.
*
            s = dat1*abs( dble( h( i, i-1 ) ) )
            t = s + h( i, i )
         ELSE IF( mod(kdefl,kexsh).EQ.0 ) THEN
*
*           Exceptional shift.
*
            s = dat1*abs( dble( h( l+1, l ) ) )
            t = s + h( l, l )
         ELSE
*
*           Wilkinson's shift.
*
            t = h( i, i )
            u = sqrt( h( i-1, i ) )*sqrt( h( i, i-1 ) )
            s = cabs1( u )
            IF( s.NE.rzero ) THEN
               x = half*( h( i-1, i-1 )-t )
               sx = cabs1( x )
               s = max( s, cabs1( x ) )
               y = s*sqrt( ( x / s )**2+( u / s )**2 )
               IF( sx.GT.rzero ) THEN
                  IF( dble( x / sx )*dble( y )+dimag( x / sx )*
     $                dimag( y ).LT.rzero )y = -y
               END IF
               t = t - u*zladiv( u, ( x+y ) )
            END IF
         END IF
*
*        Look for two consecutive small subdiagonal elements.
*
         DO 60 m = i - 1, l + 1, -1
*
*           Determine the effect of starting the single-shift QR
*           iteration at row M, and see if this would make H(M,M-1)
*           negligible.
*
            h11 = h( m, m )
            h22 = h( m+1, m+1 )
            h11s = h11 - t
            h21 = dble( h( m+1, m ) )
            s = cabs1( h11s ) + abs( h21 )
            h11s = h11s / s
            h21 = h21 / s
            v( 1 ) = h11s
            v( 2 ) = h21
            h10 = dble( h( m, m-1 ) )
            IF( abs( h10 )*abs( h21 ).LE.ulp*
     $          ( cabs1( h11s )*( cabs1( h11 )+cabs1( h22 ) ) ) )
     $          GO TO 70
   60    CONTINUE
         h11 = h( l, l )
         h22 = h( l+1, l+1 )
         h11s = h11 - t
         h21 = dble( h( l+1, l ) )
         s = cabs1( h11s ) + abs( h21 )
         h11s = h11s / s
         h21 = h21 / s
         v( 1 ) = h11s
         v( 2 ) = h21
   70    CONTINUE
*
*        Single-shift QR step
*
         DO 120 k = m, i - 1
*
*           The first iteration of this loop determines a reflection G
*           from the vector V and applies it from left and right to H,
*           thus creating a nonzero bulge below the subdiagonal.
*
*           Each subsequent iteration determines a reflection G to
*           restore the Hessenberg form in the (K-1)th column, and thus
*           chases the bulge one step toward the bottom of the active
*           submatrix.
*
*           V(2) is always real before the call to ZLARFG, and hence
*           after the call T2 ( = T1*V(2) ) is also real.
*
            IF( k.GT.m )
     $         CALL zcopy( 2, h( k, k-1 ), 1, v, 1 )
            CALL zlarfg( 2, v( 1 ), v( 2 ), 1, t1 )
            IF( k.GT.m ) THEN
               h( k, k-1 ) = v( 1 )
               h( k+1, k-1 ) = zero
            END IF
            v2 = v( 2 )
            t2 = dble( t1*v2 )
*
*           Apply G from the left to transform the rows of the matrix
*           in columns K to I2.
*
            DO 80 j = k, i2
               sum = dconjg( t1 )*h( k, j ) + t2*h( k+1, j )
               h( k, j ) = h( k, j ) - sum
               h( k+1, j ) = h( k+1, j ) - sum*v2
   80       CONTINUE
*
*           Apply G from the right to transform the columns of the
*           matrix in rows I1 to min(K+2,I).
*
            DO 90 j = i1, min( k+2, i )
               sum = t1*h( j, k ) + t2*h( j, k+1 )
               h( j, k ) = h( j, k ) - sum
               h( j, k+1 ) = h( j, k+1 ) - sum*dconjg( v2 )
   90       CONTINUE
*
            IF( wantz ) THEN
*
*              Accumulate transformations in the matrix Z
*
               DO 100 j = iloz, ihiz
                  sum = t1*z( j, k ) + t2*z( j, k+1 )
                  z( j, k ) = z( j, k ) - sum
                  z( j, k+1 ) = z( j, k+1 ) - sum*dconjg( v2 )
  100          CONTINUE
            END IF
*
            IF( k.EQ.m .AND. m.GT.l ) THEN
*
*              If the QR step was started at row M > L because two
*              consecutive small subdiagonals were found, then extra
*              scaling must be performed to ensure that H(M,M-1) remains
*              real.
*
               temp = one - t1
               temp = temp / abs( temp )
               h( m+1, m ) = h( m+1, m )*dconjg( temp )
               IF( m+2.LE.i )
     $            h( m+2, m+1 ) = h( m+2, m+1 )*temp
               DO 110 j = m, i
                  IF( j.NE.m+1 ) THEN
                     IF( i2.GT.j )
     $                  CALL zscal( i2-j, temp, h( j, j+1 ), ldh )
                     CALL zscal( j-i1, dconjg( temp ), h( i1, j ), 1 )
                     IF( wantz ) THEN
                        CALL zscal( nz, dconjg( temp ), z( iloz, j ),
     $                              1 )
                     END IF
                  END IF
  110          CONTINUE
            END IF
  120    CONTINUE
*
*        Ensure that H(I,I-1) is real.
*
         temp = h( i, i-1 )
         IF( dimag( temp ).NE.rzero ) THEN
            rtemp = abs( temp )
            h( i, i-1 ) = rtemp
            temp = temp / rtemp
            IF( i2.GT.i )
     $         CALL zscal( i2-i, dconjg( temp ), h( i, i+1 ), ldh )
            CALL zscal( i-i1, temp, h( i1, i ), 1 )
            IF( wantz ) THEN
               CALL zscal( nz, temp, z( iloz, i ), 1 )
            END IF
         END IF
*
  130 CONTINUE
*
*     Failure to converge in remaining number of iterations
*
      info = i
      RETURN
*
  140 CONTINUE
*
*     H(I,I-1) is negligible: one eigenvalue has converged.
*
      w( i ) = h( i, i )
*     reset deflation counter
      kdefl = 0
*
*     return to start of the main loop with new value of I.
*
      i = l - 1
      GO TO 30
*
  150 CONTINUE
      RETURN
*
*     End of ZLAHQR
*
      END

! ZLASSQ      
       SUBROUTINE ZLASSQ( N, X, INCX, SCALE, SUMSQ )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
      DOUBLE PRECISION   SCALE, SUMSQ
*     ..
*     .. Array Arguments ..
      COMPLEX*16         X( * )
*     ..
*
*  Purpose
*  =======
*
*  ZLASSQ returns the values scl and ssq such that
*
*     ( scl**2 )*ssq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
*
*  where x( i ) = abs( X( 1 + ( i - 1 )*INCX ) ). The value of sumsq is
*  assumed to be at least unity and the value of ssq will then satisfy
*
*     1.0 .le. ssq .le. ( sumsq + 2*n ).
*
*  scale is assumed to be non-negative and scl returns the value
*
*     scl = max( scale, abs( real( x( i ) ) ), abs( aimag( x( i ) ) ) ),
*            i
*
*  scale and sumsq must be supplied in SCALE and SUMSQ respectively.
*  SCALE and SUMSQ are overwritten by scl and ssq respectively.
*
*  The routine makes only one pass through the vector X.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of elements to be used from the vector X.
*
*  X       (input) COMPLEX*16 array, dimension (N)
*          The vector x as described above.
*             x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.
*
*  INCX    (input) INTEGER
*          The increment between successive values of the vector X.
*          INCX > 0.
*
*  SCALE   (input/output) DOUBLE PRECISION
*          On entry, the value  scale  in the equation above.
*          On exit, SCALE is overwritten with the value  scl .
*
*  SUMSQ   (input/output) DOUBLE PRECISION
*          On entry, the value  sumsq  in the equation above.
*          On exit, SUMSQ is overwritten with the value  ssq .
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            IX
      DOUBLE PRECISION   TEMP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DIMAG
*     ..
*     .. Executable Statements ..
*
      IF( N.GT.0 ) THEN
         DO 10 IX = 1, 1 + ( N-1 )*INCX, INCX
            IF( DBLE( X( IX ) ).NE.ZERO ) THEN
               TEMP1 = ABS( DBLE( X( IX ) ) )
               IF( SCALE.LT.TEMP1 ) THEN
                  SUMSQ = 1 + SUMSQ*( SCALE / TEMP1 )**2
                  SCALE = TEMP1
               ELSE
                  SUMSQ = SUMSQ + ( TEMP1 / SCALE )**2
               END IF
            END IF
            IF( DIMAG( X( IX ) ).NE.ZERO ) THEN
               TEMP1 = ABS( DIMAG( X( IX ) ) )
               IF( SCALE.LT.TEMP1 ) THEN
                  SUMSQ = 1 + SUMSQ*( SCALE / TEMP1 )**2
                  SCALE = TEMP1
               ELSE
                  SUMSQ = SUMSQ + ( TEMP1 / SCALE )**2
               END IF
            END IF
   10    CONTINUE
      END IF
*
      RETURN
*
*     End of ZLASSQ
*
      END
     
! ZLAHR2
      SUBROUTINE zlahr2( N, K, NB, A, LDA, TAU, T, LDT, Y, LDY )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            K, LDA, LDT, LDY, N, NB
*     ..
*     .. Array Arguments ..
      COMPLEX*16        A( LDA, * ), T( LDT, NB ), TAU( NB ),
     $                   Y( LDY, NB )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16        ZERO, ONE
      parameter( zero = ( 0.0d+0, 0.0d+0 ),
     $                     one = ( 1.0d+0, 0.0d+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I
      COMPLEX*16        EI
*     ..
*     .. External Subroutines ..
      EXTERNAL           zaxpy, zcopy, zgemm, zgemv, zlacpy,
     $                   zlarfg, zscal, ztrmm, ztrmv, zlacgv
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          min
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( n.LE.1 )
     $   RETURN
*
      DO 10 i = 1, nb
         IF( i.GT.1 ) THEN
*
*           Update A(K+1:N,I)
*
*           Update I-th column of A - Y * V**H
*
            CALL zlacgv( i-1, a( k+i-1, 1 ), lda )
            CALL zgemv( 'NO TRANSPOSE', n-k, i-1, -one, y(k+1,1), ldy,
     $                  a( k+i-1, 1 ), lda, one, a( k+1, i ), 1 )
            CALL zlacgv( i-1, a( k+i-1, 1 ), lda )
*
*           Apply I - V * T**H * V**H to this column (call it b) from the
*           left, using the last column of T as workspace
*
*           Let  V = ( V1 )   and   b = ( b1 )   (first I-1 rows)
*                    ( V2 )             ( b2 )
*
*           where V1 is unit lower triangular
*
*           w := V1**H * b1
*
            CALL zcopy( i-1, a( k+1, i ), 1, t( 1, nb ), 1 )
            CALL ztrmv( 'Lower', 'Conjugate transpose', 'UNIT',
     $                  i-1, a( k+1, 1 ),
     $                  lda, t( 1, nb ), 1 )
*
*           w := w + V2**H * b2
*
            CALL zgemv( 'Conjugate transpose', n-k-i+1, i-1,
     $                  one, a( k+i, 1 ),
     $                  lda, a( k+i, i ), 1, one, t( 1, nb ), 1 )
*
*           w := T**H * w
*
            CALL ztrmv( 'Upper', 'Conjugate transpose', 'NON-UNIT',
     $                  i-1, t, ldt,
     $                  t( 1, nb ), 1 )
*
*           b2 := b2 - V2*w
*
            CALL zgemv( 'NO TRANSPOSE', n-k-i+1, i-1, -one,
     $                  a( k+i, 1 ),
     $                  lda, t( 1, nb ), 1, one, a( k+i, i ), 1 )
*
*           b1 := b1 - V1*w
*
            CALL ztrmv( 'Lower', 'NO TRANSPOSE',
     $                  'UNIT', i-1,
     $                  a( k+1, 1 ), lda, t( 1, nb ), 1 )
            CALL zaxpy( i-1, -one, t( 1, nb ), 1, a( k+1, i ), 1 )
*
            a( k+i-1, i-1 ) = ei
         END IF
*
*        Generate the elementary reflector H(I) to annihilate
*        A(K+I+1:N,I)
*
         CALL zlarfg( n-k-i+1, a( k+i, i ), a( min( k+i+1, n ), i ), 1,
     $                tau( i ) )
         ei = a( k+i, i )
         a( k+i, i ) = one
*
*        Compute  Y(K+1:N,I)
*
         CALL zgemv( 'NO TRANSPOSE', n-k, n-k-i+1,
     $               one, a( k+1, i+1 ),
     $               lda, a( k+i, i ), 1, zero, y( k+1, i ), 1 )
         CALL zgemv( 'Conjugate transpose', n-k-i+1, i-1,
     $               one, a( k+i, 1 ), lda,
     $               a( k+i, i ), 1, zero, t( 1, i ), 1 )
         CALL zgemv( 'NO TRANSPOSE', n-k, i-1, -one,
     $               y( k+1, 1 ), ldy,
     $               t( 1, i ), 1, one, y( k+1, i ), 1 )
         CALL zscal( n-k, tau( i ), y( k+1, i ), 1 )
*
*        Compute T(1:I,I)
*
         CALL zscal( i-1, -tau( i ), t( 1, i ), 1 )
         CALL ztrmv( 'Upper', 'No Transpose', 'NON-UNIT',
     $               i-1, t, ldt,
     $               t( 1, i ), 1 )
         t( i, i ) = tau( i )
*
   10 CONTINUE
      a( k+nb, nb ) = ei
*
*     Compute Y(1:K,1:NB)
*
      CALL zlacpy( 'ALL', k, nb, a( 1, 2 ), lda, y, ldy )
      CALL ztrmm( 'RIGHT', 'Lower', 'NO TRANSPOSE',
     $            'UNIT', k, nb,
     $            one, a( k+1, 1 ), lda, y, ldy )
      IF( n.GT.k+nb )
     $   CALL zgemm( 'NO TRANSPOSE', 'NO TRANSPOSE', k,
     $               nb, n-k-nb, one,
     $               a( 1, 2+nb ), lda, a( k+1+nb, 1 ), lda, one, y,
     $               ldy )
      CALL ztrmm( 'RIGHT', 'Upper', 'NO TRANSPOSE',
     $            'NON-UNIT', k, nb,
     $            one, t, ldt, y, ldy )
*
      RETURN
*
*     End of ZLAHR2
*
      END
     
! ZTRMM
      SUBROUTINE ztrmm(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
*
*  -- Reference BLAS level3 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      COMPLEX*16 ALPHA
      INTEGER LDA,LDB,M,N
      CHARACTER DIAG,SIDE,TRANSA,UPLO
*     ..
*     .. Array Arguments ..
      COMPLEX*16 A(LDA,*),B(LDB,*)
*     ..
*
*  =====================================================================
*
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL lsame
*     ..
*     .. External Subroutines ..
      EXTERNAL xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC dconjg,max
*     ..
*     .. Local Scalars ..
      COMPLEX*16 TEMP
      INTEGER I,INFO,J,K,NROWA
      LOGICAL LSIDE,NOCONJ,NOUNIT,UPPER
*     ..
*     .. Parameters ..
      COMPLEX*16 ONE
      parameter(one= (1.0d+0,0.0d+0))
      COMPLEX*16 ZERO
      parameter(zero= (0.0d+0,0.0d+0))
*     ..
*
*     Test the input parameters.
*
      lside = lsame(side,'L')
      IF (lside) THEN
          nrowa = m
      ELSE
          nrowa = n
      END IF
      noconj = lsame(transa,'T')
      nounit = lsame(diag,'N')
      upper = lsame(uplo,'U')
*
      info = 0
      IF ((.NOT.lside) .AND. (.NOT.lsame(side,'R'))) THEN
          info = 1
      ELSE IF ((.NOT.upper) .AND. (.NOT.lsame(uplo,'L'))) THEN
          info = 2
      ELSE IF ((.NOT.lsame(transa,'N')) .AND.
     +         (.NOT.lsame(transa,'T')) .AND.
     +         (.NOT.lsame(transa,'C'))) THEN
          info = 3
      ELSE IF ((.NOT.lsame(diag,'U')) .AND. (.NOT.lsame(diag,'N'))) THEN
          info = 4
      ELSE IF (m.LT.0) THEN
          info = 5
      ELSE IF (n.LT.0) THEN
          info = 6
      ELSE IF (lda.LT.max(1,nrowa)) THEN
          info = 9
      ELSE IF (ldb.LT.max(1,m)) THEN
          info = 11
      END IF
      IF (info.NE.0) THEN
          CALL xerbla('ZTRMM ',info)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF (m.EQ.0 .OR. n.EQ.0) RETURN
*
*     And when  alpha.eq.zero.
*
      IF (alpha.EQ.zero) THEN
          DO 20 j = 1,n
              DO 10 i = 1,m
                  b(i,j) = zero
   10         CONTINUE
   20     CONTINUE
          RETURN
      END IF
*
*     Start the operations.
*
      IF (lside) THEN
          IF (lsame(transa,'N')) THEN
*
*           Form  B := alpha*A*B.
*
              IF (upper) THEN
                  DO 50 j = 1,n
                      DO 40 k = 1,m
                          IF (b(k,j).NE.zero) THEN
                              temp = alpha*b(k,j)
                              DO 30 i = 1,k - 1
                                  b(i,j) = b(i,j) + temp*a(i,k)
   30                         CONTINUE
                              IF (nounit) temp = temp*a(k,k)
                              b(k,j) = temp
                          END IF
   40                 CONTINUE
   50             CONTINUE
              ELSE
                  DO 80 j = 1,n
                      DO 70 k = m,1,-1
                          IF (b(k,j).NE.zero) THEN
                              temp = alpha*b(k,j)
                              b(k,j) = temp
                              IF (nounit) b(k,j) = b(k,j)*a(k,k)
                              DO 60 i = k + 1,m
                                  b(i,j) = b(i,j) + temp*a(i,k)
   60                         CONTINUE
                          END IF
   70                 CONTINUE
   80             CONTINUE
              END IF
          ELSE
*
*           Form  B := alpha*A**T*B   or   B := alpha*A**H*B.
*
              IF (upper) THEN
                  DO 120 j = 1,n
                      DO 110 i = m,1,-1
                          temp = b(i,j)
                          IF (noconj) THEN
                              IF (nounit) temp = temp*a(i,i)
                              DO 90 k = 1,i - 1
                                  temp = temp + a(k,i)*b(k,j)
   90                         CONTINUE
                          ELSE
                              IF (nounit) temp = temp*dconjg(a(i,i))
                              DO 100 k = 1,i - 1
                                  temp = temp + dconjg(a(k,i))*b(k,j)
  100                         CONTINUE
                          END IF
                          b(i,j) = alpha*temp
  110                 CONTINUE
  120             CONTINUE
              ELSE
                  DO 160 j = 1,n
                      DO 150 i = 1,m
                          temp = b(i,j)
                          IF (noconj) THEN
                              IF (nounit) temp = temp*a(i,i)
                              DO 130 k = i + 1,m
                                  temp = temp + a(k,i)*b(k,j)
  130                         CONTINUE
                          ELSE
                              IF (nounit) temp = temp*dconjg(a(i,i))
                              DO 140 k = i + 1,m
                                  temp = temp + dconjg(a(k,i))*b(k,j)
  140                         CONTINUE
                          END IF
                          b(i,j) = alpha*temp
  150                 CONTINUE
  160             CONTINUE
              END IF
          END IF
      ELSE
          IF (lsame(transa,'N')) THEN
*
*           Form  B := alpha*B*A.
*
              IF (upper) THEN
                  DO 200 j = n,1,-1
                      temp = alpha
                      IF (nounit) temp = temp*a(j,j)
                      DO 170 i = 1,m
                          b(i,j) = temp*b(i,j)
  170                 CONTINUE
                      DO 190 k = 1,j - 1
                          IF (a(k,j).NE.zero) THEN
                              temp = alpha*a(k,j)
                              DO 180 i = 1,m
                                  b(i,j) = b(i,j) + temp*b(i,k)
  180                         CONTINUE
                          END IF
  190                 CONTINUE
  200             CONTINUE
              ELSE
                  DO 240 j = 1,n
                      temp = alpha
                      IF (nounit) temp = temp*a(j,j)
                      DO 210 i = 1,m
                          b(i,j) = temp*b(i,j)
  210                 CONTINUE
                      DO 230 k = j + 1,n
                          IF (a(k,j).NE.zero) THEN
                              temp = alpha*a(k,j)
                              DO 220 i = 1,m
                                  b(i,j) = b(i,j) + temp*b(i,k)
  220                         CONTINUE
                          END IF
  230                 CONTINUE
  240             CONTINUE
              END IF
          ELSE
*
*           Form  B := alpha*B*A**T   or   B := alpha*B*A**H.
*
              IF (upper) THEN
                  DO 280 k = 1,n
                      DO 260 j = 1,k - 1
                          IF (a(j,k).NE.zero) THEN
                              IF (noconj) THEN
                                  temp = alpha*a(j,k)
                              ELSE
                                  temp = alpha*dconjg(a(j,k))
                              END IF
                              DO 250 i = 1,m
                                  b(i,j) = b(i,j) + temp*b(i,k)
  250                         CONTINUE
                          END IF
  260                 CONTINUE
                      temp = alpha
                      IF (nounit) THEN
                          IF (noconj) THEN
                              temp = temp*a(k,k)
                          ELSE
                              temp = temp*dconjg(a(k,k))
                          END IF
                      END IF
                      IF (temp.NE.one) THEN
                          DO 270 i = 1,m
                              b(i,k) = temp*b(i,k)
  270                     CONTINUE
                      END IF
  280             CONTINUE
              ELSE
                  DO 320 k = n,1,-1
                      DO 300 j = k + 1,n
                          IF (a(j,k).NE.zero) THEN
                              IF (noconj) THEN
                                  temp = alpha*a(j,k)
                              ELSE
                                  temp = alpha*dconjg(a(j,k))
                              END IF
                              DO 290 i = 1,m
                                  b(i,j) = b(i,j) + temp*b(i,k)
  290                         CONTINUE
                          END IF
  300                 CONTINUE
                      temp = alpha
                      IF (nounit) THEN
                          IF (noconj) THEN
                              temp = temp*a(k,k)
                          ELSE
                              temp = temp*dconjg(a(k,k))
                          END IF
                      END IF
                      IF (temp.NE.one) THEN
                          DO 310 i = 1,m
                              b(i,k) = temp*b(i,k)
  310                     CONTINUE
                      END IF
  320             CONTINUE
              END IF
          END IF
      END IF
*
      RETURN
*
*     End of ZTRMM
*
      END

! ZAXPY
      SUBROUTINE zaxpy(N,ZA,ZX,INCX,ZY,INCY)
*
*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      COMPLEX*16 ZA
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      COMPLEX*16 ZX(*),ZY(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER I,IX,IY
*     ..
*     .. External Functions ..
      DOUBLE PRECISION DCABS1
      EXTERNAL dcabs1
*     ..
      IF (n.LE.0) RETURN
      IF (dcabs1(za).EQ.0.0d0) RETURN
      IF (incx.EQ.1 .AND. incy.EQ.1) THEN
*
*        code for both increments equal to 1
*
         DO i = 1,n
            zy(i) = zy(i) + za*zx(i)
         END DO
      ELSE
*
*        code for unequal increments or equal increments
*          not equal to 1
*
         ix = 1
         iy = 1
         IF (incx.LT.0) ix = (-n+1)*incx + 1
         IF (incy.LT.0) iy = (-n+1)*incy + 1
         DO i = 1,n
            zy(iy) = zy(iy) + za*zx(ix)
            ix = ix + incx
            iy = iy + incy
         END DO
      END IF
*
      RETURN
*
*     End of ZAXPY
*
      END
      
! ZLARFB
      SUBROUTINE zlarfb( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV,
     $                   T, LDT, C, LDC, WORK, LDWORK )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          DIRECT, SIDE, STOREV, TRANS
      INTEGER            K, LDC, LDT, LDV, LDWORK, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         C( LDC, * ), T( LDT, * ), V( LDV, * ),
     $                   work( ldwork, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE
      parameter( one = ( 1.0d+0, 0.0d+0 ) )
*     ..
*     .. Local Scalars ..
      CHARACTER          TRANST
      INTEGER            I, J
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           lsame
*     ..
*     .. External Subroutines ..
      EXTERNAL           zcopy, zgemm, zlacgv, ztrmm
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          dconjg
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( m.LE.0 .OR. n.LE.0 )
     $   RETURN
*
      IF( lsame( trans, 'N' ) ) THEN
         transt = 'C'
      ELSE
         transt = 'N'
      END IF
*
      IF( lsame( storev, 'C' ) ) THEN
*
         IF( lsame( direct, 'F' ) ) THEN
*
*           Let  V =  ( V1 )    (first K rows)
*                     ( V2 )
*           where  V1  is unit lower triangular.
*
            IF( lsame( side, 'L' ) ) THEN
*
*              Form  H * C  or  H**H * C  where  C = ( C1 )
*                                                    ( C2 )
*
*              W := C**H * V  =  (C1**H * V1 + C2**H * V2)  (stored in WORK)
*
*              W := C1**H
*
               DO 10 j = 1, k
                  CALL zcopy( n, c( j, 1 ), ldc, work( 1, j ), 1 )
                  CALL zlacgv( n, work( 1, j ), 1 )
   10          CONTINUE
*
*              W := W * V1
*
               CALL ztrmm( 'Right', 'Lower', 'No transpose', 'Unit', n,
     $                     k, one, v, ldv, work, ldwork )
               IF( m.GT.k ) THEN
*
*                 W := W + C2**H * V2
*
                  CALL zgemm( 'Conjugate transpose', 'No transpose', n,
     $                        k, m-k, one, c( k+1, 1 ), ldc,
     $                        v( k+1, 1 ), ldv, one, work, ldwork )
               END IF
*
*              W := W * T**H  or  W * T
*
               CALL ztrmm( 'Right', 'Upper', transt, 'Non-unit', n, k,
     $                     one, t, ldt, work, ldwork )
*
*              C := C - V * W**H
*
               IF( m.GT.k ) THEN
*
*                 C2 := C2 - V2 * W**H
*
                  CALL zgemm( 'No transpose', 'Conjugate transpose',
     $                        m-k, n, k, -one, v( k+1, 1 ), ldv, work,
     $                        ldwork, one, c( k+1, 1 ), ldc )
               END IF
*
*              W := W * V1**H
*
               CALL ztrmm( 'Right', 'Lower', 'Conjugate transpose',
     $                     'Unit', n, k, one, v, ldv, work, ldwork )
*
*              C1 := C1 - W**H
*
               DO 30 j = 1, k
                  DO 20 i = 1, n
                     c( j, i ) = c( j, i ) - dconjg( work( i, j ) )
   20             CONTINUE
   30          CONTINUE
*
            ELSE IF( lsame( side, 'R' ) ) THEN
*
*              Form  C * H  or  C * H**H  where  C = ( C1  C2 )
*
*              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
*
*              W := C1
*
               DO 40 j = 1, k
                  CALL zcopy( m, c( 1, j ), 1, work( 1, j ), 1 )
   40          CONTINUE
*
*              W := W * V1
*
               CALL ztrmm( 'Right', 'Lower', 'No transpose', 'Unit', m,
     $                     k, one, v, ldv, work, ldwork )
               IF( n.GT.k ) THEN
*
*                 W := W + C2 * V2
*
                  CALL zgemm( 'No transpose', 'No transpose', m, k, n-k,
     $                        one, c( 1, k+1 ), ldc, v( k+1, 1 ), ldv,
     $                        one, work, ldwork )
               END IF
*
*              W := W * T  or  W * T**H
*
               CALL ztrmm( 'Right', 'Upper', trans, 'Non-unit', m, k,
     $                     one, t, ldt, work, ldwork )
*
*              C := C - W * V**H
*
               IF( n.GT.k ) THEN
*
*                 C2 := C2 - W * V2**H
*
                  CALL zgemm( 'No transpose', 'Conjugate transpose', m,
     $                        n-k, k, -one, work, ldwork, v( k+1, 1 ),
     $                        ldv, one, c( 1, k+1 ), ldc )
               END IF
*
*              W := W * V1**H
*
               CALL ztrmm( 'Right', 'Lower', 'Conjugate transpose',
     $                     'Unit', m, k, one, v, ldv, work, ldwork )
*
*              C1 := C1 - W
*
               DO 60 j = 1, k
                  DO 50 i = 1, m
                     c( i, j ) = c( i, j ) - work( i, j )
   50             CONTINUE
   60          CONTINUE
            END IF
*
         ELSE
*
*           Let  V =  ( V1 )
*                     ( V2 )    (last K rows)
*           where  V2  is unit upper triangular.
*
            IF( lsame( side, 'L' ) ) THEN
*
*              Form  H * C  or  H**H * C  where  C = ( C1 )
*                                                    ( C2 )
*
*              W := C**H * V  =  (C1**H * V1 + C2**H * V2)  (stored in WORK)
*
*              W := C2**H
*
               DO 70 j = 1, k
                  CALL zcopy( n, c( m-k+j, 1 ), ldc, work( 1, j ), 1 )
                  CALL zlacgv( n, work( 1, j ), 1 )
   70          CONTINUE
*
*              W := W * V2
*
               CALL ztrmm( 'Right', 'Upper', 'No transpose', 'Unit', n,
     $                     k, one, v( m-k+1, 1 ), ldv, work, ldwork )
               IF( m.GT.k ) THEN
*
*                 W := W + C1**H * V1
*
                  CALL zgemm( 'Conjugate transpose', 'No transpose', n,
     $                        k, m-k, one, c, ldc, v, ldv, one, work,
     $                        ldwork )
               END IF
*
*              W := W * T**H  or  W * T
*
               CALL ztrmm( 'Right', 'Lower', transt, 'Non-unit', n, k,
     $                     one, t, ldt, work, ldwork )
*
*              C := C - V * W**H
*
               IF( m.GT.k ) THEN
*
*                 C1 := C1 - V1 * W**H
*
                  CALL zgemm( 'No transpose', 'Conjugate transpose',
     $                        m-k, n, k, -one, v, ldv, work, ldwork,
     $                        one, c, ldc )
               END IF
*
*              W := W * V2**H
*
               CALL ztrmm( 'Right', 'Upper', 'Conjugate transpose',
     $                     'Unit', n, k, one, v( m-k+1, 1 ), ldv, work,
     $                     ldwork )
*
*              C2 := C2 - W**H
*
               DO 90 j = 1, k
                  DO 80 i = 1, n
                     c( m-k+j, i ) = c( m-k+j, i ) -
     $                               dconjg( work( i, j ) )
   80             CONTINUE
   90          CONTINUE
*
            ELSE IF( lsame( side, 'R' ) ) THEN
*
*              Form  C * H  or  C * H**H  where  C = ( C1  C2 )
*
*              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
*
*              W := C2
*
               DO 100 j = 1, k
                  CALL zcopy( m, c( 1, n-k+j ), 1, work( 1, j ), 1 )
  100          CONTINUE
*
*              W := W * V2
*
               CALL ztrmm( 'Right', 'Upper', 'No transpose', 'Unit', m,
     $                     k, one, v( n-k+1, 1 ), ldv, work, ldwork )
               IF( n.GT.k ) THEN
*
*                 W := W + C1 * V1
*
                  CALL zgemm( 'No transpose', 'No transpose', m, k, n-k,
     $                        one, c, ldc, v, ldv, one, work, ldwork )
               END IF
*
*              W := W * T  or  W * T**H
*
               CALL ztrmm( 'Right', 'Lower', trans, 'Non-unit', m, k,
     $                     one, t, ldt, work, ldwork )
*
*              C := C - W * V**H
*
               IF( n.GT.k ) THEN
*
*                 C1 := C1 - W * V1**H
*
                  CALL zgemm( 'No transpose', 'Conjugate transpose', m,
     $                        n-k, k, -one, work, ldwork, v, ldv, one,
     $                        c, ldc )
               END IF
*
*              W := W * V2**H
*
               CALL ztrmm( 'Right', 'Upper', 'Conjugate transpose',
     $                     'Unit', m, k, one, v( n-k+1, 1 ), ldv, work,
     $                     ldwork )
*
*              C2 := C2 - W
*
               DO 120 j = 1, k
                  DO 110 i = 1, m
                     c( i, n-k+j ) = c( i, n-k+j ) - work( i, j )
  110             CONTINUE
  120          CONTINUE
            END IF
         END IF
*
      ELSE IF( lsame( storev, 'R' ) ) THEN
*
         IF( lsame( direct, 'F' ) ) THEN
*
*           Let  V =  ( V1  V2 )    (V1: first K columns)
*           where  V1  is unit upper triangular.
*
            IF( lsame( side, 'L' ) ) THEN
*
*              Form  H * C  or  H**H * C  where  C = ( C1 )
*                                                    ( C2 )
*
*              W := C**H * V**H  =  (C1**H * V1**H + C2**H * V2**H) (stored in WORK)
*
*              W := C1**H
*
               DO 130 j = 1, k
                  CALL zcopy( n, c( j, 1 ), ldc, work( 1, j ), 1 )
                  CALL zlacgv( n, work( 1, j ), 1 )
  130          CONTINUE
*
*              W := W * V1**H
*
               CALL ztrmm( 'Right', 'Upper', 'Conjugate transpose',
     $                     'Unit', n, k, one, v, ldv, work, ldwork )
               IF( m.GT.k ) THEN
*
*                 W := W + C2**H * V2**H
*
                  CALL zgemm( 'Conjugate transpose',
     $                        'Conjugate transpose', n, k, m-k, one,
     $                        c( k+1, 1 ), ldc, v( 1, k+1 ), ldv, one,
     $                        work, ldwork )
               END IF
*
*              W := W * T**H  or  W * T
*
               CALL ztrmm( 'Right', 'Upper', transt, 'Non-unit', n, k,
     $                     one, t, ldt, work, ldwork )
*
*              C := C - V**H * W**H
*
               IF( m.GT.k ) THEN
*
*                 C2 := C2 - V2**H * W**H
*
                  CALL zgemm( 'Conjugate transpose',
     $                        'Conjugate transpose', m-k, n, k, -one,
     $                        v( 1, k+1 ), ldv, work, ldwork, one,
     $                        c( k+1, 1 ), ldc )
               END IF
*
*              W := W * V1
*
               CALL ztrmm( 'Right', 'Upper', 'No transpose', 'Unit', n,
     $                     k, one, v, ldv, work, ldwork )
*
*              C1 := C1 - W**H
*
               DO 150 j = 1, k
                  DO 140 i = 1, n
                     c( j, i ) = c( j, i ) - dconjg( work( i, j ) )
  140             CONTINUE
  150          CONTINUE
*
            ELSE IF( lsame( side, 'R' ) ) THEN
*
*              Form  C * H  or  C * H**H  where  C = ( C1  C2 )
*
*              W := C * V**H  =  (C1*V1**H + C2*V2**H)  (stored in WORK)
*
*              W := C1
*
               DO 160 j = 1, k
                  CALL zcopy( m, c( 1, j ), 1, work( 1, j ), 1 )
  160          CONTINUE
*
*              W := W * V1**H
*
               CALL ztrmm( 'Right', 'Upper', 'Conjugate transpose',
     $                     'Unit', m, k, one, v, ldv, work, ldwork )
               IF( n.GT.k ) THEN
*
*                 W := W + C2 * V2**H
*
                  CALL zgemm( 'No transpose', 'Conjugate transpose', m,
     $                        k, n-k, one, c( 1, k+1 ), ldc,
     $                        v( 1, k+1 ), ldv, one, work, ldwork )
               END IF
*
*              W := W * T  or  W * T**H
*
               CALL ztrmm( 'Right', 'Upper', trans, 'Non-unit', m, k,
     $                     one, t, ldt, work, ldwork )
*
*              C := C - W * V
*
               IF( n.GT.k ) THEN
*
*                 C2 := C2 - W * V2
*
                  CALL zgemm( 'No transpose', 'No transpose', m, n-k, k,
     $                        -one, work, ldwork, v( 1, k+1 ), ldv, one,
     $                        c( 1, k+1 ), ldc )
               END IF
*
*              W := W * V1
*
               CALL ztrmm( 'Right', 'Upper', 'No transpose', 'Unit', m,
     $                     k, one, v, ldv, work, ldwork )
*
*              C1 := C1 - W
*
               DO 180 j = 1, k
                  DO 170 i = 1, m
                     c( i, j ) = c( i, j ) - work( i, j )
  170             CONTINUE
  180          CONTINUE
*
            END IF
*
         ELSE
*
*           Let  V =  ( V1  V2 )    (V2: last K columns)
*           where  V2  is unit lower triangular.
*
            IF( lsame( side, 'L' ) ) THEN
*
*              Form  H * C  or  H**H * C  where  C = ( C1 )
*                                                    ( C2 )
*
*              W := C**H * V**H  =  (C1**H * V1**H + C2**H * V2**H) (stored in WORK)
*
*              W := C2**H
*
               DO 190 j = 1, k
                  CALL zcopy( n, c( m-k+j, 1 ), ldc, work( 1, j ), 1 )
                  CALL zlacgv( n, work( 1, j ), 1 )
  190          CONTINUE
*
*              W := W * V2**H
*
               CALL ztrmm( 'Right', 'Lower', 'Conjugate transpose',
     $                     'Unit', n, k, one, v( 1, m-k+1 ), ldv, work,
     $                     ldwork )
               IF( m.GT.k ) THEN
*
*                 W := W + C1**H * V1**H
*
                  CALL zgemm( 'Conjugate transpose',
     $                        'Conjugate transpose', n, k, m-k, one, c,
     $                        ldc, v, ldv, one, work, ldwork )
               END IF
*
*              W := W * T**H  or  W * T
*
               CALL ztrmm( 'Right', 'Lower', transt, 'Non-unit', n, k,
     $                     one, t, ldt, work, ldwork )
*
*              C := C - V**H * W**H
*
               IF( m.GT.k ) THEN
*
*                 C1 := C1 - V1**H * W**H
*
                  CALL zgemm( 'Conjugate transpose',
     $                        'Conjugate transpose', m-k, n, k, -one, v,
     $                        ldv, work, ldwork, one, c, ldc )
               END IF
*
*              W := W * V2
*
               CALL ztrmm( 'Right', 'Lower', 'No transpose', 'Unit', n,
     $                     k, one, v( 1, m-k+1 ), ldv, work, ldwork )
*
*              C2 := C2 - W**H
*
               DO 210 j = 1, k
                  DO 200 i = 1, n
                     c( m-k+j, i ) = c( m-k+j, i ) -
     $                               dconjg( work( i, j ) )
  200             CONTINUE
  210          CONTINUE
*
            ELSE IF( lsame( side, 'R' ) ) THEN
*
*              Form  C * H  or  C * H**H  where  C = ( C1  C2 )
*
*              W := C * V**H  =  (C1*V1**H + C2*V2**H)  (stored in WORK)
*
*              W := C2
*
               DO 220 j = 1, k
                  CALL zcopy( m, c( 1, n-k+j ), 1, work( 1, j ), 1 )
  220          CONTINUE
*
*              W := W * V2**H
*
               CALL ztrmm( 'Right', 'Lower', 'Conjugate transpose',
     $                     'Unit', m, k, one, v( 1, n-k+1 ), ldv, work,
     $                     ldwork )
               IF( n.GT.k ) THEN
*
*                 W := W + C1 * V1**H
*
                  CALL zgemm( 'No transpose', 'Conjugate transpose', m,
     $                        k, n-k, one, c, ldc, v, ldv, one, work,
     $                        ldwork )
               END IF
*
*              W := W * T  or  W * T**H
*
               CALL ztrmm( 'Right', 'Lower', trans, 'Non-unit', m, k,
     $                     one, t, ldt, work, ldwork )
*
*              C := C - W * V
*
               IF( n.GT.k ) THEN
*
*                 C1 := C1 - W * V1
*
                  CALL zgemm( 'No transpose', 'No transpose', m, n-k, k,
     $                        -one, work, ldwork, v, ldv, one, c, ldc )
               END IF
*
*              W := W * V2
*
               CALL ztrmm( 'Right', 'Lower', 'No transpose', 'Unit', m,
     $                     k, one, v( 1, n-k+1 ), ldv, work, ldwork )
*
*              C1 := C1 - W
*
               DO 240 j = 1, k
                  DO 230 i = 1, m
                     c( i, n-k+j ) = c( i, n-k+j ) - work( i, j )
  230             CONTINUE
  240          CONTINUE
*
            END IF
*
         END IF
      END IF
*
      RETURN
*
*     End of ZLARFB
*
      END
      
! ZGEHD2
      SUBROUTINE zgehd2( N, ILO, IHI, A, LDA, TAU, WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            IHI, ILO, INFO, LDA, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE
      parameter( one = ( 1.0d+0, 0.0d+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I
      COMPLEX*16         ALPHA
*     ..
*     .. External Subroutines ..
      EXTERNAL           xerbla, zlarf, zlarfg
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          dconjg, max, min
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      info = 0
      IF( n.LT.0 ) THEN
         info = -1
      ELSE IF( ilo.LT.1 .OR. ilo.GT.max( 1, n ) ) THEN
         info = -2
      ELSE IF( ihi.LT.min( ilo, n ) .OR. ihi.GT.n ) THEN
         info = -3
      ELSE IF( lda.LT.max( 1, n ) ) THEN
         info = -5
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'ZGEHD2', -info )
         RETURN
      END IF
*
      DO 10 i = ilo, ihi - 1
*
*        Compute elementary reflector H(i) to annihilate A(i+2:ihi,i)
*
         alpha = a( i+1, i )
         CALL zlarfg( ihi-i, alpha, a( min( i+2, n ), i ), 1, tau( i ) )
         a( i+1, i ) = one
*
*        Apply H(i) to A(1:ihi,i+1:ihi) from the right
*
         CALL zlarf( 'Right', ihi, ihi-i, a( i+1, i ), 1, tau( i ),
     $               a( 1, i+1 ), lda, work )
*
*        Apply H(i)**H to A(i+1:ihi,i+1:n) from the left
*
         CALL zlarf( 'Left', ihi-i, n-i, a( i+1, i ), 1,
     $               dconjg( tau( i ) ), a( i+1, i+1 ), lda, work )
*
         a( i+1, i ) = alpha
   10 CONTINUE
*
      RETURN
*
*     End of ZGEHD2
*
      END
      
! ZUNGQR
      SUBROUTINE zungqr( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ZERO
      parameter( zero = ( 0.0d+0, 0.0d+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, J, KI, KK, L, LDWORK,
     $                   LWKOPT, NB, NBMIN, NX
*     ..
*     .. External Subroutines ..
      EXTERNAL           xerbla, zlarfb, zlarft, zung2r
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          max, min
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ilaenv
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      info = 0
      nb = ilaenv( 1, 'ZUNGQR', ' ', m, n, k, -1 )
      lwkopt = max( 1, n )*nb
      work( 1 ) = lwkopt
      lquery = ( lwork.EQ.-1 )
      IF( m.LT.0 ) THEN
         info = -1
      ELSE IF( n.LT.0 .OR. n.GT.m ) THEN
         info = -2
      ELSE IF( k.LT.0 .OR. k.GT.n ) THEN
         info = -3
      ELSE IF( lda.LT.max( 1, m ) ) THEN
         info = -5
      ELSE IF( lwork.LT.max( 1, n ) .AND. .NOT.lquery ) THEN
         info = -8
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'ZUNGQR', -info )
         RETURN
      ELSE IF( lquery ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( n.LE.0 ) THEN
         work( 1 ) = 1
         RETURN
      END IF
*
      nbmin = 2
      nx = 0
      iws = n
      IF( nb.GT.1 .AND. nb.LT.k ) THEN
*
*        Determine when to cross over from blocked to unblocked code.
*
         nx = max( 0, ilaenv( 3, 'ZUNGQR', ' ', m, n, k, -1 ) )
         IF( nx.LT.k ) THEN
*
*           Determine if workspace is large enough for blocked code.
*
            ldwork = n
            iws = ldwork*nb
            IF( lwork.LT.iws ) THEN
*
*              Not enough workspace to use optimal NB:  reduce NB and
*              determine the minimum value of NB.
*
               nb = lwork / ldwork
               nbmin = max( 2, ilaenv( 2, 'ZUNGQR', ' ', m, n, k, -1 ) )
            END IF
         END IF
      END IF
*
      IF( nb.GE.nbmin .AND. nb.LT.k .AND. nx.LT.k ) THEN
*
*        Use blocked code after the last block.
*        The first kk columns are handled by the block method.
*
         ki = ( ( k-nx-1 ) / nb )*nb
         kk = min( k, ki+nb )
*
*        Set A(1:kk,kk+1:n) to zero.
*
         DO 20 j = kk + 1, n
            DO 10 i = 1, kk
               a( i, j ) = zero
   10       CONTINUE
   20    CONTINUE
      ELSE
         kk = 0
      END IF
*
*     Use unblocked code for the last or only block.
*
      IF( kk.LT.n )
     $   CALL zung2r( m-kk, n-kk, k-kk, a( kk+1, kk+1 ), lda,
     $                tau( kk+1 ), work, iinfo )
*
      IF( kk.GT.0 ) THEN
*
*        Use blocked code
*
         DO 50 i = ki + 1, 1, -nb
            ib = min( nb, k-i+1 )
            IF( i+ib.LE.n ) THEN
*
*              Form the triangular factor of the block reflector
*              H = H(i) H(i+1) . . . H(i+ib-1)
*
               CALL zlarft( 'Forward', 'Columnwise', m-i+1, ib,
     $                      a( i, i ), lda, tau( i ), work, ldwork )
*
*              Apply H to A(i:m,i+ib:n) from the left
*
               CALL zlarfb( 'Left', 'No transpose', 'Forward',
     $                      'Columnwise', m-i+1, n-i-ib+1, ib,
     $                      a( i, i ), lda, work, ldwork, a( i, i+ib ),
     $                      lda, work( ib+1 ), ldwork )
            END IF
*
*           Apply H to rows i:m of current block
*
            CALL zung2r( m-i+1, ib, ib, a( i, i ), lda, tau( i ), work,
     $                   iinfo )
*
*           Set rows 1:i-1 of current block to zero
*
            DO 40 j = i, i + ib - 1
               DO 30 l = 1, i - 1
                  a( l, j ) = zero
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
      END IF
*
      work( 1 ) = iws
      RETURN
*
*     End of ZUNGQR
*
      END

! ZTRTI2
      SUBROUTINE ztrti2( UPLO, DIAG, N, A, LDA, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE
      parameter( one = ( 1.0d+0, 0.0d+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOUNIT, UPPER
      INTEGER            J
      COMPLEX*16         AJJ
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           lsame
*     ..
*     .. External Subroutines ..
      EXTERNAL           xerbla, zscal, ztrmv
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          max
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      info = 0
      upper = lsame( uplo, 'U' )
      nounit = lsame( diag, 'N' )
      IF( .NOT.upper .AND. .NOT.lsame( uplo, 'L' ) ) THEN
         info = -1
      ELSE IF( .NOT.nounit .AND. .NOT.lsame( diag, 'U' ) ) THEN
         info = -2
      ELSE IF( n.LT.0 ) THEN
         info = -3
      ELSE IF( lda.LT.max( 1, n ) ) THEN
         info = -5
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'ZTRTI2', -info )
         RETURN
      END IF
*
      IF( upper ) THEN
*
*        Compute inverse of upper triangular matrix.
*
         DO 10 j = 1, n
            IF( nounit ) THEN
               a( j, j ) = one / a( j, j )
               ajj = -a( j, j )
            ELSE
               ajj = -one
            END IF
*
*           Compute elements 1:j-1 of j-th column.
*
            CALL ztrmv( 'Upper', 'No transpose', diag, j-1, a, lda,
     $                  a( 1, j ), 1 )
            CALL zscal( j-1, ajj, a( 1, j ), 1 )
   10    CONTINUE
      ELSE
*
*        Compute inverse of lower triangular matrix.
*
         DO 20 j = n, 1, -1
            IF( nounit ) THEN
               a( j, j ) = one / a( j, j )
               ajj = -a( j, j )
            ELSE
               ajj = -one
            END IF
            IF( j.LT.n ) THEN
*
*              Compute elements j+1:n of j-th column.
*
               CALL ztrmv( 'Lower', 'No transpose', diag, n-j,
     $                     a( j+1, j+1 ), lda, a( j+1, j ), 1 )
               CALL zscal( n-j, ajj, a( j+1, j ), 1 )
            END IF
   20    CONTINUE
      END IF
*
      RETURN
*
*     End of ZTRTI2
*
      END
      
! DCABS1
      DOUBLE PRECISION FUNCTION dcabs1(Z)
*
*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      COMPLEX*16 z
*     ..
*     ..
*  =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC abs,dble,dimag
*
      dcabs1 = abs(dble(z)) + abs(dimag(z))
      RETURN
*
*     End of DCABS1
*
      END
      
! ZTRSV
      SUBROUTINE ztrsv(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
*
*  -- Reference BLAS level2 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER INCX,LDA,N
      CHARACTER DIAG,TRANS,UPLO
*     ..
*     .. Array Arguments ..
      COMPLEX*16 A(LDA,*),X(*)
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16 ZERO
      parameter(zero= (0.0d+0,0.0d+0))
*     ..
*     .. Local Scalars ..
      COMPLEX*16 TEMP
      INTEGER I,INFO,IX,J,JX,KX
      LOGICAL NOCONJ,NOUNIT
*     ..
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL lsame
*     ..
*     .. External Subroutines ..
      EXTERNAL xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC dconjg,max
*     ..
*
*     Test the input parameters.
*
      info = 0
      IF (.NOT.lsame(uplo,'U') .AND. .NOT.lsame(uplo,'L')) THEN
          info = 1
      ELSE IF (.NOT.lsame(trans,'N') .AND. .NOT.lsame(trans,'T') .AND.
     +         .NOT.lsame(trans,'C')) THEN
          info = 2
      ELSE IF (.NOT.lsame(diag,'U') .AND. .NOT.lsame(diag,'N')) THEN
          info = 3
      ELSE IF (n.LT.0) THEN
          info = 4
      ELSE IF (lda.LT.max(1,n)) THEN
          info = 6
      ELSE IF (incx.EQ.0) THEN
          info = 8
      END IF
      IF (info.NE.0) THEN
          CALL xerbla('ZTRSV ',info)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF (n.EQ.0) RETURN
*
      noconj = lsame(trans,'T')
      nounit = lsame(diag,'N')
*
*     Set up the start point in X if the increment is not unity. This
*     will be  ( N - 1 )*INCX  too small for descending loops.
*
      IF (incx.LE.0) THEN
          kx = 1 - (n-1)*incx
      ELSE IF (incx.NE.1) THEN
          kx = 1
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF (lsame(trans,'N')) THEN
*
*        Form  x := inv( A )*x.
*
          IF (lsame(uplo,'U')) THEN
              IF (incx.EQ.1) THEN
                  DO 20 j = n,1,-1
                      IF (x(j).NE.zero) THEN
                          IF (nounit) x(j) = x(j)/a(j,j)
                          temp = x(j)
                          DO 10 i = j - 1,1,-1
                              x(i) = x(i) - temp*a(i,j)
   10                     CONTINUE
                      END IF
   20             CONTINUE
              ELSE
                  jx = kx + (n-1)*incx
                  DO 40 j = n,1,-1
                      IF (x(jx).NE.zero) THEN
                          IF (nounit) x(jx) = x(jx)/a(j,j)
                          temp = x(jx)
                          ix = jx
                          DO 30 i = j - 1,1,-1
                              ix = ix - incx
                              x(ix) = x(ix) - temp*a(i,j)
   30                     CONTINUE
                      END IF
                      jx = jx - incx
   40             CONTINUE
              END IF
          ELSE
              IF (incx.EQ.1) THEN
                  DO 60 j = 1,n
                      IF (x(j).NE.zero) THEN
                          IF (nounit) x(j) = x(j)/a(j,j)
                          temp = x(j)
                          DO 50 i = j + 1,n
                              x(i) = x(i) - temp*a(i,j)
   50                     CONTINUE
                      END IF
   60             CONTINUE
              ELSE
                  jx = kx
                  DO 80 j = 1,n
                      IF (x(jx).NE.zero) THEN
                          IF (nounit) x(jx) = x(jx)/a(j,j)
                          temp = x(jx)
                          ix = jx
                          DO 70 i = j + 1,n
                              ix = ix + incx
                              x(ix) = x(ix) - temp*a(i,j)
   70                     CONTINUE
                      END IF
                      jx = jx + incx
   80             CONTINUE
              END IF
          END IF
      ELSE
*
*        Form  x := inv( A**T )*x  or  x := inv( A**H )*x.
*
          IF (lsame(uplo,'U')) THEN
              IF (incx.EQ.1) THEN
                  DO 110 j = 1,n
                      temp = x(j)
                      IF (noconj) THEN
                          DO 90 i = 1,j - 1
                              temp = temp - a(i,j)*x(i)
   90                     CONTINUE
                          IF (nounit) temp = temp/a(j,j)
                      ELSE
                          DO 100 i = 1,j - 1
                              temp = temp - dconjg(a(i,j))*x(i)
  100                     CONTINUE
                          IF (nounit) temp = temp/dconjg(a(j,j))
                      END IF
                      x(j) = temp
  110             CONTINUE
              ELSE
                  jx = kx
                  DO 140 j = 1,n
                      ix = kx
                      temp = x(jx)
                      IF (noconj) THEN
                          DO 120 i = 1,j - 1
                              temp = temp - a(i,j)*x(ix)
                              ix = ix + incx
  120                     CONTINUE
                          IF (nounit) temp = temp/a(j,j)
                      ELSE
                          DO 130 i = 1,j - 1
                              temp = temp - dconjg(a(i,j))*x(ix)
                              ix = ix + incx
  130                     CONTINUE
                          IF (nounit) temp = temp/dconjg(a(j,j))
                      END IF
                      x(jx) = temp
                      jx = jx + incx
  140             CONTINUE
              END IF
          ELSE
              IF (incx.EQ.1) THEN
                  DO 170 j = n,1,-1
                      temp = x(j)
                      IF (noconj) THEN
                          DO 150 i = n,j + 1,-1
                              temp = temp - a(i,j)*x(i)
  150                     CONTINUE
                          IF (nounit) temp = temp/a(j,j)
                      ELSE
                          DO 160 i = n,j + 1,-1
                              temp = temp - dconjg(a(i,j))*x(i)
  160                     CONTINUE
                          IF (nounit) temp = temp/dconjg(a(j,j))
                      END IF
                      x(j) = temp
  170             CONTINUE
              ELSE
                  kx = kx + (n-1)*incx
                  jx = kx
                  DO 200 j = n,1,-1
                      ix = kx
                      temp = x(jx)
                      IF (noconj) THEN
                          DO 180 i = n,j + 1,-1
                              temp = temp - a(i,j)*x(ix)
                              ix = ix - incx
  180                     CONTINUE
                          IF (nounit) temp = temp/a(j,j)
                      ELSE
                          DO 190 i = n,j + 1,-1
                              temp = temp - dconjg(a(i,j))*x(ix)
                              ix = ix - incx
  190                     CONTINUE
                          IF (nounit) temp = temp/dconjg(a(j,j))
                      END IF
                      x(jx) = temp
                      jx = jx - incx
  200             CONTINUE
              END IF
          END IF
      END IF
*
      RETURN
*
*     End of ZTRSV
*
      END
      
! ZLADIV
      COMPLEX*16 FUNCTION zladiv( X, Y )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      COMPLEX*16         x, y
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION   zi, zr
*     ..
*     .. External Subroutines ..
      EXTERNAL           dladiv
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          dble, dcmplx, dimag
*     ..
*     .. Executable Statements ..
*
      CALL dladiv( dble( x ), dimag( x ), dble( y ), dimag( y ), zr,
     $             zi )
      zladiv = dcmplx( zr, zi )
*
      RETURN
*
*     End of ZLADIV
*
      END
      
! ZDOTU
      COMPLEX*16 FUNCTION zdotu(N,ZX,INCX,ZY,INCY)
*
*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER incx,incy,n
*     ..
*     .. Array Arguments ..
      COMPLEX*16 zx(*),zy(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      COMPLEX*16 ztemp
      INTEGER i,ix,iy
*     ..
      ztemp = (0.0d0,0.0d0)
      zdotu = (0.0d0,0.0d0)
      IF (n.LE.0) RETURN
      IF (incx.EQ.1 .AND. incy.EQ.1) THEN
*
*        code for both increments equal to 1
*
         DO i = 1,n
            ztemp = ztemp + zx(i)*zy(i)
         END DO
      ELSE
*
*        code for unequal increments or equal increments
*          not equal to 1
*
         ix = 1
         iy = 1
         IF (incx.LT.0) ix = (-n+1)*incx + 1
         IF (incy.LT.0) iy = (-n+1)*incy + 1
         DO i = 1,n
            ztemp = ztemp + zx(ix)*zy(iy)
            ix = ix + incx
            iy = iy + incy
         END DO
      END IF
      zdotu = ztemp
      RETURN
*
*     End of ZDOTU
*
      END

! ZDOTC
      COMPLEX*16 FUNCTION zdotc(N,ZX,INCX,ZY,INCY)
*
*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER incx,incy,n
*     ..
*     .. Array Arguments ..
      COMPLEX*16 zx(*),zy(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      COMPLEX*16 ztemp
      INTEGER i,ix,iy
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC dconjg
*     ..
      ztemp = (0.0d0,0.0d0)
      zdotc = (0.0d0,0.0d0)
      IF (n.LE.0) RETURN
      IF (incx.EQ.1 .AND. incy.EQ.1) THEN
*
*        code for both increments equal to 1
*
         DO i = 1,n
            ztemp = ztemp + dconjg(zx(i))*zy(i)
         END DO
      ELSE
*
*        code for unequal increments or equal increments
*          not equal to 1
*
         ix = 1
         iy = 1
         IF (incx.LT.0) ix = (-n+1)*incx + 1
         IF (incy.LT.0) iy = (-n+1)*incy + 1
         DO i = 1,n
            ztemp = ztemp + dconjg(zx(ix))*zy(iy)
            ix = ix + incx
            iy = iy + incy
         END DO
      END IF
      zdotc = ztemp
      RETURN
*
*     End of ZDOTC
*
      END
      
! ZLAQR3
      SUBROUTINE zlaqr3( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ,
     $                   IHIZ, Z, LDZ, NS, ND, SH, V, LDV, NH, T, LDT,
     $                   NV, WV, LDWV, WORK, LWORK )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            IHIZ, ILOZ, KBOT, KTOP, LDH, LDT, LDV, LDWV,
     $                   LDZ, LWORK, N, ND, NH, NS, NV, NW
      LOGICAL            WANTT, WANTZ
*     ..
*     .. Array Arguments ..
      COMPLEX*16         H( LDH, * ), SH( * ), T( LDT, * ), V( LDV, * ),
     $                   WORK( * ), WV( LDWV, * ), Z( LDZ, * )
*     ..
*
*  ================================================================
*
*     .. Parameters ..
      COMPLEX*16         ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0d0, 0.0d0 ),
     $                   one = ( 1.0d0, 0.0d0 ) )
      DOUBLE PRECISION   RZERO, RONE
      PARAMETER          ( RZERO = 0.0d0, rone = 1.0d0 )
*     ..
*     .. Local Scalars ..
      COMPLEX*16         BETA, CDUM, S, TAU
      DOUBLE PRECISION   FOO, SAFMAX, SAFMIN, SMLNUM, ULP
      INTEGER            I, IFST, ILST, INFO, INFQR, J, JW, KCOL, KLN,
     $                   knt, krow, kwtop, ltop, lwk1, lwk2, lwk3,
     $                   lwkopt, nmin
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      INTEGER            ILAENV
      EXTERNAL           dlamch, ilaenv
*     ..
*     .. External Subroutines ..
      EXTERNAL           dlabad, zcopy, zgehrd, zgemm, zlacpy, zlahqr,
     $                   zlaqr4, zlarf, zlarfg, zlaset, ztrexc, zunmhr
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, dble, dcmplx, dconjg, dimag, int, max, min
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
*     ..
*     .. Statement Function definitions ..
      cabs1( cdum ) = abs( dble( cdum ) ) + abs( dimag( cdum ) )
*     ..
*     .. Executable Statements ..
*
*     ==== Estimate optimal workspace. ====
*
      jw = min( nw, kbot-ktop+1 )
      IF( jw.LE.2 ) THEN
         lwkopt = 1
      ELSE
*
*        ==== Workspace query call to ZGEHRD ====
*
         CALL zgehrd( jw, 1, jw-1, t, ldt, work, work, -1, info )
         lwk1 = int( work( 1 ) )
*
*        ==== Workspace query call to ZUNMHR ====
*
         CALL zunmhr( 'R', 'N', jw, jw, 1, jw-1, t, ldt, work, v, ldv,
     $                work, -1, info )
         lwk2 = int( work( 1 ) )
*
*        ==== Workspace query call to ZLAQR4 ====
*
         CALL zlaqr4( .true., .true., jw, 1, jw, t, ldt, sh, 1, jw, v,
     $                ldv, work, -1, infqr )
         lwk3 = int( work( 1 ) )
*
*        ==== Optimal workspace ====
*
         lwkopt = max( jw+max( lwk1, lwk2 ), lwk3 )
      END IF
*
*     ==== Quick return in case of workspace query. ====
*
      IF( lwork.EQ.-1 ) THEN
         work( 1 ) = dcmplx( lwkopt, 0 )
         RETURN
      END IF
*
*     ==== Nothing to do ...
*     ... for an empty active block ... ====
      ns = 0
      nd = 0
      work( 1 ) = one
      IF( ktop.GT.kbot )
     $   RETURN
*     ... nor for an empty deflation window. ====
      IF( nw.LT.1 )
     $   RETURN
*
*     ==== Machine constants ====
*
      safmin = dlamch( 'SAFE MINIMUM' )
      safmax = rone / safmin
      CALL dlabad( safmin, safmax )
      ulp = dlamch( 'PRECISION' )
      smlnum = safmin*( dble( n ) / ulp )
*
*     ==== Setup deflation window ====
*
      jw = min( nw, kbot-ktop+1 )
      kwtop = kbot - jw + 1
      IF( kwtop.EQ.ktop ) THEN
         s = zero
      ELSE
         s = h( kwtop, kwtop-1 )
      END IF
*
      IF( kbot.EQ.kwtop ) THEN
*
*        ==== 1-by-1 deflation window: not much to do ====
*
         sh( kwtop ) = h( kwtop, kwtop )
         ns = 1
         nd = 0
         IF( cabs1( s ).LE.max( smlnum, ulp*cabs1( h( kwtop,
     $       kwtop ) ) ) ) THEN
            ns = 0
            nd = 1
            IF( kwtop.GT.ktop )
     $         h( kwtop, kwtop-1 ) = zero
         END IF
         work( 1 ) = one
         RETURN
      END IF
*
*     ==== Convert to spike-triangular form.  (In case of a
*     .    rare QR failure, this routine continues to do
*     .    aggressive early deflation using that part of
*     .    the deflation window that converged using INFQR
*     .    here and there to keep track.) ====
*
      CALL zlacpy( 'U', jw, jw, h( kwtop, kwtop ), ldh, t, ldt )
      CALL zcopy( jw-1, h( kwtop+1, kwtop ), ldh+1, t( 2, 1 ), ldt+1 )
*
      CALL zlaset( 'A', jw, jw, zero, one, v, ldv )
      nmin = ilaenv( 12, 'ZLAQR3', 'SV', jw, 1, jw, lwork )
      IF( jw.GT.nmin ) THEN
         CALL zlaqr4( .true., .true., jw, 1, jw, t, ldt, sh( kwtop ), 1,
     $                jw, v, ldv, work, lwork, infqr )
      ELSE
         CALL zlahqr( .true., .true., jw, 1, jw, t, ldt, sh( kwtop ), 1,
     $                jw, v, ldv, infqr )
      END IF
*
*     ==== Deflation detection loop ====
*
      ns = jw
      ilst = infqr + 1
      DO 10 knt = infqr + 1, jw
*
*        ==== Small spike tip deflation test ====
*
         foo = cabs1( t( ns, ns ) )
         IF( foo.EQ.rzero )
     $      foo = cabs1( s )
         IF( cabs1( s )*cabs1( v( 1, ns ) ).LE.max( smlnum, ulp*foo ) )
     $        THEN
*
*           ==== One more converged eigenvalue ====
*
            ns = ns - 1
         ELSE
*
*           ==== One undeflatable eigenvalue.  Move it up out of the
*           .    way.   (ZTREXC can not fail in this case.) ====
*
            ifst = ns
            CALL ztrexc( 'V', jw, t, ldt, v, ldv, ifst, ilst, info )
            ilst = ilst + 1
         END IF
   10 CONTINUE
*
*        ==== Return to Hessenberg form ====
*
      IF( ns.EQ.0 )
     $   s = zero
*
      IF( ns.LT.jw ) THEN
*
*        ==== sorting the diagonal of T improves accuracy for
*        .    graded matrices.  ====
*
         DO 30 i = infqr + 1, ns
            ifst = i
            DO 20 j = i + 1, ns
               IF( cabs1( t( j, j ) ).GT.cabs1( t( ifst, ifst ) ) )
     $            ifst = j
   20       CONTINUE
            ilst = i
            IF( ifst.NE.ilst )
     $         CALL ztrexc( 'V', jw, t, ldt, v, ldv, ifst, ilst, info )
   30    CONTINUE
      END IF
*
*     ==== Restore shift/eigenvalue array from T ====
*
      DO 40 i = infqr + 1, jw
         sh( kwtop+i-1 ) = t( i, i )
   40 CONTINUE
*
*
      IF( ns.LT.jw .OR. s.EQ.zero ) THEN
         IF( ns.GT.1 .AND. s.NE.zero ) THEN
*
*           ==== Reflect spike back into lower triangle ====
*
            CALL zcopy( ns, v, ldv, work, 1 )
            DO 50 i = 1, ns
               work( i ) = dconjg( work( i ) )
   50       CONTINUE
            beta = work( 1 )
            CALL zlarfg( ns, beta, work( 2 ), 1, tau )
            work( 1 ) = one
*
            CALL zlaset( 'L', jw-2, jw-2, zero, zero, t( 3, 1 ), ldt )
*
            CALL zlarf( 'L', ns, jw, work, 1, dconjg( tau ), t, ldt,
     $                  work( jw+1 ) )
            CALL zlarf( 'R', ns, ns, work, 1, tau, t, ldt,
     $                  work( jw+1 ) )
            CALL zlarf( 'R', jw, ns, work, 1, tau, v, ldv,
     $                  work( jw+1 ) )
*
            CALL zgehrd( jw, 1, ns, t, ldt, work, work( jw+1 ),
     $                   lwork-jw, info )
         END IF
*
*        ==== Copy updated reduced window into place ====
*
         IF( kwtop.GT.1 )
     $      h( kwtop, kwtop-1 ) = s*dconjg( v( 1, 1 ) )
         CALL zlacpy( 'U', jw, jw, t, ldt, h( kwtop, kwtop ), ldh )
         CALL zcopy( jw-1, t( 2, 1 ), ldt+1, h( kwtop+1, kwtop ),
     $               ldh+1 )
*
*        ==== Accumulate orthogonal matrix in order update
*        .    H and Z, if requested.  ====
*
         IF( ns.GT.1 .AND. s.NE.zero )
     $      CALL zunmhr( 'R', 'N', jw, ns, 1, ns, t, ldt, work, v, ldv,
     $                   work( jw+1 ), lwork-jw, info )
*
*        ==== Update vertical slab in H ====
*
         IF( wantt ) THEN
            ltop = 1
         ELSE
            ltop = ktop
         END IF
         DO 60 krow = ltop, kwtop - 1, nv
            kln = min( nv, kwtop-krow )
            CALL zgemm( 'N', 'N', kln, jw, jw, one, h( krow, kwtop ),
     $                  ldh, v, ldv, zero, wv, ldwv )
            CALL zlacpy( 'A', kln, jw, wv, ldwv, h( krow, kwtop ), ldh )
   60    CONTINUE
*
*        ==== Update horizontal slab in H ====
*
         IF( wantt ) THEN
            DO 70 kcol = kbot + 1, n, nh
               kln = min( nh, n-kcol+1 )
               CALL zgemm( 'C', 'N', jw, kln, jw, one, v, ldv,
     $                     h( kwtop, kcol ), ldh, zero, t, ldt )
               CALL zlacpy( 'A', jw, kln, t, ldt, h( kwtop, kcol ),
     $                      ldh )
   70       CONTINUE
         END IF
*
*        ==== Update vertical slab in Z ====
*
         IF( wantz ) THEN
            DO 80 krow = iloz, ihiz, nv
               kln = min( nv, ihiz-krow+1 )
               CALL zgemm( 'N', 'N', kln, jw, jw, one, z( krow, kwtop ),
     $                     ldz, v, ldv, zero, wv, ldwv )
               CALL zlacpy( 'A', kln, jw, wv, ldwv, z( krow, kwtop ),
     $                      ldz )
   80       CONTINUE
         END IF
      END IF
*
*     ==== Return the number of deflations ... ====
*
      nd = jw - ns
*
*     ==== ... and the number of shifts. (Subtracting
*     .    INFQR from the spike length takes care
*     .    of the case of a rare QR failure while
*     .    calculating eigenvalues of the deflation
*     .    window.)  ====
*
      ns = ns - infqr
*
*      ==== Return optimal workspace. ====
*
      work( 1 ) = dcmplx( lwkopt, 0 )
*
*     ==== End of ZLAQR3 ====
*
      END
      
! ZLAQR4
      SUBROUTINE zlaqr4( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ,
     $                   IHIZ, Z, LDZ, WORK, LWORK, INFO )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, LWORK, N
      LOGICAL            WANTT, WANTZ
*     ..
*     .. Array Arguments ..
      COMPLEX*16         H( LDH, * ), W( * ), WORK( * ), Z( LDZ, * )
*     ..
*
*  ================================================================
*
*     .. Parameters ..
*
*     ==== Matrices of order NTINY or smaller must be processed by
*     .    ZLAHQR because of insufficient subdiagonal scratch space.
*     .    (This is a hard limit.) ====
      INTEGER            NTINY
      parameter( ntiny = 15 )
*
*     ==== Exceptional deflation windows:  try to cure rare
*     .    slow convergence by varying the size of the
*     .    deflation window after KEXNW iterations. ====
      INTEGER            KEXNW
      parameter( kexnw = 5 )
*
*     ==== Exceptional shifts: try to cure rare slow convergence
*     .    with ad-hoc exceptional shifts every KEXSH iterations.
*     .    ====
      INTEGER            KEXSH
      parameter( kexsh = 6 )
*
*     ==== The constant WILK1 is used to form the exceptional
*     .    shifts. ====
      DOUBLE PRECISION   WILK1
      parameter( wilk1 = 0.75d0 )
      COMPLEX*16         ZERO, ONE
      parameter( zero = ( 0.0d0, 0.0d0 ),
     $                   one = ( 1.0d0, 0.0d0 ) )
      DOUBLE PRECISION   TWO
      parameter( two = 2.0d0 )
*     ..
*     .. Local Scalars ..
      COMPLEX*16         AA, BB, CC, CDUM, DD, DET, RTDISC, SWAP, TR2
      DOUBLE PRECISION   S
      INTEGER            I, INF, IT, ITMAX, K, KACC22, KBOT, KDU, KS,
     $                   kt, ktop, ku, kv, kwh, kwtop, kwv, ld, ls,
     $                   lwkopt, ndec, ndfl, nh, nho, nibble, nmin, ns,
     $                   nsmax, nsr, nve, nw, nwmax, nwr, nwupbd
      LOGICAL            SORTED
      CHARACTER          JBCMPZ*2
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ilaenv
*     ..
*     .. Local Arrays ..
      COMPLEX*16         ZDUM( 1, 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           zlacpy, zlahqr, zlaqr2, zlaqr5
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, dble, dcmplx, dimag, int, max, min, mod,
     $                   sqrt
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
*     ..
*     .. Statement Function definitions ..
      cabs1( cdum ) = abs( dble( cdum ) ) + abs( dimag( cdum ) )
*     ..
*     .. Executable Statements ..
      info = 0
*
*     ==== Quick return for N = 0: nothing to do. ====
*
      IF( n.EQ.0 ) THEN
         work( 1 ) = one
         RETURN
      END IF
*
      IF( n.LE.ntiny ) THEN
*
*        ==== Tiny matrices must use ZLAHQR. ====
*
         lwkopt = 1
         IF( lwork.NE.-1 )
     $      CALL zlahqr( wantt, wantz, n, ilo, ihi, h, ldh, w, iloz,
     $                   ihiz, z, ldz, info )
      ELSE
*
*        ==== Use small bulge multi-shift QR with aggressive early
*        .    deflation on larger-than-tiny matrices. ====
*
*        ==== Hope for the best. ====
*
         info = 0
*
*        ==== Set up job flags for ILAENV. ====
*
         IF( wantt ) THEN
            jbcmpz( 1: 1 ) = 'S'
         ELSE
            jbcmpz( 1: 1 ) = 'E'
         END IF
         IF( wantz ) THEN
            jbcmpz( 2: 2 ) = 'V'
         ELSE
            jbcmpz( 2: 2 ) = 'N'
         END IF
*
*        ==== NWR = recommended deflation window size.  At this
*        .    point,  N .GT. NTINY = 15, so there is enough
*        .    subdiagonal workspace for NWR.GE.2 as required.
*        .    (In fact, there is enough subdiagonal space for
*        .    NWR.GE.4.) ====
*
         nwr = ilaenv( 13, 'ZLAQR4', jbcmpz, n, ilo, ihi, lwork )
         nwr = max( 2, nwr )
         nwr = min( ihi-ilo+1, ( n-1 ) / 3, nwr )
*
*        ==== NSR = recommended number of simultaneous shifts.
*        .    At this point N .GT. NTINY = 15, so there is at
*        .    enough subdiagonal workspace for NSR to be even
*        .    and greater than or equal to two as required. ====
*
         nsr = ilaenv( 15, 'ZLAQR4', jbcmpz, n, ilo, ihi, lwork )
         nsr = min( nsr, ( n-3 ) / 6, ihi-ilo )
         nsr = max( 2, nsr-mod( nsr, 2 ) )
*
*        ==== Estimate optimal workspace ====
*
*        ==== Workspace query call to ZLAQR2 ====
*
         CALL zlaqr2( wantt, wantz, n, ilo, ihi, nwr+1, h, ldh, iloz,
     $                ihiz, z, ldz, ls, ld, w, h, ldh, n, h, ldh, n, h,
     $                ldh, work, -1 )
*
*        ==== Optimal workspace = MAX(ZLAQR5, ZLAQR2) ====
*
         lwkopt = max( 3*nsr / 2, int( work( 1 ) ) )
*
*        ==== Quick return in case of workspace query. ====
*
         IF( lwork.EQ.-1 ) THEN
            work( 1 ) = dcmplx( lwkopt, 0 )
            RETURN
         END IF
*
*        ==== ZLAHQR/ZLAQR0 crossover point ====
*
         nmin = ilaenv( 12, 'ZLAQR4', jbcmpz, n, ilo, ihi, lwork )
         nmin = max( ntiny, nmin )
*
*        ==== Nibble crossover point ====
*
         nibble = ilaenv( 14, 'ZLAQR4', jbcmpz, n, ilo, ihi, lwork )
         nibble = max( 0, nibble )
*
*        ==== Accumulate reflections during ttswp?  Use block
*        .    2-by-2 structure during matrix-matrix multiply? ====
*
         kacc22 = ilaenv( 16, 'ZLAQR4', jbcmpz, n, ilo, ihi, lwork )
         kacc22 = max( 0, kacc22 )
         kacc22 = min( 2, kacc22 )
*
*        ==== NWMAX = the largest possible deflation window for
*        .    which there is sufficient workspace. ====
*
         nwmax = min( ( n-1 ) / 3, lwork / 2 )
         nw = nwmax
*
*        ==== NSMAX = the Largest number of simultaneous shifts
*        .    for which there is sufficient workspace. ====
*
         nsmax = min( ( n-3 ) / 6, 2*lwork / 3 )
         nsmax = nsmax - mod( nsmax, 2 )
*
*        ==== NDFL: an iteration count restarted at deflation. ====
*
         ndfl = 1
*
*        ==== ITMAX = iteration limit ====
*
         itmax = max( 30, 2*kexsh )*max( 10, ( ihi-ilo+1 ) )
*
*        ==== Last row and column in the active block ====
*
         kbot = ihi
*
*        ==== Main Loop ====
*
         DO 70 it = 1, itmax
*
*           ==== Done when KBOT falls below ILO ====
*
            IF( kbot.LT.ilo )
     $         GO TO 80
*
*           ==== Locate active block ====
*
            DO 10 k = kbot, ilo + 1, -1
               IF( h( k, k-1 ).EQ.zero )
     $            GO TO 20
   10       CONTINUE
            k = ilo
   20       CONTINUE
            ktop = k
*
*           ==== Select deflation window size:
*           .    Typical Case:
*           .      If possible and advisable, nibble the entire
*           .      active block.  If not, use size MIN(NWR,NWMAX)
*           .      or MIN(NWR+1,NWMAX) depending upon which has
*           .      the smaller corresponding subdiagonal entry
*           .      (a heuristic).
*           .
*           .    Exceptional Case:
*           .      If there have been no deflations in KEXNW or
*           .      more iterations, then vary the deflation window
*           .      size.   At first, because, larger windows are,
*           .      in general, more powerful than smaller ones,
*           .      rapidly increase the window to the maximum possible.
*           .      Then, gradually reduce the window size. ====
*
            nh = kbot - ktop + 1
            nwupbd = min( nh, nwmax )
            IF( ndfl.LT.kexnw ) THEN
               nw = min( nwupbd, nwr )
            ELSE
               nw = min( nwupbd, 2*nw )
            END IF
            IF( nw.LT.nwmax ) THEN
               IF( nw.GE.nh-1 ) THEN
                  nw = nh
               ELSE
                  kwtop = kbot - nw + 1
                  IF( cabs1( h( kwtop, kwtop-1 ) ).GT.
     $                cabs1( h( kwtop-1, kwtop-2 ) ) )nw = nw + 1
               END IF
            END IF
            IF( ndfl.LT.kexnw ) THEN
               ndec = -1
            ELSE IF( ndec.GE.0 .OR. nw.GE.nwupbd ) THEN
               ndec = ndec + 1
               IF( nw-ndec.LT.2 )
     $            ndec = 0
               nw = nw - ndec
            END IF
*
*           ==== Aggressive early deflation:
*           .    split workspace under the subdiagonal into
*           .      - an nw-by-nw work array V in the lower
*           .        left-hand-corner,
*           .      - an NW-by-at-least-NW-but-more-is-better
*           .        (NW-by-NHO) horizontal work array along
*           .        the bottom edge,
*           .      - an at-least-NW-but-more-is-better (NHV-by-NW)
*           .        vertical work array along the left-hand-edge.
*           .        ====
*
            kv = n - nw + 1
            kt = nw + 1
            nho = ( n-nw-1 ) - kt + 1
            kwv = nw + 2
            nve = ( n-nw ) - kwv + 1
*
*           ==== Aggressive early deflation ====
*
            CALL zlaqr2( wantt, wantz, n, ktop, kbot, nw, h, ldh, iloz,
     $                   ihiz, z, ldz, ls, ld, w, h( kv, 1 ), ldh, nho,
     $                   h( kv, kt ), ldh, nve, h( kwv, 1 ), ldh, work,
     $                   lwork )
*
*           ==== Adjust KBOT accounting for new deflations. ====
*
            kbot = kbot - ld
*
*           ==== KS points to the shifts. ====
*
            ks = kbot - ls + 1
*
*           ==== Skip an expensive QR sweep if there is a (partly
*           .    heuristic) reason to expect that many eigenvalues
*           .    will deflate without it.  Here, the QR sweep is
*           .    skipped if many eigenvalues have just been deflated
*           .    or if the remaining active block is small.
*
            IF( ( ld.EQ.0 ) .OR. ( ( 100*ld.LE.nw*nibble ) .AND. ( kbot-
     $          ktop+1.GT.min( nmin, nwmax ) ) ) ) THEN
*
*              ==== NS = nominal number of simultaneous shifts.
*              .    This may be lowered (slightly) if ZLAQR2
*              .    did not provide that many shifts. ====
*
               ns = min( nsmax, nsr, max( 2, kbot-ktop ) )
               ns = ns - mod( ns, 2 )
*
*              ==== If there have been no deflations
*              .    in a multiple of KEXSH iterations,
*              .    then try exceptional shifts.
*              .    Otherwise use shifts provided by
*              .    ZLAQR2 above or from the eigenvalues
*              .    of a trailing principal submatrix. ====
*
               IF( mod( ndfl, kexsh ).EQ.0 ) THEN
                  ks = kbot - ns + 1
                  DO 30 i = kbot, ks + 1, -2
                     w( i ) = h( i, i ) + wilk1*cabs1( h( i, i-1 ) )
                     w( i-1 ) = w( i )
   30             CONTINUE
               ELSE
*
*                 ==== Got NS/2 or fewer shifts? Use ZLAHQR
*                 .    on a trailing principal submatrix to
*                 .    get more. (Since NS.LE.NSMAX.LE.(N-3)/6,
*                 .    there is enough space below the subdiagonal
*                 .    to fit an NS-by-NS scratch array.) ====
*
                  IF( kbot-ks+1.LE.ns / 2 ) THEN
                     ks = kbot - ns + 1
                     kt = n - ns + 1
                     CALL zlacpy( 'A', ns, ns, h( ks, ks ), ldh,
     $                            h( kt, 1 ), ldh )
                     CALL zlahqr( .false., .false., ns, 1, ns,
     $                            h( kt, 1 ), ldh, w( ks ), 1, 1, zdum,
     $                            1, inf )
                     ks = ks + inf
*
*                    ==== In case of a rare QR failure use
*                    .    eigenvalues of the trailing 2-by-2
*                    .    principal submatrix.  Scale to avoid
*                    .    overflows, underflows and subnormals.
*                    .    (The scale factor S can not be zero,
*                    .    because H(KBOT,KBOT-1) is nonzero.) ====
*
                     IF( ks.GE.kbot ) THEN
                        s = cabs1( h( kbot-1, kbot-1 ) ) +
     $                      cabs1( h( kbot, kbot-1 ) ) +
     $                      cabs1( h( kbot-1, kbot ) ) +
     $                      cabs1( h( kbot, kbot ) )
                        aa = h( kbot-1, kbot-1 ) / s
                        cc = h( kbot, kbot-1 ) / s
                        bb = h( kbot-1, kbot ) / s
                        dd = h( kbot, kbot ) / s
                        tr2 = ( aa+dd ) / two
                        det = ( aa-tr2 )*( dd-tr2 ) - bb*cc
                        rtdisc = sqrt( -det )
                        w( kbot-1 ) = ( tr2+rtdisc )*s
                        w( kbot ) = ( tr2-rtdisc )*s
*
                        ks = kbot - 1
                     END IF
                  END IF
*
                  IF( kbot-ks+1.GT.ns ) THEN
*
*                    ==== Sort the shifts (Helps a little) ====
*
                     sorted = .false.
                     DO 50 k = kbot, ks + 1, -1
                        IF( sorted )
     $                     GO TO 60
                        sorted = .true.
                        DO 40 i = ks, k - 1
                           IF( cabs1( w( i ) ).LT.cabs1( w( i+1 ) ) )
     $                          THEN
                              sorted = .false.
                              swap = w( i )
                              w( i ) = w( i+1 )
                              w( i+1 ) = swap
                           END IF
   40                   CONTINUE
   50                CONTINUE
   60                CONTINUE
                  END IF
               END IF
*
*              ==== If there are only two shifts, then use
*              .    only one.  ====
*
               IF( kbot-ks+1.EQ.2 ) THEN
                  IF( cabs1( w( kbot )-h( kbot, kbot ) ).LT.
     $                cabs1( w( kbot-1 )-h( kbot, kbot ) ) ) THEN
                     w( kbot-1 ) = w( kbot )
                  ELSE
                     w( kbot ) = w( kbot-1 )
                  END IF
               END IF
*
*              ==== Use up to NS of the the smallest magnitude
*              .    shifts.  If there aren't NS shifts available,
*              .    then use them all, possibly dropping one to
*              .    make the number of shifts even. ====
*
               ns = min( ns, kbot-ks+1 )
               ns = ns - mod( ns, 2 )
               ks = kbot - ns + 1
*
*              ==== Small-bulge multi-shift QR sweep:
*              .    split workspace under the subdiagonal into
*              .    - a KDU-by-KDU work array U in the lower
*              .      left-hand-corner,
*              .    - a KDU-by-at-least-KDU-but-more-is-better
*              .      (KDU-by-NHo) horizontal work array WH along
*              .      the bottom edge,
*              .    - and an at-least-KDU-but-more-is-better-by-KDU
*              .      (NVE-by-KDU) vertical work WV arrow along
*              .      the left-hand-edge. ====
*
               kdu = 2*ns
               ku = n - kdu + 1
               kwh = kdu + 1
               nho = ( n-kdu+1-4 ) - ( kdu+1 ) + 1
               kwv = kdu + 4
               nve = n - kdu - kwv + 1
*
*              ==== Small-bulge multi-shift QR sweep ====
*
               CALL zlaqr5( wantt, wantz, kacc22, n, ktop, kbot, ns,
     $                      w( ks ), h, ldh, iloz, ihiz, z, ldz, work,
     $                      3, h( ku, 1 ), ldh, nve, h( kwv, 1 ), ldh,
     $                      nho, h( ku, kwh ), ldh )
            END IF
*
*           ==== Note progress (or the lack of it). ====
*
            IF( ld.GT.0 ) THEN
               ndfl = 1
            ELSE
               ndfl = ndfl + 1
            END IF
*
*           ==== End of main loop ====
   70    CONTINUE
*
*        ==== Iteration limit exceeded.  Set INFO to show where
*        .    the problem occurred and exit. ====
*
         info = kbot
   80    CONTINUE
      END IF
*
*     ==== Return the optimal value of LWORK. ====
*
      work( 1 ) = dcmplx( lwkopt, 0 )
*
*     ==== End of ZLAQR4 ====
*
      END
      
! ZLAQR5
      SUBROUTINE zlaqr5( WANTT, WANTZ, KACC22, N, KTOP, KBOT, NSHFTS, S,
     $                   H, LDH, ILOZ, IHIZ, Z, LDZ, V, LDV, U, LDU, NV,
     $                   WV, LDWV, NH, WH, LDWH )
      IMPLICIT NONE
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            IHIZ, ILOZ, KACC22, KBOT, KTOP, LDH, LDU, LDV,
     $                   LDWH, LDWV, LDZ, N, NH, NSHFTS, NV
      LOGICAL            WANTT, WANTZ
*     ..
*     .. Array Arguments ..
      COMPLEX*16         H( LDH, * ), S( * ), U( LDU, * ), V( LDV, * ),
     $                   WH( LDWH, * ), WV( LDWV, * ), Z( LDZ, * )
*     ..
*
*  ================================================================
*     .. Parameters ..
      COMPLEX*16         ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0d0, 0.0d0 ),
     $                   one = ( 1.0d0, 0.0d0 ) )
      DOUBLE PRECISION   RZERO, RONE
      PARAMETER          ( RZERO = 0.0d0, rone = 1.0d0 )
*     ..
*     .. Local Scalars ..
      COMPLEX*16         ALPHA, BETA, CDUM, REFSUM, T1, T2, T3
      DOUBLE PRECISION   H11, H12, H21, H22, SAFMAX, SAFMIN, SCL,
     $                   smlnum, tst1, tst2, ulp
      INTEGER            I2, I4, INCOL, J, JBOT, JCOL, JLEN,
     $                   JROW, JTOP, K, K1, KDU, KMS, KRCOL,
     $                   m, m22, mbot, mtop, nbmps, ndcol,
     $                   ns, nu
      LOGICAL            ACCUM, BMP22
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. Intrinsic Functions ..
*
      INTRINSIC          abs, dble, dconjg, dimag, max, min, mod
*     ..
*     .. Local Arrays ..
      COMPLEX*16         VT( 3 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           dlabad, zgemm, zlacpy, zlaqr1, zlarfg, zlaset,
     $                   ztrmm
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
*     ..
*     .. Statement Function definitions ..
      cabs1( cdum ) = abs( dble( cdum ) ) + abs( dimag( cdum ) )
*     ..
*     .. Executable Statements ..
*
*     ==== If there are no shifts, then there is nothing to do. ====
*
      IF( nshfts.LT.2 )
     $   RETURN
*
*     ==== If the active block is empty or 1-by-1, then there
*     .    is nothing to do. ====
*
      IF( ktop.GE.kbot )
     $   RETURN
*
*     ==== NSHFTS is supposed to be even, but if it is odd,
*     .    then simply reduce it by one.  ====
*
      ns = nshfts - mod( nshfts, 2 )
*
*     ==== Machine constants for deflation ====
*
      safmin = dlamch( 'SAFE MINIMUM' )
      safmax = rone / safmin
      CALL dlabad( safmin, safmax )
      ulp = dlamch( 'PRECISION' )
      smlnum = safmin*( dble( n ) / ulp )
*
*     ==== Use accumulated reflections to update far-from-diagonal
*     .    entries ? ====
*
      accum = ( kacc22.EQ.1 ) .OR. ( kacc22.EQ.2 )
*
*     ==== clear trash ====
*
      IF( ktop+2.LE.kbot )
     $   h( ktop+2, ktop ) = zero
*
*     ==== NBMPS = number of 2-shift bulges in the chain ====
*
      nbmps = ns / 2
*
*     ==== KDU = width of slab ====
*
      kdu = 4*nbmps
*
*     ==== Create and chase chains of NBMPS bulges ====
*
      DO 180 incol = ktop - 2*nbmps + 1, kbot - 2, 2*nbmps
*
*        JTOP = Index from which updates from the right start.
*
         IF( accum ) THEN
            jtop = max( ktop, incol )
         ELSE IF( wantt ) THEN
            jtop = 1
         ELSE
            jtop = ktop
         END IF
*
         ndcol = incol + kdu
         IF( accum )
     $      CALL zlaset( 'ALL', kdu, kdu, zero, one, u, ldu )
*
*        ==== Near-the-diagonal bulge chase.  The following loop
*        .    performs the near-the-diagonal part of a small bulge
*        .    multi-shift QR sweep.  Each 4*NBMPS column diagonal
*        .    chunk extends from column INCOL to column NDCOL
*        .    (including both column INCOL and column NDCOL). The
*        .    following loop chases a 2*NBMPS+1 column long chain of
*        .    NBMPS bulges 2*NBMPS columns to the right.  (INCOL
*        .    may be less than KTOP and and NDCOL may be greater than
*        .    KBOT indicating phantom columns from which to chase
*        .    bulges before they are actually introduced or to which
*        .    to chase bulges beyond column KBOT.)  ====
*
         DO 145 krcol = incol, min( incol+2*nbmps-1, kbot-2 )
*
*           ==== Bulges number MTOP to MBOT are active double implicit
*           .    shift bulges.  There may or may not also be small
*           .    2-by-2 bulge, if there is room.  The inactive bulges
*           .    (if any) must wait until the active bulges have moved
*           .    down the diagonal to make room.  The phantom matrix
*           .    paradigm described above helps keep track.  ====
*
            mtop = max( 1, ( ktop-krcol ) / 2+1 )
            mbot = min( nbmps, ( kbot-krcol-1 ) / 2 )
            m22 = mbot + 1
            bmp22 = ( mbot.LT.nbmps ) .AND. ( krcol+2*( m22-1 ) ).EQ.
     $              ( kbot-2 )
*
*           ==== Generate reflections to chase the chain right
*           .    one column.  (The minimum value of K is KTOP-1.) ====
*
            IF ( bmp22 ) THEN
*
*              ==== Special case: 2-by-2 reflection at bottom treated
*              .    separately ====
*
               k = krcol + 2*( m22-1 )
               IF( k.EQ.ktop-1 ) THEN
                  CALL zlaqr1( 2, h( k+1, k+1 ), ldh, s( 2*m22-1 ),
     $                         s( 2*m22 ), v( 1, m22 ) )
                  beta = v( 1, m22 )
                  CALL zlarfg( 2, beta, v( 2, m22 ), 1, v( 1, m22 ) )
               ELSE
                  beta = h( k+1, k )
                  v( 2, m22 ) = h( k+2, k )
                  CALL zlarfg( 2, beta, v( 2, m22 ), 1, v( 1, m22 ) )
                  h( k+1, k ) = beta
                  h( k+2, k ) = zero
               END IF
 
*
*              ==== Perform update from right within 
*              .    computational window. ====
*
               t1 = v( 1, m22 )
               t2 = t1*dconjg( v( 2, m22 ) )
               DO 30 j = jtop, min( kbot, k+3 )
                  refsum = h( j, k+1 ) + v( 2, m22 )*h( j, k+2 )
                  h( j, k+1 ) = h( j, k+1 ) - refsum*t1
                  h( j, k+2 ) = h( j, k+2 ) - refsum*t2
   30          CONTINUE
*
*              ==== Perform update from left within 
*              .    computational window. ====
*
               IF( accum ) THEN
                  jbot = min( ndcol, kbot )
               ELSE IF( wantt ) THEN
                  jbot = n
               ELSE
                  jbot = kbot
               END IF
               t1 = dconjg( v( 1, m22 ) )
               t2 = t1*v( 2, m22 )
               DO 40 j = k+1, jbot
                  refsum = h( k+1, j ) +
     $                     dconjg( v( 2, m22 ) )*h( k+2, j )
                  h( k+1, j ) = h( k+1, j ) - refsum*t1
                  h( k+2, j ) = h( k+2, j ) - refsum*t2
   40          CONTINUE
*
*              ==== The following convergence test requires that
*              .    the tradition small-compared-to-nearby-diagonals
*              .    criterion and the Ahues & Tisseur (LAWN 122, 1997)
*              .    criteria both be satisfied.  The latter improves
*              .    accuracy in some examples. Falling back on an
*              .    alternate convergence criterion when TST1 or TST2
*              .    is zero (as done here) is traditional but probably
*              .    unnecessary. ====
*
               IF( k.GE.ktop ) THEN
                  IF( h( k+1, k ).NE.zero ) THEN
                     tst1 = cabs1( h( k, k ) ) + cabs1( h( k+1, k+1 ) )
                     IF( tst1.EQ.rzero ) THEN
                        IF( k.GE.ktop+1 )
     $                     tst1 = tst1 + cabs1( h( k, k-1 ) )
                        IF( k.GE.ktop+2 )
     $                     tst1 = tst1 + cabs1( h( k, k-2 ) )
                        IF( k.GE.ktop+3 )
     $                     tst1 = tst1 + cabs1( h( k, k-3 ) )
                        IF( k.LE.kbot-2 )
     $                     tst1 = tst1 + cabs1( h( k+2, k+1 ) )
                        IF( k.LE.kbot-3 )
     $                     tst1 = tst1 + cabs1( h( k+3, k+1 ) )
                        IF( k.LE.kbot-4 )
     $                     tst1 = tst1 + cabs1( h( k+4, k+1 ) )
                     END IF
                     IF( cabs1( h( k+1, k ) )
     $                   .LE.max( smlnum, ulp*tst1 ) ) THEN
                        h12 = max( cabs1( h( k+1, k ) ),
     $                     cabs1( h( k, k+1 ) ) )
                        h21 = min( cabs1( h( k+1, k ) ),
     $                     cabs1( h( k, k+1 ) ) )
                        h11 = max( cabs1( h( k+1, k+1 ) ),
     $                     cabs1( h( k, k )-h( k+1, k+1 ) ) )
                        h22 = min( cabs1( h( k+1, k+1 ) ),
     $                     cabs1( h( k, k )-h( k+1, k+1 ) ) )
                        scl = h11 + h12
                        tst2 = h22*( h11 / scl )
*
                        IF( tst2.EQ.rzero .OR. h21*( h12 / scl ).LE.
     $                      max( smlnum, ulp*tst2 ) )h( k+1, k ) = zero
                     END IF
                  END IF
               END IF
*
*              ==== Accumulate orthogonal transformations. ====
*
               IF( accum ) THEN
                  kms = k - incol
                  DO 50 j = max( 1, ktop-incol ), kdu
                     refsum = v( 1, m22 )*( u( j, kms+1 )+
     $                        v( 2, m22 )*u( j, kms+2 ) )
                     u( j, kms+1 ) = u( j, kms+1 ) - refsum
                     u( j, kms+2 ) = u( j, kms+2 ) -
     $                               refsum*dconjg( v( 2, m22 ) )
  50                 CONTINUE
               ELSE IF( wantz ) THEN
                  DO 60 j = iloz, ihiz
                     refsum = v( 1, m22 )*( z( j, k+1 )+v( 2, m22 )*
     $                        z( j, k+2 ) )
                     z( j, k+1 ) = z( j, k+1 ) - refsum
                     z( j, k+2 ) = z( j, k+2 ) -
     $                             refsum*dconjg( v( 2, m22 ) )
  60              CONTINUE
               END IF
            END IF
*
*           ==== Normal case: Chain of 3-by-3 reflections ====
*
            DO 80 m = mbot, mtop, -1
               k = krcol + 2*( m-1 )
               IF( k.EQ.ktop-1 ) THEN
                  CALL zlaqr1( 3, h( ktop, ktop ), ldh, s( 2*m-1 ),
     $                         s( 2*m ), v( 1, m ) )
                  alpha = v( 1, m )
                  CALL zlarfg( 3, alpha, v( 2, m ), 1, v( 1, m ) )
               ELSE
*
*                 ==== Perform delayed transformation of row below
*                 .    Mth bulge. Exploit fact that first two elements
*                 .    of row are actually zero. ====
*
                  refsum = v( 1, m )*v( 3, m )*h( k+3, k+2 )
                  h( k+3, k   ) = -refsum
                  h( k+3, k+1 ) = -refsum*dconjg( v( 2, m ) )
                  h( k+3, k+2 ) = h( k+3, k+2 ) -
     $                            refsum*dconjg( v( 3, m ) )
*
*                 ==== Calculate reflection to move
*                 .    Mth bulge one step. ====
*
                  beta      = h( k+1, k )
                  v( 2, m ) = h( k+2, k )
                  v( 3, m ) = h( k+3, k )
                  CALL zlarfg( 3, beta, v( 2, m ), 1, v( 1, m ) )
*
*                 ==== A Bulge may collapse because of vigilant
*                 .    deflation or destructive underflow.  In the
*                 .    underflow case, try the two-small-subdiagonals
*                 .    trick to try to reinflate the bulge.  ====
*
                  IF( h( k+3, k ).NE.zero .OR. h( k+3, k+1 ).NE.
     $                zero .OR. h( k+3, k+2 ).EQ.zero ) THEN
*
*                    ==== Typical case: not collapsed (yet). ====
*
                     h( k+1, k ) = beta
                     h( k+2, k ) = zero
                     h( k+3, k ) = zero
                  ELSE
*
*                    ==== Atypical case: collapsed.  Attempt to
*                    .    reintroduce ignoring H(K+1,K) and H(K+2,K).
*                    .    If the fill resulting from the new
*                    .    reflector is too large, then abandon it.
*                    .    Otherwise, use the new one. ====
*
                     CALL zlaqr1( 3, h( k+1, k+1 ), ldh, s( 2*m-1 ),
     $                            s( 2*m ), vt )
                     alpha = vt( 1 )
                     CALL zlarfg( 3, alpha, vt( 2 ), 1, vt( 1 ) )
                     refsum = dconjg( vt( 1 ) )*
     $                        ( h( k+1, k )+dconjg( vt( 2 ) )*
     $                        h( k+2, k ) )
*
                     IF( cabs1( h( k+2, k )-refsum*vt( 2 ) )+
     $                   cabs1( refsum*vt( 3 ) ).GT.ulp*
     $                   ( cabs1( h( k, k ) )+cabs1( h( k+1,
     $                   k+1 ) )+cabs1( h( k+2, k+2 ) ) ) ) THEN
*
*                       ==== Starting a new bulge here would
*                       .    create non-negligible fill.  Use
*                       .    the old one with trepidation. ====
*
                        h( k+1, k ) = beta
                        h( k+2, k ) = zero
                        h( k+3, k ) = zero
                     ELSE
*
*                       ==== Starting a new bulge here would
*                       .    create only negligible fill.
*                       .    Replace the old reflector with
*                       .    the new one. ====
*
                        h( k+1, k ) = h( k+1, k ) - refsum
                        h( k+2, k ) = zero
                        h( k+3, k ) = zero
                        v( 1, m ) = vt( 1 )
                        v( 2, m ) = vt( 2 )
                        v( 3, m ) = vt( 3 )
                     END IF
                  END IF
               END IF
*
*              ====  Apply reflection from the right and
*              .     the first column of update from the left.
*              .     These updates are required for the vigilant
*              .     deflation check. We still delay most of the
*              .     updates from the left for efficiency. ====
*
               t1 = v( 1, m )
               t2 = t1*dconjg( v( 2, m ) )
               t3 = t1*dconjg( v( 3, m ) )
               DO 70 j = jtop, min( kbot, k+3 )
                  refsum = h( j, k+1 ) + v( 2, m )*h( j, k+2 )
     $                     + v( 3, m )*h( j, k+3 )
                  h( j, k+1 ) = h( j, k+1 ) - refsum*t1
                  h( j, k+2 ) = h( j, k+2 ) - refsum*t2
                  h( j, k+3 ) = h( j, k+3 ) - refsum*t3
   70          CONTINUE
*
*              ==== Perform update from left for subsequent
*              .    column. ====
*
               t1 = dconjg( v( 1, m ) )
               t2 = t1*v( 2, m )
               t3 = t1*v( 3, m )
               refsum = h( k+1, k+1 )
     $                  + dconjg( v( 2, m ) )*h( k+2, k+1 )
     $                  + dconjg( v( 3, m ) )*h( k+3, k+1 )
               h( k+1, k+1 ) = h( k+1, k+1 ) - refsum*t1
               h( k+2, k+1 ) = h( k+2, k+1 ) - refsum*t2
               h( k+3, k+1 ) = h( k+3, k+1 ) - refsum*t3
*
*              ==== The following convergence test requires that
*              .    the tradition small-compared-to-nearby-diagonals
*              .    criterion and the Ahues & Tisseur (LAWN 122, 1997)
*              .    criteria both be satisfied.  The latter improves
*              .    accuracy in some examples. Falling back on an
*              .    alternate convergence criterion when TST1 or TST2
*              .    is zero (as done here) is traditional but probably
*              .    unnecessary. ====
*
               IF( k.LT.ktop)
     $              cycle
               IF( h( k+1, k ).NE.zero ) THEN
                  tst1 = cabs1( h( k, k ) ) + cabs1( h( k+1, k+1 ) )
                  IF( tst1.EQ.rzero ) THEN
                     IF( k.GE.ktop+1 )
     $                  tst1 = tst1 + cabs1( h( k, k-1 ) )
                     IF( k.GE.ktop+2 )
     $                  tst1 = tst1 + cabs1( h( k, k-2 ) )
                     IF( k.GE.ktop+3 )
     $                  tst1 = tst1 + cabs1( h( k, k-3 ) )
                     IF( k.LE.kbot-2 )
     $                  tst1 = tst1 + cabs1( h( k+2, k+1 ) )
                     IF( k.LE.kbot-3 )
     $                  tst1 = tst1 + cabs1( h( k+3, k+1 ) )
                     IF( k.LE.kbot-4 )
     $                  tst1 = tst1 + cabs1( h( k+4, k+1 ) )
                  END IF
                  IF( cabs1( h( k+1, k ) ).LE.max( smlnum, ulp*tst1 ) )
     $                 THEN
                     h12 = max( cabs1( h( k+1, k ) ),
     $                     cabs1( h( k, k+1 ) ) )
                     h21 = min( cabs1( h( k+1, k ) ),
     $                     cabs1( h( k, k+1 ) ) )
                     h11 = max( cabs1( h( k+1, k+1 ) ),
     $                     cabs1( h( k, k )-h( k+1, k+1 ) ) )
                     h22 = min( cabs1( h( k+1, k+1 ) ),
     $                     cabs1( h( k, k )-h( k+1, k+1 ) ) )
                     scl = h11 + h12
                     tst2 = h22*( h11 / scl )
*
                     IF( tst2.EQ.rzero .OR. h21*( h12 / scl ).LE.
     $                   max( smlnum, ulp*tst2 ) )h( k+1, k ) = zero
                  END IF
               END IF
   80       CONTINUE
*
*           ==== Multiply H by reflections from the left ====
*
            IF( accum ) THEN
               jbot = min( ndcol, kbot )
            ELSE IF( wantt ) THEN
               jbot = n
            ELSE
               jbot = kbot
            END IF
*
            DO 100 m = mbot, mtop, -1
               k = krcol + 2*( m-1 )
               t1 = dconjg( v( 1, m ) )
               t2 = t1*v( 2, m )
               t3 = t1*v( 3, m )
               DO 90 j = max( ktop, krcol + 2*m ), jbot
                  refsum = h( k+1, j ) + dconjg( v( 2, m ) )*h( k+2, j )
     $                     + dconjg( v( 3, m ) )*h( k+3, j )
                  h( k+1, j ) = h( k+1, j ) - refsum*t1
                  h( k+2, j ) = h( k+2, j ) - refsum*t2
                  h( k+3, j ) = h( k+3, j ) - refsum*t3
   90          CONTINUE
  100       CONTINUE
*
*           ==== Accumulate orthogonal transformations. ====
*
            IF( accum ) THEN
*
*              ==== Accumulate U. (If needed, update Z later
*              .    with an efficient matrix-matrix
*              .    multiply.) ====
*
               DO 120 m = mbot, mtop, -1
                  k = krcol + 2*( m-1 )
                  kms = k - incol
                  i2 = max( 1, ktop-incol )
                  i2 = max( i2, kms-(krcol-incol)+1 )
                  i4 = min( kdu, krcol + 2*( mbot-1 ) - incol + 5 )
                  t1 = v( 1, m )
                  t2 = t1*dconjg( v( 2, m ) )
                  t3 = t1*dconjg( v( 3, m ) )
                  DO 110 j = i2, i4
                     refsum = u( j, kms+1 ) + v( 2, m )*u( j, kms+2 )
     $                        + v( 3, m )*u( j, kms+3 )
                     u( j, kms+1 ) = u( j, kms+1 ) - refsum*t1
                     u( j, kms+2 ) = u( j, kms+2 ) - refsum*t2
                     u( j, kms+3 ) = u( j, kms+3 ) - refsum*t3
  110             CONTINUE
  120          CONTINUE
            ELSE IF( wantz ) THEN
*
*              ==== U is not accumulated, so update Z
*              .    now by multiplying by reflections
*              .    from the right. ====
*
               DO 140 m = mbot, mtop, -1
                  k = krcol + 2*( m-1 )
                  t1 = v( 1, m )
                  t2 = t1*dconjg( v( 2, m ) )
                  t3 = t1*dconjg( v( 3, m ) )
                  DO 130 j = iloz, ihiz
                     refsum = z( j, k+1 ) + v( 2, m )*z( j, k+2 )
     $                        + v( 3, m )*z( j, k+3 )
                     z( j, k+1 ) = z( j, k+1 ) - refsum*t1
                     z( j, k+2 ) = z( j, k+2 ) - refsum*t2
                     z( j, k+3 ) = z( j, k+3 ) - refsum*t3
  130             CONTINUE
  140          CONTINUE
            END IF
*
*           ==== End of near-the-diagonal bulge chase. ====
*
  145    CONTINUE
*
*        ==== Use U (if accumulated) to update far-from-diagonal
*        .    entries in H.  If required, use U to update Z as
*        .    well. ====
*
         IF( accum ) THEN
            IF( wantt ) THEN
               jtop = 1
               jbot = n
            ELSE
               jtop = ktop
               jbot = kbot
            END IF
            k1 = max( 1, ktop-incol )
            nu = ( kdu-max( 0, ndcol-kbot ) ) - k1 + 1
*
*           ==== Horizontal Multiply ====
*
            DO 150 jcol = min( ndcol, kbot ) + 1, jbot, nh
               jlen = min( nh, jbot-jcol+1 )
               CALL zgemm( 'C', 'N', nu, jlen, nu, one, u( k1, k1 ),
     $                     ldu, h( incol+k1, jcol ), ldh, zero, wh,
     $                     ldwh )
               CALL zlacpy( 'ALL', nu, jlen, wh, ldwh,
     $                      h( incol+k1, jcol ), ldh )
  150       CONTINUE
*
*           ==== Vertical multiply ====
*
            DO 160 jrow = jtop, max( ktop, incol ) - 1, nv
               jlen = min( nv, max( ktop, incol )-jrow )
               CALL zgemm( 'N', 'N', jlen, nu, nu, one,
     $                     h( jrow, incol+k1 ), ldh, u( k1, k1 ),
     $                     ldu, zero, wv, ldwv )
               CALL zlacpy( 'ALL', jlen, nu, wv, ldwv,
     $                      h( jrow, incol+k1 ), ldh )
  160       CONTINUE
*
*           ==== Z multiply (also vertical) ====
*
            IF( wantz ) THEN
               DO 170 jrow = iloz, ihiz, nv
                  jlen = min( nv, ihiz-jrow+1 )
                  CALL zgemm( 'N', 'N', jlen, nu, nu, one,
     $                        z( jrow, incol+k1 ), ldz, u( k1, k1 ),
     $                        ldu, zero, wv, ldwv )
                  CALL zlacpy( 'ALL', jlen, nu, wv, ldwv,
     $                         z( jrow, incol+k1 ), ldz )
  170          CONTINUE
            END IF
         END IF
  180 CONTINUE
*
*     ==== End of ZLAQR5 ====
*
      END
      
! ZLARFG
      SUBROUTINE zlarfg( N, ALPHA, X, INCX, TAU )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
      COMPLEX*16         ALPHA, TAU
*     ..
*     .. Array Arguments ..
      COMPLEX*16         X( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      parameter( one = 1.0d+0, zero = 0.0d+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            J, KNT
      DOUBLE PRECISION   ALPHI, ALPHR, BETA, RSAFMN, SAFMIN, XNORM
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLAPY3, DZNRM2
      COMPLEX*16         ZLADIV
      EXTERNAL           dlamch, dlapy3, dznrm2, zladiv
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, dble, dcmplx, dimag, sign
*     ..
*     .. External Subroutines ..
      EXTERNAL           zdscal, zscal
*     ..
*     .. Executable Statements ..
*
      IF( n.LE.0 ) THEN
         tau = zero
         RETURN
      END IF
*
      xnorm = dznrm2( n-1, x, incx )
      alphr = dble( alpha )
      alphi = dimag( alpha )
*
      IF( xnorm.EQ.zero .AND. alphi.EQ.zero ) THEN
*
*        H  =  I
*
         tau = zero
      ELSE
*
*        general case
*
         beta = -sign( dlapy3( alphr, alphi, xnorm ), alphr )
         safmin = dlamch( 'S' ) / dlamch( 'E' )
         rsafmn = one / safmin
*
         knt = 0
         IF( abs( beta ).LT.safmin ) THEN
*
*           XNORM, BETA may be inaccurate; scale X and recompute them
*
   10       CONTINUE
            knt = knt + 1
            CALL zdscal( n-1, rsafmn, x, incx )
            beta = beta*rsafmn
            alphi = alphi*rsafmn
            alphr = alphr*rsafmn
            IF( (abs( beta ).LT.safmin) .AND. (knt .LT. 20) )
     $         GO TO 10
*
*           New BETA is at most 1, at least SAFMIN
*
            xnorm = dznrm2( n-1, x, incx )
            alpha = dcmplx( alphr, alphi )
            beta = -sign( dlapy3( alphr, alphi, xnorm ), alphr )
         END IF
         tau = dcmplx( ( beta-alphr ) / beta, -alphi / beta )
         alpha = zladiv( dcmplx( one ), alpha-beta )
         CALL zscal( n-1, alpha, x, incx )
*
*        If ALPHA is subnormal, it may lose relative accuracy
*
         DO 20 j = 1, knt
            beta = beta*safmin
 20      CONTINUE
         alpha = beta
      END IF
*
      RETURN
*
*     End of ZLARFG
*
      END
      
! ZLACGV
      SUBROUTINE zlacgv( N, X, INCX )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         X( * )
*     ..
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, IOFF
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          dconjg
*     ..
*     .. Executable Statements ..
*
      IF( incx.EQ.1 ) THEN
         DO 10 i = 1, n
            x( i ) = dconjg( x( i ) )
   10    CONTINUE
      ELSE
         ioff = 1
         IF( incx.LT.0 )
     $      ioff = 1 - ( n-1 )*incx
         DO 20 i = 1, n
            x( ioff ) = dconjg( x( ioff ) )
            ioff = ioff + incx
   20    CONTINUE
      END IF
      RETURN
*
*     End of ZLACGV
*
      END
      
! ZTRMV
      SUBROUTINE ztrmv(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
*
*  -- Reference BLAS level2 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER INCX,LDA,N
      CHARACTER DIAG,TRANS,UPLO
*     ..
*     .. Array Arguments ..
      COMPLEX*16 A(LDA,*),X(*)
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16 ZERO
      parameter(zero= (0.0d+0,0.0d+0))
*     ..
*     .. Local Scalars ..
      COMPLEX*16 TEMP
      INTEGER I,INFO,IX,J,JX,KX
      LOGICAL NOCONJ,NOUNIT
*     ..
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL lsame
*     ..
*     .. External Subroutines ..
      EXTERNAL xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC dconjg,max
*     ..
*
*     Test the input parameters.
*
      info = 0
      IF (.NOT.lsame(uplo,'U') .AND. .NOT.lsame(uplo,'L')) THEN
          info = 1
      ELSE IF (.NOT.lsame(trans,'N') .AND. .NOT.lsame(trans,'T') .AND.
     +         .NOT.lsame(trans,'C')) THEN
          info = 2
      ELSE IF (.NOT.lsame(diag,'U') .AND. .NOT.lsame(diag,'N')) THEN
          info = 3
      ELSE IF (n.LT.0) THEN
          info = 4
      ELSE IF (lda.LT.max(1,n)) THEN
          info = 6
      ELSE IF (incx.EQ.0) THEN
          info = 8
      END IF
      IF (info.NE.0) THEN
          CALL xerbla('ZTRMV ',info)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF (n.EQ.0) RETURN
*
      noconj = lsame(trans,'T')
      nounit = lsame(diag,'N')
*
*     Set up the start point in X if the increment is not unity. This
*     will be  ( N - 1 )*INCX  too small for descending loops.
*
      IF (incx.LE.0) THEN
          kx = 1 - (n-1)*incx
      ELSE IF (incx.NE.1) THEN
          kx = 1
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF (lsame(trans,'N')) THEN
*
*        Form  x := A*x.
*
          IF (lsame(uplo,'U')) THEN
              IF (incx.EQ.1) THEN
                  DO 20 j = 1,n
                      IF (x(j).NE.zero) THEN
                          temp = x(j)
                          DO 10 i = 1,j - 1
                              x(i) = x(i) + temp*a(i,j)
   10                     CONTINUE
                          IF (nounit) x(j) = x(j)*a(j,j)
                      END IF
   20             CONTINUE
              ELSE
                  jx = kx
                  DO 40 j = 1,n
                      IF (x(jx).NE.zero) THEN
                          temp = x(jx)
                          ix = kx
                          DO 30 i = 1,j - 1
                              x(ix) = x(ix) + temp*a(i,j)
                              ix = ix + incx
   30                     CONTINUE
                          IF (nounit) x(jx) = x(jx)*a(j,j)
                      END IF
                      jx = jx + incx
   40             CONTINUE
              END IF
          ELSE
              IF (incx.EQ.1) THEN
                  DO 60 j = n,1,-1
                      IF (x(j).NE.zero) THEN
                          temp = x(j)
                          DO 50 i = n,j + 1,-1
                              x(i) = x(i) + temp*a(i,j)
   50                     CONTINUE
                          IF (nounit) x(j) = x(j)*a(j,j)
                      END IF
   60             CONTINUE
              ELSE
                  kx = kx + (n-1)*incx
                  jx = kx
                  DO 80 j = n,1,-1
                      IF (x(jx).NE.zero) THEN
                          temp = x(jx)
                          ix = kx
                          DO 70 i = n,j + 1,-1
                              x(ix) = x(ix) + temp*a(i,j)
                              ix = ix - incx
   70                     CONTINUE
                          IF (nounit) x(jx) = x(jx)*a(j,j)
                      END IF
                      jx = jx - incx
   80             CONTINUE
              END IF
          END IF
      ELSE
*
*        Form  x := A**T*x  or  x := A**H*x.
*
          IF (lsame(uplo,'U')) THEN
              IF (incx.EQ.1) THEN
                  DO 110 j = n,1,-1
                      temp = x(j)
                      IF (noconj) THEN
                          IF (nounit) temp = temp*a(j,j)
                          DO 90 i = j - 1,1,-1
                              temp = temp + a(i,j)*x(i)
   90                     CONTINUE
                      ELSE
                          IF (nounit) temp = temp*dconjg(a(j,j))
                          DO 100 i = j - 1,1,-1
                              temp = temp + dconjg(a(i,j))*x(i)
  100                     CONTINUE
                      END IF
                      x(j) = temp
  110             CONTINUE
              ELSE
                  jx = kx + (n-1)*incx
                  DO 140 j = n,1,-1
                      temp = x(jx)
                      ix = jx
                      IF (noconj) THEN
                          IF (nounit) temp = temp*a(j,j)
                          DO 120 i = j - 1,1,-1
                              ix = ix - incx
                              temp = temp + a(i,j)*x(ix)
  120                     CONTINUE
                      ELSE
                          IF (nounit) temp = temp*dconjg(a(j,j))
                          DO 130 i = j - 1,1,-1
                              ix = ix - incx
                              temp = temp + dconjg(a(i,j))*x(ix)
  130                     CONTINUE
                      END IF
                      x(jx) = temp
                      jx = jx - incx
  140             CONTINUE
              END IF
          ELSE
              IF (incx.EQ.1) THEN
                  DO 170 j = 1,n
                      temp = x(j)
                      IF (noconj) THEN
                          IF (nounit) temp = temp*a(j,j)
                          DO 150 i = j + 1,n
                              temp = temp + a(i,j)*x(i)
  150                     CONTINUE
                      ELSE
                          IF (nounit) temp = temp*dconjg(a(j,j))
                          DO 160 i = j + 1,n
                              temp = temp + dconjg(a(i,j))*x(i)
  160                     CONTINUE
                      END IF
                      x(j) = temp
  170             CONTINUE
              ELSE
                  jx = kx
                  DO 200 j = 1,n
                      temp = x(jx)
                      ix = jx
                      IF (noconj) THEN
                          IF (nounit) temp = temp*a(j,j)
                          DO 180 i = j + 1,n
                              ix = ix + incx
                              temp = temp + a(i,j)*x(ix)
  180                     CONTINUE
                      ELSE
                          IF (nounit) temp = temp*dconjg(a(j,j))
                          DO 190 i = j + 1,n
                              ix = ix + incx
                              temp = temp + dconjg(a(i,j))*x(ix)
  190                     CONTINUE
                      END IF
                      x(jx) = temp
                      jx = jx + incx
  200             CONTINUE
              END IF
          END IF
      END IF
*
      RETURN
*
*     End of ZTRMV
*
      END

! ZLARF
      SUBROUTINE zlarf( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE
      INTEGER            INCV, LDC, M, N
      COMPLEX*16         TAU
*     ..
*     .. Array Arguments ..
      COMPLEX*16         C( LDC, * ), V( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      parameter( one = ( 1.0d+0, 0.0d+0 ),
     $                   zero = ( 0.0d+0, 0.0d+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            APPLYLEFT
      INTEGER            I, LASTV, LASTC
*     ..
*     .. External Subroutines ..
      EXTERNAL           zgemv, zgerc
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAZLR, ILAZLC
      EXTERNAL           lsame, ilazlr, ilazlc
*     ..
*     .. Executable Statements ..
*
      applyleft = lsame( side, 'L' )
      lastv = 0
      lastc = 0
      IF( tau.NE.zero ) THEN
*     Set up variables for scanning V.  LASTV begins pointing to the end
*     of V.
         IF( applyleft ) THEN
            lastv = m
         ELSE
            lastv = n
         END IF
         IF( incv.GT.0 ) THEN
            i = 1 + (lastv-1) * incv
         ELSE
            i = 1
         END IF
*     Look for the last non-zero row in V.
         DO WHILE( lastv.GT.0 .AND. v( i ).EQ.zero )
            lastv = lastv - 1
            i = i - incv
         END DO
         IF( applyleft ) THEN
*     Scan for the last non-zero column in C(1:lastv,:).
            lastc = ilazlc(lastv, n, c, ldc)
         ELSE
*     Scan for the last non-zero row in C(:,1:lastv).
            lastc = ilazlr(m, lastv, c, ldc)
         END IF
      END IF
*     Note that lastc.eq.0 renders the BLAS operations null; no special
*     case is needed at this level.
      IF( applyleft ) THEN
*
*        Form  H * C
*
         IF( lastv.GT.0 ) THEN
*
*           w(1:lastc,1) := C(1:lastv,1:lastc)**H * v(1:lastv,1)
*
            CALL zgemv( 'Conjugate transpose', lastv, lastc, one,
     $           c, ldc, v, incv, zero, work, 1 )
*
*           C(1:lastv,1:lastc) := C(...) - v(1:lastv,1) * w(1:lastc,1)**H
*
            CALL zgerc( lastv, lastc, -tau, v, incv, work, 1, c, ldc )
         END IF
      ELSE
*
*        Form  C * H
*
         IF( lastv.GT.0 ) THEN
*
*           w(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1)
*
            CALL zgemv( 'No transpose', lastc, lastv, one, c, ldc,
     $           v, incv, zero, work, 1 )
*
*           C(1:lastc,1:lastv) := C(...) - w(1:lastc,1) * v(1:lastv,1)**H
*
            CALL zgerc( lastc, lastv, -tau, work, 1, v, incv, c, ldc )
         END IF
      END IF
      RETURN
*
*     End of ZLARF
*
      END
      
! ZUNG2R
      SUBROUTINE zung2r( M, N, K, A, LDA, TAU, WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      parameter( one = ( 1.0d+0, 0.0d+0 ),
     $                   zero = ( 0.0d+0, 0.0d+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, L
*     ..
*     .. External Subroutines ..
      EXTERNAL           xerbla, zlarf, zscal
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          max
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      info = 0
      IF( m.LT.0 ) THEN
         info = -1
      ELSE IF( n.LT.0 .OR. n.GT.m ) THEN
         info = -2
      ELSE IF( k.LT.0 .OR. k.GT.n ) THEN
         info = -3
      ELSE IF( lda.LT.max( 1, m ) ) THEN
         info = -5
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'ZUNG2R', -info )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( n.LE.0 )
     $   RETURN
*
*     Initialise columns k+1:n to columns of the unit matrix
*
      DO 20 j = k + 1, n
         DO 10 l = 1, m
            a( l, j ) = zero
   10    CONTINUE
         a( j, j ) = one
   20 CONTINUE
*
      DO 40 i = k, 1, -1
*
*        Apply H(i) to A(i:m,i:n) from the left
*
         IF( i.LT.n ) THEN
            a( i, i ) = one
            CALL zlarf( 'Left', m-i+1, n-i, a( i, i ), 1, tau( i ),
     $                  a( i, i+1 ), lda, work )
         END IF
         IF( i.LT.m )
     $      CALL zscal( m-i, -tau( i ), a( i+1, i ), 1 )
         a( i, i ) = one - tau( i )
*
*        Set A(1:i-1,i) to zero
*
         DO 30 l = 1, i - 1
            a( l, i ) = zero
   30    CONTINUE
   40 CONTINUE
      RETURN
*
*     End of ZUNG2R
*
      END
      
! ZLARFT
      SUBROUTINE zlarft( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          DIRECT, STOREV
      INTEGER            K, LDT, LDV, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         T( LDT, * ), TAU( * ), V( LDV, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      parameter( one = ( 1.0d+0, 0.0d+0 ),
     $                   zero = ( 0.0d+0, 0.0d+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, PREVLASTV, LASTV
*     ..
*     .. External Subroutines ..
      EXTERNAL           zgemv, ztrmv, zgemm
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           lsame
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( n.EQ.0 )
     $   RETURN
*
      IF( lsame( direct, 'F' ) ) THEN
         prevlastv = n
         DO i = 1, k
            prevlastv = max( prevlastv, i )
            IF( tau( i ).EQ.zero ) THEN
*
*              H(i)  =  I
*
               DO j = 1, i
                  t( j, i ) = zero
               END DO
            ELSE
*
*              general case
*
               IF( lsame( storev, 'C' ) ) THEN
*                 Skip any trailing zeros.
                  DO lastv = n, i+1, -1
                     IF( v( lastv, i ).NE.zero ) EXIT
                  END DO
                  DO j = 1, i-1
                     t( j, i ) = -tau( i ) * conjg( v( i , j ) )
                  END DO
                  j = min( lastv, prevlastv )
*
*                 T(1:i-1,i) := - tau(i) * V(i:j,1:i-1)**H * V(i:j,i)
*
                  CALL zgemv( 'Conjugate transpose', j-i, i-1,
     $                        -tau( i ), v( i+1, 1 ), ldv,
     $                        v( i+1, i ), 1, one, t( 1, i ), 1 )
               ELSE
*                 Skip any trailing zeros.
                  DO lastv = n, i+1, -1
                     IF( v( i, lastv ).NE.zero ) EXIT
                  END DO
                  DO j = 1, i-1
                     t( j, i ) = -tau( i ) * v( j , i )
                  END DO
                  j = min( lastv, prevlastv )
*
*                 T(1:i-1,i) := - tau(i) * V(1:i-1,i:j) * V(i,i:j)**H
*
                  CALL zgemm( 'N', 'C', i-1, 1, j-i, -tau( i ),
     $                        v( 1, i+1 ), ldv, v( i, i+1 ), ldv,
     $                        one, t( 1, i ), ldt )
               END IF
*
*              T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i)
*
               CALL ztrmv( 'Upper', 'No transpose', 'Non-unit', i-1, t,
     $                     ldt, t( 1, i ), 1 )
               t( i, i ) = tau( i )
               IF( i.GT.1 ) THEN
                  prevlastv = max( prevlastv, lastv )
               ELSE
                  prevlastv = lastv
               END IF
             END IF
         END DO
      ELSE
         prevlastv = 1
         DO i = k, 1, -1
            IF( tau( i ).EQ.zero ) THEN
*
*              H(i)  =  I
*
               DO j = i, k
                  t( j, i ) = zero
               END DO
            ELSE
*
*              general case
*
               IF( i.LT.k ) THEN
                  IF( lsame( storev, 'C' ) ) THEN
*                    Skip any leading zeros.
                     DO lastv = 1, i-1
                        IF( v( lastv, i ).NE.zero ) EXIT
                     END DO
                     DO j = i+1, k
                        t( j, i ) = -tau( i ) * conjg( v( n-k+i , j ) )
                     END DO
                     j = max( lastv, prevlastv )
*
*                    T(i+1:k,i) = -tau(i) * V(j:n-k+i,i+1:k)**H * V(j:n-k+i,i)
*
                     CALL zgemv( 'Conjugate transpose', n-k+i-j, k-i,
     $                           -tau( i ), v( j, i+1 ), ldv, v( j, i ),
     $                           1, one, t( i+1, i ), 1 )
                  ELSE
*                    Skip any leading zeros.
                     DO lastv = 1, i-1
                        IF( v( i, lastv ).NE.zero ) EXIT
                     END DO
                     DO j = i+1, k
                        t( j, i ) = -tau( i ) * v( j, n-k+i )
                     END DO
                     j = max( lastv, prevlastv )
*
*                    T(i+1:k,i) = -tau(i) * V(i+1:k,j:n-k+i) * V(i,j:n-k+i)**H
*
                     CALL zgemm( 'N', 'C', k-i, 1, n-k+i-j, -tau( i ),
     $                           v( i+1, j ), ldv, v( i, j ), ldv,
     $                           one, t( i+1, i ), ldt )
                  END IF
*
*                 T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i)
*
                  CALL ztrmv( 'Lower', 'No transpose', 'Non-unit', k-i,
     $                        t( i+1, i+1 ), ldt, t( i+1, i ), 1 )
                  IF( i.GT.1 ) THEN
                     prevlastv = min( prevlastv, lastv )
                  ELSE
                     prevlastv = lastv
                  END IF
               END IF
               t( i, i ) = tau( i )
            END IF
         END DO
      END IF
      RETURN
*
*     End of ZLARFT
*
      END

! ZUNMHR
      SUBROUTINE zunmhr( SIDE, TRANS, M, N, ILO, IHI, A, LDA, TAU, C,
     $                   LDC, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            IHI, ILO, INFO, LDA, LDC, LWORK, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            LEFT, LQUERY
      INTEGER            I1, I2, IINFO, LWKOPT, MI, NB, NH, NI, NQ, NW
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           lsame, ilaenv
*     ..
*     .. External Subroutines ..
      EXTERNAL           xerbla, zunmqr
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          max, min
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      info = 0
      nh = ihi - ilo
      left = lsame( side, 'L' )
      lquery = ( lwork.EQ.-1 )
*
*     NQ is the order of Q and NW is the minimum dimension of WORK
*
      IF( left ) THEN
         nq = m
         nw = max( 1, n )
      ELSE
         nq = n
         nw = max( 1, m )
      END IF
      IF( .NOT.left .AND. .NOT.lsame( side, 'R' ) ) THEN
         info = -1
      ELSE IF( .NOT.lsame( trans, 'N' ) .AND. .NOT.lsame( trans, 'C' ) )
     $          THEN
         info = -2
      ELSE IF( m.LT.0 ) THEN
         info = -3
      ELSE IF( n.LT.0 ) THEN
         info = -4
      ELSE IF( ilo.LT.1 .OR. ilo.GT.max( 1, nq ) ) THEN
         info = -5
      ELSE IF( ihi.LT.min( ilo, nq ) .OR. ihi.GT.nq ) THEN
         info = -6
      ELSE IF( lda.LT.max( 1, nq ) ) THEN
         info = -8
      ELSE IF( ldc.LT.max( 1, m ) ) THEN
         info = -11
      ELSE IF( lwork.LT.nw .AND. .NOT.lquery ) THEN
         info = -13
      END IF
*
      IF( info.EQ.0 ) THEN
         IF( left ) THEN
            nb = ilaenv( 1, 'ZUNMQR', side // trans, nh, n, nh, -1 )
         ELSE
            nb = ilaenv( 1, 'ZUNMQR', side // trans, m, nh, nh, -1 )
         END IF
         lwkopt = nw*nb
         work( 1 ) = lwkopt
      END IF
*
      IF( info.NE.0 ) THEN
         CALL xerbla( 'ZUNMHR', -info )
         RETURN
      ELSE IF( lquery ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( m.EQ.0 .OR. n.EQ.0 .OR. nh.EQ.0 ) THEN
         work( 1 ) = 1
         RETURN
      END IF
*
      IF( left ) THEN
         mi = nh
         ni = n
         i1 = ilo + 1
         i2 = 1
      ELSE
         mi = m
         ni = nh
         i1 = 1
         i2 = ilo + 1
      END IF
*
      CALL zunmqr( side, trans, mi, ni, nh, a( ilo+1, ilo ), lda,
     $             tau( ilo ), c( i1, i2 ), ldc, work, lwork, iinfo )
*
      work( 1 ) = lwkopt
      RETURN
*
*     End of ZUNMHR
*
      END
      
! ZTREXC
      SUBROUTINE ztrexc( COMPQ, N, T, LDT, Q, LDQ, IFST, ILST, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          COMPQ
      INTEGER            IFST, ILST, INFO, LDQ, LDT, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         Q( LDQ, * ), T( LDT, * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            WANTQ
      INTEGER            K, M1, M2, M3
      DOUBLE PRECISION   CS
      COMPLEX*16         SN, T11, T22, TEMP
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           lsame
*     ..
*     .. External Subroutines ..
      EXTERNAL           xerbla, zlartg, zrot
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          dconjg, max
*     ..
*     .. Executable Statements ..
*
*     Decode and test the input parameters.
*
      info = 0
      wantq = lsame( compq, 'V' )
      IF( .NOT.lsame( compq, 'N' ) .AND. .NOT.wantq ) THEN
         info = -1
      ELSE IF( n.LT.0 ) THEN
         info = -2
      ELSE IF( ldt.LT.max( 1, n ) ) THEN
         info = -4
      ELSE IF( ldq.LT.1 .OR. ( wantq .AND. ldq.LT.max( 1, n ) ) ) THEN
         info = -6
      ELSE IF(( ifst.LT.1 .OR. ifst.GT.n ).AND.( n.GT.0 )) THEN
         info = -7
      ELSE IF(( ilst.LT.1 .OR. ilst.GT.n ).AND.( n.GT.0 )) THEN
         info = -8
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'ZTREXC', -info )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( n.LE.1 .OR. ifst.EQ.ilst )
     $   RETURN
*
      IF( ifst.LT.ilst ) THEN
*
*        Move the IFST-th diagonal element forward down the diagonal.
*
         m1 = 0
         m2 = -1
         m3 = 1
      ELSE
*
*        Move the IFST-th diagonal element backward up the diagonal.
*
         m1 = -1
         m2 = 0
         m3 = -1
      END IF
*
      DO 10 k = ifst + m1, ilst + m2, m3
*
*        Interchange the k-th and (k+1)-th diagonal elements.
*
         t11 = t( k, k )
         t22 = t( k+1, k+1 )
*
*        Determine the transformation to perform the interchange.
*
         CALL zlartg( t( k, k+1 ), t22-t11, cs, sn, temp )
*
*        Apply transformation to the matrix T.
*
         IF( k+2.LE.n )
     $      CALL zrot( n-k-1, t( k, k+2 ), ldt, t( k+1, k+2 ), ldt, cs,
     $                 sn )
         CALL zrot( k-1, t( 1, k ), 1, t( 1, k+1 ), 1, cs,
     $              dconjg( sn ) )
*
         t( k, k ) = t22
         t( k+1, k+1 ) = t11
*
         IF( wantq ) THEN
*
*           Accumulate transformation in the matrix Q.
*
            CALL zrot( n, q( 1, k ), 1, q( 1, k+1 ), 1, cs,
     $                 dconjg( sn ) )
         END IF
*
   10 CONTINUE
*
      RETURN
*
*     End of ZTREXC
*
      END
      
! ZLAQR2
      SUBROUTINE zlaqr2( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ,
     $                   IHIZ, Z, LDZ, NS, ND, SH, V, LDV, NH, T, LDT,
     $                   NV, WV, LDWV, WORK, LWORK )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            IHIZ, ILOZ, KBOT, KTOP, LDH, LDT, LDV, LDWV,
     $                   LDZ, LWORK, N, ND, NH, NS, NV, NW
      LOGICAL            WANTT, WANTZ
*     ..
*     .. Array Arguments ..
      COMPLEX*16         H( LDH, * ), SH( * ), T( LDT, * ), V( LDV, * ),
     $                   WORK( * ), WV( LDWV, * ), Z( LDZ, * )
*     ..
*
*  ================================================================
*
*     .. Parameters ..
      COMPLEX*16         ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0d0, 0.0d0 ),
     $                   one = ( 1.0d0, 0.0d0 ) )
      DOUBLE PRECISION   RZERO, RONE
      PARAMETER          ( RZERO = 0.0d0, rone = 1.0d0 )
*     ..
*     .. Local Scalars ..
      COMPLEX*16         BETA, CDUM, S, TAU
      DOUBLE PRECISION   FOO, SAFMAX, SAFMIN, SMLNUM, ULP
      INTEGER            I, IFST, ILST, INFO, INFQR, J, JW, KCOL, KLN,
     $                   knt, krow, kwtop, ltop, lwk1, lwk2, lwkopt
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           dlabad, zcopy, zgehrd, zgemm, zlacpy, zlahqr,
     $                   zlarf, zlarfg, zlaset, ztrexc, zunmhr
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, dble, dcmplx, dconjg, dimag, int, max, min
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
*     ..
*     .. Statement Function definitions ..
      cabs1( cdum ) = abs( dble( cdum ) ) + abs( dimag( cdum ) )
*     ..
*     .. Executable Statements ..
*
*     ==== Estimate optimal workspace. ====
*
      jw = min( nw, kbot-ktop+1 )
      IF( jw.LE.2 ) THEN
         lwkopt = 1
      ELSE
*
*        ==== Workspace query call to ZGEHRD ====
*
         CALL zgehrd( jw, 1, jw-1, t, ldt, work, work, -1, info )
         lwk1 = int( work( 1 ) )
*
*        ==== Workspace query call to ZUNMHR ====
*
         CALL zunmhr( 'R', 'N', jw, jw, 1, jw-1, t, ldt, work, v, ldv,
     $                work, -1, info )
         lwk2 = int( work( 1 ) )
*
*        ==== Optimal workspace ====
*
         lwkopt = jw + max( lwk1, lwk2 )
      END IF
*
*     ==== Quick return in case of workspace query. ====
*
      IF( lwork.EQ.-1 ) THEN
         work( 1 ) = dcmplx( lwkopt, 0 )
         RETURN
      END IF
*
*     ==== Nothing to do ...
*     ... for an empty active block ... ====
      ns = 0
      nd = 0
      work( 1 ) = one
      IF( ktop.GT.kbot )
     $   RETURN
*     ... nor for an empty deflation window. ====
      IF( nw.LT.1 )
     $   RETURN
*
*     ==== Machine constants ====
*
      safmin = dlamch( 'SAFE MINIMUM' )
      safmax = rone / safmin
      CALL dlabad( safmin, safmax )
      ulp = dlamch( 'PRECISION' )
      smlnum = safmin*( dble( n ) / ulp )
*
*     ==== Setup deflation window ====
*
      jw = min( nw, kbot-ktop+1 )
      kwtop = kbot - jw + 1
      IF( kwtop.EQ.ktop ) THEN
         s = zero
      ELSE
         s = h( kwtop, kwtop-1 )
      END IF
*
      IF( kbot.EQ.kwtop ) THEN
*
*        ==== 1-by-1 deflation window: not much to do ====
*
         sh( kwtop ) = h( kwtop, kwtop )
         ns = 1
         nd = 0
         IF( cabs1( s ).LE.max( smlnum, ulp*cabs1( h( kwtop,
     $       kwtop ) ) ) ) THEN
            ns = 0
            nd = 1
            IF( kwtop.GT.ktop )
     $         h( kwtop, kwtop-1 ) = zero
         END IF
         work( 1 ) = one
         RETURN
      END IF
*
*     ==== Convert to spike-triangular form.  (In case of a
*     .    rare QR failure, this routine continues to do
*     .    aggressive early deflation using that part of
*     .    the deflation window that converged using INFQR
*     .    here and there to keep track.) ====
*
      CALL zlacpy( 'U', jw, jw, h( kwtop, kwtop ), ldh, t, ldt )
      CALL zcopy( jw-1, h( kwtop+1, kwtop ), ldh+1, t( 2, 1 ), ldt+1 )
*
      CALL zlaset( 'A', jw, jw, zero, one, v, ldv )
      CALL zlahqr( .true., .true., jw, 1, jw, t, ldt, sh( kwtop ), 1,
     $             jw, v, ldv, infqr )
*
*     ==== Deflation detection loop ====
*
      ns = jw
      ilst = infqr + 1
      DO 10 knt = infqr + 1, jw
*
*        ==== Small spike tip deflation test ====
*
         foo = cabs1( t( ns, ns ) )
         IF( foo.EQ.rzero )
     $      foo = cabs1( s )
         IF( cabs1( s )*cabs1( v( 1, ns ) ).LE.max( smlnum, ulp*foo ) )
     $        THEN
*
*           ==== One more converged eigenvalue ====
*
            ns = ns - 1
         ELSE
*
*           ==== One undeflatable eigenvalue.  Move it up out of the
*           .    way.   (ZTREXC can not fail in this case.) ====
*
            ifst = ns
            CALL ztrexc( 'V', jw, t, ldt, v, ldv, ifst, ilst, info )
            ilst = ilst + 1
         END IF
   10 CONTINUE
*
*        ==== Return to Hessenberg form ====
*
      IF( ns.EQ.0 )
     $   s = zero
*
      IF( ns.LT.jw ) THEN
*
*        ==== sorting the diagonal of T improves accuracy for
*        .    graded matrices.  ====
*
         DO 30 i = infqr + 1, ns
            ifst = i
            DO 20 j = i + 1, ns
               IF( cabs1( t( j, j ) ).GT.cabs1( t( ifst, ifst ) ) )
     $            ifst = j
   20       CONTINUE
            ilst = i
            IF( ifst.NE.ilst )
     $         CALL ztrexc( 'V', jw, t, ldt, v, ldv, ifst, ilst, info )
   30    CONTINUE
      END IF
*
*     ==== Restore shift/eigenvalue array from T ====
*
      DO 40 i = infqr + 1, jw
         sh( kwtop+i-1 ) = t( i, i )
   40 CONTINUE
*
*
      IF( ns.LT.jw .OR. s.EQ.zero ) THEN
         IF( ns.GT.1 .AND. s.NE.zero ) THEN
*
*           ==== Reflect spike back into lower triangle ====
*
            CALL zcopy( ns, v, ldv, work, 1 )
            DO 50 i = 1, ns
               work( i ) = dconjg( work( i ) )
   50       CONTINUE
            beta = work( 1 )
            CALL zlarfg( ns, beta, work( 2 ), 1, tau )
            work( 1 ) = one
*
            CALL zlaset( 'L', jw-2, jw-2, zero, zero, t( 3, 1 ), ldt )
*
            CALL zlarf( 'L', ns, jw, work, 1, dconjg( tau ), t, ldt,
     $                  work( jw+1 ) )
            CALL zlarf( 'R', ns, ns, work, 1, tau, t, ldt,
     $                  work( jw+1 ) )
            CALL zlarf( 'R', jw, ns, work, 1, tau, v, ldv,
     $                  work( jw+1 ) )
*
            CALL zgehrd( jw, 1, ns, t, ldt, work, work( jw+1 ),
     $                   lwork-jw, info )
         END IF
*
*        ==== Copy updated reduced window into place ====
*
         IF( kwtop.GT.1 )
     $      h( kwtop, kwtop-1 ) = s*dconjg( v( 1, 1 ) )
         CALL zlacpy( 'U', jw, jw, t, ldt, h( kwtop, kwtop ), ldh )
         CALL zcopy( jw-1, t( 2, 1 ), ldt+1, h( kwtop+1, kwtop ),
     $               ldh+1 )
*
*        ==== Accumulate orthogonal matrix in order update
*        .    H and Z, if requested.  ====
*
         IF( ns.GT.1 .AND. s.NE.zero )
     $      CALL zunmhr( 'R', 'N', jw, ns, 1, ns, t, ldt, work, v, ldv,
     $                   work( jw+1 ), lwork-jw, info )
*
*        ==== Update vertical slab in H ====
*
         IF( wantt ) THEN
            ltop = 1
         ELSE
            ltop = ktop
         END IF
         DO 60 krow = ltop, kwtop - 1, nv
            kln = min( nv, kwtop-krow )
            CALL zgemm( 'N', 'N', kln, jw, jw, one, h( krow, kwtop ),
     $                  ldh, v, ldv, zero, wv, ldwv )
            CALL zlacpy( 'A', kln, jw, wv, ldwv, h( krow, kwtop ), ldh )
   60    CONTINUE
*
*        ==== Update horizontal slab in H ====
*
         IF( wantt ) THEN
            DO 70 kcol = kbot + 1, n, nh
               kln = min( nh, n-kcol+1 )
               CALL zgemm( 'C', 'N', jw, kln, jw, one, v, ldv,
     $                     h( kwtop, kcol ), ldh, zero, t, ldt )
               CALL zlacpy( 'A', jw, kln, t, ldt, h( kwtop, kcol ),
     $                      ldh )
   70       CONTINUE
         END IF
*
*        ==== Update vertical slab in Z ====
*
         IF( wantz ) THEN
            DO 80 krow = iloz, ihiz, nv
               kln = min( nv, ihiz-krow+1 )
               CALL zgemm( 'N', 'N', kln, jw, jw, one, z( krow, kwtop ),
     $                     ldz, v, ldv, zero, wv, ldwv )
               CALL zlacpy( 'A', kln, jw, wv, ldwv, z( krow, kwtop ),
     $                      ldz )
   80       CONTINUE
         END IF
      END IF
*
*     ==== Return the number of deflations ... ====
*
      nd = jw - ns
*
*     ==== ... and the number of shifts. (Subtracting
*     .    INFQR from the spike length takes care
*     .    of the case of a rare QR failure while
*     .    calculating eigenvalues of the deflation
*     .    window.)  ====
*
      ns = ns - infqr
*
*      ==== Return optimal workspace. ====
*
      work( 1 ) = dcmplx( lwkopt, 0 )
*
*     ==== End of ZLAQR2 ====
*
      END
      
! ZLAQR1
      SUBROUTINE zlaqr1( N, H, LDH, S1, S2, V )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      COMPLEX*16         S1, S2
      INTEGER            LDH, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         H( LDH, * ), V( * )
*     ..
*
*  ================================================================
*
*     .. Parameters ..
      COMPLEX*16         ZERO
      parameter( zero = ( 0.0d0, 0.0d0 ) )
      DOUBLE PRECISION   RZERO
      parameter( rzero = 0.0d0 )
*     ..
*     .. Local Scalars ..
      COMPLEX*16         CDUM, H21S, H31S
      DOUBLE PRECISION   S
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, dble, dimag
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
*     ..
*     .. Statement Function definitions ..
      cabs1( cdum ) = abs( dble( cdum ) ) + abs( dimag( cdum ) )
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( n.NE.2 .AND. n.NE.3 ) THEN
         RETURN
      END IF
*
      IF( n.EQ.2 ) THEN
         s = cabs1( h( 1, 1 )-s2 ) + cabs1( h( 2, 1 ) )
         IF( s.EQ.rzero ) THEN
            v( 1 ) = zero
            v( 2 ) = zero
         ELSE
            h21s = h( 2, 1 ) / s
            v( 1 ) = h21s*h( 1, 2 ) + ( h( 1, 1 )-s1 )*
     $               ( ( h( 1, 1 )-s2 ) / s )
            v( 2 ) = h21s*( h( 1, 1 )+h( 2, 2 )-s1-s2 )
         END IF
      ELSE
         s = cabs1( h( 1, 1 )-s2 ) + cabs1( h( 2, 1 ) ) +
     $       cabs1( h( 3, 1 ) )
         IF( s.EQ.zero ) THEN
            v( 1 ) = zero
            v( 2 ) = zero
            v( 3 ) = zero
         ELSE
            h21s = h( 2, 1 ) / s
            h31s = h( 3, 1 ) / s
            v( 1 ) = ( h( 1, 1 )-s1 )*( ( h( 1, 1 )-s2 ) / s ) +
     $               h( 1, 2 )*h21s + h( 1, 3 )*h31s
            v( 2 ) = h21s*( h( 1, 1 )+h( 2, 2 )-s1-s2 ) + h( 2, 3 )*h31s
            v( 3 ) = h31s*( h( 1, 1 )+h( 3, 3 )-s1-s2 ) + h21s*h( 3, 2 )
         END IF
      END IF
      END
      
! ILAZLC
      INTEGER FUNCTION ilazlc( M, N, A, LDA )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            m, n, lda
*     ..
*     .. Array Arguments ..
      COMPLEX*16         a( lda, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16       zero
      parameter( zero = (0.0d+0, 0.0d+0) )
*     ..
*     .. Local Scalars ..
      INTEGER i
*     ..
*     .. Executable Statements ..
*
*     Quick test for the common case where one corner is non-zero.
      IF( n.EQ.0 ) THEN
         ilazlc = n
      ELSE IF( a(1, n).NE.zero .OR. a(m, n).NE.zero ) THEN
         ilazlc = n
      ELSE
*     Now scan each column from the end, returning with the first non-zero.
         DO ilazlc = n, 1, -1
            DO i = 1, m
               IF( a(i, ilazlc).NE.zero ) RETURN
            END DO
         END DO
      END IF
      RETURN
      END
      
! ZGERC
      SUBROUTINE zgerc(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
*
*  -- Reference BLAS level2 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      COMPLEX*16 ALPHA
      INTEGER INCX,INCY,LDA,M,N
*     ..
*     .. Array Arguments ..
      COMPLEX*16 A(LDA,*),X(*),Y(*)
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16 ZERO
      parameter(zero= (0.0d+0,0.0d+0))
*     ..
*     .. Local Scalars ..
      COMPLEX*16 TEMP
      INTEGER I,INFO,IX,J,JY,KX
*     ..
*     .. External Subroutines ..
      EXTERNAL xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC dconjg,max
*     ..
*
*     Test the input parameters.
*
      info = 0
      IF (m.LT.0) THEN
          info = 1
      ELSE IF (n.LT.0) THEN
          info = 2
      ELSE IF (incx.EQ.0) THEN
          info = 5
      ELSE IF (incy.EQ.0) THEN
          info = 7
      ELSE IF (lda.LT.max(1,m)) THEN
          info = 9
      END IF
      IF (info.NE.0) THEN
          CALL xerbla('ZGERC ',info)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((m.EQ.0) .OR. (n.EQ.0) .OR. (alpha.EQ.zero)) RETURN
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF (incy.GT.0) THEN
          jy = 1
      ELSE
          jy = 1 - (n-1)*incy
      END IF
      IF (incx.EQ.1) THEN
          DO 20 j = 1,n
              IF (y(jy).NE.zero) THEN
                  temp = alpha*dconjg(y(jy))
                  DO 10 i = 1,m
                      a(i,j) = a(i,j) + x(i)*temp
   10             CONTINUE
              END IF
              jy = jy + incy
   20     CONTINUE
      ELSE
          IF (incx.GT.0) THEN
              kx = 1
          ELSE
              kx = 1 - (m-1)*incx
          END IF
          DO 40 j = 1,n
              IF (y(jy).NE.zero) THEN
                  temp = alpha*dconjg(y(jy))
                  ix = kx
                  DO 30 i = 1,m
                      a(i,j) = a(i,j) + x(ix)*temp
                      ix = ix + incx
   30             CONTINUE
              END IF
              jy = jy + incy
   40     CONTINUE
      END IF
*
      RETURN
*
*     End of ZGERC
*
      END
      
! ZUNMQR
      SUBROUTINE zunmqr( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
     $                   WORK, LWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, LWORK, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NBMAX, LDT, TSIZE
      parameter( nbmax = 64, ldt = nbmax+1,
     $                     tsize = ldt*nbmax )
*     ..
*     .. Local Scalars ..
      LOGICAL            LEFT, LQUERY, NOTRAN
      INTEGER            I, I1, I2, I3, IB, IC, IINFO, IWT, JC, LDWORK,
     $                   lwkopt, mi, nb, nbmin, ni, nq, nw
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           lsame, ilaenv
*     ..
*     .. External Subroutines ..
      EXTERNAL           xerbla, zlarfb, zlarft, zunm2r
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          max, min
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      info = 0
      left = lsame( side, 'L' )
      notran = lsame( trans, 'N' )
      lquery = ( lwork.EQ.-1 )
*
*     NQ is the order of Q and NW is the minimum dimension of WORK
*
      IF( left ) THEN
         nq = m
         nw = max( 1, n )
      ELSE
         nq = n
         nw = max( 1, m )
      END IF
      IF( .NOT.left .AND. .NOT.lsame( side, 'R' ) ) THEN
         info = -1
      ELSE IF( .NOT.notran .AND. .NOT.lsame( trans, 'C' ) ) THEN
         info = -2
      ELSE IF( m.LT.0 ) THEN
         info = -3
      ELSE IF( n.LT.0 ) THEN
         info = -4
      ELSE IF( k.LT.0 .OR. k.GT.nq ) THEN
         info = -5
      ELSE IF( lda.LT.max( 1, nq ) ) THEN
         info = -7
      ELSE IF( ldc.LT.max( 1, m ) ) THEN
         info = -10
      ELSE IF( lwork.LT.nw .AND. .NOT.lquery ) THEN
         info = -12
      END IF
*
      IF( info.EQ.0 ) THEN
*
*        Compute the workspace requirements
*
         nb = min( nbmax, ilaenv( 1, 'ZUNMQR', side // trans, m, n, k,
     $        -1 ) )
         lwkopt = nw*nb + tsize
         work( 1 ) = lwkopt
      END IF
*
      IF( info.NE.0 ) THEN
         CALL xerbla( 'ZUNMQR', -info )
         RETURN
      ELSE IF( lquery ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( m.EQ.0 .OR. n.EQ.0 .OR. k.EQ.0 ) THEN
         work( 1 ) = 1
         RETURN
      END IF
*
      nbmin = 2
      ldwork = nw
      IF( nb.GT.1 .AND. nb.LT.k ) THEN
         IF( lwork.LT.lwkopt ) THEN
            nb = (lwork-tsize) / ldwork
            nbmin = max( 2, ilaenv( 2, 'ZUNMQR', side // trans, m, n, k,
     $              -1 ) )
         END IF
      END IF
*
      IF( nb.LT.nbmin .OR. nb.GE.k ) THEN
*
*        Use unblocked code
*
         CALL zunm2r( side, trans, m, n, k, a, lda, tau, c, ldc, work,
     $                iinfo )
      ELSE
*
*        Use blocked code
*
         iwt = 1 + nw*nb
         IF( ( left .AND. .NOT.notran ) .OR.
     $       ( .NOT.left .AND. notran ) ) THEN
            i1 = 1
            i2 = k
            i3 = nb
         ELSE
            i1 = ( ( k-1 ) / nb )*nb + 1
            i2 = 1
            i3 = -nb
         END IF
*
         IF( left ) THEN
            ni = n
            jc = 1
         ELSE
            mi = m
            ic = 1
         END IF
*
         DO 10 i = i1, i2, i3
            ib = min( nb, k-i+1 )
*
*           Form the triangular factor of the block reflector
*           H = H(i) H(i+1) . . . H(i+ib-1)
*
            CALL zlarft( 'Forward', 'Columnwise', nq-i+1, ib, a( i, i ),
     $                   lda, tau( i ), work( iwt ), ldt )
            IF( left ) THEN
*
*              H or H**H is applied to C(i:m,1:n)
*
               mi = m - i + 1
               ic = i
            ELSE
*
*              H or H**H is applied to C(1:m,i:n)
*
               ni = n - i + 1
               jc = i
            END IF
*
*           Apply H or H**H
*
            CALL zlarfb( side, trans, 'Forward', 'Columnwise', mi, ni,
     $                   ib, a( i, i ), lda, work( iwt ), ldt,
     $                   c( ic, jc ), ldc, work, ldwork )
   10    CONTINUE
      END IF
      work( 1 ) = lwkopt
      RETURN
*
*     End of ZUNMQR
*
      END
      
! ILAZLR
      INTEGER FUNCTION ilazlr( M, N, A, LDA )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            m, n, lda
*     ..
*     .. Array Arguments ..
      COMPLEX*16         a( lda, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16       zero
      parameter( zero = (0.0d+0, 0.0d+0) )
*     ..
*     .. Local Scalars ..
      INTEGER i, j
*     ..
*     .. Executable Statements ..
*
*     Quick test for the common case where one corner is non-zero.
      IF( m.EQ.0 ) THEN
         ilazlr = m
      ELSE IF( a(m, 1).NE.zero .OR. a(m, n).NE.zero ) THEN
         ilazlr = m
      ELSE
*     Scan up each column tracking the last zero row seen.
         ilazlr = 0
         DO j = 1, n
            i=m
            DO WHILE((a(max(i,1),j).EQ.zero).AND.(i.GE.1))
               i=i-1
            ENDDO
            ilazlr = max( ilazlr, i )
         END DO
      END IF
      RETURN
      END

! ZLARTG      
      SUBROUTINE ZLARTG( F, G, CS, SN, R )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   CS
      COMPLEX*16         F, G, R, SN
*     ..
*
*  Purpose
*  =======
*
*  ZLARTG generates a plane rotation so that
*
*     [  CS  SN  ]     [ F ]     [ R ]
*     [  __      ]  .  [   ]  =  [   ]   where CS**2 + |SN|**2 = 1.
*     [ -SN  CS  ]     [ G ]     [ 0 ]
*
*  This is a faster version of the BLAS1 routine ZROTG, except for
*  the following differences:
*     F and G are unchanged on return.
*     If G=0, then CS=1 and SN=0.
*     If F=0, then CS=0 and SN is chosen so that R is real.
*
*  Arguments
*  =========
*
*  F       (input) COMPLEX*16
*          The first component of vector to be rotated.
*
*  G       (input) COMPLEX*16
*          The second component of vector to be rotated.
*
*  CS      (output) DOUBLE PRECISION
*          The cosine of the rotation.
*
*  SN      (output) COMPLEX*16
*          The sine of the rotation.
*
*  R       (output) COMPLEX*16
*          The nonzero component of the rotated vector.
*
*  Further Details
*  ======= =======
*
*  3-5-96 - Modified with a new algorithm by W. Kahan and J. Demmel
*
*  This version has a few statements commented out for thread safety
*  (machine parameters are computed on each entry). 10 feb 03, SJH.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   TWO, ONE, ZERO
      PARAMETER          ( TWO = 2.0D+0, ONE = 1.0D+0, ZERO = 0.0D+0 )
      COMPLEX*16         CZERO
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
*     LOGICAL            FIRST
      INTEGER            COUNT, I
      DOUBLE PRECISION   D, DI, DR, EPS, F2, F2S, G2, G2S, SAFMIN,
     $                   SAFMN2, SAFMX2, SCALE
      COMPLEX*16         FF, FS, GS
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLAPY2
      EXTERNAL           DLAMCH, DLAPY2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCMPLX, DCONJG, DIMAG, INT, LOG,
     $                   MAX, SQRT
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   ABS1, ABSSQ
*     ..
*     .. Save statement ..
*     SAVE               FIRST, SAFMX2, SAFMIN, SAFMN2
*     ..
*     .. Data statements ..
*     DATA               FIRST / .TRUE. /
*     ..
*     .. Statement Function definitions ..
      ABS1( FF ) = MAX( ABS( DBLE( FF ) ), ABS( DIMAG( FF ) ) )
      ABSSQ( FF ) = DBLE( FF )**2 + DIMAG( FF )**2
*     ..
*     .. Executable Statements ..
*
*     IF( FIRST ) THEN
         SAFMIN = DLAMCH( 'S' )
         EPS = DLAMCH( 'E' )
         SAFMN2 = DLAMCH( 'B' )**INT( LOG( SAFMIN / EPS ) /
     $            LOG( DLAMCH( 'B' ) ) / TWO )
         SAFMX2 = ONE / SAFMN2
*        FIRST = .FALSE.
*     END IF
      SCALE = MAX( ABS1( F ), ABS1( G ) )
      FS = F
      GS = G
      COUNT = 0
      IF( SCALE.GE.SAFMX2 ) THEN
   10    CONTINUE
         COUNT = COUNT + 1
         FS = FS*SAFMN2
         GS = GS*SAFMN2
         SCALE = SCALE*SAFMN2
         IF( SCALE.GE.SAFMX2 )
     $      GO TO 10
      ELSE IF( SCALE.LE.SAFMN2 ) THEN
         IF( G.EQ.CZERO ) THEN
            CS = ONE
            SN = CZERO
            R = F
            RETURN
         END IF
   20    CONTINUE
         COUNT = COUNT - 1
         FS = FS*SAFMX2
         GS = GS*SAFMX2
         SCALE = SCALE*SAFMX2
         IF( SCALE.LE.SAFMN2 )
     $      GO TO 20
      END IF
      F2 = ABSSQ( FS )
      G2 = ABSSQ( GS )
      IF( F2.LE.MAX( G2, ONE )*SAFMIN ) THEN
*
*        This is a rare case: F is very small.
*
         IF( F.EQ.CZERO ) THEN
            CS = ZERO
            R = DLAPY2( DBLE( G ), DIMAG( G ) )
*           Do complex/real division explicitly with two real divisions
            D = DLAPY2( DBLE( GS ), DIMAG( GS ) )
            SN = DCMPLX( DBLE( GS ) / D, -DIMAG( GS ) / D )
            RETURN
         END IF
         F2S = DLAPY2( DBLE( FS ), DIMAG( FS ) )
*        G2 and G2S are accurate
*        G2 is at least SAFMIN, and G2S is at least SAFMN2
         G2S = SQRT( G2 )
*        Error in CS from underflow in F2S is at most
*        UNFL / SAFMN2 .lt. sqrt(UNFL*EPS) .lt. EPS
*        If MAX(G2,ONE)=G2, then F2 .lt. G2*SAFMIN,
*        and so CS .lt. sqrt(SAFMIN)
*        If MAX(G2,ONE)=ONE, then F2 .lt. SAFMIN
*        and so CS .lt. sqrt(SAFMIN)/SAFMN2 = sqrt(EPS)
*        Therefore, CS = F2S/G2S / sqrt( 1 + (F2S/G2S)**2 ) = F2S/G2S
         CS = F2S / G2S
*        Make sure abs(FF) = 1
*        Do complex/real division explicitly with 2 real divisions
         IF( ABS1( F ).GT.ONE ) THEN
            D = DLAPY2( DBLE( F ), DIMAG( F ) )
            FF = DCMPLX( DBLE( F ) / D, DIMAG( F ) / D )
         ELSE
            DR = SAFMX2*DBLE( F )
            DI = SAFMX2*DIMAG( F )
            D = DLAPY2( DR, DI )
            FF = DCMPLX( DR / D, DI / D )
         END IF
         SN = FF*DCMPLX( DBLE( GS ) / G2S, -DIMAG( GS ) / G2S )
         R = CS*F + SN*G
      ELSE
*
*        This is the most common case.
*        Neither F2 nor F2/G2 are less than SAFMIN
*        F2S cannot overflow, and it is accurate
*
         F2S = SQRT( ONE+G2 / F2 )
*        Do the F2S(real)*FS(complex) multiply with two real multiplies
         R = DCMPLX( F2S*DBLE( FS ), F2S*DIMAG( FS ) )
         CS = ONE / F2S
         D = F2 + G2
*        Do complex/real division explicitly with two real divisions
         SN = DCMPLX( DBLE( R ) / D, DIMAG( R ) / D )
         SN = SN*DCONJG( GS )
         IF( COUNT.NE.0 ) THEN
            IF( COUNT.GT.0 ) THEN
               DO 30 I = 1, COUNT
                  R = R*SAFMX2
   30          CONTINUE
            ELSE
               DO 40 I = 1, -COUNT
                  R = R*SAFMN2
   40          CONTINUE
            END IF
         END IF
      END IF
      RETURN
*
*     End of ZLARTG
*
      END
      
! ZROT
      SUBROUTINE zrot( N, CX, INCX, CY, INCY, C, S )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INCX, INCY, N
      DOUBLE PRECISION   C
      COMPLEX*16         S
*     ..
*     .. Array Arguments ..
      COMPLEX*16         CX( * ), CY( * )
*     ..
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, IX, IY
      COMPLEX*16         STEMP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          dconjg
*     ..
*     .. Executable Statements ..
*
      IF( n.LE.0 )
     $   RETURN
      IF( incx.EQ.1 .AND. incy.EQ.1 )
     $   GO TO 20
*
*     Code for unequal increments or equal increments not equal to 1
*
      ix = 1
      iy = 1
      IF( incx.LT.0 )
     $   ix = ( -n+1 )*incx + 1
      IF( incy.LT.0 )
     $   iy = ( -n+1 )*incy + 1
      DO 10 i = 1, n
         stemp = c*cx( ix ) + s*cy( iy )
         cy( iy ) = c*cy( iy ) - dconjg( s )*cx( ix )
         cx( ix ) = stemp
         ix = ix + incx
         iy = iy + incy
   10 CONTINUE
      RETURN
*
*     Code for both increments equal to 1
*
   20 CONTINUE
      DO 30 i = 1, n
         stemp = c*cx( i ) + s*cy( i )
         cy( i ) = c*cy( i ) - dconjg( s )*cx( i )
         cx( i ) = stemp
   30 CONTINUE
      RETURN
      END
      
! ZUNM2R
      SUBROUTINE zunm2r( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
     $                   WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE
      parameter( one = ( 1.0d+0, 0.0d+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            LEFT, NOTRAN
      INTEGER            I, I1, I2, I3, IC, JC, MI, NI, NQ
      COMPLEX*16         AII, TAUI
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           lsame
*     ..
*     .. External Subroutines ..
      EXTERNAL           xerbla, zlarf
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          dconjg, max
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      info = 0
      left = lsame( side, 'L' )
      notran = lsame( trans, 'N' )
*
*     NQ is the order of Q
*
      IF( left ) THEN
         nq = m
      ELSE
         nq = n
      END IF
      IF( .NOT.left .AND. .NOT.lsame( side, 'R' ) ) THEN
         info = -1
      ELSE IF( .NOT.notran .AND. .NOT.lsame( trans, 'C' ) ) THEN
         info = -2
      ELSE IF( m.LT.0 ) THEN
         info = -3
      ELSE IF( n.LT.0 ) THEN
         info = -4
      ELSE IF( k.LT.0 .OR. k.GT.nq ) THEN
         info = -5
      ELSE IF( lda.LT.max( 1, nq ) ) THEN
         info = -7
      ELSE IF( ldc.LT.max( 1, m ) ) THEN
         info = -10
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'ZUNM2R', -info )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( m.EQ.0 .OR. n.EQ.0 .OR. k.EQ.0 )
     $   RETURN
*
      IF( ( left .AND. .NOT.notran .OR. .NOT.left .AND. notran ) ) THEN
         i1 = 1
         i2 = k
         i3 = 1
      ELSE
         i1 = k
         i2 = 1
         i3 = -1
      END IF
*
      IF( left ) THEN
         ni = n
         jc = 1
      ELSE
         mi = m
         ic = 1
      END IF
*
      DO 10 i = i1, i2, i3
         IF( left ) THEN
*
*           H(i) or H(i)**H is applied to C(i:m,1:n)
*
            mi = m - i + 1
            ic = i
         ELSE
*
*           H(i) or H(i)**H is applied to C(1:m,i:n)
*
            ni = n - i + 1
            jc = i
         END IF
*
*        Apply H(i) or H(i)**H
*
         IF( notran ) THEN
            taui = tau( i )
         ELSE
            taui = dconjg( tau( i ) )
         END IF
         aii = a( i, i )
         a( i, i ) = one
         CALL zlarf( side, mi, ni, a( i, i ), 1, taui, c( ic, jc ), ldc,
     $               work )
         a( i, i ) = aii
   10 CONTINUE
      RETURN
*
*     End of ZUNM2R
*
      END
      
      