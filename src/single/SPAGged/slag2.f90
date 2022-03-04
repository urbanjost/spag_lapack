!*==slag2.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SLAG2 computes the eigenvalues of a 2-by-2 generalized eigenvalue problem, with scaling as necessary to avoid over-/underflow.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLAG2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slag2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slag2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slag2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLAG2( A, LDA, B, LDB, SAFMIN, SCALE1, SCALE2, WR1,
!                         WR2, WI )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LDB
!       REAL               SAFMIN, SCALE1, SCALE2, WI, WR1, WR2
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), B( LDB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLAG2 computes the eigenvalues of a 2 x 2 generalized eigenvalue
!> problem  A - w B, with scaling as necessary to avoid over-/underflow.
!>
!> The scaling factor "s" results in a modified eigenvalue equation
!>
!>     s A - w B
!>
!> where  s  is a non-negative scaling factor chosen so that  w,  w B,
!> and  s A  do not overflow and, if possible, do not underflow, either.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] A
!> \verbatim
!>          A is REAL array, dimension (LDA, 2)
!>          On entry, the 2 x 2 matrix A.  It is assumed that its 1-norm
!>          is less than 1/SAFMIN.  Entries less than
!>          sqrt(SAFMIN)*norm(A) are subject to being treated as zero.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= 2.
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is REAL array, dimension (LDB, 2)
!>          On entry, the 2 x 2 upper triangular matrix B.  It is
!>          assumed that the one-norm of B is less than 1/SAFMIN.  The
!>          diagonals should be at least sqrt(SAFMIN) times the largest
!>          element of B (in absolute value); if a diagonal is smaller
!>          than that, then  +/- sqrt(SAFMIN) will be used instead of
!>          that diagonal.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= 2.
!> \endverbatim
!>
!> \param[in] SAFMIN
!> \verbatim
!>          SAFMIN is REAL
!>          The smallest positive number s.t. 1/SAFMIN does not
!>          overflow.  (This should always be SLAMCH('S') -- it is an
!>          argument in order to avoid having to call SLAMCH frequently.)
!> \endverbatim
!>
!> \param[out] SCALE1
!> \verbatim
!>          SCALE1 is REAL
!>          A scaling factor used to avoid over-/underflow in the
!>          eigenvalue equation which defines the first eigenvalue.  If
!>          the eigenvalues are complex, then the eigenvalues are
!>          ( WR1  +/-  WI i ) / SCALE1  (which may lie outside the
!>          exponent range of the machine), SCALE1=SCALE2, and SCALE1
!>          will always be positive.  If the eigenvalues are real, then
!>          the first (real) eigenvalue is  WR1 / SCALE1 , but this may
!>          overflow or underflow, and in fact, SCALE1 may be zero or
!>          less than the underflow threshold if the exact eigenvalue
!>          is sufficiently large.
!> \endverbatim
!>
!> \param[out] SCALE2
!> \verbatim
!>          SCALE2 is REAL
!>          A scaling factor used to avoid over-/underflow in the
!>          eigenvalue equation which defines the second eigenvalue.  If
!>          the eigenvalues are complex, then SCALE2=SCALE1.  If the
!>          eigenvalues are real, then the second (real) eigenvalue is
!>          WR2 / SCALE2 , but this may overflow or underflow, and in
!>          fact, SCALE2 may be zero or less than the underflow
!>          threshold if the exact eigenvalue is sufficiently large.
!> \endverbatim
!>
!> \param[out] WR1
!> \verbatim
!>          WR1 is REAL
!>          If the eigenvalue is real, then WR1 is SCALE1 times the
!>          eigenvalue closest to the (2,2) element of A B**(-1).  If the
!>          eigenvalue is complex, then WR1=WR2 is SCALE1 times the real
!>          part of the eigenvalues.
!> \endverbatim
!>
!> \param[out] WR2
!> \verbatim
!>          WR2 is REAL
!>          If the eigenvalue is real, then WR2 is SCALE2 times the
!>          other eigenvalue.  If the eigenvalue is complex, then
!>          WR1=WR2 is SCALE1 times the real part of the eigenvalues.
!> \endverbatim
!>
!> \param[out] WI
!> \verbatim
!>          WI is REAL
!>          If the eigenvalue is real, then WI is zero.  If the
!>          eigenvalue is complex, then WI is SCALE1 times the imaginary
!>          part of the eigenvalues.  WI will always be non-negative.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \date June 2016
!
!> \ingroup realOTHERauxiliary
!
!  =====================================================================
      SUBROUTINE SLAG2(A,Lda,B,Ldb,Safmin,Scale1,Scale2,Wr1,Wr2,Wi)
      IMPLICIT NONE
!*--SLAG2159
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 ,              &
     &                      TWO = 2.0E+0 , HALF = ONE/TWO ,             &
     &                      FUZZY1 = ONE + 1.0E-5
!
! Dummy argument declarations rewritten by SPAG
!
      REAL , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(IN) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      REAL , INTENT(IN) :: Safmin
      REAL , INTENT(INOUT) :: Scale1
      REAL , INTENT(OUT) :: Scale2
      REAL , INTENT(INOUT) :: Wr1
      REAL , INTENT(INOUT) :: Wr2
      REAL , INTENT(INOUT) :: Wi
!
! Local variable declarations rewritten by SPAG
!
      REAL :: a11 , a12 , a21 , a22 , abi22 , anorm , as11 , as12 ,     &
     &        as22 , ascale , b11 , b12 , b22 , binv11 , binv22 , bmin ,&
     &        bnorm , bscale , bsize , c1 , c2 , c3 , c4 , c5 , diff ,  &
     &        discr , pp , qq , r , rtmax , rtmin , s1 , s2 , safmax ,  &
     &        shift , ss , sum , wabs , wbig , wdet , wscale , wsize ,  &
     &        wsmall
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
      rtmin = SQRT(Safmin)
      rtmax = ONE/rtmin
      safmax = ONE/Safmin
!
!     Scale A
!
      anorm = MAX(ABS(A(1,1))+ABS(A(2,1)),ABS(A(1,2))+ABS(A(2,2)),      &
     &        Safmin)
      ascale = ONE/anorm
      a11 = ascale*A(1,1)
      a21 = ascale*A(2,1)
      a12 = ascale*A(1,2)
      a22 = ascale*A(2,2)
!
!     Perturb B if necessary to insure non-singularity
!
      b11 = B(1,1)
      b12 = B(1,2)
      b22 = B(2,2)
      bmin = rtmin*MAX(ABS(b11),ABS(b12),ABS(b22),rtmin)
      IF ( ABS(b11)<bmin ) b11 = SIGN(bmin,b11)
      IF ( ABS(b22)<bmin ) b22 = SIGN(bmin,b22)
!
!     Scale B
!
      bnorm = MAX(ABS(b11),ABS(b12)+ABS(b22),Safmin)
      bsize = MAX(ABS(b11),ABS(b22))
      bscale = ONE/bsize
      b11 = b11*bscale
      b12 = b12*bscale
      b22 = b22*bscale
!
!     Compute larger eigenvalue by method described by C. van Loan
!
!     ( AS is A shifted by -SHIFT*B )
!
      binv11 = ONE/b11
      binv22 = ONE/b22
      s1 = a11*binv11
      s2 = a22*binv22
      IF ( ABS(s1)<=ABS(s2) ) THEN
         as12 = a12 - s1*b12
         as22 = a22 - s1*b22
         ss = a21*(binv11*binv22)
         abi22 = as22*binv22 - ss*b12
         pp = HALF*abi22
         shift = s1
      ELSE
         as12 = a12 - s2*b12
         as11 = a11 - s2*b11
         ss = a21*(binv11*binv22)
         abi22 = -ss*b12
         pp = HALF*(as11*binv11+abi22)
         shift = s2
      ENDIF
      qq = ss*as12
      IF ( ABS(pp*rtmin)>=ONE ) THEN
         discr = (rtmin*pp)**2 + qq*Safmin
         r = SQRT(ABS(discr))*rtmax
      ELSEIF ( pp**2+ABS(qq)<=Safmin ) THEN
         discr = (rtmax*pp)**2 + qq*safmax
         r = SQRT(ABS(discr))*rtmin
      ELSE
         discr = pp**2 + qq
         r = SQRT(ABS(discr))
      ENDIF
!
!     Note: the test of R in the following IF is to cover the case when
!           DISCR is small and negative and is flushed to zero during
!           the calculation of R.  On machines which have a consistent
!           flush-to-zero threshold and handle numbers above that
!           threshold correctly, it would not be necessary.
!
      IF ( discr>=ZERO .OR. r==ZERO ) THEN
         sum = pp + SIGN(r,pp)
         diff = pp - SIGN(r,pp)
         wbig = shift + sum
!
!        Compute smaller eigenvalue
!
         wsmall = shift + diff
         IF ( HALF*ABS(wbig)>MAX(ABS(wsmall),Safmin) ) THEN
            wdet = (a11*a22-a12*a21)*(binv11*binv22)
            wsmall = wdet/wbig
         ENDIF
!
!        Choose (real) eigenvalue closest to 2,2 element of A*B**(-1)
!        for WR1.
!
         IF ( pp>abi22 ) THEN
            Wr1 = MIN(wbig,wsmall)
            Wr2 = MAX(wbig,wsmall)
         ELSE
            Wr1 = MAX(wbig,wsmall)
            Wr2 = MIN(wbig,wsmall)
         ENDIF
         Wi = ZERO
      ELSE
!
!        Complex eigenvalues
!
         Wr1 = shift + pp
         Wr2 = Wr1
         Wi = r
      ENDIF
!
!     Further scaling to avoid underflow and overflow in computing
!     SCALE1 and overflow in computing w*B.
!
!     This scale factor (WSCALE) is bounded from above using C1 and C2,
!     and from below using C3 and C4.
!        C1 implements the condition  s A  must never overflow.
!        C2 implements the condition  w B  must never overflow.
!        C3, with C2,
!           implement the condition that s A - w B must never overflow.
!        C4 implements the condition  s    should not underflow.
!        C5 implements the condition  max(s,|w|) should be at least 2.
!
      c1 = bsize*(Safmin*MAX(ONE,ascale))
      c2 = Safmin*MAX(ONE,bnorm)
      c3 = bsize*Safmin
      IF ( ascale<=ONE .AND. bsize<=ONE ) THEN
         c4 = MIN(ONE,(ascale/Safmin)*bsize)
      ELSE
         c4 = ONE
      ENDIF
      IF ( ascale<=ONE .OR. bsize<=ONE ) THEN
         c5 = MIN(ONE,ascale*bsize)
      ELSE
         c5 = ONE
      ENDIF
!
!     Scale first eigenvalue
!
      wabs = ABS(Wr1) + ABS(Wi)
      wsize = MAX(Safmin,c1,FUZZY1*(wabs*c2+c3),                        &
     &        MIN(c4,HALF*MAX(wabs,c5)))
      IF ( wsize/=ONE ) THEN
         wscale = ONE/wsize
         IF ( wsize>ONE ) THEN
            Scale1 = (MAX(ascale,bsize)*wscale)*MIN(ascale,bsize)
         ELSE
            Scale1 = (MIN(ascale,bsize)*wscale)*MAX(ascale,bsize)
         ENDIF
         Wr1 = Wr1*wscale
         IF ( Wi/=ZERO ) THEN
            Wi = Wi*wscale
            Wr2 = Wr1
            Scale2 = Scale1
         ENDIF
      ELSE
         Scale1 = ascale*bsize
         Scale2 = Scale1
      ENDIF
!
!     Scale second eigenvalue (if real)
!
      IF ( Wi==ZERO ) THEN
         wsize = MAX(Safmin,c1,FUZZY1*(ABS(Wr2)*c2+c3),                 &
     &           MIN(c4,HALF*MAX(ABS(Wr2),c5)))
         IF ( wsize/=ONE ) THEN
            wscale = ONE/wsize
            IF ( wsize>ONE ) THEN
               Scale2 = (MAX(ascale,bsize)*wscale)*MIN(ascale,bsize)
            ELSE
               Scale2 = (MIN(ascale,bsize)*wscale)*MAX(ascale,bsize)
            ENDIF
            Wr2 = Wr2*wscale
         ELSE
            Scale2 = ascale*bsize
         ENDIF
      ENDIF
!
!     End of SLAG2
!
      END SUBROUTINE SLAG2
