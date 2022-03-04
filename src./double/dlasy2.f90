!*==dlasy2.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DLASY2 solves the Sylvester matrix equation where the matrices are of order 1 or 2.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLASY2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasy2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlasy2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlasy2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLASY2( LTRANL, LTRANR, ISGN, N1, N2, TL, LDTL, TR,
!                          LDTR, B, LDB, SCALE, X, LDX, XNORM, INFO )
!
!       .. Scalar Arguments ..
!       LOGICAL            LTRANL, LTRANR
!       INTEGER            INFO, ISGN, LDB, LDTL, LDTR, LDX, N1, N2
!       DOUBLE PRECISION   SCALE, XNORM
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   B( LDB, * ), TL( LDTL, * ), TR( LDTR, * ),
!      $                   X( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLASY2 solves for the N1 by N2 matrix X, 1 <= N1,N2 <= 2, in
!>
!>        op(TL)*X + ISGN*X*op(TR) = SCALE*B,
!>
!> where TL is N1 by N1, TR is N2 by N2, B is N1 by N2, and ISGN = 1 or
!> -1.  op(T) = T or T**T, where T**T denotes the transpose of T.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] LTRANL
!> \verbatim
!>          LTRANL is LOGICAL
!>          On entry, LTRANL specifies the op(TL):
!>             = .FALSE., op(TL) = TL,
!>             = .TRUE., op(TL) = TL**T.
!> \endverbatim
!>
!> \param[in] LTRANR
!> \verbatim
!>          LTRANR is LOGICAL
!>          On entry, LTRANR specifies the op(TR):
!>            = .FALSE., op(TR) = TR,
!>            = .TRUE., op(TR) = TR**T.
!> \endverbatim
!>
!> \param[in] ISGN
!> \verbatim
!>          ISGN is INTEGER
!>          On entry, ISGN specifies the sign of the equation
!>          as described before. ISGN may only be 1 or -1.
!> \endverbatim
!>
!> \param[in] N1
!> \verbatim
!>          N1 is INTEGER
!>          On entry, N1 specifies the order of matrix TL.
!>          N1 may only be 0, 1 or 2.
!> \endverbatim
!>
!> \param[in] N2
!> \verbatim
!>          N2 is INTEGER
!>          On entry, N2 specifies the order of matrix TR.
!>          N2 may only be 0, 1 or 2.
!> \endverbatim
!>
!> \param[in] TL
!> \verbatim
!>          TL is DOUBLE PRECISION array, dimension (LDTL,2)
!>          On entry, TL contains an N1 by N1 matrix.
!> \endverbatim
!>
!> \param[in] LDTL
!> \verbatim
!>          LDTL is INTEGER
!>          The leading dimension of the matrix TL. LDTL >= max(1,N1).
!> \endverbatim
!>
!> \param[in] TR
!> \verbatim
!>          TR is DOUBLE PRECISION array, dimension (LDTR,2)
!>          On entry, TR contains an N2 by N2 matrix.
!> \endverbatim
!>
!> \param[in] LDTR
!> \verbatim
!>          LDTR is INTEGER
!>          The leading dimension of the matrix TR. LDTR >= max(1,N2).
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (LDB,2)
!>          On entry, the N1 by N2 matrix B contains the right-hand
!>          side of the equation.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the matrix B. LDB >= max(1,N1).
!> \endverbatim
!>
!> \param[out] SCALE
!> \verbatim
!>          SCALE is DOUBLE PRECISION
!>          On exit, SCALE contains the scale factor. SCALE is chosen
!>          less than or equal to 1 to prevent the solution overflowing.
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension (LDX,2)
!>          On exit, X contains the N1 by N2 solution.
!> \endverbatim
!>
!> \param[in] LDX
!> \verbatim
!>          LDX is INTEGER
!>          The leading dimension of the matrix X. LDX >= max(1,N1).
!> \endverbatim
!>
!> \param[out] XNORM
!> \verbatim
!>          XNORM is DOUBLE PRECISION
!>          On exit, XNORM is the infinity-norm of the solution.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          On exit, INFO is set to
!>             0: successful exit.
!>             1: TL and TR have too close eigenvalues, so TL or
!>                TR is perturbed to get a nonsingular equation.
!>          NOTE: In the interests of speed, this routine does not
!>                check the inputs for errors.
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
!> \ingroup doubleSYauxiliary
!
!  =====================================================================
      SUBROUTINE DLASY2(Ltranl,Ltranr,Isgn,N1,N2,Tl,Ldtl,Tr,Ldtr,B,Ldb, &
     &                  Scale,X,Ldx,Xnorm,Info)
      USE F77KINDS                        
      USE S_DCOPY
      USE S_DLAMCH
      USE S_DSWAP
      USE S_IDAMAX
      IMPLICIT NONE
!*--DLASY2183
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              TWO = 2.0D+0 , HALF = 0.5D+0 ,      &
     &                              EIGHT = 8.0D+0
!
! Dummy argument declarations rewritten by SPAG
!
      LOGICAL , INTENT(IN) :: Ltranl
      LOGICAL , INTENT(IN) :: Ltranr
      INTEGER , INTENT(IN) :: Isgn
      INTEGER , INTENT(IN) :: N1
      INTEGER , INTENT(IN) :: N2
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Ldtl,*) :: Tl
      INTEGER , INTENT(IN) :: Ldtl
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Ldtr,*) :: Tr
      INTEGER , INTENT(IN) :: Ldtr
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      REAL(R8KIND) , INTENT(INOUT) :: Scale
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL(R8KIND) , INTENT(OUT) :: Xnorm
      INTEGER , INTENT(OUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: bet , eps , gam , l21 , sgn , smin , smlnum ,     &
     &                tau1 , temp , u11 , u12 , u22 , xmax
      LOGICAL :: bswap , xswap
      LOGICAL , DIMENSION(4) , SAVE :: bswpiv , xswpiv
      REAL(R8KIND) , DIMENSION(4) :: btmp , tmp
      INTEGER :: i , ip , ipiv , ipsv , j , jp , jpsv , k
      INTEGER , DIMENSION(4) :: jpiv
      INTEGER , DIMENSION(4) , SAVE :: locl21 , locu12 , locu22
      REAL(R8KIND) , DIMENSION(4,4) :: t16
      REAL(R8KIND) , DIMENSION(2) :: x2
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
! =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. Local Arrays ..
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Data statements ..
      DATA locu12/3 , 4 , 1 , 2/ , locl21/2 , 1 , 4 , 3/ , locu22/4 ,   &
     &     3 , 2 , 1/
      DATA xswpiv/.FALSE. , .FALSE. , .TRUE. , .TRUE./
      DATA bswpiv/.FALSE. , .TRUE. , .FALSE. , .TRUE./
!     ..
!     .. Executable Statements ..
!
!     Do not check the input parameters for errors
!
      Info = 0
!
!     Quick return if possible
!
      IF ( N1==0 .OR. N2==0 ) RETURN
!
!     Set constants to control overflow
!
      eps = DLAMCH('P')
      smlnum = DLAMCH('S')/eps
      sgn = Isgn
!
      k = N1 + N1 + N2 - 2
      IF ( k==2 ) THEN
!
!     1 by 2:
!     TL11*[X11 X12] + ISGN*[X11 X12]*op[TR11 TR12]  = [B11 B12]
!                                       [TR21 TR22]
!
!
         smin = MAX(eps*MAX(ABS(Tl(1,1)),ABS(Tr(1,1)),ABS(Tr(1,2)),     &
     &          ABS(Tr(2,1)),ABS(Tr(2,2))),smlnum)
         tmp(1) = Tl(1,1) + sgn*Tr(1,1)
         tmp(4) = Tl(1,1) + sgn*Tr(2,2)
         IF ( Ltranr ) THEN
            tmp(2) = sgn*Tr(2,1)
            tmp(3) = sgn*Tr(1,2)
         ELSE
            tmp(2) = sgn*Tr(1,2)
            tmp(3) = sgn*Tr(2,1)
         ENDIF
         btmp(1) = B(1,1)
         btmp(2) = B(1,2)
      ELSEIF ( k==3 ) THEN
!
!     2 by 1:
!          op[TL11 TL12]*[X11] + ISGN* [X11]*TR11  = [B11]
!            [TL21 TL22] [X21]         [X21]         [B21]
!
         smin = MAX(eps*MAX(ABS(Tr(1,1)),ABS(Tl(1,1)),ABS(Tl(1,2)),     &
     &          ABS(Tl(2,1)),ABS(Tl(2,2))),smlnum)
         tmp(1) = Tl(1,1) + sgn*Tr(1,1)
         tmp(4) = Tl(2,2) + sgn*Tr(1,1)
         IF ( Ltranl ) THEN
            tmp(2) = Tl(1,2)
            tmp(3) = Tl(2,1)
         ELSE
            tmp(2) = Tl(2,1)
            tmp(3) = Tl(1,2)
         ENDIF
         btmp(1) = B(1,1)
         btmp(2) = B(2,1)
      ELSEIF ( k==4 ) THEN
!
!     2 by 2:
!     op[TL11 TL12]*[X11 X12] +ISGN* [X11 X12]*op[TR11 TR12] = [B11 B12]
!       [TL21 TL22] [X21 X22]        [X21 X22]   [TR21 TR22]   [B21 B22]
!
!     Solve equivalent 4 by 4 system using complete pivoting.
!     Set pivots less than SMIN to SMIN.
!
         smin = MAX(ABS(Tr(1,1)),ABS(Tr(1,2)),ABS(Tr(2,1)),ABS(Tr(2,2)))
         smin = MAX(smin,ABS(Tl(1,1)),ABS(Tl(1,2)),ABS(Tl(2,1)),        &
     &          ABS(Tl(2,2)))
         smin = MAX(eps*smin,smlnum)
         btmp(1) = ZERO
         CALL DCOPY(16,btmp,0,t16,1)
         t16(1,1) = Tl(1,1) + sgn*Tr(1,1)
         t16(2,2) = Tl(2,2) + sgn*Tr(1,1)
         t16(3,3) = Tl(1,1) + sgn*Tr(2,2)
         t16(4,4) = Tl(2,2) + sgn*Tr(2,2)
         IF ( Ltranl ) THEN
            t16(1,2) = Tl(2,1)
            t16(2,1) = Tl(1,2)
            t16(3,4) = Tl(2,1)
            t16(4,3) = Tl(1,2)
         ELSE
            t16(1,2) = Tl(1,2)
            t16(2,1) = Tl(2,1)
            t16(3,4) = Tl(1,2)
            t16(4,3) = Tl(2,1)
         ENDIF
         IF ( Ltranr ) THEN
            t16(1,3) = sgn*Tr(1,2)
            t16(2,4) = sgn*Tr(1,2)
            t16(3,1) = sgn*Tr(2,1)
            t16(4,2) = sgn*Tr(2,1)
         ELSE
            t16(1,3) = sgn*Tr(2,1)
            t16(2,4) = sgn*Tr(2,1)
            t16(3,1) = sgn*Tr(1,2)
            t16(4,2) = sgn*Tr(1,2)
         ENDIF
         btmp(1) = B(1,1)
         btmp(2) = B(2,1)
         btmp(3) = B(1,2)
         btmp(4) = B(2,2)
!
!     Perform elimination
!
         DO i = 1 , 3
            xmax = ZERO
            DO ip = i , 4
               DO jp = i , 4
                  IF ( ABS(t16(ip,jp))>=xmax ) THEN
                     xmax = ABS(t16(ip,jp))
                     ipsv = ip
                     jpsv = jp
                  ENDIF
               ENDDO
            ENDDO
            IF ( ipsv/=i ) THEN
               CALL DSWAP(4,t16(ipsv,1),4,t16(i,1),4)
               temp = btmp(i)
               btmp(i) = btmp(ipsv)
               btmp(ipsv) = temp
            ENDIF
            IF ( jpsv/=i ) CALL DSWAP(4,t16(1,jpsv),1,t16(1,i),1)
            jpiv(i) = jpsv
            IF ( ABS(t16(i,i))<smin ) THEN
               Info = 1
               t16(i,i) = smin
            ENDIF
            DO j = i + 1 , 4
               t16(j,i) = t16(j,i)/t16(i,i)
               btmp(j) = btmp(j) - t16(j,i)*btmp(i)
               DO k = i + 1 , 4
                  t16(j,k) = t16(j,k) - t16(j,i)*t16(i,k)
               ENDDO
            ENDDO
         ENDDO
         IF ( ABS(t16(4,4))<smin ) THEN
            Info = 1
            t16(4,4) = smin
         ENDIF
         Scale = ONE
         IF ( (EIGHT*smlnum)*ABS(btmp(1))>ABS(t16(1,1)) .OR.            &
     &        (EIGHT*smlnum)*ABS(btmp(2))>ABS(t16(2,2)) .OR.            &
     &        (EIGHT*smlnum)*ABS(btmp(3))>ABS(t16(3,3)) .OR.            &
     &        (EIGHT*smlnum)*ABS(btmp(4))>ABS(t16(4,4)) ) THEN
            Scale = (ONE/EIGHT)                                         &
     &              /MAX(ABS(btmp(1)),ABS(btmp(2)),ABS(btmp(3)),        &
     &              ABS(btmp(4)))
            btmp(1) = btmp(1)*Scale
            btmp(2) = btmp(2)*Scale
            btmp(3) = btmp(3)*Scale
            btmp(4) = btmp(4)*Scale
         ENDIF
         DO i = 1 , 4
            k = 5 - i
            temp = ONE/t16(k,k)
            tmp(k) = btmp(k)*temp
            DO j = k + 1 , 4
               tmp(k) = tmp(k) - (temp*t16(k,j))*tmp(j)
            ENDDO
         ENDDO
         DO i = 1 , 3
            IF ( jpiv(4-i)/=4-i ) THEN
               temp = tmp(4-i)
               tmp(4-i) = tmp(jpiv(4-i))
               tmp(jpiv(4-i)) = temp
            ENDIF
         ENDDO
         X(1,1) = tmp(1)
         X(2,1) = tmp(2)
         X(1,2) = tmp(3)
         X(2,2) = tmp(4)
         Xnorm = MAX(ABS(tmp(1))+ABS(tmp(3)),ABS(tmp(2))+ABS(tmp(4)))
         GOTO 99999
      ELSE
!
!     1 by 1: TL11*X + SGN*X*TR11 = B11
!
         tau1 = Tl(1,1) + sgn*Tr(1,1)
         bet = ABS(tau1)
         IF ( bet<=smlnum ) THEN
            tau1 = smlnum
            bet = smlnum
            Info = 1
         ENDIF
!
         Scale = ONE
         gam = ABS(B(1,1))
         IF ( smlnum*gam>bet ) Scale = ONE/gam
!
         X(1,1) = (B(1,1)*Scale)/tau1
         Xnorm = ABS(X(1,1))
         RETURN
      ENDIF
!
!     Solve 2 by 2 system using complete pivoting.
!     Set pivots less than SMIN to SMIN.
!
      ipiv = IDAMAX(4,tmp,1)
      u11 = tmp(ipiv)
      IF ( ABS(u11)<=smin ) THEN
         Info = 1
         u11 = smin
      ENDIF
      u12 = tmp(locu12(ipiv))
      l21 = tmp(locl21(ipiv))/u11
      u22 = tmp(locu22(ipiv)) - u12*l21
      xswap = xswpiv(ipiv)
      bswap = bswpiv(ipiv)
      IF ( ABS(u22)<=smin ) THEN
         Info = 1
         u22 = smin
      ENDIF
      IF ( bswap ) THEN
         temp = btmp(2)
         btmp(2) = btmp(1) - l21*temp
         btmp(1) = temp
      ELSE
         btmp(2) = btmp(2) - l21*btmp(1)
      ENDIF
      Scale = ONE
      IF ( (TWO*smlnum)*ABS(btmp(2))>ABS(u22) .OR. (TWO*smlnum)         &
     &     *ABS(btmp(1))>ABS(u11) ) THEN
         Scale = HALF/MAX(ABS(btmp(1)),ABS(btmp(2)))
         btmp(1) = btmp(1)*Scale
         btmp(2) = btmp(2)*Scale
      ENDIF
      x2(2) = btmp(2)/u22
      x2(1) = btmp(1)/u11 - (u12/u11)*x2(2)
      IF ( xswap ) THEN
         temp = x2(2)
         x2(2) = x2(1)
         x2(1) = temp
      ENDIF
      X(1,1) = x2(1)
      IF ( N1==1 ) THEN
         X(1,2) = x2(2)
         Xnorm = ABS(X(1,1)) + ABS(X(1,2))
      ELSE
         X(2,1) = x2(2)
         Xnorm = MAX(ABS(X(1,1)),ABS(X(2,1)))
      ENDIF
!
!     End of DLASY2
!
99999 END SUBROUTINE DLASY2
