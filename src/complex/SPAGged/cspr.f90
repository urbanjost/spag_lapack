!*==cspr.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CSPR performs the symmetrical rank-1 update of a complex symmetric packed matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CSPR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cspr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cspr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cspr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CSPR( UPLO, N, ALPHA, X, INCX, AP )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INCX, N
!       COMPLEX            ALPHA
!       ..
!       .. Array Arguments ..
!       COMPLEX            AP( * ), X( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CSPR    performs the symmetric rank 1 operation
!>
!>    A := alpha*x*x**H + A,
!>
!> where alpha is a complex scalar, x is an n element vector and A is an
!> n by n symmetric matrix, supplied in packed form.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>           On entry, UPLO specifies whether the upper or lower
!>           triangular part of the matrix A is supplied in the packed
!>           array AP as follows:
!>
!>              UPLO = 'U' or 'u'   The upper triangular part of A is
!>                                  supplied in AP.
!>
!>              UPLO = 'L' or 'l'   The lower triangular part of A is
!>                                  supplied in AP.
!>
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           On entry, N specifies the order of the matrix A.
!>           N must be at least zero.
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is COMPLEX
!>           On entry, ALPHA specifies the scalar alpha.
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is COMPLEX array, dimension at least
!>           ( 1 + ( N - 1 )*abs( INCX ) ).
!>           Before entry, the incremented array X must contain the N-
!>           element vector x.
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>           On entry, INCX specifies the increment for the elements of
!>           X. INCX must not be zero.
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in,out] AP
!> \verbatim
!>          AP is COMPLEX array, dimension at least
!>           ( ( N*( N + 1 ) )/2 ).
!>           Before entry, with  UPLO = 'U' or 'u', the array AP must
!>           contain the upper triangular part of the symmetric matrix
!>           packed sequentially, column by column, so that AP( 1 )
!>           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )
!>           and a( 2, 2 ) respectively, and so on. On exit, the array
!>           AP is overwritten by the upper triangular part of the
!>           updated matrix.
!>           Before entry, with UPLO = 'L' or 'l', the array AP must
!>           contain the lower triangular part of the symmetric matrix
!>           packed sequentially, column by column, so that AP( 1 )
!>           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )
!>           and a( 3, 1 ) respectively, and so on. On exit, the array
!>           AP is overwritten by the lower triangular part of the
!>           updated matrix.
!>           Note that the imaginary parts of the diagonal elements need
!>           not be set, they are assumed to be zero, and on exit they
!>           are set to zero.
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
!> \date December 2016
!
!> \ingroup complexOTHERauxiliary
!
!  =====================================================================
      SUBROUTINE CSPR(Uplo,N,Alpha,X,Incx,Ap)
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--CSPR138
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(IN) :: Alpha
      COMPLEX , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Ap
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , info , ix , j , jx , k , kk , kx
      COMPLEX :: temp
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
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      info = 0
      IF ( .NOT.LSAME(Uplo,'U') .AND. .NOT.LSAME(Uplo,'L') ) THEN
         info = 1
      ELSEIF ( N<0 ) THEN
         info = 2
      ELSEIF ( Incx==0 ) THEN
         info = 5
      ENDIF
      IF ( info/=0 ) THEN
         CALL XERBLA('CSPR  ',info)
         RETURN
      ENDIF
!
!     Quick return if possible.
!
      IF ( (N==0) .OR. (Alpha==ZERO) ) RETURN
!
!     Set the start point in X if the increment is not unity.
!
      IF ( Incx<=0 ) THEN
         kx = 1 - (N-1)*Incx
      ELSEIF ( Incx/=1 ) THEN
         kx = 1
      ENDIF
!
!     Start the operations. In this version the elements of the array AP
!     are accessed sequentially with one pass through AP.
!
      kk = 1
      IF ( LSAME(Uplo,'U') ) THEN
!
!        Form  A  when upper triangle is stored in AP.
!
         IF ( Incx==1 ) THEN
            DO j = 1 , N
               IF ( X(j)/=ZERO ) THEN
                  temp = Alpha*X(j)
                  k = kk
                  DO i = 1 , j - 1
                     Ap(k) = Ap(k) + X(i)*temp
                     k = k + 1
                  ENDDO
                  Ap(kk+j-1) = Ap(kk+j-1) + X(j)*temp
               ELSE
                  Ap(kk+j-1) = Ap(kk+j-1)
               ENDIF
               kk = kk + j
            ENDDO
         ELSE
            jx = kx
            DO j = 1 , N
               IF ( X(jx)/=ZERO ) THEN
                  temp = Alpha*X(jx)
                  ix = kx
                  DO k = kk , kk + j - 2
                     Ap(k) = Ap(k) + X(ix)*temp
                     ix = ix + Incx
                  ENDDO
                  Ap(kk+j-1) = Ap(kk+j-1) + X(jx)*temp
               ELSE
                  Ap(kk+j-1) = Ap(kk+j-1)
               ENDIF
               jx = jx + Incx
               kk = kk + j
            ENDDO
         ENDIF
!
!        Form  A  when lower triangle is stored in AP.
!
      ELSEIF ( Incx==1 ) THEN
         DO j = 1 , N
            IF ( X(j)/=ZERO ) THEN
               temp = Alpha*X(j)
               Ap(kk) = Ap(kk) + temp*X(j)
               k = kk + 1
               DO i = j + 1 , N
                  Ap(k) = Ap(k) + X(i)*temp
                  k = k + 1
               ENDDO
            ELSE
               Ap(kk) = Ap(kk)
            ENDIF
            kk = kk + N - j + 1
         ENDDO
      ELSE
         jx = kx
         DO j = 1 , N
            IF ( X(jx)/=ZERO ) THEN
               temp = Alpha*X(jx)
               Ap(kk) = Ap(kk) + temp*X(jx)
               ix = jx
               DO k = kk + 1 , kk + N - j
                  ix = ix + Incx
                  Ap(k) = Ap(k) + X(ix)*temp
               ENDDO
            ELSE
               Ap(kk) = Ap(kk)
            ENDIF
            jx = jx + Incx
            kk = kk + N - j + 1
         ENDDO
      ENDIF
!
!
!     End of CSPR
!
      END SUBROUTINE CSPR
