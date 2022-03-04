!*==cspmv.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b CSPMV computes a matrix-vector product for complex vectors using a complex symmetric packed matrix
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CSPMV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cspmv.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cspmv.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cspmv.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CSPMV( UPLO, N, ALPHA, AP, X, INCX, BETA, Y, INCY )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INCX, INCY, N
!       COMPLEX            ALPHA, BETA
!       ..
!       .. Array Arguments ..
!       COMPLEX            AP( * ), X( * ), Y( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CSPMV  performs the matrix-vector operation
!>
!>    y := alpha*A*x + beta*y,
!>
!> where alpha and beta are scalars, x and y are n element vectors and
!> A is an n by n symmetric matrix, supplied in packed form.
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
!> \param[in] AP
!> \verbatim
!>          AP is COMPLEX array, dimension at least
!>           ( ( N*( N + 1 ) )/2 ).
!>           Before entry, with UPLO = 'U' or 'u', the array AP must
!>           contain the upper triangular part of the symmetric matrix
!>           packed sequentially, column by column, so that AP( 1 )
!>           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )
!>           and a( 2, 2 ) respectively, and so on.
!>           Before entry, with UPLO = 'L' or 'l', the array AP must
!>           contain the lower triangular part of the symmetric matrix
!>           packed sequentially, column by column, so that AP( 1 )
!>           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )
!>           and a( 3, 1 ) respectively, and so on.
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
!> \param[in] BETA
!> \verbatim
!>          BETA is COMPLEX
!>           On entry, BETA specifies the scalar beta. When BETA is
!>           supplied as zero then Y need not be set on input.
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in,out] Y
!> \verbatim
!>          Y is COMPLEX array, dimension at least
!>           ( 1 + ( N - 1 )*abs( INCY ) ).
!>           Before entry, the incremented array Y must contain the n
!>           element vector y. On exit, Y is overwritten by the updated
!>           vector y.
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>           On entry, INCY specifies the increment for the elements of
!>           Y. INCY must not be zero.
!>           Unchanged on exit.
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
      SUBROUTINE CSPMV(Uplo,N,Alpha,Ap,X,Incx,Beta,Y,Incy)
      IMPLICIT NONE
!*--CSPMV155
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER Incx , Incy , N
      COMPLEX Alpha , Beta
!     ..
!     .. Array Arguments ..
      COMPLEX Ap(*) , X(*) , Y(*)
!     ..
!
! =====================================================================
!
!     .. Parameters ..
      COMPLEX ONE
      PARAMETER (ONE=(1.0E+0,0.0E+0))
      COMPLEX ZERO
      PARAMETER (ZERO=(0.0E+0,0.0E+0))
!     ..
!     .. Local Scalars ..
      INTEGER i , info , ix , iy , j , jx , jy , k , kk , kx , ky
      COMPLEX temp1 , temp2
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA
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
         info = 6
      ELSEIF ( Incy==0 ) THEN
         info = 9
      ENDIF
      IF ( info/=0 ) THEN
         CALL XERBLA('CSPMV ',info)
         RETURN
      ENDIF
!
!     Quick return if possible.
!
      IF ( (N==0) .OR. ((Alpha==ZERO) .AND. (Beta==ONE)) ) RETURN
!
!     Set up the start points in  X  and  Y.
!
      IF ( Incx>0 ) THEN
         kx = 1
      ELSE
         kx = 1 - (N-1)*Incx
      ENDIF
      IF ( Incy>0 ) THEN
         ky = 1
      ELSE
         ky = 1 - (N-1)*Incy
      ENDIF
!
!     Start the operations. In this version the elements of the array AP
!     are accessed sequentially with one pass through AP.
!
!     First form  y := beta*y.
!
      IF ( Beta/=ONE ) THEN
         IF ( Incy/=1 ) THEN
            iy = ky
            IF ( Beta==ZERO ) THEN
               DO i = 1 , N
                  Y(iy) = ZERO
                  iy = iy + Incy
               ENDDO
            ELSE
               DO i = 1 , N
                  Y(iy) = Beta*Y(iy)
                  iy = iy + Incy
               ENDDO
            ENDIF
         ELSEIF ( Beta==ZERO ) THEN
            DO i = 1 , N
               Y(i) = ZERO
            ENDDO
         ELSE
            DO i = 1 , N
               Y(i) = Beta*Y(i)
            ENDDO
         ENDIF
      ENDIF
      IF ( Alpha==ZERO ) RETURN
      kk = 1
      IF ( LSAME(Uplo,'U') ) THEN
!
!        Form  y  when AP contains the upper triangle.
!
         IF ( (Incx==1) .AND. (Incy==1) ) THEN
            DO j = 1 , N
               temp1 = Alpha*X(j)
               temp2 = ZERO
               k = kk
               DO i = 1 , j - 1
                  Y(i) = Y(i) + temp1*Ap(k)
                  temp2 = temp2 + Ap(k)*X(i)
                  k = k + 1
               ENDDO
               Y(j) = Y(j) + temp1*Ap(kk+j-1) + Alpha*temp2
               kk = kk + j
            ENDDO
         ELSE
            jx = kx
            jy = ky
            DO j = 1 , N
               temp1 = Alpha*X(jx)
               temp2 = ZERO
               ix = kx
               iy = ky
               DO k = kk , kk + j - 2
                  Y(iy) = Y(iy) + temp1*Ap(k)
                  temp2 = temp2 + Ap(k)*X(ix)
                  ix = ix + Incx
                  iy = iy + Incy
               ENDDO
               Y(jy) = Y(jy) + temp1*Ap(kk+j-1) + Alpha*temp2
               jx = jx + Incx
               jy = jy + Incy
               kk = kk + j
            ENDDO
         ENDIF
!
!        Form  y  when AP contains the lower triangle.
!
      ELSEIF ( (Incx==1) .AND. (Incy==1) ) THEN
         DO j = 1 , N
            temp1 = Alpha*X(j)
            temp2 = ZERO
            Y(j) = Y(j) + temp1*Ap(kk)
            k = kk + 1
            DO i = j + 1 , N
               Y(i) = Y(i) + temp1*Ap(k)
               temp2 = temp2 + Ap(k)*X(i)
               k = k + 1
            ENDDO
            Y(j) = Y(j) + Alpha*temp2
            kk = kk + (N-j+1)
         ENDDO
      ELSE
         jx = kx
         jy = ky
         DO j = 1 , N
            temp1 = Alpha*X(jx)
            temp2 = ZERO
            Y(jy) = Y(jy) + temp1*Ap(kk)
            ix = jx
            iy = jy
            DO k = kk + 1 , kk + N - j
               ix = ix + Incx
               iy = iy + Incy
               Y(iy) = Y(iy) + temp1*Ap(k)
               temp2 = temp2 + Ap(k)*X(ix)
            ENDDO
            Y(jy) = Y(jy) + Alpha*temp2
            jx = jx + Incx
            jy = jy + Incy
            kk = kk + (N-j+1)
         ENDDO
      ENDIF
!
!
!     End of CSPMV
!
      END SUBROUTINE CSPMV
