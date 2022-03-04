!*==slapll.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SLAPLL measures the linear dependence of two vectors.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLAPLL + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slapll.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slapll.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slapll.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLAPLL( N, X, INCX, Y, INCY, SSMIN )
!
!       .. Scalar Arguments ..
!       INTEGER            INCX, INCY, N
!       REAL               SSMIN
!       ..
!       .. Array Arguments ..
!       REAL               X( * ), Y( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> Given two column vectors X and Y, let
!>
!>                      A = ( X Y ).
!>
!> The subroutine first computes the QR factorization of A = Q*R,
!> and then computes the SVD of the 2-by-2 upper triangular matrix R.
!> The smaller singular value of R is returned in SSMIN, which is used
!> as the measurement of the linear dependency of the vectors X and Y.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The length of the vectors X and Y.
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is REAL array,
!>                         dimension (1+(N-1)*INCX)
!>          On entry, X contains the N-vector X.
!>          On exit, X is overwritten.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>          The increment between successive elements of X. INCX > 0.
!> \endverbatim
!>
!> \param[in,out] Y
!> \verbatim
!>          Y is REAL array,
!>                         dimension (1+(N-1)*INCY)
!>          On entry, Y contains the N-vector Y.
!>          On exit, Y is overwritten.
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>          The increment between successive elements of Y. INCY > 0.
!> \endverbatim
!>
!> \param[out] SSMIN
!> \verbatim
!>          SSMIN is REAL
!>          The smallest singular value of the N-by-2 matrix A = ( X Y ).
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
!> \ingroup realOTHERauxiliary
!
!  =====================================================================
      SUBROUTINE SLAPLL(N,X,Incx,Y,Incy,Ssmin)
      USE S_SAXPY
      USE S_SDOT
      USE S_SLARFG
      USE S_SLAS2
      IMPLICIT NONE
!*--SLAPLL110
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER :: N
      REAL , INTENT(INOUT) , DIMENSION(*) :: X
      INTEGER :: Incx
      REAL , DIMENSION(*) :: Y
      INTEGER :: Incy
      REAL :: Ssmin
!
! Local variable declarations rewritten by SPAG
!
      REAL :: a11 , a12 , a22 , c , ssmax , tau
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
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF ( N<=1 ) THEN
         Ssmin = ZERO
         RETURN
      ENDIF
!
!     Compute the QR factorization of the N-by-2 matrix ( X Y )
!
      CALL SLARFG(N,X(1),X(1+Incx),Incx,tau)
      a11 = X(1)
      X(1) = ONE
!
      c = -tau*SDOT(N,X,Incx,Y,Incy)
      CALL SAXPY(N,c,X,Incx,Y,Incy)
!
      CALL SLARFG(N-1,Y(1+Incy),Y(1+2*Incy),Incy,tau)
!
      a12 = Y(1)
      a22 = Y(1+Incy)
!
!     Compute the SVD of 2-by-2 Upper triangular matrix.
!
      CALL SLAS2(a11,a12,a22,Ssmin,ssmax)
!
!
!     End of SLAPLL
!
      END SUBROUTINE SLAPLL
