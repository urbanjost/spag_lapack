!*==zlapll.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZLAPLL measures the linear dependence of two vectors.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLAPLL + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlapll.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlapll.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlapll.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLAPLL( N, X, INCX, Y, INCY, SSMIN )
!
!       .. Scalar Arguments ..
!       INTEGER            INCX, INCY, N
!       DOUBLE PRECISION   SSMIN
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         X( * ), Y( * )
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
!>          X is COMPLEX*16 array, dimension (1+(N-1)*INCX)
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
!>          Y is COMPLEX*16 array, dimension (1+(N-1)*INCY)
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
!>          SSMIN is DOUBLE PRECISION
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
!> \ingroup complex16OTHERauxiliary
!
!  =====================================================================
      SUBROUTINE ZLAPLL(N,X,Incx,Y,Incy,Ssmin)
      USE F77KINDS                        
      USE S_DLAS2
      USE S_ZAXPY
      USE S_ZDOTC
      USE S_ZLARFG
      IMPLICIT NONE
!*--ZLAPLL109
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0)
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: X
      INTEGER :: Incx
      COMPLEX(CX16KIND) , DIMENSION(*) :: Y
      INTEGER :: Incy
      REAL(R8KIND) :: Ssmin
!
! Local variable declarations rewritten by SPAG
!
      COMPLEX(CX16KIND) :: a11 , a12 , a22 , c , tau
      REAL(R8KIND) :: ssmax
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
      CALL ZLARFG(N,X(1),X(1+Incx),Incx,tau)
      a11 = X(1)
      X(1) = CONE
!
      c = -DCONJG(tau)*ZDOTC(N,X,Incx,Y,Incy)
      CALL ZAXPY(N,c,X,Incx,Y,Incy)
!
      CALL ZLARFG(N-1,Y(1+Incy),Y(1+2*Incy),Incy,tau)
!
      a12 = Y(1)
      a22 = Y(1+Incy)
!
!     Compute the SVD of 2-by-2 Upper triangular matrix.
!
      CALL DLAS2(ABS(a11),ABS(a12),ABS(a22),Ssmin,ssmax)
!
!
!     End of ZLAPLL
!
      END SUBROUTINE ZLAPLL
