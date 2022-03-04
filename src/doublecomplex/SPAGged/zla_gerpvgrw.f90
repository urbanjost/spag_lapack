!*==zla_gerpvgrw.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZLA_GERPVGRW multiplies a square real matrix by a complex matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLA_GERPVGRW + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zla_gerpvgrw.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zla_gerpvgrw.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zla_gerpvgrw.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION ZLA_GERPVGRW( N, NCOLS, A, LDA, AF,
!                LDAF )
!
!       .. Scalar Arguments ..
!       INTEGER            N, NCOLS, LDA, LDAF
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * ), AF( LDAF, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>
!> ZLA_GERPVGRW computes the reciprocal pivot growth factor
!> norm(A)/norm(U). The "max absolute element" norm is used. If this is
!> much less than 1, the stability of the LU factorization of the
!> (equilibrated) matrix A could be poor. This also means that the
!> solution X, estimated condition numbers, and error bounds could be
!> unreliable.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>     The number of linear equations, i.e., the order of the
!>     matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] NCOLS
!> \verbatim
!>          NCOLS is INTEGER
!>     The number of columns of the matrix A. NCOLS >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>     On entry, the N-by-N matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>     The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] AF
!> \verbatim
!>          AF is COMPLEX*16 array, dimension (LDAF,N)
!>     The factors L and U from the factorization
!>     A = P*L*U as computed by ZGETRF.
!> \endverbatim
!>
!> \param[in] LDAF
!> \verbatim
!>          LDAF is INTEGER
!>     The leading dimension of the array AF.  LDAF >= max(1,N).
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
!> \ingroup complex16GEcomputational
!
!  =====================================================================
      FUNCTION ZLA_GERPVGRW(N,Ncols,A,Lda,Af,Ldaf)
      USE F77KINDS                        
      IMPLICIT NONE
!*--ZLA_GERPVGRW104
      REAL(R8KIND) :: ZLA_GERPVGRW
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Ncols
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Ldaf,*) :: Af
      INTEGER , INTENT(IN) :: Ldaf
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: amax , rpvgrw , umax
      REAL(R8KIND) :: CABS1
      INTEGER :: i , j
      COMPLEX(CX16KIND) :: zdum
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Statement Functions ..
!     ..
!     .. Statement Function Definitions ..
      CABS1(zdum) = ABS(DBLE(zdum)) + ABS(DIMAG(zdum))
!     ..
!     .. Executable Statements ..
!
      rpvgrw = 1.0D+0
 
      DO j = 1 , Ncols
         amax = 0.0D+0
         umax = 0.0D+0
         DO i = 1 , N
            amax = MAX(CABS1(A(i,j)),amax)
         ENDDO
         DO i = 1 , j
            umax = MAX(CABS1(Af(i,j)),umax)
         ENDDO
         IF ( umax/=0.0D+0 ) rpvgrw = MIN(amax/umax,rpvgrw)
      ENDDO
      ZLA_GERPVGRW = rpvgrw
      END FUNCTION ZLA_GERPVGRW
