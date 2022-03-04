!*==sla_gerpvgrw.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SLA_GERPVGRW
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLA_GERPVGRW + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sla_gerpvgrw.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sla_gerpvgrw.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sla_gerpvgrw.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       REAL FUNCTION SLA_GERPVGRW( N, NCOLS, A, LDA, AF, LDAF )
!
!       .. Scalar Arguments ..
!       INTEGER            N, NCOLS, LDA, LDAF
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), AF( LDAF, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLA_GERPVGRW computes the reciprocal pivot growth factor
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
!>          A is REAL array, dimension (LDA,N)
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
!>          AF is REAL array, dimension (LDAF,N)
!>     The factors L and U from the factorization
!>     A = P*L*U as computed by SGETRF.
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
!> \date December 2016
!
!> \ingroup realGEcomputational
!
!  =====================================================================
      FUNCTION SLA_GERPVGRW(N,Ncols,A,Lda,Af,Ldaf)
      IMPLICIT NONE
!*--SLA_GERPVGRW101
      REAL :: SLA_GERPVGRW
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Ncols
      REAL , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(IN) , DIMENSION(Ldaf,*) :: Af
      INTEGER , INTENT(IN) :: Ldaf
!
! Local variable declarations rewritten by SPAG
!
      REAL :: amax , rpvgrw , umax
      INTEGER :: i , j
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
!     .. Executable Statements ..
!
      rpvgrw = 1.0
 
      DO j = 1 , Ncols
         amax = 0.0
         umax = 0.0
         DO i = 1 , N
            amax = MAX(ABS(A(i,j)),amax)
         ENDDO
         DO i = 1 , j
            umax = MAX(ABS(Af(i,j)),umax)
         ENDDO
         IF ( umax/=0.0 ) rpvgrw = MIN(amax/umax,rpvgrw)
      ENDDO
      SLA_GERPVGRW = rpvgrw
      END FUNCTION SLA_GERPVGRW
