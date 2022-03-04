!*==dlarscl2.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DLARSCL2 performs reciprocal diagonal scaling on a vector.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLARSCL2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarscl2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarscl2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarscl2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLARSCL2 ( M, N, D, X, LDX )
!
!       .. Scalar Arguments ..
!       INTEGER            M, N, LDX
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   D( * ), X( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLARSCL2 performs a reciprocal diagonal scaling on an vector:
!>   x <-- inv(D) * x
!> where the diagonal matrix D is stored as a vector.
!>
!> Eventually to be replaced by BLAS_dge_diag_scale in the new BLAS
!> standard.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>     The number of rows of D and X. M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>     The number of columns of X. N >= 0.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (M)
!>     Diagonal matrix D, stored as a vector of length M.
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension (LDX,N)
!>     On entry, the vector X to be scaled by D.
!>     On exit, the scaled vector.
!> \endverbatim
!>
!> \param[in] LDX
!> \verbatim
!>          LDX is INTEGER
!>     The leading dimension of the vector X. LDX >= M.
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
!> \ingroup doubleOTHERcomputational
!
!  =====================================================================
      SUBROUTINE DLARSCL2(M,N,D,X,Ldx)
      USE F77KINDS                        
      IMPLICIT NONE
!*--DLARSCL295
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
!
! Local variable declarations rewritten by SPAG
!
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
!     .. Executable Statements ..
!
      DO j = 1 , N
         DO i = 1 , M
            X(i,j) = X(i,j)/D(i)
         ENDDO
      ENDDO
 
      END SUBROUTINE DLARSCL2
