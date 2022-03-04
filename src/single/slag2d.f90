!*==slag2d.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SLAG2D converts a single precision matrix to a double precision matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLAG2D + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slag2d.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slag2d.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slag2d.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLAG2D( M, N, SA, LDSA, A, LDA, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LDSA, M, N
!       ..
!       .. Array Arguments ..
!       REAL               SA( LDSA, * )
!       DOUBLE PRECISION   A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLAG2D converts a SINGLE PRECISION matrix, SA, to a DOUBLE
!> PRECISION matrix, A.
!>
!> Note that while it is possible to overflow while converting
!> from double to single, it is not possible to overflow when
!> converting from single to double.
!>
!> This is an auxiliary routine so there is no argument checking.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of lines of the matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] SA
!> \verbatim
!>          SA is REAL array, dimension (LDSA,N)
!>          On entry, the M-by-N coefficient matrix SA.
!> \endverbatim
!>
!> \param[in] LDSA
!> \verbatim
!>          LDSA is INTEGER
!>          The leading dimension of the array SA.  LDSA >= max(1,M).
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On exit, the M-by-N coefficient matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
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
!> \ingroup OTHERauxiliary
!
!  =====================================================================
      SUBROUTINE SLAG2D(M,N,Sa,Ldsa,A,Lda,Info)
      IMPLICIT NONE
!*--SLAG2D108
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Info , Lda , Ldsa , M , N
!     ..
!     .. Array Arguments ..
      REAL Sa(Ldsa,*)
      DOUBLE PRECISION A(Lda,*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER i , j
!     ..
!     .. Executable Statements ..
!
      Info = 0
      DO j = 1 , N
         DO i = 1 , M
            A(i,j) = Sa(i,j)
         ENDDO
      ENDDO
!
!     End of SLAG2D
!
      END SUBROUTINE SLAG2D
