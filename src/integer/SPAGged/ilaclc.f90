!*==ilaclc.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ILACLC scans a matrix for its last non-zero column.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ILACLC + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ilaclc.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ilaclc.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ilaclc.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       INTEGER FUNCTION ILACLC( M, N, A, LDA )
!
!       .. Scalar Arguments ..
!       INTEGER            M, N, LDA
!       ..
!       .. Array Arguments ..
!       COMPLEX            A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ILACLC scans A for its last non-zero column.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          The m by n matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A. LDA >= max(1,M).
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
      FUNCTION ILACLC(M,N,A,Lda)
      IMPLICIT NONE
!*--ILACLC82
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
      INTEGER :: ILACLC
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i
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
!     .. Executable Statements ..
!
!     Quick test for the common case where one corner is non-zero.
      IF ( N==0 ) THEN
         ILACLC = N
      ELSEIF ( A(1,N)/=ZERO .OR. A(M,N)/=ZERO ) THEN
         ILACLC = N
      ELSE
!     Now scan each column from the end, returning with the first non-zero.
         DO ILACLC = N , 1 , -1
            DO i = 1 , M
               IF ( A(i,ILACLC)/=ZERO ) RETURN
            ENDDO
         ENDDO
      ENDIF
      END FUNCTION ILACLC
