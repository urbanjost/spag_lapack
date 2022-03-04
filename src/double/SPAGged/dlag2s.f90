!*==dlag2s.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DLAG2S converts a double precision matrix to a single precision matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLAG2S + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlag2s.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlag2s.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlag2s.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAG2S( M, N, A, LDA, SA, LDSA, INFO )
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
!> DLAG2S converts a DOUBLE PRECISION matrix, SA, to a SINGLE
!> PRECISION matrix, A.
!>
!> RMAX is the overflow for the SINGLE PRECISION arithmetic
!> DLAG2S checks that all the entries of A are between -RMAX and
!> RMAX. If not the conversion is aborted and a flag is raised.
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
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, the M-by-N coefficient matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[out] SA
!> \verbatim
!>          SA is REAL array, dimension (LDSA,N)
!>          On exit, if INFO=0, the M-by-N coefficient matrix SA; if
!>          INFO>0, the content of SA is unspecified.
!> \endverbatim
!>
!> \param[in] LDSA
!> \verbatim
!>          LDSA is INTEGER
!>          The leading dimension of the array SA.  LDSA >= max(1,M).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
!>          = 1:  an entry of the matrix A is greater than the SINGLE
!>                PRECISION overflow threshold, in this case, the content
!>                of SA in exit is unspecified.
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
!> \ingroup doubleOTHERauxiliary
!
!  =====================================================================
      SUBROUTINE DLAG2S(M,N,A,Lda,Sa,Ldsa,Info)
      USE F77KINDS                        
      USE S_SLAMCH
      IMPLICIT NONE
!*--DLAG2S114
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(OUT) , DIMENSION(Ldsa,*) :: Sa
      INTEGER , INTENT(IN) :: Ldsa
      INTEGER , INTENT(OUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , j
      REAL(R8KIND) :: rmax
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
!     .. External Functions ..
!     ..
!     .. Executable Statements ..
!
      rmax = SLAMCH('O')
      DO j = 1 , N
         DO i = 1 , M
            IF ( (A(i,j)<-rmax) .OR. (A(i,j)>rmax) ) THEN
               Info = 1
               GOTO 99999
            ENDIF
            Sa(i,j) = A(i,j)
         ENDDO
      ENDDO
      Info = 0
!
!     End of DLAG2S
!
99999 END SUBROUTINE DLAG2S
