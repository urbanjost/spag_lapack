!*==ilazlr.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b ILAZLR scans a matrix for its last non-zero row.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ILAZLR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ilazlr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ilazlr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ilazlr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       INTEGER FUNCTION ILAZLR( M, N, A, LDA )
!
!       .. Scalar Arguments ..
!       INTEGER            M, N, LDA
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ILAZLR scans A for its last non-zero row.
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
!>          A is COMPLEX*16 array, dimension (LDA,N)
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
!> \ingroup complex16OTHERauxiliary
!
!  =====================================================================
      INTEGER FUNCTION ILAZLR(M,N,A,Lda)
      IMPLICIT NONE
!*--ILAZLR82
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER M , N , Lda
!     ..
!     .. Array Arguments ..
      COMPLEX*16 A(Lda,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16 ZERO
      PARAMETER (ZERO=(0.0D+0,0.0D+0))
!     ..
!     .. Local Scalars ..
      INTEGER i , j
!     ..
!     .. Executable Statements ..
!
!     Quick test for the common case where one corner is non-zero.
      IF ( M==0 ) THEN
         ILAZLR = M
      ELSEIF ( A(M,1)/=ZERO .OR. A(M,N)/=ZERO ) THEN
         ILAZLR = M
      ELSE
!     Scan up each column tracking the last zero row seen.
         ILAZLR = 0
         DO j = 1 , N
            i = M
            DO WHILE ( (A(MAX(i,1),j)==ZERO) .AND. (i>=1) )
               i = i - 1
            ENDDO
            ILAZLR = MAX(ILAZLR,i)
         ENDDO
      ENDIF
      END FUNCTION ILAZLR
