!*==sqrt11.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
 
!> \brief \b SQRT11
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       REAL             FUNCTION SQRT11( M, K, A, LDA, TAU, WORK, LWORK )
!
!       .. Scalar Arguments ..
!       INTEGER            K, LDA, LWORK, M
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), TAU( * ), WORK( LWORK )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SQRT11 computes the test ratio
!>
!>       || Q'*Q - I || / (eps * m)
!>
!> where the orthogonal matrix Q is represented as a product of
!> elementary transformations.  Each transformation has the form
!>
!>    H(k) = I - tau(k) v(k) v(k)'
!>
!> where tau(k) is stored in TAU(k) and v(k) is an m-vector of the form
!> [ 0 ... 0 1 x(k) ]', where x(k) is a vector of length m-k stored
!> in A(k+1:m,k).
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
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The number of columns of A whose subdiagonal entries
!>          contain information about orthogonal transformations.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is REAL array, dimension (LDA,K)
!>          The (possibly partial) output of a QR reduction routine.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is REAL array, dimension (K)
!>          The scaling factors tau for the elementary transformations as
!>          computed by the QR factorization routine.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (LWORK)
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The length of the array WORK.  LWORK >= M*M + M.
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
!> \ingroup single_lin
!
!  =====================================================================
      REAL FUNCTION SQRT11(M,K,A,Lda,Tau,Work,Lwork)
      IMPLICIT NONE
!*--SQRT11103
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER K , Lda , Lwork , M
!     ..
!     .. Array Arguments ..
      REAL A(Lda,*) , Tau(*) , Work(Lwork)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E0,ONE=1.0E0)
!     ..
!     .. Local Scalars ..
      INTEGER info , j
!     ..
!     .. External Functions ..
      REAL SLAMCH , SLANGE
      EXTERNAL SLAMCH , SLANGE
!     ..
!     .. External Subroutines ..
      EXTERNAL SLASET , SORM2R , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC REAL
!     ..
!     .. Local Arrays ..
      REAL rdummy(1)
!     ..
!     .. Executable Statements ..
!
      SQRT11 = ZERO
!
!     Test for sufficient workspace
!
      IF ( Lwork<M*M+M ) THEN
         CALL XERBLA('SQRT11',7)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( M<=0 ) RETURN
!
      CALL SLASET('Full',M,M,ZERO,ONE,Work,M)
!
!     Form Q
!
      CALL SORM2R('Left','No transpose',M,M,K,A,Lda,Tau,Work,M,         &
     &            Work(M*M+1),info)
!
!     Form Q'*Q
!
      CALL SORM2R('Left','Transpose',M,M,K,A,Lda,Tau,Work,M,Work(M*M+1),&
     &            info)
!
      DO j = 1 , M
         Work((j-1)*M+j) = Work((j-1)*M+j) - ONE
      ENDDO
!
      SQRT11 = SLANGE('One-norm',M,M,Work,M,rdummy)                     &
     &         /(REAL(M)*SLAMCH('Epsilon'))
!
!
!     End of SQRT11
!
      END FUNCTION SQRT11
