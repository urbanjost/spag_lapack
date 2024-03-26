!*==zqrt11.f90  processed by SPAG 7.51RB at 17:35 on  4 Mar 2022
 
!> \brief \b zqrt11
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION zqrt11( M, K, A, LDA, TAU, WORK, LWORK )
!
!       .. Scalar Arguments ..
!       INTEGER            K, LDA, LWORK, M
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * ), TAU( * ), WORK( LWORK )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> zqrt11 computes the test ratio
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
!>          A is COMPLEX*16 array, dimension (LDA,K)
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
!>          TAU is COMPLEX*16 array, dimension (K)
!>          The scaling factors tau for the elementary transformations as
!>          computed by the QR factorization routine.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (LWORK)
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
!> \ingroup complex16_lin
!
!  =====================================================================
      DOUBLE PRECISION FUNCTION zqrt11(M,K,A,Lda,Tau,Work,Lwork)
      IMPLICIT NONE
!*--ZQRT11103
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
      COMPLEX*16 A(Lda,*) , Tau(*) , Work(Lwork)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
!     ..
!     .. Local Scalars ..
      INTEGER info , j
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH , ZLANGE
EXTERNAL    DLAMCH,   ZLANGE
!     ..
!     .. External Subroutines ..
EXTERNAL    ZLASET,   ZUNM2R
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DBLE , DCMPLX
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION rdummy(1)
!     ..
!     .. Executable Statements ..
!
      zqrt11 = ZERO
!
!     Test for sufficient workspace
!
      IF ( Lwork<M*M+M ) THEN
         CALL xerbla('zqrt11',7)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( M<=0 ) RETURN
!
      CALL ZLASET('Full',M,M,DCMPLX(ZERO),DCMPLX(ONE),Work,M)
!
!     Form Q
!
      CALL ZUNM2R('Left','No transpose',M,M,K,A,Lda,Tau,Work,M,         &
     &            Work(M*M+1),info)
!
!     Form Q'*Q
!
      CALL ZUNM2R('Left','Conjugate transpose',M,M,K,A,Lda,Tau,Work,M,  &
     &            Work(M*M+1),info)
!
      DO j = 1 , M
         Work((j-1)*M+j) = Work((j-1)*M+j) - ONE
      ENDDO
!
      zqrt11 = ZLANGE('One-norm',M,M,Work,M,rdummy)                     &
     &         /(DBLE(M)*DLAMCH('Epsilon'))
!
!
!     End of zqrt11
!
      END FUNCTION zqrt11
