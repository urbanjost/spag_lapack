!*==cqpt01.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b cqpt01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       REAL             FUNCTION cqpt01( M, N, K, A, AF, LDA, TAU, JPVT,
!                        WORK, LWORK )
!
!       .. Scalar Arguments ..
!       INTEGER            K, LDA, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       INTEGER            JPVT( * )
!       COMPLEX            A( LDA, * ), AF( LDA, * ), TAU( * ),
!      $                   WORK( LWORK )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> cqpt01 tests the QR-factorization with pivoting of a matrix A.  The
!> array AF contains the (possibly partial) QR-factorization of A, where
!> the upper triangle of AF(1:k,1:k) is a partial triangular factor,
!> the entries below the diagonal in the first k columns are the
!> Householder vectors, and the rest of AF contains a partially updated
!> matrix.
!>
!> This function returns ||A*P - Q*R||/(||norm(A)||*eps*M)
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrices A and AF.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrices A and AF.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The number of columns of AF that have been reduced
!>          to upper triangular form.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA, N)
!>          The original matrix A.
!> \endverbatim
!>
!> \param[in] AF
!> \verbatim
!>          AF is COMPLEX array, dimension (LDA,N)
!>          The (possibly partial) output of CGEQPF.  The upper triangle
!>          of AF(1:k,1:k) is a partial triangular factor, the entries
!>          below the diagonal in the first k columns are the Householder
!>          vectors, and the rest of AF contains a partially updated
!>          matrix.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the arrays A and AF.
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is COMPLEX array, dimension (K)
!>          Details of the Householder transformations as returned by
!>          CGEQPF.
!> \endverbatim
!>
!> \param[in] JPVT
!> \verbatim
!>          JPVT is INTEGER array, dimension (N)
!>          Pivot information as returned by CGEQPF.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (LWORK)
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The length of the array WORK.  LWORK >= M*N+N.
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
!> \ingroup complex_lin
!
!  =====================================================================
      REAL FUNCTION cqpt01(M,N,K,A,Af,Lda,Tau,Jpvt,Work,Lwork)
      IMPLICIT NONE
!*--CQPT01123
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER K , Lda , Lwork , M , N
!     ..
!     .. Array Arguments ..
      INTEGER Jpvt(*)
      COMPLEX A(Lda,*) , Af(Lda,*) , Tau(*) , Work(Lwork)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E0,ONE=1.0E0)
!     ..
!     .. Local Scalars ..
      INTEGER i , info , j
      REAL norma
!     ..
!     .. Local Arrays ..
      REAL rwork(1)
!     ..
!     .. External Functions ..
      REAL CLANGE , SLAMCH
EXTERNAL    CLANGE,   SLAMCH
!     ..
!     .. External Subroutines ..
EXTERNAL    CAXPY,   CCOPY,   CUNMQR
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC CMPLX , MAX , MIN , REAL
!     ..
!     .. Executable Statements ..
!
      cqpt01 = ZERO
!
!     Test if there is enough workspace
!
      IF ( Lwork<M*N+N ) THEN
         CALL xerbla('cqpt01',10)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( M<=0 .OR. N<=0 ) RETURN
!
      norma = CLANGE('One-norm',M,N,A,Lda,rwork)
!
      DO j = 1 , K
         DO i = 1 , MIN(j,M)
            Work((j-1)*M+i) = Af(i,j)
         ENDDO
         DO i = j + 1 , M
            Work((j-1)*M+i) = ZERO
         ENDDO
      ENDDO
      DO j = K + 1 , N
         CALL CCOPY(M,Af(1,j),1,Work((j-1)*M+1),1)
      ENDDO
!
      CALL CUNMQR('Left','No transpose',M,N,K,Af,Lda,Tau,Work,M,        &
     &            Work(M*N+1),Lwork-M*N,info)
!
      DO j = 1 , N
!
!        Compare i-th column of QR and jpvt(i)-th column of A
!
         CALL CAXPY(M,CMPLX(-ONE),A(1,Jpvt(j)),1,Work((j-1)*M+1),1)
      ENDDO
!
      cqpt01 = CLANGE('One-norm',M,N,Work,M,rwork)                      &
     &         /(REAL(MAX(M,N))*SLAMCH('Epsilon'))
      IF ( norma/=ZERO ) cqpt01 = cqpt01/norma
!
!
!     End of cqpt01
!
      END FUNCTION cqpt01
