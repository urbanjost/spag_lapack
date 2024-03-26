!*==crqt02.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b crqt02
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE crqt02( M, N, K, A, AF, Q, R, LDA, TAU, WORK, LWORK,
!                          RWORK, RESULT )
!
!       .. Scalar Arguments ..
!       INTEGER            K, LDA, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       REAL               RESULT( * ), RWORK( * )
!       COMPLEX            A( LDA, * ), AF( LDA, * ), Q( LDA, * ),
!      $                   R( LDA, * ), TAU( * ), WORK( LWORK )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> crqt02 tests CUNGRQ, which generates an m-by-n matrix Q with
!> orthonornmal rows that is defined as the product of k elementary
!> reflectors.
!>
!> Given the RQ factorization of an m-by-n matrix A, crqt02 generates
!> the orthogonal matrix Q defined by the factorization of the last k
!> rows of A; it compares R(m-k+1:m,n-m+1:n) with
!> A(m-k+1:m,1:n)*Q(n-m+1:n,1:n)', and checks that the rows of Q are
!> orthonormal.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix Q to be generated.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix Q to be generated.
!>          N >= M >= 0.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The number of elementary reflectors whose product defines the
!>          matrix Q. M >= K >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          The m-by-n matrix A which was factorized by crqt01.
!> \endverbatim
!>
!> \param[in] AF
!> \verbatim
!>          AF is COMPLEX array, dimension (LDA,N)
!>          Details of the RQ factorization of A, as returned by CGERQF.
!>          See CGERQF for further details.
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is COMPLEX array, dimension (LDA,N)
!> \endverbatim
!>
!> \param[out] R
!> \verbatim
!>          R is COMPLEX array, dimension (LDA,M)
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the arrays A, AF, Q and L. LDA >= N.
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is COMPLEX array, dimension (M)
!>          The scalar factors of the elementary reflectors corresponding
!>          to the RQ factorization in AF.
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
!>          The dimension of the array WORK.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (M)
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is REAL array, dimension (2)
!>          The test ratios:
!>          RESULT(1) = norm( R - A*Q' ) / ( N * norm(A) * EPS )
!>          RESULT(2) = norm( I - Q*Q' ) / ( N * EPS )
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
      SUBROUTINE crqt02(M,N,K,A,Af,Q,R,Lda,Tau,Work,Lwork,Rwork,Result)
      IMPLICIT NONE
!*--CRQT02139
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
      REAL Result(*) , Rwork(*)
      COMPLEX A(Lda,*) , Af(Lda,*) , Q(Lda,*) , R(Lda,*) , Tau(*) ,     &
     &        Work(Lwork)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
      COMPLEX ROGUE
      PARAMETER (ROGUE=(-1.0E+10,-1.0E+10))
!     ..
!     .. Local Scalars ..
      INTEGER info
      REAL anorm , eps , resid
!     ..
!     .. External Functions ..
      REAL CLANGE , CLANSY , SLAMCH
EXTERNAL    CLANGE,   CLANSY,   SLAMCH
!     ..
!     .. External Subroutines ..
EXTERNAL    CGEMM,   CHERK,   CLACPY,   CLASET,   CUNGRQ
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC CMPLX , MAX , REAL
!     ..
!     .. Scalars in Common ..
      CHARACTER*32 SRNamt
!     ..
!     .. Common blocks ..
      COMMON /SRNAMC/ SRNamt
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF ( M==0 .OR. N==0 .OR. K==0 ) THEN
         Result(1) = ZERO
         Result(2) = ZERO
         RETURN
      ENDIF
!
      eps = SLAMCH('Epsilon')
!
!     Copy the last k rows of the factorization to the array Q
!
      CALL CLASET('Full',M,N,ROGUE,ROGUE,Q,Lda)
      IF ( K<N ) CALL CLACPY('Full',K,N-K,Af(M-K+1,1),Lda,Q(M-K+1,1),   &
     &                       Lda)
      IF ( K>1 ) CALL CLACPY('Lower',K-1,K-1,Af(M-K+2,N-K+1),Lda,       &
     &                       Q(M-K+2,N-K+1),Lda)
!
!     Generate the last n rows of the matrix Q
!
      SRNamt = 'CUNGRQ'
      CALL CUNGRQ(M,N,K,Q,Lda,Tau(M-K+1),Work,Lwork,info)
!
!     Copy R(m-k+1:m,n-m+1:n)
!
      CALL CLASET('Full',K,M,CMPLX(ZERO),CMPLX(ZERO),R(M-K+1,N-M+1),Lda)
      CALL CLACPY('Upper',K,K,Af(M-K+1,N-K+1),Lda,R(M-K+1,N-K+1),Lda)
!
!     Compute R(m-k+1:m,n-m+1:n) - A(m-k+1:m,1:n) * Q(n-m+1:n,1:n)'
!
      CALL CGEMM('No transpose','Conjugate transpose',K,M,N,CMPLX(-ONE),&
     &           A(M-K+1,1),Lda,Q,Lda,CMPLX(ONE),R(M-K+1,N-M+1),Lda)
!
!     Compute norm( R - A*Q' ) / ( N * norm(A) * EPS ) .
!
      anorm = CLANGE('1',K,N,A(M-K+1,1),Lda,Rwork)
      resid = CLANGE('1',K,M,R(M-K+1,N-M+1),Lda,Rwork)
      IF ( anorm>ZERO ) THEN
         Result(1) = ((resid/REAL(MAX(1,N)))/anorm)/eps
      ELSE
         Result(1) = ZERO
      ENDIF
!
!     Compute I - Q*Q'
!
      CALL CLASET('Full',M,M,CMPLX(ZERO),CMPLX(ONE),R,Lda)
      CALL CHERK('Upper','No transpose',M,N,-ONE,Q,Lda,ONE,R,Lda)
!
!     Compute norm( I - Q*Q' ) / ( N * EPS ) .
!
      resid = CLANSY('1','Upper',M,R,Lda,Rwork)
!
      Result(2) = (resid/REAL(MAX(1,N)))/eps
!
!
!     End of crqt02
!
      END SUBROUTINE crqt02
