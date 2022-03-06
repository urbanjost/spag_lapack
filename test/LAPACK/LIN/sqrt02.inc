!*==sqrt02.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b sqrt02
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE sqrt02( M, N, K, A, AF, Q, R, LDA, TAU, WORK, LWORK,
!                          RWORK, RESULT )
!
!       .. Scalar Arguments ..
!       INTEGER            K, LDA, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), AF( LDA, * ), Q( LDA, * ),
!      $                   R( LDA, * ), RESULT( * ), RWORK( * ), TAU( * ),
!      $                   WORK( LWORK )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> sqrt02 tests SORGQR, which generates an m-by-n matrix Q with
!> orthonornmal columns that is defined as the product of k elementary
!> reflectors.
!>
!> Given the QR factorization of an m-by-n matrix A, sqrt02 generates
!> the orthogonal matrix Q defined by the factorization of the first k
!> columns of A; it compares R(1:n,1:k) with Q(1:m,1:n)'*A(1:m,1:k),
!> and checks that the columns of Q are orthonormal.
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
!>          M >= N >= 0.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The number of elementary reflectors whose product defines the
!>          matrix Q. N >= K >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>          The m-by-n matrix A which was factorized by sqrt01.
!> \endverbatim
!>
!> \param[in] AF
!> \verbatim
!>          AF is REAL array, dimension (LDA,N)
!>          Details of the QR factorization of A, as returned by SGEQRF.
!>          See SGEQRF for further details.
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is REAL array, dimension (LDA,N)
!> \endverbatim
!>
!> \param[out] R
!> \verbatim
!>          R is REAL array, dimension (LDA,N)
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the arrays A, AF, Q and R. LDA >= M.
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is REAL array, dimension (N)
!>          The scalar factors of the elementary reflectors corresponding
!>          to the QR factorization in AF.
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
!>          RESULT(1) = norm( R - Q'*A ) / ( M * norm(A) * EPS )
!>          RESULT(2) = norm( I - Q'*Q ) / ( M * EPS )
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
      SUBROUTINE sqrt02(M,N,K,A,Af,Q,R,Lda,Tau,Work,Lwork,Rwork,Result)
      IMPLICIT NONE
!*--SQRT02138
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
      REAL A(Lda,*) , Af(Lda,*) , Q(Lda,*) , R(Lda,*) , Result(*) ,     &
     &     Rwork(*) , Tau(*) , Work(Lwork)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
      REAL ROGUE
      PARAMETER (ROGUE=-1.0E+10)
!     ..
!     .. Local Scalars ..
      INTEGER info
      REAL anorm , eps , resid
!     ..
!     .. External Functions ..
      REAL SLAMCH , SLANGE , SLANSY
EXTERNAL    SLAMCH,   SLANGE,   SLANSY
!     ..
!     .. External Subroutines ..
EXTERNAL    SGEMM,   SLACPY,   SLASET,   SORGQR,   SSYRK
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , REAL
!     ..
!     .. Scalars in Common ..
      CHARACTER*32 SRNamt
!     ..
!     .. Common blocks ..
      COMMON /SRNAMC/ SRNamt
!     ..
!     .. Executable Statements ..
!
      eps = SLAMCH('Epsilon')
!
!     Copy the first k columns of the factorization to the array Q
!
      CALL SLASET('Full',M,N,ROGUE,ROGUE,Q,Lda)
      CALL SLACPY('Lower',M-1,K,Af(2,1),Lda,Q(2,1),Lda)
!
!     Generate the first n columns of the matrix Q
!
      SRNamt = 'SORGQR'
      CALL SORGQR(M,N,K,Q,Lda,Tau,Work,Lwork,info)
!
!     Copy R(1:n,1:k)
!
      CALL SLASET('Full',N,K,ZERO,ZERO,R,Lda)
      CALL SLACPY('Upper',N,K,Af,Lda,R,Lda)
!
!     Compute R(1:n,1:k) - Q(1:m,1:n)' * A(1:m,1:k)
!
      CALL SGEMM('Transpose','No transpose',N,K,M,-ONE,Q,Lda,A,Lda,ONE, &
     &           R,Lda)
!
!     Compute norm( R - Q'*A ) / ( M * norm(A) * EPS ) .
!
      anorm = SLANGE('1',M,K,A,Lda,Rwork)
      resid = SLANGE('1',N,K,R,Lda,Rwork)
      IF ( anorm>ZERO ) THEN
         Result(1) = ((resid/REAL(MAX(1,M)))/anorm)/eps
      ELSE
         Result(1) = ZERO
      ENDIF
!
!     Compute I - Q'*Q
!
      CALL SLASET('Full',N,N,ZERO,ONE,R,Lda)
      CALL SSYRK('Upper','Transpose',N,M,-ONE,Q,Lda,ONE,R,Lda)
!
!     Compute norm( I - Q'*Q ) / ( M * EPS ) .
!
      resid = SLANSY('1','Upper',N,R,Lda,Rwork)
!
      Result(2) = (resid/REAL(MAX(1,M)))/eps
!
!
!     End of sqrt02
!
      END SUBROUTINE sqrt02
