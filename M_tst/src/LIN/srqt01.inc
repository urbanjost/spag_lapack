!*==srqt01.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b srqt01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE srqt01( M, N, A, AF, Q, R, LDA, TAU, WORK, LWORK,
!                          RWORK, RESULT )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LWORK, M, N
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
!> srqt01 tests SGERQF, which computes the RQ factorization of an m-by-n
!> matrix A, and partially tests SORGRQ which forms the n-by-n
!> orthogonal matrix Q.
!>
!> srqt01 compares R with A*Q', and checks that Q is orthogonal.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= 0.
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
!>          A is REAL array, dimension (LDA,N)
!>          The m-by-n matrix A.
!> \endverbatim
!>
!> \param[out] AF
!> \verbatim
!>          AF is REAL array, dimension (LDA,N)
!>          Details of the RQ factorization of A, as returned by SGERQF.
!>          See SGERQF for further details.
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is REAL array, dimension (LDA,N)
!>          The n-by-n orthogonal matrix Q.
!> \endverbatim
!>
!> \param[out] R
!> \verbatim
!>          R is REAL array, dimension (LDA,max(M,N))
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the arrays A, AF, Q and L.
!>          LDA >= max(M,N).
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is REAL array, dimension (min(M,N))
!>          The scalar factors of the elementary reflectors, as returned
!>          by SGERQF.
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
!>          RWORK is REAL array, dimension (max(M,N))
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
!> \ingroup single_lin
!
!  =====================================================================
      SUBROUTINE srqt01(M,N,A,Af,Q,R,Lda,Tau,Work,Lwork,Rwork,Result)
      IMPLICIT NONE
!*--SRQT01129
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Lda , Lwork , M , N
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
      INTEGER info , minmn
      REAL anorm , eps , resid
!     ..
!     .. External Functions ..
      REAL SLAMCH , SLANGE , SLANSY
EXTERNAL    SLAMCH,   SLANGE,   SLANSY
!     ..
!     .. External Subroutines ..
EXTERNAL    SGEMM,   SGERQF,   SLACPY,   SLASET,   SORGRQ,   SSYRK
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN , REAL
!     ..
!     .. Scalars in Common ..
      CHARACTER*32 SRNamt
!     ..
!     .. Common blocks ..
      COMMON /SRNAMC/ SRNamt
!     ..
!     .. Executable Statements ..
!
      minmn = MIN(M,N)
      eps = SLAMCH('Epsilon')
!
!     Copy the matrix A to the array AF.
!
      CALL SLACPY('Full',M,N,A,Lda,Af,Lda)
!
!     Factorize the matrix A in the array AF.
!
      SRNamt = 'SGERQF'
      CALL SGERQF(M,N,Af,Lda,Tau,Work,Lwork,info)
!
!     Copy details of Q
!
      CALL SLASET('Full',N,N,ROGUE,ROGUE,Q,Lda)
      IF ( M<=N ) THEN
         IF ( M>0 .AND. M<N ) CALL SLACPY('Full',M,N-M,Af,Lda,Q(N-M+1,1)&
     &        ,Lda)
         IF ( M>1 ) CALL SLACPY('Lower',M-1,M-1,Af(2,N-M+1),Lda,        &
     &                          Q(N-M+2,N-M+1),Lda)
      ELSE
         IF ( N>1 ) CALL SLACPY('Lower',N-1,N-1,Af(M-N+2,1),Lda,Q(2,1), &
     &                          Lda)
      ENDIF
!
!     Generate the n-by-n matrix Q
!
      SRNamt = 'SORGRQ'
      CALL SORGRQ(N,N,minmn,Q,Lda,Tau,Work,Lwork,info)
!
!     Copy R
!
      CALL SLASET('Full',M,N,ZERO,ZERO,R,Lda)
      IF ( M<=N ) THEN
         IF ( M>0 ) CALL SLACPY('Upper',M,M,Af(1,N-M+1),Lda,R(1,N-M+1), &
     &                          Lda)
      ELSE
         IF ( M>N .AND. N>0 ) CALL SLACPY('Full',M-N,N,Af,Lda,R,Lda)
         IF ( N>0 ) CALL SLACPY('Upper',N,N,Af(M-N+1,1),Lda,R(M-N+1,1), &
     &                          Lda)
      ENDIF
!
!     Compute R - A*Q'
!
      CALL SGEMM('No transpose','Transpose',M,N,N,-ONE,A,Lda,Q,Lda,ONE, &
     &           R,Lda)
!
!     Compute norm( R - Q'*A ) / ( N * norm(A) * EPS ) .
!
      anorm = SLANGE('1',M,N,A,Lda,Rwork)
      resid = SLANGE('1',M,N,R,Lda,Rwork)
      IF ( anorm>ZERO ) THEN
         Result(1) = ((resid/REAL(MAX(1,N)))/anorm)/eps
      ELSE
         Result(1) = ZERO
      ENDIF
!
!     Compute I - Q*Q'
!
      CALL SLASET('Full',N,N,ZERO,ONE,R,Lda)
      CALL SSYRK('Upper','No transpose',N,N,-ONE,Q,Lda,ONE,R,Lda)
!
!     Compute norm( I - Q*Q' ) / ( N * EPS ) .
!
      resid = SLANSY('1','Upper',N,R,Lda,Rwork)
!
      Result(2) = (resid/REAL(MAX(1,N)))/eps
!
!
!     End of srqt01
!
      END SUBROUTINE srqt01
