!*==dqlt01.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b DQLT01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DQLT01( M, N, A, AF, Q, L, LDA, TAU, WORK, LWORK,
!                          RWORK, RESULT )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), AF( LDA, * ), L( LDA, * ),
!      $                   Q( LDA, * ), RESULT( * ), RWORK( * ), TAU( * ),
!      $                   WORK( LWORK )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DQLT01 tests DGEQLF, which computes the QL factorization of an m-by-n
!> matrix A, and partially tests DORGQL which forms the m-by-m
!> orthogonal matrix Q.
!>
!> DQLT01 compares L with Q'*A, and checks that Q is orthogonal.
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
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          The m-by-n matrix A.
!> \endverbatim
!>
!> \param[out] AF
!> \verbatim
!>          AF is DOUBLE PRECISION array, dimension (LDA,N)
!>          Details of the QL factorization of A, as returned by DGEQLF.
!>          See DGEQLF for further details.
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is DOUBLE PRECISION array, dimension (LDA,M)
!>          The m-by-m orthogonal matrix Q.
!> \endverbatim
!>
!> \param[out] L
!> \verbatim
!>          L is DOUBLE PRECISION array, dimension (LDA,max(M,N))
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the arrays A, AF, Q and R.
!>          LDA >= max(M,N).
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is DOUBLE PRECISION array, dimension (min(M,N))
!>          The scalar factors of the elementary reflectors, as returned
!>          by DGEQLF.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (LWORK)
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
!>          RWORK is DOUBLE PRECISION array, dimension (M)
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is DOUBLE PRECISION array, dimension (2)
!>          The test ratios:
!>          RESULT(1) = norm( L - Q'*A ) / ( M * norm(A) * EPS )
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
!> \ingroup double_lin
!
!  =====================================================================
      SUBROUTINE DQLT01(M,N,A,Af,Q,L,Lda,Tau,Work,Lwork,Rwork,Result)
      IMPLICIT NONE
!*--DQLT01129
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
      DOUBLE PRECISION A(Lda,*) , Af(Lda,*) , L(Lda,*) , Q(Lda,*) ,     &
     &                 Result(*) , Rwork(*) , Tau(*) , Work(Lwork)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
      DOUBLE PRECISION ROGUE
      PARAMETER (ROGUE=-1.0D+10)
!     ..
!     .. Local Scalars ..
      INTEGER info , minmn
      DOUBLE PRECISION anorm , eps , resid
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH , DLANGE , DLANSY
      EXTERNAL DLAMCH , DLANGE , DLANSY
!     ..
!     .. External Subroutines ..
      EXTERNAL DGEMM , DGEQLF , DLACPY , DLASET , DORGQL , DSYRK
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DBLE , MAX , MIN
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
      eps = DLAMCH('Epsilon')
!
!     Copy the matrix A to the array AF.
!
      CALL DLACPY('Full',M,N,A,Lda,Af,Lda)
!
!     Factorize the matrix A in the array AF.
!
      SRNamt = 'DGEQLF'
      CALL DGEQLF(M,N,Af,Lda,Tau,Work,Lwork,info)
!
!     Copy details of Q
!
      CALL DLASET('Full',M,M,ROGUE,ROGUE,Q,Lda)
      IF ( M>=N ) THEN
         IF ( N<M .AND. N>0 ) CALL DLACPY('Full',M-N,N,Af,Lda,Q(1,M-N+1)&
     &        ,Lda)
         IF ( N>1 ) CALL DLACPY('Upper',N-1,N-1,Af(M-N+1,2),Lda,        &
     &                          Q(M-N+1,M-N+2),Lda)
      ELSE
         IF ( M>1 ) CALL DLACPY('Upper',M-1,M-1,Af(1,N-M+2),Lda,Q(1,2), &
     &                          Lda)
      ENDIF
!
!     Generate the m-by-m matrix Q
!
      SRNamt = 'DORGQL'
      CALL DORGQL(M,M,minmn,Q,Lda,Tau,Work,Lwork,info)
!
!     Copy L
!
      CALL DLASET('Full',M,N,ZERO,ZERO,L,Lda)
      IF ( M>=N ) THEN
         IF ( N>0 ) CALL DLACPY('Lower',N,N,Af(M-N+1,1),Lda,L(M-N+1,1), &
     &                          Lda)
      ELSE
         IF ( N>M .AND. M>0 ) CALL DLACPY('Full',M,N-M,Af,Lda,L,Lda)
         IF ( M>0 ) CALL DLACPY('Lower',M,M,Af(1,N-M+1),Lda,L(1,N-M+1), &
     &                          Lda)
      ENDIF
!
!     Compute L - Q'*A
!
      CALL DGEMM('Transpose','No transpose',M,N,M,-ONE,Q,Lda,A,Lda,ONE, &
     &           L,Lda)
!
!     Compute norm( L - Q'*A ) / ( M * norm(A) * EPS ) .
!
      anorm = DLANGE('1',M,N,A,Lda,Rwork)
      resid = DLANGE('1',M,N,L,Lda,Rwork)
      IF ( anorm>ZERO ) THEN
         Result(1) = ((resid/DBLE(MAX(1,M)))/anorm)/eps
      ELSE
         Result(1) = ZERO
      ENDIF
!
!     Compute I - Q'*Q
!
      CALL DLASET('Full',M,M,ZERO,ONE,L,Lda)
      CALL DSYRK('Upper','Transpose',M,M,-ONE,Q,Lda,ONE,L,Lda)
!
!     Compute norm( I - Q'*Q ) / ( M * EPS ) .
!
      resid = DLANSY('1','Upper',M,L,Lda,Rwork)
!
      Result(2) = (resid/DBLE(MAX(1,M)))/eps
!
!
!     End of DQLT01
!
      END SUBROUTINE DQLT01
