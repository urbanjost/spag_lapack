!*==cqlt01.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b CQLT01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CQLT01( M, N, A, AF, Q, L, LDA, TAU, WORK, LWORK,
!                          RWORK, RESULT )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       REAL               RESULT( * ), RWORK( * )
!       COMPLEX            A( LDA, * ), AF( LDA, * ), L( LDA, * ),
!      $                   Q( LDA, * ), TAU( * ), WORK( LWORK )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CQLT01 tests CGEQLF, which computes the QL factorization of an m-by-n
!> matrix A, and partially tests CUNGQL which forms the m-by-m
!> orthogonal matrix Q.
!>
!> CQLT01 compares L with Q'*A, and checks that Q is orthogonal.
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
!>          A is COMPLEX array, dimension (LDA,N)
!>          The m-by-n matrix A.
!> \endverbatim
!>
!> \param[out] AF
!> \verbatim
!>          AF is COMPLEX array, dimension (LDA,N)
!>          Details of the QL factorization of A, as returned by CGEQLF.
!>          See CGEQLF for further details.
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is COMPLEX array, dimension (LDA,M)
!>          The m-by-m orthogonal matrix Q.
!> \endverbatim
!>
!> \param[out] L
!> \verbatim
!>          L is COMPLEX array, dimension (LDA,max(M,N))
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
!>          TAU is COMPLEX array, dimension (min(M,N))
!>          The scalar factors of the elementary reflectors, as returned
!>          by CGEQLF.
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
!> \ingroup complex_lin
!
!  =====================================================================
      SUBROUTINE CQLT01(M,N,A,Af,Q,L,Lda,Tau,Work,Lwork,Rwork,Result)
      IMPLICIT NONE
!*--CQLT01129
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
      REAL Result(*) , Rwork(*)
      COMPLEX A(Lda,*) , Af(Lda,*) , L(Lda,*) , Q(Lda,*) , Tau(*) ,     &
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
      INTEGER info , minmn
      REAL anorm , eps , resid
!     ..
!     .. External Functions ..
      REAL CLANGE , CLANSY , SLAMCH
      EXTERNAL CLANGE , CLANSY , SLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL CGEMM , CGEQLF , CHERK , CLACPY , CLASET , CUNGQL
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC CMPLX , MAX , MIN , REAL
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
      CALL CLACPY('Full',M,N,A,Lda,Af,Lda)
!
!     Factorize the matrix A in the array AF.
!
      SRNamt = 'CGEQLF'
      CALL CGEQLF(M,N,Af,Lda,Tau,Work,Lwork,info)
!
!     Copy details of Q
!
      CALL CLASET('Full',M,M,ROGUE,ROGUE,Q,Lda)
      IF ( M>=N ) THEN
         IF ( N<M .AND. N>0 ) CALL CLACPY('Full',M-N,N,Af,Lda,Q(1,M-N+1)&
     &        ,Lda)
         IF ( N>1 ) CALL CLACPY('Upper',N-1,N-1,Af(M-N+1,2),Lda,        &
     &                          Q(M-N+1,M-N+2),Lda)
      ELSE
         IF ( M>1 ) CALL CLACPY('Upper',M-1,M-1,Af(1,N-M+2),Lda,Q(1,2), &
     &                          Lda)
      ENDIF
!
!     Generate the m-by-m matrix Q
!
      SRNamt = 'CUNGQL'
      CALL CUNGQL(M,M,minmn,Q,Lda,Tau,Work,Lwork,info)
!
!     Copy L
!
      CALL CLASET('Full',M,N,CMPLX(ZERO),CMPLX(ZERO),L,Lda)
      IF ( M>=N ) THEN
         IF ( N>0 ) CALL CLACPY('Lower',N,N,Af(M-N+1,1),Lda,L(M-N+1,1), &
     &                          Lda)
      ELSE
         IF ( N>M .AND. M>0 ) CALL CLACPY('Full',M,N-M,Af,Lda,L,Lda)
         IF ( M>0 ) CALL CLACPY('Lower',M,M,Af(1,N-M+1),Lda,L(1,N-M+1), &
     &                          Lda)
      ENDIF
!
!     Compute L - Q'*A
!
      CALL CGEMM('Conjugate transpose','No transpose',M,N,M,CMPLX(-ONE),&
     &           Q,Lda,A,Lda,CMPLX(ONE),L,Lda)
!
!     Compute norm( L - Q'*A ) / ( M * norm(A) * EPS ) .
!
      anorm = CLANGE('1',M,N,A,Lda,Rwork)
      resid = CLANGE('1',M,N,L,Lda,Rwork)
      IF ( anorm>ZERO ) THEN
         Result(1) = ((resid/REAL(MAX(1,M)))/anorm)/eps
      ELSE
         Result(1) = ZERO
      ENDIF
!
!     Compute I - Q'*Q
!
      CALL CLASET('Full',M,M,CMPLX(ZERO),CMPLX(ONE),L,Lda)
      CALL CHERK('Upper','Conjugate transpose',M,M,-ONE,Q,Lda,ONE,L,Lda)
!
!     Compute norm( I - Q'*Q ) / ( M * EPS ) .
!
      resid = CLANSY('1','Upper',M,L,Lda,Rwork)
!
      Result(2) = (resid/REAL(MAX(1,M)))/eps
!
!
!     End of CQLT01
!
      END SUBROUTINE CQLT01
