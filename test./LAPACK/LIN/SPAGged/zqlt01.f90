!*==zqlt01.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b ZQLT01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZQLT01( M, N, A, AF, Q, L, LDA, TAU, WORK, LWORK,
!                          RWORK, RESULT )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   RESULT( * ), RWORK( * )
!       COMPLEX*16         A( LDA, * ), AF( LDA, * ), L( LDA, * ),
!      $                   Q( LDA, * ), TAU( * ), WORK( LWORK )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZQLT01 tests ZGEQLF, which computes the QL factorization of an m-by-n
!> matrix A, and partially tests ZUNGQL which forms the m-by-m
!> orthogonal matrix Q.
!>
!> ZQLT01 compares L with Q'*A, and checks that Q is orthogonal.
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
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          The m-by-n matrix A.
!> \endverbatim
!>
!> \param[out] AF
!> \verbatim
!>          AF is COMPLEX*16 array, dimension (LDA,N)
!>          Details of the QL factorization of A, as returned by ZGEQLF.
!>          See ZGEQLF for further details.
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is COMPLEX*16 array, dimension (LDA,M)
!>          The m-by-m orthogonal matrix Q.
!> \endverbatim
!>
!> \param[out] L
!> \verbatim
!>          L is COMPLEX*16 array, dimension (LDA,max(M,N))
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
!>          TAU is COMPLEX*16 array, dimension (min(M,N))
!>          The scalar factors of the elementary reflectors, as returned
!>          by ZGEQLF.
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
!> \ingroup complex16_lin
!
!  =====================================================================
      SUBROUTINE ZQLT01(M,N,A,Af,Q,L,Lda,Tau,Work,Lwork,Rwork,Result)
      IMPLICIT NONE
!*--ZQLT01129
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
      DOUBLE PRECISION Result(*) , Rwork(*)
      COMPLEX*16 A(Lda,*) , Af(Lda,*) , L(Lda,*) , Q(Lda,*) , Tau(*) ,  &
     &           Work(Lwork)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
      COMPLEX*16 ROGUE
      PARAMETER (ROGUE=(-1.0D+10,-1.0D+10))
!     ..
!     .. Local Scalars ..
      INTEGER info , minmn
      DOUBLE PRECISION anorm , eps , resid
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH , ZLANGE , ZLANSY
      EXTERNAL DLAMCH , ZLANGE , ZLANSY
!     ..
!     .. External Subroutines ..
      EXTERNAL ZGEMM , ZGEQLF , ZHERK , ZLACPY , ZLASET , ZUNGQL
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DBLE , DCMPLX , MAX , MIN
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
      CALL ZLACPY('Full',M,N,A,Lda,Af,Lda)
!
!     Factorize the matrix A in the array AF.
!
      SRNamt = 'ZGEQLF'
      CALL ZGEQLF(M,N,Af,Lda,Tau,Work,Lwork,info)
!
!     Copy details of Q
!
      CALL ZLASET('Full',M,M,ROGUE,ROGUE,Q,Lda)
      IF ( M>=N ) THEN
         IF ( N<M .AND. N>0 ) CALL ZLACPY('Full',M-N,N,Af,Lda,Q(1,M-N+1)&
     &        ,Lda)
         IF ( N>1 ) CALL ZLACPY('Upper',N-1,N-1,Af(M-N+1,2),Lda,        &
     &                          Q(M-N+1,M-N+2),Lda)
      ELSE
         IF ( M>1 ) CALL ZLACPY('Upper',M-1,M-1,Af(1,N-M+2),Lda,Q(1,2), &
     &                          Lda)
      ENDIF
!
!     Generate the m-by-m matrix Q
!
      SRNamt = 'ZUNGQL'
      CALL ZUNGQL(M,M,minmn,Q,Lda,Tau,Work,Lwork,info)
!
!     Copy L
!
      CALL ZLASET('Full',M,N,DCMPLX(ZERO),DCMPLX(ZERO),L,Lda)
      IF ( M>=N ) THEN
         IF ( N>0 ) CALL ZLACPY('Lower',N,N,Af(M-N+1,1),Lda,L(M-N+1,1), &
     &                          Lda)
      ELSE
         IF ( N>M .AND. M>0 ) CALL ZLACPY('Full',M,N-M,Af,Lda,L,Lda)
         IF ( M>0 ) CALL ZLACPY('Lower',M,M,Af(1,N-M+1),Lda,L(1,N-M+1), &
     &                          Lda)
      ENDIF
!
!     Compute L - Q'*A
!
      CALL ZGEMM('Conjugate transpose','No transpose',M,N,M,DCMPLX(-ONE)&
     &           ,Q,Lda,A,Lda,DCMPLX(ONE),L,Lda)
!
!     Compute norm( L - Q'*A ) / ( M * norm(A) * EPS ) .
!
      anorm = ZLANGE('1',M,N,A,Lda,Rwork)
      resid = ZLANGE('1',M,N,L,Lda,Rwork)
      IF ( anorm>ZERO ) THEN
         Result(1) = ((resid/DBLE(MAX(1,M)))/anorm)/eps
      ELSE
         Result(1) = ZERO
      ENDIF
!
!     Compute I - Q'*Q
!
      CALL ZLASET('Full',M,M,DCMPLX(ZERO),DCMPLX(ONE),L,Lda)
      CALL ZHERK('Upper','Conjugate transpose',M,M,-ONE,Q,Lda,ONE,L,Lda)
!
!     Compute norm( I - Q'*Q ) / ( M * EPS ) .
!
      resid = ZLANSY('1','Upper',M,L,Lda,Rwork)
!
      Result(2) = (resid/DBLE(MAX(1,M)))/eps
!
!
!     End of ZQLT01
!
      END SUBROUTINE ZQLT01
