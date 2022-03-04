!*==zhst01.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b ZHST01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZHST01( N, ILO, IHI, A, LDA, H, LDH, Q, LDQ, WORK,
!                          LWORK, RWORK, RESULT )
!
!       .. Scalar Arguments ..
!       INTEGER            IHI, ILO, LDA, LDH, LDQ, LWORK, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   RESULT( 2 ), RWORK( * )
!       COMPLEX*16         A( LDA, * ), H( LDH, * ), Q( LDQ, * ),
!      $                   WORK( LWORK )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZHST01 tests the reduction of a general matrix A to upper Hessenberg
!> form:  A = Q*H*Q'.  Two test ratios are computed;
!>
!> RESULT(1) = norm( A - Q*H*Q' ) / ( norm(A) * N * EPS )
!> RESULT(2) = norm( I - Q'*Q ) / ( N * EPS )
!>
!> The matrix Q is assumed to be given explicitly as it would be
!> following ZGEHRD + ZUNGHR.
!>
!> In this version, ILO and IHI are not used, but they could be used
!> to save some work if this is desired.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] ILO
!> \verbatim
!>          ILO is INTEGER
!> \endverbatim
!>
!> \param[in] IHI
!> \verbatim
!>          IHI is INTEGER
!>
!>          A is assumed to be upper triangular in rows and columns
!>          1:ILO-1 and IHI+1:N, so Q differs from the identity only in
!>          rows and columns ILO+1:IHI.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          The original n by n matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] H
!> \verbatim
!>          H is COMPLEX*16 array, dimension (LDH,N)
!>          The upper Hessenberg matrix H from the reduction A = Q*H*Q'
!>          as computed by ZGEHRD.  H is assumed to be zero below the
!>          first subdiagonal.
!> \endverbatim
!>
!> \param[in] LDH
!> \verbatim
!>          LDH is INTEGER
!>          The leading dimension of the array H.  LDH >= max(1,N).
!> \endverbatim
!>
!> \param[in] Q
!> \verbatim
!>          Q is COMPLEX*16 array, dimension (LDQ,N)
!>          The orthogonal matrix Q from the reduction A = Q*H*Q' as
!>          computed by ZGEHRD + ZUNGHR.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>          The leading dimension of the array Q.  LDQ >= max(1,N).
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
!>          The length of the array WORK.  LWORK >= 2*N*N.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is DOUBLE PRECISION array, dimension (2)
!>          RESULT(1) = norm( A - Q*H*Q' ) / ( norm(A) * N * EPS )
!>          RESULT(2) = norm( I - Q'*Q ) / ( N * EPS )
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
!> \ingroup complex16_eig
!
!  =====================================================================
      SUBROUTINE ZHST01(N,Ilo,Ihi,A,Lda,H,Ldh,Q,Ldq,Work,Lwork,Rwork,   &
     &                  Result)
      IMPLICIT NONE
!*--ZHST01144
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Ihi , Ilo , Lda , Ldh , Ldq , Lwork , N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION Result(2) , Rwork(*)
      COMPLEX*16 A(Lda,*) , H(Ldh,*) , Q(Ldq,*) , Work(Lwork)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ONE , ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
      INTEGER ldwork
      DOUBLE PRECISION anorm , eps , ovfl , smlnum , unfl , wnorm
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH , ZLANGE
      EXTERNAL DLAMCH , ZLANGE
!     ..
!     .. External Subroutines ..
      EXTERNAL DLABAD , ZGEMM , ZLACPY , ZUNT01
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DCMPLX , MAX , MIN
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF ( N<=0 ) THEN
         Result(1) = ZERO
         Result(2) = ZERO
         RETURN
      ENDIF
!
      unfl = DLAMCH('Safe minimum')
      eps = DLAMCH('Precision')
      ovfl = ONE/unfl
      CALL DLABAD(unfl,ovfl)
      smlnum = unfl*N/eps
!
!     Test 1:  Compute norm( A - Q*H*Q' ) / ( norm(A) * N * EPS )
!
!     Copy A to WORK
!
      ldwork = MAX(1,N)
      CALL ZLACPY(' ',N,N,A,Lda,Work,ldwork)
!
!     Compute Q*H
!
      CALL ZGEMM('No transpose','No transpose',N,N,N,DCMPLX(ONE),Q,Ldq, &
     &           H,Ldh,DCMPLX(ZERO),Work(ldwork*N+1),ldwork)
!
!     Compute A - Q*H*Q'
!
      CALL ZGEMM('No transpose','Conjugate transpose',N,N,N,DCMPLX(-ONE)&
     &           ,Work(ldwork*N+1),ldwork,Q,Ldq,DCMPLX(ONE),Work,ldwork)
!
      anorm = MAX(ZLANGE('1',N,N,A,Lda,Rwork),unfl)
      wnorm = ZLANGE('1',N,N,Work,ldwork,Rwork)
!
!     Note that RESULT(1) cannot overflow and is bounded by 1/(N*EPS)
!
      Result(1) = MIN(wnorm,anorm)/MAX(smlnum,anorm*eps)/N
!
!     Test 2:  Compute norm( I - Q'*Q ) / ( N * EPS )
!
      CALL ZUNT01('Columns',N,N,Q,Ldq,Work,Lwork,Rwork,Result(2))
!
!
!     End of ZHST01
!
      END SUBROUTINE ZHST01
