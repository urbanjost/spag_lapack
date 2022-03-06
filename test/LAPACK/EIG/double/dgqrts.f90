!*==dgqrts.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b DGQRTS
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGQRTS( N, M, P, A, AF, Q, R, LDA, TAUA, B, BF, Z, T,
!                          BWK, LDB, TAUB, WORK, LWORK, RWORK, RESULT )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LDB, LWORK, M, N, P
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), AF( LDA, * ), B( LDB, * ),
!      $                   BF( LDB, * ), BWK( LDB, * ), Q( LDA, * ),
!      $                   R( LDA, * ), RESULT( 4 ), RWORK( * ),
!      $                   T( LDB, * ), TAUA( * ), TAUB( * ),
!      $                   WORK( LWORK ), Z( LDB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGQRTS tests DGGQRF, which computes the GQR factorization of an
!> N-by-M matrix A and a N-by-P matrix B: A = Q*R and B = Q*T*Z.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of rows of the matrices A and B.  N >= 0.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of columns of the matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] P
!> \verbatim
!>          P is INTEGER
!>          The number of columns of the matrix B.  P >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,M)
!>          The N-by-M matrix A.
!> \endverbatim
!>
!> \param[out] AF
!> \verbatim
!>          AF is DOUBLE PRECISION array, dimension (LDA,N)
!>          Details of the GQR factorization of A and B, as returned
!>          by DGGQRF, see SGGQRF for further details.
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is DOUBLE PRECISION array, dimension (LDA,N)
!>          The M-by-M orthogonal matrix Q.
!> \endverbatim
!>
!> \param[out] R
!> \verbatim
!>          R is DOUBLE PRECISION array, dimension (LDA,MAX(M,N))
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the arrays A, AF, R and Q.
!>          LDA >= max(M,N).
!> \endverbatim
!>
!> \param[out] TAUA
!> \verbatim
!>          TAUA is DOUBLE PRECISION array, dimension (min(M,N))
!>          The scalar factors of the elementary reflectors, as returned
!>          by DGGQRF.
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (LDB,P)
!>          On entry, the N-by-P matrix A.
!> \endverbatim
!>
!> \param[out] BF
!> \verbatim
!>          BF is DOUBLE PRECISION array, dimension (LDB,N)
!>          Details of the GQR factorization of A and B, as returned
!>          by DGGQRF, see SGGQRF for further details.
!> \endverbatim
!>
!> \param[out] Z
!> \verbatim
!>          Z is DOUBLE PRECISION array, dimension (LDB,P)
!>          The P-by-P orthogonal matrix Z.
!> \endverbatim
!>
!> \param[out] T
!> \verbatim
!>          T is DOUBLE PRECISION array, dimension (LDB,max(P,N))
!> \endverbatim
!>
!> \param[out] BWK
!> \verbatim
!>          BWK is DOUBLE PRECISION array, dimension (LDB,N)
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the arrays B, BF, Z and T.
!>          LDB >= max(P,N).
!> \endverbatim
!>
!> \param[out] TAUB
!> \verbatim
!>          TAUB is DOUBLE PRECISION array, dimension (min(P,N))
!>          The scalar factors of the elementary reflectors, as returned
!>          by DGGRQF.
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
!>          The dimension of the array WORK, LWORK >= max(N,M,P)**2.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (max(N,M,P))
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is DOUBLE PRECISION array, dimension (4)
!>          The test ratios:
!>            RESULT(1) = norm( R - Q'*A ) / ( MAX(M,N)*norm(A)*ULP)
!>            RESULT(2) = norm( T*Z - Q'*B ) / (MAX(P,N)*norm(B)*ULP)
!>            RESULT(3) = norm( I - Q'*Q ) / ( M*ULP )
!>            RESULT(4) = norm( I - Z'*Z ) / ( P*ULP )
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
!> \ingroup double_eig
!
!  =====================================================================
      SUBROUTINE DGQRTS(N,M,P,A,Af,Q,R,Lda,Taua,B,Bf,Z,T,Bwk,Ldb,Taub,  &
     &                  Work,Lwork,Rwork,Result)
      IMPLICIT NONE
!*--DGQRTS180
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Lda , Ldb , Lwork , M , N , P
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(Lda,*) , Af(Lda,*) , B(Ldb,*) , Bf(Ldb,*) ,    &
     &                 Bwk(Ldb,*) , Q(Lda,*) , R(Lda,*) , Result(4) ,   &
     &                 Rwork(*) , T(Ldb,*) , Taua(*) , Taub(*) ,        &
     &                 Work(Lwork) , Z(Ldb,*)
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
      INTEGER info
      DOUBLE PRECISION anorm , bnorm , resid , ulp , unfl
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH , DLANGE , DLANSY
      EXTERNAL DLAMCH , DLANGE , DLANSY
!     ..
!     .. External Subroutines ..
      EXTERNAL DGEMM , DGGQRF , DLACPY , DLASET , DORGQR , DORGRQ ,     &
     &         DSYRK
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DBLE , MAX , MIN
!     ..
!     .. Executable Statements ..
!
      ulp = DLAMCH('Precision')
      unfl = DLAMCH('Safe minimum')
!
!     Copy the matrix A to the array AF.
!
      CALL DLACPY('Full',N,M,A,Lda,Af,Lda)
      CALL DLACPY('Full',N,P,B,Ldb,Bf,Ldb)
!
      anorm = MAX(DLANGE('1',N,M,A,Lda,Rwork),unfl)
      bnorm = MAX(DLANGE('1',N,P,B,Ldb,Rwork),unfl)
!
!     Factorize the matrices A and B in the arrays AF and BF.
!
      CALL DGGQRF(N,M,P,Af,Lda,Taua,Bf,Ldb,Taub,Work,Lwork,info)
!
!     Generate the N-by-N matrix Q
!
      CALL DLASET('Full',N,N,ROGUE,ROGUE,Q,Lda)
      CALL DLACPY('Lower',N-1,M,Af(2,1),Lda,Q(2,1),Lda)
      CALL DORGQR(N,N,MIN(N,M),Q,Lda,Taua,Work,Lwork,info)
!
!     Generate the P-by-P matrix Z
!
      CALL DLASET('Full',P,P,ROGUE,ROGUE,Z,Ldb)
      IF ( N<=P ) THEN
         IF ( N>0 .AND. N<P ) CALL DLACPY('Full',N,P-N,Bf,Ldb,Z(P-N+1,1)&
     &        ,Ldb)
         IF ( N>1 ) CALL DLACPY('Lower',N-1,N-1,Bf(2,P-N+1),Ldb,        &
     &                          Z(P-N+2,P-N+1),Ldb)
      ELSE
         IF ( P>1 ) CALL DLACPY('Lower',P-1,P-1,Bf(N-P+2,1),Ldb,Z(2,1), &
     &                          Ldb)
      ENDIF
      CALL DORGRQ(P,P,MIN(N,P),Z,Ldb,Taub,Work,Lwork,info)
!
!     Copy R
!
      CALL DLASET('Full',N,M,ZERO,ZERO,R,Lda)
      CALL DLACPY('Upper',N,M,Af,Lda,R,Lda)
!
!     Copy T
!
      CALL DLASET('Full',N,P,ZERO,ZERO,T,Ldb)
      IF ( N<=P ) THEN
         CALL DLACPY('Upper',N,N,Bf(1,P-N+1),Ldb,T(1,P-N+1),Ldb)
      ELSE
         CALL DLACPY('Full',N-P,P,Bf,Ldb,T,Ldb)
         CALL DLACPY('Upper',P,P,Bf(N-P+1,1),Ldb,T(N-P+1,1),Ldb)
      ENDIF
!
!     Compute R - Q'*A
!
      CALL DGEMM('Transpose','No transpose',N,M,N,-ONE,Q,Lda,A,Lda,ONE, &
     &           R,Lda)
!
!     Compute norm( R - Q'*A ) / ( MAX(M,N)*norm(A)*ULP ) .
!
      resid = DLANGE('1',N,M,R,Lda,Rwork)
      IF ( anorm>ZERO ) THEN
         Result(1) = ((resid/DBLE(MAX(1,M,N)))/anorm)/ulp
      ELSE
         Result(1) = ZERO
      ENDIF
!
!     Compute T*Z - Q'*B
!
      CALL DGEMM('No Transpose','No transpose',N,P,P,ONE,T,Ldb,Z,Ldb,   &
     &           ZERO,Bwk,Ldb)
      CALL DGEMM('Transpose','No transpose',N,P,N,-ONE,Q,Lda,B,Ldb,ONE, &
     &           Bwk,Ldb)
!
!     Compute norm( T*Z - Q'*B ) / ( MAX(P,N)*norm(A)*ULP ) .
!
      resid = DLANGE('1',N,P,Bwk,Ldb,Rwork)
      IF ( bnorm>ZERO ) THEN
         Result(2) = ((resid/DBLE(MAX(1,P,N)))/bnorm)/ulp
      ELSE
         Result(2) = ZERO
      ENDIF
!
!     Compute I - Q'*Q
!
      CALL DLASET('Full',N,N,ZERO,ONE,R,Lda)
      CALL DSYRK('Upper','Transpose',N,N,-ONE,Q,Lda,ONE,R,Lda)
!
!     Compute norm( I - Q'*Q ) / ( N * ULP ) .
!
      resid = DLANSY('1','Upper',N,R,Lda,Rwork)
      Result(3) = (resid/DBLE(MAX(1,N)))/ulp
!
!     Compute I - Z'*Z
!
      CALL DLASET('Full',P,P,ZERO,ONE,T,Ldb)
      CALL DSYRK('Upper','Transpose',P,P,-ONE,Z,Ldb,ONE,T,Ldb)
!
!     Compute norm( I - Z'*Z ) / ( P*ULP ) .
!
      resid = DLANSY('1','Upper',P,T,Ldb,Rwork)
      Result(4) = (resid/DBLE(MAX(1,P)))/ulp
!
!
!     End of DGQRTS
!
      END SUBROUTINE DGQRTS
