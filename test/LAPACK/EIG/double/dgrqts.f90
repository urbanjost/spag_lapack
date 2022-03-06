!*==dgrqts.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b DGRQTS
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGRQTS( M, P, N, A, AF, Q, R, LDA, TAUA, B, BF, Z, T,
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
!> DGRQTS tests DGGRQF, which computes the GRQ factorization of an
!> M-by-N matrix A and a P-by-N matrix B: A = R*Q and B = Z*T*Q.
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
!> \param[in] P
!> \verbatim
!>          P is INTEGER
!>          The number of rows of the matrix B.  P >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrices A and B.  N >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          The M-by-N matrix A.
!> \endverbatim
!>
!> \param[out] AF
!> \verbatim
!>          AF is DOUBLE PRECISION array, dimension (LDA,N)
!>          Details of the GRQ factorization of A and B, as returned
!>          by DGGRQF, see SGGRQF for further details.
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is DOUBLE PRECISION array, dimension (LDA,N)
!>          The N-by-N orthogonal matrix Q.
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
!>          by DGGQRC.
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (LDB,N)
!>          On entry, the P-by-N matrix A.
!> \endverbatim
!>
!> \param[out] BF
!> \verbatim
!>          BF is DOUBLE PRECISION array, dimension (LDB,N)
!>          Details of the GQR factorization of A and B, as returned
!>          by DGGRQF, see SGGRQF for further details.
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
!>          The dimension of the array WORK, LWORK >= max(M,P,N)**2.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (M)
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is DOUBLE PRECISION array, dimension (4)
!>          The test ratios:
!>            RESULT(1) = norm( R - A*Q' ) / ( MAX(M,N)*norm(A)*ULP)
!>            RESULT(2) = norm( T*Q - Z'*B ) / (MAX(P,N)*norm(B)*ULP)
!>            RESULT(3) = norm( I - Q'*Q ) / ( N*ULP )
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
      SUBROUTINE DGRQTS(M,P,N,A,Af,Q,R,Lda,Taua,B,Bf,Z,T,Bwk,Ldb,Taub,  &
     &                  Work,Lwork,Rwork,Result)
      IMPLICIT NONE
!*--DGRQTS180
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
      EXTERNAL DGEMM , DGGRQF , DLACPY , DLASET , DORGQR , DORGRQ ,     &
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
      CALL DLACPY('Full',M,N,A,Lda,Af,Lda)
      CALL DLACPY('Full',P,N,B,Ldb,Bf,Ldb)
!
      anorm = MAX(DLANGE('1',M,N,A,Lda,Rwork),unfl)
      bnorm = MAX(DLANGE('1',P,N,B,Ldb,Rwork),unfl)
!
!     Factorize the matrices A and B in the arrays AF and BF.
!
      CALL DGGRQF(M,P,N,Af,Lda,Taua,Bf,Ldb,Taub,Work,Lwork,info)
!
!     Generate the N-by-N matrix Q
!
      CALL DLASET('Full',N,N,ROGUE,ROGUE,Q,Lda)
      IF ( M<=N ) THEN
         IF ( M>0 .AND. M<N ) CALL DLACPY('Full',M,N-M,Af,Lda,Q(N-M+1,1)&
     &        ,Lda)
         IF ( M>1 ) CALL DLACPY('Lower',M-1,M-1,Af(2,N-M+1),Lda,        &
     &                          Q(N-M+2,N-M+1),Lda)
      ELSE
         IF ( N>1 ) CALL DLACPY('Lower',N-1,N-1,Af(M-N+2,1),Lda,Q(2,1), &
     &                          Lda)
      ENDIF
      CALL DORGRQ(N,N,MIN(M,N),Q,Lda,Taua,Work,Lwork,info)
!
!     Generate the P-by-P matrix Z
!
      CALL DLASET('Full',P,P,ROGUE,ROGUE,Z,Ldb)
      IF ( P>1 ) CALL DLACPY('Lower',P-1,N,Bf(2,1),Ldb,Z(2,1),Ldb)
      CALL DORGQR(P,P,MIN(P,N),Z,Ldb,Taub,Work,Lwork,info)
!
!     Copy R
!
      CALL DLASET('Full',M,N,ZERO,ZERO,R,Lda)
      IF ( M<=N ) THEN
         CALL DLACPY('Upper',M,M,Af(1,N-M+1),Lda,R(1,N-M+1),Lda)
      ELSE
         CALL DLACPY('Full',M-N,N,Af,Lda,R,Lda)
         CALL DLACPY('Upper',N,N,Af(M-N+1,1),Lda,R(M-N+1,1),Lda)
      ENDIF
!
!     Copy T
!
      CALL DLASET('Full',P,N,ZERO,ZERO,T,Ldb)
      CALL DLACPY('Upper',P,N,Bf,Ldb,T,Ldb)
!
!     Compute R - A*Q'
!
      CALL DGEMM('No transpose','Transpose',M,N,N,-ONE,A,Lda,Q,Lda,ONE, &
     &           R,Lda)
!
!     Compute norm( R - A*Q' ) / ( MAX(M,N)*norm(A)*ULP ) .
!
      resid = DLANGE('1',M,N,R,Lda,Rwork)
      IF ( anorm>ZERO ) THEN
         Result(1) = ((resid/DBLE(MAX(1,M,N)))/anorm)/ulp
      ELSE
         Result(1) = ZERO
      ENDIF
!
!     Compute T*Q - Z'*B
!
      CALL DGEMM('Transpose','No transpose',P,N,P,ONE,Z,Ldb,B,Ldb,ZERO, &
     &           Bwk,Ldb)
      CALL DGEMM('No transpose','No transpose',P,N,N,ONE,T,Ldb,Q,Lda,   &
     &           -ONE,Bwk,Ldb)
!
!     Compute norm( T*Q - Z'*B ) / ( MAX(P,N)*norm(A)*ULP ) .
!
      resid = DLANGE('1',P,N,Bwk,Ldb,Rwork)
      IF ( bnorm>ZERO ) THEN
         Result(2) = ((resid/DBLE(MAX(1,P,M)))/bnorm)/ulp
      ELSE
         Result(2) = ZERO
      ENDIF
!
!     Compute I - Q*Q'
!
      CALL DLASET('Full',N,N,ZERO,ONE,R,Lda)
      CALL DSYRK('Upper','No Transpose',N,N,-ONE,Q,Lda,ONE,R,Lda)
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
!     End of DGRQTS
!
      END SUBROUTINE DGRQTS
