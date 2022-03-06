!*==cgsvts3.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b CGSVTS3
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGSVTS3( M, P, N, A, AF, LDA, B, BF, LDB, U, LDU, V,
!                           LDV, Q, LDQ, ALPHA, BETA, R, LDR, IWORK, WORK,
!                           LWORK, RWORK, RESULT )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LDB, LDQ, LDR, LDU, LDV, LWORK, M, N, P
!       ..
!       .. Array Arguments ..
!       INTEGER            IWORK( * )
!       REAL               ALPHA( * ), BETA( * ), RESULT( 6 ), RWORK( * )
!       COMPLEX            A( LDA, * ), AF( LDA, * ), B( LDB, * ),
!      $                   BF( LDB, * ), Q( LDQ, * ), R( LDR, * ),
!      $                   U( LDU, * ), V( LDV, * ), WORK( LWORK )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CGSVTS3 tests CGGSVD3, which computes the GSVD of an M-by-N matrix A
!> and a P-by-N matrix B:
!>              U'*A*Q = D1*R and V'*B*Q = D2*R.
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
!>          A is COMPLEX array, dimension (LDA,M)
!>          The M-by-N matrix A.
!> \endverbatim
!>
!> \param[out] AF
!> \verbatim
!>          AF is COMPLEX array, dimension (LDA,N)
!>          Details of the GSVD of A and B, as returned by CGGSVD3,
!>          see CGGSVD3 for further details.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the arrays A and AF.
!>          LDA >= max( 1,M ).
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is COMPLEX array, dimension (LDB,P)
!>          On entry, the P-by-N matrix B.
!> \endverbatim
!>
!> \param[out] BF
!> \verbatim
!>          BF is COMPLEX array, dimension (LDB,N)
!>          Details of the GSVD of A and B, as returned by CGGSVD3,
!>          see CGGSVD3 for further details.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the arrays B and BF.
!>          LDB >= max(1,P).
!> \endverbatim
!>
!> \param[out] U
!> \verbatim
!>          U is COMPLEX array, dimension(LDU,M)
!>          The M by M unitary matrix U.
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER
!>          The leading dimension of the array U. LDU >= max(1,M).
!> \endverbatim
!>
!> \param[out] V
!> \verbatim
!>          V is COMPLEX array, dimension(LDV,M)
!>          The P by P unitary matrix V.
!> \endverbatim
!>
!> \param[in] LDV
!> \verbatim
!>          LDV is INTEGER
!>          The leading dimension of the array V. LDV >= max(1,P).
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is COMPLEX array, dimension(LDQ,N)
!>          The N by N unitary matrix Q.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>          The leading dimension of the array Q. LDQ >= max(1,N).
!> \endverbatim
!>
!> \param[out] ALPHA
!> \verbatim
!>          ALPHA is REAL array, dimension (N)
!> \endverbatim
!>
!> \param[out] BETA
!> \verbatim
!>          BETA is REAL array, dimension (N)
!>
!>          The generalized singular value pairs of A and B, the
!>          ``diagonal'' matrices D1 and D2 are constructed from
!>          ALPHA and BETA, see subroutine CGGSVD3 for details.
!> \endverbatim
!>
!> \param[out] R
!> \verbatim
!>          R is COMPLEX array, dimension(LDQ,N)
!>          The upper triangular matrix R.
!> \endverbatim
!>
!> \param[in] LDR
!> \verbatim
!>          LDR is INTEGER
!>          The leading dimension of the array R. LDR >= max(1,N).
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (N)
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
!>          The dimension of the array WORK,
!>          LWORK >= max(M,P,N)*max(M,P,N).
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (max(M,P,N))
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is REAL array, dimension (6)
!>          The test ratios:
!>          RESULT(1) = norm( U'*A*Q - D1*R ) / ( MAX(M,N)*norm(A)*ULP)
!>          RESULT(2) = norm( V'*B*Q - D2*R ) / ( MAX(P,N)*norm(B)*ULP)
!>          RESULT(3) = norm( I - U'*U ) / ( M*ULP )
!>          RESULT(4) = norm( I - V'*V ) / ( P*ULP )
!>          RESULT(5) = norm( I - Q'*Q ) / ( N*ULP )
!>          RESULT(6) = 0        if ALPHA is in decreasing order;
!>                    = ULPINV   otherwise.
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
!> \date August 2015
!
!> \ingroup complex_eig
!
!  =====================================================================
      SUBROUTINE CGSVTS3(M,P,N,A,Af,Lda,B,Bf,Ldb,U,Ldu,V,Ldv,Q,Ldq,     &
     &                   Alpha,Beta,R,Ldr,Iwork,Work,Lwork,Rwork,Result)
      IMPLICIT NONE
!*--CGSVTS3212
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     August 2015
!
!     .. Scalar Arguments ..
      INTEGER Lda , Ldb , Ldq , Ldr , Ldu , Ldv , Lwork , M , N , P
!     ..
!     .. Array Arguments ..
      INTEGER Iwork(*)
      REAL Alpha(*) , Beta(*) , Result(6) , Rwork(*)
      COMPLEX A(Lda,*) , Af(Lda,*) , B(Ldb,*) , Bf(Ldb,*) , Q(Ldq,*) ,  &
     &        R(Ldr,*) , U(Ldu,*) , V(Ldv,*) , Work(Lwork)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
      COMPLEX CZERO , CONE
      PARAMETER (CZERO=(0.0E+0,0.0E+0),CONE=(1.0E+0,0.0E+0))
!     ..
!     .. Local Scalars ..
      INTEGER i , info , j , k , l
      REAL anorm , bnorm , resid , temp , ulp , ulpinv , unfl
!     ..
!     .. External Functions ..
      REAL CLANGE , CLANHE , SLAMCH
      EXTERNAL CLANGE , CLANHE , SLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL CGEMM , CGGSVD3 , CHERK , CLACPY , CLASET , SCOPY
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN , REAL
!     ..
!     .. Executable Statements ..
!
      ulp = SLAMCH('Precision')
      ulpinv = ONE/ulp
      unfl = SLAMCH('Safe minimum')
!
!     Copy the matrix A to the array AF.
!
      CALL CLACPY('Full',M,N,A,Lda,Af,Lda)
      CALL CLACPY('Full',P,N,B,Ldb,Bf,Ldb)
!
      anorm = MAX(CLANGE('1',M,N,A,Lda,Rwork),unfl)
      bnorm = MAX(CLANGE('1',P,N,B,Ldb,Rwork),unfl)
!
!     Factorize the matrices A and B in the arrays AF and BF.
!
      CALL CGGSVD3('U','V','Q',M,N,P,k,l,Af,Lda,Bf,Ldb,Alpha,Beta,U,Ldu,&
     &             V,Ldv,Q,Ldq,Work,Lwork,Rwork,Iwork,info)
!
!     Copy R
!
      DO i = 1 , MIN(k+l,M)
         DO j = i , k + l
            R(i,j) = Af(i,N-k-l+j)
         ENDDO
      ENDDO
!
      IF ( M-k<l ) THEN
         DO i = M + 1 , k + l
            DO j = i , k + l
               R(i,j) = Bf(i-k,N-k-l+j)
            ENDDO
         ENDDO
      ENDIF
!
!     Compute A:= U'*A*Q - D1*R
!
      CALL CGEMM('No transpose','No transpose',M,N,N,CONE,A,Lda,Q,Ldq,  &
     &           CZERO,Work,Lda)
!
      CALL CGEMM('Conjugate transpose','No transpose',M,N,M,CONE,U,Ldu, &
     &           Work,Lda,CZERO,A,Lda)
!
      DO i = 1 , k
         DO j = i , k + l
            A(i,N-k-l+j) = A(i,N-k-l+j) - R(i,j)
         ENDDO
      ENDDO
!
      DO i = k + 1 , MIN(k+l,M)
         DO j = i , k + l
            A(i,N-k-l+j) = A(i,N-k-l+j) - Alpha(i)*R(i,j)
         ENDDO
      ENDDO
!
!     Compute norm( U'*A*Q - D1*R ) / ( MAX(1,M,N)*norm(A)*ULP ) .
!
      resid = CLANGE('1',M,N,A,Lda,Rwork)
      IF ( anorm>ZERO ) THEN
         Result(1) = ((resid/REAL(MAX(1,M,N)))/anorm)/ulp
      ELSE
         Result(1) = ZERO
      ENDIF
!
!     Compute B := V'*B*Q - D2*R
!
      CALL CGEMM('No transpose','No transpose',P,N,N,CONE,B,Ldb,Q,Ldq,  &
     &           CZERO,Work,Ldb)
!
      CALL CGEMM('Conjugate transpose','No transpose',P,N,P,CONE,V,Ldv, &
     &           Work,Ldb,CZERO,B,Ldb)
!
      DO i = 1 , l
         DO j = i , l
            B(i,N-l+j) = B(i,N-l+j) - Beta(k+i)*R(k+i,k+j)
         ENDDO
      ENDDO
!
!     Compute norm( V'*B*Q - D2*R ) / ( MAX(P,N)*norm(B)*ULP ) .
!
      resid = CLANGE('1',P,N,B,Ldb,Rwork)
      IF ( bnorm>ZERO ) THEN
         Result(2) = ((resid/REAL(MAX(1,P,N)))/bnorm)/ulp
      ELSE
         Result(2) = ZERO
      ENDIF
!
!     Compute I - U'*U
!
      CALL CLASET('Full',M,M,CZERO,CONE,Work,Ldq)
      CALL CHERK('Upper','Conjugate transpose',M,M,-ONE,U,Ldu,ONE,Work, &
     &           Ldu)
!
!     Compute norm( I - U'*U ) / ( M * ULP ) .
!
      resid = CLANHE('1','Upper',M,Work,Ldu,Rwork)
      Result(3) = (resid/REAL(MAX(1,M)))/ulp
!
!     Compute I - V'*V
!
      CALL CLASET('Full',P,P,CZERO,CONE,Work,Ldv)
      CALL CHERK('Upper','Conjugate transpose',P,P,-ONE,V,Ldv,ONE,Work, &
     &           Ldv)
!
!     Compute norm( I - V'*V ) / ( P * ULP ) .
!
      resid = CLANHE('1','Upper',P,Work,Ldv,Rwork)
      Result(4) = (resid/REAL(MAX(1,P)))/ulp
!
!     Compute I - Q'*Q
!
      CALL CLASET('Full',N,N,CZERO,CONE,Work,Ldq)
      CALL CHERK('Upper','Conjugate transpose',N,N,-ONE,Q,Ldq,ONE,Work, &
     &           Ldq)
!
!     Compute norm( I - Q'*Q ) / ( N * ULP ) .
!
      resid = CLANHE('1','Upper',N,Work,Ldq,Rwork)
      Result(5) = (resid/REAL(MAX(1,N)))/ulp
!
!     Check sorting
!
      CALL SCOPY(N,Alpha,1,Rwork,1)
      DO i = k + 1 , MIN(k+l,M)
         j = Iwork(i)
         IF ( i/=j ) THEN
            temp = Rwork(i)
            Rwork(i) = Rwork(j)
            Rwork(j) = temp
         ENDIF
      ENDDO
!
      Result(6) = ZERO
      DO i = k + 1 , MIN(k+l,M) - 1
         IF ( Rwork(i)<Rwork(i+1) ) Result(6) = ulpinv
      ENDDO
!
!
!     End of CGSVTS3
!
      END SUBROUTINE CGSVTS3
