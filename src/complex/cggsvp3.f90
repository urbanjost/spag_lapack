!*==cggsvp3.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b CGGSVP3
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CGGSVP3 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cggsvp3.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cggsvp3.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cggsvp3.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGGSVP3( JOBU, JOBV, JOBQ, M, P, N, A, LDA, B, LDB,
!                           TOLA, TOLB, K, L, U, LDU, V, LDV, Q, LDQ,
!                           IWORK, RWORK, TAU, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBQ, JOBU, JOBV
!       INTEGER            INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N, P, LWORK
!       REAL               TOLA, TOLB
!       ..
!       .. Array Arguments ..
!       INTEGER            IWORK( * )
!       REAL               RWORK( * )
!       COMPLEX            A( LDA, * ), B( LDB, * ), Q( LDQ, * ),
!      $                   TAU( * ), U( LDU, * ), V( LDV, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CGGSVP3 computes unitary matrices U, V and Q such that
!>
!>                    N-K-L  K    L
!>  U**H*A*Q =     K ( 0    A12  A13 )  if M-K-L >= 0;
!>                 L ( 0     0   A23 )
!>             M-K-L ( 0     0    0  )
!>
!>                  N-K-L  K    L
!>         =     K ( 0    A12  A13 )  if M-K-L < 0;
!>             M-K ( 0     0   A23 )
!>
!>                  N-K-L  K    L
!>  V**H*B*Q =   L ( 0     0   B13 )
!>             P-L ( 0     0    0  )
!>
!> where the K-by-K matrix A12 and L-by-L matrix B13 are nonsingular
!> upper triangular; A23 is L-by-L upper triangular if M-K-L >= 0,
!> otherwise A23 is (M-K)-by-L upper trapezoidal.  K+L = the effective
!> numerical rank of the (M+P)-by-N matrix (A**H,B**H)**H.
!>
!> This decomposition is the preprocessing step for computing the
!> Generalized Singular Value Decomposition (GSVD), see subroutine
!> CGGSVD3.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOBU
!> \verbatim
!>          JOBU is CHARACTER*1
!>          = 'U':  Unitary matrix U is computed;
!>          = 'N':  U is not computed.
!> \endverbatim
!>
!> \param[in] JOBV
!> \verbatim
!>          JOBV is CHARACTER*1
!>          = 'V':  Unitary matrix V is computed;
!>          = 'N':  V is not computed.
!> \endverbatim
!>
!> \param[in] JOBQ
!> \verbatim
!>          JOBQ is CHARACTER*1
!>          = 'Q':  Unitary matrix Q is computed;
!>          = 'N':  Q is not computed.
!> \endverbatim
!>
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
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          On entry, the M-by-N matrix A.
!>          On exit, A contains the triangular (or trapezoidal) matrix
!>          described in the Purpose section.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A. LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is COMPLEX array, dimension (LDB,N)
!>          On entry, the P-by-N matrix B.
!>          On exit, B contains the triangular matrix described in
!>          the Purpose section.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B. LDB >= max(1,P).
!> \endverbatim
!>
!> \param[in] TOLA
!> \verbatim
!>          TOLA is REAL
!> \endverbatim
!>
!> \param[in] TOLB
!> \verbatim
!>          TOLB is REAL
!>
!>          TOLA and TOLB are the thresholds to determine the effective
!>          numerical rank of matrix B and a subblock of A. Generally,
!>          they are set to
!>             TOLA = MAX(M,N)*norm(A)*MACHEPS,
!>             TOLB = MAX(P,N)*norm(B)*MACHEPS.
!>          The size of TOLA and TOLB may affect the size of backward
!>          errors of the decomposition.
!> \endverbatim
!>
!> \param[out] K
!> \verbatim
!>          K is INTEGER
!> \endverbatim
!>
!> \param[out] L
!> \verbatim
!>          L is INTEGER
!>
!>          On exit, K and L specify the dimension of the subblocks
!>          described in Purpose section.
!>          K + L = effective numerical rank of (A**H,B**H)**H.
!> \endverbatim
!>
!> \param[out] U
!> \verbatim
!>          U is COMPLEX array, dimension (LDU,M)
!>          If JOBU = 'U', U contains the unitary matrix U.
!>          If JOBU = 'N', U is not referenced.
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER
!>          The leading dimension of the array U. LDU >= max(1,M) if
!>          JOBU = 'U'; LDU >= 1 otherwise.
!> \endverbatim
!>
!> \param[out] V
!> \verbatim
!>          V is COMPLEX array, dimension (LDV,P)
!>          If JOBV = 'V', V contains the unitary matrix V.
!>          If JOBV = 'N', V is not referenced.
!> \endverbatim
!>
!> \param[in] LDV
!> \verbatim
!>          LDV is INTEGER
!>          The leading dimension of the array V. LDV >= max(1,P) if
!>          JOBV = 'V'; LDV >= 1 otherwise.
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is COMPLEX array, dimension (LDQ,N)
!>          If JOBQ = 'Q', Q contains the unitary matrix Q.
!>          If JOBQ = 'N', Q is not referenced.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>          The leading dimension of the array Q. LDQ >= max(1,N) if
!>          JOBQ = 'Q'; LDQ >= 1 otherwise.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (N)
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (2*N)
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is COMPLEX array, dimension (N)
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
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
!> \ingroup complexOTHERcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The subroutine uses LAPACK subroutine CGEQP3 for the QR factorization
!>  with column pivoting to detect the effective numerical rank of the
!>  a matrix. It may be replaced by a better rank determination strategy.
!>
!>  CGGSVP3 replaces the deprecated subroutine CGGSVP.
!>
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE CGGSVP3(Jobu,Jobv,Jobq,M,P,N,A,Lda,B,Ldb,Tola,Tolb,K,L,&
     &                   U,Ldu,V,Ldv,Q,Ldq,Iwork,Rwork,Tau,Work,Lwork,  &
     &                   Info)
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     August 2015
!
      IMPLICIT NONE
!*--CGGSVP3288
!
!     .. Scalar Arguments ..
      CHARACTER Jobq , Jobu , Jobv
      INTEGER Info , K , L , Lda , Ldb , Ldq , Ldu , Ldv , M , N , P ,  &
     &        Lwork
      REAL Tola , Tolb
!     ..
!     .. Array Arguments ..
      INTEGER Iwork(*)
      REAL Rwork(*)
      COMPLEX A(Lda,*) , B(Ldb,*) , Q(Ldq,*) , Tau(*) , U(Ldu,*) ,      &
     &        V(Ldv,*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX CZERO , CONE
      PARAMETER (CZERO=(0.0E+0,0.0E+0),CONE=(1.0E+0,0.0E+0))
!     ..
!     .. Local Scalars ..
      LOGICAL forwrd , wantq , wantu , wantv , lquery
      INTEGER i , j , lwkopt
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL CGEQP3 , CGEQR2 , CGERQ2 , CLACPY , CLAPMT , CLASET ,    &
     &         CUNG2R , CUNM2R , CUNMR2 , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , AIMAG , MAX , MIN , REAL
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      wantu = LSAME(Jobu,'U')
      wantv = LSAME(Jobv,'V')
      wantq = LSAME(Jobq,'Q')
      forwrd = .TRUE.
      lquery = (Lwork==-1)
      lwkopt = 1
!
!     Test the input arguments
!
      Info = 0
      IF ( .NOT.(wantu .OR. LSAME(Jobu,'N')) ) THEN
         Info = -1
      ELSEIF ( .NOT.(wantv .OR. LSAME(Jobv,'N')) ) THEN
         Info = -2
      ELSEIF ( .NOT.(wantq .OR. LSAME(Jobq,'N')) ) THEN
         Info = -3
      ELSEIF ( M<0 ) THEN
         Info = -4
      ELSEIF ( P<0 ) THEN
         Info = -5
      ELSEIF ( N<0 ) THEN
         Info = -6
      ELSEIF ( Lda<MAX(1,M) ) THEN
         Info = -8
      ELSEIF ( Ldb<MAX(1,P) ) THEN
         Info = -10
      ELSEIF ( Ldu<1 .OR. (wantu .AND. Ldu<M) ) THEN
         Info = -16
      ELSEIF ( Ldv<1 .OR. (wantv .AND. Ldv<P) ) THEN
         Info = -18
      ELSEIF ( Ldq<1 .OR. (wantq .AND. Ldq<N) ) THEN
         Info = -20
      ELSEIF ( Lwork<1 .AND. .NOT.lquery ) THEN
         Info = -24
      ENDIF
!
!     Compute workspace
!
      IF ( Info==0 ) THEN
         CALL CGEQP3(P,N,B,Ldb,Iwork,Tau,Work,-1,Rwork,Info)
         lwkopt = INT(Work(1))
         IF ( wantv ) lwkopt = MAX(lwkopt,P)
         lwkopt = MAX(lwkopt,MIN(N,P))
         lwkopt = MAX(lwkopt,M)
         IF ( wantq ) lwkopt = MAX(lwkopt,N)
         CALL CGEQP3(M,N,A,Lda,Iwork,Tau,Work,-1,Rwork,Info)
         lwkopt = MAX(lwkopt,INT(Work(1)))
         lwkopt = MAX(1,lwkopt)
         Work(1) = CMPLX(lwkopt)
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('CGGSVP3',-Info)
         RETURN
      ENDIF
      IF ( lquery ) RETURN
!
!     QR with column pivoting of B: B*P = V*( S11 S12 )
!                                           (  0   0  )
!
      DO i = 1 , N
         Iwork(i) = 0
      ENDDO
      CALL CGEQP3(P,N,B,Ldb,Iwork,Tau,Work,Lwork,Rwork,Info)
!
!     Update A := A*P
!
      CALL CLAPMT(forwrd,M,N,A,Lda,Iwork)
!
!     Determine the effective rank of matrix B.
!
      L = 0
      DO i = 1 , MIN(P,N)
         IF ( ABS(B(i,i))>Tolb ) L = L + 1
      ENDDO
!
      IF ( wantv ) THEN
!
!        Copy the details of V, and form V.
!
         CALL CLASET('Full',P,P,CZERO,CZERO,V,Ldv)
         IF ( P>1 ) CALL CLACPY('Lower',P-1,N,B(2,1),Ldb,V(2,1),Ldv)
         CALL CUNG2R(P,P,MIN(P,N),V,Ldv,Tau,Work,Info)
      ENDIF
!
!     Clean up B
!
      DO j = 1 , L - 1
         DO i = j + 1 , L
            B(i,j) = CZERO
         ENDDO
      ENDDO
      IF ( P>L ) CALL CLASET('Full',P-L,N,CZERO,CZERO,B(L+1,1),Ldb)
!
      IF ( wantq ) THEN
!
!        Set Q = I and Update Q := Q*P
!
         CALL CLASET('Full',N,N,CZERO,CONE,Q,Ldq)
         CALL CLAPMT(forwrd,N,N,Q,Ldq,Iwork)
      ENDIF
!
      IF ( P>=L .AND. N/=L ) THEN
!
!        RQ factorization of ( S11 S12 ) = ( 0 S12 )*Z
!
         CALL CGERQ2(L,N,B,Ldb,Tau,Work,Info)
!
!        Update A := A*Z**H
!
         CALL CUNMR2('Right','Conjugate transpose',M,N,L,B,Ldb,Tau,A,   &
     &               Lda,Work,Info)
!
!           Update Q := Q*Z**H
!
         IF ( wantq ) CALL CUNMR2('Right','Conjugate transpose',N,N,L,B,&
     &                            Ldb,Tau,Q,Ldq,Work,Info)
!
!        Clean up B
!
         CALL CLASET('Full',L,N-L,CZERO,CZERO,B,Ldb)
         DO j = N - L + 1 , N
            DO i = j - N + L + 1 , L
               B(i,j) = CZERO
            ENDDO
         ENDDO
!
      ENDIF
!
!     Let              N-L     L
!                A = ( A11    A12 ) M,
!
!     then the following does the complete QR decomposition of A11:
!
!              A11 = U*(  0  T12 )*P1**H
!                      (  0   0  )
!
      DO i = 1 , N - L
         Iwork(i) = 0
      ENDDO
      CALL CGEQP3(M,N-L,A,Lda,Iwork,Tau,Work,Lwork,Rwork,Info)
!
!     Determine the effective rank of A11
!
      K = 0
      DO i = 1 , MIN(M,N-L)
         IF ( ABS(A(i,i))>Tola ) K = K + 1
      ENDDO
!
!     Update A12 := U**H*A12, where A12 = A( 1:M, N-L+1:N )
!
      CALL CUNM2R('Left','Conjugate transpose',M,L,MIN(M,N-L),A,Lda,Tau,&
     &            A(1,N-L+1),Lda,Work,Info)
!
      IF ( wantu ) THEN
!
!        Copy the details of U, and form U
!
         CALL CLASET('Full',M,M,CZERO,CZERO,U,Ldu)
         IF ( M>1 ) CALL CLACPY('Lower',M-1,N-L,A(2,1),Lda,U(2,1),Ldu)
         CALL CUNG2R(M,M,MIN(M,N-L),U,Ldu,Tau,Work,Info)
      ENDIF
!
!
!        Update Q( 1:N, 1:N-L )  = Q( 1:N, 1:N-L )*P1
!
      IF ( wantq ) CALL CLAPMT(forwrd,N,N-L,Q,Ldq,Iwork)
!
!     Clean up A: set the strictly lower triangular part of
!     A(1:K, 1:K) = 0, and A( K+1:M, 1:N-L ) = 0.
!
      DO j = 1 , K - 1
         DO i = j + 1 , K
            A(i,j) = CZERO
         ENDDO
      ENDDO
      IF ( M>K ) CALL CLASET('Full',M-K,N-L,CZERO,CZERO,A(K+1,1),Lda)
!
      IF ( N-L>K ) THEN
!
!        RQ factorization of ( T11 T12 ) = ( 0 T12 )*Z1
!
         CALL CGERQ2(K,N-L,A,Lda,Tau,Work,Info)
!
!
!           Update Q( 1:N,1:N-L ) = Q( 1:N,1:N-L )*Z1**H
!
         IF ( wantq ) CALL CUNMR2('Right','Conjugate transpose',N,N-L,K,&
     &                            A,Lda,Tau,Q,Ldq,Work,Info)
!
!        Clean up A
!
         CALL CLASET('Full',K,N-L-K,CZERO,CZERO,A,Lda)
         DO j = N - L - K + 1 , N - L
            DO i = j - N + L + K + 1 , K
               A(i,j) = CZERO
            ENDDO
         ENDDO
!
      ENDIF
!
      IF ( M>K ) THEN
!
!        QR factorization of A( K+1:M,N-L+1:N )
!
         CALL CGEQR2(M-K,L,A(K+1,N-L+1),Lda,Tau,Work,Info)
!
!
!           Update U(:,K+1:M) := U(:,K+1:M)*U1
!
         IF ( wantu ) CALL CUNM2R('Right','No transpose',M,M-K,         &
     &                            MIN(M-K,L),A(K+1,N-L+1),Lda,Tau,      &
     &                            U(1,K+1),Ldu,Work,Info)
!
!        Clean up
!
         DO j = N - L + 1 , N
            DO i = j - N + K + L + 1 , M
               A(i,j) = CZERO
            ENDDO
         ENDDO
!
      ENDIF
!
      Work(1) = CMPLX(lwkopt)
!
!     End of CGGSVP3
!
      END SUBROUTINE CGGSVP3
