!*==dtgsja.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b DTGSJA
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DTGSJA + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtgsja.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtgsja.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtgsja.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DTGSJA( JOBU, JOBV, JOBQ, M, P, N, K, L, A, LDA, B,
!                          LDB, TOLA, TOLB, ALPHA, BETA, U, LDU, V, LDV,
!                          Q, LDQ, WORK, NCYCLE, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBQ, JOBU, JOBV
!       INTEGER            INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N,
!      $                   NCYCLE, P
!       DOUBLE PRECISION   TOLA, TOLB
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), ALPHA( * ), B( LDB, * ),
!      $                   BETA( * ), Q( LDQ, * ), U( LDU, * ),
!      $                   V( LDV, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DTGSJA computes the generalized singular value decomposition (GSVD)
!> of two real upper triangular (or trapezoidal) matrices A and B.
!>
!> On entry, it is assumed that matrices A and B have the following
!> forms, which may be obtained by the preprocessing subroutine DGGSVP
!> from a general M-by-N matrix A and P-by-N matrix B:
!>
!>              N-K-L  K    L
!>    A =    K ( 0    A12  A13 ) if M-K-L >= 0;
!>           L ( 0     0   A23 )
!>       M-K-L ( 0     0    0  )
!>
!>            N-K-L  K    L
!>    A =  K ( 0    A12  A13 ) if M-K-L < 0;
!>       M-K ( 0     0   A23 )
!>
!>            N-K-L  K    L
!>    B =  L ( 0     0   B13 )
!>       P-L ( 0     0    0  )
!>
!> where the K-by-K matrix A12 and L-by-L matrix B13 are nonsingular
!> upper triangular; A23 is L-by-L upper triangular if M-K-L >= 0,
!> otherwise A23 is (M-K)-by-L upper trapezoidal.
!>
!> On exit,
!>
!>        U**T *A*Q = D1*( 0 R ),    V**T *B*Q = D2*( 0 R ),
!>
!> where U, V and Q are orthogonal matrices.
!> R is a nonsingular upper triangular matrix, and D1 and D2 are
!> ``diagonal'' matrices, which are of the following structures:
!>
!> If M-K-L >= 0,
!>
!>                     K  L
!>        D1 =     K ( I  0 )
!>                 L ( 0  C )
!>             M-K-L ( 0  0 )
!>
!>                   K  L
!>        D2 = L   ( 0  S )
!>             P-L ( 0  0 )
!>
!>                N-K-L  K    L
!>   ( 0 R ) = K (  0   R11  R12 ) K
!>             L (  0    0   R22 ) L
!>
!> where
!>
!>   C = diag( ALPHA(K+1), ... , ALPHA(K+L) ),
!>   S = diag( BETA(K+1),  ... , BETA(K+L) ),
!>   C**2 + S**2 = I.
!>
!>   R is stored in A(1:K+L,N-K-L+1:N) on exit.
!>
!> If M-K-L < 0,
!>
!>                K M-K K+L-M
!>     D1 =   K ( I  0    0   )
!>          M-K ( 0  C    0   )
!>
!>                  K M-K K+L-M
!>     D2 =   M-K ( 0  S    0   )
!>          K+L-M ( 0  0    I   )
!>            P-L ( 0  0    0   )
!>
!>                N-K-L  K   M-K  K+L-M
!> ( 0 R ) =    K ( 0    R11  R12  R13  )
!>           M-K ( 0     0   R22  R23  )
!>         K+L-M ( 0     0    0   R33  )
!>
!> where
!> C = diag( ALPHA(K+1), ... , ALPHA(M) ),
!> S = diag( BETA(K+1),  ... , BETA(M) ),
!> C**2 + S**2 = I.
!>
!> R = ( R11 R12 R13 ) is stored in A(1:M, N-K-L+1:N) and R33 is stored
!>     (  0  R22 R23 )
!> in B(M-K+1:L,N+M-K-L+1:N) on exit.
!>
!> The computation of the orthogonal transformation matrices U, V or Q
!> is optional.  These matrices may either be formed explicitly, or they
!> may be postmultiplied into input matrices U1, V1, or Q1.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOBU
!> \verbatim
!>          JOBU is CHARACTER*1
!>          = 'U':  U must contain an orthogonal matrix U1 on entry, and
!>                  the product U1*U is returned;
!>          = 'I':  U is initialized to the unit matrix, and the
!>                  orthogonal matrix U is returned;
!>          = 'N':  U is not computed.
!> \endverbatim
!>
!> \param[in] JOBV
!> \verbatim
!>          JOBV is CHARACTER*1
!>          = 'V':  V must contain an orthogonal matrix V1 on entry, and
!>                  the product V1*V is returned;
!>          = 'I':  V is initialized to the unit matrix, and the
!>                  orthogonal matrix V is returned;
!>          = 'N':  V is not computed.
!> \endverbatim
!>
!> \param[in] JOBQ
!> \verbatim
!>          JOBQ is CHARACTER*1
!>          = 'Q':  Q must contain an orthogonal matrix Q1 on entry, and
!>                  the product Q1*Q is returned;
!>          = 'I':  Q is initialized to the unit matrix, and the
!>                  orthogonal matrix Q is returned;
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
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!> \endverbatim
!>
!> \param[in] L
!> \verbatim
!>          L is INTEGER
!>
!>          K and L specify the subblocks in the input matrices A and B:
!>          A23 = A(K+1:MIN(K+L,M),N-L+1:N) and B13 = B(1:L,N-L+1:N)
!>          of A and B, whose GSVD is going to be computed by DTGSJA.
!>          See Further Details.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, the M-by-N matrix A.
!>          On exit, A(N-K+1:N,1:MIN(K+L,M) ) contains the triangular
!>          matrix R or part of R.  See Purpose for details.
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
!>          B is DOUBLE PRECISION array, dimension (LDB,N)
!>          On entry, the P-by-N matrix B.
!>          On exit, if necessary, B(M-K+1:L,N+M-K-L+1:N) contains
!>          a part of R.  See Purpose for details.
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
!>          TOLA is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in] TOLB
!> \verbatim
!>          TOLB is DOUBLE PRECISION
!>
!>          TOLA and TOLB are the convergence criteria for the Jacobi-
!>          Kogbetliantz iteration procedure. Generally, they are the
!>          same as used in the preprocessing step, say
!>              TOLA = max(M,N)*norm(A)*MAZHEPS,
!>              TOLB = max(P,N)*norm(B)*MAZHEPS.
!> \endverbatim
!>
!> \param[out] ALPHA
!> \verbatim
!>          ALPHA is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[out] BETA
!> \verbatim
!>          BETA is DOUBLE PRECISION array, dimension (N)
!>
!>          On exit, ALPHA and BETA contain the generalized singular
!>          value pairs of A and B;
!>            ALPHA(1:K) = 1,
!>            BETA(1:K)  = 0,
!>          and if M-K-L >= 0,
!>            ALPHA(K+1:K+L) = diag(C),
!>            BETA(K+1:K+L)  = diag(S),
!>          or if M-K-L < 0,
!>            ALPHA(K+1:M)= C, ALPHA(M+1:K+L)= 0
!>            BETA(K+1:M) = S, BETA(M+1:K+L) = 1.
!>          Furthermore, if K+L < N,
!>            ALPHA(K+L+1:N) = 0 and
!>            BETA(K+L+1:N)  = 0.
!> \endverbatim
!>
!> \param[in,out] U
!> \verbatim
!>          U is DOUBLE PRECISION array, dimension (LDU,M)
!>          On entry, if JOBU = 'U', U must contain a matrix U1 (usually
!>          the orthogonal matrix returned by DGGSVP).
!>          On exit,
!>          if JOBU = 'I', U contains the orthogonal matrix U;
!>          if JOBU = 'U', U contains the product U1*U.
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
!> \param[in,out] V
!> \verbatim
!>          V is DOUBLE PRECISION array, dimension (LDV,P)
!>          On entry, if JOBV = 'V', V must contain a matrix V1 (usually
!>          the orthogonal matrix returned by DGGSVP).
!>          On exit,
!>          if JOBV = 'I', V contains the orthogonal matrix V;
!>          if JOBV = 'V', V contains the product V1*V.
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
!> \param[in,out] Q
!> \verbatim
!>          Q is DOUBLE PRECISION array, dimension (LDQ,N)
!>          On entry, if JOBQ = 'Q', Q must contain a matrix Q1 (usually
!>          the orthogonal matrix returned by DGGSVP).
!>          On exit,
!>          if JOBQ = 'I', Q contains the orthogonal matrix Q;
!>          if JOBQ = 'Q', Q contains the product Q1*Q.
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
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (2*N)
!> \endverbatim
!>
!> \param[out] NCYCLE
!> \verbatim
!>          NCYCLE is INTEGER
!>          The number of cycles required for convergence.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          = 1:  the procedure does not converge after MAXIT cycles.
!> \endverbatim
!>
!> \verbatim
!>  Internal Parameters
!>  ===================
!>
!>  MAXIT   INTEGER
!>          MAXIT specifies the total loops that the iterative procedure
!>          may take. If after MAXIT cycles, the routine fails to
!>          converge, we return INFO = 1.
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
!> \ingroup doubleOTHERcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  DTGSJA essentially uses a variant of Kogbetliantz algorithm to reduce
!>  min(L,M-K)-by-L triangular (or trapezoidal) matrix A23 and L-by-L
!>  matrix B13 to the form:
!>
!>           U1**T *A13*Q1 = C1*R1; V1**T *B13*Q1 = S1*R1,
!>
!>  where U1, V1 and Q1 are orthogonal matrix, and Z**T is the transpose
!>  of Z.  C1 and S1 are diagonal matrices satisfying
!>
!>                C1**2 + S1**2 = I,
!>
!>  and R1 is an L-by-L nonsingular upper triangular matrix.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE DTGSJA(Jobu,Jobv,Jobq,M,P,N,K,L,A,Lda,B,Ldb,Tola,Tolb, &
     &                  Alpha,Beta,U,Ldu,V,Ldv,Q,Ldq,Work,Ncycle,Info)
      IMPLICIT NONE
!*--DTGSJA381
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Jobq , Jobu , Jobv
      INTEGER Info , K , L , Lda , Ldb , Ldq , Ldu , Ldv , M , N ,      &
     &        Ncycle , P
      DOUBLE PRECISION Tola , Tolb
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(Lda,*) , Alpha(*) , B(Ldb,*) , Beta(*) ,       &
     &                 Q(Ldq,*) , U(Ldu,*) , V(Ldv,*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER MAXIT
      PARAMETER (MAXIT=40)
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
!     ..
!     .. Local Scalars ..
!
      LOGICAL initq , initu , initv , upper , wantq , wantu , wantv
      INTEGER i , j , kcycle
      DOUBLE PRECISION a1 , a2 , a3 , b1 , b2 , b3 , csq , csu , csv ,  &
     &                 error , gamma , rwk , snq , snu , snv , ssmin
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL DCOPY , DLAGS2 , DLAPLL , DLARTG , DLASET , DROT ,       &
     &         DSCAL , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , MIN
!     ..
!     .. Executable Statements ..
!
!     Decode and test the input parameters
!
      initu = LSAME(Jobu,'I')
      wantu = initu .OR. LSAME(Jobu,'U')
!
      initv = LSAME(Jobv,'I')
      wantv = initv .OR. LSAME(Jobv,'V')
!
      initq = LSAME(Jobq,'I')
      wantq = initq .OR. LSAME(Jobq,'Q')
!
      Info = 0
      IF ( .NOT.(initu .OR. wantu .OR. LSAME(Jobu,'N')) ) THEN
         Info = -1
      ELSEIF ( .NOT.(initv .OR. wantv .OR. LSAME(Jobv,'N')) ) THEN
         Info = -2
      ELSEIF ( .NOT.(initq .OR. wantq .OR. LSAME(Jobq,'N')) ) THEN
         Info = -3
      ELSEIF ( M<0 ) THEN
         Info = -4
      ELSEIF ( P<0 ) THEN
         Info = -5
      ELSEIF ( N<0 ) THEN
         Info = -6
      ELSEIF ( Lda<MAX(1,M) ) THEN
         Info = -10
      ELSEIF ( Ldb<MAX(1,P) ) THEN
         Info = -12
      ELSEIF ( Ldu<1 .OR. (wantu .AND. Ldu<M) ) THEN
         Info = -18
      ELSEIF ( Ldv<1 .OR. (wantv .AND. Ldv<P) ) THEN
         Info = -20
      ELSEIF ( Ldq<1 .OR. (wantq .AND. Ldq<N) ) THEN
         Info = -22
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('DTGSJA',-Info)
         RETURN
      ENDIF
!
!     Initialize U, V and Q, if necessary
!
      IF ( initu ) CALL DLASET('Full',M,M,ZERO,ONE,U,Ldu)
      IF ( initv ) CALL DLASET('Full',P,P,ZERO,ONE,V,Ldv)
      IF ( initq ) CALL DLASET('Full',N,N,ZERO,ONE,Q,Ldq)
!
!     Loop until convergence
!
      upper = .FALSE.
      DO kcycle = 1 , MAXIT
!
         upper = .NOT.upper
!
         DO i = 1 , L - 1
            DO j = i + 1 , L
!
               a1 = ZERO
               a2 = ZERO
               a3 = ZERO
               IF ( K+i<=M ) a1 = A(K+i,N-L+i)
               IF ( K+j<=M ) a3 = A(K+j,N-L+j)
!
               b1 = B(i,N-L+i)
               b3 = B(j,N-L+j)
!
               IF ( upper ) THEN
                  IF ( K+i<=M ) a2 = A(K+i,N-L+j)
                  b2 = B(i,N-L+j)
               ELSE
                  IF ( K+j<=M ) a2 = A(K+j,N-L+i)
                  b2 = B(j,N-L+i)
               ENDIF
!
               CALL DLAGS2(upper,a1,a2,a3,b1,b2,b3,csu,snu,csv,snv,csq, &
     &                     snq)
!
!              Update (K+I)-th and (K+J)-th rows of matrix A: U**T *A
!
               IF ( K+j<=M ) CALL DROT(L,A(K+j,N-L+1),Lda,A(K+i,N-L+1), &
     &                                 Lda,csu,snu)
!
!              Update I-th and J-th rows of matrix B: V**T *B
!
               CALL DROT(L,B(j,N-L+1),Ldb,B(i,N-L+1),Ldb,csv,snv)
!
!              Update (N-L+I)-th and (N-L+J)-th columns of matrices
!              A and B: A*Q and B*Q
!
               CALL DROT(MIN(K+L,M),A(1,N-L+j),1,A(1,N-L+i),1,csq,snq)
!
               CALL DROT(L,B(1,N-L+j),1,B(1,N-L+i),1,csq,snq)
!
               IF ( upper ) THEN
                  IF ( K+i<=M ) A(K+i,N-L+j) = ZERO
                  B(i,N-L+j) = ZERO
               ELSE
                  IF ( K+j<=M ) A(K+j,N-L+i) = ZERO
                  B(j,N-L+i) = ZERO
               ENDIF
!
!              Update orthogonal matrices U, V, Q, if desired.
!
               IF ( wantu .AND. K+j<=M )                                &
     &              CALL DROT(M,U(1,K+j),1,U(1,K+i),1,csu,snu)
!
               IF ( wantv ) CALL DROT(P,V(1,j),1,V(1,i),1,csv,snv)
!
               IF ( wantq ) CALL DROT(N,Q(1,N-L+j),1,Q(1,N-L+i),1,csq,  &
     &                                snq)
!
            ENDDO
         ENDDO
!
         IF ( .NOT.upper ) THEN
!
!           The matrices A13 and B13 were lower triangular at the start
!           of the cycle, and are now upper triangular.
!
!           Convergence test: test the parallelism of the corresponding
!           rows of A and B.
!
            error = ZERO
            DO i = 1 , MIN(L,M-K)
               CALL DCOPY(L-i+1,A(K+i,N-L+i),Lda,Work,1)
               CALL DCOPY(L-i+1,B(i,N-L+i),Ldb,Work(L+1),1)
               CALL DLAPLL(L-i+1,Work,1,Work(L+1),1,ssmin)
               error = MAX(error,ssmin)
            ENDDO
!
            IF ( ABS(error)<=MIN(Tola,Tolb) ) GOTO 100
         ENDIF
!
!        End of cycle loop
!
      ENDDO
!
!     The algorithm has not converged after MAXIT cycles.
!
      Info = 1
      GOTO 200
!
!
!     If ERROR <= MIN(TOLA,TOLB), then the algorithm has converged.
!     Compute the generalized singular value pairs (ALPHA, BETA), and
!     set the triangular matrix R to array A.
!
 100  DO i = 1 , K
         Alpha(i) = ONE
         Beta(i) = ZERO
      ENDDO
!
      DO i = 1 , MIN(L,M-K)
!
         a1 = A(K+i,N-L+i)
         b1 = B(i,N-L+i)
!
         IF ( a1/=ZERO ) THEN
            gamma = b1/a1
!
!           change sign if necessary
!
            IF ( gamma<ZERO ) THEN
               CALL DSCAL(L-i+1,-ONE,B(i,N-L+i),Ldb)
               IF ( wantv ) CALL DSCAL(P,-ONE,V(1,i),1)
            ENDIF
!
            CALL DLARTG(ABS(gamma),ONE,Beta(K+i),Alpha(K+i),rwk)
!
            IF ( Alpha(K+i)>=Beta(K+i) ) THEN
               CALL DSCAL(L-i+1,ONE/Alpha(K+i),A(K+i,N-L+i),Lda)
            ELSE
               CALL DSCAL(L-i+1,ONE/Beta(K+i),B(i,N-L+i),Ldb)
               CALL DCOPY(L-i+1,B(i,N-L+i),Ldb,A(K+i,N-L+i),Lda)
            ENDIF
!
         ELSE
!
            Alpha(K+i) = ZERO
            Beta(K+i) = ONE
            CALL DCOPY(L-i+1,B(i,N-L+i),Ldb,A(K+i,N-L+i),Lda)
!
         ENDIF
!
      ENDDO
!
!     Post-assignment
!
      DO i = M + 1 , K + L
         Alpha(i) = ZERO
         Beta(i) = ONE
      ENDDO
!
      IF ( K+L<N ) THEN
         DO i = K + L + 1 , N
            Alpha(i) = ZERO
            Beta(i) = ZERO
         ENDDO
      ENDIF
!
 200  Ncycle = kcycle
!
!     End of DTGSJA
!
      END SUBROUTINE DTGSJA
