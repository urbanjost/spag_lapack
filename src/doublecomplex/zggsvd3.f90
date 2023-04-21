!*==zggsvd3.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief <b> ZGGSVD3 computes the singular value decomposition (SVD) for OTHER matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZGGSVD3 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zggsvd3.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zggsvd3.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zggsvd3.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZGGSVD3( JOBU, JOBV, JOBQ, M, N, P, K, L, A, LDA, B,
!                           LDB, ALPHA, BETA, U, LDU, V, LDV, Q, LDQ, WORK,
!                           LWORK, RWORK, IWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBQ, JOBU, JOBV
!       INTEGER            INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N, P, LWORK
!       ..
!       .. Array Arguments ..
!       INTEGER            IWORK( * )
!       DOUBLE PRECISION   ALPHA( * ), BETA( * ), RWORK( * )
!       COMPLEX*16         A( LDA, * ), B( LDB, * ), Q( LDQ, * ),
!      $                   U( LDU, * ), V( LDV, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZGGSVD3 computes the generalized singular value decomposition (GSVD)
!> of an M-by-N complex matrix A and P-by-N complex matrix B:
!>
!>       U**H*A*Q = D1*( 0 R ),    V**H*B*Q = D2*( 0 R )
!>
!> where U, V and Q are unitary matrices.
!> Let K+L = the effective numerical rank of the
!> matrix (A**H,B**H)**H, then R is a (K+L)-by-(K+L) nonsingular upper
!> triangular matrix, D1 and D2 are M-by-(K+L) and P-by-(K+L) "diagonal"
!> matrices and of the following structures, respectively:
!>
!> If M-K-L >= 0,
!>
!>                     K  L
!>        D1 =     K ( I  0 )
!>                 L ( 0  C )
!>             M-K-L ( 0  0 )
!>
!>                   K  L
!>        D2 =   L ( 0  S )
!>             P-L ( 0  0 )
!>
!>                 N-K-L  K    L
!>   ( 0 R ) = K (  0   R11  R12 )
!>             L (  0    0   R22 )
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
!>                   K M-K K+L-M
!>        D1 =   K ( I  0    0   )
!>             M-K ( 0  C    0   )
!>
!>                     K M-K K+L-M
!>        D2 =   M-K ( 0  S    0  )
!>             K+L-M ( 0  0    I  )
!>               P-L ( 0  0    0  )
!>
!>                    N-K-L  K   M-K  K+L-M
!>   ( 0 R ) =     K ( 0    R11  R12  R13  )
!>               M-K ( 0     0   R22  R23  )
!>             K+L-M ( 0     0    0   R33  )
!>
!> where
!>
!>   C = diag( ALPHA(K+1), ... , ALPHA(M) ),
!>   S = diag( BETA(K+1),  ... , BETA(M) ),
!>   C**2 + S**2 = I.
!>
!>   (R11 R12 R13 ) is stored in A(1:M, N-K-L+1:N), and R33 is stored
!>   ( 0  R22 R23 )
!>   in B(M-K+1:L,N+M-K-L+1:N) on exit.
!>
!> The routine computes C, S, R, and optionally the unitary
!> transformation matrices U, V and Q.
!>
!> In particular, if B is an N-by-N nonsingular matrix, then the GSVD of
!> A and B implicitly gives the SVD of A*inv(B):
!>                      A*inv(B) = U*(D1*inv(D2))*V**H.
!> If ( A**H,B**H)**H has orthonormal columns, then the GSVD of A and B is also
!> equal to the CS decomposition of A and B. Furthermore, the GSVD can
!> be used to derive the solution of the eigenvalue problem:
!>                      A**H*A x = lambda* B**H*B x.
!> In some literature, the GSVD of A and B is presented in the form
!>                  U**H*A*X = ( 0 D1 ),   V**H*B*X = ( 0 D2 )
!> where U and V are orthogonal and X is nonsingular, and D1 and D2 are
!> ``diagonal''.  The former GSVD form can be converted to the latter
!> form by taking the nonsingular matrix X as
!>
!>                       X = Q*(  I   0    )
!>                             (  0 inv(R) )
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
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrices A and B.  N >= 0.
!> \endverbatim
!>
!> \param[in] P
!> \verbatim
!>          P is INTEGER
!>          The number of rows of the matrix B.  P >= 0.
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
!>          described in Purpose.
!>          K + L = effective numerical rank of (A**H,B**H)**H.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          On entry, the M-by-N matrix A.
!>          On exit, A contains the triangular matrix R, or part of R.
!>          See Purpose for details.
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
!>          B is COMPLEX*16 array, dimension (LDB,N)
!>          On entry, the P-by-N matrix B.
!>          On exit, B contains part of the triangular matrix R if
!>          M-K-L < 0.  See Purpose for details.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B. LDB >= max(1,P).
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
!>            ALPHA(K+1:K+L) = C,
!>            BETA(K+1:K+L)  = S,
!>          or if M-K-L < 0,
!>            ALPHA(K+1:M)=C, ALPHA(M+1:K+L)=0
!>            BETA(K+1:M) =S, BETA(M+1:K+L) =1
!>          and
!>            ALPHA(K+L+1:N) = 0
!>            BETA(K+L+1:N)  = 0
!> \endverbatim
!>
!> \param[out] U
!> \verbatim
!>          U is COMPLEX*16 array, dimension (LDU,M)
!>          If JOBU = 'U', U contains the M-by-M unitary matrix U.
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
!>          V is COMPLEX*16 array, dimension (LDV,P)
!>          If JOBV = 'V', V contains the P-by-P unitary matrix V.
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
!>          Q is COMPLEX*16 array, dimension (LDQ,N)
!>          If JOBQ = 'Q', Q contains the N-by-N unitary matrix Q.
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
!>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
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
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (2*N)
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (N)
!>          On exit, IWORK stores the sorting information. More
!>          precisely, the following loop will sort ALPHA
!>             for I = K+1, min(M,K+L)
!>                 swap ALPHA(I) and ALPHA(IWORK(I))
!>             endfor
!>          such that ALPHA(1) >= ALPHA(2) >= ... >= ALPHA(N).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          > 0:  if INFO = 1, the Jacobi-type procedure failed to
!>                converge.  For further details, see subroutine ZTGSJA.
!> \endverbatim
!
!> \par Internal Parameters:
!  =========================
!>
!> \verbatim
!>  TOLA    DOUBLE PRECISION
!>  TOLB    DOUBLE PRECISION
!>          TOLA and TOLB are the thresholds to determine the effective
!>          rank of (A**H,B**H)**H. Generally, they are set to
!>                   TOLA = MAX(M,N)*norm(A)*MACHEPS,
!>                   TOLB = MAX(P,N)*norm(B)*MACHEPS.
!>          The size of TOLA and TOLB may affect the size of backward
!>          errors of the decomposition.
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
!> \ingroup complex16GEsing
!
!> \par Contributors:
!  ==================
!>
!>     Ming Gu and Huan Ren, Computer Science Division, University of
!>     California at Berkeley, USA
!>
!
!> \par Further Details:
!  =====================
!>
!>  ZGGSVD3 replaces the deprecated subroutine ZGGSVD.
!>
!  =====================================================================
      SUBROUTINE ZGGSVD3(Jobu,Jobv,Jobq,M,N,P,K,L,A,Lda,B,Ldb,Alpha,    &
     &                   Beta,U,Ldu,V,Ldv,Q,Ldq,Work,Lwork,Rwork,Iwork, &
     &                   Info)
      IMPLICIT NONE
!*--ZGGSVD3357
!
!  -- LAPACK driver routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     August 2015
!
!     .. Scalar Arguments ..
      CHARACTER Jobq , Jobu , Jobv
      INTEGER Info , K , L , Lda , Ldb , Ldq , Ldu , Ldv , M , N , P ,  &
     &        Lwork
!     ..
!     .. Array Arguments ..
      INTEGER Iwork(*)
      DOUBLE PRECISION Alpha(*) , Beta(*) , Rwork(*)
      COMPLEX*16 A(Lda,*) , B(Ldb,*) , Q(Ldq,*) , U(Ldu,*) , V(Ldv,*) , &
     &           Work(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL wantq , wantu , wantv , lquery
      INTEGER i , ibnd , isub , j , ncycle , lwkopt
      DOUBLE PRECISION anorm , bnorm , smax , temp , tola , tolb , ulp ,&
     &                 unfl
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      DOUBLE PRECISION DLAMCH , ZLANGE
      EXTERNAL LSAME , DLAMCH , ZLANGE
!     ..
!     .. External Subroutines ..
      EXTERNAL DCOPY , XERBLA , ZGGSVP3 , ZTGSJA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN
!     ..
!     .. Executable Statements ..
!
!     Decode and test the input parameters
!
      wantu = LSAME(Jobu,'U')
      wantv = LSAME(Jobv,'V')
      wantq = LSAME(Jobq,'Q')
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
      ELSEIF ( N<0 ) THEN
         Info = -5
      ELSEIF ( P<0 ) THEN
         Info = -6
      ELSEIF ( Lda<MAX(1,M) ) THEN
         Info = -10
      ELSEIF ( Ldb<MAX(1,P) ) THEN
         Info = -12
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
         CALL ZGGSVP3(Jobu,Jobv,Jobq,M,P,N,A,Lda,B,Ldb,tola,tolb,K,L,U, &
     &                Ldu,V,Ldv,Q,Ldq,Iwork,Rwork,Work,Work,-1,Info)
         lwkopt = N + INT(Work(1))
         lwkopt = MAX(2*N,lwkopt)
         lwkopt = MAX(1,lwkopt)
         Work(1) = DCMPLX(lwkopt)
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZGGSVD3',-Info)
         RETURN
      ENDIF
      IF ( lquery ) RETURN
!
!     Compute the Frobenius norm of matrices A and B
!
      anorm = ZLANGE('1',M,N,A,Lda,Rwork)
      bnorm = ZLANGE('1',P,N,B,Ldb,Rwork)
!
!     Get machine precision and set up threshold for determining
!     the effective numerical rank of the matrices A and B.
!
      ulp = DLAMCH('Precision')
      unfl = DLAMCH('Safe Minimum')
      tola = MAX(M,N)*MAX(anorm,unfl)*ulp
      tolb = MAX(P,N)*MAX(bnorm,unfl)*ulp
!
      CALL ZGGSVP3(Jobu,Jobv,Jobq,M,P,N,A,Lda,B,Ldb,tola,tolb,K,L,U,Ldu,&
     &             V,Ldv,Q,Ldq,Iwork,Rwork,Work,Work(N+1),Lwork-N,Info)
!
!     Compute the GSVD of two upper "triangular" matrices
!
      CALL ZTGSJA(Jobu,Jobv,Jobq,M,P,N,K,L,A,Lda,B,Ldb,tola,tolb,Alpha, &
     &            Beta,U,Ldu,V,Ldv,Q,Ldq,Work,ncycle,Info)
!
!     Sort the singular values and store the pivot indices in IWORK
!     Copy ALPHA to RWORK, then sort ALPHA in RWORK
!
      CALL DCOPY(N,Alpha,1,Rwork,1)
      ibnd = MIN(L,M-K)
      DO i = 1 , ibnd
!
!        Scan for largest ALPHA(K+I)
!
         isub = i
         smax = Rwork(K+i)
         DO j = i + 1 , ibnd
            temp = Rwork(K+j)
            IF ( temp>smax ) THEN
               isub = j
               smax = temp
            ENDIF
         ENDDO
         IF ( isub/=i ) THEN
            Rwork(K+isub) = Rwork(K+i)
            Rwork(K+i) = smax
            Iwork(K+i) = K + isub
         ELSE
            Iwork(K+i) = K + i
         ENDIF
      ENDDO
!
      Work(1) = DCMPLX(lwkopt)
!
!     End of ZGGSVD3
!
      END SUBROUTINE ZGGSVD3
