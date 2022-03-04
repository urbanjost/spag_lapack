!*==cggglm.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CGGGLM
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CGGGLM + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cggglm.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cggglm.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cggglm.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGGGLM( N, M, P, A, LDA, B, LDB, D, X, Y, WORK, LWORK,
!                          INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LDB, LWORK, M, N, P
!       ..
!       .. Array Arguments ..
!       COMPLEX            A( LDA, * ), B( LDB, * ), D( * ), WORK( * ),
!      $                   X( * ), Y( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CGGGLM solves a general Gauss-Markov linear model (GLM) problem:
!>
!>         minimize || y ||_2   subject to   d = A*x + B*y
!>             x
!>
!> where A is an N-by-M matrix, B is an N-by-P matrix, and d is a
!> given N-vector. It is assumed that M <= N <= M+P, and
!>
!>            rank(A) = M    and    rank( A B ) = N.
!>
!> Under these assumptions, the constrained equation is always
!> consistent, and there is a unique solution x and a minimal 2-norm
!> solution y, which is obtained using a generalized QR factorization
!> of the matrices (A, B) given by
!>
!>    A = Q*(R),   B = Q*T*Z.
!>          (0)
!>
!> In particular, if matrix B is square nonsingular, then the problem
!> GLM is equivalent to the following weighted linear least squares
!> problem
!>
!>              minimize || inv(B)*(d-A*x) ||_2
!>                  x
!>
!> where inv(B) denotes the inverse of B.
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
!>          The number of columns of the matrix A.  0 <= M <= N.
!> \endverbatim
!>
!> \param[in] P
!> \verbatim
!>          P is INTEGER
!>          The number of columns of the matrix B.  P >= N-M.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,M)
!>          On entry, the N-by-M matrix A.
!>          On exit, the upper triangular part of the array A contains
!>          the M-by-M upper triangular matrix R.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A. LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is COMPLEX array, dimension (LDB,P)
!>          On entry, the N-by-P matrix B.
!>          On exit, if N <= P, the upper triangle of the subarray
!>          B(1:N,P-N+1:P) contains the N-by-N upper triangular matrix T;
!>          if N > P, the elements on and above the (N-P)th subdiagonal
!>          contain the N-by-P upper trapezoidal matrix T.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B. LDB >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is COMPLEX array, dimension (N)
!>          On entry, D is the left hand side of the GLM equation.
!>          On exit, D is destroyed.
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is COMPLEX array, dimension (M)
!> \endverbatim
!>
!> \param[out] Y
!> \verbatim
!>          Y is COMPLEX array, dimension (P)
!>
!>          On exit, X and Y are the solutions of the GLM problem.
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
!>          The dimension of the array WORK. LWORK >= max(1,N+M+P).
!>          For optimum performance, LWORK >= M+min(N,P)+max(N,P)*NB,
!>          where NB is an upper bound for the optimal blocksizes for
!>          CGEQRF, CGERQF, CUNMQR and CUNMRQ.
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
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          = 1:  the upper triangular factor R associated with A in the
!>                generalized QR factorization of the pair (A, B) is
!>                singular, so that rank(A) < M; the least squares
!>                solution could not be computed.
!>          = 2:  the bottom (N-M) by (N-M) part of the upper trapezoidal
!>                factor T associated with B in the generalized QR
!>                factorization of the pair (A, B) is singular, so that
!>                rank( A B ) < N; the least squares solution could not
!>                be computed.
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
!> \ingroup complexOTHEReigen
!
!  =====================================================================
      SUBROUTINE CGGGLM(N,M,P,A,Lda,B,Ldb,D,X,Y,Work,Lwork,Info)
      USE S_CCOPY
      USE S_CGEMV
      USE S_CGGQRF
      USE S_CTRTRS
      USE S_CUNMQR
      USE S_CUNMRQ
      USE S_ILAENV
      USE S_XERBLA
      IMPLICIT NONE
!*--CGGGLM196
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX , PARAMETER  ::  CZERO = (0.0E+0,0.0E+0) ,                &
     &                         CONE = (1.0E+0,0.0E+0)
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER :: N
      INTEGER :: M
      INTEGER :: P
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , DIMENSION(*) :: D
      COMPLEX , DIMENSION(*) :: X
      COMPLEX , DIMENSION(*) :: Y
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , lopt , lwkmin , lwkopt , nb , nb1 , nb2 , nb3 ,    &
     &           nb4 , np
      LOGICAL :: lquery
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  ===================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. External Functions ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      Info = 0
      np = MIN(N,P)
      lquery = (Lwork==-1)
      IF ( N<0 ) THEN
         Info = -1
      ELSEIF ( M<0 .OR. M>N ) THEN
         Info = -2
      ELSEIF ( P<0 .OR. P<N-M ) THEN
         Info = -3
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -5
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -7
      ENDIF
!
!     Calculate workspace
!
      IF ( Info==0 ) THEN
         IF ( N==0 ) THEN
            lwkmin = 1
            lwkopt = 1
         ELSE
            nb1 = ILAENV(1,'CGEQRF',' ',N,M,-1,-1)
            nb2 = ILAENV(1,'CGERQF',' ',N,M,-1,-1)
            nb3 = ILAENV(1,'CUNMQR',' ',N,M,P,-1)
            nb4 = ILAENV(1,'CUNMRQ',' ',N,M,P,-1)
            nb = MAX(nb1,nb2,nb3,nb4)
            lwkmin = M + N + P
            lwkopt = M + np + MAX(N,P)*nb
         ENDIF
         Work(1) = lwkopt
!
         IF ( Lwork<lwkmin .AND. .NOT.lquery ) Info = -12
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('CGGGLM',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) THEN
         DO i = 1 , M
            X(i) = CZERO
         ENDDO
         DO i = 1 , P
            Y(i) = CZERO
         ENDDO
         RETURN
      ENDIF
!
!     Compute the GQR factorization of matrices A and B:
!
!          Q**H*A = ( R11 ) M,    Q**H*B*Z**H = ( T11   T12 ) M
!                   (  0  ) N-M                 (  0    T22 ) N-M
!                      M                         M+P-N  N-M
!
!     where R11 and T22 are upper triangular, and Q and Z are
!     unitary.
!
      CALL CGGQRF(N,M,P,A,Lda,Work,B,Ldb,Work(M+1),Work(M+np+1),        &
     &            Lwork-M-np,Info)
      lopt = Work(M+np+1)
!
!     Update left-hand-side vector d = Q**H*d = ( d1 ) M
!                                               ( d2 ) N-M
!
      CALL CUNMQR('Left','Conjugate transpose',N,1,M,A,Lda,Work,D,      &
     &            MAX(1,N),Work(M+np+1),Lwork-M-np,Info)
      lopt = MAX(lopt,INT(Work(M+np+1)))
!
!     Solve T22*y2 = d2 for y2
!
      IF ( N>M ) THEN
         CALL CTRTRS('Upper','No transpose','Non unit',N-M,1,           &
     &               B(M+1,M+P-N+1),Ldb,D(M+1),N-M,Info)
!
         IF ( Info>0 ) THEN
            Info = 1
            RETURN
         ENDIF
!
         CALL CCOPY(N-M,D(M+1),1,Y(M+P-N+1),1)
      ENDIF
!
!     Set y1 = 0
!
      DO i = 1 , M + P - N
         Y(i) = CZERO
      ENDDO
!
!     Update d1 = d1 - T12*y2
!
      CALL CGEMV('No transpose',M,N-M,-CONE,B(1,M+P-N+1),Ldb,Y(M+P-N+1),&
     &           1,CONE,D,1)
!
!     Solve triangular system: R11*x = d1
!
      IF ( M>0 ) THEN
         CALL CTRTRS('Upper','No Transpose','Non unit',M,1,A,Lda,D,M,   &
     &               Info)
!
         IF ( Info>0 ) THEN
            Info = 2
            RETURN
         ENDIF
!
!        Copy D to X
!
         CALL CCOPY(M,D,1,X,1)
      ENDIF
!
!     Backward transformation y = Z**H *y
!
      CALL CUNMRQ('Left','Conjugate transpose',P,1,np,B(MAX(1,N-P+1),1),&
     &            Ldb,Work(M+1),Y,MAX(1,P),Work(M+np+1),Lwork-M-np,Info)
      Work(1) = M + np + MAX(lopt,INT(Work(M+np+1)))
!
!
!     End of CGGGLM
!
      END SUBROUTINE CGGGLM
