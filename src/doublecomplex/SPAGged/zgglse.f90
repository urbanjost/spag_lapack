!*==zgglse.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief <b> ZGGLSE solves overdetermined or underdetermined systems for OTHER matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZGGLSE + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgglse.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgglse.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgglse.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZGGLSE( M, N, P, A, LDA, B, LDB, C, D, X, WORK, LWORK,
!                          INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LDB, LWORK, M, N, P
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * ), B( LDB, * ), C( * ), D( * ),
!      $                   WORK( * ), X( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZGGLSE solves the linear equality-constrained least squares (LSE)
!> problem:
!>
!>         minimize || c - A*x ||_2   subject to   B*x = d
!>
!> where A is an M-by-N matrix, B is a P-by-N matrix, c is a given
!> M-vector, and d is a given P-vector. It is assumed that
!> P <= N <= M+P, and
!>
!>          rank(B) = P and  rank( (A) ) = N.
!>                               ( (B) )
!>
!> These conditions ensure that the LSE problem has a unique solution,
!> which is obtained using a generalized RQ factorization of the
!> matrices (B, A) given by
!>
!>    B = (0 R)*Q,   A = Z*T*Q.
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
!>          The number of columns of the matrices A and B. N >= 0.
!> \endverbatim
!>
!> \param[in] P
!> \verbatim
!>          P is INTEGER
!>          The number of rows of the matrix B. 0 <= P <= N <= M+P.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          On entry, the M-by-N matrix A.
!>          On exit, the elements on and above the diagonal of the array
!>          contain the min(M,N)-by-N upper trapezoidal matrix T.
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
!>          On exit, the upper triangle of the subarray B(1:P,N-P+1:N)
!>          contains the P-by-P upper triangular matrix R.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B. LDB >= max(1,P).
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is COMPLEX*16 array, dimension (M)
!>          On entry, C contains the right hand side vector for the
!>          least squares part of the LSE problem.
!>          On exit, the residual sum of squares for the solution
!>          is given by the sum of squares of elements N-P+1 to M of
!>          vector C.
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is COMPLEX*16 array, dimension (P)
!>          On entry, D contains the right hand side vector for the
!>          constrained equation.
!>          On exit, D is destroyed.
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is COMPLEX*16 array, dimension (N)
!>          On exit, X is the solution of the LSE problem.
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
!>          The dimension of the array WORK. LWORK >= max(1,M+N+P).
!>          For optimum performance LWORK >= P+min(M,N)+max(M,N)*NB,
!>          where NB is an upper bound for the optimal blocksizes for
!>          ZGEQRF, CGERQF, ZUNMQR and CUNMRQ.
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
!>          = 1:  the upper triangular factor R associated with B in the
!>                generalized RQ factorization of the pair (B, A) is
!>                singular, so that rank(B) < P; the least squares
!>                solution could not be computed.
!>          = 2:  the (N-P) by (N-P) part of the upper trapezoidal factor
!>                T associated with A in the generalized RQ factorization
!>                of the pair (B, A) is singular, so that
!>                rank( (A) ) < N; the least squares solution could not
!>                    ( (B) )
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
!> \ingroup complex16OTHERsolve
!
!  =====================================================================
      SUBROUTINE ZGGLSE(M,N,P,A,Lda,B,Ldb,C,D,X,Work,Lwork,Info)
      USE F77KINDS                        
      USE S_ILAENV
      USE S_XERBLA
      USE S_ZAXPY
      USE S_ZCOPY
      USE S_ZGEMV
      USE S_ZGGRQF
      USE S_ZTRMV
      USE S_ZTRTRS
      USE S_ZUNMQR
      USE S_ZUNMRQ
      IMPLICIT NONE
!*--ZGGLSE194
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0)
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER :: M
      INTEGER :: N
      INTEGER :: P
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(*) :: C
      COMPLEX(CX16KIND) , DIMENSION(*) :: D
      COMPLEX(CX16KIND) , DIMENSION(*) :: X
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: lopt , lwkmin , lwkopt , mn , nb , nb1 , nb2 , nb3 ,   &
     &           nb4 , nr
      LOGICAL :: lquery
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
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
      mn = MIN(M,N)
      lquery = (Lwork==-1)
      IF ( M<0 ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( P<0 .OR. P>N .OR. P<N-M ) THEN
         Info = -3
      ELSEIF ( Lda<MAX(1,M) ) THEN
         Info = -5
      ELSEIF ( Ldb<MAX(1,P) ) THEN
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
            nb1 = ILAENV(1,'ZGEQRF',' ',M,N,-1,-1)
            nb2 = ILAENV(1,'ZGERQF',' ',M,N,-1,-1)
            nb3 = ILAENV(1,'ZUNMQR',' ',M,N,P,-1)
            nb4 = ILAENV(1,'ZUNMRQ',' ',M,N,P,-1)
            nb = MAX(nb1,nb2,nb3,nb4)
            lwkmin = M + N + P
            lwkopt = P + mn + MAX(M,N)*nb
         ENDIF
         Work(1) = lwkopt
!
         IF ( Lwork<lwkmin .AND. .NOT.lquery ) Info = -12
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZGGLSE',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
!     Compute the GRQ factorization of matrices B and A:
!
!            B*Q**H = (  0  T12 ) P   Z**H*A*Q**H = ( R11 R12 ) N-P
!                        N-P  P                     (  0  R22 ) M+P-N
!                                                      N-P  P
!
!     where T12 and R11 are upper triangular, and Q and Z are
!     unitary.
!
      CALL ZGGRQF(P,M,N,B,Ldb,Work,A,Lda,Work(P+1),Work(P+mn+1),        &
     &            Lwork-P-mn,Info)
      lopt = Work(P+mn+1)
!
!     Update c = Z**H *c = ( c1 ) N-P
!                       ( c2 ) M+P-N
!
      CALL ZUNMQR('Left','Conjugate Transpose',M,1,mn,A,Lda,Work(P+1),C,&
     &            MAX(1,M),Work(P+mn+1),Lwork-P-mn,Info)
      lopt = MAX(lopt,INT(Work(P+mn+1)))
!
!     Solve T12*x2 = d for x2
!
      IF ( P>0 ) THEN
         CALL ZTRTRS('Upper','No transpose','Non-unit',P,1,B(1,N-P+1),  &
     &               Ldb,D,P,Info)
!
         IF ( Info>0 ) THEN
            Info = 1
            RETURN
         ENDIF
!
!        Put the solution in X
!
         CALL ZCOPY(P,D,1,X(N-P+1),1)
!
!        Update c1
!
         CALL ZGEMV('No transpose',N-P,P,-CONE,A(1,N-P+1),Lda,D,1,CONE, &
     &              C,1)
      ENDIF
!
!     Solve R11*x1 = c1 for x1
!
      IF ( N>P ) THEN
         CALL ZTRTRS('Upper','No transpose','Non-unit',N-P,1,A,Lda,C,   &
     &               N-P,Info)
!
         IF ( Info>0 ) THEN
            Info = 2
            RETURN
         ENDIF
!
!        Put the solutions in X
!
         CALL ZCOPY(N-P,C,1,X,1)
      ENDIF
!
!     Compute the residual vector:
!
      IF ( M<N ) THEN
         nr = M + P - N
         IF ( nr>0 ) CALL ZGEMV('No transpose',nr,N-M,-CONE,A(N-P+1,M+1)&
     &                          ,Lda,D(nr+1),1,CONE,C(N-P+1),1)
      ELSE
         nr = P
      ENDIF
      IF ( nr>0 ) THEN
         CALL ZTRMV('Upper','No transpose','Non unit',nr,A(N-P+1,N-P+1),&
     &              Lda,D,1)
         CALL ZAXPY(nr,-CONE,D,1,C(N-P+1),1)
      ENDIF
!
!     Backward transformation x = Q**H*x
!
      CALL ZUNMRQ('Left','Conjugate Transpose',N,1,P,B,Ldb,Work(1),X,N, &
     &            Work(P+mn+1),Lwork-P-mn,Info)
      Work(1) = P + mn + MAX(lopt,INT(Work(P+mn+1)))
!
!
!     End of ZGGLSE
!
      END SUBROUTINE ZGGLSE
