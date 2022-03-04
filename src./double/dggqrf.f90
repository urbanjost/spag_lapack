!*==dggqrf.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DGGQRF
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DGGQRF + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dggqrf.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dggqrf.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dggqrf.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGGQRF( N, M, P, A, LDA, TAUA, B, LDB, TAUB, WORK,
!                          LWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LDB, LWORK, M, N, P
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), TAUA( * ), TAUB( * ),
!      $                   WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGGQRF computes a generalized QR factorization of an N-by-M matrix A
!> and an N-by-P matrix B:
!>
!>             A = Q*R,        B = Q*T*Z,
!>
!> where Q is an N-by-N orthogonal matrix, Z is a P-by-P orthogonal
!> matrix, and R and T assume one of the forms:
!>
!> if N >= M,  R = ( R11 ) M  ,   or if N < M,  R = ( R11  R12 ) N,
!>                 (  0  ) N-M                         N   M-N
!>                    M
!>
!> where R11 is upper triangular, and
!>
!> if N <= P,  T = ( 0  T12 ) N,   or if N > P,  T = ( T11 ) N-P,
!>                  P-N  N                           ( T21 ) P
!>                                                      P
!>
!> where T12 or T21 is upper triangular.
!>
!> In particular, if B is square and nonsingular, the GQR factorization
!> of A and B implicitly gives the QR factorization of inv(B)*A:
!>
!>              inv(B)*A = Z**T*(inv(T)*R)
!>
!> where inv(B) denotes the inverse of the matrix B, and Z**T denotes the
!> transpose of the matrix Z.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of rows of the matrices A and B. N >= 0.
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
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,M)
!>          On entry, the N-by-M matrix A.
!>          On exit, the elements on and above the diagonal of the array
!>          contain the min(N,M)-by-M upper trapezoidal matrix R (R is
!>          upper triangular if N >= M); the elements below the diagonal,
!>          with the array TAUA, represent the orthogonal matrix Q as a
!>          product of min(N,M) elementary reflectors (see Further
!>          Details).
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A. LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] TAUA
!> \verbatim
!>          TAUA is DOUBLE PRECISION array, dimension (min(N,M))
!>          The scalar factors of the elementary reflectors which
!>          represent the orthogonal matrix Q (see Further Details).
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (LDB,P)
!>          On entry, the N-by-P matrix B.
!>          On exit, if N <= P, the upper triangle of the subarray
!>          B(1:N,P-N+1:P) contains the N-by-N upper triangular matrix T;
!>          if N > P, the elements on and above the (N-P)-th subdiagonal
!>          contain the N-by-P upper trapezoidal matrix T; the remaining
!>          elements, with the array TAUB, represent the orthogonal
!>          matrix Z as a product of elementary reflectors (see Further
!>          Details).
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B. LDB >= max(1,N).
!> \endverbatim
!>
!> \param[out] TAUB
!> \verbatim
!>          TAUB is DOUBLE PRECISION array, dimension (min(N,P))
!>          The scalar factors of the elementary reflectors which
!>          represent the orthogonal matrix Z (see Further Details).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK. LWORK >= max(1,N,M,P).
!>          For optimum performance LWORK >= max(N,M,P)*max(NB1,NB2,NB3),
!>          where NB1 is the optimal blocksize for the QR factorization
!>          of an N-by-M matrix, NB2 is the optimal blocksize for the
!>          RQ factorization of an N-by-P matrix, and NB3 is the optimal
!>          blocksize for a call of DORMQR.
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
!> \date December 2016
!
!> \ingroup doubleOTHERcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The matrix Q is represented as a product of elementary reflectors
!>
!>     Q = H(1) H(2) . . . H(k), where k = min(n,m).
!>
!>  Each H(i) has the form
!>
!>     H(i) = I - taua * v * v**T
!>
!>  where taua is a real scalar, and v is a real vector with
!>  v(1:i-1) = 0 and v(i) = 1; v(i+1:n) is stored on exit in A(i+1:n,i),
!>  and taua in TAUA(i).
!>  To form Q explicitly, use LAPACK subroutine DORGQR.
!>  To use Q to update another matrix, use LAPACK subroutine DORMQR.
!>
!>  The matrix Z is represented as a product of elementary reflectors
!>
!>     Z = H(1) H(2) . . . H(k), where k = min(n,p).
!>
!>  Each H(i) has the form
!>
!>     H(i) = I - taub * v * v**T
!>
!>  where taub is a real scalar, and v is a real vector with
!>  v(p-k+i+1:p) = 0 and v(p-k+i) = 1; v(1:p-k+i-1) is stored on exit in
!>  B(n-k+i,1:p-k+i-1), and taub in TAUB(i).
!>  To form Z explicitly, use LAPACK subroutine DORGRQ.
!>  To use Z to update another matrix, use LAPACK subroutine DORMRQ.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE DGGQRF(N,M,P,A,Lda,Taua,B,Ldb,Taub,Work,Lwork,Info)
      USE F77KINDS                        
      USE S_DGEQRF
      USE S_DGERQF
      USE S_DORMQR
      USE S_ILAENV
      USE S_XERBLA
      IMPLICIT NONE
!*--DGGQRF224
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER :: N
      INTEGER :: M
      INTEGER :: P
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: Taua
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(*) :: Taub
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: lopt , lwkopt , nb , nb1 , nb2 , nb3
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
      nb1 = ILAENV(1,'DGEQRF',' ',N,M,-1,-1)
      nb2 = ILAENV(1,'DGERQF',' ',N,P,-1,-1)
      nb3 = ILAENV(1,'DORMQR',' ',N,M,P,-1)
      nb = MAX(nb1,nb2,nb3)
      lwkopt = MAX(N,M,P)*nb
      Work(1) = lwkopt
      lquery = (Lwork==-1)
      IF ( N<0 ) THEN
         Info = -1
      ELSEIF ( M<0 ) THEN
         Info = -2
      ELSEIF ( P<0 ) THEN
         Info = -3
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -5
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -8
      ELSEIF ( Lwork<MAX(1,N,M,P) .AND. .NOT.lquery ) THEN
         Info = -11
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('DGGQRF',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     QR factorization of N-by-M matrix A: A = Q*R
!
      CALL DGEQRF(N,M,A,Lda,Taua,Work,Lwork,Info)
      lopt = Work(1)
!
!     Update B := Q**T*B.
!
      CALL DORMQR('Left','Transpose',N,P,MIN(N,M),A,Lda,Taua,B,Ldb,Work,&
     &            Lwork,Info)
      lopt = MAX(lopt,INT(Work(1)))
!
!     RQ factorization of N-by-P matrix B: B = T*Z.
!
      CALL DGERQF(N,P,B,Ldb,Taub,Work,Lwork,Info)
      Work(1) = MAX(lopt,INT(Work(1)))
!
!
!     End of DGGQRF
!
      END SUBROUTINE DGGQRF
