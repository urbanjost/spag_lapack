!*==cherfs.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CHERFS
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CHERFS + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cherfs.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cherfs.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cherfs.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CHERFS( UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB,
!                          X, LDX, FERR, BERR, WORK, RWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       REAL               BERR( * ), FERR( * ), RWORK( * )
!       COMPLEX            A( LDA, * ), AF( LDAF, * ), B( LDB, * ),
!      $                   WORK( * ), X( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CHERFS improves the computed solution to a system of linear
!> equations when the coefficient matrix is Hermitian indefinite, and
!> provides error bounds and backward error estimates for the solution.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  Upper triangle of A is stored;
!>          = 'L':  Lower triangle of A is stored.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of right hand sides, i.e., the number of columns
!>          of the matrices B and X.  NRHS >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          The Hermitian matrix A.  If UPLO = 'U', the leading N-by-N
!>          upper triangular part of A contains the upper triangular part
!>          of the matrix A, and the strictly lower triangular part of A
!>          is not referenced.  If UPLO = 'L', the leading N-by-N lower
!>          triangular part of A contains the lower triangular part of
!>          the matrix A, and the strictly upper triangular part of A is
!>          not referenced.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] AF
!> \verbatim
!>          AF is COMPLEX array, dimension (LDAF,N)
!>          The factored form of the matrix A.  AF contains the block
!>          diagonal matrix D and the multipliers used to obtain the
!>          factor U or L from the factorization A = U*D*U**H or
!>          A = L*D*L**H as computed by CHETRF.
!> \endverbatim
!>
!> \param[in] LDAF
!> \verbatim
!>          LDAF is INTEGER
!>          The leading dimension of the array AF.  LDAF >= max(1,N).
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          Details of the interchanges and the block structure of D
!>          as determined by CHETRF.
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is COMPLEX array, dimension (LDB,NRHS)
!>          The right hand side matrix B.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is COMPLEX array, dimension (LDX,NRHS)
!>          On entry, the solution matrix X, as computed by CHETRS.
!>          On exit, the improved solution matrix X.
!> \endverbatim
!>
!> \param[in] LDX
!> \verbatim
!>          LDX is INTEGER
!>          The leading dimension of the array X.  LDX >= max(1,N).
!> \endverbatim
!>
!> \param[out] FERR
!> \verbatim
!>          FERR is REAL array, dimension (NRHS)
!>          The estimated forward error bound for each solution vector
!>          X(j) (the j-th column of the solution matrix X).
!>          If XTRUE is the true solution corresponding to X(j), FERR(j)
!>          is an estimated upper bound for the magnitude of the largest
!>          element in (X(j) - XTRUE) divided by the magnitude of the
!>          largest element in X(j).  The estimate is as reliable as
!>          the estimate for RCOND, and is almost always a slight
!>          overestimate of the true error.
!> \endverbatim
!>
!> \param[out] BERR
!> \verbatim
!>          BERR is REAL array, dimension (NRHS)
!>          The componentwise relative backward error of each solution
!>          vector X(j) (i.e., the smallest relative change in
!>          any element of A or B that makes X(j) an exact solution).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (2*N)
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!> \endverbatim
!
!> \par Internal Parameters:
!  =========================
!>
!> \verbatim
!>  ITMAX is the maximum number of steps of iterative refinement.
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
!> \ingroup complexHEcomputational
!
!  =====================================================================
      SUBROUTINE CHERFS(Uplo,N,Nrhs,A,Lda,Af,Ldaf,Ipiv,B,Ldb,X,Ldx,Ferr,&
     &                  Berr,Work,Rwork,Info)
      USE S_CAXPY
      USE S_CCOPY
      USE S_CHEMV
      USE S_CHETRS
      USE S_CLACN2
      USE S_LSAME
      USE S_SLAMCH
      USE S_XERBLA
      IMPLICIT NONE
!*--CHERFS204
!
! PARAMETER definitions rewritten by SPAG
!
      INTEGER , PARAMETER  ::  ITMAX = 5
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0)
      REAL , PARAMETER  ::  TWO = 2.0E+0 , THREE = 3.0E+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER , INTENT(IN) :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Berr
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL :: CABS1
      INTEGER :: count , i , j , k , kase , nz
      REAL :: eps , lstres , s , safe1 , safe2 , safmin , xk
      INTEGER , DIMENSION(3) :: isave
      LOGICAL :: upper
      COMPLEX :: zdum
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
!     .. Local Arrays ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. External Functions ..
!     ..
!     .. Statement Functions ..
!     ..
!     .. Statement Function definitions ..
      CABS1(zdum) = ABS(REAL(zdum)) + ABS(AIMAG(zdum))
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
      upper = LSAME(Uplo,'U')
      IF ( .NOT.upper .AND. .NOT.LSAME(Uplo,'L') ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Nrhs<0 ) THEN
         Info = -3
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -5
      ELSEIF ( Ldaf<MAX(1,N) ) THEN
         Info = -7
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -10
      ELSEIF ( Ldx<MAX(1,N) ) THEN
         Info = -12
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('CHERFS',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 .OR. Nrhs==0 ) THEN
         DO j = 1 , Nrhs
            Ferr(j) = ZERO
            Berr(j) = ZERO
         ENDDO
         RETURN
      ENDIF
!
!     NZ = maximum number of nonzero elements in each row of A, plus 1
!
      nz = N + 1
      eps = SLAMCH('Epsilon')
      safmin = SLAMCH('Safe minimum')
      safe1 = nz*safmin
      safe2 = safe1/eps
!
!     Do for each right hand side
!
      DO j = 1 , Nrhs
!
         count = 1
         lstres = THREE
         DO
!
!        Loop until stopping criterion is satisfied.
!
!        Compute residual R = B - A * X
!
            CALL CCOPY(N,B(1,j),1,Work,1)
            CALL CHEMV(Uplo,N,-ONE,A,Lda,X(1,j),1,ONE,Work,1)
!
!        Compute componentwise relative backward error from formula
!
!        max(i) ( abs(R(i)) / ( abs(A)*abs(X) + abs(B) )(i) )
!
!        where abs(Z) is the componentwise absolute value of the matrix
!        or vector Z.  If the i-th component of the denominator is less
!        than SAFE2, then SAFE1 is added to the i-th components of the
!        numerator and denominator before dividing.
!
            DO i = 1 , N
               Rwork(i) = CABS1(B(i,j))
            ENDDO
!
!        Compute abs(A)*abs(X) + abs(B).
!
            IF ( upper ) THEN
               DO k = 1 , N
                  s = ZERO
                  xk = CABS1(X(k,j))
                  DO i = 1 , k - 1
                     Rwork(i) = Rwork(i) + CABS1(A(i,k))*xk
                     s = s + CABS1(A(i,k))*CABS1(X(i,j))
                  ENDDO
                  Rwork(k) = Rwork(k) + ABS(REAL(A(k,k)))*xk + s
               ENDDO
            ELSE
               DO k = 1 , N
                  s = ZERO
                  xk = CABS1(X(k,j))
                  Rwork(k) = Rwork(k) + ABS(REAL(A(k,k)))*xk
                  DO i = k + 1 , N
                     Rwork(i) = Rwork(i) + CABS1(A(i,k))*xk
                     s = s + CABS1(A(i,k))*CABS1(X(i,j))
                  ENDDO
                  Rwork(k) = Rwork(k) + s
               ENDDO
            ENDIF
            s = ZERO
            DO i = 1 , N
               IF ( Rwork(i)>safe2 ) THEN
                  s = MAX(s,CABS1(Work(i))/Rwork(i))
               ELSE
                  s = MAX(s,(CABS1(Work(i))+safe1)/(Rwork(i)+safe1))
               ENDIF
            ENDDO
            Berr(j) = s
!
!        Test stopping criterion. Continue iterating if
!           1) The residual BERR(J) is larger than machine epsilon, and
!           2) BERR(J) decreased by at least a factor of 2 during the
!              last iteration, and
!           3) At most ITMAX iterations tried.
!
            IF ( Berr(j)>eps .AND. TWO*Berr(j)<=lstres .AND.            &
     &           count<=ITMAX ) THEN
!
!           Update solution and try again.
!
               CALL CHETRS(Uplo,N,1,Af,Ldaf,Ipiv,Work,N,Info)
               CALL CAXPY(N,ONE,Work,1,X(1,j),1)
               lstres = Berr(j)
               count = count + 1
               CYCLE
            ENDIF
!
!        Bound error from formula
!
!        norm(X - XTRUE) / norm(X) .le. FERR =
!        norm( abs(inv(A))*
!           ( abs(R) + NZ*EPS*( abs(A)*abs(X)+abs(B) ))) / norm(X)
!
!        where
!          norm(Z) is the magnitude of the largest component of Z
!          inv(A) is the inverse of A
!          abs(Z) is the componentwise absolute value of the matrix or
!             vector Z
!          NZ is the maximum number of nonzeros in any row of A, plus 1
!          EPS is machine epsilon
!
!        The i-th component of abs(R)+NZ*EPS*(abs(A)*abs(X)+abs(B))
!        is incremented by SAFE1 if the i-th component of
!        abs(A)*abs(X) + abs(B) is less than SAFE2.
!
!        Use CLACN2 to estimate the infinity-norm of the matrix
!           inv(A) * diag(W),
!        where W = abs(R) + NZ*EPS*( abs(A)*abs(X)+abs(B) )))
!
            DO i = 1 , N
               IF ( Rwork(i)>safe2 ) THEN
                  Rwork(i) = CABS1(Work(i)) + nz*eps*Rwork(i)
               ELSE
                  Rwork(i) = CABS1(Work(i)) + nz*eps*Rwork(i) + safe1
               ENDIF
            ENDDO
!
            kase = 0
            EXIT
         ENDDO
         DO
            CALL CLACN2(N,Work(N+1),Work,Ferr(j),kase,isave)
            IF ( kase/=0 ) THEN
               IF ( kase==1 ) THEN
!
!              Multiply by diag(W)*inv(A**H).
!
                  CALL CHETRS(Uplo,N,1,Af,Ldaf,Ipiv,Work,N,Info)
                  DO i = 1 , N
                     Work(i) = Rwork(i)*Work(i)
                  ENDDO
               ELSEIF ( kase==2 ) THEN
!
!              Multiply by inv(A)*diag(W).
!
                  DO i = 1 , N
                     Work(i) = Rwork(i)*Work(i)
                  ENDDO
                  CALL CHETRS(Uplo,N,1,Af,Ldaf,Ipiv,Work,N,Info)
               ENDIF
               CYCLE
            ENDIF
!
!        Normalize error.
!
            lstres = ZERO
            DO i = 1 , N
               lstres = MAX(lstres,CABS1(X(i,j)))
            ENDDO
            IF ( lstres/=ZERO ) Ferr(j) = Ferr(j)/lstres
            EXIT
         ENDDO
!
      ENDDO
!
!
!     End of CHERFS
!
      END SUBROUTINE CHERFS
