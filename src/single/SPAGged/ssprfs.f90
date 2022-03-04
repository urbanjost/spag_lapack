!*==ssprfs.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SSPRFS
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SSPRFS + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssprfs.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssprfs.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssprfs.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SSPRFS( UPLO, N, NRHS, AP, AFP, IPIV, B, LDB, X, LDX,
!                          FERR, BERR, WORK, IWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDB, LDX, N, NRHS
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * ), IWORK( * )
!       REAL               AFP( * ), AP( * ), B( LDB, * ), BERR( * ),
!      $                   FERR( * ), WORK( * ), X( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SSPRFS improves the computed solution to a system of linear
!> equations when the coefficient matrix is symmetric indefinite
!> and packed, and provides error bounds and backward error estimates
!> for the solution.
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
!> \param[in] AP
!> \verbatim
!>          AP is REAL array, dimension (N*(N+1)/2)
!>          The upper or lower triangle of the symmetric matrix A, packed
!>          columnwise in a linear array.  The j-th column of A is stored
!>          in the array AP as follows:
!>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
!>          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.
!> \endverbatim
!>
!> \param[in] AFP
!> \verbatim
!>          AFP is REAL array, dimension (N*(N+1)/2)
!>          The factored form of the matrix A.  AFP contains the block
!>          diagonal matrix D and the multipliers used to obtain the
!>          factor U or L from the factorization A = U*D*U**T or
!>          A = L*D*L**T as computed by SSPTRF, stored as a packed
!>          triangular matrix.
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          Details of the interchanges and the block structure of D
!>          as determined by SSPTRF.
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is REAL array, dimension (LDB,NRHS)
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
!>          X is REAL array, dimension (LDX,NRHS)
!>          On entry, the solution matrix X, as computed by SSPTRS.
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
!>          WORK is REAL array, dimension (3*N)
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (N)
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
!> \ingroup realOTHERcomputational
!
!  =====================================================================
      SUBROUTINE SSPRFS(Uplo,N,Nrhs,Ap,Afp,Ipiv,B,Ldb,X,Ldx,Ferr,Berr,  &
     &                  Work,Iwork,Info)
      USE S_LSAME
      USE S_SAXPY
      USE S_SCOPY
      USE S_SLACN2
      USE S_SLAMCH
      USE S_SSPMV
      USE S_SSPTRS
      USE S_XERBLA
      IMPLICIT NONE
!*--SSPRFS191
!
! PARAMETER definitions rewritten by SPAG
!
      INTEGER , PARAMETER  ::  ITMAX = 5
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 ,              &
     &                      TWO = 2.0E+0 , THREE = 3.0E+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER , INTENT(IN) :: Nrhs
      REAL , DIMENSION(*) :: Ap
      REAL , DIMENSION(*) :: Afp
      INTEGER , DIMENSION(*) :: Ipiv
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      REAL , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Berr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: count , i , ik , j , k , kase , kk , nz
      REAL :: eps , lstres , s , safe1 , safe2 , safmin , xk
      INTEGER , DIMENSION(3) :: isave
      LOGICAL :: upper
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
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -8
      ELSEIF ( Ldx<MAX(1,N) ) THEN
         Info = -10
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('SSPRFS',-Info)
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
            CALL SCOPY(N,B(1,j),1,Work(N+1),1)
            CALL SSPMV(Uplo,N,-ONE,Ap,X(1,j),1,ONE,Work(N+1),1)
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
               Work(i) = ABS(B(i,j))
            ENDDO
!
!        Compute abs(A)*abs(X) + abs(B).
!
            kk = 1
            IF ( upper ) THEN
               DO k = 1 , N
                  s = ZERO
                  xk = ABS(X(k,j))
                  ik = kk
                  DO i = 1 , k - 1
                     Work(i) = Work(i) + ABS(Ap(ik))*xk
                     s = s + ABS(Ap(ik))*ABS(X(i,j))
                     ik = ik + 1
                  ENDDO
                  Work(k) = Work(k) + ABS(Ap(kk+k-1))*xk + s
                  kk = kk + k
               ENDDO
            ELSE
               DO k = 1 , N
                  s = ZERO
                  xk = ABS(X(k,j))
                  Work(k) = Work(k) + ABS(Ap(kk))*xk
                  ik = kk + 1
                  DO i = k + 1 , N
                     Work(i) = Work(i) + ABS(Ap(ik))*xk
                     s = s + ABS(Ap(ik))*ABS(X(i,j))
                     ik = ik + 1
                  ENDDO
                  Work(k) = Work(k) + s
                  kk = kk + (N-k+1)
               ENDDO
            ENDIF
            s = ZERO
            DO i = 1 , N
               IF ( Work(i)>safe2 ) THEN
                  s = MAX(s,ABS(Work(N+i))/Work(i))
               ELSE
                  s = MAX(s,(ABS(Work(N+i))+safe1)/(Work(i)+safe1))
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
               CALL SSPTRS(Uplo,N,1,Afp,Ipiv,Work(N+1),N,Info)
               CALL SAXPY(N,ONE,Work(N+1),1,X(1,j),1)
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
!        Use SLACN2 to estimate the infinity-norm of the matrix
!           inv(A) * diag(W),
!        where W = abs(R) + NZ*EPS*( abs(A)*abs(X)+abs(B) )))
!
            DO i = 1 , N
               IF ( Work(i)>safe2 ) THEN
                  Work(i) = ABS(Work(N+i)) + nz*eps*Work(i)
               ELSE
                  Work(i) = ABS(Work(N+i)) + nz*eps*Work(i) + safe1
               ENDIF
            ENDDO
!
            kase = 0
            EXIT
         ENDDO
         DO
            CALL SLACN2(N,Work(2*N+1),Work(N+1),Iwork,Ferr(j),kase,     &
     &                  isave)
            IF ( kase/=0 ) THEN
               IF ( kase==1 ) THEN
!
!              Multiply by diag(W)*inv(A**T).
!
                  CALL SSPTRS(Uplo,N,1,Afp,Ipiv,Work(N+1),N,Info)
                  DO i = 1 , N
                     Work(N+i) = Work(i)*Work(N+i)
                  ENDDO
               ELSEIF ( kase==2 ) THEN
!
!              Multiply by inv(A)*diag(W).
!
                  DO i = 1 , N
                     Work(N+i) = Work(i)*Work(N+i)
                  ENDDO
                  CALL SSPTRS(Uplo,N,1,Afp,Ipiv,Work(N+1),N,Info)
               ENDIF
               CYCLE
            ENDIF
!
!        Normalize error.
!
            lstres = ZERO
            DO i = 1 , N
               lstres = MAX(lstres,ABS(X(i,j)))
            ENDDO
            IF ( lstres/=ZERO ) Ferr(j) = Ferr(j)/lstres
            EXIT
         ENDDO
!
      ENDDO
!
!
!     End of SSPRFS
!
      END SUBROUTINE SSPRFS
