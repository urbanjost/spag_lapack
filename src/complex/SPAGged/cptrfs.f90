!*==cptrfs.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CPTRFS
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CPTRFS + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cptrfs.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cptrfs.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cptrfs.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CPTRFS( UPLO, N, NRHS, D, E, DF, EF, B, LDB, X, LDX,
!                          FERR, BERR, WORK, RWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDB, LDX, N, NRHS
!       ..
!       .. Array Arguments ..
!       REAL               BERR( * ), D( * ), DF( * ), FERR( * ),
!      $                   RWORK( * )
!       COMPLEX            B( LDB, * ), E( * ), EF( * ), WORK( * ),
!      $                   X( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CPTRFS improves the computed solution to a system of linear
!> equations when the coefficient matrix is Hermitian positive definite
!> and tridiagonal, and provides error bounds and backward error
!> estimates for the solution.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the superdiagonal or the subdiagonal of the
!>          tridiagonal matrix A is stored and the form of the
!>          factorization:
!>          = 'U':  E is the superdiagonal of A, and A = U**H*D*U;
!>          = 'L':  E is the subdiagonal of A, and A = L*D*L**H.
!>          (The two forms are equivalent if A is real.)
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
!>          of the matrix B.  NRHS >= 0.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is REAL array, dimension (N)
!>          The n real diagonal elements of the tridiagonal matrix A.
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is COMPLEX array, dimension (N-1)
!>          The (n-1) off-diagonal elements of the tridiagonal matrix A
!>          (see UPLO).
!> \endverbatim
!>
!> \param[in] DF
!> \verbatim
!>          DF is REAL array, dimension (N)
!>          The n diagonal elements of the diagonal matrix D from
!>          the factorization computed by CPTTRF.
!> \endverbatim
!>
!> \param[in] EF
!> \verbatim
!>          EF is COMPLEX array, dimension (N-1)
!>          The (n-1) off-diagonal elements of the unit bidiagonal
!>          factor U or L from the factorization computed by CPTTRF
!>          (see UPLO).
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
!>          On entry, the solution matrix X, as computed by CPTTRS.
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
!>          The forward error bound for each solution vector
!>          X(j) (the j-th column of the solution matrix X).
!>          If XTRUE is the true solution corresponding to X(j), FERR(j)
!>          is an estimated upper bound for the magnitude of the largest
!>          element in (X(j) - XTRUE) divided by the magnitude of the
!>          largest element in X(j).
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
!>          WORK is COMPLEX array, dimension (N)
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
!> \ingroup complexPTcomputational
!
!  =====================================================================
      SUBROUTINE CPTRFS(Uplo,N,Nrhs,D,E,Df,Ef,B,Ldb,X,Ldx,Ferr,Berr,    &
     &                  Work,Rwork,Info)
      USE S_CAXPY
      USE S_CPTTRS
      USE S_ISAMAX
      USE S_LSAME
      USE S_SLAMCH
      USE S_XERBLA
      IMPLICIT NONE
!*--CPTRFS193
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
      REAL , INTENT(IN) , DIMENSION(*) :: D
      COMPLEX , INTENT(IN) , DIMENSION(*) :: E
      REAL , DIMENSION(*) :: Df
      COMPLEX , DIMENSION(*) :: Ef
      COMPLEX , INTENT(IN) , DIMENSION(Ldb,*) :: B
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
      COMPLEX :: bi , cx , dx , ex , zdum
      REAL :: CABS1
      INTEGER :: count , i , ix , j , nz
      REAL :: eps , lstres , s , safe1 , safe2 , safmin
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
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
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
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -9
      ELSEIF ( Ldx<MAX(1,N) ) THEN
         Info = -11
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('CPTRFS',-Info)
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
      nz = 4
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
!        Compute residual R = B - A * X.  Also compute
!        abs(A)*abs(x) + abs(b) for use in the backward error bound.
!
            IF ( upper ) THEN
               IF ( N==1 ) THEN
                  bi = B(1,j)
                  dx = D(1)*X(1,j)
                  Work(1) = bi - dx
                  Rwork(1) = CABS1(bi) + CABS1(dx)
               ELSE
                  bi = B(1,j)
                  dx = D(1)*X(1,j)
                  ex = E(1)*X(2,j)
                  Work(1) = bi - dx - ex
                  Rwork(1) = CABS1(bi) + CABS1(dx) + CABS1(E(1))        &
     &                       *CABS1(X(2,j))
                  DO i = 2 , N - 1
                     bi = B(i,j)
                     cx = CONJG(E(i-1))*X(i-1,j)
                     dx = D(i)*X(i,j)
                     ex = E(i)*X(i+1,j)
                     Work(i) = bi - cx - dx - ex
                     Rwork(i) = CABS1(bi) + CABS1(E(i-1))               &
     &                          *CABS1(X(i-1,j)) + CABS1(dx)            &
     &                          + CABS1(E(i))*CABS1(X(i+1,j))
                  ENDDO
                  bi = B(N,j)
                  cx = CONJG(E(N-1))*X(N-1,j)
                  dx = D(N)*X(N,j)
                  Work(N) = bi - cx - dx
                  Rwork(N) = CABS1(bi) + CABS1(E(N-1))*CABS1(X(N-1,j))  &
     &                       + CABS1(dx)
               ENDIF
            ELSEIF ( N==1 ) THEN
               bi = B(1,j)
               dx = D(1)*X(1,j)
               Work(1) = bi - dx
               Rwork(1) = CABS1(bi) + CABS1(dx)
            ELSE
               bi = B(1,j)
               dx = D(1)*X(1,j)
               ex = CONJG(E(1))*X(2,j)
               Work(1) = bi - dx - ex
               Rwork(1) = CABS1(bi) + CABS1(dx) + CABS1(E(1))           &
     &                    *CABS1(X(2,j))
               DO i = 2 , N - 1
                  bi = B(i,j)
                  cx = E(i-1)*X(i-1,j)
                  dx = D(i)*X(i,j)
                  ex = CONJG(E(i))*X(i+1,j)
                  Work(i) = bi - cx - dx - ex
                  Rwork(i) = CABS1(bi) + CABS1(E(i-1))*CABS1(X(i-1,j))  &
     &                       + CABS1(dx) + CABS1(E(i))*CABS1(X(i+1,j))
               ENDDO
               bi = B(N,j)
               cx = E(N-1)*X(N-1,j)
               dx = D(N)*X(N,j)
               Work(N) = bi - cx - dx
               Rwork(N) = CABS1(bi) + CABS1(E(N-1))*CABS1(X(N-1,j))     &
     &                    + CABS1(dx)
            ENDIF
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
               CALL CPTTRS(Uplo,N,1,Df,Ef,Work,N,Info)
               CALL CAXPY(N,CMPLX(ONE),Work,1,X(1,j),1)
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
            DO i = 1 , N
               IF ( Rwork(i)>safe2 ) THEN
                  Rwork(i) = CABS1(Work(i)) + nz*eps*Rwork(i)
               ELSE
                  Rwork(i) = CABS1(Work(i)) + nz*eps*Rwork(i) + safe1
               ENDIF
            ENDDO
            ix = ISAMAX(N,Rwork,1)
            Ferr(j) = Rwork(ix)
!
!        Estimate the norm of inv(A).
!
!        Solve M(A) * x = e, where M(A) = (m(i,j)) is given by
!
!           m(i,j) =  abs(A(i,j)), i = j,
!           m(i,j) = -abs(A(i,j)), i .ne. j,
!
!        and e = [ 1, 1, ..., 1 ]**T.  Note M(A) = M(L)*D*M(L)**H.
!
!        Solve M(L) * x = e.
!
            Rwork(1) = ONE
            DO i = 2 , N
               Rwork(i) = ONE + Rwork(i-1)*ABS(Ef(i-1))
            ENDDO
!
!        Solve D * M(L)**H * x = b.
!
            Rwork(N) = Rwork(N)/Df(N)
            DO i = N - 1 , 1 , -1
               Rwork(i) = Rwork(i)/Df(i) + Rwork(i+1)*ABS(Ef(i))
            ENDDO
!
!        Compute norm(inv(A)) = max(x(i)), 1<=i<=n.
!
            ix = ISAMAX(N,Rwork,1)
            Ferr(j) = Ferr(j)*ABS(Rwork(ix))
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
!     End of CPTRFS
!
      END SUBROUTINE CPTRFS
