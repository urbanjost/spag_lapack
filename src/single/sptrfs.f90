!*==sptrfs.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SPTRFS
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SPTRFS + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sptrfs.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sptrfs.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sptrfs.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SPTRFS( N, NRHS, D, E, DF, EF, B, LDB, X, LDX, FERR,
!                          BERR, WORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDB, LDX, N, NRHS
!       ..
!       .. Array Arguments ..
!       REAL               B( LDB, * ), BERR( * ), D( * ), DF( * ),
!      $                   E( * ), EF( * ), FERR( * ), WORK( * ),
!      $                   X( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SPTRFS improves the computed solution to a system of linear
!> equations when the coefficient matrix is symmetric positive definite
!> and tridiagonal, and provides error bounds and backward error
!> estimates for the solution.
!> \endverbatim
!
!  Arguments:
!  ==========
!
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
!>          The n diagonal elements of the tridiagonal matrix A.
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is REAL array, dimension (N-1)
!>          The (n-1) subdiagonal elements of the tridiagonal matrix A.
!> \endverbatim
!>
!> \param[in] DF
!> \verbatim
!>          DF is REAL array, dimension (N)
!>          The n diagonal elements of the diagonal matrix D from the
!>          factorization computed by SPTTRF.
!> \endverbatim
!>
!> \param[in] EF
!> \verbatim
!>          EF is REAL array, dimension (N-1)
!>          The (n-1) subdiagonal elements of the unit bidiagonal factor
!>          L from the factorization computed by SPTTRF.
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
!>          On entry, the solution matrix X, as computed by SPTTRS.
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
!>          WORK is REAL array, dimension (2*N)
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
!> \ingroup realPTcomputational
!
!  =====================================================================
      SUBROUTINE SPTRFS(N,Nrhs,D,E,Df,Ef,B,Ldb,X,Ldx,Ferr,Berr,Work,    &
     &                  Info)
      IMPLICIT NONE
!*--SPTRFS167
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Info , Ldb , Ldx , N , Nrhs
!     ..
!     .. Array Arguments ..
      REAL B(Ldb,*) , Berr(*) , D(*) , Df(*) , E(*) , Ef(*) , Ferr(*) , &
     &     Work(*) , X(Ldx,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER ITMAX
      PARAMETER (ITMAX=5)
      REAL ZERO
      PARAMETER (ZERO=0.0E+0)
      REAL ONE
      PARAMETER (ONE=1.0E+0)
      REAL TWO
      PARAMETER (TWO=2.0E+0)
      REAL THREE
      PARAMETER (THREE=3.0E+0)
!     ..
!     .. Local Scalars ..
      INTEGER count , i , ix , j , nz
      REAL bi , cx , dx , eps , ex , lstres , s , safe1 , safe2 , safmin
!     ..
!     .. External Subroutines ..
      EXTERNAL SAXPY , SPTTRS , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX
!     ..
!     .. External Functions ..
      INTEGER ISAMAX
      REAL SLAMCH
      EXTERNAL ISAMAX , SLAMCH
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
      IF ( N<0 ) THEN
         Info = -1
      ELSEIF ( Nrhs<0 ) THEN
         Info = -2
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -8
      ELSEIF ( Ldx<MAX(1,N) ) THEN
         Info = -10
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('SPTRFS',-Info)
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
            IF ( N==1 ) THEN
               bi = B(1,j)
               dx = D(1)*X(1,j)
               Work(N+1) = bi - dx
               Work(1) = ABS(bi) + ABS(dx)
            ELSE
               bi = B(1,j)
               dx = D(1)*X(1,j)
               ex = E(1)*X(2,j)
               Work(N+1) = bi - dx - ex
               Work(1) = ABS(bi) + ABS(dx) + ABS(ex)
               DO i = 2 , N - 1
                  bi = B(i,j)
                  cx = E(i-1)*X(i-1,j)
                  dx = D(i)*X(i,j)
                  ex = E(i)*X(i+1,j)
                  Work(N+i) = bi - cx - dx - ex
                  Work(i) = ABS(bi) + ABS(cx) + ABS(dx) + ABS(ex)
               ENDDO
               bi = B(N,j)
               cx = E(N-1)*X(N-1,j)
               dx = D(N)*X(N,j)
               Work(N+N) = bi - cx - dx
               Work(N) = ABS(bi) + ABS(cx) + ABS(dx)
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
               CALL SPTTRS(N,1,Df,Ef,Work(N+1),N,Info)
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
            DO i = 1 , N
               IF ( Work(i)>safe2 ) THEN
                  Work(i) = ABS(Work(N+i)) + nz*eps*Work(i)
               ELSE
                  Work(i) = ABS(Work(N+i)) + nz*eps*Work(i) + safe1
               ENDIF
            ENDDO
            ix = ISAMAX(N,Work,1)
            Ferr(j) = Work(ix)
!
!        Estimate the norm of inv(A).
!
!        Solve M(A) * x = e, where M(A) = (m(i,j)) is given by
!
!           m(i,j) =  abs(A(i,j)), i = j,
!           m(i,j) = -abs(A(i,j)), i .ne. j,
!
!        and e = [ 1, 1, ..., 1 ]**T.  Note M(A) = M(L)*D*M(L)**T.
!
!        Solve M(L) * x = e.
!
            Work(1) = ONE
            DO i = 2 , N
               Work(i) = ONE + Work(i-1)*ABS(Ef(i-1))
            ENDDO
!
!        Solve D * M(L)**T * x = b.
!
            Work(N) = Work(N)/Df(N)
            DO i = N - 1 , 1 , -1
               Work(i) = Work(i)/Df(i) + Work(i+1)*ABS(Ef(i))
            ENDDO
!
!        Compute norm(inv(A)) = max(x(i)), 1<=i<=n.
!
            ix = ISAMAX(N,Work,1)
            Ferr(j) = Ferr(j)*ABS(Work(ix))
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
!     End of SPTRFS
!
      END SUBROUTINE SPTRFS
