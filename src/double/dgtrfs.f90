!*==dgtrfs.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b DGTRFS
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DGTRFS + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgtrfs.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgtrfs.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgtrfs.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGTRFS( TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF, DU2,
!                          IPIV, B, LDB, X, LDX, FERR, BERR, WORK, IWORK,
!                          INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANS
!       INTEGER            INFO, LDB, LDX, N, NRHS
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * ), IWORK( * )
!       DOUBLE PRECISION   B( LDB, * ), BERR( * ), D( * ), DF( * ),
!      $                   DL( * ), DLF( * ), DU( * ), DU2( * ), DUF( * ),
!      $                   FERR( * ), WORK( * ), X( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGTRFS improves the computed solution to a system of linear
!> equations when the coefficient matrix is tridiagonal, and provides
!> error bounds and backward error estimates for the solution.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          Specifies the form of the system of equations:
!>          = 'N':  A * X = B     (No transpose)
!>          = 'T':  A**T * X = B  (Transpose)
!>          = 'C':  A**H * X = B  (Conjugate transpose = Transpose)
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
!> \param[in] DL
!> \verbatim
!>          DL is DOUBLE PRECISION array, dimension (N-1)
!>          The (n-1) subdiagonal elements of A.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>          The diagonal elements of A.
!> \endverbatim
!>
!> \param[in] DU
!> \verbatim
!>          DU is DOUBLE PRECISION array, dimension (N-1)
!>          The (n-1) superdiagonal elements of A.
!> \endverbatim
!>
!> \param[in] DLF
!> \verbatim
!>          DLF is DOUBLE PRECISION array, dimension (N-1)
!>          The (n-1) multipliers that define the matrix L from the
!>          LU factorization of A as computed by DGTTRF.
!> \endverbatim
!>
!> \param[in] DF
!> \verbatim
!>          DF is DOUBLE PRECISION array, dimension (N)
!>          The n diagonal elements of the upper triangular matrix U from
!>          the LU factorization of A.
!> \endverbatim
!>
!> \param[in] DUF
!> \verbatim
!>          DUF is DOUBLE PRECISION array, dimension (N-1)
!>          The (n-1) elements of the first superdiagonal of U.
!> \endverbatim
!>
!> \param[in] DU2
!> \verbatim
!>          DU2 is DOUBLE PRECISION array, dimension (N-2)
!>          The (n-2) elements of the second superdiagonal of U.
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          The pivot indices; for 1 <= i <= n, row i of the matrix was
!>          interchanged with row IPIV(i).  IPIV(i) will always be either
!>          i or i+1; IPIV(i) = i indicates a row interchange was not
!>          required.
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)
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
!>          X is DOUBLE PRECISION array, dimension (LDX,NRHS)
!>          On entry, the solution matrix X, as computed by DGTTRS.
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
!>          FERR is DOUBLE PRECISION array, dimension (NRHS)
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
!>          BERR is DOUBLE PRECISION array, dimension (NRHS)
!>          The componentwise relative backward error of each solution
!>          vector X(j) (i.e., the smallest relative change in
!>          any element of A or B that makes X(j) an exact solution).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (3*N)
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
!> \ingroup doubleGTcomputational
!
!  =====================================================================
      SUBROUTINE DGTRFS(Trans,N,Nrhs,Dl,D,Du,Dlf,Df,Duf,Du2,Ipiv,B,Ldb, &
     &                  X,Ldx,Ferr,Berr,Work,Iwork,Info)
      IMPLICIT NONE
!*--DGTRFS212
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Trans
      INTEGER Info , Ldb , Ldx , N , Nrhs
!     ..
!     .. Array Arguments ..
      INTEGER Ipiv(*) , Iwork(*)
      DOUBLE PRECISION B(Ldb,*) , Berr(*) , D(*) , Df(*) , Dl(*) ,      &
     &                 Dlf(*) , Du(*) , Du2(*) , Duf(*) , Ferr(*) ,     &
     &                 Work(*) , X(Ldx,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER ITMAX
      PARAMETER (ITMAX=5)
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
      DOUBLE PRECISION TWO
      PARAMETER (TWO=2.0D+0)
      DOUBLE PRECISION THREE
      PARAMETER (THREE=3.0D+0)
!     ..
!     .. Local Scalars ..
      LOGICAL notran
      CHARACTER transn , transt
      INTEGER count , i , j , kase , nz
      DOUBLE PRECISION eps , lstres , s , safe1 , safe2 , safmin
!     ..
!     .. Local Arrays ..
      INTEGER isave(3)
!     ..
!     .. External Subroutines ..
      EXTERNAL DAXPY , DCOPY , DGTTRS , DLACN2 , DLAGTM , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      DOUBLE PRECISION DLAMCH
      EXTERNAL LSAME , DLAMCH
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
      notran = LSAME(Trans,'N')
      IF ( .NOT.notran .AND. .NOT.LSAME(Trans,'T') .AND.                &
     &     .NOT.LSAME(Trans,'C') ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Nrhs<0 ) THEN
         Info = -3
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -13
      ELSEIF ( Ldx<MAX(1,N) ) THEN
         Info = -15
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('DGTRFS',-Info)
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
      IF ( notran ) THEN
         transn = 'N'
         transt = 'T'
      ELSE
         transn = 'T'
         transt = 'N'
      ENDIF
!
!     NZ = maximum number of nonzero elements in each row of A, plus 1
!
      nz = 4
      eps = DLAMCH('Epsilon')
      safmin = DLAMCH('Safe minimum')
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
!        Compute residual R = B - op(A) * X,
!        where op(A) = A, A**T, or A**H, depending on TRANS.
!
            CALL DCOPY(N,B(1,j),1,Work(N+1),1)
            CALL DLAGTM(Trans,N,1,-ONE,Dl,D,Du,X(1,j),Ldx,ONE,Work(N+1),&
     &                  N)
!
!        Compute abs(op(A))*abs(x) + abs(b) for use in the backward
!        error bound.
!
            IF ( notran ) THEN
               IF ( N==1 ) THEN
                  Work(1) = ABS(B(1,j)) + ABS(D(1)*X(1,j))
               ELSE
                  Work(1) = ABS(B(1,j)) + ABS(D(1)*X(1,j))              &
     &                      + ABS(Du(1)*X(2,j))
                  DO i = 2 , N - 1
                     Work(i) = ABS(B(i,j)) + ABS(Dl(i-1)*X(i-1,j))      &
     &                         + ABS(D(i)*X(i,j)) + ABS(Du(i)*X(i+1,j))
                  ENDDO
                  Work(N) = ABS(B(N,j)) + ABS(Dl(N-1)*X(N-1,j))         &
     &                      + ABS(D(N)*X(N,j))
               ENDIF
            ELSEIF ( N==1 ) THEN
               Work(1) = ABS(B(1,j)) + ABS(D(1)*X(1,j))
            ELSE
               Work(1) = ABS(B(1,j)) + ABS(D(1)*X(1,j))                 &
     &                   + ABS(Dl(1)*X(2,j))
               DO i = 2 , N - 1
                  Work(i) = ABS(B(i,j)) + ABS(Du(i-1)*X(i-1,j))         &
     &                      + ABS(D(i)*X(i,j)) + ABS(Dl(i)*X(i+1,j))
               ENDDO
               Work(N) = ABS(B(N,j)) + ABS(Du(N-1)*X(N-1,j))            &
     &                   + ABS(D(N)*X(N,j))
            ENDIF
!
!        Compute componentwise relative backward error from formula
!
!        max(i) ( abs(R(i)) / ( abs(op(A))*abs(X) + abs(B) )(i) )
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
               CALL DGTTRS(Trans,N,1,Dlf,Df,Duf,Du2,Ipiv,Work(N+1),N,   &
     &                     Info)
               CALL DAXPY(N,ONE,Work(N+1),1,X(1,j),1)
               lstres = Berr(j)
               count = count + 1
               CYCLE
            ENDIF
!
!        Bound error from formula
!
!        norm(X - XTRUE) / norm(X) .le. FERR =
!        norm( abs(inv(op(A)))*
!           ( abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) ))) / norm(X)
!
!        where
!          norm(Z) is the magnitude of the largest component of Z
!          inv(op(A)) is the inverse of op(A)
!          abs(Z) is the componentwise absolute value of the matrix or
!             vector Z
!          NZ is the maximum number of nonzeros in any row of A, plus 1
!          EPS is machine epsilon
!
!        The i-th component of abs(R)+NZ*EPS*(abs(op(A))*abs(X)+abs(B))
!        is incremented by SAFE1 if the i-th component of
!        abs(op(A))*abs(X) + abs(B) is less than SAFE2.
!
!        Use DLACN2 to estimate the infinity-norm of the matrix
!           inv(op(A)) * diag(W),
!        where W = abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) )))
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
            CALL DLACN2(N,Work(2*N+1),Work(N+1),Iwork,Ferr(j),kase,     &
     &                  isave)
            IF ( kase/=0 ) THEN
               IF ( kase==1 ) THEN
!
!              Multiply by diag(W)*inv(op(A)**T).
!
                  CALL DGTTRS(transt,N,1,Dlf,Df,Duf,Du2,Ipiv,Work(N+1), &
     &                        N,Info)
                  DO i = 1 , N
                     Work(N+i) = Work(i)*Work(N+i)
                  ENDDO
               ELSE
!
!              Multiply by inv(op(A))*diag(W).
!
                  DO i = 1 , N
                     Work(N+i) = Work(i)*Work(N+i)
                  ENDDO
                  CALL DGTTRS(transn,N,1,Dlf,Df,Duf,Du2,Ipiv,Work(N+1), &
     &                        N,Info)
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
!     End of DGTRFS
!
      END SUBROUTINE DGTRFS
