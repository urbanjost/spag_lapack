!*==stbrfs.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b STBRFS
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download STBRFS + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stbrfs.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stbrfs.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stbrfs.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE STBRFS( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, B,
!                          LDB, X, LDX, FERR, BERR, WORK, IWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIAG, TRANS, UPLO
!       INTEGER            INFO, KD, LDAB, LDB, LDX, N, NRHS
!       ..
!       .. Array Arguments ..
!       INTEGER            IWORK( * )
!       REAL               AB( LDAB, * ), B( LDB, * ), BERR( * ),
!      $                   FERR( * ), WORK( * ), X( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> STBRFS provides error bounds and backward error estimates for the
!> solution to a system of linear equations with a triangular band
!> coefficient matrix.
!>
!> The solution matrix X must be computed by STBTRS or some other
!> means before entering this routine.  STBRFS does not do iterative
!> refinement because doing so cannot improve the backward error.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  A is upper triangular;
!>          = 'L':  A is lower triangular.
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          Specifies the form of the system of equations:
!>          = 'N':  A * X = B  (No transpose)
!>          = 'T':  A**T * X = B  (Transpose)
!>          = 'C':  A**H * X = B  (Conjugate transpose = Transpose)
!> \endverbatim
!>
!> \param[in] DIAG
!> \verbatim
!>          DIAG is CHARACTER*1
!>          = 'N':  A is non-unit triangular;
!>          = 'U':  A is unit triangular.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] KD
!> \verbatim
!>          KD is INTEGER
!>          The number of superdiagonals or subdiagonals of the
!>          triangular band matrix A.  KD >= 0.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of right hand sides, i.e., the number of columns
!>          of the matrices B and X.  NRHS >= 0.
!> \endverbatim
!>
!> \param[in] AB
!> \verbatim
!>          AB is REAL array, dimension (LDAB,N)
!>          The upper or lower triangular band matrix A, stored in the
!>          first kd+1 rows of the array. The j-th column of A is stored
!>          in the j-th column of the array AB as follows:
!>          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
!>          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
!>          If DIAG = 'U', the diagonal elements of A are not referenced
!>          and are assumed to be 1.
!> \endverbatim
!>
!> \param[in] LDAB
!> \verbatim
!>          LDAB is INTEGER
!>          The leading dimension of the array AB.  LDAB >= KD+1.
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
!> \param[in] X
!> \verbatim
!>          X is REAL array, dimension (LDX,NRHS)
!>          The solution matrix X.
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
      SUBROUTINE STBRFS(Uplo,Trans,Diag,N,Kd,Nrhs,Ab,Ldab,B,Ldb,X,Ldx,  &
     &                  Ferr,Berr,Work,Iwork,Info)
      USE S_LSAME
      USE S_SAXPY
      USE S_SCOPY
      USE S_SLACN2
      USE S_SLAMCH
      USE S_STBMV
      USE S_STBSV
      USE S_XERBLA
      IMPLICIT NONE
!*--STBRFS200
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER :: N
      INTEGER :: Kd
      INTEGER , INTENT(IN) :: Nrhs
      REAL , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      REAL , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL , INTENT(OUT) , DIMENSION(*) :: Berr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL :: eps , lstres , s , safe1 , safe2 , safmin , xk
      INTEGER :: i , j , k , kase , nz
      INTEGER , DIMENSION(3) :: isave
      LOGICAL :: notran , nounit , upper
      CHARACTER :: transt
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
      notran = LSAME(Trans,'N')
      nounit = LSAME(Diag,'N')
!
      IF ( .NOT.upper .AND. .NOT.LSAME(Uplo,'L') ) THEN
         Info = -1
      ELSEIF ( .NOT.notran .AND. .NOT.LSAME(Trans,'T') .AND.            &
     &         .NOT.LSAME(Trans,'C') ) THEN
         Info = -2
      ELSEIF ( .NOT.nounit .AND. .NOT.LSAME(Diag,'U') ) THEN
         Info = -3
      ELSEIF ( N<0 ) THEN
         Info = -4
      ELSEIF ( Kd<0 ) THEN
         Info = -5
      ELSEIF ( Nrhs<0 ) THEN
         Info = -6
      ELSEIF ( Ldab<Kd+1 ) THEN
         Info = -8
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -10
      ELSEIF ( Ldx<MAX(1,N) ) THEN
         Info = -12
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('STBRFS',-Info)
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
         transt = 'T'
      ELSE
         transt = 'N'
      ENDIF
!
!     NZ = maximum number of nonzero elements in each row of A, plus 1
!
      nz = Kd + 2
      eps = SLAMCH('Epsilon')
      safmin = SLAMCH('Safe minimum')
      safe1 = nz*safmin
      safe2 = safe1/eps
!
!     Do for each right hand side
!
      DO j = 1 , Nrhs
!
!        Compute residual R = B - op(A) * X,
!        where op(A) = A or A**T, depending on TRANS.
!
         CALL SCOPY(N,X(1,j),1,Work(N+1),1)
         CALL STBMV(Uplo,Trans,Diag,N,Kd,Ab,Ldab,Work(N+1),1)
         CALL SAXPY(N,-ONE,B(1,j),1,Work(N+1),1)
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
         DO i = 1 , N
            Work(i) = ABS(B(i,j))
         ENDDO
!
         IF ( notran ) THEN
!
!           Compute abs(A)*abs(X) + abs(B).
!
            IF ( upper ) THEN
               IF ( nounit ) THEN
                  DO k = 1 , N
                     xk = ABS(X(k,j))
                     DO i = MAX(1,k-Kd) , k
                        Work(i) = Work(i) + ABS(Ab(Kd+1+i-k,k))*xk
                     ENDDO
                  ENDDO
               ELSE
                  DO k = 1 , N
                     xk = ABS(X(k,j))
                     DO i = MAX(1,k-Kd) , k - 1
                        Work(i) = Work(i) + ABS(Ab(Kd+1+i-k,k))*xk
                     ENDDO
                     Work(k) = Work(k) + xk
                  ENDDO
               ENDIF
            ELSEIF ( nounit ) THEN
               DO k = 1 , N
                  xk = ABS(X(k,j))
                  DO i = k , MIN(N,k+Kd)
                     Work(i) = Work(i) + ABS(Ab(1+i-k,k))*xk
                  ENDDO
               ENDDO
            ELSE
               DO k = 1 , N
                  xk = ABS(X(k,j))
                  DO i = k + 1 , MIN(N,k+Kd)
                     Work(i) = Work(i) + ABS(Ab(1+i-k,k))*xk
                  ENDDO
                  Work(k) = Work(k) + xk
               ENDDO
            ENDIF
!
!           Compute abs(A**T)*abs(X) + abs(B).
!
         ELSEIF ( upper ) THEN
            IF ( nounit ) THEN
               DO k = 1 , N
                  s = ZERO
                  DO i = MAX(1,k-Kd) , k
                     s = s + ABS(Ab(Kd+1+i-k,k))*ABS(X(i,j))
                  ENDDO
                  Work(k) = Work(k) + s
               ENDDO
            ELSE
               DO k = 1 , N
                  s = ABS(X(k,j))
                  DO i = MAX(1,k-Kd) , k - 1
                     s = s + ABS(Ab(Kd+1+i-k,k))*ABS(X(i,j))
                  ENDDO
                  Work(k) = Work(k) + s
               ENDDO
            ENDIF
         ELSEIF ( nounit ) THEN
            DO k = 1 , N
               s = ZERO
               DO i = k , MIN(N,k+Kd)
                  s = s + ABS(Ab(1+i-k,k))*ABS(X(i,j))
               ENDDO
               Work(k) = Work(k) + s
            ENDDO
         ELSE
            DO k = 1 , N
               s = ABS(X(k,j))
               DO i = k + 1 , MIN(N,k+Kd)
                  s = s + ABS(Ab(1+i-k,k))*ABS(X(i,j))
               ENDDO
               Work(k) = Work(k) + s
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
!        Use SLACN2 to estimate the infinity-norm of the matrix
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
         DO
            CALL SLACN2(N,Work(2*N+1),Work(N+1),Iwork,Ferr(j),kase,     &
     &                  isave)
            IF ( kase/=0 ) THEN
               IF ( kase==1 ) THEN
!
!              Multiply by diag(W)*inv(op(A)**T).
!
                  CALL STBSV(Uplo,transt,Diag,N,Kd,Ab,Ldab,Work(N+1),1)
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
                  CALL STBSV(Uplo,Trans,Diag,N,Kd,Ab,Ldab,Work(N+1),1)
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
!     End of STBRFS
!
      END SUBROUTINE STBRFS
