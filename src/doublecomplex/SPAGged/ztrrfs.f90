!*==ztrrfs.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZTRRFS
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZTRRFS + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztrrfs.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztrrfs.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztrrfs.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZTRRFS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB, X,
!                          LDX, FERR, BERR, WORK, RWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIAG, TRANS, UPLO
!       INTEGER            INFO, LDA, LDB, LDX, N, NRHS
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   BERR( * ), FERR( * ), RWORK( * )
!       COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * ),
!      $                   X( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZTRRFS provides error bounds and backward error estimates for the
!> solution to a system of linear equations with a triangular
!> coefficient matrix.
!>
!> The solution matrix X must be computed by ZTRTRS or some other
!> means before entering this routine.  ZTRRFS does not do iterative
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
!>          = 'N':  A * X = B     (No transpose)
!>          = 'T':  A**T * X = B  (Transpose)
!>          = 'C':  A**H * X = B  (Conjugate transpose)
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
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of right hand sides, i.e., the number of columns
!>          of the matrices B and X.  NRHS >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          The triangular matrix A.  If UPLO = 'U', the leading N-by-N
!>          upper triangular part of the array A contains the upper
!>          triangular matrix, and the strictly lower triangular part of
!>          A is not referenced.  If UPLO = 'L', the leading N-by-N lower
!>          triangular part of the array A contains the lower triangular
!>          matrix, and the strictly upper triangular part of A is not
!>          referenced.  If DIAG = 'U', the diagonal elements of A are
!>          also not referenced and are assumed to be 1.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is COMPLEX*16 array, dimension (LDB,NRHS)
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
!>          X is COMPLEX*16 array, dimension (LDX,NRHS)
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
!>          WORK is COMPLEX*16 array, dimension (2*N)
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (N)
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
!> \ingroup complex16OTHERcomputational
!
!  =====================================================================
      SUBROUTINE ZTRRFS(Uplo,Trans,Diag,N,Nrhs,A,Lda,B,Ldb,X,Ldx,Ferr,  &
     &                  Berr,Work,Rwork,Info)
      USE F77KINDS                        
      USE S_DLAMCH
      USE S_LSAME
      USE S_XERBLA
      USE S_ZAXPY
      USE S_ZCOPY
      USE S_ZLACN2
      USE S_ZTRMV
      USE S_ZTRSV
      IMPLICIT NONE
!*--ZTRRFS195
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0)
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER :: N
      INTEGER , INTENT(IN) :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: Berr
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: CABS1
      REAL(R8KIND) :: eps , lstres , s , safe1 , safe2 , safmin , xk
      INTEGER :: i , j , k , kase , nz
      INTEGER , DIMENSION(3) :: isave
      LOGICAL :: notran , nounit , upper
      CHARACTER :: transn , transt
      COMPLEX(CX16KIND) :: zdum
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
      CABS1(zdum) = ABS(DBLE(zdum)) + ABS(DIMAG(zdum))
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
      ELSEIF ( Nrhs<0 ) THEN
         Info = -5
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -7
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -9
      ELSEIF ( Ldx<MAX(1,N) ) THEN
         Info = -11
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZTRRFS',-Info)
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
         transt = 'C'
      ELSE
         transn = 'C'
         transt = 'N'
      ENDIF
!
!     NZ = maximum number of nonzero elements in each row of A, plus 1
!
      nz = N + 1
      eps = DLAMCH('Epsilon')
      safmin = DLAMCH('Safe minimum')
      safe1 = nz*safmin
      safe2 = safe1/eps
!
!     Do for each right hand side
!
      DO j = 1 , Nrhs
!
!        Compute residual R = B - op(A) * X,
!        where op(A) = A, A**T, or A**H, depending on TRANS.
!
         CALL ZCOPY(N,X(1,j),1,Work,1)
         CALL ZTRMV(Uplo,Trans,Diag,N,A,Lda,Work,1)
         CALL ZAXPY(N,-ONE,B(1,j),1,Work,1)
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
            Rwork(i) = CABS1(B(i,j))
         ENDDO
!
         IF ( notran ) THEN
!
!           Compute abs(A)*abs(X) + abs(B).
!
            IF ( upper ) THEN
               IF ( nounit ) THEN
                  DO k = 1 , N
                     xk = CABS1(X(k,j))
                     DO i = 1 , k
                        Rwork(i) = Rwork(i) + CABS1(A(i,k))*xk
                     ENDDO
                  ENDDO
               ELSE
                  DO k = 1 , N
                     xk = CABS1(X(k,j))
                     DO i = 1 , k - 1
                        Rwork(i) = Rwork(i) + CABS1(A(i,k))*xk
                     ENDDO
                     Rwork(k) = Rwork(k) + xk
                  ENDDO
               ENDIF
            ELSEIF ( nounit ) THEN
               DO k = 1 , N
                  xk = CABS1(X(k,j))
                  DO i = k , N
                     Rwork(i) = Rwork(i) + CABS1(A(i,k))*xk
                  ENDDO
               ENDDO
            ELSE
               DO k = 1 , N
                  xk = CABS1(X(k,j))
                  DO i = k + 1 , N
                     Rwork(i) = Rwork(i) + CABS1(A(i,k))*xk
                  ENDDO
                  Rwork(k) = Rwork(k) + xk
               ENDDO
            ENDIF
!
!           Compute abs(A**H)*abs(X) + abs(B).
!
         ELSEIF ( upper ) THEN
            IF ( nounit ) THEN
               DO k = 1 , N
                  s = ZERO
                  DO i = 1 , k
                     s = s + CABS1(A(i,k))*CABS1(X(i,j))
                  ENDDO
                  Rwork(k) = Rwork(k) + s
               ENDDO
            ELSE
               DO k = 1 , N
                  s = CABS1(X(k,j))
                  DO i = 1 , k - 1
                     s = s + CABS1(A(i,k))*CABS1(X(i,j))
                  ENDDO
                  Rwork(k) = Rwork(k) + s
               ENDDO
            ENDIF
         ELSEIF ( nounit ) THEN
            DO k = 1 , N
               s = ZERO
               DO i = k , N
                  s = s + CABS1(A(i,k))*CABS1(X(i,j))
               ENDDO
               Rwork(k) = Rwork(k) + s
            ENDDO
         ELSE
            DO k = 1 , N
               s = CABS1(X(k,j))
               DO i = k + 1 , N
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
!        Use ZLACN2 to estimate the infinity-norm of the matrix
!           inv(op(A)) * diag(W),
!        where W = abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) )))
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
         DO
            CALL ZLACN2(N,Work(N+1),Work,Ferr(j),kase,isave)
            IF ( kase/=0 ) THEN
               IF ( kase==1 ) THEN
!
!              Multiply by diag(W)*inv(op(A)**H).
!
                  CALL ZTRSV(Uplo,transt,Diag,N,A,Lda,Work,1)
                  DO i = 1 , N
                     Work(i) = Rwork(i)*Work(i)
                  ENDDO
               ELSE
!
!              Multiply by inv(op(A))*diag(W).
!
                  DO i = 1 , N
                     Work(i) = Rwork(i)*Work(i)
                  ENDDO
                  CALL ZTRSV(Uplo,transn,Diag,N,A,Lda,Work,1)
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
!     End of ZTRRFS
!
      END SUBROUTINE ZTRRFS
