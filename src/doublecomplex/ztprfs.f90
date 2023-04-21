!*==ztprfs.f90  processed by SPAG 7.51RB at 20:09 on  3 Mar 2022
!> \brief \b ZTPRFS
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZTPRFS + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztprfs.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztprfs.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztprfs.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZTPRFS( UPLO, TRANS, DIAG, N, NRHS, AP, B, LDB, X, LDX,
!                          FERR, BERR, WORK, RWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIAG, TRANS, UPLO
!       INTEGER            INFO, LDB, LDX, N, NRHS
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   BERR( * ), FERR( * ), RWORK( * )
!       COMPLEX*16         AP( * ), B( LDB, * ), WORK( * ), X( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZTPRFS provides error bounds and backward error estimates for the
!> solution to a system of linear equations with a triangular packed
!> coefficient matrix.
!>
!> The solution matrix X must be computed by ZTPTRS or some other
!> means before entering this routine.  ZTPRFS does not do iterative
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
!> \param[in] AP
!> \verbatim
!>          AP is COMPLEX*16 array, dimension (N*(N+1)/2)
!>          The upper or lower triangular matrix A, packed columnwise in
!>          a linear array.  The j-th column of A is stored in the array
!>          AP as follows:
!>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
!>          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
!>          If DIAG = 'U', the diagonal elements of A are not referenced
!>          and are assumed to be 1.
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
      SUBROUTINE ZTPRFS(Uplo,Trans,Diag,N,Nrhs,Ap,B,Ldb,X,Ldx,Ferr,Berr,&
     &                  Work,Rwork,Info)
      IMPLICIT NONE
!*--ZTPRFS178
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Diag , Trans , Uplo
      INTEGER Info , Ldb , Ldx , N , Nrhs
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION Berr(*) , Ferr(*) , Rwork(*)
      COMPLEX*16 Ap(*) , B(Ldb,*) , Work(*) , X(Ldx,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D+0)
      COMPLEX*16 ONE
      PARAMETER (ONE=(1.0D+0,0.0D+0))
!     ..
!     .. Local Scalars ..
      LOGICAL notran , nounit , upper
      CHARACTER transn , transt
      INTEGER i , j , k , kase , kc , nz
      DOUBLE PRECISION eps , lstres , s , safe1 , safe2 , safmin , xk
      COMPLEX*16 zdum
!     ..
!     .. Local Arrays ..
      INTEGER isave(3)
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA , ZAXPY , ZCOPY , ZLACN2 , ZTPMV , ZTPSV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , DBLE , DIMAG , MAX
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      DOUBLE PRECISION DLAMCH
      EXTERNAL LSAME , DLAMCH
!     ..
!     .. Statement Functions ..
      DOUBLE PRECISION CABS1
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
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -8
      ELSEIF ( Ldx<MAX(1,N) ) THEN
         Info = -10
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZTPRFS',-Info)
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
         CALL ZTPMV(Uplo,Trans,Diag,N,Ap,Work,1)
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
               kc = 1
               IF ( nounit ) THEN
                  DO k = 1 , N
                     xk = CABS1(X(k,j))
                     DO i = 1 , k
                        Rwork(i) = Rwork(i) + CABS1(Ap(kc+i-1))*xk
                     ENDDO
                     kc = kc + k
                  ENDDO
               ELSE
                  DO k = 1 , N
                     xk = CABS1(X(k,j))
                     DO i = 1 , k - 1
                        Rwork(i) = Rwork(i) + CABS1(Ap(kc+i-1))*xk
                     ENDDO
                     Rwork(k) = Rwork(k) + xk
                     kc = kc + k
                  ENDDO
               ENDIF
            ELSE
               kc = 1
               IF ( nounit ) THEN
                  DO k = 1 , N
                     xk = CABS1(X(k,j))
                     DO i = k , N
                        Rwork(i) = Rwork(i) + CABS1(Ap(kc+i-k))*xk
                     ENDDO
                     kc = kc + N - k + 1
                  ENDDO
               ELSE
                  DO k = 1 , N
                     xk = CABS1(X(k,j))
                     DO i = k + 1 , N
                        Rwork(i) = Rwork(i) + CABS1(Ap(kc+i-k))*xk
                     ENDDO
                     Rwork(k) = Rwork(k) + xk
                     kc = kc + N - k + 1
                  ENDDO
               ENDIF
            ENDIF
!
!           Compute abs(A**H)*abs(X) + abs(B).
!
         ELSEIF ( upper ) THEN
            kc = 1
            IF ( nounit ) THEN
               DO k = 1 , N
                  s = ZERO
                  DO i = 1 , k
                     s = s + CABS1(Ap(kc+i-1))*CABS1(X(i,j))
                  ENDDO
                  Rwork(k) = Rwork(k) + s
                  kc = kc + k
               ENDDO
            ELSE
               DO k = 1 , N
                  s = CABS1(X(k,j))
                  DO i = 1 , k - 1
                     s = s + CABS1(Ap(kc+i-1))*CABS1(X(i,j))
                  ENDDO
                  Rwork(k) = Rwork(k) + s
                  kc = kc + k
               ENDDO
            ENDIF
         ELSE
            kc = 1
            IF ( nounit ) THEN
               DO k = 1 , N
                  s = ZERO
                  DO i = k , N
                     s = s + CABS1(Ap(kc+i-k))*CABS1(X(i,j))
                  ENDDO
                  Rwork(k) = Rwork(k) + s
                  kc = kc + N - k + 1
               ENDDO
            ELSE
               DO k = 1 , N
                  s = CABS1(X(k,j))
                  DO i = k + 1 , N
                     s = s + CABS1(Ap(kc+i-k))*CABS1(X(i,j))
                  ENDDO
                  Rwork(k) = Rwork(k) + s
                  kc = kc + N - k + 1
               ENDDO
            ENDIF
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
                  CALL ZTPSV(Uplo,transt,Diag,N,Ap,Work,1)
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
                  CALL ZTPSV(Uplo,transn,Diag,N,Ap,Work,1)
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
!     End of ZTPRFS
!
      END SUBROUTINE ZTPRFS
