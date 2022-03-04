!*==dgbrfs.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b DGBRFS
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DGBRFS + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgbrfs.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgbrfs.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgbrfs.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGBRFS( TRANS, N, KL, KU, NRHS, AB, LDAB, AFB, LDAFB,
!                          IPIV, B, LDB, X, LDX, FERR, BERR, WORK, IWORK,
!                          INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANS
!       INTEGER            INFO, KL, KU, LDAB, LDAFB, LDB, LDX, N, NRHS
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * ), IWORK( * )
!       DOUBLE PRECISION   AB( LDAB, * ), AFB( LDAFB, * ), B( LDB, * ),
!      $                   BERR( * ), FERR( * ), WORK( * ), X( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGBRFS improves the computed solution to a system of linear
!> equations when the coefficient matrix is banded, and provides
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
!> \param[in] KL
!> \verbatim
!>          KL is INTEGER
!>          The number of subdiagonals within the band of A.  KL >= 0.
!> \endverbatim
!>
!> \param[in] KU
!> \verbatim
!>          KU is INTEGER
!>          The number of superdiagonals within the band of A.  KU >= 0.
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
!>          AB is DOUBLE PRECISION array, dimension (LDAB,N)
!>          The original band matrix A, stored in rows 1 to KL+KU+1.
!>          The j-th column of A is stored in the j-th column of the
!>          array AB as follows:
!>          AB(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(n,j+kl).
!> \endverbatim
!>
!> \param[in] LDAB
!> \verbatim
!>          LDAB is INTEGER
!>          The leading dimension of the array AB.  LDAB >= KL+KU+1.
!> \endverbatim
!>
!> \param[in] AFB
!> \verbatim
!>          AFB is DOUBLE PRECISION array, dimension (LDAFB,N)
!>          Details of the LU factorization of the band matrix A, as
!>          computed by DGBTRF.  U is stored as an upper triangular band
!>          matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and
!>          the multipliers used during the factorization are stored in
!>          rows KL+KU+2 to 2*KL+KU+1.
!> \endverbatim
!>
!> \param[in] LDAFB
!> \verbatim
!>          LDAFB is INTEGER
!>          The leading dimension of the array AFB.  LDAFB >= 2*KL*KU+1.
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          The pivot indices from DGBTRF; for 1<=i<=N, row i of the
!>          matrix was interchanged with row IPIV(i).
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
!>          On entry, the solution matrix X, as computed by DGBTRS.
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
!> \ingroup doubleGBcomputational
!
!  =====================================================================
      SUBROUTINE DGBRFS(Trans,N,Kl,Ku,Nrhs,Ab,Ldab,Afb,Ldafb,Ipiv,B,Ldb,&
     &                  X,Ldx,Ferr,Berr,Work,Iwork,Info)
      IMPLICIT NONE
!*--DGBRFS208
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Trans
      INTEGER Info , Kl , Ku , Ldab , Ldafb , Ldb , Ldx , N , Nrhs
!     ..
!     .. Array Arguments ..
      INTEGER Ipiv(*) , Iwork(*)
      DOUBLE PRECISION Ab(Ldab,*) , Afb(Ldafb,*) , B(Ldb,*) , Berr(*) , &
     &                 Ferr(*) , Work(*) , X(Ldx,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER ITMAX
      PARAMETER (ITMAX=5)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D+0)
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D+0)
      DOUBLE PRECISION TWO
      PARAMETER (TWO=2.0D+0)
      DOUBLE PRECISION THREE
      PARAMETER (THREE=3.0D+0)
!     ..
!     .. Local Scalars ..
      LOGICAL notran
      CHARACTER transt
      INTEGER count , i , j , k , kase , kk , nz
      DOUBLE PRECISION eps , lstres , s , safe1 , safe2 , safmin , xk
!     ..
!     .. Local Arrays ..
      INTEGER isave(3)
!     ..
!     .. External Subroutines ..
      EXTERNAL DAXPY , DCOPY , DGBMV , DGBTRS , DLACN2 , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , MIN
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
      ELSEIF ( Kl<0 ) THEN
         Info = -3
      ELSEIF ( Ku<0 ) THEN
         Info = -4
      ELSEIF ( Nrhs<0 ) THEN
         Info = -5
      ELSEIF ( Ldab<Kl+Ku+1 ) THEN
         Info = -7
      ELSEIF ( Ldafb<2*Kl+Ku+1 ) THEN
         Info = -9
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -12
      ELSEIF ( Ldx<MAX(1,N) ) THEN
         Info = -14
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('DGBRFS',-Info)
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
      nz = MIN(Kl+Ku+2,N+1)
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
            CALL DGBMV(Trans,N,N,Kl,Ku,-ONE,Ab,Ldab,X(1,j),1,ONE,       &
     &                 Work(N+1),1)
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
!        Compute abs(op(A))*abs(X) + abs(B).
!
            IF ( notran ) THEN
               DO k = 1 , N
                  kk = Ku + 1 - k
                  xk = ABS(X(k,j))
                  DO i = MAX(1,k-Ku) , MIN(N,k+Kl)
                     Work(i) = Work(i) + ABS(Ab(kk+i,k))*xk
                  ENDDO
               ENDDO
            ELSE
               DO k = 1 , N
                  s = ZERO
                  kk = Ku + 1 - k
                  DO i = MAX(1,k-Ku) , MIN(N,k+Kl)
                     s = s + ABS(Ab(kk+i,k))*ABS(X(i,j))
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
               CALL DGBTRS(Trans,N,Kl,Ku,1,Afb,Ldafb,Ipiv,Work(N+1),N,  &
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
                  CALL DGBTRS(transt,N,Kl,Ku,1,Afb,Ldafb,Ipiv,Work(N+1),&
     &                        N,Info)
                  DO i = 1 , N
                     Work(N+i) = Work(N+i)*Work(i)
                  ENDDO
               ELSE
!
!              Multiply by inv(op(A))*diag(W).
!
                  DO i = 1 , N
                     Work(N+i) = Work(N+i)*Work(i)
                  ENDDO
                  CALL DGBTRS(Trans,N,Kl,Ku,1,Afb,Ldafb,Ipiv,Work(N+1), &
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
!     End of DGBRFS
!
      END SUBROUTINE DGBRFS
