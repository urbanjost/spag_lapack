!*==zcposv.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief <b> ZCPOSV computes the solution to system of linear equations A * X = B for PO matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZCPOSV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zcposv.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zcposv.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zcposv.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZCPOSV( UPLO, N, NRHS, A, LDA, B, LDB, X, LDX, WORK,
!                          SWORK, RWORK, ITER, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, ITER, LDA, LDB, LDX, N, NRHS
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   RWORK( * )
!       COMPLEX            SWORK( * )
!       COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( N, * ),
!      $                   X( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZCPOSV computes the solution to a complex system of linear equations
!>    A * X = B,
!> where A is an N-by-N Hermitian positive definite matrix and X and B
!> are N-by-NRHS matrices.
!>
!> ZCPOSV first attempts to factorize the matrix in COMPLEX and use this
!> factorization within an iterative refinement procedure to produce a
!> solution with COMPLEX*16 normwise backward error quality (see below).
!> If the approach fails the method switches to a COMPLEX*16
!> factorization and solve.
!>
!> The iterative refinement is not going to be a winning strategy if
!> the ratio COMPLEX performance over COMPLEX*16 performance is too
!> small. A reasonable strategy should take the number of right-hand
!> sides and the size of the matrix into account. This might be done
!> with a call to ILAENV in the future. Up to now, we always try
!> iterative refinement.
!>
!> The iterative refinement process is stopped if
!>     ITER > ITERMAX
!> or for all the RHS we have:
!>     RNRM < SQRT(N)*XNRM*ANRM*EPS*BWDMAX
!> where
!>     o ITER is the number of the current iteration in the iterative
!>       refinement process
!>     o RNRM is the infinity-norm of the residual
!>     o XNRM is the infinity-norm of the solution
!>     o ANRM is the infinity-operator-norm of the matrix A
!>     o EPS is the machine epsilon returned by DLAMCH('Epsilon')
!> The value ITERMAX and BWDMAX are fixed to 30 and 1.0D+00
!> respectively.
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
!>          The number of linear equations, i.e., the order of the
!>          matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of right hand sides, i.e., the number of columns
!>          of the matrix B.  NRHS >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array,
!>          dimension (LDA,N)
!>          On entry, the Hermitian matrix A. If UPLO = 'U', the leading
!>          N-by-N upper triangular part of A contains the upper
!>          triangular part of the matrix A, and the strictly lower
!>          triangular part of A is not referenced.  If UPLO = 'L', the
!>          leading N-by-N lower triangular part of A contains the lower
!>          triangular part of the matrix A, and the strictly upper
!>          triangular part of A is not referenced.
!>
!>          Note that the imaginary parts of the diagonal
!>          elements need not be set and are assumed to be zero.
!>
!>          On exit, if iterative refinement has been successfully used
!>          (INFO = 0 and ITER >= 0, see description below), then A is
!>          unchanged, if double precision factorization has been used
!>          (INFO = 0 and ITER < 0, see description below), then the
!>          array A contains the factor U or L from the Cholesky
!>          factorization A = U**H*U or A = L*L**H.
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
!>          The N-by-NRHS right hand side matrix B.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,N).
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is COMPLEX*16 array, dimension (LDX,NRHS)
!>          If INFO = 0, the N-by-NRHS solution matrix X.
!> \endverbatim
!>
!> \param[in] LDX
!> \verbatim
!>          LDX is INTEGER
!>          The leading dimension of the array X.  LDX >= max(1,N).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (N,NRHS)
!>          This array is used to hold the residual vectors.
!> \endverbatim
!>
!> \param[out] SWORK
!> \verbatim
!>          SWORK is COMPLEX array, dimension (N*(N+NRHS))
!>          This array is used to use the single precision matrix and the
!>          right-hand sides or solutions in single precision.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[out] ITER
!> \verbatim
!>          ITER is INTEGER
!>          < 0: iterative refinement has failed, COMPLEX*16
!>               factorization has been performed
!>               -1 : the routine fell back to full precision for
!>                    implementation- or machine-specific reasons
!>               -2 : narrowing the precision induced an overflow,
!>                    the routine fell back to full precision
!>               -3 : failure of CPOTRF
!>               -31: stop the iterative refinement after the 30th
!>                    iterations
!>          > 0: iterative refinement has been successfully used.
!>               Returns the number of iterations
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  if INFO = i, the leading minor of order i of
!>                (COMPLEX*16) A is not positive definite, so the
!>                factorization could not be completed, and the solution
!>                has not been computed.
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
!> \date June 2016
!
!> \ingroup complex16POsolve
!
!  =====================================================================
      SUBROUTINE ZCPOSV(Uplo,N,Nrhs,A,Lda,B,Ldb,X,Ldx,Work,Swork,Rwork, &
     &                  Iter,Info)
      USE F77KINDS                        
      USE S_CLAG2Z
      USE S_CPOTRF
      USE S_CPOTRS
      USE S_DLAMCH
      USE S_IZAMAX
      USE S_LSAME
      USE S_XERBLA
      USE S_ZAXPY
      USE S_ZHEMM
      USE S_ZLACPY
      USE S_ZLAG2C
      USE S_ZLANHE
      USE S_ZLAT2C
      USE S_ZPOTRF
      USE S_ZPOTRS
      IMPLICIT NONE
!*--ZCPOSV229
!
! PARAMETER definitions rewritten by SPAG
!
      LOGICAL , PARAMETER  ::  DOITREF = .TRUE.
      INTEGER , PARAMETER  ::  ITERMAX = 30
      REAL(R8KIND) , PARAMETER  ::  BWDMAX = 1.0E+00
      COMPLEX(CX16KIND) , PARAMETER  ::  NEGONE = (-1.0D+00,0.0D+00) ,  &
     &                 ONE = (1.0D+00,0.0D+00)
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      COMPLEX(CX16KIND) , DIMENSION(N,*) :: Work
      COMPLEX , DIMENSION(*) :: Swork
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(OUT) :: Iter
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: anrm , cte , eps , rnrm , xnrm
      REAL(R8KIND) :: CABS1
      INTEGER :: i , iiter , ptsa , ptsx
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
!
!
!
!
!     .. Local Scalars ..
!
!     .. External Subroutines ..
!     ..
!     .. External Functions ..
!     ..
!     .. Intrinsic Functions ..
!     .. Statement Functions ..
!     ..
!     .. Statement Function definitions ..
      CABS1(zdum) = ABS(DBLE(zdum)) + ABS(DIMAG(zdum))
!     ..
!     .. Executable Statements ..
!
      Info = 0
      Iter = 0
!
!     Test the input parameters.
!
      IF ( .NOT.LSAME(Uplo,'U') .AND. .NOT.LSAME(Uplo,'L') ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Nrhs<0 ) THEN
         Info = -3
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -5
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -7
      ELSEIF ( Ldx<MAX(1,N) ) THEN
         Info = -9
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZCPOSV',-Info)
         RETURN
      ENDIF
!
!     Quick return if (N.EQ.0).
!
      IF ( N==0 ) RETURN
!
!     Skip single precision iterative refinement if a priori slower
!     than double precision factorization.
!
      IF ( .NOT.DOITREF ) THEN
         Iter = -1
         GOTO 300
      ENDIF
!
!     Compute some constants.
!
      anrm = ZLANHE('I',Uplo,N,A,Lda,Rwork)
      eps = DLAMCH('Epsilon')
      cte = anrm*eps*SQRT(DBLE(N))*BWDMAX
!
!     Set the indices PTSA, PTSX for referencing SA and SX in SWORK.
!
      ptsa = 1
      ptsx = ptsa + N*N
!
!     Convert B from double precision to single precision and store the
!     result in SX.
!
      CALL ZLAG2C(N,Nrhs,B,Ldb,Swork(ptsx),N,Info)
!
      IF ( Info/=0 ) THEN
         Iter = -2
         GOTO 300
      ENDIF
!
!     Convert A from double precision to single precision and store the
!     result in SA.
!
      CALL ZLAT2C(Uplo,N,A,Lda,Swork(ptsa),N,Info)
!
      IF ( Info/=0 ) THEN
         Iter = -2
         GOTO 300
      ENDIF
!
!     Compute the Cholesky factorization of SA.
!
      CALL CPOTRF(Uplo,N,Swork(ptsa),N,Info)
!
      IF ( Info/=0 ) THEN
         Iter = -3
         GOTO 300
      ENDIF
!
!     Solve the system SA*SX = SB.
!
      CALL CPOTRS(Uplo,N,Nrhs,Swork(ptsa),N,Swork(ptsx),N,Info)
!
!     Convert SX back to COMPLEX*16
!
      CALL CLAG2Z(N,Nrhs,Swork(ptsx),N,X,Ldx,Info)
!
!     Compute R = B - AX (R is WORK).
!
      CALL ZLACPY('All',N,Nrhs,B,Ldb,Work,N)
!
      CALL ZHEMM('Left',Uplo,N,Nrhs,NEGONE,A,Lda,X,Ldx,ONE,Work,N)
!
!     Check whether the NRHS normwise backward errors satisfy the
!     stopping criterion. If yes, set ITER=0 and return.
!
      DO i = 1 , Nrhs
         xnrm = CABS1(X(IZAMAX(N,X(1,i),1),i))
         rnrm = CABS1(Work(IZAMAX(N,Work(1,i),1),i))
         IF ( rnrm>xnrm*cte ) GOTO 100
      ENDDO
!
!     If we are here, the NRHS normwise backward errors satisfy the
!     stopping criterion. We are good to exit.
!
      Iter = 0
      RETURN
!
!
 100  DO iiter = 1 , ITERMAX
!
!        Convert R (in WORK) from double precision to single precision
!        and store the result in SX.
!
         CALL ZLAG2C(N,Nrhs,Work,N,Swork(ptsx),N,Info)
!
         IF ( Info/=0 ) THEN
            Iter = -2
            GOTO 300
         ENDIF
!
!        Solve the system SA*SX = SR.
!
         CALL CPOTRS(Uplo,N,Nrhs,Swork(ptsa),N,Swork(ptsx),N,Info)
!
!        Convert SX back to double precision and update the current
!        iterate.
!
         CALL CLAG2Z(N,Nrhs,Swork(ptsx),N,Work,N,Info)
!
         DO i = 1 , Nrhs
            CALL ZAXPY(N,ONE,Work(1,i),1,X(1,i),1)
         ENDDO
!
!        Compute R = B - AX (R is WORK).
!
         CALL ZLACPY('All',N,Nrhs,B,Ldb,Work,N)
!
         CALL ZHEMM('L',Uplo,N,Nrhs,NEGONE,A,Lda,X,Ldx,ONE,Work,N)
!
!        Check whether the NRHS normwise backward errors satisfy the
!        stopping criterion. If yes, set ITER=IITER>0 and return.
!
         DO i = 1 , Nrhs
            xnrm = CABS1(X(IZAMAX(N,X(1,i),1),i))
            rnrm = CABS1(Work(IZAMAX(N,Work(1,i),1),i))
            IF ( rnrm>xnrm*cte ) GOTO 200
         ENDDO
!
!        If we are here, the NRHS normwise backward errors satisfy the
!        stopping criterion, we are good to exit.
!
         Iter = iiter
!
         RETURN
!
!
 200  ENDDO
!
!     If we are at this place of the code, this is because we have
!     performed ITER=ITERMAX iterations and never satisfied the
!     stopping criterion, set up the ITER flag accordingly and follow
!     up on double precision routine.
!
      Iter = -ITERMAX - 1
!
!
!     Single-precision iterative refinement failed to converge to a
!     satisfactory solution, so we resort to double precision.
!
 300  CALL ZPOTRF(Uplo,N,A,Lda,Info)
!
      IF ( Info/=0 ) RETURN
!
      CALL ZLACPY('All',N,Nrhs,B,Ldb,X,Ldx)
      CALL ZPOTRS(Uplo,N,Nrhs,A,Lda,X,Ldx,Info)
!
!
!     End of ZCPOSV.
!
      END SUBROUTINE ZCPOSV
