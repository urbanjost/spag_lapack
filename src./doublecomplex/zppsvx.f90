!*==zppsvx.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief <b> ZPPSVX computes the solution to system of linear equations A * X = B for OTHER matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZPPSVX + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zppsvx.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zppsvx.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zppsvx.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZPPSVX( FACT, UPLO, N, NRHS, AP, AFP, EQUED, S, B, LDB,
!                          X, LDX, RCOND, FERR, BERR, WORK, RWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          EQUED, FACT, UPLO
!       INTEGER            INFO, LDB, LDX, N, NRHS
!       DOUBLE PRECISION   RCOND
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   BERR( * ), FERR( * ), RWORK( * ), S( * )
!       COMPLEX*16         AFP( * ), AP( * ), B( LDB, * ), WORK( * ),
!      $                   X( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZPPSVX uses the Cholesky factorization A = U**H * U or A = L * L**H to
!> compute the solution to a complex system of linear equations
!>    A * X = B,
!> where A is an N-by-N Hermitian positive definite matrix stored in
!> packed format and X and B are N-by-NRHS matrices.
!>
!> Error bounds on the solution and a condition estimate are also
!> provided.
!> \endverbatim
!
!> \par Description:
!  =================
!>
!> \verbatim
!>
!> The following steps are performed:
!>
!> 1. If FACT = 'E', real scaling factors are computed to equilibrate
!>    the system:
!>       diag(S) * A * diag(S) * inv(diag(S)) * X = diag(S) * B
!>    Whether or not the system will be equilibrated depends on the
!>    scaling of the matrix A, but if equilibration is used, A is
!>    overwritten by diag(S)*A*diag(S) and B by diag(S)*B.
!>
!> 2. If FACT = 'N' or 'E', the Cholesky decomposition is used to
!>    factor the matrix A (after equilibration if FACT = 'E') as
!>       A = U**H * U ,  if UPLO = 'U', or
!>       A = L * L**H,  if UPLO = 'L',
!>    where U is an upper triangular matrix, L is a lower triangular
!>    matrix, and **H indicates conjugate transpose.
!>
!> 3. If the leading i-by-i principal minor is not positive definite,
!>    then the routine returns with INFO = i. Otherwise, the factored
!>    form of A is used to estimate the condition number of the matrix
!>    A.  If the reciprocal of the condition number is less than machine
!>    precision, INFO = N+1 is returned as a warning, but the routine
!>    still goes on to solve for X and compute error bounds as
!>    described below.
!>
!> 4. The system of equations is solved for X using the factored form
!>    of A.
!>
!> 5. Iterative refinement is applied to improve the computed solution
!>    matrix and calculate error bounds and backward error estimates
!>    for it.
!>
!> 6. If equilibration was used, the matrix X is premultiplied by
!>    diag(S) so that it solves the original system before
!>    equilibration.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] FACT
!> \verbatim
!>          FACT is CHARACTER*1
!>          Specifies whether or not the factored form of the matrix A is
!>          supplied on entry, and if not, whether the matrix A should be
!>          equilibrated before it is factored.
!>          = 'F':  On entry, AFP contains the factored form of A.
!>                  If EQUED = 'Y', the matrix A has been equilibrated
!>                  with scaling factors given by S.  AP and AFP will not
!>                  be modified.
!>          = 'N':  The matrix A will be copied to AFP and factored.
!>          = 'E':  The matrix A will be equilibrated if necessary, then
!>                  copied to AFP and factored.
!> \endverbatim
!>
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
!>          of the matrices B and X.  NRHS >= 0.
!> \endverbatim
!>
!> \param[in,out] AP
!> \verbatim
!>          AP is COMPLEX*16 array, dimension (N*(N+1)/2)
!>          On entry, the upper or lower triangle of the Hermitian matrix
!>          A, packed columnwise in a linear array, except if FACT = 'F'
!>          and EQUED = 'Y', then A must contain the equilibrated matrix
!>          diag(S)*A*diag(S).  The j-th column of A is stored in the
!>          array AP as follows:
!>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
!>          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
!>          See below for further details.  A is not modified if
!>          FACT = 'F' or 'N', or if FACT = 'E' and EQUED = 'N' on exit.
!>
!>          On exit, if FACT = 'E' and EQUED = 'Y', A is overwritten by
!>          diag(S)*A*diag(S).
!> \endverbatim
!>
!> \param[in,out] AFP
!> \verbatim
!>          AFP is COMPLEX*16 array, dimension (N*(N+1)/2)
!>          If FACT = 'F', then AFP is an input argument and on entry
!>          contains the triangular factor U or L from the Cholesky
!>          factorization A = U**H*U or A = L*L**H, in the same storage
!>          format as A.  If EQUED .ne. 'N', then AFP is the factored
!>          form of the equilibrated matrix A.
!>
!>          If FACT = 'N', then AFP is an output argument and on exit
!>          returns the triangular factor U or L from the Cholesky
!>          factorization A = U**H * U or A = L * L**H of the original
!>          matrix A.
!>
!>          If FACT = 'E', then AFP is an output argument and on exit
!>          returns the triangular factor U or L from the Cholesky
!>          factorization A = U**H * U or A = L * L**H of the equilibrated
!>          matrix A (see the description of AP for the form of the
!>          equilibrated matrix).
!> \endverbatim
!>
!> \param[in,out] EQUED
!> \verbatim
!>          EQUED is CHARACTER*1
!>          Specifies the form of equilibration that was done.
!>          = 'N':  No equilibration (always true if FACT = 'N').
!>          = 'Y':  Equilibration was done, i.e., A has been replaced by
!>                  diag(S) * A * diag(S).
!>          EQUED is an input argument if FACT = 'F'; otherwise, it is an
!>          output argument.
!> \endverbatim
!>
!> \param[in,out] S
!> \verbatim
!>          S is DOUBLE PRECISION array, dimension (N)
!>          The scale factors for A; not accessed if EQUED = 'N'.  S is
!>          an input argument if FACT = 'F'; otherwise, S is an output
!>          argument.  If FACT = 'F' and EQUED = 'Y', each element of S
!>          must be positive.
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is COMPLEX*16 array, dimension (LDB,NRHS)
!>          On entry, the N-by-NRHS right hand side matrix B.
!>          On exit, if EQUED = 'N', B is not modified; if EQUED = 'Y',
!>          B is overwritten by diag(S) * B.
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
!>          If INFO = 0 or INFO = N+1, the N-by-NRHS solution matrix X to
!>          the original system of equations.  Note that if EQUED = 'Y',
!>          A and B are modified on exit, and the solution to the
!>          equilibrated system is inv(diag(S))*X.
!> \endverbatim
!>
!> \param[in] LDX
!> \verbatim
!>          LDX is INTEGER
!>          The leading dimension of the array X.  LDX >= max(1,N).
!> \endverbatim
!>
!> \param[out] RCOND
!> \verbatim
!>          RCOND is DOUBLE PRECISION
!>          The estimate of the reciprocal condition number of the matrix
!>          A after equilibration (if done).  If RCOND is less than the
!>          machine precision (in particular, if RCOND = 0), the matrix
!>          is singular to working precision.  This condition is
!>          indicated by a return code of INFO > 0.
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
!>          > 0:  if INFO = i, and i is
!>                <= N:  the leading minor of order i of A is
!>                       not positive definite, so the factorization
!>                       could not be completed, and the solution has not
!>                       been computed. RCOND = 0 is returned.
!>                = N+1: U is nonsingular, but RCOND is less than machine
!>                       precision, meaning that the matrix is singular
!>                       to working precision.  Nevertheless, the
!>                       solution and error bounds are computed because
!>                       there are a number of situations where the
!>                       computed solution can be more accurate than the
!>                       value of RCOND would suggest.
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
!> \date April 2012
!
!> \ingroup complex16OTHERsolve
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The packed storage scheme is illustrated by the following example
!>  when N = 4, UPLO = 'U':
!>
!>  Two-dimensional storage of the Hermitian matrix A:
!>
!>     a11 a12 a13 a14
!>         a22 a23 a24
!>             a33 a34     (aij = conjg(aji))
!>                 a44
!>
!>  Packed storage of the upper triangle of A:
!>
!>  AP = [ a11, a12, a22, a13, a23, a33, a14, a24, a34, a44 ]
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE ZPPSVX(Fact,Uplo,N,Nrhs,Ap,Afp,Equed,S,B,Ldb,X,Ldx,    &
     &                  Rcond,Ferr,Berr,Work,Rwork,Info)
      USE F77KINDS                        
      USE S_DLAMCH
      USE S_LSAME
      USE S_XERBLA
      USE S_ZCOPY
      USE S_ZLACPY
      USE S_ZLANHP
      USE S_ZLAQHP
      USE S_ZPPCON
      USE S_ZPPEQU
      USE S_ZPPRFS
      USE S_ZPPTRF
      USE S_ZPPTRS
      IMPLICIT NONE
!*--ZPPSVX328
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Fact
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(*) :: Ap
      COMPLEX(CX16KIND) , DIMENSION(*) :: Afp
      CHARACTER :: Equed
      REAL(R8KIND) , DIMENSION(*) :: S
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL(R8KIND) , INTENT(INOUT) :: Rcond
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL(R8KIND) , DIMENSION(*) :: Berr
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: amax , anorm , bignum , scond , smax , smin ,     &
     &                smlnum
      LOGICAL :: equil , nofact , rcequ
      INTEGER :: i , infequ , j
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
!     .. Executable Statements ..
!
      Info = 0
      nofact = LSAME(Fact,'N')
      equil = LSAME(Fact,'E')
      IF ( nofact .OR. equil ) THEN
         Equed = 'N'
         rcequ = .FALSE.
      ELSE
         rcequ = LSAME(Equed,'Y')
         smlnum = DLAMCH('Safe minimum')
         bignum = ONE/smlnum
      ENDIF
!
!     Test the input parameters.
!
      IF ( .NOT.nofact .AND. .NOT.equil .AND. .NOT.LSAME(Fact,'F') )    &
     &     THEN
         Info = -1
      ELSEIF ( .NOT.LSAME(Uplo,'U') .AND. .NOT.LSAME(Uplo,'L') ) THEN
         Info = -2
      ELSEIF ( N<0 ) THEN
         Info = -3
      ELSEIF ( Nrhs<0 ) THEN
         Info = -4
      ELSEIF ( LSAME(Fact,'F') .AND. .NOT.(rcequ .OR. LSAME(Equed,'N')) &
     &         ) THEN
         Info = -7
      ELSE
         IF ( rcequ ) THEN
            smin = bignum
            smax = ZERO
            DO j = 1 , N
               smin = MIN(smin,S(j))
               smax = MAX(smax,S(j))
            ENDDO
            IF ( smin<=ZERO ) THEN
               Info = -8
            ELSEIF ( N>0 ) THEN
               scond = MAX(smin,smlnum)/MIN(smax,bignum)
            ELSE
               scond = ONE
            ENDIF
         ENDIF
         IF ( Info==0 ) THEN
            IF ( Ldb<MAX(1,N) ) THEN
               Info = -10
            ELSEIF ( Ldx<MAX(1,N) ) THEN
               Info = -12
            ENDIF
         ENDIF
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZPPSVX',-Info)
         RETURN
      ENDIF
!
      IF ( equil ) THEN
!
!        Compute row and column scalings to equilibrate the matrix A.
!
         CALL ZPPEQU(Uplo,N,Ap,S,scond,amax,infequ)
         IF ( infequ==0 ) THEN
!
!           Equilibrate the matrix.
!
            CALL ZLAQHP(Uplo,N,Ap,S,scond,amax,Equed)
            rcequ = LSAME(Equed,'Y')
         ENDIF
      ENDIF
!
!     Scale the right-hand side.
!
      IF ( rcequ ) THEN
         DO j = 1 , Nrhs
            DO i = 1 , N
               B(i,j) = S(i)*B(i,j)
            ENDDO
         ENDDO
      ENDIF
!
      IF ( nofact .OR. equil ) THEN
!
!        Compute the Cholesky factorization A = U**H * U or A = L * L**H.
!
         CALL ZCOPY(N*(N+1)/2,Ap,1,Afp,1)
         CALL ZPPTRF(Uplo,N,Afp,Info)
!
!        Return if INFO is non-zero.
!
         IF ( Info>0 ) THEN
            Rcond = ZERO
            RETURN
         ENDIF
      ENDIF
!
!     Compute the norm of the matrix A.
!
      anorm = ZLANHP('I',Uplo,N,Ap,Rwork)
!
!     Compute the reciprocal of the condition number of A.
!
      CALL ZPPCON(Uplo,N,Afp,anorm,Rcond,Work,Rwork,Info)
!
!     Compute the solution matrix X.
!
      CALL ZLACPY('Full',N,Nrhs,B,Ldb,X,Ldx)
      CALL ZPPTRS(Uplo,N,Nrhs,Afp,X,Ldx,Info)
!
!     Use iterative refinement to improve the computed solution and
!     compute error bounds and backward error estimates for it.
!
      CALL ZPPRFS(Uplo,N,Nrhs,Ap,Afp,B,Ldb,X,Ldx,Ferr,Berr,Work,Rwork,  &
     &            Info)
!
!     Transform the solution matrix X to a solution of the original
!     system.
!
      IF ( rcequ ) THEN
         DO j = 1 , Nrhs
            DO i = 1 , N
               X(i,j) = S(i)*X(i,j)
            ENDDO
         ENDDO
         DO j = 1 , Nrhs
            Ferr(j) = Ferr(j)/scond
         ENDDO
      ENDIF
!
!     Set INFO = N+1 if the matrix is singular to working precision.
!
      IF ( Rcond<DLAMCH('Epsilon') ) Info = N + 1
!
!
!     End of ZPPSVX
!
      END SUBROUTINE ZPPSVX
