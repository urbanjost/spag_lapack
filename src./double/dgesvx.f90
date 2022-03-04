!*==dgesvx.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief <b> DGESVX computes the solution to system of linear equations A * X = B for GE matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DGESVX + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgesvx.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgesvx.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgesvx.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGESVX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV,
!                          EQUED, R, C, B, LDB, X, LDX, RCOND, FERR, BERR,
!                          WORK, IWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          EQUED, FACT, TRANS
!       INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS
!       DOUBLE PRECISION   RCOND
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * ), IWORK( * )
!       DOUBLE PRECISION   A( LDA, * ), AF( LDAF, * ), B( LDB, * ),
!      $                   BERR( * ), C( * ), FERR( * ), R( * ),
!      $                   WORK( * ), X( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGESVX uses the LU factorization to compute the solution to a real
!> system of linear equations
!>    A * X = B,
!> where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
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
!>       TRANS = 'N':  diag(R)*A*diag(C)     *inv(diag(C))*X = diag(R)*B
!>       TRANS = 'T': (diag(R)*A*diag(C))**T *inv(diag(R))*X = diag(C)*B
!>       TRANS = 'C': (diag(R)*A*diag(C))**H *inv(diag(R))*X = diag(C)*B
!>    Whether or not the system will be equilibrated depends on the
!>    scaling of the matrix A, but if equilibration is used, A is
!>    overwritten by diag(R)*A*diag(C) and B by diag(R)*B (if TRANS='N')
!>    or diag(C)*B (if TRANS = 'T' or 'C').
!>
!> 2. If FACT = 'N' or 'E', the LU decomposition is used to factor the
!>    matrix A (after equilibration if FACT = 'E') as
!>       A = P * L * U,
!>    where P is a permutation matrix, L is a unit lower triangular
!>    matrix, and U is upper triangular.
!>
!> 3. If some U(i,i)=0, so that U is exactly singular, then the routine
!>    returns with INFO = i. Otherwise, the factored form of A is used
!>    to estimate the condition number of the matrix A.  If the
!>    reciprocal of the condition number is less than machine precision,
!>    INFO = N+1 is returned as a warning, but the routine still goes on
!>    to solve for X and compute error bounds as described below.
!>
!> 4. The system of equations is solved for X using the factored form
!>    of A.
!>
!> 5. Iterative refinement is applied to improve the computed solution
!>    matrix and calculate error bounds and backward error estimates
!>    for it.
!>
!> 6. If equilibration was used, the matrix X is premultiplied by
!>    diag(C) (if TRANS = 'N') or diag(R) (if TRANS = 'T' or 'C') so
!>    that it solves the original system before equilibration.
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
!>          = 'F':  On entry, AF and IPIV contain the factored form of A.
!>                  If EQUED is not 'N', the matrix A has been
!>                  equilibrated with scaling factors given by R and C.
!>                  A, AF, and IPIV are not modified.
!>          = 'N':  The matrix A will be copied to AF and factored.
!>          = 'E':  The matrix A will be equilibrated if necessary, then
!>                  copied to AF and factored.
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          Specifies the form of the system of equations:
!>          = 'N':  A * X = B     (No transpose)
!>          = 'T':  A**T * X = B  (Transpose)
!>          = 'C':  A**H * X = B  (Transpose)
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
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, the N-by-N matrix A.  If FACT = 'F' and EQUED is
!>          not 'N', then A must have been equilibrated by the scaling
!>          factors in R and/or C.  A is not modified if FACT = 'F' or
!>          'N', or if FACT = 'E' and EQUED = 'N' on exit.
!>
!>          On exit, if EQUED .ne. 'N', A is scaled as follows:
!>          EQUED = 'R':  A := diag(R) * A
!>          EQUED = 'C':  A := A * diag(C)
!>          EQUED = 'B':  A := diag(R) * A * diag(C).
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] AF
!> \verbatim
!>          AF is DOUBLE PRECISION array, dimension (LDAF,N)
!>          If FACT = 'F', then AF is an input argument and on entry
!>          contains the factors L and U from the factorization
!>          A = P*L*U as computed by DGETRF.  If EQUED .ne. 'N', then
!>          AF is the factored form of the equilibrated matrix A.
!>
!>          If FACT = 'N', then AF is an output argument and on exit
!>          returns the factors L and U from the factorization A = P*L*U
!>          of the original matrix A.
!>
!>          If FACT = 'E', then AF is an output argument and on exit
!>          returns the factors L and U from the factorization A = P*L*U
!>          of the equilibrated matrix A (see the description of A for
!>          the form of the equilibrated matrix).
!> \endverbatim
!>
!> \param[in] LDAF
!> \verbatim
!>          LDAF is INTEGER
!>          The leading dimension of the array AF.  LDAF >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          If FACT = 'F', then IPIV is an input argument and on entry
!>          contains the pivot indices from the factorization A = P*L*U
!>          as computed by DGETRF; row i of the matrix was interchanged
!>          with row IPIV(i).
!>
!>          If FACT = 'N', then IPIV is an output argument and on exit
!>          contains the pivot indices from the factorization A = P*L*U
!>          of the original matrix A.
!>
!>          If FACT = 'E', then IPIV is an output argument and on exit
!>          contains the pivot indices from the factorization A = P*L*U
!>          of the equilibrated matrix A.
!> \endverbatim
!>
!> \param[in,out] EQUED
!> \verbatim
!>          EQUED is CHARACTER*1
!>          Specifies the form of equilibration that was done.
!>          = 'N':  No equilibration (always true if FACT = 'N').
!>          = 'R':  Row equilibration, i.e., A has been premultiplied by
!>                  diag(R).
!>          = 'C':  Column equilibration, i.e., A has been postmultiplied
!>                  by diag(C).
!>          = 'B':  Both row and column equilibration, i.e., A has been
!>                  replaced by diag(R) * A * diag(C).
!>          EQUED is an input argument if FACT = 'F'; otherwise, it is an
!>          output argument.
!> \endverbatim
!>
!> \param[in,out] R
!> \verbatim
!>          R is DOUBLE PRECISION array, dimension (N)
!>          The row scale factors for A.  If EQUED = 'R' or 'B', A is
!>          multiplied on the left by diag(R); if EQUED = 'N' or 'C', R
!>          is not accessed.  R is an input argument if FACT = 'F';
!>          otherwise, R is an output argument.  If FACT = 'F' and
!>          EQUED = 'R' or 'B', each element of R must be positive.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is DOUBLE PRECISION array, dimension (N)
!>          The column scale factors for A.  If EQUED = 'C' or 'B', A is
!>          multiplied on the right by diag(C); if EQUED = 'N' or 'R', C
!>          is not accessed.  C is an input argument if FACT = 'F';
!>          otherwise, C is an output argument.  If FACT = 'F' and
!>          EQUED = 'C' or 'B', each element of C must be positive.
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)
!>          On entry, the N-by-NRHS right hand side matrix B.
!>          On exit,
!>          if EQUED = 'N', B is not modified;
!>          if TRANS = 'N' and EQUED = 'R' or 'B', B is overwritten by
!>          diag(R)*B;
!>          if TRANS = 'T' or 'C' and EQUED = 'C' or 'B', B is
!>          overwritten by diag(C)*B.
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
!>          X is DOUBLE PRECISION array, dimension (LDX,NRHS)
!>          If INFO = 0 or INFO = N+1, the N-by-NRHS solution matrix X
!>          to the original system of equations.  Note that A and B are
!>          modified on exit if EQUED .ne. 'N', and the solution to the
!>          equilibrated system is inv(diag(C))*X if TRANS = 'N' and
!>          EQUED = 'C' or 'B', or inv(diag(R))*X if TRANS = 'T' or 'C'
!>          and EQUED = 'R' or 'B'.
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
!>          WORK is DOUBLE PRECISION array, dimension (4*N)
!>          On exit, WORK(1) contains the reciprocal pivot growth
!>          factor norm(A)/norm(U). The "max absolute element" norm is
!>          used. If WORK(1) is much less than 1, then the stability
!>          of the LU factorization of the (equilibrated) matrix A
!>          could be poor. This also means that the solution X, condition
!>          estimator RCOND, and forward error bound FERR could be
!>          unreliable. If factorization fails with 0<INFO<=N, then
!>          WORK(1) contains the reciprocal pivot growth factor for the
!>          leading INFO columns of A.
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
!>          > 0:  if INFO = i, and i is
!>                <= N:  U(i,i) is exactly zero.  The factorization has
!>                       been completed, but the factor U is exactly
!>                       singular, so the solution and error bounds
!>                       could not be computed. RCOND = 0 is returned.
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
!> \ingroup doubleGEsolve
!
!  =====================================================================
      SUBROUTINE DGESVX(Fact,Trans,N,Nrhs,A,Lda,Af,Ldaf,Ipiv,Equed,R,C, &
     &                  B,Ldb,X,Ldx,Rcond,Ferr,Berr,Work,Iwork,Info)
      USE F77KINDS                        
      USE S_DGECON
      USE S_DGEEQU
      USE S_DGERFS
      USE S_DGETRF
      USE S_DGETRS
      USE S_DLACPY
      USE S_DLAMCH
      USE S_DLANGE
      USE S_DLANTR
      USE S_DLAQGE
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--DGESVX365
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Fact
      CHARACTER :: Trans
      INTEGER :: N
      INTEGER :: Nrhs
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      CHARACTER :: Equed
      REAL(R8KIND) , DIMENSION(*) :: R
      REAL(R8KIND) , DIMENSION(*) :: C
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL(R8KIND) , INTENT(INOUT) :: Rcond
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL(R8KIND) , DIMENSION(*) :: Berr
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: amax , anorm , bignum , colcnd , rcmax , rcmin ,  &
     &                rowcnd , rpvgrw , smlnum
      LOGICAL :: colequ , equil , nofact , notran , rowequ
      INTEGER :: i , infequ , j
      CHARACTER :: norm
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
      notran = LSAME(Trans,'N')
      IF ( nofact .OR. equil ) THEN
         Equed = 'N'
         rowequ = .FALSE.
         colequ = .FALSE.
      ELSE
         rowequ = LSAME(Equed,'R') .OR. LSAME(Equed,'B')
         colequ = LSAME(Equed,'C') .OR. LSAME(Equed,'B')
         smlnum = DLAMCH('Safe minimum')
         bignum = ONE/smlnum
      ENDIF
!
!     Test the input parameters.
!
      IF ( .NOT.nofact .AND. .NOT.equil .AND. .NOT.LSAME(Fact,'F') )    &
     &     THEN
         Info = -1
      ELSEIF ( .NOT.notran .AND. .NOT.LSAME(Trans,'T') .AND.            &
     &         .NOT.LSAME(Trans,'C') ) THEN
         Info = -2
      ELSEIF ( N<0 ) THEN
         Info = -3
      ELSEIF ( Nrhs<0 ) THEN
         Info = -4
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -6
      ELSEIF ( Ldaf<MAX(1,N) ) THEN
         Info = -8
      ELSEIF ( LSAME(Fact,'F') .AND.                                    &
     &         .NOT.(rowequ .OR. colequ .OR. LSAME(Equed,'N')) ) THEN
         Info = -10
      ELSE
         IF ( rowequ ) THEN
            rcmin = bignum
            rcmax = ZERO
            DO j = 1 , N
               rcmin = MIN(rcmin,R(j))
               rcmax = MAX(rcmax,R(j))
            ENDDO
            IF ( rcmin<=ZERO ) THEN
               Info = -11
            ELSEIF ( N>0 ) THEN
               rowcnd = MAX(rcmin,smlnum)/MIN(rcmax,bignum)
            ELSE
               rowcnd = ONE
            ENDIF
         ENDIF
         IF ( colequ .AND. Info==0 ) THEN
            rcmin = bignum
            rcmax = ZERO
            DO j = 1 , N
               rcmin = MIN(rcmin,C(j))
               rcmax = MAX(rcmax,C(j))
            ENDDO
            IF ( rcmin<=ZERO ) THEN
               Info = -12
            ELSEIF ( N>0 ) THEN
               colcnd = MAX(rcmin,smlnum)/MIN(rcmax,bignum)
            ELSE
               colcnd = ONE
            ENDIF
         ENDIF
         IF ( Info==0 ) THEN
            IF ( Ldb<MAX(1,N) ) THEN
               Info = -14
            ELSEIF ( Ldx<MAX(1,N) ) THEN
               Info = -16
            ENDIF
         ENDIF
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('DGESVX',-Info)
         RETURN
      ENDIF
!
      IF ( equil ) THEN
!
!        Compute row and column scalings to equilibrate the matrix A.
!
         CALL DGEEQU(N,N,A,Lda,R,C,rowcnd,colcnd,amax,infequ)
         IF ( infequ==0 ) THEN
!
!           Equilibrate the matrix.
!
            CALL DLAQGE(N,N,A,Lda,R,C,rowcnd,colcnd,amax,Equed)
            rowequ = LSAME(Equed,'R') .OR. LSAME(Equed,'B')
            colequ = LSAME(Equed,'C') .OR. LSAME(Equed,'B')
         ENDIF
      ENDIF
!
!     Scale the right hand side.
!
      IF ( notran ) THEN
         IF ( rowequ ) THEN
            DO j = 1 , Nrhs
               DO i = 1 , N
                  B(i,j) = R(i)*B(i,j)
               ENDDO
            ENDDO
         ENDIF
      ELSEIF ( colequ ) THEN
         DO j = 1 , Nrhs
            DO i = 1 , N
               B(i,j) = C(i)*B(i,j)
            ENDDO
         ENDDO
      ENDIF
!
      IF ( nofact .OR. equil ) THEN
!
!        Compute the LU factorization of A.
!
         CALL DLACPY('Full',N,N,A,Lda,Af,Ldaf)
         CALL DGETRF(N,N,Af,Ldaf,Ipiv,Info)
!
!        Return if INFO is non-zero.
!
         IF ( Info>0 ) THEN
!
!           Compute the reciprocal pivot growth factor of the
!           leading rank-deficient INFO columns of A.
!
            rpvgrw = DLANTR('M','U','N',Info,Info,Af,Ldaf,Work)
            IF ( rpvgrw==ZERO ) THEN
               rpvgrw = ONE
            ELSE
               rpvgrw = DLANGE('M',N,Info,A,Lda,Work)/rpvgrw
            ENDIF
            Work(1) = rpvgrw
            Rcond = ZERO
            RETURN
         ENDIF
      ENDIF
!
!     Compute the norm of the matrix A and the
!     reciprocal pivot growth factor RPVGRW.
!
      IF ( notran ) THEN
         norm = '1'
      ELSE
         norm = 'I'
      ENDIF
      anorm = DLANGE(norm,N,N,A,Lda,Work)
      rpvgrw = DLANTR('M','U','N',N,N,Af,Ldaf,Work)
      IF ( rpvgrw==ZERO ) THEN
         rpvgrw = ONE
      ELSE
         rpvgrw = DLANGE('M',N,N,A,Lda,Work)/rpvgrw
      ENDIF
!
!     Compute the reciprocal of the condition number of A.
!
      CALL DGECON(norm,N,Af,Ldaf,anorm,Rcond,Work,Iwork,Info)
!
!     Compute the solution matrix X.
!
      CALL DLACPY('Full',N,Nrhs,B,Ldb,X,Ldx)
      CALL DGETRS(Trans,N,Nrhs,Af,Ldaf,Ipiv,X,Ldx,Info)
!
!     Use iterative refinement to improve the computed solution and
!     compute error bounds and backward error estimates for it.
!
      CALL DGERFS(Trans,N,Nrhs,A,Lda,Af,Ldaf,Ipiv,B,Ldb,X,Ldx,Ferr,Berr,&
     &            Work,Iwork,Info)
!
!     Transform the solution matrix X to a solution of the original
!     system.
!
      IF ( notran ) THEN
         IF ( colequ ) THEN
            DO j = 1 , Nrhs
               DO i = 1 , N
                  X(i,j) = C(i)*X(i,j)
               ENDDO
            ENDDO
            DO j = 1 , Nrhs
               Ferr(j) = Ferr(j)/colcnd
            ENDDO
         ENDIF
      ELSEIF ( rowequ ) THEN
         DO j = 1 , Nrhs
            DO i = 1 , N
               X(i,j) = R(i)*X(i,j)
            ENDDO
         ENDDO
         DO j = 1 , Nrhs
            Ferr(j) = Ferr(j)/rowcnd
         ENDDO
      ENDIF
!
      Work(1) = rpvgrw
!
!     Set INFO = N+1 if the matrix is singular to working precision.
!
      IF ( Rcond<DLAMCH('Epsilon') ) Info = N + 1
!
!     End of DGESVX
!
      END SUBROUTINE DGESVX
