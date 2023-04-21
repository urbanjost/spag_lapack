!*==dgbsvx.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief <b> DGBSVX computes the solution to system of linear equations A * X = B for GB matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DGBSVX + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgbsvx.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgbsvx.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgbsvx.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGBSVX( FACT, TRANS, N, KL, KU, NRHS, AB, LDAB, AFB,
!                          LDAFB, IPIV, EQUED, R, C, B, LDB, X, LDX,
!                          RCOND, FERR, BERR, WORK, IWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          EQUED, FACT, TRANS
!       INTEGER            INFO, KL, KU, LDAB, LDAFB, LDB, LDX, N, NRHS
!       DOUBLE PRECISION   RCOND
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * ), IWORK( * )
!       DOUBLE PRECISION   AB( LDAB, * ), AFB( LDAFB, * ), B( LDB, * ),
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
!> DGBSVX uses the LU factorization to compute the solution to a real
!> system of linear equations A * X = B, A**T * X = B, or A**H * X = B,
!> where A is a band matrix of order N with KL subdiagonals and KU
!> superdiagonals, and X and B are N-by-NRHS matrices.
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
!> The following steps are performed by this subroutine:
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
!>       A = L * U,
!>    where L is a product of permutation and unit lower triangular
!>    matrices with KL subdiagonals, and U is upper triangular with
!>    KL+KU superdiagonals.
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
!>          = 'F':  On entry, AFB and IPIV contain the factored form of
!>                  A.  If EQUED is not 'N', the matrix A has been
!>                  equilibrated with scaling factors given by R and C.
!>                  AB, AFB, and IPIV are not modified.
!>          = 'N':  The matrix A will be copied to AFB and factored.
!>          = 'E':  The matrix A will be equilibrated if necessary, then
!>                  copied to AFB and factored.
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          Specifies the form of the system of equations.
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
!> \param[in,out] AB
!> \verbatim
!>          AB is DOUBLE PRECISION array, dimension (LDAB,N)
!>          On entry, the matrix A in band storage, in rows 1 to KL+KU+1.
!>          The j-th column of A is stored in the j-th column of the
!>          array AB as follows:
!>          AB(KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+kl)
!>
!>          If FACT = 'F' and EQUED is not 'N', then A must have been
!>          equilibrated by the scaling factors in R and/or C.  AB is not
!>          modified if FACT = 'F' or 'N', or if FACT = 'E' and
!>          EQUED = 'N' on exit.
!>
!>          On exit, if EQUED .ne. 'N', A is scaled as follows:
!>          EQUED = 'R':  A := diag(R) * A
!>          EQUED = 'C':  A := A * diag(C)
!>          EQUED = 'B':  A := diag(R) * A * diag(C).
!> \endverbatim
!>
!> \param[in] LDAB
!> \verbatim
!>          LDAB is INTEGER
!>          The leading dimension of the array AB.  LDAB >= KL+KU+1.
!> \endverbatim
!>
!> \param[in,out] AFB
!> \verbatim
!>          AFB is DOUBLE PRECISION array, dimension (LDAFB,N)
!>          If FACT = 'F', then AFB is an input argument and on entry
!>          contains details of the LU factorization of the band matrix
!>          A, as computed by DGBTRF.  U is stored as an upper triangular
!>          band matrix with KL+KU superdiagonals in rows 1 to KL+KU+1,
!>          and the multipliers used during the factorization are stored
!>          in rows KL+KU+2 to 2*KL+KU+1.  If EQUED .ne. 'N', then AFB is
!>          the factored form of the equilibrated matrix A.
!>
!>          If FACT = 'N', then AFB is an output argument and on exit
!>          returns details of the LU factorization of A.
!>
!>          If FACT = 'E', then AFB is an output argument and on exit
!>          returns details of the LU factorization of the equilibrated
!>          matrix A (see the description of AB for the form of the
!>          equilibrated matrix).
!> \endverbatim
!>
!> \param[in] LDAFB
!> \verbatim
!>          LDAFB is INTEGER
!>          The leading dimension of the array AFB.  LDAFB >= 2*KL+KU+1.
!> \endverbatim
!>
!> \param[in,out] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          If FACT = 'F', then IPIV is an input argument and on entry
!>          contains the pivot indices from the factorization A = L*U
!>          as computed by DGBTRF; row i of the matrix was interchanged
!>          with row IPIV(i).
!>
!>          If FACT = 'N', then IPIV is an output argument and on exit
!>          contains the pivot indices from the factorization A = L*U
!>          of the original matrix A.
!>
!>          If FACT = 'E', then IPIV is an output argument and on exit
!>          contains the pivot indices from the factorization A = L*U
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
!>          On entry, the right hand side matrix B.
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
!>          WORK is DOUBLE PRECISION array, dimension (3*N)
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
!>                <= N:  U(i,i) is exactly zero.  The factorization
!>                       has been completed, but the factor U is exactly
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
!> \ingroup doubleGBsolve
!
!  =====================================================================
      SUBROUTINE DGBSVX(Fact,Trans,N,Kl,Ku,Nrhs,Ab,Ldab,Afb,Ldafb,Ipiv, &
     &                  Equed,R,C,B,Ldb,X,Ldx,Rcond,Ferr,Berr,Work,     &
     &                  Iwork,Info)
      IMPLICIT NONE
!*--DGBSVX373
!
!  -- LAPACK driver routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     April 2012
!
!     .. Scalar Arguments ..
      CHARACTER Equed , Fact , Trans
      INTEGER Info , Kl , Ku , Ldab , Ldafb , Ldb , Ldx , N , Nrhs
      DOUBLE PRECISION Rcond
!     ..
!     .. Array Arguments ..
      INTEGER Ipiv(*) , Iwork(*)
      DOUBLE PRECISION Ab(Ldab,*) , Afb(Ldafb,*) , B(Ldb,*) , Berr(*) , &
     &                 C(*) , Ferr(*) , R(*) , Work(*) , X(Ldx,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
!     ..
!     .. Local Scalars ..
      LOGICAL colequ , equil , nofact , notran , rowequ
      CHARACTER norm
      INTEGER i , infequ , j , j1 , j2
      DOUBLE PRECISION amax , anorm , bignum , colcnd , rcmax , rcmin , &
     &                 rowcnd , rpvgrw , smlnum
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      DOUBLE PRECISION DLAMCH , DLANGB , DLANTB
      EXTERNAL LSAME , DLAMCH , DLANGB , DLANTB
!     ..
!     .. External Subroutines ..
      EXTERNAL DCOPY , DGBCON , DGBEQU , DGBRFS , DGBTRF , DGBTRS ,     &
     &         DLACPY , DLAQGB , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , MIN
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
      ELSEIF ( Kl<0 ) THEN
         Info = -4
      ELSEIF ( Ku<0 ) THEN
         Info = -5
      ELSEIF ( Nrhs<0 ) THEN
         Info = -6
      ELSEIF ( Ldab<Kl+Ku+1 ) THEN
         Info = -8
      ELSEIF ( Ldafb<2*Kl+Ku+1 ) THEN
         Info = -10
      ELSEIF ( LSAME(Fact,'F') .AND.                                    &
     &         .NOT.(rowequ .OR. colequ .OR. LSAME(Equed,'N')) ) THEN
         Info = -12
      ELSE
         IF ( rowequ ) THEN
            rcmin = bignum
            rcmax = ZERO
            DO j = 1 , N
               rcmin = MIN(rcmin,R(j))
               rcmax = MAX(rcmax,R(j))
            ENDDO
            IF ( rcmin<=ZERO ) THEN
               Info = -13
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
               Info = -14
            ELSEIF ( N>0 ) THEN
               colcnd = MAX(rcmin,smlnum)/MIN(rcmax,bignum)
            ELSE
               colcnd = ONE
            ENDIF
         ENDIF
         IF ( Info==0 ) THEN
            IF ( Ldb<MAX(1,N) ) THEN
               Info = -16
            ELSEIF ( Ldx<MAX(1,N) ) THEN
               Info = -18
            ENDIF
         ENDIF
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('DGBSVX',-Info)
         RETURN
      ENDIF
!
      IF ( equil ) THEN
!
!        Compute row and column scalings to equilibrate the matrix A.
!
         CALL DGBEQU(N,N,Kl,Ku,Ab,Ldab,R,C,rowcnd,colcnd,amax,infequ)
         IF ( infequ==0 ) THEN
!
!           Equilibrate the matrix.
!
            CALL DLAQGB(N,N,Kl,Ku,Ab,Ldab,R,C,rowcnd,colcnd,amax,Equed)
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
!        Compute the LU factorization of the band matrix A.
!
         DO j = 1 , N
            j1 = MAX(j-Ku,1)
            j2 = MIN(j+Kl,N)
            CALL DCOPY(j2-j1+1,Ab(Ku+1-j+j1,j),1,Afb(Kl+Ku+1-j+j1,j),1)
         ENDDO
!
         CALL DGBTRF(N,N,Kl,Ku,Afb,Ldafb,Ipiv,Info)
!
!        Return if INFO is non-zero.
!
         IF ( Info>0 ) THEN
!
!           Compute the reciprocal pivot growth factor of the
!           leading rank-deficient INFO columns of A.
!
            anorm = ZERO
            DO j = 1 , Info
               DO i = MAX(Ku+2-j,1) , MIN(N+Ku+1-j,Kl+Ku+1)
                  anorm = MAX(anorm,ABS(Ab(i,j)))
               ENDDO
            ENDDO
            rpvgrw = DLANTB('M','U','N',Info,MIN(Info-1,Kl+Ku),         &
     &               Afb(MAX(1,Kl+Ku+2-Info),1),Ldafb,Work)
            IF ( rpvgrw==ZERO ) THEN
               rpvgrw = ONE
            ELSE
               rpvgrw = anorm/rpvgrw
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
      anorm = DLANGB(norm,N,Kl,Ku,Ab,Ldab,Work)
      rpvgrw = DLANTB('M','U','N',N,Kl+Ku,Afb,Ldafb,Work)
      IF ( rpvgrw==ZERO ) THEN
         rpvgrw = ONE
      ELSE
         rpvgrw = DLANGB('M',N,Kl,Ku,Ab,Ldab,Work)/rpvgrw
      ENDIF
!
!     Compute the reciprocal of the condition number of A.
!
      CALL DGBCON(norm,N,Kl,Ku,Afb,Ldafb,Ipiv,anorm,Rcond,Work,Iwork,   &
     &            Info)
!
!     Compute the solution matrix X.
!
      CALL DLACPY('Full',N,Nrhs,B,Ldb,X,Ldx)
      CALL DGBTRS(Trans,N,Kl,Ku,Nrhs,Afb,Ldafb,Ipiv,X,Ldx,Info)
!
!     Use iterative refinement to improve the computed solution and
!     compute error bounds and backward error estimates for it.
!
      CALL DGBRFS(Trans,N,Kl,Ku,Nrhs,Ab,Ldab,Afb,Ldafb,Ipiv,B,Ldb,X,Ldx,&
     &            Ferr,Berr,Work,Iwork,Info)
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
!     Set INFO = N+1 if the matrix is singular to working precision.
!
      IF ( Rcond<DLAMCH('Epsilon') ) Info = N + 1
!
      Work(1) = rpvgrw
!
!     End of DGBSVX
!
      END SUBROUTINE DGBSVX
