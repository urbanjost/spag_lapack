!*==dsgesv.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief <b> DSGESV computes the solution to system of linear equations A * X = B for GE matrices</b> (mixed precision with iterative refinement)
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DSGESV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsgesv.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsgesv.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsgesv.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSGESV( N, NRHS, A, LDA, IPIV, B, LDB, X, LDX, WORK,
!                          SWORK, ITER, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, ITER, LDA, LDB, LDX, N, NRHS
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       REAL               SWORK( * )
!       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), WORK( N, * ),
!      $                   X( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DSGESV computes the solution to a real system of linear equations
!>    A * X = B,
!> where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
!>
!> DSGESV first attempts to factorize the matrix in SINGLE PRECISION
!> and use this factorization within an iterative refinement procedure
!> to produce a solution with DOUBLE PRECISION normwise backward error
!> quality (see below). If the approach fails the method switches to a
!> DOUBLE PRECISION factorization and solve.
!>
!> The iterative refinement is not going to be a winning strategy if
!> the ratio SINGLE PRECISION performance over DOUBLE PRECISION
!> performance is too small. A reasonable strategy should take the
!> number of right-hand sides and the size of the matrix into account.
!> This might be done with a call to ILAENV in the future. Up to now, we
!> always try iterative refinement.
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
!>          A is DOUBLE PRECISION array,
!>          dimension (LDA,N)
!>          On entry, the N-by-N coefficient matrix A.
!>          On exit, if iterative refinement has been successfully used
!>          (INFO = 0 and ITER >= 0, see description below), then A is
!>          unchanged, if double precision factorization has been used
!>          (INFO = 0 and ITER < 0, see description below), then the
!>          array A contains the factors L and U from the factorization
!>          A = P*L*U; the unit diagonal elements of L are not stored.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          The pivot indices that define the permutation matrix P;
!>          row i of the matrix was interchanged with row IPIV(i).
!>          Corresponds either to the single precision factorization
!>          (if INFO = 0 and ITER >= 0) or the double precision
!>          factorization (if INFO = 0 and ITER < 0).
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)
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
!>          X is DOUBLE PRECISION array, dimension (LDX,NRHS)
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
!>          WORK is DOUBLE PRECISION array, dimension (N,NRHS)
!>          This array is used to hold the residual vectors.
!> \endverbatim
!>
!> \param[out] SWORK
!> \verbatim
!>          SWORK is REAL array, dimension (N*(N+NRHS))
!>          This array is used to use the single precision matrix and the
!>          right-hand sides or solutions in single precision.
!> \endverbatim
!>
!> \param[out] ITER
!> \verbatim
!>          ITER is INTEGER
!>          < 0: iterative refinement has failed, double precision
!>               factorization has been performed
!>               -1 : the routine fell back to full precision for
!>                    implementation- or machine-specific reasons
!>               -2 : narrowing the precision induced an overflow,
!>                    the routine fell back to full precision
!>               -3 : failure of SGETRF
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
!>          > 0:  if INFO = i, U(i,i) computed in DOUBLE PRECISION is
!>                exactly zero.  The factorization has been completed,
!>                but the factor U is exactly singular, so the solution
!>                could not be computed.
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
!> \ingroup doubleGEsolve
!
!  =====================================================================
      SUBROUTINE DSGESV(N,Nrhs,A,Lda,Ipiv,B,Ldb,X,Ldx,Work,Swork,Iter,  &
     &                  Info)
      IMPLICIT NONE
!*--DSGESV199
!
!  -- LAPACK driver routine (version 3.8.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!
!     .. Scalar Arguments ..
      INTEGER Info , Iter , Lda , Ldb , Ldx , N , Nrhs
!     ..
!     .. Array Arguments ..
      INTEGER Ipiv(*)
      REAL Swork(*)
      DOUBLE PRECISION A(Lda,*) , B(Ldb,*) , Work(N,*) , X(Ldx,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      LOGICAL DOITREF
      PARAMETER (DOITREF=.TRUE.)
!
      INTEGER ITERMAX
      PARAMETER (ITERMAX=30)
!
      DOUBLE PRECISION BWDMAX
      PARAMETER (BWDMAX=1.0E+00)
!
      DOUBLE PRECISION NEGONE , ONE
      PARAMETER (NEGONE=-1.0D+0,ONE=1.0D+0)
!
!     .. Local Scalars ..
      INTEGER i , iiter , ptsa , ptsx
      DOUBLE PRECISION anrm , cte , eps , rnrm , xnrm
!
!     .. External Subroutines ..
      EXTERNAL DAXPY , DGEMM , DLACPY , DLAG2S , DGETRF , DGETRS ,      &
     &         SGETRF , SGETRS , SLAG2D , XERBLA
!     ..
!     .. External Functions ..
      INTEGER IDAMAX
      DOUBLE PRECISION DLAMCH , DLANGE
      EXTERNAL IDAMAX , DLAMCH , DLANGE
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , DBLE , MAX , SQRT
!     ..
!     .. Executable Statements ..
!
      Info = 0
      Iter = 0
!
!     Test the input parameters.
!
      IF ( N<0 ) THEN
         Info = -1
      ELSEIF ( Nrhs<0 ) THEN
         Info = -2
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -4
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -7
      ELSEIF ( Ldx<MAX(1,N) ) THEN
         Info = -9
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('DSGESV',-Info)
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
      anrm = DLANGE('I',N,N,A,Lda,Work)
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
      CALL DLAG2S(N,Nrhs,B,Ldb,Swork(ptsx),N,Info)
!
      IF ( Info/=0 ) THEN
         Iter = -2
         GOTO 300
      ENDIF
!
!     Convert A from double precision to single precision and store the
!     result in SA.
!
      CALL DLAG2S(N,N,A,Lda,Swork(ptsa),N,Info)
!
      IF ( Info/=0 ) THEN
         Iter = -2
         GOTO 300
      ENDIF
!
!     Compute the LU factorization of SA.
!
      CALL SGETRF(N,N,Swork(ptsa),N,Ipiv,Info)
!
      IF ( Info/=0 ) THEN
         Iter = -3
         GOTO 300
      ENDIF
!
!     Solve the system SA*SX = SB.
!
      CALL SGETRS('No transpose',N,Nrhs,Swork(ptsa),N,Ipiv,Swork(ptsx), &
     &            N,Info)
!
!     Convert SX back to double precision
!
      CALL SLAG2D(N,Nrhs,Swork(ptsx),N,X,Ldx,Info)
!
!     Compute R = B - AX (R is WORK).
!
      CALL DLACPY('All',N,Nrhs,B,Ldb,Work,N)
!
      CALL DGEMM('No Transpose','No Transpose',N,Nrhs,N,NEGONE,A,Lda,X, &
     &           Ldx,ONE,Work,N)
!
!     Check whether the NRHS normwise backward errors satisfy the
!     stopping criterion. If yes, set ITER=0 and return.
!
      DO i = 1 , Nrhs
         xnrm = ABS(X(IDAMAX(N,X(1,i),1),i))
         rnrm = ABS(Work(IDAMAX(N,Work(1,i),1),i))
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
         CALL DLAG2S(N,Nrhs,Work,N,Swork(ptsx),N,Info)
!
         IF ( Info/=0 ) THEN
            Iter = -2
            GOTO 300
         ENDIF
!
!        Solve the system SA*SX = SR.
!
         CALL SGETRS('No transpose',N,Nrhs,Swork(ptsa),N,Ipiv,          &
     &               Swork(ptsx),N,Info)
!
!        Convert SX back to double precision and update the current
!        iterate.
!
         CALL SLAG2D(N,Nrhs,Swork(ptsx),N,Work,N,Info)
!
         DO i = 1 , Nrhs
            CALL DAXPY(N,ONE,Work(1,i),1,X(1,i),1)
         ENDDO
!
!        Compute R = B - AX (R is WORK).
!
         CALL DLACPY('All',N,Nrhs,B,Ldb,Work,N)
!
         CALL DGEMM('No Transpose','No Transpose',N,Nrhs,N,NEGONE,A,Lda,&
     &              X,Ldx,ONE,Work,N)
!
!        Check whether the NRHS normwise backward errors satisfy the
!        stopping criterion. If yes, set ITER=IITER>0 and return.
!
         DO i = 1 , Nrhs
            xnrm = ABS(X(IDAMAX(N,X(1,i),1),i))
            rnrm = ABS(Work(IDAMAX(N,Work(1,i),1),i))
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
!     stopping criterion, set up the ITER flag accordingly and follow up
!     on double precision routine.
!
      Iter = -ITERMAX - 1
!
!
!     Single-precision iterative refinement failed to converge to a
!     satisfactory solution, so we resort to double precision.
!
 300  CALL DGETRF(N,N,A,Lda,Ipiv,Info)
!
      IF ( Info/=0 ) RETURN
!
      CALL DLACPY('All',N,Nrhs,B,Ldb,X,Ldx)
      CALL DGETRS('No transpose',N,Nrhs,A,Lda,Ipiv,X,Ldx,Info)
!
!
!     End of DSGESV.
!
      END SUBROUTINE DSGESV
