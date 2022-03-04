!*==cget07.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b CGET07
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGET07( TRANS, N, NRHS, A, LDA, B, LDB, X, LDX, XACT,
!                          LDXACT, FERR, CHKFERR, BERR, RESLTS )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANS
!       LOGICAL            CHKFERR
!       INTEGER            LDA, LDB, LDX, LDXACT, N, NRHS
!       ..
!       .. Array Arguments ..
!       REAL               BERR( * ), FERR( * ), RESLTS( * )
!       COMPLEX            A( LDA, * ), B( LDB, * ), X( LDX, * ),
!      $                   XACT( LDXACT, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CGET07 tests the error bounds from iterative refinement for the
!> computed solution to a system of equations op(A)*X = B, where A is a
!> general n by n matrix and op(A) = A or A**T, depending on TRANS.
!>
!> RESLTS(1) = test of the error bound
!>           = norm(X - XACT) / ( norm(X) * FERR )
!>
!> A large value is returned if this ratio is not less than one.
!>
!> RESLTS(2) = residual from the iterative refinement routine
!>           = the maximum of BERR / ( (n+1)*EPS + (*) ), where
!>             (*) = (n+1)*UNFL / (min_i (abs(op(A))*abs(X) +abs(b))_i )
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          Specifies the form of the system of equations.
!>          = 'N':  A * X = B     (No transpose)
!>          = 'T':  A**T * X = B  (Transpose)
!>          = 'C':  A**H * X = B  (Conjugate transpose = Transpose)
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of rows of the matrices X and XACT.  N >= 0.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of columns of the matrices X and XACT.  NRHS >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          The original n by n matrix A.
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
!>          B is COMPLEX array, dimension (LDB,NRHS)
!>          The right hand side vectors for the system of linear
!>          equations.
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
!>          X is COMPLEX array, dimension (LDX,NRHS)
!>          The computed solution vectors.  Each vector is stored as a
!>          column of the matrix X.
!> \endverbatim
!>
!> \param[in] LDX
!> \verbatim
!>          LDX is INTEGER
!>          The leading dimension of the array X.  LDX >= max(1,N).
!> \endverbatim
!>
!> \param[in] XACT
!> \verbatim
!>          XACT is COMPLEX array, dimension (LDX,NRHS)
!>          The exact solution vectors.  Each vector is stored as a
!>          column of the matrix XACT.
!> \endverbatim
!>
!> \param[in] LDXACT
!> \verbatim
!>          LDXACT is INTEGER
!>          The leading dimension of the array XACT.  LDXACT >= max(1,N).
!> \endverbatim
!>
!> \param[in] FERR
!> \verbatim
!>          FERR is REAL array, dimension (NRHS)
!>          The estimated forward error bounds for each solution vector
!>          X.  If XTRUE is the true solution, FERR bounds the magnitude
!>          of the largest entry in (X - XTRUE) divided by the magnitude
!>          of the largest entry in X.
!> \endverbatim
!>
!> \param[in] CHKFERR
!> \verbatim
!>          CHKFERR is LOGICAL
!>          Set to .TRUE. to check FERR, .FALSE. not to check FERR.
!>          When the test system is ill-conditioned, the "true"
!>          solution in XACT may be incorrect.
!> \endverbatim
!>
!> \param[in] BERR
!> \verbatim
!>          BERR is REAL array, dimension (NRHS)
!>          The componentwise relative backward error of each solution
!>          vector (i.e., the smallest relative change in any entry of A
!>          or B that makes X an exact solution).
!> \endverbatim
!>
!> \param[out] RESLTS
!> \verbatim
!>          RESLTS is REAL array, dimension (2)
!>          The maximum over the NRHS solution vectors of the ratios:
!>          RESLTS(1) = norm(X - XACT) / ( norm(X) * FERR )
!>          RESLTS(2) = BERR / ( (n+1)*EPS + (*) )
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
!> \ingroup complex_lin
!
!  =====================================================================
      SUBROUTINE CGET07(Trans,N,Nrhs,A,Lda,B,Ldb,X,Ldx,Xact,Ldxact,Ferr,&
     &                  Chkferr,Berr,Reslts)
      IMPLICIT NONE
!*--CGET07170
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Trans
      LOGICAL Chkferr
      INTEGER Lda , Ldb , Ldx , Ldxact , N , Nrhs
!     ..
!     .. Array Arguments ..
      REAL Berr(*) , Ferr(*) , Reslts(*)
      COMPLEX A(Lda,*) , B(Ldb,*) , X(Ldx,*) , Xact(Ldxact,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
!     ..
!     .. Local Scalars ..
      LOGICAL notran
      INTEGER i , imax , j , k
      REAL axbi , diff , eps , errbnd , ovfl , tmp , unfl , xnorm
      COMPLEX zdum
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER ICAMAX
      REAL SLAMCH
      EXTERNAL LSAME , ICAMAX , SLAMCH
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , AIMAG , MAX , MIN , REAL
!     ..
!     .. Statement Functions ..
      REAL CABS1
!     ..
!     .. Statement Function definitions ..
      CABS1(zdum) = ABS(REAL(zdum)) + ABS(AIMAG(zdum))
!     ..
!     .. Executable Statements ..
!
!     Quick exit if N = 0 or NRHS = 0.
!
      IF ( N<=0 .OR. Nrhs<=0 ) THEN
         Reslts(1) = ZERO
         Reslts(2) = ZERO
         RETURN
      ENDIF
!
      eps = SLAMCH('Epsilon')
      unfl = SLAMCH('Safe minimum')
      ovfl = ONE/unfl
      notran = LSAME(Trans,'N')
!
!     Test 1:  Compute the maximum of
!        norm(X - XACT) / ( norm(X) * FERR )
!     over all the vectors X and XACT using the infinity-norm.
!
      errbnd = ZERO
      IF ( Chkferr ) THEN
         DO j = 1 , Nrhs
            imax = ICAMAX(N,X(1,j),1)
            xnorm = MAX(CABS1(X(imax,j)),unfl)
            diff = ZERO
            DO i = 1 , N
               diff = MAX(diff,CABS1(X(i,j)-Xact(i,j)))
            ENDDO
!
            IF ( xnorm>ONE ) THEN
            ELSEIF ( diff>ovfl*xnorm ) THEN
               errbnd = ONE/eps
               CYCLE
            ENDIF
!
            IF ( diff/xnorm<=Ferr(j) ) THEN
               errbnd = MAX(errbnd,(diff/xnorm)/Ferr(j))
            ELSE
               errbnd = ONE/eps
            ENDIF
         ENDDO
      ENDIF
      Reslts(1) = errbnd
!
!     Test 2:  Compute the maximum of BERR / ( (n+1)*EPS + (*) ), where
!     (*) = (n+1)*UNFL / (min_i (abs(op(A))*abs(X) +abs(b))_i )
!
      DO k = 1 , Nrhs
         DO i = 1 , N
            tmp = CABS1(B(i,k))
            IF ( notran ) THEN
               DO j = 1 , N
                  tmp = tmp + CABS1(A(i,j))*CABS1(X(j,k))
               ENDDO
            ELSE
               DO j = 1 , N
                  tmp = tmp + CABS1(A(j,i))*CABS1(X(j,k))
               ENDDO
            ENDIF
            IF ( i==1 ) THEN
               axbi = tmp
            ELSE
               axbi = MIN(axbi,tmp)
            ENDIF
         ENDDO
         tmp = Berr(k)/((N+1)*eps+(N+1)*unfl/MAX(axbi,(N+1)*unfl))
         IF ( k==1 ) THEN
            Reslts(2) = tmp
         ELSE
            Reslts(2) = MAX(Reslts(2),tmp)
         ENDIF
      ENDDO
!
!
!     End of CGET07
!
      END SUBROUTINE CGET07
