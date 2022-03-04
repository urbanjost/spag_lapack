!*==stpt05.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b STPT05
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE STPT05( UPLO, TRANS, DIAG, N, NRHS, AP, B, LDB, X, LDX,
!                          XACT, LDXACT, FERR, BERR, RESLTS )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIAG, TRANS, UPLO
!       INTEGER            LDB, LDX, LDXACT, N, NRHS
!       ..
!       .. Array Arguments ..
!       REAL               AP( * ), B( LDB, * ), BERR( * ), FERR( * ),
!      $                   RESLTS( * ), X( LDX, * ), XACT( LDXACT, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> STPT05 tests the error bounds from iterative refinement for the
!> computed solution to a system of equations A*X = B, where A is a
!> triangular matrix in packed storage format.
!>
!> RESLTS(1) = test of the error bound
!>           = norm(X - XACT) / ( norm(X) * FERR )
!>
!> A large value is returned if this ratio is not less than one.
!>
!> RESLTS(2) = residual from the iterative refinement routine
!>           = the maximum of BERR / ( (n+1)*EPS + (*) ), where
!>             (*) = (n+1)*UNFL / (min_i (abs(A)*abs(X) +abs(b))_i )
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the matrix A is upper or lower triangular.
!>          = 'U':  Upper triangular
!>          = 'L':  Lower triangular
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          Specifies the form of the system of equations.
!>          = 'N':  A * X = B  (No transpose)
!>          = 'T':  A'* X = B  (Transpose)
!>          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
!> \endverbatim
!>
!> \param[in] DIAG
!> \verbatim
!>          DIAG is CHARACTER*1
!>          Specifies whether or not the matrix A is unit triangular.
!>          = 'N':  Non-unit triangular
!>          = 'U':  Unit triangular
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of rows of the matrices X, B, and XACT, and the
!>          order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of columns of the matrices X, B, and XACT.
!>          NRHS >= 0.
!> \endverbatim
!>
!> \param[in] AP
!> \verbatim
!>          AP is REAL array, dimension (N*(N+1)/2)
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
!>          B is REAL array, dimension (LDB,NRHS)
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
!>          X is REAL array, dimension (LDX,NRHS)
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
!>          XACT is REAL array, dimension (LDX,NRHS)
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
!> \ingroup single_lin
!
!  =====================================================================
      SUBROUTINE STPT05(Uplo,Trans,Diag,N,Nrhs,Ap,B,Ldb,X,Ldx,Xact,     &
     &                  Ldxact,Ferr,Berr,Reslts)
      IMPLICIT NONE
!*--STPT05178
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Diag , Trans , Uplo
      INTEGER Ldb , Ldx , Ldxact , N , Nrhs
!     ..
!     .. Array Arguments ..
      REAL Ap(*) , B(Ldb,*) , Berr(*) , Ferr(*) , Reslts(*) , X(Ldx,*) ,&
     &     Xact(Ldxact,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
!     ..
!     .. Local Scalars ..
      LOGICAL notran , unit , upper
      INTEGER i , ifu , imax , j , jc , k
      REAL axbi , diff , eps , errbnd , ovfl , tmp , unfl , xnorm
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER ISAMAX
      REAL SLAMCH
      EXTERNAL LSAME , ISAMAX , SLAMCH
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , MIN
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
      upper = LSAME(Uplo,'U')
      notran = LSAME(Trans,'N')
      unit = LSAME(Diag,'U')
!
!     Test 1:  Compute the maximum of
!        norm(X - XACT) / ( norm(X) * FERR )
!     over all the vectors X and XACT using the infinity-norm.
!
      errbnd = ZERO
      DO j = 1 , Nrhs
         imax = ISAMAX(N,X(1,j),1)
         xnorm = MAX(ABS(X(imax,j)),unfl)
         diff = ZERO
         DO i = 1 , N
            diff = MAX(diff,ABS(X(i,j)-Xact(i,j)))
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
      Reslts(1) = errbnd
!
!     Test 2:  Compute the maximum of BERR / ( (n+1)*EPS + (*) ), where
!     (*) = (n+1)*UNFL / (min_i (abs(A)*abs(X) +abs(b))_i )
!
      ifu = 0
      IF ( unit ) ifu = 1
      DO k = 1 , Nrhs
         DO i = 1 , N
            tmp = ABS(B(i,k))
            IF ( upper ) THEN
               jc = ((i-1)*i)/2
               IF ( .NOT.notran ) THEN
                  DO j = 1 , i - ifu
                     tmp = tmp + ABS(Ap(jc+j))*ABS(X(j,k))
                  ENDDO
                  IF ( unit ) tmp = tmp + ABS(X(i,k))
               ELSE
                  jc = jc + i
                  IF ( unit ) THEN
                     tmp = tmp + ABS(X(i,k))
                     jc = jc + i
                  ENDIF
                  DO j = i + ifu , N
                     tmp = tmp + ABS(Ap(jc))*ABS(X(j,k))
                     jc = jc + j
                  ENDDO
               ENDIF
            ELSEIF ( notran ) THEN
               jc = i
               DO j = 1 , i - ifu
                  tmp = tmp + ABS(Ap(jc))*ABS(X(j,k))
                  jc = jc + N - j
               ENDDO
               IF ( unit ) tmp = tmp + ABS(X(i,k))
            ELSE
               jc = (i-1)*(N-i) + (i*(i+1))/2
               IF ( unit ) tmp = tmp + ABS(X(i,k))
               DO j = i + ifu , N
                  tmp = tmp + ABS(Ap(jc+j-i))*ABS(X(j,k))
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
!     End of STPT05
!
      END SUBROUTINE STPT05
