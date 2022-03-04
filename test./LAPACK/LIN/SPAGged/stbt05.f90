!*==stbt05.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b STBT05
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE STBT05( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, B,
!                          LDB, X, LDX, XACT, LDXACT, FERR, BERR, RESLTS )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIAG, TRANS, UPLO
!       INTEGER            KD, LDAB, LDB, LDX, LDXACT, N, NRHS
!       ..
!       .. Array Arguments ..
!       REAL               AB( LDAB, * ), B( LDB, * ), BERR( * ),
!      $                   FERR( * ), RESLTS( * ), X( LDX, * ),
!      $                   XACT( LDXACT, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> STBT05 tests the error bounds from iterative refinement for the
!> computed solution to a system of equations A*X = B, where A is a
!> triangular band matrix.
!>
!> RESLTS(1) = test of the error bound
!>           = norm(X - XACT) / ( norm(X) * FERR )
!>
!> A large value is returned if this ratio is not less than one.
!>
!> RESLTS(2) = residual from the iterative refinement routine
!>           = the maximum of BERR / ( NZ*EPS + (*) ), where
!>             (*) = NZ*UNFL / (min_i (abs(A)*abs(X) +abs(b))_i )
!>             and NZ = max. number of nonzeros in any row of A, plus 1
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
!> \param[in] KD
!> \verbatim
!>          KD is INTEGER
!>          The number of super-diagonals of the matrix A if UPLO = 'U',
!>          or the number of sub-diagonals if UPLO = 'L'.  KD >= 0.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of columns of the matrices X, B, and XACT.
!>          NRHS >= 0.
!> \endverbatim
!>
!> \param[in] AB
!> \verbatim
!>          AB is REAL array, dimension (LDAB,N)
!>          The upper or lower triangular band matrix A, stored in the
!>          first kd+1 rows of the array. The j-th column of A is stored
!>          in the j-th column of the array AB as follows:
!>          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
!>          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
!>          If DIAG = 'U', the diagonal elements of A are not referenced
!>          and are assumed to be 1.
!> \endverbatim
!>
!> \param[in] LDAB
!> \verbatim
!>          LDAB is INTEGER
!>          The leading dimension of the array AB.  LDAB >= KD+1.
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
!>          RESLTS(2) = BERR / ( NZ*EPS + (*) )
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
      SUBROUTINE STBT05(Uplo,Trans,Diag,N,Kd,Nrhs,Ab,Ldab,B,Ldb,X,Ldx,  &
     &                  Xact,Ldxact,Ferr,Berr,Reslts)
      IMPLICIT NONE
!*--STBT05193
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Diag , Trans , Uplo
      INTEGER Kd , Ldab , Ldb , Ldx , Ldxact , N , Nrhs
!     ..
!     .. Array Arguments ..
      REAL Ab(Ldab,*) , B(Ldb,*) , Berr(*) , Ferr(*) , Reslts(*) ,      &
     &     X(Ldx,*) , Xact(Ldxact,*)
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
      INTEGER i , ifu , imax , j , k , nz
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
      nz = MIN(Kd,N-1) + 1
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
!     Test 2:  Compute the maximum of BERR / ( NZ*EPS + (*) ), where
!     (*) = NZ*UNFL / (min_i (abs(A)*abs(X) +abs(b))_i )
!
      ifu = 0
      IF ( unit ) ifu = 1
      DO k = 1 , Nrhs
         DO i = 1 , N
            tmp = ABS(B(i,k))
            IF ( upper ) THEN
               IF ( .NOT.notran ) THEN
                  DO j = MAX(i-Kd,1) , i - ifu
                     tmp = tmp + ABS(Ab(Kd+1-i+j,i))*ABS(X(j,k))
                  ENDDO
                  IF ( unit ) tmp = tmp + ABS(X(i,k))
               ELSE
                  IF ( unit ) tmp = tmp + ABS(X(i,k))
                  DO j = i + ifu , MIN(i+Kd,N)
                     tmp = tmp + ABS(Ab(Kd+1+i-j,j))*ABS(X(j,k))
                  ENDDO
               ENDIF
            ELSEIF ( notran ) THEN
               DO j = MAX(i-Kd,1) , i - ifu
                  tmp = tmp + ABS(Ab(1+i-j,j))*ABS(X(j,k))
               ENDDO
               IF ( unit ) tmp = tmp + ABS(X(i,k))
            ELSE
               IF ( unit ) tmp = tmp + ABS(X(i,k))
               DO j = i + ifu , MIN(i+Kd,N)
                  tmp = tmp + ABS(Ab(1+j-i,i))*ABS(X(j,k))
               ENDDO
            ENDIF
            IF ( i==1 ) THEN
               axbi = tmp
            ELSE
               axbi = MIN(axbi,tmp)
            ENDIF
         ENDDO
         tmp = Berr(k)/(nz*eps+nz*unfl/MAX(axbi,nz*unfl))
         IF ( k==1 ) THEN
            Reslts(2) = tmp
         ELSE
            Reslts(2) = MAX(Reslts(2),tmp)
         ENDIF
      ENDDO
!
!
!     End of STBT05
!
      END SUBROUTINE STBT05
