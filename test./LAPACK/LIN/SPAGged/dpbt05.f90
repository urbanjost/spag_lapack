!*==dpbt05.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b DPBT05
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DPBT05( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, X, LDX,
!                          XACT, LDXACT, FERR, BERR, RESLTS )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            KD, LDAB, LDB, LDX, LDXACT, N, NRHS
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   AB( LDAB, * ), B( LDB, * ), BERR( * ),
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
!> DPBT05 tests the error bounds from iterative refinement for the
!> computed solution to a system of equations A*X = B, where A is a
!> symmetric band matrix.
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
!>          Specifies whether the upper or lower triangular part of the
!>          symmetric matrix A is stored.
!>          = 'U':  Upper triangular
!>          = 'L':  Lower triangular
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
!>          AB is DOUBLE PRECISION array, dimension (LDAB,N)
!>          The upper or lower triangle of the symmetric band matrix A,
!>          stored in the first KD+1 rows of the array.  The j-th column
!>          of A is stored in the j-th column of the array AB as follows:
!>          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
!>          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
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
!>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)
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
!>          X is DOUBLE PRECISION array, dimension (LDX,NRHS)
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
!>          XACT is DOUBLE PRECISION array, dimension (LDX,NRHS)
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
!>          FERR is DOUBLE PRECISION array, dimension (NRHS)
!>          The estimated forward error bounds for each solution vector
!>          X.  If XTRUE is the true solution, FERR bounds the magnitude
!>          of the largest entry in (X - XTRUE) divided by the magnitude
!>          of the largest entry in X.
!> \endverbatim
!>
!> \param[in] BERR
!> \verbatim
!>          BERR is DOUBLE PRECISION array, dimension (NRHS)
!>          The componentwise relative backward error of each solution
!>          vector (i.e., the smallest relative change in any entry of A
!>          or B that makes X an exact solution).
!> \endverbatim
!>
!> \param[out] RESLTS
!> \verbatim
!>          RESLTS is DOUBLE PRECISION array, dimension (2)
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
!> \ingroup double_lin
!
!  =====================================================================
      SUBROUTINE DPBT05(Uplo,N,Kd,Nrhs,Ab,Ldab,B,Ldb,X,Ldx,Xact,Ldxact, &
     &                  Ferr,Berr,Reslts)
      IMPLICIT NONE
!*--DPBT05175
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER Kd , Ldab , Ldb , Ldx , Ldxact , N , Nrhs
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION Ab(Ldab,*) , B(Ldb,*) , Berr(*) , Ferr(*) ,      &
     &                 Reslts(*) , X(Ldx,*) , Xact(Ldxact,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
!     ..
!     .. Local Scalars ..
      LOGICAL upper
      INTEGER i , imax , j , k , nz
      DOUBLE PRECISION axbi , diff , eps , errbnd , ovfl , tmp , unfl , &
     &                 xnorm
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER IDAMAX
      DOUBLE PRECISION DLAMCH
      EXTERNAL LSAME , IDAMAX , DLAMCH
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
      eps = DLAMCH('Epsilon')
      unfl = DLAMCH('Safe minimum')
      ovfl = ONE/unfl
      upper = LSAME(Uplo,'U')
      nz = 2*MAX(Kd,N-1) + 1
!
!     Test 1:  Compute the maximum of
!        norm(X - XACT) / ( norm(X) * FERR )
!     over all the vectors X and XACT using the infinity-norm.
!
      errbnd = ZERO
      DO j = 1 , Nrhs
         imax = IDAMAX(N,X(1,j),1)
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
      DO k = 1 , Nrhs
         DO i = 1 , N
            tmp = ABS(B(i,k))
            IF ( upper ) THEN
               DO j = MAX(i-Kd,1) , i
                  tmp = tmp + ABS(Ab(Kd+1-i+j,i))*ABS(X(j,k))
               ENDDO
               DO j = i + 1 , MIN(i+Kd,N)
                  tmp = tmp + ABS(Ab(Kd+1+i-j,j))*ABS(X(j,k))
               ENDDO
            ELSE
               DO j = MAX(i-Kd,1) , i - 1
                  tmp = tmp + ABS(Ab(1+i-j,j))*ABS(X(j,k))
               ENDDO
               DO j = i , MIN(i+Kd,N)
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
!     End of DPBT05
!
      END SUBROUTINE DPBT05
