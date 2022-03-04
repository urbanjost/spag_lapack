!*==sgtt05.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b SGTT05
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGTT05( TRANS, N, NRHS, DL, D, DU, B, LDB, X, LDX,
!                          XACT, LDXACT, FERR, BERR, RESLTS )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANS
!       INTEGER            LDB, LDX, LDXACT, N, NRHS
!       ..
!       .. Array Arguments ..
!       REAL               B( LDB, * ), BERR( * ), D( * ), DL( * ),
!      $                   DU( * ), FERR( * ), RESLTS( * ), X( LDX, * ),
!      $                   XACT( LDXACT, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SGTT05 tests the error bounds from iterative refinement for the
!> computed solution to a system of equations A*X = B, where A is a
!> general tridiagonal matrix of order n and op(A) = A or A**T,
!> depending on TRANS.
!>
!> RESLTS(1) = test of the error bound
!>           = norm(X - XACT) / ( norm(X) * FERR )
!>
!> A large value is returned if this ratio is not less than one.
!>
!> RESLTS(2) = residual from the iterative refinement routine
!>           = the maximum of BERR / ( NZ*EPS + (*) ), where
!>             (*) = NZ*UNFL / (min_i (abs(op(A))*abs(X) +abs(b))_i )
!>             and NZ = max. number of nonzeros in any row of A, plus 1
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
!> \param[in] DL
!> \verbatim
!>          DL is REAL array, dimension (N-1)
!>          The (n-1) sub-diagonal elements of A.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is REAL array, dimension (N)
!>          The diagonal elements of A.
!> \endverbatim
!>
!> \param[in] DU
!> \verbatim
!>          DU is REAL array, dimension (N-1)
!>          The (n-1) super-diagonal elements of A.
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
      SUBROUTINE SGTT05(Trans,N,Nrhs,Dl,D,Du,B,Ldb,X,Ldx,Xact,Ldxact,   &
     &                  Ferr,Berr,Reslts)
      IMPLICIT NONE
!*--SGTT05169
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Trans
      INTEGER Ldb , Ldx , Ldxact , N , Nrhs
!     ..
!     .. Array Arguments ..
      REAL B(Ldb,*) , Berr(*) , D(*) , Dl(*) , Du(*) , Ferr(*) ,        &
     &     Reslts(*) , X(Ldx,*) , Xact(Ldxact,*)
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
      INTEGER i , imax , j , k , nz
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
      notran = LSAME(Trans,'N')
      nz = 4
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
!     (*) = NZ*UNFL / (min_i (abs(op(A))*abs(X) +abs(b))_i )
!
      DO k = 1 , Nrhs
         IF ( notran ) THEN
            IF ( N==1 ) THEN
               axbi = ABS(B(1,k)) + ABS(D(1)*X(1,k))
            ELSE
               axbi = ABS(B(1,k)) + ABS(D(1)*X(1,k)) + ABS(Du(1)*X(2,k))
               DO i = 2 , N - 1
                  tmp = ABS(B(i,k)) + ABS(Dl(i-1)*X(i-1,k))             &
     &                  + ABS(D(i)*X(i,k)) + ABS(Du(i)*X(i+1,k))
                  axbi = MIN(axbi,tmp)
               ENDDO
               tmp = ABS(B(N,k)) + ABS(Dl(N-1)*X(N-1,k))                &
     &               + ABS(D(N)*X(N,k))
               axbi = MIN(axbi,tmp)
            ENDIF
         ELSEIF ( N==1 ) THEN
            axbi = ABS(B(1,k)) + ABS(D(1)*X(1,k))
         ELSE
            axbi = ABS(B(1,k)) + ABS(D(1)*X(1,k)) + ABS(Dl(1)*X(2,k))
            DO i = 2 , N - 1
               tmp = ABS(B(i,k)) + ABS(Du(i-1)*X(i-1,k))                &
     &               + ABS(D(i)*X(i,k)) + ABS(Dl(i)*X(i+1,k))
               axbi = MIN(axbi,tmp)
            ENDDO
            tmp = ABS(B(N,k)) + ABS(Du(N-1)*X(N-1,k)) + ABS(D(N)*X(N,k))
            axbi = MIN(axbi,tmp)
         ENDIF
         tmp = Berr(k)/(nz*eps+nz*unfl/MAX(axbi,nz*unfl))
         IF ( k==1 ) THEN
            Reslts(2) = tmp
         ELSE
            Reslts(2) = MAX(Reslts(2),tmp)
         ENDIF
      ENDDO
!
!
!     End of SGTT05
!
      END SUBROUTINE SGTT05
