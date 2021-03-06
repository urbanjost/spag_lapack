!*==sgtsv.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief <b> SGTSV computes the solution to system of linear equations A * X = B for GT matrices </b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SGTSV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgtsv.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgtsv.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgtsv.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGTSV( N, NRHS, DL, D, DU, B, LDB, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDB, N, NRHS
!       ..
!       .. Array Arguments ..
!       REAL               B( LDB, * ), D( * ), DL( * ), DU( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SGTSV  solves the equation
!>
!>    A*X = B,
!>
!> where A is an n by n tridiagonal matrix, by Gaussian elimination with
!> partial pivoting.
!>
!> Note that the equation  A**T*X = B  may be solved by interchanging the
!> order of the arguments DU and DL.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of right hand sides, i.e., the number of columns
!>          of the matrix B.  NRHS >= 0.
!> \endverbatim
!>
!> \param[in,out] DL
!> \verbatim
!>          DL is REAL array, dimension (N-1)
!>          On entry, DL must contain the (n-1) sub-diagonal elements of
!>          A.
!>
!>          On exit, DL is overwritten by the (n-2) elements of the
!>          second super-diagonal of the upper triangular matrix U from
!>          the LU factorization of A, in DL(1), ..., DL(n-2).
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is REAL array, dimension (N)
!>          On entry, D must contain the diagonal elements of A.
!>
!>          On exit, D is overwritten by the n diagonal elements of U.
!> \endverbatim
!>
!> \param[in,out] DU
!> \verbatim
!>          DU is REAL array, dimension (N-1)
!>          On entry, DU must contain the (n-1) super-diagonal elements
!>          of A.
!>
!>          On exit, DU is overwritten by the (n-1) elements of the first
!>          super-diagonal of U.
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is REAL array, dimension (LDB,NRHS)
!>          On entry, the N by NRHS matrix of right hand side matrix B.
!>          On exit, if INFO = 0, the N by NRHS solution matrix X.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,N).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
!>          > 0: if INFO = i, U(i,i) is exactly zero, and the solution
!>               has not been computed.  The factorization has not been
!>               completed unless i = N.
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
!> \ingroup realGTsolve
!
!  =====================================================================
      SUBROUTINE SGTSV(N,Nrhs,Dl,D,Du,B,Ldb,Info)
      IMPLICIT NONE
!*--SGTSV131
!
!  -- LAPACK driver routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Info , Ldb , N , Nrhs
!     ..
!     .. Array Arguments ..
      REAL B(Ldb,*) , D(*) , Dl(*) , Du(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO
      PARAMETER (ZERO=0.0E+0)
!     ..
!     .. Local Scalars ..
      INTEGER i , j
      REAL fact , temp
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA
!     ..
!     .. Executable Statements ..
!
      Info = 0
      IF ( N<0 ) THEN
         Info = -1
      ELSEIF ( Nrhs<0 ) THEN
         Info = -2
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -7
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('SGTSV ',-Info)
         RETURN
      ENDIF
!
      IF ( N==0 ) RETURN
!
      IF ( Nrhs==1 ) THEN
         DO i = 1 , N - 2
            IF ( ABS(D(i))>=ABS(Dl(i)) ) THEN
!
!              No row interchange required
!
               IF ( D(i)/=ZERO ) THEN
                  fact = Dl(i)/D(i)
                  D(i+1) = D(i+1) - fact*Du(i)
                  B(i+1,1) = B(i+1,1) - fact*B(i,1)
               ELSE
                  Info = i
                  RETURN
               ENDIF
               Dl(i) = ZERO
            ELSE
!
!              Interchange rows I and I+1
!
               fact = D(i)/Dl(i)
               D(i) = Dl(i)
               temp = D(i+1)
               D(i+1) = Du(i) - fact*temp
               Dl(i) = Du(i+1)
               Du(i+1) = -fact*Dl(i)
               Du(i) = temp
               temp = B(i,1)
               B(i,1) = B(i+1,1)
               B(i+1,1) = temp - fact*B(i+1,1)
            ENDIF
         ENDDO
         IF ( N>1 ) THEN
            i = N - 1
            IF ( ABS(D(i))<ABS(Dl(i)) ) THEN
               fact = D(i)/Dl(i)
               D(i) = Dl(i)
               temp = D(i+1)
               D(i+1) = Du(i) - fact*temp
               Du(i) = temp
               temp = B(i,1)
               B(i,1) = B(i+1,1)
               B(i+1,1) = temp - fact*B(i+1,1)
            ELSEIF ( D(i)/=ZERO ) THEN
               fact = Dl(i)/D(i)
               D(i+1) = D(i+1) - fact*Du(i)
               B(i+1,1) = B(i+1,1) - fact*B(i,1)
            ELSE
               Info = i
               RETURN
            ENDIF
         ENDIF
         IF ( D(N)==ZERO ) THEN
            Info = N
            RETURN
         ENDIF
      ELSE
         DO i = 1 , N - 2
            IF ( ABS(D(i))>=ABS(Dl(i)) ) THEN
!
!              No row interchange required
!
               IF ( D(i)/=ZERO ) THEN
                  fact = Dl(i)/D(i)
                  D(i+1) = D(i+1) - fact*Du(i)
                  DO j = 1 , Nrhs
                     B(i+1,j) = B(i+1,j) - fact*B(i,j)
                  ENDDO
               ELSE
                  Info = i
                  RETURN
               ENDIF
               Dl(i) = ZERO
            ELSE
!
!              Interchange rows I and I+1
!
               fact = D(i)/Dl(i)
               D(i) = Dl(i)
               temp = D(i+1)
               D(i+1) = Du(i) - fact*temp
               Dl(i) = Du(i+1)
               Du(i+1) = -fact*Dl(i)
               Du(i) = temp
               DO j = 1 , Nrhs
                  temp = B(i,j)
                  B(i,j) = B(i+1,j)
                  B(i+1,j) = temp - fact*B(i+1,j)
               ENDDO
            ENDIF
         ENDDO
         IF ( N>1 ) THEN
            i = N - 1
            IF ( ABS(D(i))<ABS(Dl(i)) ) THEN
               fact = D(i)/Dl(i)
               D(i) = Dl(i)
               temp = D(i+1)
               D(i+1) = Du(i) - fact*temp
               Du(i) = temp
               DO j = 1 , Nrhs
                  temp = B(i,j)
                  B(i,j) = B(i+1,j)
                  B(i+1,j) = temp - fact*B(i+1,j)
               ENDDO
            ELSEIF ( D(i)/=ZERO ) THEN
               fact = Dl(i)/D(i)
               D(i+1) = D(i+1) - fact*Du(i)
               DO j = 1 , Nrhs
                  B(i+1,j) = B(i+1,j) - fact*B(i,j)
               ENDDO
            ELSE
               Info = i
               RETURN
            ENDIF
         ENDIF
         IF ( D(N)==ZERO ) THEN
            Info = N
            RETURN
         ENDIF
      ENDIF
!
!     Back solve with the matrix U from the factorization.
!
      IF ( Nrhs<=2 ) THEN
         j = 1
         DO
            B(N,j) = B(N,j)/D(N)
            IF ( N>1 ) B(N-1,j) = (B(N-1,j)-Du(N-1)*B(N,j))/D(N-1)
            DO i = N - 2 , 1 , -1
               B(i,j) = (B(i,j)-Du(i)*B(i+1,j)-Dl(i)*B(i+2,j))/D(i)
            ENDDO
            IF ( j<Nrhs ) THEN
               j = j + 1
               CYCLE
            ENDIF
            EXIT
         ENDDO
      ELSE
         DO j = 1 , Nrhs
            B(N,j) = B(N,j)/D(N)
            IF ( N>1 ) B(N-1,j) = (B(N-1,j)-Du(N-1)*B(N,j))/D(N-1)
            DO i = N - 2 , 1 , -1
               B(i,j) = (B(i,j)-Du(i)*B(i+1,j)-Dl(i)*B(i+2,j))/D(i)
            ENDDO
         ENDDO
      ENDIF
!
!
!     End of SGTSV
!
      END SUBROUTINE SGTSV
