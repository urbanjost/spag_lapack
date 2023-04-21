!*==cgtsv.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief <b> CGTSV computes the solution to system of linear equations A * X = B for GT matrices </b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CGTSV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgtsv.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgtsv.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgtsv.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGTSV( N, NRHS, DL, D, DU, B, LDB, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDB, N, NRHS
!       ..
!       .. Array Arguments ..
!       COMPLEX            B( LDB, * ), D( * ), DL( * ), DU( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CGTSV  solves the equation
!>
!>    A*X = B,
!>
!> where A is an N-by-N tridiagonal matrix, by Gaussian elimination with
!> partial pivoting.
!>
!> Note that the equation  A**T *X = B  may be solved by interchanging the
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
!>          DL is COMPLEX array, dimension (N-1)
!>          On entry, DL must contain the (n-1) subdiagonal elements of
!>          A.
!>          On exit, DL is overwritten by the (n-2) elements of the
!>          second superdiagonal of the upper triangular matrix U from
!>          the LU factorization of A, in DL(1), ..., DL(n-2).
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is COMPLEX array, dimension (N)
!>          On entry, D must contain the diagonal elements of A.
!>          On exit, D is overwritten by the n diagonal elements of U.
!> \endverbatim
!>
!> \param[in,out] DU
!> \verbatim
!>          DU is COMPLEX array, dimension (N-1)
!>          On entry, DU must contain the (n-1) superdiagonal elements
!>          of A.
!>          On exit, DU is overwritten by the (n-1) elements of the first
!>          superdiagonal of U.
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is COMPLEX array, dimension (LDB,NRHS)
!>          On entry, the N-by-NRHS right hand side matrix B.
!>          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
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
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  if INFO = i, U(i,i) is exactly zero, and the solution
!>                has not been computed.  The factorization has not been
!>                completed unless i = N.
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
!> \ingroup complexGTsolve
!
!  =====================================================================
      SUBROUTINE CGTSV(N,Nrhs,Dl,D,Du,B,Ldb,Info)
      IMPLICIT NONE
!*--CGTSV128
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
      COMPLEX B(Ldb,*) , D(*) , Dl(*) , Du(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX ZERO
      PARAMETER (ZERO=(0.0E+0,0.0E+0))
!     ..
!     .. Local Scalars ..
      INTEGER j , k
      COMPLEX mult , temp , zdum
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , AIMAG , MAX , REAL
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA
!     ..
!     .. Statement Functions ..
      REAL CABS1
!     ..
!     .. Statement Function definitions ..
      CABS1(zdum) = ABS(REAL(zdum)) + ABS(AIMAG(zdum))
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
         CALL XERBLA('CGTSV ',-Info)
         RETURN
      ENDIF
!
      IF ( N==0 ) RETURN
!
      DO k = 1 , N - 1
         IF ( Dl(k)==ZERO ) THEN
!
!           Subdiagonal is zero, no elimination is required.
!
            IF ( D(k)==ZERO ) THEN
!
!              Diagonal is zero: set INFO = K and return; a unique
!              solution can not be found.
!
               Info = k
               RETURN
            ENDIF
         ELSEIF ( CABS1(D(k))>=CABS1(Dl(k)) ) THEN
!
!           No row interchange required
!
            mult = Dl(k)/D(k)
            D(k+1) = D(k+1) - mult*Du(k)
            DO j = 1 , Nrhs
               B(k+1,j) = B(k+1,j) - mult*B(k,j)
            ENDDO
            IF ( k<(N-1) ) Dl(k) = ZERO
         ELSE
!
!           Interchange rows K and K+1
!
            mult = D(k)/Dl(k)
            D(k) = Dl(k)
            temp = D(k+1)
            D(k+1) = Du(k) - mult*temp
            IF ( k<(N-1) ) THEN
               Dl(k) = Du(k+1)
               Du(k+1) = -mult*Dl(k)
            ENDIF
            Du(k) = temp
            DO j = 1 , Nrhs
               temp = B(k,j)
               B(k,j) = B(k+1,j)
               B(k+1,j) = temp - mult*B(k+1,j)
            ENDDO
         ENDIF
      ENDDO
      IF ( D(N)==ZERO ) THEN
         Info = N
         RETURN
      ENDIF
!
!     Back solve with the matrix U from the factorization.
!
      DO j = 1 , Nrhs
         B(N,j) = B(N,j)/D(N)
         IF ( N>1 ) B(N-1,j) = (B(N-1,j)-Du(N-1)*B(N,j))/D(N-1)
         DO k = N - 2 , 1 , -1
            B(k,j) = (B(k,j)-Du(k)*B(k+1,j)-Dl(k)*B(k+2,j))/D(k)
         ENDDO
      ENDDO
!
!
!     End of CGTSV
!
      END SUBROUTINE CGTSV
