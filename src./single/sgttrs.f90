!*==sgttrs.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SGTTRS
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SGTTRS + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgttrs.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgttrs.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgttrs.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGTTRS( TRANS, N, NRHS, DL, D, DU, DU2, IPIV, B, LDB,
!                          INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANS
!       INTEGER            INFO, LDB, N, NRHS
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       REAL               B( LDB, * ), D( * ), DL( * ), DU( * ), DU2( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SGTTRS solves one of the systems of equations
!>    A*X = B  or  A**T*X = B,
!> with a tridiagonal matrix A using the LU factorization computed
!> by SGTTRF.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          Specifies the form of the system of equations.
!>          = 'N':  A * X = B  (No transpose)
!>          = 'T':  A**T* X = B  (Transpose)
!>          = 'C':  A**T* X = B  (Conjugate transpose = Transpose)
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of right hand sides, i.e., the number of columns
!>          of the matrix B.  NRHS >= 0.
!> \endverbatim
!>
!> \param[in] DL
!> \verbatim
!>          DL is REAL array, dimension (N-1)
!>          The (n-1) multipliers that define the matrix L from the
!>          LU factorization of A.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is REAL array, dimension (N)
!>          The n diagonal elements of the upper triangular matrix U from
!>          the LU factorization of A.
!> \endverbatim
!>
!> \param[in] DU
!> \verbatim
!>          DU is REAL array, dimension (N-1)
!>          The (n-1) elements of the first super-diagonal of U.
!> \endverbatim
!>
!> \param[in] DU2
!> \verbatim
!>          DU2 is REAL array, dimension (N-2)
!>          The (n-2) elements of the second super-diagonal of U.
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          The pivot indices; for 1 <= i <= n, row i of the matrix was
!>          interchanged with row IPIV(i).  IPIV(i) will always be either
!>          i or i+1; IPIV(i) = i indicates a row interchange was not
!>          required.
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is REAL array, dimension (LDB,NRHS)
!>          On entry, the matrix of right hand side vectors B.
!>          On exit, B is overwritten by the solution vectors X.
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
!> \ingroup realGTcomputational
!
!  =====================================================================
      SUBROUTINE SGTTRS(Trans,N,Nrhs,Dl,D,Du,Du2,Ipiv,B,Ldb,Info)
      USE S_ILAENV
      USE S_SGTTS2
      USE S_XERBLA
      IMPLICIT NONE
!*--SGTTRS144
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Trans
      INTEGER :: N
      INTEGER :: Nrhs
      REAL , DIMENSION(*) :: Dl
      REAL , DIMENSION(*) :: D
      REAL , DIMENSION(*) :: Du
      REAL , DIMENSION(*) :: Du2
      INTEGER , DIMENSION(*) :: Ipiv
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: itrans , j , jb , nb
      LOGICAL :: notran
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
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
      notran = (Trans=='N' .OR. Trans=='n')
      IF ( .NOT.notran .AND. .NOT.(Trans=='T' .OR. Trans=='t') .AND.    &
     &     .NOT.(Trans=='C' .OR. Trans=='c') ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Nrhs<0 ) THEN
         Info = -3
      ELSEIF ( Ldb<MAX(N,1) ) THEN
         Info = -10
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('SGTTRS',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 .OR. Nrhs==0 ) RETURN
!
!     Decode TRANS
!
      IF ( notran ) THEN
         itrans = 0
      ELSE
         itrans = 1
      ENDIF
!
!     Determine the number of right-hand sides to solve at a time.
!
      IF ( Nrhs==1 ) THEN
         nb = 1
      ELSE
         nb = MAX(1,ILAENV(1,'SGTTRS',Trans,N,Nrhs,-1,-1))
      ENDIF
!
      IF ( nb>=Nrhs ) THEN
         CALL SGTTS2(itrans,N,Nrhs,Dl,D,Du,Du2,Ipiv,B,Ldb)
      ELSE
         DO j = 1 , Nrhs , nb
            jb = MIN(Nrhs-j+1,nb)
            CALL SGTTS2(itrans,N,jb,Dl,D,Du,Du2,Ipiv,B(1,j),Ldb)
         ENDDO
      ENDIF
!
!     End of SGTTRS
!
      END SUBROUTINE SGTTRS
