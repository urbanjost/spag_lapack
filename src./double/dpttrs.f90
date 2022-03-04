!*==dpttrs.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DPTTRS
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DPTTRS + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dpttrs.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dpttrs.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dpttrs.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DPTTRS( N, NRHS, D, E, B, LDB, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDB, N, NRHS
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   B( LDB, * ), D( * ), E( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DPTTRS solves a tridiagonal system of the form
!>    A * X = B
!> using the L*D*L**T factorization of A computed by DPTTRF.  D is a
!> diagonal matrix specified in the vector D, L is a unit bidiagonal
!> matrix whose subdiagonal is specified in the vector E, and X and B
!> are N by NRHS matrices.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the tridiagonal matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of right hand sides, i.e., the number of columns
!>          of the matrix B.  NRHS >= 0.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>          The n diagonal elements of the diagonal matrix D from the
!>          L*D*L**T factorization of A.
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is DOUBLE PRECISION array, dimension (N-1)
!>          The (n-1) subdiagonal elements of the unit bidiagonal factor
!>          L from the L*D*L**T factorization of A.  E can also be regarded
!>          as the superdiagonal of the unit bidiagonal factor U from the
!>          factorization A = U**T*D*U.
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)
!>          On entry, the right hand side vectors B for the system of
!>          linear equations.
!>          On exit, the solution vectors, X.
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
!>          < 0: if INFO = -k, the k-th argument had an illegal value
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
!> \ingroup doublePTcomputational
!
!  =====================================================================
      SUBROUTINE DPTTRS(N,Nrhs,D,E,B,Ldb,Info)
      USE F77KINDS                        
      USE S_DPTTS2
      USE S_ILAENV
      USE S_XERBLA
      IMPLICIT NONE
!*--DPTTRS117
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER :: N
      INTEGER :: Nrhs
      REAL(R8KIND) , DIMENSION(*) :: D
      REAL(R8KIND) , DIMENSION(*) :: E
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: j , jb , nb
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
!     Test the input arguments.
!
      Info = 0
      IF ( N<0 ) THEN
         Info = -1
      ELSEIF ( Nrhs<0 ) THEN
         Info = -2
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -6
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('DPTTRS',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 .OR. Nrhs==0 ) RETURN
!
!     Determine the number of right-hand sides to solve at a time.
!
      IF ( Nrhs==1 ) THEN
         nb = 1
      ELSE
         nb = MAX(1,ILAENV(1,'DPTTRS',' ',N,Nrhs,-1,-1))
      ENDIF
!
      IF ( nb>=Nrhs ) THEN
         CALL DPTTS2(N,Nrhs,D,E,B,Ldb)
      ELSE
         DO j = 1 , Nrhs , nb
            jb = MIN(Nrhs-j+1,nb)
            CALL DPTTS2(N,jb,D,E,B(1,j),Ldb)
         ENDDO
      ENDIF
!
!
!     End of DPTTRS
!
      END SUBROUTINE DPTTRS
