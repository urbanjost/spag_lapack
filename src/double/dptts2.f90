!*==dptts2.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b DPTTS2 solves a tridiagonal system of the form AX=B using the L D LH factorization computed by spttrf.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DPTTS2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dptts2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dptts2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dptts2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DPTTS2( N, NRHS, D, E, B, LDB )
!
!       .. Scalar Arguments ..
!       INTEGER            LDB, N, NRHS
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
!> DPTTS2 solves a tridiagonal system of the form
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
      SUBROUTINE DPTTS2(N,Nrhs,D,E,B,Ldb)
      IMPLICIT NONE
!*--DPTTS2106
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Ldb , N , Nrhs
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION B(Ldb,*) , D(*) , E(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER i , j
!     ..
!     .. External Subroutines ..
      EXTERNAL DSCAL
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF ( N<=1 ) THEN
         IF ( N==1 ) CALL DSCAL(Nrhs,1.D0/D(1),B,Ldb)
         RETURN
      ENDIF
!
!     Solve A * X = B using the factorization A = L*D*L**T,
!     overwriting each right hand side vector with its solution.
!
      DO j = 1 , Nrhs
!
!           Solve L * x = b.
!
         DO i = 2 , N
            B(i,j) = B(i,j) - B(i-1,j)*E(i-1)
         ENDDO
!
!           Solve D * L**T * x = b.
!
         B(N,j) = B(N,j)/D(N)
         DO i = N - 1 , 1 , -1
            B(i,j) = B(i,j)/D(i) - B(i+1,j)*E(i)
         ENDDO
      ENDDO
!
!
!     End of DPTTS2
!
      END SUBROUTINE DPTTS2
