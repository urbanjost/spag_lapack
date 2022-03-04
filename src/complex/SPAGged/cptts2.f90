!*==cptts2.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CPTTS2 solves a tridiagonal system of the form AX=B using the L D LH factorization computed by spttrf.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CPTTS2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cptts2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cptts2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cptts2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CPTTS2( IUPLO, N, NRHS, D, E, B, LDB )
!
!       .. Scalar Arguments ..
!       INTEGER            IUPLO, LDB, N, NRHS
!       ..
!       .. Array Arguments ..
!       REAL               D( * )
!       COMPLEX            B( LDB, * ), E( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CPTTS2 solves a tridiagonal system of the form
!>    A * X = B
!> using the factorization A = U**H*D*U or A = L*D*L**H computed by CPTTRF.
!> D is a diagonal matrix specified in the vector D, U (or L) is a unit
!> bidiagonal matrix whose superdiagonal (subdiagonal) is specified in
!> the vector E, and X and B are N by NRHS matrices.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] IUPLO
!> \verbatim
!>          IUPLO is INTEGER
!>          Specifies the form of the factorization and whether the
!>          vector E is the superdiagonal of the upper bidiagonal factor
!>          U or the subdiagonal of the lower bidiagonal factor L.
!>          = 1:  A = U**H *D*U, E is the superdiagonal of U
!>          = 0:  A = L*D*L**H, E is the subdiagonal of L
!> \endverbatim
!>
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
!>          D is REAL array, dimension (N)
!>          The n diagonal elements of the diagonal matrix D from the
!>          factorization A = U**H *D*U or A = L*D*L**H.
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is COMPLEX array, dimension (N-1)
!>          If IUPLO = 1, the (n-1) superdiagonal elements of the unit
!>          bidiagonal factor U from the factorization A = U**H*D*U.
!>          If IUPLO = 0, the (n-1) subdiagonal elements of the unit
!>          bidiagonal factor L from the factorization A = L*D*L**H.
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is COMPLEX array, dimension (LDB,NRHS)
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
!> \date June 2016
!
!> \ingroup complexPTcomputational
!
!  =====================================================================
      SUBROUTINE CPTTS2(Iuplo,N,Nrhs,D,E,B,Ldb)
      USE S_CSSCAL
      IMPLICIT NONE
!*--CPTTS2118
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: Iuplo
      INTEGER , INTENT(IN) :: N
      INTEGER :: Nrhs
      REAL , DIMENSION(*) :: D
      COMPLEX , INTENT(IN) , DIMENSION(*) :: E
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , j
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
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF ( N<=1 ) THEN
         IF ( N==1 ) CALL CSSCAL(Nrhs,1./D(1),B,Ldb)
         RETURN
      ENDIF
!
      IF ( Iuplo==1 ) THEN
!
!        Solve A * X = B using the factorization A = U**H *D*U,
!        overwriting each right hand side vector with its solution.
!
         IF ( Nrhs<=2 ) THEN
            j = 1
            DO
!
!           Solve U**H * x = b.
!
               DO i = 2 , N
                  B(i,j) = B(i,j) - B(i-1,j)*CONJG(E(i-1))
               ENDDO
!
!           Solve D * U * x = b.
!
               DO i = 1 , N
                  B(i,j) = B(i,j)/D(i)
               ENDDO
               DO i = N - 1 , 1 , -1
                  B(i,j) = B(i,j) - B(i+1,j)*E(i)
               ENDDO
               IF ( j<Nrhs ) THEN
                  j = j + 1
                  CYCLE
               ENDIF
               EXIT
            ENDDO
         ELSE
            DO j = 1 , Nrhs
!
!              Solve U**H * x = b.
!
               DO i = 2 , N
                  B(i,j) = B(i,j) - B(i-1,j)*CONJG(E(i-1))
               ENDDO
!
!              Solve D * U * x = b.
!
               B(N,j) = B(N,j)/D(N)
               DO i = N - 1 , 1 , -1
                  B(i,j) = B(i,j)/D(i) - B(i+1,j)*E(i)
               ENDDO
            ENDDO
         ENDIF
!
!        Solve A * X = B using the factorization A = L*D*L**H,
!        overwriting each right hand side vector with its solution.
!
      ELSEIF ( Nrhs<=2 ) THEN
         j = 1
         DO
!
!           Solve L * x = b.
!
            DO i = 2 , N
               B(i,j) = B(i,j) - B(i-1,j)*E(i-1)
            ENDDO
!
!           Solve D * L**H * x = b.
!
            DO i = 1 , N
               B(i,j) = B(i,j)/D(i)
            ENDDO
            DO i = N - 1 , 1 , -1
               B(i,j) = B(i,j) - B(i+1,j)*CONJG(E(i))
            ENDDO
            IF ( j<Nrhs ) THEN
               j = j + 1
               CYCLE
            ENDIF
            EXIT
         ENDDO
      ELSE
         DO j = 1 , Nrhs
!
!              Solve L * x = b.
!
            DO i = 2 , N
               B(i,j) = B(i,j) - B(i-1,j)*E(i-1)
            ENDDO
!
!              Solve D * L**H * x = b.
!
            B(N,j) = B(N,j)/D(N)
            DO i = N - 1 , 1 , -1
               B(i,j) = B(i,j)/D(i) - B(i+1,j)*CONJG(E(i))
            ENDDO
         ENDDO
      ENDIF
!
!
!     End of CPTTS2
!
      END SUBROUTINE CPTTS2
