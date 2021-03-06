!*==cgtts2.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b CGTTS2 solves a system of linear equations with a tridiagonal matrix using the LU factorization computed by sgttrf.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CGTTS2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgtts2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgtts2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgtts2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGTTS2( ITRANS, N, NRHS, DL, D, DU, DU2, IPIV, B, LDB )
!
!       .. Scalar Arguments ..
!       INTEGER            ITRANS, LDB, N, NRHS
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX            B( LDB, * ), D( * ), DL( * ), DU( * ), DU2( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CGTTS2 solves one of the systems of equations
!>    A * X = B,  A**T * X = B,  or  A**H * X = B,
!> with a tridiagonal matrix A using the LU factorization computed
!> by CGTTRF.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ITRANS
!> \verbatim
!>          ITRANS is INTEGER
!>          Specifies the form of the system of equations.
!>          = 0:  A * X = B     (No transpose)
!>          = 1:  A**T * X = B  (Transpose)
!>          = 2:  A**H * X = B  (Conjugate transpose)
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
!>          DL is COMPLEX array, dimension (N-1)
!>          The (n-1) multipliers that define the matrix L from the
!>          LU factorization of A.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is COMPLEX array, dimension (N)
!>          The n diagonal elements of the upper triangular matrix U from
!>          the LU factorization of A.
!> \endverbatim
!>
!> \param[in] DU
!> \verbatim
!>          DU is COMPLEX array, dimension (N-1)
!>          The (n-1) elements of the first super-diagonal of U.
!> \endverbatim
!>
!> \param[in] DU2
!> \verbatim
!>          DU2 is COMPLEX array, dimension (N-2)
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
!>          B is COMPLEX array, dimension (LDB,NRHS)
!>          On entry, the matrix of right hand side vectors B.
!>          On exit, B is overwritten by the solution vectors X.
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
!> \ingroup complexGTcomputational
!
!  =====================================================================
      SUBROUTINE CGTTS2(Itrans,N,Nrhs,Dl,D,Du,Du2,Ipiv,B,Ldb)
      IMPLICIT NONE
!*--CGTTS2132
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Itrans , Ldb , N , Nrhs
!     ..
!     .. Array Arguments ..
      INTEGER Ipiv(*)
      COMPLEX B(Ldb,*) , D(*) , Dl(*) , Du(*) , Du2(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER i , j
      COMPLEX temp
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC CONJG
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF ( N==0 .OR. Nrhs==0 ) RETURN
!
      IF ( Itrans==0 ) THEN
!
!        Solve A*X = B using the LU factorization of A,
!        overwriting each right hand side vector with its solution.
!
         IF ( Nrhs<=1 ) THEN
            j = 1
            DO
!
!           Solve L*x = b.
!
               DO i = 1 , N - 1
                  IF ( Ipiv(i)==i ) THEN
                     B(i+1,j) = B(i+1,j) - Dl(i)*B(i,j)
                  ELSE
                     temp = B(i,j)
                     B(i,j) = B(i+1,j)
                     B(i+1,j) = temp - Dl(i)*B(i,j)
                  ENDIF
               ENDDO
!
!           Solve U*x = b.
!
               B(N,j) = B(N,j)/D(N)
               IF ( N>1 ) B(N-1,j) = (B(N-1,j)-Du(N-1)*B(N,j))/D(N-1)
               DO i = N - 2 , 1 , -1
                  B(i,j) = (B(i,j)-Du(i)*B(i+1,j)-Du2(i)*B(i+2,j))/D(i)
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
!           Solve L*x = b.
!
               DO i = 1 , N - 1
                  IF ( Ipiv(i)==i ) THEN
                     B(i+1,j) = B(i+1,j) - Dl(i)*B(i,j)
                  ELSE
                     temp = B(i,j)
                     B(i,j) = B(i+1,j)
                     B(i+1,j) = temp - Dl(i)*B(i,j)
                  ENDIF
               ENDDO
!
!           Solve U*x = b.
!
               B(N,j) = B(N,j)/D(N)
               IF ( N>1 ) B(N-1,j) = (B(N-1,j)-Du(N-1)*B(N,j))/D(N-1)
               DO i = N - 2 , 1 , -1
                  B(i,j) = (B(i,j)-Du(i)*B(i+1,j)-Du2(i)*B(i+2,j))/D(i)
               ENDDO
            ENDDO
         ENDIF
      ELSEIF ( Itrans==1 ) THEN
!
!        Solve A**T * X = B.
!
         IF ( Nrhs<=1 ) THEN
            j = 1
            DO
!
!           Solve U**T * x = b.
!
               B(1,j) = B(1,j)/D(1)
               IF ( N>1 ) B(2,j) = (B(2,j)-Du(1)*B(1,j))/D(2)
               DO i = 3 , N
                  B(i,j) = (B(i,j)-Du(i-1)*B(i-1,j)-Du2(i-2)*B(i-2,j))  &
     &                     /D(i)
               ENDDO
!
!           Solve L**T * x = b.
!
               DO i = N - 1 , 1 , -1
                  IF ( Ipiv(i)==i ) THEN
                     B(i,j) = B(i,j) - Dl(i)*B(i+1,j)
                  ELSE
                     temp = B(i+1,j)
                     B(i+1,j) = B(i,j) - Dl(i)*temp
                     B(i,j) = temp
                  ENDIF
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
!           Solve U**T * x = b.
!
               B(1,j) = B(1,j)/D(1)
               IF ( N>1 ) B(2,j) = (B(2,j)-Du(1)*B(1,j))/D(2)
               DO i = 3 , N
                  B(i,j) = (B(i,j)-Du(i-1)*B(i-1,j)-Du2(i-2)*B(i-2,j))  &
     &                     /D(i)
               ENDDO
!
!           Solve L**T * x = b.
!
               DO i = N - 1 , 1 , -1
                  IF ( Ipiv(i)==i ) THEN
                     B(i,j) = B(i,j) - Dl(i)*B(i+1,j)
                  ELSE
                     temp = B(i+1,j)
                     B(i+1,j) = B(i,j) - Dl(i)*temp
                     B(i,j) = temp
                  ENDIF
               ENDDO
            ENDDO
         ENDIF
!
!        Solve A**H * X = B.
!
      ELSEIF ( Nrhs<=1 ) THEN
         j = 1
         DO
!
!           Solve U**H * x = b.
!
            B(1,j) = B(1,j)/CONJG(D(1))
            IF ( N>1 ) B(2,j) = (B(2,j)-CONJG(Du(1))*B(1,j))/CONJG(D(2))
            DO i = 3 , N
               B(i,j) = (B(i,j)-CONJG(Du(i-1))*B(i-1,j)-CONJG(Du2(i-2)) &
     &                  *B(i-2,j))/CONJG(D(i))
            ENDDO
!
!           Solve L**H * x = b.
!
            DO i = N - 1 , 1 , -1
               IF ( Ipiv(i)==i ) THEN
                  B(i,j) = B(i,j) - CONJG(Dl(i))*B(i+1,j)
               ELSE
                  temp = B(i+1,j)
                  B(i+1,j) = B(i,j) - CONJG(Dl(i))*temp
                  B(i,j) = temp
               ENDIF
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
!           Solve U**H * x = b.
!
            B(1,j) = B(1,j)/CONJG(D(1))
            IF ( N>1 ) B(2,j) = (B(2,j)-CONJG(Du(1))*B(1,j))/CONJG(D(2))
            DO i = 3 , N
               B(i,j) = (B(i,j)-CONJG(Du(i-1))*B(i-1,j)-CONJG(Du2(i-2)) &
     &                  *B(i-2,j))/CONJG(D(i))
            ENDDO
!
!           Solve L**H * x = b.
!
            DO i = N - 1 , 1 , -1
               IF ( Ipiv(i)==i ) THEN
                  B(i,j) = B(i,j) - CONJG(Dl(i))*B(i+1,j)
               ELSE
                  temp = B(i+1,j)
                  B(i+1,j) = B(i,j) - CONJG(Dl(i))*temp
                  B(i,j) = temp
               ENDIF
            ENDDO
         ENDDO
      ENDIF
!
!     End of CGTTS2
!
      END SUBROUTINE CGTTS2
