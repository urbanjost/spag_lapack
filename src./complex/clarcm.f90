!*==clarcm.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CLARCM copies all or part of a real two-dimensional array to a complex array.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLARCM + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clarcm.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clarcm.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clarcm.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLARCM( M, N, A, LDA, B, LDB, C, LDC, RWORK )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LDB, LDC, M, N
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), RWORK( * )
!       COMPLEX            B( LDB, * ), C( LDC, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLARCM performs a very simple matrix-matrix multiplication:
!>          C := A * B,
!> where A is M by M and real; B is M by N and complex;
!> C is M by N and complex.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A and of the matrix C.
!>          M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns and rows of the matrix B and
!>          the number of columns of the matrix C.
!>          N >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is REAL array, dimension (LDA, M)
!>          On entry, A contains the M by M matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A. LDA >=max(1,M).
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is COMPLEX array, dimension (LDB, N)
!>          On entry, B contains the M by N matrix B.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B. LDB >=max(1,M).
!> \endverbatim
!>
!> \param[out] C
!> \verbatim
!>          C is COMPLEX array, dimension (LDC, N)
!>          On exit, C contains the M by N matrix C.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of the array C. LDC >=max(1,M).
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (2*M*N)
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
!> \ingroup complexOTHERauxiliary
!
!  =====================================================================
      SUBROUTINE CLARCM(M,N,A,Lda,B,Ldb,C,Ldc,Rwork)
      USE S_SGEMM
      IMPLICIT NONE
!*--CLARCM119
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ONE = 1.0E0 , ZERO = 0.0E0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER :: M
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , INTENT(IN) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER , INTENT(IN) :: Ldc
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , j , l
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible.
!
      IF ( (M==0) .OR. (N==0) ) RETURN
!
      DO j = 1 , N
         DO i = 1 , M
            Rwork((j-1)*M+i) = REAL(B(i,j))
         ENDDO
      ENDDO
!
      l = M*N + 1
      CALL SGEMM('N','N',M,N,M,ONE,A,Lda,Rwork,M,ZERO,Rwork(l),M)
      DO j = 1 , N
         DO i = 1 , M
            C(i,j) = Rwork(l+(j-1)*M+i-1)
         ENDDO
      ENDDO
!
      DO j = 1 , N
         DO i = 1 , M
            Rwork((j-1)*M+i) = AIMAG(B(i,j))
         ENDDO
      ENDDO
      CALL SGEMM('N','N',M,N,M,ONE,A,Lda,Rwork,M,ZERO,Rwork(l),M)
      DO j = 1 , N
         DO i = 1 , M
            C(i,j) = CMPLX(REAL(C(i,j)),Rwork(l+(j-1)*M+i-1))
         ENDDO
      ENDDO
!
!
!     End of CLARCM
!
      END SUBROUTINE CLARCM
