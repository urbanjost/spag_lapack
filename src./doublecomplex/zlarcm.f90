!*==zlarcm.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZLARCM copies all or part of a real two-dimensional array to a complex array.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLARCM + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlarcm.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlarcm.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlarcm.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLARCM( M, N, A, LDA, B, LDB, C, LDC, RWORK )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LDB, LDC, M, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), RWORK( * )
!       COMPLEX*16         B( LDB, * ), C( LDC, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLARCM performs a very simple matrix-matrix multiplication:
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
!>          A is DOUBLE PRECISION array, dimension (LDA, M)
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
!>          B is COMPLEX*16 array, dimension (LDB, N)
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
!>          C is COMPLEX*16 array, dimension (LDC, N)
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
!>          RWORK is DOUBLE PRECISION array, dimension (2*M*N)
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
!> \ingroup complex16OTHERauxiliary
!
!  =====================================================================
      SUBROUTINE ZLARCM(M,N,A,Lda,B,Ldb,C,Ldc,Rwork)
      USE F77KINDS                        
      USE S_DGEMM
      IMPLICIT NONE
!*--ZLARCM120
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D0 , ZERO = 0.0D0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER :: M
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER , INTENT(IN) :: Ldc
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
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
            Rwork((j-1)*M+i) = DBLE(B(i,j))
         ENDDO
      ENDDO
!
      l = M*N + 1
      CALL DGEMM('N','N',M,N,M,ONE,A,Lda,Rwork,M,ZERO,Rwork(l),M)
      DO j = 1 , N
         DO i = 1 , M
            C(i,j) = Rwork(l+(j-1)*M+i-1)
         ENDDO
      ENDDO
!
      DO j = 1 , N
         DO i = 1 , M
            Rwork((j-1)*M+i) = DIMAG(B(i,j))
         ENDDO
      ENDDO
      CALL DGEMM('N','N',M,N,M,ONE,A,Lda,Rwork,M,ZERO,Rwork(l),M)
      DO j = 1 , N
         DO i = 1 , M
            C(i,j) = DCMPLX(DBLE(C(i,j)),Rwork(l+(j-1)*M+i-1))
         ENDDO
      ENDDO
!
!
!     End of ZLARCM
!
      END SUBROUTINE ZLARCM
