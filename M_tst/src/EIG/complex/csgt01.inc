!*==csgt01.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b csgt01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE csgt01( ITYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D,
!                          WORK, RWORK, RESULT )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            ITYPE, LDA, LDB, LDZ, M, N
!       ..
!       .. Array Arguments ..
!       REAL               D( * ), RESULT( * ), RWORK( * )
!       COMPLEX            A( LDA, * ), B( LDB, * ), WORK( * ),
!      $                   Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> csgt01 checks a decomposition of the form
!>
!>    A Z   =  B Z D or
!>    A B Z =  Z D or
!>    B A Z =  Z D
!>
!> where A is a Hermitian matrix, B is Hermitian positive definite,
!> Z is unitary, and D is diagonal.
!>
!> One of the following test ratios is computed:
!>
!> ITYPE = 1:  RESULT(1) = | A Z - B Z D | / ( |A| |Z| n ulp )
!>
!> ITYPE = 2:  RESULT(1) = | A B Z - Z D | / ( |A| |Z| n ulp )
!>
!> ITYPE = 3:  RESULT(1) = | B A Z - Z D | / ( |A| |Z| n ulp )
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ITYPE
!> \verbatim
!>          ITYPE is INTEGER
!>          The form of the Hermitian generalized eigenproblem.
!>          = 1:  A*z = (lambda)*B*z
!>          = 2:  A*B*z = (lambda)*z
!>          = 3:  B*A*z = (lambda)*z
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the upper or lower triangular part of the
!>          Hermitian matrices A and B is stored.
!>          = 'U':  Upper triangular
!>          = 'L':  Lower triangular
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of eigenvalues found.  M >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA, N)
!>          The original Hermitian matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is COMPLEX array, dimension (LDB, N)
!>          The original Hermitian positive definite matrix B.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,N).
!> \endverbatim
!>
!> \param[in] Z
!> \verbatim
!>          Z is COMPLEX array, dimension (LDZ, M)
!>          The computed eigenvectors of the generalized eigenproblem.
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>          The leading dimension of the array Z.  LDZ >= max(1,N).
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is REAL array, dimension (M)
!>          The computed eigenvalues of the generalized eigenproblem.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (N*N)
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (N)
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is REAL array, dimension (1)
!>          The test ratio as described above.
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
!> \ingroup complex_eig
!
!  =====================================================================
      SUBROUTINE csgt01(Itype,Uplo,N,M,A,Lda,B,Ldb,Z,Ldz,D,Work,Rwork,  &
     &                  Result)
      IMPLICIT NONE
!*--CSGT01156
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER Itype , Lda , Ldb , Ldz , M , N
!     ..
!     .. Array Arguments ..
      REAL D(*) , Result(*) , Rwork(*)
      COMPLEX A(Lda,*) , B(Ldb,*) , Work(*) , Z(Ldz,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
      COMPLEX CZERO , CONE
      PARAMETER (CZERO=(0.0E+0,0.0E+0),CONE=(1.0E+0,0.0E+0))
!     ..
!     .. Local Scalars ..
      INTEGER i
      REAL anorm , ulp
!     ..
!     .. External Functions ..
      REAL CLANGE , CLANHE , SLAMCH
      EXTERNAL CLANGE , CLANHE , SLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL CHEMM , CSSCAL
!     ..
!     .. Executable Statements ..
!
      Result(1) = ZERO
      IF ( N<=0 ) RETURN
!
      ulp = SLAMCH('Epsilon')
!
!     Compute product of 1-norms of A and Z.
!
      anorm = CLANHE('1',Uplo,N,A,Lda,Rwork)*CLANGE('1',N,M,Z,Ldz,Rwork)
      IF ( anorm==ZERO ) anorm = ONE
!
      IF ( Itype==1 ) THEN
!
!        Norm of AZ - BZD
!
         CALL CHEMM('Left',Uplo,N,M,CONE,A,Lda,Z,Ldz,CZERO,Work,N)
         DO i = 1 , M
            CALL CSSCAL(N,D(i),Z(1,i),1)
         ENDDO
         CALL CHEMM('Left',Uplo,N,M,CONE,B,Ldb,Z,Ldz,-CONE,Work,N)
!
         Result(1) = (CLANGE('1',N,M,Work,N,Rwork)/anorm)/(N*ulp)
!
      ELSEIF ( Itype==2 ) THEN
!
!        Norm of ABZ - ZD
!
         CALL CHEMM('Left',Uplo,N,M,CONE,B,Ldb,Z,Ldz,CZERO,Work,N)
         DO i = 1 , M
            CALL CSSCAL(N,D(i),Z(1,i),1)
         ENDDO
         CALL CHEMM('Left',Uplo,N,M,CONE,A,Lda,Work,N,-CONE,Z,Ldz)
!
         Result(1) = (CLANGE('1',N,M,Z,Ldz,Rwork)/anorm)/(N*ulp)
!
      ELSEIF ( Itype==3 ) THEN
!
!        Norm of BAZ - ZD
!
         CALL CHEMM('Left',Uplo,N,M,CONE,A,Lda,Z,Ldz,CZERO,Work,N)
         DO i = 1 , M
            CALL CSSCAL(N,D(i),Z(1,i),1)
         ENDDO
         CALL CHEMM('Left',Uplo,N,M,CONE,B,Ldb,Work,N,-CONE,Z,Ldz)
!
         Result(1) = (CLANGE('1',N,M,Z,Ldz,Rwork)/anorm)/(N*ulp)
      ENDIF
!
!
!     End of csgt01
!
      END SUBROUTINE csgt01
