!*==zsgt01.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b ZSGT01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZSGT01( ITYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D,
!                          WORK, RWORK, RESULT )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            ITYPE, LDA, LDB, LDZ, M, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   D( * ), RESULT( * ), RWORK( * )
!       COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * ),
!      $                   Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CDGT01 checks a decomposition of the form
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
!>          A is COMPLEX*16 array, dimension (LDA, N)
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
!>          B is COMPLEX*16 array, dimension (LDB, N)
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
!>          Z is COMPLEX*16 array, dimension (LDZ, M)
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
!>          D is DOUBLE PRECISION array, dimension (M)
!>          The computed eigenvalues of the generalized eigenproblem.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (N*N)
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is DOUBLE PRECISION array, dimension (1)
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
!> \ingroup complex16_eig
!
!  =====================================================================
      SUBROUTINE ZSGT01(Itype,Uplo,N,M,A,Lda,B,Ldb,Z,Ldz,D,Work,Rwork,  &
     &                  Result)
      IMPLICIT NONE
!*--ZSGT01156
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
      DOUBLE PRECISION D(*) , Result(*) , Rwork(*)
      COMPLEX*16 A(Lda,*) , B(Ldb,*) , Work(*) , Z(Ldz,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
      COMPLEX*16 CZERO , CONE
      PARAMETER (CZERO=(0.0D+0,0.0D+0),CONE=(1.0D+0,0.0D+0))
!     ..
!     .. Local Scalars ..
      INTEGER i
      DOUBLE PRECISION anorm , ulp
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH , ZLANGE , ZLANHE
      EXTERNAL DLAMCH , ZLANGE , ZLANHE
!     ..
!     .. External Subroutines ..
      EXTERNAL ZDSCAL , ZHEMM
!     ..
!     .. Executable Statements ..
!
      Result(1) = ZERO
      IF ( N<=0 ) RETURN
!
      ulp = DLAMCH('Epsilon')
!
!     Compute product of 1-norms of A and Z.
!
      anorm = ZLANHE('1',Uplo,N,A,Lda,Rwork)*ZLANGE('1',N,M,Z,Ldz,Rwork)
      IF ( anorm==ZERO ) anorm = ONE
!
      IF ( Itype==1 ) THEN
!
!        Norm of AZ - BZD
!
         CALL ZHEMM('Left',Uplo,N,M,CONE,A,Lda,Z,Ldz,CZERO,Work,N)
         DO i = 1 , M
            CALL ZDSCAL(N,D(i),Z(1,i),1)
         ENDDO
         CALL ZHEMM('Left',Uplo,N,M,CONE,B,Ldb,Z,Ldz,-CONE,Work,N)
!
         Result(1) = (ZLANGE('1',N,M,Work,N,Rwork)/anorm)/(N*ulp)
!
      ELSEIF ( Itype==2 ) THEN
!
!        Norm of ABZ - ZD
!
         CALL ZHEMM('Left',Uplo,N,M,CONE,B,Ldb,Z,Ldz,CZERO,Work,N)
         DO i = 1 , M
            CALL ZDSCAL(N,D(i),Z(1,i),1)
         ENDDO
         CALL ZHEMM('Left',Uplo,N,M,CONE,A,Lda,Work,N,-CONE,Z,Ldz)
!
         Result(1) = (ZLANGE('1',N,M,Z,Ldz,Rwork)/anorm)/(N*ulp)
!
      ELSEIF ( Itype==3 ) THEN
!
!        Norm of BAZ - ZD
!
         CALL ZHEMM('Left',Uplo,N,M,CONE,A,Lda,Z,Ldz,CZERO,Work,N)
         DO i = 1 , M
            CALL ZDSCAL(N,D(i),Z(1,i),1)
         ENDDO
         CALL ZHEMM('Left',Uplo,N,M,CONE,B,Ldb,Work,N,-CONE,Z,Ldz)
!
         Result(1) = (ZLANGE('1',N,M,Z,Ldz,Rwork)/anorm)/(N*ulp)
      ENDIF
!
!
!     End of CDGT01
!
      END SUBROUTINE ZSGT01
