!*==dsgt01.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b DSGT01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSGT01( ITYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D,
!                          WORK, RESULT )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            ITYPE, LDA, LDB, LDZ, M, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), D( * ), RESULT( * ),
!      $                   WORK( * ), Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DDGT01 checks a decomposition of the form
!>
!>    A Z   =  B Z D or
!>    A B Z =  Z D or
!>    B A Z =  Z D
!>
!> where A is a symmetric matrix, B is
!> symmetric positive definite, Z is orthogonal, and D is diagonal.
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
!>          The form of the symmetric generalized eigenproblem.
!>          = 1:  A*z = (lambda)*B*z
!>          = 2:  A*B*z = (lambda)*z
!>          = 3:  B*A*z = (lambda)*z
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the upper or lower triangular part of the
!>          symmetric matrices A and B is stored.
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
!>          The number of eigenvalues found.  0 <= M <= N.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA, N)
!>          The original symmetric matrix A.
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
!>          B is DOUBLE PRECISION array, dimension (LDB, N)
!>          The original symmetric positive definite matrix B.
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
!>          Z is DOUBLE PRECISION array, dimension (LDZ, M)
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
!>          WORK is DOUBLE PRECISION array, dimension (N*N)
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
!> \ingroup double_eig
!
!  =====================================================================
      SUBROUTINE DSGT01(Itype,Uplo,N,M,A,Lda,B,Ldb,Z,Ldz,D,Work,Result)
      IMPLICIT NONE
!*--DSGT01149
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
      DOUBLE PRECISION A(Lda,*) , B(Ldb,*) , D(*) , Result(*) , Work(*) &
     &                 , Z(Ldz,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
!     ..
!     .. Local Scalars ..
      INTEGER i
      DOUBLE PRECISION anorm , ulp
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH , DLANGE , DLANSY
      EXTERNAL DLAMCH , DLANGE , DLANSY
!     ..
!     .. External Subroutines ..
      EXTERNAL DSCAL , DSYMM
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
      anorm = DLANSY('1',Uplo,N,A,Lda,Work)*DLANGE('1',N,M,Z,Ldz,Work)
      IF ( anorm==ZERO ) anorm = ONE
!
      IF ( Itype==1 ) THEN
!
!        Norm of AZ - BZD
!
         CALL DSYMM('Left',Uplo,N,M,ONE,A,Lda,Z,Ldz,ZERO,Work,N)
         DO i = 1 , M
            CALL DSCAL(N,D(i),Z(1,i),1)
         ENDDO
         CALL DSYMM('Left',Uplo,N,M,ONE,B,Ldb,Z,Ldz,-ONE,Work,N)
!
         Result(1) = (DLANGE('1',N,M,Work,N,Work)/anorm)/(N*ulp)
!
      ELSEIF ( Itype==2 ) THEN
!
!        Norm of ABZ - ZD
!
         CALL DSYMM('Left',Uplo,N,M,ONE,B,Ldb,Z,Ldz,ZERO,Work,N)
         DO i = 1 , M
            CALL DSCAL(N,D(i),Z(1,i),1)
         ENDDO
         CALL DSYMM('Left',Uplo,N,M,ONE,A,Lda,Work,N,-ONE,Z,Ldz)
!
         Result(1) = (DLANGE('1',N,M,Z,Ldz,Work)/anorm)/(N*ulp)
!
      ELSEIF ( Itype==3 ) THEN
!
!        Norm of BAZ - ZD
!
         CALL DSYMM('Left',Uplo,N,M,ONE,A,Lda,Z,Ldz,ZERO,Work,N)
         DO i = 1 , M
            CALL DSCAL(N,D(i),Z(1,i),1)
         ENDDO
         CALL DSYMM('Left',Uplo,N,M,ONE,B,Ldb,Work,N,-ONE,Z,Ldz)
!
         Result(1) = (DLANGE('1',N,M,Z,Ldz,Work)/anorm)/(N*ulp)
      ENDIF
!
!
!     End of DDGT01
!
      END SUBROUTINE DSGT01
