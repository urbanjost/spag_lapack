!*==zhet22.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b zhet22
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE zhet22( ITYPE, UPLO, N, M, KBAND, A, LDA, D, E, U, LDU,
!                          V, LDV, TAU, WORK, RWORK, RESULT )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            ITYPE, KBAND, LDA, LDU, LDV, M, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   D( * ), E( * ), RESULT( 2 ), RWORK( * )
!       COMPLEX*16         A( LDA, * ), TAU( * ), U( LDU, * ),
!      $                   V( LDV, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>      zhet22  generally checks a decomposition of the form
!>
!>              A U = U S
!>
!>      where A is complex Hermitian, the columns of U are orthonormal,
!>      and S is diagonal (if KBAND=0) or symmetric tridiagonal (if
!>      KBAND=1).  If ITYPE=1, then U is represented as a dense matrix,
!>      otherwise the U is expressed as a product of Householder
!>      transformations, whose vectors are stored in the array "V" and
!>      whose scaling constants are in "TAU"; we shall use the letter
!>      "V" to refer to the product of Householder transformations
!>      (which should be equal to U).
!>
!>      Specifically, if ITYPE=1, then:
!>
!>              RESULT(1) = | U**H A U - S | / ( |A| m ulp ) and
!>              RESULT(2) = | I - U**H U | / ( m ulp )
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \verbatim
!>  ITYPE   INTEGER
!>          Specifies the type of tests to be performed.
!>          1: U expressed as a dense orthogonal matrix:
!>             RESULT(1) = | A - U S U**H | / ( |A| n ulp )   *and
!>             RESULT(2) = | I - U U**H | / ( n ulp )
!>
!>  UPLO    CHARACTER
!>          If UPLO='U', the upper triangle of A will be used and the
!>          (strictly) lower triangle will not be referenced.  If
!>          UPLO='L', the lower triangle of A will be used and the
!>          (strictly) upper triangle will not be referenced.
!>          Not modified.
!>
!>  N       INTEGER
!>          The size of the matrix.  If it is zero, zhet22 does nothing.
!>          It must be at least zero.
!>          Not modified.
!>
!>  M       INTEGER
!>          The number of columns of U.  If it is zero, zhet22 does
!>          nothing.  It must be at least zero.
!>          Not modified.
!>
!>  KBAND   INTEGER
!>          The bandwidth of the matrix.  It may only be zero or one.
!>          If zero, then S is diagonal, and E is not referenced.  If
!>          one, then S is symmetric tri-diagonal.
!>          Not modified.
!>
!>  A       COMPLEX*16 array, dimension (LDA , N)
!>          The original (unfactored) matrix.  It is assumed to be
!>          symmetric, and only the upper (UPLO='U') or only the lower
!>          (UPLO='L') will be referenced.
!>          Not modified.
!>
!>  LDA     INTEGER
!>          The leading dimension of A.  It must be at least 1
!>          and at least N.
!>          Not modified.
!>
!>  D       DOUBLE PRECISION array, dimension (N)
!>          The diagonal of the (symmetric tri-) diagonal matrix.
!>          Not modified.
!>
!>  E       DOUBLE PRECISION array, dimension (N)
!>          The off-diagonal of the (symmetric tri-) diagonal matrix.
!>          E(1) is ignored, E(2) is the (1,2) and (2,1) element, etc.
!>          Not referenced if KBAND=0.
!>          Not modified.
!>
!>  U       COMPLEX*16 array, dimension (LDU, N)
!>          If ITYPE=1, this contains the orthogonal matrix in
!>          the decomposition, expressed as a dense matrix.
!>          Not modified.
!>
!>  LDU     INTEGER
!>          The leading dimension of U.  LDU must be at least N and
!>          at least 1.
!>          Not modified.
!>
!>  V       COMPLEX*16 array, dimension (LDV, N)
!>          If ITYPE=2 or 3, the lower triangle of this array contains
!>          the Householder vectors used to describe the orthogonal
!>          matrix in the decomposition.  If ITYPE=1, then it is not
!>          referenced.
!>          Not modified.
!>
!>  LDV     INTEGER
!>          The leading dimension of V.  LDV must be at least N and
!>          at least 1.
!>          Not modified.
!>
!>  TAU     COMPLEX*16 array, dimension (N)
!>          If ITYPE >= 2, then TAU(j) is the scalar factor of
!>          v(j) v(j)**H in the Householder transformation H(j) of
!>          the product  U = H(1)...H(n-2)
!>          If ITYPE < 2, then TAU is not referenced.
!>          Not modified.
!>
!>  WORK    COMPLEX*16 array, dimension (2*N**2)
!>          Workspace.
!>          Modified.
!>
!>  RWORK   DOUBLE PRECISION array, dimension (N)
!>          Workspace.
!>          Modified.
!>
!>  RESULT  DOUBLE PRECISION array, dimension (2)
!>          The values computed by the two tests described above.  The
!>          values are currently limited to 1/ulp, to avoid overflow.
!>          RESULT(1) is always modified.  RESULT(2) is modified only
!>          if LDU is at least N.
!>          Modified.
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
      SUBROUTINE zhet22(Itype,Uplo,N,M,Kband,A,Lda,D,E,U,Ldu,V,Ldv,Tau, &
     &                  Work,Rwork,Result)
      IMPLICIT NONE
!*--ZHET22165
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER Itype , Kband , Lda , Ldu , Ldv , M , N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION D(*) , E(*) , Result(2) , Rwork(*)
      COMPLEX*16 A(Lda,*) , Tau(*) , U(Ldu,*) , V(Ldv,*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
      COMPLEX*16 CZERO , CONE
      PARAMETER (CZERO=(0.0D0,0.0D0),CONE=(1.0D0,0.0D0))
!     ..
!     .. Local Scalars ..
      INTEGER j , jj , jj1 , jj2 , nn , nnp1
      DOUBLE PRECISION anorm , ulp , unfl , wnorm
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH , ZLANHE
      EXTERNAL DLAMCH , ZLANHE
!     ..
!     .. External Subroutines ..
      EXTERNAL ZGEMM , ZHEMM 
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DBLE , MAX , MIN
!     ..
!     .. Executable Statements ..
!
      Result(1) = ZERO
      Result(2) = ZERO
      IF ( N<=0 .OR. M<=0 ) RETURN
!
      unfl = DLAMCH('Safe minimum')
      ulp = DLAMCH('Precision')
!
!     Do Test 1
!
!     Norm of A:
!
      anorm = MAX(ZLANHE('1',Uplo,N,A,Lda,Rwork),unfl)
!
!     Compute error matrix:
!
!     ITYPE=1: error = U**H A U - S
!
      CALL ZHEMM('L',Uplo,N,M,CONE,A,Lda,U,Ldu,CZERO,Work,N)
      nn = N*N
      nnp1 = nn + 1
      CALL ZGEMM('C','N',M,M,N,CONE,U,Ldu,Work,N,CZERO,Work(nnp1),N)
      DO j = 1 , M
         jj = nn + (j-1)*N + j
         Work(jj) = Work(jj) - D(j)
      ENDDO
      IF ( Kband==1 .AND. N>1 ) THEN
         DO j = 2 , M
            jj1 = nn + (j-1)*N + j - 1
            jj2 = nn + (j-2)*N + j
            Work(jj1) = Work(jj1) - E(j-1)
            Work(jj2) = Work(jj2) - E(j-1)
         ENDDO
      ENDIF
      wnorm = ZLANHE('1',Uplo,M,Work(nnp1),N,Rwork)
!
      IF ( anorm>wnorm ) THEN
         Result(1) = (wnorm/anorm)/(M*ulp)
      ELSEIF ( anorm<ONE ) THEN
         Result(1) = (MIN(wnorm,M*anorm)/anorm)/(M*ulp)
      ELSE
         Result(1) = MIN(wnorm/anorm,DBLE(M))/(M*ulp)
      ENDIF
!
!     Do Test 2
!
!     Compute  U**H U - I
!
      IF ( Itype==1 ) CALL zunt01('Columns',N,M,U,Ldu,Work,2*N*N,Rwork, &
     &                            Result(2))
!
!
!     End of zhet22
!
      END SUBROUTINE zhet22
