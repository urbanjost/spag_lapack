!*==zrzt01.f90  processed by SPAG 7.51RB at 17:35 on  4 Mar 2022
!> \brief \b zrzt01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION zrzt01( M, N, A, AF, LDA, TAU, WORK,
!                        LWORK )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * ), AF( LDA, * ), TAU( * ),
!      $                   WORK( LWORK )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> zrzt01 returns
!>      || A - R*Q || / ( M * eps * ||A|| )
!> for an upper trapezoidal A that was factored with ZTZRZF.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrices A and AF.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrices A and AF.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          The original upper trapezoidal M by N matrix A.
!> \endverbatim
!>
!> \param[in] AF
!> \verbatim
!>          AF is COMPLEX*16 array, dimension (LDA,N)
!>          The output of ZTZRZF for input matrix A.
!>          The lower triangle is not referenced.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the arrays A and AF.
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is COMPLEX*16 array, dimension (M)
!>          Details of the  Householder transformations as returned by
!>          ZTZRZF.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (LWORK)
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The length of the array WORK.  LWORK >= m*n + m.
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
!> \ingroup complex16_lin
!
!  =====================================================================
      DOUBLE PRECISION FUNCTION zrzt01(M,N,A,Af,Lda,Tau,Work,Lwork)
      IMPLICIT NONE
!*--ZRZT01101
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Lda , Lwork , M , N
!     ..
!     .. Array Arguments ..
      COMPLEX*16 A(Lda,*) , Af(Lda,*) , Tau(*) , Work(Lwork)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
!     ..
!     .. Local Scalars ..
      INTEGER i , info , j
      DOUBLE PRECISION norma
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION rwork(1)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH , ZLANGE
EXTERNAL    DLAMCH,   ZLANGE
!     ..
!     .. External Subroutines ..
EXTERNAL    ZAXPY,   ZLASET,   ZUNMRZ
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DBLE , DCMPLX , MAX
!     ..
!     .. Executable Statements ..
!
      zrzt01 = ZERO
!
      IF ( Lwork<M*N+M ) THEN
         CALL xerbla('zrzt01',8)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( M<=0 .OR. N<=0 ) RETURN
!
      norma = ZLANGE('One-norm',M,N,A,Lda,rwork)
!
!     Copy upper triangle R
!
      CALL ZLASET('Full',M,N,DCMPLX(ZERO),DCMPLX(ZERO),Work,M)
      DO j = 1 , M
         DO i = 1 , j
            Work((j-1)*M+i) = Af(i,j)
         ENDDO
      ENDDO
!
!     R = R * P(1) * ... *P(m)
!
      CALL ZUNMRZ('Right','No tranpose',M,N,M,N-M,Af,Lda,Tau,Work,M,    &
     &            Work(M*N+1),Lwork-M*N,info)
!
!     R = R - A
!
      DO i = 1 , N
         CALL ZAXPY(M,DCMPLX(-ONE),A(1,i),1,Work((i-1)*M+1),1)
      ENDDO
!
      zrzt01 = ZLANGE('One-norm',M,N,Work,M,rwork)
!
      zrzt01 = zrzt01/(DLAMCH('Epsilon')*DBLE(MAX(M,N)))
      IF ( norma/=ZERO ) zrzt01 = zrzt01/norma
!
!
!     End of zrzt01
!
      END FUNCTION zrzt01
