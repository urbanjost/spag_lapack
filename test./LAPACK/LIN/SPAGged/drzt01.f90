!*==drzt01.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b DRZT01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION DRZT01( M, N, A, AF, LDA, TAU, WORK,
!                        LWORK )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), AF( LDA, * ), TAU( * ),
!      $                   WORK( LWORK )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DRZT01 returns
!>      || A - R*Q || / ( M * eps * ||A|| )
!> for an upper trapezoidal A that was factored with DTZRZF.
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
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          The original upper trapezoidal M by N matrix A.
!> \endverbatim
!>
!> \param[in] AF
!> \verbatim
!>          AF is DOUBLE PRECISION array, dimension (LDA,N)
!>          The output of DTZRZF for input matrix A.
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
!>          TAU is DOUBLE PRECISION array, dimension (M)
!>          Details of the Householder transformations as returned by
!>          DTZRZF.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (LWORK)
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The length of the array WORK.  LWORK >= m*n + m*nb.
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
!> \ingroup double_lin
!
!  =====================================================================
      DOUBLE PRECISION FUNCTION DRZT01(M,N,A,Af,Lda,Tau,Work,Lwork)
      IMPLICIT NONE
!*--DRZT01101
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
      DOUBLE PRECISION A(Lda,*) , Af(Lda,*) , Tau(*) , Work(Lwork)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
!     ..
!     .. Local Scalars ..
      INTEGER i , info , j
      DOUBLE PRECISION norma
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION rwork(1)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH , DLANGE
      EXTERNAL DLAMCH , DLANGE
!     ..
!     .. External Subroutines ..
      EXTERNAL DAXPY , DLASET , DORMRZ , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DBLE , MAX
!     ..
!     .. Executable Statements ..
!
      DRZT01 = ZERO
!
      IF ( Lwork<M*N+M ) THEN
         CALL XERBLA('DRZT01',8)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( M<=0 .OR. N<=0 ) RETURN
!
      norma = DLANGE('One-norm',M,N,A,Lda,rwork)
!
!     Copy upper triangle R
!
      CALL DLASET('Full',M,N,ZERO,ZERO,Work,M)
      DO j = 1 , M
         DO i = 1 , j
            Work((j-1)*M+i) = Af(i,j)
         ENDDO
      ENDDO
!
!     R = R * P(1) * ... *P(m)
!
      CALL DORMRZ('Right','No tranpose',M,N,M,N-M,Af,Lda,Tau,Work,M,    &
     &            Work(M*N+1),Lwork-M*N,info)
!
!     R = R - A
!
      DO i = 1 , N
         CALL DAXPY(M,-ONE,A(1,i),1,Work((i-1)*M+1),1)
      ENDDO
!
      DRZT01 = DLANGE('One-norm',M,N,Work,M,rwork)
!
      DRZT01 = DRZT01/(DLAMCH('Epsilon')*DBLE(MAX(M,N)))
      IF ( norma/=ZERO ) DRZT01 = DRZT01/norma
!
!
!     End of DRZT01
!
      END FUNCTION DRZT01
