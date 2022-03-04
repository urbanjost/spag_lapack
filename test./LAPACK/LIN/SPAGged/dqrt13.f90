!*==dqrt13.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b DQRT13
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DQRT13( SCALE, M, N, A, LDA, NORMA, ISEED )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, M, N, SCALE
!       DOUBLE PRECISION   NORMA
!       ..
!       .. Array Arguments ..
!       INTEGER            ISEED( 4 )
!       DOUBLE PRECISION   A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DQRT13 generates a full-rank matrix that may be scaled to have large
!> or small norm.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SCALE
!> \verbatim
!>          SCALE is INTEGER
!>          SCALE = 1: normally scaled matrix
!>          SCALE = 2: matrix scaled up
!>          SCALE = 3: matrix scaled down
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of A.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          The M-by-N matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.
!> \endverbatim
!>
!> \param[out] NORMA
!> \verbatim
!>          NORMA is DOUBLE PRECISION
!>          The one-norm of A.
!> \endverbatim
!>
!> \param[in,out] ISEED
!> \verbatim
!>          ISEED is integer array, dimension (4)
!>          Seed for random number generator
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
      SUBROUTINE DQRT13(Scale,M,N,A,Lda,Norma,Iseed)
      IMPLICIT NONE
!*--DQRT1395
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Lda , M , N , Scale
      DOUBLE PRECISION Norma
!     ..
!     .. Array Arguments ..
      INTEGER Iseed(4)
      DOUBLE PRECISION A(Lda,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
!     ..
!     .. Local Scalars ..
      INTEGER info , j
      DOUBLE PRECISION bignum , smlnum
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DASUM , DLAMCH , DLANGE
      EXTERNAL DASUM , DLAMCH , DLANGE
!     ..
!     .. External Subroutines ..
      EXTERNAL DLABAD , DLARNV , DLASCL
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC SIGN
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION dummy(1)
!     ..
!     .. Executable Statements ..
!
      IF ( M<=0 .OR. N<=0 ) RETURN
!
!     benign matrix
!
      DO j = 1 , N
         CALL DLARNV(2,Iseed,M,A(1,j))
         IF ( j<=M ) A(j,j) = A(j,j) + SIGN(DASUM(M,A(1,j),1),A(j,j))
      ENDDO
!
!     scaled versions
!
      IF ( Scale/=1 ) THEN
         Norma = DLANGE('Max',M,N,A,Lda,dummy)
         smlnum = DLAMCH('Safe minimum')
         bignum = ONE/smlnum
         CALL DLABAD(smlnum,bignum)
         smlnum = smlnum/DLAMCH('Epsilon')
         bignum = ONE/smlnum
!
         IF ( Scale==2 ) THEN
!
!           matrix scaled up
!
            CALL DLASCL('General',0,0,Norma,bignum,M,N,A,Lda,info)
         ELSEIF ( Scale==3 ) THEN
!
!           matrix scaled down
!
            CALL DLASCL('General',0,0,Norma,smlnum,M,N,A,Lda,info)
         ENDIF
      ENDIF
!
      Norma = DLANGE('One-norm',M,N,A,Lda,dummy)
!
!     End of DQRT13
!
      END SUBROUTINE DQRT13
