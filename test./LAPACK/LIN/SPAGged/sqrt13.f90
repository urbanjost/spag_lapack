!*==sqrt13.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b SQRT13
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SQRT13( SCALE, M, N, A, LDA, NORMA, ISEED )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, M, N, SCALE
!       REAL               NORMA
!       ..
!       .. Array Arguments ..
!       INTEGER            ISEED( 4 )
!       REAL               A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SQRT13 generates a full-rank matrix that may be scaled to have large
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
!>          A is REAL array, dimension (LDA,N)
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
!>          NORMA is REAL
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
!> \ingroup single_lin
!
!  =====================================================================
      SUBROUTINE SQRT13(Scale,M,N,A,Lda,Norma,Iseed)
      IMPLICIT NONE
!*--SQRT1395
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Lda , M , N , Scale
      REAL Norma
!     ..
!     .. Array Arguments ..
      INTEGER Iseed(4)
      REAL A(Lda,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ONE
      PARAMETER (ONE=1.0E0)
!     ..
!     .. Local Scalars ..
      INTEGER info , j
      REAL bignum , smlnum
!     ..
!     .. External Functions ..
      REAL SASUM , SLAMCH , SLANGE
      EXTERNAL SASUM , SLAMCH , SLANGE
!     ..
!     .. External Subroutines ..
      EXTERNAL SLABAD , SLARNV , SLASCL
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC SIGN
!     ..
!     .. Local Arrays ..
      REAL dummy(1)
!     ..
!     .. Executable Statements ..
!
      IF ( M<=0 .OR. N<=0 ) RETURN
!
!     benign matrix
!
      DO j = 1 , N
         CALL SLARNV(2,Iseed,M,A(1,j))
         IF ( j<=M ) A(j,j) = A(j,j) + SIGN(SASUM(M,A(1,j),1),A(j,j))
      ENDDO
!
!     scaled versions
!
      IF ( Scale/=1 ) THEN
         Norma = SLANGE('Max',M,N,A,Lda,dummy)
         smlnum = SLAMCH('Safe minimum')
         bignum = ONE/smlnum
         CALL SLABAD(smlnum,bignum)
         smlnum = smlnum/SLAMCH('Epsilon')
         bignum = ONE/smlnum
!
         IF ( Scale==2 ) THEN
!
!           matrix scaled up
!
            CALL SLASCL('General',0,0,Norma,bignum,M,N,A,Lda,info)
         ELSEIF ( Scale==3 ) THEN
!
!           matrix scaled down
!
            CALL SLASCL('General',0,0,Norma,smlnum,M,N,A,Lda,info)
         ENDIF
      ENDIF
!
      Norma = SLANGE('One-norm',M,N,A,Lda,dummy)
!
!     End of SQRT13
!
      END SUBROUTINE SQRT13
