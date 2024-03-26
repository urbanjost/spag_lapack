!*==zqrt13.f90  processed by SPAG 7.51RB at 17:35 on  4 Mar 2022
!> \brief \b zqrt13
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE zqrt13( SCALE, M, N, A, LDA, NORMA, ISEED )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, M, N, SCALE
!       DOUBLE PRECISION   NORMA
!       ..
!       .. Array Arguments ..
!       INTEGER            ISEED( 4 )
!       COMPLEX*16         A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> zqrt13 generates a full-rank matrix that may be scaled to have large
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
!>          A is COMPLEX*16 array, dimension (LDA,N)
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
!> \ingroup complex16_lin
!
!  =====================================================================
      SUBROUTINE zqrt13(Scale,M,N,A,Lda,Norma,Iseed)
      IMPLICIT NONE
!*--ZQRT1395
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
      COMPLEX*16 A(Lda,*)
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
      DOUBLE PRECISION DLAMCH , DZASUM , ZLANGE
EXTERNAL    DLAMCH,   DZASUM,   ZLANGE
!     ..
!     .. External Subroutines ..
EXTERNAL    DLABAD,   ZLARNV,   ZLASCL
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DBLE , DCMPLX , SIGN
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
         CALL ZLARNV(2,Iseed,M,A(1,j))
         IF ( j<=M ) A(j,j) = A(j,j)                                    &
     &                        + DCMPLX(SIGN(DZASUM(M,A(1,j),1),DBLE     &
     &                        (A(j,j))))
      ENDDO
!
!     scaled versions
!
      IF ( Scale/=1 ) THEN
         Norma = ZLANGE('Max',M,N,A,Lda,dummy)
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
            CALL ZLASCL('General',0,0,Norma,bignum,M,N,A,Lda,info)
         ELSEIF ( Scale==3 ) THEN
!
!           matrix scaled down
!
            CALL ZLASCL('General',0,0,Norma,smlnum,M,N,A,Lda,info)
         ENDIF
      ENDIF
!
      Norma = ZLANGE('One-norm',M,N,A,Lda,dummy)
!
!     End of zqrt13
!
      END SUBROUTINE zqrt13
