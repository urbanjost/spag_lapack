!*==dlatb5.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b dlatb5
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE dlatb5( PATH, IMAT, N, TYPE, KL, KU, ANORM, MODE,
!                          CNDNUM, DIST )
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION   ANORM, CNDNUM
!       INTEGER            IMAT, KL, KU, MODE, N
!       CHARACTER          DIST, TYPE
!       CHARACTER*3        PATH
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> dlatb5 sets parameters for the matrix generator based on the type
!> of matrix to be generated.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] PATH
!> \verbatim
!>          PATH is CHARACTER*3
!>          The LAPACK path name.
!> \endverbatim
!>
!> \param[in] IMAT
!> \verbatim
!>          IMAT is INTEGER
!>          An integer key describing which matrix to generate for this
!>          path.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of rows and columns in the matrix to be generated.
!> \endverbatim
!>
!> \param[out] TYPE
!> \verbatim
!>          TYPE is CHARACTER*1
!>          The type of the matrix to be generated:
!>          = 'S':  symmetric matrix
!>          = 'P':  symmetric positive (semi)definite matrix
!>          = 'N':  nonsymmetric matrix
!> \endverbatim
!>
!> \param[out] KL
!> \verbatim
!>          KL is INTEGER
!>          The lower band width of the matrix to be generated.
!> \endverbatim
!>
!> \param[out] KU
!> \verbatim
!>          KU is INTEGER
!>          The upper band width of the matrix to be generated.
!> \endverbatim
!>
!> \param[out] ANORM
!> \verbatim
!>          ANORM is DOUBLE PRECISION
!>          The desired norm of the matrix to be generated.  The diagonal
!>          matrix of singular values or eigenvalues is scaled by this
!>          value.
!> \endverbatim
!>
!> \param[out] MODE
!> \verbatim
!>          MODE is INTEGER
!>          A key indicating how to choose the vector of eigenvalues.
!> \endverbatim
!>
!> \param[out] CNDNUM
!> \verbatim
!>          CNDNUM is DOUBLE PRECISION
!>          The desired condition number.
!> \endverbatim
!>
!> \param[out] DIST
!> \verbatim
!>          DIST is CHARACTER*1
!>          The type of distribution to be used by the random number
!>          generator.
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
      SUBROUTINE dlatb5(Path,Imat,N,Type,Kl,Ku,Anorm,Mode,Cndnum,Dist)
      IMPLICIT NONE
!*--DLATB5117
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION Anorm , Cndnum
      INTEGER Imat , Kl , Ku , Mode , N
      CHARACTER Dist , Type
      CHARACTER*3 Path
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION SHRINK , TENTH
      PARAMETER (SHRINK=0.25D0,TENTH=0.1D+0)
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D+0)
      DOUBLE PRECISION TWO
      PARAMETER (TWO=2.0D+0)
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION badc1 , badc2 , eps , large , small
      LOGICAL first
      CHARACTER*2 c2
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH
EXTERNAL    DLAMCH
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , SQRT
!     ..
!     .. External Subroutines ..
EXTERNAL    DLABAD
!     ..
!     .. Save statement ..
      SAVE eps , small , large , badc1 , badc2 , first
!     ..
!     .. Data statements ..
      DATA first/.TRUE./
!     ..
!     .. Executable Statements ..
!
!     Set some constants for use in the subroutine.
!
      IF ( first ) THEN
         first = .FALSE.
         eps = DLAMCH('Precision')
         badc2 = TENTH/eps
         badc1 = SQRT(badc2)
         small = DLAMCH('Safe minimum')
         large = ONE/small
!
!        If it looks like we're on a Cray, take the square root of
!        SMALL and LARGE to avoid overflow and underflow problems.
!
         CALL DLABAD(small,large)
         small = SHRINK*(small/eps)
         large = ONE/small
      ENDIF
!
      c2 = Path(2:3)
!
!     Set some parameters
!
      Dist = 'S'
      Mode = 3
!
!     Set TYPE, the type of matrix to be generated.
!
      Type = c2(1:1)
!
!     Set the lower and upper bandwidths.
!
      IF ( Imat==1 ) THEN
         Kl = 0
      ELSE
         Kl = MAX(N-1,0)
      ENDIF
      Ku = Kl
!
!     Set the condition number and norm.etc
!
      IF ( Imat==3 ) THEN
         Cndnum = 1.0D12
         Mode = 2
      ELSEIF ( Imat==4 ) THEN
         Cndnum = 1.0D12
         Mode = 1
      ELSEIF ( Imat==5 ) THEN
         Cndnum = 1.0D12
         Mode = 3
      ELSEIF ( Imat==6 ) THEN
         Cndnum = badc1
      ELSEIF ( Imat==7 ) THEN
         Cndnum = badc2
      ELSE
         Cndnum = TWO
      ENDIF
!
      IF ( Imat==8 ) THEN
         Anorm = small
      ELSEIF ( Imat==9 ) THEN
         Anorm = large
      ELSE
         Anorm = ONE
      ENDIF
!
      IF ( N<=1 ) Cndnum = ONE
!
!
!     End of dlatb5
!
      END SUBROUTINE dlatb5
