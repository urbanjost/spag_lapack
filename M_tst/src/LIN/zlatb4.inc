!*==zlatb4.f90  processed by SPAG 7.51RB at 17:35 on  4 Mar 2022
!> \brief \b zlatb4
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE zlatb4( PATH, IMAT, M, N, TYPE, KL, KU, ANORM, MODE,
!                          CNDNUM, DIST )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIST, TYPE
!       CHARACTER*3        PATH
!       INTEGER            IMAT, KL, KU, M, MODE, N
!       DOUBLE PRECISION   ANORM, CNDNUM
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> zlatb4 sets parameters for the matrix generator based on the type of
!> matrix to be generated.
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
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows in the matrix to be generated.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns in the matrix to be generated.
!> \endverbatim
!>
!> \param[out] TYPE
!> \verbatim
!>          TYPE is CHARACTER*1
!>          The type of the matrix to be generated:
!>          = 'S':  symmetric matrix
!>          = 'H':  Hermitian matrix
!>          = 'P':  Hermitian positive (semi)definite matrix
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
!> \ingroup complex16_lin
!
!  =====================================================================
      SUBROUTINE zlatb4(Path,Imat,M,N,Type,Kl,Ku,Anorm,Mode,Cndnum,Dist)
      IMPLICIT NONE
!*--ZLATB4124
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Dist , Type
      CHARACTER*3 Path
      INTEGER Imat , Kl , Ku , M , Mode , N
      DOUBLE PRECISION Anorm , Cndnum
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
      LOGICAL first
      CHARACTER*2 c2
      INTEGER mat
      DOUBLE PRECISION badc1 , badc2 , eps , large , small
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN
      DOUBLE PRECISION DLAMCH
EXTERNAL    LSAMEN,   DLAMCH
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , SQRT
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
!     Set some parameters we don't plan to change.
!
      Dist = 'S'
      Mode = 3
!
!     xQR, xLQ, xQL, xRQ:  Set parameters to generate a general
!                          M x N matrix.
!
      IF ( LSAMEN(2,c2,'QR') .OR. LSAMEN(2,c2,'LQ') .OR.                &
     &     LSAMEN(2,c2,'QL') .OR. LSAMEN(2,c2,'RQ') ) THEN
!
!        Set TYPE, the type of matrix to be generated.
!
         Type = 'N'
!
!        Set the lower and upper bandwidths.
!
         IF ( Imat==1 ) THEN
            Kl = 0
            Ku = 0
         ELSEIF ( Imat==2 ) THEN
            Kl = 0
            Ku = MAX(N-1,0)
         ELSEIF ( Imat==3 ) THEN
            Kl = MAX(M-1,0)
            Ku = 0
         ELSE
            Kl = MAX(M-1,0)
            Ku = MAX(N-1,0)
         ENDIF
!
!        Set the condition number and norm.
!
         IF ( Imat==5 ) THEN
            Cndnum = badc1
         ELSEIF ( Imat==6 ) THEN
            Cndnum = badc2
         ELSE
            Cndnum = TWO
         ENDIF
!
         IF ( Imat==7 ) THEN
            Anorm = small
         ELSEIF ( Imat==8 ) THEN
            Anorm = large
         ELSE
            Anorm = ONE
         ENDIF
!
      ELSEIF ( LSAMEN(2,c2,'GE') ) THEN
!
!        xGE:  Set parameters to generate a general M x N matrix.
!
!        Set TYPE, the type of matrix to be generated.
!
         Type = 'N'
!
!        Set the lower and upper bandwidths.
!
         IF ( Imat==1 ) THEN
            Kl = 0
            Ku = 0
         ELSEIF ( Imat==2 ) THEN
            Kl = 0
            Ku = MAX(N-1,0)
         ELSEIF ( Imat==3 ) THEN
            Kl = MAX(M-1,0)
            Ku = 0
         ELSE
            Kl = MAX(M-1,0)
            Ku = MAX(N-1,0)
         ENDIF
!
!        Set the condition number and norm.
!
         IF ( Imat==8 ) THEN
            Cndnum = badc1
         ELSEIF ( Imat==9 ) THEN
            Cndnum = badc2
         ELSE
            Cndnum = TWO
         ENDIF
!
         IF ( Imat==10 ) THEN
            Anorm = small
         ELSEIF ( Imat==11 ) THEN
            Anorm = large
         ELSE
            Anorm = ONE
         ENDIF
!
      ELSEIF ( LSAMEN(2,c2,'GB') ) THEN
!
!        xGB:  Set parameters to generate a general banded matrix.
!
!        Set TYPE, the type of matrix to be generated.
!
         Type = 'N'
!
!        Set the condition number and norm.
!
         IF ( Imat==5 ) THEN
            Cndnum = badc1
         ELSEIF ( Imat==6 ) THEN
            Cndnum = TENTH*badc2
         ELSE
            Cndnum = TWO
         ENDIF
!
         IF ( Imat==7 ) THEN
            Anorm = small
         ELSEIF ( Imat==8 ) THEN
            Anorm = large
         ELSE
            Anorm = ONE
         ENDIF
!
      ELSEIF ( LSAMEN(2,c2,'GT') ) THEN
!
!        xGT:  Set parameters to generate a general tridiagonal matrix.
!
!        Set TYPE, the type of matrix to be generated.
!
         Type = 'N'
!
!        Set the lower and upper bandwidths.
!
         IF ( Imat==1 ) THEN
            Kl = 0
         ELSE
            Kl = 1
         ENDIF
         Ku = Kl
!
!        Set the condition number and norm.
!
         IF ( Imat==3 ) THEN
            Cndnum = badc1
         ELSEIF ( Imat==4 ) THEN
            Cndnum = badc2
         ELSE
            Cndnum = TWO
         ENDIF
!
         IF ( Imat==5 .OR. Imat==11 ) THEN
            Anorm = small
         ELSEIF ( Imat==6 .OR. Imat==12 ) THEN
            Anorm = large
         ELSE
            Anorm = ONE
         ENDIF
!
      ELSEIF ( LSAMEN(2,c2,'PO') .OR. LSAMEN(2,c2,'PP') ) THEN
!
!        xPO, xPP: Set parameters to generate a
!        symmetric or Hermitian positive definite matrix.
!
!        Set TYPE, the type of matrix to be generated.
!
         Type = c2(1:1)
!
!        Set the lower and upper bandwidths.
!
         IF ( Imat==1 ) THEN
            Kl = 0
         ELSE
            Kl = MAX(N-1,0)
         ENDIF
         Ku = Kl
!
!        Set the condition number and norm.
!
         IF ( Imat==6 ) THEN
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
      ELSEIF ( LSAMEN(2,c2,'HE') .OR. LSAMEN(2,c2,'HP') .OR.            &
     &         LSAMEN(2,c2,'SY') .OR. LSAMEN(2,c2,'SP') ) THEN
!
!        xHE, xHP, xSY, xSP: Set parameters to generate a
!        symmetric or Hermitian matrix.
!
!        Set TYPE, the type of matrix to be generated.
!
         Type = c2(1:1)
!
!        Set the lower and upper bandwidths.
!
         IF ( Imat==1 ) THEN
            Kl = 0
         ELSE
            Kl = MAX(N-1,0)
         ENDIF
         Ku = Kl
!
!        Set the condition number and norm.
!
         IF ( Imat==7 ) THEN
            Cndnum = badc1
         ELSEIF ( Imat==8 ) THEN
            Cndnum = badc2
         ELSE
            Cndnum = TWO
         ENDIF
!
         IF ( Imat==9 ) THEN
            Anorm = small
         ELSEIF ( Imat==10 ) THEN
            Anorm = large
         ELSE
            Anorm = ONE
         ENDIF
!
      ELSEIF ( LSAMEN(2,c2,'PB') ) THEN
!
!        xPB:  Set parameters to generate a symmetric band matrix.
!
!        Set TYPE, the type of matrix to be generated.
!
         Type = 'P'
!
!        Set the norm and condition number.
!
         IF ( Imat==5 ) THEN
            Cndnum = badc1
         ELSEIF ( Imat==6 ) THEN
            Cndnum = badc2
         ELSE
            Cndnum = TWO
         ENDIF
!
         IF ( Imat==7 ) THEN
            Anorm = small
         ELSEIF ( Imat==8 ) THEN
            Anorm = large
         ELSE
            Anorm = ONE
         ENDIF
!
      ELSEIF ( LSAMEN(2,c2,'PT') ) THEN
!
!        xPT:  Set parameters to generate a symmetric positive definite
!        tridiagonal matrix.
!
         Type = 'P'
         IF ( Imat==1 ) THEN
            Kl = 0
         ELSE
            Kl = 1
         ENDIF
         Ku = Kl
!
!        Set the condition number and norm.
!
         IF ( Imat==3 ) THEN
            Cndnum = badc1
         ELSEIF ( Imat==4 ) THEN
            Cndnum = badc2
         ELSE
            Cndnum = TWO
         ENDIF
!
         IF ( Imat==5 .OR. Imat==11 ) THEN
            Anorm = small
         ELSEIF ( Imat==6 .OR. Imat==12 ) THEN
            Anorm = large
         ELSE
            Anorm = ONE
         ENDIF
!
      ELSEIF ( LSAMEN(2,c2,'TR') .OR. LSAMEN(2,c2,'TP') ) THEN
!
!        xTR, xTP:  Set parameters to generate a triangular matrix
!
!        Set TYPE, the type of matrix to be generated.
!
         Type = 'N'
!
!        Set the lower and upper bandwidths.
!
         mat = ABS(Imat)
         IF ( mat==1 .OR. mat==7 ) THEN
            Kl = 0
            Ku = 0
         ELSEIF ( Imat<0 ) THEN
            Kl = MAX(N-1,0)
            Ku = 0
         ELSE
            Kl = 0
            Ku = MAX(N-1,0)
         ENDIF
!
!        Set the condition number and norm.
!
         IF ( mat==3 .OR. mat==9 ) THEN
            Cndnum = badc1
         ELSEIF ( mat==4 .OR. mat==10 ) THEN
            Cndnum = badc2
         ELSE
            Cndnum = TWO
         ENDIF
!
         IF ( mat==5 ) THEN
            Anorm = small
         ELSEIF ( mat==6 ) THEN
            Anorm = large
         ELSE
            Anorm = ONE
         ENDIF
!
      ELSEIF ( LSAMEN(2,c2,'TB') ) THEN
!
!        xTB:  Set parameters to generate a triangular band matrix.
!
!        Set TYPE, the type of matrix to be generated.
!
         Type = 'N'
!
!        Set the norm and condition number.
!
         IF ( Imat==2 .OR. Imat==8 ) THEN
            Cndnum = badc1
         ELSEIF ( Imat==3 .OR. Imat==9 ) THEN
            Cndnum = badc2
         ELSE
            Cndnum = TWO
         ENDIF
!
         IF ( Imat==4 ) THEN
            Anorm = small
         ELSEIF ( Imat==5 ) THEN
            Anorm = large
         ELSE
            Anorm = ONE
         ENDIF
      ENDIF
      IF ( N<=1 ) Cndnum = ONE
!
!
!     End of zlatb4
!
      END SUBROUTINE zlatb4
