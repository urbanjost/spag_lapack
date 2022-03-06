!*==dlatb9.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b DLATB9
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLATB9( PATH, IMAT, M, P, N, TYPE, KLA, KUA, KLB, KUB,
!                          ANORM, BNORM, MODEA, MODEB, CNDNMA, CNDNMB,
!                          DISTA, DISTB )
!
!       .. Scalar Arguments ..
!       CHARACTER          DISTA, DISTB, TYPE
!       CHARACTER*3        PATH
!       INTEGER            IMAT, KLA, KLB, KUA, KUB, M, MODEA, MODEB, N, P
!       DOUBLE PRECISION   ANORM, BNORM, CNDNMA, CNDNMB
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLATB9 sets parameters for the matrix generator based on the type of
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
!>          = 1:   A: diagonal, B: upper triangular
!>          = 2:   A: upper triangular, B: upper triangular
!>          = 3:   A: lower triangular, B: upper triangular
!>          Else:  A: general dense, B: general dense
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows in the matrix to be generated.
!> \endverbatim
!>
!> \param[in] P
!> \verbatim
!>          P is INTEGER
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
!>          = 'S':  symmetric matrix;
!>          = 'P':  symmetric positive (semi)definite matrix;
!>          = 'N':  nonsymmetric matrix.
!> \endverbatim
!>
!> \param[out] KLA
!> \verbatim
!>          KLA is INTEGER
!>          The lower band width of the matrix to be generated.
!> \endverbatim
!>
!> \param[out] KUA
!> \verbatim
!>          KUA is INTEGER
!>          The upper band width of the matrix to be generated.
!> \endverbatim
!>
!> \param[out] KLB
!> \verbatim
!>          KLB is INTEGER
!>          The lower band width of the matrix to be generated.
!> \endverbatim
!>
!> \param[out] KUB
!> \verbatim
!>          KUA is INTEGER
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
!> \param[out] BNORM
!> \verbatim
!>          BNORM is DOUBLE PRECISION
!>          The desired norm of the matrix to be generated.  The diagonal
!>          matrix of singular values or eigenvalues is scaled by this
!>          value.
!> \endverbatim
!>
!> \param[out] MODEA
!> \verbatim
!>          MODEA is INTEGER
!>          A key indicating how to choose the vector of eigenvalues.
!> \endverbatim
!>
!> \param[out] MODEB
!> \verbatim
!>          MODEB is INTEGER
!>          A key indicating how to choose the vector of eigenvalues.
!> \endverbatim
!>
!> \param[out] CNDNMA
!> \verbatim
!>          CNDNMA is DOUBLE PRECISION
!>          The desired condition number.
!> \endverbatim
!>
!> \param[out] CNDNMB
!> \verbatim
!>          CNDNMB is DOUBLE PRECISION
!>          The desired condition number.
!> \endverbatim
!>
!> \param[out] DISTA
!> \verbatim
!>          DISTA is CHARACTER*1
!>          The type of distribution to be used by the random number
!>          generator.
!> \endverbatim
!>
!> \param[out] DISTB
!> \verbatim
!>          DISTB is CHARACTER*1
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
!> \ingroup double_eig
!
!  =====================================================================
      SUBROUTINE DLATB9(Path,Imat,M,P,N,Type,Kla,Kua,Klb,Kub,Anorm,     &
     &                  Bnorm,Modea,Modeb,Cndnma,Cndnmb,Dista,Distb)
      IMPLICIT NONE
!*--DLATB9173
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Dista , Distb , Type
      CHARACTER*3 Path
      INTEGER Imat , Kla , Klb , Kua , Kub , M , Modea , Modeb , N , P
      DOUBLE PRECISION Anorm , Bnorm , Cndnma , Cndnmb
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION SHRINK , TENTH
      PARAMETER (SHRINK=0.25D0,TENTH=0.1D+0)
      DOUBLE PRECISION ONE , TEN
      PARAMETER (ONE=1.0D+0,TEN=1.0D+1)
!     ..
!     .. Local Scalars ..
      LOGICAL first
      DOUBLE PRECISION badc1 , badc2 , eps , large , small
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN
      DOUBLE PRECISION DLAMCH
      EXTERNAL LSAMEN , DLAMCH
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , SQRT
!     ..
!     .. External Subroutines ..
      EXTERNAL DLABAD
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
!     Set some parameters we don't plan to change.
!
      Type = 'N'
      Dista = 'S'
      Distb = 'S'
      Modea = 3
      Modeb = 4
!
!     Set the lower and upper bandwidths.
!
      IF ( LSAMEN(3,Path,'GRQ') .OR. LSAMEN(3,Path,'LSE') .OR.          &
     &     LSAMEN(3,Path,'GSV') ) THEN
!
!        A: M by N, B: P by N
!
         IF ( Imat==1 ) THEN
!
!           A: diagonal, B: upper triangular
!
            Kla = 0
            Kua = 0
            Klb = 0
            Kub = MAX(N-1,0)
!
         ELSEIF ( Imat==2 ) THEN
!
!           A: upper triangular, B: upper triangular
!
            Kla = 0
            Kua = MAX(N-1,0)
            Klb = 0
            Kub = MAX(N-1,0)
!
         ELSEIF ( Imat==3 ) THEN
!
!           A: lower triangular, B: upper triangular
!
            Kla = MAX(M-1,0)
            Kua = 0
            Klb = 0
            Kub = MAX(N-1,0)
!
         ELSE
!
!           A: general dense, B: general dense
!
            Kla = MAX(M-1,0)
            Kua = MAX(N-1,0)
            Klb = MAX(P-1,0)
            Kub = MAX(N-1,0)
!
         ENDIF
!
      ELSEIF ( LSAMEN(3,Path,'GQR') .OR. LSAMEN(3,Path,'GLM') ) THEN
!
!        A: N by M, B: N by P
!
         IF ( Imat==1 ) THEN
!
!           A: diagonal, B: lower triangular
!
            Kla = 0
            Kua = 0
            Klb = MAX(N-1,0)
            Kub = 0
         ELSEIF ( Imat==2 ) THEN
!
!           A: lower triangular, B: diagonal
!
            Kla = MAX(N-1,0)
            Kua = 0
            Klb = 0
            Kub = 0
!
         ELSEIF ( Imat==3 ) THEN
!
!           A: lower triangular, B: upper triangular
!
            Kla = MAX(N-1,0)
            Kua = 0
            Klb = 0
            Kub = MAX(P-1,0)
!
         ELSE
!
!           A: general dense, B: general dense
!
            Kla = MAX(N-1,0)
            Kua = MAX(M-1,0)
            Klb = MAX(N-1,0)
            Kub = MAX(P-1,0)
         ENDIF
!
      ENDIF
!
!     Set the condition number and norm.
!
      Cndnma = TEN*TEN
      Cndnmb = TEN
      IF ( LSAMEN(3,Path,'GQR') .OR. LSAMEN(3,Path,'GRQ') .OR.          &
     &     LSAMEN(3,Path,'GSV') ) THEN
         IF ( Imat==5 ) THEN
            Cndnma = badc1
            Cndnmb = badc1
         ELSEIF ( Imat==6 ) THEN
            Cndnma = badc2
            Cndnmb = badc2
         ELSEIF ( Imat==7 ) THEN
            Cndnma = badc1
            Cndnmb = badc2
         ELSEIF ( Imat==8 ) THEN
            Cndnma = badc2
            Cndnmb = badc1
         ENDIF
      ENDIF
!
      Anorm = TEN
      Bnorm = TEN*TEN*TEN
      IF ( LSAMEN(3,Path,'GQR') .OR. LSAMEN(3,Path,'GRQ') ) THEN
         IF ( Imat==7 ) THEN
            Anorm = small
            Bnorm = large
         ELSEIF ( Imat==8 ) THEN
            Anorm = large
            Bnorm = small
         ENDIF
      ENDIF
!
      IF ( N<=1 ) THEN
         Cndnma = ONE
         Cndnmb = ONE
      ENDIF
!
!
!     End of DLATB9
!
      END SUBROUTINE DLATB9
