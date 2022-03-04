!*==srotmg.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
!> \brief \b SROTMG
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SROTMG(SD1,SD2,SX1,SY1,SPARAM)
!
!       .. Scalar Arguments ..
!       REAL SD1,SD2,SX1,SY1
!       ..
!       .. Array Arguments ..
!       REAL SPARAM(5)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    CONSTRUCT THE MODIFIED GIVENS TRANSFORMATION MATRIX H WHICH ZEROS
!>    THE SECOND COMPONENT OF THE 2-VECTOR  (SQRT(SD1)*SX1,SQRT(SD2)*>    SY2)**T.
!>    WITH SPARAM(1)=SFLAG, H HAS ONE OF THE FOLLOWING FORMS..
!>
!>    SFLAG=-1.E0     SFLAG=0.E0        SFLAG=1.E0     SFLAG=-2.E0
!>
!>      (SH11  SH12)    (1.E0  SH12)    (SH11  1.E0)    (1.E0  0.E0)
!>    H=(          )    (          )    (          )    (          )
!>      (SH21  SH22),   (SH21  1.E0),   (-1.E0 SH22),   (0.E0  1.E0).
!>    LOCATIONS 2-4 OF SPARAM CONTAIN SH11,SH21,SH12, AND SH22
!>    RESPECTIVELY. (VALUES OF 1.E0, -1.E0, OR 0.E0 IMPLIED BY THE
!>    VALUE OF SPARAM(1) ARE NOT STORED IN SPARAM.)
!>
!>    THE VALUES OF GAMSQ AND RGAMSQ SET IN THE DATA STATEMENT MAY BE
!>    INEXACT.  THIS IS OK AS THEY ARE ONLY USED FOR TESTING THE SIZE
!>    OF SD1 AND SD2.  ALL ACTUAL SCALING OF DATA IS DONE USING GAM.
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in,out] SD1
!> \verbatim
!>          SD1 is REAL
!> \endverbatim
!>
!> \param[in,out] SD2
!> \verbatim
!>          SD2 is REAL
!> \endverbatim
!>
!> \param[in,out] SX1
!> \verbatim
!>          SX1 is REAL
!> \endverbatim
!>
!> \param[in] SY1
!> \verbatim
!>          SY1 is REAL
!> \endverbatim
!>
!> \param[out] SPARAM
!> \verbatim
!>          SPARAM is REAL array, dimension (5)
!>     SPARAM(1)=SFLAG
!>     SPARAM(2)=SH11
!>     SPARAM(3)=SH21
!>     SPARAM(4)=SH12
!>     SPARAM(5)=SH22
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
!> \date November 2017
!
!> \ingroup single_blas_level1
!
!  =====================================================================
      SUBROUTINE SROTMG(Sd1,Sd2,Sx1,Sy1,Sparam)
      IMPLICIT NONE
!*--SROTMG94
!
! Dummy argument declarations rewritten by SPAG
!
      REAL , INTENT(INOUT) :: Sd1
      REAL , INTENT(INOUT) :: Sd2
      REAL , INTENT(INOUT) :: Sx1
      REAL , INTENT(IN) :: Sy1
      REAL , INTENT(OUT) , DIMENSION(5) :: Sparam
!
! Local variable declarations rewritten by SPAG
!
      REAL , SAVE :: gam , gamsq , one , rgamsq , two , zero
      REAL :: sflag , sh11 , sh12 , sh21 , sh22 , sp1 , sp2 , sq1 ,     &
     &        sq2 , stemp , su
!
! End of declarations rewritten by SPAG
!
!
! Local variable declarations rewritten by SPAG
!
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Data statements ..
!
      DATA zero , one , two/0.E0 , 1.E0 , 2.E0/
      DATA gam , gamsq , rgamsq/4096.E0 , 1.67772E7 , 5.96046E-8/
!     ..
 
      IF ( Sd1<zero ) THEN
!        GO ZERO-H-D-AND-SX1..
         sflag = -one
         sh11 = zero
         sh12 = zero
         sh21 = zero
         sh22 = zero
!
         Sd1 = zero
         Sd2 = zero
         Sx1 = zero
      ELSE
!        CASE-SD1-NONNEGATIVE
         sp2 = Sd2*Sy1
         IF ( sp2==zero ) THEN
            sflag = -two
            Sparam(1) = sflag
            RETURN
         ENDIF
!        REGULAR-CASE..
         sp1 = Sd1*Sx1
         sq2 = sp2*Sy1
         sq1 = sp1*Sx1
!
         IF ( ABS(sq1)>ABS(sq2) ) THEN
            sh21 = -Sy1/Sx1
            sh12 = sp2/sp1
!
            su = one - sh12*sh21
!
            IF ( su>zero ) THEN
               sflag = zero
               Sd1 = Sd1/su
               Sd2 = Sd2/su
               Sx1 = Sx1*su
            ELSE
!            This code path if here for safety. We do not expect this
!            condition to ever hold except in edge cases with rounding
!            errors. See DOI: 10.1145/355841.355847
               sflag = -one
               sh11 = zero
               sh12 = zero
               sh21 = zero
               sh22 = zero
!
               Sd1 = zero
               Sd2 = zero
               Sx1 = zero
            ENDIF
 
         ELSEIF ( sq2<zero ) THEN
!              GO ZERO-H-D-AND-SX1..
            sflag = -one
            sh11 = zero
            sh12 = zero
            sh21 = zero
            sh22 = zero
!
            Sd1 = zero
            Sd2 = zero
            Sx1 = zero
         ELSE
            sflag = one
            sh11 = sp1/sp2
            sh22 = Sx1/Sy1
            su = one + sh11*sh22
            stemp = Sd2/su
            Sd2 = Sd1/su
            Sd1 = stemp
            Sx1 = Sy1*su
         ENDIF
 
!     PROCEDURE..SCALE-CHECK
         IF ( Sd1/=zero ) THEN
            DO WHILE ( (Sd1<=rgamsq) .OR. (Sd1>=gamsq) )
               IF ( sflag==zero ) THEN
                  sh11 = one
                  sh22 = one
                  sflag = -one
               ELSE
                  sh21 = -one
                  sh12 = one
                  sflag = -one
               ENDIF
               IF ( Sd1<=rgamsq ) THEN
                  Sd1 = Sd1*gam**2
                  Sx1 = Sx1/gam
                  sh11 = sh11/gam
                  sh12 = sh12/gam
               ELSE
                  Sd1 = Sd1/gam**2
                  Sx1 = Sx1*gam
                  sh11 = sh11*gam
                  sh12 = sh12*gam
               ENDIF
            ENDDO
         ENDIF
 
         IF ( Sd2/=zero ) THEN
            DO WHILE ( (ABS(Sd2)<=rgamsq) .OR. (ABS(Sd2)>=gamsq) )
               IF ( sflag==zero ) THEN
                  sh11 = one
                  sh22 = one
                  sflag = -one
               ELSE
                  sh21 = -one
                  sh12 = one
                  sflag = -one
               ENDIF
               IF ( ABS(Sd2)<=rgamsq ) THEN
                  Sd2 = Sd2*gam**2
                  sh21 = sh21/gam
                  sh22 = sh22/gam
               ELSE
                  Sd2 = Sd2/gam**2
                  sh21 = sh21*gam
                  sh22 = sh22*gam
               ENDIF
            ENDDO
         ENDIF
 
      ENDIF
 
      IF ( sflag<zero ) THEN
         Sparam(2) = sh11
         Sparam(3) = sh21
         Sparam(4) = sh12
         Sparam(5) = sh22
      ELSEIF ( sflag==zero ) THEN
         Sparam(3) = sh21
         Sparam(4) = sh12
      ELSE
         Sparam(2) = sh11
         Sparam(5) = sh22
      ENDIF
 
      Sparam(1) = sflag
      END SUBROUTINE SROTMG
