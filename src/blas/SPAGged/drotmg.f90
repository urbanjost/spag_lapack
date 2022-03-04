!*==drotmg.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DROTMG
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DROTMG(DD1,DD2,DX1,DY1,DPARAM)
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION DD1,DD2,DX1,DY1
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION DPARAM(5)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    CONSTRUCT THE MODIFIED GIVENS TRANSFORMATION MATRIX H WHICH ZEROS
!>    THE SECOND COMPONENT OF THE 2-VECTOR  (DSQRT(DD1)*DX1,DSQRT(DD2)*>    DY2)**T.
!>    WITH DPARAM(1)=DFLAG, H HAS ONE OF THE FOLLOWING FORMS..
!>
!>    DFLAG=-1.D0     DFLAG=0.D0        DFLAG=1.D0     DFLAG=-2.D0
!>
!>      (DH11  DH12)    (1.D0  DH12)    (DH11  1.D0)    (1.D0  0.D0)
!>    H=(          )    (          )    (          )    (          )
!>      (DH21  DH22),   (DH21  1.D0),   (-1.D0 DH22),   (0.D0  1.D0).
!>    LOCATIONS 2-4 OF DPARAM CONTAIN DH11, DH21, DH12, AND DH22
!>    RESPECTIVELY. (VALUES OF 1.D0, -1.D0, OR 0.D0 IMPLIED BY THE
!>    VALUE OF DPARAM(1) ARE NOT STORED IN DPARAM.)
!>
!>    THE VALUES OF GAMSQ AND RGAMSQ SET IN THE DATA STATEMENT MAY BE
!>    INEXACT.  THIS IS OK AS THEY ARE ONLY USED FOR TESTING THE SIZE
!>    OF DD1 AND DD2.  ALL ACTUAL SCALING OF DATA IS DONE USING GAM.
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in,out] DD1
!> \verbatim
!>          DD1 is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in,out] DD2
!> \verbatim
!>          DD2 is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in,out] DX1
!> \verbatim
!>          DX1 is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in] DY1
!> \verbatim
!>          DY1 is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[out] DPARAM
!> \verbatim
!>          DPARAM is DOUBLE PRECISION array, dimension (5)
!>     DPARAM(1)=DFLAG
!>     DPARAM(2)=DH11
!>     DPARAM(3)=DH21
!>     DPARAM(4)=DH12
!>     DPARAM(5)=DH22
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
!> \ingroup double_blas_level1
!
!  =====================================================================
      SUBROUTINE DROTMG(Dd1,Dd2,Dx1,Dy1,Dparam)
      USE F77KINDS                        
      IMPLICIT NONE
!*--DROTMG95
!
! Dummy argument declarations rewritten by SPAG
!
      REAL(R8KIND) , INTENT(INOUT) :: Dd1
      REAL(R8KIND) , INTENT(INOUT) :: Dd2
      REAL(R8KIND) , INTENT(INOUT) :: Dx1
      REAL(R8KIND) , INTENT(IN) :: Dy1
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(5) :: Dparam
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: dflag , dh11 , dh12 , dh21 , dh22 , dp1 , dp2 ,   &
     &                dq1 , dq2 , dtemp , du
      REAL(R8KIND) , SAVE :: gam , gamsq , one , rgamsq , two , zero
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
      DATA zero , one , two/0.D0 , 1.D0 , 2.D0/
      DATA gam , gamsq , rgamsq/4096.D0 , 16777216.D0 , 5.9604645D-8/
!     ..
 
      IF ( Dd1<zero ) THEN
!        GO ZERO-H-D-AND-DX1..
         dflag = -one
         dh11 = zero
         dh12 = zero
         dh21 = zero
         dh22 = zero
!
         Dd1 = zero
         Dd2 = zero
         Dx1 = zero
      ELSE
!        CASE-DD1-NONNEGATIVE
         dp2 = Dd2*Dy1
         IF ( dp2==zero ) THEN
            dflag = -two
            Dparam(1) = dflag
            RETURN
         ENDIF
!        REGULAR-CASE..
         dp1 = Dd1*Dx1
         dq2 = dp2*Dy1
         dq1 = dp1*Dx1
!
         IF ( DABS(dq1)>DABS(dq2) ) THEN
            dh21 = -Dy1/Dx1
            dh12 = dp2/dp1
!
            du = one - dh12*dh21
!
            IF ( du>zero ) THEN
               dflag = zero
               Dd1 = Dd1/du
               Dd2 = Dd2/du
               Dx1 = Dx1*du
            ELSE
!            This code path if here for safety. We do not expect this
!            condition to ever hold except in edge cases with rounding
!            errors. See DOI: 10.1145/355841.355847
               dflag = -one
               dh11 = zero
               dh12 = zero
               dh21 = zero
               dh22 = zero
!
               Dd1 = zero
               Dd2 = zero
               Dx1 = zero
            ENDIF
 
         ELSEIF ( dq2<zero ) THEN
!              GO ZERO-H-D-AND-DX1..
            dflag = -one
            dh11 = zero
            dh12 = zero
            dh21 = zero
            dh22 = zero
!
            Dd1 = zero
            Dd2 = zero
            Dx1 = zero
         ELSE
            dflag = one
            dh11 = dp1/dp2
            dh22 = Dx1/Dy1
            du = one + dh11*dh22
            dtemp = Dd2/du
            Dd2 = Dd1/du
            Dd1 = dtemp
            Dx1 = Dy1*du
         ENDIF
 
!     PROCEDURE..SCALE-CHECK
         IF ( Dd1/=zero ) THEN
            DO WHILE ( (Dd1<=rgamsq) .OR. (Dd1>=gamsq) )
               IF ( dflag==zero ) THEN
                  dh11 = one
                  dh22 = one
                  dflag = -one
               ELSE
                  dh21 = -one
                  dh12 = one
                  dflag = -one
               ENDIF
               IF ( Dd1<=rgamsq ) THEN
                  Dd1 = Dd1*gam**2
                  Dx1 = Dx1/gam
                  dh11 = dh11/gam
                  dh12 = dh12/gam
               ELSE
                  Dd1 = Dd1/gam**2
                  Dx1 = Dx1*gam
                  dh11 = dh11*gam
                  dh12 = dh12*gam
               ENDIF
            ENDDO
         ENDIF
 
         IF ( Dd2/=zero ) THEN
            DO WHILE ( (DABS(Dd2)<=rgamsq) .OR. (DABS(Dd2)>=gamsq) )
               IF ( dflag==zero ) THEN
                  dh11 = one
                  dh22 = one
                  dflag = -one
               ELSE
                  dh21 = -one
                  dh12 = one
                  dflag = -one
               ENDIF
               IF ( DABS(Dd2)<=rgamsq ) THEN
                  Dd2 = Dd2*gam**2
                  dh21 = dh21/gam
                  dh22 = dh22/gam
               ELSE
                  Dd2 = Dd2/gam**2
                  dh21 = dh21*gam
                  dh22 = dh22*gam
               ENDIF
            ENDDO
         ENDIF
 
      ENDIF
 
      IF ( dflag<zero ) THEN
         Dparam(2) = dh11
         Dparam(3) = dh21
         Dparam(4) = dh12
         Dparam(5) = dh22
      ELSEIF ( dflag==zero ) THEN
         Dparam(3) = dh21
         Dparam(4) = dh12
      ELSE
         Dparam(2) = dh11
         Dparam(5) = dh22
      ENDIF
 
      Dparam(1) = dflag
      END SUBROUTINE DROTMG
