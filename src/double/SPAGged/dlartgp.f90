!*==dlartgp.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DLARTGP generates a plane rotation so that the diagonal is nonnegative.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLARTGP + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlartgp.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlartgp.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlartgp.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLARTGP( F, G, CS, SN, R )
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION   CS, F, G, R, SN
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLARTGP generates a plane rotation so that
!>
!>    [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.
!>    [ -SN  CS  ]     [ G ]     [ 0 ]
!>
!> This is a slower, more accurate version of the Level 1 BLAS routine DROTG,
!> with the following other differences:
!>    F and G are unchanged on return.
!>    If G=0, then CS=(+/-)1 and SN=0.
!>    If F=0 and (G .ne. 0), then CS=0 and SN=(+/-)1.
!>
!> The sign is chosen so that R >= 0.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] F
!> \verbatim
!>          F is DOUBLE PRECISION
!>          The first component of vector to be rotated.
!> \endverbatim
!>
!> \param[in] G
!> \verbatim
!>          G is DOUBLE PRECISION
!>          The second component of vector to be rotated.
!> \endverbatim
!>
!> \param[out] CS
!> \verbatim
!>          CS is DOUBLE PRECISION
!>          The cosine of the rotation.
!> \endverbatim
!>
!> \param[out] SN
!> \verbatim
!>          SN is DOUBLE PRECISION
!>          The sine of the rotation.
!> \endverbatim
!>
!> \param[out] R
!> \verbatim
!>          R is DOUBLE PRECISION
!>          The nonzero component of the rotated vector.
!>
!>  This version has a few statements commented out for thread safety
!>  (machine parameters are computed on each entry). 10 feb 03, SJH.
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
!> \ingroup OTHERauxiliary
!
!  =====================================================================
      SUBROUTINE DLARTGP(F,G,Cs,Sn,R)
      USE F77KINDS                        
      USE S_DLAMCH
      IMPLICIT NONE
!*--DLARTGP101
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0 ,        &
     &                              TWO = 2.0D0
!
! Dummy argument declarations rewritten by SPAG
!
      REAL(R8KIND) , INTENT(IN) :: F
      REAL(R8KIND) , INTENT(IN) :: G
      REAL(R8KIND) , INTENT(INOUT) :: Cs
      REAL(R8KIND) , INTENT(INOUT) :: Sn
      REAL(R8KIND) , INTENT(INOUT) :: R
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: count , i
      REAL(R8KIND) :: eps , f1 , g1 , safmin , safmn2 , safmx2 , scale
!
! End of declarations rewritten by SPAG
!
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     LOGICAL            FIRST
!     ..
!     .. External Functions ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Save statement ..
!     SAVE               FIRST, SAFMX2, SAFMIN, SAFMN2
!     ..
!     .. Data statements ..
!     DATA               FIRST / .TRUE. /
!     ..
!     .. Executable Statements ..
!
!     IF( FIRST ) THEN
      safmin = DLAMCH('S')
      eps = DLAMCH('E')
      safmn2 = DLAMCH('B')**INT(LOG(safmin/eps)/LOG(DLAMCH('B'))/TWO)
      safmx2 = ONE/safmn2
!        FIRST = .FALSE.
!     END IF
      IF ( G==ZERO ) THEN
         Cs = SIGN(ONE,F)
         Sn = ZERO
         R = ABS(F)
      ELSEIF ( F==ZERO ) THEN
         Cs = ZERO
         Sn = SIGN(ONE,G)
         R = ABS(G)
      ELSE
         f1 = F
         g1 = G
         scale = MAX(ABS(f1),ABS(g1))
         IF ( scale>=safmx2 ) THEN
            count = 0
            DO
               count = count + 1
               f1 = f1*safmn2
               g1 = g1*safmn2
               scale = MAX(ABS(f1),ABS(g1))
               IF ( scale<safmx2 .OR. count>=20 ) THEN
                  R = SQRT(f1**2+g1**2)
                  Cs = f1/R
                  Sn = g1/R
                  DO i = 1 , count
                     R = R*safmx2
                  ENDDO
                  EXIT
               ENDIF
            ENDDO
         ELSEIF ( scale<=safmn2 ) THEN
            count = 0
            DO
               count = count + 1
               f1 = f1*safmx2
               g1 = g1*safmx2
               scale = MAX(ABS(f1),ABS(g1))
               IF ( scale>safmn2 ) THEN
                  R = SQRT(f1**2+g1**2)
                  Cs = f1/R
                  Sn = g1/R
                  DO i = 1 , count
                     R = R*safmn2
                  ENDDO
                  EXIT
               ENDIF
            ENDDO
         ELSE
            R = SQRT(f1**2+g1**2)
            Cs = f1/R
            Sn = g1/R
         ENDIF
         IF ( R<ZERO ) THEN
            Cs = -Cs
            Sn = -Sn
            R = -R
         ENDIF
      ENDIF
!
!     End of DLARTGP
!
      END SUBROUTINE DLARTGP
