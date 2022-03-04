!*==slartg.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SLARTG generates a plane rotation with real cosine and real sine.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLARTG + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slartg.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slartg.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slartg.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLARTG( F, G, CS, SN, R )
!
!       .. Scalar Arguments ..
!       REAL               CS, F, G, R, SN
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLARTG generate a plane rotation so that
!>
!>    [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.
!>    [ -SN  CS  ]     [ G ]     [ 0 ]
!>
!> This is a slower, more accurate version of the BLAS1 routine SROTG,
!> with the following other differences:
!>    F and G are unchanged on return.
!>    If G=0, then CS=1 and SN=0.
!>    If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any
!>       floating point operations (saves work in SBDSQR when
!>       there are zeros on the diagonal).
!>
!> If F exceeds G in magnitude, CS will be positive.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] F
!> \verbatim
!>          F is REAL
!>          The first component of vector to be rotated.
!> \endverbatim
!>
!> \param[in] G
!> \verbatim
!>          G is REAL
!>          The second component of vector to be rotated.
!> \endverbatim
!>
!> \param[out] CS
!> \verbatim
!>          CS is REAL
!>          The cosine of the rotation.
!> \endverbatim
!>
!> \param[out] SN
!> \verbatim
!>          SN is REAL
!>          The sine of the rotation.
!> \endverbatim
!>
!> \param[out] R
!> \verbatim
!>          R is REAL
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
      SUBROUTINE SLARTG(F,G,Cs,Sn,R)
      USE S_SLAMCH
      IMPLICIT NONE
!*--SLARTG102
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0 , TWO = 2.0E0
!
! Dummy argument declarations rewritten by SPAG
!
      REAL , INTENT(IN) :: F
      REAL , INTENT(IN) :: G
      REAL , INTENT(INOUT) :: Cs
      REAL , INTENT(INOUT) :: Sn
      REAL , INTENT(INOUT) :: R
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: count , i
      REAL :: eps , f1 , g1 , safmin , safmn2 , safmx2 , scale
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
      safmin = SLAMCH('S')
      eps = SLAMCH('E')
      safmn2 = SLAMCH('B')**INT(LOG(safmin/eps)/LOG(SLAMCH('B'))/TWO)
      safmx2 = ONE/safmn2
!        FIRST = .FALSE.
!     END IF
      IF ( G==ZERO ) THEN
         Cs = ONE
         Sn = ZERO
         R = F
      ELSEIF ( F==ZERO ) THEN
         Cs = ZERO
         Sn = ONE
         R = G
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
         IF ( ABS(F)>ABS(G) .AND. Cs<ZERO ) THEN
            Cs = -Cs
            Sn = -Sn
            R = -R
         ENDIF
      ENDIF
!
!     End of SLARTG
!
      END SUBROUTINE SLARTG
