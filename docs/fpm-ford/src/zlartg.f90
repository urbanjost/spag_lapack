!*==zlartg.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
 
!> \brief \b ZLARTG generates a plane rotation with real cosine and complex sine.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLARTG + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlartg.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlartg.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlartg.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLARTG( F, G, CS, SN, R )
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION   CS
!       COMPLEX*16         F, G, R, SN
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLARTG generates a plane rotation so that
!>
!>    [  CS  SN  ]     [ F ]     [ R ]
!>    [  __      ]  .  [   ]  =  [   ]   where CS**2 + |SN|**2 = 1.
!>    [ -SN  CS  ]     [ G ]     [ 0 ]
!>
!> This is a faster version of the BLAS1 routine ZROTG, except for
!> the following differences:
!>    F and G are unchanged on return.
!>    If G=0, then CS=1 and SN=0.
!>    If F=0, then CS=0 and SN is chosen so that R is real.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] F
!> \verbatim
!>          F is COMPLEX*16
!>          The first component of vector to be rotated.
!> \endverbatim
!>
!> \param[in] G
!> \verbatim
!>          G is COMPLEX*16
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
!>          SN is COMPLEX*16
!>          The sine of the rotation.
!> \endverbatim
!>
!> \param[out] R
!> \verbatim
!>          R is COMPLEX*16
!>          The nonzero component of the rotated vector.
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
!> \ingroup complex16OTHERauxiliary
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  3-5-96 - Modified with a new algorithm by W. Kahan and J. Demmel
!>
!>  This version has a few statements commented out for thread safety
!>  (machine parameters are computed on each entry). 10 feb 03, SJH.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE ZLARTG(F,G,Cs,Sn,R)
      IMPLICIT NONE
!*--ZLARTG108
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION Cs
      COMPLEX*16 F , G , R , Sn
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION TWO , ONE , ZERO
      PARAMETER (TWO=2.0D+0,ONE=1.0D+0,ZERO=0.0D+0)
      COMPLEX*16 CZERO
      PARAMETER (CZERO=(0.0D+0,0.0D+0))
!     ..
!     .. Local Scalars ..
!     LOGICAL            FIRST
      INTEGER count , i
      DOUBLE PRECISION d , di , dr , eps , f2 , f2s , g2 , g2s ,        &
     &                 safmin , safmn2 , safmx2 , scale
      COMPLEX*16 ff , fs , gs
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH , DLAPY2
      LOGICAL DISNAN
      EXTERNAL DLAMCH , DLAPY2 , DISNAN
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , DBLE , DCMPLX , DCONJG , DIMAG , INT , LOG , MAX ,&
     &          SQRT
!     ..
!     .. Statement Functions ..
      DOUBLE PRECISION ABS1 , ABSSQ
!     ..
!     .. Statement Function definitions ..
      ABS1(ff) = MAX(ABS(DBLE(ff)),ABS(DIMAG(ff)))
      ABSSQ(ff) = DBLE(ff)**2 + DIMAG(ff)**2
!     ..
!     .. Executable Statements ..
!
      safmin = DLAMCH('S')
      eps = DLAMCH('E')
      safmn2 = DLAMCH('B')**INT(LOG(safmin/eps)/LOG(DLAMCH('B'))/TWO)
      safmx2 = ONE/safmn2
      scale = MAX(ABS1(F),ABS1(G))
      fs = F
      gs = G
      count = 0
      IF ( scale>=safmx2 ) THEN
         DO
            count = count + 1
            fs = fs*safmn2
            gs = gs*safmn2
            scale = scale*safmn2
            IF ( scale<safmx2 .OR. count>=20 ) EXIT
         ENDDO
      ELSEIF ( scale<=safmn2 ) THEN
         IF ( G==CZERO .OR. DISNAN(ABS(G)) ) THEN
            Cs = ONE
            Sn = CZERO
            R = F
            RETURN
         ENDIF
         DO
            count = count - 1
            fs = fs*safmx2
            gs = gs*safmx2
            scale = scale*safmx2
            IF ( scale>safmn2 ) EXIT
         ENDDO
      ENDIF
      f2 = ABSSQ(fs)
      g2 = ABSSQ(gs)
      IF ( f2<=MAX(g2,ONE)*safmin ) THEN
!
!        This is a rare case: F is very small.
!
         IF ( F==CZERO ) THEN
            Cs = ZERO
            R = DLAPY2(DBLE(G),DIMAG(G))
!           Do complex/real division explicitly with two real divisions
            d = DLAPY2(DBLE(gs),DIMAG(gs))
            Sn = DCMPLX(DBLE(gs)/d,-DIMAG(gs)/d)
            RETURN
         ENDIF
         f2s = DLAPY2(DBLE(fs),DIMAG(fs))
!        G2 and G2S are accurate
!        G2 is at least SAFMIN, and G2S is at least SAFMN2
         g2s = SQRT(g2)
!        Error in CS from underflow in F2S is at most
!        UNFL / SAFMN2 .lt. sqrt(UNFL*EPS) .lt. EPS
!        If MAX(G2,ONE)=G2, then F2 .lt. G2*SAFMIN,
!        and so CS .lt. sqrt(SAFMIN)
!        If MAX(G2,ONE)=ONE, then F2 .lt. SAFMIN
!        and so CS .lt. sqrt(SAFMIN)/SAFMN2 = sqrt(EPS)
!        Therefore, CS = F2S/G2S / sqrt( 1 + (F2S/G2S)**2 ) = F2S/G2S
         Cs = f2s/g2s
!        Make sure abs(FF) = 1
!        Do complex/real division explicitly with 2 real divisions
         IF ( ABS1(F)>ONE ) THEN
            d = DLAPY2(DBLE(F),DIMAG(F))
            ff = DCMPLX(DBLE(F)/d,DIMAG(F)/d)
         ELSE
            dr = safmx2*DBLE(F)
            di = safmx2*DIMAG(F)
            d = DLAPY2(dr,di)
            ff = DCMPLX(dr/d,di/d)
         ENDIF
         Sn = ff*DCMPLX(DBLE(gs)/g2s,-DIMAG(gs)/g2s)
         R = Cs*F + Sn*G
      ELSE
!
!        This is the most common case.
!        Neither F2 nor F2/G2 are less than SAFMIN
!        F2S cannot overflow, and it is accurate
!
         f2s = SQRT(ONE+g2/f2)
!        Do the F2S(real)*FS(complex) multiply with two real multiplies
         R = DCMPLX(f2s*DBLE(fs),f2s*DIMAG(fs))
         Cs = ONE/f2s
         d = f2 + g2
!        Do complex/real division explicitly with two real divisions
         Sn = DCMPLX(DBLE(R)/d,DIMAG(R)/d)
         Sn = Sn*DCONJG(gs)
         IF ( count/=0 ) THEN
            IF ( count>0 ) THEN
               DO i = 1 , count
                  R = R*safmx2
               ENDDO
            ELSE
               DO i = 1 , -count
                  R = R*safmn2
               ENDDO
            ENDIF
         ENDIF
      ENDIF
!
!     End of ZLARTG
!
      END SUBROUTINE ZLARTG
