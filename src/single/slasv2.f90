!*==slasv2.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SLASV2 computes the singular value decomposition of a 2-by-2 triangular matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLASV2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasv2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasv2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasv2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLASV2( F, G, H, SSMIN, SSMAX, SNR, CSR, SNL, CSL )
!
!       .. Scalar Arguments ..
!       REAL               CSL, CSR, F, G, H, SNL, SNR, SSMAX, SSMIN
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLASV2 computes the singular value decomposition of a 2-by-2
!> triangular matrix
!>    [  F   G  ]
!>    [  0   H  ].
!> On return, abs(SSMAX) is the larger singular value, abs(SSMIN) is the
!> smaller singular value, and (CSL,SNL) and (CSR,SNR) are the left and
!> right singular vectors for abs(SSMAX), giving the decomposition
!>
!>    [ CSL  SNL ] [  F   G  ] [ CSR -SNR ]  =  [ SSMAX   0   ]
!>    [-SNL  CSL ] [  0   H  ] [ SNR  CSR ]     [  0    SSMIN ].
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] F
!> \verbatim
!>          F is REAL
!>          The (1,1) element of the 2-by-2 matrix.
!> \endverbatim
!>
!> \param[in] G
!> \verbatim
!>          G is REAL
!>          The (1,2) element of the 2-by-2 matrix.
!> \endverbatim
!>
!> \param[in] H
!> \verbatim
!>          H is REAL
!>          The (2,2) element of the 2-by-2 matrix.
!> \endverbatim
!>
!> \param[out] SSMIN
!> \verbatim
!>          SSMIN is REAL
!>          abs(SSMIN) is the smaller singular value.
!> \endverbatim
!>
!> \param[out] SSMAX
!> \verbatim
!>          SSMAX is REAL
!>          abs(SSMAX) is the larger singular value.
!> \endverbatim
!>
!> \param[out] SNL
!> \verbatim
!>          SNL is REAL
!> \endverbatim
!>
!> \param[out] CSL
!> \verbatim
!>          CSL is REAL
!>          The vector (CSL, SNL) is a unit left singular vector for the
!>          singular value abs(SSMAX).
!> \endverbatim
!>
!> \param[out] SNR
!> \verbatim
!>          SNR is REAL
!> \endverbatim
!>
!> \param[out] CSR
!> \verbatim
!>          CSR is REAL
!>          The vector (CSR, SNR) is a unit right singular vector for the
!>          singular value abs(SSMAX).
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
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Any input parameter may be aliased with any output parameter.
!>
!>  Barring over/underflow and assuming a guard digit in subtraction, all
!>  output quantities are correct to within a few units in the last
!>  place (ulps).
!>
!>  In IEEE arithmetic, the code works correctly if one matrix element is
!>  infinite.
!>
!>  Overflow will not occur unless the largest singular value itself
!>  overflows or is within a few ulps of overflow. (On machines with
!>  partial overflow, like the Cray, overflow may occur if the largest
!>  singular value is within a factor of 2 of overflow.)
!>
!>  Underflow is harmless if underflow is gradual. Otherwise, results
!>  may correspond to a matrix modified by perturbations of size near
!>  the underflow threshold.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE SLASV2(F,G,H,Ssmin,Ssmax,Snr,Csr,Snl,Csl)
      IMPLICIT NONE
!*--SLASV2142
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      REAL Csl , Csr , F , G , H , Snl , Snr , Ssmax , Ssmin
!     ..
!
! =====================================================================
!
!     .. Parameters ..
      REAL ZERO
      PARAMETER (ZERO=0.0E0)
      REAL HALF
      PARAMETER (HALF=0.5E0)
      REAL ONE
      PARAMETER (ONE=1.0E0)
      REAL TWO
      PARAMETER (TWO=2.0E0)
      REAL FOUR
      PARAMETER (FOUR=4.0E0)
!     ..
!     .. Local Scalars ..
      LOGICAL gasmal , swap
      INTEGER pmax
      REAL a , clt , crt , d , fa , ft , ga , gt , ha , ht , l , m ,    &
     &     mm , r , s , slt , srt , t , temp , tsign , tt
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , SIGN , SQRT
!     ..
!     .. External Functions ..
      REAL SLAMCH
      EXTERNAL SLAMCH
!     ..
!     .. Executable Statements ..
!
      ft = F
      fa = ABS(ft)
      ht = H
      ha = ABS(H)
!
!     PMAX points to the maximum absolute element of matrix
!       PMAX = 1 if F largest in absolute values
!       PMAX = 2 if G largest in absolute values
!       PMAX = 3 if H largest in absolute values
!
      pmax = 1
      swap = (ha>fa)
      IF ( swap ) THEN
         pmax = 3
         temp = ft
         ft = ht
         ht = temp
         temp = fa
         fa = ha
         ha = temp
!
!        Now FA .ge. HA
!
      ENDIF
      gt = G
      ga = ABS(gt)
      IF ( ga==ZERO ) THEN
!
!        Diagonal matrix
!
         Ssmin = ha
         Ssmax = fa
         clt = ONE
         crt = ONE
         slt = ZERO
         srt = ZERO
      ELSE
         gasmal = .TRUE.
         IF ( ga>fa ) THEN
            pmax = 2
            IF ( (fa/ga)<SLAMCH('EPS') ) THEN
!
!              Case of very large GA
!
               gasmal = .FALSE.
               Ssmax = ga
               IF ( ha>ONE ) THEN
                  Ssmin = fa/(ga/ha)
               ELSE
                  Ssmin = (fa/ga)*ha
               ENDIF
               clt = ONE
               slt = ht/gt
               srt = ONE
               crt = ft/gt
            ENDIF
         ENDIF
         IF ( gasmal ) THEN
!
!           Normal case
!
            d = fa - ha
            IF ( d==fa ) THEN
!
!              Copes with infinite F or H
!
               l = ONE
            ELSE
               l = d/fa
            ENDIF
!
!           Note that 0 .le. L .le. 1
!
            m = gt/ft
!
!           Note that abs(M) .le. 1/macheps
!
            t = TWO - l
!
!           Note that T .ge. 1
!
            mm = m*m
            tt = t*t
            s = SQRT(tt+mm)
!
!           Note that 1 .le. S .le. 1 + 1/macheps
!
            IF ( l==ZERO ) THEN
               r = ABS(m)
            ELSE
               r = SQRT(l*l+mm)
            ENDIF
!
!           Note that 0 .le. R .le. 1 + 1/macheps
!
            a = HALF*(s+r)
!
!           Note that 1 .le. A .le. 1 + abs(M)
!
            Ssmin = ha/a
            Ssmax = fa*a
            IF ( mm/=ZERO ) THEN
               t = (m/(s+t)+m/(r+l))*(ONE+a)
!
!              Note that M is very tiny
!
            ELSEIF ( l==ZERO ) THEN
               t = SIGN(TWO,ft)*SIGN(ONE,gt)
            ELSE
               t = gt/SIGN(d,ft) + m/t
            ENDIF
            l = SQRT(t*t+FOUR)
            crt = TWO/l
            srt = t/l
            clt = (crt+srt*m)/a
            slt = (ht/ft)*srt/a
         ENDIF
      ENDIF
      IF ( swap ) THEN
         Csl = srt
         Snl = crt
         Csr = slt
         Snr = clt
      ELSE
         Csl = clt
         Snl = slt
         Csr = crt
         Snr = srt
      ENDIF
!
!     Correct signs of SSMAX and SSMIN
!
      IF ( pmax==1 ) tsign = SIGN(ONE,Csr)*SIGN(ONE,Csl)*SIGN(ONE,F)
      IF ( pmax==2 ) tsign = SIGN(ONE,Snr)*SIGN(ONE,Csl)*SIGN(ONE,G)
      IF ( pmax==3 ) tsign = SIGN(ONE,Snr)*SIGN(ONE,Snl)*SIGN(ONE,H)
      Ssmax = SIGN(Ssmax,tsign)
      Ssmin = SIGN(Ssmin,tsign*SIGN(ONE,F)*SIGN(ONE,H))
!
!     End of SLASV2
!
      END SUBROUTINE SLASV2
