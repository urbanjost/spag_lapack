!*==slasq5.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief <b> SLASQ5 computes one dqds transform in ping-pong form. Used by sbdsqr and sstegr. </b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLASQ5 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasq5.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasq5.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasq5.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLASQ5( I0, N0, Z, PP, TAU, SIGMA, DMIN, DMIN1, DMIN2, DN,
!                          DNM1, DNM2, IEEE, EPS )
!
!       .. Scalar Arguments ..
!       LOGICAL            IEEE
!       INTEGER            I0, N0, PP
!       REAL               EPS, DMIN, DMIN1, DMIN2, DN, DNM1, DNM2, SIGMA, TAU
!       ..
!       .. Array Arguments ..
!       REAL               Z( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLASQ5 computes one dqds transform in ping-pong form, one
!> version for IEEE machines another for non IEEE machines.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] I0
!> \verbatim
!>          I0 is INTEGER
!>        First index.
!> \endverbatim
!>
!> \param[in] N0
!> \verbatim
!>          N0 is INTEGER
!>        Last index.
!> \endverbatim
!>
!> \param[in] Z
!> \verbatim
!>          Z is REAL array, dimension ( 4*N )
!>        Z holds the qd array. EMIN is stored in Z(4*N0) to avoid
!>        an extra argument.
!> \endverbatim
!>
!> \param[in] PP
!> \verbatim
!>          PP is INTEGER
!>        PP=0 for ping, PP=1 for pong.
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is REAL
!>        This is the shift.
!> \endverbatim
!>
!> \param[in] SIGMA
!> \verbatim
!>          SIGMA is REAL
!>        This is the accumulated shift up to this step.
!> \endverbatim
!>
!> \param[out] DMIN
!> \verbatim
!>          DMIN is REAL
!>        Minimum value of d.
!> \endverbatim
!>
!> \param[out] DMIN1
!> \verbatim
!>          DMIN1 is REAL
!>        Minimum value of d, excluding D( N0 ).
!> \endverbatim
!>
!> \param[out] DMIN2
!> \verbatim
!>          DMIN2 is REAL
!>        Minimum value of d, excluding D( N0 ) and D( N0-1 ).
!> \endverbatim
!>
!> \param[out] DN
!> \verbatim
!>          DN is REAL
!>        d(N0), the last value of d.
!> \endverbatim
!>
!> \param[out] DNM1
!> \verbatim
!>          DNM1 is REAL
!>        d(N0-1).
!> \endverbatim
!>
!> \param[out] DNM2
!> \verbatim
!>          DNM2 is REAL
!>        d(N0-2).
!> \endverbatim
!>
!> \param[in] IEEE
!> \verbatim
!>          IEEE is LOGICAL
!>        Flag for IEEE or non IEEE arithmetic.
!> \endverbatim
!>
!> \param[in] EPS
!> \verbatim
!>         EPS is REAL
!>        This is the value of epsilon used.
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
!> \ingroup auxOTHERcomputational
!
!  =====================================================================
      SUBROUTINE SLASQ5(I0,N0,Z,Pp,Tau,Sigma,Dmin,Dmin1,Dmin2,Dn,Dnm1,  &
     &                  Dnm2,Ieee,Eps)
      IMPLICIT NONE
!*--SLASQ5148
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      LOGICAL Ieee
      INTEGER I0 , N0 , Pp
      REAL Dmin , Dmin1 , Dmin2 , Dn , Dnm1 , Dnm2 , Tau , Sigma , Eps
!     ..
!     .. Array Arguments ..
      REAL Z(*)
!     ..
!
!  =====================================================================
!
!     .. Parameter ..
      REAL ZERO , HALF
      PARAMETER (ZERO=0.0E0,HALF=0.5)
!     ..
!     .. Local Scalars ..
      INTEGER j4 , j4p2
      REAL d , emin , temp , dthresh
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MIN
!     ..
!     .. Executable Statements ..
!
      IF ( (N0-I0-1)<=0 ) RETURN
!
      dthresh = Eps*(Sigma+Tau)
      IF ( Tau<dthresh*HALF ) Tau = ZERO
      IF ( Tau/=ZERO ) THEN
         j4 = 4*I0 + Pp - 3
         emin = Z(j4+4)
         d = Z(j4) - Tau
         Dmin = d
         Dmin1 = -Z(j4)
!
         IF ( Ieee ) THEN
!
!     Code for IEEE arithmetic.
!
            IF ( Pp==0 ) THEN
               DO j4 = 4*I0 , 4*(N0-3) , 4
                  Z(j4-2) = d + Z(j4-1)
                  temp = Z(j4+1)/Z(j4-2)
                  d = d*temp - Tau
                  Dmin = MIN(Dmin,d)
                  Z(j4) = Z(j4-1)*temp
                  emin = MIN(Z(j4),emin)
               ENDDO
            ELSE
               DO j4 = 4*I0 , 4*(N0-3) , 4
                  Z(j4-3) = d + Z(j4)
                  temp = Z(j4+2)/Z(j4-3)
                  d = d*temp - Tau
                  Dmin = MIN(Dmin,d)
                  Z(j4-1) = Z(j4)*temp
                  emin = MIN(Z(j4-1),emin)
               ENDDO
            ENDIF
!
!     Unroll last two steps.
!
            Dnm2 = d
            Dmin2 = Dmin
            j4 = 4*(N0-2) - Pp
            j4p2 = j4 + 2*Pp - 1
            Z(j4-2) = Dnm2 + Z(j4p2)
            Z(j4) = Z(j4p2+2)*(Z(j4p2)/Z(j4-2))
            Dnm1 = Z(j4p2+2)*(Dnm2/Z(j4-2)) - Tau
            Dmin = MIN(Dmin,Dnm1)
!
            Dmin1 = Dmin
            j4 = j4 + 4
            j4p2 = j4 + 2*Pp - 1
            Z(j4-2) = Dnm1 + Z(j4p2)
            Z(j4) = Z(j4p2+2)*(Z(j4p2)/Z(j4-2))
            Dn = Z(j4p2+2)*(Dnm1/Z(j4-2)) - Tau
            Dmin = MIN(Dmin,Dn)
!
         ELSE
!
!     Code for non IEEE arithmetic.
!
            IF ( Pp==0 ) THEN
               DO j4 = 4*I0 , 4*(N0-3) , 4
                  Z(j4-2) = d + Z(j4-1)
                  IF ( d<ZERO ) THEN
                     RETURN
                  ELSE
                     Z(j4) = Z(j4+1)*(Z(j4-1)/Z(j4-2))
                     d = Z(j4+1)*(d/Z(j4-2)) - Tau
                  ENDIF
                  Dmin = MIN(Dmin,d)
                  emin = MIN(emin,Z(j4))
               ENDDO
            ELSE
               DO j4 = 4*I0 , 4*(N0-3) , 4
                  Z(j4-3) = d + Z(j4)
                  IF ( d<ZERO ) THEN
                     RETURN
                  ELSE
                     Z(j4-1) = Z(j4+2)*(Z(j4)/Z(j4-3))
                     d = Z(j4+2)*(d/Z(j4-3)) - Tau
                  ENDIF
                  Dmin = MIN(Dmin,d)
                  emin = MIN(emin,Z(j4-1))
               ENDDO
            ENDIF
!
!     Unroll last two steps.
!
            Dnm2 = d
            Dmin2 = Dmin
            j4 = 4*(N0-2) - Pp
            j4p2 = j4 + 2*Pp - 1
            Z(j4-2) = Dnm2 + Z(j4p2)
            IF ( Dnm2<ZERO ) THEN
               RETURN
            ELSE
               Z(j4) = Z(j4p2+2)*(Z(j4p2)/Z(j4-2))
               Dnm1 = Z(j4p2+2)*(Dnm2/Z(j4-2)) - Tau
            ENDIF
            Dmin = MIN(Dmin,Dnm1)
!
            Dmin1 = Dmin
            j4 = j4 + 4
            j4p2 = j4 + 2*Pp - 1
            Z(j4-2) = Dnm1 + Z(j4p2)
            IF ( Dnm1<ZERO ) THEN
               RETURN
            ELSE
               Z(j4) = Z(j4p2+2)*(Z(j4p2)/Z(j4-2))
               Dn = Z(j4p2+2)*(Dnm1/Z(j4-2)) - Tau
            ENDIF
            Dmin = MIN(Dmin,Dn)
!
         ENDIF
!
      ELSE
!     This is the version that sets d's to zero if they are small enough
         j4 = 4*I0 + Pp - 3
         emin = Z(j4+4)
         d = Z(j4) - Tau
         Dmin = d
         Dmin1 = -Z(j4)
         IF ( Ieee ) THEN
!
!     Code for IEEE arithmetic.
!
            IF ( Pp==0 ) THEN
               DO j4 = 4*I0 , 4*(N0-3) , 4
                  Z(j4-2) = d + Z(j4-1)
                  temp = Z(j4+1)/Z(j4-2)
                  d = d*temp - Tau
                  IF ( d<dthresh ) d = ZERO
                  Dmin = MIN(Dmin,d)
                  Z(j4) = Z(j4-1)*temp
                  emin = MIN(Z(j4),emin)
               ENDDO
            ELSE
               DO j4 = 4*I0 , 4*(N0-3) , 4
                  Z(j4-3) = d + Z(j4)
                  temp = Z(j4+2)/Z(j4-3)
                  d = d*temp - Tau
                  IF ( d<dthresh ) d = ZERO
                  Dmin = MIN(Dmin,d)
                  Z(j4-1) = Z(j4)*temp
                  emin = MIN(Z(j4-1),emin)
               ENDDO
            ENDIF
!
!     Unroll last two steps.
!
            Dnm2 = d
            Dmin2 = Dmin
            j4 = 4*(N0-2) - Pp
            j4p2 = j4 + 2*Pp - 1
            Z(j4-2) = Dnm2 + Z(j4p2)
            Z(j4) = Z(j4p2+2)*(Z(j4p2)/Z(j4-2))
            Dnm1 = Z(j4p2+2)*(Dnm2/Z(j4-2)) - Tau
            Dmin = MIN(Dmin,Dnm1)
!
            Dmin1 = Dmin
            j4 = j4 + 4
            j4p2 = j4 + 2*Pp - 1
            Z(j4-2) = Dnm1 + Z(j4p2)
            Z(j4) = Z(j4p2+2)*(Z(j4p2)/Z(j4-2))
            Dn = Z(j4p2+2)*(Dnm1/Z(j4-2)) - Tau
            Dmin = MIN(Dmin,Dn)
!
         ELSE
!
!     Code for non IEEE arithmetic.
!
            IF ( Pp==0 ) THEN
               DO j4 = 4*I0 , 4*(N0-3) , 4
                  Z(j4-2) = d + Z(j4-1)
                  IF ( d<ZERO ) THEN
                     RETURN
                  ELSE
                     Z(j4) = Z(j4+1)*(Z(j4-1)/Z(j4-2))
                     d = Z(j4+1)*(d/Z(j4-2)) - Tau
                  ENDIF
                  IF ( d<dthresh ) d = ZERO
                  Dmin = MIN(Dmin,d)
                  emin = MIN(emin,Z(j4))
               ENDDO
            ELSE
               DO j4 = 4*I0 , 4*(N0-3) , 4
                  Z(j4-3) = d + Z(j4)
                  IF ( d<ZERO ) THEN
                     RETURN
                  ELSE
                     Z(j4-1) = Z(j4+2)*(Z(j4)/Z(j4-3))
                     d = Z(j4+2)*(d/Z(j4-3)) - Tau
                  ENDIF
                  IF ( d<dthresh ) d = ZERO
                  Dmin = MIN(Dmin,d)
                  emin = MIN(emin,Z(j4-1))
               ENDDO
            ENDIF
!
!     Unroll last two steps.
!
            Dnm2 = d
            Dmin2 = Dmin
            j4 = 4*(N0-2) - Pp
            j4p2 = j4 + 2*Pp - 1
            Z(j4-2) = Dnm2 + Z(j4p2)
            IF ( Dnm2<ZERO ) THEN
               RETURN
            ELSE
               Z(j4) = Z(j4p2+2)*(Z(j4p2)/Z(j4-2))
               Dnm1 = Z(j4p2+2)*(Dnm2/Z(j4-2)) - Tau
            ENDIF
            Dmin = MIN(Dmin,Dnm1)
!
            Dmin1 = Dmin
            j4 = j4 + 4
            j4p2 = j4 + 2*Pp - 1
            Z(j4-2) = Dnm1 + Z(j4p2)
            IF ( Dnm1<ZERO ) THEN
               RETURN
            ELSE
               Z(j4) = Z(j4p2+2)*(Z(j4p2)/Z(j4-2))
               Dn = Z(j4p2+2)*(Dnm1/Z(j4-2)) - Tau
            ENDIF
            Dmin = MIN(Dmin,Dn)
!
         ENDIF
!
      ENDIF
      Z(j4+2) = Dn
      Z(4*N0-Pp) = emin
!
!     End of SLASQ5
!
      END SUBROUTINE SLASQ5
