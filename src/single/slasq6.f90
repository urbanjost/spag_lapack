!*==slasq6.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SLASQ6 computes one dqd transform in ping-pong form. Used by sbdsqr and sstegr.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLASQ6 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasq6.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasq6.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasq6.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLASQ6( I0, N0, Z, PP, DMIN, DMIN1, DMIN2, DN,
!                          DNM1, DNM2 )
!
!       .. Scalar Arguments ..
!       INTEGER            I0, N0, PP
!       REAL               DMIN, DMIN1, DMIN2, DN, DNM1, DNM2
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
!> SLASQ6 computes one dqd (shift equal to zero) transform in
!> ping-pong form, with protection against underflow and overflow.
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
      SUBROUTINE SLASQ6(I0,N0,Z,Pp,Dmin,Dmin1,Dmin2,Dn,Dnm1,Dnm2)
      IMPLICIT NONE
!*--SLASQ6122
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER I0 , N0 , Pp
      REAL Dmin , Dmin1 , Dmin2 , Dn , Dnm1 , Dnm2
!     ..
!     .. Array Arguments ..
      REAL Z(*)
!     ..
!
!  =====================================================================
!
!     .. Parameter ..
      REAL ZERO
      PARAMETER (ZERO=0.0E0)
!     ..
!     .. Local Scalars ..
      INTEGER j4 , j4p2
      REAL d , emin , safmin , temp
!     ..
!     .. External Function ..
      REAL SLAMCH
      EXTERNAL SLAMCH
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MIN
!     ..
!     .. Executable Statements ..
!
      IF ( (N0-I0-1)<=0 ) RETURN
!
      safmin = SLAMCH('Safe minimum')
      j4 = 4*I0 + Pp - 3
      emin = Z(j4+4)
      d = Z(j4)
      Dmin = d
!
      IF ( Pp==0 ) THEN
         DO j4 = 4*I0 , 4*(N0-3) , 4
            Z(j4-2) = d + Z(j4-1)
            IF ( Z(j4-2)==ZERO ) THEN
               Z(j4) = ZERO
               d = Z(j4+1)
               Dmin = d
               emin = ZERO
            ELSEIF ( safmin*Z(j4+1)<Z(j4-2) .AND. safmin*Z(j4-2)<Z(j4+1)&
     &               ) THEN
               temp = Z(j4+1)/Z(j4-2)
               Z(j4) = Z(j4-1)*temp
               d = d*temp
            ELSE
               Z(j4) = Z(j4+1)*(Z(j4-1)/Z(j4-2))
               d = Z(j4+1)*(d/Z(j4-2))
            ENDIF
            Dmin = MIN(Dmin,d)
            emin = MIN(emin,Z(j4))
         ENDDO
      ELSE
         DO j4 = 4*I0 , 4*(N0-3) , 4
            Z(j4-3) = d + Z(j4)
            IF ( Z(j4-3)==ZERO ) THEN
               Z(j4-1) = ZERO
               d = Z(j4+2)
               Dmin = d
               emin = ZERO
            ELSEIF ( safmin*Z(j4+2)<Z(j4-3) .AND. safmin*Z(j4-3)<Z(j4+2)&
     &               ) THEN
               temp = Z(j4+2)/Z(j4-3)
               Z(j4-1) = Z(j4)*temp
               d = d*temp
            ELSE
               Z(j4-1) = Z(j4+2)*(Z(j4)/Z(j4-3))
               d = Z(j4+2)*(d/Z(j4-3))
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
      IF ( Z(j4-2)==ZERO ) THEN
         Z(j4) = ZERO
         Dnm1 = Z(j4p2+2)
         Dmin = Dnm1
         emin = ZERO
      ELSEIF ( safmin*Z(j4p2+2)<Z(j4-2) .AND. safmin*Z(j4-2)<Z(j4p2+2) )&
     &         THEN
         temp = Z(j4p2+2)/Z(j4-2)
         Z(j4) = Z(j4p2)*temp
         Dnm1 = Dnm2*temp
      ELSE
         Z(j4) = Z(j4p2+2)*(Z(j4p2)/Z(j4-2))
         Dnm1 = Z(j4p2+2)*(Dnm2/Z(j4-2))
      ENDIF
      Dmin = MIN(Dmin,Dnm1)
!
      Dmin1 = Dmin
      j4 = j4 + 4
      j4p2 = j4 + 2*Pp - 1
      Z(j4-2) = Dnm1 + Z(j4p2)
      IF ( Z(j4-2)==ZERO ) THEN
         Z(j4) = ZERO
         Dn = Z(j4p2+2)
         Dmin = Dn
         emin = ZERO
      ELSEIF ( safmin*Z(j4p2+2)<Z(j4-2) .AND. safmin*Z(j4-2)<Z(j4p2+2) )&
     &         THEN
         temp = Z(j4p2+2)/Z(j4-2)
         Z(j4) = Z(j4p2)*temp
         Dn = Dnm1*temp
      ELSE
         Z(j4) = Z(j4p2+2)*(Z(j4p2)/Z(j4-2))
         Dn = Z(j4p2+2)*(Dnm1/Z(j4-2))
      ENDIF
      Dmin = MIN(Dmin,Dn)
!
      Z(j4+2) = Dn
      Z(4*N0-Pp) = emin
!
!     End of SLASQ6
!
      END SUBROUTINE SLASQ6
