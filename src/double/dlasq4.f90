!*==dlasq4.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b DLASQ4 computes an approximation to the smallest eigenvalue using values of d from the previous transform. Used by sbdsqr.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLASQ4 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasq4.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlasq4.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlasq4.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLASQ4( I0, N0, Z, PP, N0IN, DMIN, DMIN1, DMIN2, DN,
!                          DN1, DN2, TAU, TTYPE, G )
!
!       .. Scalar Arguments ..
!       INTEGER            I0, N0, N0IN, PP, TTYPE
!       DOUBLE PRECISION   DMIN, DMIN1, DMIN2, DN, DN1, DN2, G, TAU
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   Z( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLASQ4 computes an approximation TAU to the smallest eigenvalue
!> using values of d from the previous transform.
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
!>          Z is DOUBLE PRECISION array, dimension ( 4*N0 )
!>        Z holds the qd array.
!> \endverbatim
!>
!> \param[in] PP
!> \verbatim
!>          PP is INTEGER
!>        PP=0 for ping, PP=1 for pong.
!> \endverbatim
!>
!> \param[in] N0IN
!> \verbatim
!>          N0IN is INTEGER
!>        The value of N0 at start of EIGTEST.
!> \endverbatim
!>
!> \param[in] DMIN
!> \verbatim
!>          DMIN is DOUBLE PRECISION
!>        Minimum value of d.
!> \endverbatim
!>
!> \param[in] DMIN1
!> \verbatim
!>          DMIN1 is DOUBLE PRECISION
!>        Minimum value of d, excluding D( N0 ).
!> \endverbatim
!>
!> \param[in] DMIN2
!> \verbatim
!>          DMIN2 is DOUBLE PRECISION
!>        Minimum value of d, excluding D( N0 ) and D( N0-1 ).
!> \endverbatim
!>
!> \param[in] DN
!> \verbatim
!>          DN is DOUBLE PRECISION
!>        d(N)
!> \endverbatim
!>
!> \param[in] DN1
!> \verbatim
!>          DN1 is DOUBLE PRECISION
!>        d(N-1)
!> \endverbatim
!>
!> \param[in] DN2
!> \verbatim
!>          DN2 is DOUBLE PRECISION
!>        d(N-2)
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is DOUBLE PRECISION
!>        This is the shift.
!> \endverbatim
!>
!> \param[out] TTYPE
!> \verbatim
!>          TTYPE is INTEGER
!>        Shift type.
!> \endverbatim
!>
!> \param[in,out] G
!> \verbatim
!>          G is DOUBLE PRECISION
!>        G is passed as an argument in order to save its value between
!>        calls to DLASQ4.
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
!> \date June 2016
!
!> \ingroup auxOTHERcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  CNST1 = 9/16
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE DLASQ4(I0,N0,Z,Pp,N0in,Dmin,Dmin1,Dmin2,Dn,Dn1,Dn2,Tau,&
     &                  Ttype,G)
      IMPLICIT NONE
!*--DLASQ4155
!
!  -- LAPACK computational routine (version 3.7.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!
!     .. Scalar Arguments ..
      INTEGER I0 , N0 , N0in , Pp , Ttype
      DOUBLE PRECISION Dmin , Dmin1 , Dmin2 , Dn , Dn1 , Dn2 , G , Tau
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION Z(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION CNST1 , CNST2 , CNST3
      PARAMETER (CNST1=0.5630D0,CNST2=1.010D0,CNST3=1.050D0)
      DOUBLE PRECISION QURTR , THIRD , HALF , ZERO , ONE , TWO , HUNDRD
      PARAMETER (QURTR=0.250D0,THIRD=0.3330D0,HALF=0.50D0,ZERO=0.0D0,   &
     &           ONE=1.0D0,TWO=2.0D0,HUNDRD=100.0D0)
!     ..
!     .. Local Scalars ..
      INTEGER i4 , nn , np
      DOUBLE PRECISION a2 , b1 , b2 , gam , gap1 , gap2 , s
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN , SQRT
!     ..
!     .. Executable Statements ..
!
!     A negative DMIN forces the shift to take that absolute value
!     TTYPE records the type of shift.
!
      IF ( Dmin<=ZERO ) THEN
         Tau = -Dmin
         Ttype = -1
         RETURN
      ENDIF
!
      nn = 4*N0 + Pp
      IF ( N0in==N0 ) THEN
!
!        No eigenvalues deflated.
!
         IF ( Dmin==Dn .OR. Dmin==Dn1 ) THEN
!
            b1 = SQRT(Z(nn-3))*SQRT(Z(nn-5))
            b2 = SQRT(Z(nn-7))*SQRT(Z(nn-9))
            a2 = Z(nn-7) + Z(nn-5)
!
!           Cases 2 and 3.
!
            IF ( Dmin==Dn .AND. Dmin1==Dn1 ) THEN
               gap2 = Dmin2 - a2 - Dmin2*QURTR
               IF ( gap2>ZERO .AND. gap2>b2 ) THEN
                  gap1 = a2 - Dn - (b2/gap2)*b2
               ELSE
                  gap1 = a2 - Dn - (b1+b2)
               ENDIF
               IF ( gap1>ZERO .AND. gap1>b1 ) THEN
                  s = MAX(Dn-(b1/gap1)*b1,HALF*Dmin)
                  Ttype = -2
               ELSE
                  s = ZERO
                  IF ( Dn>b1 ) s = Dn - b1
                  IF ( a2>(b1+b2) ) s = MIN(s,a2-(b1+b2))
                  s = MAX(s,THIRD*Dmin)
                  Ttype = -3
               ENDIF
            ELSE
!
!              Case 4.
!
               Ttype = -4
               s = QURTR*Dmin
               IF ( Dmin==Dn ) THEN
                  gam = Dn
                  a2 = ZERO
                  IF ( Z(nn-5)>Z(nn-7) ) RETURN
                  b2 = Z(nn-5)/Z(nn-7)
                  np = nn - 9
               ELSE
                  np = nn - 2*Pp
                  gam = Dn1
                  IF ( Z(np-4)>Z(np-2) ) RETURN
                  a2 = Z(np-4)/Z(np-2)
                  IF ( Z(nn-9)>Z(nn-11) ) RETURN
                  b2 = Z(nn-9)/Z(nn-11)
                  np = nn - 13
               ENDIF
!
!              Approximate contribution to norm squared from I < NN-1.
!
               a2 = a2 + b2
               DO i4 = np , 4*I0 - 1 + Pp , -4
                  IF ( b2==ZERO ) EXIT
                  b1 = b2
                  IF ( Z(i4)>Z(i4-2) ) RETURN
                  b2 = b2*(Z(i4)/Z(i4-2))
                  a2 = a2 + b2
                  IF ( HUNDRD*MAX(b2,b1)<a2 .OR. CNST1<a2 ) EXIT
               ENDDO
               a2 = CNST3*a2
!
!              Rayleigh quotient residual bound.
!
               IF ( a2<CNST1 ) s = gam*(ONE-SQRT(a2))/(ONE+a2)
            ENDIF
         ELSEIF ( Dmin==Dn2 ) THEN
!
!           Case 5.
!
            Ttype = -5
            s = QURTR*Dmin
!
!           Compute contribution to norm squared from I > NN-2.
!
            np = nn - 2*Pp
            b1 = Z(np-2)
            b2 = Z(np-6)
            gam = Dn2
            IF ( Z(np-8)>b2 .OR. Z(np-4)>b1 ) RETURN
            a2 = (Z(np-8)/b2)*(ONE+Z(np-4)/b1)
!
!           Approximate contribution to norm squared from I < NN-2.
!
            IF ( N0-I0>2 ) THEN
               b2 = Z(nn-13)/Z(nn-15)
               a2 = a2 + b2
               DO i4 = nn - 17 , 4*I0 - 1 + Pp , -4
                  IF ( b2==ZERO ) EXIT
                  b1 = b2
                  IF ( Z(i4)>Z(i4-2) ) RETURN
                  b2 = b2*(Z(i4)/Z(i4-2))
                  a2 = a2 + b2
                  IF ( HUNDRD*MAX(b2,b1)<a2 .OR. CNST1<a2 ) EXIT
               ENDDO
               a2 = CNST3*a2
            ENDIF
!
            IF ( a2<CNST1 ) s = gam*(ONE-SQRT(a2))/(ONE+a2)
         ELSE
!
!           Case 6, no information to guide us.
!
            IF ( Ttype==-6 ) THEN
               G = G + THIRD*(ONE-G)
            ELSEIF ( Ttype==-18 ) THEN
               G = QURTR*THIRD
            ELSE
               G = QURTR
            ENDIF
            s = G*Dmin
            Ttype = -6
         ENDIF
!
      ELSEIF ( N0in==(N0+1) ) THEN
!
!        One eigenvalue just deflated. Use DMIN1, DN1 for DMIN and DN.
!
         IF ( Dmin1==Dn1 .AND. Dmin2==Dn2 ) THEN
!
!           Cases 7 and 8.
!
            Ttype = -7
            s = THIRD*Dmin1
            IF ( Z(nn-5)>Z(nn-7) ) RETURN
            b1 = Z(nn-5)/Z(nn-7)
            b2 = b1
            IF ( b2/=ZERO ) THEN
               DO i4 = 4*N0 - 9 + Pp , 4*I0 - 1 + Pp , -4
                  a2 = b1
                  IF ( Z(i4)>Z(i4-2) ) RETURN
                  b1 = b1*(Z(i4)/Z(i4-2))
                  b2 = b2 + b1
                  IF ( HUNDRD*MAX(b1,a2)<b2 ) EXIT
               ENDDO
            ENDIF
            b2 = SQRT(CNST3*b2)
            a2 = Dmin1/(ONE+b2**2)
            gap2 = HALF*Dmin2 - a2
            IF ( gap2>ZERO .AND. gap2>b2*a2 ) THEN
               s = MAX(s,a2*(ONE-CNST2*a2*(b2/gap2)*b2))
            ELSE
               s = MAX(s,a2*(ONE-CNST2*b2))
               Ttype = -8
            ENDIF
         ELSE
!
!           Case 9.
!
            s = QURTR*Dmin1
            IF ( Dmin1==Dn1 ) s = HALF*Dmin1
            Ttype = -9
         ENDIF
!
      ELSEIF ( N0in==(N0+2) ) THEN
!
!        Two eigenvalues deflated. Use DMIN2, DN2 for DMIN and DN.
!
!        Cases 10 and 11.
!
         IF ( Dmin2==Dn2 .AND. TWO*Z(nn-5)<Z(nn-7) ) THEN
            Ttype = -10
            s = THIRD*Dmin2
            IF ( Z(nn-5)>Z(nn-7) ) RETURN
            b1 = Z(nn-5)/Z(nn-7)
            b2 = b1
            IF ( b2/=ZERO ) THEN
               DO i4 = 4*N0 - 9 + Pp , 4*I0 - 1 + Pp , -4
                  IF ( Z(i4)>Z(i4-2) ) RETURN
                  b1 = b1*(Z(i4)/Z(i4-2))
                  b2 = b2 + b1
                  IF ( HUNDRD*b1<b2 ) EXIT
               ENDDO
            ENDIF
            b2 = SQRT(CNST3*b2)
            a2 = Dmin2/(ONE+b2**2)
            gap2 = Z(nn-7) + Z(nn-9) - SQRT(Z(nn-11))*SQRT(Z(nn-9)) - a2
            IF ( gap2>ZERO .AND. gap2>b2*a2 ) THEN
               s = MAX(s,a2*(ONE-CNST2*a2*(b2/gap2)*b2))
            ELSE
               s = MAX(s,a2*(ONE-CNST2*b2))
            ENDIF
         ELSE
            s = QURTR*Dmin2
            Ttype = -11
         ENDIF
      ELSEIF ( N0in>(N0+2) ) THEN
!
!        Case 12, more than two eigenvalues deflated. No information.
!
         s = ZERO
         Ttype = -12
      ENDIF
!
      Tau = s
!
!     End of DLASQ4
!
      END SUBROUTINE DLASQ4
