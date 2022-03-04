!*==lsame.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b LSAME
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!      LOGICAL FUNCTION LSAME( CA, CB )
!
!     .. Scalar Arguments ..
!      CHARACTER          CA, CB
!     ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> LSAME returns .TRUE. if CA is the same letter as CB regardless of
!> case.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] CA
!> \verbatim
!> \endverbatim
!>
!> \param[in] CB
!> \verbatim
!>          CA and CB specify the single characters to be compared.
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
!> \ingroup auxOTHERauxiliary
!
!  =====================================================================
      LOGICAL FUNCTION LSAME(Ca,Cb)
      IMPLICIT NONE
!*--LSAME55
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Ca , Cb
!     ..
!
! =====================================================================
!
!     .. Intrinsic Functions ..
      INTRINSIC ICHAR
!     ..
!     .. Local Scalars ..
      INTEGER inta , intb , zcode
!     ..
!     .. Executable Statements ..
!
!     Test if the characters are equal
!
      LSAME = Ca==Cb
      IF ( LSAME ) RETURN
!
!     Now test for equivalence if both characters are alphabetic.
!
      zcode = ICHAR('Z')
!
!     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
!     machines, on which ICHAR returns a value with bit 8 set.
!     ICHAR('A') on Prime machines returns 193 which is the same as
!     ICHAR('A') on an EBCDIC machine.
!
      inta = ICHAR(Ca)
      intb = ICHAR(Cb)
!
      IF ( zcode==90 .OR. zcode==122 ) THEN
!
!        ASCII is assumed - ZCODE is the ASCII code of either lower or
!        upper case 'Z'.
!
         IF ( inta>=97 .AND. inta<=122 ) inta = inta - 32
         IF ( intb>=97 .AND. intb<=122 ) intb = intb - 32
!
      ELSEIF ( zcode==233 .OR. zcode==169 ) THEN
!
!        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
!        upper case 'Z'.
!
         IF ( inta>=129 .AND. inta<=137 .OR. inta>=145 .AND.            &
     &        inta<=153 .OR. inta>=162 .AND. inta<=169 ) inta = inta +  &
     &        64
         IF ( intb>=129 .AND. intb<=137 .OR. intb>=145 .AND.            &
     &        intb<=153 .OR. intb>=162 .AND. intb<=169 ) intb = intb +  &
     &        64
!
      ELSEIF ( zcode==218 .OR. zcode==250 ) THEN
!
!        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
!        plus 128 of either lower or upper case 'Z'.
!
         IF ( inta>=225 .AND. inta<=250 ) inta = inta - 32
         IF ( intb>=225 .AND. intb<=250 ) intb = intb - 32
      ENDIF
      LSAME = inta==intb
!
!     RETURN
!
!     End of LSAME
!
      END FUNCTION LSAME
