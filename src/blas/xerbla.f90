!*==xerbla.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b XERBLA
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE XERBLA( SRNAME, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER*(*)      SRNAME
!       INTEGER            INFO
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> XERBLA  is an error handler for the LAPACK routines.
!> It is called by an LAPACK routine if an input parameter has an
!> invalid value.  A message is printed and execution stops.
!>
!> Installers may consider modifying the STOP statement in order to
!> call system-specific exception-handling facilities.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SRNAME
!> \verbatim
!>          SRNAME is CHARACTER*(*)
!>          The name of the routine which called XERBLA.
!> \endverbatim
!>
!> \param[in] INFO
!> \verbatim
!>          INFO is INTEGER
!>          The position of the invalid parameter in the parameter list
!>          of the calling routine.
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
!> \ingroup aux_blas
!
!  =====================================================================
      SUBROUTINE XERBLA(Srname,Info)
      IMPLICIT NONE
!*--XERBLA64
!
!  -- Reference BLAS level1 routine (version 3.7.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER*(*) Srname
      INTEGER Info
!     ..
!
! =====================================================================
!
!     .. Intrinsic Functions ..
      INTRINSIC LEN_TRIM
!     ..
!     .. Executable Statements ..
!
      WRITE (*,FMT=99001) Srname(1:LEN_TRIM(Srname)) , Info
!
      STOP
!
99001 FORMAT (' ** On entry to ',A,' parameter number ',I2,' had ',     &
     &        'an illegal value')
!
!     End of XERBLA
!
      END SUBROUTINE XERBLA
