!*==tstiee.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b TSTIEE
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
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
      PROGRAM TSTIEE
      IMPLICIT NONE
!*--TSTIEE24
!
!  -- LAPACK test routine (version 3.7.0) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. External Functions ..
      INTEGER ILAENV
      EXTERNAL ILAENV
!     ..
!     .. Local Scalars ..
      INTEGER ieeeok
!     ..
!     .. Executable Statements ..
!
      WRITE (6,FMT=*)                                                   &
     &               'We are about to check whether infinity arithmetic'
      WRITE (6,FMT=*) 'can be trusted.  If this test hangs, set'
      WRITE (6,FMT=*) 'ILAENV = 0 for ISPEC = 11 in LAPACK/SRC/ilaenv.f'
!
      ieeeok = ILAENV(11,'ILAENV','N',1,2,3,4)
      WRITE (6,FMT=*)
!
      IF ( ieeeok==0 ) THEN
         WRITE (6,FMT=*)                                                &
     &           'Infinity arithmetic did not perform per the ieee spec'
      ELSE
         WRITE (6,FMT=*)                                                &
     &             'Infinity arithmetic performed as per the ieee spec.'
         WRITE (6,FMT=*)                                                &
     &            'However, this is not an exhaustive test and does not'
         WRITE (6,FMT=*) 'guarantee that infinity arithmetic meets the' &
     &                   , ' ieee spec.'
      ENDIF
!
      WRITE (6,FMT=*)
!     ILAENV( 10, ...) checks both infinity and NaN arithmetic
!     infinity has already been checked so checking NaN now
      WRITE (6,FMT=*) 'We are about to check whether NaN arithmetic'
      WRITE (6,FMT=*) 'can be trusted.  If this test hangs, set'
      WRITE (6,FMT=*) 'ILAENV = 0 for ISPEC = 10 in LAPACK/SRC/ilaenv.f'
      ieeeok = ILAENV(10,'ILAENV','N',1,2,3,4)
!
      WRITE (6,FMT=*)
      IF ( ieeeok==0 ) THEN
         WRITE (6,FMT=*)                                                &
     &                'NaN arithmetic did not perform per the ieee spec'
      ELSE
         WRITE (6,FMT=*) 'NaN arithmetic performed as per the ieee' ,   &
     &                   ' spec.'
         WRITE (6,FMT=*)                                                &
     &            'However, this is not an exhaustive test and does not'
         WRITE (6,FMT=*) 'guarantee that NaN arithmetic meets the' ,    &
     &                   ' ieee spec.'
      ENDIF
      WRITE (6,FMT=*)
!
      END PROGRAM TSTIEE
