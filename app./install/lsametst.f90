PROGRAM LSAMETST
!*==aa0003.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
      IMPLICIT NONE
!*--AA00033
!> \brief \b LSAMETST
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!      PROGRAM LSAMETST
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
!  =====================================================================      PROGRAM LSAMETST
!
!  -- LAPACK test routine (version 3.7.0) --
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!  =====================================================================
!     .. Local Scalars ..
      INTEGER i1 , i2
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ICHAR
!     ..
!     .. Executable Statements ..
!
!
!     Determine the character set.
!
      i1 = ICHAR('A')
      i2 = ICHAR('a')
      IF ( i2-i1==32 ) THEN
         WRITE (*,*) ' ASCII character set'
      ELSE
         WRITE (*,*) ' Non-ASCII character set, IOFF should be ' ,      &
     &               i2 - i1
      ENDIF
!
!     Test LSAME.
!
      IF ( .NOT.LSAME('A','A') ) WRITE (*,99001) 'A' , 'A'
      IF ( .NOT.LSAME('A','a') ) WRITE (*,99001) 'A' , 'a'
      IF ( .NOT.LSAME('a','A') ) WRITE (*,99001) 'a' , 'A'
      IF ( .NOT.LSAME('a','a') ) WRITE (*,99001) 'a' , 'a'
      IF ( LSAME('A','B') ) WRITE (*,99002) 'A' , 'B'
      IF ( LSAME('A','b') ) WRITE (*,99002) 'A' , 'b'
      IF ( LSAME('a','B') ) WRITE (*,99002) 'a' , 'B'
      IF ( LSAME('a','b') ) WRITE (*,99002) 'a' , 'b'
      IF ( LSAME('O','/') ) WRITE (*,99002) 'O' , '/'
      IF ( LSAME('/','O') ) WRITE (*,99002) '/' , 'O'
      IF ( LSAME('o','/') ) WRITE (*,99002) 'o' , '/'
      IF ( LSAME('/','o') ) WRITE (*,99002) '/' , 'o'
      WRITE (*,*) ' Tests completed'
!
99001 FORMAT (' *** Error:  LSAME( ',A1,', ',A1,') is .FALSE.')
99002 FORMAT (' *** Error:  LSAME( ',A1,', ',A1,') is .TRUE.')
END PROGRAM LSAMETST
