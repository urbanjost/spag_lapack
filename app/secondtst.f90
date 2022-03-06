!*==aa0004.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!*--AA00043
!> \brief \b SECONDTST
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \date November 2017
!
!> \ingroup auxOTHERcomputational
!
!  =====================================================================      
PROGRAM SECONDTST
      IMPLICIT NONE
!
!  -- LAPACK test routine (version 3.8.0) --
!
!  -- LAPACK computational routine (version 3.8.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2017
!
! =====================================================================
!
!     .. Parameters ..
      INTEGER NMAX , ITS
      PARAMETER (NMAX=1000,ITS=50000)
!     ..
!     .. Local Scalars ..
      INTEGER i , j
      REAL alpha , avg , t1 , t2 , tnosec , total
!     ..
!     .. Local Arrays ..
      REAL x(NMAX) , y(NMAX)
!     ..
!     .. External Functions ..
      REAL SECOND
      EXTERNAL SECOND
!     ..
!     .. External Subroutines ..
      EXTERNAL MYSUB
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC REAL
!     ..
!     .. Executable Statements ..
!
!    .. Figure TOTAL flops ..
      total = REAL(NMAX)*REAL(ITS)*2.0
!
!     Initialize X and Y
!
      DO i = 1 , NMAX
         x(i) = REAL(1)/REAL(i)
         y(i) = REAL(NMAX-i)/REAL(NMAX)
      ENDDO
      alpha = 0.315
!
!     Time TOTAL SAXPY operations
!
      t1 = SECOND()
      DO j = 1 , ITS
         DO i = 1 , NMAX
            y(i) = y(i) + alpha*x(i)
         ENDDO
         alpha = -alpha
      ENDDO
      t2 = SECOND()
      tnosec = t2 - t1
      WRITE (6,99001) total , tnosec
!
99001 FORMAT (' Time for ',G10.3,' SAXPY ops = ',G10.3,' seconds')
      IF ( tnosec>0.0 ) THEN
         WRITE (6,99002) (total/1.0E6)/tnosec
99002    FORMAT (' SAXPY performance rate        = ',G10.3,' mflops ')
      ELSE
         WRITE (6,99003)
99003    FORMAT (' *** Warning:  Time for operations was less or equal',&
     &           ' than zero => timing in TESTING might be dubious')
      ENDIF
!
!     Time TOTAL SAXPY operations with SECOND in the outer loop
!
      t1 = SECOND()
      DO j = 1 , ITS
         DO i = 1 , NMAX
            y(i) = y(i) + alpha*x(i)
         ENDDO
         alpha = -alpha
         t2 = SECOND()
      ENDDO
!
!     Compute the time used in milliseconds used by an average call
!     to SECOND.
!
      WRITE (6,99004) t2 - t1
99004 FORMAT (' Including SECOND, time        = ',G10.3,' seconds')
      avg = ((t2-t1)-tnosec)*1000.0E+00/REAL(ITS)
      IF ( avg>0.0 ) WRITE (6,99005) avg
99005 FORMAT (' Average time for SECOND       = ',G10.3,' milliseconds')
!
!     Compute the equivalent number of floating point operations used
!     by an average call to SECOND.
!
      IF ( (avg>0.0) .AND. (tnosec>0.0) ) WRITE (6,99006) (avg/1000)    &
     &     *total/tnosec
99006 FORMAT (' Equivalent floating point ops = ',G10.3,' ops')
      CALL MYSUB(NMAX,x,y)
END PROGRAM SECONDTST
!*==mysub.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
      SUBROUTINE MYSUB(N,X,Y)
      IMPLICIT NONE
!*--MYSUB123
      INTEGER N
      REAL X(N) , Y(N)
      END SUBROUTINE MYSUB
