!*==aa0002.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
      IMPLICIT NONE
!*--AA00023
!> \brief \b DSECNDTST
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!      PROGRAM DSECNDTST
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
!> \ingroup auxOTHERauxiliary
!
!  =====================================================================      PROGRAM DSECNDTST
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
      DOUBLE PRECISION alpha , avg , t1 , t2 , tnosec , total
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION x(NMAX) , y(NMAX)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DSECND
      EXTERNAL DSECND
!     ..
!     .. External Subroutines ..
      EXTERNAL MYSUB
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DBLE
!     ..
!     .. Executable Statements ..
!
!    .. Figure TOTAL flops ..
      total = DBLE(NMAX)*DBLE(ITS)*2.0
!
!     Initialize X and Y
!
      DO i = 1 , NMAX
         x(i) = DBLE(1)/DBLE(i)
         y(i) = DBLE(NMAX-i)/DBLE(NMAX)
      ENDDO
      alpha = 0.315D0
!
!     Time TOTAL SAXPY operations
!
      t1 = DSECND()
      DO j = 1 , ITS
         DO i = 1 , NMAX
            y(i) = y(i) + alpha*x(i)
         ENDDO
         alpha = -alpha
      ENDDO
      t2 = DSECND()
      tnosec = t2 - t1
      WRITE (6,99001) total , tnosec
!
99001 FORMAT (' Time for ',G10.3,' DAXPY ops = ',G10.3,' seconds')
      IF ( tnosec>0.0 ) THEN
         WRITE (6,99002) (total/1.0D6)/tnosec
99002    FORMAT (' DAXPY performance rate        = ',G10.3,' mflops ')
      ELSE
         WRITE (6,99003)
99003    FORMAT (' *** Warning:  Time for operations was less or equal',&
     &           ' than zero => timing in TESTING might be dubious')
      ENDIF
!
!     Time TOTAL DAXPY operations with DSECND in the outer loop
!
      t1 = DSECND()
      DO j = 1 , ITS
         DO i = 1 , NMAX
            y(i) = y(i) + alpha*x(i)
         ENDDO
         alpha = -alpha
         t2 = DSECND()
      ENDDO
!
!     Compute the time used in milliseconds used by an average call
!     to DSECND.
!
      WRITE (6,99004) t2 - t1
99004 FORMAT (' Including DSECND, time        = ',G10.3,' seconds')
      avg = ((t2-t1)-tnosec)*1000.0D+00/DBLE(ITS)
      IF ( avg>0.0 ) WRITE (6,99005) avg
99005 FORMAT (' Average time for DSECND       = ',G10.3,' milliseconds')
!
!     Compute the equivalent number of floating point operations used
!     by an average call to DSECND.
!
      IF ( (avg>0.0) .AND. (tnosec>0.0) ) WRITE (6,99006) (avg/1000)    &
     &     *total/tnosec
99006 FORMAT (' Equivalent floating point ops = ',G10.3,' ops')
      CALL MYSUB(NMAX,x,y)
      END PROGRAM AA0002
!*==mysub.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
      SUBROUTINE MYSUB(N,X,Y)
      IMPLICIT NONE
!*--MYSUB127
      INTEGER N
      DOUBLE PRECISION X(N) , Y(N)
      END SUBROUTINE MYSUB
