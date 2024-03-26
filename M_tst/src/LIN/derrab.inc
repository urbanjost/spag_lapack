!*==derrab.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b derrab
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE derrab( NUNIT )
!
!       .. Scalar Arguments ..
!       INTEGER            NUNIT
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> derrab tests the error exits for DSGESV.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] NUNIT
!> \verbatim
!>          NUNIT is INTEGER
!>          The unit number for output.
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
!> \ingroup double_lin
!
!  =====================================================================
      SUBROUTINE derrab(Nunit)
      IMPLICIT NONE
!*--DERRAB51
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Nunit
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER NMAX
      PARAMETER (NMAX=4)
!     ..
!     .. Local Scalars ..
      INTEGER i , info , iter , j
!     ..
!     .. Local Arrays ..
      INTEGER ip(NMAX)
      DOUBLE PRECISION a(NMAX,NMAX) , af(NMAX,NMAX) , b(NMAX) , c(NMAX) &
     &                 , r(NMAX) , r1(NMAX) , r2(NMAX) , w(2*NMAX) ,    &
     &                 x(NMAX)
      DOUBLE PRECISION work(1)
      REAL swork(1)
!     ..
!     .. External Subroutines ..
EXTERNAL    DSGESV
!     ..
!     .. Scalars in Common ..
      LOGICAL LERr , OK
      CHARACTER*32 SRNamt
      INTEGER INFot , NOUt
!     ..
!     .. Common blocks ..
      COMMON /INFOC / INFot , NOUt , OK , LERr
      COMMON /SRNAMC/ SRNamt
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DBLE
!     ..
!     .. Executable Statements ..
!
      NOUt = Nunit
      WRITE (NOUt,FMT=*)
!
!     Set the variables to innocuous values.
!
      DO j = 1 , NMAX
         DO i = 1 , NMAX
            a(i,j) = 1.D0/DBLE(i+j)
            af(i,j) = 1.D0/DBLE(i+j)
         ENDDO
         b(j) = 0.D0
         r1(j) = 0.D0
         r2(j) = 0.D0
         w(j) = 0.D0
         x(j) = 0.D0
         c(j) = 0.D0
         r(j) = 0.D0
         ip(j) = j
      ENDDO
      OK = .TRUE.
!
      SRNamt = 'DSGESV'
      INFot = 1
      CALL DSGESV(-1,0,a,1,ip,b,1,x,1,work,swork,iter,info)
      CALL CHKXER('DSGESV',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DSGESV(0,-1,a,1,ip,b,1,x,1,work,swork,iter,info)
      CALL CHKXER('DSGESV',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL DSGESV(2,1,a,1,ip,b,2,x,2,work,swork,iter,info)
      CALL CHKXER('DSGESV',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL DSGESV(2,1,a,2,ip,b,1,x,2,work,swork,iter,info)
      CALL CHKXER('DSGESV',INFot,NOUt,LERr,OK)
      INFot = 9
      CALL DSGESV(2,1,a,2,ip,b,2,x,1,work,swork,iter,info)
      CALL CHKXER('DSGESV',INFot,NOUt,LERr,OK)
!
!     Print a summary line.
!
      IF ( OK ) THEN
         WRITE (NOUt,FMT=99001) 'DSGESV'
      ELSE
         WRITE (NOUt,FMT=99002) 'DSGESV'
      ENDIF
!
99001 FORMAT (1X,A6,' drivers passed the tests of the error exits')
99002 FORMAT (' *** ',A6,' drivers failed the tests of the error ',     &
     &        'exits ***')
!
!
!     End of derrab
!
      END SUBROUTINE derrab
