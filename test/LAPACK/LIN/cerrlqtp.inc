!*==cerrlqtp.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b zerrlqtp
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE cerrlqtp( PATH, NUNIT )
!
!       .. Scalar Arguments ..
!       CHARACTER*3        PATH
!       INTEGER            NUNIT
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> cerrlqtp tests the error exits for the complex routines
!> that use the LQT decomposition of a triangular-pentagonal matrix.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] PATH
!> \verbatim
!>          PATH is CHARACTER*3
!>          The LAPACK path name for the routines to be tested.
!> \endverbatim
!>
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
      SUBROUTINE cerrlqtp(Path,Nunit)
      IMPLICIT NONE
!*--CERRLQTP59
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER*3 Path
      INTEGER Nunit
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER NMAX
      PARAMETER (NMAX=2)
!     ..
!     .. Local Scalars ..
      INTEGER i , info , j
!     ..
!     .. Local Arrays ..
      COMPLEX a(NMAX,NMAX) , t(NMAX,NMAX) , w(NMAX) , b(NMAX,NMAX) ,    &
     &        c(NMAX,NMAX)
!     ..
!     .. External Subroutines ..
EXTERNAL    CHKXER,   CTPLQT2,   CTPLQT,   CTPMLQT
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
      INTRINSIC REAL , CMPLX
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
            a(i,j) = 1.E0/CMPLX(REAL(i+j),0.E0)
            c(i,j) = 1.E0/CMPLX(REAL(i+j),0.E0)
            t(i,j) = 1.E0/CMPLX(REAL(i+j),0.E0)
         ENDDO
         w(j) = 0.E0
      ENDDO
      OK = .TRUE.
!
!     Error exits for TPLQT factorization
!
!     CTPLQT
!
      SRNamt = 'CTPLQT'
      INFot = 1
      CALL CTPLQT(-1,1,0,1,a,1,b,1,t,1,w,info)
      CALL CHKXER('CTPLQT',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CTPLQT(1,-1,0,1,a,1,b,1,t,1,w,info)
      CALL CHKXER('CTPLQT',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL CTPLQT(0,1,-1,1,a,1,b,1,t,1,w,info)
      CALL CHKXER('CTPLQT',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL CTPLQT(0,1,1,1,a,1,b,1,t,1,w,info)
      CALL CHKXER('CTPLQT',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL CTPLQT(0,1,0,0,a,1,b,1,t,1,w,info)
      CALL CHKXER('CTPLQT',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL CTPLQT(1,1,0,2,a,1,b,1,t,1,w,info)
      CALL CHKXER('CTPLQT',INFot,NOUt,LERr,OK)
      INFot = 6
      CALL CTPLQT(2,1,0,2,a,1,b,1,t,1,w,info)
      CALL CHKXER('CTPLQT',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL CTPLQT(2,1,0,1,a,2,b,1,t,1,w,info)
      CALL CHKXER('CTPLQT',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL CTPLQT(2,2,1,2,a,2,b,2,t,1,w,info)
      CALL CHKXER('CTPLQT',INFot,NOUt,LERr,OK)
!
!     CTPLQT2
!
      SRNamt = 'CTPLQT2'
      INFot = 1
      CALL CTPLQT2(-1,0,0,a,1,b,1,t,1,info)
      CALL CHKXER('CTPLQT2',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CTPLQT2(0,-1,0,a,1,b,1,t,1,info)
      CALL CHKXER('CTPLQT2',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL CTPLQT2(0,0,-1,a,1,b,1,t,1,info)
      CALL CHKXER('CTPLQT2',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL CTPLQT2(2,2,0,a,1,b,2,t,2,info)
      CALL CHKXER('CTPLQT2',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL CTPLQT2(2,2,0,a,2,b,1,t,2,info)
      CALL CHKXER('CTPLQT2',INFot,NOUt,LERr,OK)
      INFot = 9
      CALL CTPLQT2(2,2,0,a,2,b,2,t,1,info)
      CALL CHKXER('CTPLQT2',INFot,NOUt,LERr,OK)
!
!     CTPMLQT
!
      SRNamt = 'CTPMLQT'
      INFot = 1
      CALL CTPMLQT('/','N',0,0,0,0,1,a,1,t,1,b,1,c,1,w,info)
      CALL CHKXER('CTPMLQT',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CTPMLQT('L','/',0,0,0,0,1,a,1,t,1,b,1,c,1,w,info)
      CALL CHKXER('CTPMLQT',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL CTPMLQT('L','N',-1,0,0,0,1,a,1,t,1,b,1,c,1,w,info)
      CALL CHKXER('CTPMLQT',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL CTPMLQT('L','N',0,-1,0,0,1,a,1,t,1,b,1,c,1,w,info)
      CALL CHKXER('CTPMLQT',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL CTPMLQT('L','N',0,0,-1,0,1,a,1,t,1,b,1,c,1,w,info)
      INFot = 6
      CALL CTPMLQT('L','N',0,0,0,-1,1,a,1,t,1,b,1,c,1,w,info)
      CALL CHKXER('CTPMLQT',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL CTPMLQT('L','N',0,0,0,0,0,a,1,t,1,b,1,c,1,w,info)
      CALL CHKXER('CTPMLQT',INFot,NOUt,LERr,OK)
      INFot = 9
      CALL CTPMLQT('R','N',2,2,2,1,1,a,1,t,1,b,1,c,1,w,info)
      CALL CHKXER('CTPMLQT',INFot,NOUt,LERr,OK)
      INFot = 11
      CALL CTPMLQT('R','N',1,1,1,1,1,a,1,t,0,b,1,c,1,w,info)
      CALL CHKXER('CTPMLQT',INFot,NOUt,LERr,OK)
      INFot = 13
      CALL CTPMLQT('L','N',1,1,1,1,1,a,1,t,1,b,0,c,1,w,info)
      CALL CHKXER('CTPMLQT',INFot,NOUt,LERr,OK)
      INFot = 15
      CALL CTPMLQT('L','N',1,1,1,1,1,a,1,t,1,b,1,c,0,w,info)
      CALL CHKXER('CTPMLQT',INFot,NOUt,LERr,OK)
!
!     Print a summary line.
!
      CALL alaesm(Path,OK,NOUt)
!
!
!     End of cerrlqt
!
      END SUBROUTINE cerrlqtp
