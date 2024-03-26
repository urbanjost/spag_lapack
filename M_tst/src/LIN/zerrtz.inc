!*==zerrtz.f90  processed by SPAG 7.51RB at 17:35 on  4 Mar 2022
!> \brief \b zerrtz
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE zerrtz( PATH, NUNIT )
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
!> zerrtz tests the error exits for ZTZRZF.
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
!> \ingroup complex16_lin
!
!  =====================================================================
      SUBROUTINE zerrtz(Path,Nunit)
      IMPLICIT NONE
!*--ZERRTZ58
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
      CHARACTER*2 c2
      INTEGER info
!     ..
!     .. Local Arrays ..
      COMPLEX*16 a(NMAX,NMAX) , tau(NMAX) , w(NMAX)
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN
EXTERNAL    LSAMEN
!     ..
!     .. External Subroutines ..
EXTERNAL    ZTZRZF
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
      INTRINSIC DCMPLX
!     ..
!     .. Executable Statements ..
!
      NOUt = Nunit
      c2 = Path(2:3)
      a(1,1) = DCMPLX(1.D+0,-1.D+0)
      a(1,2) = DCMPLX(2.D+0,-2.D+0)
      a(2,2) = DCMPLX(3.D+0,-3.D+0)
      a(2,1) = DCMPLX(4.D+0,-4.D+0)
      w(1) = DCMPLX(0.D+0,0.D+0)
      w(2) = DCMPLX(0.D+0,0.D+0)
      OK = .TRUE.
!
!     Test error exits for the trapezoidal routines.
      WRITE (NOUt,FMT=*)
      IF ( LSAMEN(2,c2,'TZ') ) THEN
!
!
!        ZTZRZF
!
         SRNamt = 'ZTZRZF'
         INFot = 1
         CALL ZTZRZF(-1,0,a,1,tau,w,1,info)
         CALL CHKXER('ZTZRZF',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZTZRZF(1,0,a,1,tau,w,1,info)
         CALL CHKXER('ZTZRZF',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZTZRZF(2,2,a,1,tau,w,1,info)
         CALL CHKXER('ZTZRZF',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZTZRZF(2,2,a,2,tau,w,0,info)
         CALL CHKXER('ZTZRZF',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZTZRZF(2,3,a,2,tau,w,1,info)
         CALL CHKXER('ZTZRZF',INFot,NOUt,LERr,OK)
      ENDIF
!
!     Print a summary line.
!
      CALL alaesm(Path,OK,NOUt)
!
!
!     End of zerrtz
!
      END SUBROUTINE zerrtz
