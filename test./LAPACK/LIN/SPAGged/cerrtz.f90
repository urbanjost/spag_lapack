!*==cerrtz.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b CERRTZ
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CERRTZ( PATH, NUNIT )
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
!> CERRTZ tests the error exits for CTZRZF.
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
!> \ingroup complex_lin
!
!  =====================================================================
      SUBROUTINE CERRTZ(Path,Nunit)
      IMPLICIT NONE
!*--CERRTZ58
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
      COMPLEX a(NMAX,NMAX) , tau(NMAX) , w(NMAX)
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN
      EXTERNAL LSAMEN
!     ..
!     .. External Subroutines ..
      EXTERNAL ALAESM , CHKXER , CTZRZF
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
      INTRINSIC CMPLX
!     ..
!     .. Executable Statements ..
!
      NOUt = Nunit
      c2 = Path(2:3)
      a(1,1) = CMPLX(1.E+0,-1.E+0)
      a(1,2) = CMPLX(2.E+0,-2.E+0)
      a(2,2) = CMPLX(3.E+0,-3.E+0)
      a(2,1) = CMPLX(4.E+0,-4.E+0)
      w(1) = CMPLX(0.E+0,0.E+0)
      w(2) = CMPLX(0.E+0,0.E+0)
      OK = .TRUE.
!
!     Test error exits for the trapezoidal routines.
!
      WRITE (NOUt,FMT=*)
      IF ( LSAMEN(2,c2,'TZ') ) THEN
!
!        CTZRZF
!
         SRNamt = 'CTZRZF'
         INFot = 1
         CALL CTZRZF(-1,0,a,1,tau,w,1,info)
         CALL CHKXER('CTZRZF',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CTZRZF(1,0,a,1,tau,w,1,info)
         CALL CHKXER('CTZRZF',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CTZRZF(2,2,a,1,tau,w,1,info)
         CALL CHKXER('CTZRZF',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CTZRZF(2,2,a,2,tau,w,0,info)
         CALL CHKXER('CTZRZF',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CTZRZF(2,3,a,2,tau,w,1,info)
         CALL CHKXER('CTZRZF',INFot,NOUt,LERr,OK)
      ENDIF
!
!     Print a summary line.
!
      CALL ALAESM(Path,OK,NOUt)
!
!
!     End of CERRTZ
!
      END SUBROUTINE CERRTZ
