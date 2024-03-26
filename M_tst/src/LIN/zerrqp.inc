!*==zerrqp.f90  processed by SPAG 7.51RB at 17:35 on  4 Mar 2022
!> \brief \b zerrqp
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE zerrqp( PATH, NUNIT )
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
!> zerrqp tests the error exits for ZGEQP3.
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
      SUBROUTINE zerrqp(Path,Nunit)
      IMPLICIT NONE
!*--ZERRQP58
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
      PARAMETER (NMAX=3)
!     ..
!     .. Local Scalars ..
      CHARACTER*2 c2
      INTEGER info , lw
!     ..
!     .. Local Arrays ..
      INTEGER ip(NMAX)
      DOUBLE PRECISION rw(2*NMAX)
      COMPLEX*16 a(NMAX,NMAX) , tau(NMAX) , w(2*NMAX+3*NMAX)
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN
EXTERNAL    LSAMEN
!     ..
!     .. External Subroutines ..
EXTERNAL    ZGEQP3
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
      lw = NMAX + 1
      a(1,1) = DCMPLX(1.0D+0,-1.0D+0)
      a(1,2) = DCMPLX(2.0D+0,-2.0D+0)
      a(2,2) = DCMPLX(3.0D+0,-3.0D+0)
      a(2,1) = DCMPLX(4.0D+0,-4.0D+0)
      OK = .TRUE.
      WRITE (NOUt,FMT=*)
!
!     Test error exits for QR factorization with pivoting
!
      IF ( LSAMEN(2,c2,'QP') ) THEN
!
!        ZGEQP3
!
         SRNamt = 'ZGEQP3'
         INFot = 1
         CALL ZGEQP3(-1,0,a,1,ip,tau,w,lw,rw,info)
         CALL CHKXER('ZGEQP3',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGEQP3(1,-1,a,1,ip,tau,w,lw,rw,info)
         CALL CHKXER('ZGEQP3',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZGEQP3(2,3,a,1,ip,tau,w,lw,rw,info)
         CALL CHKXER('ZGEQP3',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZGEQP3(2,2,a,2,ip,tau,w,lw-10,rw,info)
         CALL CHKXER('ZGEQP3',INFot,NOUt,LERr,OK)
      ENDIF
!
!     Print a summary line.
!
      CALL alaesm(Path,OK,NOUt)
!
!
!     End of zerrqp
!
      END SUBROUTINE zerrqp
