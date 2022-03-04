!*==derrqp.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b DERRQP
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DERRQP( PATH, NUNIT )
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
!> DERRQP tests the error exits for DGEQP3.
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
      SUBROUTINE DERRQP(Path,Nunit)
      IMPLICIT NONE
!*--DERRQP58
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
      DOUBLE PRECISION a(NMAX,NMAX) , tau(NMAX) , w(3*NMAX+1)
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN
      EXTERNAL LSAMEN
!     ..
!     .. External Subroutines ..
      EXTERNAL ALAESM , CHKXER , DGEQP3
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
!     .. Executable Statements ..
!
      NOUt = Nunit
      WRITE (NOUt,FMT=*)
      c2 = Path(2:3)
      lw = 3*NMAX + 1
      a(1,1) = 1.0D+0
      a(1,2) = 2.0D+0
      a(2,2) = 3.0D+0
      a(2,1) = 4.0D+0
      OK = .TRUE.
!
      IF ( LSAMEN(2,c2,'QP') ) THEN
!
!        Test error exits for QR factorization with pivoting
!
!        DGEQP3
!
         SRNamt = 'DGEQP3'
         INFot = 1
         CALL DGEQP3(-1,0,a,1,ip,tau,w,lw,info)
         CALL CHKXER('DGEQP3',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGEQP3(1,-1,a,1,ip,tau,w,lw,info)
         CALL CHKXER('DGEQP3',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DGEQP3(2,3,a,1,ip,tau,w,lw,info)
         CALL CHKXER('DGEQP3',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DGEQP3(2,2,a,2,ip,tau,w,lw-10,info)
         CALL CHKXER('DGEQP3',INFot,NOUt,LERr,OK)
      ENDIF
!
!     Print a summary line.
!
      CALL ALAESM(Path,OK,NOUt)
!
!
!     End of DERRQP
!
      END SUBROUTINE DERRQP