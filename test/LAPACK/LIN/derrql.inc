!*==derrql.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b derrql
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE derrql( PATH, NUNIT )
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
!> derrql tests the error exits for the DOUBLE PRECISION routines
!> that use the QL decomposition of a general matrix.
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
      SUBROUTINE derrql(Path,Nunit)
      IMPLICIT NONE
!*--DERRQL59
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
      DOUBLE PRECISION a(NMAX,NMAX) , af(NMAX,NMAX) , b(NMAX) , w(NMAX) &
     &                 , x(NMAX)
!     ..
!     .. External Subroutines ..
EXTERNAL    CHKXER,   DGEQL2,   DGEQLF,   DORG2L
EXTERNAL    DORGQL,   DORM2L,   DORMQL
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
         w(j) = 0.D0
         x(j) = 0.D0
      ENDDO
      OK = .TRUE.
!
!     Error exits for QL factorization
!
!     DGEQLF
!
      SRNamt = 'DGEQLF'
      INFot = 1
      CALL DGEQLF(-1,0,a,1,b,w,1,info)
      CALL CHKXER('DGEQLF',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DGEQLF(0,-1,a,1,b,w,1,info)
      CALL CHKXER('DGEQLF',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL DGEQLF(2,1,a,1,b,w,1,info)
      CALL CHKXER('DGEQLF',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL DGEQLF(1,2,a,1,b,w,1,info)
      CALL CHKXER('DGEQLF',INFot,NOUt,LERr,OK)
!
!     DGEQL2
!
      SRNamt = 'DGEQL2'
      INFot = 1
      CALL DGEQL2(-1,0,a,1,b,w,info)
      CALL CHKXER('DGEQL2',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DGEQL2(0,-1,a,1,b,w,info)
      CALL CHKXER('DGEQL2',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL DGEQL2(2,1,a,1,b,w,info)
      CALL CHKXER('DGEQL2',INFot,NOUt,LERr,OK)
!
!     dgeqls
!
      SRNamt = 'dgeqls'
      INFot = 1
      CALL dgeqls(-1,0,0,a,1,x,b,1,w,1,info)
      CALL CHKXER('dgeqls',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL dgeqls(0,-1,0,a,1,x,b,1,w,1,info)
      CALL CHKXER('dgeqls',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL dgeqls(1,2,0,a,1,x,b,1,w,1,info)
      CALL CHKXER('dgeqls',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL dgeqls(0,0,-1,a,1,x,b,1,w,1,info)
      CALL CHKXER('dgeqls',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL dgeqls(2,1,0,a,1,x,b,2,w,1,info)
      CALL CHKXER('dgeqls',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL dgeqls(2,1,0,a,2,x,b,1,w,1,info)
      CALL CHKXER('dgeqls',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL dgeqls(1,1,2,a,1,x,b,1,w,1,info)
      CALL CHKXER('dgeqls',INFot,NOUt,LERr,OK)
!
!     DORGQL
!
      SRNamt = 'DORGQL'
      INFot = 1
      CALL DORGQL(-1,0,0,a,1,x,w,1,info)
      CALL CHKXER('DORGQL',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DORGQL(0,-1,0,a,1,x,w,1,info)
      CALL CHKXER('DORGQL',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DORGQL(1,2,0,a,1,x,w,2,info)
      CALL CHKXER('DORGQL',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL DORGQL(0,0,-1,a,1,x,w,1,info)
      CALL CHKXER('DORGQL',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL DORGQL(1,1,2,a,1,x,w,1,info)
      CALL CHKXER('DORGQL',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL DORGQL(2,1,0,a,1,x,w,1,info)
      CALL CHKXER('DORGQL',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL DORGQL(2,2,0,a,2,x,w,1,info)
      CALL CHKXER('DORGQL',INFot,NOUt,LERr,OK)
!
!     DORG2L
!
      SRNamt = 'DORG2L'
      INFot = 1
      CALL DORG2L(-1,0,0,a,1,x,w,info)
      CALL CHKXER('DORG2L',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DORG2L(0,-1,0,a,1,x,w,info)
      CALL CHKXER('DORG2L',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DORG2L(1,2,0,a,1,x,w,info)
      CALL CHKXER('DORG2L',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL DORG2L(0,0,-1,a,1,x,w,info)
      CALL CHKXER('DORG2L',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL DORG2L(2,1,2,a,2,x,w,info)
      CALL CHKXER('DORG2L',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL DORG2L(2,1,0,a,1,x,w,info)
      CALL CHKXER('DORG2L',INFot,NOUt,LERr,OK)
!
!     DORMQL
!
      SRNamt = 'DORMQL'
      INFot = 1
      CALL DORMQL('/','N',0,0,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('DORMQL',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DORMQL('L','/',0,0,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('DORMQL',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL DORMQL('L','N',-1,0,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('DORMQL',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL DORMQL('L','N',0,-1,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('DORMQL',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL DORMQL('L','N',0,0,-1,a,1,x,af,1,w,1,info)
      CALL CHKXER('DORMQL',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL DORMQL('L','N',0,1,1,a,1,x,af,1,w,1,info)
      CALL CHKXER('DORMQL',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL DORMQL('R','N',1,0,1,a,1,x,af,1,w,1,info)
      CALL CHKXER('DORMQL',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL DORMQL('L','N',2,1,0,a,1,x,af,2,w,1,info)
      CALL CHKXER('DORMQL',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL DORMQL('R','N',1,2,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('DORMQL',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL DORMQL('L','N',2,1,0,a,2,x,af,1,w,1,info)
      CALL CHKXER('DORMQL',INFot,NOUt,LERr,OK)
      INFot = 12
      CALL DORMQL('L','N',1,2,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('DORMQL',INFot,NOUt,LERr,OK)
      INFot = 12
      CALL DORMQL('R','N',2,1,0,a,1,x,af,2,w,1,info)
      CALL CHKXER('DORMQL',INFot,NOUt,LERr,OK)
!
!     DORM2L
!
      SRNamt = 'DORM2L'
      INFot = 1
      CALL DORM2L('/','N',0,0,0,a,1,x,af,1,w,info)
      CALL CHKXER('DORM2L',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DORM2L('L','/',0,0,0,a,1,x,af,1,w,info)
      CALL CHKXER('DORM2L',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL DORM2L('L','N',-1,0,0,a,1,x,af,1,w,info)
      CALL CHKXER('DORM2L',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL DORM2L('L','N',0,-1,0,a,1,x,af,1,w,info)
      CALL CHKXER('DORM2L',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL DORM2L('L','N',0,0,-1,a,1,x,af,1,w,info)
      CALL CHKXER('DORM2L',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL DORM2L('L','N',0,1,1,a,1,x,af,1,w,info)
      CALL CHKXER('DORM2L',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL DORM2L('R','N',1,0,1,a,1,x,af,1,w,info)
      CALL CHKXER('DORM2L',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL DORM2L('L','N',2,1,0,a,1,x,af,2,w,info)
      CALL CHKXER('DORM2L',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL DORM2L('R','N',1,2,0,a,1,x,af,1,w,info)
      CALL CHKXER('DORM2L',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL DORM2L('L','N',2,1,0,a,2,x,af,1,w,info)
      CALL CHKXER('DORM2L',INFot,NOUt,LERr,OK)
!
!     Print a summary line.
!
      CALL alaesm(Path,OK,NOUt)
!
!
!     End of derrql
!
      END SUBROUTINE derrql
