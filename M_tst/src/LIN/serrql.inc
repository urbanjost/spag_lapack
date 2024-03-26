!*==serrql.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b serrql
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE serrql( PATH, NUNIT )
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
!> serrql tests the error exits for the REAL routines
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
!> \ingroup single_lin
!
!  =====================================================================
      SUBROUTINE serrql(Path,Nunit)
      IMPLICIT NONE
!*--SERRQL59
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
      REAL a(NMAX,NMAX) , af(NMAX,NMAX) , b(NMAX) , w(NMAX) , x(NMAX)
!     ..
!     .. External Subroutines ..
EXTERNAL    SGEQL2,   SGEQLF,   SORG2L
EXTERNAL    SORGQL,   SORM2L,   SORMQL
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
      INTRINSIC REAL
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
            a(i,j) = 1./REAL(i+j)
            af(i,j) = 1./REAL(i+j)
         ENDDO
         b(j) = 0.
         w(j) = 0.
         x(j) = 0.
      ENDDO
      OK = .TRUE.
!
!     Error exits for QL factorization
!
!     SGEQLF
!
      SRNamt = 'SGEQLF'
      INFot = 1
      CALL SGEQLF(-1,0,a,1,b,w,1,info)
      CALL CHKXER('SGEQLF',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL SGEQLF(0,-1,a,1,b,w,1,info)
      CALL CHKXER('SGEQLF',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL SGEQLF(2,1,a,1,b,w,1,info)
      CALL CHKXER('SGEQLF',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL SGEQLF(1,2,a,1,b,w,1,info)
      CALL CHKXER('SGEQLF',INFot,NOUt,LERr,OK)
!
!     SGEQL2
!
      SRNamt = 'SGEQL2'
      INFot = 1
      CALL SGEQL2(-1,0,a,1,b,w,info)
      CALL CHKXER('SGEQL2',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL SGEQL2(0,-1,a,1,b,w,info)
      CALL CHKXER('SGEQL2',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL SGEQL2(2,1,a,1,b,w,info)
      CALL CHKXER('SGEQL2',INFot,NOUt,LERr,OK)
!
!     sgeqls
!
      SRNamt = 'sgeqls'
      INFot = 1
      CALL sgeqls(-1,0,0,a,1,x,b,1,w,1,info)
      CALL CHKXER('sgeqls',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL sgeqls(0,-1,0,a,1,x,b,1,w,1,info)
      CALL CHKXER('sgeqls',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL sgeqls(1,2,0,a,1,x,b,1,w,1,info)
      CALL CHKXER('sgeqls',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL sgeqls(0,0,-1,a,1,x,b,1,w,1,info)
      CALL CHKXER('sgeqls',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL sgeqls(2,1,0,a,1,x,b,2,w,1,info)
      CALL CHKXER('sgeqls',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL sgeqls(2,1,0,a,2,x,b,1,w,1,info)
      CALL CHKXER('sgeqls',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL sgeqls(1,1,2,a,1,x,b,1,w,1,info)
      CALL CHKXER('sgeqls',INFot,NOUt,LERr,OK)
!
!     SORGQL
!
      SRNamt = 'SORGQL'
      INFot = 1
      CALL SORGQL(-1,0,0,a,1,x,w,1,info)
      CALL CHKXER('SORGQL',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL SORGQL(0,-1,0,a,1,x,w,1,info)
      CALL CHKXER('SORGQL',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL SORGQL(1,2,0,a,1,x,w,2,info)
      CALL CHKXER('SORGQL',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL SORGQL(0,0,-1,a,1,x,w,1,info)
      CALL CHKXER('SORGQL',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL SORGQL(1,1,2,a,1,x,w,1,info)
      CALL CHKXER('SORGQL',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL SORGQL(2,1,0,a,1,x,w,1,info)
      CALL CHKXER('SORGQL',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL SORGQL(2,2,0,a,2,x,w,1,info)
      CALL CHKXER('SORGQL',INFot,NOUt,LERr,OK)
!
!     SORG2L
!
      SRNamt = 'SORG2L'
      INFot = 1
      CALL SORG2L(-1,0,0,a,1,x,w,info)
      CALL CHKXER('SORG2L',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL SORG2L(0,-1,0,a,1,x,w,info)
      CALL CHKXER('SORG2L',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL SORG2L(1,2,0,a,1,x,w,info)
      CALL CHKXER('SORG2L',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL SORG2L(0,0,-1,a,1,x,w,info)
      CALL CHKXER('SORG2L',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL SORG2L(2,1,2,a,2,x,w,info)
      CALL CHKXER('SORG2L',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL SORG2L(2,1,0,a,1,x,w,info)
      CALL CHKXER('SORG2L',INFot,NOUt,LERr,OK)
!
!     SORMQL
!
      SRNamt = 'SORMQL'
      INFot = 1
      CALL SORMQL('/','N',0,0,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('SORMQL',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL SORMQL('L','/',0,0,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('SORMQL',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL SORMQL('L','N',-1,0,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('SORMQL',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL SORMQL('L','N',0,-1,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('SORMQL',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL SORMQL('L','N',0,0,-1,a,1,x,af,1,w,1,info)
      CALL CHKXER('SORMQL',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL SORMQL('L','N',0,1,1,a,1,x,af,1,w,1,info)
      CALL CHKXER('SORMQL',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL SORMQL('R','N',1,0,1,a,1,x,af,1,w,1,info)
      CALL CHKXER('SORMQL',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL SORMQL('L','N',2,1,0,a,1,x,af,2,w,1,info)
      CALL CHKXER('SORMQL',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL SORMQL('R','N',1,2,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('SORMQL',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL SORMQL('L','N',2,1,0,a,2,x,af,1,w,1,info)
      CALL CHKXER('SORMQL',INFot,NOUt,LERr,OK)
      INFot = 12
      CALL SORMQL('L','N',1,2,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('SORMQL',INFot,NOUt,LERr,OK)
      INFot = 12
      CALL SORMQL('R','N',2,1,0,a,1,x,af,2,w,1,info)
      CALL CHKXER('SORMQL',INFot,NOUt,LERr,OK)
!
!     SORM2L
!
      SRNamt = 'SORM2L'
      INFot = 1
      CALL SORM2L('/','N',0,0,0,a,1,x,af,1,w,info)
      CALL CHKXER('SORM2L',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL SORM2L('L','/',0,0,0,a,1,x,af,1,w,info)
      CALL CHKXER('SORM2L',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL SORM2L('L','N',-1,0,0,a,1,x,af,1,w,info)
      CALL CHKXER('SORM2L',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL SORM2L('L','N',0,-1,0,a,1,x,af,1,w,info)
      CALL CHKXER('SORM2L',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL SORM2L('L','N',0,0,-1,a,1,x,af,1,w,info)
      CALL CHKXER('SORM2L',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL SORM2L('L','N',0,1,1,a,1,x,af,1,w,info)
      CALL CHKXER('SORM2L',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL SORM2L('R','N',1,0,1,a,1,x,af,1,w,info)
      CALL CHKXER('SORM2L',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL SORM2L('L','N',2,1,0,a,1,x,af,2,w,info)
      CALL CHKXER('SORM2L',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL SORM2L('R','N',1,2,0,a,1,x,af,1,w,info)
      CALL CHKXER('SORM2L',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL SORM2L('L','N',2,1,0,a,2,x,af,1,w,info)
      CALL CHKXER('SORM2L',INFot,NOUt,LERr,OK)
!
!     Print a summary line.
!
      CALL alaesm(Path,OK,NOUt)
!
!
!     End of serrql
!
      END SUBROUTINE serrql
