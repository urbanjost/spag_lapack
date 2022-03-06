!*==zerrql.f90  processed by SPAG 7.51RB at 17:35 on  4 Mar 2022
!> \brief \b zerrql
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE zerrql( PATH, NUNIT )
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
!> zerrql tests the error exits for the COMPLEX*16 routines
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
!> \ingroup complex16_lin
!
!  =====================================================================
      SUBROUTINE zerrql(Path,Nunit)
      IMPLICIT NONE
!*--ZERRQL59
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
      COMPLEX*16 a(NMAX,NMAX) , af(NMAX,NMAX) , b(NMAX) , w(NMAX) ,     &
     &           x(NMAX)
!     ..
!     .. External Subroutines ..
EXTERNAL    CHKXER,   ZGEQL2,   ZGEQLF,   ZUNG2L
EXTERNAL    ZUNGQL,   ZUNM2L,   ZUNMQL
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
      INTRINSIC DBLE , DCMPLX
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
            a(i,j) = DCMPLX(1.D0/DBLE(i+j),-1.D0/DBLE(i+j))
            af(i,j) = DCMPLX(1.D0/DBLE(i+j),-1.D0/DBLE(i+j))
         ENDDO
         b(j) = 0.D0
         w(j) = 0.D0
         x(j) = 0.D0
      ENDDO
      OK = .TRUE.
!
!     Error exits for QL factorization
!
!     ZGEQLF
!
      SRNamt = 'ZGEQLF'
      INFot = 1
      CALL ZGEQLF(-1,0,a,1,b,w,1,info)
      CALL CHKXER('ZGEQLF',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZGEQLF(0,-1,a,1,b,w,1,info)
      CALL CHKXER('ZGEQLF',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL ZGEQLF(2,1,a,1,b,w,1,info)
      CALL CHKXER('ZGEQLF',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL ZGEQLF(1,2,a,1,b,w,1,info)
      CALL CHKXER('ZGEQLF',INFot,NOUt,LERr,OK)
!
!     ZGEQL2
!
      SRNamt = 'ZGEQL2'
      INFot = 1
      CALL ZGEQL2(-1,0,a,1,b,w,info)
      CALL CHKXER('ZGEQL2',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZGEQL2(0,-1,a,1,b,w,info)
      CALL CHKXER('ZGEQL2',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL ZGEQL2(2,1,a,1,b,w,info)
      CALL CHKXER('ZGEQL2',INFot,NOUt,LERr,OK)
!
!     zgeqls
!
      SRNamt = 'zgeqls'
      INFot = 1
      CALL zgeqls(-1,0,0,a,1,x,b,1,w,1,info)
      CALL CHKXER('zgeqls',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL zgeqls(0,-1,0,a,1,x,b,1,w,1,info)
      CALL CHKXER('zgeqls',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL zgeqls(1,2,0,a,1,x,b,1,w,1,info)
      CALL CHKXER('zgeqls',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL zgeqls(0,0,-1,a,1,x,b,1,w,1,info)
      CALL CHKXER('zgeqls',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL zgeqls(2,1,0,a,1,x,b,2,w,1,info)
      CALL CHKXER('zgeqls',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL zgeqls(2,1,0,a,2,x,b,1,w,1,info)
      CALL CHKXER('zgeqls',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL zgeqls(1,1,2,a,1,x,b,1,w,1,info)
      CALL CHKXER('zgeqls',INFot,NOUt,LERr,OK)
!
!     ZUNGQL
!
      SRNamt = 'ZUNGQL'
      INFot = 1
      CALL ZUNGQL(-1,0,0,a,1,x,w,1,info)
      CALL CHKXER('ZUNGQL',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZUNGQL(0,-1,0,a,1,x,w,1,info)
      CALL CHKXER('ZUNGQL',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZUNGQL(1,2,0,a,1,x,w,2,info)
      CALL CHKXER('ZUNGQL',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL ZUNGQL(0,0,-1,a,1,x,w,1,info)
      CALL CHKXER('ZUNGQL',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL ZUNGQL(1,1,2,a,1,x,w,1,info)
      CALL CHKXER('ZUNGQL',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL ZUNGQL(2,1,0,a,1,x,w,1,info)
      CALL CHKXER('ZUNGQL',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL ZUNGQL(2,2,0,a,2,x,w,1,info)
      CALL CHKXER('ZUNGQL',INFot,NOUt,LERr,OK)
!
!     ZUNG2L
!
      SRNamt = 'ZUNG2L'
      INFot = 1
      CALL ZUNG2L(-1,0,0,a,1,x,w,info)
      CALL CHKXER('ZUNG2L',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZUNG2L(0,-1,0,a,1,x,w,info)
      CALL CHKXER('ZUNG2L',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZUNG2L(1,2,0,a,1,x,w,info)
      CALL CHKXER('ZUNG2L',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL ZUNG2L(0,0,-1,a,1,x,w,info)
      CALL CHKXER('ZUNG2L',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL ZUNG2L(2,1,2,a,2,x,w,info)
      CALL CHKXER('ZUNG2L',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL ZUNG2L(2,1,0,a,1,x,w,info)
      CALL CHKXER('ZUNG2L',INFot,NOUt,LERr,OK)
!
!     ZUNMQL
!
      SRNamt = 'ZUNMQL'
      INFot = 1
      CALL ZUNMQL('/','N',0,0,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('ZUNMQL',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZUNMQL('L','/',0,0,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('ZUNMQL',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL ZUNMQL('L','N',-1,0,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('ZUNMQL',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL ZUNMQL('L','N',0,-1,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('ZUNMQL',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL ZUNMQL('L','N',0,0,-1,a,1,x,af,1,w,1,info)
      CALL CHKXER('ZUNMQL',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL ZUNMQL('L','N',0,1,1,a,1,x,af,1,w,1,info)
      CALL CHKXER('ZUNMQL',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL ZUNMQL('R','N',1,0,1,a,1,x,af,1,w,1,info)
      CALL CHKXER('ZUNMQL',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL ZUNMQL('L','N',2,1,0,a,1,x,af,2,w,1,info)
      CALL CHKXER('ZUNMQL',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL ZUNMQL('R','N',1,2,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('ZUNMQL',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL ZUNMQL('L','N',2,1,0,a,2,x,af,1,w,1,info)
      CALL CHKXER('ZUNMQL',INFot,NOUt,LERr,OK)
      INFot = 12
      CALL ZUNMQL('L','N',1,2,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('ZUNMQL',INFot,NOUt,LERr,OK)
      INFot = 12
      CALL ZUNMQL('R','N',2,1,0,a,1,x,af,2,w,1,info)
      CALL CHKXER('ZUNMQL',INFot,NOUt,LERr,OK)
!
!     ZUNM2L
!
      SRNamt = 'ZUNM2L'
      INFot = 1
      CALL ZUNM2L('/','N',0,0,0,a,1,x,af,1,w,info)
      CALL CHKXER('ZUNM2L',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZUNM2L('L','/',0,0,0,a,1,x,af,1,w,info)
      CALL CHKXER('ZUNM2L',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL ZUNM2L('L','N',-1,0,0,a,1,x,af,1,w,info)
      CALL CHKXER('ZUNM2L',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL ZUNM2L('L','N',0,-1,0,a,1,x,af,1,w,info)
      CALL CHKXER('ZUNM2L',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL ZUNM2L('L','N',0,0,-1,a,1,x,af,1,w,info)
      CALL CHKXER('ZUNM2L',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL ZUNM2L('L','N',0,1,1,a,1,x,af,1,w,info)
      CALL CHKXER('ZUNM2L',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL ZUNM2L('R','N',1,0,1,a,1,x,af,1,w,info)
      CALL CHKXER('ZUNM2L',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL ZUNM2L('L','N',2,1,0,a,1,x,af,2,w,info)
      CALL CHKXER('ZUNM2L',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL ZUNM2L('R','N',1,2,0,a,1,x,af,1,w,info)
      CALL CHKXER('ZUNM2L',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL ZUNM2L('L','N',2,1,0,a,2,x,af,1,w,info)
      CALL CHKXER('ZUNM2L',INFot,NOUt,LERr,OK)
!
!     Print a summary line.
!
      CALL alaesm(Path,OK,NOUt)
!
!
!     End of zerrql
!
      END SUBROUTINE zerrql
