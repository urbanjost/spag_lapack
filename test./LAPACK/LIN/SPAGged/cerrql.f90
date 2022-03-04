!*==cerrql.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b CERRQL
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CERRQL( PATH, NUNIT )
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
!> CERRQL tests the error exits for the COMPLEX routines
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
!> \ingroup complex_lin
!
!  =====================================================================
      SUBROUTINE CERRQL(Path,Nunit)
      IMPLICIT NONE
!*--CERRQL59
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
      COMPLEX a(NMAX,NMAX) , af(NMAX,NMAX) , b(NMAX) , w(NMAX) , x(NMAX)
!     ..
!     .. External Subroutines ..
      EXTERNAL ALAESM , CGEQL2 , CGEQLF , CGEQLS , CHKXER , CUNG2L ,    &
     &         CUNGQL , CUNM2L , CUNMQL
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
      INTRINSIC CMPLX , REAL
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
            a(i,j) = CMPLX(1./REAL(i+j),-1./REAL(i+j))
            af(i,j) = CMPLX(1./REAL(i+j),-1./REAL(i+j))
         ENDDO
         b(j) = 0.
         w(j) = 0.
         x(j) = 0.
      ENDDO
      OK = .TRUE.
!
!     Error exits for QL factorization
!
!     CGEQLF
!
      SRNamt = 'CGEQLF'
      INFot = 1
      CALL CGEQLF(-1,0,a,1,b,w,1,info)
      CALL CHKXER('CGEQLF',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CGEQLF(0,-1,a,1,b,w,1,info)
      CALL CHKXER('CGEQLF',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL CGEQLF(2,1,a,1,b,w,1,info)
      CALL CHKXER('CGEQLF',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL CGEQLF(1,2,a,1,b,w,1,info)
      CALL CHKXER('CGEQLF',INFot,NOUt,LERr,OK)
!
!     CGEQL2
!
      SRNamt = 'CGEQL2'
      INFot = 1
      CALL CGEQL2(-1,0,a,1,b,w,info)
      CALL CHKXER('CGEQL2',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CGEQL2(0,-1,a,1,b,w,info)
      CALL CHKXER('CGEQL2',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL CGEQL2(2,1,a,1,b,w,info)
      CALL CHKXER('CGEQL2',INFot,NOUt,LERr,OK)
!
!     CGEQLS
!
      SRNamt = 'CGEQLS'
      INFot = 1
      CALL CGEQLS(-1,0,0,a,1,x,b,1,w,1,info)
      CALL CHKXER('CGEQLS',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CGEQLS(0,-1,0,a,1,x,b,1,w,1,info)
      CALL CHKXER('CGEQLS',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CGEQLS(1,2,0,a,1,x,b,1,w,1,info)
      CALL CHKXER('CGEQLS',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL CGEQLS(0,0,-1,a,1,x,b,1,w,1,info)
      CALL CHKXER('CGEQLS',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL CGEQLS(2,1,0,a,1,x,b,2,w,1,info)
      CALL CHKXER('CGEQLS',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL CGEQLS(2,1,0,a,2,x,b,1,w,1,info)
      CALL CHKXER('CGEQLS',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL CGEQLS(1,1,2,a,1,x,b,1,w,1,info)
      CALL CHKXER('CGEQLS',INFot,NOUt,LERr,OK)
!
!     CUNGQL
!
      SRNamt = 'CUNGQL'
      INFot = 1
      CALL CUNGQL(-1,0,0,a,1,x,w,1,info)
      CALL CHKXER('CUNGQL',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CUNGQL(0,-1,0,a,1,x,w,1,info)
      CALL CHKXER('CUNGQL',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CUNGQL(1,2,0,a,1,x,w,2,info)
      CALL CHKXER('CUNGQL',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL CUNGQL(0,0,-1,a,1,x,w,1,info)
      CALL CHKXER('CUNGQL',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL CUNGQL(1,1,2,a,1,x,w,1,info)
      CALL CHKXER('CUNGQL',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL CUNGQL(2,1,0,a,1,x,w,1,info)
      CALL CHKXER('CUNGQL',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL CUNGQL(2,2,0,a,2,x,w,1,info)
      CALL CHKXER('CUNGQL',INFot,NOUt,LERr,OK)
!
!     CUNG2L
!
      SRNamt = 'CUNG2L'
      INFot = 1
      CALL CUNG2L(-1,0,0,a,1,x,w,info)
      CALL CHKXER('CUNG2L',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CUNG2L(0,-1,0,a,1,x,w,info)
      CALL CHKXER('CUNG2L',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CUNG2L(1,2,0,a,1,x,w,info)
      CALL CHKXER('CUNG2L',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL CUNG2L(0,0,-1,a,1,x,w,info)
      CALL CHKXER('CUNG2L',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL CUNG2L(2,1,2,a,2,x,w,info)
      CALL CHKXER('CUNG2L',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL CUNG2L(2,1,0,a,1,x,w,info)
      CALL CHKXER('CUNG2L',INFot,NOUt,LERr,OK)
!
!     CUNMQL
!
      SRNamt = 'CUNMQL'
      INFot = 1
      CALL CUNMQL('/','N',0,0,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('CUNMQL',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CUNMQL('L','/',0,0,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('CUNMQL',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL CUNMQL('L','N',-1,0,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('CUNMQL',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL CUNMQL('L','N',0,-1,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('CUNMQL',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL CUNMQL('L','N',0,0,-1,a,1,x,af,1,w,1,info)
      CALL CHKXER('CUNMQL',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL CUNMQL('L','N',0,1,1,a,1,x,af,1,w,1,info)
      CALL CHKXER('CUNMQL',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL CUNMQL('R','N',1,0,1,a,1,x,af,1,w,1,info)
      CALL CHKXER('CUNMQL',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL CUNMQL('L','N',2,1,0,a,1,x,af,2,w,1,info)
      CALL CHKXER('CUNMQL',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL CUNMQL('R','N',1,2,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('CUNMQL',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL CUNMQL('L','N',2,1,0,a,2,x,af,1,w,1,info)
      CALL CHKXER('CUNMQL',INFot,NOUt,LERr,OK)
      INFot = 12
      CALL CUNMQL('L','N',1,2,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('CUNMQL',INFot,NOUt,LERr,OK)
      INFot = 12
      CALL CUNMQL('R','N',2,1,0,a,1,x,af,2,w,1,info)
      CALL CHKXER('CUNMQL',INFot,NOUt,LERr,OK)
!
!     CUNM2L
!
      SRNamt = 'CUNM2L'
      INFot = 1
      CALL CUNM2L('/','N',0,0,0,a,1,x,af,1,w,info)
      CALL CHKXER('CUNM2L',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CUNM2L('L','/',0,0,0,a,1,x,af,1,w,info)
      CALL CHKXER('CUNM2L',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL CUNM2L('L','N',-1,0,0,a,1,x,af,1,w,info)
      CALL CHKXER('CUNM2L',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL CUNM2L('L','N',0,-1,0,a,1,x,af,1,w,info)
      CALL CHKXER('CUNM2L',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL CUNM2L('L','N',0,0,-1,a,1,x,af,1,w,info)
      CALL CHKXER('CUNM2L',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL CUNM2L('L','N',0,1,1,a,1,x,af,1,w,info)
      CALL CHKXER('CUNM2L',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL CUNM2L('R','N',1,0,1,a,1,x,af,1,w,info)
      CALL CHKXER('CUNM2L',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL CUNM2L('L','N',2,1,0,a,1,x,af,2,w,info)
      CALL CHKXER('CUNM2L',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL CUNM2L('R','N',1,2,0,a,1,x,af,1,w,info)
      CALL CHKXER('CUNM2L',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL CUNM2L('L','N',2,1,0,a,2,x,af,1,w,info)
      CALL CHKXER('CUNM2L',INFot,NOUt,LERr,OK)
!
!     Print a summary line.
!
      CALL ALAESM(Path,OK,NOUt)
!
!
!     End of CERRQL
!
      END SUBROUTINE CERRQL
