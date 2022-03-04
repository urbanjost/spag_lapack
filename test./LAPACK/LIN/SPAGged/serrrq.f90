!*==serrrq.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b SERRRQ
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SERRRQ( PATH, NUNIT )
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
!> SERRRQ tests the error exits for the REAL routines
!> that use the RQ decomposition of a general matrix.
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
      SUBROUTINE SERRRQ(Path,Nunit)
      IMPLICIT NONE
!*--SERRRQ59
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
      EXTERNAL ALAESM , CHKXER , SGERQ2 , SGERQF , SGERQS , SORGR2 ,    &
     &         SORGRQ , SORMR2 , SORMRQ
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
!     Error exits for RQ factorization
!
!     SGERQF
!
      SRNamt = 'SGERQF'
      INFot = 1
      CALL SGERQF(-1,0,a,1,b,w,1,info)
      CALL CHKXER('SGERQF',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL SGERQF(0,-1,a,1,b,w,1,info)
      CALL CHKXER('SGERQF',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL SGERQF(2,1,a,1,b,w,2,info)
      CALL CHKXER('SGERQF',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL SGERQF(2,1,a,2,b,w,1,info)
      CALL CHKXER('SGERQF',INFot,NOUt,LERr,OK)
!
!     SGERQ2
!
      SRNamt = 'SGERQ2'
      INFot = 1
      CALL SGERQ2(-1,0,a,1,b,w,info)
      CALL CHKXER('SGERQ2',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL SGERQ2(0,-1,a,1,b,w,info)
      CALL CHKXER('SGERQ2',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL SGERQ2(2,1,a,1,b,w,info)
      CALL CHKXER('SGERQ2',INFot,NOUt,LERr,OK)
!
!     SGERQS
!
      SRNamt = 'SGERQS'
      INFot = 1
      CALL SGERQS(-1,0,0,a,1,x,b,1,w,1,info)
      CALL CHKXER('SGERQS',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL SGERQS(0,-1,0,a,1,x,b,1,w,1,info)
      CALL CHKXER('SGERQS',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL SGERQS(2,1,0,a,2,x,b,1,w,1,info)
      CALL CHKXER('SGERQS',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL SGERQS(0,0,-1,a,1,x,b,1,w,1,info)
      CALL CHKXER('SGERQS',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL SGERQS(2,2,0,a,1,x,b,2,w,1,info)
      CALL CHKXER('SGERQS',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL SGERQS(2,2,0,a,2,x,b,1,w,1,info)
      CALL CHKXER('SGERQS',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL SGERQS(1,1,2,a,1,x,b,1,w,1,info)
      CALL CHKXER('SGERQS',INFot,NOUt,LERr,OK)
!
!     SORGRQ
!
      SRNamt = 'SORGRQ'
      INFot = 1
      CALL SORGRQ(-1,0,0,a,1,x,w,1,info)
      CALL CHKXER('SORGRQ',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL SORGRQ(0,-1,0,a,1,x,w,1,info)
      CALL CHKXER('SORGRQ',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL SORGRQ(2,1,0,a,2,x,w,2,info)
      CALL CHKXER('SORGRQ',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL SORGRQ(0,0,-1,a,1,x,w,1,info)
      CALL CHKXER('SORGRQ',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL SORGRQ(1,2,2,a,1,x,w,1,info)
      CALL CHKXER('SORGRQ',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL SORGRQ(2,2,0,a,1,x,w,2,info)
      CALL CHKXER('SORGRQ',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL SORGRQ(2,2,0,a,2,x,w,1,info)
      CALL CHKXER('SORGRQ',INFot,NOUt,LERr,OK)
!
!     SORGR2
!
      SRNamt = 'SORGR2'
      INFot = 1
      CALL SORGR2(-1,0,0,a,1,x,w,info)
      CALL CHKXER('SORGR2',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL SORGR2(0,-1,0,a,1,x,w,info)
      CALL CHKXER('SORGR2',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL SORGR2(2,1,0,a,2,x,w,info)
      CALL CHKXER('SORGR2',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL SORGR2(0,0,-1,a,1,x,w,info)
      CALL CHKXER('SORGR2',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL SORGR2(1,2,2,a,2,x,w,info)
      CALL CHKXER('SORGR2',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL SORGR2(2,2,0,a,1,x,w,info)
      CALL CHKXER('SORGR2',INFot,NOUt,LERr,OK)
!
!     SORMRQ
!
      SRNamt = 'SORMRQ'
      INFot = 1
      CALL SORMRQ('/','N',0,0,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('SORMRQ',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL SORMRQ('L','/',0,0,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('SORMRQ',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL SORMRQ('L','N',-1,0,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('SORMRQ',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL SORMRQ('L','N',0,-1,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('SORMRQ',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL SORMRQ('L','N',0,0,-1,a,1,x,af,1,w,1,info)
      CALL CHKXER('SORMRQ',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL SORMRQ('L','N',0,1,1,a,1,x,af,1,w,1,info)
      CALL CHKXER('SORMRQ',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL SORMRQ('R','N',1,0,1,a,1,x,af,1,w,1,info)
      CALL CHKXER('SORMRQ',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL SORMRQ('L','N',2,1,2,a,1,x,af,2,w,1,info)
      CALL CHKXER('SORMRQ',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL SORMRQ('R','N',1,2,2,a,1,x,af,1,w,1,info)
      CALL CHKXER('SORMRQ',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL SORMRQ('L','N',2,1,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('SORMRQ',INFot,NOUt,LERr,OK)
      INFot = 12
      CALL SORMRQ('L','N',1,2,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('SORMRQ',INFot,NOUt,LERr,OK)
      INFot = 12
      CALL SORMRQ('R','N',2,1,0,a,1,x,af,2,w,1,info)
      CALL CHKXER('SORMRQ',INFot,NOUt,LERr,OK)
!
!     SORMR2
!
      SRNamt = 'SORMR2'
      INFot = 1
      CALL SORMR2('/','N',0,0,0,a,1,x,af,1,w,info)
      CALL CHKXER('SORMR2',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL SORMR2('L','/',0,0,0,a,1,x,af,1,w,info)
      CALL CHKXER('SORMR2',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL SORMR2('L','N',-1,0,0,a,1,x,af,1,w,info)
      CALL CHKXER('SORMR2',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL SORMR2('L','N',0,-1,0,a,1,x,af,1,w,info)
      CALL CHKXER('SORMR2',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL SORMR2('L','N',0,0,-1,a,1,x,af,1,w,info)
      CALL CHKXER('SORMR2',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL SORMR2('L','N',0,1,1,a,1,x,af,1,w,info)
      CALL CHKXER('SORMR2',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL SORMR2('R','N',1,0,1,a,1,x,af,1,w,info)
      CALL CHKXER('SORMR2',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL SORMR2('L','N',2,1,2,a,1,x,af,2,w,info)
      CALL CHKXER('SORMR2',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL SORMR2('R','N',1,2,2,a,1,x,af,1,w,info)
      CALL CHKXER('SORMR2',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL SORMR2('L','N',2,1,0,a,1,x,af,1,w,info)
      CALL CHKXER('SORMR2',INFot,NOUt,LERr,OK)
!
!     Print a summary line.
!
      CALL ALAESM(Path,OK,NOUt)
!
!
!     End of SERRRQ
!
      END SUBROUTINE SERRRQ
