!*==serrlq.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b serrlq
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE serrlq( PATH, NUNIT )
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
!> serrlq tests the error exits for the REAL routines
!> that use the LQ decomposition of a general matrix.
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
      SUBROUTINE serrlq(Path,Nunit)
      IMPLICIT NONE
!*--SERRLQ59
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
EXTERNAL    SGELQ2,   SGELQF,   SORGL2
EXTERNAL    SORGLQ,   SORML2,   SORMLQ
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
!     Error exits for LQ factorization
!
!     SGELQF
!
      SRNamt = 'SGELQF'
      INFot = 1
      CALL SGELQF(-1,0,a,1,b,w,1,info)
      CALL CHKXER('SGELQF',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL SGELQF(0,-1,a,1,b,w,1,info)
      CALL CHKXER('SGELQF',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL SGELQF(2,1,a,1,b,w,2,info)
      CALL CHKXER('SGELQF',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL SGELQF(2,1,a,2,b,w,1,info)
      CALL CHKXER('SGELQF',INFot,NOUt,LERr,OK)
!
!     SGELQ2
!
      SRNamt = 'SGELQ2'
      INFot = 1
      CALL SGELQ2(-1,0,a,1,b,w,info)
      CALL CHKXER('SGELQ2',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL SGELQ2(0,-1,a,1,b,w,info)
      CALL CHKXER('SGELQ2',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL SGELQ2(2,1,a,1,b,w,info)
      CALL CHKXER('SGELQ2',INFot,NOUt,LERr,OK)
!
!     sgelqs
!
      SRNamt = 'sgelqs'
      INFot = 1
      CALL sgelqs(-1,0,0,a,1,x,b,1,w,1,info)
      CALL CHKXER('sgelqs',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL sgelqs(0,-1,0,a,1,x,b,1,w,1,info)
      CALL CHKXER('sgelqs',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL sgelqs(2,1,0,a,2,x,b,1,w,1,info)
      CALL CHKXER('sgelqs',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL sgelqs(0,0,-1,a,1,x,b,1,w,1,info)
      CALL CHKXER('sgelqs',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL sgelqs(2,2,0,a,1,x,b,2,w,1,info)
      CALL CHKXER('sgelqs',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL sgelqs(1,2,0,a,1,x,b,1,w,1,info)
      CALL CHKXER('sgelqs',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL sgelqs(1,1,2,a,1,x,b,1,w,1,info)
      CALL CHKXER('sgelqs',INFot,NOUt,LERr,OK)
!
!     SORGLQ
!
      SRNamt = 'SORGLQ'
      INFot = 1
      CALL SORGLQ(-1,0,0,a,1,x,w,1,info)
      CALL CHKXER('SORGLQ',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL SORGLQ(0,-1,0,a,1,x,w,1,info)
      CALL CHKXER('SORGLQ',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL SORGLQ(2,1,0,a,2,x,w,2,info)
      CALL CHKXER('SORGLQ',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL SORGLQ(0,0,-1,a,1,x,w,1,info)
      CALL CHKXER('SORGLQ',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL SORGLQ(1,1,2,a,1,x,w,1,info)
      CALL CHKXER('SORGLQ',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL SORGLQ(2,2,0,a,1,x,w,2,info)
      CALL CHKXER('SORGLQ',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL SORGLQ(2,2,0,a,2,x,w,1,info)
      CALL CHKXER('SORGLQ',INFot,NOUt,LERr,OK)
!
!     SORGL2
!
      SRNamt = 'SORGL2'
      INFot = 1
      CALL SORGL2(-1,0,0,a,1,x,w,info)
      CALL CHKXER('SORGL2',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL SORGL2(0,-1,0,a,1,x,w,info)
      CALL CHKXER('SORGL2',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL SORGL2(2,1,0,a,2,x,w,info)
      CALL CHKXER('SORGL2',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL SORGL2(0,0,-1,a,1,x,w,info)
      CALL CHKXER('SORGL2',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL SORGL2(1,1,2,a,1,x,w,info)
      CALL CHKXER('SORGL2',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL SORGL2(2,2,0,a,1,x,w,info)
      CALL CHKXER('SORGL2',INFot,NOUt,LERr,OK)
!
!     SORMLQ
!
      SRNamt = 'SORMLQ'
      INFot = 1
      CALL SORMLQ('/','N',0,0,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('SORMLQ',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL SORMLQ('L','/',0,0,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('SORMLQ',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL SORMLQ('L','N',-1,0,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('SORMLQ',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL SORMLQ('L','N',0,-1,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('SORMLQ',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL SORMLQ('L','N',0,0,-1,a,1,x,af,1,w,1,info)
      CALL CHKXER('SORMLQ',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL SORMLQ('L','N',0,1,1,a,1,x,af,1,w,1,info)
      CALL CHKXER('SORMLQ',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL SORMLQ('R','N',1,0,1,a,1,x,af,1,w,1,info)
      CALL CHKXER('SORMLQ',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL SORMLQ('L','N',2,0,2,a,1,x,af,2,w,1,info)
      CALL CHKXER('SORMLQ',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL SORMLQ('R','N',0,2,2,a,1,x,af,1,w,1,info)
      CALL CHKXER('SORMLQ',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL SORMLQ('L','N',2,1,0,a,2,x,af,1,w,1,info)
      CALL CHKXER('SORMLQ',INFot,NOUt,LERr,OK)
      INFot = 12
      CALL SORMLQ('L','N',1,2,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('SORMLQ',INFot,NOUt,LERr,OK)
      INFot = 12
      CALL SORMLQ('R','N',2,1,0,a,1,x,af,2,w,1,info)
      CALL CHKXER('SORMLQ',INFot,NOUt,LERr,OK)
!
!     SORML2
!
      SRNamt = 'SORML2'
      INFot = 1
      CALL SORML2('/','N',0,0,0,a,1,x,af,1,w,info)
      CALL CHKXER('SORML2',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL SORML2('L','/',0,0,0,a,1,x,af,1,w,info)
      CALL CHKXER('SORML2',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL SORML2('L','N',-1,0,0,a,1,x,af,1,w,info)
      CALL CHKXER('SORML2',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL SORML2('L','N',0,-1,0,a,1,x,af,1,w,info)
      CALL CHKXER('SORML2',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL SORML2('L','N',0,0,-1,a,1,x,af,1,w,info)
      CALL CHKXER('SORML2',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL SORML2('L','N',0,1,1,a,1,x,af,1,w,info)
      CALL CHKXER('SORML2',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL SORML2('R','N',1,0,1,a,1,x,af,1,w,info)
      CALL CHKXER('SORML2',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL SORML2('L','N',2,1,2,a,1,x,af,2,w,info)
      CALL CHKXER('SORML2',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL SORML2('R','N',1,2,2,a,1,x,af,1,w,info)
      CALL CHKXER('SORML2',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL SORML2('L','N',2,1,0,a,2,x,af,1,w,info)
      CALL CHKXER('SORML2',INFot,NOUt,LERr,OK)
!
!     Print a summary line.
!
      CALL alaesm(Path,OK,NOUt)
!
!
!     End of serrlq
!
      END SUBROUTINE serrlq
