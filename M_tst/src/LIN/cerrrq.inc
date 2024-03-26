!*==cerrrq.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b cerrrq
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE cerrrq( PATH, NUNIT )
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
!> cerrrq tests the error exits for the COMPLEX routines
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
!> \ingroup complex_lin
!
!  =====================================================================
      SUBROUTINE cerrrq(Path,Nunit)
      IMPLICIT NONE
!*--CERRRQ59
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
EXTERNAL    CGERQ2,   CGERQF,   CUNGR2
EXTERNAL    CUNGRQ,   CUNMR2,   CUNMRQ
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
!     Error exits for RQ factorization
!
!     CGERQF
!
      SRNamt = 'CGERQF'
      INFot = 1
      CALL CGERQF(-1,0,a,1,b,w,1,info)
      CALL CHKXER('CGERQF',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CGERQF(0,-1,a,1,b,w,1,info)
      CALL CHKXER('CGERQF',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL CGERQF(2,1,a,1,b,w,2,info)
      CALL CHKXER('CGERQF',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL CGERQF(2,1,a,2,b,w,1,info)
      CALL CHKXER('CGERQF',INFot,NOUt,LERr,OK)
!
!     CGERQ2
!
      SRNamt = 'CGERQ2'
      INFot = 1
      CALL CGERQ2(-1,0,a,1,b,w,info)
      CALL CHKXER('CGERQ2',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CGERQ2(0,-1,a,1,b,w,info)
      CALL CHKXER('CGERQ2',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL CGERQ2(2,1,a,1,b,w,info)
      CALL CHKXER('CGERQ2',INFot,NOUt,LERr,OK)
!
!     cgerqs
!
      SRNamt = 'cgerqs'
      INFot = 1
      CALL cgerqs(-1,0,0,a,1,x,b,1,w,1,info)
      CALL CHKXER('cgerqs',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL cgerqs(0,-1,0,a,1,x,b,1,w,1,info)
      CALL CHKXER('cgerqs',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL cgerqs(2,1,0,a,2,x,b,1,w,1,info)
      CALL CHKXER('cgerqs',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL cgerqs(0,0,-1,a,1,x,b,1,w,1,info)
      CALL CHKXER('cgerqs',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL cgerqs(2,2,0,a,1,x,b,2,w,1,info)
      CALL CHKXER('cgerqs',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL cgerqs(2,2,0,a,2,x,b,1,w,1,info)
      CALL CHKXER('cgerqs',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL cgerqs(1,1,2,a,1,x,b,1,w,1,info)
      CALL CHKXER('cgerqs',INFot,NOUt,LERr,OK)
!
!     CUNGRQ
!
      SRNamt = 'CUNGRQ'
      INFot = 1
      CALL CUNGRQ(-1,0,0,a,1,x,w,1,info)
      CALL CHKXER('CUNGRQ',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CUNGRQ(0,-1,0,a,1,x,w,1,info)
      CALL CHKXER('CUNGRQ',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CUNGRQ(2,1,0,a,2,x,w,2,info)
      CALL CHKXER('CUNGRQ',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL CUNGRQ(0,0,-1,a,1,x,w,1,info)
      CALL CHKXER('CUNGRQ',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL CUNGRQ(1,2,2,a,1,x,w,1,info)
      CALL CHKXER('CUNGRQ',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL CUNGRQ(2,2,0,a,1,x,w,2,info)
      CALL CHKXER('CUNGRQ',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL CUNGRQ(2,2,0,a,2,x,w,1,info)
      CALL CHKXER('CUNGRQ',INFot,NOUt,LERr,OK)
!
!     CUNGR2
!
      SRNamt = 'CUNGR2'
      INFot = 1
      CALL CUNGR2(-1,0,0,a,1,x,w,info)
      CALL CHKXER('CUNGR2',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CUNGR2(0,-1,0,a,1,x,w,info)
      CALL CHKXER('CUNGR2',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CUNGR2(2,1,0,a,2,x,w,info)
      CALL CHKXER('CUNGR2',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL CUNGR2(0,0,-1,a,1,x,w,info)
      CALL CHKXER('CUNGR2',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL CUNGR2(1,2,2,a,2,x,w,info)
      CALL CHKXER('CUNGR2',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL CUNGR2(2,2,0,a,1,x,w,info)
      CALL CHKXER('CUNGR2',INFot,NOUt,LERr,OK)
!
!     CUNMRQ
!
      SRNamt = 'CUNMRQ'
      INFot = 1
      CALL CUNMRQ('/','N',0,0,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('CUNMRQ',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CUNMRQ('L','/',0,0,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('CUNMRQ',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL CUNMRQ('L','N',-1,0,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('CUNMRQ',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL CUNMRQ('L','N',0,-1,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('CUNMRQ',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL CUNMRQ('L','N',0,0,-1,a,1,x,af,1,w,1,info)
      CALL CHKXER('CUNMRQ',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL CUNMRQ('L','N',0,1,1,a,1,x,af,1,w,1,info)
      CALL CHKXER('CUNMRQ',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL CUNMRQ('R','N',1,0,1,a,1,x,af,1,w,1,info)
      CALL CHKXER('CUNMRQ',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL CUNMRQ('L','N',2,1,2,a,1,x,af,2,w,1,info)
      CALL CHKXER('CUNMRQ',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL CUNMRQ('R','N',1,2,2,a,1,x,af,1,w,1,info)
      CALL CHKXER('CUNMRQ',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL CUNMRQ('L','N',2,1,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('CUNMRQ',INFot,NOUt,LERr,OK)
      INFot = 12
      CALL CUNMRQ('L','N',1,2,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('CUNMRQ',INFot,NOUt,LERr,OK)
      INFot = 12
      CALL CUNMRQ('R','N',2,1,0,a,1,x,af,2,w,1,info)
      CALL CHKXER('CUNMRQ',INFot,NOUt,LERr,OK)
!
!     CUNMR2
!
      SRNamt = 'CUNMR2'
      INFot = 1
      CALL CUNMR2('/','N',0,0,0,a,1,x,af,1,w,info)
      CALL CHKXER('CUNMR2',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CUNMR2('L','/',0,0,0,a,1,x,af,1,w,info)
      CALL CHKXER('CUNMR2',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL CUNMR2('L','N',-1,0,0,a,1,x,af,1,w,info)
      CALL CHKXER('CUNMR2',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL CUNMR2('L','N',0,-1,0,a,1,x,af,1,w,info)
      CALL CHKXER('CUNMR2',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL CUNMR2('L','N',0,0,-1,a,1,x,af,1,w,info)
      CALL CHKXER('CUNMR2',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL CUNMR2('L','N',0,1,1,a,1,x,af,1,w,info)
      CALL CHKXER('CUNMR2',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL CUNMR2('R','N',1,0,1,a,1,x,af,1,w,info)
      CALL CHKXER('CUNMR2',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL CUNMR2('L','N',2,1,2,a,1,x,af,2,w,info)
      CALL CHKXER('CUNMR2',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL CUNMR2('R','N',1,2,2,a,1,x,af,1,w,info)
      CALL CHKXER('CUNMR2',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL CUNMR2('L','N',2,1,0,a,1,x,af,1,w,info)
      CALL CHKXER('CUNMR2',INFot,NOUt,LERr,OK)
!
!     Print a summary line.
!
      CALL alaesm(Path,OK,NOUt)
!
!
!     End of cerrrq
!
      END SUBROUTINE cerrrq
