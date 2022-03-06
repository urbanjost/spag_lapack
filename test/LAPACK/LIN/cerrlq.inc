!*==cerrlq.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b cerrlq
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE cerrlq( PATH, NUNIT )
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
!> cerrlq tests the error exits for the COMPLEX routines
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
!> \ingroup complex_lin
!
!  =====================================================================
      SUBROUTINE cerrlq(Path,Nunit)
      IMPLICIT NONE
!*--CERRLQ59
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
EXTERNAL    CGELQ2,   CGELQF,   CHKXER,   CUNGL2
EXTERNAL    CUNGLQ,   CUNML2,   CUNMLQ
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
!     Error exits for LQ factorization
!
!     CGELQF
!
      SRNamt = 'CGELQF'
      INFot = 1
      CALL CGELQF(-1,0,a,1,b,w,1,info)
      CALL CHKXER('CGELQF',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CGELQF(0,-1,a,1,b,w,1,info)
      CALL CHKXER('CGELQF',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL CGELQF(2,1,a,1,b,w,2,info)
      CALL CHKXER('CGELQF',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL CGELQF(2,1,a,2,b,w,1,info)
      CALL CHKXER('CGELQF',INFot,NOUt,LERr,OK)
!
!     CGELQ2
!
      SRNamt = 'CGELQ2'
      INFot = 1
      CALL CGELQ2(-1,0,a,1,b,w,info)
      CALL CHKXER('CGELQ2',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CGELQ2(0,-1,a,1,b,w,info)
      CALL CHKXER('CGELQ2',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL CGELQ2(2,1,a,1,b,w,info)
      CALL CHKXER('CGELQ2',INFot,NOUt,LERr,OK)
!
!     cgelqs
!
      SRNamt = 'cgelqs'
      INFot = 1
      CALL cgelqs(-1,0,0,a,1,x,b,1,w,1,info)
      CALL CHKXER('cgelqs',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL cgelqs(0,-1,0,a,1,x,b,1,w,1,info)
      CALL CHKXER('cgelqs',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL cgelqs(2,1,0,a,2,x,b,1,w,1,info)
      CALL CHKXER('cgelqs',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL cgelqs(0,0,-1,a,1,x,b,1,w,1,info)
      CALL CHKXER('cgelqs',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL cgelqs(2,2,0,a,1,x,b,2,w,1,info)
      CALL CHKXER('cgelqs',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL cgelqs(1,2,0,a,1,x,b,1,w,1,info)
      CALL CHKXER('cgelqs',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL cgelqs(1,1,2,a,1,x,b,1,w,1,info)
      CALL CHKXER('cgelqs',INFot,NOUt,LERr,OK)
!
!     CUNGLQ
!
      SRNamt = 'CUNGLQ'
      INFot = 1
      CALL CUNGLQ(-1,0,0,a,1,x,w,1,info)
      CALL CHKXER('CUNGLQ',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CUNGLQ(0,-1,0,a,1,x,w,1,info)
      CALL CHKXER('CUNGLQ',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CUNGLQ(2,1,0,a,2,x,w,2,info)
      CALL CHKXER('CUNGLQ',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL CUNGLQ(0,0,-1,a,1,x,w,1,info)
      CALL CHKXER('CUNGLQ',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL CUNGLQ(1,1,2,a,1,x,w,1,info)
      CALL CHKXER('CUNGLQ',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL CUNGLQ(2,2,0,a,1,x,w,2,info)
      CALL CHKXER('CUNGLQ',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL CUNGLQ(2,2,0,a,2,x,w,1,info)
      CALL CHKXER('CUNGLQ',INFot,NOUt,LERr,OK)
!
!     CUNGL2
!
      SRNamt = 'CUNGL2'
      INFot = 1
      CALL CUNGL2(-1,0,0,a,1,x,w,info)
      CALL CHKXER('CUNGL2',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CUNGL2(0,-1,0,a,1,x,w,info)
      CALL CHKXER('CUNGL2',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CUNGL2(2,1,0,a,2,x,w,info)
      CALL CHKXER('CUNGL2',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL CUNGL2(0,0,-1,a,1,x,w,info)
      CALL CHKXER('CUNGL2',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL CUNGL2(1,1,2,a,1,x,w,info)
      CALL CHKXER('CUNGL2',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL CUNGL2(2,2,0,a,1,x,w,info)
      CALL CHKXER('CUNGL2',INFot,NOUt,LERr,OK)
!
!     CUNMLQ
!
      SRNamt = 'CUNMLQ'
      INFot = 1
      CALL CUNMLQ('/','N',0,0,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('CUNMLQ',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CUNMLQ('L','/',0,0,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('CUNMLQ',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL CUNMLQ('L','N',-1,0,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('CUNMLQ',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL CUNMLQ('L','N',0,-1,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('CUNMLQ',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL CUNMLQ('L','N',0,0,-1,a,1,x,af,1,w,1,info)
      CALL CHKXER('CUNMLQ',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL CUNMLQ('L','N',0,1,1,a,1,x,af,1,w,1,info)
      CALL CHKXER('CUNMLQ',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL CUNMLQ('R','N',1,0,1,a,1,x,af,1,w,1,info)
      CALL CHKXER('CUNMLQ',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL CUNMLQ('L','N',2,0,2,a,1,x,af,2,w,1,info)
      CALL CHKXER('CUNMLQ',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL CUNMLQ('R','N',0,2,2,a,1,x,af,1,w,1,info)
      CALL CHKXER('CUNMLQ',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL CUNMLQ('L','N',2,1,0,a,2,x,af,1,w,1,info)
      CALL CHKXER('CUNMLQ',INFot,NOUt,LERr,OK)
      INFot = 12
      CALL CUNMLQ('L','N',1,2,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('CUNMLQ',INFot,NOUt,LERr,OK)
      INFot = 12
      CALL CUNMLQ('R','N',2,1,0,a,1,x,af,2,w,1,info)
      CALL CHKXER('CUNMLQ',INFot,NOUt,LERr,OK)
!
!     CUNML2
!
      SRNamt = 'CUNML2'
      INFot = 1
      CALL CUNML2('/','N',0,0,0,a,1,x,af,1,w,info)
      CALL CHKXER('CUNML2',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CUNML2('L','/',0,0,0,a,1,x,af,1,w,info)
      CALL CHKXER('CUNML2',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL CUNML2('L','N',-1,0,0,a,1,x,af,1,w,info)
      CALL CHKXER('CUNML2',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL CUNML2('L','N',0,-1,0,a,1,x,af,1,w,info)
      CALL CHKXER('CUNML2',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL CUNML2('L','N',0,0,-1,a,1,x,af,1,w,info)
      CALL CHKXER('CUNML2',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL CUNML2('L','N',0,1,1,a,1,x,af,1,w,info)
      CALL CHKXER('CUNML2',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL CUNML2('R','N',1,0,1,a,1,x,af,1,w,info)
      CALL CHKXER('CUNML2',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL CUNML2('L','N',2,1,2,a,1,x,af,2,w,info)
      CALL CHKXER('CUNML2',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL CUNML2('R','N',1,2,2,a,1,x,af,1,w,info)
      CALL CHKXER('CUNML2',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL CUNML2('L','N',2,1,0,a,2,x,af,1,w,info)
      CALL CHKXER('CUNML2',INFot,NOUt,LERr,OK)
!
!     Print a summary line.
!
      CALL alaesm(Path,OK,NOUt)
!
!
!     End of cerrlq
!
      END SUBROUTINE cerrlq
