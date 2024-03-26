!*==zerrlq.f90  processed by SPAG 7.51RB at 17:35 on  4 Mar 2022
!> \brief \b zerrlq
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE zerrlq( PATH, NUNIT )
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
!> zerrlq tests the error exits for the COMPLEX*16 routines
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
!> \ingroup complex16_lin
!
!  =====================================================================
      SUBROUTINE zerrlq(Path,Nunit)
      IMPLICIT NONE
!*--ZERRLQ59
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
EXTERNAL    ZGELQ2,   ZGELQF,   ZUNGL2
EXTERNAL    ZUNGLQ,   ZUNML2,   ZUNMLQ
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
!     Error exits for LQ factorization
!
!     ZGELQF
!
      SRNamt = 'ZGELQF'
      INFot = 1
      CALL ZGELQF(-1,0,a,1,b,w,1,info)
      CALL CHKXER('ZGELQF',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZGELQF(0,-1,a,1,b,w,1,info)
      CALL CHKXER('ZGELQF',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL ZGELQF(2,1,a,1,b,w,2,info)
      CALL CHKXER('ZGELQF',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL ZGELQF(2,1,a,2,b,w,1,info)
      CALL CHKXER('ZGELQF',INFot,NOUt,LERr,OK)
!
!     ZGELQ2
!
      SRNamt = 'ZGELQ2'
      INFot = 1
      CALL ZGELQ2(-1,0,a,1,b,w,info)
      CALL CHKXER('ZGELQ2',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZGELQ2(0,-1,a,1,b,w,info)
      CALL CHKXER('ZGELQ2',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL ZGELQ2(2,1,a,1,b,w,info)
      CALL CHKXER('ZGELQ2',INFot,NOUt,LERr,OK)
!
!     zgelqs
!
      SRNamt = 'zgelqs'
      INFot = 1
      CALL zgelqs(-1,0,0,a,1,x,b,1,w,1,info)
      CALL CHKXER('zgelqs',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL zgelqs(0,-1,0,a,1,x,b,1,w,1,info)
      CALL CHKXER('zgelqs',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL zgelqs(2,1,0,a,2,x,b,1,w,1,info)
      CALL CHKXER('zgelqs',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL zgelqs(0,0,-1,a,1,x,b,1,w,1,info)
      CALL CHKXER('zgelqs',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL zgelqs(2,2,0,a,1,x,b,2,w,1,info)
      CALL CHKXER('zgelqs',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL zgelqs(1,2,0,a,1,x,b,1,w,1,info)
      CALL CHKXER('zgelqs',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL zgelqs(1,1,2,a,1,x,b,1,w,1,info)
      CALL CHKXER('zgelqs',INFot,NOUt,LERr,OK)
!
!     ZUNGLQ
!
      SRNamt = 'ZUNGLQ'
      INFot = 1
      CALL ZUNGLQ(-1,0,0,a,1,x,w,1,info)
      CALL CHKXER('ZUNGLQ',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZUNGLQ(0,-1,0,a,1,x,w,1,info)
      CALL CHKXER('ZUNGLQ',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZUNGLQ(2,1,0,a,2,x,w,2,info)
      CALL CHKXER('ZUNGLQ',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL ZUNGLQ(0,0,-1,a,1,x,w,1,info)
      CALL CHKXER('ZUNGLQ',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL ZUNGLQ(1,1,2,a,1,x,w,1,info)
      CALL CHKXER('ZUNGLQ',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL ZUNGLQ(2,2,0,a,1,x,w,2,info)
      CALL CHKXER('ZUNGLQ',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL ZUNGLQ(2,2,0,a,2,x,w,1,info)
      CALL CHKXER('ZUNGLQ',INFot,NOUt,LERr,OK)
!
!     ZUNGL2
!
      SRNamt = 'ZUNGL2'
      INFot = 1
      CALL ZUNGL2(-1,0,0,a,1,x,w,info)
      CALL CHKXER('ZUNGL2',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZUNGL2(0,-1,0,a,1,x,w,info)
      CALL CHKXER('ZUNGL2',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZUNGL2(2,1,0,a,2,x,w,info)
      CALL CHKXER('ZUNGL2',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL ZUNGL2(0,0,-1,a,1,x,w,info)
      CALL CHKXER('ZUNGL2',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL ZUNGL2(1,1,2,a,1,x,w,info)
      CALL CHKXER('ZUNGL2',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL ZUNGL2(2,2,0,a,1,x,w,info)
      CALL CHKXER('ZUNGL2',INFot,NOUt,LERr,OK)
!
!     ZUNMLQ
!
      SRNamt = 'ZUNMLQ'
      INFot = 1
      CALL ZUNMLQ('/','N',0,0,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('ZUNMLQ',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZUNMLQ('L','/',0,0,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('ZUNMLQ',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL ZUNMLQ('L','N',-1,0,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('ZUNMLQ',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL ZUNMLQ('L','N',0,-1,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('ZUNMLQ',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL ZUNMLQ('L','N',0,0,-1,a,1,x,af,1,w,1,info)
      CALL CHKXER('ZUNMLQ',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL ZUNMLQ('L','N',0,1,1,a,1,x,af,1,w,1,info)
      CALL CHKXER('ZUNMLQ',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL ZUNMLQ('R','N',1,0,1,a,1,x,af,1,w,1,info)
      CALL CHKXER('ZUNMLQ',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL ZUNMLQ('L','N',2,0,2,a,1,x,af,2,w,1,info)
      CALL CHKXER('ZUNMLQ',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL ZUNMLQ('R','N',0,2,2,a,1,x,af,1,w,1,info)
      CALL CHKXER('ZUNMLQ',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL ZUNMLQ('L','N',2,1,0,a,2,x,af,1,w,1,info)
      CALL CHKXER('ZUNMLQ',INFot,NOUt,LERr,OK)
      INFot = 12
      CALL ZUNMLQ('L','N',1,2,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('ZUNMLQ',INFot,NOUt,LERr,OK)
      INFot = 12
      CALL ZUNMLQ('R','N',2,1,0,a,1,x,af,2,w,1,info)
      CALL CHKXER('ZUNMLQ',INFot,NOUt,LERr,OK)
!
!     ZUNML2
!
      SRNamt = 'ZUNML2'
      INFot = 1
      CALL ZUNML2('/','N',0,0,0,a,1,x,af,1,w,info)
      CALL CHKXER('ZUNML2',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZUNML2('L','/',0,0,0,a,1,x,af,1,w,info)
      CALL CHKXER('ZUNML2',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL ZUNML2('L','N',-1,0,0,a,1,x,af,1,w,info)
      CALL CHKXER('ZUNML2',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL ZUNML2('L','N',0,-1,0,a,1,x,af,1,w,info)
      CALL CHKXER('ZUNML2',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL ZUNML2('L','N',0,0,-1,a,1,x,af,1,w,info)
      CALL CHKXER('ZUNML2',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL ZUNML2('L','N',0,1,1,a,1,x,af,1,w,info)
      CALL CHKXER('ZUNML2',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL ZUNML2('R','N',1,0,1,a,1,x,af,1,w,info)
      CALL CHKXER('ZUNML2',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL ZUNML2('L','N',2,1,2,a,1,x,af,2,w,info)
      CALL CHKXER('ZUNML2',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL ZUNML2('R','N',1,2,2,a,1,x,af,1,w,info)
      CALL CHKXER('ZUNML2',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL ZUNML2('L','N',2,1,0,a,2,x,af,1,w,info)
      CALL CHKXER('ZUNML2',INFot,NOUt,LERr,OK)
!
!     Print a summary line.
!
      CALL alaesm(Path,OK,NOUt)
!
!
!     End of zerrlq
!
      END SUBROUTINE zerrlq
