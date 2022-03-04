!*==zerrrq.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b ZERRRQ
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZERRRQ( PATH, NUNIT )
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
!> ZERRRQ tests the error exits for the COMPLEX*16 routines
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
!> \ingroup complex16_lin
!
!  =====================================================================
      SUBROUTINE ZERRRQ(Path,Nunit)
      IMPLICIT NONE
!*--ZERRRQ59
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
      EXTERNAL ALAESM , CHKXER , ZGERQ2 , ZGERQF , ZGERQS , ZUNGR2 ,    &
     &         ZUNGRQ , ZUNMR2 , ZUNMRQ
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
!     Error exits for RQ factorization
!
!     ZGERQF
!
      SRNamt = 'ZGERQF'
      INFot = 1
      CALL ZGERQF(-1,0,a,1,b,w,1,info)
      CALL CHKXER('ZGERQF',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZGERQF(0,-1,a,1,b,w,1,info)
      CALL CHKXER('ZGERQF',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL ZGERQF(2,1,a,1,b,w,2,info)
      CALL CHKXER('ZGERQF',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL ZGERQF(2,1,a,2,b,w,1,info)
      CALL CHKXER('ZGERQF',INFot,NOUt,LERr,OK)
!
!     ZGERQ2
!
      SRNamt = 'ZGERQ2'
      INFot = 1
      CALL ZGERQ2(-1,0,a,1,b,w,info)
      CALL CHKXER('ZGERQ2',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZGERQ2(0,-1,a,1,b,w,info)
      CALL CHKXER('ZGERQ2',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL ZGERQ2(2,1,a,1,b,w,info)
      CALL CHKXER('ZGERQ2',INFot,NOUt,LERr,OK)
!
!     ZGERQS
!
      SRNamt = 'ZGERQS'
      INFot = 1
      CALL ZGERQS(-1,0,0,a,1,x,b,1,w,1,info)
      CALL CHKXER('ZGERQS',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZGERQS(0,-1,0,a,1,x,b,1,w,1,info)
      CALL CHKXER('ZGERQS',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZGERQS(2,1,0,a,2,x,b,1,w,1,info)
      CALL CHKXER('ZGERQS',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL ZGERQS(0,0,-1,a,1,x,b,1,w,1,info)
      CALL CHKXER('ZGERQS',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL ZGERQS(2,2,0,a,1,x,b,2,w,1,info)
      CALL CHKXER('ZGERQS',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL ZGERQS(2,2,0,a,2,x,b,1,w,1,info)
      CALL CHKXER('ZGERQS',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL ZGERQS(1,1,2,a,1,x,b,1,w,1,info)
      CALL CHKXER('ZGERQS',INFot,NOUt,LERr,OK)
!
!     ZUNGRQ
!
      SRNamt = 'ZUNGRQ'
      INFot = 1
      CALL ZUNGRQ(-1,0,0,a,1,x,w,1,info)
      CALL CHKXER('ZUNGRQ',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZUNGRQ(0,-1,0,a,1,x,w,1,info)
      CALL CHKXER('ZUNGRQ',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZUNGRQ(2,1,0,a,2,x,w,2,info)
      CALL CHKXER('ZUNGRQ',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL ZUNGRQ(0,0,-1,a,1,x,w,1,info)
      CALL CHKXER('ZUNGRQ',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL ZUNGRQ(1,2,2,a,1,x,w,1,info)
      CALL CHKXER('ZUNGRQ',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL ZUNGRQ(2,2,0,a,1,x,w,2,info)
      CALL CHKXER('ZUNGRQ',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL ZUNGRQ(2,2,0,a,2,x,w,1,info)
      CALL CHKXER('ZUNGRQ',INFot,NOUt,LERr,OK)
!
!     ZUNGR2
!
      SRNamt = 'ZUNGR2'
      INFot = 1
      CALL ZUNGR2(-1,0,0,a,1,x,w,info)
      CALL CHKXER('ZUNGR2',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZUNGR2(0,-1,0,a,1,x,w,info)
      CALL CHKXER('ZUNGR2',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZUNGR2(2,1,0,a,2,x,w,info)
      CALL CHKXER('ZUNGR2',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL ZUNGR2(0,0,-1,a,1,x,w,info)
      CALL CHKXER('ZUNGR2',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL ZUNGR2(1,2,2,a,2,x,w,info)
      CALL CHKXER('ZUNGR2',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL ZUNGR2(2,2,0,a,1,x,w,info)
      CALL CHKXER('ZUNGR2',INFot,NOUt,LERr,OK)
!
!     ZUNMRQ
!
      SRNamt = 'ZUNMRQ'
      INFot = 1
      CALL ZUNMRQ('/','N',0,0,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('ZUNMRQ',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZUNMRQ('L','/',0,0,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('ZUNMRQ',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL ZUNMRQ('L','N',-1,0,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('ZUNMRQ',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL ZUNMRQ('L','N',0,-1,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('ZUNMRQ',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL ZUNMRQ('L','N',0,0,-1,a,1,x,af,1,w,1,info)
      CALL CHKXER('ZUNMRQ',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL ZUNMRQ('L','N',0,1,1,a,1,x,af,1,w,1,info)
      CALL CHKXER('ZUNMRQ',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL ZUNMRQ('R','N',1,0,1,a,1,x,af,1,w,1,info)
      CALL CHKXER('ZUNMRQ',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL ZUNMRQ('L','N',2,1,2,a,1,x,af,2,w,1,info)
      CALL CHKXER('ZUNMRQ',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL ZUNMRQ('R','N',1,2,2,a,1,x,af,1,w,1,info)
      CALL CHKXER('ZUNMRQ',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL ZUNMRQ('L','N',2,1,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('ZUNMRQ',INFot,NOUt,LERr,OK)
      INFot = 12
      CALL ZUNMRQ('L','N',1,2,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('ZUNMRQ',INFot,NOUt,LERr,OK)
      INFot = 12
      CALL ZUNMRQ('R','N',2,1,0,a,1,x,af,2,w,1,info)
      CALL CHKXER('ZUNMRQ',INFot,NOUt,LERr,OK)
!
!     ZUNMR2
!
      SRNamt = 'ZUNMR2'
      INFot = 1
      CALL ZUNMR2('/','N',0,0,0,a,1,x,af,1,w,info)
      CALL CHKXER('ZUNMR2',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZUNMR2('L','/',0,0,0,a,1,x,af,1,w,info)
      CALL CHKXER('ZUNMR2',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL ZUNMR2('L','N',-1,0,0,a,1,x,af,1,w,info)
      CALL CHKXER('ZUNMR2',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL ZUNMR2('L','N',0,-1,0,a,1,x,af,1,w,info)
      CALL CHKXER('ZUNMR2',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL ZUNMR2('L','N',0,0,-1,a,1,x,af,1,w,info)
      CALL CHKXER('ZUNMR2',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL ZUNMR2('L','N',0,1,1,a,1,x,af,1,w,info)
      CALL CHKXER('ZUNMR2',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL ZUNMR2('R','N',1,0,1,a,1,x,af,1,w,info)
      CALL CHKXER('ZUNMR2',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL ZUNMR2('L','N',2,1,2,a,1,x,af,2,w,info)
      CALL CHKXER('ZUNMR2',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL ZUNMR2('R','N',1,2,2,a,1,x,af,1,w,info)
      CALL CHKXER('ZUNMR2',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL ZUNMR2('L','N',2,1,0,a,1,x,af,1,w,info)
      CALL CHKXER('ZUNMR2',INFot,NOUt,LERr,OK)
!
!     Print a summary line.
!
      CALL ALAESM(Path,OK,NOUt)
!
!
!     End of ZERRRQ
!
      END SUBROUTINE ZERRRQ