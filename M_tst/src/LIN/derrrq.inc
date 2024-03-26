!*==derrrq.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b derrrq
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE derrrq( PATH, NUNIT )
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
!> derrrq tests the error exits for the DOUBLE PRECISION routines
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
!> \ingroup double_lin
!
!  =====================================================================
      SUBROUTINE derrrq(Path,Nunit)
      IMPLICIT NONE
!*--DERRRQ59
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
EXTERNAL    DGERQ2,   DGERQF,   DORGR2
EXTERNAL    DORGRQ,   DORMR2,   DORMRQ
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
!     Error exits for RQ factorization
!
!     DGERQF
!
      SRNamt = 'DGERQF'
      INFot = 1
      CALL DGERQF(-1,0,a,1,b,w,1,info)
      CALL CHKXER('DGERQF',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DGERQF(0,-1,a,1,b,w,1,info)
      CALL CHKXER('DGERQF',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL DGERQF(2,1,a,1,b,w,2,info)
      CALL CHKXER('DGERQF',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL DGERQF(2,1,a,2,b,w,1,info)
      CALL CHKXER('DGERQF',INFot,NOUt,LERr,OK)
!
!     DGERQ2
!
      SRNamt = 'DGERQ2'
      INFot = 1
      CALL DGERQ2(-1,0,a,1,b,w,info)
      CALL CHKXER('DGERQ2',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DGERQ2(0,-1,a,1,b,w,info)
      CALL CHKXER('DGERQ2',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL DGERQ2(2,1,a,1,b,w,info)
      CALL CHKXER('DGERQ2',INFot,NOUt,LERr,OK)
!
!     dgerqs
!
      SRNamt = 'dgerqs'
      INFot = 1
      CALL dgerqs(-1,0,0,a,1,x,b,1,w,1,info)
      CALL CHKXER('dgerqs',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL dgerqs(0,-1,0,a,1,x,b,1,w,1,info)
      CALL CHKXER('dgerqs',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL dgerqs(2,1,0,a,2,x,b,1,w,1,info)
      CALL CHKXER('dgerqs',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL dgerqs(0,0,-1,a,1,x,b,1,w,1,info)
      CALL CHKXER('dgerqs',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL dgerqs(2,2,0,a,1,x,b,2,w,1,info)
      CALL CHKXER('dgerqs',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL dgerqs(2,2,0,a,2,x,b,1,w,1,info)
      CALL CHKXER('dgerqs',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL dgerqs(1,1,2,a,1,x,b,1,w,1,info)
      CALL CHKXER('dgerqs',INFot,NOUt,LERr,OK)
!
!     DORGRQ
!
      SRNamt = 'DORGRQ'
      INFot = 1
      CALL DORGRQ(-1,0,0,a,1,x,w,1,info)
      CALL CHKXER('DORGRQ',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DORGRQ(0,-1,0,a,1,x,w,1,info)
      CALL CHKXER('DORGRQ',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DORGRQ(2,1,0,a,2,x,w,2,info)
      CALL CHKXER('DORGRQ',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL DORGRQ(0,0,-1,a,1,x,w,1,info)
      CALL CHKXER('DORGRQ',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL DORGRQ(1,2,2,a,1,x,w,1,info)
      CALL CHKXER('DORGRQ',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL DORGRQ(2,2,0,a,1,x,w,2,info)
      CALL CHKXER('DORGRQ',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL DORGRQ(2,2,0,a,2,x,w,1,info)
      CALL CHKXER('DORGRQ',INFot,NOUt,LERr,OK)
!
!     DORGR2
!
      SRNamt = 'DORGR2'
      INFot = 1
      CALL DORGR2(-1,0,0,a,1,x,w,info)
      CALL CHKXER('DORGR2',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DORGR2(0,-1,0,a,1,x,w,info)
      CALL CHKXER('DORGR2',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DORGR2(2,1,0,a,2,x,w,info)
      CALL CHKXER('DORGR2',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL DORGR2(0,0,-1,a,1,x,w,info)
      CALL CHKXER('DORGR2',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL DORGR2(1,2,2,a,2,x,w,info)
      CALL CHKXER('DORGR2',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL DORGR2(2,2,0,a,1,x,w,info)
      CALL CHKXER('DORGR2',INFot,NOUt,LERr,OK)
!
!     DORMRQ
!
      SRNamt = 'DORMRQ'
      INFot = 1
      CALL DORMRQ('/','N',0,0,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('DORMRQ',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DORMRQ('L','/',0,0,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('DORMRQ',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL DORMRQ('L','N',-1,0,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('DORMRQ',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL DORMRQ('L','N',0,-1,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('DORMRQ',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL DORMRQ('L','N',0,0,-1,a,1,x,af,1,w,1,info)
      CALL CHKXER('DORMRQ',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL DORMRQ('L','N',0,1,1,a,1,x,af,1,w,1,info)
      CALL CHKXER('DORMRQ',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL DORMRQ('R','N',1,0,1,a,1,x,af,1,w,1,info)
      CALL CHKXER('DORMRQ',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL DORMRQ('L','N',2,1,2,a,1,x,af,2,w,1,info)
      CALL CHKXER('DORMRQ',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL DORMRQ('R','N',1,2,2,a,1,x,af,1,w,1,info)
      CALL CHKXER('DORMRQ',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL DORMRQ('L','N',2,1,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('DORMRQ',INFot,NOUt,LERr,OK)
      INFot = 12
      CALL DORMRQ('L','N',1,2,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('DORMRQ',INFot,NOUt,LERr,OK)
      INFot = 12
      CALL DORMRQ('R','N',2,1,0,a,1,x,af,2,w,1,info)
      CALL CHKXER('DORMRQ',INFot,NOUt,LERr,OK)
!
!     DORMR2
!
      SRNamt = 'DORMR2'
      INFot = 1
      CALL DORMR2('/','N',0,0,0,a,1,x,af,1,w,info)
      CALL CHKXER('DORMR2',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DORMR2('L','/',0,0,0,a,1,x,af,1,w,info)
      CALL CHKXER('DORMR2',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL DORMR2('L','N',-1,0,0,a,1,x,af,1,w,info)
      CALL CHKXER('DORMR2',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL DORMR2('L','N',0,-1,0,a,1,x,af,1,w,info)
      CALL CHKXER('DORMR2',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL DORMR2('L','N',0,0,-1,a,1,x,af,1,w,info)
      CALL CHKXER('DORMR2',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL DORMR2('L','N',0,1,1,a,1,x,af,1,w,info)
      CALL CHKXER('DORMR2',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL DORMR2('R','N',1,0,1,a,1,x,af,1,w,info)
      CALL CHKXER('DORMR2',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL DORMR2('L','N',2,1,2,a,1,x,af,2,w,info)
      CALL CHKXER('DORMR2',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL DORMR2('R','N',1,2,2,a,1,x,af,1,w,info)
      CALL CHKXER('DORMR2',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL DORMR2('L','N',2,1,0,a,1,x,af,1,w,info)
      CALL CHKXER('DORMR2',INFot,NOUt,LERr,OK)
!
!     Print a summary line.
!
      CALL alaesm(Path,OK,NOUt)
!
!
!     End of derrrq
!
      END SUBROUTINE derrrq
