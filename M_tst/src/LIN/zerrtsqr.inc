!*==zerrtsqr.f90  processed by SPAG 7.51RB at 17:35 on  4 Mar 2022
!> \brief \b zerrtsqr
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE zerrtsqr( PATH, NUNIT )
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
!> zerrtsqr tests the error exits for the ZOUBLE PRECISION routines
!> that use the TSQR decomposition of a general matrix.
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
!> \author Univ. of Colorado Zenver
!> \author NAG Ltd.
!
!> \date December 2016
!
!> \ingroup double_lin
!
!  =====================================================================
      SUBROUTINE zerrtsqr(Path,Nunit)
      IMPLICIT NONE
!*--ZERRTSQR59
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
      INTEGER i , info , j , nb
!     ..
!     .. Local Arrays ..
      COMPLEX*16 a(NMAX,NMAX) , t(NMAX,NMAX) , w(NMAX) , c(NMAX,NMAX) , &
     &           tau(NMAX)
!     ..
!     .. External Subroutines ..
EXTERNAL    ZGEQR,   ZGEMQR,   ZGELQ,   ZGEMLQ
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
            c(i,j) = 1.D0/DBLE(i+j)
            t(i,j) = 1.D0/DBLE(i+j)
         ENDDO
         w(j) = 0.D0
      ENDDO
      OK = .TRUE.
!
!     Error exits for TS factorization
!
!     ZGEQR
!
      SRNamt = 'ZGEQR'
      INFot = 1
      CALL ZGEQR(-1,0,a,1,tau,1,w,1,info)
      CALL CHKXER('ZGEQR',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZGEQR(0,-1,a,1,tau,1,w,1,info)
      CALL CHKXER('ZGEQR',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL ZGEQR(1,1,a,0,tau,1,w,1,info)
      CALL CHKXER('ZGEQR',INFot,NOUt,LERr,OK)
      INFot = 6
      CALL ZGEQR(3,2,a,3,tau,1,w,1,info)
      CALL CHKXER('ZGEQR',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL ZGEQR(3,2,a,3,tau,8,w,0,info)
      CALL CHKXER('ZGEQR',INFot,NOUt,LERr,OK)
!
!     ZGEMQR
!
      tau(1) = 1
      tau(2) = 1
      SRNamt = 'ZGEMQR'
      nb = 1
      INFot = 1
      CALL ZGEMQR('/','N',0,0,0,a,1,tau,1,c,1,w,1,info)
      CALL CHKXER('ZGEMQR',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZGEMQR('L','/',0,0,0,a,1,tau,1,c,1,w,1,info)
      CALL CHKXER('ZGEMQR',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL ZGEMQR('L','N',-1,0,0,a,1,tau,1,c,1,w,1,info)
      CALL CHKXER('ZGEMQR',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL ZGEMQR('L','N',0,-1,0,a,1,tau,1,c,1,w,1,info)
      CALL CHKXER('ZGEMQR',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL ZGEMQR('L','N',0,0,-1,a,1,tau,1,c,1,w,1,info)
      CALL CHKXER('ZGEMQR',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL ZGEMQR('R','N',0,0,-1,a,1,tau,1,c,1,w,1,info)
      CALL CHKXER('ZGEMQR',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL ZGEMQR('L','N',2,1,0,a,0,tau,1,c,1,w,1,info)
      CALL CHKXER('ZGEMQR',INFot,NOUt,LERr,OK)
      INFot = 9
      CALL ZGEMQR('R','N',2,2,1,a,2,tau,0,c,1,w,1,info)
      CALL CHKXER('ZGEMQR',INFot,NOUt,LERr,OK)
      INFot = 9
      CALL ZGEMQR('L','N',2,2,1,a,2,tau,0,c,1,w,1,info)
      CALL CHKXER('ZGEMQR',INFot,NOUt,LERr,OK)
      INFot = 11
      CALL ZGEMQR('L','N',2,1,1,a,2,tau,6,c,0,w,1,info)
      CALL CHKXER('ZGEMQR',INFot,NOUt,LERr,OK)
      INFot = 13
      CALL ZGEMQR('L','N',2,2,1,a,2,tau,6,c,2,w,0,info)
      CALL CHKXER('ZGEMQR',INFot,NOUt,LERr,OK)
!
!     ZGELQ
!
      SRNamt = 'ZGELQ'
      INFot = 1
      CALL ZGELQ(-1,0,a,1,tau,1,w,1,info)
      CALL CHKXER('ZGELQ',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZGELQ(0,-1,a,1,tau,1,w,1,info)
      CALL CHKXER('ZGELQ',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL ZGELQ(1,1,a,0,tau,1,w,1,info)
      CALL CHKXER('ZGELQ',INFot,NOUt,LERr,OK)
      INFot = 6
      CALL ZGELQ(2,3,a,3,tau,1,w,1,info)
      CALL CHKXER('ZGELQ',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL ZGELQ(2,3,a,3,tau,8,w,0,info)
      CALL CHKXER('ZGELQ',INFot,NOUt,LERr,OK)
!
!     ZGEMLQ
!
      tau(1) = 1
      tau(2) = 1
      SRNamt = 'ZGEMLQ'
      nb = 1
      INFot = 1
      CALL ZGEMLQ('/','N',0,0,0,a,1,tau,1,c,1,w,1,info)
      CALL CHKXER('ZGEMLQ',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZGEMLQ('L','/',0,0,0,a,1,tau,1,c,1,w,1,info)
      CALL CHKXER('ZGEMLQ',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL ZGEMLQ('L','N',-1,0,0,a,1,tau,1,c,1,w,1,info)
      CALL CHKXER('ZGEMLQ',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL ZGEMLQ('L','N',0,-1,0,a,1,tau,1,c,1,w,1,info)
      CALL CHKXER('ZGEMLQ',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL ZGEMLQ('L','N',0,0,-1,a,1,tau,1,c,1,w,1,info)
      CALL CHKXER('ZGEMLQ',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL ZGEMLQ('R','N',0,0,-1,a,1,tau,1,c,1,w,1,info)
      CALL CHKXER('ZGEMLQ',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL ZGEMLQ('L','N',1,2,0,a,0,tau,1,c,1,w,1,info)
      CALL CHKXER('ZGEMLQ',INFot,NOUt,LERr,OK)
      INFot = 9
      CALL ZGEMLQ('R','N',2,2,1,a,1,tau,0,c,1,w,1,info)
      CALL CHKXER('ZGEMLQ',INFot,NOUt,LERr,OK)
      INFot = 9
      CALL ZGEMLQ('L','N',2,2,1,a,1,tau,0,c,1,w,1,info)
      CALL CHKXER('ZGEMLQ',INFot,NOUt,LERr,OK)
      INFot = 11
      CALL ZGEMLQ('L','N',1,2,1,a,1,tau,6,c,0,w,1,info)
      CALL CHKXER('ZGEMLQ',INFot,NOUt,LERr,OK)
      INFot = 13
      CALL ZGEMLQ('L','N',2,2,1,a,2,tau,6,c,2,w,0,info)
      CALL CHKXER('ZGEMLQ',INFot,NOUt,LERr,OK)
!
!     Print a summary line.
!
      CALL alaesm(Path,OK,NOUt)
!
!
!     End of derrtsqr
!
      END SUBROUTINE zerrtsqr
