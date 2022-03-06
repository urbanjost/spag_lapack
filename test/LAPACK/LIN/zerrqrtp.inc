!*==zerrqrtp.f90  processed by SPAG 7.51RB at 17:35 on  4 Mar 2022
!> \brief \b zerrqrtp
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE zerrqrtp( PATH, NUNIT )
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
!> zerrqrtp tests the error exits for the COMPLEX*16 routines
!> that use the QRT decomposition of a triangular-pentagonal matrix.
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
      SUBROUTINE zerrqrtp(Path,Nunit)
      IMPLICIT NONE
!*--ZERRQRTP59
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
      COMPLEX*16 a(NMAX,NMAX) , t(NMAX,NMAX) , w(NMAX) , b(NMAX,NMAX) , &
     &           c(NMAX,NMAX)
!     ..
!     .. External Subroutines ..
EXTERNAL    CHKXER,   ZTPQRT2,   ZTPQRT,   ZTPMQRT
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
            a(i,j) = 1.D0/DCMPLX(DBLE(i+j),0.D0)
            c(i,j) = 1.D0/DCMPLX(DBLE(i+j),0.D0)
            t(i,j) = 1.D0/DCMPLX(DBLE(i+j),0.D0)
         ENDDO
         w(j) = DCMPLX(0.D0,0.D0)
      ENDDO
      OK = .TRUE.
!
!     Error exits for TPQRT factorization
!
!     ZTPQRT
!
      SRNamt = 'ZTPQRT'
      INFot = 1
      CALL ZTPQRT(-1,1,0,1,a,1,b,1,t,1,w,info)
      CALL CHKXER('ZTPQRT',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZTPQRT(1,-1,0,1,a,1,b,1,t,1,w,info)
      CALL CHKXER('ZTPQRT',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL ZTPQRT(0,1,-1,1,a,1,b,1,t,1,w,info)
      CALL CHKXER('ZTPQRT',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL ZTPQRT(0,1,1,1,a,1,b,1,t,1,w,info)
      CALL CHKXER('ZTPQRT',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL ZTPQRT(0,1,0,0,a,1,b,1,t,1,w,info)
      CALL CHKXER('ZTPQRT',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL ZTPQRT(0,1,0,2,a,1,b,1,t,1,w,info)
      CALL CHKXER('ZTPQRT',INFot,NOUt,LERr,OK)
      INFot = 6
      CALL ZTPQRT(1,2,0,2,a,1,b,1,t,1,w,info)
      CALL CHKXER('ZTPQRT',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL ZTPQRT(2,1,0,1,a,1,b,1,t,1,w,info)
      CALL CHKXER('ZTPQRT',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL ZTPQRT(2,2,1,2,a,2,b,2,t,1,w,info)
      CALL CHKXER('ZTPQRT',INFot,NOUt,LERr,OK)
!
!     ZTPQRT2
!
      SRNamt = 'ZTPQRT2'
      INFot = 1
      CALL ZTPQRT2(-1,0,0,a,1,b,1,t,1,info)
      CALL CHKXER('ZTPQRT2',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZTPQRT2(0,-1,0,a,1,b,1,t,1,info)
      CALL CHKXER('ZTPQRT2',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL ZTPQRT2(0,0,-1,a,1,b,1,t,1,info)
      CALL CHKXER('ZTPQRT2',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL ZTPQRT2(2,2,0,a,1,b,2,t,2,info)
      CALL CHKXER('ZTPQRT2',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL ZTPQRT2(2,2,0,a,2,b,1,t,2,info)
      CALL CHKXER('ZTPQRT2',INFot,NOUt,LERr,OK)
      INFot = 9
      CALL ZTPQRT2(2,2,0,a,2,b,2,t,1,info)
      CALL CHKXER('ZTPQRT2',INFot,NOUt,LERr,OK)
!
!     ZTPMQRT
!
      SRNamt = 'ZTPMQRT'
      INFot = 1
      CALL ZTPMQRT('/','N',0,0,0,0,1,a,1,t,1,b,1,c,1,w,info)
      CALL CHKXER('ZTPMQRT',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZTPMQRT('L','/',0,0,0,0,1,a,1,t,1,b,1,c,1,w,info)
      CALL CHKXER('ZTPMQRT',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL ZTPMQRT('L','N',-1,0,0,0,1,a,1,t,1,b,1,c,1,w,info)
      CALL CHKXER('ZTPMQRT',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL ZTPMQRT('L','N',0,-1,0,0,1,a,1,t,1,b,1,c,1,w,info)
      CALL CHKXER('ZTPMQRT',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL ZTPMQRT('L','N',0,0,-1,0,1,a,1,t,1,b,1,c,1,w,info)
      INFot = 6
      CALL ZTPMQRT('L','N',0,0,0,-1,1,a,1,t,1,b,1,c,1,w,info)
      CALL CHKXER('ZTPMQRT',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL ZTPMQRT('L','N',0,0,0,0,0,a,1,t,1,b,1,c,1,w,info)
      CALL CHKXER('ZTPMQRT',INFot,NOUt,LERr,OK)
      INFot = 9
      CALL ZTPMQRT('R','N',1,2,1,1,1,a,1,t,1,b,1,c,1,w,info)
      CALL CHKXER('ZTPMQRT',INFot,NOUt,LERr,OK)
      INFot = 9
      CALL ZTPMQRT('L','N',2,1,1,1,1,a,1,t,1,b,1,c,1,w,info)
      CALL CHKXER('ZTPMQRT',INFot,NOUt,LERr,OK)
      INFot = 11
      CALL ZTPMQRT('R','N',1,1,1,1,1,a,1,t,0,b,1,c,1,w,info)
      CALL CHKXER('ZTPMQRT',INFot,NOUt,LERr,OK)
      INFot = 13
      CALL ZTPMQRT('L','N',1,1,1,1,1,a,1,t,1,b,0,c,1,w,info)
      CALL CHKXER('ZTPMQRT',INFot,NOUt,LERr,OK)
      INFot = 15
      CALL ZTPMQRT('L','N',1,1,1,1,1,a,1,t,1,b,1,c,0,w,info)
      CALL CHKXER('ZTPMQRT',INFot,NOUt,LERr,OK)
!
!     Print a summary line.
!
      CALL alaesm(Path,OK,NOUt)
!
!
!     End of zerrqrt
!
      END SUBROUTINE zerrqrtp
