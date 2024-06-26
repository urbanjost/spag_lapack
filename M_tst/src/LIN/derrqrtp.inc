!*==derrqrtp.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b derrqrtp
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE derrqrtp( PATH, NUNIT )
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
!> derrqrtp tests the error exits for the REAL routines
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
!> \ingroup double_lin
!
!  =====================================================================
      SUBROUTINE derrqrtp(Path,Nunit)
      IMPLICIT NONE
!*--DERRQRTP59
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
      DOUBLE PRECISION a(NMAX,NMAX) , t(NMAX,NMAX) , w(NMAX) ,          &
     &                 b(NMAX,NMAX) , c(NMAX,NMAX)
!     ..
!     .. External Subroutines ..
EXTERNAL    DTPQRT2,   DTPQRT,   DTPMQRT
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
         w(j) = 0.0
      ENDDO
      OK = .TRUE.
!
!     Error exits for TPQRT factorization
!
!     DTPQRT
!
      SRNamt = 'DTPQRT'
      INFot = 1
      CALL DTPQRT(-1,1,0,1,a,1,b,1,t,1,w,info)
      CALL CHKXER('DTPQRT',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DTPQRT(1,-1,0,1,a,1,b,1,t,1,w,info)
      CALL CHKXER('DTPQRT',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL DTPQRT(0,1,-1,1,a,1,b,1,t,1,w,info)
      CALL CHKXER('DTPQRT',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL DTPQRT(0,1,1,1,a,1,b,1,t,1,w,info)
      CALL CHKXER('DTPQRT',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL DTPQRT(0,1,0,0,a,1,b,1,t,1,w,info)
      CALL CHKXER('DTPQRT',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL DTPQRT(0,1,0,2,a,1,b,1,t,1,w,info)
      CALL CHKXER('DTPQRT',INFot,NOUt,LERr,OK)
      INFot = 6
      CALL DTPQRT(1,2,0,2,a,1,b,1,t,1,w,info)
      CALL CHKXER('DTPQRT',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL DTPQRT(2,1,0,1,a,1,b,1,t,1,w,info)
      CALL CHKXER('DTPQRT',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL DTPQRT(2,2,1,2,a,2,b,2,t,1,w,info)
      CALL CHKXER('DTPQRT',INFot,NOUt,LERr,OK)
!
!     DTPQRT2
!
      SRNamt = 'DTPQRT2'
      INFot = 1
      CALL DTPQRT2(-1,0,0,a,1,b,1,t,1,info)
      CALL CHKXER('DTPQRT2',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DTPQRT2(0,-1,0,a,1,b,1,t,1,info)
      CALL CHKXER('DTPQRT2',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL DTPQRT2(0,0,-1,a,1,b,1,t,1,info)
      CALL CHKXER('DTPQRT2',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL DTPQRT2(2,2,0,a,1,b,2,t,2,info)
      CALL CHKXER('DTPQRT2',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL DTPQRT2(2,2,0,a,2,b,1,t,2,info)
      CALL CHKXER('DTPQRT2',INFot,NOUt,LERr,OK)
      INFot = 9
      CALL DTPQRT2(2,2,0,a,2,b,2,t,1,info)
      CALL CHKXER('DTPQRT2',INFot,NOUt,LERr,OK)
!
!     DTPMQRT
!
      SRNamt = 'DTPMQRT'
      INFot = 1
      CALL DTPMQRT('/','N',0,0,0,0,1,a,1,t,1,b,1,c,1,w,info)
      CALL CHKXER('DTPMQRT',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DTPMQRT('L','/',0,0,0,0,1,a,1,t,1,b,1,c,1,w,info)
      CALL CHKXER('DTPMQRT',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL DTPMQRT('L','N',-1,0,0,0,1,a,1,t,1,b,1,c,1,w,info)
      CALL CHKXER('DTPMQRT',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL DTPMQRT('L','N',0,-1,0,0,1,a,1,t,1,b,1,c,1,w,info)
      CALL CHKXER('DTPMQRT',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL DTPMQRT('L','N',0,0,-1,0,1,a,1,t,1,b,1,c,1,w,info)
      INFot = 6
      CALL DTPMQRT('L','N',0,0,0,-1,1,a,1,t,1,b,1,c,1,w,info)
      CALL CHKXER('DTPMQRT',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL DTPMQRT('L','N',0,0,0,0,0,a,1,t,1,b,1,c,1,w,info)
      CALL CHKXER('DTPMQRT',INFot,NOUt,LERr,OK)
      INFot = 9
      CALL DTPMQRT('R','N',1,2,1,1,1,a,1,t,1,b,1,c,1,w,info)
      CALL CHKXER('DTPMQRT',INFot,NOUt,LERr,OK)
      INFot = 9
      CALL DTPMQRT('L','N',2,1,1,1,1,a,1,t,1,b,1,c,1,w,info)
      CALL CHKXER('DTPMQRT',INFot,NOUt,LERr,OK)
      INFot = 11
      CALL DTPMQRT('R','N',1,1,1,1,1,a,1,t,0,b,1,c,1,w,info)
      CALL CHKXER('DTPMQRT',INFot,NOUt,LERr,OK)
      INFot = 13
      CALL DTPMQRT('L','N',1,1,1,1,1,a,1,t,1,b,0,c,1,w,info)
      CALL CHKXER('DTPMQRT',INFot,NOUt,LERr,OK)
      INFot = 15
      CALL DTPMQRT('L','N',1,1,1,1,1,a,1,t,1,b,1,c,0,w,info)
      CALL CHKXER('DTPMQRT',INFot,NOUt,LERr,OK)
!
!     Print a summary line.
!
      CALL alaesm(Path,OK,NOUt)
!
!
!     End of derrqrt
!
      END SUBROUTINE derrqrtp
