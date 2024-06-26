!*==zerrqrt.f90  processed by SPAG 7.51RB at 17:35 on  4 Mar 2022
!> \brief \b zerrqrt
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE zerrqrt( PATH, NUNIT )
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
!> zerrqrt tests the error exits for the COMPLEX*16 routines
!> that use the QRT decomposition of a general matrix.
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
      SUBROUTINE zerrqrt(Path,Nunit)
      IMPLICIT NONE
!*--ZERRQRT59
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
      COMPLEX*16 a(NMAX,NMAX) , t(NMAX,NMAX) , w(NMAX) , c(NMAX,NMAX)
!     ..
!     .. External Subroutines ..
EXTERNAL    ZGEQRT2,   ZGEQRT3,   ZGEQRT,   ZGEMQRT
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
         w(j) = 0.D0
      ENDDO
      OK = .TRUE.
!
!     Error exits for QRT factorization
!
!     ZGEQRT
!
      SRNamt = 'ZGEQRT'
      INFot = 1
      CALL ZGEQRT(-1,0,1,a,1,t,1,w,info)
      CALL CHKXER('ZGEQRT',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZGEQRT(0,-1,1,a,1,t,1,w,info)
      CALL CHKXER('ZGEQRT',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL ZGEQRT(0,0,0,a,1,t,1,w,info)
      CALL CHKXER('ZGEQRT',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL ZGEQRT(2,1,1,a,1,t,1,w,info)
      CALL CHKXER('ZGEQRT',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL ZGEQRT(2,2,2,a,2,t,1,w,info)
      CALL CHKXER('ZGEQRT',INFot,NOUt,LERr,OK)
!
!     ZGEQRT2
!
      SRNamt = 'ZGEQRT2'
      INFot = 1
      CALL ZGEQRT2(-1,0,a,1,t,1,info)
      CALL CHKXER('ZGEQRT2',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZGEQRT2(0,-1,a,1,t,1,info)
      CALL CHKXER('ZGEQRT2',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL ZGEQRT2(2,1,a,1,t,1,info)
      CALL CHKXER('ZGEQRT2',INFot,NOUt,LERr,OK)
      INFot = 6
      CALL ZGEQRT2(2,2,a,2,t,1,info)
      CALL CHKXER('ZGEQRT2',INFot,NOUt,LERr,OK)
!
!     ZGEQRT3
!
      SRNamt = 'ZGEQRT3'
      INFot = 1
      CALL ZGEQRT3(-1,0,a,1,t,1,info)
      CALL CHKXER('ZGEQRT3',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZGEQRT3(0,-1,a,1,t,1,info)
      CALL CHKXER('ZGEQRT3',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL ZGEQRT3(2,1,a,1,t,1,info)
      CALL CHKXER('ZGEQRT3',INFot,NOUt,LERr,OK)
      INFot = 6
      CALL ZGEQRT3(2,2,a,2,t,1,info)
      CALL CHKXER('ZGEQRT3',INFot,NOUt,LERr,OK)
!
!     ZGEMQRT
!
      SRNamt = 'ZGEMQRT'
      INFot = 1
      CALL ZGEMQRT('/','N',0,0,0,1,a,1,t,1,c,1,w,info)
      CALL CHKXER('ZGEMQRT',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZGEMQRT('L','/',0,0,0,1,a,1,t,1,c,1,w,info)
      CALL CHKXER('ZGEMQRT',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL ZGEMQRT('L','N',-1,0,0,1,a,1,t,1,c,1,w,info)
      CALL CHKXER('ZGEMQRT',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL ZGEMQRT('L','N',0,-1,0,1,a,1,t,1,c,1,w,info)
      CALL CHKXER('ZGEMQRT',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL ZGEMQRT('L','N',0,0,-1,1,a,1,t,1,c,1,w,info)
      CALL CHKXER('ZGEMQRT',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL ZGEMQRT('R','N',0,0,-1,1,a,1,t,1,c,1,w,info)
      CALL CHKXER('ZGEMQRT',INFot,NOUt,LERr,OK)
      INFot = 6
      CALL ZGEMQRT('L','N',0,0,0,0,a,1,t,1,c,1,w,info)
      CALL CHKXER('ZGEMQRT',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL ZGEMQRT('R','N',1,2,1,1,a,1,t,1,c,1,w,info)
      CALL CHKXER('ZGEMQRT',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL ZGEMQRT('L','N',2,1,1,1,a,1,t,1,c,1,w,info)
      CALL CHKXER('ZGEMQRT',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL ZGEMQRT('R','N',1,1,1,1,a,1,t,0,c,1,w,info)
      CALL CHKXER('ZGEMQRT',INFot,NOUt,LERr,OK)
      INFot = 12
      CALL ZGEMQRT('L','N',1,1,1,1,a,1,t,1,c,0,w,info)
      CALL CHKXER('ZGEMQRT',INFot,NOUt,LERr,OK)
!
!     Print a summary line.
!
      CALL alaesm(Path,OK,NOUt)
!
!
!     End of zerrqrt
!
      END SUBROUTINE zerrqrt
