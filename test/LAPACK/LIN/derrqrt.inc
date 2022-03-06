!*==derrqrt.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b derrqrt
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE derrqrt( PATH, NUNIT )
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
!> derrqrt tests the error exits for the DOUBLE PRECISION routines
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
!> \ingroup double_lin
!
!  =====================================================================
      SUBROUTINE derrqrt(Path,Nunit)
      IMPLICIT NONE
!*--DERRQRT59
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
     &                 c(NMAX,NMAX)
!     ..
!     .. External Subroutines ..
EXTERNAL    CHKXER,   DGEQRT2,   DGEQRT3,   DGEQRT,   DGEMQRT
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
!     Error exits for QRT factorization
!
!     DGEQRT
!
      SRNamt = 'DGEQRT'
      INFot = 1
      CALL DGEQRT(-1,0,1,a,1,t,1,w,info)
      CALL CHKXER('DGEQRT',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DGEQRT(0,-1,1,a,1,t,1,w,info)
      CALL CHKXER('DGEQRT',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL DGEQRT(0,0,0,a,1,t,1,w,info)
      CALL CHKXER('DGEQRT',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL DGEQRT(2,1,1,a,1,t,1,w,info)
      CALL CHKXER('DGEQRT',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL DGEQRT(2,2,2,a,2,t,1,w,info)
      CALL CHKXER('DGEQRT',INFot,NOUt,LERr,OK)
!
!     DGEQRT2
!
      SRNamt = 'DGEQRT2'
      INFot = 1
      CALL DGEQRT2(-1,0,a,1,t,1,info)
      CALL CHKXER('DGEQRT2',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DGEQRT2(0,-1,a,1,t,1,info)
      CALL CHKXER('DGEQRT2',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL DGEQRT2(2,1,a,1,t,1,info)
      CALL CHKXER('DGEQRT2',INFot,NOUt,LERr,OK)
      INFot = 6
      CALL DGEQRT2(2,2,a,2,t,1,info)
      CALL CHKXER('DGEQRT2',INFot,NOUt,LERr,OK)
!
!     DGEQRT3
!
      SRNamt = 'DGEQRT3'
      INFot = 1
      CALL DGEQRT3(-1,0,a,1,t,1,info)
      CALL CHKXER('DGEQRT3',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DGEQRT3(0,-1,a,1,t,1,info)
      CALL CHKXER('DGEQRT3',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL DGEQRT3(2,1,a,1,t,1,info)
      CALL CHKXER('DGEQRT3',INFot,NOUt,LERr,OK)
      INFot = 6
      CALL DGEQRT3(2,2,a,2,t,1,info)
      CALL CHKXER('DGEQRT3',INFot,NOUt,LERr,OK)
!
!     DGEMQRT
!
      SRNamt = 'DGEMQRT'
      INFot = 1
      CALL DGEMQRT('/','N',0,0,0,1,a,1,t,1,c,1,w,info)
      CALL CHKXER('DGEMQRT',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DGEMQRT('L','/',0,0,0,1,a,1,t,1,c,1,w,info)
      CALL CHKXER('DGEMQRT',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL DGEMQRT('L','N',-1,0,0,1,a,1,t,1,c,1,w,info)
      CALL CHKXER('DGEMQRT',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL DGEMQRT('L','N',0,-1,0,1,a,1,t,1,c,1,w,info)
      CALL CHKXER('DGEMQRT',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL DGEMQRT('L','N',0,0,-1,1,a,1,t,1,c,1,w,info)
      CALL CHKXER('DGEMQRT',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL DGEMQRT('R','N',0,0,-1,1,a,1,t,1,c,1,w,info)
      CALL CHKXER('DGEMQRT',INFot,NOUt,LERr,OK)
      INFot = 6
      CALL DGEMQRT('L','N',0,0,0,0,a,1,t,1,c,1,w,info)
      CALL CHKXER('DGEMQRT',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL DGEMQRT('R','N',1,2,1,1,a,1,t,1,c,1,w,info)
      CALL CHKXER('DGEMQRT',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL DGEMQRT('L','N',2,1,1,1,a,1,t,1,c,1,w,info)
      CALL CHKXER('DGEMQRT',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL DGEMQRT('R','N',1,1,1,1,a,1,t,0,c,1,w,info)
      CALL CHKXER('DGEMQRT',INFot,NOUt,LERr,OK)
      INFot = 12
      CALL DGEMQRT('L','N',1,1,1,1,a,1,t,1,c,0,w,info)
      CALL CHKXER('DGEMQRT',INFot,NOUt,LERr,OK)
!
!     Print a summary line.
!
      CALL alaesm(Path,OK,NOUt)
!
!
!     End of derrqrt
!
      END SUBROUTINE derrqrt
