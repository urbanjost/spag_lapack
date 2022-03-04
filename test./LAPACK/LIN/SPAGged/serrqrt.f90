!*==serrqrt.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b SERRQRT
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SERRQRT( PATH, NUNIT )
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
!> SERRQRT tests the error exits for the REAL routines
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
!> \ingroup single_lin
!
!  =====================================================================
      SUBROUTINE SERRQRT(Path,Nunit)
      IMPLICIT NONE
!*--SERRQRT59
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
      REAL a(NMAX,NMAX) , t(NMAX,NMAX) , w(NMAX) , c(NMAX,NMAX)
!     ..
!     .. External Subroutines ..
      EXTERNAL ALAESM , CHKXER , SGEQRT2 , SGEQRT3 , SGEQRT , SGEMQRT
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
      INTRINSIC FLOAT
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
            a(i,j) = 1.0/FLOAT(i+j)
            c(i,j) = 1.0/FLOAT(i+j)
            t(i,j) = 1.0/FLOAT(i+j)
         ENDDO
         w(j) = 0.0
      ENDDO
      OK = .TRUE.
!
!     Error exits for QRT factorization
!
!     SGEQRT
!
      SRNamt = 'SGEQRT'
      INFot = 1
      CALL SGEQRT(-1,0,1,a,1,t,1,w,info)
      CALL CHKXER('SGEQRT',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL SGEQRT(0,-1,1,a,1,t,1,w,info)
      CALL CHKXER('SGEQRT',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL SGEQRT(0,0,0,a,1,t,1,w,info)
      CALL CHKXER('SGEQRT',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL SGEQRT(2,1,1,a,1,t,1,w,info)
      CALL CHKXER('SGEQRT',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL SGEQRT(2,2,2,a,2,t,1,w,info)
      CALL CHKXER('SGEQRT',INFot,NOUt,LERr,OK)
!
!     SGEQRT2
!
      SRNamt = 'SGEQRT2'
      INFot = 1
      CALL SGEQRT2(-1,0,a,1,t,1,info)
      CALL CHKXER('SGEQRT2',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL SGEQRT2(0,-1,a,1,t,1,info)
      CALL CHKXER('SGEQRT2',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL SGEQRT2(2,1,a,1,t,1,info)
      CALL CHKXER('SGEQRT2',INFot,NOUt,LERr,OK)
      INFot = 6
      CALL SGEQRT2(2,2,a,2,t,1,info)
      CALL CHKXER('SGEQRT2',INFot,NOUt,LERr,OK)
!
!     SGEQRT3
!
      SRNamt = 'SGEQRT3'
      INFot = 1
      CALL SGEQRT3(-1,0,a,1,t,1,info)
      CALL CHKXER('SGEQRT3',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL SGEQRT3(0,-1,a,1,t,1,info)
      CALL CHKXER('SGEQRT3',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL SGEQRT3(2,1,a,1,t,1,info)
      CALL CHKXER('SGEQRT3',INFot,NOUt,LERr,OK)
      INFot = 6
      CALL SGEQRT3(2,2,a,2,t,1,info)
      CALL CHKXER('SGEQRT3',INFot,NOUt,LERr,OK)
!
!     SGEMQRT
!
      SRNamt = 'SGEMQRT'
      INFot = 1
      CALL SGEMQRT('/','N',0,0,0,1,a,1,t,1,c,1,w,info)
      CALL CHKXER('SGEMQRT',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL SGEMQRT('L','/',0,0,0,1,a,1,t,1,c,1,w,info)
      CALL CHKXER('SGEMQRT',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL SGEMQRT('L','N',-1,0,0,1,a,1,t,1,c,1,w,info)
      CALL CHKXER('SGEMQRT',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL SGEMQRT('L','N',0,-1,0,1,a,1,t,1,c,1,w,info)
      CALL CHKXER('SGEMQRT',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL SGEMQRT('L','N',0,0,-1,1,a,1,t,1,c,1,w,info)
      CALL CHKXER('SGEMQRT',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL SGEMQRT('R','N',0,0,-1,1,a,1,t,1,c,1,w,info)
      CALL CHKXER('SGEMQRT',INFot,NOUt,LERr,OK)
      INFot = 6
      CALL SGEMQRT('L','N',0,0,0,0,a,1,t,1,c,1,w,info)
      CALL CHKXER('SGEMQRT',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL SGEMQRT('R','N',1,2,1,1,a,1,t,1,c,1,w,info)
      CALL CHKXER('SGEMQRT',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL SGEMQRT('L','N',2,1,1,1,a,1,t,1,c,1,w,info)
      CALL CHKXER('SGEMQRT',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL SGEMQRT('R','N',1,1,1,1,a,1,t,0,c,1,w,info)
      CALL CHKXER('SGEMQRT',INFot,NOUt,LERr,OK)
      INFot = 12
      CALL SGEMQRT('L','N',1,1,1,1,a,1,t,1,c,0,w,info)
      CALL CHKXER('SGEMQRT',INFot,NOUt,LERr,OK)
!
!     Print a summary line.
!
      CALL ALAESM(Path,OK,NOUt)
!
!
!     End of SERRQRT
!
      END SUBROUTINE SERRQRT
