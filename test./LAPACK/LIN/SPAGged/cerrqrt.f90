!*==cerrqrt.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b CERRQRT
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CERRQRT( PATH, NUNIT )
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
!> CERRQRT tests the error exits for the COMPLEX routines
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
!> \ingroup complex_lin
!
!  =====================================================================
      SUBROUTINE CERRQRT(Path,Nunit)
      IMPLICIT NONE
!*--CERRQRT59
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
      COMPLEX a(NMAX,NMAX) , t(NMAX,NMAX) , w(NMAX) , c(NMAX,NMAX)
!     ..
!     .. External Subroutines ..
      EXTERNAL ALAESM , CHKXER , CGEQRT2 , CGEQRT3 , CGEQRT , CGEMQRT
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
      INTRINSIC FLOAT , CMPLX
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
            a(i,j) = 1.0/CMPLX(FLOAT(i+j),0.0)
            c(i,j) = 1.0/CMPLX(FLOAT(i+j),0.0)
            t(i,j) = 1.0/CMPLX(FLOAT(i+j),0.0)
         ENDDO
         w(j) = 0.0
      ENDDO
      OK = .TRUE.
!
!     Error exits for QRT factorization
!
!     CGEQRT
!
      SRNamt = 'CGEQRT'
      INFot = 1
      CALL CGEQRT(-1,0,1,a,1,t,1,w,info)
      CALL CHKXER('CGEQRT',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CGEQRT(0,-1,1,a,1,t,1,w,info)
      CALL CHKXER('CGEQRT',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL CGEQRT(0,0,0,a,1,t,1,w,info)
      CALL CHKXER('CGEQRT',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL CGEQRT(2,1,1,a,1,t,1,w,info)
      CALL CHKXER('CGEQRT',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL CGEQRT(2,2,2,a,2,t,1,w,info)
      CALL CHKXER('CGEQRT',INFot,NOUt,LERr,OK)
!
!     CGEQRT2
!
      SRNamt = 'CGEQRT2'
      INFot = 1
      CALL CGEQRT2(-1,0,a,1,t,1,info)
      CALL CHKXER('CGEQRT2',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CGEQRT2(0,-1,a,1,t,1,info)
      CALL CHKXER('CGEQRT2',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL CGEQRT2(2,1,a,1,t,1,info)
      CALL CHKXER('CGEQRT2',INFot,NOUt,LERr,OK)
      INFot = 6
      CALL CGEQRT2(2,2,a,2,t,1,info)
      CALL CHKXER('CGEQRT2',INFot,NOUt,LERr,OK)
!
!     CGEQRT3
!
      SRNamt = 'CGEQRT3'
      INFot = 1
      CALL CGEQRT3(-1,0,a,1,t,1,info)
      CALL CHKXER('CGEQRT3',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CGEQRT3(0,-1,a,1,t,1,info)
      CALL CHKXER('CGEQRT3',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL CGEQRT3(2,1,a,1,t,1,info)
      CALL CHKXER('CGEQRT3',INFot,NOUt,LERr,OK)
      INFot = 6
      CALL CGEQRT3(2,2,a,2,t,1,info)
      CALL CHKXER('CGEQRT3',INFot,NOUt,LERr,OK)
!
!     CGEMQRT
!
      SRNamt = 'CGEMQRT'
      INFot = 1
      CALL CGEMQRT('/','N',0,0,0,1,a,1,t,1,c,1,w,info)
      CALL CHKXER('CGEMQRT',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CGEMQRT('L','/',0,0,0,1,a,1,t,1,c,1,w,info)
      CALL CHKXER('CGEMQRT',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL CGEMQRT('L','N',-1,0,0,1,a,1,t,1,c,1,w,info)
      CALL CHKXER('CGEMQRT',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL CGEMQRT('L','N',0,-1,0,1,a,1,t,1,c,1,w,info)
      CALL CHKXER('CGEMQRT',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL CGEMQRT('L','N',0,0,-1,1,a,1,t,1,c,1,w,info)
      CALL CHKXER('CGEMQRT',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL CGEMQRT('R','N',0,0,-1,1,a,1,t,1,c,1,w,info)
      CALL CHKXER('CGEMQRT',INFot,NOUt,LERr,OK)
      INFot = 6
      CALL CGEMQRT('L','N',0,0,0,0,a,1,t,1,c,1,w,info)
      CALL CHKXER('CGEMQRT',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL CGEMQRT('R','N',1,2,1,1,a,1,t,1,c,1,w,info)
      CALL CHKXER('CGEMQRT',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL CGEMQRT('L','N',2,1,1,1,a,1,t,1,c,1,w,info)
      CALL CHKXER('CGEMQRT',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL CGEMQRT('R','N',1,1,1,1,a,1,t,0,c,1,w,info)
      CALL CHKXER('CGEMQRT',INFot,NOUt,LERr,OK)
      INFot = 12
      CALL CGEMQRT('L','N',1,1,1,1,a,1,t,1,c,0,w,info)
      CALL CHKXER('CGEMQRT',INFot,NOUt,LERr,OK)
!
!     Print a summary line.
!
      CALL ALAESM(Path,OK,NOUt)
!
!
!     End of CERRQRT
!
      END SUBROUTINE CERRQRT
