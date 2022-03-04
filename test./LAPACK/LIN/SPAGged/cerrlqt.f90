!*==cerrlqt.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b CERRLQT
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CERRLQT( PATH, NUNIT )
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
!> CERRLQT tests the error exits for the COMPLEX routines
!> that use the LQT decomposition of a general matrix.
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
      SUBROUTINE CERRLQT(Path,Nunit)
      IMPLICIT NONE
!*--CERRLQT59
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
      EXTERNAL ALAESM , CHKXER , CGELQT3 , CGELQT , CGEMLQT
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
      INTRINSIC REAL , CMPLX
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
            a(i,j) = 1.E0/CMPLX(REAL(i+j),0.E0)
            c(i,j) = 1.E0/CMPLX(REAL(i+j),0.E0)
            t(i,j) = 1.E0/CMPLX(REAL(i+j),0.E0)
         ENDDO
         w(j) = 0.E0
      ENDDO
      OK = .TRUE.
!
!     Error exits for LQT factorization
!
!     CGELQT
!
      SRNamt = 'CGELQT'
      INFot = 1
      CALL CGELQT(-1,0,1,a,1,t,1,w,info)
      CALL CHKXER('CGELQT',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CGELQT(0,-1,1,a,1,t,1,w,info)
      CALL CHKXER('CGELQT',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL CGELQT(0,0,0,a,1,t,1,w,info)
      CALL CHKXER('CGELQT',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL CGELQT(2,1,1,a,1,t,1,w,info)
      CALL CHKXER('CGELQT',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL CGELQT(2,2,2,a,2,t,1,w,info)
      CALL CHKXER('CGELQT',INFot,NOUt,LERr,OK)
!
!     CGELQT3
!
      SRNamt = 'CGELQT3'
      INFot = 1
      CALL CGELQT3(-1,0,a,1,t,1,info)
      CALL CHKXER('CGELQT3',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CGELQT3(0,-1,a,1,t,1,info)
      CALL CHKXER('CGELQT3',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL CGELQT3(2,2,a,1,t,1,info)
      CALL CHKXER('CGELQT3',INFot,NOUt,LERr,OK)
      INFot = 6
      CALL CGELQT3(2,2,a,2,t,1,info)
      CALL CHKXER('CGELQT3',INFot,NOUt,LERr,OK)
!
!     CGEMLQT
!
      SRNamt = 'CGEMLQT'
      INFot = 1
      CALL CGEMLQT('/','N',0,0,0,1,a,1,t,1,c,1,w,info)
      CALL CHKXER('CGEMLQT',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CGEMLQT('L','/',0,0,0,1,a,1,t,1,c,1,w,info)
      CALL CHKXER('CGEMLQT',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL CGEMLQT('L','N',-1,0,0,1,a,1,t,1,c,1,w,info)
      CALL CHKXER('CGEMLQT',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL CGEMLQT('L','N',0,-1,0,1,a,1,t,1,c,1,w,info)
      CALL CHKXER('CGEMLQT',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL CGEMLQT('L','N',0,0,-1,1,a,1,t,1,c,1,w,info)
      CALL CHKXER('CGEMLQT',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL CGEMLQT('R','N',0,0,-1,1,a,1,t,1,c,1,w,info)
      CALL CHKXER('CGEMLQT',INFot,NOUt,LERr,OK)
      INFot = 6
      CALL CGEMLQT('L','N',0,0,0,0,a,1,t,1,c,1,w,info)
      CALL CHKXER('CGEMLQT',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL CGEMLQT('R','N',2,2,2,1,a,1,t,1,c,1,w,info)
      CALL CHKXER('CGEMLQT',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL CGEMLQT('L','N',2,2,2,1,a,1,t,1,c,1,w,info)
      CALL CHKXER('CGEMLQT',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL CGEMLQT('R','N',1,1,1,1,a,1,t,0,c,1,w,info)
      CALL CHKXER('CGEMLQT',INFot,NOUt,LERr,OK)
      INFot = 12
      CALL CGEMLQT('L','N',1,1,1,1,a,1,t,1,c,0,w,info)
      CALL CHKXER('CGEMLQT',INFot,NOUt,LERr,OK)
!
!     Print a summary line.
!
      CALL ALAESM(Path,OK,NOUt)
!
!
!     End of CERRLQT
!
      END SUBROUTINE CERRLQT
