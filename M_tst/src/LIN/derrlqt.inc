!*==derrlqt.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b DERLQT
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE derrlqt( PATH, NUNIT )
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
!> derrlqt tests the error exits for the DOUBLE PRECISION routines
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
      SUBROUTINE derrlqt(Path,Nunit)
      IMPLICIT NONE
!*--DERRLQT59
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
EXTERNAL    DGELQT3,   DGELQT,   DGEMLQT
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
!     Error exits for LQT factorization
!
!     DGELQT
!
      SRNamt = 'DGELQT'
      INFot = 1
      CALL DGELQT(-1,0,1,a,1,t,1,w,info)
      CALL CHKXER('DGELQT',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DGELQT(0,-1,1,a,1,t,1,w,info)
      CALL CHKXER('DGELQT',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL DGELQT(0,0,0,a,1,t,1,w,info)
      CALL CHKXER('DGELQT',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL DGELQT(2,1,1,a,1,t,1,w,info)
      CALL CHKXER('DGELQT',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL DGELQT(2,2,2,a,2,t,1,w,info)
      CALL CHKXER('DGELQT',INFot,NOUt,LERr,OK)
!
!     DGELQT3
!
      SRNamt = 'DGELQT3'
      INFot = 1
      CALL DGELQT3(-1,0,a,1,t,1,info)
      CALL CHKXER('DGELQT3',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DGELQT3(0,-1,a,1,t,1,info)
      CALL CHKXER('DGELQT3',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL DGELQT3(2,2,a,1,t,1,info)
      CALL CHKXER('DGELQT3',INFot,NOUt,LERr,OK)
      INFot = 6
      CALL DGELQT3(2,2,a,2,t,1,info)
      CALL CHKXER('DGELQT3',INFot,NOUt,LERr,OK)
!
!     DGEMLQT
!
      SRNamt = 'DGEMLQT'
      INFot = 1
      CALL DGEMLQT('/','N',0,0,0,1,a,1,t,1,c,1,w,info)
      CALL CHKXER('DGEMLQT',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DGEMLQT('L','/',0,0,0,1,a,1,t,1,c,1,w,info)
      CALL CHKXER('DGEMLQT',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL DGEMLQT('L','N',-1,0,0,1,a,1,t,1,c,1,w,info)
      CALL CHKXER('DGEMLQT',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL DGEMLQT('L','N',0,-1,0,1,a,1,t,1,c,1,w,info)
      CALL CHKXER('DGEMLQT',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL DGEMLQT('L','N',0,0,-1,1,a,1,t,1,c,1,w,info)
      CALL CHKXER('DGEMLQT',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL DGEMLQT('R','N',0,0,-1,1,a,1,t,1,c,1,w,info)
      CALL CHKXER('DGEMLQT',INFot,NOUt,LERr,OK)
      INFot = 6
      CALL DGEMLQT('L','N',0,0,0,0,a,1,t,1,c,1,w,info)
      CALL CHKXER('DGEMLQT',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL DGEMLQT('R','N',2,2,2,1,a,1,t,1,c,1,w,info)
      CALL CHKXER('DGEMLQT',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL DGEMLQT('L','N',2,2,2,1,a,1,t,1,c,1,w,info)
      CALL CHKXER('DGEMLQT',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL DGEMLQT('R','N',1,1,1,1,a,1,t,0,c,1,w,info)
      CALL CHKXER('DGEMLQT',INFot,NOUt,LERr,OK)
      INFot = 12
      CALL DGEMLQT('L','N',1,1,1,1,a,1,t,1,c,0,w,info)
      CALL CHKXER('DGEMLQT',INFot,NOUt,LERr,OK)
!
!     Print a summary line.
!
      CALL alaesm(Path,OK,NOUt)
!
!
!     End of derrlqt
!
      END SUBROUTINE derrlqt
