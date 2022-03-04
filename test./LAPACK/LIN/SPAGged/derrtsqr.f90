!*==derrtsqr.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b DERRTSQR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DERRTSQR( PATH, NUNIT )
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
!> DERRTSQR tests the error exits for the DOUBLE PRECISION routines
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
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \date December 2016
!
!> \ingroup double_lin
!
!  =====================================================================
      SUBROUTINE DERRTSQR(Path,Nunit)
      IMPLICIT NONE
!*--DERRTSQR59
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
      DOUBLE PRECISION a(NMAX,NMAX) , t(NMAX,NMAX) , w(NMAX) ,          &
     &                 c(NMAX,NMAX) , tau(NMAX*2)
!     ..
!     .. External Subroutines ..
      EXTERNAL ALAESM , CHKXER , DGEQR , DGEMQR , DGELQ , DGEMLQ
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
!     DGEQR
!
      SRNamt = 'DGEQR'
      INFot = 1
      CALL DGEQR(-1,0,a,1,tau,1,w,1,info)
      CALL CHKXER('DGEQR',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DGEQR(0,-1,a,1,tau,1,w,1,info)
      CALL CHKXER('DGEQR',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL DGEQR(1,1,a,0,tau,1,w,1,info)
      CALL CHKXER('DGEQR',INFot,NOUt,LERr,OK)
      INFot = 6
      CALL DGEQR(3,2,a,3,tau,1,w,1,info)
      CALL CHKXER('DGEQR',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL DGEQR(3,2,a,3,tau,7,w,0,info)
      CALL CHKXER('DGEQR',INFot,NOUt,LERr,OK)
!
!     DGEMQR
!
      tau(1) = 1
      tau(2) = 1
      tau(3) = 1
      tau(4) = 1
      SRNamt = 'DGEMQR'
      nb = 1
      INFot = 1
      CALL DGEMQR('/','N',0,0,0,a,1,tau,1,c,1,w,1,info)
      CALL CHKXER('DGEMQR',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DGEMQR('L','/',0,0,0,a,1,tau,1,c,1,w,1,info)
      CALL CHKXER('DGEMQR',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL DGEMQR('L','N',-1,0,0,a,1,tau,1,c,1,w,1,info)
      CALL CHKXER('DGEMQR',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL DGEMQR('L','N',0,-1,0,a,1,tau,1,c,1,w,1,info)
      CALL CHKXER('DGEMQR',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL DGEMQR('L','N',0,0,-1,a,1,tau,1,c,1,w,1,info)
      CALL CHKXER('DGEMQR',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL DGEMQR('R','N',0,0,-1,a,1,tau,1,c,1,w,1,info)
      CALL CHKXER('DGEMQR',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL DGEMQR('L','N',2,1,0,a,0,tau,1,c,1,w,1,info)
      CALL CHKXER('DGEMQR',INFot,NOUt,LERr,OK)
      INFot = 9
      CALL DGEMQR('R','N',2,2,1,a,2,tau,0,c,1,w,1,info)
      CALL CHKXER('DGEMQR',INFot,NOUt,LERr,OK)
      INFot = 9
      CALL DGEMQR('L','N',2,2,1,a,2,tau,0,c,1,w,1,info)
      CALL CHKXER('DGEMQR',INFot,NOUt,LERr,OK)
      INFot = 11
      CALL DGEMQR('L','N',2,1,1,a,2,tau,6,c,0,w,1,info)
      CALL CHKXER('DGEMQR',INFot,NOUt,LERr,OK)
      INFot = 13
      CALL DGEMQR('L','N',2,2,1,a,2,tau,6,c,2,w,0,info)
      CALL CHKXER('DGEMQR',INFot,NOUt,LERr,OK)
!
!     DGELQ
!
      SRNamt = 'DGELQ'
      INFot = 1
      CALL DGELQ(-1,0,a,1,tau,1,w,1,info)
      CALL CHKXER('DGELQ',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DGELQ(0,-1,a,1,tau,1,w,1,info)
      CALL CHKXER('DGELQ',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL DGELQ(1,1,a,0,tau,1,w,1,info)
      CALL CHKXER('DGELQ',INFot,NOUt,LERr,OK)
      INFot = 6
      CALL DGELQ(2,3,a,3,tau,1,w,1,info)
      CALL CHKXER('DGELQ',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL DGELQ(2,3,a,3,tau,7,w,0,info)
      CALL CHKXER('DGELQ',INFot,NOUt,LERr,OK)
!
!     DGEMLQ
!
      tau(1) = 1
      tau(2) = 1
      SRNamt = 'DGEMLQ'
      nb = 1
      INFot = 1
      CALL DGEMLQ('/','N',0,0,0,a,1,tau,1,c,1,w,1,info)
      CALL CHKXER('DGEMLQ',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DGEMLQ('L','/',0,0,0,a,1,tau,1,c,1,w,1,info)
      CALL CHKXER('DGEMLQ',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL DGEMLQ('L','N',-1,0,0,a,1,tau,1,c,1,w,1,info)
      CALL CHKXER('DGEMLQ',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL DGEMLQ('L','N',0,-1,0,a,1,tau,1,c,1,w,1,info)
      CALL CHKXER('DGEMLQ',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL DGEMLQ('L','N',0,0,-1,a,1,tau,1,c,1,w,1,info)
      CALL CHKXER('DGEMLQ',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL DGEMLQ('R','N',0,0,-1,a,1,tau,1,c,1,w,1,info)
      CALL CHKXER('DGEMLQ',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL DGEMLQ('L','N',1,2,0,a,0,tau,1,c,1,w,1,info)
      CALL CHKXER('DGEMLQ',INFot,NOUt,LERr,OK)
      INFot = 9
      CALL DGEMLQ('R','N',2,2,1,a,1,tau,0,c,1,w,1,info)
      CALL CHKXER('DGEMLQ',INFot,NOUt,LERr,OK)
      INFot = 9
      CALL DGEMLQ('L','N',2,2,1,a,1,tau,0,c,1,w,1,info)
      CALL CHKXER('DGEMLQ',INFot,NOUt,LERr,OK)
      INFot = 11
      CALL DGEMLQ('L','N',1,2,1,a,1,tau,6,c,0,w,1,info)
      CALL CHKXER('DGEMLQ',INFot,NOUt,LERr,OK)
      INFot = 13
      CALL DGEMLQ('L','N',2,2,1,a,2,tau,6,c,2,w,0,info)
      CALL CHKXER('DGEMLQ',INFot,NOUt,LERr,OK)
!
!     Print a summary line.
!
      CALL ALAESM(Path,OK,NOUt)
!
!
!     End of DERRTSQR
!
      END SUBROUTINE DERRTSQR
