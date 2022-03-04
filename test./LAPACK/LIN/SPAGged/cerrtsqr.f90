!*==cerrtsqr.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b CERRTSQR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CERRTSQR( PATH, NUNIT )
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
!> CERRTSQR tests the error exits for the COMPLEX routines
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
      SUBROUTINE CERRTSQR(Path,Nunit)
      IMPLICIT NONE
!*--CERRTSQR59
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
      COMPLEX a(NMAX,NMAX) , t(NMAX,NMAX) , w(NMAX) , c(NMAX,NMAX) ,    &
     &        tau(NMAX)
!     ..
!     .. External Subroutines ..
      EXTERNAL ALAESM , CHKXER , CGEQR , CGEMQR , CGELQ , CGEMLQ
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
      INTRINSIC REAL
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
!     Error exits for TS factorization
!
!     CGEQR
!
      SRNamt = 'CGEQR'
      INFot = 1
      CALL CGEQR(-1,0,a,1,tau,1,w,1,info)
      CALL CHKXER('CGEQR',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CGEQR(0,-1,a,1,tau,1,w,1,info)
      CALL CHKXER('CGEQR',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL CGEQR(1,1,a,0,tau,1,w,1,info)
      CALL CHKXER('CGEQR',INFot,NOUt,LERr,OK)
      INFot = 6
      CALL CGEQR(3,2,a,3,tau,1,w,1,info)
      CALL CHKXER('CGEQR',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL CGEQR(3,2,a,3,tau,8,w,0,info)
      CALL CHKXER('CGEQR',INFot,NOUt,LERr,OK)
!
!     CGEMQR
!
      tau(1) = 1
      tau(2) = 1
      SRNamt = 'CGEMQR'
      nb = 1
      INFot = 1
      CALL CGEMQR('/','N',0,0,0,a,1,tau,1,c,1,w,1,info)
      CALL CHKXER('CGEMQR',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CGEMQR('L','/',0,0,0,a,1,tau,1,c,1,w,1,info)
      CALL CHKXER('CGEMQR',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL CGEMQR('L','N',-1,0,0,a,1,tau,1,c,1,w,1,info)
      CALL CHKXER('CGEMQR',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL CGEMQR('L','N',0,-1,0,a,1,tau,1,c,1,w,1,info)
      CALL CHKXER('CGEMQR',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL CGEMQR('L','N',0,0,-1,a,1,tau,1,c,1,w,1,info)
      CALL CHKXER('CGEMQR',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL CGEMQR('R','N',0,0,-1,a,1,tau,1,c,1,w,1,info)
      CALL CHKXER('CGEMQR',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL CGEMQR('L','N',2,1,0,a,0,tau,1,c,1,w,1,info)
      CALL CHKXER('CGEMQR',INFot,NOUt,LERr,OK)
      INFot = 9
      CALL CGEMQR('R','N',2,2,1,a,2,tau,0,c,1,w,1,info)
      CALL CHKXER('CGEMQR',INFot,NOUt,LERr,OK)
      INFot = 9
      CALL CGEMQR('L','N',2,2,1,a,2,tau,0,c,1,w,1,info)
      CALL CHKXER('CGEMQR',INFot,NOUt,LERr,OK)
      INFot = 11
      CALL CGEMQR('L','N',2,1,1,a,2,tau,6,c,0,w,1,info)
      CALL CHKXER('CGEMQR',INFot,NOUt,LERr,OK)
      INFot = 13
      CALL CGEMQR('L','N',2,2,1,a,2,tau,6,c,2,w,0,info)
      CALL CHKXER('CGEMQR',INFot,NOUt,LERr,OK)
!
!     CGELQ
!
      SRNamt = 'CGELQ'
      INFot = 1
      CALL CGELQ(-1,0,a,1,tau,1,w,1,info)
      CALL CHKXER('CGELQ',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CGELQ(0,-1,a,1,tau,1,w,1,info)
      CALL CHKXER('CGELQ',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL CGELQ(1,1,a,0,tau,1,w,1,info)
      CALL CHKXER('CGELQ',INFot,NOUt,LERr,OK)
      INFot = 6
      CALL CGELQ(2,3,a,3,tau,1,w,1,info)
      CALL CHKXER('CGELQ',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL CGELQ(2,3,a,3,tau,8,w,0,info)
      CALL CHKXER('CGELQ',INFot,NOUt,LERr,OK)
!
!     CGEMLQ
!
      tau(1) = 1
      tau(2) = 1
      SRNamt = 'CGEMLQ'
      nb = 1
      INFot = 1
      CALL CGEMLQ('/','N',0,0,0,a,1,tau,1,c,1,w,1,info)
      CALL CHKXER('CGEMLQ',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CGEMLQ('L','/',0,0,0,a,1,tau,1,c,1,w,1,info)
      CALL CHKXER('CGEMLQ',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL CGEMLQ('L','N',-1,0,0,a,1,tau,1,c,1,w,1,info)
      CALL CHKXER('CGEMLQ',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL CGEMLQ('L','N',0,-1,0,a,1,tau,1,c,1,w,1,info)
      CALL CHKXER('CGEMLQ',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL CGEMLQ('L','N',0,0,-1,a,1,tau,1,c,1,w,1,info)
      CALL CHKXER('CGEMLQ',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL CGEMLQ('R','N',0,0,-1,a,1,tau,1,c,1,w,1,info)
      CALL CHKXER('CGEMLQ',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL CGEMLQ('L','N',1,2,0,a,0,tau,1,c,1,w,1,info)
      CALL CHKXER('CGEMLQ',INFot,NOUt,LERr,OK)
      INFot = 9
      CALL CGEMLQ('R','N',2,2,1,a,1,tau,0,c,1,w,1,info)
      CALL CHKXER('CGEMLQ',INFot,NOUt,LERr,OK)
      INFot = 9
      CALL CGEMLQ('L','N',2,2,1,a,1,tau,0,c,1,w,1,info)
      CALL CHKXER('CGEMLQ',INFot,NOUt,LERr,OK)
      INFot = 11
      CALL CGEMLQ('L','N',1,2,1,a,1,tau,6,c,0,w,1,info)
      CALL CHKXER('CGEMLQ',INFot,NOUt,LERr,OK)
      INFot = 13
      CALL CGEMLQ('L','N',2,2,1,a,2,tau,6,c,2,w,0,info)
      CALL CHKXER('CGEMLQ',INFot,NOUt,LERr,OK)
!
!     Print a summary line.
!
      CALL ALAESM(Path,OK,NOUt)
!
!
!     End of CERRTSQR
!
      END SUBROUTINE CERRTSQR
