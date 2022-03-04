!*==zerrps.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b ZERRPS
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZERRPS( PATH, NUNIT )
!
!       .. Scalar Arguments ..
!       INTEGER            NUNIT
!       CHARACTER*3        PATH
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZERRPS tests the error exits for the COMPLEX routines
!> for ZPSTRF.
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
      SUBROUTINE ZERRPS(Path,Nunit)
      IMPLICIT NONE
!*--ZERRPS59
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Nunit
      CHARACTER*3 Path
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER NMAX
      PARAMETER (NMAX=4)
!     ..
!     .. Local Scalars ..
      INTEGER i , info , j , rank
!     ..
!     .. Local Arrays ..
      COMPLEX*16 a(NMAX,NMAX)
      DOUBLE PRECISION rwork(2*NMAX)
      INTEGER piv(NMAX)
!     ..
!     .. External Subroutines ..
      EXTERNAL ALAESM , CHKXER , ZPSTF2 , ZPSTRF
!     ..
!     .. Scalars in Common ..
      INTEGER INFot , NOUt
      LOGICAL LERr , OK
      CHARACTER*32 SRNamt
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
!
         ENDDO
         piv(j) = j
         rwork(j) = 0.D0
         rwork(NMAX+j) = 0.D0
!
      ENDDO
      OK = .TRUE.
!
!
!        Test error exits of the routines that use the Cholesky
!        decomposition of an Hermitian positive semidefinite matrix.
!
!        ZPSTRF
!
      SRNamt = 'ZPSTRF'
      INFot = 1
      CALL ZPSTRF('/',0,a,1,piv,rank,-1.D0,rwork,info)
      CALL CHKXER('ZPSTRF',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZPSTRF('U',-1,a,1,piv,rank,-1.D0,rwork,info)
      CALL CHKXER('ZPSTRF',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL ZPSTRF('U',2,a,1,piv,rank,-1.D0,rwork,info)
      CALL CHKXER('ZPSTRF',INFot,NOUt,LERr,OK)
!
!        ZPSTF2
!
      SRNamt = 'ZPSTF2'
      INFot = 1
      CALL ZPSTF2('/',0,a,1,piv,rank,-1.D0,rwork,info)
      CALL CHKXER('ZPSTF2',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZPSTF2('U',-1,a,1,piv,rank,-1.D0,rwork,info)
      CALL CHKXER('ZPSTF2',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL ZPSTF2('U',2,a,1,piv,rank,-1.D0,rwork,info)
      CALL CHKXER('ZPSTF2',INFot,NOUt,LERr,OK)
!
!
!     Print a summary line.
!
      CALL ALAESM(Path,OK,NOUt)
!
!
!     End of ZERRPS
!
      END SUBROUTINE ZERRPS
