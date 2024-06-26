!*==cerrps.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b cerrps
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE cerrps( PATH, NUNIT )
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
!> cerrps tests the error exits for the COMPLEX routines
!> for CPSTRF..
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
      SUBROUTINE cerrps(Path,Nunit)
      IMPLICIT NONE
!*--CERRPS59
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
      COMPLEX a(NMAX,NMAX)
      REAL rwork(2*NMAX)
      INTEGER piv(NMAX)
!     ..
!     .. External Subroutines ..
EXTERNAL    CPSTF2,   CPSTRF
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
            a(i,j) = 1.0/REAL(i+j)
!
         ENDDO
         piv(j) = j
         rwork(j) = 0.
         rwork(NMAX+j) = 0.
!
      ENDDO
      OK = .TRUE.
!
!
!        Test error exits of the routines that use the Cholesky
!        decomposition of an Hermitian positive semidefinite matrix.
!
!        CPSTRF
!
      SRNamt = 'CPSTRF'
      INFot = 1
      CALL CPSTRF('/',0,a,1,piv,rank,-1.0,rwork,info)
      CALL CHKXER('CPSTRF',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CPSTRF('U',-1,a,1,piv,rank,-1.0,rwork,info)
      CALL CHKXER('CPSTRF',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL CPSTRF('U',2,a,1,piv,rank,-1.0,rwork,info)
      CALL CHKXER('CPSTRF',INFot,NOUt,LERr,OK)
!
!        CPSTF2
!
      SRNamt = 'CPSTF2'
      INFot = 1
      CALL CPSTF2('/',0,a,1,piv,rank,-1.0,rwork,info)
      CALL CHKXER('CPSTF2',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CPSTF2('U',-1,a,1,piv,rank,-1.0,rwork,info)
      CALL CHKXER('CPSTF2',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL CPSTF2('U',2,a,1,piv,rank,-1.0,rwork,info)
      CALL CHKXER('CPSTF2',INFot,NOUt,LERr,OK)
!
!
!     Print a summary line.
!
      CALL alaesm(Path,OK,NOUt)
!
!
!     End of cerrps
!
      END SUBROUTINE cerrps
