!*==derrps.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b derrps
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE derrps( PATH, NUNIT )
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
!> derrps tests the error exits for the DOUBLE PRECISION routines
!> for DPSTRF.
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
      SUBROUTINE derrps(Path,Nunit)
      IMPLICIT NONE
!*--DERRPS59
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
      DOUBLE PRECISION a(NMAX,NMAX) , work(2*NMAX)
      INTEGER piv(NMAX)
!     ..
!     .. External Subroutines ..
EXTERNAL    DPSTF2,   DPSTRF
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
         work(j) = 0.D0
         work(NMAX+j) = 0.D0
!
      ENDDO
      OK = .TRUE.
!
!
!        Test error exits of the routines that use the Cholesky
!        decomposition of a symmetric positive semidefinite matrix.
!
!        DPSTRF
!
      SRNamt = 'DPSTRF'
      INFot = 1
      CALL DPSTRF('/',0,a,1,piv,rank,-1.D0,work,info)
      CALL CHKXER('DPSTRF',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DPSTRF('U',-1,a,1,piv,rank,-1.D0,work,info)
      CALL CHKXER('DPSTRF',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL DPSTRF('U',2,a,1,piv,rank,-1.D0,work,info)
      CALL CHKXER('DPSTRF',INFot,NOUt,LERr,OK)
!
!        DPSTF2
!
      SRNamt = 'DPSTF2'
      INFot = 1
      CALL DPSTF2('/',0,a,1,piv,rank,-1.D0,work,info)
      CALL CHKXER('DPSTF2',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DPSTF2('U',-1,a,1,piv,rank,-1.D0,work,info)
      CALL CHKXER('DPSTF2',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL DPSTF2('U',2,a,1,piv,rank,-1.D0,work,info)
      CALL CHKXER('DPSTF2',INFot,NOUt,LERr,OK)
!
!
!     Print a summary line.
!
      CALL alaesm(Path,OK,NOUt)
!
!
!     End of derrps
!
      END SUBROUTINE derrps
