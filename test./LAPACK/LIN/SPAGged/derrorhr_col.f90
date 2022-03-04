!*==derrorhr_col.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b DERRORHR_COL
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DERRORHR_COL( PATH, NUNIT )
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
!> DERRORHR_COL tests the error exits for DORHR_COL that does
!> Householder reconstruction from the output of tall-skinny
!> factorization DLATSQR.
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
!> \date November 2019
!
!> \ingroup double_lin
!
!  =====================================================================
      SUBROUTINE DERRORHR_COL(Path,Nunit)
      IMPLICIT NONE
!*--DERRORHR_COL60
!
!  -- LAPACK test routine (version 3.9.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2019
!
!     .. Scalar Arguments ..
      CHARACTER(LEN=3) Path
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
      DOUBLE PRECISION a(NMAX,NMAX) , t(NMAX,NMAX) , d(NMAX)
!     ..
!     .. External Subroutines ..
      EXTERNAL ALAESM , CHKXER , DORHR_COL
!     ..
!     .. Scalars in Common ..
      LOGICAL LERr , OK
      CHARACTER(LEN=32) SRNamt
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
            a(i,j) = 1.D+0/DBLE(i+j)
            t(i,j) = 1.D+0/DBLE(i+j)
         ENDDO
         d(j) = 0.D+0
      ENDDO
      OK = .TRUE.
!
!     Error exits for Householder reconstruction
!
!     DORHR_COL
!
      SRNamt = 'DORHR_COL'
!
      INFot = 1
      CALL DORHR_COL(-1,0,1,a,1,t,1,d,info)
      CALL CHKXER('DORHR_COL',INFot,NOUt,LERr,OK)
!
      INFot = 2
      CALL DORHR_COL(0,-1,1,a,1,t,1,d,info)
      CALL CHKXER('DORHR_COL',INFot,NOUt,LERr,OK)
      CALL DORHR_COL(1,2,1,a,1,t,1,d,info)
      CALL CHKXER('DORHR_COL',INFot,NOUt,LERr,OK)
!
      INFot = 3
      CALL DORHR_COL(0,0,-1,a,1,t,1,d,info)
      CALL CHKXER('DORHR_COL',INFot,NOUt,LERr,OK)
!
      CALL DORHR_COL(0,0,0,a,1,t,1,d,info)
      CALL CHKXER('DORHR_COL',INFot,NOUt,LERr,OK)
!
      INFot = 5
      CALL DORHR_COL(0,0,1,a,-1,t,1,d,info)
      CALL CHKXER('DORHR_COL',INFot,NOUt,LERr,OK)
!
      CALL DORHR_COL(0,0,1,a,0,t,1,d,info)
      CALL CHKXER('DORHR_COL',INFot,NOUt,LERr,OK)
!
      CALL DORHR_COL(2,0,1,a,1,t,1,d,info)
      CALL CHKXER('DORHR_COL',INFot,NOUt,LERr,OK)
!
      INFot = 7
      CALL DORHR_COL(0,0,1,a,1,t,-1,d,info)
      CALL CHKXER('DORHR_COL',INFot,NOUt,LERr,OK)
!
      CALL DORHR_COL(0,0,1,a,1,t,0,d,info)
      CALL CHKXER('DORHR_COL',INFot,NOUt,LERr,OK)
!
      CALL DORHR_COL(4,3,2,a,4,t,1,d,info)
      CALL CHKXER('DORHR_COL',INFot,NOUt,LERr,OK)
!
!     Print a summary line.
!
      CALL ALAESM(Path,OK,NOUt)
!
!
!     End of DERRORHR_COL
!
      END SUBROUTINE DERRORHR_COL
