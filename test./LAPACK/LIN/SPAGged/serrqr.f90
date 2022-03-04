!*==serrqr.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b SERRQR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SERRQR( PATH, NUNIT )
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
!> SERRQR tests the error exits for the REAL routines
!> that use the QR decomposition of a general matrix.
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
      SUBROUTINE SERRQR(Path,Nunit)
      IMPLICIT NONE
!*--SERRQR59
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
      REAL a(NMAX,NMAX) , af(NMAX,NMAX) , b(NMAX) , w(NMAX) , x(NMAX)
!     ..
!     .. External Subroutines ..
      EXTERNAL ALAESM , CHKXER , SGEQR2 , SGEQR2P , SGEQRF , SGEQRFP ,  &
     &         SGEQRS , SORG2R , SORGQR , SORM2R , SORMQR
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
            a(i,j) = 1./REAL(i+j)
            af(i,j) = 1./REAL(i+j)
         ENDDO
         b(j) = 0.
         w(j) = 0.
         x(j) = 0.
      ENDDO
      OK = .TRUE.
!
!     Error exits for QR factorization
!
!     SGEQRF
!
      SRNamt = 'SGEQRF'
      INFot = 1
      CALL SGEQRF(-1,0,a,1,b,w,1,info)
      CALL CHKXER('SGEQRF',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL SGEQRF(0,-1,a,1,b,w,1,info)
      CALL CHKXER('SGEQRF',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL SGEQRF(2,1,a,1,b,w,1,info)
      CALL CHKXER('SGEQRF',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL SGEQRF(1,2,a,1,b,w,1,info)
      CALL CHKXER('SGEQRF',INFot,NOUt,LERr,OK)
!
!     SGEQRFP
!
      SRNamt = 'SGEQRFP'
      INFot = 1
      CALL SGEQRFP(-1,0,a,1,b,w,1,info)
      CALL CHKXER('SGEQRFP',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL SGEQRFP(0,-1,a,1,b,w,1,info)
      CALL CHKXER('SGEQRFP',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL SGEQRFP(2,1,a,1,b,w,1,info)
      CALL CHKXER('SGEQRFP',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL SGEQRFP(1,2,a,1,b,w,1,info)
      CALL CHKXER('SGEQRFP',INFot,NOUt,LERr,OK)
!
!     SGEQR2
!
      SRNamt = 'SGEQR2'
      INFot = 1
      CALL SGEQR2(-1,0,a,1,b,w,info)
      CALL CHKXER('SGEQR2',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL SGEQR2(0,-1,a,1,b,w,info)
      CALL CHKXER('SGEQR2',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL SGEQR2(2,1,a,1,b,w,info)
      CALL CHKXER('SGEQR2',INFot,NOUt,LERr,OK)
!
!     SGEQR2P
!
      SRNamt = 'SGEQR2P'
      INFot = 1
      CALL SGEQR2P(-1,0,a,1,b,w,info)
      CALL CHKXER('SGEQR2P',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL SGEQR2P(0,-1,a,1,b,w,info)
      CALL CHKXER('SGEQR2P',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL SGEQR2P(2,1,a,1,b,w,info)
      CALL CHKXER('SGEQR2P',INFot,NOUt,LERr,OK)
!
!     SGEQRS
!
      SRNamt = 'SGEQRS'
      INFot = 1
      CALL SGEQRS(-1,0,0,a,1,x,b,1,w,1,info)
      CALL CHKXER('SGEQRS',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL SGEQRS(0,-1,0,a,1,x,b,1,w,1,info)
      CALL CHKXER('SGEQRS',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL SGEQRS(1,2,0,a,2,x,b,2,w,1,info)
      CALL CHKXER('SGEQRS',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL SGEQRS(0,0,-1,a,1,x,b,1,w,1,info)
      CALL CHKXER('SGEQRS',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL SGEQRS(2,1,0,a,1,x,b,2,w,1,info)
      CALL CHKXER('SGEQRS',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL SGEQRS(2,1,0,a,2,x,b,1,w,1,info)
      CALL CHKXER('SGEQRS',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL SGEQRS(1,1,2,a,1,x,b,1,w,1,info)
      CALL CHKXER('SGEQRS',INFot,NOUt,LERr,OK)
!
!     SORGQR
!
      SRNamt = 'SORGQR'
      INFot = 1
      CALL SORGQR(-1,0,0,a,1,x,w,1,info)
      CALL CHKXER('SORGQR',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL SORGQR(0,-1,0,a,1,x,w,1,info)
      CALL CHKXER('SORGQR',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL SORGQR(1,2,0,a,1,x,w,2,info)
      CALL CHKXER('SORGQR',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL SORGQR(0,0,-1,a,1,x,w,1,info)
      CALL CHKXER('SORGQR',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL SORGQR(1,1,2,a,1,x,w,1,info)
      CALL CHKXER('SORGQR',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL SORGQR(2,2,0,a,1,x,w,2,info)
      CALL CHKXER('SORGQR',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL SORGQR(2,2,0,a,2,x,w,1,info)
      CALL CHKXER('SORGQR',INFot,NOUt,LERr,OK)
!
!     SORG2R
!
      SRNamt = 'SORG2R'
      INFot = 1
      CALL SORG2R(-1,0,0,a,1,x,w,info)
      CALL CHKXER('SORG2R',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL SORG2R(0,-1,0,a,1,x,w,info)
      CALL CHKXER('SORG2R',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL SORG2R(1,2,0,a,1,x,w,info)
      CALL CHKXER('SORG2R',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL SORG2R(0,0,-1,a,1,x,w,info)
      CALL CHKXER('SORG2R',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL SORG2R(2,1,2,a,2,x,w,info)
      CALL CHKXER('SORG2R',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL SORG2R(2,1,0,a,1,x,w,info)
      CALL CHKXER('SORG2R',INFot,NOUt,LERr,OK)
!
!     SORMQR
!
      SRNamt = 'SORMQR'
      INFot = 1
      CALL SORMQR('/','N',0,0,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('SORMQR',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL SORMQR('L','/',0,0,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('SORMQR',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL SORMQR('L','N',-1,0,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('SORMQR',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL SORMQR('L','N',0,-1,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('SORMQR',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL SORMQR('L','N',0,0,-1,a,1,x,af,1,w,1,info)
      CALL CHKXER('SORMQR',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL SORMQR('L','N',0,1,1,a,1,x,af,1,w,1,info)
      CALL CHKXER('SORMQR',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL SORMQR('R','N',1,0,1,a,1,x,af,1,w,1,info)
      CALL CHKXER('SORMQR',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL SORMQR('L','N',2,1,0,a,1,x,af,2,w,1,info)
      CALL CHKXER('SORMQR',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL SORMQR('R','N',1,2,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('SORMQR',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL SORMQR('L','N',2,1,0,a,2,x,af,1,w,1,info)
      CALL CHKXER('SORMQR',INFot,NOUt,LERr,OK)
      INFot = 12
      CALL SORMQR('L','N',1,2,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('SORMQR',INFot,NOUt,LERr,OK)
      INFot = 12
      CALL SORMQR('R','N',2,1,0,a,1,x,af,2,w,1,info)
      CALL CHKXER('SORMQR',INFot,NOUt,LERr,OK)
!
!     SORM2R
!
      SRNamt = 'SORM2R'
      INFot = 1
      CALL SORM2R('/','N',0,0,0,a,1,x,af,1,w,info)
      CALL CHKXER('SORM2R',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL SORM2R('L','/',0,0,0,a,1,x,af,1,w,info)
      CALL CHKXER('SORM2R',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL SORM2R('L','N',-1,0,0,a,1,x,af,1,w,info)
      CALL CHKXER('SORM2R',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL SORM2R('L','N',0,-1,0,a,1,x,af,1,w,info)
      CALL CHKXER('SORM2R',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL SORM2R('L','N',0,0,-1,a,1,x,af,1,w,info)
      CALL CHKXER('SORM2R',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL SORM2R('L','N',0,1,1,a,1,x,af,1,w,info)
      CALL CHKXER('SORM2R',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL SORM2R('R','N',1,0,1,a,1,x,af,1,w,info)
      CALL CHKXER('SORM2R',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL SORM2R('L','N',2,1,0,a,1,x,af,2,w,info)
      CALL CHKXER('SORM2R',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL SORM2R('R','N',1,2,0,a,1,x,af,1,w,info)
      CALL CHKXER('SORM2R',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL SORM2R('L','N',2,1,0,a,2,x,af,1,w,info)
      CALL CHKXER('SORM2R',INFot,NOUt,LERr,OK)
!
!     Print a summary line.
!
      CALL ALAESM(Path,OK,NOUt)
!
!
!     End of SERRQR
!
      END SUBROUTINE SERRQR
