!*==derrqr.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b derrqr
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE derrqr( PATH, NUNIT )
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
!> derrqr tests the error exits for the DOUBLE PRECISION routines
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
!> \ingroup double_lin
!
!  =====================================================================
      SUBROUTINE derrqr(Path,Nunit)
      IMPLICIT NONE
!*--DERRQR59
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
      DOUBLE PRECISION a(NMAX,NMAX) , af(NMAX,NMAX) , b(NMAX) , w(NMAX) &
     &                 , x(NMAX)
!     ..
!     .. External Subroutines ..
EXTERNAL    DGEQR2,   DGEQR2P,   DGEQRF,   DGEQRFP
EXTERNAL    DORG2R,   DORGQR,   DORM2R,   DORMQR
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
            af(i,j) = 1.D0/DBLE(i+j)
         ENDDO
         b(j) = 0.D0
         w(j) = 0.D0
         x(j) = 0.D0
      ENDDO
      OK = .TRUE.
!
!     Error exits for QR factorization
!
!     DGEQRF
!
      SRNamt = 'DGEQRF'
      INFot = 1
      CALL DGEQRF(-1,0,a,1,b,w,1,info)
      CALL CHKXER('DGEQRF',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DGEQRF(0,-1,a,1,b,w,1,info)
      CALL CHKXER('DGEQRF',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL DGEQRF(2,1,a,1,b,w,1,info)
      CALL CHKXER('DGEQRF',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL DGEQRF(1,2,a,1,b,w,1,info)
      CALL CHKXER('DGEQRF',INFot,NOUt,LERr,OK)
!
!     DGEQRFP
!
      SRNamt = 'DGEQRFP'
      INFot = 1
      CALL DGEQRFP(-1,0,a,1,b,w,1,info)
      CALL CHKXER('DGEQRFP',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DGEQRFP(0,-1,a,1,b,w,1,info)
      CALL CHKXER('DGEQRFP',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL DGEQRFP(2,1,a,1,b,w,1,info)
      CALL CHKXER('DGEQRFP',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL DGEQRFP(1,2,a,1,b,w,1,info)
      CALL CHKXER('DGEQRFP',INFot,NOUt,LERr,OK)
!
!     DGEQR2
!
      SRNamt = 'DGEQR2'
      INFot = 1
      CALL DGEQR2(-1,0,a,1,b,w,info)
      CALL CHKXER('DGEQR2',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DGEQR2(0,-1,a,1,b,w,info)
      CALL CHKXER('DGEQR2',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL DGEQR2(2,1,a,1,b,w,info)
      CALL CHKXER('DGEQR2',INFot,NOUt,LERr,OK)
!
!     DGEQR2P
!
      SRNamt = 'DGEQR2P'
      INFot = 1
      CALL DGEQR2P(-1,0,a,1,b,w,info)
      CALL CHKXER('DGEQR2P',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DGEQR2P(0,-1,a,1,b,w,info)
      CALL CHKXER('DGEQR2P',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL DGEQR2P(2,1,a,1,b,w,info)
      CALL CHKXER('DGEQR2P',INFot,NOUt,LERr,OK)
!
!     dgeqrs
!
      SRNamt = 'dgeqrs'
      INFot = 1
      CALL dgeqrs(-1,0,0,a,1,x,b,1,w,1,info)
      CALL CHKXER('dgeqrs',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL dgeqrs(0,-1,0,a,1,x,b,1,w,1,info)
      CALL CHKXER('dgeqrs',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL dgeqrs(1,2,0,a,2,x,b,2,w,1,info)
      CALL CHKXER('dgeqrs',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL dgeqrs(0,0,-1,a,1,x,b,1,w,1,info)
      CALL CHKXER('dgeqrs',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL dgeqrs(2,1,0,a,1,x,b,2,w,1,info)
      CALL CHKXER('dgeqrs',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL dgeqrs(2,1,0,a,2,x,b,1,w,1,info)
      CALL CHKXER('dgeqrs',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL dgeqrs(1,1,2,a,1,x,b,1,w,1,info)
      CALL CHKXER('dgeqrs',INFot,NOUt,LERr,OK)
!
!     DORGQR
!
      SRNamt = 'DORGQR'
      INFot = 1
      CALL DORGQR(-1,0,0,a,1,x,w,1,info)
      CALL CHKXER('DORGQR',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DORGQR(0,-1,0,a,1,x,w,1,info)
      CALL CHKXER('DORGQR',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DORGQR(1,2,0,a,1,x,w,2,info)
      CALL CHKXER('DORGQR',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL DORGQR(0,0,-1,a,1,x,w,1,info)
      CALL CHKXER('DORGQR',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL DORGQR(1,1,2,a,1,x,w,1,info)
      CALL CHKXER('DORGQR',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL DORGQR(2,2,0,a,1,x,w,2,info)
      CALL CHKXER('DORGQR',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL DORGQR(2,2,0,a,2,x,w,1,info)
      CALL CHKXER('DORGQR',INFot,NOUt,LERr,OK)
!
!     DORG2R
!
      SRNamt = 'DORG2R'
      INFot = 1
      CALL DORG2R(-1,0,0,a,1,x,w,info)
      CALL CHKXER('DORG2R',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DORG2R(0,-1,0,a,1,x,w,info)
      CALL CHKXER('DORG2R',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DORG2R(1,2,0,a,1,x,w,info)
      CALL CHKXER('DORG2R',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL DORG2R(0,0,-1,a,1,x,w,info)
      CALL CHKXER('DORG2R',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL DORG2R(2,1,2,a,2,x,w,info)
      CALL CHKXER('DORG2R',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL DORG2R(2,1,0,a,1,x,w,info)
      CALL CHKXER('DORG2R',INFot,NOUt,LERr,OK)
!
!     DORMQR
!
      SRNamt = 'DORMQR'
      INFot = 1
      CALL DORMQR('/','N',0,0,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('DORMQR',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DORMQR('L','/',0,0,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('DORMQR',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL DORMQR('L','N',-1,0,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('DORMQR',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL DORMQR('L','N',0,-1,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('DORMQR',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL DORMQR('L','N',0,0,-1,a,1,x,af,1,w,1,info)
      CALL CHKXER('DORMQR',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL DORMQR('L','N',0,1,1,a,1,x,af,1,w,1,info)
      CALL CHKXER('DORMQR',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL DORMQR('R','N',1,0,1,a,1,x,af,1,w,1,info)
      CALL CHKXER('DORMQR',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL DORMQR('L','N',2,1,0,a,1,x,af,2,w,1,info)
      CALL CHKXER('DORMQR',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL DORMQR('R','N',1,2,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('DORMQR',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL DORMQR('L','N',2,1,0,a,2,x,af,1,w,1,info)
      CALL CHKXER('DORMQR',INFot,NOUt,LERr,OK)
      INFot = 12
      CALL DORMQR('L','N',1,2,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('DORMQR',INFot,NOUt,LERr,OK)
      INFot = 12
      CALL DORMQR('R','N',2,1,0,a,1,x,af,2,w,1,info)
      CALL CHKXER('DORMQR',INFot,NOUt,LERr,OK)
!
!     DORM2R
!
      SRNamt = 'DORM2R'
      INFot = 1
      CALL DORM2R('/','N',0,0,0,a,1,x,af,1,w,info)
      CALL CHKXER('DORM2R',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DORM2R('L','/',0,0,0,a,1,x,af,1,w,info)
      CALL CHKXER('DORM2R',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL DORM2R('L','N',-1,0,0,a,1,x,af,1,w,info)
      CALL CHKXER('DORM2R',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL DORM2R('L','N',0,-1,0,a,1,x,af,1,w,info)
      CALL CHKXER('DORM2R',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL DORM2R('L','N',0,0,-1,a,1,x,af,1,w,info)
      CALL CHKXER('DORM2R',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL DORM2R('L','N',0,1,1,a,1,x,af,1,w,info)
      CALL CHKXER('DORM2R',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL DORM2R('R','N',1,0,1,a,1,x,af,1,w,info)
      CALL CHKXER('DORM2R',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL DORM2R('L','N',2,1,0,a,1,x,af,2,w,info)
      CALL CHKXER('DORM2R',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL DORM2R('R','N',1,2,0,a,1,x,af,1,w,info)
      CALL CHKXER('DORM2R',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL DORM2R('L','N',2,1,0,a,2,x,af,1,w,info)
      CALL CHKXER('DORM2R',INFot,NOUt,LERr,OK)
!
!     Print a summary line.
!
      CALL alaesm(Path,OK,NOUt)
!
!
!     End of derrqr
!
      END SUBROUTINE derrqr
