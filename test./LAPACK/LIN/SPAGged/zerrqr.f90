!*==zerrqr.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b ZERRQR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZERRQR( PATH, NUNIT )
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
!> ZERRQR tests the error exits for the COMPLEX*16 routines
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
!> \ingroup complex16_lin
!
!  =====================================================================
      SUBROUTINE ZERRQR(Path,Nunit)
      IMPLICIT NONE
!*--ZERRQR59
!
!  -- LAPACK test routine ((version 3.7.0) --
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
      COMPLEX*16 a(NMAX,NMAX) , af(NMAX,NMAX) , b(NMAX) , w(NMAX) ,     &
     &           x(NMAX)
!     ..
!     .. External Subroutines ..
      EXTERNAL ALAESM , CHKXER , ZGEQR2 , ZGEQR2P , ZGEQRF , ZGEQRFP ,  &
     &         ZGEQRS , ZUNG2R , ZUNGQR , ZUNM2R , ZUNMQR
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
      INTRINSIC DBLE , DCMPLX
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
            a(i,j) = DCMPLX(1.D0/DBLE(i+j),-1.D0/DBLE(i+j))
            af(i,j) = DCMPLX(1.D0/DBLE(i+j),-1.D0/DBLE(i+j))
         ENDDO
         b(j) = 0.D0
         w(j) = 0.D0
         x(j) = 0.D0
      ENDDO
      OK = .TRUE.
!
!     Error exits for QR factorization
!
!     ZGEQRF
!
      SRNamt = 'ZGEQRF'
      INFot = 1
      CALL ZGEQRF(-1,0,a,1,b,w,1,info)
      CALL CHKXER('ZGEQRF',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZGEQRF(0,-1,a,1,b,w,1,info)
      CALL CHKXER('ZGEQRF',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL ZGEQRF(2,1,a,1,b,w,1,info)
      CALL CHKXER('ZGEQRF',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL ZGEQRF(1,2,a,1,b,w,1,info)
      CALL CHKXER('ZGEQRF',INFot,NOUt,LERr,OK)
!
!     ZGEQRFP
!
      SRNamt = 'ZGEQRFP'
      INFot = 1
      CALL ZGEQRFP(-1,0,a,1,b,w,1,info)
      CALL CHKXER('ZGEQRFP',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZGEQRFP(0,-1,a,1,b,w,1,info)
      CALL CHKXER('ZGEQRFP',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL ZGEQRFP(2,1,a,1,b,w,1,info)
      CALL CHKXER('ZGEQRFP',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL ZGEQRFP(1,2,a,1,b,w,1,info)
      CALL CHKXER('ZGEQRFP',INFot,NOUt,LERr,OK)
!
!     ZGEQR2
!
      SRNamt = 'ZGEQR2'
      INFot = 1
      CALL ZGEQR2(-1,0,a,1,b,w,info)
      CALL CHKXER('ZGEQR2',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZGEQR2(0,-1,a,1,b,w,info)
      CALL CHKXER('ZGEQR2',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL ZGEQR2(2,1,a,1,b,w,info)
      CALL CHKXER('ZGEQR2',INFot,NOUt,LERr,OK)
!
!     ZGEQR2P
!
      SRNamt = 'ZGEQR2P'
      INFot = 1
      CALL ZGEQR2P(-1,0,a,1,b,w,info)
      CALL CHKXER('ZGEQR2P',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZGEQR2P(0,-1,a,1,b,w,info)
      CALL CHKXER('ZGEQR2P',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL ZGEQR2P(2,1,a,1,b,w,info)
      CALL CHKXER('ZGEQR2P',INFot,NOUt,LERr,OK)
!
!     ZGEQRS
!
      SRNamt = 'ZGEQRS'
      INFot = 1
      CALL ZGEQRS(-1,0,0,a,1,x,b,1,w,1,info)
      CALL CHKXER('ZGEQRS',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZGEQRS(0,-1,0,a,1,x,b,1,w,1,info)
      CALL CHKXER('ZGEQRS',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZGEQRS(1,2,0,a,2,x,b,2,w,1,info)
      CALL CHKXER('ZGEQRS',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL ZGEQRS(0,0,-1,a,1,x,b,1,w,1,info)
      CALL CHKXER('ZGEQRS',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL ZGEQRS(2,1,0,a,1,x,b,2,w,1,info)
      CALL CHKXER('ZGEQRS',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL ZGEQRS(2,1,0,a,2,x,b,1,w,1,info)
      CALL CHKXER('ZGEQRS',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL ZGEQRS(1,1,2,a,1,x,b,1,w,1,info)
      CALL CHKXER('ZGEQRS',INFot,NOUt,LERr,OK)
!
!     ZUNGQR
!
      SRNamt = 'ZUNGQR'
      INFot = 1
      CALL ZUNGQR(-1,0,0,a,1,x,w,1,info)
      CALL CHKXER('ZUNGQR',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZUNGQR(0,-1,0,a,1,x,w,1,info)
      CALL CHKXER('ZUNGQR',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZUNGQR(1,2,0,a,1,x,w,2,info)
      CALL CHKXER('ZUNGQR',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL ZUNGQR(0,0,-1,a,1,x,w,1,info)
      CALL CHKXER('ZUNGQR',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL ZUNGQR(1,1,2,a,1,x,w,1,info)
      CALL CHKXER('ZUNGQR',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL ZUNGQR(2,2,0,a,1,x,w,2,info)
      CALL CHKXER('ZUNGQR',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL ZUNGQR(2,2,0,a,2,x,w,1,info)
      CALL CHKXER('ZUNGQR',INFot,NOUt,LERr,OK)
!
!     ZUNG2R
!
      SRNamt = 'ZUNG2R'
      INFot = 1
      CALL ZUNG2R(-1,0,0,a,1,x,w,info)
      CALL CHKXER('ZUNG2R',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZUNG2R(0,-1,0,a,1,x,w,info)
      CALL CHKXER('ZUNG2R',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZUNG2R(1,2,0,a,1,x,w,info)
      CALL CHKXER('ZUNG2R',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL ZUNG2R(0,0,-1,a,1,x,w,info)
      CALL CHKXER('ZUNG2R',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL ZUNG2R(2,1,2,a,2,x,w,info)
      CALL CHKXER('ZUNG2R',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL ZUNG2R(2,1,0,a,1,x,w,info)
      CALL CHKXER('ZUNG2R',INFot,NOUt,LERr,OK)
!
!     ZUNMQR
!
      SRNamt = 'ZUNMQR'
      INFot = 1
      CALL ZUNMQR('/','N',0,0,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('ZUNMQR',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZUNMQR('L','/',0,0,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('ZUNMQR',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL ZUNMQR('L','N',-1,0,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('ZUNMQR',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL ZUNMQR('L','N',0,-1,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('ZUNMQR',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL ZUNMQR('L','N',0,0,-1,a,1,x,af,1,w,1,info)
      CALL CHKXER('ZUNMQR',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL ZUNMQR('L','N',0,1,1,a,1,x,af,1,w,1,info)
      CALL CHKXER('ZUNMQR',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL ZUNMQR('R','N',1,0,1,a,1,x,af,1,w,1,info)
      CALL CHKXER('ZUNMQR',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL ZUNMQR('L','N',2,1,0,a,1,x,af,2,w,1,info)
      CALL CHKXER('ZUNMQR',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL ZUNMQR('R','N',1,2,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('ZUNMQR',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL ZUNMQR('L','N',2,1,0,a,2,x,af,1,w,1,info)
      CALL CHKXER('ZUNMQR',INFot,NOUt,LERr,OK)
      INFot = 12
      CALL ZUNMQR('L','N',1,2,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('ZUNMQR',INFot,NOUt,LERr,OK)
      INFot = 12
      CALL ZUNMQR('R','N',2,1,0,a,1,x,af,2,w,1,info)
      CALL CHKXER('ZUNMQR',INFot,NOUt,LERr,OK)
!
!     ZUNM2R
!
      SRNamt = 'ZUNM2R'
      INFot = 1
      CALL ZUNM2R('/','N',0,0,0,a,1,x,af,1,w,info)
      CALL CHKXER('ZUNM2R',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZUNM2R('L','/',0,0,0,a,1,x,af,1,w,info)
      CALL CHKXER('ZUNM2R',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL ZUNM2R('L','N',-1,0,0,a,1,x,af,1,w,info)
      CALL CHKXER('ZUNM2R',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL ZUNM2R('L','N',0,-1,0,a,1,x,af,1,w,info)
      CALL CHKXER('ZUNM2R',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL ZUNM2R('L','N',0,0,-1,a,1,x,af,1,w,info)
      CALL CHKXER('ZUNM2R',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL ZUNM2R('L','N',0,1,1,a,1,x,af,1,w,info)
      CALL CHKXER('ZUNM2R',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL ZUNM2R('R','N',1,0,1,a,1,x,af,1,w,info)
      CALL CHKXER('ZUNM2R',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL ZUNM2R('L','N',2,1,0,a,1,x,af,2,w,info)
      CALL CHKXER('ZUNM2R',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL ZUNM2R('R','N',1,2,0,a,1,x,af,1,w,info)
      CALL CHKXER('ZUNM2R',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL ZUNM2R('L','N',2,1,0,a,2,x,af,1,w,info)
      CALL CHKXER('ZUNM2R',INFot,NOUt,LERr,OK)
!
!     Print a summary line.
!
      CALL ALAESM(Path,OK,NOUt)
!
!
!     End of ZERRQR
!
      END SUBROUTINE ZERRQR
