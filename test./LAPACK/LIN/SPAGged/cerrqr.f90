!*==cerrqr.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b CERRQR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CERRQR( PATH, NUNIT )
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
!> CERRQR tests the error exits for the COMPLEX routines
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
!> \ingroup complex_lin
!
!  =====================================================================
      SUBROUTINE CERRQR(Path,Nunit)
      IMPLICIT NONE
!*--CERRQR59
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
      COMPLEX a(NMAX,NMAX) , af(NMAX,NMAX) , b(NMAX) , w(NMAX) , x(NMAX)
!     ..
!     .. External Subroutines ..
      EXTERNAL ALAESM , CGEQR2 , CGEQR2P , CGEQRF , CGEQRFP , CGEQRS ,  &
     &         CHKXER , CUNG2R , CUNGQR , CUNM2R , CUNMQR
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
      INTRINSIC CMPLX , REAL
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
            a(i,j) = CMPLX(1./REAL(i+j),-1./REAL(i+j))
            af(i,j) = CMPLX(1./REAL(i+j),-1./REAL(i+j))
         ENDDO
         b(j) = 0.
         w(j) = 0.
         x(j) = 0.
      ENDDO
      OK = .TRUE.
!
!     Error exits for QR factorization
!
!     CGEQRF
!
      SRNamt = 'CGEQRF'
      INFot = 1
      CALL CGEQRF(-1,0,a,1,b,w,1,info)
      CALL CHKXER('CGEQRF',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CGEQRF(0,-1,a,1,b,w,1,info)
      CALL CHKXER('CGEQRF',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL CGEQRF(2,1,a,1,b,w,1,info)
      CALL CHKXER('CGEQRF',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL CGEQRF(1,2,a,1,b,w,1,info)
      CALL CHKXER('CGEQRF',INFot,NOUt,LERr,OK)
!
!     CGEQRFP
!
      SRNamt = 'CGEQRFP'
      INFot = 1
      CALL CGEQRFP(-1,0,a,1,b,w,1,info)
      CALL CHKXER('CGEQRFP',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CGEQRFP(0,-1,a,1,b,w,1,info)
      CALL CHKXER('CGEQRFP',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL CGEQRFP(2,1,a,1,b,w,1,info)
      CALL CHKXER('CGEQRFP',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL CGEQRFP(1,2,a,1,b,w,1,info)
      CALL CHKXER('CGEQRFP',INFot,NOUt,LERr,OK)
!
!     CGEQR2
!
      SRNamt = 'CGEQR2'
      INFot = 1
      CALL CGEQR2(-1,0,a,1,b,w,info)
      CALL CHKXER('CGEQR2',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CGEQR2(0,-1,a,1,b,w,info)
      CALL CHKXER('CGEQR2',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL CGEQR2(2,1,a,1,b,w,info)
      CALL CHKXER('CGEQR2',INFot,NOUt,LERr,OK)
!
!     CGEQR2P
!
      SRNamt = 'CGEQR2P'
      INFot = 1
      CALL CGEQR2P(-1,0,a,1,b,w,info)
      CALL CHKXER('CGEQR2P',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CGEQR2P(0,-1,a,1,b,w,info)
      CALL CHKXER('CGEQR2P',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL CGEQR2P(2,1,a,1,b,w,info)
      CALL CHKXER('CGEQR2P',INFot,NOUt,LERr,OK)
!
!     CGEQRS
!
      SRNamt = 'CGEQRS'
      INFot = 1
      CALL CGEQRS(-1,0,0,a,1,x,b,1,w,1,info)
      CALL CHKXER('CGEQRS',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CGEQRS(0,-1,0,a,1,x,b,1,w,1,info)
      CALL CHKXER('CGEQRS',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CGEQRS(1,2,0,a,2,x,b,2,w,1,info)
      CALL CHKXER('CGEQRS',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL CGEQRS(0,0,-1,a,1,x,b,1,w,1,info)
      CALL CHKXER('CGEQRS',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL CGEQRS(2,1,0,a,1,x,b,2,w,1,info)
      CALL CHKXER('CGEQRS',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL CGEQRS(2,1,0,a,2,x,b,1,w,1,info)
      CALL CHKXER('CGEQRS',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL CGEQRS(1,1,2,a,1,x,b,1,w,1,info)
      CALL CHKXER('CGEQRS',INFot,NOUt,LERr,OK)
!
!     CUNGQR
!
      SRNamt = 'CUNGQR'
      INFot = 1
      CALL CUNGQR(-1,0,0,a,1,x,w,1,info)
      CALL CHKXER('CUNGQR',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CUNGQR(0,-1,0,a,1,x,w,1,info)
      CALL CHKXER('CUNGQR',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CUNGQR(1,2,0,a,1,x,w,2,info)
      CALL CHKXER('CUNGQR',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL CUNGQR(0,0,-1,a,1,x,w,1,info)
      CALL CHKXER('CUNGQR',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL CUNGQR(1,1,2,a,1,x,w,1,info)
      CALL CHKXER('CUNGQR',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL CUNGQR(2,2,0,a,1,x,w,2,info)
      CALL CHKXER('CUNGQR',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL CUNGQR(2,2,0,a,2,x,w,1,info)
      CALL CHKXER('CUNGQR',INFot,NOUt,LERr,OK)
!
!     CUNG2R
!
      SRNamt = 'CUNG2R'
      INFot = 1
      CALL CUNG2R(-1,0,0,a,1,x,w,info)
      CALL CHKXER('CUNG2R',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CUNG2R(0,-1,0,a,1,x,w,info)
      CALL CHKXER('CUNG2R',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CUNG2R(1,2,0,a,1,x,w,info)
      CALL CHKXER('CUNG2R',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL CUNG2R(0,0,-1,a,1,x,w,info)
      CALL CHKXER('CUNG2R',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL CUNG2R(2,1,2,a,2,x,w,info)
      CALL CHKXER('CUNG2R',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL CUNG2R(2,1,0,a,1,x,w,info)
      CALL CHKXER('CUNG2R',INFot,NOUt,LERr,OK)
!
!     CUNMQR
!
      SRNamt = 'CUNMQR'
      INFot = 1
      CALL CUNMQR('/','N',0,0,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('CUNMQR',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CUNMQR('L','/',0,0,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('CUNMQR',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL CUNMQR('L','N',-1,0,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('CUNMQR',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL CUNMQR('L','N',0,-1,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('CUNMQR',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL CUNMQR('L','N',0,0,-1,a,1,x,af,1,w,1,info)
      CALL CHKXER('CUNMQR',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL CUNMQR('L','N',0,1,1,a,1,x,af,1,w,1,info)
      CALL CHKXER('CUNMQR',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL CUNMQR('R','N',1,0,1,a,1,x,af,1,w,1,info)
      CALL CHKXER('CUNMQR',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL CUNMQR('L','N',2,1,0,a,1,x,af,2,w,1,info)
      CALL CHKXER('CUNMQR',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL CUNMQR('R','N',1,2,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('CUNMQR',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL CUNMQR('L','N',2,1,0,a,2,x,af,1,w,1,info)
      CALL CHKXER('CUNMQR',INFot,NOUt,LERr,OK)
      INFot = 12
      CALL CUNMQR('L','N',1,2,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('CUNMQR',INFot,NOUt,LERr,OK)
      INFot = 12
      CALL CUNMQR('R','N',2,1,0,a,1,x,af,2,w,1,info)
      CALL CHKXER('CUNMQR',INFot,NOUt,LERr,OK)
!
!     CUNM2R
!
      SRNamt = 'CUNM2R'
      INFot = 1
      CALL CUNM2R('/','N',0,0,0,a,1,x,af,1,w,info)
      CALL CHKXER('CUNM2R',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CUNM2R('L','/',0,0,0,a,1,x,af,1,w,info)
      CALL CHKXER('CUNM2R',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL CUNM2R('L','N',-1,0,0,a,1,x,af,1,w,info)
      CALL CHKXER('CUNM2R',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL CUNM2R('L','N',0,-1,0,a,1,x,af,1,w,info)
      CALL CHKXER('CUNM2R',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL CUNM2R('L','N',0,0,-1,a,1,x,af,1,w,info)
      CALL CHKXER('CUNM2R',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL CUNM2R('L','N',0,1,1,a,1,x,af,1,w,info)
      CALL CHKXER('CUNM2R',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL CUNM2R('R','N',1,0,1,a,1,x,af,1,w,info)
      CALL CHKXER('CUNM2R',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL CUNM2R('L','N',2,1,0,a,1,x,af,2,w,info)
      CALL CHKXER('CUNM2R',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL CUNM2R('R','N',1,2,0,a,1,x,af,1,w,info)
      CALL CHKXER('CUNM2R',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL CUNM2R('L','N',2,1,0,a,2,x,af,1,w,info)
      CALL CHKXER('CUNM2R',INFot,NOUt,LERr,OK)
!
!     Print a summary line.
!
      CALL ALAESM(Path,OK,NOUt)
!
!
!     End of CERRQR
!
      END SUBROUTINE CERRQR
