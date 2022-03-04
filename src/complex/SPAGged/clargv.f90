!*==clargv.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CLARGV generates a vector of plane rotations with real cosines and complex sines.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLARGV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clargv.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clargv.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clargv.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLARGV( N, X, INCX, Y, INCY, C, INCC )
!
!       .. Scalar Arguments ..
!       INTEGER            INCC, INCX, INCY, N
!       ..
!       .. Array Arguments ..
!       REAL               C( * )
!       COMPLEX            X( * ), Y( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLARGV generates a vector of complex plane rotations with real
!> cosines, determined by elements of the complex vectors x and y.
!> For i = 1,2,...,n
!>
!>    (        c(i)   s(i) ) ( x(i) ) = ( r(i) )
!>    ( -conjg(s(i))  c(i) ) ( y(i) ) = (   0  )
!>
!>    where c(i)**2 + ABS(s(i))**2 = 1
!>
!> The following conventions are used (these are the same as in CLARTG,
!> but differ from the BLAS1 routine CROTG):
!>    If y(i)=0, then c(i)=1 and s(i)=0.
!>    If x(i)=0, then c(i)=0 and s(i) is chosen so that r(i) is real.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of plane rotations to be generated.
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is COMPLEX array, dimension (1+(N-1)*INCX)
!>          On entry, the vector x.
!>          On exit, x(i) is overwritten by r(i), for i = 1,...,n.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>          The increment between elements of X. INCX > 0.
!> \endverbatim
!>
!> \param[in,out] Y
!> \verbatim
!>          Y is COMPLEX array, dimension (1+(N-1)*INCY)
!>          On entry, the vector y.
!>          On exit, the sines of the plane rotations.
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>          The increment between elements of Y. INCY > 0.
!> \endverbatim
!>
!> \param[out] C
!> \verbatim
!>          C is REAL array, dimension (1+(N-1)*INCC)
!>          The cosines of the plane rotations.
!> \endverbatim
!>
!> \param[in] INCC
!> \verbatim
!>          INCC is INTEGER
!>          The increment between elements of C. INCC > 0.
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
!> \ingroup complexOTHERauxiliary
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  6-6-96 - Modified with a new algorithm by W. Kahan and J. Demmel
!>
!>  This version has a few statements commented out for thread safety
!>  (machine parameters are computed on each entry). 10 feb 03, SJH.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE CLARGV(N,X,Incx,Y,Incy,C,Incc)
      USE S_SLAMCH
      USE S_SLAPY2
      IMPLICIT NONE
!*--CLARGV128
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  TWO = 2.0E+0 , ONE = 1.0E+0 , ZERO = 0.0E+0
      COMPLEX , PARAMETER  ::  CZERO = (0.0E+0,0.0E+0)
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      REAL , INTENT(OUT) , DIMENSION(*) :: C
      INTEGER , INTENT(IN) :: Incc
!
! Local variable declarations rewritten by SPAG
!
      REAL :: ABS1 , ABSSQ
      INTEGER :: count , i , ic , ix , iy , j
      REAL :: cs , d , di , dr , eps , f2 , f2s , g2 , g2s , safmin ,   &
     &        safmn2 , safmx2 , scale
      COMPLEX :: f , ff , fs , g , gs , r , sn
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     LOGICAL            FIRST
!     ..
!     .. External Functions ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Statement Functions ..
!     ..
!     .. Save statement ..
!     SAVE               FIRST, SAFMX2, SAFMIN, SAFMN2
!     ..
!     .. Data statements ..
!     DATA               FIRST / .TRUE. /
!     ..
!     .. Statement Function definitions ..
      ABS1(ff) = MAX(ABS(REAL(ff)),ABS(AIMAG(ff)))
      ABSSQ(ff) = REAL(ff)**2 + AIMAG(ff)**2
!     ..
!     .. Executable Statements ..
!
!     IF( FIRST ) THEN
!        FIRST = .FALSE.
      safmin = SLAMCH('S')
      eps = SLAMCH('E')
      safmn2 = SLAMCH('B')**INT(LOG(safmin/eps)/LOG(SLAMCH('B'))/TWO)
      safmx2 = ONE/safmn2
!     END IF
      ix = 1
      iy = 1
      ic = 1
      DO i = 1 , N
         f = X(ix)
         g = Y(iy)
!
!        Use identical algorithm as in CLARTG
!
         scale = MAX(ABS1(f),ABS1(g))
         fs = f
         gs = g
         count = 0
         IF ( scale>=safmx2 ) THEN
            DO
               count = count + 1
               fs = fs*safmn2
               gs = gs*safmn2
               scale = scale*safmn2
               IF ( scale<safmx2 .OR. count>=20 ) EXIT
            ENDDO
         ELSEIF ( scale<=safmn2 ) THEN
            IF ( g==CZERO ) THEN
               cs = ONE
               sn = CZERO
               r = f
               GOTO 50
            ENDIF
            DO
               count = count - 1
               fs = fs*safmx2
               gs = gs*safmx2
               scale = scale*safmx2
               IF ( scale>safmn2 ) EXIT
            ENDDO
         ENDIF
         f2 = ABSSQ(fs)
         g2 = ABSSQ(gs)
         IF ( f2<=MAX(g2,ONE)*safmin ) THEN
!
!           This is a rare case: F is very small.
!
            IF ( f==CZERO ) THEN
               cs = ZERO
               r = SLAPY2(REAL(g),AIMAG(g))
!              Do complex/real division explicitly with two real
!              divisions
               d = SLAPY2(REAL(gs),AIMAG(gs))
               sn = CMPLX(REAL(gs)/d,-AIMAG(gs)/d)
               GOTO 50
            ENDIF
            f2s = SLAPY2(REAL(fs),AIMAG(fs))
!           G2 and G2S are accurate
!           G2 is at least SAFMIN, and G2S is at least SAFMN2
            g2s = SQRT(g2)
!           Error in CS from underflow in F2S is at most
!           UNFL / SAFMN2 .lt. sqrt(UNFL*EPS) .lt. EPS
!           If MAX(G2,ONE)=G2, then F2 .lt. G2*SAFMIN,
!           and so CS .lt. sqrt(SAFMIN)
!           If MAX(G2,ONE)=ONE, then F2 .lt. SAFMIN
!           and so CS .lt. sqrt(SAFMIN)/SAFMN2 = sqrt(EPS)
!           Therefore, CS = F2S/G2S / sqrt( 1 + (F2S/G2S)**2 ) = F2S/G2S
            cs = f2s/g2s
!           Make sure abs(FF) = 1
!           Do complex/real division explicitly with 2 real divisions
            IF ( ABS1(f)>ONE ) THEN
               d = SLAPY2(REAL(f),AIMAG(f))
               ff = CMPLX(REAL(f)/d,AIMAG(f)/d)
            ELSE
               dr = safmx2*REAL(f)
               di = safmx2*AIMAG(f)
               d = SLAPY2(dr,di)
               ff = CMPLX(dr/d,di/d)
            ENDIF
            sn = ff*CMPLX(REAL(gs)/g2s,-AIMAG(gs)/g2s)
            r = cs*f + sn*g
         ELSE
!
!           This is the most common case.
!           Neither F2 nor F2/G2 are less than SAFMIN
!           F2S cannot overflow, and it is accurate
!
            f2s = SQRT(ONE+g2/f2)
!           Do the F2S(real)*FS(complex) multiply with two real
!           multiplies
            r = CMPLX(f2s*REAL(fs),f2s*AIMAG(fs))
            cs = ONE/f2s
            d = f2 + g2
!           Do complex/real division explicitly with two real divisions
            sn = CMPLX(REAL(r)/d,AIMAG(r)/d)
            sn = sn*CONJG(gs)
            IF ( count/=0 ) THEN
               IF ( count>0 ) THEN
                  DO j = 1 , count
                     r = r*safmx2
                  ENDDO
               ELSE
                  DO j = 1 , -count
                     r = r*safmn2
                  ENDDO
               ENDIF
            ENDIF
         ENDIF
 50      C(ic) = cs
         Y(iy) = sn
         X(ix) = r
         ic = ic + Incc
         iy = iy + Incy
         ix = ix + Incx
      ENDDO
!
!     End of CLARGV
!
      END SUBROUTINE CLARGV
