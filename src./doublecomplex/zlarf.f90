!*==zlarf.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZLARF applies an elementary reflector to a general rectangular matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLARF + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlarf.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlarf.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlarf.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          SIDE
!       INTEGER            INCV, LDC, M, N
!       COMPLEX*16         TAU
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         C( LDC, * ), V( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLARF applies a complex elementary reflector H to a complex M-by-N
!> matrix C, from either the left or the right. H is represented in the
!> form
!>
!>       H = I - tau * v * v**H
!>
!> where tau is a complex scalar and v is a complex vector.
!>
!> If tau = 0, then H is taken to be the unit matrix.
!>
!> To apply H**H, supply conjg(tau) instead
!> tau.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'L': form  H * C
!>          = 'R': form  C * H
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix C.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix C.
!> \endverbatim
!>
!> \param[in] V
!> \verbatim
!>          V is COMPLEX*16 array, dimension
!>                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'
!>                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'
!>          The vector v in the representation of H. V is not used if
!>          TAU = 0.
!> \endverbatim
!>
!> \param[in] INCV
!> \verbatim
!>          INCV is INTEGER
!>          The increment between elements of v. INCV <> 0.
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is COMPLEX*16
!>          The value tau in the representation of H.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is COMPLEX*16 array, dimension (LDC,N)
!>          On entry, the M-by-N matrix C.
!>          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
!>          or C * H if SIDE = 'R'.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of the array C. LDC >= max(1,M).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension
!>                         (N) if SIDE = 'L'
!>                      or (M) if SIDE = 'R'
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
!> \ingroup complex16OTHERauxiliary
!
!  =====================================================================
      SUBROUTINE ZLARF(Side,M,N,V,Incv,Tau,C,Ldc,Work)
      USE F77KINDS                        
      USE S_ILAZLC
      USE S_ILAZLR
      USE S_LSAME
      USE S_ZGEMV
      USE S_ZGERC
      IMPLICIT NONE
!*--ZLARF138
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &                 ZERO = (0.0D+0,0.0D+0)
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Side
      INTEGER :: M
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(*) :: V
      INTEGER :: Incv
      COMPLEX(CX16KIND) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
!
! Local variable declarations rewritten by SPAG
!
      LOGICAL :: applyleft
      INTEGER :: i , lastc , lastv
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
!     ..
!     .. External Subroutines ..
!     ..
!     .. External Functions ..
!     ..
!     .. Executable Statements ..
!
      applyleft = LSAME(Side,'L')
      lastv = 0
      lastc = 0
      IF ( Tau/=ZERO ) THEN
!     Set up variables for scanning V.  LASTV begins pointing to the end
!     of V.
         IF ( applyleft ) THEN
            lastv = M
         ELSE
            lastv = N
         ENDIF
         IF ( Incv>0 ) THEN
            i = 1 + (lastv-1)*Incv
         ELSE
            i = 1
         ENDIF
!     Look for the last non-zero row in V.
         DO WHILE ( lastv>0 .AND. V(i)==ZERO )
            lastv = lastv - 1
            i = i - Incv
         ENDDO
         IF ( applyleft ) THEN
!     Scan for the last non-zero column in C(1:lastv,:).
            lastc = ILAZLC(lastv,N,C,Ldc)
         ELSE
!     Scan for the last non-zero row in C(:,1:lastv).
            lastc = ILAZLR(M,lastv,C,Ldc)
         ENDIF
      ENDIF
!     Note that lastc.eq.0 renders the BLAS operations null; no special
!     case is needed at this level.
      IF ( applyleft ) THEN
!
!        Form  H * C
!
         IF ( lastv>0 ) THEN
!
!           w(1:lastc,1) := C(1:lastv,1:lastc)**H * v(1:lastv,1)
!
            CALL ZGEMV('Conjugate transpose',lastv,lastc,ONE,C,Ldc,V,   &
     &                 Incv,ZERO,Work,1)
!
!           C(1:lastv,1:lastc) := C(...) - v(1:lastv,1) * w(1:lastc,1)**H
!
            CALL ZGERC(lastv,lastc,-Tau,V,Incv,Work,1,C,Ldc)
         ENDIF
!
!        Form  C * H
!
      ELSEIF ( lastv>0 ) THEN
!
!           w(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1)
!
         CALL ZGEMV('No transpose',lastc,lastv,ONE,C,Ldc,V,Incv,ZERO,   &
     &              Work,1)
!
!           C(1:lastc,1:lastv) := C(...) - w(1:lastc,1) * v(1:lastv,1)**H
!
         CALL ZGERC(lastc,lastv,-Tau,Work,1,V,Incv,C,Ldc)
      ENDIF
!
!     End of ZLARF
!
      END SUBROUTINE ZLARF
