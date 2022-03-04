!*==dlarz.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DLARZ applies an elementary reflector (as returned by stzrzf) to a general matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLARZ + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarz.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarz.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarz.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLARZ( SIDE, M, N, L, V, INCV, TAU, C, LDC, WORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          SIDE
!       INTEGER            INCV, L, LDC, M, N
!       DOUBLE PRECISION   TAU
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   C( LDC, * ), V( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLARZ applies a real elementary reflector H to a real M-by-N
!> matrix C, from either the left or the right. H is represented in the
!> form
!>
!>       H = I - tau * v * v**T
!>
!> where tau is a real scalar and v is a real vector.
!>
!> If tau = 0, then H is taken to be the unit matrix.
!>
!>
!> H is a product of k elementary reflectors as returned by DTZRZF.
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
!> \param[in] L
!> \verbatim
!>          L is INTEGER
!>          The number of entries of the vector V containing
!>          the meaningful part of the Householder vectors.
!>          If SIDE = 'L', M >= L >= 0, if SIDE = 'R', N >= L >= 0.
!> \endverbatim
!>
!> \param[in] V
!> \verbatim
!>          V is DOUBLE PRECISION array, dimension (1+(L-1)*abs(INCV))
!>          The vector v in the representation of H as returned by
!>          DTZRZF. V is not used if TAU = 0.
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
!>          TAU is DOUBLE PRECISION
!>          The value tau in the representation of H.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is DOUBLE PRECISION array, dimension (LDC,N)
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
!>          WORK is DOUBLE PRECISION array, dimension
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
!> \ingroup doubleOTHERcomputational
!
!> \par Contributors:
!  ==================
!>
!>    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE DLARZ(Side,M,N,L,V,Incv,Tau,C,Ldc,Work)
      USE F77KINDS                        
      USE S_DAXPY
      USE S_DCOPY
      USE S_DGEMV
      USE S_DGER
      USE S_LSAME
      IMPLICIT NONE
!*--DLARZ155
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Side
      INTEGER :: M
      INTEGER :: N
      INTEGER :: L
      REAL(R8KIND) , DIMENSION(*) :: V
      INTEGER :: Incv
      REAL(R8KIND) :: Tau
      REAL(R8KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL(R8KIND) , DIMENSION(*) :: Work
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
!     .. External Subroutines ..
!     ..
!     .. External Functions ..
!     ..
!     .. Executable Statements ..
!
      IF ( LSAME(Side,'L') ) THEN
!
!        Form  H * C
!
         IF ( Tau/=ZERO ) THEN
!
!           w( 1:n ) = C( 1, 1:n )
!
            CALL DCOPY(N,C,Ldc,Work,1)
!
!           w( 1:n ) = w( 1:n ) + C( m-l+1:m, 1:n )**T * v( 1:l )
!
            CALL DGEMV('Transpose',L,N,ONE,C(M-L+1,1),Ldc,V,Incv,ONE,   &
     &                 Work,1)
!
!           C( 1, 1:n ) = C( 1, 1:n ) - tau * w( 1:n )
!
            CALL DAXPY(N,-Tau,Work,1,C,Ldc)
!
!           C( m-l+1:m, 1:n ) = C( m-l+1:m, 1:n ) - ...
!                               tau * v( 1:l ) * w( 1:n )**T
!
            CALL DGER(L,N,-Tau,V,Incv,Work,1,C(M-L+1,1),Ldc)
         ENDIF
!
!
!        Form  C * H
!
      ELSEIF ( Tau/=ZERO ) THEN
!
!           w( 1:m ) = C( 1:m, 1 )
!
         CALL DCOPY(M,C,1,Work,1)
!
!           w( 1:m ) = w( 1:m ) + C( 1:m, n-l+1:n, 1:n ) * v( 1:l )
!
         CALL DGEMV('No transpose',M,L,ONE,C(1,N-L+1),Ldc,V,Incv,ONE,   &
     &              Work,1)
!
!           C( 1:m, 1 ) = C( 1:m, 1 ) - tau * w( 1:m )
!
         CALL DAXPY(M,-Tau,Work,1,C,1)
!
!           C( 1:m, n-l+1:n ) = C( 1:m, n-l+1:n ) - ...
!                               tau * w( 1:m ) * v( 1:l )**T
!
         CALL DGER(M,L,-Tau,Work,1,V,Incv,C(1,N-L+1),Ldc)
!
!
      ENDIF
!
!
!     End of DLARZ
!
      END SUBROUTINE DLARZ
