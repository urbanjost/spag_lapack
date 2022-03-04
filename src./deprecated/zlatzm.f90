!*==zlatzm.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZLATZM
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLATZM + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlatzm.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlatzm.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlatzm.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLATZM( SIDE, M, N, V, INCV, TAU, C1, C2, LDC, WORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          SIDE
!       INTEGER            INCV, LDC, M, N
!       COMPLEX*16         TAU
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         C1( LDC, * ), C2( LDC, * ), V( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> This routine is deprecated and has been replaced by routine ZUNMRZ.
!>
!> ZLATZM applies a Householder matrix generated by ZTZRQF to a matrix.
!>
!> Let P = I - tau*u*u**H,   u = ( 1 ),
!>                               ( v )
!> where v is an (m-1) vector if SIDE = 'L', or a (n-1) vector if
!> SIDE = 'R'.
!>
!> If SIDE equals 'L', let
!>        C = [ C1 ] 1
!>            [ C2 ] m-1
!>              n
!> Then C is overwritten by P*C.
!>
!> If SIDE equals 'R', let
!>        C = [ C1, C2 ] m
!>               1  n-1
!> Then C is overwritten by C*P.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'L': form P * C
!>          = 'R': form C * P
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
!>                  (1 + (M-1)*abs(INCV)) if SIDE = 'L'
!>                  (1 + (N-1)*abs(INCV)) if SIDE = 'R'
!>          The vector v in the representation of P. V is not used
!>          if TAU = 0.
!> \endverbatim
!>
!> \param[in] INCV
!> \verbatim
!>          INCV is INTEGER
!>          The increment between elements of v. INCV <> 0
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is COMPLEX*16
!>          The value tau in the representation of P.
!> \endverbatim
!>
!> \param[in,out] C1
!> \verbatim
!>          C1 is COMPLEX*16 array, dimension
!>                         (LDC,N) if SIDE = 'L'
!>                         (M,1)   if SIDE = 'R'
!>          On entry, the n-vector C1 if SIDE = 'L', or the m-vector C1
!>          if SIDE = 'R'.
!>
!>          On exit, the first row of P*C if SIDE = 'L', or the first
!>          column of C*P if SIDE = 'R'.
!> \endverbatim
!>
!> \param[in,out] C2
!> \verbatim
!>          C2 is COMPLEX*16 array, dimension
!>                         (LDC, N)   if SIDE = 'L'
!>                         (LDC, N-1) if SIDE = 'R'
!>          On entry, the (m - 1) x n matrix C2 if SIDE = 'L', or the
!>          m x (n - 1) matrix C2 if SIDE = 'R'.
!>
!>          On exit, rows 2:m of P*C if SIDE = 'L', or columns 2:m of C*P
!>          if SIDE = 'R'.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of the arrays C1 and C2.
!>          LDC >= max(1,M).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension
!>                      (N) if SIDE = 'L'
!>                      (M) if SIDE = 'R'
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
!> \ingroup complex16OTHERcomputational
!
!  =====================================================================
      SUBROUTINE ZLATZM(Side,M,N,V,Incv,Tau,C1,C2,Ldc,Work)
      USE F77KINDS                        
      USE S_LSAME
      USE S_ZAXPY
      USE S_ZCOPY
      USE S_ZGEMV
      USE S_ZGERC
      USE S_ZGERU
      USE S_ZLACGV
      IMPLICIT NONE
!*--ZLATZM164
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
      COMPLEX(CX16KIND) , DIMENSION(Ldc,*) :: C1
      COMPLEX(CX16KIND) , DIMENSION(Ldc,*) :: C2
      INTEGER :: Ldc
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
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
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
      IF ( (MIN(M,N)==0) .OR. (Tau==ZERO) ) RETURN
!
      IF ( LSAME(Side,'L') ) THEN
!
!        w :=  ( C1 + v**H * C2 )**H
!
         CALL ZCOPY(N,C1,Ldc,Work,1)
         CALL ZLACGV(N,Work,1)
         CALL ZGEMV('Conjugate transpose',M-1,N,ONE,C2,Ldc,V,Incv,ONE,  &
     &              Work,1)
!
!        [ C1 ] := [ C1 ] - tau* [ 1 ] * w**H
!        [ C2 ]    [ C2 ]        [ v ]
!
         CALL ZLACGV(N,Work,1)
         CALL ZAXPY(N,-Tau,Work,1,C1,Ldc)
         CALL ZGERU(M-1,N,-Tau,V,Incv,Work,1,C2,Ldc)
!
      ELSEIF ( LSAME(Side,'R') ) THEN
!
!        w := C1 + C2 * v
!
         CALL ZCOPY(M,C1,1,Work,1)
         CALL ZGEMV('No transpose',M,N-1,ONE,C2,Ldc,V,Incv,ONE,Work,1)
!
!        [ C1, C2 ] := [ C1, C2 ] - tau* w * [ 1 , v**H]
!
         CALL ZAXPY(M,-Tau,Work,1,C1,1)
         CALL ZGERC(M,N-1,-Tau,Work,1,V,Incv,C2,Ldc)
      ENDIF
!
!
!     End of ZLATZM
!
      END SUBROUTINE ZLATZM