!*==clarfx.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CLARFX applies an elementary reflector to a general rectangular matrix, with loop unrolling when the reflector has order ≤ 10.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLARFX + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clarfx.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clarfx.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clarfx.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLARFX( SIDE, M, N, V, TAU, C, LDC, WORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          SIDE
!       INTEGER            LDC, M, N
!       COMPLEX            TAU
!       ..
!       .. Array Arguments ..
!       COMPLEX            C( LDC, * ), V( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLARFX applies a complex elementary reflector H to a complex m by n
!> matrix C, from either the left or the right. H is represented in the
!> form
!>
!>       H = I - tau * v * v**H
!>
!> where tau is a complex scalar and v is a complex vector.
!>
!> If tau = 0, then H is taken to be the unit matrix
!>
!> This version uses inline code if H has order < 11.
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
!>          V is COMPLEX array, dimension (M) if SIDE = 'L'
!>                                        or (N) if SIDE = 'R'
!>          The vector v in the representation of H.
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is COMPLEX
!>          The value tau in the representation of H.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is COMPLEX array, dimension (LDC,N)
!>          On entry, the m by n matrix C.
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
!>          WORK is COMPLEX array, dimension (N) if SIDE = 'L'
!>                                            or (M) if SIDE = 'R'
!>          WORK is not referenced if H has order < 11.
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
!  =====================================================================
      SUBROUTINE CLARFX(Side,M,N,V,Tau,C,Ldc,Work)
      USE S_CLARF
      USE S_LSAME
      IMPLICIT NONE
!*--CLARFX125
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0) ,                 &
     &                         ONE = (1.0E+0,0.0E+0)
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Side
      INTEGER :: M
      INTEGER :: N
      COMPLEX , DIMENSION(*) :: V
      COMPLEX :: Tau
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX , DIMENSION(*) :: Work
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: j
      COMPLEX :: sum , t1 , t10 , t2 , t3 , t4 , t5 , t6 , t7 , t8 ,    &
     &           t9 , v1 , v10 , v2 , v3 , v4 , v5 , v6 , v7 , v8 , v9
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
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
      IF ( Tau==ZERO ) RETURN
      IF ( LSAME(Side,'L') ) THEN
!
!        Form  H * C, where H has order m.
!
         IF ( M==1 ) THEN
!
!        Special code for 1 x 1 Householder
!
            t1 = ONE - Tau*V(1)*CONJG(V(1))
            DO j = 1 , N
               C(1,j) = t1*C(1,j)
            ENDDO
         ELSEIF ( M==2 ) THEN
!
!        Special code for 2 x 2 Householder
!
            v1 = CONJG(V(1))
            t1 = Tau*CONJG(v1)
            v2 = CONJG(V(2))
            t2 = Tau*CONJG(v2)
            DO j = 1 , N
               sum = v1*C(1,j) + v2*C(2,j)
               C(1,j) = C(1,j) - sum*t1
               C(2,j) = C(2,j) - sum*t2
            ENDDO
         ELSEIF ( M==3 ) THEN
!
!        Special code for 3 x 3 Householder
!
            v1 = CONJG(V(1))
            t1 = Tau*CONJG(v1)
            v2 = CONJG(V(2))
            t2 = Tau*CONJG(v2)
            v3 = CONJG(V(3))
            t3 = Tau*CONJG(v3)
            DO j = 1 , N
               sum = v1*C(1,j) + v2*C(2,j) + v3*C(3,j)
               C(1,j) = C(1,j) - sum*t1
               C(2,j) = C(2,j) - sum*t2
               C(3,j) = C(3,j) - sum*t3
            ENDDO
         ELSEIF ( M==4 ) THEN
!
!        Special code for 4 x 4 Householder
!
            v1 = CONJG(V(1))
            t1 = Tau*CONJG(v1)
            v2 = CONJG(V(2))
            t2 = Tau*CONJG(v2)
            v3 = CONJG(V(3))
            t3 = Tau*CONJG(v3)
            v4 = CONJG(V(4))
            t4 = Tau*CONJG(v4)
            DO j = 1 , N
               sum = v1*C(1,j) + v2*C(2,j) + v3*C(3,j) + v4*C(4,j)
               C(1,j) = C(1,j) - sum*t1
               C(2,j) = C(2,j) - sum*t2
               C(3,j) = C(3,j) - sum*t3
               C(4,j) = C(4,j) - sum*t4
            ENDDO
         ELSEIF ( M==5 ) THEN
!
!        Special code for 5 x 5 Householder
!
            v1 = CONJG(V(1))
            t1 = Tau*CONJG(v1)
            v2 = CONJG(V(2))
            t2 = Tau*CONJG(v2)
            v3 = CONJG(V(3))
            t3 = Tau*CONJG(v3)
            v4 = CONJG(V(4))
            t4 = Tau*CONJG(v4)
            v5 = CONJG(V(5))
            t5 = Tau*CONJG(v5)
            DO j = 1 , N
               sum = v1*C(1,j) + v2*C(2,j) + v3*C(3,j) + v4*C(4,j)      &
     &               + v5*C(5,j)
               C(1,j) = C(1,j) - sum*t1
               C(2,j) = C(2,j) - sum*t2
               C(3,j) = C(3,j) - sum*t3
               C(4,j) = C(4,j) - sum*t4
               C(5,j) = C(5,j) - sum*t5
            ENDDO
         ELSEIF ( M==6 ) THEN
!
!        Special code for 6 x 6 Householder
!
            v1 = CONJG(V(1))
            t1 = Tau*CONJG(v1)
            v2 = CONJG(V(2))
            t2 = Tau*CONJG(v2)
            v3 = CONJG(V(3))
            t3 = Tau*CONJG(v3)
            v4 = CONJG(V(4))
            t4 = Tau*CONJG(v4)
            v5 = CONJG(V(5))
            t5 = Tau*CONJG(v5)
            v6 = CONJG(V(6))
            t6 = Tau*CONJG(v6)
            DO j = 1 , N
               sum = v1*C(1,j) + v2*C(2,j) + v3*C(3,j) + v4*C(4,j)      &
     &               + v5*C(5,j) + v6*C(6,j)
               C(1,j) = C(1,j) - sum*t1
               C(2,j) = C(2,j) - sum*t2
               C(3,j) = C(3,j) - sum*t3
               C(4,j) = C(4,j) - sum*t4
               C(5,j) = C(5,j) - sum*t5
               C(6,j) = C(6,j) - sum*t6
            ENDDO
         ELSEIF ( M==7 ) THEN
!
!        Special code for 7 x 7 Householder
!
            v1 = CONJG(V(1))
            t1 = Tau*CONJG(v1)
            v2 = CONJG(V(2))
            t2 = Tau*CONJG(v2)
            v3 = CONJG(V(3))
            t3 = Tau*CONJG(v3)
            v4 = CONJG(V(4))
            t4 = Tau*CONJG(v4)
            v5 = CONJG(V(5))
            t5 = Tau*CONJG(v5)
            v6 = CONJG(V(6))
            t6 = Tau*CONJG(v6)
            v7 = CONJG(V(7))
            t7 = Tau*CONJG(v7)
            DO j = 1 , N
               sum = v1*C(1,j) + v2*C(2,j) + v3*C(3,j) + v4*C(4,j)      &
     &               + v5*C(5,j) + v6*C(6,j) + v7*C(7,j)
               C(1,j) = C(1,j) - sum*t1
               C(2,j) = C(2,j) - sum*t2
               C(3,j) = C(3,j) - sum*t3
               C(4,j) = C(4,j) - sum*t4
               C(5,j) = C(5,j) - sum*t5
               C(6,j) = C(6,j) - sum*t6
               C(7,j) = C(7,j) - sum*t7
            ENDDO
         ELSEIF ( M==8 ) THEN
!
!        Special code for 8 x 8 Householder
!
            v1 = CONJG(V(1))
            t1 = Tau*CONJG(v1)
            v2 = CONJG(V(2))
            t2 = Tau*CONJG(v2)
            v3 = CONJG(V(3))
            t3 = Tau*CONJG(v3)
            v4 = CONJG(V(4))
            t4 = Tau*CONJG(v4)
            v5 = CONJG(V(5))
            t5 = Tau*CONJG(v5)
            v6 = CONJG(V(6))
            t6 = Tau*CONJG(v6)
            v7 = CONJG(V(7))
            t7 = Tau*CONJG(v7)
            v8 = CONJG(V(8))
            t8 = Tau*CONJG(v8)
            DO j = 1 , N
               sum = v1*C(1,j) + v2*C(2,j) + v3*C(3,j) + v4*C(4,j)      &
     &               + v5*C(5,j) + v6*C(6,j) + v7*C(7,j) + v8*C(8,j)
               C(1,j) = C(1,j) - sum*t1
               C(2,j) = C(2,j) - sum*t2
               C(3,j) = C(3,j) - sum*t3
               C(4,j) = C(4,j) - sum*t4
               C(5,j) = C(5,j) - sum*t5
               C(6,j) = C(6,j) - sum*t6
               C(7,j) = C(7,j) - sum*t7
               C(8,j) = C(8,j) - sum*t8
            ENDDO
         ELSEIF ( M==9 ) THEN
!
!        Special code for 9 x 9 Householder
!
            v1 = CONJG(V(1))
            t1 = Tau*CONJG(v1)
            v2 = CONJG(V(2))
            t2 = Tau*CONJG(v2)
            v3 = CONJG(V(3))
            t3 = Tau*CONJG(v3)
            v4 = CONJG(V(4))
            t4 = Tau*CONJG(v4)
            v5 = CONJG(V(5))
            t5 = Tau*CONJG(v5)
            v6 = CONJG(V(6))
            t6 = Tau*CONJG(v6)
            v7 = CONJG(V(7))
            t7 = Tau*CONJG(v7)
            v8 = CONJG(V(8))
            t8 = Tau*CONJG(v8)
            v9 = CONJG(V(9))
            t9 = Tau*CONJG(v9)
            DO j = 1 , N
               sum = v1*C(1,j) + v2*C(2,j) + v3*C(3,j) + v4*C(4,j)      &
     &               + v5*C(5,j) + v6*C(6,j) + v7*C(7,j) + v8*C(8,j)    &
     &               + v9*C(9,j)
               C(1,j) = C(1,j) - sum*t1
               C(2,j) = C(2,j) - sum*t2
               C(3,j) = C(3,j) - sum*t3
               C(4,j) = C(4,j) - sum*t4
               C(5,j) = C(5,j) - sum*t5
               C(6,j) = C(6,j) - sum*t6
               C(7,j) = C(7,j) - sum*t7
               C(8,j) = C(8,j) - sum*t8
               C(9,j) = C(9,j) - sum*t9
            ENDDO
         ELSEIF ( M==10 ) THEN
!
!        Special code for 10 x 10 Householder
!
            v1 = CONJG(V(1))
            t1 = Tau*CONJG(v1)
            v2 = CONJG(V(2))
            t2 = Tau*CONJG(v2)
            v3 = CONJG(V(3))
            t3 = Tau*CONJG(v3)
            v4 = CONJG(V(4))
            t4 = Tau*CONJG(v4)
            v5 = CONJG(V(5))
            t5 = Tau*CONJG(v5)
            v6 = CONJG(V(6))
            t6 = Tau*CONJG(v6)
            v7 = CONJG(V(7))
            t7 = Tau*CONJG(v7)
            v8 = CONJG(V(8))
            t8 = Tau*CONJG(v8)
            v9 = CONJG(V(9))
            t9 = Tau*CONJG(v9)
            v10 = CONJG(V(10))
            t10 = Tau*CONJG(v10)
            DO j = 1 , N
               sum = v1*C(1,j) + v2*C(2,j) + v3*C(3,j) + v4*C(4,j)      &
     &               + v5*C(5,j) + v6*C(6,j) + v7*C(7,j) + v8*C(8,j)    &
     &               + v9*C(9,j) + v10*C(10,j)
               C(1,j) = C(1,j) - sum*t1
               C(2,j) = C(2,j) - sum*t2
               C(3,j) = C(3,j) - sum*t3
               C(4,j) = C(4,j) - sum*t4
               C(5,j) = C(5,j) - sum*t5
               C(6,j) = C(6,j) - sum*t6
               C(7,j) = C(7,j) - sum*t7
               C(8,j) = C(8,j) - sum*t8
               C(9,j) = C(9,j) - sum*t9
               C(10,j) = C(10,j) - sum*t10
            ENDDO
         ELSE
!
!        Code for general M
!
            CALL CLARF(Side,M,N,V,1,Tau,C,Ldc,Work)
         ENDIF
!
!        Form  C * H, where H has order n.
!
      ELSEIF ( N==1 ) THEN
!
!        Special code for 1 x 1 Householder
!
         t1 = ONE - Tau*V(1)*CONJG(V(1))
         DO j = 1 , M
            C(j,1) = t1*C(j,1)
         ENDDO
      ELSEIF ( N==2 ) THEN
!
!        Special code for 2 x 2 Householder
!
         v1 = V(1)
         t1 = Tau*CONJG(v1)
         v2 = V(2)
         t2 = Tau*CONJG(v2)
         DO j = 1 , M
            sum = v1*C(j,1) + v2*C(j,2)
            C(j,1) = C(j,1) - sum*t1
            C(j,2) = C(j,2) - sum*t2
         ENDDO
      ELSEIF ( N==3 ) THEN
!
!        Special code for 3 x 3 Householder
!
         v1 = V(1)
         t1 = Tau*CONJG(v1)
         v2 = V(2)
         t2 = Tau*CONJG(v2)
         v3 = V(3)
         t3 = Tau*CONJG(v3)
         DO j = 1 , M
            sum = v1*C(j,1) + v2*C(j,2) + v3*C(j,3)
            C(j,1) = C(j,1) - sum*t1
            C(j,2) = C(j,2) - sum*t2
            C(j,3) = C(j,3) - sum*t3
         ENDDO
      ELSEIF ( N==4 ) THEN
!
!        Special code for 4 x 4 Householder
!
         v1 = V(1)
         t1 = Tau*CONJG(v1)
         v2 = V(2)
         t2 = Tau*CONJG(v2)
         v3 = V(3)
         t3 = Tau*CONJG(v3)
         v4 = V(4)
         t4 = Tau*CONJG(v4)
         DO j = 1 , M
            sum = v1*C(j,1) + v2*C(j,2) + v3*C(j,3) + v4*C(j,4)
            C(j,1) = C(j,1) - sum*t1
            C(j,2) = C(j,2) - sum*t2
            C(j,3) = C(j,3) - sum*t3
            C(j,4) = C(j,4) - sum*t4
         ENDDO
      ELSEIF ( N==5 ) THEN
!
!        Special code for 5 x 5 Householder
!
         v1 = V(1)
         t1 = Tau*CONJG(v1)
         v2 = V(2)
         t2 = Tau*CONJG(v2)
         v3 = V(3)
         t3 = Tau*CONJG(v3)
         v4 = V(4)
         t4 = Tau*CONJG(v4)
         v5 = V(5)
         t5 = Tau*CONJG(v5)
         DO j = 1 , M
            sum = v1*C(j,1) + v2*C(j,2) + v3*C(j,3) + v4*C(j,4)         &
     &            + v5*C(j,5)
            C(j,1) = C(j,1) - sum*t1
            C(j,2) = C(j,2) - sum*t2
            C(j,3) = C(j,3) - sum*t3
            C(j,4) = C(j,4) - sum*t4
            C(j,5) = C(j,5) - sum*t5
         ENDDO
      ELSEIF ( N==6 ) THEN
!
!        Special code for 6 x 6 Householder
!
         v1 = V(1)
         t1 = Tau*CONJG(v1)
         v2 = V(2)
         t2 = Tau*CONJG(v2)
         v3 = V(3)
         t3 = Tau*CONJG(v3)
         v4 = V(4)
         t4 = Tau*CONJG(v4)
         v5 = V(5)
         t5 = Tau*CONJG(v5)
         v6 = V(6)
         t6 = Tau*CONJG(v6)
         DO j = 1 , M
            sum = v1*C(j,1) + v2*C(j,2) + v3*C(j,3) + v4*C(j,4)         &
     &            + v5*C(j,5) + v6*C(j,6)
            C(j,1) = C(j,1) - sum*t1
            C(j,2) = C(j,2) - sum*t2
            C(j,3) = C(j,3) - sum*t3
            C(j,4) = C(j,4) - sum*t4
            C(j,5) = C(j,5) - sum*t5
            C(j,6) = C(j,6) - sum*t6
         ENDDO
      ELSEIF ( N==7 ) THEN
!
!        Special code for 7 x 7 Householder
!
         v1 = V(1)
         t1 = Tau*CONJG(v1)
         v2 = V(2)
         t2 = Tau*CONJG(v2)
         v3 = V(3)
         t3 = Tau*CONJG(v3)
         v4 = V(4)
         t4 = Tau*CONJG(v4)
         v5 = V(5)
         t5 = Tau*CONJG(v5)
         v6 = V(6)
         t6 = Tau*CONJG(v6)
         v7 = V(7)
         t7 = Tau*CONJG(v7)
         DO j = 1 , M
            sum = v1*C(j,1) + v2*C(j,2) + v3*C(j,3) + v4*C(j,4)         &
     &            + v5*C(j,5) + v6*C(j,6) + v7*C(j,7)
            C(j,1) = C(j,1) - sum*t1
            C(j,2) = C(j,2) - sum*t2
            C(j,3) = C(j,3) - sum*t3
            C(j,4) = C(j,4) - sum*t4
            C(j,5) = C(j,5) - sum*t5
            C(j,6) = C(j,6) - sum*t6
            C(j,7) = C(j,7) - sum*t7
         ENDDO
      ELSEIF ( N==8 ) THEN
!
!        Special code for 8 x 8 Householder
!
         v1 = V(1)
         t1 = Tau*CONJG(v1)
         v2 = V(2)
         t2 = Tau*CONJG(v2)
         v3 = V(3)
         t3 = Tau*CONJG(v3)
         v4 = V(4)
         t4 = Tau*CONJG(v4)
         v5 = V(5)
         t5 = Tau*CONJG(v5)
         v6 = V(6)
         t6 = Tau*CONJG(v6)
         v7 = V(7)
         t7 = Tau*CONJG(v7)
         v8 = V(8)
         t8 = Tau*CONJG(v8)
         DO j = 1 , M
            sum = v1*C(j,1) + v2*C(j,2) + v3*C(j,3) + v4*C(j,4)         &
     &            + v5*C(j,5) + v6*C(j,6) + v7*C(j,7) + v8*C(j,8)
            C(j,1) = C(j,1) - sum*t1
            C(j,2) = C(j,2) - sum*t2
            C(j,3) = C(j,3) - sum*t3
            C(j,4) = C(j,4) - sum*t4
            C(j,5) = C(j,5) - sum*t5
            C(j,6) = C(j,6) - sum*t6
            C(j,7) = C(j,7) - sum*t7
            C(j,8) = C(j,8) - sum*t8
         ENDDO
      ELSEIF ( N==9 ) THEN
!
!        Special code for 9 x 9 Householder
!
         v1 = V(1)
         t1 = Tau*CONJG(v1)
         v2 = V(2)
         t2 = Tau*CONJG(v2)
         v3 = V(3)
         t3 = Tau*CONJG(v3)
         v4 = V(4)
         t4 = Tau*CONJG(v4)
         v5 = V(5)
         t5 = Tau*CONJG(v5)
         v6 = V(6)
         t6 = Tau*CONJG(v6)
         v7 = V(7)
         t7 = Tau*CONJG(v7)
         v8 = V(8)
         t8 = Tau*CONJG(v8)
         v9 = V(9)
         t9 = Tau*CONJG(v9)
         DO j = 1 , M
            sum = v1*C(j,1) + v2*C(j,2) + v3*C(j,3) + v4*C(j,4)         &
     &            + v5*C(j,5) + v6*C(j,6) + v7*C(j,7) + v8*C(j,8)       &
     &            + v9*C(j,9)
            C(j,1) = C(j,1) - sum*t1
            C(j,2) = C(j,2) - sum*t2
            C(j,3) = C(j,3) - sum*t3
            C(j,4) = C(j,4) - sum*t4
            C(j,5) = C(j,5) - sum*t5
            C(j,6) = C(j,6) - sum*t6
            C(j,7) = C(j,7) - sum*t7
            C(j,8) = C(j,8) - sum*t8
            C(j,9) = C(j,9) - sum*t9
         ENDDO
      ELSEIF ( N==10 ) THEN
!
!        Special code for 10 x 10 Householder
!
         v1 = V(1)
         t1 = Tau*CONJG(v1)
         v2 = V(2)
         t2 = Tau*CONJG(v2)
         v3 = V(3)
         t3 = Tau*CONJG(v3)
         v4 = V(4)
         t4 = Tau*CONJG(v4)
         v5 = V(5)
         t5 = Tau*CONJG(v5)
         v6 = V(6)
         t6 = Tau*CONJG(v6)
         v7 = V(7)
         t7 = Tau*CONJG(v7)
         v8 = V(8)
         t8 = Tau*CONJG(v8)
         v9 = V(9)
         t9 = Tau*CONJG(v9)
         v10 = V(10)
         t10 = Tau*CONJG(v10)
         DO j = 1 , M
            sum = v1*C(j,1) + v2*C(j,2) + v3*C(j,3) + v4*C(j,4)         &
     &            + v5*C(j,5) + v6*C(j,6) + v7*C(j,7) + v8*C(j,8)       &
     &            + v9*C(j,9) + v10*C(j,10)
            C(j,1) = C(j,1) - sum*t1
            C(j,2) = C(j,2) - sum*t2
            C(j,3) = C(j,3) - sum*t3
            C(j,4) = C(j,4) - sum*t4
            C(j,5) = C(j,5) - sum*t5
            C(j,6) = C(j,6) - sum*t6
            C(j,7) = C(j,7) - sum*t7
            C(j,8) = C(j,8) - sum*t8
            C(j,9) = C(j,9) - sum*t9
            C(j,10) = C(j,10) - sum*t10
         ENDDO
      ELSE
!
!        Code for general N
!
         CALL CLARF(Side,M,N,V,1,Tau,C,Ldc,Work)
      ENDIF
!
!     End of CLARFX
!
      END SUBROUTINE CLARFX
