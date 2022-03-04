!*==dget34.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b DGET34
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGET34( RMAX, LMAX, NINFO, KNT )
!
!       .. Scalar Arguments ..
!       INTEGER            KNT, LMAX
!       DOUBLE PRECISION   RMAX
!       ..
!       .. Array Arguments ..
!       INTEGER            NINFO( 2 )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGET34 tests DLAEXC, a routine for swapping adjacent blocks (either
!> 1 by 1 or 2 by 2) on the diagonal of a matrix in real Schur form.
!> Thus, DLAEXC computes an orthogonal matrix Q such that
!>
!>     Q' * [ A B ] * Q  = [ C1 B1 ]
!>          [ 0 C ]        [ 0  A1 ]
!>
!> where C1 is similar to C and A1 is similar to A.  Both A and C are
!> assumed to be in standard form (equal diagonal entries and
!> offdiagonal with differing signs) and A1 and C1 are returned with the
!> same properties.
!>
!> The test code verifies these last assertions, as well as that
!> the residual in the above equation is small.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[out] RMAX
!> \verbatim
!>          RMAX is DOUBLE PRECISION
!>          Value of the largest test ratio.
!> \endverbatim
!>
!> \param[out] LMAX
!> \verbatim
!>          LMAX is INTEGER
!>          Example number where largest test ratio achieved.
!> \endverbatim
!>
!> \param[out] NINFO
!> \verbatim
!>          NINFO is INTEGER array, dimension (2)
!>          NINFO(J) is the number of examples where INFO=J occurred.
!> \endverbatim
!>
!> \param[out] KNT
!> \verbatim
!>          KNT is INTEGER
!>          Total number of examples tested.
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
!> \ingroup double_eig
!
!  =====================================================================
      SUBROUTINE DGET34(Rmax,Lmax,Ninfo,Knt)
      IMPLICIT NONE
!*--DGET3486
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Knt , Lmax
      DOUBLE PRECISION Rmax
!     ..
!     .. Array Arguments ..
      INTEGER Ninfo(2)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , HALF , ONE
      PARAMETER (ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0)
      DOUBLE PRECISION TWO , THREE
      PARAMETER (TWO=2.0D0,THREE=3.0D0)
      INTEGER LWORK
      PARAMETER (LWORK=32)
!     ..
!     .. Local Scalars ..
      INTEGER i , ia , ia11 , ia12 , ia21 , ia22 , iam , ib , ic ,      &
     &        ic11 , ic12 , ic21 , ic22 , icm , info , j
      DOUBLE PRECISION bignum , eps , res , smlnum , tnrm
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION q(4,4) , result(2) , t(4,4) , t1(4,4) , val(9) , &
     &                 vm(2) , work(LWORK)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH
      EXTERNAL DLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL DCOPY , DHST01 , DLABAD , DLAEXC
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , DBLE , MAX , SIGN , SQRT
!     ..
!     .. Executable Statements ..
!
!     Get machine parameters
!
      eps = DLAMCH('P')
      smlnum = DLAMCH('S')/eps
      bignum = ONE/smlnum
      CALL DLABAD(smlnum,bignum)
!
!     Set up test case parameters
!
      val(1) = ZERO
      val(2) = SQRT(smlnum)
      val(3) = ONE
      val(4) = TWO
      val(5) = SQRT(bignum)
      val(6) = -SQRT(smlnum)
      val(7) = -ONE
      val(8) = -TWO
      val(9) = -SQRT(bignum)
      vm(1) = ONE
      vm(2) = ONE + TWO*eps
      CALL DCOPY(16,val(4),0,t(1,1),1)
!
      Ninfo(1) = 0
      Ninfo(2) = 0
      Knt = 0
      Lmax = 0
      Rmax = ZERO
!
!     Begin test loop
!
      DO ia = 1 , 9
         DO iam = 1 , 2
            DO ib = 1 , 9
               DO ic = 1 , 9
                  t(1,1) = val(ia)*vm(iam)
                  t(2,2) = val(ic)
                  t(1,2) = val(ib)
                  t(2,1) = ZERO
                  tnrm = MAX(ABS(t(1,1)),ABS(t(2,2)),ABS(t(1,2)))
                  CALL DCOPY(16,t,1,t1,1)
                  CALL DCOPY(16,val(1),0,q,1)
                  CALL DCOPY(4,val(3),0,q,5)
                  CALL DLAEXC(.TRUE.,2,t,4,q,4,1,1,1,work,info)
                  IF ( info/=0 ) Ninfo(info) = Ninfo(info) + 1
                  CALL DHST01(2,1,2,t1,4,t,4,q,4,work,LWORK,result)
                  res = result(1) + result(2)
                  IF ( info/=0 ) res = res + ONE/eps
                  IF ( t(1,1)/=t1(2,2) ) res = res + ONE/eps
                  IF ( t(2,2)/=t1(1,1) ) res = res + ONE/eps
                  IF ( t(2,1)/=ZERO ) res = res + ONE/eps
                  Knt = Knt + 1
                  IF ( res>Rmax ) THEN
                     Lmax = Knt
                     Rmax = res
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO
!
      DO ia = 1 , 5
         DO iam = 1 , 2
            DO ib = 1 , 5
               DO ic11 = 1 , 5
                  DO ic12 = 2 , 5
                     DO ic21 = 2 , 4
                        DO ic22 = -1 , 1 , 2
                           t(1,1) = val(ia)*vm(iam)
                           t(1,2) = val(ib)
                           t(1,3) = -TWO*val(ib)
                           t(2,1) = ZERO
                           t(2,2) = val(ic11)
                           t(2,3) = val(ic12)
                           t(3,1) = ZERO
                           t(3,2) = -val(ic21)
                           t(3,3) = val(ic11)*DBLE(ic22)
                           tnrm = MAX(ABS(t(1,1)),ABS(t(1,2)),          &
     &                            ABS(t(1,3)),ABS(t(2,2)),ABS(t(2,3)),  &
     &                            ABS(t(3,2)),ABS(t(3,3)))
                           CALL DCOPY(16,t,1,t1,1)
                           CALL DCOPY(16,val(1),0,q,1)
                           CALL DCOPY(4,val(3),0,q,5)
                           CALL DLAEXC(.TRUE.,3,t,4,q,4,1,1,2,work,info)
                           IF ( info/=0 ) Ninfo(info) = Ninfo(info) + 1
                           CALL DHST01(3,1,3,t1,4,t,4,q,4,work,LWORK,   &
     &                                 result)
                           res = result(1) + result(2)
                           IF ( info==0 ) THEN
                              IF ( t1(1,1)/=t(3,3) ) res = res + ONE/eps
                              IF ( t(3,1)/=ZERO ) res = res + ONE/eps
                              IF ( t(3,2)/=ZERO ) res = res + ONE/eps
                              IF ( t(2,1)/=0 .AND.                      &
     &                             (t(1,1)/=t(2,2) .OR. SIGN(ONE,t(1,2))&
     &                             ==SIGN(ONE,t(2,1))) ) res = res +    &
     &                             ONE/eps
                           ENDIF
                           Knt = Knt + 1
                           IF ( res>Rmax ) THEN
                              Lmax = Knt
                              Rmax = res
                           ENDIF
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
!
      DO ia11 = 1 , 5
         DO ia12 = 2 , 5
            DO ia21 = 2 , 4
               DO ia22 = -1 , 1 , 2
                  DO icm = 1 , 2
                     DO ib = 1 , 5
                        DO ic = 1 , 5
                           t(1,1) = val(ia11)
                           t(1,2) = val(ia12)
                           t(1,3) = -TWO*val(ib)
                           t(2,1) = -val(ia21)
                           t(2,2) = val(ia11)*DBLE(ia22)
                           t(2,3) = val(ib)
                           t(3,1) = ZERO
                           t(3,2) = ZERO
                           t(3,3) = val(ic)*vm(icm)
                           tnrm = MAX(ABS(t(1,1)),ABS(t(1,2)),          &
     &                            ABS(t(1,3)),ABS(t(2,2)),ABS(t(2,3)),  &
     &                            ABS(t(3,2)),ABS(t(3,3)))
                           CALL DCOPY(16,t,1,t1,1)
                           CALL DCOPY(16,val(1),0,q,1)
                           CALL DCOPY(4,val(3),0,q,5)
                           CALL DLAEXC(.TRUE.,3,t,4,q,4,1,2,1,work,info)
                           IF ( info/=0 ) Ninfo(info) = Ninfo(info) + 1
                           CALL DHST01(3,1,3,t1,4,t,4,q,4,work,LWORK,   &
     &                                 result)
                           res = result(1) + result(2)
                           IF ( info==0 ) THEN
                              IF ( t1(3,3)/=t(1,1) ) res = res + ONE/eps
                              IF ( t(2,1)/=ZERO ) res = res + ONE/eps
                              IF ( t(3,1)/=ZERO ) res = res + ONE/eps
                              IF ( t(3,2)/=0 .AND.                      &
     &                             (t(2,2)/=t(3,3) .OR. SIGN(ONE,t(2,3))&
     &                             ==SIGN(ONE,t(3,2))) ) res = res +    &
     &                             ONE/eps
                           ENDIF
                           Knt = Knt + 1
                           IF ( res>Rmax ) THEN
                              Lmax = Knt
                              Rmax = res
                           ENDIF
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
!
      DO ia11 = 1 , 5
         DO ia12 = 2 , 5
            DO ia21 = 2 , 4
               DO ia22 = -1 , 1 , 2
                  DO ib = 1 , 5
                     DO ic11 = 3 , 4
                        DO ic12 = 3 , 4
                           DO ic21 = 3 , 4
                              DO ic22 = -1 , 1 , 2
                                 DO icm = 5 , 7
                                    iam = 1
                                    t(1,1) = val(ia11)*vm(iam)
                                    t(1,2) = val(ia12)*vm(iam)
                                    t(1,3) = -TWO*val(ib)
                                    t(1,4) = HALF*val(ib)
                                    t(2,1) = -t(1,2)*val(ia21)
                                    t(2,2) = val(ia11)*DBLE(ia22)       &
     &                                 *vm(iam)
                                    t(2,3) = val(ib)
                                    t(2,4) = THREE*val(ib)
                                    t(3,1) = ZERO
                                    t(3,2) = ZERO
                                    t(3,3) = val(ic11)*ABS(val(icm))
                                    t(3,4) = val(ic12)*ABS(val(icm))
                                    t(4,1) = ZERO
                                    t(4,2) = ZERO
                                    t(4,3) = -t(3,4)*val(ic21)          &
     &                                 *ABS(val(icm))
                                    t(4,4) = val(ic11)*DBLE(ic22)       &
     &                                 *ABS(val(icm))
                                    tnrm = ZERO
                                    DO i = 1 , 4
                                       DO j = 1 , 4
                                         tnrm = MAX(tnrm,ABS(t(i,j)))
                                       ENDDO
                                    ENDDO
                                    CALL DCOPY(16,t,1,t1,1)
                                    CALL DCOPY(16,val(1),0,q,1)
                                    CALL DCOPY(4,val(3),0,q,5)
                                    CALL DLAEXC(.TRUE.,4,t,4,q,4,1,2,2, &
     &                                 work,info)
                                    IF ( info/=0 ) Ninfo(info)          &
     &                                 = Ninfo(info) + 1
                                    CALL DHST01(4,1,4,t1,4,t,4,q,4,work,&
     &                                 LWORK,result)
                                    res = result(1) + result(2)
                                    IF ( info==0 ) THEN
                                       IF ( t(3,1)/=ZERO ) res = res +  &
     &                                    ONE/eps
                                       IF ( t(4,1)/=ZERO ) res = res +  &
     &                                    ONE/eps
                                       IF ( t(3,2)/=ZERO ) res = res +  &
     &                                    ONE/eps
                                       IF ( t(4,2)/=ZERO ) res = res +  &
     &                                    ONE/eps
                                       IF ( t(2,1)/=0 .AND.             &
     &                                    (t(1,1)/=t(2,2) .OR.          &
     &                                    SIGN(ONE,t(1,2))              &
     &                                    ==SIGN(ONE,t(2,1))) )         &
     &                                    res = res + ONE/eps
                                       IF ( t(4,3)/=0 .AND.             &
     &                                    (t(3,3)/=t(4,4) .OR.          &
     &                                    SIGN(ONE,t(3,4))              &
     &                                    ==SIGN(ONE,t(4,3))) )         &
     &                                    res = res + ONE/eps
                                    ENDIF
                                    Knt = Knt + 1
                                    IF ( res>Rmax ) THEN
                                       Lmax = Knt
                                       Rmax = res
                                    ENDIF
                                 ENDDO
                              ENDDO
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
!
!
!     End of DGET34
!
      END SUBROUTINE DGET34
