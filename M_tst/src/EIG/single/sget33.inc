!*==sget33.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b sget33
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE sget33( RMAX, LMAX, NINFO, KNT )
!
!       .. Scalar Arguments ..
!       INTEGER            KNT, LMAX, NINFO
!       REAL               RMAX
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> sget33 tests SLANV2, a routine for putting 2 by 2 blocks into
!> standard form.  In other words, it computes a two by two rotation
!> [[C,S];[-S,C]] where in
!>
!>    [ C S ][T(1,1) T(1,2)][ C -S ] = [ T11 T12 ]
!>    [-S C ][T(2,1) T(2,2)][ S  C ]   [ T21 T22 ]
!>
!> either
!>    1) T21=0 (real eigenvalues), or
!>    2) T11=T22 and T21*T12<0 (complex conjugate eigenvalues).
!> We also  verify that the residual is small.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[out] RMAX
!> \verbatim
!>          RMAX is REAL
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
!>          NINFO is INTEGER
!>          Number of examples returned with INFO .NE. 0.
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
!> \ingroup single_eig
!
!  =====================================================================
      SUBROUTINE sget33(Rmax,Lmax,Ninfo,Knt)
      IMPLICIT NONE
!*--SGET3380
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Knt , Lmax , Ninfo
      REAL Rmax
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E0,ONE=1.0E0)
      REAL TWO , FOUR
      PARAMETER (TWO=2.0E0,FOUR=4.0E0)
!     ..
!     .. Local Scalars ..
      INTEGER i1 , i2 , i3 , i4 , im1 , im2 , im3 , im4 , j1 , j2 , j3
      REAL bignum , cs , eps , res , smlnum , sn , sum , tnrm , wi1 ,   &
     &     wi2 , wr1 , wr2
!     ..
!     .. Local Arrays ..
      REAL q(2,2) , t(2,2) , t1(2,2) , t2(2,2) , val(4) , vm(3)
!     ..
!     .. External Functions ..
      REAL SLAMCH
      EXTERNAL SLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL SLABAD , SLANV2
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , SIGN
!     ..
!     .. Executable Statements ..
!
!     Get machine parameters
!
      eps = SLAMCH('P')
      smlnum = SLAMCH('S')/eps
      bignum = ONE/smlnum
      CALL SLABAD(smlnum,bignum)
!
!     Set up test case parameters
!
      val(1) = ONE
      val(2) = ONE + TWO*eps
      val(3) = TWO
      val(4) = TWO - FOUR*eps
      vm(1) = smlnum
      vm(2) = ONE
      vm(3) = bignum
!
      Knt = 0
      Ninfo = 0
      Lmax = 0
      Rmax = ZERO
!
!     Begin test loop
!
      DO i1 = 1 , 4
         DO i2 = 1 , 4
            DO i3 = 1 , 4
               DO i4 = 1 , 4
                  DO im1 = 1 , 3
                     DO im2 = 1 , 3
                        DO im3 = 1 , 3
                           DO im4 = 1 , 3
                              t(1,1) = val(i1)*vm(im1)
                              t(1,2) = val(i2)*vm(im2)
                              t(2,1) = -val(i3)*vm(im3)
                              t(2,2) = val(i4)*vm(im4)
                              tnrm = MAX(ABS(t(1,1)),ABS(t(1,2)),       &
     &                               ABS(t(2,1)),ABS(t(2,2)))
                              t1(1,1) = t(1,1)
                              t1(1,2) = t(1,2)
                              t1(2,1) = t(2,1)
                              t1(2,2) = t(2,2)
                              q(1,1) = ONE
                              q(1,2) = ZERO
                              q(2,1) = ZERO
                              q(2,2) = ONE
!
                              CALL SLANV2(t(1,1),t(1,2),t(2,1),t(2,2),  &
     &                           wr1,wi1,wr2,wi2,cs,sn)
                              DO j1 = 1 , 2
                                 res = q(j1,1)*cs + q(j1,2)*sn
                                 q(j1,2) = -q(j1,1)*sn + q(j1,2)*cs
                                 q(j1,1) = res
                              ENDDO
!
                              res = ZERO
                              res = res + ABS(q(1,1)**2+q(1,2)**2-ONE)  &
     &                              /eps
                              res = res + ABS(q(2,2)**2+q(2,1)**2-ONE)  &
     &                              /eps
                              res = res +                               &
     &                              ABS(q(1,1)*q(2,1)+q(1,2)*q(2,2))/eps
                              DO j1 = 1 , 2
                                 DO j2 = 1 , 2
                                    t2(j1,j2) = ZERO
                                    DO j3 = 1 , 2
                                       t2(j1,j2) = t2(j1,j2) + t1(j1,j3)&
     &                                    *q(j3,j2)
                                    ENDDO
                                 ENDDO
                              ENDDO
                              DO j1 = 1 , 2
                                 DO j2 = 1 , 2
                                    sum = t(j1,j2)
                                    DO j3 = 1 , 2
                                       sum = sum - q(j3,j1)*t2(j3,j2)
                                    ENDDO
                                    res = res + ABS(sum)/eps/tnrm
                                 ENDDO
                              ENDDO
                              IF ( t(2,1)/=ZERO .AND.                   &
     &                             (t(1,1)/=t(2,2) .OR. SIGN(ONE,t(1,2))&
     &                             *SIGN(ONE,t(2,1))>ZERO) ) res = res +&
     &                             ONE/eps
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
!
!
!     End of sget33
!
      END SUBROUTINE sget33
