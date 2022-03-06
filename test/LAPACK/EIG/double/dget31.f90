!*==dget31.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b DGET31
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGET31( RMAX, LMAX, NINFO, KNT )
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
!> DGET31 tests DLALN2, a routine for solving
!>
!>    (ca A - w D)X = sB
!>
!> where A is an NA by NA matrix (NA=1 or 2 only), w is a real (NW=1) or
!> complex (NW=2) constant, ca is a real constant, D is an NA by NA real
!> diagonal matrix, and B is an NA by NW matrix (when NW=2 the second
!> column of B contains the imaginary part of the solution).  The code
!> returns X and s, where s is a scale factor, less than or equal to 1,
!> which is chosen to avoid overflow in X.
!>
!> If any singular values of ca A-w D are less than another input
!> parameter SMIN, they are perturbed up to SMIN.
!>
!> The test condition is that the scaled residual
!>
!>     norm( (ca A-w D)*X - s*B ) /
!>           ( max( ulp*norm(ca A-w D), SMIN )*norm(X) )
!>
!> should be on the order of 1.  Here, ulp is the machine precision.
!> Also, it is verified that SCALE is less than or equal to 1, and that
!> XNORM = infinity-norm(X).
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
!>          NINFO is INTEGER array, dimension (3)
!>          NINFO(1) = number of examples with INFO less than 0
!>          NINFO(2) = number of examples with INFO greater than 0
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
      SUBROUTINE DGET31(Rmax,Lmax,Ninfo,Knt)
      IMPLICIT NONE
!*--DGET3195
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
      DOUBLE PRECISION TWO , THREE , FOUR
      PARAMETER (TWO=2.0D0,THREE=3.0D0,FOUR=4.0D0)
      DOUBLE PRECISION SEVEN , TEN
      PARAMETER (SEVEN=7.0D0,TEN=10.0D0)
      DOUBLE PRECISION TWNONE
      PARAMETER (TWNONE=21.0D0)
!     ..
!     .. Local Scalars ..
      INTEGER ia , ib , ica , id1 , id2 , info , ismin , itrans , iwi , &
     &        iwr , na , nw
      DOUBLE PRECISION bignum , ca , d1 , d2 , den , eps , res , scale ,&
     &                 smin , smlnum , tmp , unfl , wi , wr , xnorm
!     ..
!     .. Local Arrays ..
      LOGICAL ltrans(0:1)
      DOUBLE PRECISION a(2,2) , b(2,2) , vab(3) , vca(5) , vdd(4) ,     &
     &                 vsmin(4) , vwi(4) , vwr(4) , x(2,2)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH
      EXTERNAL DLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL DLABAD , DLALN2
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , SQRT
!     ..
!     .. Data statements ..
      DATA ltrans/.FALSE. , .TRUE./
!     ..
!     .. Executable Statements ..
!
!     Get machine parameters
!
      eps = DLAMCH('P')
      unfl = DLAMCH('U')
      smlnum = DLAMCH('S')/eps
      bignum = ONE/smlnum
      CALL DLABAD(smlnum,bignum)
!
!     Set up test case parameters
!
      vsmin(1) = smlnum
      vsmin(2) = eps
      vsmin(3) = ONE/(TEN*TEN)
      vsmin(4) = ONE/eps
      vab(1) = SQRT(smlnum)
      vab(2) = ONE
      vab(3) = SQRT(bignum)
      vwr(1) = ZERO
      vwr(2) = HALF
      vwr(3) = TWO
      vwr(4) = ONE
      vwi(1) = smlnum
      vwi(2) = eps
      vwi(3) = ONE
      vwi(4) = TWO
      vdd(1) = SQRT(smlnum)
      vdd(2) = ONE
      vdd(3) = TWO
      vdd(4) = SQRT(bignum)
      vca(1) = ZERO
      vca(2) = SQRT(smlnum)
      vca(3) = eps
      vca(4) = HALF
      vca(5) = ONE
!
      Knt = 0
      Ninfo(1) = 0
      Ninfo(2) = 0
      Lmax = 0
      Rmax = ZERO
!
!     Begin test loop
!
      DO id1 = 1 , 4
         d1 = vdd(id1)
         DO id2 = 1 , 4
            d2 = vdd(id2)
            DO ica = 1 , 5
               ca = vca(ica)
               DO itrans = 0 , 1
                  DO ismin = 1 , 4
                     smin = vsmin(ismin)
!
                     na = 1
                     nw = 1
                     DO ia = 1 , 3
                        a(1,1) = vab(ia)
                        DO ib = 1 , 3
                           b(1,1) = vab(ib)
                           DO iwr = 1 , 4
                              IF ( d1==ONE .AND. d2==ONE .AND. ca==ONE )&
     &                             THEN
                                 wr = vwr(iwr)*a(1,1)
                              ELSE
                                 wr = vwr(iwr)
                              ENDIF
                              wi = ZERO
                              CALL DLALN2(ltrans(itrans),na,nw,smin,ca, &
     &                           a,2,d1,d2,b,2,wr,wi,x,2,scale,xnorm,   &
     &                           info)
                              IF ( info<0 ) Ninfo(1) = Ninfo(1) + 1
                              IF ( info>0 ) Ninfo(2) = Ninfo(2) + 1
                              res = ABS((ca*a(1,1)-wr*d1)*x(1,1)        &
     &                              -scale*b(1,1))
                              IF ( info==0 ) THEN
                                 den = MAX                              &
     &                                 (eps*(ABS((ca*a(1,1)-wr*d1)*x(1, &
     &                                 1))),smlnum)
                              ELSE
                                 den = MAX(smin*ABS(x(1,1)),smlnum)
                              ENDIF
                              res = res/den
                              IF ( ABS(x(1,1))<unfl .AND. ABS(b(1,1))   &
     &                             <=smlnum*ABS(ca*a(1,1)-wr*d1) )      &
     &                             res = ZERO
                              IF ( scale>ONE ) res = res + ONE/eps
                              res = res + ABS(xnorm-ABS(x(1,1)))        &
     &                              /MAX(smlnum,xnorm)/eps
                              IF ( info/=0 .AND. info/=1 ) res = res +  &
     &                             ONE/eps
                              Knt = Knt + 1
                              IF ( res>Rmax ) THEN
                                 Lmax = Knt
                                 Rmax = res
                              ENDIF
                           ENDDO
                        ENDDO
                     ENDDO
!
                     na = 1
                     nw = 2
                     DO ia = 1 , 3
                        a(1,1) = vab(ia)
                        DO ib = 1 , 3
                           b(1,1) = vab(ib)
                           b(1,2) = -HALF*vab(ib)
                           DO iwr = 1 , 4
                              IF ( d1==ONE .AND. d2==ONE .AND. ca==ONE )&
     &                             THEN
                                 wr = vwr(iwr)*a(1,1)
                              ELSE
                                 wr = vwr(iwr)
                              ENDIF
                              DO iwi = 1 , 4
                                 IF ( d1==ONE .AND. d2==ONE .AND.       &
     &                                ca==ONE ) THEN
                                    wi = vwi(iwi)*a(1,1)
                                 ELSE
                                    wi = vwi(iwi)
                                 ENDIF
                                 CALL DLALN2(ltrans(itrans),na,nw,smin, &
     &                              ca,a,2,d1,d2,b,2,wr,wi,x,2,scale,   &
     &                              xnorm,info)
                                 IF ( info<0 ) Ninfo(1) = Ninfo(1) + 1
                                 IF ( info>0 ) Ninfo(2) = Ninfo(2) + 1
                                 res = ABS((ca*a(1,1)-wr*d1)*x(1,1)     &
     &                                 +(wi*d1)*x(1,2)-scale*b(1,1))
                                 res = res +                            &
     &                                 ABS((-wi*d1)*x(1,1)+(ca*a(1,1)   &
     &                                 -wr*d1)*x(1,2)-scale*b(1,2))
                                 IF ( info==0 ) THEN
                                    den = MAX                           &
     &                                 (eps*(MAX(ABS(ca*a(1,1)-wr*d1),  &
     &                                 ABS(d1*wi))                      &
     &                                 *(ABS(x(1,1))+ABS(x(1,2)))),     &
     &                                 smlnum)
                                 ELSE
                                    den = MAX                           &
     &                                 (smin*(ABS(x(1,1))+ABS(x(1,2))), &
     &                                 smlnum)
                                 ENDIF
                                 res = res/den
                                 IF ( ABS(x(1,1))<unfl .AND. ABS(x(1,2))&
     &                                <unfl .AND. ABS(b(1,1))           &
     &                                <=smlnum*ABS(ca*a(1,1)-wr*d1) )   &
     &                                res = ZERO
                                 IF ( scale>ONE ) res = res + ONE/eps
                                 res = res +                            &
     &                                 ABS(xnorm-ABS(x(1,1))-ABS(x(1,2))&
     &                                 )/MAX(smlnum,xnorm)/eps
                                 IF ( info/=0 .AND. info/=1 )           &
     &                                res = res + ONE/eps
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
                     na = 2
                     nw = 1
                     DO ia = 1 , 3
                        a(1,1) = vab(ia)
                        a(1,2) = -THREE*vab(ia)
                        a(2,1) = -SEVEN*vab(ia)
                        a(2,2) = TWNONE*vab(ia)
                        DO ib = 1 , 3
                           b(1,1) = vab(ib)
                           b(2,1) = -TWO*vab(ib)
                           DO iwr = 1 , 4
                              IF ( d1==ONE .AND. d2==ONE .AND. ca==ONE )&
     &                             THEN
                                 wr = vwr(iwr)*a(1,1)
                              ELSE
                                 wr = vwr(iwr)
                              ENDIF
                              wi = ZERO
                              CALL DLALN2(ltrans(itrans),na,nw,smin,ca, &
     &                           a,2,d1,d2,b,2,wr,wi,x,2,scale,xnorm,   &
     &                           info)
                              IF ( info<0 ) Ninfo(1) = Ninfo(1) + 1
                              IF ( info>0 ) Ninfo(2) = Ninfo(2) + 1
                              IF ( itrans==1 ) THEN
                                 tmp = a(1,2)
                                 a(1,2) = a(2,1)
                                 a(2,1) = tmp
                              ENDIF
                              res = ABS((ca*a(1,1)-wr*d1)*x(1,1)        &
     &                              +(ca*a(1,2))*x(2,1)-scale*b(1,1))
                              res = res +                               &
     &                              ABS((ca*a(2,1))*x(1,1)+(ca*a(2,2)   &
     &                              -wr*d2)*x(2,1)-scale*b(2,1))
                              IF ( info==0 ) THEN
                                 den = MAX                              &
     &                                 (eps*(MAX(ABS(ca*a(1,1)-wr*d1)+  &
     &                                 ABS(ca*a(1,2)),ABS(ca*a(2,1))    &
     &                                 +ABS(ca*a(2,2)-wr*d2))           &
     &                                 *MAX(ABS(x(1,1)),ABS(x(2,1)))),  &
     &                                 smlnum)
                              ELSE
                                 den = MAX                              &
     &                                 (eps*(MAX(smin/eps,MAX(ABS(ca*a  &
     &                                 (1,1)-wr*d1)+ABS(ca*a(1,2)),     &
     &                                 ABS(ca*a(2,1))                   &
     &                                 +ABS(ca*a(2,2)-wr*d2)))          &
     &                                 *MAX(ABS(x(1,1)),ABS(x(2,1)))),  &
     &                                 smlnum)
                              ENDIF
                              res = res/den
                              IF ( ABS(x(1,1))<unfl .AND. ABS(x(2,1))   &
     &                             <unfl .AND. ABS(b(1,1))+ABS(b(2,1))  &
     &                             <=smlnum*(ABS(ca*a(1,1)-wr*d1)       &
     &                             +ABS(ca*a(1,2))+ABS(ca*a(2,1))       &
     &                             +ABS(ca*a(2,2)-wr*d2)) ) res = ZERO
                              IF ( scale>ONE ) res = res + ONE/eps
                              res = res +                               &
     &                              ABS(xnorm-MAX(ABS(x(1,1)),ABS(x(2,1)&
     &                              )))/MAX(smlnum,xnorm)/eps
                              IF ( info/=0 .AND. info/=1 ) res = res +  &
     &                             ONE/eps
                              Knt = Knt + 1
                              IF ( res>Rmax ) THEN
                                 Lmax = Knt
                                 Rmax = res
                              ENDIF
                           ENDDO
                        ENDDO
                     ENDDO
!
                     na = 2
                     nw = 2
                     DO ia = 1 , 3
                        a(1,1) = vab(ia)*TWO
                        a(1,2) = -THREE*vab(ia)
                        a(2,1) = -SEVEN*vab(ia)
                        a(2,2) = TWNONE*vab(ia)
                        DO ib = 1 , 3
                           b(1,1) = vab(ib)
                           b(2,1) = -TWO*vab(ib)
                           b(1,2) = FOUR*vab(ib)
                           b(2,2) = -SEVEN*vab(ib)
                           DO iwr = 1 , 4
                              IF ( d1==ONE .AND. d2==ONE .AND. ca==ONE )&
     &                             THEN
                                 wr = vwr(iwr)*a(1,1)
                              ELSE
                                 wr = vwr(iwr)
                              ENDIF
                              DO iwi = 1 , 4
                                 IF ( d1==ONE .AND. d2==ONE .AND.       &
     &                                ca==ONE ) THEN
                                    wi = vwi(iwi)*a(1,1)
                                 ELSE
                                    wi = vwi(iwi)
                                 ENDIF
                                 CALL DLALN2(ltrans(itrans),na,nw,smin, &
     &                              ca,a,2,d1,d2,b,2,wr,wi,x,2,scale,   &
     &                              xnorm,info)
                                 IF ( info<0 ) Ninfo(1) = Ninfo(1) + 1
                                 IF ( info>0 ) Ninfo(2) = Ninfo(2) + 1
                                 IF ( itrans==1 ) THEN
                                    tmp = a(1,2)
                                    a(1,2) = a(2,1)
                                    a(2,1) = tmp
                                 ENDIF
                                 res = ABS((ca*a(1,1)-wr*d1)*x(1,1)     &
     &                                 +(ca*a(1,2))*x(2,1)+(wi*d1)      &
     &                                 *x(1,2)-scale*b(1,1))
                                 res = res +                            &
     &                                 ABS((ca*a(1,1)-wr*d1)*x(1,2)     &
     &                                 +(ca*a(1,2))*x(2,2)-(wi*d1)      &
     &                                 *x(1,1)-scale*b(1,2))
                                 res = res +                            &
     &                                 ABS((ca*a(2,1))*x(1,1)+(ca*a(2,2)&
     &                                 -wr*d2)*x(2,1)+(wi*d2)*x(2,2)    &
     &                                 -scale*b(2,1))
                                 res = res +                            &
     &                                 ABS((ca*a(2,1))*x(1,2)+(ca*a(2,2)&
     &                                 -wr*d2)*x(2,2)-(wi*d2)*x(2,1)    &
     &                                 -scale*b(2,2))
                                 IF ( info==0 ) THEN
                                    den = MAX                           &
     &                                 (eps*(MAX(ABS(ca*a(1,1)-wr*d1)   &
     &                                 +ABS(ca*a(1,2))+ABS(wi*d1),      &
     &                                 ABS(ca*a(2,1))                   &
     &                                 +ABS(ca*a(2,2)-wr*d2)+ABS(wi*d2))&
     &                                 *MAX(ABS(x(1,1))+ABS(x(2,1)),    &
     &                                 ABS(x(1,2))+ABS(x(2,2)))),smlnum)
                                 ELSE
                                    den = MAX                           &
     &                                 (eps*(MAX(smin/eps,MAX(ABS(ca*a  &
     &                                 (1,1)-wr*d1)+ABS(ca*a(1,2))      &
     &                                 +ABS(wi*d1),ABS(ca*a(2,1))       &
     &                                 +ABS(ca*a(2,2)-wr*d2)+ABS(wi*d2))&
     &                                 )                                &
     &                                 *MAX(ABS(x(1,1))+ABS(x(2,1)),ABS(&
     &                                 x(1,2))+ABS(x(2,2)))),smlnum)
                                 ENDIF
                                 res = res/den
                                 IF ( ABS(x(1,1))<unfl .AND. ABS(x(2,1))&
     &                                <unfl .AND. ABS(x(1,2))<unfl .AND.&
     &                                ABS(x(2,2))<unfl .AND. ABS(b(1,1))&
     &                                +ABS(b(2,1))                      &
     &                                <=smlnum*(ABS(ca*a(1,1)-wr*d1)    &
     &                                +ABS(ca*a(1,2))+ABS(ca*a(2,1))    &
     &                                +ABS(ca*a(2,2)-wr*d2)+ABS(wi*d2)  &
     &                                +ABS(wi*d1)) ) res = ZERO
                                 IF ( scale>ONE ) res = res + ONE/eps
                                 res = res +                            &
     &                                 ABS(xnorm-MAX(ABS(x(1,1))+ABS    &
     &                                 (x(1,2)),ABS(x(2,1))+ABS(x(2,2)))&
     &                                 )/MAX(smlnum,xnorm)/eps
                                 IF ( info/=0 .AND. info/=1 )           &
     &                                res = res + ONE/eps
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
!
!
!     End of DGET31
!
      END SUBROUTINE DGET31
