!*==dget32.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b DGET32
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGET32( RMAX, LMAX, NINFO, KNT )
!
!       .. Scalar Arguments ..
!       INTEGER            KNT, LMAX, NINFO
!       DOUBLE PRECISION   RMAX
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGET32 tests DLASY2, a routine for solving
!>
!>         op(TL)*X + ISGN*X*op(TR) = SCALE*B
!>
!> where TL is N1 by N1, TR is N2 by N2, and N1,N2 =1 or 2 only.
!> X and B are N1 by N2, op() is an optional transpose, an
!> ISGN = 1 or -1. SCALE is chosen less than or equal to 1 to
!> avoid overflow in X.
!>
!> The test condition is that the scaled residual
!>
!> norm( op(TL)*X + ISGN*X*op(TR) = SCALE*B )
!>      / ( max( ulp*norm(TL), ulp*norm(TR)) * norm(X), SMLNUM )
!>
!> should be on the order of 1. Here, ulp is the machine precision.
!> Also, it is verified that SCALE is less than or equal to 1, and
!> that XNORM = infinity-norm(X).
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
!>          NINFO is INTEGER
!>          Number of examples returned with INFO.NE.0.
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
      SUBROUTINE DGET32(Rmax,Lmax,Ninfo,Knt)
      IMPLICIT NONE
!*--DGET3286
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Knt , Lmax , Ninfo
      DOUBLE PRECISION Rmax
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
      DOUBLE PRECISION TWO , FOUR , EIGHT
      PARAMETER (TWO=2.0D0,FOUR=4.0D0,EIGHT=8.0D0)
!     ..
!     .. Local Scalars ..
      LOGICAL ltranl , ltranr
      INTEGER ib , ib1 , ib2 , ib3 , info , isgn , itl , itlscl , itr , &
     &        itranl , itranr , itrscl , n1 , n2
      DOUBLE PRECISION bignum , den , eps , res , scale , sgn , smlnum ,&
     &                 tmp , tnrm , xnorm , xnrm
!     ..
!     .. Local Arrays ..
      INTEGER itval(2,2,8)
      DOUBLE PRECISION b(2,2) , tl(2,2) , tr(2,2) , val(3) , x(2,2)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH
      EXTERNAL DLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL DLABAD , DLASY2
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , MIN , SQRT
!     ..
!     .. Data statements ..
      DATA itval/8 , 4 , 2 , 1 , 4 , 8 , 1 , 2 , 2 , 1 , 8 , 4 , 1 , 2 ,&
     &     4 , 8 , 9 , 4 , 2 , 1 , 4 , 9 , 1 , 2 , 2 , 1 , 9 , 4 , 1 ,  &
     &     2 , 4 , 9/
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
      val(1) = SQRT(smlnum)
      val(2) = ONE
      val(3) = SQRT(bignum)
!
      Knt = 0
      Ninfo = 0
      Lmax = 0
      Rmax = ZERO
!
!     Begin test loop
!
      DO itranl = 0 , 1
         DO itranr = 0 , 1
            DO isgn = -1 , 1 , 2
               sgn = isgn
               ltranl = itranl==1
               ltranr = itranr==1
!
               n1 = 1
               n2 = 1
               DO itl = 1 , 3
                  DO itr = 1 , 3
                     DO ib = 1 , 3
                        tl(1,1) = val(itl)
                        tr(1,1) = val(itr)
                        b(1,1) = val(ib)
                        Knt = Knt + 1
                        CALL DLASY2(ltranl,ltranr,isgn,n1,n2,tl,2,tr,2, &
     &                              b,2,scale,x,2,xnorm,info)
                        IF ( info/=0 ) Ninfo = Ninfo + 1
                        res = ABS((tl(1,1)+sgn*tr(1,1))*x(1,1)          &
     &                        -scale*b(1,1))
                        IF ( info==0 ) THEN
                           den = MAX                                    &
     &                           (eps*((ABS(tr(1,1))+ABS(tl(1,1)))*ABS  &
     &                           (x(1,1))),smlnum)
                        ELSE
                           den = smlnum*MAX(ABS(x(1,1)),ONE)
                        ENDIF
                        res = res/den
                        IF ( scale>ONE ) res = res + ONE/eps
                        res = res + ABS(xnorm-ABS(x(1,1)))              &
     &                        /MAX(smlnum,xnorm)/eps
                        IF ( info/=0 .AND. info/=1 ) res = res + ONE/eps
                        IF ( res>Rmax ) THEN
                           Lmax = Knt
                           Rmax = res
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO
!
               n1 = 2
               n2 = 1
               DO itl = 1 , 8
                  DO itlscl = 1 , 3
                     DO itr = 1 , 3
                        DO ib1 = 1 , 3
                           DO ib2 = 1 , 3
                              b(1,1) = val(ib1)
                              b(2,1) = -FOUR*val(ib2)
                              tl(1,1) = itval(1,1,itl)*val(itlscl)
                              tl(2,1) = itval(2,1,itl)*val(itlscl)
                              tl(1,2) = itval(1,2,itl)*val(itlscl)
                              tl(2,2) = itval(2,2,itl)*val(itlscl)
                              tr(1,1) = val(itr)
                              Knt = Knt + 1
                              CALL DLASY2(ltranl,ltranr,isgn,n1,n2,tl,2,&
     &                           tr,2,b,2,scale,x,2,xnorm,info)
                              IF ( info/=0 ) Ninfo = Ninfo + 1
                              IF ( ltranl ) THEN
                                 tmp = tl(1,2)
                                 tl(1,2) = tl(2,1)
                                 tl(2,1) = tmp
                              ENDIF
                              res = ABS((tl(1,1)+sgn*tr(1,1))*x(1,1)    &
     &                              +tl(1,2)*x(2,1)-scale*b(1,1))
                              res = res +                               &
     &                              ABS((tl(2,2)+sgn*tr(1,1))*x(2,1)    &
     &                              +tl(2,1)*x(1,1)-scale*b(2,1))
                              tnrm = ABS(tr(1,1)) + ABS(tl(1,1))        &
     &                               + ABS(tl(1,2)) + ABS(tl(2,1))      &
     &                               + ABS(tl(2,2))
                              xnrm = MAX(ABS(x(1,1)),ABS(x(2,1)))
                              den = MAX(smlnum,smlnum*xnrm,(tnrm*eps)   &
     &                              *xnrm)
                              res = res/den
                              IF ( scale>ONE ) res = res + ONE/eps
                              res = res + ABS(xnorm-xnrm)               &
     &                              /MAX(smlnum,xnorm)/eps
                              IF ( res>Rmax ) THEN
                                 Lmax = Knt
                                 Rmax = res
                              ENDIF
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
!
               n1 = 1
               n2 = 2
               DO itr = 1 , 8
                  DO itrscl = 1 , 3
                     DO itl = 1 , 3
                        DO ib1 = 1 , 3
                           DO ib2 = 1 , 3
                              b(1,1) = val(ib1)
                              b(1,2) = -TWO*val(ib2)
                              tr(1,1) = itval(1,1,itr)*val(itrscl)
                              tr(2,1) = itval(2,1,itr)*val(itrscl)
                              tr(1,2) = itval(1,2,itr)*val(itrscl)
                              tr(2,2) = itval(2,2,itr)*val(itrscl)
                              tl(1,1) = val(itl)
                              Knt = Knt + 1
                              CALL DLASY2(ltranl,ltranr,isgn,n1,n2,tl,2,&
     &                           tr,2,b,2,scale,x,2,xnorm,info)
                              IF ( info/=0 ) Ninfo = Ninfo + 1
                              IF ( ltranr ) THEN
                                 tmp = tr(1,2)
                                 tr(1,2) = tr(2,1)
                                 tr(2,1) = tmp
                              ENDIF
                              tnrm = ABS(tl(1,1)) + ABS(tr(1,1))        &
     &                               + ABS(tr(1,2)) + ABS(tr(2,2))      &
     &                               + ABS(tr(2,1))
                              xnrm = ABS(x(1,1)) + ABS(x(1,2))
                              res = ABS(((tl(1,1)+sgn*tr(1,1)))*(x(1,1))&
     &                              +(sgn*tr(2,1))*(x(1,2))             &
     &                              -(scale*b(1,1)))
                              res = res +                               &
     &                              ABS(((tl(1,1)+sgn*tr(2,2)))*(x(1,2))&
     &                              +(sgn*tr(1,2))*(x(1,1))             &
     &                              -(scale*b(1,2)))
                              den = MAX(smlnum,smlnum*xnrm,(tnrm*eps)   &
     &                              *xnrm)
                              res = res/den
                              IF ( scale>ONE ) res = res + ONE/eps
                              res = res + ABS(xnorm-xnrm)               &
     &                              /MAX(smlnum,xnorm)/eps
                              IF ( res>Rmax ) THEN
                                 Lmax = Knt
                                 Rmax = res
                              ENDIF
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
!
               n1 = 2
               n2 = 2
               DO itr = 1 , 8
                  DO itrscl = 1 , 3
                     DO itl = 1 , 8
                        DO itlscl = 1 , 3
                           DO ib1 = 1 , 3
                              DO ib2 = 1 , 3
                                 DO ib3 = 1 , 3
                                    b(1,1) = val(ib1)
                                    b(2,1) = -FOUR*val(ib2)
                                    b(1,2) = -TWO*val(ib3)
                                    b(2,2)                              &
     &                                 = EIGHT*MIN(val(ib1),val(ib2),   &
     &                                 val(ib3))
                                    tr(1,1) = itval(1,1,itr)*val(itrscl)
                                    tr(2,1) = itval(2,1,itr)*val(itrscl)
                                    tr(1,2) = itval(1,2,itr)*val(itrscl)
                                    tr(2,2) = itval(2,2,itr)*val(itrscl)
                                    tl(1,1) = itval(1,1,itl)*val(itlscl)
                                    tl(2,1) = itval(2,1,itl)*val(itlscl)
                                    tl(1,2) = itval(1,2,itl)*val(itlscl)
                                    tl(2,2) = itval(2,2,itl)*val(itlscl)
                                    Knt = Knt + 1
                                    CALL DLASY2(ltranl,ltranr,isgn,n1,  &
     &                                 n2,tl,2,tr,2,b,2,scale,x,2,xnorm,&
     &                                 info)
                                    IF ( info/=0 ) Ninfo = Ninfo + 1
                                    IF ( ltranr ) THEN
                                       tmp = tr(1,2)
                                       tr(1,2) = tr(2,1)
                                       tr(2,1) = tmp
                                    ENDIF
                                    IF ( ltranl ) THEN
                                       tmp = tl(1,2)
                                       tl(1,2) = tl(2,1)
                                       tl(2,1) = tmp
                                    ENDIF
                                    tnrm = ABS(tr(1,1)) + ABS(tr(2,1))  &
     &                                 + ABS(tr(1,2)) + ABS(tr(2,2))    &
     &                                 + ABS(tl(1,1)) + ABS(tl(2,1))    &
     &                                 + ABS(tl(1,2)) + ABS(tl(2,2))
                                    xnrm = MAX(ABS(x(1,1))+ABS(x(1,2)), &
     &                                 ABS(x(2,1))+ABS(x(2,2)))
                                    res = ABS(((tl(1,1)+sgn*tr(1,1)))   &
     &                                 *(x(1,1))+(sgn*tr(2,1))*(x(1,2)) &
     &                                 +(tl(1,2))*(x(2,1))              &
     &                                 -(scale*b(1,1)))
                                    res = res +                         &
     &                                 ABS((tl(1,1))*(x(1,2))+(sgn*tr(1,&
     &                                 2))*(x(1,1))+(sgn*tr(2,2))       &
     &                                 *(x(1,2))+(tl(1,2))*(x(2,2))     &
     &                                 -(scale*b(1,2)))
                                    res = res +                         &
     &                                 ABS((tl(2,1))*(x(1,1))+(sgn*tr(1,&
     &                                 1))*(x(2,1))+(sgn*tr(2,1))       &
     &                                 *(x(2,2))+(tl(2,2))*(x(2,1))     &
     &                                 -(scale*b(2,1)))
                                    res = res +                         &
     &                                 ABS(((tl(2,2)+sgn*tr(2,2)))      &
     &                                 *(x(2,2))+(sgn*tr(1,2))*(x(2,1)) &
     &                                 +(tl(2,1))*(x(1,2))              &
     &                                 -(scale*b(2,2)))
                                    den = MAX(smlnum,smlnum*xnrm,       &
     &                                 (tnrm*eps)*xnrm)
                                    res = res/den
                                    IF ( scale>ONE ) res = res + ONE/eps
                                    res = res + ABS(xnorm-xnrm)         &
     &                                 /MAX(smlnum,xnorm)/eps
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
!     End of DGET32
!
      END SUBROUTINE DGET32
