!*==dget39.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b DGET39
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGET39( RMAX, LMAX, NINFO, KNT )
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
!> DGET39 tests DLAQTR, a routine for solving the real or
!> special complex quasi upper triangular system
!>
!>      op(T)*p = scale*c,
!> or
!>      op(T + iB)*(p+iq) = scale*(c+id),
!>
!> in real arithmetic. T is upper quasi-triangular.
!> If it is complex, then the first diagonal block of T must be
!> 1 by 1, B has the special structure
!>
!>                B = [ b(1) b(2) ... b(n) ]
!>                    [       w            ]
!>                    [           w        ]
!>                    [              .     ]
!>                    [                 w  ]
!>
!> op(A) = A or A', where A' denotes the conjugate transpose of
!> the matrix A.
!>
!> On input, X = [ c ].  On output, X = [ p ].
!>               [ d ]                  [ q ]
!>
!> Scale is an output less than or equal to 1, chosen to avoid
!> overflow in X.
!> This subroutine is specially designed for the condition number
!> estimation in the eigenproblem routine DTRSNA.
!>
!> The test code verifies that the following residual is order 1:
!>
!>      ||(T+i*B)*(x1+i*x2) - scale*(d1+i*d2)||
!>    -----------------------------------------
!>        max(ulp*(||T||+||B||)*(||x1||+||x2||),
!>            (||T||+||B||)*smlnum/ulp,
!>            smlnum)
!>
!> (The (||T||+||B||)*smlnum/ulp term accounts for possible
!>  (gradual or nongradual) underflow in x1 and x2.)
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
!>          Number of examples where INFO is nonzero.
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
      SUBROUTINE DGET39(Rmax,Lmax,Ninfo,Knt)
      IMPLICIT NONE
!*--DGET39107
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
      INTEGER LDT , LDT2
      PARAMETER (LDT=10,LDT2=2*LDT)
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
!     ..
!     .. Local Scalars ..
      INTEGER i , info , ivm1 , ivm2 , ivm3 , ivm4 , ivm5 , j , k , n , &
     &        ndim
      DOUBLE PRECISION bignum , domin , dumm , eps , norm , normtb ,    &
     &                 resid , scale , smlnum , w , xnorm
!     ..
!     .. External Functions ..
      INTEGER IDAMAX
      DOUBLE PRECISION DASUM , DDOT , DLAMCH , DLANGE
      EXTERNAL IDAMAX , DASUM , DDOT , DLAMCH , DLANGE
!     ..
!     .. External Subroutines ..
      EXTERNAL DCOPY , DGEMV , DLABAD , DLAQTR
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , COS , DBLE , MAX , SIN , SQRT
!     ..
!     .. Local Arrays ..
      INTEGER idim(6) , ival(5,5,6)
      DOUBLE PRECISION b(LDT) , d(LDT2) , dum(1) , t(LDT,LDT) , vm1(5) ,&
     &                 vm2(5) , vm3(5) , vm4(5) , vm5(3) , work(LDT) ,  &
     &                 x(LDT2) , y(LDT2)
!     ..
!     .. Data statements ..
      DATA idim/4 , 5*5/
      DATA ival/3 , 4*0 , 1 , 1 , -1 , 0 , 0 , 3 , 2 , 1 , 0 , 0 , 4 ,  &
     &     3 , 2 , 2 , 0 , 5*0 , 1 , 4*0 , 2 , 2 , 3*0 , 3 , 3 , 4 , 0 ,&
     &     0 , 4 , 2 , 2 , 3 , 0 , 4*1 , 5 , 1 , 4*0 , 2 , 4 , -2 , 0 , &
     &     0 , 3 , 3 , 4 , 0 , 0 , 4 , 2 , 2 , 3 , 0 , 5*1 , 1 , 4*0 ,  &
     &     2 , 1 , -1 , 0 , 0 , 9 , 8 , 1 , 0 , 0 , 4 , 9 , 1 , 2 , -1 ,&
     &     5*2 , 9 , 4*0 , 6 , 4 , 0 , 0 , 0 , 3 , 2 , 1 , 1 , 0 , 5 ,  &
     &     1 , -1 , 1 , 0 , 5*2 , 4 , 4*0 , 2 , 2 , 0 , 0 , 0 , 1 , 4 , &
     &     4 , 0 , 0 , 2 , 4 , 2 , 2 , -1 , 5*2/
!     ..
!     .. Executable Statements ..
!
!     Get machine parameters
!
      eps = DLAMCH('P')
      smlnum = DLAMCH('S')
      bignum = ONE/smlnum
      CALL DLABAD(smlnum,bignum)
!
!     Set up test case parameters
!
      vm1(1) = ONE
      vm1(2) = SQRT(smlnum)
      vm1(3) = SQRT(vm1(2))
      vm1(4) = SQRT(bignum)
      vm1(5) = SQRT(vm1(4))
!
      vm2(1) = ONE
      vm2(2) = SQRT(smlnum)
      vm2(3) = SQRT(vm2(2))
      vm2(4) = SQRT(bignum)
      vm2(5) = SQRT(vm2(4))
!
      vm3(1) = ONE
      vm3(2) = SQRT(smlnum)
      vm3(3) = SQRT(vm3(2))
      vm3(4) = SQRT(bignum)
      vm3(5) = SQRT(vm3(4))
!
      vm4(1) = ONE
      vm4(2) = SQRT(smlnum)
      vm4(3) = SQRT(vm4(2))
      vm4(4) = SQRT(bignum)
      vm4(5) = SQRT(vm4(4))
!
      vm5(1) = ONE
      vm5(2) = eps
      vm5(3) = SQRT(smlnum)
!
!     Initialization
!
      Knt = 0
      Rmax = ZERO
      Ninfo = 0
      smlnum = smlnum/eps
!
!     Begin test loop
!
      DO ivm5 = 1 , 3
         DO ivm4 = 1 , 5
            DO ivm3 = 1 , 5
               DO ivm2 = 1 , 5
                  DO ivm1 = 1 , 5
                     DO ndim = 1 , 6
!
                        n = idim(ndim)
                        DO i = 1 , n
                           DO j = 1 , n
                              t(i,j) = DBLE(ival(i,j,ndim))*vm1(ivm1)
                              IF ( i>=j ) t(i,j) = t(i,j)*vm5(ivm5)
                           ENDDO
                        ENDDO
!
                        w = ONE*vm2(ivm2)
!
                        DO i = 1 , n
                           b(i) = COS(DBLE(i))*vm3(ivm3)
                        ENDDO
!
                        DO i = 1 , 2*n
                           d(i) = SIN(DBLE(i))*vm4(ivm4)
                        ENDDO
!
                        norm = DLANGE('1',n,n,t,LDT,work)
                        k = IDAMAX(n,b,1)
                        normtb = norm + ABS(b(k)) + ABS(w)
!
                        CALL DCOPY(n,d,1,x,1)
                        Knt = Knt + 1
                        CALL DLAQTR(.FALSE.,.TRUE.,n,t,LDT,dum,dumm,    &
     &                              scale,x,work,info)
                        IF ( info/=0 ) Ninfo = Ninfo + 1
!
!                       || T*x - scale*d || /
!                         max(ulp*||T||*||x||,smlnum/ulp*||T||,smlnum)
!
                        CALL DCOPY(n,d,1,y,1)
                        CALL DGEMV('No transpose',n,n,ONE,t,LDT,x,1,    &
     &                             -scale,y,1)
                        xnorm = DASUM(n,x,1)
                        resid = DASUM(n,y,1)
                        domin = MAX(smlnum,(smlnum/eps)*norm,(norm*eps) &
     &                          *xnorm)
                        resid = resid/domin
                        IF ( resid>Rmax ) THEN
                           Rmax = resid
                           Lmax = Knt
                        ENDIF
!
                        CALL DCOPY(n,d,1,x,1)
                        Knt = Knt + 1
                        CALL DLAQTR(.TRUE.,.TRUE.,n,t,LDT,dum,dumm,     &
     &                              scale,x,work,info)
                        IF ( info/=0 ) Ninfo = Ninfo + 1
!
!                       || T*x - scale*d || /
!                         max(ulp*||T||*||x||,smlnum/ulp*||T||,smlnum)
!
                        CALL DCOPY(n,d,1,y,1)
                        CALL DGEMV('Transpose',n,n,ONE,t,LDT,x,1,-scale,&
     &                             y,1)
                        xnorm = DASUM(n,x,1)
                        resid = DASUM(n,y,1)
                        domin = MAX(smlnum,(smlnum/eps)*norm,(norm*eps) &
     &                          *xnorm)
                        resid = resid/domin
                        IF ( resid>Rmax ) THEN
                           Rmax = resid
                           Lmax = Knt
                        ENDIF
!
                        CALL DCOPY(2*n,d,1,x,1)
                        Knt = Knt + 1
                        CALL DLAQTR(.FALSE.,.FALSE.,n,t,LDT,b,w,scale,x,&
     &                              work,info)
                        IF ( info/=0 ) Ninfo = Ninfo + 1
!
!                       ||(T+i*B)*(x1+i*x2) - scale*(d1+i*d2)|| /
!                          max(ulp*(||T||+||B||)*(||x1||+||x2||),
!                                  smlnum/ulp * (||T||+||B||), smlnum )
!
!
                        CALL DCOPY(2*n,d,1,y,1)
                        y(1) = DDOT(n,b,1,x(1+n),1) + scale*y(1)
                        DO i = 2 , n
                           y(i) = w*x(i+n) + scale*y(i)
                        ENDDO
                        CALL DGEMV('No transpose',n,n,ONE,t,LDT,x,1,    &
     &                             -ONE,y,1)
!
                        y(1+n) = DDOT(n,b,1,x,1) - scale*y(1+n)
                        DO i = 2 , n
                           y(i+n) = w*x(i) - scale*y(i+n)
                        ENDDO
                        CALL DGEMV('No transpose',n,n,ONE,t,LDT,x(1+n), &
     &                             1,ONE,y(1+n),1)
!
                        resid = DASUM(2*n,y,1)
                        domin = MAX(smlnum,(smlnum/eps)*normtb,         &
     &                          eps*(normtb*DASUM(2*n,x,1)))
                        resid = resid/domin
                        IF ( resid>Rmax ) THEN
                           Rmax = resid
                           Lmax = Knt
                        ENDIF
!
                        CALL DCOPY(2*n,d,1,x,1)
                        Knt = Knt + 1
                        CALL DLAQTR(.TRUE.,.FALSE.,n,t,LDT,b,w,scale,x, &
     &                              work,info)
                        IF ( info/=0 ) Ninfo = Ninfo + 1
!
!                       ||(T+i*B)*(x1+i*x2) - scale*(d1+i*d2)|| /
!                          max(ulp*(||T||+||B||)*(||x1||+||x2||),
!                                  smlnum/ulp * (||T||+||B||), smlnum )
!
                        CALL DCOPY(2*n,d,1,y,1)
                        y(1) = b(1)*x(1+n) - scale*y(1)
                        DO i = 2 , n
                           y(i) = b(i)*x(1+n) + w*x(i+n) - scale*y(i)
                        ENDDO
                        CALL DGEMV('Transpose',n,n,ONE,t,LDT,x,1,ONE,y, &
     &                             1)
!
                        y(1+n) = b(1)*x(1) + scale*y(1+n)
                        DO i = 2 , n
                           y(i+n) = b(i)*x(1) + w*x(i) + scale*y(i+n)
                        ENDDO
                        CALL DGEMV('Transpose',n,n,ONE,t,LDT,x(1+n),1,  &
     &                             -ONE,y(1+n),1)
!
                        resid = DASUM(2*n,y,1)
                        domin = MAX(smlnum,(smlnum/eps)*normtb,         &
     &                          eps*(normtb*DASUM(2*n,x,1)))
                        resid = resid/domin
                        IF ( resid>Rmax ) THEN
                           Rmax = resid
                           Lmax = Knt
                        ENDIF
!
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
!
!
!     End of DGET39
!
      END SUBROUTINE DGET39
