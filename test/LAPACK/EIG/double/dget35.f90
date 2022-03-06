!*==dget35.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b DGET35
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGET35( RMAX, LMAX, NINFO, KNT )
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
!> DGET35 tests DTRSYL, a routine for solving the Sylvester matrix
!> equation
!>
!>    op(A)*X + ISGN*X*op(B) = scale*C,
!>
!> A and B are assumed to be in Schur canonical form, op() represents an
!> optional transpose, and ISGN can be -1 or +1.  Scale is an output
!> less than or equal to 1, chosen to avoid overflow in X.
!>
!> The test code verifies that the following residual is order 1:
!>
!>    norm(op(A)*X + ISGN*X*op(B) - scale*C) /
!>        (EPS*max(norm(A),norm(B))*norm(X))
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
      SUBROUTINE DGET35(Rmax,Lmax,Ninfo,Knt)
      IMPLICIT NONE
!*--DGET3582
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
      DOUBLE PRECISION TWO , FOUR
      PARAMETER (TWO=2.0D0,FOUR=4.0D0)
!     ..
!     .. Local Scalars ..
      CHARACTER trana , tranb
      INTEGER i , ima , imb , imlda1 , imlda2 , imldb1 , imloff , info ,&
     &        isgn , itrana , itranb , j , m , n
      DOUBLE PRECISION bignum , cnrm , eps , res , res1 , rmul , scale ,&
     &                 smlnum , tnrm , xnrm
!     ..
!     .. Local Arrays ..
      INTEGER idim(8) , ival(6,6,8)
      DOUBLE PRECISION a(6,6) , b(6,6) , c(6,6) , cc(6,6) , dum(1) ,    &
     &                 vm1(3) , vm2(3)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH , DLANGE
      EXTERNAL DLAMCH , DLANGE
!     ..
!     .. External Subroutines ..
      EXTERNAL DGEMM , DLABAD , DTRSYL
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , DBLE , MAX , SIN , SQRT
!     ..
!     .. Data statements ..
      DATA idim/1 , 2 , 3 , 4 , 3 , 3 , 6 , 4/
      DATA ival/1 , 35*0 , 1 , 2 , 4*0 , -2 , 0 , 28*0 , 1 , 5*0 , 5 ,  &
     &     1 , 2 , 3*0 , -8 , -2 , 1 , 21*0 , 3 , 4 , 4*0 , -5 , 3 ,    &
     &     4*0 , 1 , 2 , 1 , 4 , 2*0 , -3 , -9 , -1 , 1 , 14*0 , 1 ,    &
     &     5*0 , 2 , 3 , 4*0 , 5 , 6 , 7 , 21*0 , 1 , 5*0 , 1 , 3 , -4 ,&
     &     3*0 , 2 , 5 , 2 , 21*0 , 1 , 2 , 4*0 , -2 , 0 , 4*0 , 5 , 6 ,&
     &     3 , 4 , 2*0 , -1 , -9 , -5 , 2 , 2*0 , 4*8 , 5 , 6 , 4*9 ,   &
     &     -7 , 5 , 1 , 5*0 , 1 , 5 , 2 , 3*0 , 2 , -21 , 5 , 3*0 , 1 , &
     &     2 , 3 , 4 , 14*0/
!     ..
!     .. Executable Statements ..
!
!     Get machine parameters
!
      eps = DLAMCH('P')
      smlnum = DLAMCH('S')*FOUR/eps
      bignum = ONE/smlnum
      CALL DLABAD(smlnum,bignum)
!
!     Set up test case parameters
!
      vm1(1) = SQRT(smlnum)
      vm1(2) = ONE
      vm1(3) = SQRT(bignum)
      vm2(1) = ONE
      vm2(2) = ONE + TWO*eps
      vm2(3) = TWO
!
      Knt = 0
      Ninfo = 0
      Lmax = 0
      Rmax = ZERO
!
!     Begin test loop
!
      DO itrana = 1 , 2
         DO itranb = 1 , 2
            DO isgn = -1 , 1 , 2
               DO ima = 1 , 8
                  DO imlda1 = 1 , 3
                     DO imlda2 = 1 , 3
                        DO imloff = 1 , 2
                           DO imb = 1 , 8
                              DO imldb1 = 1 , 3
                                 IF ( itrana==1 ) trana = 'N'
                                 IF ( itrana==2 ) trana = 'T'
                                 IF ( itranb==1 ) tranb = 'N'
                                 IF ( itranb==2 ) tranb = 'T'
                                 m = idim(ima)
                                 n = idim(imb)
                                 tnrm = ZERO
                                 DO i = 1 , m
                                    DO j = 1 , m
                                       a(i,j) = ival(i,j,ima)
                                       IF ( ABS(i-j)<=1 ) THEN
                                         a(i,j) = a(i,j)*vm1(imlda1)
                                         a(i,j) = a(i,j)*vm2(imlda2)
                                       ELSE
                                         a(i,j) = a(i,j)*vm1(imloff)
                                       ENDIF
                                       tnrm = MAX(tnrm,ABS(a(i,j)))
                                    ENDDO
                                 ENDDO
                                 DO i = 1 , n
                                    DO j = 1 , n
                                       b(i,j) = ival(i,j,imb)
                                       IF ( ABS(i-j)<=1 ) THEN
                                         b(i,j) = b(i,j)*vm1(imldb1)
                                       ELSE
                                         b(i,j) = b(i,j)*vm1(imloff)
                                       ENDIF
                                       tnrm = MAX(tnrm,ABS(b(i,j)))
                                    ENDDO
                                 ENDDO
                                 cnrm = ZERO
                                 DO i = 1 , m
                                    DO j = 1 , n
                                       c(i,j) = SIN(DBLE(i*j))
                                       cnrm = MAX(cnrm,c(i,j))
                                       cc(i,j) = c(i,j)
                                    ENDDO
                                 ENDDO
                                 Knt = Knt + 1
                                 CALL DTRSYL(trana,tranb,isgn,m,n,a,6,b,&
     &                              6,c,6,scale,info)
                                 IF ( info/=0 ) Ninfo = Ninfo + 1
                                 xnrm = DLANGE('M',m,n,c,6,dum)
                                 rmul = ONE
                                 IF ( xnrm>ONE .AND. tnrm>ONE ) THEN
                                    IF ( xnrm>bignum/tnrm )             &
     &                                 rmul = ONE/MAX(xnrm,tnrm)
                                 ENDIF
                                 CALL DGEMM(trana,'N',m,n,m,rmul,a,6,c, &
     &                              6,-scale*rmul,cc,6)
                                 CALL DGEMM('N',tranb,m,n,n,DBLE(isgn)  &
     &                              *rmul,c,6,b,6,ONE,cc,6)
                                 res1 = DLANGE('M',m,n,cc,6,dum)
                                 res = res1/MAX(smlnum,smlnum*xnrm,     &
     &                                 ((rmul*tnrm)*eps)*xnrm)
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
!     End of DGET35
!
      END SUBROUTINE DGET35
