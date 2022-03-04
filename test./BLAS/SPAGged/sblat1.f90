!*==sblat1.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
 
!> \brief \b SBLAT1
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       PROGRAM SBLAT1
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    Test program for the REAL Level 1 BLAS.
!>
!>    Based upon the original BLAS test routine together with:
!>    F06EAF Example Program Text
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
!> \date April 2012
!
!> \ingroup single_blas_testing
!
!  =====================================================================
      PROGRAM SBLAT1
      IMPLICIT NONE
!*--SBLAT142
!
!  -- Reference BLAS test routine (version 3.8.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     April 2012
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER NOUT
      PARAMETER (NOUT=6)
!     .. Scalars in Common ..
      INTEGER ICAse , INCx , INCy , N
      LOGICAL PASs
!     .. Local Scalars ..
      REAL sfac
      INTEGER ic
!     .. External Subroutines ..
      EXTERNAL CHECK0 , CHECK1 , CHECK2 , CHECK3 , HEADER
!     .. Common blocks ..
      COMMON /COMBLA/ ICAse , N , INCx , INCy , PASs
!     .. Data statements ..
      DATA sfac/9.765625E-4/
!     .. Executable Statements ..
      WRITE (NOUT,99001)
!
99001 FORMAT (' Real BLAS Test Program Results',/1X)
      DO ic = 1 , 13
         ICAse = ic
         CALL HEADER
!
!        .. Initialize  PASS,  INCX,  and INCY for a new case. ..
!        .. the value 9999 for INCX or INCY will appear in the ..
!        .. detailed  output, if any, for cases  that do not involve ..
!        .. these parameters ..
!
         PASs = .TRUE.
         INCx = 9999
         INCy = 9999
         IF ( ICAse==3 .OR. ICAse==11 ) THEN
            CALL CHECK0(sfac)
         ELSEIF ( ICAse==7 .OR. ICAse==8 .OR. ICAse==9 .OR. ICAse==10 ) &
     &            THEN
            CALL CHECK1(sfac)
         ELSEIF ( ICAse==1 .OR. ICAse==2 .OR. ICAse==5 .OR.             &
     &            ICAse==6 .OR. ICAse==12 .OR. ICAse==13 ) THEN
            CALL CHECK2(sfac)
         ELSEIF ( ICAse==4 ) THEN
            CALL CHECK3(sfac)
         ENDIF
!        -- Print
         IF ( PASs ) WRITE (NOUT,99002)
99002    FORMAT ('                                    ----- PASS -----')
      ENDDO
      STOP
      END PROGRAM SBLAT1
!*==header.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
      SUBROUTINE HEADER
      IMPLICIT NONE
!*--HEADER102
!     .. Parameters ..
      INTEGER NOUT
      PARAMETER (NOUT=6)
!     .. Scalars in Common ..
      INTEGER ICAse , INCx , INCy , N
      LOGICAL PASs
!     .. Local Arrays ..
      CHARACTER*6 l(13)
!     .. Common blocks ..
      COMMON /COMBLA/ ICAse , N , INCx , INCy , PASs
!     .. Data statements ..
      DATA l(1)/' SDOT '/
      DATA l(2)/'SAXPY '/
      DATA l(3)/'SROTG '/
      DATA l(4)/' SROT '/
      DATA l(5)/'SCOPY '/
      DATA l(6)/'SSWAP '/
      DATA l(7)/'SNRM2 '/
      DATA l(8)/'SASUM '/
      DATA l(9)/'SSCAL '/
      DATA l(10)/'ISAMAX'/
      DATA l(11)/'SROTMG'/
      DATA l(12)/'SROTM '/
      DATA l(13)/'SDSDOT'/
!     .. Executable Statements ..
      WRITE (NOUT,99001) ICAse , l(ICAse)
!
99001 FORMAT (/' Test of subprogram number',I3,12X,A6)
      RETURN
      END SUBROUTINE HEADER
!*==check0.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
      SUBROUTINE CHECK0(Sfac)
      IMPLICIT NONE
!*--CHECK0136
!     .. Parameters ..
      INTEGER NOUT
      PARAMETER (NOUT=6)
!     .. Scalar Arguments ..
      REAL Sfac
!     .. Scalars in Common ..
      INTEGER ICAse , INCx , INCy , N
      LOGICAL PASs
!     .. Local Scalars ..
      REAL d12 , sa , sb , sc , ss
      INTEGER i , k
!     .. Local Arrays ..
      REAL da1(8) , datrue(8) , db1(8) , dbtrue(8) , dc1(8) , ds1(8) ,  &
     &     dab(4,9) , dtemp(9) , dtrue(9,9)
!     .. External Subroutines ..
      EXTERNAL SROTG , SROTMG , STEST , STEST1
!     .. Common blocks ..
      COMMON /COMBLA/ ICAse , N , INCx , INCy , PASs
!     .. Data statements ..
      DATA da1/0.3E0 , 0.4E0 , -0.3E0 , -0.4E0 , -0.3E0 , 0.0E0 ,       &
     &     0.0E0 , 1.0E0/
      DATA db1/0.4E0 , 0.3E0 , 0.4E0 , 0.3E0 , -0.4E0 , 0.0E0 , 1.0E0 , &
     &     0.0E0/
      DATA dc1/0.6E0 , 0.8E0 , -0.6E0 , 0.8E0 , 0.6E0 , 1.0E0 , 0.0E0 , &
     &     1.0E0/
      DATA ds1/0.8E0 , 0.6E0 , 0.8E0 , -0.6E0 , 0.8E0 , 0.0E0 , 1.0E0 , &
     &     0.0E0/
      DATA datrue/0.5E0 , 0.5E0 , 0.5E0 , -0.5E0 , -0.5E0 , 0.0E0 ,     &
     &     1.0E0 , 1.0E0/
      DATA dbtrue/0.0E0 , 0.6E0 , 0.0E0 , -0.6E0 , 0.0E0 , 0.0E0 ,      &
     &     1.0E0 , 0.0E0/
!     INPUT FOR MODIFIED GIVENS
      DATA dab/.1E0 , .3E0 , 1.2E0 , .2E0 , .7E0 , .2E0 , .6E0 , 4.2E0 ,&
     &     0.E0 , 0.E0 , 0.E0 , 0.E0 , 4.E0 , -1.E0 , 2.E0 , 4.E0 ,     &
     &     6.E-10 , 2.E-2 , 1.E5 , 10.E0 , 4.E10 , 2.E-2 , 1.E-5 ,      &
     &     10.E0 , 2.E-10 , 4.E-2 , 1.E5 , 10.E0 , 2.E10 , 4.E-2 ,      &
     &     1.E-5 , 10.E0 , 4.E0 , -2.E0 , 8.E0 , 4.E0/
!    TRUE RESULTS FOR MODIFIED GIVENS
      DATA dtrue/0.E0 , 0.E0 , 1.3E0 , .2E0 , 0.E0 , 0.E0 , 0.E0 ,      &
     &     .5E0 , 0.E0 , 0.E0 , 0.E0 , 4.5E0 , 4.2E0 , 1.E0 , .5E0 ,    &
     &     0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , -2.E0 ,     &
     &     0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 4.E0 ,      &
     &     -1.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 15.E-3 , 0.E0 ,   &
     &     10.E0 , -1.E0 , 0.E0 , -1.E-4 , 0.E0 , 1.E0 , 0.E0 , 0.E0 ,  &
     &     6144.E-5 , 10.E0 , -1.E0 , 4096.E0 , -1.E6 , 0.E0 , 1.E0 ,   &
     &     0.E0 , 0.E0 , 15.E0 , 10.E0 , -1.E0 , 5.E-5 , 0.E0 , 1.E0 ,  &
     &     0.E0 , 0.E0 , 0.E0 , 15.E0 , 10.E0 , -1.E0 , 5.E5 ,          &
     &     -4096.E0 , 1.E0 , 4096.E-6 , 0.E0 , 0.E0 , 7.E0 , 4.E0 ,     &
     &     0.E0 , 0.E0 , -.5E0 , -.25E0 , 0.E0/
!                   4096 = 2 ** 12
      DATA d12/4096.E0/
      dtrue(1,1) = 12.E0/130.E0
      dtrue(2,1) = 36.E0/130.E0
      dtrue(7,1) = -1.E0/6.E0
      dtrue(1,2) = 14.E0/75.E0
      dtrue(2,2) = 49.E0/75.E0
      dtrue(9,2) = 1.E0/7.E0
      dtrue(1,5) = 45.E-11*(d12*d12)
      dtrue(3,5) = 4.E5/(3.E0*d12)
      dtrue(6,5) = 1.E0/d12
      dtrue(8,5) = 1.E4/(3.E0*d12)
      dtrue(1,6) = 4.E10/(1.5E0*d12*d12)
      dtrue(2,6) = 2.E-2/1.5E0
      dtrue(8,6) = 5.E-7*d12
      dtrue(1,7) = 4.E0/150.E0
      dtrue(2,7) = (2.E-10/1.5E0)*(d12*d12)
      dtrue(7,7) = -dtrue(6,5)
      dtrue(9,7) = 1.E4/d12
      dtrue(1,8) = dtrue(1,7)
      dtrue(2,8) = 2.E10/(1.5E0*d12*d12)
      dtrue(1,9) = 32.E0/7.E0
      dtrue(2,9) = -16.E0/7.E0
!     .. Executable Statements ..
!
!     Compute true values which cannot be prestored
!     in decimal notation
!
      dbtrue(1) = 1.0E0/0.6E0
      dbtrue(3) = -1.0E0/0.6E0
      dbtrue(5) = 1.0E0/0.6E0
!
      DO k = 1 , 8
!        .. Set N=K for identification in output if any ..
         N = k
         IF ( ICAse==3 ) THEN
!           .. SROTG ..
            IF ( k>8 ) EXIT
            sa = da1(k)
            sb = db1(k)
            CALL SROTG(sa,sb,sc,ss)
            CALL STEST1(sa,datrue(k),datrue(k),Sfac)
            CALL STEST1(sb,dbtrue(k),dbtrue(k),Sfac)
            CALL STEST1(sc,dc1(k),dc1(k),Sfac)
            CALL STEST1(ss,ds1(k),ds1(k),Sfac)
         ELSEIF ( ICAse==11 ) THEN
!           .. SROTMG ..
            DO i = 1 , 4
               dtemp(i) = dab(i,k)
               dtemp(i+4) = 0.0
            ENDDO
            dtemp(9) = 0.0
            CALL SROTMG(dtemp(1),dtemp(2),dtemp(3),dtemp(4),dtemp(5))
            CALL STEST(9,dtemp,dtrue(1,k),dtrue(1,k),Sfac)
         ELSE
            WRITE (NOUT,*) ' Shouldn''t be here in CHECK0'
            STOP
         ENDIF
      ENDDO
      END SUBROUTINE CHECK0
!*==check1.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
      SUBROUTINE CHECK1(Sfac)
      IMPLICIT NONE
!*--CHECK1249
!     .. Parameters ..
      INTEGER NOUT
      PARAMETER (NOUT=6)
!     .. Scalar Arguments ..
      REAL Sfac
!     .. Scalars in Common ..
      INTEGER ICAse , INCx , INCy , N
      LOGICAL PASs
!     .. Local Scalars ..
      INTEGER i , ix , len , np1
!     .. Local Arrays ..
      REAL dtrue1(5) , dtrue3(5) , dtrue5(8,5,2) , dv(8,5,2) , dvr(8) , &
     &     sa(10) , stemp(1) , strue(8) , sx(8) , sxr(15)
      INTEGER itrue2(5) , itruec(5)
!     .. External Functions ..
      REAL SASUM , SNRM2
      INTEGER ISAMAX
      EXTERNAL SASUM , SNRM2 , ISAMAX
!     .. External Subroutines ..
      EXTERNAL ITEST1 , SSCAL , STEST , STEST1
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     .. Common blocks ..
      COMMON /COMBLA/ ICAse , N , INCx , INCy , PASs
!     .. Data statements ..
      DATA sa/0.3E0 , -1.0E0 , 0.0E0 , 1.0E0 , 0.3E0 , 0.3E0 , 0.3E0 ,  &
     &     0.3E0 , 0.3E0 , 0.3E0/
      DATA dv/0.1E0 , 2.0E0 , 2.0E0 , 2.0E0 , 2.0E0 , 2.0E0 , 2.0E0 ,   &
     &     2.0E0 , 0.3E0 , 3.0E0 , 3.0E0 , 3.0E0 , 3.0E0 , 3.0E0 ,      &
     &     3.0E0 , 3.0E0 , 0.3E0 , -0.4E0 , 4.0E0 , 4.0E0 , 4.0E0 ,     &
     &     4.0E0 , 4.0E0 , 4.0E0 , 0.2E0 , -0.6E0 , 0.3E0 , 5.0E0 ,     &
     &     5.0E0 , 5.0E0 , 5.0E0 , 5.0E0 , 0.1E0 , -0.3E0 , 0.5E0 ,     &
     &     -0.1E0 , 6.0E0 , 6.0E0 , 6.0E0 , 6.0E0 , 0.1E0 , 8.0E0 ,     &
     &     8.0E0 , 8.0E0 , 8.0E0 , 8.0E0 , 8.0E0 , 8.0E0 , 0.3E0 ,      &
     &     9.0E0 , 9.0E0 , 9.0E0 , 9.0E0 , 9.0E0 , 9.0E0 , 9.0E0 ,      &
     &     0.3E0 , 2.0E0 , -0.4E0 , 2.0E0 , 2.0E0 , 2.0E0 , 2.0E0 ,     &
     &     2.0E0 , 0.2E0 , 3.0E0 , -0.6E0 , 5.0E0 , 0.3E0 , 2.0E0 ,     &
     &     2.0E0 , 2.0E0 , 0.1E0 , 4.0E0 , -0.3E0 , 6.0E0 , -0.5E0 ,    &
     &     7.0E0 , -0.1E0 , 3.0E0/
      DATA dvr/8.0E0 , -7.0E0 , 9.0E0 , 5.0E0 , 9.0E0 , 8.0E0 , 7.0E0 , &
     &     7.0E0/
      DATA dtrue1/0.0E0 , 0.3E0 , 0.5E0 , 0.7E0 , 0.6E0/
      DATA dtrue3/0.0E0 , 0.3E0 , 0.7E0 , 1.1E0 , 1.0E0/
      DATA dtrue5/0.10E0 , 2.0E0 , 2.0E0 , 2.0E0 , 2.0E0 , 2.0E0 ,      &
     &     2.0E0 , 2.0E0 , -0.3E0 , 3.0E0 , 3.0E0 , 3.0E0 , 3.0E0 ,     &
     &     3.0E0 , 3.0E0 , 3.0E0 , 0.0E0 , 0.0E0 , 4.0E0 , 4.0E0 ,      &
     &     4.0E0 , 4.0E0 , 4.0E0 , 4.0E0 , 0.20E0 , -0.60E0 , 0.30E0 ,  &
     &     5.0E0 , 5.0E0 , 5.0E0 , 5.0E0 , 5.0E0 , 0.03E0 , -0.09E0 ,   &
     &     0.15E0 , -0.03E0 , 6.0E0 , 6.0E0 , 6.0E0 , 6.0E0 , 0.10E0 ,  &
     &     8.0E0 , 8.0E0 , 8.0E0 , 8.0E0 , 8.0E0 , 8.0E0 , 8.0E0 ,      &
     &     0.09E0 , 9.0E0 , 9.0E0 , 9.0E0 , 9.0E0 , 9.0E0 , 9.0E0 ,     &
     &     9.0E0 , 0.09E0 , 2.0E0 , -0.12E0 , 2.0E0 , 2.0E0 , 2.0E0 ,   &
     &     2.0E0 , 2.0E0 , 0.06E0 , 3.0E0 , -0.18E0 , 5.0E0 , 0.09E0 ,  &
     &     2.0E0 , 2.0E0 , 2.0E0 , 0.03E0 , 4.0E0 , -0.09E0 , 6.0E0 ,   &
     &     -0.15E0 , 7.0E0 , -0.03E0 , 3.0E0/
      DATA itrue2/0 , 1 , 2 , 2 , 3/
      DATA itruec/0 , 1 , 1 , 1 , 1/
!     .. Executable Statements ..
      DO INCx = 1 , 2
         DO np1 = 1 , 5
            N = np1 - 1
            len = 2*MAX(N,1)
!           .. Set vector arguments ..
            DO i = 1 , len
               sx(i) = dv(i,np1,INCx)
            ENDDO
!
            IF ( ICAse==7 ) THEN
!              .. SNRM2 ..
               stemp(1) = dtrue1(np1)
               CALL STEST1(SNRM2(N,sx,INCx),stemp(1),stemp,Sfac)
            ELSEIF ( ICAse==8 ) THEN
!              .. SASUM ..
               stemp(1) = dtrue3(np1)
               CALL STEST1(SASUM(N,sx,INCx),stemp(1),stemp,Sfac)
            ELSEIF ( ICAse==9 ) THEN
!              .. SSCAL ..
               CALL SSCAL(N,sa((INCx-1)*5+np1),sx,INCx)
               DO i = 1 , len
                  strue(i) = dtrue5(i,np1,INCx)
               ENDDO
               CALL STEST(len,sx,strue,strue,Sfac)
            ELSEIF ( ICAse==10 ) THEN
!              .. ISAMAX ..
               CALL ITEST1(ISAMAX(N,sx,INCx),itrue2(np1))
               DO i = 1 , len
                  sx(i) = 42.0E0
               ENDDO
               CALL ITEST1(ISAMAX(N,sx,INCx),itruec(np1))
            ELSE
               WRITE (NOUT,*) ' Shouldn''t be here in CHECK1'
               STOP
            ENDIF
         ENDDO
         IF ( ICAse==10 ) THEN
            N = 8
            ix = 1
            DO i = 1 , N
               sxr(ix) = dvr(i)
               ix = ix + INCx
            ENDDO
            CALL ITEST1(ISAMAX(N,sxr,INCx),3)
         ENDIF
      ENDDO
      END SUBROUTINE CHECK1
!*==check2.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
      SUBROUTINE CHECK2(Sfac)
      IMPLICIT NONE
!*--CHECK2358
!     .. Parameters ..
      INTEGER NOUT
      PARAMETER (NOUT=6)
!     .. Scalar Arguments ..
      REAL Sfac
!     .. Scalars in Common ..
      INTEGER ICAse , INCx , INCy , N
      LOGICAL PASs
!     .. Local Scalars ..
      REAL sa
      INTEGER i , j , ki , kn , kni , kpar , ksize , lenx , leny ,      &
     &        lincx , lincy , mx , my
!     .. Local Arrays ..
      REAL dt10x(7,4,4) , dt10y(7,4,4) , dt7(4,4) , dt8(7,4,4) , dx1(7) &
     &     , dy1(7) , ssize1(4) , ssize2(14,2) , ssize3(4) , ssize(7) , &
     &     stx(7) , sty(7) , sx(7) , sy(7) , dpar(5,4) , dt19x(7,4,16) ,&
     &     dt19xa(7,4,4) , dt19xb(7,4,4) , dt19xc(7,4,4) , dt19xd(7,4,4)&
     &     , dt19y(7,4,16) , dt19ya(7,4,4) , dt19yb(7,4,4) ,            &
     &     dt19yc(7,4,4) , dt19yd(7,4,4) , dtemp(5) , st7b(4,4) ,       &
     &     sty0(1) , sx0(1) , sy0(1)
      INTEGER incxs(4) , incys(4) , lens(4,2) , ns(4)
!     .. External Functions ..
      REAL SDOT , SDSDOT
      EXTERNAL SDOT , SDSDOT
!     .. External Subroutines ..
      EXTERNAL SAXPY , SCOPY , SROTM , SSWAP , STEST , STEST1
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MIN
!     .. Common blocks ..
      COMMON /COMBLA/ ICAse , N , INCx , INCy , PASs
!     .. Data statements ..
      EQUIVALENCE (dt19x(1,1,1),dt19xa(1,1,1))
      EQUIVALENCE (dt19x(1,1,5),dt19xb(1,1,1))
      EQUIVALENCE (dt19x(1,1,9),dt19xc(1,1,1))
      EQUIVALENCE (dt19x(1,1,13),dt19xd(1,1,1))
      EQUIVALENCE (dt19y(1,1,1),dt19ya(1,1,1))
      EQUIVALENCE (dt19y(1,1,5),dt19yb(1,1,1))
      EQUIVALENCE (dt19y(1,1,9),dt19yc(1,1,1))
      EQUIVALENCE (dt19y(1,1,13),dt19yd(1,1,1))
 
      DATA sa/0.3E0/
      DATA incxs/1 , 2 , -2 , -1/
      DATA incys/1 , -2 , 1 , -2/
      DATA lens/1 , 1 , 2 , 4 , 1 , 1 , 3 , 7/
      DATA ns/0 , 1 , 2 , 4/
      DATA dx1/0.6E0 , 0.1E0 , -0.5E0 , 0.8E0 , 0.9E0 , -0.3E0 , -0.4E0/
      DATA dy1/0.5E0 , -0.9E0 , 0.3E0 , 0.7E0 , -0.6E0 , 0.2E0 , 0.8E0/
      DATA dt7/0.0E0 , 0.30E0 , 0.21E0 , 0.62E0 , 0.0E0 , 0.30E0 ,      &
     &     -0.07E0 , 0.85E0 , 0.0E0 , 0.30E0 , -0.79E0 , -0.74E0 ,      &
     &     0.0E0 , 0.30E0 , 0.33E0 , 1.27E0/
      DATA st7b/.1 , .4 , .31 , .72 , .1 , .4 , .03 , .95 , .1 , .4 ,   &
     &     -.69 , -.64 , .1 , .4 , .43 , 1.37/
      DATA dt8/0.5E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,  &
     &     0.68E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,     &
     &     0.68E0 , -0.87E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,   &
     &     0.68E0 , -0.87E0 , 0.15E0 , 0.94E0 , 0.0E0 , 0.0E0 , 0.0E0 , &
     &     0.5E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,      &
     &     0.68E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,     &
     &     0.35E0 , -0.9E0 , 0.48E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,   &
     &     0.38E0 , -0.9E0 , 0.57E0 , 0.7E0 , -0.75E0 , 0.2E0 , 0.98E0 ,&
     &     0.5E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,      &
     &     0.68E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,     &
     &     0.35E0 , -0.72E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,   &
     &     0.38E0 , -0.63E0 , 0.15E0 , 0.88E0 , 0.0E0 , 0.0E0 , 0.0E0 , &
     &     0.5E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,      &
     &     0.68E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,     &
     &     0.68E0 , -0.9E0 , 0.33E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,   &
     &     0.68E0 , -0.9E0 , 0.33E0 , 0.7E0 , -0.75E0 , 0.2E0 , 1.04E0/
      DATA dt10x/0.6E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,&
     &     0.5E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,      &
     &     0.5E0 , -0.9E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,     &
     &     0.5E0 , -0.9E0 , 0.3E0 , 0.7E0 , 0.0E0 , 0.0E0 , 0.0E0 ,     &
     &     0.6E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,      &
     &     0.5E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,      &
     &     0.3E0 , 0.1E0 , 0.5E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,      &
     &     0.8E0 , 0.1E0 , -0.6E0 , 0.8E0 , 0.3E0 , -0.3E0 , 0.5E0 ,    &
     &     0.6E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,      &
     &     0.5E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,      &
     &     -0.9E0 , 0.1E0 , 0.5E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,     &
     &     0.7E0 , 0.1E0 , 0.3E0 , 0.8E0 , -0.9E0 , -0.3E0 , 0.5E0 ,    &
     &     0.6E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,      &
     &     0.5E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,      &
     &     0.5E0 , 0.3E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,      &
     &     0.5E0 , 0.3E0 , -0.6E0 , 0.8E0 , 0.0E0 , 0.0E0 , 0.0E0/
      DATA dt10y/0.5E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,&
     &     0.6E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,      &
     &     0.6E0 , 0.1E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,      &
     &     0.6E0 , 0.1E0 , -0.5E0 , 0.8E0 , 0.0E0 , 0.0E0 , 0.0E0 ,     &
     &     0.5E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,      &
     &     0.6E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,      &
     &     -0.5E0 , -0.9E0 , 0.6E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,    &
     &     -0.4E0 , -0.9E0 , 0.9E0 , 0.7E0 , -0.5E0 , 0.2E0 , 0.6E0 ,   &
     &     0.5E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,      &
     &     0.6E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,      &
     &     -0.5E0 , 0.6E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,     &
     &     -0.4E0 , 0.9E0 , -0.5E0 , 0.6E0 , 0.0E0 , 0.0E0 , 0.0E0 ,    &
     &     0.5E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,      &
     &     0.6E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,      &
     &     0.6E0 , -0.9E0 , 0.1E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,     &
     &     0.6E0 , -0.9E0 , 0.1E0 , 0.7E0 , -0.5E0 , 0.2E0 , 0.8E0/
      DATA ssize1/0.0E0 , 0.3E0 , 1.6E0 , 3.2E0/
      DATA ssize2/0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,       &
     &     0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,      &
     &     0.0E0 , 1.17E0 , 1.17E0 , 1.17E0 , 1.17E0 , 1.17E0 , 1.17E0 ,&
     &     1.17E0 , 1.17E0 , 1.17E0 , 1.17E0 , 1.17E0 , 1.17E0 ,        &
     &     1.17E0 , 1.17E0/
      DATA ssize3/.1 , .4 , 1.7 , 3.3/
!
!                         FOR DROTM
!
      DATA dpar/ - 2.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , -1.E0 , 2.E0 ,    &
     &     -3.E0 , -4.E0 , 5.E0 , 0.E0 , 0.E0 , 2.E0 , -3.E0 , 0.E0 ,   &
     &     1.E0 , 5.E0 , 2.E0 , 0.E0 , -4.E0/
!                        TRUE X RESULTS F0R ROTATIONS DROTM
      DATA dt19xa/.6E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 ,      &
     &     .6E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , .6E0 ,      &
     &     0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , .6E0 , 0.E0 ,      &
     &     0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , .6E0 , 0.E0 , 0.E0 ,      &
     &     0.E0 , 0.E0 , 0.E0 , 0.E0 , -.8E0 , 0.E0 , 0.E0 , 0.E0 ,     &
     &     0.E0 , 0.E0 , 0.E0 , -.9E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 ,     &
     &     0.E0 , 0.E0 , 3.5E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 ,     &
     &     0.E0 , .6E0 , .1E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 ,      &
     &     -.8E0 , 3.8E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , -.9E0 ,   &
     &     2.8E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 3.5E0 , -.4E0 ,   &
     &     0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , .6E0 , .1E0 , -.5E0 ,     &
     &     .8E0 , 0.E0 , 0.E0 , 0.E0 , -.8E0 , 3.8E0 , -2.2E0 , -1.2E0 ,&
     &     0.E0 , 0.E0 , 0.E0 , -.9E0 , 2.8E0 , -1.4E0 , -1.3E0 , 0.E0 ,&
     &     0.E0 , 0.E0 , 3.5E0 , -.4E0 , -2.2E0 , 4.7E0 , 0.E0 , 0.E0 , &
     &     0.E0/
!
      DATA dt19xb/.6E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 ,      &
     &     .6E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , .6E0 ,      &
     &     0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , .6E0 , 0.E0 ,      &
     &     0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , .6E0 , 0.E0 , 0.E0 ,      &
     &     0.E0 , 0.E0 , 0.E0 , 0.E0 , -.8E0 , 0.E0 , 0.E0 , 0.E0 ,     &
     &     0.E0 , 0.E0 , 0.E0 , -.9E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 ,     &
     &     0.E0 , 0.E0 , 3.5E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 ,     &
     &     0.E0 , .6E0 , .1E0 , -.5E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 ,     &
     &     0.E0 , .1E0 , -3.0E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , -.3E0 ,   &
     &     .1E0 , -2.0E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 3.3E0 , .1E0 ,   &
     &     -2.0E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , .6E0 , .1E0 , -.5E0 ,   &
     &     .8E0 , .9E0 , -.3E0 , -.4E0 , -2.0E0 , .1E0 , 1.4E0 , .8E0 , &
     &     .6E0 , -.3E0 , -2.8E0 , -1.8E0 , .1E0 , 1.3E0 , .8E0 , 0.E0 ,&
     &     -.3E0 , -1.9E0 , 3.8E0 , .1E0 , -3.1E0 , .8E0 , 4.8E0 ,      &
     &     -.3E0 , -1.5E0/
!
      DATA dt19xc/.6E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 ,      &
     &     .6E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , .6E0 ,      &
     &     0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , .6E0 , 0.E0 ,      &
     &     0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , .6E0 , 0.E0 , 0.E0 ,      &
     &     0.E0 , 0.E0 , 0.E0 , 0.E0 , -.8E0 , 0.E0 , 0.E0 , 0.E0 ,     &
     &     0.E0 , 0.E0 , 0.E0 , -.9E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 ,     &
     &     0.E0 , 0.E0 , 3.5E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 ,     &
     &     0.E0 , .6E0 , .1E0 , -.5E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 ,     &
     &     4.8E0 , .1E0 , -3.0E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 3.3E0 ,  &
     &     .1E0 , -2.0E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 2.1E0 , .1E0 ,   &
     &     -2.0E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , .6E0 , .1E0 , -.5E0 ,   &
     &     .8E0 , .9E0 , -.3E0 , -.4E0 , -1.6E0 , .1E0 , -2.2E0 , .8E0 ,&
     &     5.4E0 , -.3E0 , -2.8E0 , -1.5E0 , .1E0 , -1.4E0 , .8E0 ,     &
     &     3.6E0 , -.3E0 , -1.9E0 , 3.7E0 , .1E0 , -2.2E0 , .8E0 ,      &
     &     3.6E0 , -.3E0 , -1.5E0/
!
      DATA dt19xd/.6E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 ,      &
     &     .6E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , .6E0 ,      &
     &     0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , .6E0 , 0.E0 ,      &
     &     0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , .6E0 , 0.E0 , 0.E0 ,      &
     &     0.E0 , 0.E0 , 0.E0 , 0.E0 , -.8E0 , 0.E0 , 0.E0 , 0.E0 ,     &
     &     0.E0 , 0.E0 , 0.E0 , -.9E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 ,     &
     &     0.E0 , 0.E0 , 3.5E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 ,     &
     &     0.E0 , .6E0 , .1E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 ,      &
     &     -.8E0 , -1.0E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , -.9E0 ,  &
     &     -.8E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 3.5E0 , .8E0 ,    &
     &     0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , .6E0 , .1E0 , -.5E0 ,     &
     &     .8E0 , 0.E0 , 0.E0 , 0.E0 , -.8E0 , -1.0E0 , 1.4E0 , -1.6E0 ,&
     &     0.E0 , 0.E0 , 0.E0 , -.9E0 , -.8E0 , 1.3E0 , -1.6E0 , 0.E0 , &
     &     0.E0 , 0.E0 , 3.5E0 , .8E0 , -3.1E0 , 4.8E0 , 0.E0 , 0.E0 ,  &
     &     0.E0/
!                        TRUE Y RESULTS FOR ROTATIONS DROTM
      DATA dt19ya/.5E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 ,      &
     &     .5E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , .5E0 ,      &
     &     0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , .5E0 , 0.E0 ,      &
     &     0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , .5E0 , 0.E0 , 0.E0 ,      &
     &     0.E0 , 0.E0 , 0.E0 , 0.E0 , .7E0 , 0.E0 , 0.E0 , 0.E0 ,      &
     &     0.E0 , 0.E0 , 0.E0 , 1.7E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 ,     &
     &     0.E0 , 0.E0 , -2.6E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 ,    &
     &     0.E0 , .5E0 , -.9E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 ,     &
     &     .7E0 , -4.8E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 1.7E0 ,   &
     &     -.7E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , -2.6E0 , 3.5E0 ,  &
     &     0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , .5E0 , -.9E0 , .3E0 ,     &
     &     .7E0 , 0.E0 , 0.E0 , 0.E0 , .7E0 , -4.8E0 , 3.0E0 , 1.1E0 ,  &
     &     0.E0 , 0.E0 , 0.E0 , 1.7E0 , -.7E0 , -.7E0 , 2.3E0 , 0.E0 ,  &
     &     0.E0 , 0.E0 , -2.6E0 , 3.5E0 , -.7E0 , -3.6E0 , 0.E0 , 0.E0 ,&
     &     0.E0/
!
      DATA dt19yb/.5E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 ,      &
     &     .5E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , .5E0 ,      &
     &     0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , .5E0 , 0.E0 ,      &
     &     0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , .5E0 , 0.E0 , 0.E0 ,      &
     &     0.E0 , 0.E0 , 0.E0 , 0.E0 , .7E0 , 0.E0 , 0.E0 , 0.E0 ,      &
     &     0.E0 , 0.E0 , 0.E0 , 1.7E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 ,     &
     &     0.E0 , 0.E0 , -2.6E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 ,    &
     &     0.E0 , .5E0 , -.9E0 , .3E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 ,     &
     &     4.0E0 , -.9E0 , -.3E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , -.5E0 ,  &
     &     -.9E0 , 1.5E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , -1.5E0 , -.9E0 , &
     &     -1.8E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , .5E0 , -.9E0 , .3E0 ,   &
     &     .7E0 , -.6E0 , .2E0 , .8E0 , 3.7E0 , -.9E0 , -1.2E0 , .7E0 , &
     &     -1.5E0 , .2E0 , 2.2E0 , -.3E0 , -.9E0 , 2.1E0 , .7E0 ,       &
     &     -1.6E0 , .2E0 , 2.0E0 , -1.6E0 , -.9E0 , -2.1E0 , .7E0 ,     &
     &     2.9E0 , .2E0 , -3.8E0/
!
      DATA dt19yc/.5E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 ,      &
     &     .5E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , .5E0 ,      &
     &     0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , .5E0 , 0.E0 ,      &
     &     0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , .5E0 , 0.E0 , 0.E0 ,      &
     &     0.E0 , 0.E0 , 0.E0 , 0.E0 , .7E0 , 0.E0 , 0.E0 , 0.E0 ,      &
     &     0.E0 , 0.E0 , 0.E0 , 1.7E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 ,     &
     &     0.E0 , 0.E0 , -2.6E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 ,    &
     &     0.E0 , .5E0 , -.9E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 ,     &
     &     4.0E0 , -6.3E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , -.5E0 ,  &
     &     .3E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , -1.5E0 , 3.0E0 ,   &
     &     0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , .5E0 , -.9E0 , .3E0 ,     &
     &     .7E0 , 0.E0 , 0.E0 , 0.E0 , 3.7E0 , -7.2E0 , 3.0E0 , 1.7E0 , &
     &     0.E0 , 0.E0 , 0.E0 , -.3E0 , .9E0 , -.7E0 , 1.9E0 , 0.E0 ,   &
     &     0.E0 , 0.E0 , -1.6E0 , 2.7E0 , -.7E0 , -3.4E0 , 0.E0 , 0.E0 ,&
     &     0.E0/
!
      DATA dt19yd/.5E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 ,      &
     &     .5E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , .5E0 ,      &
     &     0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , .5E0 , 0.E0 ,      &
     &     0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , .5E0 , 0.E0 , 0.E0 ,      &
     &     0.E0 , 0.E0 , 0.E0 , 0.E0 , .7E0 , 0.E0 , 0.E0 , 0.E0 ,      &
     &     0.E0 , 0.E0 , 0.E0 , 1.7E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 ,     &
     &     0.E0 , 0.E0 , -2.6E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 ,    &
     &     0.E0 , .5E0 , -.9E0 , .3E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 ,     &
     &     .7E0 , -.9E0 , 1.2E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , 1.7E0 ,   &
     &     -.9E0 , .5E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , -2.6E0 , -.9E0 ,  &
     &     -1.3E0 , 0.E0 , 0.E0 , 0.E0 , 0.E0 , .5E0 , -.9E0 , .3E0 ,   &
     &     .7E0 , -.6E0 , .2E0 , .8E0 , .7E0 , -.9E0 , 1.2E0 , .7E0 ,   &
     &     -1.5E0 , .2E0 , 1.6E0 , 1.7E0 , -.9E0 , .5E0 , .7E0 ,        &
     &     -1.6E0 , .2E0 , 2.4E0 , -2.6E0 , -.9E0 , -1.3E0 , .7E0 ,     &
     &     2.9E0 , .2E0 , -4.0E0/
!
!     .. Executable Statements ..
!
      DO ki = 1 , 4
         INCx = incxs(ki)
         INCy = incys(ki)
         mx = ABS(INCx)
         my = ABS(INCy)
!
         DO kn = 1 , 4
            N = ns(kn)
            ksize = MIN(2,kn)
            lenx = lens(kn,mx)
            leny = lens(kn,my)
!           .. Initialize all argument arrays ..
            DO i = 1 , 7
               sx(i) = dx1(i)
               sy(i) = dy1(i)
            ENDDO
!
            IF ( ICAse==1 ) THEN
!              .. SDOT ..
               CALL STEST1(SDOT(N,sx,INCx,sy,INCy),dt7(kn,ki),ssize1(kn)&
     &                     ,Sfac)
            ELSEIF ( ICAse==2 ) THEN
!              .. SAXPY ..
               CALL SAXPY(N,sa,sx,INCx,sy,INCy)
               DO j = 1 , leny
                  sty(j) = dt8(j,kn,ki)
               ENDDO
               CALL STEST(leny,sy,sty,ssize2(1,ksize),Sfac)
            ELSEIF ( ICAse==5 ) THEN
!              .. SCOPY ..
               DO i = 1 , 7
                  sty(i) = dt10y(i,kn,ki)
               ENDDO
               CALL SCOPY(N,sx,INCx,sy,INCy)
               CALL STEST(leny,sy,sty,ssize2(1,1),1.0E0)
               IF ( ki==1 ) THEN
                  sx0(1) = 42.0E0
                  sy0(1) = 43.0E0
                  IF ( N==0 ) THEN
                     sty0(1) = sy0(1)
                  ELSE
                     sty0(1) = sx0(1)
                  ENDIF
                  lincx = INCx
                  INCx = 0
                  lincy = INCy
                  INCy = 0
                  CALL SCOPY(N,sx0,INCx,sy0,INCy)
                  CALL STEST(1,sy0,sty0,ssize2(1,1),1.0E0)
                  INCx = lincx
                  INCy = lincy
               ENDIF
            ELSEIF ( ICAse==6 ) THEN
!              .. SSWAP ..
               CALL SSWAP(N,sx,INCx,sy,INCy)
               DO i = 1 , 7
                  stx(i) = dt10x(i,kn,ki)
                  sty(i) = dt10y(i,kn,ki)
               ENDDO
               CALL STEST(lenx,sx,stx,ssize2(1,1),1.0E0)
               CALL STEST(leny,sy,sty,ssize2(1,1),1.0E0)
            ELSEIF ( ICAse==12 ) THEN
!              .. SROTM ..
               kni = kn + 4*(ki-1)
               DO kpar = 1 , 4
                  DO i = 1 , 7
                     sx(i) = dx1(i)
                     sy(i) = dy1(i)
                     stx(i) = dt19x(i,kpar,kni)
                     sty(i) = dt19y(i,kpar,kni)
                  ENDDO
!
                  DO i = 1 , 5
                     dtemp(i) = dpar(i,kpar)
                  ENDDO
!
                  DO i = 1 , lenx
                     ssize(i) = stx(i)
                  ENDDO
!                   SEE REMARK ABOVE ABOUT DT11X(1,2,7)
!                       AND DT11X(5,3,8).
                  IF ( (kpar==2) .AND. (kni==7) ) ssize(1) = 2.4E0
                  IF ( (kpar==3) .AND. (kni==8) ) ssize(5) = 1.8E0
!
                  CALL SROTM(N,sx,INCx,sy,INCy,dtemp)
                  CALL STEST(lenx,sx,stx,ssize,Sfac)
                  CALL STEST(leny,sy,sty,sty,Sfac)
               ENDDO
            ELSEIF ( ICAse==13 ) THEN
!              .. SDSROT ..
               CALL STEST1(SDSDOT(N,.1,sx,INCx,sy,INCy),st7b(kn,ki),    &
     &                     ssize3(kn),Sfac)
            ELSE
               WRITE (NOUT,*) ' Shouldn''t be here in CHECK2'
               STOP
            ENDIF
         ENDDO
      ENDDO
      END SUBROUTINE CHECK2
!*==check3.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
      SUBROUTINE CHECK3(Sfac)
      IMPLICIT NONE
!*--CHECK3705
!     .. Parameters ..
      INTEGER NOUT
      PARAMETER (NOUT=6)
!     .. Scalar Arguments ..
      REAL Sfac
!     .. Scalars in Common ..
      INTEGER ICAse , INCx , INCy , N
      LOGICAL PASs
!     .. Local Scalars ..
      REAL sc , ss
      INTEGER i , k , ki , kn , ksize , lenx , leny , mx , my
!     .. Local Arrays ..
      REAL copyx(5) , copyy(5) , dt9x(7,4,4) , dt9y(7,4,4) , dx1(7) ,   &
     &     dy1(7) , mwpc(11) , mwps(11) , mwpstx(5) , mwpsty(5) ,       &
     &     mwptx(11,5) , mwpty(11,5) , mwpx(5) , mwpy(5) , ssize2(14,2) &
     &     , stx(7) , sty(7) , sx(7) , sy(7)
      INTEGER incxs(4) , incys(4) , lens(4,2) , mwpinx(11) , mwpiny(11) &
     &        , mwpn(11) , ns(4)
!     .. External Subroutines ..
      EXTERNAL SROT , STEST
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MIN
!     .. Common blocks ..
      COMMON /COMBLA/ ICAse , N , INCx , INCy , PASs
!     .. Data statements ..
      DATA incxs/1 , 2 , -2 , -1/
      DATA incys/1 , -2 , 1 , -2/
      DATA lens/1 , 1 , 2 , 4 , 1 , 1 , 3 , 7/
      DATA ns/0 , 1 , 2 , 4/
      DATA dx1/0.6E0 , 0.1E0 , -0.5E0 , 0.8E0 , 0.9E0 , -0.3E0 , -0.4E0/
      DATA dy1/0.5E0 , -0.9E0 , 0.3E0 , 0.7E0 , -0.6E0 , 0.2E0 , 0.8E0/
      DATA sc , ss/0.8E0 , 0.6E0/
      DATA dt9x/0.6E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , &
     &     0.78E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,     &
     &     0.78E0 , -0.46E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,   &
     &     0.78E0 , -0.46E0 , -0.22E0 , 1.06E0 , 0.0E0 , 0.0E0 , 0.0E0 ,&
     &     0.6E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,      &
     &     0.78E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,     &
     &     0.66E0 , 0.1E0 , -0.1E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,    &
     &     0.96E0 , 0.1E0 , -0.76E0 , 0.8E0 , 0.90E0 , -0.3E0 ,         &
     &     -0.02E0 , 0.6E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,    &
     &     0.0E0 , 0.78E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,     &
     &     0.0E0 , -0.06E0 , 0.1E0 , -0.1E0 , 0.0E0 , 0.0E0 , 0.0E0 ,   &
     &     0.0E0 , 0.90E0 , 0.1E0 , -0.22E0 , 0.8E0 , 0.18E0 , -0.3E0 , &
     &     -0.02E0 , 0.6E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,    &
     &     0.0E0 , 0.78E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,     &
     &     0.0E0 , 0.78E0 , 0.26E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,    &
     &     0.0E0 , 0.78E0 , 0.26E0 , -0.76E0 , 1.12E0 , 0.0E0 , 0.0E0 , &
     &     0.0E0/
      DATA dt9y/0.5E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , &
     &     0.04E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,     &
     &     0.04E0 , -0.78E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,   &
     &     0.04E0 , -0.78E0 , 0.54E0 , 0.08E0 , 0.0E0 , 0.0E0 , 0.0E0 , &
     &     0.5E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,      &
     &     0.04E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,     &
     &     0.7E0 , -0.9E0 , -0.12E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,   &
     &     0.64E0 , -0.9E0 , -0.30E0 , 0.7E0 , -0.18E0 , 0.2E0 ,        &
     &     0.28E0 , 0.5E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,     &
     &     0.0E0 , 0.04E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,     &
     &     0.0E0 , 0.7E0 , -1.08E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,    &
     &     0.0E0 , 0.64E0 , -1.26E0 , 0.54E0 , 0.20E0 , 0.0E0 , 0.0E0 , &
     &     0.0E0 , 0.5E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,      &
     &     0.0E0 , 0.04E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,     &
     &     0.0E0 , 0.04E0 , -0.9E0 , 0.18E0 , 0.0E0 , 0.0E0 , 0.0E0 ,   &
     &     0.0E0 , 0.04E0 , -0.9E0 , 0.18E0 , 0.7E0 , -0.18E0 , 0.2E0 , &
     &     0.16E0/
      DATA ssize2/0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,       &
     &     0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 , 0.0E0 ,      &
     &     0.0E0 , 1.17E0 , 1.17E0 , 1.17E0 , 1.17E0 , 1.17E0 , 1.17E0 ,&
     &     1.17E0 , 1.17E0 , 1.17E0 , 1.17E0 , 1.17E0 , 1.17E0 ,        &
     &     1.17E0 , 1.17E0/
!     .. Executable Statements ..
!
      DO ki = 1 , 4
         INCx = incxs(ki)
         INCy = incys(ki)
         mx = ABS(INCx)
         my = ABS(INCy)
!
         DO kn = 1 , 4
            N = ns(kn)
            ksize = MIN(2,kn)
            lenx = lens(kn,mx)
            leny = lens(kn,my)
!
            IF ( ICAse==4 ) THEN
!              .. SROT ..
               DO i = 1 , 7
                  sx(i) = dx1(i)
                  sy(i) = dy1(i)
                  stx(i) = dt9x(i,kn,ki)
                  sty(i) = dt9y(i,kn,ki)
               ENDDO
               CALL SROT(N,sx,INCx,sy,INCy,sc,ss)
               CALL STEST(lenx,sx,stx,ssize2(1,ksize),Sfac)
               CALL STEST(leny,sy,sty,ssize2(1,ksize),Sfac)
            ELSE
               WRITE (NOUT,*) ' Shouldn''t be here in CHECK3'
               STOP
            ENDIF
         ENDDO
      ENDDO
!
      mwpc(1) = 1
      DO i = 2 , 11
         mwpc(i) = 0
      ENDDO
      mwps(1) = 0
      DO i = 2 , 6
         mwps(i) = 1
      ENDDO
      DO i = 7 , 11
         mwps(i) = -1
      ENDDO
      mwpinx(1) = 1
      mwpinx(2) = 1
      mwpinx(3) = 1
      mwpinx(4) = -1
      mwpinx(5) = 1
      mwpinx(6) = -1
      mwpinx(7) = 1
      mwpinx(8) = 1
      mwpinx(9) = -1
      mwpinx(10) = 1
      mwpinx(11) = -1
      mwpiny(1) = 1
      mwpiny(2) = 1
      mwpiny(3) = -1
      mwpiny(4) = -1
      mwpiny(5) = 2
      mwpiny(6) = 1
      mwpiny(7) = 1
      mwpiny(8) = -1
      mwpiny(9) = -1
      mwpiny(10) = 2
      mwpiny(11) = 1
      DO i = 1 , 11
         mwpn(i) = 5
      ENDDO
      mwpn(5) = 3
      mwpn(10) = 3
      DO i = 1 , 5
         mwpx(i) = i
         mwpy(i) = i
         mwptx(1,i) = i
         mwpty(1,i) = i
         mwptx(2,i) = i
         mwpty(2,i) = -i
         mwptx(3,i) = 6 - i
         mwpty(3,i) = i - 6
         mwptx(4,i) = i
         mwpty(4,i) = -i
         mwptx(6,i) = 6 - i
         mwpty(6,i) = i - 6
         mwptx(7,i) = -i
         mwpty(7,i) = i
         mwptx(8,i) = i - 6
         mwpty(8,i) = 6 - i
         mwptx(9,i) = -i
         mwpty(9,i) = i
         mwptx(11,i) = i - 6
         mwpty(11,i) = 6 - i
      ENDDO
      mwptx(5,1) = 1
      mwptx(5,2) = 3
      mwptx(5,3) = 5
      mwptx(5,4) = 4
      mwptx(5,5) = 5
      mwpty(5,1) = -1
      mwpty(5,2) = 2
      mwpty(5,3) = -2
      mwpty(5,4) = 4
      mwpty(5,5) = -3
      mwptx(10,1) = -1
      mwptx(10,2) = -3
      mwptx(10,3) = -5
      mwptx(10,4) = 4
      mwptx(10,5) = 5
      mwpty(10,1) = 1
      mwpty(10,2) = 2
      mwpty(10,3) = 2
      mwpty(10,4) = 4
      mwpty(10,5) = 3
      DO i = 1 , 11
         INCx = mwpinx(i)
         INCy = mwpiny(i)
         DO k = 1 , 5
            copyx(k) = mwpx(k)
            copyy(k) = mwpy(k)
            mwpstx(k) = mwptx(i,k)
            mwpsty(k) = mwpty(i,k)
         ENDDO
         CALL SROT(mwpn(i),copyx,INCx,copyy,INCy,mwpc(i),mwps(i))
         CALL STEST(5,copyx,mwpstx,mwpstx,Sfac)
         CALL STEST(5,copyy,mwpsty,mwpsty,Sfac)
      ENDDO
      END SUBROUTINE CHECK3
!*==stest.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
      SUBROUTINE STEST(Len,Scomp,Strue,Ssize,Sfac)
      IMPLICIT NONE
!*--STEST906
!     ********************************* STEST **************************
!
!     THIS SUBR COMPARES ARRAYS  SCOMP() AND STRUE() OF LENGTH LEN TO
!     SEE IF THE TERM BY TERM DIFFERENCES, MULTIPLIED BY SFAC, ARE
!     NEGLIGIBLE.
!
!     C. L. LAWSON, JPL, 1974 DEC 10
!
!     .. Parameters ..
      INTEGER NOUT
      REAL ZERO
      PARAMETER (NOUT=6,ZERO=0.0E0)
!     .. Scalar Arguments ..
      REAL Sfac
      INTEGER Len
!     .. Array Arguments ..
      REAL Scomp(Len) , Ssize(Len) , Strue(Len)
!     .. Scalars in Common ..
      INTEGER ICAse , INCx , INCy , N
      LOGICAL PASs
!     .. Local Scalars ..
      REAL sd
      INTEGER i
!     .. External Functions ..
      REAL SDIFF
      EXTERNAL SDIFF
!     .. Intrinsic Functions ..
      INTRINSIC ABS
!     .. Common blocks ..
      COMMON /COMBLA/ ICAse , N , INCx , INCy , PASs
!     .. Executable Statements ..
!
      DO i = 1 , Len
         sd = Scomp(i) - Strue(i)
         IF ( ABS(Sfac*sd)>ABS(Ssize(i))*EPSILON(ZERO) ) THEN
!
!                             HERE    SCOMP(I) IS NOT CLOSE TO STRUE(I).
!
            IF ( PASs ) THEN
!                             PRINT FAIL MESSAGE AND HEADER.
               PASs = .FALSE.
               WRITE (NOUT,99001)
!
99001          FORMAT ('                                       FAIL')
               WRITE (NOUT,99002)
99002          FORMAT (/                                                &
     &               ' CASE  N INCX INCY  I                            '&
     &               ,                                                  &
     &        ' COMP(I)                             TRUE(I)  DIFFERENCE'&
     &        ,'     SIZE(I)',/1X)
            ENDIF
            WRITE (NOUT,99003) ICAse , N , INCx , INCy , i , Scomp(i) , &
     &                         Strue(i) , sd , Ssize(i)
99003       FORMAT (1X,I4,I3,2I5,I3,2E36.8,2E12.4)
         ENDIF
      ENDDO
      RETURN
      END SUBROUTINE STEST
!*==stest1.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
      SUBROUTINE STEST1(Scomp1,Strue1,Ssize,Sfac)
      IMPLICIT NONE
!*--STEST1968
!     ************************* STEST1 *****************************
!
!     THIS IS AN INTERFACE SUBROUTINE TO ACCOMMODATE THE FORTRAN
!     REQUIREMENT THAT WHEN A DUMMY ARGUMENT IS AN ARRAY, THE
!     ACTUAL ARGUMENT MUST ALSO BE AN ARRAY OR AN ARRAY ELEMENT.
!
!     C.L. LAWSON, JPL, 1978 DEC 6
!
!     .. Scalar Arguments ..
      REAL Scomp1 , Sfac , Strue1
!     .. Array Arguments ..
      REAL Ssize(*)
!     .. Local Arrays ..
      REAL scomp(1) , strue(1)
!     .. External Subroutines ..
      EXTERNAL STEST
!     .. Executable Statements ..
!
      scomp(1) = Scomp1
      strue(1) = Strue1
      CALL STEST(1,scomp,strue,Ssize,Sfac)
!
      END SUBROUTINE STEST1
!*==sdiff.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
      REAL FUNCTION SDIFF(Sa,Sb)
      IMPLICIT NONE
!*--SDIFF995
!     ********************************* SDIFF **************************
!     COMPUTES DIFFERENCE OF TWO NUMBERS.  C. L. LAWSON, JPL 1974 FEB 15
!
!     .. Scalar Arguments ..
      REAL Sa , Sb
!     .. Executable Statements ..
      SDIFF = Sa - Sb
      END FUNCTION SDIFF
!*==itest1.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
      SUBROUTINE ITEST1(Icomp,Itrue)
      IMPLICIT NONE
!*--ITEST11007
!     ********************************* ITEST1 *************************
!
!     THIS SUBROUTINE COMPARES THE VARIABLES ICOMP AND ITRUE FOR
!     EQUALITY.
!     C. L. LAWSON, JPL, 1974 DEC 10
!
!     .. Parameters ..
      INTEGER NOUT
      PARAMETER (NOUT=6)
!     .. Scalar Arguments ..
      INTEGER Icomp , Itrue
!     .. Scalars in Common ..
      INTEGER ICAse , INCx , INCy , N
      LOGICAL PASs
!     .. Local Scalars ..
      INTEGER id
!     .. Common blocks ..
      COMMON /COMBLA/ ICAse , N , INCx , INCy , PASs
!     .. Executable Statements ..
!
      IF ( Icomp/=Itrue ) THEN
!
!                            HERE ICOMP IS NOT EQUAL TO ITRUE.
!
         IF ( PASs ) THEN
!                             PRINT FAIL MESSAGE AND HEADER.
            PASs = .FALSE.
            WRITE (NOUT,99001)
!
99001       FORMAT ('                                       FAIL')
            WRITE (NOUT,99002)
99002       FORMAT (/' CASE  N INCX INCY                               '&
     &              ,                                                   &
     &        ' COMP                                TRUE     DIFFERENCE'&
     &        ,/1X)
         ENDIF
         id = Icomp - Itrue
         WRITE (NOUT,99003) ICAse , N , INCx , INCy , Icomp , Itrue , id
99003    FORMAT (1X,I4,I3,2I5,2I36,I12)
      ENDIF
      RETURN
      END SUBROUTINE ITEST1
