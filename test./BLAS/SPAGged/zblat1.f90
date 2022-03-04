!*==zblat1.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
 
!> \brief \b ZBLAT1
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       PROGRAM ZBLAT1
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    Test program for the COMPLEX*16 Level 1 BLAS.
!>
!>    Based upon the original BLAS test routine together with:
!>    F06GAF Example Program Text
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
!> \ingroup complex16_blas_testing
!
!  =====================================================================
      PROGRAM ZBLAT1
      IMPLICIT NONE
!*--ZBLAT142
!
!  -- Reference BLAS test routine (version 3.7.0) --
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
      INTEGER ICAse , INCx , INCy , MODe , N
      LOGICAL PASs
!     .. Local Scalars ..
      DOUBLE PRECISION sfac
      INTEGER ic
!     .. External Subroutines ..
      EXTERNAL CHECK1 , CHECK2 , HEADER
!     .. Common blocks ..
      COMMON /COMBLA/ ICAse , N , INCx , INCy , MODe , PASs
!     .. Data statements ..
      DATA sfac/9.765625D-4/
!     .. Executable Statements ..
      WRITE (NOUT,99001)
!
99001 FORMAT (' Complex BLAS Test Program Results',/1X)
      DO ic = 1 , 10
         ICAse = ic
         CALL HEADER
!
!        Initialize PASS, INCX, INCY, and MODE for a new case.
!        The value 9999 for INCX, INCY or MODE will appear in the
!        detailed  output, if any, for cases that do not involve
!        these parameters.
!
         PASs = .TRUE.
         INCx = 9999
         INCy = 9999
         MODe = 9999
         IF ( ICAse<=5 ) THEN
            CALL CHECK2(sfac)
         ELSEIF ( ICAse>=6 ) THEN
            CALL CHECK1(sfac)
         ENDIF
!        -- Print
         IF ( PASs ) WRITE (NOUT,99002)
99002    FORMAT ('                                    ----- PASS -----')
      ENDDO
      STOP
      END PROGRAM ZBLAT1
!*==header.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
      SUBROUTINE HEADER
      IMPLICIT NONE
!*--HEADER97
!     .. Parameters ..
      INTEGER NOUT
      PARAMETER (NOUT=6)
!     .. Scalars in Common ..
      INTEGER ICAse , INCx , INCy , MODe , N
      LOGICAL PASs
!     .. Local Arrays ..
      CHARACTER*6 l(10)
!     .. Common blocks ..
      COMMON /COMBLA/ ICAse , N , INCx , INCy , MODe , PASs
!     .. Data statements ..
      DATA l(1)/'ZDOTC '/
      DATA l(2)/'ZDOTU '/
      DATA l(3)/'ZAXPY '/
      DATA l(4)/'ZCOPY '/
      DATA l(5)/'ZSWAP '/
      DATA l(6)/'DZNRM2'/
      DATA l(7)/'DZASUM'/
      DATA l(8)/'ZSCAL '/
      DATA l(9)/'ZDSCAL'/
      DATA l(10)/'IZAMAX'/
!     .. Executable Statements ..
      WRITE (NOUT,99001) ICAse , l(ICAse)
!
99001 FORMAT (/' Test of subprogram number',I3,12X,A6)
      RETURN
      END SUBROUTINE HEADER
!*==check1.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
      SUBROUTINE CHECK1(Sfac)
      IMPLICIT NONE
!*--CHECK1128
!     .. Parameters ..
      INTEGER NOUT
      PARAMETER (NOUT=6)
!     .. Scalar Arguments ..
      DOUBLE PRECISION Sfac
!     .. Scalars in Common ..
      INTEGER ICAse , INCx , INCy , MODe , N
      LOGICAL PASs
!     .. Local Scalars ..
      COMPLEX*16 ca
      DOUBLE PRECISION sa
      INTEGER i , ix , j , len , np1
!     .. Local Arrays ..
      COMPLEX*16 ctrue5(8,5,2) , ctrue6(8,5,2) , cv(8,5,2) , cvr(8) ,   &
     &           cx(8) , cxr(15) , mwpcs(5) , mwpct(5)
      DOUBLE PRECISION strue2(5) , strue4(5)
      INTEGER itrue3(5) , itruec(5)
!     .. External Functions ..
      DOUBLE PRECISION DZASUM , DZNRM2
      INTEGER IZAMAX
      EXTERNAL DZASUM , DZNRM2 , IZAMAX
!     .. External Subroutines ..
      EXTERNAL ZSCAL , ZDSCAL , CTEST , ITEST1 , STEST1
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     .. Common blocks ..
      COMMON /COMBLA/ ICAse , N , INCx , INCy , MODe , PASs
!     .. Data statements ..
      DATA sa , ca/0.3D0 , (0.4D0,-0.7D0)/
      DATA ((cv(i,j,1),i=1,8),j=1,5)/(0.1D0,0.1D0) , (1.0D0,2.0D0) ,    &
     &      (1.0D0,2.0D0) , (1.0D0,2.0D0) , (1.0D0,2.0D0) ,             &
     &      (1.0D0,2.0D0) , (1.0D0,2.0D0) , (1.0D0,2.0D0) ,             &
     &      (0.3D0,-0.4D0) , (3.0D0,4.0D0) , (3.0D0,4.0D0) ,            &
     &      (3.0D0,4.0D0) , (3.0D0,4.0D0) , (3.0D0,4.0D0) ,             &
     &      (3.0D0,4.0D0) , (3.0D0,4.0D0) , (0.1D0,-0.3D0) ,            &
     &      (0.5D0,-0.1D0) , (5.0D0,6.0D0) , (5.0D0,6.0D0) ,            &
     &      (5.0D0,6.0D0) , (5.0D0,6.0D0) , (5.0D0,6.0D0) ,             &
     &      (5.0D0,6.0D0) , (0.1D0,0.1D0) , (-0.6D0,0.1D0) ,            &
     &      (0.1D0,-0.3D0) , (7.0D0,8.0D0) , (7.0D0,8.0D0) ,            &
     &      (7.0D0,8.0D0) , (7.0D0,8.0D0) , (7.0D0,8.0D0) ,             &
     &      (0.3D0,0.1D0) , (0.5D0,0.0D0) , (0.0D0,0.5D0) ,             &
     &      (0.0D0,0.2D0) , (2.0D0,3.0D0) , (2.0D0,3.0D0) ,             &
     &      (2.0D0,3.0D0) , (2.0D0,3.0D0)/
      DATA ((cv(i,j,2),i=1,8),j=1,5)/(0.1D0,0.1D0) , (4.0D0,5.0D0) ,    &
     &      (4.0D0,5.0D0) , (4.0D0,5.0D0) , (4.0D0,5.0D0) ,             &
     &      (4.0D0,5.0D0) , (4.0D0,5.0D0) , (4.0D0,5.0D0) ,             &
     &      (0.3D0,-0.4D0) , (6.0D0,7.0D0) , (6.0D0,7.0D0) ,            &
     &      (6.0D0,7.0D0) , (6.0D0,7.0D0) , (6.0D0,7.0D0) ,             &
     &      (6.0D0,7.0D0) , (6.0D0,7.0D0) , (0.1D0,-0.3D0) ,            &
     &      (8.0D0,9.0D0) , (0.5D0,-0.1D0) , (2.0D0,5.0D0) ,            &
     &      (2.0D0,5.0D0) , (2.0D0,5.0D0) , (2.0D0,5.0D0) ,             &
     &      (2.0D0,5.0D0) , (0.1D0,0.1D0) , (3.0D0,6.0D0) ,             &
     &      (-0.6D0,0.1D0) , (4.0D0,7.0D0) , (0.1D0,-0.3D0) ,           &
     &      (7.0D0,2.0D0) , (7.0D0,2.0D0) , (7.0D0,2.0D0) ,             &
     &      (0.3D0,0.1D0) , (5.0D0,8.0D0) , (0.5D0,0.0D0) ,             &
     &      (6.0D0,9.0D0) , (0.0D0,0.5D0) , (8.0D0,3.0D0) ,             &
     &      (0.0D0,0.2D0) , (9.0D0,4.0D0)/
      DATA cvr/(8.0D0,8.0D0) , (-7.0D0,-7.0D0) , (9.0D0,9.0D0) ,        &
     &     (5.0D0,5.0D0) , (9.0D0,9.0D0) , (8.0D0,8.0D0) , (7.0D0,7.0D0)&
     &     , (7.0D0,7.0D0)/
      DATA strue2/0.0D0 , 0.5D0 , 0.6D0 , 0.7D0 , 0.8D0/
      DATA strue4/0.0D0 , 0.7D0 , 1.0D0 , 1.3D0 , 1.6D0/
      DATA ((ctrue5(i,j,1),i=1,8),j=1,5)/(0.1D0,0.1D0) , (1.0D0,2.0D0) ,&
     &      (1.0D0,2.0D0) , (1.0D0,2.0D0) , (1.0D0,2.0D0) ,             &
     &      (1.0D0,2.0D0) , (1.0D0,2.0D0) , (1.0D0,2.0D0) ,             &
     &      (-0.16D0,-0.37D0) , (3.0D0,4.0D0) , (3.0D0,4.0D0) ,         &
     &      (3.0D0,4.0D0) , (3.0D0,4.0D0) , (3.0D0,4.0D0) ,             &
     &      (3.0D0,4.0D0) , (3.0D0,4.0D0) , (-0.17D0,-0.19D0) ,         &
     &      (0.13D0,-0.39D0) , (5.0D0,6.0D0) , (5.0D0,6.0D0) ,          &
     &      (5.0D0,6.0D0) , (5.0D0,6.0D0) , (5.0D0,6.0D0) ,             &
     &      (5.0D0,6.0D0) , (0.11D0,-0.03D0) , (-0.17D0,0.46D0) ,       &
     &      (-0.17D0,-0.19D0) , (7.0D0,8.0D0) , (7.0D0,8.0D0) ,         &
     &      (7.0D0,8.0D0) , (7.0D0,8.0D0) , (7.0D0,8.0D0) ,             &
     &      (0.19D0,-0.17D0) , (0.20D0,-0.35D0) , (0.35D0,0.20D0) ,     &
     &      (0.14D0,0.08D0) , (2.0D0,3.0D0) , (2.0D0,3.0D0) ,           &
     &      (2.0D0,3.0D0) , (2.0D0,3.0D0)/
      DATA ((ctrue5(i,j,2),i=1,8),j=1,5)/(0.1D0,0.1D0) , (4.0D0,5.0D0) ,&
     &      (4.0D0,5.0D0) , (4.0D0,5.0D0) , (4.0D0,5.0D0) ,             &
     &      (4.0D0,5.0D0) , (4.0D0,5.0D0) , (4.0D0,5.0D0) ,             &
     &      (-0.16D0,-0.37D0) , (6.0D0,7.0D0) , (6.0D0,7.0D0) ,         &
     &      (6.0D0,7.0D0) , (6.0D0,7.0D0) , (6.0D0,7.0D0) ,             &
     &      (6.0D0,7.0D0) , (6.0D0,7.0D0) , (-0.17D0,-0.19D0) ,         &
     &      (8.0D0,9.0D0) , (0.13D0,-0.39D0) , (2.0D0,5.0D0) ,          &
     &      (2.0D0,5.0D0) , (2.0D0,5.0D0) , (2.0D0,5.0D0) ,             &
     &      (2.0D0,5.0D0) , (0.11D0,-0.03D0) , (3.0D0,6.0D0) ,          &
     &      (-0.17D0,0.46D0) , (4.0D0,7.0D0) , (-0.17D0,-0.19D0) ,      &
     &      (7.0D0,2.0D0) , (7.0D0,2.0D0) , (7.0D0,2.0D0) ,             &
     &      (0.19D0,-0.17D0) , (5.0D0,8.0D0) , (0.20D0,-0.35D0) ,       &
     &      (6.0D0,9.0D0) , (0.35D0,0.20D0) , (8.0D0,3.0D0) ,           &
     &      (0.14D0,0.08D0) , (9.0D0,4.0D0)/
      DATA ((ctrue6(i,j,1),i=1,8),j=1,5)/(0.1D0,0.1D0) , (1.0D0,2.0D0) ,&
     &      (1.0D0,2.0D0) , (1.0D0,2.0D0) , (1.0D0,2.0D0) ,             &
     &      (1.0D0,2.0D0) , (1.0D0,2.0D0) , (1.0D0,2.0D0) ,             &
     &      (0.09D0,-0.12D0) , (3.0D0,4.0D0) , (3.0D0,4.0D0) ,          &
     &      (3.0D0,4.0D0) , (3.0D0,4.0D0) , (3.0D0,4.0D0) ,             &
     &      (3.0D0,4.0D0) , (3.0D0,4.0D0) , (0.03D0,-0.09D0) ,          &
     &      (0.15D0,-0.03D0) , (5.0D0,6.0D0) , (5.0D0,6.0D0) ,          &
     &      (5.0D0,6.0D0) , (5.0D0,6.0D0) , (5.0D0,6.0D0) ,             &
     &      (5.0D0,6.0D0) , (0.03D0,0.03D0) , (-0.18D0,0.03D0) ,        &
     &      (0.03D0,-0.09D0) , (7.0D0,8.0D0) , (7.0D0,8.0D0) ,          &
     &      (7.0D0,8.0D0) , (7.0D0,8.0D0) , (7.0D0,8.0D0) ,             &
     &      (0.09D0,0.03D0) , (0.15D0,0.00D0) , (0.00D0,0.15D0) ,       &
     &      (0.00D0,0.06D0) , (2.0D0,3.0D0) , (2.0D0,3.0D0) ,           &
     &      (2.0D0,3.0D0) , (2.0D0,3.0D0)/
      DATA ((ctrue6(i,j,2),i=1,8),j=1,5)/(0.1D0,0.1D0) , (4.0D0,5.0D0) ,&
     &      (4.0D0,5.0D0) , (4.0D0,5.0D0) , (4.0D0,5.0D0) ,             &
     &      (4.0D0,5.0D0) , (4.0D0,5.0D0) , (4.0D0,5.0D0) ,             &
     &      (0.09D0,-0.12D0) , (6.0D0,7.0D0) , (6.0D0,7.0D0) ,          &
     &      (6.0D0,7.0D0) , (6.0D0,7.0D0) , (6.0D0,7.0D0) ,             &
     &      (6.0D0,7.0D0) , (6.0D0,7.0D0) , (0.03D0,-0.09D0) ,          &
     &      (8.0D0,9.0D0) , (0.15D0,-0.03D0) , (2.0D0,5.0D0) ,          &
     &      (2.0D0,5.0D0) , (2.0D0,5.0D0) , (2.0D0,5.0D0) ,             &
     &      (2.0D0,5.0D0) , (0.03D0,0.03D0) , (3.0D0,6.0D0) ,           &
     &      (-0.18D0,0.03D0) , (4.0D0,7.0D0) , (0.03D0,-0.09D0) ,       &
     &      (7.0D0,2.0D0) , (7.0D0,2.0D0) , (7.0D0,2.0D0) ,             &
     &      (0.09D0,0.03D0) , (5.0D0,8.0D0) , (0.15D0,0.00D0) ,         &
     &      (6.0D0,9.0D0) , (0.00D0,0.15D0) , (8.0D0,3.0D0) ,           &
     &      (0.00D0,0.06D0) , (9.0D0,4.0D0)/
      DATA itrue3/0 , 1 , 2 , 2 , 2/
      DATA itruec/0 , 1 , 1 , 1 , 1/
!     .. Executable Statements ..
      DO INCx = 1 , 2
         DO np1 = 1 , 5
            N = np1 - 1
            len = 2*MAX(N,1)
!           .. Set vector arguments ..
            DO i = 1 , len
               cx(i) = cv(i,np1,INCx)
            ENDDO
            IF ( ICAse==6 ) THEN
!              .. DZNRM2 ..
               CALL STEST1(DZNRM2(N,cx,INCx),strue2(np1),strue2(np1),   &
     &                     Sfac)
            ELSEIF ( ICAse==7 ) THEN
!              .. DZASUM ..
               CALL STEST1(DZASUM(N,cx,INCx),strue4(np1),strue4(np1),   &
     &                     Sfac)
            ELSEIF ( ICAse==8 ) THEN
!              .. ZSCAL ..
               CALL ZSCAL(N,ca,cx,INCx)
               CALL CTEST(len,cx,ctrue5(1,np1,INCx),ctrue5(1,np1,INCx), &
     &                    Sfac)
            ELSEIF ( ICAse==9 ) THEN
!              .. ZDSCAL ..
               CALL ZDSCAL(N,sa,cx,INCx)
               CALL CTEST(len,cx,ctrue6(1,np1,INCx),ctrue6(1,np1,INCx), &
     &                    Sfac)
            ELSEIF ( ICAse==10 ) THEN
!              .. IZAMAX ..
               CALL ITEST1(IZAMAX(N,cx,INCx),itrue3(np1))
               DO i = 1 , len
                  cx(i) = (42.0D0,43.0D0)
               ENDDO
               CALL ITEST1(IZAMAX(N,cx,INCx),itruec(np1))
            ELSE
               WRITE (NOUT,*) ' Shouldn''t be here in CHECK1'
               STOP
            ENDIF
!
         ENDDO
         IF ( ICAse==10 ) THEN
            N = 8
            ix = 1
            DO i = 1 , N
               cxr(ix) = cvr(i)
               ix = ix + INCx
            ENDDO
            CALL ITEST1(IZAMAX(N,cxr,INCx),3)
         ENDIF
      ENDDO
!
      INCx = 1
      IF ( ICAse==8 ) THEN
!        ZSCAL
!        Add a test for alpha equal to zero.
         ca = (0.0D0,0.0D0)
         DO i = 1 , 5
            mwpct(i) = (0.0D0,0.0D0)
            mwpcs(i) = (1.0D0,1.0D0)
         ENDDO
         CALL ZSCAL(5,ca,cx,INCx)
         CALL CTEST(5,cx,mwpct,mwpcs,Sfac)
      ELSEIF ( ICAse==9 ) THEN
!        ZDSCAL
!        Add a test for alpha equal to zero.
         sa = 0.0D0
         DO i = 1 , 5
            mwpct(i) = (0.0D0,0.0D0)
            mwpcs(i) = (1.0D0,1.0D0)
         ENDDO
         CALL ZDSCAL(5,sa,cx,INCx)
         CALL CTEST(5,cx,mwpct,mwpcs,Sfac)
!        Add a test for alpha equal to one.
         sa = 1.0D0
         DO i = 1 , 5
            mwpct(i) = cx(i)
            mwpcs(i) = cx(i)
         ENDDO
         CALL ZDSCAL(5,sa,cx,INCx)
         CALL CTEST(5,cx,mwpct,mwpcs,Sfac)
!        Add a test for alpha equal to minus one.
         sa = -1.0D0
         DO i = 1 , 5
            mwpct(i) = -cx(i)
            mwpcs(i) = -cx(i)
         ENDDO
         CALL ZDSCAL(5,sa,cx,INCx)
         CALL CTEST(5,cx,mwpct,mwpcs,Sfac)
      ENDIF
      END SUBROUTINE CHECK1
!*==check2.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
      SUBROUTINE CHECK2(Sfac)
      IMPLICIT NONE
!*--CHECK2342
!     .. Parameters ..
      INTEGER NOUT
      PARAMETER (NOUT=6)
!     .. Scalar Arguments ..
      DOUBLE PRECISION Sfac
!     .. Scalars in Common ..
      INTEGER ICAse , INCx , INCy , MODe , N
      LOGICAL PASs
!     .. Local Scalars ..
      COMPLEX*16 ca
      INTEGER i , j , ki , kn , ksize , lenx , leny , lincx , lincy ,   &
     &        mx , my
!     .. Local Arrays ..
      COMPLEX*16 cdot(1) , csize1(4) , csize2(7,2) , csize3(14) ,       &
     &           ct10x(7,4,4) , ct10y(7,4,4) , ct6(4,4) , ct7(4,4) ,    &
     &           ct8(7,4,4) , cty0(1) , cx(7) , cx0(1) , cx1(7) ,       &
     &           cy(7) , cy0(1) , cy1(7)
      INTEGER incxs(4) , incys(4) , lens(4,2) , ns(4)
!     .. External Functions ..
      COMPLEX*16 ZDOTC , ZDOTU
      EXTERNAL ZDOTC , ZDOTU
!     .. External Subroutines ..
      EXTERNAL ZAXPY , ZCOPY , ZSWAP , CTEST
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MIN
!     .. Common blocks ..
      COMMON /COMBLA/ ICAse , N , INCx , INCy , MODe , PASs
!     .. Data statements ..
      DATA ca/(0.4D0,-0.7D0)/
      DATA incxs/1 , 2 , -2 , -1/
      DATA incys/1 , -2 , 1 , -2/
      DATA lens/1 , 1 , 2 , 4 , 1 , 1 , 3 , 7/
      DATA ns/0 , 1 , 2 , 4/
      DATA cx1/(0.7D0,-0.8D0) , (-0.4D0,-0.7D0) , (-0.1D0,-0.9D0) ,     &
     &     (0.2D0,-0.8D0) , (-0.9D0,-0.4D0) , (0.1D0,0.4D0) ,           &
     &     (-0.6D0,0.6D0)/
      DATA cy1/(0.6D0,-0.6D0) , (-0.9D0,0.5D0) , (0.7D0,-0.6D0) ,       &
     &     (0.1D0,-0.5D0) , (-0.1D0,-0.2D0) , (-0.5D0,-0.3D0) ,         &
     &     (0.8D0,-0.7D0)/
      DATA ((ct8(i,j,1),i=1,7),j=1,4)/(0.6D0,-0.6D0) , (0.0D0,0.0D0) ,  &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0) ,             &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.32D0,-1.41D0) ,          &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0) ,             &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0) ,             &
     &      (0.32D0,-1.41D0) , (-1.55D0,0.5D0) , (0.0D0,0.0D0) ,        &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0) ,             &
     &      (0.0D0,0.0D0) , (0.32D0,-1.41D0) , (-1.55D0,0.5D0) ,        &
     &      (0.03D0,-0.89D0) , (-0.38D0,-0.96D0) , (0.0D0,0.0D0) ,      &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0)/
      DATA ((ct8(i,j,2),i=1,7),j=1,4)/(0.6D0,-0.6D0) , (0.0D0,0.0D0) ,  &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0) ,             &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.32D0,-1.41D0) ,          &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0) ,             &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0) ,             &
     &      (-0.07D0,-0.89D0) , (-0.9D0,0.5D0) , (0.42D0,-1.41D0) ,     &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0) ,             &
     &      (0.0D0,0.0D0) , (0.78D0,0.06D0) , (-0.9D0,0.5D0) ,          &
     &      (0.06D0,-0.13D0) , (0.1D0,-0.5D0) , (-0.77D0,-0.49D0) ,     &
     &      (-0.5D0,-0.3D0) , (0.52D0,-1.51D0)/
      DATA ((ct8(i,j,3),i=1,7),j=1,4)/(0.6D0,-0.6D0) , (0.0D0,0.0D0) ,  &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0) ,             &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.32D0,-1.41D0) ,          &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0) ,             &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0) ,             &
     &      (-0.07D0,-0.89D0) , (-1.18D0,-0.31D0) , (0.0D0,0.0D0) ,     &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0) ,             &
     &      (0.0D0,0.0D0) , (0.78D0,0.06D0) , (-1.54D0,0.97D0) ,        &
     &      (0.03D0,-0.89D0) , (-0.18D0,-1.31D0) , (0.0D0,0.0D0) ,      &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0)/
      DATA ((ct8(i,j,4),i=1,7),j=1,4)/(0.6D0,-0.6D0) , (0.0D0,0.0D0) ,  &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0) ,             &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.32D0,-1.41D0) ,          &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0) ,             &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0) ,             &
     &      (0.32D0,-1.41D0) , (-0.9D0,0.5D0) , (0.05D0,-0.6D0) ,       &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0) ,             &
     &      (0.0D0,0.0D0) , (0.32D0,-1.41D0) , (-0.9D0,0.5D0) ,         &
     &      (0.05D0,-0.6D0) , (0.1D0,-0.5D0) , (-0.77D0,-0.49D0) ,      &
     &      (-0.5D0,-0.3D0) , (0.32D0,-1.16D0)/
      DATA ct7/(0.0D0,0.0D0) , (-0.06D0,-0.90D0) , (0.65D0,-0.47D0) ,   &
     &     (-0.34D0,-1.22D0) , (0.0D0,0.0D0) , (-0.06D0,-0.90D0) ,      &
     &     (-0.59D0,-1.46D0) , (-1.04D0,-0.04D0) , (0.0D0,0.0D0) ,      &
     &     (-0.06D0,-0.90D0) , (-0.83D0,0.59D0) , (0.07D0,-0.37D0) ,    &
     &     (0.0D0,0.0D0) , (-0.06D0,-0.90D0) , (-0.76D0,-1.15D0) ,      &
     &     (-1.33D0,-1.82D0)/
      DATA ct6/(0.0D0,0.0D0) , (0.90D0,0.06D0) , (0.91D0,-0.77D0) ,     &
     &     (1.80D0,-0.10D0) , (0.0D0,0.0D0) , (0.90D0,0.06D0) ,         &
     &     (1.45D0,0.74D0) , (0.20D0,0.90D0) , (0.0D0,0.0D0) ,          &
     &     (0.90D0,0.06D0) , (-0.55D0,0.23D0) , (0.83D0,-0.39D0) ,      &
     &     (0.0D0,0.0D0) , (0.90D0,0.06D0) , (1.04D0,0.79D0) ,          &
     &     (1.95D0,1.22D0)/
      DATA ((ct10x(i,j,1),i=1,7),j=1,4)/(0.7D0,-0.8D0) , (0.0D0,0.0D0) ,&
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0) ,             &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.6D0,-0.6D0) ,            &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0) ,             &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0) ,             &
     &      (0.6D0,-0.6D0) , (-0.9D0,0.5D0) , (0.0D0,0.0D0) ,           &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0) ,             &
     &      (0.0D0,0.0D0) , (0.6D0,-0.6D0) , (-0.9D0,0.5D0) ,           &
     &      (0.7D0,-0.6D0) , (0.1D0,-0.5D0) , (0.0D0,0.0D0) ,           &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0)/
      DATA ((ct10x(i,j,2),i=1,7),j=1,4)/(0.7D0,-0.8D0) , (0.0D0,0.0D0) ,&
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0) ,             &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.6D0,-0.6D0) ,            &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0) ,             &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0) ,             &
     &      (0.7D0,-0.6D0) , (-0.4D0,-0.7D0) , (0.6D0,-0.6D0) ,         &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0) ,             &
     &      (0.0D0,0.0D0) , (0.8D0,-0.7D0) , (-0.4D0,-0.7D0) ,          &
     &      (-0.1D0,-0.2D0) , (0.2D0,-0.8D0) , (0.7D0,-0.6D0) ,         &
     &      (0.1D0,0.4D0) , (0.6D0,-0.6D0)/
      DATA ((ct10x(i,j,3),i=1,7),j=1,4)/(0.7D0,-0.8D0) , (0.0D0,0.0D0) ,&
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0) ,             &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.6D0,-0.6D0) ,            &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0) ,             &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0) ,             &
     &      (-0.9D0,0.5D0) , (-0.4D0,-0.7D0) , (0.6D0,-0.6D0) ,         &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0) ,             &
     &      (0.0D0,0.0D0) , (0.1D0,-0.5D0) , (-0.4D0,-0.7D0) ,          &
     &      (0.7D0,-0.6D0) , (0.2D0,-0.8D0) , (-0.9D0,0.5D0) ,          &
     &      (0.1D0,0.4D0) , (0.6D0,-0.6D0)/
      DATA ((ct10x(i,j,4),i=1,7),j=1,4)/(0.7D0,-0.8D0) , (0.0D0,0.0D0) ,&
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0) ,             &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.6D0,-0.6D0) ,            &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0) ,             &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0) ,             &
     &      (0.6D0,-0.6D0) , (0.7D0,-0.6D0) , (0.0D0,0.0D0) ,           &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0) ,             &
     &      (0.0D0,0.0D0) , (0.6D0,-0.6D0) , (0.7D0,-0.6D0) ,           &
     &      (-0.1D0,-0.2D0) , (0.8D0,-0.7D0) , (0.0D0,0.0D0) ,          &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0)/
      DATA ((ct10y(i,j,1),i=1,7),j=1,4)/(0.6D0,-0.6D0) , (0.0D0,0.0D0) ,&
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0) ,             &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.7D0,-0.8D0) ,            &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0) ,             &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0) ,             &
     &      (0.7D0,-0.8D0) , (-0.4D0,-0.7D0) , (0.0D0,0.0D0) ,          &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0) ,             &
     &      (0.0D0,0.0D0) , (0.7D0,-0.8D0) , (-0.4D0,-0.7D0) ,          &
     &      (-0.1D0,-0.9D0) , (0.2D0,-0.8D0) , (0.0D0,0.0D0) ,          &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0)/
      DATA ((ct10y(i,j,2),i=1,7),j=1,4)/(0.6D0,-0.6D0) , (0.0D0,0.0D0) ,&
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0) ,             &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.7D0,-0.8D0) ,            &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0) ,             &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0) ,             &
     &      (-0.1D0,-0.9D0) , (-0.9D0,0.5D0) , (0.7D0,-0.8D0) ,         &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0) ,             &
     &      (0.0D0,0.0D0) , (-0.6D0,0.6D0) , (-0.9D0,0.5D0) ,           &
     &      (-0.9D0,-0.4D0) , (0.1D0,-0.5D0) , (-0.1D0,-0.9D0) ,        &
     &      (-0.5D0,-0.3D0) , (0.7D0,-0.8D0)/
      DATA ((ct10y(i,j,3),i=1,7),j=1,4)/(0.6D0,-0.6D0) , (0.0D0,0.0D0) ,&
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0) ,             &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.7D0,-0.8D0) ,            &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0) ,             &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0) ,             &
     &      (-0.1D0,-0.9D0) , (0.7D0,-0.8D0) , (0.0D0,0.0D0) ,          &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0) ,             &
     &      (0.0D0,0.0D0) , (-0.6D0,0.6D0) , (-0.9D0,-0.4D0) ,          &
     &      (-0.1D0,-0.9D0) , (0.7D0,-0.8D0) , (0.0D0,0.0D0) ,          &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0)/
      DATA ((ct10y(i,j,4),i=1,7),j=1,4)/(0.6D0,-0.6D0) , (0.0D0,0.0D0) ,&
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0) ,             &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.7D0,-0.8D0) ,            &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0) ,             &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0) ,             &
     &      (0.7D0,-0.8D0) , (-0.9D0,0.5D0) , (-0.4D0,-0.7D0) ,         &
     &      (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0) ,             &
     &      (0.0D0,0.0D0) , (0.7D0,-0.8D0) , (-0.9D0,0.5D0) ,           &
     &      (-0.4D0,-0.7D0) , (0.1D0,-0.5D0) , (-0.1D0,-0.9D0) ,        &
     &      (-0.5D0,-0.3D0) , (0.2D0,-0.8D0)/
      DATA csize1/(0.0D0,0.0D0) , (0.9D0,0.9D0) , (1.63D0,1.73D0) ,     &
     &     (2.90D0,2.78D0)/
      DATA csize3/(0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0) ,       &
     &     (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0)&
     &     , (1.17D0,1.17D0) , (1.17D0,1.17D0) , (1.17D0,1.17D0) ,      &
     &     (1.17D0,1.17D0) , (1.17D0,1.17D0) , (1.17D0,1.17D0) ,        &
     &     (1.17D0,1.17D0)/
      DATA csize2/(0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0) ,       &
     &     (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0) , (0.0D0,0.0D0)&
     &     , (1.54D0,1.54D0) , (1.54D0,1.54D0) , (1.54D0,1.54D0) ,      &
     &     (1.54D0,1.54D0) , (1.54D0,1.54D0) , (1.54D0,1.54D0) ,        &
     &     (1.54D0,1.54D0)/
!     .. Executable Statements ..
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
!           .. initialize all argument arrays ..
            DO i = 1 , 7
               cx(i) = cx1(i)
               cy(i) = cy1(i)
            ENDDO
            IF ( ICAse==1 ) THEN
!              .. ZDOTC ..
               cdot(1) = ZDOTC(N,cx,INCx,cy,INCy)
               CALL CTEST(1,cdot,ct6(kn,ki),csize1(kn),Sfac)
            ELSEIF ( ICAse==2 ) THEN
!              .. ZDOTU ..
               cdot(1) = ZDOTU(N,cx,INCx,cy,INCy)
               CALL CTEST(1,cdot,ct7(kn,ki),csize1(kn),Sfac)
            ELSEIF ( ICAse==3 ) THEN
!              .. ZAXPY ..
               CALL ZAXPY(N,ca,cx,INCx,cy,INCy)
               CALL CTEST(leny,cy,ct8(1,kn,ki),csize2(1,ksize),Sfac)
            ELSEIF ( ICAse==4 ) THEN
!              .. ZCOPY ..
               CALL ZCOPY(N,cx,INCx,cy,INCy)
               CALL CTEST(leny,cy,ct10y(1,kn,ki),csize3,1.0D0)
               IF ( ki==1 ) THEN
                  cx0(1) = (42.0D0,43.0D0)
                  cy0(1) = (44.0D0,45.0D0)
                  IF ( N==0 ) THEN
                     cty0(1) = cy0(1)
                  ELSE
                     cty0(1) = cx0(1)
                  ENDIF
                  lincx = INCx
                  INCx = 0
                  lincy = INCy
                  INCy = 0
                  CALL ZCOPY(N,cx0,INCx,cy0,INCy)
                  CALL CTEST(1,cy0,cty0,csize3,1.0D0)
                  INCx = lincx
                  INCy = lincy
               ENDIF
            ELSEIF ( ICAse==5 ) THEN
!              .. ZSWAP ..
               CALL ZSWAP(N,cx,INCx,cy,INCy)
               CALL CTEST(lenx,cx,ct10x(1,kn,ki),csize3,1.0D0)
               CALL CTEST(leny,cy,ct10y(1,kn,ki),csize3,1.0D0)
            ELSE
               WRITE (NOUT,*) ' Shouldn''t be here in CHECK2'
               STOP
            ENDIF
!
         ENDDO
      ENDDO
      END SUBROUTINE CHECK2
!*==stest.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
      SUBROUTINE STEST(Len,Scomp,Strue,Ssize,Sfac)
      IMPLICIT NONE
!*--STEST592
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
      DOUBLE PRECISION ZERO
      PARAMETER (NOUT=6,ZERO=0.0D0)
!     .. Scalar Arguments ..
      DOUBLE PRECISION Sfac
      INTEGER Len
!     .. Array Arguments ..
      DOUBLE PRECISION Scomp(Len) , Ssize(Len) , Strue(Len)
!     .. Scalars in Common ..
      INTEGER ICAse , INCx , INCy , MODe , N
      LOGICAL PASs
!     .. Local Scalars ..
      DOUBLE PRECISION sd
      INTEGER i
!     .. External Functions ..
      DOUBLE PRECISION SDIFF
      EXTERNAL SDIFF
!     .. Intrinsic Functions ..
      INTRINSIC ABS
!     .. Common blocks ..
      COMMON /COMBLA/ ICAse , N , INCx , INCy , MODe , PASs
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
     &          ' CASE  N INCX INCY MODE  I                            '&
     &          ,                                                       &
     &        ' COMP(I)                             TRUE(I)  DIFFERENCE'&
     &        ,'     SIZE(I)',/1X)
            ENDIF
            WRITE (NOUT,99003) ICAse , N , INCx , INCy , MODe , i ,     &
     &                         Scomp(i) , Strue(i) , sd , Ssize(i)
99003       FORMAT (1X,I4,I3,3I5,I3,2D36.8,2D12.4)
         ENDIF
      ENDDO
      RETURN
      END SUBROUTINE STEST
!*==stest1.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
      SUBROUTINE STEST1(Scomp1,Strue1,Ssize,Sfac)
      IMPLICIT NONE
!*--STEST1654
!     ************************* STEST1 *****************************
!
!     THIS IS AN INTERFACE SUBROUTINE TO ACCOMMODATE THE FORTRAN
!     REQUIREMENT THAT WHEN A DUMMY ARGUMENT IS AN ARRAY, THE
!     ACTUAL ARGUMENT MUST ALSO BE AN ARRAY OR AN ARRAY ELEMENT.
!
!     C.L. LAWSON, JPL, 1978 DEC 6
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION Scomp1 , Sfac , Strue1
!     .. Array Arguments ..
      DOUBLE PRECISION Ssize(*)
!     .. Local Arrays ..
      DOUBLE PRECISION scomp(1) , strue(1)
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
      DOUBLE PRECISION FUNCTION SDIFF(Sa,Sb)
      IMPLICIT NONE
!*--SDIFF681
!     ********************************* SDIFF **************************
!     COMPUTES DIFFERENCE OF TWO NUMBERS.  C. L. LAWSON, JPL 1974 FEB 15
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION Sa , Sb
!     .. Executable Statements ..
      SDIFF = Sa - Sb
      END FUNCTION SDIFF
!*==ctest.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
      SUBROUTINE CTEST(Len,Ccomp,Ctrue,Csize,Sfac)
      IMPLICIT NONE
!*--CTEST693
!     **************************** CTEST *****************************
!
!     C.L. LAWSON, JPL, 1978 DEC 6
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION Sfac
      INTEGER Len
!     .. Array Arguments ..
      COMPLEX*16 Ccomp(Len) , Csize(Len) , Ctrue(Len)
!     .. Local Scalars ..
      INTEGER i
!     .. Local Arrays ..
      DOUBLE PRECISION scomp(20) , ssize(20) , strue(20)
!     .. External Subroutines ..
      EXTERNAL STEST
!     .. Intrinsic Functions ..
      INTRINSIC DIMAG , DBLE
!     .. Executable Statements ..
      DO i = 1 , Len
         scomp(2*i-1) = DBLE(Ccomp(i))
         scomp(2*i) = DIMAG(Ccomp(i))
         strue(2*i-1) = DBLE(Ctrue(i))
         strue(2*i) = DIMAG(Ctrue(i))
         ssize(2*i-1) = DBLE(Csize(i))
         ssize(2*i) = DIMAG(Csize(i))
      ENDDO
!
      CALL STEST(2*Len,scomp,strue,ssize,Sfac)
      END SUBROUTINE CTEST
!*==itest1.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
      SUBROUTINE ITEST1(Icomp,Itrue)
      IMPLICIT NONE
!*--ITEST1726
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
      INTEGER ICAse , INCx , INCy , MODe , N
      LOGICAL PASs
!     .. Local Scalars ..
      INTEGER id
!     .. Common blocks ..
      COMMON /COMBLA/ ICAse , N , INCx , INCy , MODe , PASs
!     .. Executable Statements ..
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
99002       FORMAT (/                                                   &
     &          ' CASE  N INCX INCY MODE                               '&
     &          ,                                                       &
     &        ' COMP                                TRUE     DIFFERENCE'&
     &        ,/1X)
         ENDIF
         id = Icomp - Itrue
         WRITE (NOUT,99003) ICAse , N , INCx , INCy , MODe , Icomp ,    &
     &                      Itrue , id
99003    FORMAT (1X,I4,I3,3I5,2I36,I12)
      ENDIF
      RETURN
      END SUBROUTINE ITEST1
