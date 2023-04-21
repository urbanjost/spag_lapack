!*==cblat1.f90  processed by SPAG 7.51RB at 17:52 on  4 Mar 2022
!> \brief \b CBLAT1
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       PROGRAM CBLAT1
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    Test program for the COMPLEX Level 1 BLAS.
!>    Based upon the original BLAS test routine together with:
!>
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
!> \ingroup complex_blas_testing
!
!  =====================================================================
      PROGRAM CBLAT1
      IMPLICIT NONE
!*--CBLAT141
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
      REAL sfac
      INTEGER ic
!     .. External Subroutines ..
      EXTERNAL CHECK1 , CHECK2 , HEADER
!     .. Common blocks ..
      COMMON /COMBLA/ ICAse , N , INCx , INCy , MODe , PASs
!     .. Data statements ..
      DATA sfac/9.765625E-4/
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
      END PROGRAM CBLAT1
!*==header.f90  processed by SPAG 7.51RB at 17:52 on  4 Mar 2022
      SUBROUTINE HEADER
      IMPLICIT NONE
!*--HEADER96
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
      DATA l(1)/'CDOTC '/
      DATA l(2)/'CDOTU '/
      DATA l(3)/'CAXPY '/
      DATA l(4)/'CCOPY '/
      DATA l(5)/'CSWAP '/
      DATA l(6)/'SCNRM2'/
      DATA l(7)/'SCASUM'/
      DATA l(8)/'CSCAL '/
      DATA l(9)/'CSSCAL'/
      DATA l(10)/'ICAMAX'/
!     .. Executable Statements ..
      WRITE (NOUT,99001) ICAse , l(ICAse)
!
99001 FORMAT (/' Test of subprogram number',I3,12X,A6)
      RETURN
      END SUBROUTINE HEADER
!*==check1.f90  processed by SPAG 7.51RB at 17:52 on  4 Mar 2022
      SUBROUTINE CHECK1(Sfac)
      IMPLICIT NONE
!*--CHECK1127
!     .. Parameters ..
      INTEGER NOUT
      PARAMETER (NOUT=6)
!     .. Scalar Arguments ..
      REAL Sfac
!     .. Scalars in Common ..
      INTEGER ICAse , INCx , INCy , MODe , N
      LOGICAL PASs
!     .. Local Scalars ..
      COMPLEX ca
      REAL sa
      INTEGER i , ix , j , len , np1
!     .. Local Arrays ..
      COMPLEX ctrue5(8,5,2) , ctrue6(8,5,2) , cv(8,5,2) , cvr(8) ,      &
     &        cx(8) , cxr(15) , mwpcs(5) , mwpct(5)
      REAL strue2(5) , strue4(5)
      INTEGER itrue3(5) , itruec(5)
!     .. External Functions ..
      REAL SCASUM , SCNRM2
      INTEGER ICAMAX
      EXTERNAL SCASUM , SCNRM2 , ICAMAX
!     .. External Subroutines ..
      EXTERNAL CSCAL , CSSCAL , CTEST , ITEST1 , STEST1
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     .. Common blocks ..
      COMMON /COMBLA/ ICAse , N , INCx , INCy , MODe , PASs
!     .. Data statements ..
      DATA sa , ca/0.3E0 , (0.4E0,-0.7E0)/
      DATA ((cv(i,j,1),i=1,8),j=1,5)/(0.1E0,0.1E0) , (1.0E0,2.0E0) ,    &
     &      (1.0E0,2.0E0) , (1.0E0,2.0E0) , (1.0E0,2.0E0) ,             &
     &      (1.0E0,2.0E0) , (1.0E0,2.0E0) , (1.0E0,2.0E0) ,             &
     &      (0.3E0,-0.4E0) , (3.0E0,4.0E0) , (3.0E0,4.0E0) ,            &
     &      (3.0E0,4.0E0) , (3.0E0,4.0E0) , (3.0E0,4.0E0) ,             &
     &      (3.0E0,4.0E0) , (3.0E0,4.0E0) , (0.1E0,-0.3E0) ,            &
     &      (0.5E0,-0.1E0) , (5.0E0,6.0E0) , (5.0E0,6.0E0) ,            &
     &      (5.0E0,6.0E0) , (5.0E0,6.0E0) , (5.0E0,6.0E0) ,             &
     &      (5.0E0,6.0E0) , (0.1E0,0.1E0) , (-0.6E0,0.1E0) ,            &
     &      (0.1E0,-0.3E0) , (7.0E0,8.0E0) , (7.0E0,8.0E0) ,            &
     &      (7.0E0,8.0E0) , (7.0E0,8.0E0) , (7.0E0,8.0E0) ,             &
     &      (0.3E0,0.1E0) , (0.5E0,0.0E0) , (0.0E0,0.5E0) ,             &
     &      (0.0E0,0.2E0) , (2.0E0,3.0E0) , (2.0E0,3.0E0) ,             &
     &      (2.0E0,3.0E0) , (2.0E0,3.0E0)/
      DATA ((cv(i,j,2),i=1,8),j=1,5)/(0.1E0,0.1E0) , (4.0E0,5.0E0) ,    &
     &      (4.0E0,5.0E0) , (4.0E0,5.0E0) , (4.0E0,5.0E0) ,             &
     &      (4.0E0,5.0E0) , (4.0E0,5.0E0) , (4.0E0,5.0E0) ,             &
     &      (0.3E0,-0.4E0) , (6.0E0,7.0E0) , (6.0E0,7.0E0) ,            &
     &      (6.0E0,7.0E0) , (6.0E0,7.0E0) , (6.0E0,7.0E0) ,             &
     &      (6.0E0,7.0E0) , (6.0E0,7.0E0) , (0.1E0,-0.3E0) ,            &
     &      (8.0E0,9.0E0) , (0.5E0,-0.1E0) , (2.0E0,5.0E0) ,            &
     &      (2.0E0,5.0E0) , (2.0E0,5.0E0) , (2.0E0,5.0E0) ,             &
     &      (2.0E0,5.0E0) , (0.1E0,0.1E0) , (3.0E0,6.0E0) ,             &
     &      (-0.6E0,0.1E0) , (4.0E0,7.0E0) , (0.1E0,-0.3E0) ,           &
     &      (7.0E0,2.0E0) , (7.0E0,2.0E0) , (7.0E0,2.0E0) ,             &
     &      (0.3E0,0.1E0) , (5.0E0,8.0E0) , (0.5E0,0.0E0) ,             &
     &      (6.0E0,9.0E0) , (0.0E0,0.5E0) , (8.0E0,3.0E0) ,             &
     &      (0.0E0,0.2E0) , (9.0E0,4.0E0)/
      DATA cvr/(8.0E0,8.0E0) , (-7.0E0,-7.0E0) , (9.0E0,9.0E0) ,        &
     &     (5.0E0,5.0E0) , (9.0E0,9.0E0) , (8.0E0,8.0E0) , (7.0E0,7.0E0)&
     &     , (7.0E0,7.0E0)/
      DATA strue2/0.0E0 , 0.5E0 , 0.6E0 , 0.7E0 , 0.8E0/
      DATA strue4/0.0E0 , 0.7E0 , 1.0E0 , 1.3E0 , 1.6E0/
      DATA ((ctrue5(i,j,1),i=1,8),j=1,5)/(0.1E0,0.1E0) , (1.0E0,2.0E0) ,&
     &      (1.0E0,2.0E0) , (1.0E0,2.0E0) , (1.0E0,2.0E0) ,             &
     &      (1.0E0,2.0E0) , (1.0E0,2.0E0) , (1.0E0,2.0E0) ,             &
     &      (-0.16E0,-0.37E0) , (3.0E0,4.0E0) , (3.0E0,4.0E0) ,         &
     &      (3.0E0,4.0E0) , (3.0E0,4.0E0) , (3.0E0,4.0E0) ,             &
     &      (3.0E0,4.0E0) , (3.0E0,4.0E0) , (-0.17E0,-0.19E0) ,         &
     &      (0.13E0,-0.39E0) , (5.0E0,6.0E0) , (5.0E0,6.0E0) ,          &
     &      (5.0E0,6.0E0) , (5.0E0,6.0E0) , (5.0E0,6.0E0) ,             &
     &      (5.0E0,6.0E0) , (0.11E0,-0.03E0) , (-0.17E0,0.46E0) ,       &
     &      (-0.17E0,-0.19E0) , (7.0E0,8.0E0) , (7.0E0,8.0E0) ,         &
     &      (7.0E0,8.0E0) , (7.0E0,8.0E0) , (7.0E0,8.0E0) ,             &
     &      (0.19E0,-0.17E0) , (0.20E0,-0.35E0) , (0.35E0,0.20E0) ,     &
     &      (0.14E0,0.08E0) , (2.0E0,3.0E0) , (2.0E0,3.0E0) ,           &
     &      (2.0E0,3.0E0) , (2.0E0,3.0E0)/
      DATA ((ctrue5(i,j,2),i=1,8),j=1,5)/(0.1E0,0.1E0) , (4.0E0,5.0E0) ,&
     &      (4.0E0,5.0E0) , (4.0E0,5.0E0) , (4.0E0,5.0E0) ,             &
     &      (4.0E0,5.0E0) , (4.0E0,5.0E0) , (4.0E0,5.0E0) ,             &
     &      (-0.16E0,-0.37E0) , (6.0E0,7.0E0) , (6.0E0,7.0E0) ,         &
     &      (6.0E0,7.0E0) , (6.0E0,7.0E0) , (6.0E0,7.0E0) ,             &
     &      (6.0E0,7.0E0) , (6.0E0,7.0E0) , (-0.17E0,-0.19E0) ,         &
     &      (8.0E0,9.0E0) , (0.13E0,-0.39E0) , (2.0E0,5.0E0) ,          &
     &      (2.0E0,5.0E0) , (2.0E0,5.0E0) , (2.0E0,5.0E0) ,             &
     &      (2.0E0,5.0E0) , (0.11E0,-0.03E0) , (3.0E0,6.0E0) ,          &
     &      (-0.17E0,0.46E0) , (4.0E0,7.0E0) , (-0.17E0,-0.19E0) ,      &
     &      (7.0E0,2.0E0) , (7.0E0,2.0E0) , (7.0E0,2.0E0) ,             &
     &      (0.19E0,-0.17E0) , (5.0E0,8.0E0) , (0.20E0,-0.35E0) ,       &
     &      (6.0E0,9.0E0) , (0.35E0,0.20E0) , (8.0E0,3.0E0) ,           &
     &      (0.14E0,0.08E0) , (9.0E0,4.0E0)/
      DATA ((ctrue6(i,j,1),i=1,8),j=1,5)/(0.1E0,0.1E0) , (1.0E0,2.0E0) ,&
     &      (1.0E0,2.0E0) , (1.0E0,2.0E0) , (1.0E0,2.0E0) ,             &
     &      (1.0E0,2.0E0) , (1.0E0,2.0E0) , (1.0E0,2.0E0) ,             &
     &      (0.09E0,-0.12E0) , (3.0E0,4.0E0) , (3.0E0,4.0E0) ,          &
     &      (3.0E0,4.0E0) , (3.0E0,4.0E0) , (3.0E0,4.0E0) ,             &
     &      (3.0E0,4.0E0) , (3.0E0,4.0E0) , (0.03E0,-0.09E0) ,          &
     &      (0.15E0,-0.03E0) , (5.0E0,6.0E0) , (5.0E0,6.0E0) ,          &
     &      (5.0E0,6.0E0) , (5.0E0,6.0E0) , (5.0E0,6.0E0) ,             &
     &      (5.0E0,6.0E0) , (0.03E0,0.03E0) , (-0.18E0,0.03E0) ,        &
     &      (0.03E0,-0.09E0) , (7.0E0,8.0E0) , (7.0E0,8.0E0) ,          &
     &      (7.0E0,8.0E0) , (7.0E0,8.0E0) , (7.0E0,8.0E0) ,             &
     &      (0.09E0,0.03E0) , (0.15E0,0.00E0) , (0.00E0,0.15E0) ,       &
     &      (0.00E0,0.06E0) , (2.0E0,3.0E0) , (2.0E0,3.0E0) ,           &
     &      (2.0E0,3.0E0) , (2.0E0,3.0E0)/
      DATA ((ctrue6(i,j,2),i=1,8),j=1,5)/(0.1E0,0.1E0) , (4.0E0,5.0E0) ,&
     &      (4.0E0,5.0E0) , (4.0E0,5.0E0) , (4.0E0,5.0E0) ,             &
     &      (4.0E0,5.0E0) , (4.0E0,5.0E0) , (4.0E0,5.0E0) ,             &
     &      (0.09E0,-0.12E0) , (6.0E0,7.0E0) , (6.0E0,7.0E0) ,          &
     &      (6.0E0,7.0E0) , (6.0E0,7.0E0) , (6.0E0,7.0E0) ,             &
     &      (6.0E0,7.0E0) , (6.0E0,7.0E0) , (0.03E0,-0.09E0) ,          &
     &      (8.0E0,9.0E0) , (0.15E0,-0.03E0) , (2.0E0,5.0E0) ,          &
     &      (2.0E0,5.0E0) , (2.0E0,5.0E0) , (2.0E0,5.0E0) ,             &
     &      (2.0E0,5.0E0) , (0.03E0,0.03E0) , (3.0E0,6.0E0) ,           &
     &      (-0.18E0,0.03E0) , (4.0E0,7.0E0) , (0.03E0,-0.09E0) ,       &
     &      (7.0E0,2.0E0) , (7.0E0,2.0E0) , (7.0E0,2.0E0) ,             &
     &      (0.09E0,0.03E0) , (5.0E0,8.0E0) , (0.15E0,0.00E0) ,         &
     &      (6.0E0,9.0E0) , (0.00E0,0.15E0) , (8.0E0,3.0E0) ,           &
     &      (0.00E0,0.06E0) , (9.0E0,4.0E0)/
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
!              .. SCNRM2 ..
               CALL STEST1(SCNRM2(N,cx,INCx),strue2(np1),strue2(np1),   &
     &                     Sfac)
            ELSEIF ( ICAse==7 ) THEN
!              .. SCASUM ..
               CALL STEST1(SCASUM(N,cx,INCx),strue4(np1),strue4(np1),   &
     &                     Sfac)
            ELSEIF ( ICAse==8 ) THEN
!              .. CSCAL ..
               CALL CSCAL(N,ca,cx,INCx)
               CALL CTEST(len,cx,ctrue5(1,np1,INCx),ctrue5(1,np1,INCx), &
     &                    Sfac)
            ELSEIF ( ICAse==9 ) THEN
!              .. CSSCAL ..
               CALL CSSCAL(N,sa,cx,INCx)
               CALL CTEST(len,cx,ctrue6(1,np1,INCx),ctrue6(1,np1,INCx), &
     &                    Sfac)
            ELSEIF ( ICAse==10 ) THEN
!              .. ICAMAX ..
               CALL ITEST1(ICAMAX(N,cx,INCx),itrue3(np1))
               DO i = 1 , len
                  cx(i) = (42.0E0,43.0E0)
               ENDDO
               CALL ITEST1(ICAMAX(N,cx,INCx),itruec(np1))
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
            CALL ITEST1(ICAMAX(N,cxr,INCx),3)
         ENDIF
      ENDDO
!
      INCx = 1
      IF ( ICAse==8 ) THEN
!        CSCAL
!        Add a test for alpha equal to zero.
         ca = (0.0E0,0.0E0)
         DO i = 1 , 5
            mwpct(i) = (0.0E0,0.0E0)
            mwpcs(i) = (1.0E0,1.0E0)
         ENDDO
         CALL CSCAL(5,ca,cx,INCx)
         CALL CTEST(5,cx,mwpct,mwpcs,Sfac)
      ELSEIF ( ICAse==9 ) THEN
!        CSSCAL
!        Add a test for alpha equal to zero.
         sa = 0.0E0
         DO i = 1 , 5
            mwpct(i) = (0.0E0,0.0E0)
            mwpcs(i) = (1.0E0,1.0E0)
         ENDDO
         CALL CSSCAL(5,sa,cx,INCx)
         CALL CTEST(5,cx,mwpct,mwpcs,Sfac)
!        Add a test for alpha equal to one.
         sa = 1.0E0
         DO i = 1 , 5
            mwpct(i) = cx(i)
            mwpcs(i) = cx(i)
         ENDDO
         CALL CSSCAL(5,sa,cx,INCx)
         CALL CTEST(5,cx,mwpct,mwpcs,Sfac)
!        Add a test for alpha equal to minus one.
         sa = -1.0E0
         DO i = 1 , 5
            mwpct(i) = -cx(i)
            mwpcs(i) = -cx(i)
         ENDDO
         CALL CSSCAL(5,sa,cx,INCx)
         CALL CTEST(5,cx,mwpct,mwpcs,Sfac)
      ENDIF
      END SUBROUTINE CHECK1
!*==check2.f90  processed by SPAG 7.51RB at 17:52 on  4 Mar 2022
      SUBROUTINE CHECK2(Sfac)
      IMPLICIT NONE
!*--CHECK2341
!     .. Parameters ..
      INTEGER NOUT
      PARAMETER (NOUT=6)
!     .. Scalar Arguments ..
      REAL Sfac
!     .. Scalars in Common ..
      INTEGER ICAse , INCx , INCy , MODe , N
      LOGICAL PASs
!     .. Local Scalars ..
      COMPLEX ca
      INTEGER i , j , ki , kn , ksize , lenx , leny , lincx , lincy ,   &
     &        mx , my
!     .. Local Arrays ..
      COMPLEX cdot(1) , csize1(4) , csize2(7,2) , csize3(14) ,          &
     &        ct10x(7,4,4) , ct10y(7,4,4) , ct6(4,4) , ct7(4,4) ,       &
     &        ct8(7,4,4) , cty0(1) , cx(7) , cx0(1) , cx1(7) , cy(7) ,  &
     &        cy0(1) , cy1(7)
      INTEGER incxs(4) , incys(4) , lens(4,2) , ns(4)
!     .. External Functions ..
      COMPLEX CDOTC , CDOTU
      EXTERNAL CDOTC , CDOTU
!     .. External Subroutines ..
      EXTERNAL CAXPY , CCOPY , CSWAP , CTEST
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MIN
!     .. Common blocks ..
      COMMON /COMBLA/ ICAse , N , INCx , INCy , MODe , PASs
!     .. Data statements ..
      DATA ca/(0.4E0,-0.7E0)/
      DATA incxs/1 , 2 , -2 , -1/
      DATA incys/1 , -2 , 1 , -2/
      DATA lens/1 , 1 , 2 , 4 , 1 , 1 , 3 , 7/
      DATA ns/0 , 1 , 2 , 4/
      DATA cx1/(0.7E0,-0.8E0) , (-0.4E0,-0.7E0) , (-0.1E0,-0.9E0) ,     &
     &     (0.2E0,-0.8E0) , (-0.9E0,-0.4E0) , (0.1E0,0.4E0) ,           &
     &     (-0.6E0,0.6E0)/
      DATA cy1/(0.6E0,-0.6E0) , (-0.9E0,0.5E0) , (0.7E0,-0.6E0) ,       &
     &     (0.1E0,-0.5E0) , (-0.1E0,-0.2E0) , (-0.5E0,-0.3E0) ,         &
     &     (0.8E0,-0.7E0)/
      DATA ((ct8(i,j,1),i=1,7),j=1,4)/(0.6E0,-0.6E0) , (0.0E0,0.0E0) ,  &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0) ,             &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.32E0,-1.41E0) ,          &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0) ,             &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0) ,             &
     &      (0.32E0,-1.41E0) , (-1.55E0,0.5E0) , (0.0E0,0.0E0) ,        &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0) ,             &
     &      (0.0E0,0.0E0) , (0.32E0,-1.41E0) , (-1.55E0,0.5E0) ,        &
     &      (0.03E0,-0.89E0) , (-0.38E0,-0.96E0) , (0.0E0,0.0E0) ,      &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0)/
      DATA ((ct8(i,j,2),i=1,7),j=1,4)/(0.6E0,-0.6E0) , (0.0E0,0.0E0) ,  &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0) ,             &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.32E0,-1.41E0) ,          &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0) ,             &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0) ,             &
     &      (-0.07E0,-0.89E0) , (-0.9E0,0.5E0) , (0.42E0,-1.41E0) ,     &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0) ,             &
     &      (0.0E0,0.0E0) , (0.78E0,0.06E0) , (-0.9E0,0.5E0) ,          &
     &      (0.06E0,-0.13E0) , (0.1E0,-0.5E0) , (-0.77E0,-0.49E0) ,     &
     &      (-0.5E0,-0.3E0) , (0.52E0,-1.51E0)/
      DATA ((ct8(i,j,3),i=1,7),j=1,4)/(0.6E0,-0.6E0) , (0.0E0,0.0E0) ,  &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0) ,             &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.32E0,-1.41E0) ,          &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0) ,             &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0) ,             &
     &      (-0.07E0,-0.89E0) , (-1.18E0,-0.31E0) , (0.0E0,0.0E0) ,     &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0) ,             &
     &      (0.0E0,0.0E0) , (0.78E0,0.06E0) , (-1.54E0,0.97E0) ,        &
     &      (0.03E0,-0.89E0) , (-0.18E0,-1.31E0) , (0.0E0,0.0E0) ,      &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0)/
      DATA ((ct8(i,j,4),i=1,7),j=1,4)/(0.6E0,-0.6E0) , (0.0E0,0.0E0) ,  &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0) ,             &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.32E0,-1.41E0) ,          &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0) ,             &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0) ,             &
     &      (0.32E0,-1.41E0) , (-0.9E0,0.5E0) , (0.05E0,-0.6E0) ,       &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0) ,             &
     &      (0.0E0,0.0E0) , (0.32E0,-1.41E0) , (-0.9E0,0.5E0) ,         &
     &      (0.05E0,-0.6E0) , (0.1E0,-0.5E0) , (-0.77E0,-0.49E0) ,      &
     &      (-0.5E0,-0.3E0) , (0.32E0,-1.16E0)/
      DATA ct7/(0.0E0,0.0E0) , (-0.06E0,-0.90E0) , (0.65E0,-0.47E0) ,   &
     &     (-0.34E0,-1.22E0) , (0.0E0,0.0E0) , (-0.06E0,-0.90E0) ,      &
     &     (-0.59E0,-1.46E0) , (-1.04E0,-0.04E0) , (0.0E0,0.0E0) ,      &
     &     (-0.06E0,-0.90E0) , (-0.83E0,0.59E0) , (0.07E0,-0.37E0) ,    &
     &     (0.0E0,0.0E0) , (-0.06E0,-0.90E0) , (-0.76E0,-1.15E0) ,      &
     &     (-1.33E0,-1.82E0)/
      DATA ct6/(0.0E0,0.0E0) , (0.90E0,0.06E0) , (0.91E0,-0.77E0) ,     &
     &     (1.80E0,-0.10E0) , (0.0E0,0.0E0) , (0.90E0,0.06E0) ,         &
     &     (1.45E0,0.74E0) , (0.20E0,0.90E0) , (0.0E0,0.0E0) ,          &
     &     (0.90E0,0.06E0) , (-0.55E0,0.23E0) , (0.83E0,-0.39E0) ,      &
     &     (0.0E0,0.0E0) , (0.90E0,0.06E0) , (1.04E0,0.79E0) ,          &
     &     (1.95E0,1.22E0)/
      DATA ((ct10x(i,j,1),i=1,7),j=1,4)/(0.7E0,-0.8E0) , (0.0E0,0.0E0) ,&
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0) ,             &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.6E0,-0.6E0) ,            &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0) ,             &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0) ,             &
     &      (0.6E0,-0.6E0) , (-0.9E0,0.5E0) , (0.0E0,0.0E0) ,           &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0) ,             &
     &      (0.0E0,0.0E0) , (0.6E0,-0.6E0) , (-0.9E0,0.5E0) ,           &
     &      (0.7E0,-0.6E0) , (0.1E0,-0.5E0) , (0.0E0,0.0E0) ,           &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0)/
      DATA ((ct10x(i,j,2),i=1,7),j=1,4)/(0.7E0,-0.8E0) , (0.0E0,0.0E0) ,&
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0) ,             &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.6E0,-0.6E0) ,            &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0) ,             &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0) ,             &
     &      (0.7E0,-0.6E0) , (-0.4E0,-0.7E0) , (0.6E0,-0.6E0) ,         &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0) ,             &
     &      (0.0E0,0.0E0) , (0.8E0,-0.7E0) , (-0.4E0,-0.7E0) ,          &
     &      (-0.1E0,-0.2E0) , (0.2E0,-0.8E0) , (0.7E0,-0.6E0) ,         &
     &      (0.1E0,0.4E0) , (0.6E0,-0.6E0)/
      DATA ((ct10x(i,j,3),i=1,7),j=1,4)/(0.7E0,-0.8E0) , (0.0E0,0.0E0) ,&
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0) ,             &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.6E0,-0.6E0) ,            &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0) ,             &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0) ,             &
     &      (-0.9E0,0.5E0) , (-0.4E0,-0.7E0) , (0.6E0,-0.6E0) ,         &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0) ,             &
     &      (0.0E0,0.0E0) , (0.1E0,-0.5E0) , (-0.4E0,-0.7E0) ,          &
     &      (0.7E0,-0.6E0) , (0.2E0,-0.8E0) , (-0.9E0,0.5E0) ,          &
     &      (0.1E0,0.4E0) , (0.6E0,-0.6E0)/
      DATA ((ct10x(i,j,4),i=1,7),j=1,4)/(0.7E0,-0.8E0) , (0.0E0,0.0E0) ,&
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0) ,             &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.6E0,-0.6E0) ,            &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0) ,             &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0) ,             &
     &      (0.6E0,-0.6E0) , (0.7E0,-0.6E0) , (0.0E0,0.0E0) ,           &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0) ,             &
     &      (0.0E0,0.0E0) , (0.6E0,-0.6E0) , (0.7E0,-0.6E0) ,           &
     &      (-0.1E0,-0.2E0) , (0.8E0,-0.7E0) , (0.0E0,0.0E0) ,          &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0)/
      DATA ((ct10y(i,j,1),i=1,7),j=1,4)/(0.6E0,-0.6E0) , (0.0E0,0.0E0) ,&
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0) ,             &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.7E0,-0.8E0) ,            &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0) ,             &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0) ,             &
     &      (0.7E0,-0.8E0) , (-0.4E0,-0.7E0) , (0.0E0,0.0E0) ,          &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0) ,             &
     &      (0.0E0,0.0E0) , (0.7E0,-0.8E0) , (-0.4E0,-0.7E0) ,          &
     &      (-0.1E0,-0.9E0) , (0.2E0,-0.8E0) , (0.0E0,0.0E0) ,          &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0)/
      DATA ((ct10y(i,j,2),i=1,7),j=1,4)/(0.6E0,-0.6E0) , (0.0E0,0.0E0) ,&
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0) ,             &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.7E0,-0.8E0) ,            &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0) ,             &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0) ,             &
     &      (-0.1E0,-0.9E0) , (-0.9E0,0.5E0) , (0.7E0,-0.8E0) ,         &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0) ,             &
     &      (0.0E0,0.0E0) , (-0.6E0,0.6E0) , (-0.9E0,0.5E0) ,           &
     &      (-0.9E0,-0.4E0) , (0.1E0,-0.5E0) , (-0.1E0,-0.9E0) ,        &
     &      (-0.5E0,-0.3E0) , (0.7E0,-0.8E0)/
      DATA ((ct10y(i,j,3),i=1,7),j=1,4)/(0.6E0,-0.6E0) , (0.0E0,0.0E0) ,&
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0) ,             &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.7E0,-0.8E0) ,            &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0) ,             &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0) ,             &
     &      (-0.1E0,-0.9E0) , (0.7E0,-0.8E0) , (0.0E0,0.0E0) ,          &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0) ,             &
     &      (0.0E0,0.0E0) , (-0.6E0,0.6E0) , (-0.9E0,-0.4E0) ,          &
     &      (-0.1E0,-0.9E0) , (0.7E0,-0.8E0) , (0.0E0,0.0E0) ,          &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0)/
      DATA ((ct10y(i,j,4),i=1,7),j=1,4)/(0.6E0,-0.6E0) , (0.0E0,0.0E0) ,&
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0) ,             &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.7E0,-0.8E0) ,            &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0) ,             &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0) ,             &
     &      (0.7E0,-0.8E0) , (-0.9E0,0.5E0) , (-0.4E0,-0.7E0) ,         &
     &      (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0) ,             &
     &      (0.0E0,0.0E0) , (0.7E0,-0.8E0) , (-0.9E0,0.5E0) ,           &
     &      (-0.4E0,-0.7E0) , (0.1E0,-0.5E0) , (-0.1E0,-0.9E0) ,        &
     &      (-0.5E0,-0.3E0) , (0.2E0,-0.8E0)/
      DATA csize1/(0.0E0,0.0E0) , (0.9E0,0.9E0) , (1.63E0,1.73E0) ,     &
     &     (2.90E0,2.78E0)/
      DATA csize3/(0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0) ,       &
     &     (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0)&
     &     , (1.17E0,1.17E0) , (1.17E0,1.17E0) , (1.17E0,1.17E0) ,      &
     &     (1.17E0,1.17E0) , (1.17E0,1.17E0) , (1.17E0,1.17E0) ,        &
     &     (1.17E0,1.17E0)/
      DATA csize2/(0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0) ,       &
     &     (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0) , (0.0E0,0.0E0)&
     &     , (1.54E0,1.54E0) , (1.54E0,1.54E0) , (1.54E0,1.54E0) ,      &
     &     (1.54E0,1.54E0) , (1.54E0,1.54E0) , (1.54E0,1.54E0) ,        &
     &     (1.54E0,1.54E0)/
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
!              .. CDOTC ..
               cdot(1) = CDOTC(N,cx,INCx,cy,INCy)
               CALL CTEST(1,cdot,ct6(kn,ki),csize1(kn),Sfac)
            ELSEIF ( ICAse==2 ) THEN
!              .. CDOTU ..
               cdot(1) = CDOTU(N,cx,INCx,cy,INCy)
               CALL CTEST(1,cdot,ct7(kn,ki),csize1(kn),Sfac)
            ELSEIF ( ICAse==3 ) THEN
!              .. CAXPY ..
               CALL CAXPY(N,ca,cx,INCx,cy,INCy)
               CALL CTEST(leny,cy,ct8(1,kn,ki),csize2(1,ksize),Sfac)
            ELSEIF ( ICAse==4 ) THEN
!              .. CCOPY ..
               CALL CCOPY(N,cx,INCx,cy,INCy)
               CALL CTEST(leny,cy,ct10y(1,kn,ki),csize3,1.0E0)
               IF ( ki==1 ) THEN
                  cx0(1) = (42.0E0,43.0E0)
                  cy0(1) = (44.0E0,45.0E0)
                  IF ( N==0 ) THEN
                     cty0(1) = cy0(1)
                  ELSE
                     cty0(1) = cx0(1)
                  ENDIF
                  lincx = INCx
                  INCx = 0
                  lincy = INCy
                  INCy = 0
                  CALL CCOPY(N,cx0,INCx,cy0,INCy)
                  CALL CTEST(1,cy0,cty0,csize3,1.0E0)
                  INCx = lincx
                  INCy = lincy
               ENDIF
            ELSEIF ( ICAse==5 ) THEN
!              .. CSWAP ..
               CALL CSWAP(N,cx,INCx,cy,INCy)
               CALL CTEST(lenx,cx,ct10x(1,kn,ki),csize3,1.0E0)
               CALL CTEST(leny,cy,ct10y(1,kn,ki),csize3,1.0E0)
            ELSE
               WRITE (NOUT,*) ' Shouldn''t be here in CHECK2'
               STOP
            ENDIF
!
         ENDDO
      ENDDO
      END SUBROUTINE CHECK2
!*==stest.f90  processed by SPAG 7.51RB at 17:52 on  4 Mar 2022
      SUBROUTINE STEST(Len,Scomp,Strue,Ssize,Sfac)
      IMPLICIT NONE
!*--STEST591
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
      INTEGER ICAse , INCx , INCy , MODe , N
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
99003       FORMAT (1X,I4,I3,3I5,I3,2E36.8,2E12.4)
         ENDIF
      ENDDO
      RETURN
      END SUBROUTINE STEST
!*==stest1.f90  processed by SPAG 7.51RB at 17:52 on  4 Mar 2022
      SUBROUTINE STEST1(Scomp1,Strue1,Ssize,Sfac)
      IMPLICIT NONE
!*--STEST1653
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
!*==sdiff.f90  processed by SPAG 7.51RB at 17:52 on  4 Mar 2022
      REAL FUNCTION SDIFF(Sa,Sb)
      IMPLICIT NONE
!*--SDIFF680
!     ********************************* SDIFF **************************
!     COMPUTES DIFFERENCE OF TWO NUMBERS.  C. L. LAWSON, JPL 1974 FEB 15
!
!     .. Scalar Arguments ..
      REAL Sa , Sb
!     .. Executable Statements ..
      SDIFF = Sa - Sb
      END FUNCTION SDIFF
!*==ctest.f90  processed by SPAG 7.51RB at 17:52 on  4 Mar 2022
      SUBROUTINE CTEST(Len,Ccomp,Ctrue,Csize,Sfac)
      IMPLICIT NONE
!*--CTEST692
!     **************************** CTEST *****************************
!
!     C.L. LAWSON, JPL, 1978 DEC 6
!
!     .. Scalar Arguments ..
      REAL Sfac
      INTEGER Len
!     .. Array Arguments ..
      COMPLEX Ccomp(Len) , Csize(Len) , Ctrue(Len)
!     .. Local Scalars ..
      INTEGER i
!     .. Local Arrays ..
      REAL scomp(20) , ssize(20) , strue(20)
!     .. External Subroutines ..
      EXTERNAL STEST
!     .. Intrinsic Functions ..
      INTRINSIC AIMAG , REAL
!     .. Executable Statements ..
      DO i = 1 , Len
         scomp(2*i-1) = REAL(Ccomp(i))
         scomp(2*i) = AIMAG(Ccomp(i))
         strue(2*i-1) = REAL(Ctrue(i))
         strue(2*i) = AIMAG(Ctrue(i))
         ssize(2*i-1) = REAL(Csize(i))
         ssize(2*i) = AIMAG(Csize(i))
      ENDDO
!
      CALL STEST(2*Len,scomp,strue,ssize,Sfac)
      END SUBROUTINE CTEST
!*==itest1.f90  processed by SPAG 7.51RB at 17:52 on  4 Mar 2022
      SUBROUTINE ITEST1(Icomp,Itrue)
      IMPLICIT NONE
!*--ITEST1725
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
