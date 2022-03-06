!*==dchkeq.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b dchkeq
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE dchkeq( THRESH, NOUT )
!
!       .. Scalar Arguments ..
!       INTEGER            NOUT
!       DOUBLE PRECISION   THRESH
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> dchkeq tests DGEEQU, DGBEQU, DPOEQU, DPPEQU and DPBEQU
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] THRESH
!> \verbatim
!>          THRESH is DOUBLE PRECISION
!>          Threshold for testing routines. Should be between 2 and 10.
!> \endverbatim
!>
!> \param[in] NOUT
!> \verbatim
!>          NOUT is INTEGER
!>          The unit number for output.
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
!> \ingroup double_lin
!
!  =====================================================================
      SUBROUTINE dchkeq(Thresh,Nout)
      IMPLICIT NONE
!*--DCHKEQ58
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Nout
      DOUBLE PRECISION Thresh
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE , TEN
      PARAMETER (ZERO=0.0D0,ONE=1.0D+0,TEN=1.0D1)
      INTEGER NSZ , NSZB
      PARAMETER (NSZ=5,NSZB=3*NSZ-2)
      INTEGER NSZP , NPOW
      PARAMETER (NSZP=(NSZ*(NSZ+1))/2,NPOW=2*NSZ+1)
!     ..
!     .. Local Scalars ..
      LOGICAL ok
      CHARACTER*3 path
      INTEGER i , info , j , kl , ku , m , n
      DOUBLE PRECISION ccond , eps , norm , ratio , rcmax , rcmin ,     &
     &                 rcond
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION a(NSZ,NSZ) , ab(NSZB,NSZ) , ap(NSZP) , c(NSZ) ,  &
     &                 pow(NPOW) , r(NSZ) , reslts(5) , rpow(NPOW)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH
EXTERNAL    DLAMCH
!     ..
!     .. External Subroutines ..
EXTERNAL    DGBEQU,   DGEEQU,   DPBEQU,   DPOEQU,   DPPEQU
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , MIN
!     ..
!     .. Executable Statements ..
!
      path(1:1) = 'Double precision'
      path(2:3) = 'EQ'
!
      eps = DLAMCH('P')
      DO i = 1 , 5
         reslts(i) = ZERO
      ENDDO
      DO i = 1 , NPOW
         pow(i) = TEN**(i-1)
         rpow(i) = ONE/pow(i)
      ENDDO
!
!     Test DGEEQU
!
      DO n = 0 , NSZ
         DO m = 0 , NSZ
!
            DO j = 1 , NSZ
               DO i = 1 , NSZ
                  IF ( i<=m .AND. j<=n ) THEN
                     a(i,j) = pow(i+j+1)*(-1)**(i+j)
                  ELSE
                     a(i,j) = ZERO
                  ENDIF
               ENDDO
            ENDDO
!
            CALL DGEEQU(m,n,a,NSZ,r,c,rcond,ccond,norm,info)
!
            IF ( info/=0 ) THEN
               reslts(1) = ONE
            ELSEIF ( n/=0 .AND. m/=0 ) THEN
               reslts(1) = MAX(reslts(1),ABS((rcond-rpow(m))/rpow(m)))
               reslts(1) = MAX(reslts(1),ABS((ccond-rpow(n))/rpow(n)))
               reslts(1) = MAX(reslts(1),ABS((norm-pow(n+m+1))/pow(n+m+1&
     &                     )))
               DO i = 1 , m
                  reslts(1) = MAX(reslts(1),ABS((r(i)-rpow(i+n+1))/rpow(&
     &                        i+n+1)))
               ENDDO
               DO j = 1 , n
                  reslts(1) = MAX(reslts(1),ABS((c(j)-pow(n-j+1))/pow(n-&
     &                        j+1)))
               ENDDO
            ENDIF
!
         ENDDO
      ENDDO
!
!     Test with zero rows and columns
!
      DO j = 1 , NSZ
         a(MAX(NSZ-1,1),j) = ZERO
      ENDDO
      CALL DGEEQU(NSZ,NSZ,a,NSZ,r,c,rcond,ccond,norm,info)
      IF ( info/=MAX(NSZ-1,1) ) reslts(1) = ONE
!
      DO j = 1 , NSZ
         a(MAX(NSZ-1,1),j) = ONE
      ENDDO
      DO i = 1 , NSZ
         a(i,MAX(NSZ-1,1)) = ZERO
      ENDDO
      CALL DGEEQU(NSZ,NSZ,a,NSZ,r,c,rcond,ccond,norm,info)
      IF ( info/=NSZ+MAX(NSZ-1,1) ) reslts(1) = ONE
      reslts(1) = reslts(1)/eps
!
!     Test DGBEQU
!
      DO n = 0 , NSZ
         DO m = 0 , NSZ
            DO kl = 0 , MAX(m-1,0)
               DO ku = 0 , MAX(n-1,0)
!
                  DO j = 1 , NSZ
                     DO i = 1 , NSZB
                        ab(i,j) = ZERO
                     ENDDO
                  ENDDO
                  DO j = 1 , n
                     DO i = 1 , m
                        IF ( i<=MIN(m,j+kl) .AND. i>=MAX(1,j-ku) .AND.  &
     &                       j<=n ) ab(ku+1+i-j,j) = pow(i+j+1)*(-1)    &
     &                       **(i+j)
                     ENDDO
                  ENDDO
!
                  CALL DGBEQU(m,n,kl,ku,ab,NSZB,r,c,rcond,ccond,norm,   &
     &                        info)
!
                  IF ( info/=0 ) THEN
                     IF ( .NOT.((n+kl<m .AND. info==n+kl+1) .OR.        &
     &                    (m+ku<n .AND. info==2*m+ku+1)) ) reslts(2)    &
     &                    = ONE
                  ELSEIF ( n/=0 .AND. m/=0 ) THEN
!
                     rcmin = r(1)
                     rcmax = r(1)
                     DO i = 1 , m
                        rcmin = MIN(rcmin,r(i))
                        rcmax = MAX(rcmax,r(i))
                     ENDDO
                     ratio = rcmin/rcmax
                     reslts(2) = MAX(reslts(2),ABS((rcond-ratio)/ratio))
!
                     rcmin = c(1)
                     rcmax = c(1)
                     DO j = 1 , n
                        rcmin = MIN(rcmin,c(j))
                        rcmax = MAX(rcmax,c(j))
                     ENDDO
                     ratio = rcmin/rcmax
                     reslts(2) = MAX(reslts(2),ABS((ccond-ratio)/ratio))
!
                     reslts(2) = MAX(reslts(2),ABS((norm-pow(n+m+1))/pow&
     &                           (n+m+1)))
                     DO i = 1 , m
                        rcmax = ZERO
                        DO j = 1 , n
                           IF ( i<=j+kl .AND. i>=j-ku ) THEN
                              ratio = ABS(r(i)*pow(i+j+1)*c(j))
                              rcmax = MAX(rcmax,ratio)
                           ENDIF
                        ENDDO
                        reslts(2) = MAX(reslts(2),ABS(ONE-rcmax))
                     ENDDO
!
                     DO j = 1 , n
                        rcmax = ZERO
                        DO i = 1 , m
                           IF ( i<=j+kl .AND. i>=j-ku ) THEN
                              ratio = ABS(r(i)*pow(i+j+1)*c(j))
                              rcmax = MAX(rcmax,ratio)
                           ENDIF
                        ENDDO
                        reslts(2) = MAX(reslts(2),ABS(ONE-rcmax))
                     ENDDO
                  ENDIF
!
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      reslts(2) = reslts(2)/eps
!
!     Test DPOEQU
!
      DO n = 0 , NSZ
!
         DO i = 1 , NSZ
            DO j = 1 , NSZ
               IF ( i<=n .AND. j==i ) THEN
                  a(i,j) = pow(i+j+1)*(-1)**(i+j)
               ELSE
                  a(i,j) = ZERO
               ENDIF
            ENDDO
         ENDDO
!
         CALL DPOEQU(n,a,NSZ,r,rcond,norm,info)
!
         IF ( info/=0 ) THEN
            reslts(3) = ONE
         ELSEIF ( n/=0 ) THEN
            reslts(3) = MAX(reslts(3),ABS((rcond-rpow(n))/rpow(n)))
            reslts(3) = MAX(reslts(3),ABS((norm-pow(2*n+1))/pow(2*n+1)))
            DO i = 1 , n
               reslts(3) = MAX(reslts(3),ABS((r(i)-rpow(i+1))/rpow(i+1))&
     &                     )
            ENDDO
         ENDIF
      ENDDO
      a(MAX(NSZ-1,1),MAX(NSZ-1,1)) = -ONE
      CALL DPOEQU(NSZ,a,NSZ,r,rcond,norm,info)
      IF ( info/=MAX(NSZ-1,1) ) reslts(3) = ONE
      reslts(3) = reslts(3)/eps
!
!     Test DPPEQU
!
      DO n = 0 , NSZ
!
!        Upper triangular packed storage
!
         DO i = 1 , (n*(n+1))/2
            ap(i) = ZERO
         ENDDO
         DO i = 1 , n
            ap((i*(i+1))/2) = pow(2*i+1)
         ENDDO
!
         CALL DPPEQU('U',n,ap,r,rcond,norm,info)
!
         IF ( info/=0 ) THEN
            reslts(4) = ONE
         ELSEIF ( n/=0 ) THEN
            reslts(4) = MAX(reslts(4),ABS((rcond-rpow(n))/rpow(n)))
            reslts(4) = MAX(reslts(4),ABS((norm-pow(2*n+1))/pow(2*n+1)))
            DO i = 1 , n
               reslts(4) = MAX(reslts(4),ABS((r(i)-rpow(i+1))/rpow(i+1))&
     &                     )
            ENDDO
         ENDIF
!
!        Lower triangular packed storage
!
         DO i = 1 , (n*(n+1))/2
            ap(i) = ZERO
         ENDDO
         j = 1
         DO i = 1 , n
            ap(j) = pow(2*i+1)
            j = j + (n-i+1)
         ENDDO
!
         CALL DPPEQU('L',n,ap,r,rcond,norm,info)
!
         IF ( info/=0 ) THEN
            reslts(4) = ONE
         ELSEIF ( n/=0 ) THEN
            reslts(4) = MAX(reslts(4),ABS((rcond-rpow(n))/rpow(n)))
            reslts(4) = MAX(reslts(4),ABS((norm-pow(2*n+1))/pow(2*n+1)))
            DO i = 1 , n
               reslts(4) = MAX(reslts(4),ABS((r(i)-rpow(i+1))/rpow(i+1))&
     &                     )
            ENDDO
         ENDIF
!
      ENDDO
      i = (NSZ*(NSZ+1))/2 - 2
      ap(i) = -ONE
      CALL DPPEQU('L',NSZ,ap,r,rcond,norm,info)
      IF ( info/=MAX(NSZ-1,1) ) reslts(4) = ONE
      reslts(4) = reslts(4)/eps
!
!     Test DPBEQU
!
      DO n = 0 , NSZ
         DO kl = 0 , MAX(n-1,0)
!
!           Test upper triangular storage
!
            DO j = 1 , NSZ
               DO i = 1 , NSZB
                  ab(i,j) = ZERO
               ENDDO
            ENDDO
            DO j = 1 , n
               ab(kl+1,j) = pow(2*j+1)
            ENDDO
!
            CALL DPBEQU('U',n,kl,ab,NSZB,r,rcond,norm,info)
!
            IF ( info/=0 ) THEN
               reslts(5) = ONE
            ELSEIF ( n/=0 ) THEN
               reslts(5) = MAX(reslts(5),ABS((rcond-rpow(n))/rpow(n)))
               reslts(5) = MAX(reslts(5),ABS((norm-pow(2*n+1))/pow(2*n+1&
     &                     )))
               DO i = 1 , n
                  reslts(5) = MAX(reslts(5),ABS((r(i)-rpow(i+1))/rpow(i+&
     &                        1)))
               ENDDO
            ENDIF
            IF ( n/=0 ) THEN
               ab(kl+1,MAX(n-1,1)) = -ONE
               CALL DPBEQU('U',n,kl,ab,NSZB,r,rcond,norm,info)
               IF ( info/=MAX(n-1,1) ) reslts(5) = ONE
            ENDIF
!
!           Test lower triangular storage
!
            DO j = 1 , NSZ
               DO i = 1 , NSZB
                  ab(i,j) = ZERO
               ENDDO
            ENDDO
            DO j = 1 , n
               ab(1,j) = pow(2*j+1)
            ENDDO
!
            CALL DPBEQU('L',n,kl,ab,NSZB,r,rcond,norm,info)
!
            IF ( info/=0 ) THEN
               reslts(5) = ONE
            ELSEIF ( n/=0 ) THEN
               reslts(5) = MAX(reslts(5),ABS((rcond-rpow(n))/rpow(n)))
               reslts(5) = MAX(reslts(5),ABS((norm-pow(2*n+1))/pow(2*n+1&
     &                     )))
               DO i = 1 , n
                  reslts(5) = MAX(reslts(5),ABS((r(i)-rpow(i+1))/rpow(i+&
     &                        1)))
               ENDDO
            ENDIF
            IF ( n/=0 ) THEN
               ab(1,MAX(n-1,1)) = -ONE
               CALL DPBEQU('L',n,kl,ab,NSZB,r,rcond,norm,info)
               IF ( info/=MAX(n-1,1) ) reslts(5) = ONE
            ENDIF
         ENDDO
      ENDDO
      reslts(5) = reslts(5)/eps
      ok = (reslts(1)<=Thresh) .AND. (reslts(2)<=Thresh) .AND.          &
     &     (reslts(3)<=Thresh) .AND. (reslts(4)<=Thresh) .AND.          &
     &     (reslts(5)<=Thresh)
      WRITE (Nout,FMT=*)
      IF ( ok ) THEN
         WRITE (Nout,FMT=99001) path
      ELSE
         IF ( reslts(1)>Thresh ) WRITE (Nout,FMT=99002) reslts(1) ,     &
     &                                  Thresh
         IF ( reslts(2)>Thresh ) WRITE (Nout,FMT=99003) reslts(2) ,     &
     &                                  Thresh
         IF ( reslts(3)>Thresh ) WRITE (Nout,FMT=99004) reslts(3) ,     &
     &                                  Thresh
         IF ( reslts(4)>Thresh ) WRITE (Nout,FMT=99005) reslts(4) ,     &
     &                                  Thresh
         IF ( reslts(5)>Thresh ) WRITE (Nout,FMT=99006) reslts(5) ,     &
     &                                  Thresh
      ENDIF
99001 FORMAT (1X,'All tests for ',A3,' routines passed the threshold')
99002 FORMAT (' DGEEQU failed test with value ',D10.3,' exceeding',     &
     &        ' threshold ',D10.3)
99003 FORMAT (' DGBEQU failed test with value ',D10.3,' exceeding',     &
     &        ' threshold ',D10.3)
99004 FORMAT (' DPOEQU failed test with value ',D10.3,' exceeding',     &
     &        ' threshold ',D10.3)
99005 FORMAT (' DPPEQU failed test with value ',D10.3,' exceeding',     &
     &        ' threshold ',D10.3)
99006 FORMAT (' DPBEQU failed test with value ',D10.3,' exceeding',     &
     &        ' threshold ',D10.3)
!
!     End of dchkeq
!
      END SUBROUTINE dchkeq
