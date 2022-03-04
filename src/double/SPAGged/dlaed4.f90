!*==dlaed4.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DLAED4 used by sstedc. Finds a single root of the secular equation.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLAED4 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaed4.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaed4.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaed4.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAED4( N, I, D, Z, DELTA, RHO, DLAM, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            I, INFO, N
!       DOUBLE PRECISION   DLAM, RHO
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   D( * ), DELTA( * ), Z( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> This subroutine computes the I-th updated eigenvalue of a symmetric
!> rank-one modification to a diagonal matrix whose elements are
!> given in the array d, and that
!>
!>            D(i) < D(j)  for  i < j
!>
!> and that RHO > 0.  This is arranged by the calling routine, and is
!> no loss in generality.  The rank-one modified system is thus
!>
!>            diag( D )  +  RHO * Z * Z_transpose.
!>
!> where we assume the Euclidean norm of Z is 1.
!>
!> The method consists of approximating the rational functions in the
!> secular equation by simpler interpolating rational functions.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>         The length of all arrays.
!> \endverbatim
!>
!> \param[in] I
!> \verbatim
!>          I is INTEGER
!>         The index of the eigenvalue to be computed.  1 <= I <= N.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>         The original eigenvalues.  It is assumed that they are in
!>         order, D(I) < D(J)  for I < J.
!> \endverbatim
!>
!> \param[in] Z
!> \verbatim
!>          Z is DOUBLE PRECISION array, dimension (N)
!>         The components of the updating vector.
!> \endverbatim
!>
!> \param[out] DELTA
!> \verbatim
!>          DELTA is DOUBLE PRECISION array, dimension (N)
!>         If N > 2, DELTA contains (D(j) - lambda_I) in its  j-th
!>         component.  If N = 1, then DELTA(1) = 1. If N = 2, see DLAED5
!>         for detail. The vector DELTA contains the information necessary
!>         to construct the eigenvectors by DLAED3 and DLAED9.
!> \endverbatim
!>
!> \param[in] RHO
!> \verbatim
!>          RHO is DOUBLE PRECISION
!>         The scalar in the symmetric updating formula.
!> \endverbatim
!>
!> \param[out] DLAM
!> \verbatim
!>          DLAM is DOUBLE PRECISION
!>         The computed lambda_I, the I-th updated eigenvalue.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>         = 0:  successful exit
!>         > 0:  if INFO = 1, the updating process failed.
!> \endverbatim
!
!> \par Internal Parameters:
!  =========================
!>
!> \verbatim
!>  Logical variable ORGATI (origin-at-i?) is used for distinguishing
!>  whether D(i) or D(i+1) is treated as the origin.
!>
!>            ORGATI = .true.    origin at i
!>            ORGATI = .false.   origin at i+1
!>
!>   Logical variable SWTCH3 (switch-for-3-poles?) is for noting
!>   if we are working with THREE poles!
!>
!>   MAXIT is the maximum number of iterations allowed for each
!>   eigenvalue.
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
!> \ingroup auxOTHERcomputational
!
!> \par Contributors:
!  ==================
!>
!>     Ren-Cang Li, Computer Science Division, University of California
!>     at Berkeley, USA
!>
!  =====================================================================
      SUBROUTINE DLAED4(N,I,D,Z,Delta,Rho,Dlam,Info)
      USE F77KINDS                        
      USE S_DLAED5
      USE S_DLAED6
      USE S_DLAMCH
      IMPLICIT NONE
!*--DLAED4153
!
! PARAMETER definitions rewritten by SPAG
!
      INTEGER , PARAMETER  ::  MAXIT = 30
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0 ,        &
     &                              TWO = 2.0D0 , THREE = 3.0D0 ,       &
     &                              FOUR = 4.0D0 , EIGHT = 8.0D0 ,      &
     &                              TEN = 10.0D0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: N
      INTEGER :: I
      REAL(R8KIND) , DIMENSION(*) :: D
      REAL(R8KIND) , DIMENSION(*) :: Z
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Delta
      REAL(R8KIND) :: Rho
      REAL(R8KIND) :: Dlam
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: a , b , c , del , dltlb , dltub , dphi , dpsi ,   &
     &                dw , eps , erretm , eta , midpt , phi , prew ,    &
     &                psi , rhoinv , tau , temp , temp1 , w
      INTEGER :: ii , iim1 , iip1 , ip1 , iter , j , niter
      LOGICAL :: orgati , swtch , swtch3
      REAL(R8KIND) , DIMENSION(3) :: zz
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. Local Arrays ..
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Since this routine is called in an inner loop, we do no argument
!     checking.
!
!     Quick return for N=1 and 2.
!
      Info = 0
      IF ( N==1 ) THEN
!
!         Presumably, I=1 upon entry
!
         Dlam = D(1) + Rho*Z(1)*Z(1)
         Delta(1) = ONE
         RETURN
      ENDIF
      IF ( N==2 ) THEN
         CALL DLAED5(I,D,Z,Delta,Rho,Dlam)
         RETURN
      ENDIF
!
!     Compute machine epsilon
!
      eps = DLAMCH('Epsilon')
      rhoinv = ONE/Rho
!
!     The case I = N
!
      IF ( I==N ) THEN
!
!        Initialize some basic variables
!
         ii = N - 1
         niter = 1
!
!        Calculate initial guess
!
         midpt = Rho/TWO
!
!        If ||Z||_2 is not one, then TEMP should be set to
!        RHO * ||Z||_2^2 / TWO
!
         DO j = 1 , N
            Delta(j) = (D(j)-D(I)) - midpt
         ENDDO
!
         psi = ZERO
         DO j = 1 , N - 2
            psi = psi + Z(j)*Z(j)/Delta(j)
         ENDDO
!
         c = rhoinv + psi
         w = c + Z(ii)*Z(ii)/Delta(ii) + Z(N)*Z(N)/Delta(N)
!
         IF ( w<=ZERO ) THEN
            temp = Z(N-1)*Z(N-1)/(D(N)-D(N-1)+Rho) + Z(N)*Z(N)/Rho
            IF ( c<=temp ) THEN
               tau = Rho
            ELSE
               del = D(N) - D(N-1)
               a = -c*del + Z(N-1)*Z(N-1) + Z(N)*Z(N)
               b = Z(N)*Z(N)*del
               IF ( a<ZERO ) THEN
                  tau = TWO*b/(SQRT(a*a+FOUR*b*c)-a)
               ELSE
                  tau = (a+SQRT(a*a+FOUR*b*c))/(TWO*c)
               ENDIF
            ENDIF
!
!           It can be proved that
!               D(N)+RHO/2 <= LAMBDA(N) < D(N)+TAU <= D(N)+RHO
!
            dltlb = midpt
            dltub = Rho
         ELSE
            del = D(N) - D(N-1)
            a = -c*del + Z(N-1)*Z(N-1) + Z(N)*Z(N)
            b = Z(N)*Z(N)*del
            IF ( a<ZERO ) THEN
               tau = TWO*b/(SQRT(a*a+FOUR*b*c)-a)
            ELSE
               tau = (a+SQRT(a*a+FOUR*b*c))/(TWO*c)
            ENDIF
!
!           It can be proved that
!               D(N) < D(N)+TAU < LAMBDA(N) < D(N)+RHO/2
!
            dltlb = ZERO
            dltub = midpt
         ENDIF
!
         DO j = 1 , N
            Delta(j) = (D(j)-D(I)) - tau
         ENDDO
!
!        Evaluate PSI and the derivative DPSI
!
         dpsi = ZERO
         psi = ZERO
         erretm = ZERO
         DO j = 1 , ii
            temp = Z(j)/Delta(j)
            psi = psi + Z(j)*temp
            dpsi = dpsi + temp*temp
            erretm = erretm + psi
         ENDDO
         erretm = ABS(erretm)
!
!        Evaluate PHI and the derivative DPHI
!
         temp = Z(N)/Delta(N)
         phi = Z(N)*temp
         dphi = temp*temp
         erretm = EIGHT*(-phi-psi) + erretm - phi + rhoinv + ABS(tau)   &
     &            *(dpsi+dphi)
!
         w = rhoinv + phi + psi
!
!        Test for convergence
!
         IF ( ABS(w)<=eps*erretm ) THEN
            Dlam = D(I) + tau
            GOTO 99999
         ENDIF
!
         IF ( w<=ZERO ) THEN
            dltlb = MAX(dltlb,tau)
         ELSE
            dltub = MIN(dltub,tau)
         ENDIF
!
!        Calculate the new step
!
         niter = niter + 1
         c = w - Delta(N-1)*dpsi - Delta(N)*dphi
         a = (Delta(N-1)+Delta(N))*w - Delta(N-1)*Delta(N)*(dpsi+dphi)
         b = Delta(N-1)*Delta(N)*w
         IF ( c<ZERO ) c = ABS(c)
         IF ( c==ZERO ) THEN
!          ETA = B/A
!           ETA = RHO - TAU
            eta = dltub - tau
         ELSEIF ( a>=ZERO ) THEN
            eta = (a+SQRT(ABS(a*a-FOUR*b*c)))/(TWO*c)
         ELSE
            eta = TWO*b/(a-SQRT(ABS(a*a-FOUR*b*c)))
         ENDIF
!
!        Note, eta should be positive if w is negative, and
!        eta should be negative otherwise. However,
!        if for some reason caused by roundoff, eta*w > 0,
!        we simply use one Newton step instead. This way
!        will guarantee eta*w < 0.
!
         IF ( w*eta>ZERO ) eta = -w/(dpsi+dphi)
         temp = tau + eta
         IF ( temp>dltub .OR. temp<dltlb ) THEN
            IF ( w<ZERO ) THEN
               eta = (dltub-tau)/TWO
            ELSE
               eta = (dltlb-tau)/TWO
            ENDIF
         ENDIF
         DO j = 1 , N
            Delta(j) = Delta(j) - eta
         ENDDO
!
         tau = tau + eta
!
!        Evaluate PSI and the derivative DPSI
!
         dpsi = ZERO
         psi = ZERO
         erretm = ZERO
         DO j = 1 , ii
            temp = Z(j)/Delta(j)
            psi = psi + Z(j)*temp
            dpsi = dpsi + temp*temp
            erretm = erretm + psi
         ENDDO
         erretm = ABS(erretm)
!
!        Evaluate PHI and the derivative DPHI
!
         temp = Z(N)/Delta(N)
         phi = Z(N)*temp
         dphi = temp*temp
         erretm = EIGHT*(-phi-psi) + erretm - phi + rhoinv + ABS(tau)   &
     &            *(dpsi+dphi)
!
         w = rhoinv + phi + psi
!
!        Main loop to update the values of the array   DELTA
!
         iter = niter + 1
!
         DO niter = iter , MAXIT
!
!           Test for convergence
!
            IF ( ABS(w)<=eps*erretm ) THEN
               Dlam = D(I) + tau
               GOTO 99999
            ENDIF
!
            IF ( w<=ZERO ) THEN
               dltlb = MAX(dltlb,tau)
            ELSE
               dltub = MIN(dltub,tau)
            ENDIF
!
!           Calculate the new step
!
            c = w - Delta(N-1)*dpsi - Delta(N)*dphi
            a = (Delta(N-1)+Delta(N))*w - Delta(N-1)*Delta(N)           &
     &          *(dpsi+dphi)
            b = Delta(N-1)*Delta(N)*w
            IF ( a>=ZERO ) THEN
               eta = (a+SQRT(ABS(a*a-FOUR*b*c)))/(TWO*c)
            ELSE
               eta = TWO*b/(a-SQRT(ABS(a*a-FOUR*b*c)))
            ENDIF
!
!           Note, eta should be positive if w is negative, and
!           eta should be negative otherwise. However,
!           if for some reason caused by roundoff, eta*w > 0,
!           we simply use one Newton step instead. This way
!           will guarantee eta*w < 0.
!
            IF ( w*eta>ZERO ) eta = -w/(dpsi+dphi)
            temp = tau + eta
            IF ( temp>dltub .OR. temp<dltlb ) THEN
               IF ( w<ZERO ) THEN
                  eta = (dltub-tau)/TWO
               ELSE
                  eta = (dltlb-tau)/TWO
               ENDIF
            ENDIF
            DO j = 1 , N
               Delta(j) = Delta(j) - eta
            ENDDO
!
            tau = tau + eta
!
!           Evaluate PSI and the derivative DPSI
!
            dpsi = ZERO
            psi = ZERO
            erretm = ZERO
            DO j = 1 , ii
               temp = Z(j)/Delta(j)
               psi = psi + Z(j)*temp
               dpsi = dpsi + temp*temp
               erretm = erretm + psi
            ENDDO
            erretm = ABS(erretm)
!
!           Evaluate PHI and the derivative DPHI
!
            temp = Z(N)/Delta(N)
            phi = Z(N)*temp
            dphi = temp*temp
            erretm = EIGHT*(-phi-psi) + erretm - phi + rhoinv + ABS(tau)&
     &               *(dpsi+dphi)
!
            w = rhoinv + phi + psi
         ENDDO
!
!        Return with INFO = 1, NITER = MAXIT and not converged
!
         Info = 1
         Dlam = D(I) + tau
!
!        End for the case I = N
!
      ELSE
!
!        The case for I < N
!
         niter = 1
         ip1 = I + 1
!
!        Calculate initial guess
!
         del = D(ip1) - D(I)
         midpt = del/TWO
         DO j = 1 , N
            Delta(j) = (D(j)-D(I)) - midpt
         ENDDO
!
         psi = ZERO
         DO j = 1 , I - 1
            psi = psi + Z(j)*Z(j)/Delta(j)
         ENDDO
!
         phi = ZERO
         DO j = N , I + 2 , -1
            phi = phi + Z(j)*Z(j)/Delta(j)
         ENDDO
         c = rhoinv + psi + phi
         w = c + Z(I)*Z(I)/Delta(I) + Z(ip1)*Z(ip1)/Delta(ip1)
!
         IF ( w>ZERO ) THEN
!
!           d(i)< the ith eigenvalue < (d(i)+d(i+1))/2
!
!           We choose d(i) as origin.
!
            orgati = .TRUE.
            a = c*del + Z(I)*Z(I) + Z(ip1)*Z(ip1)
            b = Z(I)*Z(I)*del
            IF ( a>ZERO ) THEN
               tau = TWO*b/(a+SQRT(ABS(a*a-FOUR*b*c)))
            ELSE
               tau = (a-SQRT(ABS(a*a-FOUR*b*c)))/(TWO*c)
            ENDIF
            dltlb = ZERO
            dltub = midpt
         ELSE
!
!           (d(i)+d(i+1))/2 <= the ith eigenvalue < d(i+1)
!
!           We choose d(i+1) as origin.
!
            orgati = .FALSE.
            a = c*del - Z(I)*Z(I) - Z(ip1)*Z(ip1)
            b = Z(ip1)*Z(ip1)*del
            IF ( a<ZERO ) THEN
               tau = TWO*b/(a-SQRT(ABS(a*a+FOUR*b*c)))
            ELSE
               tau = -(a+SQRT(ABS(a*a+FOUR*b*c)))/(TWO*c)
            ENDIF
            dltlb = -midpt
            dltub = ZERO
         ENDIF
!
         IF ( orgati ) THEN
            DO j = 1 , N
               Delta(j) = (D(j)-D(I)) - tau
            ENDDO
         ELSE
            DO j = 1 , N
               Delta(j) = (D(j)-D(ip1)) - tau
            ENDDO
         ENDIF
         IF ( orgati ) THEN
            ii = I
         ELSE
            ii = I + 1
         ENDIF
         iim1 = ii - 1
         iip1 = ii + 1
!
!        Evaluate PSI and the derivative DPSI
!
         dpsi = ZERO
         psi = ZERO
         erretm = ZERO
         DO j = 1 , iim1
            temp = Z(j)/Delta(j)
            psi = psi + Z(j)*temp
            dpsi = dpsi + temp*temp
            erretm = erretm + psi
         ENDDO
         erretm = ABS(erretm)
!
!        Evaluate PHI and the derivative DPHI
!
         dphi = ZERO
         phi = ZERO
         DO j = N , iip1 , -1
            temp = Z(j)/Delta(j)
            phi = phi + Z(j)*temp
            dphi = dphi + temp*temp
            erretm = erretm + phi
         ENDDO
!
         w = rhoinv + phi + psi
!
!        W is the value of the secular function with
!        its ii-th element removed.
!
         swtch3 = .FALSE.
         IF ( orgati ) THEN
            IF ( w<ZERO ) swtch3 = .TRUE.
         ELSE
            IF ( w>ZERO ) swtch3 = .TRUE.
         ENDIF
         IF ( ii==1 .OR. ii==N ) swtch3 = .FALSE.
!
         temp = Z(ii)/Delta(ii)
         dw = dpsi + dphi + temp*temp
         temp = Z(ii)*temp
         w = w + temp
         erretm = EIGHT*(phi-psi) + erretm + TWO*rhoinv +               &
     &            THREE*ABS(temp) + ABS(tau)*dw
!
!        Test for convergence
!
         IF ( ABS(w)<=eps*erretm ) THEN
            IF ( orgati ) THEN
               Dlam = D(I) + tau
            ELSE
               Dlam = D(ip1) + tau
            ENDIF
            GOTO 99999
         ENDIF
!
         IF ( w<=ZERO ) THEN
            dltlb = MAX(dltlb,tau)
         ELSE
            dltub = MIN(dltub,tau)
         ENDIF
!
!        Calculate the new step
!
         niter = niter + 1
         IF ( .NOT.swtch3 ) THEN
            IF ( orgati ) THEN
               c = w - Delta(ip1)*dw - (D(I)-D(ip1))*(Z(I)/Delta(I))**2
            ELSE
               c = w - Delta(I)*dw - (D(ip1)-D(I))*(Z(ip1)/Delta(ip1))  &
     &             **2
            ENDIF
            a = (Delta(I)+Delta(ip1))*w - Delta(I)*Delta(ip1)*dw
            b = Delta(I)*Delta(ip1)*w
            IF ( c==ZERO ) THEN
               IF ( a==ZERO ) THEN
                  IF ( orgati ) THEN
                     a = Z(I)*Z(I) + Delta(ip1)*Delta(ip1)*(dpsi+dphi)
                  ELSE
                     a = Z(ip1)*Z(ip1) + Delta(I)*Delta(I)*(dpsi+dphi)
                  ENDIF
               ENDIF
               eta = b/a
            ELSEIF ( a<=ZERO ) THEN
               eta = (a-SQRT(ABS(a*a-FOUR*b*c)))/(TWO*c)
            ELSE
               eta = TWO*b/(a+SQRT(ABS(a*a-FOUR*b*c)))
            ENDIF
         ELSE
!
!           Interpolation using THREE most relevant poles
!
            temp = rhoinv + psi + phi
            IF ( orgati ) THEN
               temp1 = Z(iim1)/Delta(iim1)
               temp1 = temp1*temp1
               c = temp - Delta(iip1)*(dpsi+dphi) - (D(iim1)-D(iip1))   &
     &             *temp1
               zz(1) = Z(iim1)*Z(iim1)
               zz(3) = Delta(iip1)*Delta(iip1)*((dpsi-temp1)+dphi)
            ELSE
               temp1 = Z(iip1)/Delta(iip1)
               temp1 = temp1*temp1
               c = temp - Delta(iim1)*(dpsi+dphi) - (D(iip1)-D(iim1))   &
     &             *temp1
               zz(1) = Delta(iim1)*Delta(iim1)*(dpsi+(dphi-temp1))
               zz(3) = Z(iip1)*Z(iip1)
            ENDIF
            zz(2) = Z(ii)*Z(ii)
            CALL DLAED6(niter,orgati,c,Delta(iim1),zz,w,eta,Info)
            IF ( Info/=0 ) GOTO 99999
         ENDIF
!
!        Note, eta should be positive if w is negative, and
!        eta should be negative otherwise. However,
!        if for some reason caused by roundoff, eta*w > 0,
!        we simply use one Newton step instead. This way
!        will guarantee eta*w < 0.
!
         IF ( w*eta>=ZERO ) eta = -w/dw
         temp = tau + eta
         IF ( temp>dltub .OR. temp<dltlb ) THEN
            IF ( w<ZERO ) THEN
               eta = (dltub-tau)/TWO
            ELSE
               eta = (dltlb-tau)/TWO
            ENDIF
         ENDIF
!
         prew = w
!
         DO j = 1 , N
            Delta(j) = Delta(j) - eta
         ENDDO
!
!        Evaluate PSI and the derivative DPSI
!
         dpsi = ZERO
         psi = ZERO
         erretm = ZERO
         DO j = 1 , iim1
            temp = Z(j)/Delta(j)
            psi = psi + Z(j)*temp
            dpsi = dpsi + temp*temp
            erretm = erretm + psi
         ENDDO
         erretm = ABS(erretm)
!
!        Evaluate PHI and the derivative DPHI
!
         dphi = ZERO
         phi = ZERO
         DO j = N , iip1 , -1
            temp = Z(j)/Delta(j)
            phi = phi + Z(j)*temp
            dphi = dphi + temp*temp
            erretm = erretm + phi
         ENDDO
!
         temp = Z(ii)/Delta(ii)
         dw = dpsi + dphi + temp*temp
         temp = Z(ii)*temp
         w = rhoinv + phi + psi + temp
         erretm = EIGHT*(phi-psi) + erretm + TWO*rhoinv +               &
     &            THREE*ABS(temp) + ABS(tau+eta)*dw
!
         swtch = .FALSE.
         IF ( orgati ) THEN
            IF ( -w>ABS(prew)/TEN ) swtch = .TRUE.
         ELSE
            IF ( w>ABS(prew)/TEN ) swtch = .TRUE.
         ENDIF
!
         tau = tau + eta
!
!        Main loop to update the values of the array   DELTA
!
         iter = niter + 1
!
         DO niter = iter , MAXIT
!
!           Test for convergence
!
            IF ( ABS(w)<=eps*erretm ) THEN
               IF ( orgati ) THEN
                  Dlam = D(I) + tau
               ELSE
                  Dlam = D(ip1) + tau
               ENDIF
               GOTO 99999
            ENDIF
!
            IF ( w<=ZERO ) THEN
               dltlb = MAX(dltlb,tau)
            ELSE
               dltub = MIN(dltub,tau)
            ENDIF
!
!           Calculate the new step
!
            IF ( .NOT.swtch3 ) THEN
               IF ( swtch ) THEN
                  temp = Z(ii)/Delta(ii)
                  IF ( orgati ) THEN
                     dpsi = dpsi + temp*temp
                  ELSE
                     dphi = dphi + temp*temp
                  ENDIF
                  c = w - Delta(I)*dpsi - Delta(ip1)*dphi
               ELSEIF ( orgati ) THEN
                  c = w - Delta(ip1)*dw - (D(I)-D(ip1))*(Z(I)/Delta(I)) &
     &                **2
               ELSE
                  c = w - Delta(I)*dw - (D(ip1)-D(I))                   &
     &                *(Z(ip1)/Delta(ip1))**2
               ENDIF
               a = (Delta(I)+Delta(ip1))*w - Delta(I)*Delta(ip1)*dw
               b = Delta(I)*Delta(ip1)*w
               IF ( c==ZERO ) THEN
                  IF ( a==ZERO ) THEN
                     IF ( swtch ) THEN
                        a = Delta(I)*Delta(I)*dpsi + Delta(ip1)         &
     &                      *Delta(ip1)*dphi
                     ELSEIF ( orgati ) THEN
                        a = Z(I)*Z(I) + Delta(ip1)*Delta(ip1)           &
     &                      *(dpsi+dphi)
                     ELSE
                        a = Z(ip1)*Z(ip1) + Delta(I)*Delta(I)           &
     &                      *(dpsi+dphi)
                     ENDIF
                  ENDIF
                  eta = b/a
               ELSEIF ( a<=ZERO ) THEN
                  eta = (a-SQRT(ABS(a*a-FOUR*b*c)))/(TWO*c)
               ELSE
                  eta = TWO*b/(a+SQRT(ABS(a*a-FOUR*b*c)))
               ENDIF
            ELSE
!
!              Interpolation using THREE most relevant poles
!
               temp = rhoinv + psi + phi
               IF ( swtch ) THEN
                  c = temp - Delta(iim1)*dpsi - Delta(iip1)*dphi
                  zz(1) = Delta(iim1)*Delta(iim1)*dpsi
                  zz(3) = Delta(iip1)*Delta(iip1)*dphi
               ELSEIF ( orgati ) THEN
                  temp1 = Z(iim1)/Delta(iim1)
                  temp1 = temp1*temp1
                  c = temp - Delta(iip1)*(dpsi+dphi) - (D(iim1)-D(iip1))&
     &                *temp1
                  zz(1) = Z(iim1)*Z(iim1)
                  zz(3) = Delta(iip1)*Delta(iip1)*((dpsi-temp1)+dphi)
               ELSE
                  temp1 = Z(iip1)/Delta(iip1)
                  temp1 = temp1*temp1
                  c = temp - Delta(iim1)*(dpsi+dphi) - (D(iip1)-D(iim1))&
     &                *temp1
                  zz(1) = Delta(iim1)*Delta(iim1)*(dpsi+(dphi-temp1))
                  zz(3) = Z(iip1)*Z(iip1)
               ENDIF
               CALL DLAED6(niter,orgati,c,Delta(iim1),zz,w,eta,Info)
               IF ( Info/=0 ) GOTO 99999
            ENDIF
!
!           Note, eta should be positive if w is negative, and
!           eta should be negative otherwise. However,
!           if for some reason caused by roundoff, eta*w > 0,
!           we simply use one Newton step instead. This way
!           will guarantee eta*w < 0.
!
            IF ( w*eta>=ZERO ) eta = -w/dw
            temp = tau + eta
            IF ( temp>dltub .OR. temp<dltlb ) THEN
               IF ( w<ZERO ) THEN
                  eta = (dltub-tau)/TWO
               ELSE
                  eta = (dltlb-tau)/TWO
               ENDIF
            ENDIF
!
            DO j = 1 , N
               Delta(j) = Delta(j) - eta
            ENDDO
!
            tau = tau + eta
            prew = w
!
!           Evaluate PSI and the derivative DPSI
!
            dpsi = ZERO
            psi = ZERO
            erretm = ZERO
            DO j = 1 , iim1
               temp = Z(j)/Delta(j)
               psi = psi + Z(j)*temp
               dpsi = dpsi + temp*temp
               erretm = erretm + psi
            ENDDO
            erretm = ABS(erretm)
!
!           Evaluate PHI and the derivative DPHI
!
            dphi = ZERO
            phi = ZERO
            DO j = N , iip1 , -1
               temp = Z(j)/Delta(j)
               phi = phi + Z(j)*temp
               dphi = dphi + temp*temp
               erretm = erretm + phi
            ENDDO
!
            temp = Z(ii)/Delta(ii)
            dw = dpsi + dphi + temp*temp
            temp = Z(ii)*temp
            w = rhoinv + phi + psi + temp
            erretm = EIGHT*(phi-psi) + erretm + TWO*rhoinv +            &
     &               THREE*ABS(temp) + ABS(tau)*dw
            IF ( w*prew>ZERO .AND. ABS(w)>ABS(prew)/TEN )               &
     &           swtch = .NOT.swtch
!
         ENDDO
!
!        Return with INFO = 1, NITER = MAXIT and not converged
!
         Info = 1
         IF ( orgati ) THEN
            Dlam = D(I) + tau
         ELSE
            Dlam = D(ip1) + tau
         ENDIF
!
      ENDIF
!
!
!
!     End of DLAED4
!
99999 END SUBROUTINE DLAED4
