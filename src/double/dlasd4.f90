!*==dlasd4.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b DLASD4 computes the square root of the i-th updated eigenvalue of a positive symmetric rank-one modification to a positive diagonal matrix. Used by dbdsdc.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLASD4 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasd4.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlasd4.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlasd4.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLASD4( N, I, D, Z, DELTA, RHO, SIGMA, WORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            I, INFO, N
!       DOUBLE PRECISION   RHO, SIGMA
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   D( * ), DELTA( * ), WORK( * ), Z( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> This subroutine computes the square root of the I-th updated
!> eigenvalue of a positive symmetric rank-one modification to
!> a positive diagonal matrix whose entries are given as the squares
!> of the corresponding entries in the array d, and that
!>
!>        0 <= D(i) < D(j)  for  i < j
!>
!> and that RHO > 0. This is arranged by the calling routine, and is
!> no loss in generality.  The rank-one modified system is thus
!>
!>        diag( D ) * diag( D ) +  RHO * Z * Z_transpose.
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
!>          D is DOUBLE PRECISION array, dimension ( N )
!>         The original eigenvalues.  It is assumed that they are in
!>         order, 0 <= D(I) < D(J)  for I < J.
!> \endverbatim
!>
!> \param[in] Z
!> \verbatim
!>          Z is DOUBLE PRECISION array, dimension ( N )
!>         The components of the updating vector.
!> \endverbatim
!>
!> \param[out] DELTA
!> \verbatim
!>          DELTA is DOUBLE PRECISION array, dimension ( N )
!>         If N .ne. 1, DELTA contains (D(j) - sigma_I) in its  j-th
!>         component.  If N = 1, then DELTA(1) = 1.  The vector DELTA
!>         contains the information necessary to construct the
!>         (singular) eigenvectors.
!> \endverbatim
!>
!> \param[in] RHO
!> \verbatim
!>          RHO is DOUBLE PRECISION
!>         The scalar in the symmetric updating formula.
!> \endverbatim
!>
!> \param[out] SIGMA
!> \verbatim
!>          SIGMA is DOUBLE PRECISION
!>         The computed sigma_I, the I-th updated eigenvalue.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension ( N )
!>         If N .ne. 1, WORK contains (D(j) + sigma_I) in its  j-th
!>         component.  If N = 1, then WORK( 1 ) = 1.
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
!>  Logical variable SWTCH3 (switch-for-3-poles?) is for noting
!>  if we are working with THREE poles!
!>
!>  MAXIT is the maximum number of iterations allowed for each
!>  eigenvalue.
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
!> \ingroup OTHERauxiliary
!
!> \par Contributors:
!  ==================
!>
!>     Ren-Cang Li, Computer Science Division, University of California
!>     at Berkeley, USA
!>
!  =====================================================================
      SUBROUTINE DLASD4(N,I,D,Z,Delta,Rho,Sigma,Work,Info)
      IMPLICIT NONE
!*--DLASD4157
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER I , Info , N
      DOUBLE PRECISION Rho , Sigma
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION D(*) , Delta(*) , Work(*) , Z(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER MAXIT
      PARAMETER (MAXIT=400)
      DOUBLE PRECISION ZERO , ONE , TWO , THREE , FOUR , EIGHT , TEN
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0,TWO=2.0D+0,THREE=3.0D+0,        &
     &           FOUR=4.0D+0,EIGHT=8.0D+0,TEN=10.0D+0)
!     ..
!     .. Local Scalars ..
      LOGICAL orgati , swtch , swtch3 , geomavg
      INTEGER ii , iim1 , iip1 , ip1 , iter , j , niter
      DOUBLE PRECISION a , b , c , delsq , delsq2 , sq2 , dphi , dpsi , &
     &                 dtiim , dtiip , dtipsq , dtisq , dtnsq , dtnsq1 ,&
     &                 dw , eps , erretm , eta , phi , prew , psi ,     &
     &                 rhoinv , sglb , sgub , tau , tau2 , temp ,       &
     &                 temp1 , temp2 , w
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION dd(3) , zz(3)
!     ..
!     .. External Subroutines ..
      EXTERNAL DLAED6 , DLASD5
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH
      EXTERNAL DLAMCH
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , MIN , SQRT
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
!        Presumably, I=1 upon entry
!
         Sigma = SQRT(D(1)*D(1)+Rho*Z(1)*Z(1))
         Delta(1) = ONE
         Work(1) = ONE
         RETURN
      ENDIF
      IF ( N==2 ) THEN
         CALL DLASD5(I,D,Z,Delta,Rho,Sigma,Work)
         RETURN
      ENDIF
!
!     Compute machine epsilon
!
      eps = DLAMCH('Epsilon')
      rhoinv = ONE/Rho
      tau2 = ZERO
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
         temp = Rho/TWO
!
!        If ||Z||_2 is not one, then TEMP should be set to
!        RHO * ||Z||_2^2 / TWO
!
         temp1 = temp/(D(N)+SQRT(D(N)*D(N)+temp))
         DO j = 1 , N
            Work(j) = D(j) + D(N) + temp1
            Delta(j) = (D(j)-D(N)) - temp1
         ENDDO
!
         psi = ZERO
         DO j = 1 , N - 2
            psi = psi + Z(j)*Z(j)/(Delta(j)*Work(j))
         ENDDO
!
         c = rhoinv + psi
         w = c + Z(ii)*Z(ii)/(Delta(ii)*Work(ii)) + Z(N)*Z(N)           &
     &       /(Delta(N)*Work(N))
!
         IF ( w<=ZERO ) THEN
            temp1 = SQRT(D(N)*D(N)+Rho)
            temp = Z(N-1)*Z(N-1)                                        &
     &             /((D(N-1)+temp1)*(D(N)-D(N-1)+Rho/(D(N)+temp1)))     &
     &             + Z(N)*Z(N)/Rho
!
!           The following TAU2 is to approximate
!           SIGMA_n^2 - D( N )*D( N )
!
            IF ( c<=temp ) THEN
               tau = Rho
            ELSE
               delsq = (D(N)-D(N-1))*(D(N)+D(N-1))
               a = -c*delsq + Z(N-1)*Z(N-1) + Z(N)*Z(N)
               b = Z(N)*Z(N)*delsq
               IF ( a<ZERO ) THEN
                  tau2 = TWO*b/(SQRT(a*a+FOUR*b*c)-a)
               ELSE
                  tau2 = (a+SQRT(a*a+FOUR*b*c))/(TWO*c)
               ENDIF
               tau = tau2/(D(N)+SQRT(D(N)*D(N)+tau2))
            ENDIF
!
!           It can be proved that
!               D(N)^2+RHO/2 <= SIGMA_n^2 < D(N)^2+TAU2 <= D(N)^2+RHO
!
         ELSE
            delsq = (D(N)-D(N-1))*(D(N)+D(N-1))
            a = -c*delsq + Z(N-1)*Z(N-1) + Z(N)*Z(N)
            b = Z(N)*Z(N)*delsq
!
!           The following TAU2 is to approximate
!           SIGMA_n^2 - D( N )*D( N )
!
            IF ( a<ZERO ) THEN
               tau2 = TWO*b/(SQRT(a*a+FOUR*b*c)-a)
            ELSE
               tau2 = (a+SQRT(a*a+FOUR*b*c))/(TWO*c)
            ENDIF
            tau = tau2/(D(N)+SQRT(D(N)*D(N)+tau2))
 
!
!           It can be proved that
!           D(N)^2 < D(N)^2+TAU2 < SIGMA(N)^2 < D(N)^2+RHO/2
!
         ENDIF
!
!        The following TAU is to approximate SIGMA_n - D( N )
!
!         TAU = TAU2 / ( D( N )+SQRT( D( N )*D( N )+TAU2 ) )
!
         Sigma = D(N) + tau
         DO j = 1 , N
            Delta(j) = (D(j)-D(N)) - tau
            Work(j) = D(j) + D(N) + tau
         ENDDO
!
!        Evaluate PSI and the derivative DPSI
!
         dpsi = ZERO
         psi = ZERO
         erretm = ZERO
         DO j = 1 , ii
            temp = Z(j)/(Delta(j)*Work(j))
            psi = psi + Z(j)*temp
            dpsi = dpsi + temp*temp
            erretm = erretm + psi
         ENDDO
         erretm = ABS(erretm)
!
!        Evaluate PHI and the derivative DPHI
!
         temp = Z(N)/(Delta(N)*Work(N))
         phi = Z(N)*temp
         dphi = temp*temp
         erretm = EIGHT*(-phi-psi) + erretm - phi + rhoinv
!    $          + ABS( TAU2 )*( DPSI+DPHI )
!
         w = rhoinv + phi + psi
!
!        Test for convergence
!
         IF ( ABS(w)<=eps*erretm ) GOTO 99999
!
!        Calculate the new step
!
         niter = niter + 1
         dtnsq1 = Work(N-1)*Delta(N-1)
         dtnsq = Work(N)*Delta(N)
         c = w - dtnsq1*dpsi - dtnsq*dphi
         a = (dtnsq+dtnsq1)*w - dtnsq*dtnsq1*(dpsi+dphi)
         b = dtnsq*dtnsq1*w
         IF ( c<ZERO ) c = ABS(c)
         IF ( c==ZERO ) THEN
            eta = Rho - Sigma*Sigma
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
         temp = eta - dtnsq
         IF ( temp>Rho ) eta = Rho + dtnsq
!
         eta = eta/(Sigma+SQRT(eta+Sigma*Sigma))
         tau = tau + eta
         Sigma = Sigma + eta
!
         DO j = 1 , N
            Delta(j) = Delta(j) - eta
            Work(j) = Work(j) + eta
         ENDDO
!
!        Evaluate PSI and the derivative DPSI
!
         dpsi = ZERO
         psi = ZERO
         erretm = ZERO
         DO j = 1 , ii
            temp = Z(j)/(Work(j)*Delta(j))
            psi = psi + Z(j)*temp
            dpsi = dpsi + temp*temp
            erretm = erretm + psi
         ENDDO
         erretm = ABS(erretm)
!
!        Evaluate PHI and the derivative DPHI
!
         tau2 = Work(N)*Delta(N)
         temp = Z(N)/tau2
         phi = Z(N)*temp
         dphi = temp*temp
         erretm = EIGHT*(-phi-psi) + erretm - phi + rhoinv
!    $          + ABS( TAU2 )*( DPSI+DPHI )
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
            IF ( ABS(w)<=eps*erretm ) GOTO 99999
!
!           Calculate the new step
!
            dtnsq1 = Work(N-1)*Delta(N-1)
            dtnsq = Work(N)*Delta(N)
            c = w - dtnsq1*dpsi - dtnsq*dphi
            a = (dtnsq+dtnsq1)*w - dtnsq1*dtnsq*(dpsi+dphi)
            b = dtnsq1*dtnsq*w
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
            temp = eta - dtnsq
            IF ( temp<=ZERO ) eta = eta/TWO
!
            eta = eta/(Sigma+SQRT(eta+Sigma*Sigma))
            tau = tau + eta
            Sigma = Sigma + eta
!
            DO j = 1 , N
               Delta(j) = Delta(j) - eta
               Work(j) = Work(j) + eta
            ENDDO
!
!           Evaluate PSI and the derivative DPSI
!
            dpsi = ZERO
            psi = ZERO
            erretm = ZERO
            DO j = 1 , ii
               temp = Z(j)/(Work(j)*Delta(j))
               psi = psi + Z(j)*temp
               dpsi = dpsi + temp*temp
               erretm = erretm + psi
            ENDDO
            erretm = ABS(erretm)
!
!           Evaluate PHI and the derivative DPHI
!
            tau2 = Work(N)*Delta(N)
            temp = Z(N)/tau2
            phi = Z(N)*temp
            dphi = temp*temp
            erretm = EIGHT*(-phi-psi) + erretm - phi + rhoinv
!    $             + ABS( TAU2 )*( DPSI+DPHI )
!
            w = rhoinv + phi + psi
         ENDDO
!
!        Return with INFO = 1, NITER = MAXIT and not converged
!
         Info = 1
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
         delsq = (D(ip1)-D(I))*(D(ip1)+D(I))
         delsq2 = delsq/TWO
         sq2 = SQRT((D(I)*D(I)+D(ip1)*D(ip1))/TWO)
         temp = delsq2/(D(I)+sq2)
         DO j = 1 , N
            Work(j) = D(j) + D(I) + temp
            Delta(j) = (D(j)-D(I)) - temp
         ENDDO
!
         psi = ZERO
         DO j = 1 , I - 1
            psi = psi + Z(j)*Z(j)/(Work(j)*Delta(j))
         ENDDO
!
         phi = ZERO
         DO j = N , I + 2 , -1
            phi = phi + Z(j)*Z(j)/(Work(j)*Delta(j))
         ENDDO
         c = rhoinv + psi + phi
         w = c + Z(I)*Z(I)/(Work(I)*Delta(I)) + Z(ip1)*Z(ip1)           &
     &       /(Work(ip1)*Delta(ip1))
!
         geomavg = .FALSE.
         IF ( w>ZERO ) THEN
!
!           d(i)^2 < the ith sigma^2 < (d(i)^2+d(i+1)^2)/2
!
!           We choose d(i) as origin.
!
            orgati = .TRUE.
            ii = I
            sglb = ZERO
            sgub = delsq2/(D(I)+sq2)
            a = c*delsq + Z(I)*Z(I) + Z(ip1)*Z(ip1)
            b = Z(I)*Z(I)*delsq
            IF ( a>ZERO ) THEN
               tau2 = TWO*b/(a+SQRT(ABS(a*a-FOUR*b*c)))
            ELSE
               tau2 = (a-SQRT(ABS(a*a-FOUR*b*c)))/(TWO*c)
            ENDIF
!
!           TAU2 now is an estimation of SIGMA^2 - D( I )^2. The
!           following, however, is the corresponding estimation of
!           SIGMA - D( I ).
!
            tau = tau2/(D(I)+SQRT(D(I)*D(I)+tau2))
            temp = SQRT(eps)
            IF ( (D(I)<=temp*D(ip1)) .AND. (ABS(Z(I))<=temp) .AND.      &
     &           (D(I)>ZERO) ) THEN
               tau = MIN(TEN*D(I),sgub)
               geomavg = .TRUE.
            ENDIF
         ELSE
!
!           (d(i)^2+d(i+1)^2)/2 <= the ith sigma^2 < d(i+1)^2/2
!
!           We choose d(i+1) as origin.
!
            orgati = .FALSE.
            ii = ip1
            sglb = -delsq2/(D(ii)+sq2)
            sgub = ZERO
            a = c*delsq - Z(I)*Z(I) - Z(ip1)*Z(ip1)
            b = Z(ip1)*Z(ip1)*delsq
            IF ( a<ZERO ) THEN
               tau2 = TWO*b/(a-SQRT(ABS(a*a+FOUR*b*c)))
            ELSE
               tau2 = -(a+SQRT(ABS(a*a+FOUR*b*c)))/(TWO*c)
            ENDIF
!
!           TAU2 now is an estimation of SIGMA^2 - D( IP1 )^2. The
!           following, however, is the corresponding estimation of
!           SIGMA - D( IP1 ).
!
            tau = tau2/(D(ip1)+SQRT(ABS(D(ip1)*D(ip1)+tau2)))
         ENDIF
!
         Sigma = D(ii) + tau
         DO j = 1 , N
            Work(j) = D(j) + D(ii) + tau
            Delta(j) = (D(j)-D(ii)) - tau
         ENDDO
         iim1 = ii - 1
         iip1 = ii + 1
!
!        Evaluate PSI and the derivative DPSI
!
         dpsi = ZERO
         psi = ZERO
         erretm = ZERO
         DO j = 1 , iim1
            temp = Z(j)/(Work(j)*Delta(j))
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
            temp = Z(j)/(Work(j)*Delta(j))
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
         temp = Z(ii)/(Work(ii)*Delta(ii))
         dw = dpsi + dphi + temp*temp
         temp = Z(ii)*temp
         w = w + temp
         erretm = EIGHT*(phi-psi) + erretm + TWO*rhoinv +               &
     &            THREE*ABS(temp)
!    $          + ABS( TAU2 )*DW
!
!        Test for convergence
!
         IF ( ABS(w)<=eps*erretm ) GOTO 99999
!
         IF ( w<=ZERO ) THEN
            sglb = MAX(sglb,tau)
         ELSE
            sgub = MIN(sgub,tau)
         ENDIF
!
!        Calculate the new step
!
         niter = niter + 1
         IF ( .NOT.swtch3 ) THEN
            dtipsq = Work(ip1)*Delta(ip1)
            dtisq = Work(I)*Delta(I)
            IF ( orgati ) THEN
               c = w - dtipsq*dw + delsq*(Z(I)/dtisq)**2
            ELSE
               c = w - dtisq*dw - delsq*(Z(ip1)/dtipsq)**2
            ENDIF
            a = (dtipsq+dtisq)*w - dtipsq*dtisq*dw
            b = dtipsq*dtisq*w
            IF ( c==ZERO ) THEN
               IF ( a==ZERO ) THEN
                  IF ( orgati ) THEN
                     a = Z(I)*Z(I) + dtipsq*dtipsq*(dpsi+dphi)
                  ELSE
                     a = Z(ip1)*Z(ip1) + dtisq*dtisq*(dpsi+dphi)
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
            dtiim = Work(iim1)*Delta(iim1)
            dtiip = Work(iip1)*Delta(iip1)
            temp = rhoinv + psi + phi
            IF ( orgati ) THEN
               temp1 = Z(iim1)/dtiim
               temp1 = temp1*temp1
               c = (temp-dtiip*(dpsi+dphi)) - (D(iim1)-D(iip1))         &
     &             *(D(iim1)+D(iip1))*temp1
               zz(1) = Z(iim1)*Z(iim1)
               IF ( dpsi<temp1 ) THEN
                  zz(3) = dtiip*dtiip*dphi
               ELSE
                  zz(3) = dtiip*dtiip*((dpsi-temp1)+dphi)
               ENDIF
            ELSE
               temp1 = Z(iip1)/dtiip
               temp1 = temp1*temp1
               c = (temp-dtiim*(dpsi+dphi)) - (D(iip1)-D(iim1))         &
     &             *(D(iim1)+D(iip1))*temp1
               IF ( dphi<temp1 ) THEN
                  zz(1) = dtiim*dtiim*dpsi
               ELSE
                  zz(1) = dtiim*dtiim*(dpsi+(dphi-temp1))
               ENDIF
               zz(3) = Z(iip1)*Z(iip1)
            ENDIF
            zz(2) = Z(ii)*Z(ii)
            dd(1) = dtiim
            dd(2) = Delta(ii)*Work(ii)
            dd(3) = dtiip
            CALL DLAED6(niter,orgati,c,dd,zz,w,eta,Info)
!
            IF ( Info/=0 ) THEN
!
!              If INFO is not 0, i.e., DLAED6 failed, switch back
!              to 2 pole interpolation.
!
               swtch3 = .FALSE.
               Info = 0
               dtipsq = Work(ip1)*Delta(ip1)
               dtisq = Work(I)*Delta(I)
               IF ( orgati ) THEN
                  c = w - dtipsq*dw + delsq*(Z(I)/dtisq)**2
               ELSE
                  c = w - dtisq*dw - delsq*(Z(ip1)/dtipsq)**2
               ENDIF
               a = (dtipsq+dtisq)*w - dtipsq*dtisq*dw
               b = dtipsq*dtisq*w
               IF ( c==ZERO ) THEN
                  IF ( a==ZERO ) THEN
                     IF ( orgati ) THEN
                        a = Z(I)*Z(I) + dtipsq*dtipsq*(dpsi+dphi)
                     ELSE
                        a = Z(ip1)*Z(ip1) + dtisq*dtisq*(dpsi+dphi)
                     ENDIF
                  ENDIF
                  eta = b/a
               ELSEIF ( a<=ZERO ) THEN
                  eta = (a-SQRT(ABS(a*a-FOUR*b*c)))/(TWO*c)
               ELSE
                  eta = TWO*b/(a+SQRT(ABS(a*a-FOUR*b*c)))
               ENDIF
            ENDIF
         ENDIF
!
!        Note, eta should be positive if w is negative, and
!        eta should be negative otherwise. However,
!        if for some reason caused by roundoff, eta*w > 0,
!        we simply use one Newton step instead. This way
!        will guarantee eta*w < 0.
!
         IF ( w*eta>=ZERO ) eta = -w/dw
!
         eta = eta/(Sigma+SQRT(Sigma*Sigma+eta))
         temp = tau + eta
         IF ( temp>sgub .OR. temp<sglb ) THEN
            IF ( w<ZERO ) THEN
               eta = (sgub-tau)/TWO
            ELSE
               eta = (sglb-tau)/TWO
            ENDIF
            IF ( geomavg ) THEN
               IF ( w<ZERO ) THEN
                  IF ( tau>ZERO ) eta = SQRT(sgub*tau) - tau
               ELSEIF ( sglb>ZERO ) THEN
                  eta = SQRT(sglb*tau) - tau
               ENDIF
            ENDIF
         ENDIF
!
         prew = w
!
         tau = tau + eta
         Sigma = Sigma + eta
!
         DO j = 1 , N
            Work(j) = Work(j) + eta
            Delta(j) = Delta(j) - eta
         ENDDO
!
!        Evaluate PSI and the derivative DPSI
!
         dpsi = ZERO
         psi = ZERO
         erretm = ZERO
         DO j = 1 , iim1
            temp = Z(j)/(Work(j)*Delta(j))
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
            temp = Z(j)/(Work(j)*Delta(j))
            phi = phi + Z(j)*temp
            dphi = dphi + temp*temp
            erretm = erretm + phi
         ENDDO
!
         tau2 = Work(ii)*Delta(ii)
         temp = Z(ii)/tau2
         dw = dpsi + dphi + temp*temp
         temp = Z(ii)*temp
         w = rhoinv + phi + psi + temp
         erretm = EIGHT*(phi-psi) + erretm + TWO*rhoinv +               &
     &            THREE*ABS(temp)
!    $          + ABS( TAU2 )*DW
!
         swtch = .FALSE.
         IF ( orgati ) THEN
            IF ( -w>ABS(prew)/TEN ) swtch = .TRUE.
         ELSE
            IF ( w>ABS(prew)/TEN ) swtch = .TRUE.
         ENDIF
!
!        Main loop to update the values of the array   DELTA and WORK
!
         iter = niter + 1
!
         DO niter = iter , MAXIT
!
!           Test for convergence
!
!     $          .OR. (SGUB-SGLB).LE.EIGHT*ABS(SGUB+SGLB) ) THEN
            IF ( ABS(w)<=eps*erretm ) GOTO 99999
!
            IF ( w<=ZERO ) THEN
               sglb = MAX(sglb,tau)
            ELSE
               sgub = MIN(sgub,tau)
            ENDIF
!
!           Calculate the new step
!
            IF ( .NOT.swtch3 ) THEN
               dtipsq = Work(ip1)*Delta(ip1)
               dtisq = Work(I)*Delta(I)
               IF ( swtch ) THEN
                  temp = Z(ii)/(Work(ii)*Delta(ii))
                  IF ( orgati ) THEN
                     dpsi = dpsi + temp*temp
                  ELSE
                     dphi = dphi + temp*temp
                  ENDIF
                  c = w - dtisq*dpsi - dtipsq*dphi
               ELSEIF ( orgati ) THEN
                  c = w - dtipsq*dw + delsq*(Z(I)/dtisq)**2
               ELSE
                  c = w - dtisq*dw - delsq*(Z(ip1)/dtipsq)**2
               ENDIF
               a = (dtipsq+dtisq)*w - dtipsq*dtisq*dw
               b = dtipsq*dtisq*w
               IF ( c==ZERO ) THEN
                  IF ( a==ZERO ) THEN
                     IF ( swtch ) THEN
                        a = dtisq*dtisq*dpsi + dtipsq*dtipsq*dphi
                     ELSEIF ( orgati ) THEN
                        a = Z(I)*Z(I) + dtipsq*dtipsq*(dpsi+dphi)
                     ELSE
                        a = Z(ip1)*Z(ip1) + dtisq*dtisq*(dpsi+dphi)
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
               dtiim = Work(iim1)*Delta(iim1)
               dtiip = Work(iip1)*Delta(iip1)
               temp = rhoinv + psi + phi
               IF ( swtch ) THEN
                  c = temp - dtiim*dpsi - dtiip*dphi
                  zz(1) = dtiim*dtiim*dpsi
                  zz(3) = dtiip*dtiip*dphi
               ELSEIF ( orgati ) THEN
                  temp1 = Z(iim1)/dtiim
                  temp1 = temp1*temp1
                  temp2 = (D(iim1)-D(iip1))*(D(iim1)+D(iip1))*temp1
                  c = temp - dtiip*(dpsi+dphi) - temp2
                  zz(1) = Z(iim1)*Z(iim1)
                  IF ( dpsi<temp1 ) THEN
                     zz(3) = dtiip*dtiip*dphi
                  ELSE
                     zz(3) = dtiip*dtiip*((dpsi-temp1)+dphi)
                  ENDIF
               ELSE
                  temp1 = Z(iip1)/dtiip
                  temp1 = temp1*temp1
                  temp2 = (D(iip1)-D(iim1))*(D(iim1)+D(iip1))*temp1
                  c = temp - dtiim*(dpsi+dphi) - temp2
                  IF ( dphi<temp1 ) THEN
                     zz(1) = dtiim*dtiim*dpsi
                  ELSE
                     zz(1) = dtiim*dtiim*(dpsi+(dphi-temp1))
                  ENDIF
                  zz(3) = Z(iip1)*Z(iip1)
               ENDIF
               dd(1) = dtiim
               dd(2) = Delta(ii)*Work(ii)
               dd(3) = dtiip
               CALL DLAED6(niter,orgati,c,dd,zz,w,eta,Info)
!
               IF ( Info/=0 ) THEN
!
!                 If INFO is not 0, i.e., DLAED6 failed, switch
!                 back to two pole interpolation
!
                  swtch3 = .FALSE.
                  Info = 0
                  dtipsq = Work(ip1)*Delta(ip1)
                  dtisq = Work(I)*Delta(I)
                  IF ( swtch ) THEN
                     temp = Z(ii)/(Work(ii)*Delta(ii))
                     IF ( orgati ) THEN
                        dpsi = dpsi + temp*temp
                     ELSE
                        dphi = dphi + temp*temp
                     ENDIF
                     c = w - dtisq*dpsi - dtipsq*dphi
                  ELSEIF ( orgati ) THEN
                     c = w - dtipsq*dw + delsq*(Z(I)/dtisq)**2
                  ELSE
                     c = w - dtisq*dw - delsq*(Z(ip1)/dtipsq)**2
                  ENDIF
                  a = (dtipsq+dtisq)*w - dtipsq*dtisq*dw
                  b = dtipsq*dtisq*w
                  IF ( c==ZERO ) THEN
                     IF ( a==ZERO ) THEN
                        IF ( swtch ) THEN
                           a = dtisq*dtisq*dpsi + dtipsq*dtipsq*dphi
                        ELSEIF ( orgati ) THEN
                           a = Z(I)*Z(I) + dtipsq*dtipsq*(dpsi+dphi)
                        ELSE
                           a = Z(ip1)*Z(ip1) + dtisq*dtisq*(dpsi+dphi)
                        ENDIF
                     ENDIF
                     eta = b/a
                  ELSEIF ( a<=ZERO ) THEN
                     eta = (a-SQRT(ABS(a*a-FOUR*b*c)))/(TWO*c)
                  ELSE
                     eta = TWO*b/(a+SQRT(ABS(a*a-FOUR*b*c)))
                  ENDIF
               ENDIF
            ENDIF
!
!           Note, eta should be positive if w is negative, and
!           eta should be negative otherwise. However,
!           if for some reason caused by roundoff, eta*w > 0,
!           we simply use one Newton step instead. This way
!           will guarantee eta*w < 0.
!
            IF ( w*eta>=ZERO ) eta = -w/dw
!
            eta = eta/(Sigma+SQRT(Sigma*Sigma+eta))
            temp = tau + eta
            IF ( temp>sgub .OR. temp<sglb ) THEN
               IF ( w<ZERO ) THEN
                  eta = (sgub-tau)/TWO
               ELSE
                  eta = (sglb-tau)/TWO
               ENDIF
               IF ( geomavg ) THEN
                  IF ( w<ZERO ) THEN
                     IF ( tau>ZERO ) eta = SQRT(sgub*tau) - tau
                  ELSEIF ( sglb>ZERO ) THEN
                     eta = SQRT(sglb*tau) - tau
                  ENDIF
               ENDIF
            ENDIF
!
            prew = w
!
            tau = tau + eta
            Sigma = Sigma + eta
!
            DO j = 1 , N
               Work(j) = Work(j) + eta
               Delta(j) = Delta(j) - eta
            ENDDO
!
!           Evaluate PSI and the derivative DPSI
!
            dpsi = ZERO
            psi = ZERO
            erretm = ZERO
            DO j = 1 , iim1
               temp = Z(j)/(Work(j)*Delta(j))
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
               temp = Z(j)/(Work(j)*Delta(j))
               phi = phi + Z(j)*temp
               dphi = dphi + temp*temp
               erretm = erretm + phi
            ENDDO
!
            tau2 = Work(ii)*Delta(ii)
            temp = Z(ii)/tau2
            dw = dpsi + dphi + temp*temp
            temp = Z(ii)*temp
            w = rhoinv + phi + psi + temp
            erretm = EIGHT*(phi-psi) + erretm + TWO*rhoinv +            &
     &               THREE*ABS(temp)
!    $             + ABS( TAU2 )*DW
!
            IF ( w*prew>ZERO .AND. ABS(w)>ABS(prew)/TEN )               &
     &           swtch = .NOT.swtch
!
         ENDDO
!
!        Return with INFO = 1, NITER = MAXIT and not converged
!
         Info = 1
!
      ENDIF
!
!
!     End of DLASD4
!
99999 END SUBROUTINE DLASD4
