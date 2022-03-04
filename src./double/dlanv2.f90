!*==dlanv2.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DLANV2 computes the Schur factorization of a real 2-by-2 nonsymmetric matrix in standard form.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLANV2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlanv2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlanv2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlanv2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLANV2( A, B, C, D, RT1R, RT1I, RT2R, RT2I, CS, SN )
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION   A, B, C, CS, D, RT1I, RT1R, RT2I, RT2R, SN
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLANV2 computes the Schur factorization of a real 2-by-2 nonsymmetric
!> matrix in standard form:
!>
!>      [ A  B ] = [ CS -SN ] [ AA  BB ] [ CS  SN ]
!>      [ C  D ]   [ SN  CS ] [ CC  DD ] [-SN  CS ]
!>
!> where either
!> 1) CC = 0 so that AA and DD are real eigenvalues of the matrix, or
!> 2) AA = DD and BB*CC < 0, so that AA + or - sqrt(BB*CC) are complex
!> conjugate eigenvalues.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is DOUBLE PRECISION
!>          On entry, the elements of the input matrix.
!>          On exit, they are overwritten by the elements of the
!>          standardised Schur form.
!> \endverbatim
!>
!> \param[out] RT1R
!> \verbatim
!>          RT1R is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[out] RT1I
!> \verbatim
!>          RT1I is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[out] RT2R
!> \verbatim
!>          RT2R is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[out] RT2I
!> \verbatim
!>          RT2I is DOUBLE PRECISION
!>          The real and imaginary parts of the eigenvalues. If the
!>          eigenvalues are a complex conjugate pair, RT1I > 0.
!> \endverbatim
!>
!> \param[out] CS
!> \verbatim
!>          CS is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[out] SN
!> \verbatim
!>          SN is DOUBLE PRECISION
!>          Parameters of the rotation matrix.
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
!> \ingroup doubleOTHERauxiliary
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Modified by V. Sima, Research Institute for Informatics, Bucharest,
!>  Romania, to reduce the risk of cancellation errors,
!>  when computing real eigenvalues, and to ensure, if possible, that
!>  abs(RT1R) >= abs(RT2R).
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE DLANV2(A,B,C,D,Rt1r,Rt1i,Rt2r,Rt2i,Cs,Sn)
      USE F77KINDS                        
      USE S_DLAMCH
      USE S_DLAPY2
      IMPLICIT NONE
!*--DLANV2134
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , HALF = 0.5D+0 ,     &
     &                              ONE = 1.0D+0 , TWO = 2.0D0 ,        &
     &                              MULTPL = 4.0D+0
!
! Dummy argument declarations rewritten by SPAG
!
      REAL(R8KIND) , INTENT(INOUT) :: A
      REAL(R8KIND) , INTENT(INOUT) :: B
      REAL(R8KIND) , INTENT(INOUT) :: C
      REAL(R8KIND) , INTENT(INOUT) :: D
      REAL(R8KIND) , INTENT(OUT) :: Rt1r
      REAL(R8KIND) , INTENT(INOUT) :: Rt1i
      REAL(R8KIND) , INTENT(OUT) :: Rt2r
      REAL(R8KIND) , INTENT(OUT) :: Rt2i
      REAL(R8KIND) , INTENT(INOUT) :: Cs
      REAL(R8KIND) , INTENT(INOUT) :: Sn
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: aa , bb , bcmax , bcmis , cc , cs1 , dd , eps ,   &
     &                p , sab , sac , safmin , safmn2 , safmx2 , scale ,&
     &                sigma , sn1 , tau , temp , z
      INTEGER :: count
!
! End of declarations rewritten by SPAG
!
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. External Functions ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
      safmin = DLAMCH('S')
      eps = DLAMCH('P')
      safmn2 = DLAMCH('B')**INT(LOG(safmin/eps)/LOG(DLAMCH('B'))/TWO)
      safmx2 = ONE/safmn2
      IF ( C==ZERO ) THEN
         Cs = ONE
         Sn = ZERO
!
      ELSEIF ( B==ZERO ) THEN
!
!        Swap rows and columns
!
         Cs = ZERO
         Sn = ONE
         temp = D
         D = A
         A = temp
         B = -C
         C = ZERO
!
      ELSEIF ( (A-D)==ZERO .AND. SIGN(ONE,B)/=SIGN(ONE,C) ) THEN
         Cs = ONE
         Sn = ZERO
!
      ELSE
!
         temp = A - D
         p = HALF*temp
         bcmax = MAX(ABS(B),ABS(C))
         bcmis = MIN(ABS(B),ABS(C))*SIGN(ONE,B)*SIGN(ONE,C)
         scale = MAX(ABS(p),bcmax)
         z = (p/scale)*p + (bcmax/scale)*bcmis
!
!        If Z is of the order of the machine accuracy, postpone the
!        decision on the nature of eigenvalues
!
         IF ( z>=MULTPL*eps ) THEN
!
!           Real eigenvalues. Compute A and D.
!
            z = p + SIGN(SQRT(scale)*SQRT(z),p)
            A = D + z
            D = D - (bcmax/z)*bcmis
!
!           Compute B and the rotation matrix
!
            tau = DLAPY2(C,z)
            Cs = z/tau
            Sn = C/tau
            B = B - C
            C = ZERO
!
         ELSE
!
!           Complex eigenvalues, or real (almost) equal eigenvalues.
!           Make diagonal elements equal.
!
            count = 0
            sigma = B + C
 20         DO
               count = count + 1
               scale = MAX(ABS(temp),ABS(sigma))
               IF ( scale>=safmx2 ) THEN
                  sigma = sigma*safmn2
                  temp = temp*safmn2
                  IF ( count<=20 ) CYCLE
               ENDIF
               EXIT
            ENDDO
            IF ( scale<=safmn2 ) THEN
               sigma = sigma*safmx2
               temp = temp*safmx2
               IF ( count<=20 ) GOTO 20
            ENDIF
            p = HALF*temp
            tau = DLAPY2(sigma,temp)
            Cs = SQRT(HALF*(ONE+ABS(sigma)/tau))
            Sn = -(p/(tau*Cs))*SIGN(ONE,sigma)
!
!           Compute [ AA  BB ] = [ A  B ] [ CS -SN ]
!                   [ CC  DD ]   [ C  D ] [ SN  CS ]
!
            aa = A*Cs + B*Sn
            bb = -A*Sn + B*Cs
            cc = C*Cs + D*Sn
            dd = -C*Sn + D*Cs
!
!           Compute [ A  B ] = [ CS  SN ] [ AA  BB ]
!                   [ C  D ]   [-SN  CS ] [ CC  DD ]
!
            A = aa*Cs + cc*Sn
            B = bb*Cs + dd*Sn
            C = -aa*Sn + cc*Cs
            D = -bb*Sn + dd*Cs
!
            temp = HALF*(A+D)
            A = temp
            D = temp
!
            IF ( C/=ZERO ) THEN
               IF ( B==ZERO ) THEN
                  B = -C
                  C = ZERO
                  temp = Cs
                  Cs = -Sn
                  Sn = temp
               ELSEIF ( SIGN(ONE,B)==SIGN(ONE,C) ) THEN
!
!                    Real eigenvalues: reduce to upper triangular form
!
                  sab = SQRT(ABS(B))
                  sac = SQRT(ABS(C))
                  p = SIGN(sab*sac,C)
                  tau = ONE/SQRT(ABS(B+C))
                  A = temp + p
                  D = temp - p
                  B = B - C
                  C = ZERO
                  cs1 = sab*tau
                  sn1 = sac*tau
                  temp = Cs*cs1 - Sn*sn1
                  Sn = Cs*sn1 + Sn*cs1
                  Cs = temp
               ENDIF
            ENDIF
         ENDIF
!
      ENDIF
!
!     Store eigenvalues in (RT1R,RT1I) and (RT2R,RT2I).
!
      Rt1r = A
      Rt2r = D
      IF ( C==ZERO ) THEN
         Rt1i = ZERO
         Rt2i = ZERO
      ELSE
         Rt1i = SQRT(ABS(B))*SQRT(ABS(C))
         Rt2i = -Rt1i
      ENDIF
!
!     End of DLANV2
!
      END SUBROUTINE DLANV2
