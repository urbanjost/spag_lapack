!*==dlaev2.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DLAEV2 computes the eigenvalues and eigenvectors of a 2-by-2 symmetric/Hermitian matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLAEV2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaev2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaev2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaev2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAEV2( A, B, C, RT1, RT2, CS1, SN1 )
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION   A, B, C, CS1, RT1, RT2, SN1
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAEV2 computes the eigendecomposition of a 2-by-2 symmetric matrix
!>    [  A   B  ]
!>    [  B   C  ].
!> On return, RT1 is the eigenvalue of larger absolute value, RT2 is the
!> eigenvalue of smaller absolute value, and (CS1,SN1) is the unit right
!> eigenvector for RT1, giving the decomposition
!>
!>    [ CS1  SN1 ] [  A   B  ] [ CS1 -SN1 ]  =  [ RT1  0  ]
!>    [-SN1  CS1 ] [  B   C  ] [ SN1  CS1 ]     [  0  RT2 ].
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION
!>          The (1,1) element of the 2-by-2 matrix.
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is DOUBLE PRECISION
!>          The (1,2) element and the conjugate of the (2,1) element of
!>          the 2-by-2 matrix.
!> \endverbatim
!>
!> \param[in] C
!> \verbatim
!>          C is DOUBLE PRECISION
!>          The (2,2) element of the 2-by-2 matrix.
!> \endverbatim
!>
!> \param[out] RT1
!> \verbatim
!>          RT1 is DOUBLE PRECISION
!>          The eigenvalue of larger absolute value.
!> \endverbatim
!>
!> \param[out] RT2
!> \verbatim
!>          RT2 is DOUBLE PRECISION
!>          The eigenvalue of smaller absolute value.
!> \endverbatim
!>
!> \param[out] CS1
!> \verbatim
!>          CS1 is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[out] SN1
!> \verbatim
!>          SN1 is DOUBLE PRECISION
!>          The vector (CS1, SN1) is a unit right eigenvector for RT1.
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
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  RT1 is accurate to a few ulps barring over/underflow.
!>
!>  RT2 may be inaccurate if there is massive cancellation in the
!>  determinant A*C-B*B; higher precision or correctly rounded or
!>  correctly truncated arithmetic would be needed to compute RT2
!>  accurately in all cases.
!>
!>  CS1 and SN1 are accurate to a few ulps barring over/underflow.
!>
!>  Overflow is possible only if RT1 is within a factor of 5 of overflow.
!>  Underflow is harmless if the input data is 0 or exceeds
!>     underflow_threshold / macheps.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE DLAEV2(A,B,C,Rt1,Rt2,Cs1,Sn1)
      USE F77KINDS                        
      IMPLICIT NONE
!*--DLAEV2125
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D0 , TWO = 2.0D0 ,         &
     &                              ZERO = 0.0D0 , HALF = 0.5D0
!
! Dummy argument declarations rewritten by SPAG
!
      REAL(R8KIND) , INTENT(IN) :: A
      REAL(R8KIND) , INTENT(IN) :: B
      REAL(R8KIND) , INTENT(IN) :: C
      REAL(R8KIND) , INTENT(INOUT) :: Rt1
      REAL(R8KIND) , INTENT(OUT) :: Rt2
      REAL(R8KIND) , INTENT(INOUT) :: Cs1
      REAL(R8KIND) , INTENT(INOUT) :: Sn1
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: ab , acmn , acmx , acs , adf , cs , ct , df , rt ,&
     &                sm , tb , tn
      INTEGER :: sgn1 , sgn2
!
! End of declarations rewritten by SPAG
!
!     ..
!
! =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Compute the eigenvalues
!
      sm = A + C
      df = A - C
      adf = ABS(df)
      tb = B + B
      ab = ABS(tb)
      IF ( ABS(A)>ABS(C) ) THEN
         acmx = A
         acmn = C
      ELSE
         acmx = C
         acmn = A
      ENDIF
      IF ( adf>ab ) THEN
         rt = adf*SQRT(ONE+(ab/adf)**2)
      ELSEIF ( adf<ab ) THEN
         rt = ab*SQRT(ONE+(adf/ab)**2)
      ELSE
!
!        Includes case AB=ADF=0
!
         rt = ab*SQRT(TWO)
      ENDIF
      IF ( sm<ZERO ) THEN
         Rt1 = HALF*(sm-rt)
         sgn1 = -1
!
!        Order of execution important.
!        To get fully accurate smaller eigenvalue,
!        next line needs to be executed in higher precision.
!
         Rt2 = (acmx/Rt1)*acmn - (B/Rt1)*B
      ELSEIF ( sm>ZERO ) THEN
         Rt1 = HALF*(sm+rt)
         sgn1 = 1
!
!        Order of execution important.
!        To get fully accurate smaller eigenvalue,
!        next line needs to be executed in higher precision.
!
         Rt2 = (acmx/Rt1)*acmn - (B/Rt1)*B
      ELSE
!
!        Includes case RT1 = RT2 = 0
!
         Rt1 = HALF*rt
         Rt2 = -HALF*rt
         sgn1 = 1
      ENDIF
!
!     Compute the eigenvector
!
      IF ( df>=ZERO ) THEN
         cs = df + rt
         sgn2 = 1
      ELSE
         cs = df - rt
         sgn2 = -1
      ENDIF
      acs = ABS(cs)
      IF ( acs>ab ) THEN
         ct = -tb/cs
         Sn1 = ONE/SQRT(ONE+ct*ct)
         Cs1 = ct*Sn1
      ELSEIF ( ab==ZERO ) THEN
         Cs1 = ONE
         Sn1 = ZERO
      ELSE
         tn = -cs/tb
         Cs1 = ONE/SQRT(ONE+tn*tn)
         Sn1 = tn*Cs1
      ENDIF
      IF ( sgn1==sgn2 ) THEN
         tn = Cs1
         Cs1 = -Sn1
         Sn1 = tn
      ENDIF
!
!     End of DLAEV2
!
      END SUBROUTINE DLAEV2
