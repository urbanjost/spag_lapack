!*==zlaev2.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZLAEV2 computes the eigenvalues and eigenvectors of a 2-by-2 symmetric/Hermitian matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLAEV2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaev2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaev2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaev2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLAEV2( A, B, C, RT1, RT2, CS1, SN1 )
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION   CS1, RT1, RT2
!       COMPLEX*16         A, B, C, SN1
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLAEV2 computes the eigendecomposition of a 2-by-2 Hermitian matrix
!>    [  A         B  ]
!>    [  CONJG(B)  C  ].
!> On return, RT1 is the eigenvalue of larger absolute value, RT2 is the
!> eigenvalue of smaller absolute value, and (CS1,SN1) is the unit right
!> eigenvector for RT1, giving the decomposition
!>
!> [ CS1  CONJG(SN1) ] [    A     B ] [ CS1 -CONJG(SN1) ] = [ RT1  0  ]
!> [-SN1     CS1     ] [ CONJG(B) C ] [ SN1     CS1     ]   [  0  RT2 ].
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16
!>         The (1,1) element of the 2-by-2 matrix.
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is COMPLEX*16
!>         The (1,2) element and the conjugate of the (2,1) element of
!>         the 2-by-2 matrix.
!> \endverbatim
!>
!> \param[in] C
!> \verbatim
!>          C is COMPLEX*16
!>         The (2,2) element of the 2-by-2 matrix.
!> \endverbatim
!>
!> \param[out] RT1
!> \verbatim
!>          RT1 is DOUBLE PRECISION
!>         The eigenvalue of larger absolute value.
!> \endverbatim
!>
!> \param[out] RT2
!> \verbatim
!>          RT2 is DOUBLE PRECISION
!>         The eigenvalue of smaller absolute value.
!> \endverbatim
!>
!> \param[out] CS1
!> \verbatim
!>          CS1 is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[out] SN1
!> \verbatim
!>          SN1 is COMPLEX*16
!>         The vector (CS1, SN1) is a unit right eigenvector for RT1.
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
!> \ingroup complex16OTHERauxiliary
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
      SUBROUTINE ZLAEV2(A,B,C,Rt1,Rt2,Cs1,Sn1)
      USE F77KINDS                        
      USE S_DLAEV2
      IMPLICIT NONE
!*--ZLAEV2127
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
!
! Dummy argument declarations rewritten by SPAG
!
      COMPLEX(CX16KIND) , INTENT(IN) :: A
      COMPLEX(CX16KIND) , INTENT(IN) :: B
      COMPLEX(CX16KIND) , INTENT(IN) :: C
      REAL(R8KIND) :: Rt1
      REAL(R8KIND) :: Rt2
      REAL(R8KIND) :: Cs1
      COMPLEX(CX16KIND) , INTENT(OUT) :: Sn1
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: t
      COMPLEX(CX16KIND) :: w
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
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
      IF ( ABS(B)==ZERO ) THEN
         w = ONE
      ELSE
         w = DCONJG(B)/ABS(B)
      ENDIF
      CALL DLAEV2(DBLE(A),ABS(B),DBLE(C),Rt1,Rt2,Cs1,t)
      Sn1 = w*t
!
!     End of ZLAEV2
!
      END SUBROUTINE ZLAEV2
