!*==zlaesy.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZLAESY computes the eigenvalues and eigenvectors of a 2-by-2 complex symmetric matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLAESY + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaesy.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaesy.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaesy.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLAESY( A, B, C, RT1, RT2, EVSCAL, CS1, SN1 )
!
!       .. Scalar Arguments ..
!       COMPLEX*16         A, B, C, CS1, EVSCAL, RT1, RT2, SN1
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLAESY computes the eigendecomposition of a 2-by-2 symmetric matrix
!>    ( ( A, B );( B, C ) )
!> provided the norm of the matrix of eigenvectors is larger than
!> some threshold value.
!>
!> RT1 is the eigenvalue of larger absolute value, and RT2 of
!> smaller absolute value.  If the eigenvectors are computed, then
!> on return ( CS1, SN1 ) is the unit eigenvector for RT1, hence
!>
!> [  CS1     SN1   ] . [ A  B ] . [ CS1    -SN1   ] = [ RT1  0  ]
!> [ -SN1     CS1   ]   [ B  C ]   [ SN1     CS1   ]   [  0  RT2 ]
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16
!>          The ( 1, 1 ) element of input matrix.
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is COMPLEX*16
!>          The ( 1, 2 ) element of input matrix.  The ( 2, 1 ) element
!>          is also given by B, since the 2-by-2 matrix is symmetric.
!> \endverbatim
!>
!> \param[in] C
!> \verbatim
!>          C is COMPLEX*16
!>          The ( 2, 2 ) element of input matrix.
!> \endverbatim
!>
!> \param[out] RT1
!> \verbatim
!>          RT1 is COMPLEX*16
!>          The eigenvalue of larger modulus.
!> \endverbatim
!>
!> \param[out] RT2
!> \verbatim
!>          RT2 is COMPLEX*16
!>          The eigenvalue of smaller modulus.
!> \endverbatim
!>
!> \param[out] EVSCAL
!> \verbatim
!>          EVSCAL is COMPLEX*16
!>          The complex value by which the eigenvector matrix was scaled
!>          to make it orthonormal.  If EVSCAL is zero, the eigenvectors
!>          were not computed.  This means one of two things:  the 2-by-2
!>          matrix could not be diagonalized, or the norm of the matrix
!>          of eigenvectors before scaling was larger than the threshold
!>          value THRESH (set below).
!> \endverbatim
!>
!> \param[out] CS1
!> \verbatim
!>          CS1 is COMPLEX*16
!> \endverbatim
!>
!> \param[out] SN1
!> \verbatim
!>          SN1 is COMPLEX*16
!>          If EVSCAL .NE. 0,  ( CS1, SN1 ) is the unit right eigenvector
!>          for RT1.
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
!> \ingroup complex16SYauxiliary
!
!  =====================================================================
      SUBROUTINE ZLAESY(A,B,C,Rt1,Rt2,Evscal,Cs1,Sn1)
      USE F77KINDS                        
      IMPLICIT NONE
!*--ZLAESY120
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D0,0.0D0)
      REAL(R8KIND) , PARAMETER  ::  HALF = 0.5D0 , THRESH = 0.1D0
!
! Dummy argument declarations rewritten by SPAG
!
      COMPLEX(CX16KIND) , INTENT(IN) :: A
      COMPLEX(CX16KIND) , INTENT(IN) :: B
      COMPLEX(CX16KIND) , INTENT(IN) :: C
      COMPLEX(CX16KIND) , INTENT(INOUT) :: Rt1
      COMPLEX(CX16KIND) , INTENT(INOUT) :: Rt2
      COMPLEX(CX16KIND) , INTENT(INOUT) :: Evscal
      COMPLEX(CX16KIND) , INTENT(OUT) :: Cs1
      COMPLEX(CX16KIND) , INTENT(INOUT) :: Sn1
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: babs , evnorm , tabs , z
      COMPLEX(CX16KIND) :: s , t , tmp
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
!
!     Special case:  The matrix is actually diagonal.
!     To avoid divide by zero later, we treat this case separately.
!
      IF ( ABS(B)==ZERO ) THEN
         Rt1 = A
         Rt2 = C
         IF ( ABS(Rt1)<ABS(Rt2) ) THEN
            tmp = Rt1
            Rt1 = Rt2
            Rt2 = tmp
            Cs1 = ZERO
            Sn1 = ONE
         ELSE
            Cs1 = ONE
            Sn1 = ZERO
         ENDIF
      ELSE
!
!        Compute the eigenvalues and eigenvectors.
!        The characteristic equation is
!           lambda **2 - (A+C) lambda + (A*C - B*B)
!        and we solve it using the quadratic formula.
!
         s = (A+C)*HALF
         t = (A-C)*HALF
!
!        Take the square root carefully to avoid over/under flow.
!
         babs = ABS(B)
         tabs = ABS(t)
         z = MAX(babs,tabs)
         IF ( z>ZERO ) t = z*SQRT((t/z)**2+(B/z)**2)
!
!        Compute the two eigenvalues.  RT1 and RT2 are exchanged
!        if necessary so that RT1 will have the greater magnitude.
!
         Rt1 = s + t
         Rt2 = s - t
         IF ( ABS(Rt1)<ABS(Rt2) ) THEN
            tmp = Rt1
            Rt1 = Rt2
            Rt2 = tmp
         ENDIF
!
!        Choose CS1 = 1 and SN1 to satisfy the first equation, then
!        scale the components of this eigenvector so that the matrix
!        of eigenvectors X satisfies  X * X**T = I .  (No scaling is
!        done if the norm of the eigenvalue matrix is less than THRESH.)
!
         Sn1 = (Rt1-A)/B
         tabs = ABS(Sn1)
         IF ( tabs>ONE ) THEN
            t = tabs*SQRT((ONE/tabs)**2+(Sn1/tabs)**2)
         ELSE
            t = SQRT(CONE+Sn1*Sn1)
         ENDIF
         evnorm = ABS(t)
         IF ( evnorm>=THRESH ) THEN
            Evscal = CONE/t
            Cs1 = Evscal
            Sn1 = Sn1*Evscal
         ELSE
            Evscal = ZERO
         ENDIF
      ENDIF
!
!     End of ZLAESY
!
      END SUBROUTINE ZLAESY
