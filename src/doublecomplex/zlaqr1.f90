!*==zlaqr1.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b ZLAQR1 sets a scalar multiple of the first column of the product of 2-by-2 or 3-by-3 matrix H and specified shifts.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLAQR1 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaqr1.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaqr1.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaqr1.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLAQR1( N, H, LDH, S1, S2, V )
!
!       .. Scalar Arguments ..
!       COMPLEX*16         S1, S2
!       INTEGER            LDH, N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         H( LDH, * ), V( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>      Given a 2-by-2 or 3-by-3 matrix H, ZLAQR1 sets v to a
!>      scalar multiple of the first column of the product
!>
!>      (*)  K = (H - s1*I)*(H - s2*I)
!>
!>      scaling to avoid overflows and most underflows.
!>
!>      This is useful for starting double implicit shift bulges
!>      in the QR algorithm.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>              Order of the matrix H. N must be either 2 or 3.
!> \endverbatim
!>
!> \param[in] H
!> \verbatim
!>          H is COMPLEX*16 array, dimension (LDH,N)
!>              The 2-by-2 or 3-by-3 matrix H in (*).
!> \endverbatim
!>
!> \param[in] LDH
!> \verbatim
!>          LDH is INTEGER
!>              The leading dimension of H as declared in
!>              the calling procedure.  LDH >= N
!> \endverbatim
!>
!> \param[in] S1
!> \verbatim
!>          S1 is COMPLEX*16
!> \endverbatim
!>
!> \param[in] S2
!> \verbatim
!>          S2 is COMPLEX*16
!>
!>          S1 and S2 are the shifts defining K in (*) above.
!> \endverbatim
!>
!> \param[out] V
!> \verbatim
!>          V is COMPLEX*16 array, dimension (N)
!>              A scalar multiple of the first column of the
!>              matrix K in (*).
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
!> \date June 2017
!
!> \ingroup complex16OTHERauxiliary
!
!> \par Contributors:
!  ==================
!>
!>       Karen Braman and Ralph Byers, Department of Mathematics,
!>       University of Kansas, USA
!>
!  =====================================================================
      SUBROUTINE ZLAQR1(N,H,Ldh,S1,S2,V)
      IMPLICIT NONE
!*--ZLAQR1111
!
!  -- LAPACK auxiliary routine (version 3.7.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2017
!
!     .. Scalar Arguments ..
      COMPLEX*16 S1 , S2
      INTEGER Ldh , N
!     ..
!     .. Array Arguments ..
      COMPLEX*16 H(Ldh,*) , V(*)
!     ..
!
!  ================================================================
!
!     .. Parameters ..
      COMPLEX*16 ZERO
      PARAMETER (ZERO=(0.0D0,0.0D0))
      DOUBLE PRECISION RZERO
      PARAMETER (RZERO=0.0D0)
!     ..
!     .. Local Scalars ..
      COMPLEX*16 cdum , h21s , h31s
      DOUBLE PRECISION s
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , DBLE , DIMAG
!     ..
!     .. Statement Functions ..
      DOUBLE PRECISION CABS1
!     ..
!     .. Statement Function definitions ..
      CABS1(cdum) = ABS(DBLE(cdum)) + ABS(DIMAG(cdum))
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF ( N/=2 .AND. N/=3 ) RETURN
!
      IF ( N==2 ) THEN
         s = CABS1(H(1,1)-S2) + CABS1(H(2,1))
         IF ( s==RZERO ) THEN
            V(1) = ZERO
            V(2) = ZERO
         ELSE
            h21s = H(2,1)/s
            V(1) = h21s*H(1,2) + (H(1,1)-S1)*((H(1,1)-S2)/s)
            V(2) = h21s*(H(1,1)+H(2,2)-S1-S2)
         ENDIF
      ELSE
         s = CABS1(H(1,1)-S2) + CABS1(H(2,1)) + CABS1(H(3,1))
         IF ( s==ZERO ) THEN
            V(1) = ZERO
            V(2) = ZERO
            V(3) = ZERO
         ELSE
            h21s = H(2,1)/s
            h31s = H(3,1)/s
            V(1) = (H(1,1)-S1)*((H(1,1)-S2)/s) + H(1,2)*h21s + H(1,3)   &
     &             *h31s
            V(2) = h21s*(H(1,1)+H(2,2)-S1-S2) + H(2,3)*h31s
            V(3) = h31s*(H(1,1)+H(3,3)-S1-S2) + h21s*H(3,2)
         ENDIF
      ENDIF
      END SUBROUTINE ZLAQR1
