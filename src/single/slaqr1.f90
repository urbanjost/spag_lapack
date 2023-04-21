!*==slaqr1.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SLAQR1 sets a scalar multiple of the first column of the product of 2-by-2 or 3-by-3 matrix H and specified shifts.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLAQR1 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaqr1.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaqr1.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaqr1.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLAQR1( N, H, LDH, SR1, SI1, SR2, SI2, V )
!
!       .. Scalar Arguments ..
!       REAL               SI1, SI2, SR1, SR2
!       INTEGER            LDH, N
!       ..
!       .. Array Arguments ..
!       REAL               H( LDH, * ), V( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>      Given a 2-by-2 or 3-by-3 matrix H, SLAQR1 sets v to a
!>      scalar multiple of the first column of the product
!>
!>      (*)  K = (H - (sr1 + i*si1)*I)*(H - (sr2 + i*si2)*I)
!>
!>      scaling to avoid overflows and most underflows. It
!>      is assumed that either
!>
!>              1) sr1 = sr2 and si1 = -si2
!>          or
!>              2) si1 = si2 = 0.
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
!>          H is REAL array, dimension (LDH,N)
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
!> \param[in] SR1
!> \verbatim
!>          SR1 is REAL
!> \endverbatim
!>
!> \param[in] SI1
!> \verbatim
!>          SI1 is REAL
!> \endverbatim
!>
!> \param[in] SR2
!> \verbatim
!>          SR2 is REAL
!> \endverbatim
!>
!> \param[in] SI2
!> \verbatim
!>          SI2 is REAL
!>              The shifts in (*).
!> \endverbatim
!>
!> \param[out] V
!> \verbatim
!>          V is REAL array, dimension (N)
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
!> \ingroup realOTHERauxiliary
!
!> \par Contributors:
!  ==================
!>
!>       Karen Braman and Ralph Byers, Department of Mathematics,
!>       University of Kansas, USA
!>
!  =====================================================================
      SUBROUTINE SLAQR1(N,H,Ldh,Sr1,Si1,Sr2,Si2,V)
      IMPLICIT NONE
!*--SLAQR1125
!
!  -- LAPACK auxiliary routine (version 3.7.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2017
!
!     .. Scalar Arguments ..
      REAL Si1 , Si2 , Sr1 , Sr2
      INTEGER Ldh , N
!     ..
!     .. Array Arguments ..
      REAL H(Ldh,*) , V(*)
!     ..
!
!  ================================================================
!
!     .. Parameters ..
      REAL ZERO
      PARAMETER (ZERO=0.0E0)
!     ..
!     .. Local Scalars ..
      REAL h21s , h31s , s
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF ( N/=2 .AND. N/=3 ) RETURN
!
      IF ( N==2 ) THEN
         s = ABS(H(1,1)-Sr2) + ABS(Si2) + ABS(H(2,1))
         IF ( s==ZERO ) THEN
            V(1) = ZERO
            V(2) = ZERO
         ELSE
            h21s = H(2,1)/s
            V(1) = h21s*H(1,2) + (H(1,1)-Sr1)*((H(1,1)-Sr2)/s)          &
     &             - Si1*(Si2/s)
            V(2) = h21s*(H(1,1)+H(2,2)-Sr1-Sr2)
         ENDIF
      ELSE
         s = ABS(H(1,1)-Sr2) + ABS(Si2) + ABS(H(2,1)) + ABS(H(3,1))
         IF ( s==ZERO ) THEN
            V(1) = ZERO
            V(2) = ZERO
            V(3) = ZERO
         ELSE
            h21s = H(2,1)/s
            h31s = H(3,1)/s
            V(1) = (H(1,1)-Sr1)*((H(1,1)-Sr2)/s) - Si1*(Si2/s) + H(1,2) &
     &             *h21s + H(1,3)*h31s
            V(2) = h21s*(H(1,1)+H(2,2)-Sr1-Sr2) + H(2,3)*h31s
            V(3) = h31s*(H(1,1)+H(3,3)-Sr1-Sr2) + h21s*H(3,2)
         ENDIF
      ENDIF
      END SUBROUTINE SLAQR1
