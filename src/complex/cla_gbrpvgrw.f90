!*==cla_gbrpvgrw.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b CLA_GBRPVGRW computes the reciprocal pivot growth factor norm(A)/norm(U) for a general banded matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLA_GBRPVGRW + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_gbrpvgrw.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_gbrpvgrw.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_gbrpvgrw.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       REAL FUNCTION CLA_GBRPVGRW( N, KL, KU, NCOLS, AB, LDAB, AFB,
!                                   LDAFB )
!
!       .. Scalar Arguments ..
!       INTEGER            N, KL, KU, NCOLS, LDAB, LDAFB
!       ..
!       .. Array Arguments ..
!       COMPLEX            AB( LDAB, * ), AFB( LDAFB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLA_GBRPVGRW computes the reciprocal pivot growth factor
!> norm(A)/norm(U). The "max absolute element" norm is used. If this is
!> much less than 1, the stability of the LU factorization of the
!> (equilibrated) matrix A could be poor. This also means that the
!> solution X, estimated condition numbers, and error bounds could be
!> unreliable.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>     The number of linear equations, i.e., the order of the
!>     matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] KL
!> \verbatim
!>          KL is INTEGER
!>     The number of subdiagonals within the band of A.  KL >= 0.
!> \endverbatim
!>
!> \param[in] KU
!> \verbatim
!>          KU is INTEGER
!>     The number of superdiagonals within the band of A.  KU >= 0.
!> \endverbatim
!>
!> \param[in] NCOLS
!> \verbatim
!>          NCOLS is INTEGER
!>     The number of columns of the matrix A.  NCOLS >= 0.
!> \endverbatim
!>
!> \param[in] AB
!> \verbatim
!>          AB is COMPLEX array, dimension (LDAB,N)
!>     On entry, the matrix A in band storage, in rows 1 to KL+KU+1.
!>     The j-th column of A is stored in the j-th column of the
!>     array AB as follows:
!>     AB(KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+kl)
!> \endverbatim
!>
!> \param[in] LDAB
!> \verbatim
!>          LDAB is INTEGER
!>     The leading dimension of the array AB.  LDAB >= KL+KU+1.
!> \endverbatim
!>
!> \param[in] AFB
!> \verbatim
!>          AFB is COMPLEX array, dimension (LDAFB,N)
!>     Details of the LU factorization of the band matrix A, as
!>     computed by CGBTRF.  U is stored as an upper triangular
!>     band matrix with KL+KU superdiagonals in rows 1 to KL+KU+1,
!>     and the multipliers used during the factorization are stored
!>     in rows KL+KU+2 to 2*KL+KU+1.
!> \endverbatim
!>
!> \param[in] LDAFB
!> \verbatim
!>          LDAFB is INTEGER
!>     The leading dimension of the array AFB.  LDAFB >= 2*KL+KU+1.
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
!> \ingroup complexGBcomputational
!
!  =====================================================================
      REAL FUNCTION CLA_GBRPVGRW(N,Kl,Ku,Ncols,Ab,Ldab,Afb,Ldafb)
      IMPLICIT NONE
!*--CLA_GBRPVGRW120
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER N , Kl , Ku , Ncols , Ldab , Ldafb
!     ..
!     .. Array Arguments ..
      COMPLEX Ab(Ldab,*) , Afb(Ldafb,*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER i , j , kd
      REAL amax , umax , rpvgrw
      COMPLEX zdum
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , MIN , REAL , AIMAG
!     ..
!     .. Statement Functions ..
      REAL CABS1
!     ..
!     .. Statement Function Definitions ..
      CABS1(zdum) = ABS(REAL(zdum)) + ABS(AIMAG(zdum))
!     ..
!     .. Executable Statements ..
!
      rpvgrw = 1.0
 
      kd = Ku + 1
      DO j = 1 , Ncols
         amax = 0.0
         umax = 0.0
         DO i = MAX(j-Ku,1) , MIN(j+Kl,N)
            amax = MAX(CABS1(Ab(kd+i-j,j)),amax)
         ENDDO
         DO i = MAX(j-Ku,1) , j
            umax = MAX(CABS1(Afb(kd+i-j,j)),umax)
         ENDDO
         IF ( umax/=0.0 ) rpvgrw = MIN(amax/umax,rpvgrw)
      ENDDO
      CLA_GBRPVGRW = rpvgrw
      END FUNCTION CLA_GBRPVGRW
