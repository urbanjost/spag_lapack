!*==cpttrf.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b CPTTRF
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CPTTRF + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpttrf.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpttrf.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpttrf.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CPTTRF( N, D, E, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, N
!       ..
!       .. Array Arguments ..
!       REAL               D( * )
!       COMPLEX            E( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CPTTRF computes the L*D*L**H factorization of a complex Hermitian
!> positive definite tridiagonal matrix A.  The factorization may also
!> be regarded as having the form A = U**H *D*U.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is REAL array, dimension (N)
!>          On entry, the n diagonal elements of the tridiagonal matrix
!>          A.  On exit, the n diagonal elements of the diagonal matrix
!>          D from the L*D*L**H factorization of A.
!> \endverbatim
!>
!> \param[in,out] E
!> \verbatim
!>          E is COMPLEX array, dimension (N-1)
!>          On entry, the (n-1) subdiagonal elements of the tridiagonal
!>          matrix A.  On exit, the (n-1) subdiagonal elements of the
!>          unit bidiagonal factor L from the L*D*L**H factorization of A.
!>          E can also be regarded as the superdiagonal of the unit
!>          bidiagonal factor U from the U**H *D*U factorization of A.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -k, the k-th argument had an illegal value
!>          > 0: if INFO = k, the leading minor of order k is not
!>               positive definite; if k < N, the factorization could not
!>               be completed, while if k = N, the factorization was
!>               completed, but D(N) <= 0.
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
!> \ingroup complexPTcomputational
!
!  =====================================================================
      SUBROUTINE CPTTRF(N,D,E,Info)
      IMPLICIT NONE
!*--CPTTRF96
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Info , N
!     ..
!     .. Array Arguments ..
      REAL D(*)
      COMPLEX E(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO
      PARAMETER (ZERO=0.0E+0)
!     ..
!     .. Local Scalars ..
      INTEGER i , i4
      REAL eii , eir , f , g
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC AIMAG , CMPLX , MOD , REAL
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
      IF ( N<0 ) THEN
         Info = -1
         CALL XERBLA('CPTTRF',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
!     Compute the L*D*L**H (or U**H *D*U) factorization of A.
!
      i4 = MOD(N-1,4)
      DO i = 1 , i4
         IF ( D(i)<=ZERO ) THEN
            Info = i
            GOTO 99999
         ENDIF
         eir = REAL(E(i))
         eii = AIMAG(E(i))
         f = eir/D(i)
         g = eii/D(i)
         E(i) = CMPLX(f,g)
         D(i+1) = D(i+1) - f*eir - g*eii
      ENDDO
!
      DO i = i4 + 1 , N - 4 , 4
!
!        Drop out of the loop if d(i) <= 0: the matrix is not positive
!        definite.
!
         IF ( D(i)<=ZERO ) THEN
            Info = i
            GOTO 99999
         ENDIF
!
!        Solve for e(i) and d(i+1).
!
         eir = REAL(E(i))
         eii = AIMAG(E(i))
         f = eir/D(i)
         g = eii/D(i)
         E(i) = CMPLX(f,g)
         D(i+1) = D(i+1) - f*eir - g*eii
!
         IF ( D(i+1)<=ZERO ) THEN
            Info = i + 1
            GOTO 99999
         ENDIF
!
!        Solve for e(i+1) and d(i+2).
!
         eir = REAL(E(i+1))
         eii = AIMAG(E(i+1))
         f = eir/D(i+1)
         g = eii/D(i+1)
         E(i+1) = CMPLX(f,g)
         D(i+2) = D(i+2) - f*eir - g*eii
!
         IF ( D(i+2)<=ZERO ) THEN
            Info = i + 2
            GOTO 99999
         ENDIF
!
!        Solve for e(i+2) and d(i+3).
!
         eir = REAL(E(i+2))
         eii = AIMAG(E(i+2))
         f = eir/D(i+2)
         g = eii/D(i+2)
         E(i+2) = CMPLX(f,g)
         D(i+3) = D(i+3) - f*eir - g*eii
!
         IF ( D(i+3)<=ZERO ) THEN
            Info = i + 3
            GOTO 99999
         ENDIF
!
!        Solve for e(i+3) and d(i+4).
!
         eir = REAL(E(i+3))
         eii = AIMAG(E(i+3))
         f = eir/D(i+3)
         g = eii/D(i+3)
         E(i+3) = CMPLX(f,g)
         D(i+4) = D(i+4) - f*eir - g*eii
      ENDDO
!
!     Check d(n) for positive definiteness.
!
      IF ( D(N)<=ZERO ) Info = N
!
!
!     End of CPTTRF
!
99999 END SUBROUTINE CPTTRF
