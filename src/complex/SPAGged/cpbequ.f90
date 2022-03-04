!*==cpbequ.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CPBEQU
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CPBEQU + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpbequ.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpbequ.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpbequ.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CPBEQU( UPLO, N, KD, AB, LDAB, S, SCOND, AMAX, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, KD, LDAB, N
!       REAL               AMAX, SCOND
!       ..
!       .. Array Arguments ..
!       REAL               S( * )
!       COMPLEX            AB( LDAB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CPBEQU computes row and column scalings intended to equilibrate a
!> Hermitian positive definite band matrix A and reduce its condition
!> number (with respect to the two-norm).  S contains the scale factors,
!> S(i) = 1/sqrt(A(i,i)), chosen so that the scaled matrix B with
!> elements B(i,j) = S(i)*A(i,j)*S(j) has ones on the diagonal.  This
!> choice of S puts the condition number of B within a factor N of the
!> smallest possible condition number over all possible diagonal
!> scalings.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  Upper triangular of A is stored;
!>          = 'L':  Lower triangular of A is stored.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] KD
!> \verbatim
!>          KD is INTEGER
!>          The number of superdiagonals of the matrix A if UPLO = 'U',
!>          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
!> \endverbatim
!>
!> \param[in] AB
!> \verbatim
!>          AB is COMPLEX array, dimension (LDAB,N)
!>          The upper or lower triangle of the Hermitian band matrix A,
!>          stored in the first KD+1 rows of the array.  The j-th column
!>          of A is stored in the j-th column of the array AB as follows:
!>          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
!>          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
!> \endverbatim
!>
!> \param[in] LDAB
!> \verbatim
!>          LDAB is INTEGER
!>          The leading dimension of the array A.  LDAB >= KD+1.
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is REAL array, dimension (N)
!>          If INFO = 0, S contains the scale factors for A.
!> \endverbatim
!>
!> \param[out] SCOND
!> \verbatim
!>          SCOND is REAL
!>          If INFO = 0, S contains the ratio of the smallest S(i) to
!>          the largest S(i).  If SCOND >= 0.1 and AMAX is neither too
!>          large nor too small, it is not worth scaling by S.
!> \endverbatim
!>
!> \param[out] AMAX
!> \verbatim
!>          AMAX is REAL
!>          Absolute value of largest matrix element.  If AMAX is very
!>          close to overflow or very close to underflow, the matrix
!>          should be scaled.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          > 0:  if INFO = i, the i-th diagonal element is nonpositive.
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
!> \ingroup complexOTHERcomputational
!
!  =====================================================================
      SUBROUTINE CPBEQU(Uplo,N,Kd,Ab,Ldab,S,Scond,Amax,Info)
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--CPBEQU136
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Kd
      COMPLEX , INTENT(IN) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL , INTENT(INOUT) , DIMENSION(*) :: S
      REAL , INTENT(OUT) :: Scond
      REAL , INTENT(INOUT) :: Amax
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , j
      REAL :: smin
      LOGICAL :: upper
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
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
      upper = LSAME(Uplo,'U')
      IF ( .NOT.upper .AND. .NOT.LSAME(Uplo,'L') ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Kd<0 ) THEN
         Info = -3
      ELSEIF ( Ldab<Kd+1 ) THEN
         Info = -5
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('CPBEQU',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) THEN
         Scond = ONE
         Amax = ZERO
         RETURN
      ENDIF
!
      IF ( upper ) THEN
         j = Kd + 1
      ELSE
         j = 1
      ENDIF
!
!     Initialize SMIN and AMAX.
!
      S(1) = REAL(Ab(j,1))
      smin = S(1)
      Amax = S(1)
!
!     Find the minimum and maximum diagonal elements.
!
      DO i = 2 , N
         S(i) = REAL(Ab(j,i))
         smin = MIN(smin,S(i))
         Amax = MAX(Amax,S(i))
      ENDDO
!
      IF ( smin<=ZERO ) THEN
!
!        Find the first non-positive diagonal element and return.
!
         DO i = 1 , N
            IF ( S(i)<=ZERO ) THEN
               Info = i
               RETURN
            ENDIF
         ENDDO
      ELSE
!
!        Set the scale factors to the reciprocals
!        of the diagonal elements.
!
         DO i = 1 , N
            S(i) = ONE/SQRT(S(i))
         ENDDO
!
!        Compute SCOND = min(S(I)) / max(S(I))
!
         Scond = SQRT(smin)/SQRT(Amax)
      ENDIF
!
!     End of CPBEQU
!
      END SUBROUTINE CPBEQU
