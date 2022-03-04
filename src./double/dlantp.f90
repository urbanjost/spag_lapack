!*==dlantp.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DLANTP returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the element of largest absolute value of a triangular matrix supplied in packed form.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLANTP + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlantp.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlantp.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlantp.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION DLANTP( NORM, UPLO, DIAG, N, AP, WORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIAG, NORM, UPLO
!       INTEGER            N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   AP( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLANTP  returns the value of the one norm,  or the Frobenius norm, or
!> the  infinity norm,  or the  element of  largest absolute value  of a
!> triangular matrix A, supplied in packed form.
!> \endverbatim
!>
!> \return DLANTP
!> \verbatim
!>
!>    DLANTP = ( max(abs(A(i,j))), NORM = 'M' or 'm'
!>             (
!>             ( norm1(A),         NORM = '1', 'O' or 'o'
!>             (
!>             ( normI(A),         NORM = 'I' or 'i'
!>             (
!>             ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
!>
!> where  norm1  denotes the  one norm of a matrix (maximum column sum),
!> normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
!> normF  denotes the  Frobenius norm of a matrix (square root of sum of
!> squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] NORM
!> \verbatim
!>          NORM is CHARACTER*1
!>          Specifies the value to be returned in DLANTP as described
!>          above.
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the matrix A is upper or lower triangular.
!>          = 'U':  Upper triangular
!>          = 'L':  Lower triangular
!> \endverbatim
!>
!> \param[in] DIAG
!> \verbatim
!>          DIAG is CHARACTER*1
!>          Specifies whether or not the matrix A is unit triangular.
!>          = 'N':  Non-unit triangular
!>          = 'U':  Unit triangular
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.  When N = 0, DLANTP is
!>          set to zero.
!> \endverbatim
!>
!> \param[in] AP
!> \verbatim
!>          AP is DOUBLE PRECISION array, dimension (N*(N+1)/2)
!>          The upper or lower triangular matrix A, packed columnwise in
!>          a linear array.  The j-th column of A is stored in the array
!>          AP as follows:
!>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
!>          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
!>          Note that when DIAG = 'U', the elements of the array AP
!>          corresponding to the diagonal elements of the matrix A are
!>          not referenced, but are assumed to be one.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)),
!>          where LWORK >= N when NORM = 'I'; otherwise, WORK is not
!>          referenced.
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
!  =====================================================================
      FUNCTION DLANTP(Norm,Uplo,Diag,N,Ap,Work)
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
      USE F77KINDS                        
      USE S_DCOMBSSQ
      USE S_DISNAN
      USE S_DLASSQ
      USE S_LSAME
      IMPLICIT NONE
!*--DLANTP139
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      REAL(R8KIND) :: DLANTP
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Norm
      CHARACTER :: Uplo
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , DIMENSION(*) :: Ap
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) , DIMENSION(2) :: colssq , ssq
      INTEGER :: i , j , k
      REAL(R8KIND) :: sum , value
      LOGICAL :: udiag
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
! =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. Local Arrays ..
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
      IF ( N==0 ) THEN
         value = ZERO
      ELSEIF ( LSAME(Norm,'M') ) THEN
!
!        Find max(abs(A(i,j))).
!
         k = 1
         IF ( LSAME(Diag,'U') ) THEN
            value = ONE
            IF ( LSAME(Uplo,'U') ) THEN
               DO j = 1 , N
                  DO i = k , k + j - 2
                     sum = ABS(Ap(i))
                     IF ( value<sum .OR. DISNAN(sum) ) value = sum
                  ENDDO
                  k = k + j
               ENDDO
            ELSE
               DO j = 1 , N
                  DO i = k + 1 , k + N - j
                     sum = ABS(Ap(i))
                     IF ( value<sum .OR. DISNAN(sum) ) value = sum
                  ENDDO
                  k = k + N - j + 1
               ENDDO
            ENDIF
         ELSE
            value = ZERO
            IF ( LSAME(Uplo,'U') ) THEN
               DO j = 1 , N
                  DO i = k , k + j - 1
                     sum = ABS(Ap(i))
                     IF ( value<sum .OR. DISNAN(sum) ) value = sum
                  ENDDO
                  k = k + j
               ENDDO
            ELSE
               DO j = 1 , N
                  DO i = k , k + N - j
                     sum = ABS(Ap(i))
                     IF ( value<sum .OR. DISNAN(sum) ) value = sum
                  ENDDO
                  k = k + N - j + 1
               ENDDO
            ENDIF
         ENDIF
      ELSEIF ( (LSAME(Norm,'O')) .OR. (Norm=='1') ) THEN
!
!        Find norm1(A).
!
         value = ZERO
         k = 1
         udiag = LSAME(Diag,'U')
         IF ( LSAME(Uplo,'U') ) THEN
            DO j = 1 , N
               IF ( udiag ) THEN
                  sum = ONE
                  DO i = k , k + j - 2
                     sum = sum + ABS(Ap(i))
                  ENDDO
               ELSE
                  sum = ZERO
                  DO i = k , k + j - 1
                     sum = sum + ABS(Ap(i))
                  ENDDO
               ENDIF
               k = k + j
               IF ( value<sum .OR. DISNAN(sum) ) value = sum
            ENDDO
         ELSE
            DO j = 1 , N
               IF ( udiag ) THEN
                  sum = ONE
                  DO i = k + 1 , k + N - j
                     sum = sum + ABS(Ap(i))
                  ENDDO
               ELSE
                  sum = ZERO
                  DO i = k , k + N - j
                     sum = sum + ABS(Ap(i))
                  ENDDO
               ENDIF
               k = k + N - j + 1
               IF ( value<sum .OR. DISNAN(sum) ) value = sum
            ENDDO
         ENDIF
      ELSEIF ( LSAME(Norm,'I') ) THEN
!
!        Find normI(A).
!
         k = 1
         IF ( LSAME(Uplo,'U') ) THEN
            IF ( LSAME(Diag,'U') ) THEN
               DO i = 1 , N
                  Work(i) = ONE
               ENDDO
               DO j = 1 , N
                  DO i = 1 , j - 1
                     Work(i) = Work(i) + ABS(Ap(k))
                     k = k + 1
                  ENDDO
                  k = k + 1
               ENDDO
            ELSE
               DO i = 1 , N
                  Work(i) = ZERO
               ENDDO
               DO j = 1 , N
                  DO i = 1 , j
                     Work(i) = Work(i) + ABS(Ap(k))
                     k = k + 1
                  ENDDO
               ENDDO
            ENDIF
         ELSEIF ( LSAME(Diag,'U') ) THEN
            DO i = 1 , N
               Work(i) = ONE
            ENDDO
            DO j = 1 , N
               k = k + 1
               DO i = j + 1 , N
                  Work(i) = Work(i) + ABS(Ap(k))
                  k = k + 1
               ENDDO
            ENDDO
         ELSE
            DO i = 1 , N
               Work(i) = ZERO
            ENDDO
            DO j = 1 , N
               DO i = j , N
                  Work(i) = Work(i) + ABS(Ap(k))
                  k = k + 1
               ENDDO
            ENDDO
         ENDIF
         value = ZERO
         DO i = 1 , N
            sum = Work(i)
            IF ( value<sum .OR. DISNAN(sum) ) value = sum
         ENDDO
      ELSEIF ( (LSAME(Norm,'F')) .OR. (LSAME(Norm,'E')) ) THEN
!
!        Find normF(A).
!        SSQ(1) is scale
!        SSQ(2) is sum-of-squares
!        For better accuracy, sum each column separately.
!
         IF ( LSAME(Uplo,'U') ) THEN
            IF ( LSAME(Diag,'U') ) THEN
               ssq(1) = ONE
               ssq(2) = N
               k = 2
               DO j = 2 , N
                  colssq(1) = ZERO
                  colssq(2) = ONE
                  CALL DLASSQ(j-1,Ap(k),1,colssq(1),colssq(2))
                  CALL DCOMBSSQ(ssq,colssq)
                  k = k + j
               ENDDO
            ELSE
               ssq(1) = ZERO
               ssq(2) = ONE
               k = 1
               DO j = 1 , N
                  colssq(1) = ZERO
                  colssq(2) = ONE
                  CALL DLASSQ(j,Ap(k),1,colssq(1),colssq(2))
                  CALL DCOMBSSQ(ssq,colssq)
                  k = k + j
               ENDDO
            ENDIF
         ELSEIF ( LSAME(Diag,'U') ) THEN
            ssq(1) = ONE
            ssq(2) = N
            k = 2
            DO j = 1 , N - 1
               colssq(1) = ZERO
               colssq(2) = ONE
               CALL DLASSQ(N-j,Ap(k),1,colssq(1),colssq(2))
               CALL DCOMBSSQ(ssq,colssq)
               k = k + N - j + 1
            ENDDO
         ELSE
            ssq(1) = ZERO
            ssq(2) = ONE
            k = 1
            DO j = 1 , N
               colssq(1) = ZERO
               colssq(2) = ONE
               CALL DLASSQ(N-j+1,Ap(k),1,colssq(1),colssq(2))
               CALL DCOMBSSQ(ssq,colssq)
               k = k + N - j + 1
            ENDDO
         ENDIF
         value = ssq(1)*SQRT(ssq(2))
      ENDIF
!
      DLANTP = value
!
!     End of DLANTP
!
      END FUNCTION DLANTP
