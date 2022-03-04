!*==clansp.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CLANSP returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the element of largest absolute value of a symmetric matrix supplied in packed form.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLANSP + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clansp.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clansp.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clansp.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       REAL             FUNCTION CLANSP( NORM, UPLO, N, AP, WORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          NORM, UPLO
!       INTEGER            N
!       ..
!       .. Array Arguments ..
!       REAL               WORK( * )
!       COMPLEX            AP( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLANSP  returns the value of the one norm,  or the Frobenius norm, or
!> the  infinity norm,  or the  element of  largest absolute value  of a
!> complex symmetric matrix A,  supplied in packed form.
!> \endverbatim
!>
!> \return CLANSP
!> \verbatim
!>
!>    CLANSP = ( max(abs(A(i,j))), NORM = 'M' or 'm'
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
!>          Specifies the value to be returned in CLANSP as described
!>          above.
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the upper or lower triangular part of the
!>          symmetric matrix A is supplied.
!>          = 'U':  Upper triangular part of A is supplied
!>          = 'L':  Lower triangular part of A is supplied
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.  When N = 0, CLANSP is
!>          set to zero.
!> \endverbatim
!>
!> \param[in] AP
!> \verbatim
!>          AP is COMPLEX array, dimension (N*(N+1)/2)
!>          The upper or lower triangle of the symmetric matrix A, packed
!>          columnwise in a linear array.  The j-th column of A is stored
!>          in the array AP as follows:
!>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
!>          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (MAX(1,LWORK)),
!>          where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,
!>          WORK is not referenced.
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
!> \ingroup complexOTHERauxiliary
!
!  =====================================================================
      FUNCTION CLANSP(Norm,Uplo,N,Ap,Work)
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
      USE S_CLASSQ
      USE S_LSAME
      USE S_SCOMBSSQ
      USE S_SISNAN
      IMPLICIT NONE
!*--CLANSP129
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      REAL :: CLANSP
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Norm
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX , DIMENSION(*) :: Ap
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
!
! Local variable declarations rewritten by SPAG
!
      REAL :: absa , sum , value
      REAL , DIMENSION(2) :: colssq , ssq
      INTEGER :: i , j , k
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
         value = ZERO
         IF ( LSAME(Uplo,'U') ) THEN
            k = 1
            DO j = 1 , N
               DO i = k , k + j - 1
                  sum = ABS(Ap(i))
                  IF ( value<sum .OR. SISNAN(sum) ) value = sum
               ENDDO
               k = k + j
            ENDDO
         ELSE
            k = 1
            DO j = 1 , N
               DO i = k , k + N - j
                  sum = ABS(Ap(i))
                  IF ( value<sum .OR. SISNAN(sum) ) value = sum
               ENDDO
               k = k + N - j + 1
            ENDDO
         ENDIF
      ELSEIF ( (LSAME(Norm,'I')) .OR. (LSAME(Norm,'O')) .OR. (Norm=='1')&
     &         ) THEN
!
!        Find normI(A) ( = norm1(A), since A is symmetric).
!
         value = ZERO
         k = 1
         IF ( LSAME(Uplo,'U') ) THEN
            DO j = 1 , N
               sum = ZERO
               DO i = 1 , j - 1
                  absa = ABS(Ap(k))
                  sum = sum + absa
                  Work(i) = Work(i) + absa
                  k = k + 1
               ENDDO
               Work(j) = sum + ABS(Ap(k))
               k = k + 1
            ENDDO
            DO i = 1 , N
               sum = Work(i)
               IF ( value<sum .OR. SISNAN(sum) ) value = sum
            ENDDO
         ELSE
            DO i = 1 , N
               Work(i) = ZERO
            ENDDO
            DO j = 1 , N
               sum = Work(j) + ABS(Ap(k))
               k = k + 1
               DO i = j + 1 , N
                  absa = ABS(Ap(k))
                  sum = sum + absa
                  Work(i) = Work(i) + absa
                  k = k + 1
               ENDDO
               IF ( value<sum .OR. SISNAN(sum) ) value = sum
            ENDDO
         ENDIF
      ELSEIF ( (LSAME(Norm,'F')) .OR. (LSAME(Norm,'E')) ) THEN
!
!        Find normF(A).
!        SSQ(1) is scale
!        SSQ(2) is sum-of-squares
!        For better accuracy, sum each column separately.
!
         ssq(1) = ZERO
         ssq(2) = ONE
!
!        Sum off-diagonals
!
         k = 2
         IF ( LSAME(Uplo,'U') ) THEN
            DO j = 2 , N
               colssq(1) = ZERO
               colssq(2) = ONE
               CALL CLASSQ(j-1,Ap(k),1,colssq(1),colssq(2))
               CALL SCOMBSSQ(ssq,colssq)
               k = k + j
            ENDDO
         ELSE
            DO j = 1 , N - 1
               colssq(1) = ZERO
               colssq(2) = ONE
               CALL CLASSQ(N-j,Ap(k),1,colssq(1),colssq(2))
               CALL SCOMBSSQ(ssq,colssq)
               k = k + N - j + 1
            ENDDO
         ENDIF
         ssq(2) = 2*ssq(2)
!
!        Sum diagonal
!
         k = 1
         colssq(1) = ZERO
         colssq(2) = ONE
         DO i = 1 , N
            IF ( REAL(Ap(k))/=ZERO ) THEN
               absa = ABS(REAL(Ap(k)))
               IF ( colssq(1)<absa ) THEN
                  colssq(2) = ONE + colssq(2)*(colssq(1)/absa)**2
                  colssq(1) = absa
               ELSE
                  colssq(2) = colssq(2) + (absa/colssq(1))**2
               ENDIF
            ENDIF
            IF ( AIMAG(Ap(k))/=ZERO ) THEN
               absa = ABS(AIMAG(Ap(k)))
               IF ( colssq(1)<absa ) THEN
                  colssq(2) = ONE + colssq(2)*(colssq(1)/absa)**2
                  colssq(1) = absa
               ELSE
                  colssq(2) = colssq(2) + (absa/colssq(1))**2
               ENDIF
            ENDIF
            IF ( LSAME(Uplo,'U') ) THEN
               k = k + i + 1
            ELSE
               k = k + N - i + 1
            ENDIF
         ENDDO
         CALL SCOMBSSQ(ssq,colssq)
         value = ssq(1)*SQRT(ssq(2))
      ENDIF
!
      CLANSP = value
!
!     End of CLANSP
!
      END FUNCTION CLANSP
