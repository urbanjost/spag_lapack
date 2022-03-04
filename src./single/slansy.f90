!*==slansy.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SLANSY returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the element of largest absolute value of a real symmetric matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLANSY + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slansy.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slansy.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slansy.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       REAL             FUNCTION SLANSY( NORM, UPLO, N, A, LDA, WORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          NORM, UPLO
!       INTEGER            LDA, N
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLANSY  returns the value of the one norm,  or the Frobenius norm, or
!> the  infinity norm,  or the  element of  largest absolute value  of a
!> real symmetric matrix A.
!> \endverbatim
!>
!> \return SLANSY
!> \verbatim
!>
!>    SLANSY = ( max(abs(A(i,j))), NORM = 'M' or 'm'
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
!>          Specifies the value to be returned in SLANSY as described
!>          above.
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the upper or lower triangular part of the
!>          symmetric matrix A is to be referenced.
!>          = 'U':  Upper triangular part of A is referenced
!>          = 'L':  Lower triangular part of A is referenced
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.  When N = 0, SLANSY is
!>          set to zero.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>          The symmetric matrix A.  If UPLO = 'U', the leading n by n
!>          upper triangular part of A contains the upper triangular part
!>          of the matrix A, and the strictly lower triangular part of A
!>          is not referenced.  If UPLO = 'L', the leading n by n lower
!>          triangular part of A contains the lower triangular part of
!>          the matrix A, and the strictly upper triangular part of A is
!>          not referenced.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(N,1).
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
!> \ingroup realSYauxiliary
!
!  =====================================================================
      FUNCTION SLANSY(Norm,Uplo,N,A,Lda,Work)
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
      USE S_LSAME
      USE S_SCOMBSSQ
      USE S_SISNAN
      USE S_SLASSQ
      IMPLICIT NONE
!*--SLANSY136
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      REAL :: SLANSY
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Norm
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
!
! Local variable declarations rewritten by SPAG
!
      REAL :: absa , sum , value
      REAL , DIMENSION(2) :: colssq , ssq
      INTEGER :: i , j
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
            DO j = 1 , N
               DO i = 1 , j
                  sum = ABS(A(i,j))
                  IF ( value<sum .OR. SISNAN(sum) ) value = sum
               ENDDO
            ENDDO
         ELSE
            DO j = 1 , N
               DO i = j , N
                  sum = ABS(A(i,j))
                  IF ( value<sum .OR. SISNAN(sum) ) value = sum
               ENDDO
            ENDDO
         ENDIF
      ELSEIF ( (LSAME(Norm,'I')) .OR. (LSAME(Norm,'O')) .OR. (Norm=='1')&
     &         ) THEN
!
!        Find normI(A) ( = norm1(A), since A is symmetric).
!
         value = ZERO
         IF ( LSAME(Uplo,'U') ) THEN
            DO j = 1 , N
               sum = ZERO
               DO i = 1 , j - 1
                  absa = ABS(A(i,j))
                  sum = sum + absa
                  Work(i) = Work(i) + absa
               ENDDO
               Work(j) = sum + ABS(A(j,j))
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
               sum = Work(j) + ABS(A(j,j))
               DO i = j + 1 , N
                  absa = ABS(A(i,j))
                  sum = sum + absa
                  Work(i) = Work(i) + absa
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
         IF ( LSAME(Uplo,'U') ) THEN
            DO j = 2 , N
               colssq(1) = ZERO
               colssq(2) = ONE
               CALL SLASSQ(j-1,A(1,j),1,colssq(1),colssq(2))
               CALL SCOMBSSQ(ssq,colssq)
            ENDDO
         ELSE
            DO j = 1 , N - 1
               colssq(1) = ZERO
               colssq(2) = ONE
               CALL SLASSQ(N-j,A(j+1,j),1,colssq(1),colssq(2))
               CALL SCOMBSSQ(ssq,colssq)
            ENDDO
         ENDIF
         ssq(2) = 2*ssq(2)
!
!        Sum diagonal
!
         colssq(1) = ZERO
         colssq(2) = ONE
         CALL SLASSQ(N,A,Lda+1,colssq(1),colssq(2))
         CALL SCOMBSSQ(ssq,colssq)
         value = ssq(1)*SQRT(ssq(2))
      ENDIF
!
      SLANSY = value
!
!     End of SLANSY
!
      END FUNCTION SLANSY
