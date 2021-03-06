!*==zlantr.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b ZLANTR returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the element of largest absolute value of a trapezoidal or triangular matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLANTR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlantr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlantr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlantr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION ZLANTR( NORM, UPLO, DIAG, M, N, A, LDA,
!                        WORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIAG, NORM, UPLO
!       INTEGER            LDA, M, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   WORK( * )
!       COMPLEX*16         A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLANTR  returns the value of the one norm,  or the Frobenius norm, or
!> the  infinity norm,  or the  element of  largest absolute value  of a
!> trapezoidal or triangular matrix A.
!> \endverbatim
!>
!> \return ZLANTR
!> \verbatim
!>
!>    ZLANTR = ( max(abs(A(i,j))), NORM = 'M' or 'm'
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
!>          Specifies the value to be returned in ZLANTR as described
!>          above.
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the matrix A is upper or lower trapezoidal.
!>          = 'U':  Upper trapezoidal
!>          = 'L':  Lower trapezoidal
!>          Note that A is triangular instead of trapezoidal if M = N.
!> \endverbatim
!>
!> \param[in] DIAG
!> \verbatim
!>          DIAG is CHARACTER*1
!>          Specifies whether or not the matrix A has unit diagonal.
!>          = 'N':  Non-unit diagonal
!>          = 'U':  Unit diagonal
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= 0, and if
!>          UPLO = 'U', M <= N.  When M = 0, ZLANTR is set to zero.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  N >= 0, and if
!>          UPLO = 'L', N <= M.  When N = 0, ZLANTR is set to zero.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          The trapezoidal matrix A (A is triangular if M = N).
!>          If UPLO = 'U', the leading m by n upper trapezoidal part of
!>          the array A contains the upper trapezoidal matrix, and the
!>          strictly lower triangular part of A is not referenced.
!>          If UPLO = 'L', the leading m by n lower trapezoidal part of
!>          the array A contains the lower trapezoidal matrix, and the
!>          strictly upper triangular part of A is not referenced.  Note
!>          that when DIAG = 'U', the diagonal elements of A are not
!>          referenced and are assumed to be one.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(M,1).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)),
!>          where LWORK >= M when NORM = 'I'; otherwise, WORK is not
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
!> \ingroup complex16OTHERauxiliary
!
!  =====================================================================
      DOUBLE PRECISION FUNCTION ZLANTR(Norm,Uplo,Diag,M,N,A,Lda,Work)
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
      IMPLICIT NONE
!*--ZLANTR151
!     .. Scalar Arguments ..
      CHARACTER Diag , Norm , Uplo
      INTEGER Lda , M , N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION Work(*)
      COMPLEX*16 A(Lda,*)
!     ..
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ONE , ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
      LOGICAL udiag
      INTEGER i , j
      DOUBLE PRECISION sum , value
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION ssq(2) , colssq(2)
!     ..
!     .. External Functions ..
      LOGICAL LSAME , DISNAN
      EXTERNAL LSAME , DISNAN
!     ..
!     .. External Subroutines ..
      EXTERNAL ZLASSQ , DCOMBSSQ
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MIN , SQRT
!     ..
!     .. Executable Statements ..
!
      IF ( MIN(M,N)==0 ) THEN
         value = ZERO
      ELSEIF ( LSAME(Norm,'M') ) THEN
!
!        Find max(abs(A(i,j))).
!
         IF ( LSAME(Diag,'U') ) THEN
            value = ONE
            IF ( LSAME(Uplo,'U') ) THEN
               DO j = 1 , N
                  DO i = 1 , MIN(M,j-1)
                     sum = ABS(A(i,j))
                     IF ( value<sum .OR. DISNAN(sum) ) value = sum
                  ENDDO
               ENDDO
            ELSE
               DO j = 1 , N
                  DO i = j + 1 , M
                     sum = ABS(A(i,j))
                     IF ( value<sum .OR. DISNAN(sum) ) value = sum
                  ENDDO
               ENDDO
            ENDIF
         ELSE
            value = ZERO
            IF ( LSAME(Uplo,'U') ) THEN
               DO j = 1 , N
                  DO i = 1 , MIN(M,j)
                     sum = ABS(A(i,j))
                     IF ( value<sum .OR. DISNAN(sum) ) value = sum
                  ENDDO
               ENDDO
            ELSE
               DO j = 1 , N
                  DO i = j , M
                     sum = ABS(A(i,j))
                     IF ( value<sum .OR. DISNAN(sum) ) value = sum
                  ENDDO
               ENDDO
            ENDIF
         ENDIF
      ELSEIF ( (LSAME(Norm,'O')) .OR. (Norm=='1') ) THEN
!
!        Find norm1(A).
!
         value = ZERO
         udiag = LSAME(Diag,'U')
         IF ( LSAME(Uplo,'U') ) THEN
            DO j = 1 , N
               IF ( (udiag) .AND. (j<=M) ) THEN
                  sum = ONE
                  DO i = 1 , j - 1
                     sum = sum + ABS(A(i,j))
                  ENDDO
               ELSE
                  sum = ZERO
                  DO i = 1 , MIN(M,j)
                     sum = sum + ABS(A(i,j))
                  ENDDO
               ENDIF
               IF ( value<sum .OR. DISNAN(sum) ) value = sum
            ENDDO
         ELSE
            DO j = 1 , N
               IF ( udiag ) THEN
                  sum = ONE
                  DO i = j + 1 , M
                     sum = sum + ABS(A(i,j))
                  ENDDO
               ELSE
                  sum = ZERO
                  DO i = j , M
                     sum = sum + ABS(A(i,j))
                  ENDDO
               ENDIF
               IF ( value<sum .OR. DISNAN(sum) ) value = sum
            ENDDO
         ENDIF
      ELSEIF ( LSAME(Norm,'I') ) THEN
!
!        Find normI(A).
!
         IF ( LSAME(Uplo,'U') ) THEN
            IF ( LSAME(Diag,'U') ) THEN
               DO i = 1 , M
                  Work(i) = ONE
               ENDDO
               DO j = 1 , N
                  DO i = 1 , MIN(M,j-1)
                     Work(i) = Work(i) + ABS(A(i,j))
                  ENDDO
               ENDDO
            ELSE
               DO i = 1 , M
                  Work(i) = ZERO
               ENDDO
               DO j = 1 , N
                  DO i = 1 , MIN(M,j)
                     Work(i) = Work(i) + ABS(A(i,j))
                  ENDDO
               ENDDO
            ENDIF
         ELSEIF ( LSAME(Diag,'U') ) THEN
            DO i = 1 , MIN(M,N)
               Work(i) = ONE
            ENDDO
            DO i = N + 1 , M
               Work(i) = ZERO
            ENDDO
            DO j = 1 , N
               DO i = j + 1 , M
                  Work(i) = Work(i) + ABS(A(i,j))
               ENDDO
            ENDDO
         ELSE
            DO i = 1 , M
               Work(i) = ZERO
            ENDDO
            DO j = 1 , N
               DO i = j , M
                  Work(i) = Work(i) + ABS(A(i,j))
               ENDDO
            ENDDO
         ENDIF
         value = ZERO
         DO i = 1 , M
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
               ssq(2) = MIN(M,N)
               DO j = 2 , N
                  colssq(1) = ZERO
                  colssq(2) = ONE
                  CALL ZLASSQ(MIN(M,j-1),A(1,j),1,colssq(1),colssq(2))
                  CALL DCOMBSSQ(ssq,colssq)
               ENDDO
            ELSE
               ssq(1) = ZERO
               ssq(2) = ONE
               DO j = 1 , N
                  colssq(1) = ZERO
                  colssq(2) = ONE
                  CALL ZLASSQ(MIN(M,j),A(1,j),1,colssq(1),colssq(2))
                  CALL DCOMBSSQ(ssq,colssq)
               ENDDO
            ENDIF
         ELSEIF ( LSAME(Diag,'U') ) THEN
            ssq(1) = ONE
            ssq(2) = MIN(M,N)
            DO j = 1 , N
               colssq(1) = ZERO
               colssq(2) = ONE
               CALL ZLASSQ(M-j,A(MIN(M,j+1),j),1,colssq(1),colssq(2))
               CALL DCOMBSSQ(ssq,colssq)
            ENDDO
         ELSE
            ssq(1) = ZERO
            ssq(2) = ONE
            DO j = 1 , N
               colssq(1) = ZERO
               colssq(2) = ONE
               CALL ZLASSQ(M-j+1,A(j,j),1,colssq(1),colssq(2))
               CALL DCOMBSSQ(ssq,colssq)
            ENDDO
         ENDIF
         value = ssq(1)*SQRT(ssq(2))
      ENDIF
!
      ZLANTR = value
!
!     End of ZLANTR
!
      END FUNCTION ZLANTR
