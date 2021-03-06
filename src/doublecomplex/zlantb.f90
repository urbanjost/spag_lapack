!*==zlantb.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b ZLANTB returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the element of largest absolute value of a triangular band matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLANTB + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlantb.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlantb.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlantb.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION ZLANTB( NORM, UPLO, DIAG, N, K, AB,
!                        LDAB, WORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIAG, NORM, UPLO
!       INTEGER            K, LDAB, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   WORK( * )
!       COMPLEX*16         AB( LDAB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLANTB  returns the value of the one norm,  or the Frobenius norm, or
!> the  infinity norm,  or the element of  largest absolute value  of an
!> n by n triangular band matrix A,  with ( k + 1 ) diagonals.
!> \endverbatim
!>
!> \return ZLANTB
!> \verbatim
!>
!>    ZLANTB = ( max(abs(A(i,j))), NORM = 'M' or 'm'
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
!>          Specifies the value to be returned in ZLANTB as described
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
!>          The order of the matrix A.  N >= 0.  When N = 0, ZLANTB is
!>          set to zero.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The number of super-diagonals of the matrix A if UPLO = 'U',
!>          or the number of sub-diagonals of the matrix A if UPLO = 'L'.
!>          K >= 0.
!> \endverbatim
!>
!> \param[in] AB
!> \verbatim
!>          AB is COMPLEX*16 array, dimension (LDAB,N)
!>          The upper or lower triangular band matrix A, stored in the
!>          first k+1 rows of AB.  The j-th column of A is stored
!>          in the j-th column of the array AB as follows:
!>          if UPLO = 'U', AB(k+1+i-j,j) = A(i,j) for max(1,j-k)<=i<=j;
!>          if UPLO = 'L', AB(1+i-j,j)   = A(i,j) for j<=i<=min(n,j+k).
!>          Note that when DIAG = 'U', the elements of the array AB
!>          corresponding to the diagonal elements of the matrix A are
!>          not referenced, but are assumed to be one.
!> \endverbatim
!>
!> \param[in] LDAB
!> \verbatim
!>          LDAB is INTEGER
!>          The leading dimension of the array AB.  LDAB >= K+1.
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
!> \ingroup complex16OTHERauxiliary
!
!  =====================================================================
      DOUBLE PRECISION FUNCTION ZLANTB(Norm,Uplo,Diag,N,K,Ab,Ldab,Work)
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
      IMPLICIT NONE
!*--ZLANTB150
!     .. Scalar Arguments ..
      CHARACTER Diag , Norm , Uplo
      INTEGER K , Ldab , N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION Work(*)
      COMPLEX*16 Ab(Ldab,*)
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
      INTEGER i , j , l
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
      INTRINSIC ABS , MAX , MIN , SQRT
!     ..
!     .. Executable Statements ..
!
      IF ( N==0 ) THEN
         value = ZERO
      ELSEIF ( LSAME(Norm,'M') ) THEN
!
!        Find max(abs(A(i,j))).
!
         IF ( LSAME(Diag,'U') ) THEN
            value = ONE
            IF ( LSAME(Uplo,'U') ) THEN
               DO j = 1 , N
                  DO i = MAX(K+2-j,1) , K
                     sum = ABS(Ab(i,j))
                     IF ( value<sum .OR. DISNAN(sum) ) value = sum
                  ENDDO
               ENDDO
            ELSE
               DO j = 1 , N
                  DO i = 2 , MIN(N+1-j,K+1)
                     sum = ABS(Ab(i,j))
                     IF ( value<sum .OR. DISNAN(sum) ) value = sum
                  ENDDO
               ENDDO
            ENDIF
         ELSE
            value = ZERO
            IF ( LSAME(Uplo,'U') ) THEN
               DO j = 1 , N
                  DO i = MAX(K+2-j,1) , K + 1
                     sum = ABS(Ab(i,j))
                     IF ( value<sum .OR. DISNAN(sum) ) value = sum
                  ENDDO
               ENDDO
            ELSE
               DO j = 1 , N
                  DO i = 1 , MIN(N+1-j,K+1)
                     sum = ABS(Ab(i,j))
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
               IF ( udiag ) THEN
                  sum = ONE
                  DO i = MAX(K+2-j,1) , K
                     sum = sum + ABS(Ab(i,j))
                  ENDDO
               ELSE
                  sum = ZERO
                  DO i = MAX(K+2-j,1) , K + 1
                     sum = sum + ABS(Ab(i,j))
                  ENDDO
               ENDIF
               IF ( value<sum .OR. DISNAN(sum) ) value = sum
            ENDDO
         ELSE
            DO j = 1 , N
               IF ( udiag ) THEN
                  sum = ONE
                  DO i = 2 , MIN(N+1-j,K+1)
                     sum = sum + ABS(Ab(i,j))
                  ENDDO
               ELSE
                  sum = ZERO
                  DO i = 1 , MIN(N+1-j,K+1)
                     sum = sum + ABS(Ab(i,j))
                  ENDDO
               ENDIF
               IF ( value<sum .OR. DISNAN(sum) ) value = sum
            ENDDO
         ENDIF
      ELSEIF ( LSAME(Norm,'I') ) THEN
!
!        Find normI(A).
!
         value = ZERO
         IF ( LSAME(Uplo,'U') ) THEN
            IF ( LSAME(Diag,'U') ) THEN
               DO i = 1 , N
                  Work(i) = ONE
               ENDDO
               DO j = 1 , N
                  l = K + 1 - j
                  DO i = MAX(1,j-K) , j - 1
                     Work(i) = Work(i) + ABS(Ab(l+i,j))
                  ENDDO
               ENDDO
            ELSE
               DO i = 1 , N
                  Work(i) = ZERO
               ENDDO
               DO j = 1 , N
                  l = K + 1 - j
                  DO i = MAX(1,j-K) , j
                     Work(i) = Work(i) + ABS(Ab(l+i,j))
                  ENDDO
               ENDDO
            ENDIF
         ELSEIF ( LSAME(Diag,'U') ) THEN
            DO i = 1 , N
               Work(i) = ONE
            ENDDO
            DO j = 1 , N
               l = 1 - j
               DO i = j + 1 , MIN(N,j+K)
                  Work(i) = Work(i) + ABS(Ab(l+i,j))
               ENDDO
            ENDDO
         ELSE
            DO i = 1 , N
               Work(i) = ZERO
            ENDDO
            DO j = 1 , N
               l = 1 - j
               DO i = j , MIN(N,j+K)
                  Work(i) = Work(i) + ABS(Ab(l+i,j))
               ENDDO
            ENDDO
         ENDIF
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
               IF ( K>0 ) THEN
                  DO j = 2 , N
                     colssq(1) = ZERO
                     colssq(2) = ONE
                     CALL ZLASSQ(MIN(j-1,K),Ab(MAX(K+2-j,1),j),1,       &
     &                           colssq(1),colssq(2))
                     CALL DCOMBSSQ(ssq,colssq)
                  ENDDO
               ENDIF
            ELSE
               ssq(1) = ZERO
               ssq(2) = ONE
               DO j = 1 , N
                  colssq(1) = ZERO
                  colssq(2) = ONE
                  CALL ZLASSQ(MIN(j,K+1),Ab(MAX(K+2-j,1),j),1,colssq(1),&
     &                        colssq(2))
                  CALL DCOMBSSQ(ssq,colssq)
               ENDDO
            ENDIF
         ELSEIF ( LSAME(Diag,'U') ) THEN
            ssq(1) = ONE
            ssq(2) = N
            IF ( K>0 ) THEN
               DO j = 1 , N - 1
                  colssq(1) = ZERO
                  colssq(2) = ONE
                  CALL ZLASSQ(MIN(N-j,K),Ab(2,j),1,colssq(1),colssq(2))
                  CALL DCOMBSSQ(ssq,colssq)
               ENDDO
            ENDIF
         ELSE
            ssq(1) = ZERO
            ssq(2) = ONE
            DO j = 1 , N
               colssq(1) = ZERO
               colssq(2) = ONE
               CALL ZLASSQ(MIN(N-j+1,K+1),Ab(1,j),1,colssq(1),colssq(2))
               CALL DCOMBSSQ(ssq,colssq)
            ENDDO
         ENDIF
         value = ssq(1)*SQRT(ssq(2))
      ENDIF
!
      ZLANTB = value
!
!     End of ZLANTB
!
      END FUNCTION ZLANTB
