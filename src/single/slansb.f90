!*==slansb.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SLANSB returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the element of largest absolute value of a symmetric band matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLANSB + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slansb.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slansb.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slansb.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       REAL             FUNCTION SLANSB( NORM, UPLO, N, K, AB, LDAB,
!                        WORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          NORM, UPLO
!       INTEGER            K, LDAB, N
!       ..
!       .. Array Arguments ..
!       REAL               AB( LDAB, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLANSB  returns the value of the one norm,  or the Frobenius norm, or
!> the  infinity norm,  or the element of  largest absolute value  of an
!> n by n symmetric band matrix A,  with k super-diagonals.
!> \endverbatim
!>
!> \return SLANSB
!> \verbatim
!>
!>    SLANSB = ( max(abs(A(i,j))), NORM = 'M' or 'm'
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
!>          Specifies the value to be returned in SLANSB as described
!>          above.
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the upper or lower triangular part of the
!>          band matrix A is supplied.
!>          = 'U':  Upper triangular part is supplied
!>          = 'L':  Lower triangular part is supplied
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.  When N = 0, SLANSB is
!>          set to zero.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The number of super-diagonals or sub-diagonals of the
!>          band matrix A.  K >= 0.
!> \endverbatim
!>
!> \param[in] AB
!> \verbatim
!>          AB is REAL array, dimension (LDAB,N)
!>          The upper or lower triangle of the symmetric band matrix A,
!>          stored in the first K+1 rows of AB.  The j-th column of A is
!>          stored in the j-th column of the array AB as follows:
!>          if UPLO = 'U', AB(k+1+i-j,j) = A(i,j) for max(1,j-k)<=i<=j;
!>          if UPLO = 'L', AB(1+i-j,j)   = A(i,j) for j<=i<=min(n,j+k).
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
!> \ingroup realOTHERauxiliary
!
!  =====================================================================
      REAL FUNCTION SLANSB(Norm,Uplo,N,K,Ab,Ldab,Work)
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
      IMPLICIT NONE
!*--SLANSB138
!     .. Scalar Arguments ..
      CHARACTER Norm , Uplo
      INTEGER K , Ldab , N
!     ..
!     .. Array Arguments ..
      REAL Ab(Ldab,*) , Work(*)
!     ..
!
! =====================================================================
!
!     .. Parameters ..
      REAL ONE , ZERO
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
!     ..
!     .. Local Scalars ..
      INTEGER i , j , l
      REAL absa , sum , value
!     ..
!     .. Local Arrays ..
      REAL ssq(2) , colssq(2)
!     ..
!     .. External Functions ..
      LOGICAL LSAME , SISNAN
      EXTERNAL LSAME , SISNAN
!     ..
!     .. External Subroutines ..
      EXTERNAL SLASSQ , SCOMBSSQ
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
         value = ZERO
         IF ( LSAME(Uplo,'U') ) THEN
            DO j = 1 , N
               DO i = MAX(K+2-j,1) , K + 1
                  sum = ABS(Ab(i,j))
                  IF ( value<sum .OR. SISNAN(sum) ) value = sum
               ENDDO
            ENDDO
         ELSE
            DO j = 1 , N
               DO i = 1 , MIN(N+1-j,K+1)
                  sum = ABS(Ab(i,j))
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
               l = K + 1 - j
               DO i = MAX(1,j-K) , j - 1
                  absa = ABS(Ab(l+i,j))
                  sum = sum + absa
                  Work(i) = Work(i) + absa
               ENDDO
               Work(j) = sum + ABS(Ab(K+1,j))
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
               sum = Work(j) + ABS(Ab(1,j))
               l = 1 - j
               DO i = j + 1 , MIN(N,j+K)
                  absa = ABS(Ab(l+i,j))
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
         IF ( K>0 ) THEN
            IF ( LSAME(Uplo,'U') ) THEN
               DO j = 2 , N
                  colssq(1) = ZERO
                  colssq(2) = ONE
                  CALL SLASSQ(MIN(j-1,K),Ab(MAX(K+2-j,1),j),1,colssq(1),&
     &                        colssq(2))
                  CALL SCOMBSSQ(ssq,colssq)
               ENDDO
               l = K + 1
            ELSE
               DO j = 1 , N - 1
                  colssq(1) = ZERO
                  colssq(2) = ONE
                  CALL SLASSQ(MIN(N-j,K),Ab(2,j),1,colssq(1),colssq(2))
                  CALL SCOMBSSQ(ssq,colssq)
               ENDDO
               l = 1
            ENDIF
            ssq(2) = 2*ssq(2)
         ELSE
            l = 1
         ENDIF
!
!        Sum diagonal
!
         colssq(1) = ZERO
         colssq(2) = ONE
         CALL SLASSQ(N,Ab(l,1),Ldab,colssq(1),colssq(2))
         CALL SCOMBSSQ(ssq,colssq)
         value = ssq(1)*SQRT(ssq(2))
      ENDIF
!
      SLANSB = value
!
!     End of SLANSB
!
      END FUNCTION SLANSB
