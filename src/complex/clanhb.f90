!*==clanhb.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b CLANHB returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the element of largest absolute value of a Hermitian band matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLANHB + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clanhb.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clanhb.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clanhb.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       REAL             FUNCTION CLANHB( NORM, UPLO, N, K, AB, LDAB,
!                        WORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          NORM, UPLO
!       INTEGER            K, LDAB, N
!       ..
!       .. Array Arguments ..
!       REAL               WORK( * )
!       COMPLEX            AB( LDAB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLANHB  returns the value of the one norm,  or the Frobenius norm, or
!> the  infinity norm,  or the element of  largest absolute value  of an
!> n by n hermitian band matrix A,  with k super-diagonals.
!> \endverbatim
!>
!> \return CLANHB
!> \verbatim
!>
!>    CLANHB = ( max(abs(A(i,j))), NORM = 'M' or 'm'
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
!>          Specifies the value to be returned in CLANHB as described
!>          above.
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the upper or lower triangular part of the
!>          band matrix A is supplied.
!>          = 'U':  Upper triangular
!>          = 'L':  Lower triangular
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.  When N = 0, CLANHB is
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
!>          AB is COMPLEX array, dimension (LDAB,N)
!>          The upper or lower triangle of the hermitian band matrix A,
!>          stored in the first K+1 rows of AB.  The j-th column of A is
!>          stored in the j-th column of the array AB as follows:
!>          if UPLO = 'U', AB(k+1+i-j,j) = A(i,j) for max(1,j-k)<=i<=j;
!>          if UPLO = 'L', AB(1+i-j,j)   = A(i,j) for j<=i<=min(n,j+k).
!>          Note that the imaginary parts of the diagonal elements need
!>          not be set and are assumed to be zero.
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
!> \ingroup complexOTHERauxiliary
!
!  =====================================================================
      REAL FUNCTION CLANHB(Norm,Uplo,N,K,Ab,Ldab,Work)
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
      IMPLICIT NONE
!*--CLANHB141
!     .. Scalar Arguments ..
      CHARACTER Norm , Uplo
      INTEGER K , Ldab , N
!     ..
!     .. Array Arguments ..
      REAL Work(*)
      COMPLEX Ab(Ldab,*)
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
      EXTERNAL CLASSQ , SCOMBSSQ
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , MIN , REAL , SQRT
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
               DO i = MAX(K+2-j,1) , K
                  sum = ABS(Ab(i,j))
                  IF ( value<sum .OR. SISNAN(sum) ) value = sum
               ENDDO
               sum = ABS(REAL(Ab(K+1,j)))
               IF ( value<sum .OR. SISNAN(sum) ) value = sum
            ENDDO
         ELSE
            DO j = 1 , N
               sum = ABS(REAL(Ab(1,j)))
               IF ( value<sum .OR. SISNAN(sum) ) value = sum
               DO i = 2 , MIN(N+1-j,K+1)
                  sum = ABS(Ab(i,j))
                  IF ( value<sum .OR. SISNAN(sum) ) value = sum
               ENDDO
            ENDDO
         ENDIF
      ELSEIF ( (LSAME(Norm,'I')) .OR. (LSAME(Norm,'O')) .OR. (Norm=='1')&
     &         ) THEN
!
!        Find normI(A) ( = norm1(A), since A is hermitian).
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
               Work(j) = sum + ABS(REAL(Ab(K+1,j)))
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
               sum = Work(j) + ABS(REAL(Ab(1,j)))
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
                  CALL CLASSQ(MIN(j-1,K),Ab(MAX(K+2-j,1),j),1,colssq(1),&
     &                        colssq(2))
                  CALL SCOMBSSQ(ssq,colssq)
               ENDDO
               l = K + 1
            ELSE
               DO j = 1 , N - 1
                  colssq(1) = ZERO
                  colssq(2) = ONE
                  CALL CLASSQ(MIN(N-j,K),Ab(2,j),1,colssq(1),colssq(2))
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
         DO j = 1 , N
            IF ( REAL(Ab(l,j))/=ZERO ) THEN
               absa = ABS(REAL(Ab(l,j)))
               IF ( colssq(1)<absa ) THEN
                  colssq(2) = ONE + colssq(2)*(colssq(1)/absa)**2
                  colssq(1) = absa
               ELSE
                  colssq(2) = colssq(2) + (absa/colssq(1))**2
               ENDIF
            ENDIF
         ENDDO
         CALL SCOMBSSQ(ssq,colssq)
         value = ssq(1)*SQRT(ssq(2))
      ENDIF
!
      CLANHB = value
!
!     End of CLANHB
!
      END FUNCTION CLANHB
