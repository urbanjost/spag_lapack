!*==zlange.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b ZLANGE returns the value of the 1-norm, Frobenius norm, infinity-norm, or the largest absolute value of any element of a general rectangular matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLANGE + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlange.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlange.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlange.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION ZLANGE( NORM, M, N, A, LDA, WORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          NORM
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
!> ZLANGE  returns the value of the one norm,  or the Frobenius norm, or
!> the  infinity norm,  or the  element of  largest absolute value  of a
!> complex matrix A.
!> \endverbatim
!>
!> \return ZLANGE
!> \verbatim
!>
!>    ZLANGE = ( max(abs(A(i,j))), NORM = 'M' or 'm'
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
!>          Specifies the value to be returned in ZLANGE as described
!>          above.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= 0.  When M = 0,
!>          ZLANGE is set to zero.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  N >= 0.  When N = 0,
!>          ZLANGE is set to zero.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          The m by n matrix A.
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
!> \ingroup complex16GEauxiliary
!
!  =====================================================================
      DOUBLE PRECISION FUNCTION ZLANGE(Norm,M,N,A,Lda,Work)
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
      IMPLICIT NONE
!*--ZLANGE125
!     .. Scalar Arguments ..
      CHARACTER Norm
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
      INTEGER i , j
      DOUBLE PRECISION sum , value , temp
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
         value = ZERO
         DO j = 1 , N
            DO i = 1 , M
               temp = ABS(A(i,j))
               IF ( value<temp .OR. DISNAN(temp) ) value = temp
            ENDDO
         ENDDO
      ELSEIF ( (LSAME(Norm,'O')) .OR. (Norm=='1') ) THEN
!
!        Find norm1(A).
!
         value = ZERO
         DO j = 1 , N
            sum = ZERO
            DO i = 1 , M
               sum = sum + ABS(A(i,j))
            ENDDO
            IF ( value<sum .OR. DISNAN(sum) ) value = sum
         ENDDO
      ELSEIF ( LSAME(Norm,'I') ) THEN
!
!        Find normI(A).
!
         DO i = 1 , M
            Work(i) = ZERO
         ENDDO
         DO j = 1 , N
            DO i = 1 , M
               Work(i) = Work(i) + ABS(A(i,j))
            ENDDO
         ENDDO
         value = ZERO
         DO i = 1 , M
            temp = Work(i)
            IF ( value<temp .OR. DISNAN(temp) ) value = temp
         ENDDO
      ELSEIF ( (LSAME(Norm,'F')) .OR. (LSAME(Norm,'E')) ) THEN
!
!        Find normF(A).
!        SSQ(1) is scale
!        SSQ(2) is sum-of-squares
!        For better accuracy, sum each column separately.
!
         ssq(1) = ZERO
         ssq(2) = ONE
         DO j = 1 , N
            colssq(1) = ZERO
            colssq(2) = ONE
            CALL ZLASSQ(M,A(1,j),1,colssq(1),colssq(2))
            CALL DCOMBSSQ(ssq,colssq)
         ENDDO
         value = ssq(1)*SQRT(ssq(2))
      ENDIF
!
      ZLANGE = value
!
!     End of ZLANGE
!
      END FUNCTION ZLANGE
