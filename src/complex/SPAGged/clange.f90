!*==clange.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CLANGE returns the value of the 1-norm, Frobenius norm, infinity-norm, or the largest absolute value of any element of a general rectangular matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLANGE + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clange.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clange.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clange.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       REAL             FUNCTION CLANGE( NORM, M, N, A, LDA, WORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          NORM
!       INTEGER            LDA, M, N
!       ..
!       .. Array Arguments ..
!       REAL               WORK( * )
!       COMPLEX            A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLANGE  returns the value of the one norm,  or the Frobenius norm, or
!> the  infinity norm,  or the  element of  largest absolute value  of a
!> complex matrix A.
!> \endverbatim
!>
!> \return CLANGE
!> \verbatim
!>
!>    CLANGE = ( max(abs(A(i,j))), NORM = 'M' or 'm'
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
!>          Specifies the value to be returned in CLANGE as described
!>          above.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= 0.  When M = 0,
!>          CLANGE is set to zero.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  N >= 0.  When N = 0,
!>          CLANGE is set to zero.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
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
!>          WORK is REAL array, dimension (MAX(1,LWORK)),
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
!> \ingroup complexGEauxiliary
!
!  =====================================================================
      FUNCTION CLANGE(Norm,M,N,A,Lda,Work)
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
!*--CLANGE129
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      REAL :: CLANGE
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Norm
      INTEGER :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
!
! Local variable declarations rewritten by SPAG
!
      REAL , DIMENSION(2) :: colssq , ssq
      INTEGER :: i , j
      REAL :: sum , temp , value
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
               IF ( value<temp .OR. SISNAN(temp) ) value = temp
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
            IF ( value<sum .OR. SISNAN(sum) ) value = sum
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
            IF ( value<temp .OR. SISNAN(temp) ) value = temp
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
            CALL CLASSQ(M,A(1,j),1,colssq(1),colssq(2))
            CALL SCOMBSSQ(ssq,colssq)
         ENDDO
         value = ssq(1)*SQRT(ssq(2))
      ENDIF
!
      CLANGE = value
!
!     End of CLANGE
!
      END FUNCTION CLANGE
