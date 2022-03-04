!*==dlanhs.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DLANHS returns the value of the 1-norm, Frobenius norm, infinity-norm, or the largest absolute value of any element of an upper Hessenberg matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLANHS + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlanhs.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlanhs.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlanhs.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION DLANHS( NORM, N, A, LDA, WORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          NORM
!       INTEGER            LDA, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLANHS  returns the value of the one norm,  or the Frobenius norm, or
!> the  infinity norm,  or the  element of  largest absolute value  of a
!> Hessenberg matrix A.
!> \endverbatim
!>
!> \return DLANHS
!> \verbatim
!>
!>    DLANHS = ( max(abs(A(i,j))), NORM = 'M' or 'm'
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
!>          Specifies the value to be returned in DLANHS as described
!>          above.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.  When N = 0, DLANHS is
!>          set to zero.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          The n by n upper Hessenberg matrix A; the part of A below the
!>          first sub-diagonal is not referenced.
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
      FUNCTION DLANHS(Norm,N,A,Lda,Work)
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
!*--DLANHS123
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      REAL(R8KIND) :: DLANHS
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Norm
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) , DIMENSION(2) :: colssq , ssq
      INTEGER :: i , j
      REAL(R8KIND) :: sum , value
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
         DO j = 1 , N
            DO i = 1 , MIN(N,j+1)
               sum = ABS(A(i,j))
               IF ( value<sum .OR. DISNAN(sum) ) value = sum
            ENDDO
         ENDDO
      ELSEIF ( (LSAME(Norm,'O')) .OR. (Norm=='1') ) THEN
!
!        Find norm1(A).
!
         value = ZERO
         DO j = 1 , N
            sum = ZERO
            DO i = 1 , MIN(N,j+1)
               sum = sum + ABS(A(i,j))
            ENDDO
            IF ( value<sum .OR. DISNAN(sum) ) value = sum
         ENDDO
      ELSEIF ( LSAME(Norm,'I') ) THEN
!
!        Find normI(A).
!
         DO i = 1 , N
            Work(i) = ZERO
         ENDDO
         DO j = 1 , N
            DO i = 1 , MIN(N,j+1)
               Work(i) = Work(i) + ABS(A(i,j))
            ENDDO
         ENDDO
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
         ssq(1) = ZERO
         ssq(2) = ONE
         DO j = 1 , N
            colssq(1) = ZERO
            colssq(2) = ONE
            CALL DLASSQ(MIN(N,j+1),A(1,j),1,colssq(1),colssq(2))
            CALL DCOMBSSQ(ssq,colssq)
         ENDDO
         value = ssq(1)*SQRT(ssq(2))
      ENDIF
!
      DLANHS = value
!
!     End of DLANHS
!
      END FUNCTION DLANHS
