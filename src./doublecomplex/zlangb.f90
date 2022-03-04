!*==zlangb.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZLANGB returns the value of the 1-norm, Frobenius norm, infinity-norm, or the largest absolute value of any element of general band matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLANGB + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlangb.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlangb.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlangb.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION ZLANGB( NORM, N, KL, KU, AB, LDAB,
!                        WORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          NORM
!       INTEGER            KL, KU, LDAB, N
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
!> ZLANGB  returns the value of the one norm,  or the Frobenius norm, or
!> the  infinity norm,  or the element of  largest absolute value  of an
!> n by n band matrix  A,  with kl sub-diagonals and ku super-diagonals.
!> \endverbatim
!>
!> \return ZLANGB
!> \verbatim
!>
!>    ZLANGB = ( max(abs(A(i,j))), NORM = 'M' or 'm'
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
!>          Specifies the value to be returned in ZLANGB as described
!>          above.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.  When N = 0, ZLANGB is
!>          set to zero.
!> \endverbatim
!>
!> \param[in] KL
!> \verbatim
!>          KL is INTEGER
!>          The number of sub-diagonals of the matrix A.  KL >= 0.
!> \endverbatim
!>
!> \param[in] KU
!> \verbatim
!>          KU is INTEGER
!>          The number of super-diagonals of the matrix A.  KU >= 0.
!> \endverbatim
!>
!> \param[in] AB
!> \verbatim
!>          AB is COMPLEX*16 array, dimension (LDAB,N)
!>          The band matrix A, stored in rows 1 to KL+KU+1.  The j-th
!>          column of A is stored in the j-th column of the array AB as
!>          follows:
!>          AB(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(n,j+kl).
!> \endverbatim
!>
!> \param[in] LDAB
!> \verbatim
!>          LDAB is INTEGER
!>          The leading dimension of the array AB.  LDAB >= KL+KU+1.
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
!> \ingroup complex16GBauxiliary
!
!  =====================================================================
      FUNCTION ZLANGB(Norm,N,Kl,Ku,Ab,Ldab,Work)
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
      USE F77KINDS                        
      USE S_DCOMBSSQ
      USE S_DISNAN
      USE S_LSAME
      USE S_ZLASSQ
      IMPLICIT NONE
!*--ZLANGB139
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      REAL(R8KIND) :: ZLANGB
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Norm
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Kl
      INTEGER , INTENT(IN) :: Ku
      COMPLEX(CX16KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) , DIMENSION(2) :: colssq , ssq
      INTEGER :: i , j , k , l
      REAL(R8KIND) :: sum , temp , value
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
            DO i = MAX(Ku+2-j,1) , MIN(N+Ku+1-j,Kl+Ku+1)
               temp = ABS(Ab(i,j))
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
            DO i = MAX(Ku+2-j,1) , MIN(N+Ku+1-j,Kl+Ku+1)
               sum = sum + ABS(Ab(i,j))
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
            k = Ku + 1 - j
            DO i = MAX(1,j-Ku) , MIN(N,j+Kl)
               Work(i) = Work(i) + ABS(Ab(k+i,j))
            ENDDO
         ENDDO
         value = ZERO
         DO i = 1 , N
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
            l = MAX(1,j-Ku)
            k = Ku + 1 - j + l
            colssq(1) = ZERO
            colssq(2) = ONE
            CALL ZLASSQ(MIN(N,j+Kl)-l+1,Ab(k,j),1,colssq(1),colssq(2))
            CALL DCOMBSSQ(ssq,colssq)
         ENDDO
         value = ssq(1)*SQRT(ssq(2))
      ENDIF
!
      ZLANGB = value
!
!     End of ZLANGB
!
      END FUNCTION ZLANGB
