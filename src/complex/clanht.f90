!*==clanht.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b CLANHT returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the element of largest absolute value of a complex Hermitian tridiagonal matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLANHT + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clanht.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clanht.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clanht.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       REAL             FUNCTION CLANHT( NORM, N, D, E )
!
!       .. Scalar Arguments ..
!       CHARACTER          NORM
!       INTEGER            N
!       ..
!       .. Array Arguments ..
!       REAL               D( * )
!       COMPLEX            E( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLANHT  returns the value of the one norm,  or the Frobenius norm, or
!> the  infinity norm,  or the  element of  largest absolute value  of a
!> complex Hermitian tridiagonal matrix A.
!> \endverbatim
!>
!> \return CLANHT
!> \verbatim
!>
!>    CLANHT = ( max(abs(A(i,j))), NORM = 'M' or 'm'
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
!>          Specifies the value to be returned in CLANHT as described
!>          above.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.  When N = 0, CLANHT is
!>          set to zero.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is REAL array, dimension (N)
!>          The diagonal elements of A.
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is COMPLEX array, dimension (N-1)
!>          The (n-1) sub-diagonal or super-diagonal elements of A.
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
      REAL FUNCTION CLANHT(Norm,N,D,E)
      IMPLICIT NONE
!*--CLANHT105
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Norm
      INTEGER N
!     ..
!     .. Array Arguments ..
      REAL D(*)
      COMPLEX E(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ONE , ZERO
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
!     ..
!     .. Local Scalars ..
      INTEGER i
      REAL anorm , scale , sum
!     ..
!     .. External Functions ..
      LOGICAL LSAME , SISNAN
      EXTERNAL LSAME , SISNAN
!     ..
!     .. External Subroutines ..
      EXTERNAL CLASSQ , SLASSQ
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , SQRT
!     ..
!     .. Executable Statements ..
!
      IF ( N<=0 ) THEN
         anorm = ZERO
      ELSEIF ( LSAME(Norm,'M') ) THEN
!
!        Find max(abs(A(i,j))).
!
         anorm = ABS(D(N))
         DO i = 1 , N - 1
            sum = ABS(D(i))
            IF ( anorm<sum .OR. SISNAN(sum) ) anorm = sum
            sum = ABS(E(i))
            IF ( anorm<sum .OR. SISNAN(sum) ) anorm = sum
         ENDDO
      ELSEIF ( LSAME(Norm,'O') .OR. Norm=='1' .OR. LSAME(Norm,'I') )    &
     &         THEN
!
!        Find norm1(A).
!
         IF ( N==1 ) THEN
            anorm = ABS(D(1))
         ELSE
            anorm = ABS(D(1)) + ABS(E(1))
            sum = ABS(E(N-1)) + ABS(D(N))
            IF ( anorm<sum .OR. SISNAN(sum) ) anorm = sum
            DO i = 2 , N - 1
               sum = ABS(D(i)) + ABS(E(i)) + ABS(E(i-1))
               IF ( anorm<sum .OR. SISNAN(sum) ) anorm = sum
            ENDDO
         ENDIF
      ELSEIF ( (LSAME(Norm,'F')) .OR. (LSAME(Norm,'E')) ) THEN
!
!        Find normF(A).
!
         scale = ZERO
         sum = ONE
         IF ( N>1 ) THEN
            CALL CLASSQ(N-1,E,1,scale,sum)
            sum = 2*sum
         ENDIF
         CALL SLASSQ(N,D,1,scale,sum)
         anorm = scale*SQRT(sum)
      ENDIF
!
      CLANHT = anorm
!
!     End of CLANHT
!
      END FUNCTION CLANHT
