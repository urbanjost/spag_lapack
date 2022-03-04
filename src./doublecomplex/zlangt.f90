!*==zlangt.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZLANGT returns the value of the 1-norm, Frobenius norm, infinity-norm, or the largest absolute value of any element of a general tridiagonal matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLANGT + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlangt.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlangt.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlangt.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION ZLANGT( NORM, N, DL, D, DU )
!
!       .. Scalar Arguments ..
!       CHARACTER          NORM
!       INTEGER            N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         D( * ), DL( * ), DU( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLANGT  returns the value of the one norm,  or the Frobenius norm, or
!> the  infinity norm,  or the  element of  largest absolute value  of a
!> complex tridiagonal matrix A.
!> \endverbatim
!>
!> \return ZLANGT
!> \verbatim
!>
!>    ZLANGT = ( max(abs(A(i,j))), NORM = 'M' or 'm'
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
!>          Specifies the value to be returned in ZLANGT as described
!>          above.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.  When N = 0, ZLANGT is
!>          set to zero.
!> \endverbatim
!>
!> \param[in] DL
!> \verbatim
!>          DL is COMPLEX*16 array, dimension (N-1)
!>          The (n-1) sub-diagonal elements of A.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is COMPLEX*16 array, dimension (N)
!>          The diagonal elements of A.
!> \endverbatim
!>
!> \param[in] DU
!> \verbatim
!>          DU is COMPLEX*16 array, dimension (N-1)
!>          The (n-1) super-diagonal elements of A.
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
      FUNCTION ZLANGT(Norm,N,Dl,D,Du)
      USE F77KINDS                        
      USE S_DISNAN
      USE S_LSAME
      USE S_ZLASSQ
      IMPLICIT NONE
!*--ZLANGT114
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      REAL(R8KIND) :: ZLANGT
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Norm
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(*) :: Dl
      COMPLEX(CX16KIND) , DIMENSION(*) :: D
      COMPLEX(CX16KIND) , DIMENSION(*) :: Du
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: anorm , scale , sum , temp
      INTEGER :: i
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
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
            IF ( anorm<ABS(Dl(i)) .OR. DISNAN(ABS(Dl(i))) )             &
     &           anorm = ABS(Dl(i))
            IF ( anorm<ABS(D(i)) .OR. DISNAN(ABS(D(i))) )               &
     &           anorm = ABS(D(i))
            IF ( anorm<ABS(Du(i)) .OR. DISNAN(ABS(Du(i))) )             &
     &           anorm = ABS(Du(i))
         ENDDO
      ELSEIF ( LSAME(Norm,'O') .OR. Norm=='1' ) THEN
!
!        Find norm1(A).
!
         IF ( N==1 ) THEN
            anorm = ABS(D(1))
         ELSE
            anorm = ABS(D(1)) + ABS(Dl(1))
            temp = ABS(D(N)) + ABS(Du(N-1))
            IF ( anorm<temp .OR. DISNAN(temp) ) anorm = temp
            DO i = 2 , N - 1
               temp = ABS(D(i)) + ABS(Dl(i)) + ABS(Du(i-1))
               IF ( anorm<temp .OR. DISNAN(temp) ) anorm = temp
            ENDDO
         ENDIF
      ELSEIF ( LSAME(Norm,'I') ) THEN
!
!        Find normI(A).
!
         IF ( N==1 ) THEN
            anorm = ABS(D(1))
         ELSE
            anorm = ABS(D(1)) + ABS(Du(1))
            temp = ABS(D(N)) + ABS(Dl(N-1))
            IF ( anorm<temp .OR. DISNAN(temp) ) anorm = temp
            DO i = 2 , N - 1
               temp = ABS(D(i)) + ABS(Du(i)) + ABS(Dl(i-1))
               IF ( anorm<temp .OR. DISNAN(temp) ) anorm = temp
            ENDDO
         ENDIF
      ELSEIF ( (LSAME(Norm,'F')) .OR. (LSAME(Norm,'E')) ) THEN
!
!        Find normF(A).
!
         scale = ZERO
         sum = ONE
         CALL ZLASSQ(N,D,1,scale,sum)
         IF ( N>1 ) THEN
            CALL ZLASSQ(N-1,Dl,1,scale,sum)
            CALL ZLASSQ(N-1,Du,1,scale,sum)
         ENDIF
         anorm = scale*SQRT(sum)
      ENDIF
!
      ZLANGT = anorm
!
!     End of ZLANGT
!
      END FUNCTION ZLANGT
