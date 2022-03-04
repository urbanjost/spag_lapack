!*==dgeequ.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DGEEQU
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DGEEQU + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgeequ.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgeequ.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgeequ.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGEEQU( M, N, A, LDA, R, C, ROWCND, COLCND, AMAX,
!                          INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, M, N
!       DOUBLE PRECISION   AMAX, COLCND, ROWCND
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), C( * ), R( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGEEQU computes row and column scalings intended to equilibrate an
!> M-by-N matrix A and reduce its condition number.  R returns the row
!> scale factors and C the column scale factors, chosen to try to make
!> the largest element in each row and column of the matrix B with
!> elements B(i,j)=R(i)*A(i,j)*C(j) have absolute value 1.
!>
!> R(i) and C(j) are restricted to be between SMLNUM = smallest safe
!> number and BIGNUM = largest safe number.  Use of these scaling
!> factors is not guaranteed to reduce the condition number of A but
!> works well in practice.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          The M-by-N matrix whose equilibration factors are
!>          to be computed.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[out] R
!> \verbatim
!>          R is DOUBLE PRECISION array, dimension (M)
!>          If INFO = 0 or INFO > M, R contains the row scale factors
!>          for A.
!> \endverbatim
!>
!> \param[out] C
!> \verbatim
!>          C is DOUBLE PRECISION array, dimension (N)
!>          If INFO = 0,  C contains the column scale factors for A.
!> \endverbatim
!>
!> \param[out] ROWCND
!> \verbatim
!>          ROWCND is DOUBLE PRECISION
!>          If INFO = 0 or INFO > M, ROWCND contains the ratio of the
!>          smallest R(i) to the largest R(i).  If ROWCND >= 0.1 and
!>          AMAX is neither too large nor too small, it is not worth
!>          scaling by R.
!> \endverbatim
!>
!> \param[out] COLCND
!> \verbatim
!>          COLCND is DOUBLE PRECISION
!>          If INFO = 0, COLCND contains the ratio of the smallest
!>          C(i) to the largest C(i).  If COLCND >= 0.1, it is not
!>          worth scaling by C.
!> \endverbatim
!>
!> \param[out] AMAX
!> \verbatim
!>          AMAX is DOUBLE PRECISION
!>          Absolute value of largest matrix element.  If AMAX is very
!>          close to overflow or very close to underflow, the matrix
!>          should be scaled.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  if INFO = i,  and i is
!>                <= M:  the i-th row of A is exactly zero
!>                >  M:  the (i-M)-th column of A is exactly zero
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
!> \ingroup doubleGEcomputational
!
!  =====================================================================
      SUBROUTINE DGEEQU(M,N,A,Lda,R,C,Rowcnd,Colcnd,Amax,Info)
      USE F77KINDS                        
      USE S_DLAMCH
      USE S_XERBLA
      IMPLICIT NONE
!*--DGEEQU145
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: R
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: C
      REAL(R8KIND) , INTENT(OUT) :: Rowcnd
      REAL(R8KIND) , INTENT(OUT) :: Colcnd
      REAL(R8KIND) , INTENT(OUT) :: Amax
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: bignum , rcmax , rcmin , smlnum
      INTEGER :: i , j
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
!     Test the input parameters.
!
      Info = 0
      IF ( M<0 ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Lda<MAX(1,M) ) THEN
         Info = -4
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('DGEEQU',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( M==0 .OR. N==0 ) THEN
         Rowcnd = ONE
         Colcnd = ONE
         Amax = ZERO
         RETURN
      ENDIF
!
!     Get machine constants.
!
      smlnum = DLAMCH('S')
      bignum = ONE/smlnum
!
!     Compute row scale factors.
!
      DO i = 1 , M
         R(i) = ZERO
      ENDDO
!
!     Find the maximum element in each row.
!
      DO j = 1 , N
         DO i = 1 , M
            R(i) = MAX(R(i),ABS(A(i,j)))
         ENDDO
      ENDDO
!
!     Find the maximum and minimum scale factors.
!
      rcmin = bignum
      rcmax = ZERO
      DO i = 1 , M
         rcmax = MAX(rcmax,R(i))
         rcmin = MIN(rcmin,R(i))
      ENDDO
      Amax = rcmax
!
      IF ( rcmin==ZERO ) THEN
!
!        Find the first zero scale factor and return an error code.
!
         DO i = 1 , M
            IF ( R(i)==ZERO ) THEN
               Info = i
               RETURN
            ENDIF
         ENDDO
      ELSE
!
!        Invert the scale factors.
!
         DO i = 1 , M
            R(i) = ONE/MIN(MAX(R(i),smlnum),bignum)
         ENDDO
!
!        Compute ROWCND = min(R(I)) / max(R(I))
!
         Rowcnd = MAX(rcmin,smlnum)/MIN(rcmax,bignum)
      ENDIF
!
!     Compute column scale factors
!
      DO j = 1 , N
         C(j) = ZERO
      ENDDO
!
!     Find the maximum element in each column,
!     assuming the row scaling computed above.
!
      DO j = 1 , N
         DO i = 1 , M
            C(j) = MAX(C(j),ABS(A(i,j))*R(i))
         ENDDO
      ENDDO
!
!     Find the maximum and minimum scale factors.
!
      rcmin = bignum
      rcmax = ZERO
      DO j = 1 , N
         rcmin = MIN(rcmin,C(j))
         rcmax = MAX(rcmax,C(j))
      ENDDO
!
      IF ( rcmin==ZERO ) THEN
!
!        Find the first zero scale factor and return an error code.
!
         DO j = 1 , N
            IF ( C(j)==ZERO ) THEN
               Info = M + j
               RETURN
            ENDIF
         ENDDO
      ELSE
!
!        Invert the scale factors.
!
         DO j = 1 , N
            C(j) = ONE/MIN(MAX(C(j),smlnum),bignum)
         ENDDO
!
!        Compute COLCND = min(C(J)) / max(C(J))
!
         Colcnd = MAX(rcmin,smlnum)/MIN(rcmax,bignum)
      ENDIF
!
!
!     End of DGEEQU
!
      END SUBROUTINE DGEEQU
