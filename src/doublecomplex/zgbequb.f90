!*==zgbequb.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b ZGBEQUB
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZGBEQUB + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgbequb.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgbequb.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgbequb.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZGBEQUB( M, N, KL, KU, AB, LDAB, R, C, ROWCND, COLCND,
!                           AMAX, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, KL, KU, LDAB, M, N
!       DOUBLE PRECISION   AMAX, COLCND, ROWCND
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   C( * ), R( * )
!       COMPLEX*16         AB( LDAB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZGBEQUB computes row and column scalings intended to equilibrate an
!> M-by-N matrix A and reduce its condition number.  R returns the row
!> scale factors and C the column scale factors, chosen to try to make
!> the largest element in each row and column of the matrix B with
!> elements B(i,j)=R(i)*A(i,j)*C(j) have an absolute value of at most
!> the radix.
!>
!> R(i) and C(j) are restricted to be a power of the radix between
!> SMLNUM = smallest safe number and BIGNUM = largest safe number.  Use
!> of these scaling factors is not guaranteed to reduce the condition
!> number of A but works well in practice.
!>
!> This routine differs from ZGEEQU by restricting the scaling factors
!> to a power of the radix.  Barring over- and underflow, scaling by
!> these factors introduces no additional rounding errors.  However, the
!> scaled entries' magnitudes are no longer approximately 1 but lie
!> between sqrt(radix) and 1/sqrt(radix).
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
!> \param[in] KL
!> \verbatim
!>          KL is INTEGER
!>          The number of subdiagonals within the band of A.  KL >= 0.
!> \endverbatim
!>
!> \param[in] KU
!> \verbatim
!>          KU is INTEGER
!>          The number of superdiagonals within the band of A.  KU >= 0.
!> \endverbatim
!>
!> \param[in] AB
!> \verbatim
!>          AB is COMPLEX*16 array, dimension (LDAB,N)
!>          On entry, the matrix A in band storage, in rows 1 to KL+KU+1.
!>          The j-th column of A is stored in the j-th column of the
!>          array AB as follows:
!>          AB(KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+kl)
!> \endverbatim
!>
!> \param[in] LDAB
!> \verbatim
!>          LDAB is INTEGER
!>          The leading dimension of the array A.  LDAB >= max(1,M).
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
!> \date June 2016
!
!> \ingroup complex16GBcomputational
!
!  =====================================================================
      SUBROUTINE ZGBEQUB(M,N,Kl,Ku,Ab,Ldab,R,C,Rowcnd,Colcnd,Amax,Info)
      IMPLICIT NONE
!*--ZGBEQUB164
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!
!     .. Scalar Arguments ..
      INTEGER Info , Kl , Ku , Ldab , M , N
      DOUBLE PRECISION Amax , Colcnd , Rowcnd
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION C(*) , R(*)
      COMPLEX*16 Ab(Ldab,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ONE , ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
      INTEGER i , j , kd
      DOUBLE PRECISION bignum , rcmax , rcmin , smlnum , radix , logrdx
      COMPLEX*16 zdum
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH
      EXTERNAL DLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , MIN , LOG , REAL , DIMAG
!     ..
!     .. Statement Functions ..
      DOUBLE PRECISION CABS1
!     ..
!     .. Statement Function definitions ..
      CABS1(zdum) = ABS(DBLE(zdum)) + ABS(DIMAG(zdum))
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
      ELSEIF ( Kl<0 ) THEN
         Info = -3
      ELSEIF ( Ku<0 ) THEN
         Info = -4
      ELSEIF ( Ldab<Kl+Ku+1 ) THEN
         Info = -6
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZGBEQUB',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible.
!
      IF ( M==0 .OR. N==0 ) THEN
         Rowcnd = ONE
         Colcnd = ONE
         Amax = ZERO
         RETURN
      ENDIF
!
!     Get machine constants.  Assume SMLNUM is a power of the radix.
!
      smlnum = DLAMCH('S')
      bignum = ONE/smlnum
      radix = DLAMCH('B')
      logrdx = LOG(radix)
!
!     Compute row scale factors.
!
      DO i = 1 , M
         R(i) = ZERO
      ENDDO
!
!     Find the maximum element in each row.
!
      kd = Ku + 1
      DO j = 1 , N
         DO i = MAX(j-Ku,1) , MIN(j+Kl,M)
            R(i) = MAX(R(i),CABS1(Ab(kd+i-j,j)))
         ENDDO
      ENDDO
      DO i = 1 , M
         IF ( R(i)>ZERO ) R(i) = radix**INT(LOG(R(i))/logrdx)
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
!        Compute ROWCND = min(R(I)) / max(R(I)).
!
         Rowcnd = MAX(rcmin,smlnum)/MIN(rcmax,bignum)
      ENDIF
!
!     Compute column scale factors.
!
      DO j = 1 , N
         C(j) = ZERO
      ENDDO
!
!     Find the maximum element in each column,
!     assuming the row scaling computed above.
!
      DO j = 1 , N
         DO i = MAX(j-Ku,1) , MIN(j+Kl,M)
            C(j) = MAX(C(j),CABS1(Ab(kd+i-j,j))*R(i))
         ENDDO
         IF ( C(j)>ZERO ) C(j) = radix**INT(LOG(C(j))/logrdx)
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
!        Compute COLCND = min(C(J)) / max(C(J)).
!
         Colcnd = MAX(rcmin,smlnum)/MIN(rcmax,bignum)
      ENDIF
!
!
!     End of ZGBEQUB
!
      END SUBROUTINE ZGBEQUB
