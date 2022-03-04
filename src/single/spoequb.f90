!*==spoequb.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SPOEQUB
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SPOEQUB + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/spoequb.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/spoequb.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/spoequb.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SPOEQUB( N, A, LDA, S, SCOND, AMAX, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, N
!       REAL               AMAX, SCOND
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), S( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SPOEQUB computes row and column scalings intended to equilibrate a
!> symmetric positive definite matrix A and reduce its condition number
!> (with respect to the two-norm).  S contains the scale factors,
!> S(i) = 1/sqrt(A(i,i)), chosen so that the scaled matrix B with
!> elements B(i,j) = S(i)*A(i,j)*S(j) has ones on the diagonal.  This
!> choice of S puts the condition number of B within a factor N of the
!> smallest possible condition number over all possible diagonal
!> scalings.
!>
!> This routine differs from SPOEQU by restricting the scaling factors
!> to a power of the radix.  Barring over- and underflow, scaling by
!> these factors introduces no additional rounding errors.  However, the
!> scaled diagonal entries are no longer approximately 1 but lie
!> between sqrt(radix) and 1/sqrt(radix).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>          The N-by-N symmetric positive definite matrix whose scaling
!>          factors are to be computed.  Only the diagonal elements of A
!>          are referenced.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is REAL array, dimension (N)
!>          If INFO = 0, S contains the scale factors for A.
!> \endverbatim
!>
!> \param[out] SCOND
!> \verbatim
!>          SCOND is REAL
!>          If INFO = 0, S contains the ratio of the smallest S(i) to
!>          the largest S(i).  If SCOND >= 0.1 and AMAX is neither too
!>          large nor too small, it is not worth scaling by S.
!> \endverbatim
!>
!> \param[out] AMAX
!> \verbatim
!>          AMAX is REAL
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
!>          > 0:  if INFO = i, the i-th diagonal element is nonpositive.
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
!> \ingroup realPOcomputational
!
!  =====================================================================
      SUBROUTINE SPOEQUB(N,A,Lda,S,Scond,Amax,Info)
      IMPLICIT NONE
!*--SPOEQUB122
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Info , Lda , N
      REAL Amax , Scond
!     ..
!     .. Array Arguments ..
      REAL A(Lda,*) , S(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
!     ..
!     .. Local Scalars ..
      INTEGER i
      REAL smin , base , tmp
!     ..
!     .. External Functions ..
      REAL SLAMCH
      EXTERNAL SLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN , SQRT , LOG , INT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
!     Positive definite only performs 1 pass of equilibration.
!
      Info = 0
      IF ( N<0 ) THEN
         Info = -1
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -3
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('SPOEQUB',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible.
!
      IF ( N==0 ) THEN
         Scond = ONE
         Amax = ZERO
         RETURN
      ENDIF
 
      base = SLAMCH('B')
      tmp = -0.5/LOG(base)
!
!     Find the minimum and maximum diagonal elements.
!
      S(1) = A(1,1)
      smin = S(1)
      Amax = S(1)
      DO i = 2 , N
         S(i) = A(i,i)
         smin = MIN(smin,S(i))
         Amax = MAX(Amax,S(i))
      ENDDO
!
      IF ( smin<=ZERO ) THEN
!
!        Find the first non-positive diagonal element and return.
!
         DO i = 1 , N
            IF ( S(i)<=ZERO ) THEN
               Info = i
               RETURN
            ENDIF
         ENDDO
      ELSE
!
!        Set the scale factors to the reciprocals
!        of the diagonal elements.
!
         DO i = 1 , N
            S(i) = base**INT(tmp*LOG(S(i)))
         ENDDO
!
!        Compute SCOND = min(S(I)) / max(S(I)).
!
         Scond = SQRT(smin)/SQRT(Amax)
      ENDIF
!
!
!     End of SPOEQUB
!
      END SUBROUTINE SPOEQUB
