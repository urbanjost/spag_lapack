!*==cpptrf.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b CPPTRF
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CPPTRF + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpptrf.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpptrf.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpptrf.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CPPTRF( UPLO, N, AP, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, N
!       ..
!       .. Array Arguments ..
!       COMPLEX            AP( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CPPTRF computes the Cholesky factorization of a complex Hermitian
!> positive definite matrix A stored in packed format.
!>
!> The factorization has the form
!>    A = U**H * U,  if UPLO = 'U', or
!>    A = L  * L**H,  if UPLO = 'L',
!> where U is an upper triangular matrix and L is lower triangular.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  Upper triangle of A is stored;
!>          = 'L':  Lower triangle of A is stored.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] AP
!> \verbatim
!>          AP is COMPLEX array, dimension (N*(N+1)/2)
!>          On entry, the upper or lower triangle of the Hermitian matrix
!>          A, packed columnwise in a linear array.  The j-th column of A
!>          is stored in the array AP as follows:
!>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
!>          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
!>          See below for further details.
!>
!>          On exit, if INFO = 0, the triangular factor U or L from the
!>          Cholesky factorization A = U**H*U or A = L*L**H, in the same
!>          storage format as A.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  if INFO = i, the leading minor of order i is not
!>                positive definite, and the factorization could not be
!>                completed.
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
!> \ingroup complexOTHERcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The packed storage scheme is illustrated by the following example
!>  when N = 4, UPLO = 'U':
!>
!>  Two-dimensional storage of the Hermitian matrix A:
!>
!>     a11 a12 a13 a14
!>         a22 a23 a24
!>             a33 a34     (aij = conjg(aji))
!>                 a44
!>
!>  Packed storage of the upper triangle of A:
!>
!>  AP = [ a11, a12, a22, a13, a23, a33, a14, a24, a34, a44 ]
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE CPPTRF(Uplo,N,Ap,Info)
      IMPLICIT NONE
!*--CPPTRF123
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER Info , N
!     ..
!     .. Array Arguments ..
      COMPLEX Ap(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
!     ..
!     .. Local Scalars ..
      LOGICAL upper
      INTEGER j , jc , jj
      REAL ajj
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      COMPLEX CDOTC
      EXTERNAL LSAME , CDOTC
!     ..
!     .. External Subroutines ..
      EXTERNAL CHPR , CSSCAL , CTPSV , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC REAL , SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
      upper = LSAME(Uplo,'U')
      IF ( .NOT.upper .AND. .NOT.LSAME(Uplo,'L') ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('CPPTRF',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
      IF ( upper ) THEN
!
!        Compute the Cholesky factorization A = U**H * U.
!
         jj = 0
         DO j = 1 , N
            jc = jj + 1
            jj = jj + j
!
!           Compute elements 1:J-1 of column J.
!
            IF ( j>1 ) CALL CTPSV('Upper','Conjugate transpose',        &
     &                            'Non-unit',j-1,Ap,Ap(jc),1)
!
!           Compute U(J,J) and test for non-positive-definiteness.
!
            ajj = REAL(Ap(jj)) - CDOTC(j-1,Ap(jc),1,Ap(jc),1)
            IF ( ajj<=ZERO ) THEN
               Ap(jj) = ajj
               GOTO 100
            ENDIF
            Ap(jj) = SQRT(ajj)
         ENDDO
      ELSE
!
!        Compute the Cholesky factorization A = L * L**H.
!
         jj = 1
         DO j = 1 , N
!
!           Compute L(J,J) and test for non-positive-definiteness.
!
            ajj = REAL(Ap(jj))
            IF ( ajj<=ZERO ) THEN
               Ap(jj) = ajj
               GOTO 100
            ENDIF
            ajj = SQRT(ajj)
            Ap(jj) = ajj
!
!           Compute elements J+1:N of column J and update the trailing
!           submatrix.
!
            IF ( j<N ) THEN
               CALL CSSCAL(N-j,ONE/ajj,Ap(jj+1),1)
               CALL CHPR('Lower',N-j,-ONE,Ap(jj+1),1,Ap(jj+N-j+1))
               jj = jj + N - j + 1
            ENDIF
         ENDDO
      ENDIF
      GOTO 99999
!
 100  Info = j
!
!
!     End of CPPTRF
!
99999 END SUBROUTINE CPPTRF
