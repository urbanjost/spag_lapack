!*==spptri.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SPPTRI
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SPPTRI + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/spptri.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/spptri.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/spptri.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SPPTRI( UPLO, N, AP, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, N
!       ..
!       .. Array Arguments ..
!       REAL               AP( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SPPTRI computes the inverse of a real symmetric positive definite
!> matrix A using the Cholesky factorization A = U**T*U or A = L*L**T
!> computed by SPPTRF.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  Upper triangular factor is stored in AP;
!>          = 'L':  Lower triangular factor is stored in AP.
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
!>          AP is REAL array, dimension (N*(N+1)/2)
!>          On entry, the triangular factor U or L from the Cholesky
!>          factorization A = U**T*U or A = L*L**T, packed columnwise as
!>          a linear array.  The j-th column of U or L is stored in the
!>          array AP as follows:
!>          if UPLO = 'U', AP(i + (j-1)*j/2) = U(i,j) for 1<=i<=j;
!>          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = L(i,j) for j<=i<=n.
!>
!>          On exit, the upper or lower triangle of the (symmetric)
!>          inverse of A, overwriting the input factor U or L.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  if INFO = i, the (i,i) element of the factor U or L is
!>                zero, and the inverse could not be computed.
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
!> \ingroup realOTHERcomputational
!
!  =====================================================================
      SUBROUTINE SPPTRI(Uplo,N,Ap,Info)
      IMPLICIT NONE
!*--SPPTRI97
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
      REAL Ap(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ONE
      PARAMETER (ONE=1.0E+0)
!     ..
!     .. Local Scalars ..
      LOGICAL upper
      INTEGER j , jc , jj , jjn
      REAL ajj
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      REAL SDOT
      EXTERNAL LSAME , SDOT
!     ..
!     .. External Subroutines ..
      EXTERNAL SSCAL , SSPR , STPMV , STPTRI , XERBLA
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
         CALL XERBLA('SPPTRI',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
!     Invert the triangular Cholesky factor U or L.
!
      CALL STPTRI(Uplo,'Non-unit',N,Ap,Info)
      IF ( Info>0 ) RETURN
!
      IF ( upper ) THEN
!
!        Compute the product inv(U) * inv(U)**T.
!
         jj = 0
         DO j = 1 , N
            jc = jj + 1
            jj = jj + j
            IF ( j>1 ) CALL SSPR('Upper',j-1,ONE,Ap(jc),1,Ap)
            ajj = Ap(jj)
            CALL SSCAL(j,ajj,Ap(jc),1)
         ENDDO
!
      ELSE
!
!        Compute the product inv(L)**T * inv(L).
!
         jj = 1
         DO j = 1 , N
            jjn = jj + N - j + 1
            Ap(jj) = SDOT(N-j+1,Ap(jj),1,Ap(jj),1)
            IF ( j<N ) CALL STPMV('Lower','Transpose','Non-unit',N-j,   &
     &                            Ap(jjn),Ap(jj+1),1)
            jj = jjn
         ENDDO
      ENDIF
!
!
!     End of SPPTRI
!
      END SUBROUTINE SPPTRI
