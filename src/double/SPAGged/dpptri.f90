!*==dpptri.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DPPTRI
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DPPTRI + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dpptri.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dpptri.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dpptri.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DPPTRI( UPLO, N, AP, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   AP( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DPPTRI computes the inverse of a real symmetric positive definite
!> matrix A using the Cholesky factorization A = U**T*U or A = L*L**T
!> computed by DPPTRF.
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
!>          AP is DOUBLE PRECISION array, dimension (N*(N+1)/2)
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
!> \ingroup doubleOTHERcomputational
!
!  =====================================================================
      SUBROUTINE DPPTRI(Uplo,N,Ap,Info)
      USE F77KINDS                        
      USE S_DDOT
      USE S_DSCAL
      USE S_DSPR
      USE S_DTPMV
      USE S_DTPTRI
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--DPPTRI105
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Ap
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: ajj
      INTEGER :: j , jc , jj , jjn
      LOGICAL :: upper
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
         CALL XERBLA('DPPTRI',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
!     Invert the triangular Cholesky factor U or L.
!
      CALL DTPTRI(Uplo,'Non-unit',N,Ap,Info)
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
            IF ( j>1 ) CALL DSPR('Upper',j-1,ONE,Ap(jc),1,Ap)
            ajj = Ap(jj)
            CALL DSCAL(j,ajj,Ap(jc),1)
         ENDDO
!
      ELSE
!
!        Compute the product inv(L)**T * inv(L).
!
         jj = 1
         DO j = 1 , N
            jjn = jj + N - j + 1
            Ap(jj) = DDOT(N-j+1,Ap(jj),1,Ap(jj),1)
            IF ( j<N ) CALL DTPMV('Lower','Transpose','Non-unit',N-j,   &
     &                            Ap(jjn),Ap(jj+1),1)
            jj = jjn
         ENDDO
      ENDIF
!
!
!     End of DPPTRI
!
      END SUBROUTINE DPPTRI
