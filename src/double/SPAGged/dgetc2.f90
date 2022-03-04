!*==dgetc2.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DGETC2 computes the LU factorization with complete pivoting of the general n-by-n matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DGETC2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgetc2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgetc2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgetc2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGETC2( N, A, LDA, IPIV, JPIV, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * ), JPIV( * )
!       DOUBLE PRECISION   A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGETC2 computes an LU factorization with complete pivoting of the
!> n-by-n matrix A. The factorization has the form A = P * L * U * Q,
!> where P and Q are permutation matrices, L is lower triangular with
!> unit diagonal elements and U is upper triangular.
!>
!> This is the Level 2 BLAS algorithm.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A. N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA, N)
!>          On entry, the n-by-n matrix A to be factored.
!>          On exit, the factors L and U from the factorization
!>          A = P*L*U*Q; the unit diagonal elements of L are not stored.
!>          If U(k, k) appears to be less than SMIN, U(k, k) is given the
!>          value of SMIN, i.e., giving a nonsingular perturbed system.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension(N).
!>          The pivot indices; for 1 <= i <= N, row i of the
!>          matrix has been interchanged with row IPIV(i).
!> \endverbatim
!>
!> \param[out] JPIV
!> \verbatim
!>          JPIV is INTEGER array, dimension(N).
!>          The pivot indices; for 1 <= j <= N, column j of the
!>          matrix has been interchanged with column JPIV(j).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>           = 0: successful exit
!>           > 0: if INFO = k, U(k, k) is likely to produce overflow if
!>                we try to solve for x in Ax = b. So U is perturbed to
!>                avoid the overflow.
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
!> \ingroup doubleGEauxiliary
!
!> \par Contributors:
!  ==================
!>
!>     Bo Kagstrom and Peter Poromaa, Department of Computing Science,
!>     Umea University, S-901 87 Umea, Sweden.
!
!  =====================================================================
      SUBROUTINE DGETC2(N,A,Lda,Ipiv,Jpiv,Info)
      USE F77KINDS                        
      USE S_DGER
      USE S_DLABAD
      USE S_DLAMCH
      USE S_DSWAP
      IMPLICIT NONE
!*--DGETC2120
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Jpiv
      INTEGER , INTENT(OUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: bignum , eps , smin , smlnum , xmax
      INTEGER :: i , ip , ipv , j , jp , jpv
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
!     .. External Subroutines ..
!     ..
!     .. External Functions ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
      Info = 0
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
!     Set constants to control overflow
!
      eps = DLAMCH('P')
      smlnum = DLAMCH('S')/eps
      bignum = ONE/smlnum
      CALL DLABAD(smlnum,bignum)
!
!     Handle the case N=1 by itself
!
      IF ( N==1 ) THEN
         Ipiv(1) = 1
         Jpiv(1) = 1
         IF ( ABS(A(1,1))<smlnum ) THEN
            Info = 1
            A(1,1) = smlnum
         ENDIF
         RETURN
      ENDIF
!
!     Factorize A using complete pivoting.
!     Set pivots less than SMIN to SMIN.
!
      DO i = 1 , N - 1
!
!        Find max element in matrix A
!
         xmax = ZERO
         DO ip = i , N
            DO jp = i , N
               IF ( ABS(A(ip,jp))>=xmax ) THEN
                  xmax = ABS(A(ip,jp))
                  ipv = ip
                  jpv = jp
               ENDIF
            ENDDO
         ENDDO
         IF ( i==1 ) smin = MAX(eps*xmax,smlnum)
!
!        Swap rows
!
         IF ( ipv/=i ) CALL DSWAP(N,A(ipv,1),Lda,A(i,1),Lda)
         Ipiv(i) = ipv
!
!        Swap columns
!
         IF ( jpv/=i ) CALL DSWAP(N,A(1,jpv),1,A(1,i),1)
         Jpiv(i) = jpv
!
!        Check for singularity
!
         IF ( ABS(A(i,i))<smin ) THEN
            Info = i
            A(i,i) = smin
         ENDIF
         DO j = i + 1 , N
            A(j,i) = A(j,i)/A(i,i)
         ENDDO
         CALL DGER(N-i,N-i,-ONE,A(i+1,i),1,A(i,i+1),Lda,A(i+1,i+1),Lda)
      ENDDO
!
      IF ( ABS(A(N,N))<smin ) THEN
         Info = N
         A(N,N) = smin
      ENDIF
!
!     Set last pivots to N
!
      Ipiv(N) = N
      Jpiv(N) = N
!
!
!     End of DGETC2
!
      END SUBROUTINE DGETC2
