!*==zgetc2.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b ZGETC2 computes the LU factorization with complete pivoting of the general n-by-n matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZGETC2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgetc2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgetc2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgetc2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZGETC2( N, A, LDA, IPIV, JPIV, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * ), JPIV( * )
!       COMPLEX*16         A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZGETC2 computes an LU factorization, using complete pivoting, of the
!> n-by-n matrix A. The factorization has the form A = P * L * U * Q,
!> where P and Q are permutation matrices, L is lower triangular with
!> unit diagonal elements and U is upper triangular.
!>
!> This is a level 1 BLAS version of the algorithm.
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
!>          A is COMPLEX*16 array, dimension (LDA, N)
!>          On entry, the n-by-n matrix to be factored.
!>          On exit, the factors L and U from the factorization
!>          A = P*L*U*Q; the unit diagonal elements of L are not stored.
!>          If U(k, k) appears to be less than SMIN, U(k, k) is given the
!>          value of SMIN, giving a nonsingular perturbed system.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1, N).
!> \endverbatim
!>
!> \param[out] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N).
!>          The pivot indices; for 1 <= i <= N, row i of the
!>          matrix has been interchanged with row IPIV(i).
!> \endverbatim
!>
!> \param[out] JPIV
!> \verbatim
!>          JPIV is INTEGER array, dimension (N).
!>          The pivot indices; for 1 <= j <= N, column j of the
!>          matrix has been interchanged with column JPIV(j).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>           = 0: successful exit
!>           > 0: if INFO = k, U(k, k) is likely to produce overflow if
!>                one tries to solve for x in Ax = b. So U is perturbed
!>                to avoid the overflow.
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
!> \ingroup complex16GEauxiliary
!
!> \par Contributors:
!  ==================
!>
!>     Bo Kagstrom and Peter Poromaa, Department of Computing Science,
!>     Umea University, S-901 87 Umea, Sweden.
!
!  =====================================================================
      SUBROUTINE ZGETC2(N,A,Lda,Ipiv,Jpiv,Info)
      IMPLICIT NONE
!*--ZGETC2115
!
!  -- LAPACK auxiliary routine (version 3.8.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!
!     .. Scalar Arguments ..
      INTEGER Info , Lda , N
!     ..
!     .. Array Arguments ..
      INTEGER Ipiv(*) , Jpiv(*)
      COMPLEX*16 A(Lda,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
!     ..
!     .. Local Scalars ..
      INTEGER i , ip , ipv , j , jp , jpv
      DOUBLE PRECISION bignum , eps , smin , smlnum , xmax
!     ..
!     .. External Subroutines ..
      EXTERNAL ZGERU , ZSWAP , DLABAD
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH
      EXTERNAL DLAMCH
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , DCMPLX , MAX
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
            A(1,1) = DCMPLX(smlnum,ZERO)
         ENDIF
         RETURN
      ENDIF
!
!     Factorize A using complete pivoting.
!     Set pivots less than SMIN to SMIN
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
         IF ( ipv/=i ) CALL ZSWAP(N,A(ipv,1),Lda,A(i,1),Lda)
         Ipiv(i) = ipv
!
!        Swap columns
!
         IF ( jpv/=i ) CALL ZSWAP(N,A(1,jpv),1,A(1,i),1)
         Jpiv(i) = jpv
!
!        Check for singularity
!
         IF ( ABS(A(i,i))<smin ) THEN
            Info = i
            A(i,i) = DCMPLX(smin,ZERO)
         ENDIF
         DO j = i + 1 , N
            A(j,i) = A(j,i)/A(i,i)
         ENDDO
         CALL ZGERU(N-i,N-i,-DCMPLX(ONE),A(i+1,i),1,A(i,i+1),Lda,       &
     &              A(i+1,i+1),Lda)
      ENDDO
!
      IF ( ABS(A(N,N))<smin ) THEN
         Info = N
         A(N,N) = DCMPLX(smin,ZERO)
      ENDIF
!
!     Set last pivots to N
!
      Ipiv(N) = N
      Jpiv(N) = N
!
!
!     End of ZGETC2
!
      END SUBROUTINE ZGETC2
