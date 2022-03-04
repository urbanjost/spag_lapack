!*==ssytri_rook.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SSYTRI_ROOK
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SSYTRI_ROOK + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssytri_rook.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssytri_rook.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssytri_rook.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SSYTRI_ROOK( UPLO, N, A, LDA, IPIV, WORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDA, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       REAL               A( LDA, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SSYTRI_ROOK computes the inverse of a real symmetric
!> matrix A using the factorization A = U*D*U**T or A = L*D*L**T
!> computed by SSYTRF_ROOK.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the details of the factorization are stored
!>          as an upper or lower triangular matrix.
!>          = 'U':  Upper triangular, form is A = U*D*U**T;
!>          = 'L':  Lower triangular, form is A = L*D*L**T.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>          On entry, the block diagonal matrix D and the multipliers
!>          used to obtain the factor U or L as computed by SSYTRF_ROOK.
!>
!>          On exit, if INFO = 0, the (symmetric) inverse of the original
!>          matrix.  If UPLO = 'U', the upper triangular part of the
!>          inverse is formed and the part of A below the diagonal is not
!>          referenced; if UPLO = 'L' the lower triangular part of the
!>          inverse is formed and the part of A above the diagonal is
!>          not referenced.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          Details of the interchanges and the block structure of D
!>          as determined by SSYTRF_ROOK.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
!>          > 0: if INFO = i, D(i,i) = 0; the matrix is singular and its
!>               inverse could not be computed.
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
!> \date April 2012
!
!> \ingroup realSYcomputational
!
!> \par Contributors:
!  ==================
!>
!> \verbatim
!>
!>   April 2012, Igor Kozachenko,
!>                  Computer Science Division,
!>                  University of California, Berkeley
!>
!>  September 2007, Sven Hammarling, Nicholas J. Higham, Craig Lucas,
!>                  School of Mathematics,
!>                  University of Manchester
!>
!> \endverbatim
!
!  =====================================================================
      SUBROUTINE SSYTRI_ROOK(Uplo,N,A,Lda,Ipiv,Work,Info)
      IMPLICIT NONE
!*--SSYTRI_ROOK133
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     April 2012
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER Info , Lda , N
!     ..
!     .. Array Arguments ..
      INTEGER Ipiv(*)
      REAL A(Lda,*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ONE , ZERO
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
!     ..
!     .. Local Scalars ..
      LOGICAL upper
      INTEGER k , kp , kstep
      REAL ak , akkp1 , akp1 , d , t , temp
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      REAL SDOT
      EXTERNAL LSAME , SDOT
!     ..
!     .. External Subroutines ..
      EXTERNAL SCOPY , SSWAP , SSYMV , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX
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
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -4
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('SSYTRI_ROOK',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
!     Check that the diagonal matrix D is nonsingular.
!
      IF ( upper ) THEN
!
!        Upper triangular storage: examine D from bottom to top
!
         DO Info = N , 1 , -1
            IF ( Ipiv(Info)>0 .AND. A(Info,Info)==ZERO ) RETURN
         ENDDO
      ELSE
!
!        Lower triangular storage: examine D from top to bottom.
!
         DO Info = 1 , N
            IF ( Ipiv(Info)>0 .AND. A(Info,Info)==ZERO ) RETURN
         ENDDO
      ENDIF
      Info = 0
!
      IF ( upper ) THEN
!
!        Compute inv(A) from the factorization A = U*D*U**T.
!
!        K is the main loop index, increasing from 1 to N in steps of
!        1 or 2, depending on the size of the diagonal blocks.
!
         k = 1
!
!        If K > N, exit from loop.
!
         DO WHILE ( k<=N )
!
            IF ( Ipiv(k)>0 ) THEN
!
!           1 x 1 diagonal block
!
!           Invert the diagonal block.
!
               A(k,k) = ONE/A(k,k)
!
!           Compute column K of the inverse.
!
               IF ( k>1 ) THEN
                  CALL SCOPY(k-1,A(1,k),1,Work,1)
                  CALL SSYMV(Uplo,k-1,-ONE,A,Lda,Work,1,ZERO,A(1,k),1)
                  A(k,k) = A(k,k) - SDOT(k-1,Work,1,A(1,k),1)
               ENDIF
               kstep = 1
            ELSE
!
!           2 x 2 diagonal block
!
!           Invert the diagonal block.
!
               t = ABS(A(k,k+1))
               ak = A(k,k)/t
               akp1 = A(k+1,k+1)/t
               akkp1 = A(k,k+1)/t
               d = t*(ak*akp1-ONE)
               A(k,k) = akp1/d
               A(k+1,k+1) = ak/d
               A(k,k+1) = -akkp1/d
!
!           Compute columns K and K+1 of the inverse.
!
               IF ( k>1 ) THEN
                  CALL SCOPY(k-1,A(1,k),1,Work,1)
                  CALL SSYMV(Uplo,k-1,-ONE,A,Lda,Work,1,ZERO,A(1,k),1)
                  A(k,k) = A(k,k) - SDOT(k-1,Work,1,A(1,k),1)
                  A(k,k+1) = A(k,k+1) - SDOT(k-1,A(1,k),1,A(1,k+1),1)
                  CALL SCOPY(k-1,A(1,k+1),1,Work,1)
                  CALL SSYMV(Uplo,k-1,-ONE,A,Lda,Work,1,ZERO,A(1,k+1),1)
                  A(k+1,k+1) = A(k+1,k+1) - SDOT(k-1,Work,1,A(1,k+1),1)
               ENDIF
               kstep = 2
            ENDIF
!
            IF ( kstep==1 ) THEN
!
!           Interchange rows and columns K and IPIV(K) in the leading
!           submatrix A(1:k+1,1:k+1)
!
               kp = Ipiv(k)
               IF ( kp/=k ) THEN
                  IF ( kp>1 ) CALL SSWAP(kp-1,A(1,k),1,A(1,kp),1)
                  CALL SSWAP(k-kp-1,A(kp+1,k),1,A(kp,kp+1),Lda)
                  temp = A(k,k)
                  A(k,k) = A(kp,kp)
                  A(kp,kp) = temp
               ENDIF
            ELSE
!
!           Interchange rows and columns K and K+1 with -IPIV(K) and
!           -IPIV(K+1)in the leading submatrix A(1:k+1,1:k+1)
!
               kp = -Ipiv(k)
               IF ( kp/=k ) THEN
                  IF ( kp>1 ) CALL SSWAP(kp-1,A(1,k),1,A(1,kp),1)
                  CALL SSWAP(k-kp-1,A(kp+1,k),1,A(kp,kp+1),Lda)
!
                  temp = A(k,k)
                  A(k,k) = A(kp,kp)
                  A(kp,kp) = temp
                  temp = A(k,k+1)
                  A(k,k+1) = A(kp,k+1)
                  A(kp,k+1) = temp
               ENDIF
!
               k = k + 1
               kp = -Ipiv(k)
               IF ( kp/=k ) THEN
                  IF ( kp>1 ) CALL SSWAP(kp-1,A(1,k),1,A(1,kp),1)
                  CALL SSWAP(k-kp-1,A(kp+1,k),1,A(kp,kp+1),Lda)
                  temp = A(k,k)
                  A(k,k) = A(kp,kp)
                  A(kp,kp) = temp
               ENDIF
            ENDIF
!
            k = k + 1
         ENDDO
!
      ELSE
!
!        Compute inv(A) from the factorization A = L*D*L**T.
!
!        K is the main loop index, increasing from 1 to N in steps of
!        1 or 2, depending on the size of the diagonal blocks.
!
         k = N
!
!        If K < 1, exit from loop.
!
         DO WHILE ( k>=1 )
!
            IF ( Ipiv(k)>0 ) THEN
!
!           1 x 1 diagonal block
!
!           Invert the diagonal block.
!
               A(k,k) = ONE/A(k,k)
!
!           Compute column K of the inverse.
!
               IF ( k<N ) THEN
                  CALL SCOPY(N-k,A(k+1,k),1,Work,1)
                  CALL SSYMV(Uplo,N-k,-ONE,A(k+1,k+1),Lda,Work,1,ZERO,  &
     &                       A(k+1,k),1)
                  A(k,k) = A(k,k) - SDOT(N-k,Work,1,A(k+1,k),1)
               ENDIF
               kstep = 1
            ELSE
!
!           2 x 2 diagonal block
!
!           Invert the diagonal block.
!
               t = ABS(A(k,k-1))
               ak = A(k-1,k-1)/t
               akp1 = A(k,k)/t
               akkp1 = A(k,k-1)/t
               d = t*(ak*akp1-ONE)
               A(k-1,k-1) = akp1/d
               A(k,k) = ak/d
               A(k,k-1) = -akkp1/d
!
!           Compute columns K-1 and K of the inverse.
!
               IF ( k<N ) THEN
                  CALL SCOPY(N-k,A(k+1,k),1,Work,1)
                  CALL SSYMV(Uplo,N-k,-ONE,A(k+1,k+1),Lda,Work,1,ZERO,  &
     &                       A(k+1,k),1)
                  A(k,k) = A(k,k) - SDOT(N-k,Work,1,A(k+1,k),1)
                  A(k,k-1) = A(k,k-1)                                   &
     &                       - SDOT(N-k,A(k+1,k),1,A(k+1,k-1),1)
                  CALL SCOPY(N-k,A(k+1,k-1),1,Work,1)
                  CALL SSYMV(Uplo,N-k,-ONE,A(k+1,k+1),Lda,Work,1,ZERO,  &
     &                       A(k+1,k-1),1)
                  A(k-1,k-1) = A(k-1,k-1)                               &
     &                         - SDOT(N-k,Work,1,A(k+1,k-1),1)
               ENDIF
               kstep = 2
            ENDIF
!
            IF ( kstep==1 ) THEN
!
!           Interchange rows and columns K and IPIV(K) in the trailing
!           submatrix A(k-1:n,k-1:n)
!
               kp = Ipiv(k)
               IF ( kp/=k ) THEN
                  IF ( kp<N ) CALL SSWAP(N-kp,A(kp+1,k),1,A(kp+1,kp),1)
                  CALL SSWAP(kp-k-1,A(k+1,k),1,A(kp,k+1),Lda)
                  temp = A(k,k)
                  A(k,k) = A(kp,kp)
                  A(kp,kp) = temp
               ENDIF
            ELSE
!
!           Interchange rows and columns K and K-1 with -IPIV(K) and
!           -IPIV(K-1) in the trailing submatrix A(k-1:n,k-1:n)
!
               kp = -Ipiv(k)
               IF ( kp/=k ) THEN
                  IF ( kp<N ) CALL SSWAP(N-kp,A(kp+1,k),1,A(kp+1,kp),1)
                  CALL SSWAP(kp-k-1,A(k+1,k),1,A(kp,k+1),Lda)
!
                  temp = A(k,k)
                  A(k,k) = A(kp,kp)
                  A(kp,kp) = temp
                  temp = A(k,k-1)
                  A(k,k-1) = A(kp,k-1)
                  A(kp,k-1) = temp
               ENDIF
!
               k = k - 1
               kp = -Ipiv(k)
               IF ( kp/=k ) THEN
                  IF ( kp<N ) CALL SSWAP(N-kp,A(kp+1,k),1,A(kp+1,kp),1)
                  CALL SSWAP(kp-k-1,A(k+1,k),1,A(kp,k+1),Lda)
                  temp = A(k,k)
                  A(k,k) = A(kp,kp)
                  A(kp,kp) = temp
               ENDIF
            ENDIF
!
            k = k - 1
         ENDDO
      ENDIF
!
!
!     End of SSYTRI_ROOK
!
      END SUBROUTINE SSYTRI_ROOK