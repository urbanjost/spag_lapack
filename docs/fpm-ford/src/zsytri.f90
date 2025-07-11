!*==zsytri.f90  processed by SPAG 7.51RB at 20:09 on  3 Mar 2022
 
!> \brief \b ZSYTRI
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZSYTRI + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsytri.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsytri.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsytri.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZSYTRI( UPLO, N, A, LDA, IPIV, WORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDA, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX*16         A( LDA, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZSYTRI computes the inverse of a complex symmetric indefinite matrix
!> A using the factorization A = U*D*U**T or A = L*D*L**T computed by
!> ZSYTRF.
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
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          On entry, the block diagonal matrix D and the multipliers
!>          used to obtain the factor U or L as computed by ZSYTRF.
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
!>          as determined by ZSYTRF.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (2*N)
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
!> \date December 2016
!
!> \ingroup complex16SYcomputational
!
!  =====================================================================
      SUBROUTINE ZSYTRI(Uplo,N,A,Lda,Ipiv,Work,Info)
      IMPLICIT NONE
!*--ZSYTRI119
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER Info , Lda , N
!     ..
!     .. Array Arguments ..
      INTEGER Ipiv(*)
      COMPLEX*16 A(Lda,*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16 ONE , ZERO
      PARAMETER (ONE=(1.0D+0,0.0D+0),ZERO=(0.0D+0,0.0D+0))
!     ..
!     .. Local Scalars ..
      LOGICAL upper
      INTEGER k , kp , kstep
      COMPLEX*16 ak , akkp1 , akp1 , d , t , temp
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      COMPLEX*16 ZDOTU
      EXTERNAL LSAME , ZDOTU
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA , ZCOPY , ZSWAP , ZSYMV
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
         CALL XERBLA('ZSYTRI',-Info)
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
                  CALL ZCOPY(k-1,A(1,k),1,Work,1)
                  CALL ZSYMV(Uplo,k-1,-ONE,A,Lda,Work,1,ZERO,A(1,k),1)
                  A(k,k) = A(k,k) - ZDOTU(k-1,Work,1,A(1,k),1)
               ENDIF
               kstep = 1
            ELSE
!
!           2 x 2 diagonal block
!
!           Invert the diagonal block.
!
               t = A(k,k+1)
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
                  CALL ZCOPY(k-1,A(1,k),1,Work,1)
                  CALL ZSYMV(Uplo,k-1,-ONE,A,Lda,Work,1,ZERO,A(1,k),1)
                  A(k,k) = A(k,k) - ZDOTU(k-1,Work,1,A(1,k),1)
                  A(k,k+1) = A(k,k+1) - ZDOTU(k-1,A(1,k),1,A(1,k+1),1)
                  CALL ZCOPY(k-1,A(1,k+1),1,Work,1)
                  CALL ZSYMV(Uplo,k-1,-ONE,A,Lda,Work,1,ZERO,A(1,k+1),1)
                  A(k+1,k+1) = A(k+1,k+1) - ZDOTU(k-1,Work,1,A(1,k+1),1)
               ENDIF
               kstep = 2
            ENDIF
!
            kp = ABS(Ipiv(k))
            IF ( kp/=k ) THEN
!
!           Interchange rows and columns K and KP in the leading
!           submatrix A(1:k+1,1:k+1)
!
               CALL ZSWAP(kp-1,A(1,k),1,A(1,kp),1)
               CALL ZSWAP(k-kp-1,A(kp+1,k),1,A(kp,kp+1),Lda)
               temp = A(k,k)
               A(k,k) = A(kp,kp)
               A(kp,kp) = temp
               IF ( kstep==2 ) THEN
                  temp = A(k,k+1)
                  A(k,k+1) = A(kp,k+1)
                  A(kp,k+1) = temp
               ENDIF
            ENDIF
!
            k = k + kstep
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
                  CALL ZCOPY(N-k,A(k+1,k),1,Work,1)
                  CALL ZSYMV(Uplo,N-k,-ONE,A(k+1,k+1),Lda,Work,1,ZERO,  &
     &                       A(k+1,k),1)
                  A(k,k) = A(k,k) - ZDOTU(N-k,Work,1,A(k+1,k),1)
               ENDIF
               kstep = 1
            ELSE
!
!           2 x 2 diagonal block
!
!           Invert the diagonal block.
!
               t = A(k,k-1)
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
                  CALL ZCOPY(N-k,A(k+1,k),1,Work,1)
                  CALL ZSYMV(Uplo,N-k,-ONE,A(k+1,k+1),Lda,Work,1,ZERO,  &
     &                       A(k+1,k),1)
                  A(k,k) = A(k,k) - ZDOTU(N-k,Work,1,A(k+1,k),1)
                  A(k,k-1) = A(k,k-1)                                   &
     &                       - ZDOTU(N-k,A(k+1,k),1,A(k+1,k-1),1)
                  CALL ZCOPY(N-k,A(k+1,k-1),1,Work,1)
                  CALL ZSYMV(Uplo,N-k,-ONE,A(k+1,k+1),Lda,Work,1,ZERO,  &
     &                       A(k+1,k-1),1)
                  A(k-1,k-1) = A(k-1,k-1)                               &
     &                         - ZDOTU(N-k,Work,1,A(k+1,k-1),1)
               ENDIF
               kstep = 2
            ENDIF
!
            kp = ABS(Ipiv(k))
            IF ( kp/=k ) THEN
!
!           Interchange rows and columns K and KP in the trailing
!           submatrix A(k-1:n,k-1:n)
!
               IF ( kp<N ) CALL ZSWAP(N-kp,A(kp+1,k),1,A(kp+1,kp),1)
               CALL ZSWAP(kp-k-1,A(k+1,k),1,A(kp,k+1),Lda)
               temp = A(k,k)
               A(k,k) = A(kp,kp)
               A(kp,kp) = temp
               IF ( kstep==2 ) THEN
                  temp = A(k,k-1)
                  A(k,k-1) = A(kp,k-1)
                  A(kp,k-1) = temp
               ENDIF
            ENDIF
!
            k = k - kstep
         ENDDO
      ENDIF
!
!
!     End of ZSYTRI
!
      END SUBROUTINE ZSYTRI
