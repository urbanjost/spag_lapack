!*==zhetri.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZHETRI
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZHETRI + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhetri.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhetri.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhetri.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZHETRI( UPLO, N, A, LDA, IPIV, WORK, INFO )
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
!> ZHETRI computes the inverse of a complex Hermitian indefinite matrix
!> A using the factorization A = U*D*U**H or A = L*D*L**H computed by
!> ZHETRF.
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
!>          = 'U':  Upper triangular, form is A = U*D*U**H;
!>          = 'L':  Lower triangular, form is A = L*D*L**H.
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
!>          used to obtain the factor U or L as computed by ZHETRF.
!>
!>          On exit, if INFO = 0, the (Hermitian) inverse of the original
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
!>          as determined by ZHETRF.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (N)
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
!> \ingroup complex16HEcomputational
!
!  =====================================================================
      SUBROUTINE ZHETRI(Uplo,N,A,Lda,Ipiv,Work,Info)
      USE F77KINDS                        
      USE S_LSAME
      USE S_XERBLA
      USE S_ZCOPY
      USE S_ZDOTC
      USE S_ZHEMV
      USE S_ZSWAP
      IMPLICIT NONE
!*--ZHETRI125
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0) ,       &
     &                 ZERO = (0.0D+0,0.0D+0)
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: ak , akp1 , d , t
      COMPLEX(CX16KIND) :: akkp1 , temp
      INTEGER :: j , k , kp , kstep
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
!     .. Intrinsic Functions ..
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
         CALL XERBLA('ZHETRI',-Info)
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
!        Compute inv(A) from the factorization A = U*D*U**H.
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
               A(k,k) = ONE/DBLE(A(k,k))
!
!           Compute column K of the inverse.
!
               IF ( k>1 ) THEN
                  CALL ZCOPY(k-1,A(1,k),1,Work,1)
                  CALL ZHEMV(Uplo,k-1,-CONE,A,Lda,Work,1,ZERO,A(1,k),1)
                  A(k,k) = A(k,k) - DBLE(ZDOTC(k-1,Work,1,A(1,k),1))
               ENDIF
               kstep = 1
            ELSE
!
!           2 x 2 diagonal block
!
!           Invert the diagonal block.
!
               t = ABS(A(k,k+1))
               ak = DBLE(A(k,k))/t
               akp1 = DBLE(A(k+1,k+1))/t
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
                  CALL ZHEMV(Uplo,k-1,-CONE,A,Lda,Work,1,ZERO,A(1,k),1)
                  A(k,k) = A(k,k) - DBLE(ZDOTC(k-1,Work,1,A(1,k),1))
                  A(k,k+1) = A(k,k+1) - ZDOTC(k-1,A(1,k),1,A(1,k+1),1)
                  CALL ZCOPY(k-1,A(1,k+1),1,Work,1)
                  CALL ZHEMV(Uplo,k-1,-CONE,A,Lda,Work,1,ZERO,A(1,k+1), &
     &                       1)
                  A(k+1,k+1) = A(k+1,k+1)                               &
     &                         - DBLE(ZDOTC(k-1,Work,1,A(1,k+1),1))
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
               DO j = kp + 1 , k - 1
                  temp = DCONJG(A(j,k))
                  A(j,k) = DCONJG(A(kp,j))
                  A(kp,j) = temp
               ENDDO
               A(kp,k) = DCONJG(A(kp,k))
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
!        Compute inv(A) from the factorization A = L*D*L**H.
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
               A(k,k) = ONE/DBLE(A(k,k))
!
!           Compute column K of the inverse.
!
               IF ( k<N ) THEN
                  CALL ZCOPY(N-k,A(k+1,k),1,Work,1)
                  CALL ZHEMV(Uplo,N-k,-CONE,A(k+1,k+1),Lda,Work,1,ZERO, &
     &                       A(k+1,k),1)
                  A(k,k) = A(k,k) - DBLE(ZDOTC(N-k,Work,1,A(k+1,k),1))
               ENDIF
               kstep = 1
            ELSE
!
!           2 x 2 diagonal block
!
!           Invert the diagonal block.
!
               t = ABS(A(k,k-1))
               ak = DBLE(A(k-1,k-1))/t
               akp1 = DBLE(A(k,k))/t
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
                  CALL ZHEMV(Uplo,N-k,-CONE,A(k+1,k+1),Lda,Work,1,ZERO, &
     &                       A(k+1,k),1)
                  A(k,k) = A(k,k) - DBLE(ZDOTC(N-k,Work,1,A(k+1,k),1))
                  A(k,k-1) = A(k,k-1)                                   &
     &                       - ZDOTC(N-k,A(k+1,k),1,A(k+1,k-1),1)
                  CALL ZCOPY(N-k,A(k+1,k-1),1,Work,1)
                  CALL ZHEMV(Uplo,N-k,-CONE,A(k+1,k+1),Lda,Work,1,ZERO, &
     &                       A(k+1,k-1),1)
                  A(k-1,k-1) = A(k-1,k-1)                               &
     &                         - DBLE(ZDOTC(N-k,Work,1,A(k+1,k-1),1))
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
               DO j = k + 1 , kp - 1
                  temp = DCONJG(A(j,k))
                  A(j,k) = DCONJG(A(kp,j))
                  A(kp,j) = temp
               ENDDO
               A(kp,k) = DCONJG(A(kp,k))
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
!     End of ZHETRI
!
      END SUBROUTINE ZHETRI