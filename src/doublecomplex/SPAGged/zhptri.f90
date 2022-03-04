!*==zhptri.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZHPTRI
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZHPTRI + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhptri.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhptri.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhptri.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZHPTRI( UPLO, N, AP, IPIV, WORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX*16         AP( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZHPTRI computes the inverse of a complex Hermitian indefinite matrix
!> A in packed storage using the factorization A = U*D*U**H or
!> A = L*D*L**H computed by ZHPTRF.
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
!> \param[in,out] AP
!> \verbatim
!>          AP is COMPLEX*16 array, dimension (N*(N+1)/2)
!>          On entry, the block diagonal matrix D and the multipliers
!>          used to obtain the factor U or L as computed by ZHPTRF,
!>          stored as a packed triangular matrix.
!>
!>          On exit, if INFO = 0, the (Hermitian) inverse of the original
!>          matrix, stored as a packed triangular matrix. The j-th column
!>          of inv(A) is stored in the array AP as follows:
!>          if UPLO = 'U', AP(i + (j-1)*j/2) = inv(A)(i,j) for 1<=i<=j;
!>          if UPLO = 'L',
!>             AP(i + (j-1)*(2n-j)/2) = inv(A)(i,j) for j<=i<=n.
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          Details of the interchanges and the block structure of D
!>          as determined by ZHPTRF.
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
!> \ingroup complex16OTHERcomputational
!
!  =====================================================================
      SUBROUTINE ZHPTRI(Uplo,N,Ap,Ipiv,Work,Info)
      USE F77KINDS                        
      USE S_LSAME
      USE S_XERBLA
      USE S_ZCOPY
      USE S_ZDOTC
      USE S_ZHPMV
      USE S_ZSWAP
      IMPLICIT NONE
!*--ZHPTRI120
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
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Ap
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: ak , akp1 , d , t
      COMPLEX(CX16KIND) :: akkp1 , temp
      INTEGER :: j , k , kc , kcnext , kp , kpc , kstep , kx , npp
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
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZHPTRI',-Info)
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
         kp = N*(N+1)/2
         DO Info = N , 1 , -1
            IF ( Ipiv(Info)>0 .AND. Ap(kp)==ZERO ) RETURN
            kp = kp - Info
         ENDDO
      ELSE
!
!        Lower triangular storage: examine D from top to bottom.
!
         kp = 1
         DO Info = 1 , N
            IF ( Ipiv(Info)>0 .AND. Ap(kp)==ZERO ) RETURN
            kp = kp + N - Info + 1
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
         kc = 1
!
!        If K > N, exit from loop.
!
         DO WHILE ( k<=N )
!
            kcnext = kc + k
            IF ( Ipiv(k)>0 ) THEN
!
!           1 x 1 diagonal block
!
!           Invert the diagonal block.
!
               Ap(kc+k-1) = ONE/DBLE(Ap(kc+k-1))
!
!           Compute column K of the inverse.
!
               IF ( k>1 ) THEN
                  CALL ZCOPY(k-1,Ap(kc),1,Work,1)
                  CALL ZHPMV(Uplo,k-1,-CONE,Ap,Work,1,ZERO,Ap(kc),1)
                  Ap(kc+k-1) = Ap(kc+k-1)                               &
     &                         - DBLE(ZDOTC(k-1,Work,1,Ap(kc),1))
               ENDIF
               kstep = 1
            ELSE
!
!           2 x 2 diagonal block
!
!           Invert the diagonal block.
!
               t = ABS(Ap(kcnext+k-1))
               ak = DBLE(Ap(kc+k-1))/t
               akp1 = DBLE(Ap(kcnext+k))/t
               akkp1 = Ap(kcnext+k-1)/t
               d = t*(ak*akp1-ONE)
               Ap(kc+k-1) = akp1/d
               Ap(kcnext+k) = ak/d
               Ap(kcnext+k-1) = -akkp1/d
!
!           Compute columns K and K+1 of the inverse.
!
               IF ( k>1 ) THEN
                  CALL ZCOPY(k-1,Ap(kc),1,Work,1)
                  CALL ZHPMV(Uplo,k-1,-CONE,Ap,Work,1,ZERO,Ap(kc),1)
                  Ap(kc+k-1) = Ap(kc+k-1)                               &
     &                         - DBLE(ZDOTC(k-1,Work,1,Ap(kc),1))
                  Ap(kcnext+k-1) = Ap(kcnext+k-1)                       &
     &                             - ZDOTC(k-1,Ap(kc),1,Ap(kcnext),1)
                  CALL ZCOPY(k-1,Ap(kcnext),1,Work,1)
                  CALL ZHPMV(Uplo,k-1,-CONE,Ap,Work,1,ZERO,Ap(kcnext),1)
                  Ap(kcnext+k) = Ap(kcnext+k)                           &
     &                           - DBLE(ZDOTC(k-1,Work,1,Ap(kcnext),1))
               ENDIF
               kstep = 2
               kcnext = kcnext + k + 1
            ENDIF
!
            kp = ABS(Ipiv(k))
            IF ( kp/=k ) THEN
!
!           Interchange rows and columns K and KP in the leading
!           submatrix A(1:k+1,1:k+1)
!
               kpc = (kp-1)*kp/2 + 1
               CALL ZSWAP(kp-1,Ap(kc),1,Ap(kpc),1)
               kx = kpc + kp - 1
               DO j = kp + 1 , k - 1
                  kx = kx + j - 1
                  temp = DCONJG(Ap(kc+j-1))
                  Ap(kc+j-1) = DCONJG(Ap(kx))
                  Ap(kx) = temp
               ENDDO
               Ap(kc+kp-1) = DCONJG(Ap(kc+kp-1))
               temp = Ap(kc+k-1)
               Ap(kc+k-1) = Ap(kpc+kp-1)
               Ap(kpc+kp-1) = temp
               IF ( kstep==2 ) THEN
                  temp = Ap(kc+k+k-1)
                  Ap(kc+k+k-1) = Ap(kc+k+kp-1)
                  Ap(kc+k+kp-1) = temp
               ENDIF
            ENDIF
!
            k = k + kstep
            kc = kcnext
         ENDDO
!
      ELSE
!
!        Compute inv(A) from the factorization A = L*D*L**H.
!
!        K is the main loop index, increasing from 1 to N in steps of
!        1 or 2, depending on the size of the diagonal blocks.
!
         npp = N*(N+1)/2
         k = N
         kc = npp
!
!        If K < 1, exit from loop.
!
         DO WHILE ( k>=1 )
!
            kcnext = kc - (N-k+2)
            IF ( Ipiv(k)>0 ) THEN
!
!           1 x 1 diagonal block
!
!           Invert the diagonal block.
!
               Ap(kc) = ONE/DBLE(Ap(kc))
!
!           Compute column K of the inverse.
!
               IF ( k<N ) THEN
                  CALL ZCOPY(N-k,Ap(kc+1),1,Work,1)
                  CALL ZHPMV(Uplo,N-k,-CONE,Ap(kc+N-k+1),Work,1,ZERO,   &
     &                       Ap(kc+1),1)
                  Ap(kc) = Ap(kc) - DBLE(ZDOTC(N-k,Work,1,Ap(kc+1),1))
               ENDIF
               kstep = 1
            ELSE
!
!           2 x 2 diagonal block
!
!           Invert the diagonal block.
!
               t = ABS(Ap(kcnext+1))
               ak = DBLE(Ap(kcnext))/t
               akp1 = DBLE(Ap(kc))/t
               akkp1 = Ap(kcnext+1)/t
               d = t*(ak*akp1-ONE)
               Ap(kcnext) = akp1/d
               Ap(kc) = ak/d
               Ap(kcnext+1) = -akkp1/d
!
!           Compute columns K-1 and K of the inverse.
!
               IF ( k<N ) THEN
                  CALL ZCOPY(N-k,Ap(kc+1),1,Work,1)
                  CALL ZHPMV(Uplo,N-k,-CONE,Ap(kc+(N-k+1)),Work,1,ZERO, &
     &                       Ap(kc+1),1)
                  Ap(kc) = Ap(kc) - DBLE(ZDOTC(N-k,Work,1,Ap(kc+1),1))
                  Ap(kcnext+1) = Ap(kcnext+1)                           &
     &                           - ZDOTC(N-k,Ap(kc+1),1,Ap(kcnext+2),1)
                  CALL ZCOPY(N-k,Ap(kcnext+2),1,Work,1)
                  CALL ZHPMV(Uplo,N-k,-CONE,Ap(kc+(N-k+1)),Work,1,ZERO, &
     &                       Ap(kcnext+2),1)
                  Ap(kcnext) = Ap(kcnext)                               &
     &                         - DBLE(ZDOTC(N-k,Work,1,Ap(kcnext+2),1))
               ENDIF
               kstep = 2
               kcnext = kcnext - (N-k+3)
            ENDIF
!
            kp = ABS(Ipiv(k))
            IF ( kp/=k ) THEN
!
!           Interchange rows and columns K and KP in the trailing
!           submatrix A(k-1:n,k-1:n)
!
               kpc = npp - (N-kp+1)*(N-kp+2)/2 + 1
               IF ( kp<N ) CALL ZSWAP(N-kp,Ap(kc+kp-k+1),1,Ap(kpc+1),1)
               kx = kc + kp - k
               DO j = k + 1 , kp - 1
                  kx = kx + N - j + 1
                  temp = DCONJG(Ap(kc+j-k))
                  Ap(kc+j-k) = DCONJG(Ap(kx))
                  Ap(kx) = temp
               ENDDO
               Ap(kc+kp-k) = DCONJG(Ap(kc+kp-k))
               temp = Ap(kc)
               Ap(kc) = Ap(kpc)
               Ap(kpc) = temp
               IF ( kstep==2 ) THEN
                  temp = Ap(kc-N+k-1)
                  Ap(kc-N+k-1) = Ap(kc-N+kp-1)
                  Ap(kc-N+kp-1) = temp
               ENDIF
            ENDIF
!
            k = k - kstep
            kc = kcnext
         ENDDO
      ENDIF
!
!
!     End of ZHPTRI
!
      END SUBROUTINE ZHPTRI
