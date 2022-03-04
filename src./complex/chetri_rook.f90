!*==chetri_rook.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CHETRI_ROOK computes the inverse of HE matrix using the factorization obtained with the bounded Bunch-Kaufman ("rook") diagonal pivoting method.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CHETRI_ROOK + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chetri_rook.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chetri_rook.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chetri_rook.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CHETRI_ROOK( UPLO, N, A, LDA, IPIV, WORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDA, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX            A( LDA, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CHETRI_ROOK computes the inverse of a complex Hermitian indefinite matrix
!> A using the factorization A = U*D*U**H or A = L*D*L**H computed by
!> CHETRF_ROOK.
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
!>          A is COMPLEX array, dimension (LDA,N)
!>          On entry, the block diagonal matrix D and the multipliers
!>          used to obtain the factor U or L as computed by CHETRF_ROOK.
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
!>          as determined by CHETRF_ROOK.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (N)
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
!> \date November 2013
!
!> \ingroup complexHEcomputational
!
!> \par Contributors:
!  ==================
!>
!> \verbatim
!>
!>  November 2013,  Igor Kozachenko,
!>                  Computer Science Division,
!>                  University of California, Berkeley
!>
!>  September 2007, Sven Hammarling, Nicholas J. Higham, Craig Lucas,
!>                  School of Mathematics,
!>                  University of Manchester
!> \endverbatim
!
!  =====================================================================
      SUBROUTINE CHETRI_ROOK(Uplo,N,A,Lda,Ipiv,Work,Info)
      USE S_CCOPY
      USE S_CDOTC
      USE S_CHEMV
      USE S_CSWAP
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--CHETRI_ROOK138
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0) ,                 &
     &                         CZERO = (0.0E+0,0.0E+0)
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL :: ak , akp1 , d , t
      COMPLEX :: akkp1 , temp
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
         CALL XERBLA('CHETRI_ROOK',-Info)
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
            IF ( Ipiv(Info)>0 .AND. A(Info,Info)==CZERO ) RETURN
         ENDDO
      ELSE
!
!        Lower triangular storage: examine D from top to bottom.
!
         DO Info = 1 , N
            IF ( Ipiv(Info)>0 .AND. A(Info,Info)==CZERO ) RETURN
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
               A(k,k) = ONE/REAL(A(k,k))
!
!           Compute column K of the inverse.
!
               IF ( k>1 ) THEN
                  CALL CCOPY(k-1,A(1,k),1,Work,1)
                  CALL CHEMV(Uplo,k-1,-CONE,A,Lda,Work,1,CZERO,A(1,k),1)
                  A(k,k) = A(k,k) - REAL(CDOTC(k-1,Work,1,A(1,k),1))
               ENDIF
               kstep = 1
            ELSE
!
!           2 x 2 diagonal block
!
!           Invert the diagonal block.
!
               t = ABS(A(k,k+1))
               ak = REAL(A(k,k))/t
               akp1 = REAL(A(k+1,k+1))/t
               akkp1 = A(k,k+1)/t
               d = t*(ak*akp1-ONE)
               A(k,k) = akp1/d
               A(k+1,k+1) = ak/d
               A(k,k+1) = -akkp1/d
!
!           Compute columns K and K+1 of the inverse.
!
               IF ( k>1 ) THEN
                  CALL CCOPY(k-1,A(1,k),1,Work,1)
                  CALL CHEMV(Uplo,k-1,-CONE,A,Lda,Work,1,CZERO,A(1,k),1)
                  A(k,k) = A(k,k) - REAL(CDOTC(k-1,Work,1,A(1,k),1))
                  A(k,k+1) = A(k,k+1) - CDOTC(k-1,A(1,k),1,A(1,k+1),1)
                  CALL CCOPY(k-1,A(1,k+1),1,Work,1)
                  CALL CHEMV(Uplo,k-1,-CONE,A,Lda,Work,1,CZERO,A(1,k+1),&
     &                       1)
                  A(k+1,k+1) = A(k+1,k+1)                               &
     &                         - REAL(CDOTC(k-1,Work,1,A(1,k+1),1))
               ENDIF
               kstep = 2
            ENDIF
!
            IF ( kstep==1 ) THEN
!
!           Interchange rows and columns K and IPIV(K) in the leading
!           submatrix A(1:k,1:k)
!
               kp = Ipiv(k)
               IF ( kp/=k ) THEN
!
                  IF ( kp>1 ) CALL CSWAP(kp-1,A(1,k),1,A(1,kp),1)
!
                  DO j = kp + 1 , k - 1
                     temp = CONJG(A(j,k))
                     A(j,k) = CONJG(A(kp,j))
                     A(kp,j) = temp
                  ENDDO
!
                  A(kp,k) = CONJG(A(kp,k))
!
                  temp = A(k,k)
                  A(k,k) = A(kp,kp)
                  A(kp,kp) = temp
               ENDIF
            ELSE
!
!           Interchange rows and columns K and K+1 with -IPIV(K) and
!           -IPIV(K+1) in the leading submatrix A(k+1:n,k+1:n)
!
!           (1) Interchange rows and columns K and -IPIV(K)
!
               kp = -Ipiv(k)
               IF ( kp/=k ) THEN
!
                  IF ( kp>1 ) CALL CSWAP(kp-1,A(1,k),1,A(1,kp),1)
!
                  DO j = kp + 1 , k - 1
                     temp = CONJG(A(j,k))
                     A(j,k) = CONJG(A(kp,j))
                     A(kp,j) = temp
                  ENDDO
!
                  A(kp,k) = CONJG(A(kp,k))
!
                  temp = A(k,k)
                  A(k,k) = A(kp,kp)
                  A(kp,kp) = temp
!
                  temp = A(k,k+1)
                  A(k,k+1) = A(kp,k+1)
                  A(kp,k+1) = temp
               ENDIF
!
!           (2) Interchange rows and columns K+1 and -IPIV(K+1)
!
               k = k + 1
               kp = -Ipiv(k)
               IF ( kp/=k ) THEN
!
                  IF ( kp>1 ) CALL CSWAP(kp-1,A(1,k),1,A(1,kp),1)
!
                  DO j = kp + 1 , k - 1
                     temp = CONJG(A(j,k))
                     A(j,k) = CONJG(A(kp,j))
                     A(kp,j) = temp
                  ENDDO
!
                  A(kp,k) = CONJG(A(kp,k))
!
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
!        Compute inv(A) from the factorization A = L*D*L**H.
!
!        K is the main loop index, decreasing from N to 1 in steps of
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
               A(k,k) = ONE/REAL(A(k,k))
!
!           Compute column K of the inverse.
!
               IF ( k<N ) THEN
                  CALL CCOPY(N-k,A(k+1,k),1,Work,1)
                  CALL CHEMV(Uplo,N-k,-CONE,A(k+1,k+1),Lda,Work,1,CZERO,&
     &                       A(k+1,k),1)
                  A(k,k) = A(k,k) - REAL(CDOTC(N-k,Work,1,A(k+1,k),1))
               ENDIF
               kstep = 1
            ELSE
!
!           2 x 2 diagonal block
!
!           Invert the diagonal block.
!
               t = ABS(A(k,k-1))
               ak = REAL(A(k-1,k-1))/t
               akp1 = REAL(A(k,k))/t
               akkp1 = A(k,k-1)/t
               d = t*(ak*akp1-ONE)
               A(k-1,k-1) = akp1/d
               A(k,k) = ak/d
               A(k,k-1) = -akkp1/d
!
!           Compute columns K-1 and K of the inverse.
!
               IF ( k<N ) THEN
                  CALL CCOPY(N-k,A(k+1,k),1,Work,1)
                  CALL CHEMV(Uplo,N-k,-CONE,A(k+1,k+1),Lda,Work,1,CZERO,&
     &                       A(k+1,k),1)
                  A(k,k) = A(k,k) - REAL(CDOTC(N-k,Work,1,A(k+1,k),1))
                  A(k,k-1) = A(k,k-1)                                   &
     &                       - CDOTC(N-k,A(k+1,k),1,A(k+1,k-1),1)
                  CALL CCOPY(N-k,A(k+1,k-1),1,Work,1)
                  CALL CHEMV(Uplo,N-k,-CONE,A(k+1,k+1),Lda,Work,1,CZERO,&
     &                       A(k+1,k-1),1)
                  A(k-1,k-1) = A(k-1,k-1)                               &
     &                         - REAL(CDOTC(N-k,Work,1,A(k+1,k-1),1))
               ENDIF
               kstep = 2
            ENDIF
!
            IF ( kstep==1 ) THEN
!
!           Interchange rows and columns K and IPIV(K) in the trailing
!           submatrix A(k:n,k:n)
!
               kp = Ipiv(k)
               IF ( kp/=k ) THEN
!
                  IF ( kp<N ) CALL CSWAP(N-kp,A(kp+1,k),1,A(kp+1,kp),1)
!
                  DO j = k + 1 , kp - 1
                     temp = CONJG(A(j,k))
                     A(j,k) = CONJG(A(kp,j))
                     A(kp,j) = temp
                  ENDDO
!
                  A(kp,k) = CONJG(A(kp,k))
!
                  temp = A(k,k)
                  A(k,k) = A(kp,kp)
                  A(kp,kp) = temp
               ENDIF
            ELSE
!
!           Interchange rows and columns K and K-1 with -IPIV(K) and
!           -IPIV(K-1) in the trailing submatrix A(k-1:n,k-1:n)
!
!           (1) Interchange rows and columns K and -IPIV(K)
!
               kp = -Ipiv(k)
               IF ( kp/=k ) THEN
!
                  IF ( kp<N ) CALL CSWAP(N-kp,A(kp+1,k),1,A(kp+1,kp),1)
!
                  DO j = k + 1 , kp - 1
                     temp = CONJG(A(j,k))
                     A(j,k) = CONJG(A(kp,j))
                     A(kp,j) = temp
                  ENDDO
!
                  A(kp,k) = CONJG(A(kp,k))
!
                  temp = A(k,k)
                  A(k,k) = A(kp,kp)
                  A(kp,kp) = temp
!
                  temp = A(k,k-1)
                  A(k,k-1) = A(kp,k-1)
                  A(kp,k-1) = temp
               ENDIF
!
!           (2) Interchange rows and columns K-1 and -IPIV(K-1)
!
               k = k - 1
               kp = -Ipiv(k)
               IF ( kp/=k ) THEN
!
                  IF ( kp<N ) CALL CSWAP(N-kp,A(kp+1,k),1,A(kp+1,kp),1)
!
                  DO j = k + 1 , kp - 1
                     temp = CONJG(A(j,k))
                     A(j,k) = CONJG(A(kp,j))
                     A(kp,j) = temp
                  ENDDO
!
                  A(kp,k) = CONJG(A(kp,k))
!
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
!     End of CHETRI_ROOK
!
      END SUBROUTINE CHETRI_ROOK
