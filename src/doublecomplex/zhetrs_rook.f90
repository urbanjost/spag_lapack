!*==zhetrs_rook.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b ZHETRS_ROOK computes the solution to a system of linear equations A * X = B for HE matrices using factorization obtained with one of the bounded diagonal pivoting methods (max 2 interchanges)
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZHETRS_ROOK + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhetrs_rook.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhetrs_rook.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhetrs_rook.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZHETRS_ROOK( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDA, LDB, N, NRHS
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX            A( LDA, * ), B( LDB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZHETRS_ROOK solves a system of linear equations A*X = B with a complex
!> Hermitian matrix A using the factorization A = U*D*U**H or
!> A = L*D*L**H computed by ZHETRF_ROOK.
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
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of right hand sides, i.e., the number of columns
!>          of the matrix B.  NRHS >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          The block diagonal matrix D and the multipliers used to
!>          obtain the factor U or L as computed by ZHETRF_ROOK.
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
!>          as determined by ZHETRF_ROOK.
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is COMPLEX*16 array, dimension (LDB,NRHS)
!>          On entry, the right hand side matrix B.
!>          On exit, the solution matrix X.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,N).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
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
!> \ingroup complex16HEcomputational
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
!>
!> \endverbatim
!
!  =====================================================================
      SUBROUTINE ZHETRS_ROOK(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Info)
      IMPLICIT NONE
!*--ZHETRS_ROOK139
!
!  -- LAPACK computational routine (version 3.5.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2013
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER Info , Lda , Ldb , N , Nrhs
!     ..
!     .. Array Arguments ..
      INTEGER Ipiv(*)
      COMPLEX*16 A(Lda,*) , B(Ldb,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16 ONE
      PARAMETER (ONE=(1.0D+0,0.0D+0))
!     ..
!     .. Local Scalars ..
      LOGICAL upper
      INTEGER j , k , kp
      DOUBLE PRECISION s
      COMPLEX*16 ak , akm1 , akm1k , bk , bkm1 , denom
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL ZGEMV , ZGERU , ZLACGV , ZDSCAL , ZSWAP , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DCONJG , MAX , DBLE
!     ..
!     .. Executable Statements ..
!
      Info = 0
      upper = LSAME(Uplo,'U')
      IF ( .NOT.upper .AND. .NOT.LSAME(Uplo,'L') ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Nrhs<0 ) THEN
         Info = -3
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -5
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -8
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZHETRS_ROOK',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 .OR. Nrhs==0 ) RETURN
!
      IF ( upper ) THEN
!
!        Solve A*X = B, where A = U*D*U**H.
!
!        First solve U*D*X = B, overwriting B with X.
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
!           Interchange rows K and IPIV(K).
!
               kp = Ipiv(k)
               IF ( kp/=k ) CALL ZSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
!
!           Multiply by inv(U(K)), where U(K) is the transformation
!           stored in column K of A.
!
               CALL ZGERU(k-1,Nrhs,-ONE,A(1,k),1,B(k,1),Ldb,B(1,1),Ldb)
!
!           Multiply by the inverse of the diagonal block.
!
               s = DBLE(ONE)/DBLE(A(k,k))
               CALL ZDSCAL(Nrhs,s,B(k,1),Ldb)
               k = k - 1
            ELSE
!
!           2 x 2 diagonal block
!
!           Interchange rows K and -IPIV(K), then K-1 and -IPIV(K-1)
!
               kp = -Ipiv(k)
               IF ( kp/=k ) CALL ZSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
!
               kp = -Ipiv(k-1)
               IF ( kp/=k-1 ) CALL ZSWAP(Nrhs,B(k-1,1),Ldb,B(kp,1),Ldb)
!
!           Multiply by inv(U(K)), where U(K) is the transformation
!           stored in columns K-1 and K of A.
!
               CALL ZGERU(k-2,Nrhs,-ONE,A(1,k),1,B(k,1),Ldb,B(1,1),Ldb)
               CALL ZGERU(k-2,Nrhs,-ONE,A(1,k-1),1,B(k-1,1),Ldb,B(1,1), &
     &                    Ldb)
!
!           Multiply by the inverse of the diagonal block.
!
               akm1k = A(k-1,k)
               akm1 = A(k-1,k-1)/akm1k
               ak = A(k,k)/DCONJG(akm1k)
               denom = akm1*ak - ONE
               DO j = 1 , Nrhs
                  bkm1 = B(k-1,j)/akm1k
                  bk = B(k,j)/DCONJG(akm1k)
                  B(k-1,j) = (ak*bkm1-bk)/denom
                  B(k,j) = (akm1*bk-bkm1)/denom
               ENDDO
               k = k - 2
!
            ENDIF
         ENDDO
!
!        Next solve U**H *X = B, overwriting B with X.
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
!           Multiply by inv(U**H(K)), where U(K) is the transformation
!           stored in column K of A.
!
               IF ( k>1 ) THEN
                  CALL ZLACGV(Nrhs,B(k,1),Ldb)
                  CALL ZGEMV('Conjugate transpose',k-1,Nrhs,-ONE,B,Ldb, &
     &                       A(1,k),1,ONE,B(k,1),Ldb)
                  CALL ZLACGV(Nrhs,B(k,1),Ldb)
               ENDIF
!
!           Interchange rows K and IPIV(K).
!
               kp = Ipiv(k)
               IF ( kp/=k ) CALL ZSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
               k = k + 1
            ELSE
!
!           2 x 2 diagonal block
!
!           Multiply by inv(U**H(K+1)), where U(K+1) is the transformation
!           stored in columns K and K+1 of A.
!
               IF ( k>1 ) THEN
                  CALL ZLACGV(Nrhs,B(k,1),Ldb)
                  CALL ZGEMV('Conjugate transpose',k-1,Nrhs,-ONE,B,Ldb, &
     &                       A(1,k),1,ONE,B(k,1),Ldb)
                  CALL ZLACGV(Nrhs,B(k,1),Ldb)
!
                  CALL ZLACGV(Nrhs,B(k+1,1),Ldb)
                  CALL ZGEMV('Conjugate transpose',k-1,Nrhs,-ONE,B,Ldb, &
     &                       A(1,k+1),1,ONE,B(k+1,1),Ldb)
                  CALL ZLACGV(Nrhs,B(k+1,1),Ldb)
               ENDIF
!
!           Interchange rows K and -IPIV(K), then K+1 and -IPIV(K+1)
!
               kp = -Ipiv(k)
               IF ( kp/=k ) CALL ZSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
!
               kp = -Ipiv(k+1)
               IF ( kp/=k+1 ) CALL ZSWAP(Nrhs,B(k+1,1),Ldb,B(kp,1),Ldb)
!
               k = k + 2
!
            ENDIF
         ENDDO
!
      ELSE
!
!        Solve A*X = B, where A = L*D*L**H.
!
!        First solve L*D*X = B, overwriting B with X.
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
!           Interchange rows K and IPIV(K).
!
               kp = Ipiv(k)
               IF ( kp/=k ) CALL ZSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
!
!           Multiply by inv(L(K)), where L(K) is the transformation
!           stored in column K of A.
!
               IF ( k<N ) CALL ZGERU(N-k,Nrhs,-ONE,A(k+1,k),1,B(k,1),   &
     &                               Ldb,B(k+1,1),Ldb)
!
!           Multiply by the inverse of the diagonal block.
!
               s = DBLE(ONE)/DBLE(A(k,k))
               CALL ZDSCAL(Nrhs,s,B(k,1),Ldb)
               k = k + 1
            ELSE
!
!           2 x 2 diagonal block
!
!           Interchange rows K and -IPIV(K), then K+1 and -IPIV(K+1)
!
               kp = -Ipiv(k)
               IF ( kp/=k ) CALL ZSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
!
               kp = -Ipiv(k+1)
               IF ( kp/=k+1 ) CALL ZSWAP(Nrhs,B(k+1,1),Ldb,B(kp,1),Ldb)
!
!           Multiply by inv(L(K)), where L(K) is the transformation
!           stored in columns K and K+1 of A.
!
               IF ( k<N-1 ) THEN
                  CALL ZGERU(N-k-1,Nrhs,-ONE,A(k+2,k),1,B(k,1),Ldb,     &
     &                       B(k+2,1),Ldb)
                  CALL ZGERU(N-k-1,Nrhs,-ONE,A(k+2,k+1),1,B(k+1,1),Ldb, &
     &                       B(k+2,1),Ldb)
               ENDIF
!
!           Multiply by the inverse of the diagonal block.
!
               akm1k = A(k+1,k)
               akm1 = A(k,k)/DCONJG(akm1k)
               ak = A(k+1,k+1)/akm1k
               denom = akm1*ak - ONE
               DO j = 1 , Nrhs
                  bkm1 = B(k,j)/DCONJG(akm1k)
                  bk = B(k+1,j)/akm1k
                  B(k,j) = (ak*bkm1-bk)/denom
                  B(k+1,j) = (akm1*bk-bkm1)/denom
               ENDDO
               k = k + 2
!
            ENDIF
         ENDDO
!
!        Next solve L**H *X = B, overwriting B with X.
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
!           Multiply by inv(L**H(K)), where L(K) is the transformation
!           stored in column K of A.
!
               IF ( k<N ) THEN
                  CALL ZLACGV(Nrhs,B(k,1),Ldb)
                  CALL ZGEMV('Conjugate transpose',N-k,Nrhs,-ONE,       &
     &                       B(k+1,1),Ldb,A(k+1,k),1,ONE,B(k,1),Ldb)
                  CALL ZLACGV(Nrhs,B(k,1),Ldb)
               ENDIF
!
!           Interchange rows K and IPIV(K).
!
               kp = Ipiv(k)
               IF ( kp/=k ) CALL ZSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
               k = k - 1
            ELSE
!
!           2 x 2 diagonal block
!
!           Multiply by inv(L**H(K-1)), where L(K-1) is the transformation
!           stored in columns K-1 and K of A.
!
               IF ( k<N ) THEN
                  CALL ZLACGV(Nrhs,B(k,1),Ldb)
                  CALL ZGEMV('Conjugate transpose',N-k,Nrhs,-ONE,       &
     &                       B(k+1,1),Ldb,A(k+1,k),1,ONE,B(k,1),Ldb)
                  CALL ZLACGV(Nrhs,B(k,1),Ldb)
!
                  CALL ZLACGV(Nrhs,B(k-1,1),Ldb)
                  CALL ZGEMV('Conjugate transpose',N-k,Nrhs,-ONE,       &
     &                       B(k+1,1),Ldb,A(k+1,k-1),1,ONE,B(k-1,1),Ldb)
                  CALL ZLACGV(Nrhs,B(k-1,1),Ldb)
               ENDIF
!
!           Interchange rows K and -IPIV(K), then K-1 and -IPIV(K-1)
!
               kp = -Ipiv(k)
               IF ( kp/=k ) CALL ZSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
!
               kp = -Ipiv(k-1)
               IF ( kp/=k-1 ) CALL ZSWAP(Nrhs,B(k-1,1),Ldb,B(kp,1),Ldb)
!
               k = k - 2
!
            ENDIF
         ENDDO
      ENDIF
!
!
!     End of ZHETRS_ROOK
!
      END SUBROUTINE ZHETRS_ROOK
