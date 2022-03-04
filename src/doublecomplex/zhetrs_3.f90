!*==zhetrs_3.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b ZHETRS_3
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZHETRS_3 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhetrs_3.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhetrs_3.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhetrs_3.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZHETRS_3( UPLO, N, NRHS, A, LDA, E, IPIV, B, LDB,
!                            INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDA, LDB, N, NRHS
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX*16         A( LDA, * ), B( LDB, * ), E( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!> ZHETRS_3 solves a system of linear equations A * X = B with a complex
!> Hermitian matrix A using the factorization computed
!> by ZHETRF_RK or ZHETRF_BK:
!>
!>    A = P*U*D*(U**H)*(P**T) or A = P*L*D*(L**H)*(P**T),
!>
!> where U (or L) is unit upper (or lower) triangular matrix,
!> U**H (or L**H) is the conjugate of U (or L), P is a permutation
!> matrix, P**T is the transpose of P, and D is Hermitian and block
!> diagonal with 1-by-1 and 2-by-2 diagonal blocks.
!>
!> This algorithm is using Level 3 BLAS.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the details of the factorization are
!>          stored as an upper or lower triangular matrix:
!>          = 'U':  Upper triangular, form is A = P*U*D*(U**H)*(P**T);
!>          = 'L':  Lower triangular, form is A = P*L*D*(L**H)*(P**T).
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
!>          Diagonal of the block diagonal matrix D and factors U or L
!>          as computed by ZHETRF_RK and ZHETRF_BK:
!>            a) ONLY diagonal elements of the Hermitian block diagonal
!>               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k);
!>               (superdiagonal (or subdiagonal) elements of D
!>                should be provided on entry in array E), and
!>            b) If UPLO = 'U': factor U in the superdiagonal part of A.
!>               If UPLO = 'L': factor L in the subdiagonal part of A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is COMPLEX*16 array, dimension (N)
!>          On entry, contains the superdiagonal (or subdiagonal)
!>          elements of the Hermitian block diagonal matrix D
!>          with 1-by-1 or 2-by-2 diagonal blocks, where
!>          If UPLO = 'U': E(i) = D(i-1,i),i=2:N, E(1) not referenced;
!>          If UPLO = 'L': E(i) = D(i+1,i),i=1:N-1, E(N) not referenced.
!>
!>          NOTE: For 1-by-1 diagonal block D(k), where
!>          1 <= k <= N, the element E(k) is not referenced in both
!>          UPLO = 'U' or UPLO = 'L' cases.
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          Details of the interchanges and the block structure of D
!>          as determined by ZHETRF_RK or ZHETRF_BK.
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
!> \date June 2017
!
!> \ingroup complex16HEcomputational
!
!> \par Contributors:
!  ==================
!>
!> \verbatim
!>
!>  June 2017,  Igor Kozachenko,
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
      SUBROUTINE ZHETRS_3(Uplo,N,Nrhs,A,Lda,E,Ipiv,B,Ldb,Info)
      IMPLICIT NONE
!*--ZHETRS_3168
!
!  -- LAPACK computational routine (version 3.7.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2017
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER Info , Lda , Ldb , N , Nrhs
!     ..
!     .. Array Arguments ..
      INTEGER Ipiv(*)
      COMPLEX*16 A(Lda,*) , B(Ldb,*) , E(*)
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
      INTEGER i , j , k , kp
      DOUBLE PRECISION s
      COMPLEX*16 ak , akm1 , akm1k , bk , bkm1 , denom
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL ZDSCAL , ZSWAP , ZTRSM , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , DBLE , DCONJG , MAX
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
         Info = -9
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZHETRS_3',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 .OR. Nrhs==0 ) RETURN
!
      IF ( upper ) THEN
!
!        Begin Upper
!
!        Solve A*X = B, where A = U*D*U**H.
!
!        P**T * B
!
!        Interchange rows K and IPIV(K) of matrix B in the same order
!        that the formation order of IPIV(I) vector for Upper case.
!
!        (We can do the simple loop over IPIV with decrement -1,
!        since the ABS value of IPIV(I) represents the row index
!        of the interchange with row i in both 1x1 and 2x2 pivot cases)
!
         DO k = N , 1 , -1
            kp = ABS(Ipiv(k))
            IF ( kp/=k ) CALL ZSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
         ENDDO
!
!        Compute (U \P**T * B) -> B    [ (U \P**T * B) ]
!
         CALL ZTRSM('L','U','N','U',N,Nrhs,ONE,A,Lda,B,Ldb)
!
!        Compute D \ B -> B   [ D \ (U \P**T * B) ]
!
         i = N
         DO WHILE ( i>=1 )
            IF ( Ipiv(i)>0 ) THEN
               s = DBLE(ONE)/DBLE(A(i,i))
               CALL ZDSCAL(Nrhs,s,B(i,1),Ldb)
            ELSEIF ( i>1 ) THEN
               akm1k = E(i)
               akm1 = A(i-1,i-1)/akm1k
               ak = A(i,i)/DCONJG(akm1k)
               denom = akm1*ak - ONE
               DO j = 1 , Nrhs
                  bkm1 = B(i-1,j)/akm1k
                  bk = B(i,j)/DCONJG(akm1k)
                  B(i-1,j) = (ak*bkm1-bk)/denom
                  B(i,j) = (akm1*bk-bkm1)/denom
               ENDDO
               i = i - 1
            ENDIF
            i = i - 1
         ENDDO
!
!        Compute (U**H \ B) -> B   [ U**H \ (D \ (U \P**T * B) ) ]
!
         CALL ZTRSM('L','U','C','U',N,Nrhs,ONE,A,Lda,B,Ldb)
!
!        P * B  [ P * (U**H \ (D \ (U \P**T * B) )) ]
!
!        Interchange rows K and IPIV(K) of matrix B in reverse order
!        from the formation order of IPIV(I) vector for Upper case.
!
!        (We can do the simple loop over IPIV with increment 1,
!        since the ABS value of IPIV(I) represents the row index
!        of the interchange with row i in both 1x1 and 2x2 pivot cases)
!
         DO k = 1 , N , 1
            kp = ABS(Ipiv(k))
            IF ( kp/=k ) CALL ZSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
         ENDDO
!
      ELSE
!
!        Begin Lower
!
!        Solve A*X = B, where A = L*D*L**H.
!
!        P**T * B
!        Interchange rows K and IPIV(K) of matrix B in the same order
!        that the formation order of IPIV(I) vector for Lower case.
!
!        (We can do the simple loop over IPIV with increment 1,
!        since the ABS value of IPIV(I) represents the row index
!        of the interchange with row i in both 1x1 and 2x2 pivot cases)
!
         DO k = 1 , N , 1
            kp = ABS(Ipiv(k))
            IF ( kp/=k ) CALL ZSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
         ENDDO
!
!        Compute (L \P**T * B) -> B    [ (L \P**T * B) ]
!
         CALL ZTRSM('L','L','N','U',N,Nrhs,ONE,A,Lda,B,Ldb)
!
!        Compute D \ B -> B   [ D \ (L \P**T * B) ]
!
         i = 1
         DO WHILE ( i<=N )
            IF ( Ipiv(i)>0 ) THEN
               s = DBLE(ONE)/DBLE(A(i,i))
               CALL ZDSCAL(Nrhs,s,B(i,1),Ldb)
            ELSEIF ( i<N ) THEN
               akm1k = E(i)
               akm1 = A(i,i)/DCONJG(akm1k)
               ak = A(i+1,i+1)/akm1k
               denom = akm1*ak - ONE
               DO j = 1 , Nrhs
                  bkm1 = B(i,j)/DCONJG(akm1k)
                  bk = B(i+1,j)/akm1k
                  B(i,j) = (ak*bkm1-bk)/denom
                  B(i+1,j) = (akm1*bk-bkm1)/denom
               ENDDO
               i = i + 1
            ENDIF
            i = i + 1
         ENDDO
!
!        Compute (L**H \ B) -> B   [ L**H \ (D \ (L \P**T * B) ) ]
!
         CALL ZTRSM('L','L','C','U',N,Nrhs,ONE,A,Lda,B,Ldb)
!
!        P * B  [ P * (L**H \ (D \ (L \P**T * B) )) ]
!
!        Interchange rows K and IPIV(K) of matrix B in reverse order
!        from the formation order of IPIV(I) vector for Lower case.
!
!        (We can do the simple loop over IPIV with decrement -1,
!        since the ABS value of IPIV(I) represents the row index
!        of the interchange with row i in both 1x1 and 2x2 pivot cases)
!
         DO k = N , 1 , -1
            kp = ABS(Ipiv(k))
            IF ( kp/=k ) CALL ZSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
         ENDDO
!
!        END Lower
!
      ENDIF
!
!
!     End of ZHETRS_3
!
      END SUBROUTINE ZHETRS_3
