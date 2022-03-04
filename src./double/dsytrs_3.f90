!*==dsytrs_3.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DSYTRS_3
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DSYTRS_3 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsytrs_3.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsytrs_3.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsytrs_3.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSYTRS_3( UPLO, N, NRHS, A, LDA, E, IPIV, B, LDB,
!                            INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDA, LDB, N, NRHS
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), E( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!> DSYTRS_3 solves a system of linear equations A * X = B with a real
!> symmetric matrix A using the factorization computed
!> by DSYTRF_RK or DSYTRF_BK:
!>
!>    A = P*U*D*(U**T)*(P**T) or A = P*L*D*(L**T)*(P**T),
!>
!> where U (or L) is unit upper (or lower) triangular matrix,
!> U**T (or L**T) is the transpose of U (or L), P is a permutation
!> matrix, P**T is the transpose of P, and D is symmetric and block
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
!>          = 'U':  Upper triangular, form is A = P*U*D*(U**T)*(P**T);
!>          = 'L':  Lower triangular, form is A = P*L*D*(L**T)*(P**T).
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
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          Diagonal of the block diagonal matrix D and factors U or L
!>          as computed by DSYTRF_RK and DSYTRF_BK:
!>            a) ONLY diagonal elements of the symmetric block diagonal
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
!>          E is DOUBLE PRECISION array, dimension (N)
!>          On entry, contains the superdiagonal (or subdiagonal)
!>          elements of the symmetric block diagonal matrix D
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
!>          as determined by DSYTRF_RK or DSYTRF_BK.
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)
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
!> \ingroup doubleSYcomputational
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
      SUBROUTINE DSYTRS_3(Uplo,N,Nrhs,A,Lda,E,Ipiv,B,Ldb,Info)
      USE F77KINDS                        
      USE S_DSCAL
      USE S_DSWAP
      USE S_DTRSM
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--DSYTRS_3174
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: E
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: ak , akm1 , akm1k , bk , bkm1 , denom
      INTEGER :: i , j , k , kp
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
         CALL XERBLA('DSYTRS_3',-Info)
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
!        Solve A*X = B, where A = U*D*U**T.
!
!        P**T * B
!
!        Interchange rows K and IPIV(K) of matrix B in the same order
!        that the formation order of IPIV(I) vector for Upper case.
!
!        (We can do the simple loop over IPIV with decrement -1,
!        since the ABS value of IPIV( I ) represents the row index
!        of the interchange with row i in both 1x1 and 2x2 pivot cases)
!
         DO k = N , 1 , -1
            kp = ABS(Ipiv(k))
            IF ( kp/=k ) CALL DSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
         ENDDO
!
!        Compute (U \P**T * B) -> B    [ (U \P**T * B) ]
!
         CALL DTRSM('L','U','N','U',N,Nrhs,ONE,A,Lda,B,Ldb)
!
!        Compute D \ B -> B   [ D \ (U \P**T * B) ]
!
         i = N
         DO WHILE ( i>=1 )
            IF ( Ipiv(i)>0 ) THEN
               CALL DSCAL(Nrhs,ONE/A(i,i),B(i,1),Ldb)
            ELSEIF ( i>1 ) THEN
               akm1k = E(i)
               akm1 = A(i-1,i-1)/akm1k
               ak = A(i,i)/akm1k
               denom = akm1*ak - ONE
               DO j = 1 , Nrhs
                  bkm1 = B(i-1,j)/akm1k
                  bk = B(i,j)/akm1k
                  B(i-1,j) = (ak*bkm1-bk)/denom
                  B(i,j) = (akm1*bk-bkm1)/denom
               ENDDO
               i = i - 1
            ENDIF
            i = i - 1
         ENDDO
!
!        Compute (U**T \ B) -> B   [ U**T \ (D \ (U \P**T * B) ) ]
!
         CALL DTRSM('L','U','T','U',N,Nrhs,ONE,A,Lda,B,Ldb)
!
!        P * B  [ P * (U**T \ (D \ (U \P**T * B) )) ]
!
!        Interchange rows K and IPIV(K) of matrix B in reverse order
!        from the formation order of IPIV(I) vector for Upper case.
!
!        (We can do the simple loop over IPIV with increment 1,
!        since the ABS value of IPIV(I) represents the row index
!        of the interchange with row i in both 1x1 and 2x2 pivot cases)
!
         DO k = 1 , N
            kp = ABS(Ipiv(k))
            IF ( kp/=k ) CALL DSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
         ENDDO
!
      ELSE
!
!        Begin Lower
!
!        Solve A*X = B, where A = L*D*L**T.
!
!        P**T * B
!        Interchange rows K and IPIV(K) of matrix B in the same order
!        that the formation order of IPIV(I) vector for Lower case.
!
!        (We can do the simple loop over IPIV with increment 1,
!        since the ABS value of IPIV(I) represents the row index
!        of the interchange with row i in both 1x1 and 2x2 pivot cases)
!
         DO k = 1 , N
            kp = ABS(Ipiv(k))
            IF ( kp/=k ) CALL DSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
         ENDDO
!
!        Compute (L \P**T * B) -> B    [ (L \P**T * B) ]
!
         CALL DTRSM('L','L','N','U',N,Nrhs,ONE,A,Lda,B,Ldb)
!
!        Compute D \ B -> B   [ D \ (L \P**T * B) ]
!
         i = 1
         DO WHILE ( i<=N )
            IF ( Ipiv(i)>0 ) THEN
               CALL DSCAL(Nrhs,ONE/A(i,i),B(i,1),Ldb)
            ELSEIF ( i<N ) THEN
               akm1k = E(i)
               akm1 = A(i,i)/akm1k
               ak = A(i+1,i+1)/akm1k
               denom = akm1*ak - ONE
               DO j = 1 , Nrhs
                  bkm1 = B(i,j)/akm1k
                  bk = B(i+1,j)/akm1k
                  B(i,j) = (ak*bkm1-bk)/denom
                  B(i+1,j) = (akm1*bk-bkm1)/denom
               ENDDO
               i = i + 1
            ENDIF
            i = i + 1
         ENDDO
!
!        Compute (L**T \ B) -> B   [ L**T \ (D \ (L \P**T * B) ) ]
!
         CALL DTRSM('L','L','T','U',N,Nrhs,ONE,A,Lda,B,Ldb)
!
!        P * B  [ P * (L**T \ (D \ (L \P**T * B) )) ]
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
            IF ( kp/=k ) CALL DSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
         ENDDO
!
!        END Lower
!
      ENDIF
!
!
!     End of DSYTRS_3
!
      END SUBROUTINE DSYTRS_3
