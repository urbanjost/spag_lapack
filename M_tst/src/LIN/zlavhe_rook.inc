!*==zlavhe_rook.f90  processed by SPAG 7.51RB at 17:35 on  4 Mar 2022
!> \brief \b zlavhe_rook
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE zlavhe_rook( UPLO, TRANS, DIAG, N, NRHS, A, LDA, IPIV, B,
!                               LDB, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIAG, TRANS, UPLO
!       INTEGER            INFO, LDA, LDB, N, NRHS
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX*16         A( LDA, * ), B( LDB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> zlavhe_rook performs one of the matrix-vector operations
!>    x := A*x  or  x := A^H*x,
!> where x is an N element vector and  A is one of the factors
!> from the block U*D*U' or L*D*L' factorization computed by ZHETRF_ROOK.
!>
!> If TRANS = 'N', multiplies by U  or U * D  (or L  or L * D)
!> If TRANS = 'C', multiplies by U' or D * U' (or L' or D * L')
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the factor stored in A is upper or lower
!>          triangular.
!>          = 'U':  Upper triangular
!>          = 'L':  Lower triangular
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          Specifies the operation to be performed:
!>          = 'N':  x := A*x
!>          = 'C':   x := A^H*x
!> \endverbatim
!>
!> \param[in] DIAG
!> \verbatim
!>          DIAG is CHARACTER*1
!>          Specifies whether or not the diagonal blocks are unit
!>          matrices.  If the diagonal blocks are assumed to be unit,
!>          then A = U or A = L, otherwise A = U*D or A = L*D.
!>          = 'U':  Diagonal blocks are assumed to be unit matrices.
!>          = 'N':  Diagonal blocks are assumed to be non-unit matrices.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of rows and columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of right hand sides, i.e., the number of vectors
!>          x to be multiplied by A.  NRHS >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          The block diagonal matrix D and the multipliers used to
!>          obtain the factor U or L as computed by ZHETRF_ROOK.
!>          Stored as a 2-D triangular matrix.
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
!>          IPIV is INTEGER array, dimension (N)
!>          Details of the interchanges and the block structure of D,
!>          as determined by ZHETRF_ROOK.
!>          If UPLO = 'U':
!>             Only the last KB elements of IPIV are set.
!>
!>             If IPIV(k) > 0, then rows and columns k and IPIV(k) were
!>             interchanged and D(k,k) is a 1-by-1 diagonal block.
!>
!>             If IPIV(k) < 0 and IPIV(k-1) < 0, then rows and
!>             columns k and -IPIV(k) were interchanged and rows and
!>             columns k-1 and -IPIV(k-1) were inerchaged,
!>             D(k-1:k,k-1:k) is a 2-by-2 diagonal block.
!>
!>          If UPLO = 'L':
!>             Only the first KB elements of IPIV are set.
!>
!>             If IPIV(k) > 0, then rows and columns k and IPIV(k)
!>             were interchanged and D(k,k) is a 1-by-1 diagonal block.
!>
!>             If IPIV(k) < 0 and IPIV(k+1) < 0, then rows and
!>             columns k and -IPIV(k) were interchanged and rows and
!>             columns k+1 and -IPIV(k+1) were inerchaged,
!>             D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is COMPLEX array, dimension (LDB,NRHS)
!>          On entry, B contains NRHS vectors of length N.
!>          On exit, B is overwritten with the product A * B.
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
!>          = 0: successful exit
!>          < 0: if INFO = -k, the k-th argument had an illegal value
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
!> \ingroup complex16_lin
!
!  =====================================================================
      SUBROUTINE zlavhe_rook(Uplo,Trans,Diag,N,Nrhs,A,Lda,Ipiv,B,Ldb,   &
     &                       Info)
      IMPLICIT NONE
!*--ZLAVHE_ROOK157
!
!  -- LAPACK test routine (version 3.5.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2013
!
!     .. Scalar Arguments ..
      CHARACTER Diag , Trans , Uplo
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
      COMPLEX*16 CONE
      PARAMETER (CONE=(1.0D+0,0.0D+0))
!     ..
!     .. Local Scalars ..
      LOGICAL nounit
      INTEGER j , k , kp
      COMPLEX*16 d11 , d12 , d21 , d22 , t1 , t2
!     ..
!     .. External Functions ..
      LOGICAL LSAME
EXTERNAL    LSAME
!     ..
!     .. External Subroutines ..
EXTERNAL    ZGEMV,   ZGERU,   ZLACGV,   ZSCAL,   ZSWAP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , DCONJG , MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
      IF ( .NOT.LSAME(Uplo,'U') .AND. .NOT.LSAME(Uplo,'L') ) THEN
         Info = -1
      ELSEIF ( .NOT.LSAME(Trans,'N') .AND. .NOT.LSAME(Trans,'C') ) THEN
         Info = -2
      ELSEIF ( .NOT.LSAME(Diag,'U') .AND. .NOT.LSAME(Diag,'N') ) THEN
         Info = -3
      ELSEIF ( N<0 ) THEN
         Info = -4
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -6
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -9
      ENDIF
      IF ( Info/=0 ) THEN
         CALL xerbla('zlavhe_rook ',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible.
!
      IF ( N==0 ) RETURN
!
      nounit = LSAME(Diag,'N')
!------------------------------------------
!
!     Compute  B := A * B  (No transpose)
!
!------------------------------------------
      IF ( LSAME(Trans,'N') ) THEN
!
!        Compute  B := U*B
!        where U = P(m)*inv(U(m))* ... *P(1)*inv(U(1))
!
         IF ( LSAME(Uplo,'U') ) THEN
!
!        Loop forward applying the transformations.
!
            k = 1
            DO WHILE ( k<=N )
               IF ( Ipiv(k)>0 ) THEN
!
!              1 x 1 pivot block
!
!              Multiply by the diagonal element if forming U * D.
!
                  IF ( nounit ) CALL ZSCAL(Nrhs,A(k,k),B(k,1),Ldb)
!
!              Multiply by  P(K) * inv(U(K))  if K > 1.
!
                  IF ( k>1 ) THEN
!
!                 Apply the transformation.
!
                     CALL ZGERU(k-1,Nrhs,CONE,A(1,k),1,B(k,1),Ldb,B(1,1)&
     &                          ,Ldb)
!
!                 Interchange if P(K) != I.
!
                     kp = Ipiv(k)
                     IF ( kp/=k ) CALL ZSWAP(Nrhs,B(k,1),Ldb,B(kp,1),   &
     &                    Ldb)
                  ENDIF
                  k = k + 1
               ELSE
!
!              2 x 2 pivot block
!
!              Multiply by the diagonal block if forming U * D.
!
                  IF ( nounit ) THEN
                     d11 = A(k,k)
                     d22 = A(k+1,k+1)
                     d12 = A(k,k+1)
                     d21 = DCONJG(d12)
                     DO j = 1 , Nrhs
                        t1 = B(k,j)
                        t2 = B(k+1,j)
                        B(k,j) = d11*t1 + d12*t2
                        B(k+1,j) = d21*t1 + d22*t2
                     ENDDO
                  ENDIF
!
!              Multiply by  P(K) * inv(U(K))  if K > 1.
!
                  IF ( k>1 ) THEN
!
!                 Apply the transformations.
!
                     CALL ZGERU(k-1,Nrhs,CONE,A(1,k),1,B(k,1),Ldb,B(1,1)&
     &                          ,Ldb)
                     CALL ZGERU(k-1,Nrhs,CONE,A(1,k+1),1,B(k+1,1),Ldb,  &
     &                          B(1,1),Ldb)
!
!                 Interchange if a permutation was applied at the
!                 K-th step of the factorization.
!
!                 Swap the first of pair with IMAXth
!
                     kp = ABS(Ipiv(k))
                     IF ( kp/=k ) CALL ZSWAP(Nrhs,B(k,1),Ldb,B(kp,1),   &
     &                    Ldb)
!
!                 NOW swap the first of pair with Pth
!
                     kp = ABS(Ipiv(k+1))
                     IF ( kp/=k+1 )                                     &
     &                    CALL ZSWAP(Nrhs,B(k+1,1),Ldb,B(kp,1),Ldb)
                  ENDIF
                  k = k + 2
               ENDIF
            ENDDO
!
!        Compute  B := L*B
!        where L = P(1)*inv(L(1))* ... *P(m)*inv(L(m)) .
!
         ELSE
!
!           Loop backward applying the transformations to B.
!
            k = N
            DO WHILE ( k>=1 )
!
!           Test the pivot index.  If greater than zero, a 1 x 1
!           pivot was used, otherwise a 2 x 2 pivot was used.
!
               IF ( Ipiv(k)>0 ) THEN
!
!              1 x 1 pivot block:
!
!              Multiply by the diagonal element if forming L * D.
!
                  IF ( nounit ) CALL ZSCAL(Nrhs,A(k,k),B(k,1),Ldb)
!
!              Multiply by  P(K) * inv(L(K))  if K < N.
!
                  IF ( k/=N ) THEN
                     kp = Ipiv(k)
!
!                 Apply the transformation.
!
                     CALL ZGERU(N-k,Nrhs,CONE,A(k+1,k),1,B(k,1),Ldb,    &
     &                          B(k+1,1),Ldb)
!
!                 Interchange if a permutation was applied at the
!                 K-th step of the factorization.
!
                     IF ( kp/=k ) CALL ZSWAP(Nrhs,B(k,1),Ldb,B(kp,1),   &
     &                    Ldb)
                  ENDIF
                  k = k - 1
!
               ELSE
!
!              2 x 2 pivot block:
!
!              Multiply by the diagonal block if forming L * D.
!
                  IF ( nounit ) THEN
                     d11 = A(k-1,k-1)
                     d22 = A(k,k)
                     d21 = A(k,k-1)
                     d12 = DCONJG(d21)
                     DO j = 1 , Nrhs
                        t1 = B(k-1,j)
                        t2 = B(k,j)
                        B(k-1,j) = d11*t1 + d12*t2
                        B(k,j) = d21*t1 + d22*t2
                     ENDDO
                  ENDIF
!
!              Multiply by  P(K) * inv(L(K))  if K < N.
!
                  IF ( k/=N ) THEN
!
!                 Apply the transformation.
!
                     CALL ZGERU(N-k,Nrhs,CONE,A(k+1,k),1,B(k,1),Ldb,    &
     &                          B(k+1,1),Ldb)
                     CALL ZGERU(N-k,Nrhs,CONE,A(k+1,k-1),1,B(k-1,1),Ldb,&
     &                          B(k+1,1),Ldb)
!
!                 Interchange if a permutation was applied at the
!                 K-th step of the factorization.
!
!
!                 Swap the second of pair with IMAXth
!
                     kp = ABS(Ipiv(k))
                     IF ( kp/=k ) CALL ZSWAP(Nrhs,B(k,1),Ldb,B(kp,1),   &
     &                    Ldb)
!
!                 NOW swap the first of pair with Pth
!
                     kp = ABS(Ipiv(k-1))
                     IF ( kp/=k-1 )                                     &
     &                    CALL ZSWAP(Nrhs,B(k-1,1),Ldb,B(kp,1),Ldb)
!
                  ENDIF
                  k = k - 2
               ENDIF
            ENDDO
         ENDIF
!--------------------------------------------------
!
!     Compute  B := A^H * B  (conjugate transpose)
!
!--------------------------------------------------
!
!        Form  B := U^H*B
!        where U  = P(m)*inv(U(m))* ... *P(1)*inv(U(1))
!        and   U^H = inv(U^H(1))*P(1)* ... *inv(U^H(m))*P(m)
!
      ELSEIF ( LSAME(Uplo,'U') ) THEN
!
!           Loop backward applying the transformations.
!
         k = N
         DO WHILE ( k>=1 )
!
!           1 x 1 pivot block.
!
            IF ( Ipiv(k)>0 ) THEN
               IF ( k>1 ) THEN
!
!                 Interchange if P(K) != I.
!
                  kp = Ipiv(k)
                  IF ( kp/=k ) CALL ZSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
!
!                 Apply the transformation
!                    y = y - B' DCONJG(x),
!                 where x is a column of A and y is a row of B.
!
                  CALL ZLACGV(Nrhs,B(k,1),Ldb)
                  CALL ZGEMV('Conjugate',k-1,Nrhs,CONE,B,Ldb,A(1,k),1,  &
     &                       CONE,B(k,1),Ldb)
                  CALL ZLACGV(Nrhs,B(k,1),Ldb)
               ENDIF
               IF ( nounit ) CALL ZSCAL(Nrhs,A(k,k),B(k,1),Ldb)
               k = k - 1
!
!           2 x 2 pivot block.
!
            ELSE
               IF ( k>2 ) THEN
!
!                 Swap the second of pair with Pth
!
                  kp = ABS(Ipiv(k))
                  IF ( kp/=k ) CALL ZSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
!
!                 Now swap the first of pair with IMAX(r)th
!
                  kp = ABS(Ipiv(k-1))
                  IF ( kp/=k-1 ) CALL ZSWAP(Nrhs,B(k-1,1),Ldb,B(kp,1),  &
     &                 Ldb)
!
!                 Apply the transformations
!                    y = y - B' DCONJG(x),
!                 where x is a block column of A and y is a block
!                 row of B.
!
                  CALL ZLACGV(Nrhs,B(k,1),Ldb)
                  CALL ZGEMV('Conjugate',k-2,Nrhs,CONE,B,Ldb,A(1,k),1,  &
     &                       CONE,B(k,1),Ldb)
                  CALL ZLACGV(Nrhs,B(k,1),Ldb)
!
                  CALL ZLACGV(Nrhs,B(k-1,1),Ldb)
                  CALL ZGEMV('Conjugate',k-2,Nrhs,CONE,B,Ldb,A(1,k-1),1,&
     &                       CONE,B(k-1,1),Ldb)
                  CALL ZLACGV(Nrhs,B(k-1,1),Ldb)
               ENDIF
!
!              Multiply by the diagonal block if non-unit.
!
               IF ( nounit ) THEN
                  d11 = A(k-1,k-1)
                  d22 = A(k,k)
                  d12 = A(k-1,k)
                  d21 = DCONJG(d12)
                  DO j = 1 , Nrhs
                     t1 = B(k-1,j)
                     t2 = B(k,j)
                     B(k-1,j) = d11*t1 + d12*t2
                     B(k,j) = d21*t1 + d22*t2
                  ENDDO
               ENDIF
               k = k - 2
            ENDIF
         ENDDO
!
!        Form  B := L^H*B
!        where L  = P(1)*inv(L(1))* ... *P(m)*inv(L(m))
!        and   L^H = inv(L^H(m))*P(m)* ... *inv(L^H(1))*P(1)
!
      ELSE
!
!           Loop forward applying the L-transformations.
!
         k = 1
         DO WHILE ( k<=N )
!
!           1 x 1 pivot block
!
            IF ( Ipiv(k)>0 ) THEN
               IF ( k<N ) THEN
!
!                 Interchange if P(K) != I.
!
                  kp = Ipiv(k)
                  IF ( kp/=k ) CALL ZSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
!
!                 Apply the transformation
!
                  CALL ZLACGV(Nrhs,B(k,1),Ldb)
                  CALL ZGEMV('Conjugate',N-k,Nrhs,CONE,B(k+1,1),Ldb,    &
     &                       A(k+1,k),1,CONE,B(k,1),Ldb)
                  CALL ZLACGV(Nrhs,B(k,1),Ldb)
               ENDIF
               IF ( nounit ) CALL ZSCAL(Nrhs,A(k,k),B(k,1),Ldb)
               k = k + 1
!
!           2 x 2 pivot block.
!
            ELSE
               IF ( k<N-1 ) THEN
!
!                 Swap the first of pair with Pth
!
                  kp = ABS(Ipiv(k))
                  IF ( kp/=k ) CALL ZSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
!
!                 Now swap the second of pair with IMAX(r)th
!
                  kp = ABS(Ipiv(k+1))
                  IF ( kp/=k+1 ) CALL ZSWAP(Nrhs,B(k+1,1),Ldb,B(kp,1),  &
     &                 Ldb)
!
!                 Apply the transformation
!
                  CALL ZLACGV(Nrhs,B(k+1,1),Ldb)
                  CALL ZGEMV('Conjugate',N-k-1,Nrhs,CONE,B(k+2,1),Ldb,  &
     &                       A(k+2,k+1),1,CONE,B(k+1,1),Ldb)
                  CALL ZLACGV(Nrhs,B(k+1,1),Ldb)
!
                  CALL ZLACGV(Nrhs,B(k,1),Ldb)
                  CALL ZGEMV('Conjugate',N-k-1,Nrhs,CONE,B(k+2,1),Ldb,  &
     &                       A(k+2,k),1,CONE,B(k,1),Ldb)
                  CALL ZLACGV(Nrhs,B(k,1),Ldb)
               ENDIF
!
!              Multiply by the diagonal block if non-unit.
!
               IF ( nounit ) THEN
                  d11 = A(k,k)
                  d22 = A(k+1,k+1)
                  d21 = A(k+1,k)
                  d12 = DCONJG(d21)
                  DO j = 1 , Nrhs
                     t1 = B(k,j)
                     t2 = B(k+1,j)
                     B(k,j) = d11*t1 + d12*t2
                     B(k+1,j) = d21*t1 + d22*t2
                  ENDDO
               ENDIF
               k = k + 2
            ENDIF
         ENDDO
!
      ENDIF
!
!     End of zlavhe_rook
!
      END SUBROUTINE zlavhe_rook
