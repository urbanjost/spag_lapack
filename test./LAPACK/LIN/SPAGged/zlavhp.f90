!*==zlavhp.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b ZLAVHP
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLAVHP( UPLO, TRANS, DIAG, N, NRHS, A, IPIV, B, LDB,
!                          INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIAG, TRANS, UPLO
!       INTEGER            INFO, LDB, N, NRHS
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX*16         A( * ), B( LDB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    ZLAVHP  performs one of the matrix-vector operations
!>       x := A*x  or  x := A^H*x,
!>    where x is an N element vector and  A is one of the factors
!>    from the symmetric factorization computed by ZHPTRF.
!>    ZHPTRF produces a factorization of the form
!>         U * D * U^H     or     L * D * L^H,
!>    where U (or L) is a product of permutation and unit upper (lower)
!>    triangular matrices, U^H (or L^H) is the conjugate transpose of
!>    U (or L), and D is Hermitian and block diagonal with 1 x 1 and
!>    2 x 2 diagonal blocks.  The multipliers for the transformations
!>    and the upper or lower triangular parts of the diagonal blocks
!>    are stored columnwise in packed format in the linear array A.
!>
!>    If TRANS = 'N' or 'n', ZLAVHP multiplies either by U or U * D
!>    (or L or L * D).
!>    If TRANS = 'C' or 'c', ZLAVHP multiplies either by U^H or D * U^H
!>    (or L^H or D * L^H ).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \verbatim
!>  UPLO   - CHARACTER*1
!>           On entry, UPLO specifies whether the triangular matrix
!>           stored in A is upper or lower triangular.
!>              UPLO = 'U' or 'u'   The matrix is upper triangular.
!>              UPLO = 'L' or 'l'   The matrix is lower triangular.
!>           Unchanged on exit.
!>
!>  TRANS  - CHARACTER*1
!>           On entry, TRANS specifies the operation to be performed as
!>           follows:
!>              TRANS = 'N' or 'n'   x := A*x.
!>              TRANS = 'C' or 'c'   x := A^H*x.
!>           Unchanged on exit.
!>
!>  DIAG   - CHARACTER*1
!>           On entry, DIAG specifies whether the diagonal blocks are
!>           assumed to be unit matrices, as follows:
!>              DIAG = 'U' or 'u'   Diagonal blocks are unit matrices.
!>              DIAG = 'N' or 'n'   Diagonal blocks are non-unit.
!>           Unchanged on exit.
!>
!>  N      - INTEGER
!>           On entry, N specifies the order of the matrix A.
!>           N must be at least zero.
!>           Unchanged on exit.
!>
!>  NRHS   - INTEGER
!>           On entry, NRHS specifies the number of right hand sides,
!>           i.e., the number of vectors x to be multiplied by A.
!>           NRHS must be at least zero.
!>           Unchanged on exit.
!>
!>  A      - COMPLEX*16 array, dimension( N*(N+1)/2 )
!>           On entry, A contains a block diagonal matrix and the
!>           multipliers of the transformations used to obtain it,
!>           stored as a packed triangular matrix.
!>           Unchanged on exit.
!>
!>  IPIV   - INTEGER array, dimension( N )
!>           On entry, IPIV contains the vector of pivot indices as
!>           determined by ZSPTRF or ZHPTRF.
!>           If IPIV( K ) = K, no interchange was done.
!>           If IPIV( K ) <> K but IPIV( K ) > 0, then row K was inter-
!>           changed with row IPIV( K ) and a 1 x 1 pivot block was used.
!>           If IPIV( K ) < 0 and UPLO = 'U', then row K-1 was exchanged
!>           with row | IPIV( K ) | and a 2 x 2 pivot block was used.
!>           If IPIV( K ) < 0 and UPLO = 'L', then row K+1 was exchanged
!>           with row | IPIV( K ) | and a 2 x 2 pivot block was used.
!>
!>  B      - COMPLEX*16 array, dimension( LDB, NRHS )
!>           On entry, B contains NRHS vectors of length N.
!>           On exit, B is overwritten with the product A * B.
!>
!>  LDB    - INTEGER
!>           On entry, LDB contains the leading dimension of B as
!>           declared in the calling program.  LDB must be at least
!>           max( 1, N ).
!>           Unchanged on exit.
!>
!>  INFO   - INTEGER
!>           INFO is the error flag.
!>           On exit, a value of 0 indicates a successful exit.
!>           A negative value, say -K, indicates that the K-th argument
!>           has an illegal value.
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
!> \ingroup complex16_lin
!
!  =====================================================================
      SUBROUTINE ZLAVHP(Uplo,Trans,Diag,N,Nrhs,A,Ipiv,B,Ldb,Info)
      IMPLICIT NONE
!*--ZLAVHP134
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Diag , Trans , Uplo
      INTEGER Info , Ldb , N , Nrhs
!     ..
!     .. Array Arguments ..
      INTEGER Ipiv(*)
      COMPLEX*16 A(*) , B(Ldb,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16 ONE
      PARAMETER (ONE=(1.0D+0,0.0D+0))
!     ..
!     .. Local Scalars ..
      LOGICAL nounit
      INTEGER j , k , kc , kcnext , kp
      COMPLEX*16 d11 , d12 , d21 , d22 , t1 , t2
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA , ZGEMV , ZGERU , ZLACGV , ZSCAL , ZSWAP
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
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -8
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZLAVHP ',-Info)
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
            kc = 1
            DO WHILE ( k<=N )
!
!           1 x 1 pivot block
!
               IF ( Ipiv(k)>0 ) THEN
!
!              Multiply by the diagonal element if forming U * D.
!
                  IF ( nounit ) CALL ZSCAL(Nrhs,A(kc+k-1),B(k,1),Ldb)
!
!              Multiply by P(K) * inv(U(K))  if K > 1.
!
                  IF ( k>1 ) THEN
!
!                 Apply the transformation.
!
                     CALL ZGERU(k-1,Nrhs,ONE,A(kc),1,B(k,1),Ldb,B(1,1), &
     &                          Ldb)
!
!                 Interchange if P(K) != I.
!
                     kp = Ipiv(k)
                     IF ( kp/=k ) CALL ZSWAP(Nrhs,B(k,1),Ldb,B(kp,1),   &
     &                    Ldb)
                  ENDIF
                  kc = kc + k
                  k = k + 1
               ELSE
!
!              2 x 2 pivot block
!
                  kcnext = kc + k
!
!              Multiply by the diagonal block if forming U * D.
!
                  IF ( nounit ) THEN
                     d11 = A(kcnext-1)
                     d22 = A(kcnext+k)
                     d12 = A(kcnext+k-1)
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
                     CALL ZGERU(k-1,Nrhs,ONE,A(kc),1,B(k,1),Ldb,B(1,1), &
     &                          Ldb)
                     CALL ZGERU(k-1,Nrhs,ONE,A(kcnext),1,B(k+1,1),Ldb,  &
     &                          B(1,1),Ldb)
!
!                 Interchange if P(K) != I.
!
                     kp = ABS(Ipiv(k))
                     IF ( kp/=k ) CALL ZSWAP(Nrhs,B(k,1),Ldb,B(kp,1),   &
     &                    Ldb)
                  ENDIF
                  kc = kcnext + k + 1
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
            kc = N*(N+1)/2 + 1
            DO WHILE ( k>=1 )
               kc = kc - (N-k+1)
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
                  IF ( nounit ) CALL ZSCAL(Nrhs,A(kc),B(k,1),Ldb)
!
!              Multiply by  P(K) * inv(L(K))  if K < N.
!
                  IF ( k/=N ) THEN
                     kp = Ipiv(k)
!
!                 Apply the transformation.
!
                     CALL ZGERU(N-k,Nrhs,ONE,A(kc+1),1,B(k,1),Ldb,      &
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
                  kcnext = kc - (N-k+2)
!
!              Multiply by the diagonal block if forming L * D.
!
                  IF ( nounit ) THEN
                     d11 = A(kcnext)
                     d22 = A(kc)
                     d21 = A(kcnext+1)
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
                     CALL ZGERU(N-k,Nrhs,ONE,A(kc+1),1,B(k,1),Ldb,      &
     &                          B(k+1,1),Ldb)
                     CALL ZGERU(N-k,Nrhs,ONE,A(kcnext+2),1,B(k-1,1),Ldb,&
     &                          B(k+1,1),Ldb)
!
!                 Interchange if a permutation was applied at the
!                 K-th step of the factorization.
!
                     kp = ABS(Ipiv(k))
                     IF ( kp/=k ) CALL ZSWAP(Nrhs,B(k,1),Ldb,B(kp,1),   &
     &                    Ldb)
                  ENDIF
                  kc = kcnext
                  k = k - 2
               ENDIF
            ENDDO
         ENDIF
!-------------------------------------------------
!
!     Compute  B := A^H * B  (conjugate transpose)
!
!-------------------------------------------------
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
         kc = N*(N+1)/2 + 1
         DO WHILE ( k>=1 )
            kc = kc - k
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
!                 Apply the transformation:
!                    y := y - B' * conjg(x)
!                 where x is a column of A and y is a row of B.
!
                  CALL ZLACGV(Nrhs,B(k,1),Ldb)
                  CALL ZGEMV('Conjugate',k-1,Nrhs,ONE,B,Ldb,A(kc),1,ONE,&
     &                       B(k,1),Ldb)
                  CALL ZLACGV(Nrhs,B(k,1),Ldb)
               ENDIF
               IF ( nounit ) CALL ZSCAL(Nrhs,A(kc+k-1),B(k,1),Ldb)
               k = k - 1
!
!           2 x 2 pivot block.
!
            ELSE
               kcnext = kc - (k-1)
               IF ( k>2 ) THEN
!
!                 Interchange if P(K) != I.
!
                  kp = ABS(Ipiv(k))
                  IF ( kp/=k-1 ) CALL ZSWAP(Nrhs,B(k-1,1),Ldb,B(kp,1),  &
     &                 Ldb)
!
!                 Apply the transformations.
!
                  CALL ZLACGV(Nrhs,B(k,1),Ldb)
                  CALL ZGEMV('Conjugate',k-2,Nrhs,ONE,B,Ldb,A(kc),1,ONE,&
     &                       B(k,1),Ldb)
                  CALL ZLACGV(Nrhs,B(k,1),Ldb)
!
                  CALL ZLACGV(Nrhs,B(k-1,1),Ldb)
                  CALL ZGEMV('Conjugate',k-2,Nrhs,ONE,B,Ldb,A(kcnext),1,&
     &                       ONE,B(k-1,1),Ldb)
                  CALL ZLACGV(Nrhs,B(k-1,1),Ldb)
               ENDIF
!
!              Multiply by the diagonal block if non-unit.
!
               IF ( nounit ) THEN
                  d11 = A(kc-1)
                  d22 = A(kc+k-1)
                  d12 = A(kc+k-2)
                  d21 = DCONJG(d12)
                  DO j = 1 , Nrhs
                     t1 = B(k-1,j)
                     t2 = B(k,j)
                     B(k-1,j) = d11*t1 + d12*t2
                     B(k,j) = d21*t1 + d22*t2
                  ENDDO
               ENDIF
               kc = kcnext
               k = k - 2
            ENDIF
         ENDDO
!
!        Form  B := L^H*B
!        where L  = P(1)*inv(L(1))* ... *P(m)*inv(L(m))
!        and   L^H = inv(L(m))*P(m)* ... *inv(L(1))*P(1)
!
      ELSE
!
!           Loop forward applying the L-transformations.
!
         k = 1
         kc = 1
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
                  CALL ZGEMV('Conjugate',N-k,Nrhs,ONE,B(k+1,1),Ldb,     &
     &                       A(kc+1),1,ONE,B(k,1),Ldb)
                  CALL ZLACGV(Nrhs,B(k,1),Ldb)
               ENDIF
               IF ( nounit ) CALL ZSCAL(Nrhs,A(kc),B(k,1),Ldb)
               kc = kc + N - k + 1
               k = k + 1
!
!           2 x 2 pivot block.
!
            ELSE
               kcnext = kc + N - k + 1
               IF ( k<N-1 ) THEN
!
!              Interchange if P(K) != I.
!
                  kp = ABS(Ipiv(k))
                  IF ( kp/=k+1 ) CALL ZSWAP(Nrhs,B(k+1,1),Ldb,B(kp,1),  &
     &                 Ldb)
!
!                 Apply the transformation
!
                  CALL ZLACGV(Nrhs,B(k+1,1),Ldb)
                  CALL ZGEMV('Conjugate',N-k-1,Nrhs,ONE,B(k+2,1),Ldb,   &
     &                       A(kcnext+1),1,ONE,B(k+1,1),Ldb)
                  CALL ZLACGV(Nrhs,B(k+1,1),Ldb)
!
                  CALL ZLACGV(Nrhs,B(k,1),Ldb)
                  CALL ZGEMV('Conjugate',N-k-1,Nrhs,ONE,B(k+2,1),Ldb,   &
     &                       A(kc+2),1,ONE,B(k,1),Ldb)
                  CALL ZLACGV(Nrhs,B(k,1),Ldb)
               ENDIF
!
!              Multiply by the diagonal block if non-unit.
!
               IF ( nounit ) THEN
                  d11 = A(kc)
                  d22 = A(kcnext)
                  d21 = A(kc+1)
                  d12 = DCONJG(d21)
                  DO j = 1 , Nrhs
                     t1 = B(k,j)
                     t2 = B(k+1,j)
                     B(k,j) = d11*t1 + d12*t2
                     B(k+1,j) = d21*t1 + d22*t2
                  ENDDO
               ENDIF
               kc = kcnext + (N-k)
               k = k + 2
            ENDIF
         ENDDO
!
      ENDIF
!
!     End of ZLAVHP
!
      END SUBROUTINE ZLAVHP
