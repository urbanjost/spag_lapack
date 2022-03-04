!*==slavsp.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b SLAVSP
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLAVSP( UPLO, TRANS, DIAG, N, NRHS, A, IPIV, B, LDB,
!                          INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIAG, TRANS, UPLO
!       INTEGER            INFO, LDB, N, NRHS
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       REAL               A( * ), B( LDB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLAVSP  performs one of the matrix-vector operations
!>    x := A*x  or  x := A'*x,
!> where x is an N element vector and  A is one of the factors
!> from the block U*D*U' or L*D*L' factorization computed by SSPTRF.
!>
!> If TRANS = 'N', multiplies by U  or U * D  (or L  or L * D)
!> If TRANS = 'T', multiplies by U' or D * U' (or L' or D * L' )
!> If TRANS = 'C', multiplies by U' or D * U' (or L' or D * L' )
!> \endverbatim
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
!>          = 'T':  x := A'*x
!>          = 'C':  x := A'*x
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
!>          A is REAL array, dimension (N*(N+1)/2)
!>          The block diagonal matrix D and the multipliers used to
!>          obtain the factor U or L, stored as a packed triangular
!>          matrix as computed by SSPTRF.
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          The pivot indices from SSPTRF.
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is REAL array, dimension (LDB,NRHS)
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
!> \date December 2016
!
!> \ingroup single_lin
!
!  =====================================================================
      SUBROUTINE SLAVSP(Uplo,Trans,Diag,N,Nrhs,A,Ipiv,B,Ldb,Info)
      IMPLICIT NONE
!*--SLAVSP133
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
      REAL A(*) , B(Ldb,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ONE
      PARAMETER (ONE=1.0E+0)
!     ..
!     .. Local Scalars ..
      LOGICAL nounit
      INTEGER j , k , kc , kcnext , kp
      REAL d11 , d12 , d21 , d22 , t1 , t2
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL SGEMV , SGER , SSCAL , SSWAP , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
      IF ( .NOT.LSAME(Uplo,'U') .AND. .NOT.LSAME(Uplo,'L') ) THEN
         Info = -1
      ELSEIF ( .NOT.LSAME(Trans,'N') .AND. .NOT.LSAME(Trans,'T') .AND.  &
     &         .NOT.LSAME(Trans,'C') ) THEN
         Info = -2
      ELSEIF ( .NOT.LSAME(Diag,'U') .AND. .NOT.LSAME(Diag,'N') ) THEN
         Info = -3
      ELSEIF ( N<0 ) THEN
         Info = -4
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -8
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('SLAVSP ',-Info)
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
                  IF ( nounit ) CALL SSCAL(Nrhs,A(kc+k-1),B(k,1),Ldb)
!
!              Multiply by P(K) * inv(U(K))  if K > 1.
!
                  IF ( k>1 ) THEN
!
!                 Apply the transformation.
!
                     CALL SGER(k-1,Nrhs,ONE,A(kc),1,B(k,1),Ldb,B(1,1),  &
     &                         Ldb)
!
!                 Interchange if P(K) != I.
!
                     kp = Ipiv(k)
                     IF ( kp/=k ) CALL SSWAP(Nrhs,B(k,1),Ldb,B(kp,1),   &
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
                     d21 = d12
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
                     CALL SGER(k-1,Nrhs,ONE,A(kc),1,B(k,1),Ldb,B(1,1),  &
     &                         Ldb)
                     CALL SGER(k-1,Nrhs,ONE,A(kcnext),1,B(k+1,1),Ldb,   &
     &                         B(1,1),Ldb)
!
!                 Interchange if P(K) != I.
!
                     kp = ABS(Ipiv(k))
                     IF ( kp/=k ) CALL SSWAP(Nrhs,B(k,1),Ldb,B(kp,1),   &
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
                  IF ( nounit ) CALL SSCAL(Nrhs,A(kc),B(k,1),Ldb)
!
!              Multiply by  P(K) * inv(L(K))  if K < N.
!
                  IF ( k/=N ) THEN
                     kp = Ipiv(k)
!
!                 Apply the transformation.
!
                     CALL SGER(N-k,Nrhs,ONE,A(kc+1),1,B(k,1),Ldb,       &
     &                         B(k+1,1),Ldb)
!
!                 Interchange if a permutation was applied at the
!                 K-th step of the factorization.
!
                     IF ( kp/=k ) CALL SSWAP(Nrhs,B(k,1),Ldb,B(kp,1),   &
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
                     d12 = d21
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
                     CALL SGER(N-k,Nrhs,ONE,A(kc+1),1,B(k,1),Ldb,       &
     &                         B(k+1,1),Ldb)
                     CALL SGER(N-k,Nrhs,ONE,A(kcnext+2),1,B(k-1,1),Ldb, &
     &                         B(k+1,1),Ldb)
!
!                 Interchange if a permutation was applied at the
!                 K-th step of the factorization.
!
                     kp = ABS(Ipiv(k))
                     IF ( kp/=k ) CALL SSWAP(Nrhs,B(k,1),Ldb,B(kp,1),   &
     &                    Ldb)
                  ENDIF
                  kc = kcnext
                  k = k - 2
               ENDIF
            ENDDO
         ENDIF
!----------------------------------------
!
!     Compute  B := A' * B  (transpose)
!
!----------------------------------------
!
!        Form  B := U'*B
!        where U  = P(m)*inv(U(m))* ... *P(1)*inv(U(1))
!        and   U' = inv(U'(1))*P(1)* ... *inv(U'(m))*P(m)
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
                  IF ( kp/=k ) CALL SSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
!
!                 Apply the transformation
!
                  CALL SGEMV('Transpose',k-1,Nrhs,ONE,B,Ldb,A(kc),1,ONE,&
     &                       B(k,1),Ldb)
               ENDIF
               IF ( nounit ) CALL SSCAL(Nrhs,A(kc+k-1),B(k,1),Ldb)
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
                  IF ( kp/=k-1 ) CALL SSWAP(Nrhs,B(k-1,1),Ldb,B(kp,1),  &
     &                 Ldb)
!
!                 Apply the transformations
!
                  CALL SGEMV('Transpose',k-2,Nrhs,ONE,B,Ldb,A(kc),1,ONE,&
     &                       B(k,1),Ldb)
                  CALL SGEMV('Transpose',k-2,Nrhs,ONE,B,Ldb,A(kcnext),1,&
     &                       ONE,B(k-1,1),Ldb)
               ENDIF
!
!              Multiply by the diagonal block if non-unit.
!
               IF ( nounit ) THEN
                  d11 = A(kc-1)
                  d22 = A(kc+k-1)
                  d12 = A(kc+k-2)
                  d21 = d12
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
!        Form  B := L'*B
!        where L  = P(1)*inv(L(1))* ... *P(m)*inv(L(m))
!        and   L' = inv(L(m))*P(m)* ... *inv(L(1))*P(1)
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
                  IF ( kp/=k ) CALL SSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
!
!                 Apply the transformation
!
                  CALL SGEMV('Transpose',N-k,Nrhs,ONE,B(k+1,1),Ldb,     &
     &                       A(kc+1),1,ONE,B(k,1),Ldb)
               ENDIF
               IF ( nounit ) CALL SSCAL(Nrhs,A(kc),B(k,1),Ldb)
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
                  IF ( kp/=k+1 ) CALL SSWAP(Nrhs,B(k+1,1),Ldb,B(kp,1),  &
     &                 Ldb)
!
!                 Apply the transformation
!
                  CALL SGEMV('Transpose',N-k-1,Nrhs,ONE,B(k+2,1),Ldb,   &
     &                       A(kcnext+1),1,ONE,B(k+1,1),Ldb)
                  CALL SGEMV('Transpose',N-k-1,Nrhs,ONE,B(k+2,1),Ldb,   &
     &                       A(kc+2),1,ONE,B(k,1),Ldb)
               ENDIF
!
!              Multiply by the diagonal block if non-unit.
!
               IF ( nounit ) THEN
                  d11 = A(kc)
                  d22 = A(kcnext)
                  d21 = A(kc+1)
                  d12 = d21
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
!     End of SLAVSP
!
      END SUBROUTINE SLAVSP
