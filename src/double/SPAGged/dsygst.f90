!*==dsygst.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DSYGST
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DSYGST + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsygst.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsygst.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsygst.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSYGST( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, ITYPE, LDA, LDB, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DSYGST reduces a real symmetric-definite generalized eigenproblem
!> to standard form.
!>
!> If ITYPE = 1, the problem is A*x = lambda*B*x,
!> and A is overwritten by inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T)
!>
!> If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or
!> B*A*x = lambda*x, and A is overwritten by U*A*U**T or L**T*A*L.
!>
!> B must have been previously factorized as U**T*U or L*L**T by DPOTRF.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ITYPE
!> \verbatim
!>          ITYPE is INTEGER
!>          = 1: compute inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T);
!>          = 2 or 3: compute U*A*U**T or L**T*A*L.
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  Upper triangle of A is stored and B is factored as
!>                  U**T*U;
!>          = 'L':  Lower triangle of A is stored and B is factored as
!>                  L*L**T.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrices A and B.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!>          N-by-N upper triangular part of A contains the upper
!>          triangular part of the matrix A, and the strictly lower
!>          triangular part of A is not referenced.  If UPLO = 'L', the
!>          leading N-by-N lower triangular part of A contains the lower
!>          triangular part of the matrix A, and the strictly upper
!>          triangular part of A is not referenced.
!>
!>          On exit, if INFO = 0, the transformed matrix, stored in the
!>          same format as A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (LDB,N)
!>          The triangular factor from the Cholesky factorization of B,
!>          as returned by DPOTRF.
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
!> \date December 2016
!
!> \ingroup doubleSYcomputational
!
!  =====================================================================
      SUBROUTINE DSYGST(Itype,Uplo,N,A,Lda,B,Ldb,Info)
      USE F77KINDS                        
      USE S_DSYGS2
      USE S_DSYMM
      USE S_DSYR2K
      USE S_DTRMM
      USE S_DTRSM
      USE S_ILAENV
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--DSYGST140
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D0 , HALF = 0.5D0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER :: Itype
      CHARACTER :: Uplo
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: k , kb , nb
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
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. External Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
      upper = LSAME(Uplo,'U')
      IF ( Itype<1 .OR. Itype>3 ) THEN
         Info = -1
      ELSEIF ( .NOT.upper .AND. .NOT.LSAME(Uplo,'L') ) THEN
         Info = -2
      ELSEIF ( N<0 ) THEN
         Info = -3
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -5
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -7
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('DSYGST',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
!     Determine the block size for this environment.
!
      nb = ILAENV(1,'DSYGST',Uplo,N,-1,-1,-1)
!
      IF ( nb<=1 .OR. nb>=N ) THEN
!
!        Use unblocked code
!
         CALL DSYGS2(Itype,Uplo,N,A,Lda,B,Ldb,Info)
!
!        Use blocked code
!
      ELSEIF ( Itype==1 ) THEN
         IF ( upper ) THEN
!
!              Compute inv(U**T)*A*inv(U)
!
            DO k = 1 , N , nb
               kb = MIN(N-k+1,nb)
!
!                 Update the upper triangle of A(k:n,k:n)
!
               CALL DSYGS2(Itype,Uplo,kb,A(k,k),Lda,B(k,k),Ldb,Info)
               IF ( k+kb<=N ) THEN
                  CALL DTRSM('Left',Uplo,'Transpose','Non-unit',kb,     &
     &                       N-k-kb+1,ONE,B(k,k),Ldb,A(k,k+kb),Lda)
                  CALL DSYMM('Left',Uplo,kb,N-k-kb+1,-HALF,A(k,k),Lda,  &
     &                       B(k,k+kb),Ldb,ONE,A(k,k+kb),Lda)
                  CALL DSYR2K(Uplo,'Transpose',N-k-kb+1,kb,-ONE,        &
     &                        A(k,k+kb),Lda,B(k,k+kb),Ldb,ONE,          &
     &                        A(k+kb,k+kb),Lda)
                  CALL DSYMM('Left',Uplo,kb,N-k-kb+1,-HALF,A(k,k),Lda,  &
     &                       B(k,k+kb),Ldb,ONE,A(k,k+kb),Lda)
                  CALL DTRSM('Right',Uplo,'No transpose','Non-unit',kb, &
     &                       N-k-kb+1,ONE,B(k+kb,k+kb),Ldb,A(k,k+kb),   &
     &                       Lda)
               ENDIF
            ENDDO
         ELSE
!
!              Compute inv(L)*A*inv(L**T)
!
            DO k = 1 , N , nb
               kb = MIN(N-k+1,nb)
!
!                 Update the lower triangle of A(k:n,k:n)
!
               CALL DSYGS2(Itype,Uplo,kb,A(k,k),Lda,B(k,k),Ldb,Info)
               IF ( k+kb<=N ) THEN
                  CALL DTRSM('Right',Uplo,'Transpose','Non-unit',       &
     &                       N-k-kb+1,kb,ONE,B(k,k),Ldb,A(k+kb,k),Lda)
                  CALL DSYMM('Right',Uplo,N-k-kb+1,kb,-HALF,A(k,k),Lda, &
     &                       B(k+kb,k),Ldb,ONE,A(k+kb,k),Lda)
                  CALL DSYR2K(Uplo,'No transpose',N-k-kb+1,kb,-ONE,     &
     &                        A(k+kb,k),Lda,B(k+kb,k),Ldb,ONE,          &
     &                        A(k+kb,k+kb),Lda)
                  CALL DSYMM('Right',Uplo,N-k-kb+1,kb,-HALF,A(k,k),Lda, &
     &                       B(k+kb,k),Ldb,ONE,A(k+kb,k),Lda)
                  CALL DTRSM('Left',Uplo,'No transpose','Non-unit',     &
     &                       N-k-kb+1,kb,ONE,B(k+kb,k+kb),Ldb,A(k+kb,k),&
     &                       Lda)
               ENDIF
            ENDDO
         ENDIF
      ELSEIF ( upper ) THEN
!
!              Compute U*A*U**T
!
         DO k = 1 , N , nb
            kb = MIN(N-k+1,nb)
!
!                 Update the upper triangle of A(1:k+kb-1,1:k+kb-1)
!
            CALL DTRMM('Left',Uplo,'No transpose','Non-unit',k-1,kb,ONE,&
     &                 B,Ldb,A(1,k),Lda)
            CALL DSYMM('Right',Uplo,k-1,kb,HALF,A(k,k),Lda,B(1,k),Ldb,  &
     &                 ONE,A(1,k),Lda)
            CALL DSYR2K(Uplo,'No transpose',k-1,kb,ONE,A(1,k),Lda,B(1,k)&
     &                  ,Ldb,ONE,A,Lda)
            CALL DSYMM('Right',Uplo,k-1,kb,HALF,A(k,k),Lda,B(1,k),Ldb,  &
     &                 ONE,A(1,k),Lda)
            CALL DTRMM('Right',Uplo,'Transpose','Non-unit',k-1,kb,ONE,  &
     &                 B(k,k),Ldb,A(1,k),Lda)
            CALL DSYGS2(Itype,Uplo,kb,A(k,k),Lda,B(k,k),Ldb,Info)
         ENDDO
      ELSE
!
!              Compute L**T*A*L
!
         DO k = 1 , N , nb
            kb = MIN(N-k+1,nb)
!
!                 Update the lower triangle of A(1:k+kb-1,1:k+kb-1)
!
            CALL DTRMM('Right',Uplo,'No transpose','Non-unit',kb,k-1,   &
     &                 ONE,B,Ldb,A(k,1),Lda)
            CALL DSYMM('Left',Uplo,kb,k-1,HALF,A(k,k),Lda,B(k,1),Ldb,   &
     &                 ONE,A(k,1),Lda)
            CALL DSYR2K(Uplo,'Transpose',k-1,kb,ONE,A(k,1),Lda,B(k,1),  &
     &                  Ldb,ONE,A,Lda)
            CALL DSYMM('Left',Uplo,kb,k-1,HALF,A(k,k),Lda,B(k,1),Ldb,   &
     &                 ONE,A(k,1),Lda)
            CALL DTRMM('Left',Uplo,'Transpose','Non-unit',kb,k-1,ONE,   &
     &                 B(k,k),Ldb,A(k,1),Lda)
            CALL DSYGS2(Itype,Uplo,kb,A(k,k),Lda,B(k,k),Ldb,Info)
         ENDDO
      ENDIF
!
!     End of DSYGST
!
      END SUBROUTINE DSYGST
