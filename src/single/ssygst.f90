!*==ssygst.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SSYGST
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SSYGST + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssygst.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssygst.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssygst.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SSYGST( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, ITYPE, LDA, LDB, N
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), B( LDB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SSYGST reduces a real symmetric-definite generalized eigenproblem
!> to standard form.
!>
!> If ITYPE = 1, the problem is A*x = lambda*B*x,
!> and A is overwritten by inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T)
!>
!> If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or
!> B*A*x = lambda*x, and A is overwritten by U*A*U**T or L**T*A*L.
!>
!> B must have been previously factorized as U**T*U or L*L**T by SPOTRF.
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
!>          A is REAL array, dimension (LDA,N)
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
!>          B is REAL array, dimension (LDB,N)
!>          The triangular factor from the Cholesky factorization of B,
!>          as returned by SPOTRF.
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
!> \ingroup realSYcomputational
!
!  =====================================================================
      SUBROUTINE SSYGST(Itype,Uplo,N,A,Lda,B,Ldb,Info)
      IMPLICIT NONE
!*--SSYGST131
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER Info , Itype , Lda , Ldb , N
!     ..
!     .. Array Arguments ..
      REAL A(Lda,*) , B(Ldb,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ONE , HALF
      PARAMETER (ONE=1.0,HALF=0.5)
!     ..
!     .. Local Scalars ..
      LOGICAL upper
      INTEGER k , kb , nb
!     ..
!     .. External Subroutines ..
      EXTERNAL SSYGS2 , SSYMM , SSYR2K , STRMM , STRSM , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER ILAENV
      EXTERNAL LSAME , ILAENV
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
         CALL XERBLA('SSYGST',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
!     Determine the block size for this environment.
!
      nb = ILAENV(1,'SSYGST',Uplo,N,-1,-1,-1)
!
      IF ( nb<=1 .OR. nb>=N ) THEN
!
!        Use unblocked code
!
         CALL SSYGS2(Itype,Uplo,N,A,Lda,B,Ldb,Info)
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
               CALL SSYGS2(Itype,Uplo,kb,A(k,k),Lda,B(k,k),Ldb,Info)
               IF ( k+kb<=N ) THEN
                  CALL STRSM('Left',Uplo,'Transpose','Non-unit',kb,     &
     &                       N-k-kb+1,ONE,B(k,k),Ldb,A(k,k+kb),Lda)
                  CALL SSYMM('Left',Uplo,kb,N-k-kb+1,-HALF,A(k,k),Lda,  &
     &                       B(k,k+kb),Ldb,ONE,A(k,k+kb),Lda)
                  CALL SSYR2K(Uplo,'Transpose',N-k-kb+1,kb,-ONE,        &
     &                        A(k,k+kb),Lda,B(k,k+kb),Ldb,ONE,          &
     &                        A(k+kb,k+kb),Lda)
                  CALL SSYMM('Left',Uplo,kb,N-k-kb+1,-HALF,A(k,k),Lda,  &
     &                       B(k,k+kb),Ldb,ONE,A(k,k+kb),Lda)
                  CALL STRSM('Right',Uplo,'No transpose','Non-unit',kb, &
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
               CALL SSYGS2(Itype,Uplo,kb,A(k,k),Lda,B(k,k),Ldb,Info)
               IF ( k+kb<=N ) THEN
                  CALL STRSM('Right',Uplo,'Transpose','Non-unit',       &
     &                       N-k-kb+1,kb,ONE,B(k,k),Ldb,A(k+kb,k),Lda)
                  CALL SSYMM('Right',Uplo,N-k-kb+1,kb,-HALF,A(k,k),Lda, &
     &                       B(k+kb,k),Ldb,ONE,A(k+kb,k),Lda)
                  CALL SSYR2K(Uplo,'No transpose',N-k-kb+1,kb,-ONE,     &
     &                        A(k+kb,k),Lda,B(k+kb,k),Ldb,ONE,          &
     &                        A(k+kb,k+kb),Lda)
                  CALL SSYMM('Right',Uplo,N-k-kb+1,kb,-HALF,A(k,k),Lda, &
     &                       B(k+kb,k),Ldb,ONE,A(k+kb,k),Lda)
                  CALL STRSM('Left',Uplo,'No transpose','Non-unit',     &
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
            CALL STRMM('Left',Uplo,'No transpose','Non-unit',k-1,kb,ONE,&
     &                 B,Ldb,A(1,k),Lda)
            CALL SSYMM('Right',Uplo,k-1,kb,HALF,A(k,k),Lda,B(1,k),Ldb,  &
     &                 ONE,A(1,k),Lda)
            CALL SSYR2K(Uplo,'No transpose',k-1,kb,ONE,A(1,k),Lda,B(1,k)&
     &                  ,Ldb,ONE,A,Lda)
            CALL SSYMM('Right',Uplo,k-1,kb,HALF,A(k,k),Lda,B(1,k),Ldb,  &
     &                 ONE,A(1,k),Lda)
            CALL STRMM('Right',Uplo,'Transpose','Non-unit',k-1,kb,ONE,  &
     &                 B(k,k),Ldb,A(1,k),Lda)
            CALL SSYGS2(Itype,Uplo,kb,A(k,k),Lda,B(k,k),Ldb,Info)
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
            CALL STRMM('Right',Uplo,'No transpose','Non-unit',kb,k-1,   &
     &                 ONE,B,Ldb,A(k,1),Lda)
            CALL SSYMM('Left',Uplo,kb,k-1,HALF,A(k,k),Lda,B(k,1),Ldb,   &
     &                 ONE,A(k,1),Lda)
            CALL SSYR2K(Uplo,'Transpose',k-1,kb,ONE,A(k,1),Lda,B(k,1),  &
     &                  Ldb,ONE,A,Lda)
            CALL SSYMM('Left',Uplo,kb,k-1,HALF,A(k,k),Lda,B(k,1),Ldb,   &
     &                 ONE,A(k,1),Lda)
            CALL STRMM('Left',Uplo,'Transpose','Non-unit',kb,k-1,ONE,   &
     &                 B(k,k),Ldb,A(k,1),Lda)
            CALL SSYGS2(Itype,Uplo,kb,A(k,k),Lda,B(k,k),Ldb,Info)
         ENDDO
      ENDIF
!
!     End of SSYGST
!
      END SUBROUTINE SSYGST
