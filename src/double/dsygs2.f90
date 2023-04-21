!*==dsygs2.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b DSYGS2 reduces a symmetric definite generalized eigenproblem to standard form, using the factorization results obtained from spotrf (unblocked algorithm).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DSYGS2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsygs2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsygs2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsygs2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSYGS2( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
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
!> DSYGS2 reduces a real symmetric-definite generalized eigenproblem
!> to standard form.
!>
!> If ITYPE = 1, the problem is A*x = lambda*B*x,
!> and A is overwritten by inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T)
!>
!> If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or
!> B*A*x = lambda*x, and A is overwritten by U*A*U**T or L**T *A*L.
!>
!> B must have been previously factorized as U**T *U or L*L**T by DPOTRF.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ITYPE
!> \verbatim
!>          ITYPE is INTEGER
!>          = 1: compute inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T);
!>          = 2 or 3: compute U*A*U**T or L**T *A*L.
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the upper or lower triangular part of the
!>          symmetric matrix A is stored, and how B has been factorized.
!>          = 'U':  Upper triangular
!>          = 'L':  Lower triangular
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
!>          n by n upper triangular part of A contains the upper
!>          triangular part of the matrix A, and the strictly lower
!>          triangular part of A is not referenced.  If UPLO = 'L', the
!>          leading n by n lower triangular part of A contains the lower
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
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
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
      SUBROUTINE DSYGS2(Itype,Uplo,N,A,Lda,B,Ldb,Info)
      IMPLICIT NONE
!*--DSYGS2131
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
      DOUBLE PRECISION A(Lda,*) , B(Ldb,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ONE , HALF
      PARAMETER (ONE=1.0D0,HALF=0.5D0)
!     ..
!     .. Local Scalars ..
      LOGICAL upper
      INTEGER k
      DOUBLE PRECISION akk , bkk , ct
!     ..
!     .. External Subroutines ..
      EXTERNAL DAXPY , DSCAL , DSYR2 , DTRMV , DTRSV , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
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
         CALL XERBLA('DSYGS2',-Info)
         RETURN
      ENDIF
!
      IF ( Itype==1 ) THEN
         IF ( upper ) THEN
!
!           Compute inv(U**T)*A*inv(U)
!
            DO k = 1 , N
!
!              Update the upper triangle of A(k:n,k:n)
!
               akk = A(k,k)
               bkk = B(k,k)
               akk = akk/bkk**2
               A(k,k) = akk
               IF ( k<N ) THEN
                  CALL DSCAL(N-k,ONE/bkk,A(k,k+1),Lda)
                  ct = -HALF*akk
                  CALL DAXPY(N-k,ct,B(k,k+1),Ldb,A(k,k+1),Lda)
                  CALL DSYR2(Uplo,N-k,-ONE,A(k,k+1),Lda,B(k,k+1),Ldb,   &
     &                       A(k+1,k+1),Lda)
                  CALL DAXPY(N-k,ct,B(k,k+1),Ldb,A(k,k+1),Lda)
                  CALL DTRSV(Uplo,'Transpose','Non-unit',N-k,B(k+1,k+1),&
     &                       Ldb,A(k,k+1),Lda)
               ENDIF
            ENDDO
         ELSE
!
!           Compute inv(L)*A*inv(L**T)
!
            DO k = 1 , N
!
!              Update the lower triangle of A(k:n,k:n)
!
               akk = A(k,k)
               bkk = B(k,k)
               akk = akk/bkk**2
               A(k,k) = akk
               IF ( k<N ) THEN
                  CALL DSCAL(N-k,ONE/bkk,A(k+1,k),1)
                  ct = -HALF*akk
                  CALL DAXPY(N-k,ct,B(k+1,k),1,A(k+1,k),1)
                  CALL DSYR2(Uplo,N-k,-ONE,A(k+1,k),1,B(k+1,k),1,       &
     &                       A(k+1,k+1),Lda)
                  CALL DAXPY(N-k,ct,B(k+1,k),1,A(k+1,k),1)
                  CALL DTRSV(Uplo,'No transpose','Non-unit',N-k,        &
     &                       B(k+1,k+1),Ldb,A(k+1,k),1)
               ENDIF
            ENDDO
         ENDIF
      ELSEIF ( upper ) THEN
!
!           Compute U*A*U**T
!
         DO k = 1 , N
!
!              Update the upper triangle of A(1:k,1:k)
!
            akk = A(k,k)
            bkk = B(k,k)
            CALL DTRMV(Uplo,'No transpose','Non-unit',k-1,B,Ldb,A(1,k), &
     &                 1)
            ct = HALF*akk
            CALL DAXPY(k-1,ct,B(1,k),1,A(1,k),1)
            CALL DSYR2(Uplo,k-1,ONE,A(1,k),1,B(1,k),1,A,Lda)
            CALL DAXPY(k-1,ct,B(1,k),1,A(1,k),1)
            CALL DSCAL(k-1,bkk,A(1,k),1)
            A(k,k) = akk*bkk**2
         ENDDO
      ELSE
!
!           Compute L**T *A*L
!
         DO k = 1 , N
!
!              Update the lower triangle of A(1:k,1:k)
!
            akk = A(k,k)
            bkk = B(k,k)
            CALL DTRMV(Uplo,'Transpose','Non-unit',k-1,B,Ldb,A(k,1),Lda)
            ct = HALF*akk
            CALL DAXPY(k-1,ct,B(k,1),Ldb,A(k,1),Lda)
            CALL DSYR2(Uplo,k-1,ONE,A(k,1),Lda,B(k,1),Ldb,A,Lda)
            CALL DAXPY(k-1,ct,B(k,1),Ldb,A(k,1),Lda)
            CALL DSCAL(k-1,bkk,A(k,1),Lda)
            A(k,k) = akk*bkk**2
         ENDDO
      ENDIF
!
!     End of DSYGS2
!
      END SUBROUTINE DSYGS2
