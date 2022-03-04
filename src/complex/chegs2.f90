!*==chegs2.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b CHEGS2 reduces a Hermitian definite generalized eigenproblem to standard form, using the factorization results obtained from cpotrf (unblocked algorithm).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CHEGS2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chegs2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chegs2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chegs2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CHEGS2( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, ITYPE, LDA, LDB, N
!       ..
!       .. Array Arguments ..
!       COMPLEX            A( LDA, * ), B( LDB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CHEGS2 reduces a complex Hermitian-definite generalized
!> eigenproblem to standard form.
!>
!> If ITYPE = 1, the problem is A*x = lambda*B*x,
!> and A is overwritten by inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H)
!>
!> If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or
!> B*A*x = lambda*x, and A is overwritten by U*A*U**H or L**H *A*L.
!>
!> B must have been previously factorized as U**H *U or L*L**H by ZPOTRF.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ITYPE
!> \verbatim
!>          ITYPE is INTEGER
!>          = 1: compute inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H);
!>          = 2 or 3: compute U*A*U**H or L**H *A*L.
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the upper or lower triangular part of the
!>          Hermitian matrix A is stored, and how B has been factorized.
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
!>          A is COMPLEX array, dimension (LDA,N)
!>          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading
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
!> \param[in,out] B
!> \verbatim
!>          B is COMPLEX array, dimension (LDB,N)
!>          The triangular factor from the Cholesky factorization of B,
!>          as returned by CPOTRF.
!>          B is modified by the routine but restored on exit.
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
!> \ingroup complexHEcomputational
!
!  =====================================================================
      SUBROUTINE CHEGS2(Itype,Uplo,N,A,Lda,B,Ldb,Info)
      IMPLICIT NONE
!*--CHEGS2132
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
      COMPLEX A(Lda,*) , B(Ldb,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ONE , HALF
      PARAMETER (ONE=1.0E+0,HALF=0.5E+0)
      COMPLEX CONE
      PARAMETER (CONE=(1.0E+0,0.0E+0))
!     ..
!     .. Local Scalars ..
      LOGICAL upper
      INTEGER k
      REAL akk , bkk
      COMPLEX ct
!     ..
!     .. External Subroutines ..
      EXTERNAL CAXPY , CHER2 , CLACGV , CSSCAL , CTRMV , CTRSV , XERBLA
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
         CALL XERBLA('CHEGS2',-Info)
         RETURN
      ENDIF
!
      IF ( Itype==1 ) THEN
         IF ( upper ) THEN
!
!           Compute inv(U**H)*A*inv(U)
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
                  CALL CSSCAL(N-k,ONE/bkk,A(k,k+1),Lda)
                  ct = -HALF*akk
                  CALL CLACGV(N-k,A(k,k+1),Lda)
                  CALL CLACGV(N-k,B(k,k+1),Ldb)
                  CALL CAXPY(N-k,ct,B(k,k+1),Ldb,A(k,k+1),Lda)
                  CALL CHER2(Uplo,N-k,-CONE,A(k,k+1),Lda,B(k,k+1),Ldb,  &
     &                       A(k+1,k+1),Lda)
                  CALL CAXPY(N-k,ct,B(k,k+1),Ldb,A(k,k+1),Lda)
                  CALL CLACGV(N-k,B(k,k+1),Ldb)
                  CALL CTRSV(Uplo,'Conjugate transpose','Non-unit',N-k, &
     &                       B(k+1,k+1),Ldb,A(k,k+1),Lda)
                  CALL CLACGV(N-k,A(k,k+1),Lda)
               ENDIF
            ENDDO
         ELSE
!
!           Compute inv(L)*A*inv(L**H)
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
                  CALL CSSCAL(N-k,ONE/bkk,A(k+1,k),1)
                  ct = -HALF*akk
                  CALL CAXPY(N-k,ct,B(k+1,k),1,A(k+1,k),1)
                  CALL CHER2(Uplo,N-k,-CONE,A(k+1,k),1,B(k+1,k),1,      &
     &                       A(k+1,k+1),Lda)
                  CALL CAXPY(N-k,ct,B(k+1,k),1,A(k+1,k),1)
                  CALL CTRSV(Uplo,'No transpose','Non-unit',N-k,        &
     &                       B(k+1,k+1),Ldb,A(k+1,k),1)
               ENDIF
            ENDDO
         ENDIF
      ELSEIF ( upper ) THEN
!
!           Compute U*A*U**H
!
         DO k = 1 , N
!
!              Update the upper triangle of A(1:k,1:k)
!
            akk = A(k,k)
            bkk = B(k,k)
            CALL CTRMV(Uplo,'No transpose','Non-unit',k-1,B,Ldb,A(1,k), &
     &                 1)
            ct = HALF*akk
            CALL CAXPY(k-1,ct,B(1,k),1,A(1,k),1)
            CALL CHER2(Uplo,k-1,CONE,A(1,k),1,B(1,k),1,A,Lda)
            CALL CAXPY(k-1,ct,B(1,k),1,A(1,k),1)
            CALL CSSCAL(k-1,bkk,A(1,k),1)
            A(k,k) = akk*bkk**2
         ENDDO
      ELSE
!
!           Compute L**H *A*L
!
         DO k = 1 , N
!
!              Update the lower triangle of A(1:k,1:k)
!
            akk = A(k,k)
            bkk = B(k,k)
            CALL CLACGV(k-1,A(k,1),Lda)
            CALL CTRMV(Uplo,'Conjugate transpose','Non-unit',k-1,B,Ldb, &
     &                 A(k,1),Lda)
            ct = HALF*akk
            CALL CLACGV(k-1,B(k,1),Ldb)
            CALL CAXPY(k-1,ct,B(k,1),Ldb,A(k,1),Lda)
            CALL CHER2(Uplo,k-1,CONE,A(k,1),Lda,B(k,1),Ldb,A,Lda)
            CALL CAXPY(k-1,ct,B(k,1),Ldb,A(k,1),Lda)
            CALL CLACGV(k-1,B(k,1),Ldb)
            CALL CSSCAL(k-1,bkk,A(k,1),Lda)
            CALL CLACGV(k-1,A(k,1),Lda)
            A(k,k) = akk*bkk**2
         ENDDO
      ENDIF
!
!     End of CHEGS2
!
      END SUBROUTINE CHEGS2
