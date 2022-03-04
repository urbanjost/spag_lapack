!*==zhegs2.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZHEGS2 reduces a Hermitian definite generalized eigenproblem to standard form, using the factorization results obtained from cpotrf (unblocked algorithm).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZHEGS2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhegs2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhegs2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhegs2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZHEGS2( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, ITYPE, LDA, LDB, N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * ), B( LDB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZHEGS2 reduces a complex Hermitian-definite generalized
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
!>          A is COMPLEX*16 array, dimension (LDA,N)
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
!>          B is COMPLEX*16 array, dimension (LDB,N)
!>          The triangular factor from the Cholesky factorization of B,
!>          as returned by ZPOTRF.
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
!> \ingroup complex16HEcomputational
!
!  =====================================================================
      SUBROUTINE ZHEGS2(Itype,Uplo,N,A,Lda,B,Ldb,Info)
      USE F77KINDS                        
      USE S_LSAME
      USE S_XERBLA
      USE S_ZAXPY
      USE S_ZDSCAL
      USE S_ZHER2
      USE S_ZLACGV
      USE S_ZTRMV
      USE S_ZTRSV
      IMPLICIT NONE
!*--ZHEGS2141
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , HALF = 0.5D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0)
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: Itype
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: akk , bkk
      COMPLEX(CX16KIND) :: ct
      INTEGER :: k
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
         CALL XERBLA('ZHEGS2',-Info)
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
                  CALL ZDSCAL(N-k,ONE/bkk,A(k,k+1),Lda)
                  ct = -HALF*akk
                  CALL ZLACGV(N-k,A(k,k+1),Lda)
                  CALL ZLACGV(N-k,B(k,k+1),Ldb)
                  CALL ZAXPY(N-k,ct,B(k,k+1),Ldb,A(k,k+1),Lda)
                  CALL ZHER2(Uplo,N-k,-CONE,A(k,k+1),Lda,B(k,k+1),Ldb,  &
     &                       A(k+1,k+1),Lda)
                  CALL ZAXPY(N-k,ct,B(k,k+1),Ldb,A(k,k+1),Lda)
                  CALL ZLACGV(N-k,B(k,k+1),Ldb)
                  CALL ZTRSV(Uplo,'Conjugate transpose','Non-unit',N-k, &
     &                       B(k+1,k+1),Ldb,A(k,k+1),Lda)
                  CALL ZLACGV(N-k,A(k,k+1),Lda)
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
                  CALL ZDSCAL(N-k,ONE/bkk,A(k+1,k),1)
                  ct = -HALF*akk
                  CALL ZAXPY(N-k,ct,B(k+1,k),1,A(k+1,k),1)
                  CALL ZHER2(Uplo,N-k,-CONE,A(k+1,k),1,B(k+1,k),1,      &
     &                       A(k+1,k+1),Lda)
                  CALL ZAXPY(N-k,ct,B(k+1,k),1,A(k+1,k),1)
                  CALL ZTRSV(Uplo,'No transpose','Non-unit',N-k,        &
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
            CALL ZTRMV(Uplo,'No transpose','Non-unit',k-1,B,Ldb,A(1,k), &
     &                 1)
            ct = HALF*akk
            CALL ZAXPY(k-1,ct,B(1,k),1,A(1,k),1)
            CALL ZHER2(Uplo,k-1,CONE,A(1,k),1,B(1,k),1,A,Lda)
            CALL ZAXPY(k-1,ct,B(1,k),1,A(1,k),1)
            CALL ZDSCAL(k-1,bkk,A(1,k),1)
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
            CALL ZLACGV(k-1,A(k,1),Lda)
            CALL ZTRMV(Uplo,'Conjugate transpose','Non-unit',k-1,B,Ldb, &
     &                 A(k,1),Lda)
            ct = HALF*akk
            CALL ZLACGV(k-1,B(k,1),Ldb)
            CALL ZAXPY(k-1,ct,B(k,1),Ldb,A(k,1),Lda)
            CALL ZHER2(Uplo,k-1,CONE,A(k,1),Lda,B(k,1),Ldb,A,Lda)
            CALL ZAXPY(k-1,ct,B(k,1),Ldb,A(k,1),Lda)
            CALL ZLACGV(k-1,B(k,1),Ldb)
            CALL ZDSCAL(k-1,bkk,A(k,1),Lda)
            CALL ZLACGV(k-1,A(k,1),Lda)
            A(k,k) = akk*bkk**2
         ENDDO
      ENDIF
!
!     End of ZHEGS2
!
      END SUBROUTINE ZHEGS2
