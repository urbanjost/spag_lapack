!*==zhetrs_aa.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b ZHETRS_AA
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZHETRS_AA + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhetrs_aa.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhetrs_aa.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhetrs_aa.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZHETRS_AA( UPLO, N, NRHS, A, LDA, IPIV, B, LDB,
!                             WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            N, NRHS, LDA, LDB, LWORK, INFO
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * )
!       ..
!
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZHETRS_AA solves a system of linear equations A*X = B with a complex
!> hermitian matrix A using the factorization A = U**H*T*U or
!> A = L*T*L**H computed by ZHETRF_AA.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the details of the factorization are stored
!>          as an upper or lower triangular matrix.
!>          = 'U':  Upper triangular, form is A = U**H*T*U;
!>          = 'L':  Lower triangular, form is A = L*T*L**H.
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
!>          Details of factors computed by ZHETRF_AA.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          Details of the interchanges as computed by ZHETRF_AA.
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
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK. LWORK >= max(1,3*N-2).
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
!> \date November 2017
!
!> \ingroup complex16HEcomputational
!
!  =====================================================================
      SUBROUTINE ZHETRS_AA(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Work,Lwork,Info)
!
!  -- LAPACK computational routine (version 3.8.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2017
!
      IMPLICIT NONE
!*--ZHETRS_AA141
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER N , Nrhs , Lda , Ldb , Lwork , Info
!     ..
!     .. Array Arguments ..
      INTEGER Ipiv(*)
      COMPLEX*16 A(Lda,*) , B(Ldb,*) , Work(*)
!     ..
!
!  =====================================================================
!
      COMPLEX*16 ONE
      PARAMETER (ONE=1.0D+0)
!     ..
!     .. Local Scalars ..
      LOGICAL lquery , upper
      INTEGER k , kp , lwkopt
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL ZGTSV , ZSWAP , ZTRSM , ZLACGV , ZLACPY , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     ..
!     .. Executable Statements ..
!
      Info = 0
      upper = LSAME(Uplo,'U')
      lquery = (Lwork==-1)
      IF ( .NOT.upper .AND. .NOT.LSAME(Uplo,'L') ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Nrhs<0 ) THEN
         Info = -3
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -5
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -8
      ELSEIF ( Lwork<MAX(1,3*N-2) .AND. .NOT.lquery ) THEN
         Info = -10
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZHETRS_AA',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         lwkopt = (3*N-2)
         Work(1) = lwkopt
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 .OR. Nrhs==0 ) RETURN
!
      IF ( upper ) THEN
!
!        Solve A*X = B, where A = U**H*T*U.
!
!        1) Forward substitution with U**H
!
         IF ( N>1 ) THEN
!
!           Pivot, P**T * B -> B
!
            DO k = 1 , N
               kp = Ipiv(k)
               IF ( kp/=k ) CALL ZSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
            ENDDO
!
!           Compute U**H \ B -> B    [ (U**H \P**T * B) ]
!
            CALL ZTRSM('L','U','C','U',N-1,Nrhs,ONE,A(1,2),Lda,B(2,1),  &
     &                 Ldb)
         ENDIF
!
!        2) Solve with triangular matrix T
!
!        Compute T \ B -> B   [ T \ (U**H \P**T * B) ]
!
         CALL ZLACPY('F',1,N,A(1,1),Lda+1,Work(N),1)
         IF ( N>1 ) THEN
            CALL ZLACPY('F',1,N-1,A(1,2),Lda+1,Work(2*N),1)
            CALL ZLACPY('F',1,N-1,A(1,2),Lda+1,Work(1),1)
            CALL ZLACGV(N-1,Work(1),1)
         ENDIF
         CALL ZGTSV(N,Nrhs,Work(1),Work(N),Work(2*N),B,Ldb,Info)
!
!        3) Backward substitution with U
!
         IF ( N>1 ) THEN
!
!           Compute U \ B -> B   [ U \ (T \ (U**H \P**T * B) ) ]
!
            CALL ZTRSM('L','U','N','U',N-1,Nrhs,ONE,A(1,2),Lda,B(2,1),  &
     &                 Ldb)
!
!           Pivot, P * B  [ P * (U**H \ (T \ (U \P**T * B) )) ]
!
            DO k = N , 1 , -1
               kp = Ipiv(k)
               IF ( kp/=k ) CALL ZSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
            ENDDO
         ENDIF
!
      ELSE
!
!        Solve A*X = B, where A = L*T*L**H.
!
!        1) Forward substitution with L
!
         IF ( N>1 ) THEN
!
!           Pivot, P**T * B -> B
!
            DO k = 1 , N
               kp = Ipiv(k)
               IF ( kp/=k ) CALL ZSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
            ENDDO
!
!           Compute L \ B -> B    [ (L \P**T * B) ]
!
            CALL ZTRSM('L','L','N','U',N-1,Nrhs,ONE,A(2,1),Lda,B(2,1),  &
     &                 Ldb)
         ENDIF
!
!        2) Solve with triangular matrix T
!
!        Compute T \ B -> B   [ T \ (L \P**T * B) ]
!
         CALL ZLACPY('F',1,N,A(1,1),Lda+1,Work(N),1)
         IF ( N>1 ) THEN
            CALL ZLACPY('F',1,N-1,A(2,1),Lda+1,Work(1),1)
            CALL ZLACPY('F',1,N-1,A(2,1),Lda+1,Work(2*N),1)
            CALL ZLACGV(N-1,Work(2*N),1)
         ENDIF
         CALL ZGTSV(N,Nrhs,Work(1),Work(N),Work(2*N),B,Ldb,Info)
!
!        3) Backward substitution with L**H
!
         IF ( N>1 ) THEN
!
!           Compute L**H \ B -> B   [ L**H \ (T \ (L \P**T * B) ) ]
!
            CALL ZTRSM('L','L','C','U',N-1,Nrhs,ONE,A(2,1),Lda,B(2,1),  &
     &                 Ldb)
!
!           Pivot, P * B  [ P * (L**H \ (T \ (L \P**T * B) )) ]
!
            DO k = N , 1 , -1
               kp = Ipiv(k)
               IF ( kp/=k ) CALL ZSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
            ENDDO
         ENDIF
!
      ENDIF
!
!
!     End of ZHETRS_AA
!
      END SUBROUTINE ZHETRS_AA
