!*==zsytrs_aa.f90  processed by SPAG 7.51RB at 20:09 on  3 Mar 2022
!> \brief \b ZSYTRS_AA
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZSYTRS_AA + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsytrs_aa.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsytrs_aa.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsytrs_aa.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZSYTRS_AA( UPLO, N, NRHS, A, LDA, IPIV, B, LDB,
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
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZSYTRS_AA solves a system of linear equations A*X = B with a complex
!> symmetric matrix A using the factorization A = U**T*T*U or
!> A = L*T*L**T computed by ZSYTRF_AA.
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
!>          = 'U':  Upper triangular, form is A = U**T*T*U;
!>          = 'L':  Lower triangular, form is A = L*T*L**T.
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
!>          Details of factors computed by ZSYTRF_AA.
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
!>          Details of the interchanges as computed by ZSYTRF_AA.
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
!> \ingroup complex16SYcomputational
!
!  =====================================================================
      SUBROUTINE ZSYTRS_AA(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Work,Lwork,Info)
!
!  -- LAPACK computational routine (version 3.8.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2017
!
      IMPLICIT NONE
!*--ZSYTRS_AA140
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
      EXTERNAL ZGTSV , ZSWAP , ZLACPY , ZTRSM , XERBLA
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
         CALL XERBLA('ZSYTRS_AA',-Info)
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
!        Solve A*X = B, where A = U**T*T*U.
!
!        1) Forward substitution with U**T
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
!           Compute U**T \ B -> B    [ (U**T \P**T * B) ]
!
            CALL ZTRSM('L','U','T','U',N-1,Nrhs,ONE,A(1,2),Lda,B(2,1),  &
     &                 Ldb)
         ENDIF
!
!        2) Solve with triangular matrix T
!
!        Compute T \ B -> B   [ T \ (U**T \P**T * B) ]
!
         CALL ZLACPY('F',1,N,A(1,1),Lda+1,Work(N),1)
         IF ( N>1 ) THEN
            CALL ZLACPY('F',1,N-1,A(1,2),Lda+1,Work(1),1)
            CALL ZLACPY('F',1,N-1,A(1,2),Lda+1,Work(2*N),1)
         ENDIF
         CALL ZGTSV(N,Nrhs,Work(1),Work(N),Work(2*N),B,Ldb,Info)
!
!        3) Backward substitution with U
!
         IF ( N>1 ) THEN
!
!           Compute U \ B -> B   [ U \ (T \ (U**T \P**T * B) ) ]
!
            CALL ZTRSM('L','U','N','U',N-1,Nrhs,ONE,A(1,2),Lda,B(2,1),  &
     &                 Ldb)
!
!           Pivot, P * B -> B  [ P * (U \ (T \ (U**T \P**T * B) )) ]
!
            DO k = N , 1 , -1
               kp = Ipiv(k)
               IF ( kp/=k ) CALL ZSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
            ENDDO
         ENDIF
!
      ELSE
!
!        Solve A*X = B, where A = L*T*L**T.
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
         ENDIF
         CALL ZGTSV(N,Nrhs,Work(1),Work(N),Work(2*N),B,Ldb,Info)
!
!        3) Backward substitution with L**T
!
         IF ( N>1 ) THEN
!
!           Compute (L**T \ B) -> B   [ L**T \ (T \ (L \P**T * B) ) ]
!
            CALL ZTRSM('L','L','T','U',N-1,Nrhs,ONE,A(2,1),Lda,B(2,1),  &
     &                 Ldb)
!
!           Pivot, P * B -> B  [ P * (L**T \ (T \ (L \P**T * B) )) ]
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
!     End of ZSYTRS_AA
!
      END SUBROUTINE ZSYTRS_AA
