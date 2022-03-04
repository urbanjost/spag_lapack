!*==chetrs_aa.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CHETRS_AA
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CHETRS_AA + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chetrs_aa.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chetrs_aa.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chetrs_aa.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CHETRS_AA( UPLO, N, NRHS, A, LDA, IPIV, B, LDB,
!                             WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            N, NRHS, LDA, LDB, LWORK, INFO
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX            A( LDA, * ), B( LDB, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CHETRS_AA solves a system of linear equations A*X = B with a complex
!> hermitian matrix A using the factorization A = U**H*T*U or
!> A = L*T*L**H computed by CHETRF_AA.
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
!>          A is COMPLEX array, dimension (LDA,N)
!>          Details of factors computed by CHETRF_AA.
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
!>          Details of the interchanges as computed by CHETRF_AA.
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is COMPLEX array, dimension (LDB,NRHS)
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
!>          WORK is COMPLEX array, dimension (MAX(1,LWORK))
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
!> \ingroup complexHEcomputational
!
!  =====================================================================
      SUBROUTINE CHETRS_AA(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Work,Lwork,Info)
!
!  -- LAPACK computational routine (version 3.8.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2017
!
      USE S_CGTSV
      USE S_CLACGV
      USE S_CLACPY
      USE S_CSWAP
      USE S_CTRSM
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--CHETRS_AA147
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX , PARAMETER  ::  ONE = 1.0E+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: k , kp , lwkopt
      LOGICAL :: lquery , upper
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
!     ..
!     .. Local Scalars ..
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
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
         CALL XERBLA('CHETRS_AA',-Info)
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
            k = 1
            DO WHILE ( k<=N )
               kp = Ipiv(k)
               IF ( kp/=k ) CALL CSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
               k = k + 1
            ENDDO
!
!           Compute U**H \ B -> B    [ (U**H \P**T * B) ]
!
            CALL CTRSM('L','U','C','U',N-1,Nrhs,ONE,A(1,2),Lda,B(2,1),  &
     &                 Ldb)
         ENDIF
!
!        2) Solve with triangular matrix T
!
!        Compute T \ B -> B   [ T \ (U**H \P**T * B) ]
!
         CALL CLACPY('F',1,N,A(1,1),Lda+1,Work(N),1)
         IF ( N>1 ) THEN
            CALL CLACPY('F',1,N-1,A(1,2),Lda+1,Work(2*N),1)
            CALL CLACPY('F',1,N-1,A(1,2),Lda+1,Work(1),1)
            CALL CLACGV(N-1,Work(1),1)
         ENDIF
         CALL CGTSV(N,Nrhs,Work(1),Work(N),Work(2*N),B,Ldb,Info)
!
!        3) Backward substitution with U
!
         IF ( N>1 ) THEN
!
!           Compute U \ B -> B   [ U \ (T \ (U**H \P**T * B) ) ]
!
            CALL CTRSM('L','U','N','U',N-1,Nrhs,ONE,A(1,2),Lda,B(2,1),  &
     &                 Ldb)
!
!           Pivot, P * B  -> B [ P * (U \ (T \ (U**H \P**T * B) )) ]
!
            k = N
            DO WHILE ( k>=1 )
               kp = Ipiv(k)
               IF ( kp/=k ) CALL CSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
               k = k - 1
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
            k = 1
            DO WHILE ( k<=N )
               kp = Ipiv(k)
               IF ( kp/=k ) CALL CSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
               k = k + 1
            ENDDO
!
!           Compute L \ B -> B    [ (L \P**T * B) ]
!
            CALL CTRSM('L','L','N','U',N-1,Nrhs,ONE,A(2,1),Lda,B(2,1),  &
     &                 Ldb)
         ENDIF
!
!        2) Solve with triangular matrix T
!
!        Compute T \ B -> B   [ T \ (L \P**T * B) ]
!
         CALL CLACPY('F',1,N,A(1,1),Lda+1,Work(N),1)
         IF ( N>1 ) THEN
            CALL CLACPY('F',1,N-1,A(2,1),Lda+1,Work(1),1)
            CALL CLACPY('F',1,N-1,A(2,1),Lda+1,Work(2*N),1)
            CALL CLACGV(N-1,Work(2*N),1)
         ENDIF
         CALL CGTSV(N,Nrhs,Work(1),Work(N),Work(2*N),B,Ldb,Info)
!
!        3) Backward substitution with L**H
!
         IF ( N>1 ) THEN
!
!           Compute (L**H \ B) -> B   [ L**H \ (T \ (L \P**T * B) ) ]
!
            CALL CTRSM('L','L','C','U',N-1,Nrhs,ONE,A(2,1),Lda,B(2,1),  &
     &                 Ldb)
!
!           Pivot, P * B -> B  [ P * (L**H \ (T \ (L \P**T * B) )) ]
!
            k = N
            DO WHILE ( k>=1 )
               kp = Ipiv(k)
               IF ( kp/=k ) CALL CSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
               k = k - 1
            ENDDO
         ENDIF
!
      ENDIF
!
!
!     End of CHETRS_AA
!
      END SUBROUTINE CHETRS_AA
