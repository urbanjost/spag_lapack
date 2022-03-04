!*==ssytrs_rook.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SSYTRS_ROOK
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SSYTRS_ROOK + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssytrs_rook.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssytrs_rook.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssytrs_rook.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SSYTRS_ROOK( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDA, LDB, N, NRHS
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       REAL               A( LDA, * ), B( LDB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SSYTRS_ROOK solves a system of linear equations A*X = B with
!> a real symmetric matrix A using the factorization A = U*D*U**T or
!> A = L*D*L**T computed by SSYTRF_ROOK.
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
!>          = 'U':  Upper triangular, form is A = U*D*U**T;
!>          = 'L':  Lower triangular, form is A = L*D*L**T.
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
!>          A is REAL array, dimension (LDA,N)
!>          The block diagonal matrix D and the multipliers used to
!>          obtain the factor U or L as computed by SSYTRF_ROOK.
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
!>          Details of the interchanges and the block structure of D
!>          as determined by SSYTRF_ROOK.
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is REAL array, dimension (LDB,NRHS)
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
!> \date April 2012
!
!> \ingroup realSYcomputational
!
!> \par Contributors:
!  ==================
!>
!> \verbatim
!>
!>   April 2012, Igor Kozachenko,
!>                  Computer Science Division,
!>                  University of California, Berkeley
!>
!>  September 2007, Sven Hammarling, Nicholas J. Higham, Craig Lucas,
!>                  School of Mathematics,
!>                  University of Manchester
!>
!> \endverbatim
!
!  =====================================================================
      SUBROUTINE SSYTRS_ROOK(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Info)
      USE S_LSAME
      USE S_SGEMV
      USE S_SGER
      USE S_SSCAL
      USE S_SSWAP
      USE S_XERBLA
      IMPLICIT NONE
!*--SSYTRS_ROOK145
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER :: Nrhs
      REAL , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      REAL , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL :: ak , akm1 , akm1k , bk , bkm1 , denom
      INTEGER :: j , k , kp
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
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('SSYTRS_ROOK',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 .OR. Nrhs==0 ) RETURN
!
      IF ( upper ) THEN
!
!        Solve A*X = B, where A = U*D*U**T.
!
!        First solve U*D*X = B, overwriting B with X.
!
!        K is the main loop index, decreasing from N to 1 in steps of
!        1 or 2, depending on the size of the diagonal blocks.
!
         k = N
!
!        If K < 1, exit from loop.
!
         DO WHILE ( k>=1 )
!
            IF ( Ipiv(k)>0 ) THEN
!
!           1 x 1 diagonal block
!
!           Interchange rows K and IPIV(K).
!
               kp = Ipiv(k)
               IF ( kp/=k ) CALL SSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
!
!           Multiply by inv(U(K)), where U(K) is the transformation
!           stored in column K of A.
!
               CALL SGER(k-1,Nrhs,-ONE,A(1,k),1,B(k,1),Ldb,B(1,1),Ldb)
!
!           Multiply by the inverse of the diagonal block.
!
               CALL SSCAL(Nrhs,ONE/A(k,k),B(k,1),Ldb)
               k = k - 1
            ELSE
!
!           2 x 2 diagonal block
!
!           Interchange rows K and -IPIV(K) THEN K-1 and -IPIV(K-1)
!
               kp = -Ipiv(k)
               IF ( kp/=k ) CALL SSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
!
               kp = -Ipiv(k-1)
               IF ( kp/=k-1 ) CALL SSWAP(Nrhs,B(k-1,1),Ldb,B(kp,1),Ldb)
!
!           Multiply by inv(U(K)), where U(K) is the transformation
!           stored in columns K-1 and K of A.
!
               IF ( k>2 ) THEN
                  CALL SGER(k-2,Nrhs,-ONE,A(1,k),1,B(k,1),Ldb,B(1,1),   &
     &                      Ldb)
                  CALL SGER(k-2,Nrhs,-ONE,A(1,k-1),1,B(k-1,1),Ldb,B(1,1)&
     &                      ,Ldb)
               ENDIF
!
!           Multiply by the inverse of the diagonal block.
!
               akm1k = A(k-1,k)
               akm1 = A(k-1,k-1)/akm1k
               ak = A(k,k)/akm1k
               denom = akm1*ak - ONE
               DO j = 1 , Nrhs
                  bkm1 = B(k-1,j)/akm1k
                  bk = B(k,j)/akm1k
                  B(k-1,j) = (ak*bkm1-bk)/denom
                  B(k,j) = (akm1*bk-bkm1)/denom
               ENDDO
               k = k - 2
!
            ENDIF
         ENDDO
!
!        Next solve U**T *X = B, overwriting B with X.
!
!        K is the main loop index, increasing from 1 to N in steps of
!        1 or 2, depending on the size of the diagonal blocks.
!
         k = 1
!
!        If K > N, exit from loop.
!
         DO WHILE ( k<=N )
!
            IF ( Ipiv(k)>0 ) THEN
!
!           1 x 1 diagonal block
!
!           Multiply by inv(U**T(K)), where U(K) is the transformation
!           stored in column K of A.
!
               IF ( k>1 ) CALL SGEMV('Transpose',k-1,Nrhs,-ONE,B,Ldb,   &
     &                               A(1,k),1,ONE,B(k,1),Ldb)
!
!           Interchange rows K and IPIV(K).
!
               kp = Ipiv(k)
               IF ( kp/=k ) CALL SSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
               k = k + 1
            ELSE
!
!           2 x 2 diagonal block
!
!           Multiply by inv(U**T(K+1)), where U(K+1) is the transformation
!           stored in columns K and K+1 of A.
!
               IF ( k>1 ) THEN
                  CALL SGEMV('Transpose',k-1,Nrhs,-ONE,B,Ldb,A(1,k),1,  &
     &                       ONE,B(k,1),Ldb)
                  CALL SGEMV('Transpose',k-1,Nrhs,-ONE,B,Ldb,A(1,k+1),1,&
     &                       ONE,B(k+1,1),Ldb)
               ENDIF
!
!           Interchange rows K and -IPIV(K) THEN K+1 and -IPIV(K+1).
!
               kp = -Ipiv(k)
               IF ( kp/=k ) CALL SSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
!
               kp = -Ipiv(k+1)
               IF ( kp/=k+1 ) CALL SSWAP(Nrhs,B(k+1,1),Ldb,B(kp,1),Ldb)
!
               k = k + 2
!
            ENDIF
         ENDDO
!
      ELSE
!
!        Solve A*X = B, where A = L*D*L**T.
!
!        First solve L*D*X = B, overwriting B with X.
!
!        K is the main loop index, increasing from 1 to N in steps of
!        1 or 2, depending on the size of the diagonal blocks.
!
         k = 1
!
!        If K > N, exit from loop.
!
         DO WHILE ( k<=N )
!
            IF ( Ipiv(k)>0 ) THEN
!
!           1 x 1 diagonal block
!
!           Interchange rows K and IPIV(K).
!
               kp = Ipiv(k)
               IF ( kp/=k ) CALL SSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
!
!           Multiply by inv(L(K)), where L(K) is the transformation
!           stored in column K of A.
!
               IF ( k<N ) CALL SGER(N-k,Nrhs,-ONE,A(k+1,k),1,B(k,1),Ldb,&
     &                              B(k+1,1),Ldb)
!
!           Multiply by the inverse of the diagonal block.
!
               CALL SSCAL(Nrhs,ONE/A(k,k),B(k,1),Ldb)
               k = k + 1
            ELSE
!
!           2 x 2 diagonal block
!
!           Interchange rows K and -IPIV(K) THEN K+1 and -IPIV(K+1)
!
               kp = -Ipiv(k)
               IF ( kp/=k ) CALL SSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
!
               kp = -Ipiv(k+1)
               IF ( kp/=k+1 ) CALL SSWAP(Nrhs,B(k+1,1),Ldb,B(kp,1),Ldb)
!
!           Multiply by inv(L(K)), where L(K) is the transformation
!           stored in columns K and K+1 of A.
!
               IF ( k<N-1 ) THEN
                  CALL SGER(N-k-1,Nrhs,-ONE,A(k+2,k),1,B(k,1),Ldb,      &
     &                      B(k+2,1),Ldb)
                  CALL SGER(N-k-1,Nrhs,-ONE,A(k+2,k+1),1,B(k+1,1),Ldb,  &
     &                      B(k+2,1),Ldb)
               ENDIF
!
!           Multiply by the inverse of the diagonal block.
!
               akm1k = A(k+1,k)
               akm1 = A(k,k)/akm1k
               ak = A(k+1,k+1)/akm1k
               denom = akm1*ak - ONE
               DO j = 1 , Nrhs
                  bkm1 = B(k,j)/akm1k
                  bk = B(k+1,j)/akm1k
                  B(k,j) = (ak*bkm1-bk)/denom
                  B(k+1,j) = (akm1*bk-bkm1)/denom
               ENDDO
               k = k + 2
!
            ENDIF
         ENDDO
!
!        Next solve L**T *X = B, overwriting B with X.
!
!        K is the main loop index, decreasing from N to 1 in steps of
!        1 or 2, depending on the size of the diagonal blocks.
!
         k = N
!
!        If K < 1, exit from loop.
!
         DO WHILE ( k>=1 )
!
            IF ( Ipiv(k)>0 ) THEN
!
!           1 x 1 diagonal block
!
!           Multiply by inv(L**T(K)), where L(K) is the transformation
!           stored in column K of A.
!
               IF ( k<N ) CALL SGEMV('Transpose',N-k,Nrhs,-ONE,B(k+1,1),&
     &                               Ldb,A(k+1,k),1,ONE,B(k,1),Ldb)
!
!           Interchange rows K and IPIV(K).
!
               kp = Ipiv(k)
               IF ( kp/=k ) CALL SSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
               k = k - 1
            ELSE
!
!           2 x 2 diagonal block
!
!           Multiply by inv(L**T(K-1)), where L(K-1) is the transformation
!           stored in columns K-1 and K of A.
!
               IF ( k<N ) THEN
                  CALL SGEMV('Transpose',N-k,Nrhs,-ONE,B(k+1,1),Ldb,    &
     &                       A(k+1,k),1,ONE,B(k,1),Ldb)
                  CALL SGEMV('Transpose',N-k,Nrhs,-ONE,B(k+1,1),Ldb,    &
     &                       A(k+1,k-1),1,ONE,B(k-1,1),Ldb)
               ENDIF
!
!           Interchange rows K and -IPIV(K) THEN K-1 and -IPIV(K-1)
!
               kp = -Ipiv(k)
               IF ( kp/=k ) CALL SSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
!
               kp = -Ipiv(k-1)
               IF ( kp/=k-1 ) CALL SSWAP(Nrhs,B(k-1,1),Ldb,B(kp,1),Ldb)
!
               k = k - 2
!
            ENDIF
         ENDDO
      ENDIF
!
!
!     End of SSYTRS_ROOK
!
      END SUBROUTINE SSYTRS_ROOK
