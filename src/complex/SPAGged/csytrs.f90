!*==csytrs.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CSYTRS
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CSYTRS + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csytrs.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csytrs.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csytrs.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CSYTRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDA, LDB, N, NRHS
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX            A( LDA, * ), B( LDB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CSYTRS solves a system of linear equations A*X = B with a complex
!> symmetric matrix A using the factorization A = U*D*U**T or
!> A = L*D*L**T computed by CSYTRF.
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
!>          A is COMPLEX array, dimension (LDA,N)
!>          The block diagonal matrix D and the multipliers used to
!>          obtain the factor U or L as computed by CSYTRF.
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
!>          as determined by CSYTRF.
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
!> \ingroup complexSYcomputational
!
!  =====================================================================
      SUBROUTINE CSYTRS(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Info)
      USE S_CGEMV
      USE S_CGERU
      USE S_CSCAL
      USE S_CSWAP
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--CSYTRS130
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0)
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      COMPLEX :: ak , akm1 , akm1k , bk , bkm1 , denom
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
         CALL XERBLA('CSYTRS',-Info)
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
               IF ( kp/=k ) CALL CSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
!
!           Multiply by inv(U(K)), where U(K) is the transformation
!           stored in column K of A.
!
               CALL CGERU(k-1,Nrhs,-ONE,A(1,k),1,B(k,1),Ldb,B(1,1),Ldb)
!
!           Multiply by the inverse of the diagonal block.
!
               CALL CSCAL(Nrhs,ONE/A(k,k),B(k,1),Ldb)
               k = k - 1
            ELSE
!
!           2 x 2 diagonal block
!
!           Interchange rows K-1 and -IPIV(K).
!
               kp = -Ipiv(k)
               IF ( kp/=k-1 ) CALL CSWAP(Nrhs,B(k-1,1),Ldb,B(kp,1),Ldb)
!
!           Multiply by inv(U(K)), where U(K) is the transformation
!           stored in columns K-1 and K of A.
!
               CALL CGERU(k-2,Nrhs,-ONE,A(1,k),1,B(k,1),Ldb,B(1,1),Ldb)
               CALL CGERU(k-2,Nrhs,-ONE,A(1,k-1),1,B(k-1,1),Ldb,B(1,1), &
     &                    Ldb)
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
               CALL CGEMV('Transpose',k-1,Nrhs,-ONE,B,Ldb,A(1,k),1,ONE, &
     &                    B(k,1),Ldb)
!
!           Interchange rows K and IPIV(K).
!
               kp = Ipiv(k)
               IF ( kp/=k ) CALL CSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
               k = k + 1
            ELSE
!
!           2 x 2 diagonal block
!
!           Multiply by inv(U**T(K+1)), where U(K+1) is the transformation
!           stored in columns K and K+1 of A.
!
               CALL CGEMV('Transpose',k-1,Nrhs,-ONE,B,Ldb,A(1,k),1,ONE, &
     &                    B(k,1),Ldb)
               CALL CGEMV('Transpose',k-1,Nrhs,-ONE,B,Ldb,A(1,k+1),1,   &
     &                    ONE,B(k+1,1),Ldb)
!
!           Interchange rows K and -IPIV(K).
!
               kp = -Ipiv(k)
               IF ( kp/=k ) CALL CSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
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
               IF ( kp/=k ) CALL CSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
!
!           Multiply by inv(L(K)), where L(K) is the transformation
!           stored in column K of A.
!
               IF ( k<N ) CALL CGERU(N-k,Nrhs,-ONE,A(k+1,k),1,B(k,1),   &
     &                               Ldb,B(k+1,1),Ldb)
!
!           Multiply by the inverse of the diagonal block.
!
               CALL CSCAL(Nrhs,ONE/A(k,k),B(k,1),Ldb)
               k = k + 1
            ELSE
!
!           2 x 2 diagonal block
!
!           Interchange rows K+1 and -IPIV(K).
!
               kp = -Ipiv(k)
               IF ( kp/=k+1 ) CALL CSWAP(Nrhs,B(k+1,1),Ldb,B(kp,1),Ldb)
!
!           Multiply by inv(L(K)), where L(K) is the transformation
!           stored in columns K and K+1 of A.
!
               IF ( k<N-1 ) THEN
                  CALL CGERU(N-k-1,Nrhs,-ONE,A(k+2,k),1,B(k,1),Ldb,     &
     &                       B(k+2,1),Ldb)
                  CALL CGERU(N-k-1,Nrhs,-ONE,A(k+2,k+1),1,B(k+1,1),Ldb, &
     &                       B(k+2,1),Ldb)
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
               IF ( k<N ) CALL CGEMV('Transpose',N-k,Nrhs,-ONE,B(k+1,1),&
     &                               Ldb,A(k+1,k),1,ONE,B(k,1),Ldb)
!
!           Interchange rows K and IPIV(K).
!
               kp = Ipiv(k)
               IF ( kp/=k ) CALL CSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
               k = k - 1
            ELSE
!
!           2 x 2 diagonal block
!
!           Multiply by inv(L**T(K-1)), where L(K-1) is the transformation
!           stored in columns K-1 and K of A.
!
               IF ( k<N ) THEN
                  CALL CGEMV('Transpose',N-k,Nrhs,-ONE,B(k+1,1),Ldb,    &
     &                       A(k+1,k),1,ONE,B(k,1),Ldb)
                  CALL CGEMV('Transpose',N-k,Nrhs,-ONE,B(k+1,1),Ldb,    &
     &                       A(k+1,k-1),1,ONE,B(k-1,1),Ldb)
               ENDIF
!
!           Interchange rows K and -IPIV(K).
!
               kp = -Ipiv(k)
               IF ( kp/=k ) CALL CSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
               k = k - 2
!
            ENDIF
         ENDDO
      ENDIF
!
!
!     End of CSYTRS
!
      END SUBROUTINE CSYTRS
