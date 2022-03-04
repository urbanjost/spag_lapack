!*==chetrs2.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CHETRS2
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CHETRS2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chetrs2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chetrs2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chetrs2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CHETRS2( UPLO, N, NRHS, A, LDA, IPIV, B, LDB,
!                           WORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDA, LDB, N, NRHS
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
!> CHETRS2 solves a system of linear equations A*X = B with a complex
!> Hermitian matrix A using the factorization A = U*D*U**H or
!> A = L*D*L**H computed by CHETRF and converted by CSYCONV.
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
!>          = 'U':  Upper triangular, form is A = U*D*U**H;
!>          = 'L':  Lower triangular, form is A = L*D*L**H.
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
!>          obtain the factor U or L as computed by CHETRF.
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
!>          as determined by CHETRF.
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
!>          WORK is COMPLEX array, dimension (N)
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
!> \ingroup complexHEcomputational
!
!  =====================================================================
      SUBROUTINE CHETRS2(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Work,Info)
      USE S_CSSCAL
      USE S_CSWAP
      USE S_CSYCONV
      USE S_CTRSM
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--CHETRS2136
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0)
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      COMPLEX :: ak , akm1 , akm1k , bk , bkm1 , denom
      INTEGER :: i , iinfo , j , k , kp
      REAL :: s
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
         CALL XERBLA('CHETRS2',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 .OR. Nrhs==0 ) RETURN
!
!     Convert A
!
      CALL CSYCONV(Uplo,'C',N,A,Lda,Ipiv,Work,iinfo)
!
      IF ( upper ) THEN
!
!        Solve A*X = B, where A = U*D*U**H.
!
!       P**T * B
         k = N
         DO WHILE ( k>=1 )
            IF ( Ipiv(k)>0 ) THEN
!           1 x 1 diagonal block
!           Interchange rows K and IPIV(K).
               kp = Ipiv(k)
               IF ( kp/=k ) CALL CSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
               k = k - 1
            ELSE
!           2 x 2 diagonal block
!           Interchange rows K-1 and -IPIV(K).
               kp = -Ipiv(k)
               IF ( kp==-Ipiv(k-1) )                                    &
     &              CALL CSWAP(Nrhs,B(k-1,1),Ldb,B(kp,1),Ldb)
               k = k - 2
            ENDIF
         ENDDO
!
!  Compute (U \P**T * B) -> B    [ (U \P**T * B) ]
!
         CALL CTRSM('L','U','N','U',N,Nrhs,ONE,A,Lda,B,Ldb)
!
!  Compute D \ B -> B   [ D \ (U \P**T * B) ]
!
         i = N
         DO WHILE ( i>=1 )
            IF ( Ipiv(i)>0 ) THEN
               s = REAL(ONE)/REAL(A(i,i))
               CALL CSSCAL(Nrhs,s,B(i,1),Ldb)
            ELSEIF ( i>1 ) THEN
               IF ( Ipiv(i-1)==Ipiv(i) ) THEN
                  akm1k = Work(i)
                  akm1 = A(i-1,i-1)/akm1k
                  ak = A(i,i)/CONJG(akm1k)
                  denom = akm1*ak - ONE
                  DO j = 1 , Nrhs
                     bkm1 = B(i-1,j)/akm1k
                     bk = B(i,j)/CONJG(akm1k)
                     B(i-1,j) = (ak*bkm1-bk)/denom
                     B(i,j) = (akm1*bk-bkm1)/denom
                  ENDDO
                  i = i - 1
               ENDIF
            ENDIF
            i = i - 1
         ENDDO
!
!      Compute (U**H \ B) -> B   [ U**H \ (D \ (U \P**T * B) ) ]
!
         CALL CTRSM('L','U','C','U',N,Nrhs,ONE,A,Lda,B,Ldb)
!
!       P * B  [ P * (U**H \ (D \ (U \P**T * B) )) ]
!
         k = 1
         DO WHILE ( k<=N )
            IF ( Ipiv(k)>0 ) THEN
!           1 x 1 diagonal block
!           Interchange rows K and IPIV(K).
               kp = Ipiv(k)
               IF ( kp/=k ) CALL CSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
               k = k + 1
            ELSE
!           2 x 2 diagonal block
!           Interchange rows K-1 and -IPIV(K).
               kp = -Ipiv(k)
               IF ( k<N .AND. kp==-Ipiv(k+1) )                          &
     &              CALL CSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
               k = k + 2
            ENDIF
         ENDDO
!
      ELSE
!
!        Solve A*X = B, where A = L*D*L**H.
!
!       P**T * B
         k = 1
         DO WHILE ( k<=N )
            IF ( Ipiv(k)>0 ) THEN
!           1 x 1 diagonal block
!           Interchange rows K and IPIV(K).
               kp = Ipiv(k)
               IF ( kp/=k ) CALL CSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
               k = k + 1
            ELSE
!           2 x 2 diagonal block
!           Interchange rows K and -IPIV(K+1).
               kp = -Ipiv(k+1)
               IF ( kp==-Ipiv(k) ) CALL CSWAP(Nrhs,B(k+1,1),Ldb,B(kp,1),&
     &              Ldb)
               k = k + 2
            ENDIF
         ENDDO
!
!  Compute (L \P**T * B) -> B    [ (L \P**T * B) ]
!
         CALL CTRSM('L','L','N','U',N,Nrhs,ONE,A,Lda,B,Ldb)
!
!  Compute D \ B -> B   [ D \ (L \P**T * B) ]
!
         i = 1
         DO WHILE ( i<=N )
            IF ( Ipiv(i)>0 ) THEN
               s = REAL(ONE)/REAL(A(i,i))
               CALL CSSCAL(Nrhs,s,B(i,1),Ldb)
            ELSE
               akm1k = Work(i)
               akm1 = A(i,i)/CONJG(akm1k)
               ak = A(i+1,i+1)/akm1k
               denom = akm1*ak - ONE
               DO j = 1 , Nrhs
                  bkm1 = B(i,j)/CONJG(akm1k)
                  bk = B(i+1,j)/akm1k
                  B(i,j) = (ak*bkm1-bk)/denom
                  B(i+1,j) = (akm1*bk-bkm1)/denom
               ENDDO
               i = i + 1
            ENDIF
            i = i + 1
         ENDDO
!
!  Compute (L**H \ B) -> B   [ L**H \ (D \ (L \P**T * B) ) ]
!
         CALL CTRSM('L','L','C','U',N,Nrhs,ONE,A,Lda,B,Ldb)
!
!       P * B  [ P * (L**H \ (D \ (L \P**T * B) )) ]
!
         k = N
         DO WHILE ( k>=1 )
            IF ( Ipiv(k)>0 ) THEN
!           1 x 1 diagonal block
!           Interchange rows K and IPIV(K).
               kp = Ipiv(k)
               IF ( kp/=k ) CALL CSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
               k = k - 1
            ELSE
!           2 x 2 diagonal block
!           Interchange rows K-1 and -IPIV(K).
               kp = -Ipiv(k)
               IF ( k>1 .AND. kp==-Ipiv(k-1) )                          &
     &              CALL CSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
               k = k - 2
            ENDIF
         ENDDO
!
      ENDIF
!
!     Revert A
!
      CALL CSYCONV(Uplo,'R',N,A,Lda,Ipiv,Work,iinfo)
!
!
!     End of CHETRS2
!
      END SUBROUTINE CHETRS2
