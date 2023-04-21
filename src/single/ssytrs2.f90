!*==ssytrs2.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SSYTRS2
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SSYTRS2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssytrs2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssytrs2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssytrs2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SSYTRS2( UPLO, N, NRHS, A, LDA, IPIV, B, LDB,
!                           WORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDA, LDB, N, NRHS
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       REAL               A( LDA, * ), B( LDB, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SSYTRS2 solves a system of linear equations A*X = B with a real
!> symmetric matrix A using the factorization A = U*D*U**T or
!> A = L*D*L**T computed by SSYTRF and converted by SSYCONV.
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
!> \param[in,out] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>          The block diagonal matrix D and the multipliers used to
!>          obtain the factor U or L as computed by SSYTRF.
!>          Note that A is input / output. This might be counter-intuitive,
!>          and one may think that A is input only. A is input / output. This
!>          is because, at the start of the subroutine, we permute A in a
!>          "better" form and then we permute A back to its original form at
!>          the end.
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
!>          as determined by SSYTRF.
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
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (N)
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
      SUBROUTINE SSYTRS2(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Work,Info)
      IMPLICIT NONE
!*--SSYTRS2135
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER Info , Lda , Ldb , N , Nrhs
!     ..
!     .. Array Arguments ..
      INTEGER Ipiv(*)
      REAL A(Lda,*) , B(Ldb,*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ONE
      PARAMETER (ONE=1.0E+0)
!     ..
!     .. Local Scalars ..
      LOGICAL upper
      INTEGER i , iinfo , j , k , kp
      REAL ak , akm1 , akm1k , bk , bkm1 , denom
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL SSCAL , SSYCONV , SSWAP , STRSM , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX
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
         CALL XERBLA('SSYTRS2',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 .OR. Nrhs==0 ) RETURN
!
!     Convert A
!
      CALL SSYCONV(Uplo,'C',N,A,Lda,Ipiv,Work,iinfo)
!
      IF ( upper ) THEN
!
!        Solve A*X = B, where A = U*D*U**T.
!
!       P**T * B
         k = N
         DO WHILE ( k>=1 )
            IF ( Ipiv(k)>0 ) THEN
!           1 x 1 diagonal block
!           Interchange rows K and IPIV(K).
               kp = Ipiv(k)
               IF ( kp/=k ) CALL SSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
               k = k - 1
            ELSE
!           2 x 2 diagonal block
!           Interchange rows K-1 and -IPIV(K).
               kp = -Ipiv(k)
               IF ( kp==-Ipiv(k-1) )                                    &
     &              CALL SSWAP(Nrhs,B(k-1,1),Ldb,B(kp,1),Ldb)
               k = k - 2
            ENDIF
         ENDDO
!
!  Compute (U \P**T * B) -> B    [ (U \P**T * B) ]
!
         CALL STRSM('L','U','N','U',N,Nrhs,ONE,A,Lda,B,Ldb)
!
!  Compute D \ B -> B   [ D \ (U \P**T * B) ]
!
         i = N
         DO WHILE ( i>=1 )
            IF ( Ipiv(i)>0 ) THEN
               CALL SSCAL(Nrhs,ONE/A(i,i),B(i,1),Ldb)
            ELSEIF ( i>1 ) THEN
               IF ( Ipiv(i-1)==Ipiv(i) ) THEN
                  akm1k = Work(i)
                  akm1 = A(i-1,i-1)/akm1k
                  ak = A(i,i)/akm1k
                  denom = akm1*ak - ONE
                  DO j = 1 , Nrhs
                     bkm1 = B(i-1,j)/akm1k
                     bk = B(i,j)/akm1k
                     B(i-1,j) = (ak*bkm1-bk)/denom
                     B(i,j) = (akm1*bk-bkm1)/denom
                  ENDDO
                  i = i - 1
               ENDIF
            ENDIF
            i = i - 1
         ENDDO
!
!      Compute (U**T \ B) -> B   [ U**T \ (D \ (U \P**T * B) ) ]
!
         CALL STRSM('L','U','T','U',N,Nrhs,ONE,A,Lda,B,Ldb)
!
!       P * B  [ P * (U**T \ (D \ (U \P**T * B) )) ]
!
         k = 1
         DO WHILE ( k<=N )
            IF ( Ipiv(k)>0 ) THEN
!           1 x 1 diagonal block
!           Interchange rows K and IPIV(K).
               kp = Ipiv(k)
               IF ( kp/=k ) CALL SSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
               k = k + 1
            ELSE
!           2 x 2 diagonal block
!           Interchange rows K-1 and -IPIV(K).
               kp = -Ipiv(k)
               IF ( k<N .AND. kp==-Ipiv(k+1) )                          &
     &              CALL SSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
               k = k + 2
            ENDIF
         ENDDO
!
      ELSE
!
!        Solve A*X = B, where A = L*D*L**T.
!
!       P**T * B
         k = 1
         DO WHILE ( k<=N )
            IF ( Ipiv(k)>0 ) THEN
!           1 x 1 diagonal block
!           Interchange rows K and IPIV(K).
               kp = Ipiv(k)
               IF ( kp/=k ) CALL SSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
               k = k + 1
            ELSE
!           2 x 2 diagonal block
!           Interchange rows K and -IPIV(K+1).
               kp = -Ipiv(k+1)
               IF ( kp==-Ipiv(k) ) CALL SSWAP(Nrhs,B(k+1,1),Ldb,B(kp,1),&
     &              Ldb)
               k = k + 2
            ENDIF
         ENDDO
!
!  Compute (L \P**T * B) -> B    [ (L \P**T * B) ]
!
         CALL STRSM('L','L','N','U',N,Nrhs,ONE,A,Lda,B,Ldb)
!
!  Compute D \ B -> B   [ D \ (L \P**T * B) ]
!
         i = 1
         DO WHILE ( i<=N )
            IF ( Ipiv(i)>0 ) THEN
               CALL SSCAL(Nrhs,ONE/A(i,i),B(i,1),Ldb)
            ELSE
               akm1k = Work(i)
               akm1 = A(i,i)/akm1k
               ak = A(i+1,i+1)/akm1k
               denom = akm1*ak - ONE
               DO j = 1 , Nrhs
                  bkm1 = B(i,j)/akm1k
                  bk = B(i+1,j)/akm1k
                  B(i,j) = (ak*bkm1-bk)/denom
                  B(i+1,j) = (akm1*bk-bkm1)/denom
               ENDDO
               i = i + 1
            ENDIF
            i = i + 1
         ENDDO
!
!  Compute (L**T \ B) -> B   [ L**T \ (D \ (L \P**T * B) ) ]
!
         CALL STRSM('L','L','T','U',N,Nrhs,ONE,A,Lda,B,Ldb)
!
!       P * B  [ P * (L**T \ (D \ (L \P**T * B) )) ]
!
         k = N
         DO WHILE ( k>=1 )
            IF ( Ipiv(k)>0 ) THEN
!           1 x 1 diagonal block
!           Interchange rows K and IPIV(K).
               kp = Ipiv(k)
               IF ( kp/=k ) CALL SSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
               k = k - 1
            ELSE
!           2 x 2 diagonal block
!           Interchange rows K-1 and -IPIV(K).
               kp = -Ipiv(k)
               IF ( k>1 .AND. kp==-Ipiv(k-1) )                          &
     &              CALL SSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
               k = k - 2
            ENDIF
         ENDDO
!
      ENDIF
!
!     Revert A
!
      CALL SSYCONV(Uplo,'R',N,A,Lda,Ipiv,Work,iinfo)
!
!
!     End of SSYTRS2
!
      END SUBROUTINE SSYTRS2
