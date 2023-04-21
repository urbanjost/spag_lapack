!*==chptrs.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b CHPTRS
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CHPTRS + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chptrs.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chptrs.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chptrs.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CHPTRS( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDB, N, NRHS
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX            AP( * ), B( LDB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CHPTRS solves a system of linear equations A*X = B with a complex
!> Hermitian matrix A stored in packed format using the factorization
!> A = U*D*U**H or A = L*D*L**H computed by CHPTRF.
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
!> \param[in] AP
!> \verbatim
!>          AP is COMPLEX array, dimension (N*(N+1)/2)
!>          The block diagonal matrix D and the multipliers used to
!>          obtain the factor U or L as computed by CHPTRF, stored as a
!>          packed triangular matrix.
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          Details of the interchanges and the block structure of D
!>          as determined by CHPTRF.
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
!>          < 0: if INFO = -i, the i-th argument had an illegal value
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
!> \ingroup complexOTHERcomputational
!
!  =====================================================================
      SUBROUTINE CHPTRS(Uplo,N,Nrhs,Ap,Ipiv,B,Ldb,Info)
      IMPLICIT NONE
!*--CHPTRS119
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER Info , Ldb , N , Nrhs
!     ..
!     .. Array Arguments ..
      INTEGER Ipiv(*)
      COMPLEX Ap(*) , B(Ldb,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX ONE
      PARAMETER (ONE=(1.0E+0,0.0E+0))
!     ..
!     .. Local Scalars ..
      LOGICAL upper
      INTEGER j , k , kc , kp
      REAL s
      COMPLEX ak , akm1 , akm1k , bk , bkm1 , denom
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL CGEMV , CGERU , CLACGV , CSSCAL , CSWAP , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC CONJG , MAX , REAL
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
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -7
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('CHPTRS',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 .OR. Nrhs==0 ) RETURN
!
      IF ( upper ) THEN
!
!        Solve A*X = B, where A = U*D*U**H.
!
!        First solve U*D*X = B, overwriting B with X.
!
!        K is the main loop index, decreasing from N to 1 in steps of
!        1 or 2, depending on the size of the diagonal blocks.
!
         k = N
         kc = N*(N+1)/2 + 1
!
!        If K < 1, exit from loop.
!
         DO WHILE ( k>=1 )
!
            kc = kc - k
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
               CALL CGERU(k-1,Nrhs,-ONE,Ap(kc),1,B(k,1),Ldb,B(1,1),Ldb)
!
!           Multiply by the inverse of the diagonal block.
!
               s = REAL(ONE)/REAL(Ap(kc+k-1))
               CALL CSSCAL(Nrhs,s,B(k,1),Ldb)
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
               CALL CGERU(k-2,Nrhs,-ONE,Ap(kc),1,B(k,1),Ldb,B(1,1),Ldb)
               CALL CGERU(k-2,Nrhs,-ONE,Ap(kc-(k-1)),1,B(k-1,1),Ldb,    &
     &                    B(1,1),Ldb)
!
!           Multiply by the inverse of the diagonal block.
!
               akm1k = Ap(kc+k-2)
               akm1 = Ap(kc-1)/akm1k
               ak = Ap(kc+k-1)/CONJG(akm1k)
               denom = akm1*ak - ONE
               DO j = 1 , Nrhs
                  bkm1 = B(k-1,j)/akm1k
                  bk = B(k,j)/CONJG(akm1k)
                  B(k-1,j) = (ak*bkm1-bk)/denom
                  B(k,j) = (akm1*bk-bkm1)/denom
               ENDDO
               kc = kc - k + 1
               k = k - 2
!
            ENDIF
         ENDDO
!
!        Next solve U**H *X = B, overwriting B with X.
!
!        K is the main loop index, increasing from 1 to N in steps of
!        1 or 2, depending on the size of the diagonal blocks.
!
         k = 1
         kc = 1
!
!        If K > N, exit from loop.
!
         DO WHILE ( k<=N )
!
            IF ( Ipiv(k)>0 ) THEN
!
!           1 x 1 diagonal block
!
!           Multiply by inv(U**H(K)), where U(K) is the transformation
!           stored in column K of A.
!
               IF ( k>1 ) THEN
                  CALL CLACGV(Nrhs,B(k,1),Ldb)
                  CALL CGEMV('Conjugate transpose',k-1,Nrhs,-ONE,B,Ldb, &
     &                       Ap(kc),1,ONE,B(k,1),Ldb)
                  CALL CLACGV(Nrhs,B(k,1),Ldb)
               ENDIF
!
!           Interchange rows K and IPIV(K).
!
               kp = Ipiv(k)
               IF ( kp/=k ) CALL CSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
               kc = kc + k
               k = k + 1
            ELSE
!
!           2 x 2 diagonal block
!
!           Multiply by inv(U**H(K+1)), where U(K+1) is the transformation
!           stored in columns K and K+1 of A.
!
               IF ( k>1 ) THEN
                  CALL CLACGV(Nrhs,B(k,1),Ldb)
                  CALL CGEMV('Conjugate transpose',k-1,Nrhs,-ONE,B,Ldb, &
     &                       Ap(kc),1,ONE,B(k,1),Ldb)
                  CALL CLACGV(Nrhs,B(k,1),Ldb)
!
                  CALL CLACGV(Nrhs,B(k+1,1),Ldb)
                  CALL CGEMV('Conjugate transpose',k-1,Nrhs,-ONE,B,Ldb, &
     &                       Ap(kc+k),1,ONE,B(k+1,1),Ldb)
                  CALL CLACGV(Nrhs,B(k+1,1),Ldb)
               ENDIF
!
!           Interchange rows K and -IPIV(K).
!
               kp = -Ipiv(k)
               IF ( kp/=k ) CALL CSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
               kc = kc + 2*k + 1
               k = k + 2
!
            ENDIF
         ENDDO
!
      ELSE
!
!        Solve A*X = B, where A = L*D*L**H.
!
!        First solve L*D*X = B, overwriting B with X.
!
!        K is the main loop index, increasing from 1 to N in steps of
!        1 or 2, depending on the size of the diagonal blocks.
!
         k = 1
         kc = 1
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
               IF ( k<N ) CALL CGERU(N-k,Nrhs,-ONE,Ap(kc+1),1,B(k,1),   &
     &                               Ldb,B(k+1,1),Ldb)
!
!           Multiply by the inverse of the diagonal block.
!
               s = REAL(ONE)/REAL(Ap(kc))
               CALL CSSCAL(Nrhs,s,B(k,1),Ldb)
               kc = kc + N - k + 1
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
                  CALL CGERU(N-k-1,Nrhs,-ONE,Ap(kc+2),1,B(k,1),Ldb,     &
     &                       B(k+2,1),Ldb)
                  CALL CGERU(N-k-1,Nrhs,-ONE,Ap(kc+N-k+2),1,B(k+1,1),   &
     &                       Ldb,B(k+2,1),Ldb)
               ENDIF
!
!           Multiply by the inverse of the diagonal block.
!
               akm1k = Ap(kc+1)
               akm1 = Ap(kc)/CONJG(akm1k)
               ak = Ap(kc+N-k+1)/akm1k
               denom = akm1*ak - ONE
               DO j = 1 , Nrhs
                  bkm1 = B(k,j)/CONJG(akm1k)
                  bk = B(k+1,j)/akm1k
                  B(k,j) = (ak*bkm1-bk)/denom
                  B(k+1,j) = (akm1*bk-bkm1)/denom
               ENDDO
               kc = kc + 2*(N-k) + 1
               k = k + 2
!
            ENDIF
         ENDDO
!
!        Next solve L**H *X = B, overwriting B with X.
!
!        K is the main loop index, decreasing from N to 1 in steps of
!        1 or 2, depending on the size of the diagonal blocks.
!
         k = N
         kc = N*(N+1)/2 + 1
!
!        If K < 1, exit from loop.
!
         DO WHILE ( k>=1 )
!
            kc = kc - (N-k+1)
            IF ( Ipiv(k)>0 ) THEN
!
!           1 x 1 diagonal block
!
!           Multiply by inv(L**H(K)), where L(K) is the transformation
!           stored in column K of A.
!
               IF ( k<N ) THEN
                  CALL CLACGV(Nrhs,B(k,1),Ldb)
                  CALL CGEMV('Conjugate transpose',N-k,Nrhs,-ONE,       &
     &                       B(k+1,1),Ldb,Ap(kc+1),1,ONE,B(k,1),Ldb)
                  CALL CLACGV(Nrhs,B(k,1),Ldb)
               ENDIF
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
!           Multiply by inv(L**H(K-1)), where L(K-1) is the transformation
!           stored in columns K-1 and K of A.
!
               IF ( k<N ) THEN
                  CALL CLACGV(Nrhs,B(k,1),Ldb)
                  CALL CGEMV('Conjugate transpose',N-k,Nrhs,-ONE,       &
     &                       B(k+1,1),Ldb,Ap(kc+1),1,ONE,B(k,1),Ldb)
                  CALL CLACGV(Nrhs,B(k,1),Ldb)
!
                  CALL CLACGV(Nrhs,B(k-1,1),Ldb)
                  CALL CGEMV('Conjugate transpose',N-k,Nrhs,-ONE,       &
     &                       B(k+1,1),Ldb,Ap(kc-(N-k)),1,ONE,B(k-1,1),  &
     &                       Ldb)
                  CALL CLACGV(Nrhs,B(k-1,1),Ldb)
               ENDIF
!
!           Interchange rows K and -IPIV(K).
!
               kp = -Ipiv(k)
               IF ( kp/=k ) CALL CSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
               kc = kc - (N-k+2)
               k = k - 2
!
            ENDIF
         ENDDO
      ENDIF
!
!
!     End of CHPTRS
!
      END SUBROUTINE CHPTRS
