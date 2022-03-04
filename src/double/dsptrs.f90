!*==dsptrs.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b DSPTRS
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DSPTRS + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsptrs.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsptrs.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsptrs.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSPTRS( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDB, N, NRHS
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       DOUBLE PRECISION   AP( * ), B( LDB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DSPTRS solves a system of linear equations A*X = B with a real
!> symmetric matrix A stored in packed format using the factorization
!> A = U*D*U**T or A = L*D*L**T computed by DSPTRF.
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
!> \param[in] AP
!> \verbatim
!>          AP is DOUBLE PRECISION array, dimension (N*(N+1)/2)
!>          The block diagonal matrix D and the multipliers used to
!>          obtain the factor U or L as computed by DSPTRF, stored as a
!>          packed triangular matrix.
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          Details of the interchanges and the block structure of D
!>          as determined by DSPTRF.
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)
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
!> \ingroup doubleOTHERcomputational
!
!  =====================================================================
      SUBROUTINE DSPTRS(Uplo,N,Nrhs,Ap,Ipiv,B,Ldb,Info)
      IMPLICIT NONE
!*--DSPTRS119
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
      DOUBLE PRECISION Ap(*) , B(Ldb,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D+0)
!     ..
!     .. Local Scalars ..
      LOGICAL upper
      INTEGER j , k , kc , kp
      DOUBLE PRECISION ak , akm1 , akm1k , bk , bkm1 , denom
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL DGEMV , DGER , DSCAL , DSWAP , XERBLA
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
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -7
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('DSPTRS',-Info)
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
               IF ( kp/=k ) CALL DSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
!
!           Multiply by inv(U(K)), where U(K) is the transformation
!           stored in column K of A.
!
               CALL DGER(k-1,Nrhs,-ONE,Ap(kc),1,B(k,1),Ldb,B(1,1),Ldb)
!
!           Multiply by the inverse of the diagonal block.
!
               CALL DSCAL(Nrhs,ONE/Ap(kc+k-1),B(k,1),Ldb)
               k = k - 1
            ELSE
!
!           2 x 2 diagonal block
!
!           Interchange rows K-1 and -IPIV(K).
!
               kp = -Ipiv(k)
               IF ( kp/=k-1 ) CALL DSWAP(Nrhs,B(k-1,1),Ldb,B(kp,1),Ldb)
!
!           Multiply by inv(U(K)), where U(K) is the transformation
!           stored in columns K-1 and K of A.
!
               CALL DGER(k-2,Nrhs,-ONE,Ap(kc),1,B(k,1),Ldb,B(1,1),Ldb)
               CALL DGER(k-2,Nrhs,-ONE,Ap(kc-(k-1)),1,B(k-1,1),Ldb,     &
     &                   B(1,1),Ldb)
!
!           Multiply by the inverse of the diagonal block.
!
               akm1k = Ap(kc+k-2)
               akm1 = Ap(kc-1)/akm1k
               ak = Ap(kc+k-1)/akm1k
               denom = akm1*ak - ONE
               DO j = 1 , Nrhs
                  bkm1 = B(k-1,j)/akm1k
                  bk = B(k,j)/akm1k
                  B(k-1,j) = (ak*bkm1-bk)/denom
                  B(k,j) = (akm1*bk-bkm1)/denom
               ENDDO
               kc = kc - k + 1
               k = k - 2
!
            ENDIF
         ENDDO
!
!        Next solve U**T*X = B, overwriting B with X.
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
!           Multiply by inv(U**T(K)), where U(K) is the transformation
!           stored in column K of A.
!
               CALL DGEMV('Transpose',k-1,Nrhs,-ONE,B,Ldb,Ap(kc),1,ONE, &
     &                    B(k,1),Ldb)
!
!           Interchange rows K and IPIV(K).
!
               kp = Ipiv(k)
               IF ( kp/=k ) CALL DSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
               kc = kc + k
               k = k + 1
            ELSE
!
!           2 x 2 diagonal block
!
!           Multiply by inv(U**T(K+1)), where U(K+1) is the transformation
!           stored in columns K and K+1 of A.
!
               CALL DGEMV('Transpose',k-1,Nrhs,-ONE,B,Ldb,Ap(kc),1,ONE, &
     &                    B(k,1),Ldb)
               CALL DGEMV('Transpose',k-1,Nrhs,-ONE,B,Ldb,Ap(kc+k),1,   &
     &                    ONE,B(k+1,1),Ldb)
!
!           Interchange rows K and -IPIV(K).
!
               kp = -Ipiv(k)
               IF ( kp/=k ) CALL DSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
               kc = kc + 2*k + 1
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
               IF ( kp/=k ) CALL DSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
!
!           Multiply by inv(L(K)), where L(K) is the transformation
!           stored in column K of A.
!
               IF ( k<N ) CALL DGER(N-k,Nrhs,-ONE,Ap(kc+1),1,B(k,1),Ldb,&
     &                              B(k+1,1),Ldb)
!
!           Multiply by the inverse of the diagonal block.
!
               CALL DSCAL(Nrhs,ONE/Ap(kc),B(k,1),Ldb)
               kc = kc + N - k + 1
               k = k + 1
            ELSE
!
!           2 x 2 diagonal block
!
!           Interchange rows K+1 and -IPIV(K).
!
               kp = -Ipiv(k)
               IF ( kp/=k+1 ) CALL DSWAP(Nrhs,B(k+1,1),Ldb,B(kp,1),Ldb)
!
!           Multiply by inv(L(K)), where L(K) is the transformation
!           stored in columns K and K+1 of A.
!
               IF ( k<N-1 ) THEN
                  CALL DGER(N-k-1,Nrhs,-ONE,Ap(kc+2),1,B(k,1),Ldb,      &
     &                      B(k+2,1),Ldb)
                  CALL DGER(N-k-1,Nrhs,-ONE,Ap(kc+N-k+2),1,B(k+1,1),Ldb,&
     &                      B(k+2,1),Ldb)
               ENDIF
!
!           Multiply by the inverse of the diagonal block.
!
               akm1k = Ap(kc+1)
               akm1 = Ap(kc)/akm1k
               ak = Ap(kc+N-k+1)/akm1k
               denom = akm1*ak - ONE
               DO j = 1 , Nrhs
                  bkm1 = B(k,j)/akm1k
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
!        Next solve L**T*X = B, overwriting B with X.
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
!           Multiply by inv(L**T(K)), where L(K) is the transformation
!           stored in column K of A.
!
               IF ( k<N ) CALL DGEMV('Transpose',N-k,Nrhs,-ONE,B(k+1,1),&
     &                               Ldb,Ap(kc+1),1,ONE,B(k,1),Ldb)
!
!           Interchange rows K and IPIV(K).
!
               kp = Ipiv(k)
               IF ( kp/=k ) CALL DSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
               k = k - 1
            ELSE
!
!           2 x 2 diagonal block
!
!           Multiply by inv(L**T(K-1)), where L(K-1) is the transformation
!           stored in columns K-1 and K of A.
!
               IF ( k<N ) THEN
                  CALL DGEMV('Transpose',N-k,Nrhs,-ONE,B(k+1,1),Ldb,    &
     &                       Ap(kc+1),1,ONE,B(k,1),Ldb)
                  CALL DGEMV('Transpose',N-k,Nrhs,-ONE,B(k+1,1),Ldb,    &
     &                       Ap(kc-(N-k)),1,ONE,B(k-1,1),Ldb)
               ENDIF
!
!           Interchange rows K and -IPIV(K).
!
               kp = -Ipiv(k)
               IF ( kp/=k ) CALL DSWAP(Nrhs,B(k,1),Ldb,B(kp,1),Ldb)
               kc = kc - (N-k+2)
               k = k - 2
!
            ENDIF
         ENDDO
      ENDIF
!
!
!     End of DSPTRS
!
      END SUBROUTINE DSPTRS
