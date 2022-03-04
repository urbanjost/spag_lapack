!*==zhetri_3x.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZHETRI_3X
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZHETRI_3X + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhetri_3x.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhetri_3x.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhetri_3x.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZHETRI_3X( UPLO, N, A, LDA, E, IPIV, WORK, NB, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDA, N, NB
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX*16         A( LDA, * ),  E( * ), WORK( N+NB+1, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!> ZHETRI_3X computes the inverse of a complex Hermitian indefinite
!> matrix A using the factorization computed by ZHETRF_RK or ZHETRF_BK:
!>
!>     A = P*U*D*(U**H)*(P**T) or A = P*L*D*(L**H)*(P**T),
!>
!> where U (or L) is unit upper (or lower) triangular matrix,
!> U**H (or L**H) is the conjugate of U (or L), P is a permutation
!> matrix, P**T is the transpose of P, and D is Hermitian and block
!> diagonal with 1-by-1 and 2-by-2 diagonal blocks.
!>
!> This is the blocked version of the algorithm, calling Level 3 BLAS.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the details of the factorization are
!>          stored as an upper or lower triangular matrix.
!>          = 'U':  Upper triangle of A is stored;
!>          = 'L':  Lower triangle of A is stored.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          On entry, diagonal of the block diagonal matrix D and
!>          factors U or L as computed by ZHETRF_RK and ZHETRF_BK:
!>            a) ONLY diagonal elements of the Hermitian block diagonal
!>               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k);
!>               (superdiagonal (or subdiagonal) elements of D
!>                should be provided on entry in array E), and
!>            b) If UPLO = 'U': factor U in the superdiagonal part of A.
!>               If UPLO = 'L': factor L in the subdiagonal part of A.
!>
!>          On exit, if INFO = 0, the Hermitian inverse of the original
!>          matrix.
!>             If UPLO = 'U': the upper triangular part of the inverse
!>             is formed and the part of A below the diagonal is not
!>             referenced;
!>             If UPLO = 'L': the lower triangular part of the inverse
!>             is formed and the part of A above the diagonal is not
!>             referenced.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is COMPLEX*16 array, dimension (N)
!>          On entry, contains the superdiagonal (or subdiagonal)
!>          elements of the Hermitian block diagonal matrix D
!>          with 1-by-1 or 2-by-2 diagonal blocks, where
!>          If UPLO = 'U': E(i) = D(i-1,i), i=2:N, E(1) not referenced;
!>          If UPLO = 'L': E(i) = D(i+1,i), i=1:N-1, E(N) not referenced.
!>
!>          NOTE: For 1-by-1 diagonal block D(k), where
!>          1 <= k <= N, the element E(k) is not referenced in both
!>          UPLO = 'U' or UPLO = 'L' cases.
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          Details of the interchanges and the block structure of D
!>          as determined by ZHETRF_RK or ZHETRF_BK.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (N+NB+1,NB+3).
!> \endverbatim
!>
!> \param[in] NB
!> \verbatim
!>          NB is INTEGER
!>          Block size.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
!>          > 0: if INFO = i, D(i,i) = 0; the matrix is singular and its
!>               inverse could not be computed.
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
!> \date June 2017
!
!> \ingroup complex16HEcomputational
!
!> \par Contributors:
!  ==================
!> \verbatim
!>
!>  June 2017,  Igor Kozachenko,
!>                  Computer Science Division,
!>                  University of California, Berkeley
!>
!> \endverbatim
!
!  =====================================================================
      SUBROUTINE ZHETRI_3X(Uplo,N,A,Lda,E,Ipiv,Work,Nb,Info)
      USE F77KINDS                        
      USE S_LSAME
      USE S_XERBLA
      USE S_ZGEMM
      USE S_ZHESWAPR
      USE S_ZTRMM
      USE S_ZTRTRI
      IMPLICIT NONE
!*--ZHETRI_3X170
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0) ,       &
     &                 CZERO = (0.0D+0,0.0D+0)
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: E
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(N+Nb+1,*) :: Work
      INTEGER , INTENT(IN) :: Nb
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: ak , akp1 , t
      COMPLEX(CX16KIND) :: akkp1 , d , u01_ip1_j , u01_i_j , u11_ip1_j ,&
     &                     u11_i_j
      INTEGER :: cut , i , icount , invd , ip , j , k , nnb , u11
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
!     Test the input parameters.
!
      Info = 0
      upper = LSAME(Uplo,'U')
      IF ( .NOT.upper .AND. .NOT.LSAME(Uplo,'L') ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -4
      ENDIF
!
!     Quick return if possible
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZHETRI_3X',-Info)
         RETURN
      ENDIF
      IF ( N==0 ) RETURN
!
!     Workspace got Non-diag elements of D
!
      DO k = 1 , N
         Work(k,1) = E(k)
      ENDDO
!
!     Check that the diagonal matrix D is nonsingular.
!
      IF ( upper ) THEN
!
!        Upper triangular storage: examine D from bottom to top
!
         DO Info = N , 1 , -1
            IF ( Ipiv(Info)>0 .AND. A(Info,Info)==CZERO ) RETURN
         ENDDO
      ELSE
!
!        Lower triangular storage: examine D from top to bottom.
!
         DO Info = 1 , N
            IF ( Ipiv(Info)>0 .AND. A(Info,Info)==CZERO ) RETURN
         ENDDO
      ENDIF
!
      Info = 0
!
!     Splitting Workspace
!     U01 is a block ( N, NB+1 )
!     The first element of U01 is in WORK( 1, 1 )
!     U11 is a block ( NB+1, NB+1 )
!     The first element of U11 is in WORK( N+1, 1 )
!
      u11 = N
!
!     INVD is a block ( N, 2 )
!     The first element of INVD is in WORK( 1, INVD )
!
      invd = Nb + 2
 
      IF ( upper ) THEN
!
!        Begin Upper
!
!        invA = P * inv(U**H) * inv(D) * inv(U) * P**T.
!
         CALL ZTRTRI(Uplo,'U',N,A,Lda,Info)
!
!        inv(D) and inv(D) * inv(U)
!
         k = 1
         DO WHILE ( k<=N )
            IF ( Ipiv(k)>0 ) THEN
!              1 x 1 diagonal NNB
               Work(k,invd) = ONE/DBLE(A(k,k))
               Work(k,invd+1) = CZERO
            ELSE
!              2 x 2 diagonal NNB
               t = ABS(Work(k+1,1))
               ak = DBLE(A(k,k))/t
               akp1 = DBLE(A(k+1,k+1))/t
               akkp1 = Work(k+1,1)/t
               d = t*(ak*akp1-CONE)
               Work(k,invd) = akp1/d
               Work(k+1,invd+1) = ak/d
               Work(k,invd+1) = -akkp1/d
               Work(k+1,invd) = DCONJG(Work(k,invd+1))
               k = k + 1
            ENDIF
            k = k + 1
         ENDDO
!
!        inv(U**H) = (inv(U))**H
!
!        inv(U**H) * inv(D) * inv(U)
!
         cut = N
         DO WHILE ( cut>0 )
            nnb = Nb
            IF ( cut<=nnb ) THEN
               nnb = cut
            ELSE
               icount = 0
!              count negative elements,
               DO i = cut + 1 - nnb , cut
                  IF ( Ipiv(i)<0 ) icount = icount + 1
               ENDDO
!              need a even number for a clear cut
               IF ( MOD(icount,2)==1 ) nnb = nnb + 1
            ENDIF
 
            cut = cut - nnb
!
!           U01 Block
!
            DO i = 1 , cut
               DO j = 1 , nnb
                  Work(i,j) = A(i,cut+j)
               ENDDO
            ENDDO
!
!           U11 Block
!
            DO i = 1 , nnb
               Work(u11+i,i) = CONE
               DO j = 1 , i - 1
                  Work(u11+i,j) = CZERO
               ENDDO
               DO j = i + 1 , nnb
                  Work(u11+i,j) = A(cut+i,cut+j)
               ENDDO
            ENDDO
!
!           invD * U01
!
            i = 1
            DO WHILE ( i<=cut )
               IF ( Ipiv(i)>0 ) THEN
                  DO j = 1 , nnb
                     Work(i,j) = Work(i,invd)*Work(i,j)
                  ENDDO
               ELSE
                  DO j = 1 , nnb
                     u01_i_j = Work(i,j)
                     u01_ip1_j = Work(i+1,j)
                     Work(i,j) = Work(i,invd)*u01_i_j + Work(i,invd+1)  &
     &                           *u01_ip1_j
                     Work(i+1,j) = Work(i+1,invd)*u01_i_j +             &
     &                             Work(i+1,invd+1)*u01_ip1_j
                  ENDDO
                  i = i + 1
               ENDIF
               i = i + 1
            ENDDO
!
!           invD1 * U11
!
            i = 1
            DO WHILE ( i<=nnb )
               IF ( Ipiv(cut+i)>0 ) THEN
                  DO j = i , nnb
                     Work(u11+i,j) = Work(cut+i,invd)*Work(u11+i,j)
                  ENDDO
               ELSE
                  DO j = i , nnb
                     u11_i_j = Work(u11+i,j)
                     u11_ip1_j = Work(u11+i+1,j)
                     Work(u11+i,j) = Work(cut+i,invd)*Work(u11+i,j)     &
     &                               + Work(cut+i,invd+1)               &
     &                               *Work(u11+i+1,j)
                     Work(u11+i+1,j) = Work(cut+i+1,invd)*u11_i_j +     &
     &                                 Work(cut+i+1,invd+1)*u11_ip1_j
                  ENDDO
                  i = i + 1
               ENDIF
               i = i + 1
            ENDDO
!
!           U11**H * invD1 * U11 -> U11
!
            CALL ZTRMM('L','U','C','U',nnb,nnb,CONE,A(cut+1,cut+1),Lda, &
     &                 Work(u11+1,1),N+Nb+1)
!
            DO i = 1 , nnb
               DO j = i , nnb
                  A(cut+i,cut+j) = Work(u11+i,j)
               ENDDO
            ENDDO
!
!           U01**H * invD * U01 -> A( CUT+I, CUT+J )
!
            CALL ZGEMM('C','N',nnb,nnb,cut,CONE,A(1,cut+1),Lda,Work,    &
     &                 N+Nb+1,CZERO,Work(u11+1,1),N+Nb+1)
 
!
!           U11 =  U11**H * invD1 * U11 + U01**H * invD * U01
!
            DO i = 1 , nnb
               DO j = i , nnb
                  A(cut+i,cut+j) = A(cut+i,cut+j) + Work(u11+i,j)
               ENDDO
            ENDDO
!
!           U01 =  U00**H * invD0 * U01
!
            CALL ZTRMM('L',Uplo,'C','U',cut,nnb,CONE,A,Lda,Work,N+Nb+1)
 
!
!           Update U01
!
            DO i = 1 , cut
               DO j = 1 , nnb
                  A(i,cut+j) = Work(i,j)
               ENDDO
            ENDDO
!
!           Next Block
!
         ENDDO
!
!        Apply PERMUTATIONS P and P**T:
!        P * inv(U**H) * inv(D) * inv(U) * P**T.
!        Interchange rows and columns I and IPIV(I) in reverse order
!        from the formation order of IPIV vector for Upper case.
!
!        ( We can use a loop over IPIV with increment 1,
!        since the ABS value of IPIV(I) represents the row (column)
!        index of the interchange with row (column) i in both 1x1
!        and 2x2 pivot cases, i.e. we don't need separate code branches
!        for 1x1 and 2x2 pivot cases )
!
         DO i = 1 , N
            ip = ABS(Ipiv(i))
            IF ( ip/=i ) THEN
               IF ( i<ip ) CALL ZHESWAPR(Uplo,N,A,Lda,i,ip)
               IF ( i>ip ) CALL ZHESWAPR(Uplo,N,A,Lda,ip,i)
            ENDIF
         ENDDO
!
      ELSE
!
!        Begin Lower
!
!        inv A = P * inv(L**H) * inv(D) * inv(L) * P**T.
!
         CALL ZTRTRI(Uplo,'U',N,A,Lda,Info)
!
!        inv(D) and inv(D) * inv(L)
!
         k = N
         DO WHILE ( k>=1 )
            IF ( Ipiv(k)>0 ) THEN
!              1 x 1 diagonal NNB
               Work(k,invd) = ONE/DBLE(A(k,k))
               Work(k,invd+1) = CZERO
            ELSE
!              2 x 2 diagonal NNB
               t = ABS(Work(k-1,1))
               ak = DBLE(A(k-1,k-1))/t
               akp1 = DBLE(A(k,k))/t
               akkp1 = Work(k-1,1)/t
               d = t*(ak*akp1-CONE)
               Work(k-1,invd) = akp1/d
               Work(k,invd) = ak/d
               Work(k,invd+1) = -akkp1/d
               Work(k-1,invd+1) = DCONJG(Work(k,invd+1))
               k = k - 1
            ENDIF
            k = k - 1
         ENDDO
!
!        inv(L**H) = (inv(L))**H
!
!        inv(L**H) * inv(D) * inv(L)
!
         cut = 0
         DO WHILE ( cut<N )
            nnb = Nb
            IF ( (cut+nnb)>N ) THEN
               nnb = N - cut
            ELSE
               icount = 0
!              count negative elements,
               DO i = cut + 1 , cut + nnb
                  IF ( Ipiv(i)<0 ) icount = icount + 1
               ENDDO
!              need a even number for a clear cut
               IF ( MOD(icount,2)==1 ) nnb = nnb + 1
            ENDIF
!
!           L21 Block
!
            DO i = 1 , N - cut - nnb
               DO j = 1 , nnb
                  Work(i,j) = A(cut+nnb+i,cut+j)
               ENDDO
            ENDDO
!
!           L11 Block
!
            DO i = 1 , nnb
               Work(u11+i,i) = CONE
               DO j = i + 1 , nnb
                  Work(u11+i,j) = CZERO
               ENDDO
               DO j = 1 , i - 1
                  Work(u11+i,j) = A(cut+i,cut+j)
               ENDDO
            ENDDO
!
!           invD*L21
!
            i = N - cut - nnb
            DO WHILE ( i>=1 )
               IF ( Ipiv(cut+nnb+i)>0 ) THEN
                  DO j = 1 , nnb
                     Work(i,j) = Work(cut+nnb+i,invd)*Work(i,j)
                  ENDDO
               ELSE
                  DO j = 1 , nnb
                     u01_i_j = Work(i,j)
                     u01_ip1_j = Work(i-1,j)
                     Work(i,j) = Work(cut+nnb+i,invd)*u01_i_j +         &
     &                           Work(cut+nnb+i,invd+1)*u01_ip1_j
                     Work(i-1,j) = Work(cut+nnb+i-1,invd+1)*u01_i_j +   &
     &                             Work(cut+nnb+i-1,invd)*u01_ip1_j
                  ENDDO
                  i = i - 1
               ENDIF
               i = i - 1
            ENDDO
!
!           invD1*L11
!
            i = nnb
            DO WHILE ( i>=1 )
               IF ( Ipiv(cut+i)>0 ) THEN
                  DO j = 1 , nnb
                     Work(u11+i,j) = Work(cut+i,invd)*Work(u11+i,j)
                  ENDDO
 
               ELSE
                  DO j = 1 , nnb
                     u11_i_j = Work(u11+i,j)
                     u11_ip1_j = Work(u11+i-1,j)
                     Work(u11+i,j) = Work(cut+i,invd)*Work(u11+i,j)     &
     &                               + Work(cut+i,invd+1)*u11_ip1_j
                     Work(u11+i-1,j) = Work(cut+i-1,invd+1)*u11_i_j +   &
     &                                 Work(cut+i-1,invd)*u11_ip1_j
                  ENDDO
                  i = i - 1
               ENDIF
               i = i - 1
            ENDDO
!
!           L11**H * invD1 * L11 -> L11
!
            CALL ZTRMM('L',Uplo,'C','U',nnb,nnb,CONE,A(cut+1,cut+1),Lda,&
     &                 Work(u11+1,1),N+Nb+1)
 
!
            DO i = 1 , nnb
               DO j = 1 , i
                  A(cut+i,cut+j) = Work(u11+i,j)
               ENDDO
            ENDDO
!
            IF ( (cut+nnb)<N ) THEN
!
!              L21**H * invD2*L21 -> A( CUT+I, CUT+J )
!
               CALL ZGEMM('C','N',nnb,nnb,N-nnb-cut,CONE,               &
     &                    A(cut+nnb+1,cut+1),Lda,Work,N+Nb+1,CZERO,     &
     &                    Work(u11+1,1),N+Nb+1)
 
!
!              L11 =  L11**H * invD1 * L11 + U01**H * invD * U01
!
               DO i = 1 , nnb
                  DO j = 1 , i
                     A(cut+i,cut+j) = A(cut+i,cut+j) + Work(u11+i,j)
                  ENDDO
               ENDDO
!
!              L01 =  L22**H * invD2 * L21
!
               CALL ZTRMM('L',Uplo,'C','U',N-nnb-cut,nnb,CONE,          &
     &                    A(cut+nnb+1,cut+nnb+1),Lda,Work,N+Nb+1)
!
!              Update L21
!
               DO i = 1 , N - cut - nnb
                  DO j = 1 , nnb
                     A(cut+nnb+i,cut+j) = Work(i,j)
                  ENDDO
               ENDDO
!
            ELSE
!
!              L11 =  L11**H * invD1 * L11
!
               DO i = 1 , nnb
                  DO j = 1 , i
                     A(cut+i,cut+j) = Work(u11+i,j)
                  ENDDO
               ENDDO
            ENDIF
!
!           Next Block
!
            cut = cut + nnb
!
         ENDDO
!
!        Apply PERMUTATIONS P and P**T:
!        P * inv(L**H) * inv(D) * inv(L) * P**T.
!        Interchange rows and columns I and IPIV(I) in reverse order
!        from the formation order of IPIV vector for Lower case.
!
!        ( We can use a loop over IPIV with increment -1,
!        since the ABS value of IPIV(I) represents the row (column)
!        index of the interchange with row (column) i in both 1x1
!        and 2x2 pivot cases, i.e. we don't need separate code branches
!        for 1x1 and 2x2 pivot cases )
!
         DO i = N , 1 , -1
            ip = ABS(Ipiv(i))
            IF ( ip/=i ) THEN
               IF ( i<ip ) CALL ZHESWAPR(Uplo,N,A,Lda,i,ip)
               IF ( i>ip ) CALL ZHESWAPR(Uplo,N,A,Lda,ip,i)
            ENDIF
         ENDDO
!
      ENDIF
!
!
!     End of ZHETRI_3X
!
      END SUBROUTINE ZHETRI_3X
