!*==zhetri2x.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZHETRI2X
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZHETRI2X + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhetri2x.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhetri2x.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhetri2x.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZHETRI2X( UPLO, N, A, LDA, IPIV, WORK, NB, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDA, N, NB
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX*16            A( LDA, * ), WORK( N+NB+1,* )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZHETRI2X computes the inverse of a COMPLEX*16 Hermitian indefinite matrix
!> A using the factorization A = U*D*U**H or A = L*D*L**H computed by
!> ZHETRF.
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
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          On entry, the NNB diagonal matrix D and the multipliers
!>          used to obtain the factor U or L as computed by ZHETRF.
!>
!>          On exit, if INFO = 0, the (symmetric) inverse of the original
!>          matrix.  If UPLO = 'U', the upper triangular part of the
!>          inverse is formed and the part of A below the diagonal is not
!>          referenced; if UPLO = 'L' the lower triangular part of the
!>          inverse is formed and the part of A above the diagonal is
!>          not referenced.
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
!>          Details of the interchanges and the NNB structure of D
!>          as determined by ZHETRF.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (N+NB+1,NB+3)
!> \endverbatim
!>
!> \param[in] NB
!> \verbatim
!>          NB is INTEGER
!>          Block size
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
!  =====================================================================
      SUBROUTINE ZHETRI2X(Uplo,N,A,Lda,Ipiv,Work,Nb,Info)
      USE F77KINDS                        
      USE S_LSAME
      USE S_XERBLA
      USE S_ZGEMM
      USE S_ZHESWAPR
      USE S_ZSYCONV
      USE S_ZTRMM
      USE S_ZTRTRI
      IMPLICIT NONE
!*--ZHETRI2X132
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0) ,       &
     &                 ZERO = (0.0D+0,0.0D+0)
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(N+Nb+1,*) :: Work
      INTEGER , INTENT(IN) :: Nb
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      COMPLEX(CX16KIND) :: ak , akkp1 , akp1 , d , t , u01_ip1_j ,      &
     &                     u01_i_j , u11_ip1_j , u11_i_j
      INTEGER :: count , cut , i , iinfo , invd , ip , j , k , nnb , u11
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
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZHETRI2X',-Info)
         RETURN
      ENDIF
      IF ( N==0 ) RETURN
!
!     Convert A
!     Workspace got Non-diag elements of D
!
      CALL ZSYCONV(Uplo,'C',N,A,Lda,Ipiv,Work,iinfo)
!
!     Check that the diagonal matrix D is nonsingular.
!
      IF ( upper ) THEN
!
!        Upper triangular storage: examine D from bottom to top
!
         DO Info = N , 1 , -1
            IF ( Ipiv(Info)>0 .AND. A(Info,Info)==ZERO ) RETURN
         ENDDO
      ELSE
!
!        Lower triangular storage: examine D from top to bottom.
!
         DO Info = 1 , N
            IF ( Ipiv(Info)>0 .AND. A(Info,Info)==ZERO ) RETURN
         ENDDO
      ENDIF
      Info = 0
!
!  Splitting Workspace
!     U01 is a block (N,NB+1)
!     The first element of U01 is in WORK(1,1)
!     U11 is a block (NB+1,NB+1)
!     The first element of U11 is in WORK(N+1,1)
      u11 = N
!     INVD is a block (N,2)
!     The first element of INVD is in WORK(1,INVD)
      invd = Nb + 2
 
      IF ( upper ) THEN
!
!        invA = P * inv(U**H)*inv(D)*inv(U)*P**H.
!
         CALL ZTRTRI(Uplo,'U',N,A,Lda,Info)
!
!       inv(D) and inv(D)*inv(U)
!
         k = 1
         DO WHILE ( k<=N )
            IF ( Ipiv(k)>0 ) THEN
!           1 x 1 diagonal NNB
               Work(k,invd) = ONE/REAL(A(k,k))
               Work(k,invd+1) = 0
               k = k + 1
            ELSE
!           2 x 2 diagonal NNB
               t = ABS(Work(k+1,1))
               ak = REAL(A(k,k))/t
               akp1 = REAL(A(k+1,k+1))/t
               akkp1 = Work(k+1,1)/t
               d = t*(ak*akp1-ONE)
               Work(k,invd) = akp1/d
               Work(k+1,invd+1) = ak/d
               Work(k,invd+1) = -akkp1/d
               Work(k+1,invd) = DCONJG(Work(k,invd+1))
               k = k + 2
            ENDIF
         ENDDO
!
!       inv(U**H) = (inv(U))**H
!
!       inv(U**H)*inv(D)*inv(U)
!
         cut = N
         DO WHILE ( cut>0 )
            nnb = Nb
            IF ( cut<=nnb ) THEN
               nnb = cut
            ELSE
               count = 0
!             count negative elements,
               DO i = cut + 1 - nnb , cut
                  IF ( Ipiv(i)<0 ) count = count + 1
               ENDDO
!             need a even number for a clear cut
               IF ( MOD(count,2)==1 ) nnb = nnb + 1
            ENDIF
 
            cut = cut - nnb
!
!          U01 Block
!
            DO i = 1 , cut
               DO j = 1 , nnb
                  Work(i,j) = A(i,cut+j)
               ENDDO
            ENDDO
!
!          U11 Block
!
            DO i = 1 , nnb
               Work(u11+i,i) = CONE
               DO j = 1 , i - 1
                  Work(u11+i,j) = ZERO
               ENDDO
               DO j = i + 1 , nnb
                  Work(u11+i,j) = A(cut+i,cut+j)
               ENDDO
            ENDDO
!
!          invD*U01
!
            i = 1
            DO WHILE ( i<=cut )
               IF ( Ipiv(i)>0 ) THEN
                  DO j = 1 , nnb
                     Work(i,j) = Work(i,invd)*Work(i,j)
                  ENDDO
                  i = i + 1
               ELSE
                  DO j = 1 , nnb
                     u01_i_j = Work(i,j)
                     u01_ip1_j = Work(i+1,j)
                     Work(i,j) = Work(i,invd)*u01_i_j + Work(i,invd+1)  &
     &                           *u01_ip1_j
                     Work(i+1,j) = Work(i+1,invd)*u01_i_j +             &
     &                             Work(i+1,invd+1)*u01_ip1_j
                  ENDDO
                  i = i + 2
               ENDIF
            ENDDO
!
!        invD1*U11
!
            i = 1
            DO WHILE ( i<=nnb )
               IF ( Ipiv(cut+i)>0 ) THEN
                  DO j = i , nnb
                     Work(u11+i,j) = Work(cut+i,invd)*Work(u11+i,j)
                  ENDDO
                  i = i + 1
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
                  i = i + 2
               ENDIF
            ENDDO
!
!       U11**H*invD1*U11->U11
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
!          U01**H*invD*U01->A(CUT+I,CUT+J)
!
            CALL ZGEMM('C','N',nnb,nnb,cut,CONE,A(1,cut+1),Lda,Work,    &
     &                 N+Nb+1,ZERO,Work(u11+1,1),N+Nb+1)
!
!        U11 =  U11**H*invD1*U11 + U01**H*invD*U01
!
            DO i = 1 , nnb
               DO j = i , nnb
                  A(cut+i,cut+j) = A(cut+i,cut+j) + Work(u11+i,j)
               ENDDO
            ENDDO
!
!        U01 =  U00**H*invD0*U01
!
            CALL ZTRMM('L',Uplo,'C','U',cut,nnb,CONE,A,Lda,Work,N+Nb+1)
 
!
!        Update U01
!
            DO i = 1 , cut
               DO j = 1 , nnb
                  A(i,cut+j) = Work(i,j)
               ENDDO
            ENDDO
!
!      Next Block
!
         ENDDO
!
!        Apply PERMUTATIONS P and P**H: P * inv(U**H)*inv(D)*inv(U) *P**H
!
         i = 1
         DO WHILE ( i<=N )
            IF ( Ipiv(i)>0 ) THEN
               ip = Ipiv(i)
               IF ( i<ip ) CALL ZHESWAPR(Uplo,N,A,Lda,i,ip)
               IF ( i>ip ) CALL ZHESWAPR(Uplo,N,A,Lda,ip,i)
            ELSE
               ip = -Ipiv(i)
               i = i + 1
               IF ( (i-1)<ip ) CALL ZHESWAPR(Uplo,N,A,Lda,i-1,ip)
               IF ( (i-1)>ip ) CALL ZHESWAPR(Uplo,N,A,Lda,ip,i-1)
            ENDIF
            i = i + 1
         ENDDO
      ELSE
!
!        LOWER...
!
!        invA = P * inv(U**H)*inv(D)*inv(U)*P**H.
!
         CALL ZTRTRI(Uplo,'U',N,A,Lda,Info)
!
!       inv(D) and inv(D)*inv(U)
!
         k = N
         DO WHILE ( k>=1 )
            IF ( Ipiv(k)>0 ) THEN
!           1 x 1 diagonal NNB
               Work(k,invd) = ONE/REAL(A(k,k))
               Work(k,invd+1) = 0
               k = k - 1
            ELSE
!           2 x 2 diagonal NNB
               t = ABS(Work(k-1,1))
               ak = REAL(A(k-1,k-1))/t
               akp1 = REAL(A(k,k))/t
               akkp1 = Work(k-1,1)/t
               d = t*(ak*akp1-ONE)
               Work(k-1,invd) = akp1/d
               Work(k,invd) = ak/d
               Work(k,invd+1) = -akkp1/d
               Work(k-1,invd+1) = DCONJG(Work(k,invd+1))
               k = k - 2
            ENDIF
         ENDDO
!
!       inv(U**H) = (inv(U))**H
!
!       inv(U**H)*inv(D)*inv(U)
!
         cut = 0
         DO WHILE ( cut<N )
            nnb = Nb
            IF ( cut+nnb>=N ) THEN
               nnb = N - cut
            ELSE
               count = 0
!             count negative elements,
               DO i = cut + 1 , cut + nnb
                  IF ( Ipiv(i)<0 ) count = count + 1
               ENDDO
!             need a even number for a clear cut
               IF ( MOD(count,2)==1 ) nnb = nnb + 1
            ENDIF
!      L21 Block
            DO i = 1 , N - cut - nnb
               DO j = 1 , nnb
                  Work(i,j) = A(cut+nnb+i,cut+j)
               ENDDO
            ENDDO
!     L11 Block
            DO i = 1 , nnb
               Work(u11+i,i) = CONE
               DO j = i + 1 , nnb
                  Work(u11+i,j) = ZERO
               ENDDO
               DO j = 1 , i - 1
                  Work(u11+i,j) = A(cut+i,cut+j)
               ENDDO
            ENDDO
!
!          invD*L21
!
            i = N - cut - nnb
            DO WHILE ( i>=1 )
               IF ( Ipiv(cut+nnb+i)>0 ) THEN
                  DO j = 1 , nnb
                     Work(i,j) = Work(cut+nnb+i,invd)*Work(i,j)
                  ENDDO
                  i = i - 1
               ELSE
                  DO j = 1 , nnb
                     u01_i_j = Work(i,j)
                     u01_ip1_j = Work(i-1,j)
                     Work(i,j) = Work(cut+nnb+i,invd)*u01_i_j +         &
     &                           Work(cut+nnb+i,invd+1)*u01_ip1_j
                     Work(i-1,j) = Work(cut+nnb+i-1,invd+1)*u01_i_j +   &
     &                             Work(cut+nnb+i-1,invd)*u01_ip1_j
                  ENDDO
                  i = i - 2
               ENDIF
            ENDDO
!
!        invD1*L11
!
            i = nnb
            DO WHILE ( i>=1 )
               IF ( Ipiv(cut+i)>0 ) THEN
                  DO j = 1 , nnb
                     Work(u11+i,j) = Work(cut+i,invd)*Work(u11+i,j)
                  ENDDO
                  i = i - 1
               ELSE
                  DO j = 1 , nnb
                     u11_i_j = Work(u11+i,j)
                     u11_ip1_j = Work(u11+i-1,j)
                     Work(u11+i,j) = Work(cut+i,invd)*Work(u11+i,j)     &
     &                               + Work(cut+i,invd+1)*u11_ip1_j
                     Work(u11+i-1,j) = Work(cut+i-1,invd+1)*u11_i_j +   &
     &                                 Work(cut+i-1,invd)*u11_ip1_j
                  ENDDO
                  i = i - 2
               ENDIF
            ENDDO
!
!       L11**H*invD1*L11->L11
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
!          L21**H*invD2*L21->A(CUT+I,CUT+J)
!
               CALL ZGEMM('C','N',nnb,nnb,N-nnb-cut,CONE,               &
     &                    A(cut+nnb+1,cut+1),Lda,Work,N+Nb+1,ZERO,      &
     &                    Work(u11+1,1),N+Nb+1)
 
!
!        L11 =  L11**H*invD1*L11 + U01**H*invD*U01
!
               DO i = 1 , nnb
                  DO j = 1 , i
                     A(cut+i,cut+j) = A(cut+i,cut+j) + Work(u11+i,j)
                  ENDDO
               ENDDO
!
!        L01 =  L22**H*invD2*L21
!
               CALL ZTRMM('L',Uplo,'C','U',N-nnb-cut,nnb,CONE,          &
     &                    A(cut+nnb+1,cut+nnb+1),Lda,Work,N+Nb+1)
 
!      Update L21
               DO i = 1 , N - cut - nnb
                  DO j = 1 , nnb
                     A(cut+nnb+i,cut+j) = Work(i,j)
                  ENDDO
               ENDDO
            ELSE
!
!        L11 =  L11**H*invD1*L11
!
               DO i = 1 , nnb
                  DO j = 1 , i
                     A(cut+i,cut+j) = Work(u11+i,j)
                  ENDDO
               ENDDO
            ENDIF
!
!      Next Block
!
            cut = cut + nnb
         ENDDO
!
!        Apply PERMUTATIONS P and P**H: P * inv(U**H)*inv(D)*inv(U) *P**H
!
         i = N
         DO WHILE ( i>=1 )
            IF ( Ipiv(i)>0 ) THEN
               ip = Ipiv(i)
               IF ( i<ip ) CALL ZHESWAPR(Uplo,N,A,Lda,i,ip)
               IF ( i>ip ) CALL ZHESWAPR(Uplo,N,A,Lda,ip,i)
            ELSE
               ip = -Ipiv(i)
               IF ( i<ip ) CALL ZHESWAPR(Uplo,N,A,Lda,i,ip)
               IF ( i>ip ) CALL ZHESWAPR(Uplo,N,A,Lda,ip,i)
               i = i - 1
            ENDIF
            i = i - 1
         ENDDO
      ENDIF
!
!
!     End of ZHETRI2X
!
      END SUBROUTINE ZHETRI2X
