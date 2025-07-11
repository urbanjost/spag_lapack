!*==dsytri2x.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b DSYTRI2X
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DSYTRI2X + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsytri2x.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsytri2x.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsytri2x.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSYTRI2X( UPLO, N, A, LDA, IPIV, WORK, NB, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDA, N, NB
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       DOUBLE PRECISION   A( LDA, * ), WORK( N+NB+1,* )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DSYTRI2X computes the inverse of a real symmetric indefinite matrix
!> A using the factorization A = U*D*U**T or A = L*D*L**T computed by
!> DSYTRF.
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
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, the NNB diagonal matrix D and the multipliers
!>          used to obtain the factor U or L as computed by DSYTRF.
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
!>          as determined by DSYTRF.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (N+NB+1,NB+3)
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
!> \ingroup doubleSYcomputational
!
!  =====================================================================
      SUBROUTINE DSYTRI2X(Uplo,N,A,Lda,Ipiv,Work,Nb,Info)
      IMPLICIT NONE
!*--DSYTRI2X124
!
!  -- LAPACK computational routine (version 3.7.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2017
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER Info , Lda , N , Nb
!     ..
!     .. Array Arguments ..
      INTEGER Ipiv(*)
      DOUBLE PRECISION A(Lda,*) , Work(N+Nb+1,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ONE , ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
      LOGICAL upper
      INTEGER i , iinfo , ip , k , cut , nnb
      INTEGER count
      INTEGER j , u11 , invd
 
      DOUBLE PRECISION ak , akkp1 , akp1 , d , t
      DOUBLE PRECISION u01_i_j , u01_ip1_j
      DOUBLE PRECISION u11_i_j , u11_ip1_j
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL DSYCONV , XERBLA , DTRTRI
      EXTERNAL DGEMM , DTRMM , DSYSWAPR
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX
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
         CALL XERBLA('DSYTRI2X',-Info)
         RETURN
      ENDIF
      IF ( N==0 ) RETURN
!
!     Convert A
!     Workspace got Non-diag elements of D
!
      CALL DSYCONV(Uplo,'C',N,A,Lda,Ipiv,Work,iinfo)
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
!        invA = P * inv(U**T)*inv(D)*inv(U)*P**T.
!
         CALL DTRTRI(Uplo,'U',N,A,Lda,Info)
!
!       inv(D) and inv(D)*inv(U)
!
         k = 1
         DO WHILE ( k<=N )
            IF ( Ipiv(k)>0 ) THEN
!           1 x 1 diagonal NNB
               Work(k,invd) = ONE/A(k,k)
               Work(k,invd+1) = 0
               k = k + 1
            ELSE
!           2 x 2 diagonal NNB
               t = Work(k+1,1)
               ak = A(k,k)/t
               akp1 = A(k+1,k+1)/t
               akkp1 = Work(k+1,1)/t
               d = t*(ak*akp1-ONE)
               Work(k,invd) = akp1/d
               Work(k+1,invd+1) = ak/d
               Work(k,invd+1) = -akkp1/d
               Work(k+1,invd) = -akkp1/d
               k = k + 2
            ENDIF
         ENDDO
!
!       inv(U**T) = (inv(U))**T
!
!       inv(U**T)*inv(D)*inv(U)
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
               Work(u11+i,i) = ONE
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
!       U11**T*invD1*U11->U11
!
            CALL DTRMM('L','U','T','U',nnb,nnb,ONE,A(cut+1,cut+1),Lda,  &
     &                 Work(u11+1,1),N+Nb+1)
!
            DO i = 1 , nnb
               DO j = i , nnb
                  A(cut+i,cut+j) = Work(u11+i,j)
               ENDDO
            ENDDO
!
!          U01**T*invD*U01->A(CUT+I,CUT+J)
!
            CALL DGEMM('T','N',nnb,nnb,cut,ONE,A(1,cut+1),Lda,Work,     &
     &                 N+Nb+1,ZERO,Work(u11+1,1),N+Nb+1)
 
!
!        U11 =  U11**T*invD1*U11 + U01**T*invD*U01
!
            DO i = 1 , nnb
               DO j = i , nnb
                  A(cut+i,cut+j) = A(cut+i,cut+j) + Work(u11+i,j)
               ENDDO
            ENDDO
!
!        U01 =  U00**T*invD0*U01
!
            CALL DTRMM('L',Uplo,'T','U',cut,nnb,ONE,A,Lda,Work,N+Nb+1)
 
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
!        Apply PERMUTATIONS P and P**T: P * inv(U**T)*inv(D)*inv(U) *P**T
!
         i = 1
         DO WHILE ( i<=N )
            IF ( Ipiv(i)>0 ) THEN
               ip = Ipiv(i)
               IF ( i<ip ) CALL DSYSWAPR(Uplo,N,A,Lda,i,ip)
               IF ( i>ip ) CALL DSYSWAPR(Uplo,N,A,Lda,ip,i)
            ELSE
               ip = -Ipiv(i)
               i = i + 1
               IF ( (i-1)<ip ) CALL DSYSWAPR(Uplo,N,A,Lda,i-1,ip)
               IF ( (i-1)>ip ) CALL DSYSWAPR(Uplo,N,A,Lda,ip,i-1)
            ENDIF
            i = i + 1
         ENDDO
      ELSE
!
!        LOWER...
!
!        invA = P * inv(U**T)*inv(D)*inv(U)*P**T.
!
         CALL DTRTRI(Uplo,'U',N,A,Lda,Info)
!
!       inv(D) and inv(D)*inv(U)
!
         k = N
         DO WHILE ( k>=1 )
            IF ( Ipiv(k)>0 ) THEN
!           1 x 1 diagonal NNB
               Work(k,invd) = ONE/A(k,k)
               Work(k,invd+1) = 0
               k = k - 1
            ELSE
!           2 x 2 diagonal NNB
               t = Work(k-1,1)
               ak = A(k-1,k-1)/t
               akp1 = A(k,k)/t
               akkp1 = Work(k-1,1)/t
               d = t*(ak*akp1-ONE)
               Work(k-1,invd) = akp1/d
               Work(k,invd) = ak/d
               Work(k,invd+1) = -akkp1/d
               Work(k-1,invd+1) = -akkp1/d
               k = k - 2
            ENDIF
         ENDDO
!
!       inv(U**T) = (inv(U))**T
!
!       inv(U**T)*inv(D)*inv(U)
!
         cut = 0
         DO WHILE ( cut<N )
            nnb = Nb
            IF ( cut+nnb>N ) THEN
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
!     L21 Block
            DO i = 1 , N - cut - nnb
               DO j = 1 , nnb
                  Work(i,j) = A(cut+nnb+i,cut+j)
               ENDDO
            ENDDO
!     L11 Block
            DO i = 1 , nnb
               Work(u11+i,i) = ONE
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
!       L11**T*invD1*L11->L11
!
            CALL DTRMM('L',Uplo,'T','U',nnb,nnb,ONE,A(cut+1,cut+1),Lda, &
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
!          L21**T*invD2*L21->A(CUT+I,CUT+J)
!
               CALL DGEMM('T','N',nnb,nnb,N-nnb-cut,ONE,                &
     &                    A(cut+nnb+1,cut+1),Lda,Work,N+Nb+1,ZERO,      &
     &                    Work(u11+1,1),N+Nb+1)
 
!
!        L11 =  L11**T*invD1*L11 + U01**T*invD*U01
!
               DO i = 1 , nnb
                  DO j = 1 , i
                     A(cut+i,cut+j) = A(cut+i,cut+j) + Work(u11+i,j)
                  ENDDO
               ENDDO
!
!        L01 =  L22**T*invD2*L21
!
               CALL DTRMM('L',Uplo,'T','U',N-nnb-cut,nnb,ONE,           &
     &                    A(cut+nnb+1,cut+nnb+1),Lda,Work,N+Nb+1)
!
!      Update L21
!
               DO i = 1 , N - cut - nnb
                  DO j = 1 , nnb
                     A(cut+nnb+i,cut+j) = Work(i,j)
                  ENDDO
               ENDDO
 
            ELSE
!
!        L11 =  L11**T*invD1*L11
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
!        Apply PERMUTATIONS P and P**T: P * inv(U**T)*inv(D)*inv(U) *P**T
!
         i = N
         DO WHILE ( i>=1 )
            IF ( Ipiv(i)>0 ) THEN
               ip = Ipiv(i)
               IF ( i<ip ) CALL DSYSWAPR(Uplo,N,A,Lda,i,ip)
               IF ( i>ip ) CALL DSYSWAPR(Uplo,N,A,Lda,ip,i)
            ELSE
               ip = -Ipiv(i)
               IF ( i<ip ) CALL DSYSWAPR(Uplo,N,A,Lda,i,ip)
               IF ( i>ip ) CALL DSYSWAPR(Uplo,N,A,Lda,ip,i)
               i = i - 1
            ENDIF
            i = i - 1
         ENDDO
      ENDIF
!
!
!     End of DSYTRI2X
!
      END SUBROUTINE DSYTRI2X
