!*==ssytrf_aa_2stage.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SSYTRF_AA_2STAGE
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SSYTRF_AA_2STAGE + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssytrf_aa_2stage.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssytrf_aa_2stage.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssytrf_aa_2stage.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!      SUBROUTINE SSYTRF_AA_2STAGE( UPLO, N, A, LDA, TB, LTB, IPIV,
!                                   IPIV2, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            N, LDA, LTB, LWORK, INFO
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * ), IPIV2( * )
!       REAL               A( LDA, * ), TB( * ), WORK( * )
!       ..
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SSYTRF_AA_2STAGE computes the factorization of a real symmetric matrix A
!> using the Aasen's algorithm.  The form of the factorization is
!>
!>    A = U**T*T*U  or  A = L*T*L**T
!>
!> where U (or L) is a product of permutation and unit upper (lower)
!> triangular matrices, and T is a symmetric band matrix with the
!> bandwidth of NB (NB is internally selected and stored in TB( 1 ), and T is
!> LU factorized with partial pivoting).
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
!>          A is REAL array, dimension (LDA,N)
!>          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!>          N-by-N upper triangular part of A contains the upper
!>          triangular part of the matrix A, and the strictly lower
!>          triangular part of A is not referenced.  If UPLO = 'L', the
!>          leading N-by-N lower triangular part of A contains the lower
!>          triangular part of the matrix A, and the strictly upper
!>          triangular part of A is not referenced.
!>
!>          On exit, L is stored below (or above) the subdiaonal blocks,
!>          when UPLO  is 'L' (or 'U').
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] TB
!> \verbatim
!>          TB is REAL array, dimension (LTB)
!>          On exit, details of the LU factorization of the band matrix.
!> \endverbatim
!>
!> \param[in] LTB
!> \verbatim
!>          LTB is INTEGER
!>          The size of the array TB. LTB >= 4*N, internally
!>          used to select NB such that LTB >= (3*NB+1)*N.
!>
!>          If LTB = -1, then a workspace query is assumed; the
!>          routine only calculates the optimal size of LTB,
!>          returns this value as the first entry of TB, and
!>          no error message related to LTB is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          On exit, it contains the details of the interchanges, i.e.,
!>          the row and column k of A were interchanged with the
!>          row and column IPIV(k).
!> \endverbatim
!>
!> \param[out] IPIV2
!> \verbatim
!>          IPIV2 is INTEGER array, dimension (N)
!>          On exit, it contains the details of the interchanges, i.e.,
!>          the row and column k of T were interchanged with the
!>          row and column IPIV(k).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL workspace of size LWORK
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The size of WORK. LWORK >= N, internally used to select NB
!>          such that LWORK >= N*NB.
!>
!>          If LWORK = -1, then a workspace query is assumed; the
!>          routine only calculates the optimal size of the WORK array,
!>          returns this value as the first entry of the WORK array, and
!>          no error message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          > 0:  if INFO = i, band LU factorization failed on i-th column
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
!> \ingroup realSYcomputational
!
!  =====================================================================
      SUBROUTINE SSYTRF_AA_2STAGE(Uplo,N,A,Lda,Tb,Ltb,Ipiv,Ipiv2,Work,  &
     &                            Lwork,Info)
!
!  -- LAPACK computational routine (version 3.8.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2017
!
      USE S_ILAENV
      USE S_LSAME
      USE S_SCOPY
      USE S_SGBTRF
      USE S_SGEMM
      USE S_SGETRF
      USE S_SLACPY
      USE S_SLASET
      USE S_SSWAP
      USE S_SSYGST
      USE S_STRSM
      USE S_XERBLA
      IMPLICIT NONE
!*--SSYTRF_AA_2STAGE182
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(INOUT) , DIMENSION(*) :: Tb
      INTEGER , INTENT(IN) :: Ltb
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      INTEGER , DIMENSION(*) :: Ipiv2
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , i1 , i2 , iinfo , j , jb , k , kb , ldtb , nb ,    &
     &           nt , td
      REAL :: piv
      LOGICAL :: tquery , upper , wquery
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!     .. Parameters ..
!
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
      wquery = (Lwork==-1)
      tquery = (Ltb==-1)
      IF ( .NOT.upper .AND. .NOT.LSAME(Uplo,'L') ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -4
      ELSEIF ( Ltb<4*N .AND. .NOT.tquery ) THEN
         Info = -6
      ELSEIF ( Lwork<N .AND. .NOT.wquery ) THEN
         Info = -10
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('SSYTRF_AA_2STAGE',-Info)
         RETURN
      ENDIF
!
!     Answer the query
!
      nb = ILAENV(1,'SSYTRF_AA_2STAGE',Uplo,N,-1,-1,-1)
      IF ( Info==0 ) THEN
         IF ( tquery ) Tb(1) = (3*nb+1)*N
         IF ( wquery ) Work(1) = N*nb
      ENDIF
      IF ( tquery .OR. wquery ) RETURN
!
!     Quick return
!
      IF ( N==0 ) RETURN
!
!     Determine the number of the block size
!
      ldtb = Ltb/N
      IF ( ldtb<3*nb+1 ) nb = (ldtb-1)/3
      IF ( Lwork<nb*N ) nb = Lwork/N
!
!     Determine the number of the block columns
!
      nt = (N+nb-1)/nb
      td = 2*nb
      kb = MIN(nb,N)
!
!     Initialize vectors/matrices
!
      DO j = 1 , kb
         Ipiv(j) = j
      ENDDO
!
!     Save NB
!
      Tb(1) = nb
!
      IF ( upper ) THEN
!
!        .....................................................
!        Factorize A as U**T*D*U using the upper triangle of A
!        .....................................................
!
         DO j = 0 , nt - 1
!
!           Generate Jth column of W and H
!
            kb = MIN(nb,N-j*nb)
            DO i = 1 , j - 1
               IF ( i==1 ) THEN
!                 H(I,J) = T(I,I)*U(I,J) + T(I+1,I)*U(I+1,J)
                  IF ( i==(j-1) ) THEN
                     jb = nb + kb
                  ELSE
                     jb = 2*nb
                  ENDIF
                  CALL SGEMM('NoTranspose','NoTranspose',nb,kb,jb,ONE,  &
     &                       Tb(td+1+(i*nb)*ldtb),ldtb-1,               &
     &                       A((i-1)*nb+1,j*nb+1),Lda,ZERO,Work(i*nb+1),&
     &                       N)
               ELSE
!                 H(I,J) = T(I,I-1)*U(I-1,J) + T(I,I)*U(I,J) + T(I,I+1)*U(I+1,J)
                  IF ( i==j-1 ) THEN
                     jb = 2*nb + kb
                  ELSE
                     jb = 3*nb
                  ENDIF
                  CALL SGEMM('NoTranspose','NoTranspose',nb,kb,jb,ONE,  &
     &                       Tb(td+nb+1+((i-1)*nb)*ldtb),ldtb-1,        &
     &                       A((i-2)*nb+1,j*nb+1),Lda,ZERO,Work(i*nb+1),&
     &                       N)
               ENDIF
            ENDDO
!
!           Compute T(J,J)
!
            CALL SLACPY('Upper',kb,kb,A(j*nb+1,j*nb+1),Lda,             &
     &                  Tb(td+1+(j*nb)*ldtb),ldtb-1)
            IF ( j>1 ) THEN
!              T(J,J) = U(1:J,J)'*H(1:J)
               CALL SGEMM('Transpose','NoTranspose',kb,kb,(j-1)*nb,-ONE,&
     &                    A(1,j*nb+1),Lda,Work(nb+1),N,ONE,             &
     &                    Tb(td+1+(j*nb)*ldtb),ldtb-1)
!              T(J,J) += U(J,J)'*T(J,J-1)*U(J-1,J)
               CALL SGEMM('Transpose','NoTranspose',kb,nb,kb,ONE,       &
     &                    A((j-1)*nb+1,j*nb+1),Lda,                     &
     &                    Tb(td+nb+1+((j-1)*nb)*ldtb),ldtb-1,ZERO,      &
     &                    Work(1),N)
               CALL SGEMM('NoTranspose','NoTranspose',kb,kb,nb,-ONE,    &
     &                    Work(1),N,A((j-2)*nb+1,j*nb+1),Lda,ONE,       &
     &                    Tb(td+1+(j*nb)*ldtb),ldtb-1)
            ENDIF
            IF ( j>0 ) CALL SSYGST(1,'Upper',kb,Tb(td+1+(j*nb)*ldtb),   &
     &                             ldtb-1,A((j-1)*nb+1,j*nb+1),Lda,     &
     &                             iinfo)
!
!           Expand T(J,J) into full format
!
            DO i = 1 , kb
               DO k = i + 1 , kb
                  Tb(td+(k-i)+1+(j*nb+i-1)*ldtb)                        &
     &               = Tb(td-(k-(i+1))+(j*nb+k-1)*ldtb)
               ENDDO
            ENDDO
!
            IF ( j<nt-1 ) THEN
               IF ( j>0 ) THEN
!
!                 Compute H(J,J)
!
                  IF ( j==1 ) THEN
                     CALL SGEMM('NoTranspose','NoTranspose',kb,kb,kb,   &
     &                          ONE,Tb(td+1+(j*nb)*ldtb),ldtb-1,        &
     &                          A((j-1)*nb+1,j*nb+1),Lda,ZERO,          &
     &                          Work(j*nb+1),N)
                  ELSE
                     CALL SGEMM('NoTranspose','NoTranspose',kb,kb,nb+kb,&
     &                          ONE,Tb(td+nb+1+((j-1)*nb)*ldtb),ldtb-1, &
     &                          A((j-2)*nb+1,j*nb+1),Lda,ZERO,          &
     &                          Work(j*nb+1),N)
                  ENDIF
!
!                 Update with the previous column
!
                  CALL SGEMM('Transpose','NoTranspose',nb,N-(j+1)*nb,   &
     &                       j*nb,-ONE,Work(nb+1),N,A(1,(j+1)*nb+1),Lda,&
     &                       ONE,A(j*nb+1,(j+1)*nb+1),Lda)
               ENDIF
!
!              Copy panel to workspace to call SGETRF
!
               DO k = 1 , nb
                  CALL SCOPY(N-(j+1)*nb,A(j*nb+k,(j+1)*nb+1),Lda,       &
     &                       Work(1+(k-1)*N),1)
               ENDDO
!
!              Factorize panel
!
               CALL SGETRF(N-(j+1)*nb,nb,Work,N,Ipiv((j+1)*nb+1),iinfo)
!               IF (IINFO.NE.0 .AND. INFO.EQ.0) THEN
!                  INFO = IINFO+(J+1)*NB
!               END IF
!
!              Copy panel back
!
               DO k = 1 , nb
                  CALL SCOPY(N-(j+1)*nb,Work(1+(k-1)*N),1,              &
     &                       A(j*nb+k,(j+1)*nb+1),Lda)
               ENDDO
!
!              Compute T(J+1, J), zero out for GEMM update
!
               kb = MIN(nb,N-(j+1)*nb)
               CALL SLASET('Full',kb,nb,ZERO,ZERO,                      &
     &                     Tb(td+nb+1+(j*nb)*ldtb),ldtb-1)
               CALL SLACPY('Upper',kb,nb,Work,N,Tb(td+nb+1+(j*nb)*ldtb),&
     &                     ldtb-1)
               IF ( j>0 ) CALL STRSM('R','U','N','U',kb,nb,ONE,         &
     &                               A((j-1)*nb+1,j*nb+1),Lda,          &
     &                               Tb(td+nb+1+(j*nb)*ldtb),ldtb-1)
!
!              Copy T(J,J+1) into T(J+1, J), both upper/lower for GEMM
!              updates
!
               DO k = 1 , nb
                  DO i = 1 , kb
                     Tb(td-nb+k-i+1+(j*nb+nb+i-1)*ldtb)                 &
     &                  = Tb(td+nb+i-k+1+(j*nb+k-1)*ldtb)
                  ENDDO
               ENDDO
               CALL SLASET('Lower',kb,nb,ZERO,ONE,A(j*nb+1,(j+1)*nb+1), &
     &                     Lda)
!
!              Apply pivots to trailing submatrix of A
!
               DO k = 1 , kb
!                 > Adjust ipiv
                  Ipiv((j+1)*nb+k) = Ipiv((j+1)*nb+k) + (j+1)*nb
!
                  i1 = (j+1)*nb + k
                  i2 = Ipiv((j+1)*nb+k)
                  IF ( i1/=i2 ) THEN
!                    > Apply pivots to previous columns of L
                     CALL SSWAP(k-1,A((j+1)*nb+1,i1),1,A((j+1)*nb+1,i2),&
     &                          1)
!                    > Swap A(I1+1:M, I1) with A(I2, I1+1:M)
                     IF ( i2>(i1+1) )                                   &
     &                    CALL SSWAP(i2-i1-1,A(i1,i1+1),Lda,A(i1+1,i2), &
     &                    1)
!                    > Swap A(I2+1:M, I1) with A(I2+1:M, I2)
                     IF ( i2<N )                                        &
     &                    CALL SSWAP(N-i2,A(i1,i2+1),Lda,A(i2,i2+1),Lda)
!                    > Swap A(I1, I1) with A(I2, I2)
                     piv = A(i1,i1)
                     A(i1,i1) = A(i2,i2)
                     A(i2,i2) = piv
!                    > Apply pivots to previous columns of L
                     IF ( j>0 ) CALL SSWAP(j*nb,A(1,i1),1,A(1,i2),1)
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ELSE
!
!        .....................................................
!        Factorize A as L*D*L**T using the lower triangle of A
!        .....................................................
!
         DO j = 0 , nt - 1
!
!           Generate Jth column of W and H
!
            kb = MIN(nb,N-j*nb)
            DO i = 1 , j - 1
               IF ( i==1 ) THEN
!                  H(I,J) = T(I,I)*L(J,I)' + T(I+1,I)'*L(J,I+1)'
                  IF ( i==(j-1) ) THEN
                     jb = nb + kb
                  ELSE
                     jb = 2*nb
                  ENDIF
                  CALL SGEMM('NoTranspose','Transpose',nb,kb,jb,ONE,    &
     &                       Tb(td+1+(i*nb)*ldtb),ldtb-1,               &
     &                       A(j*nb+1,(i-1)*nb+1),Lda,ZERO,Work(i*nb+1),&
     &                       N)
               ELSE
!                 H(I,J) = T(I,I-1)*L(J,I-1)' + T(I,I)*L(J,I)' + T(I,I+1)*L(J,I+1)'
                  IF ( i==j-1 ) THEN
                     jb = 2*nb + kb
                  ELSE
                     jb = 3*nb
                  ENDIF
                  CALL SGEMM('NoTranspose','Transpose',nb,kb,jb,ONE,    &
     &                       Tb(td+nb+1+((i-1)*nb)*ldtb),ldtb-1,        &
     &                       A(j*nb+1,(i-2)*nb+1),Lda,ZERO,Work(i*nb+1),&
     &                       N)
               ENDIF
            ENDDO
!
!           Compute T(J,J)
!
            CALL SLACPY('Lower',kb,kb,A(j*nb+1,j*nb+1),Lda,             &
     &                  Tb(td+1+(j*nb)*ldtb),ldtb-1)
            IF ( j>1 ) THEN
!              T(J,J) = L(J,1:J)*H(1:J)
               CALL SGEMM('NoTranspose','NoTranspose',kb,kb,(j-1)*nb,   &
     &                    -ONE,A(j*nb+1,1),Lda,Work(nb+1),N,ONE,        &
     &                    Tb(td+1+(j*nb)*ldtb),ldtb-1)
!              T(J,J) += L(J,J)*T(J,J-1)*L(J,J-1)'
               CALL SGEMM('NoTranspose','NoTranspose',kb,nb,kb,ONE,     &
     &                    A(j*nb+1,(j-1)*nb+1),Lda,                     &
     &                    Tb(td+nb+1+((j-1)*nb)*ldtb),ldtb-1,ZERO,      &
     &                    Work(1),N)
               CALL SGEMM('NoTranspose','Transpose',kb,kb,nb,-ONE,      &
     &                    Work(1),N,A(j*nb+1,(j-2)*nb+1),Lda,ONE,       &
     &                    Tb(td+1+(j*nb)*ldtb),ldtb-1)
            ENDIF
            IF ( j>0 ) CALL SSYGST(1,'Lower',kb,Tb(td+1+(j*nb)*ldtb),   &
     &                             ldtb-1,A(j*nb+1,(j-1)*nb+1),Lda,     &
     &                             iinfo)
!
!           Expand T(J,J) into full format
!
            DO i = 1 , kb
               DO k = i + 1 , kb
                  Tb(td-(k-(i+1))+(j*nb+k-1)*ldtb)                      &
     &               = Tb(td+(k-i)+1+(j*nb+i-1)*ldtb)
               ENDDO
            ENDDO
!
            IF ( j<nt-1 ) THEN
               IF ( j>0 ) THEN
!
!                 Compute H(J,J)
!
                  IF ( j==1 ) THEN
                     CALL SGEMM('NoTranspose','Transpose',kb,kb,kb,ONE, &
     &                          Tb(td+1+(j*nb)*ldtb),ldtb-1,            &
     &                          A(j*nb+1,(j-1)*nb+1),Lda,ZERO,          &
     &                          Work(j*nb+1),N)
                  ELSE
                     CALL SGEMM('NoTranspose','Transpose',kb,kb,nb+kb,  &
     &                          ONE,Tb(td+nb+1+((j-1)*nb)*ldtb),ldtb-1, &
     &                          A(j*nb+1,(j-2)*nb+1),Lda,ZERO,          &
     &                          Work(j*nb+1),N)
                  ENDIF
!
!                 Update with the previous column
!
                  CALL SGEMM('NoTranspose','NoTranspose',N-(j+1)*nb,nb, &
     &                       j*nb,-ONE,A((j+1)*nb+1,1),Lda,Work(nb+1),N,&
     &                       ONE,A((j+1)*nb+1,j*nb+1),Lda)
               ENDIF
!
!              Factorize panel
!
               CALL SGETRF(N-(j+1)*nb,nb,A((j+1)*nb+1,j*nb+1),Lda,      &
     &                     Ipiv((j+1)*nb+1),iinfo)
!               IF (IINFO.NE.0 .AND. INFO.EQ.0) THEN
!                  INFO = IINFO+(J+1)*NB
!               END IF
!
!              Compute T(J+1, J), zero out for GEMM update
!
               kb = MIN(nb,N-(j+1)*nb)
               CALL SLASET('Full',kb,nb,ZERO,ZERO,                      &
     &                     Tb(td+nb+1+(j*nb)*ldtb),ldtb-1)
               CALL SLACPY('Upper',kb,nb,A((j+1)*nb+1,j*nb+1),Lda,      &
     &                     Tb(td+nb+1+(j*nb)*ldtb),ldtb-1)
               IF ( j>0 ) CALL STRSM('R','L','T','U',kb,nb,ONE,         &
     &                               A(j*nb+1,(j-1)*nb+1),Lda,          &
     &                               Tb(td+nb+1+(j*nb)*ldtb),ldtb-1)
!
!              Copy T(J+1,J) into T(J, J+1), both upper/lower for GEMM
!              updates
!
               DO k = 1 , nb
                  DO i = 1 , kb
                     Tb(td-nb+k-i+1+(j*nb+nb+i-1)*ldtb)                 &
     &                  = Tb(td+nb+i-k+1+(j*nb+k-1)*ldtb)
                  ENDDO
               ENDDO
               CALL SLASET('Upper',kb,nb,ZERO,ONE,A((j+1)*nb+1,j*nb+1), &
     &                     Lda)
!
!              Apply pivots to trailing submatrix of A
!
               DO k = 1 , kb
!                 > Adjust ipiv
                  Ipiv((j+1)*nb+k) = Ipiv((j+1)*nb+k) + (j+1)*nb
!
                  i1 = (j+1)*nb + k
                  i2 = Ipiv((j+1)*nb+k)
                  IF ( i1/=i2 ) THEN
!                    > Apply pivots to previous columns of L
                     CALL SSWAP(k-1,A(i1,(j+1)*nb+1),Lda,               &
     &                          A(i2,(j+1)*nb+1),Lda)
!                    > Swap A(I1+1:M, I1) with A(I2, I1+1:M)
                     IF ( i2>(i1+1) )                                   &
     &                    CALL SSWAP(i2-i1-1,A(i1+1,i1),1,A(i2,i1+1),   &
     &                    Lda)
!                    > Swap A(I2+1:M, I1) with A(I2+1:M, I2)
                     IF ( i2<N )                                        &
     &                    CALL SSWAP(N-i2,A(i2+1,i1),1,A(i2+1,i2),1)
!                    > Swap A(I1, I1) with A(I2, I2)
                     piv = A(i1,i1)
                     A(i1,i1) = A(i2,i2)
                     A(i2,i2) = piv
!                    > Apply pivots to previous columns of L
                     IF ( j>0 ) CALL SSWAP(j*nb,A(i1,1),Lda,A(i2,1),Lda)
                  ENDIF
               ENDDO
!
!              Apply pivots to previous columns of L
!
!               CALL SLASWP( J*NB, A( 1, 1 ), LDA,
!     $                     (J+1)*NB+1, (J+1)*NB+KB, IPIV, 1 )
            ENDIF
         ENDDO
      ENDIF
!
!     Factor the band matrix
      CALL SGBTRF(N,N,nb,nb,Tb,ldtb,Ipiv2,Info)
!
!
!     End of SSYTRF_AA_2STAGE
!
      END SUBROUTINE SSYTRF_AA_2STAGE
