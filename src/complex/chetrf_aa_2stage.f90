!*==chetrf_aa_2stage.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b CHETRF_AA_2STAGE
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CHETRF_AA_2STAGE + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chetrf_aa_2stage.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chetrf_aa_2stage.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chetrf_aa_2stage.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!      SUBROUTINE CHETRF_AA_2STAGE( UPLO, N, A, LDA, TB, LTB, IPIV,
!                                   IPIV2, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            N, LDA, LTB, LWORK, INFO
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * ), IPIV2( * )
!       COMPLEX            A( LDA, * ), TB( * ), WORK( * )
!       ..
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CHETRF_AA_2STAGE computes the factorization of a real hermitian matrix A
!> using the Aasen's algorithm.  The form of the factorization is
!>
!>    A = U**T*T*U  or  A = L*T*L**T
!>
!> where U (or L) is a product of permutation and unit upper (lower)
!> triangular matrices, and T is a hermitian band matrix with the
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
!>          A is COMPLEX array, dimension (LDA,N)
!>          On entry, the hermitian matrix A.  If UPLO = 'U', the leading
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
!>          TB is COMPLEX array, dimension (LTB)
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
!>          WORK is COMPLEX workspace of size LWORK
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
!> \ingroup complexSYcomputational
!
!  =====================================================================
      SUBROUTINE CHETRF_AA_2STAGE(Uplo,N,A,Lda,Tb,Ltb,Ipiv,Ipiv2,Work,  &
     &                            Lwork,Info)
!
!  -- LAPACK computational routine (version 3.8.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2017
!
      IMPLICIT NONE
!*--CHETRF_AA_2STAGE170
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER N , Lda , Ltb , Lwork , Info
!     ..
!     .. Array Arguments ..
      INTEGER Ipiv(*) , Ipiv2(*)
      COMPLEX A(Lda,*) , Tb(*) , Work(*)
!     ..
!
!  =====================================================================
!     .. Parameters ..
      COMPLEX ZERO , ONE
      PARAMETER (ZERO=(0.0E+0,0.0E+0),ONE=(1.0E+0,0.0E+0))
!
!     .. Local Scalars ..
      LOGICAL upper , tquery , wquery
      INTEGER i , j , k , i1 , i2 , td
      INTEGER ldtb , nb , kb , jb , nt , iinfo
      COMPLEX piv
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER ILAENV
      EXTERNAL LSAME , ILAENV
 
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA , CCOPY , CLACGV , CLACPY , CLASET , CGBTRF ,     &
     &         CGEMM , CGETRF , CHEGST , CSWAP , CTRSM
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC CONJG , MIN , MAX
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
         CALL XERBLA('CHETRF_AA_2STAGE',-Info)
         RETURN
      ENDIF
!
!     Answer the query
!
      nb = ILAENV(1,'CHETRF_AA_2STAGE',Uplo,N,-1,-1,-1)
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
!                  H(I,J) = T(I,I)*U(I,J) + T(I+1,I)*U(I+1,J)
                  IF ( i==(j-1) ) THEN
                     jb = nb + kb
                  ELSE
                     jb = 2*nb
                  ENDIF
                  CALL CGEMM('NoTranspose','NoTranspose',nb,kb,jb,ONE,  &
     &                       Tb(td+1+(i*nb)*ldtb),ldtb-1,               &
     &                       A((i-1)*nb+1,j*nb+1),Lda,ZERO,Work(i*nb+1),&
     &                       N)
               ELSE
!                 H(I,J) = T(I,I-1)*U(I-1,J) + T(I,I)*U(I,J) + T(I,I+1)*U(I+1,J)
                  IF ( i==(j-1) ) THEN
                     jb = 2*nb + kb
                  ELSE
                     jb = 3*nb
                  ENDIF
                  CALL CGEMM('NoTranspose','NoTranspose',nb,kb,jb,ONE,  &
     &                       Tb(td+nb+1+((i-1)*nb)*ldtb),ldtb-1,        &
     &                       A((i-2)*nb+1,j*nb+1),Lda,ZERO,Work(i*nb+1),&
     &                       N)
               ENDIF
            ENDDO
!
!           Compute T(J,J)
!
            CALL CLACPY('Upper',kb,kb,A(j*nb+1,j*nb+1),Lda,             &
     &                  Tb(td+1+(j*nb)*ldtb),ldtb-1)
            IF ( j>1 ) THEN
!              T(J,J) = U(1:J,J)'*H(1:J)
               CALL CGEMM('Conjugate transpose','NoTranspose',kb,kb,    &
     &                    (j-1)*nb,-ONE,A(1,j*nb+1),Lda,Work(nb+1),N,   &
     &                    ONE,Tb(td+1+(j*nb)*ldtb),ldtb-1)
!              T(J,J) += U(J,J)'*T(J,J-1)*U(J-1,J)
               CALL CGEMM('Conjugate transpose','NoTranspose',kb,nb,kb, &
     &                    ONE,A((j-1)*nb+1,j*nb+1),Lda,                 &
     &                    Tb(td+nb+1+((j-1)*nb)*ldtb),ldtb-1,ZERO,      &
     &                    Work(1),N)
               CALL CGEMM('NoTranspose','NoTranspose',kb,kb,nb,-ONE,    &
     &                    Work(1),N,A((j-2)*nb+1,j*nb+1),Lda,ONE,       &
     &                    Tb(td+1+(j*nb)*ldtb),ldtb-1)
            ENDIF
            IF ( j>0 ) CALL CHEGST(1,'Upper',kb,Tb(td+1+(j*nb)*ldtb),   &
     &                             ldtb-1,A((j-1)*nb+1,j*nb+1),Lda,     &
     &                             iinfo)
!
!           Expand T(J,J) into full format
!
            DO i = 1 , kb
               Tb(td+1+(j*nb+i-1)*ldtb) = REAL(Tb(td+1+(j*nb+i-1)*ldtb))
               DO k = i + 1 , kb
                  Tb(td+(k-i)+1+(j*nb+i-1)*ldtb)                        &
     &               = CONJG(Tb(td-(k-(i+1))+(j*nb+k-1)*ldtb))
               ENDDO
            ENDDO
!
            IF ( j<nt-1 ) THEN
               IF ( j>0 ) THEN
!
!                 Compute H(J,J)
!
                  IF ( j==1 ) THEN
                     CALL CGEMM('NoTranspose','NoTranspose',kb,kb,kb,   &
     &                          ONE,Tb(td+1+(j*nb)*ldtb),ldtb-1,        &
     &                          A((j-1)*nb+1,j*nb+1),Lda,ZERO,          &
     &                          Work(j*nb+1),N)
                  ELSE
                     CALL CGEMM('NoTranspose','NoTranspose',kb,kb,nb+kb,&
     &                          ONE,Tb(td+nb+1+((j-1)*nb)*ldtb),ldtb-1, &
     &                          A((j-2)*nb+1,j*nb+1),Lda,ZERO,          &
     &                          Work(j*nb+1),N)
                  ENDIF
!
!                 Update with the previous column
!
                  CALL CGEMM('Conjugate transpose','NoTranspose',nb,    &
     &                       N-(j+1)*nb,j*nb,-ONE,Work(nb+1),N,         &
     &                       A(1,(j+1)*nb+1),Lda,ONE,                   &
     &                       A(j*nb+1,(j+1)*nb+1),Lda)
               ENDIF
!
!              Copy panel to workspace to call CGETRF
!
               DO k = 1 , nb
                  CALL CCOPY(N-(j+1)*nb,A(j*nb+k,(j+1)*nb+1),Lda,       &
     &                       Work(1+(k-1)*N),1)
               ENDDO
!
!              Factorize panel
!
               CALL CGETRF(N-(j+1)*nb,nb,Work,N,Ipiv((j+1)*nb+1),iinfo)
!               IF (IINFO.NE.0 .AND. INFO.EQ.0) THEN
!                  INFO = IINFO+(J+1)*NB
!               END IF
!
!              Copy panel back
!
               DO k = 1 , nb
!
!                  Copy only L-factor
!
                  CALL CCOPY(N-k-(j+1)*nb,Work(k+1+(k-1)*N),1,          &
     &                       A(j*nb+k,(j+1)*nb+k+1),Lda)
!
!                  Transpose U-factor to be copied back into T(J+1, J)
!
                  CALL CLACGV(k,Work(1+(k-1)*N),1)
               ENDDO
!
!              Compute T(J+1, J), zero out for GEMM update
!
               kb = MIN(nb,N-(j+1)*nb)
               CALL CLASET('Full',kb,nb,ZERO,ZERO,                      &
     &                     Tb(td+nb+1+(j*nb)*ldtb),ldtb-1)
               CALL CLACPY('Upper',kb,nb,Work,N,Tb(td+nb+1+(j*nb)*ldtb),&
     &                     ldtb-1)
               IF ( j>0 ) CALL CTRSM('R','U','N','U',kb,nb,ONE,         &
     &                               A((j-1)*nb+1,j*nb+1),Lda,          &
     &                               Tb(td+nb+1+(j*nb)*ldtb),ldtb-1)
!
!              Copy T(J,J+1) into T(J+1, J), both upper/lower for GEMM
!              updates
!
               DO k = 1 , nb
                  DO i = 1 , kb
                     Tb(td-nb+k-i+1+(j*nb+nb+i-1)*ldtb)                 &
     &                  = CONJG(Tb(td+nb+i-k+1+(j*nb+k-1)*ldtb))
                  ENDDO
               ENDDO
               CALL CLASET('Lower',kb,nb,ZERO,ONE,A(j*nb+1,(j+1)*nb+1), &
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
                     CALL CSWAP(k-1,A((j+1)*nb+1,i1),1,A((j+1)*nb+1,i2),&
     &                          1)
!                    > Swap A(I1+1:M, I1) with A(I2, I1+1:M)
                     IF ( i2>(i1+1) ) THEN
                        CALL CSWAP(i2-i1-1,A(i1,i1+1),Lda,A(i1+1,i2),1)
                        CALL CLACGV(i2-i1-1,A(i1+1,i2),1)
                     ENDIF
                     CALL CLACGV(i2-i1,A(i1,i1+1),Lda)
!                    > Swap A(I2+1:M, I1) with A(I2+1:M, I2)
                     IF ( i2<N )                                        &
     &                    CALL CSWAP(N-i2,A(i1,i2+1),Lda,A(i2,i2+1),Lda)
!                    > Swap A(I1, I1) with A(I2, I2)
                     piv = A(i1,i1)
                     A(i1,i1) = A(i2,i2)
                     A(i2,i2) = piv
!                    > Apply pivots to previous columns of L
                     IF ( j>0 ) CALL CSWAP(j*nb,A(1,i1),1,A(1,i2),1)
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
                  CALL CGEMM('NoTranspose','Conjugate transpose',nb,kb, &
     &                       jb,ONE,Tb(td+1+(i*nb)*ldtb),ldtb-1,        &
     &                       A(j*nb+1,(i-1)*nb+1),Lda,ZERO,Work(i*nb+1),&
     &                       N)
               ELSE
!                 H(I,J) = T(I,I-1)*L(J,I-1)' + T(I,I)*L(J,I)' + T(I,I+1)*L(J,I+1)'
                  IF ( i==(j-1) ) THEN
                     jb = 2*nb + kb
                  ELSE
                     jb = 3*nb
                  ENDIF
                  CALL CGEMM('NoTranspose','Conjugate transpose',nb,kb, &
     &                       jb,ONE,Tb(td+nb+1+((i-1)*nb)*ldtb),ldtb-1, &
     &                       A(j*nb+1,(i-2)*nb+1),Lda,ZERO,Work(i*nb+1),&
     &                       N)
               ENDIF
            ENDDO
!
!           Compute T(J,J)
!
            CALL CLACPY('Lower',kb,kb,A(j*nb+1,j*nb+1),Lda,             &
     &                  Tb(td+1+(j*nb)*ldtb),ldtb-1)
            IF ( j>1 ) THEN
!              T(J,J) = L(J,1:J)*H(1:J)
               CALL CGEMM('NoTranspose','NoTranspose',kb,kb,(j-1)*nb,   &
     &                    -ONE,A(j*nb+1,1),Lda,Work(nb+1),N,ONE,        &
     &                    Tb(td+1+(j*nb)*ldtb),ldtb-1)
!              T(J,J) += L(J,J)*T(J,J-1)*L(J,J-1)'
               CALL CGEMM('NoTranspose','NoTranspose',kb,nb,kb,ONE,     &
     &                    A(j*nb+1,(j-1)*nb+1),Lda,                     &
     &                    Tb(td+nb+1+((j-1)*nb)*ldtb),ldtb-1,ZERO,      &
     &                    Work(1),N)
               CALL CGEMM('NoTranspose','Conjugate transpose',kb,kb,nb, &
     &                    -ONE,Work(1),N,A(j*nb+1,(j-2)*nb+1),Lda,ONE,  &
     &                    Tb(td+1+(j*nb)*ldtb),ldtb-1)
            ENDIF
            IF ( j>0 ) CALL CHEGST(1,'Lower',kb,Tb(td+1+(j*nb)*ldtb),   &
     &                             ldtb-1,A(j*nb+1,(j-1)*nb+1),Lda,     &
     &                             iinfo)
!
!           Expand T(J,J) into full format
!
            DO i = 1 , kb
               Tb(td+1+(j*nb+i-1)*ldtb) = REAL(Tb(td+1+(j*nb+i-1)*ldtb))
               DO k = i + 1 , kb
                  Tb(td-(k-(i+1))+(j*nb+k-1)*ldtb)                      &
     &               = CONJG(Tb(td+(k-i)+1+(j*nb+i-1)*ldtb))
               ENDDO
            ENDDO
!
            IF ( j<nt-1 ) THEN
               IF ( j>0 ) THEN
!
!                 Compute H(J,J)
!
                  IF ( j==1 ) THEN
                     CALL CGEMM('NoTranspose','Conjugate transpose',kb, &
     &                          kb,kb,ONE,Tb(td+1+(j*nb)*ldtb),ldtb-1,  &
     &                          A(j*nb+1,(j-1)*nb+1),Lda,ZERO,          &
     &                          Work(j*nb+1),N)
                  ELSE
                     CALL CGEMM('NoTranspose','Conjugate transpose',kb, &
     &                          kb,nb+kb,ONE,Tb(td+nb+1+((j-1)*nb)*ldtb)&
     &                          ,ldtb-1,A(j*nb+1,(j-2)*nb+1),Lda,ZERO,  &
     &                          Work(j*nb+1),N)
                  ENDIF
!
!                 Update with the previous column
!
                  CALL CGEMM('NoTranspose','NoTranspose',N-(j+1)*nb,nb, &
     &                       j*nb,-ONE,A((j+1)*nb+1,1),Lda,Work(nb+1),N,&
     &                       ONE,A((j+1)*nb+1,j*nb+1),Lda)
               ENDIF
!
!              Factorize panel
!
               CALL CGETRF(N-(j+1)*nb,nb,A((j+1)*nb+1,j*nb+1),Lda,      &
     &                     Ipiv((j+1)*nb+1),iinfo)
!               IF (IINFO.NE.0 .AND. INFO.EQ.0) THEN
!                  INFO = IINFO+(J+1)*NB
!               END IF
!
!              Compute T(J+1, J), zero out for GEMM update
!
               kb = MIN(nb,N-(j+1)*nb)
               CALL CLASET('Full',kb,nb,ZERO,ZERO,                      &
     &                     Tb(td+nb+1+(j*nb)*ldtb),ldtb-1)
               CALL CLACPY('Upper',kb,nb,A((j+1)*nb+1,j*nb+1),Lda,      &
     &                     Tb(td+nb+1+(j*nb)*ldtb),ldtb-1)
               IF ( j>0 ) CALL CTRSM('R','L','C','U',kb,nb,ONE,         &
     &                               A(j*nb+1,(j-1)*nb+1),Lda,          &
     &                               Tb(td+nb+1+(j*nb)*ldtb),ldtb-1)
!
!              Copy T(J+1,J) into T(J, J+1), both upper/lower for GEMM
!              updates
!
               DO k = 1 , nb
                  DO i = 1 , kb
                     Tb(td-nb+k-i+1+(j*nb+nb+i-1)*ldtb)                 &
     &                  = CONJG(Tb(td+nb+i-k+1+(j*nb+k-1)*ldtb))
                  ENDDO
               ENDDO
               CALL CLASET('Upper',kb,nb,ZERO,ONE,A((j+1)*nb+1,j*nb+1), &
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
                     CALL CSWAP(k-1,A(i1,(j+1)*nb+1),Lda,               &
     &                          A(i2,(j+1)*nb+1),Lda)
!                    > Swap A(I1+1:M, I1) with A(I2, I1+1:M)
                     IF ( i2>(i1+1) ) THEN
                        CALL CSWAP(i2-i1-1,A(i1+1,i1),1,A(i2,i1+1),Lda)
                        CALL CLACGV(i2-i1-1,A(i2,i1+1),Lda)
                     ENDIF
                     CALL CLACGV(i2-i1,A(i1+1,i1),1)
!                    > Swap A(I2+1:M, I1) with A(I2+1:M, I2)
                     IF ( i2<N )                                        &
     &                    CALL CSWAP(N-i2,A(i2+1,i1),1,A(i2+1,i2),1)
!                    > Swap A(I1, I1) with A(I2, I2)
                     piv = A(i1,i1)
                     A(i1,i1) = A(i2,i2)
                     A(i2,i2) = piv
!                    > Apply pivots to previous columns of L
                     IF ( j>0 ) CALL CSWAP(j*nb,A(i1,1),Lda,A(i2,1),Lda)
                  ENDIF
               ENDDO
!
!              Apply pivots to previous columns of L
!
!               CALL CLASWP( J*NB, A( 1, 1 ), LDA,
!     $                     (J+1)*NB+1, (J+1)*NB+KB, IPIV, 1 )
            ENDIF
         ENDDO
      ENDIF
!
!     Factor the band matrix
      CALL CGBTRF(N,N,nb,nb,Tb,ldtb,Ipiv2,Info)
!
!
!     End of CHETRF_AA_2STAGE
!
      END SUBROUTINE CHETRF_AA_2STAGE
