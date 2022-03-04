!*==ssytrf_aa.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SSYTRF_AA
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SSYTRF_AA + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssytrf_aa.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssytrf_aa.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssytrf_aa.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SSYTRF_AA( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            N, LDA, LWORK, INFO
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       REAL   A( LDA, * ), WORK( * )
!       ..
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SSYTRF_AA computes the factorization of a real symmetric matrix A
!> using the Aasen's algorithm.  The form of the factorization is
!>
!>    A = U**T*T*U  or  A = L*T*L**T
!>
!> where U (or L) is a product of permutation and unit upper (lower)
!> triangular matrices, and T is a symmetric tridiagonal matrix.
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
!>          On exit, the tridiagonal matrix is stored in the diagonals
!>          and the subdiagonals of A just below (or above) the diagonals,
!>          and L is stored below (or above) the subdiaonals, when UPLO
!>          is 'L' (or 'U').
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
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
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The length of WORK.  LWORK >= MAX(1,2*N). For optimum performance
!>          LWORK >= N*(1+NB), where NB is the optimal blocksize.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
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
      SUBROUTINE SSYTRF_AA(Uplo,N,A,Lda,Ipiv,Work,Lwork,Info)
!
!  -- LAPACK computational routine (version 3.8.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2017
!
      USE S_ILAENV
      USE S_LSAME
      USE S_SCOPY
      USE S_SGEMM
      USE S_SGEMV
      USE S_SLASYF_AA
      USE S_SSCAL
      USE S_SSWAP
      USE S_XERBLA
      IMPLICIT NONE
!*--SSYTRF_AA151
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL :: alpha
      INTEGER :: j , j1 , j2 , j3 , jb , k1 , k2 , lwkopt , mj , nb , nj
      LOGICAL :: lquery , upper
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
!     Determine the block size
!
      nb = ILAENV(1,'SSYTRF_AA',Uplo,N,-1,-1,-1)
!
!     Test the input parameters.
!
      Info = 0
      upper = LSAME(Uplo,'U')
      lquery = (Lwork==-1)
      IF ( .NOT.upper .AND. .NOT.LSAME(Uplo,'L') ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -4
      ELSEIF ( Lwork<MAX(1,2*N) .AND. .NOT.lquery ) THEN
         Info = -7
      ENDIF
!
      IF ( Info==0 ) THEN
         lwkopt = (nb+1)*N
         Work(1) = lwkopt
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('SSYTRF_AA',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return
!
      IF ( N==0 ) RETURN
      Ipiv(1) = 1
      IF ( N==1 ) RETURN
!
!     Adjust block size based on the workspace size
!
      IF ( Lwork<((1+nb)*N) ) nb = (Lwork-N)/N
!
      IF ( upper ) THEN
!
!        .....................................................
!        Factorize A as U**T*D*U using the upper triangle of A
!        .....................................................
!
!        Copy first row A(1, 1:N) into H(1:n) (stored in WORK(1:N))
!
         CALL SCOPY(N,A(1,1),Lda,Work(1),1)
!
!        J is the main loop index, increasing from 1 to N in steps of
!        JB, where JB is the number of columns factorized by SLASYF;
!        JB is either NB, or N-J+1 for the last block
!
         j = 0
         DO WHILE ( j<N )
!
!        each step of the main loop
!         J is the last column of the previous panel
!         J1 is the first column of the current panel
!         K1 identifies if the previous column of the panel has been
!          explicitly stored, e.g., K1=1 for the first panel, and
!          K1=0 for the rest
!
            j1 = j + 1
            jb = MIN(N-j1+1,nb)
            k1 = MAX(1,j) - j
!
!        Panel factorization
!
            CALL SLASYF_AA(Uplo,2-k1,N-j,jb,A(MAX(1,j),j+1),Lda,        &
     &                     Ipiv(j+1),Work,N,Work(N*nb+1))
!
!        Adjust IPIV and apply it back (J-th step picks (J+1)-th pivot)
!
            DO j2 = j + 2 , MIN(N,j+jb+1)
               Ipiv(j2) = Ipiv(j2) + j
               IF ( (j2/=Ipiv(j2)) .AND. ((j1-k1)>2) )                  &
     &              CALL SSWAP(j1-k1-2,A(1,j2),1,A(1,Ipiv(j2)),1)
            ENDDO
            j = j + jb
!
!        Trailing submatrix update, where
!         the row A(J1-1, J2-1:N) stores U(J1, J2+1:N) and
!         WORK stores the current block of the auxiriarly matrix H
!
            IF ( j<N ) THEN
!
!           If first panel and JB=1 (NB=1), then nothing to do
!
               IF ( j1>1 .OR. jb>1 ) THEN
!
!              Merge rank-1 update with BLAS-3 update
!
                  alpha = A(j,j+1)
                  A(j,j+1) = ONE
                  CALL SCOPY(N-j,A(j-1,j+1),Lda,Work((j+1-j1+1)+jb*N),1)
                  CALL SSCAL(N-j,alpha,Work((j+1-j1+1)+jb*N),1)
!
!              K1 identifies if the previous column of the panel has been
!               explicitly stored, e.g., K1=1 and K2= 0 for the first panel,
!               while K1=0 and K2=1 for the rest
!
                  IF ( j1>1 ) THEN
!
!                 Not first panel
!
                     k2 = 1
                  ELSE
!
!                 First panel
!
                     k2 = 0
!
!                 First update skips the first column
!
                     jb = jb - 1
                  ENDIF
!
                  DO j2 = j + 1 , N , nb
                     nj = MIN(nb,N-j2+1)
!
!                 Update (J2, J2) diagonal block with SGEMV
!
                     j3 = j2
                     DO mj = nj - 1 , 1 , -1
                        CALL SGEMV('No transpose',mj,jb+1,-ONE,         &
     &                             Work(j3-j1+1+k1*N),N,A(j1-k2,j3),1,  &
     &                             ONE,A(j3,j3),Lda)
                        j3 = j3 + 1
                     ENDDO
!
!                 Update off-diagonal block of J2-th block row with SGEMM
!
                     CALL SGEMM('Transpose','Transpose',nj,N-j3+1,jb+1, &
     &                          -ONE,A(j1-k2,j2),Lda,Work(j3-j1+1+k1*N),&
     &                          N,ONE,A(j2,j3),Lda)
                  ENDDO
!
!              Recover T( J, J+1 )
!
                  A(j,j+1) = alpha
               ENDIF
!
!           WORK(J+1, 1) stores H(J+1, 1)
!
               CALL SCOPY(N-j,A(j+1,j+1),Lda,Work(1),1)
            ENDIF
         ENDDO
      ELSE
!
!        .....................................................
!        Factorize A as L*D*L**T using the lower triangle of A
!        .....................................................
!
!        copy first column A(1:N, 1) into H(1:N, 1)
!         (stored in WORK(1:N))
!
         CALL SCOPY(N,A(1,1),1,Work(1),1)
!
!        J is the main loop index, increasing from 1 to N in steps of
!        JB, where JB is the number of columns factorized by SLASYF;
!        JB is either NB, or N-J+1 for the last block
!
         j = 0
         DO WHILE ( j<N )
!
!        each step of the main loop
!         J is the last column of the previous panel
!         J1 is the first column of the current panel
!         K1 identifies if the previous column of the panel has been
!          explicitly stored, e.g., K1=1 for the first panel, and
!          K1=0 for the rest
!
            j1 = j + 1
            jb = MIN(N-j1+1,nb)
            k1 = MAX(1,j) - j
!
!        Panel factorization
!
            CALL SLASYF_AA(Uplo,2-k1,N-j,jb,A(j+1,MAX(1,j)),Lda,        &
     &                     Ipiv(j+1),Work,N,Work(N*nb+1))
!
!        Adjust IPIV and apply it back (J-th step picks (J+1)-th pivot)
!
            DO j2 = j + 2 , MIN(N,j+jb+1)
               Ipiv(j2) = Ipiv(j2) + j
               IF ( (j2/=Ipiv(j2)) .AND. ((j1-k1)>2) )                  &
     &              CALL SSWAP(j1-k1-2,A(j2,1),Lda,A(Ipiv(j2),1),Lda)
            ENDDO
            j = j + jb
!
!        Trailing submatrix update, where
!          A(J2+1, J1-1) stores L(J2+1, J1) and
!          WORK(J2+1, 1) stores H(J2+1, 1)
!
            IF ( j<N ) THEN
!
!           if first panel and JB=1 (NB=1), then nothing to do
!
               IF ( j1>1 .OR. jb>1 ) THEN
!
!              Merge rank-1 update with BLAS-3 update
!
                  alpha = A(j+1,j)
                  A(j+1,j) = ONE
                  CALL SCOPY(N-j,A(j+1,j-1),1,Work((j+1-j1+1)+jb*N),1)
                  CALL SSCAL(N-j,alpha,Work((j+1-j1+1)+jb*N),1)
!
!              K1 identifies if the previous column of the panel has been
!               explicitly stored, e.g., K1=1 and K2= 0 for the first panel,
!               while K1=0 and K2=1 for the rest
!
                  IF ( j1>1 ) THEN
!
!                 Not first panel
!
                     k2 = 1
                  ELSE
!
!                 First panel
!
                     k2 = 0
!
!                 First update skips the first column
!
                     jb = jb - 1
                  ENDIF
!
                  DO j2 = j + 1 , N , nb
                     nj = MIN(nb,N-j2+1)
!
!                 Update (J2, J2) diagonal block with SGEMV
!
                     j3 = j2
                     DO mj = nj - 1 , 1 , -1
                        CALL SGEMV('No transpose',mj,jb+1,-ONE,         &
     &                             Work(j3-j1+1+k1*N),N,A(j3,j1-k2),Lda,&
     &                             ONE,A(j3,j3),1)
                        j3 = j3 + 1
                     ENDDO
!
!                 Update off-diagonal block in J2-th block column with SGEMM
!
                     CALL SGEMM('No transpose','Transpose',N-j3+1,nj,   &
     &                          jb+1,-ONE,Work(j3-j1+1+k1*N),N,         &
     &                          A(j2,j1-k2),Lda,ONE,A(j3,j2),Lda)
                  ENDDO
!
!              Recover T( J+1, J )
!
                  A(j+1,j) = alpha
               ENDIF
!
!           WORK(J+1, 1) stores H(J+1, 1)
!
               CALL SCOPY(N-j,A(j+1,j+1),1,Work(1),1)
            ENDIF
         ENDDO
      ENDIF
!
!
!     End of SSYTRF_AA
!
      END SUBROUTINE SSYTRF_AA
