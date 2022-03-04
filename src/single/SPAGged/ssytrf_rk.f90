!*==ssytrf_rk.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SSYTRF_RK computes the factorization of a real symmetric indefinite matrix using the bounded Bunch-Kaufman (rook) diagonal pivoting method (BLAS3 blocked algorithm).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SSYTRF_RK + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssytrf_rk.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssytrf_rk.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssytrf_rk.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SSYTRF_RK( UPLO, N, A, LDA, E, IPIV, WORK, LWORK,
!                             INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDA, LWORK, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       REAL               A( LDA, * ), E ( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!> SSYTRF_RK computes the factorization of a real symmetric matrix A
!> using the bounded Bunch-Kaufman (rook) diagonal pivoting method:
!>
!>    A = P*U*D*(U**T)*(P**T) or A = P*L*D*(L**T)*(P**T),
!>
!> where U (or L) is unit upper (or lower) triangular matrix,
!> U**T (or L**T) is the transpose of U (or L), P is a permutation
!> matrix, P**T is the transpose of P, and D is symmetric and block
!> diagonal with 1-by-1 and 2-by-2 diagonal blocks.
!>
!> This is the blocked version of the algorithm, calling Level 3 BLAS.
!> For more information see Further Details section.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the upper or lower triangular part of the
!>          symmetric matrix A is stored:
!>          = 'U':  Upper triangular
!>          = 'L':  Lower triangular
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
!>          On entry, the symmetric matrix A.
!>            If UPLO = 'U': the leading N-by-N upper triangular part
!>            of A contains the upper triangular part of the matrix A,
!>            and the strictly lower triangular part of A is not
!>            referenced.
!>
!>            If UPLO = 'L': the leading N-by-N lower triangular part
!>            of A contains the lower triangular part of the matrix A,
!>            and the strictly upper triangular part of A is not
!>            referenced.
!>
!>          On exit, contains:
!>            a) ONLY diagonal elements of the symmetric block diagonal
!>               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k);
!>               (superdiagonal (or subdiagonal) elements of D
!>                are stored on exit in array E), and
!>            b) If UPLO = 'U': factor U in the superdiagonal part of A.
!>               If UPLO = 'L': factor L in the subdiagonal part of A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] E
!> \verbatim
!>          E is REAL array, dimension (N)
!>          On exit, contains the superdiagonal (or subdiagonal)
!>          elements of the symmetric block diagonal matrix D
!>          with 1-by-1 or 2-by-2 diagonal blocks, where
!>          If UPLO = 'U': E(i) = D(i-1,i), i=2:N, E(1) is set to 0;
!>          If UPLO = 'L': E(i) = D(i+1,i), i=1:N-1, E(N) is set to 0.
!>
!>          NOTE: For 1-by-1 diagonal block D(k), where
!>          1 <= k <= N, the element E(k) is set to 0 in both
!>          UPLO = 'U' or UPLO = 'L' cases.
!> \endverbatim
!>
!> \param[out] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          IPIV describes the permutation matrix P in the factorization
!>          of matrix A as follows. The absolute value of IPIV(k)
!>          represents the index of row and column that were
!>          interchanged with the k-th row and column. The value of UPLO
!>          describes the order in which the interchanges were applied.
!>          Also, the sign of IPIV represents the block structure of
!>          the symmetric block diagonal matrix D with 1-by-1 or 2-by-2
!>          diagonal blocks which correspond to 1 or 2 interchanges
!>          at each factorization step. For more info see Further
!>          Details section.
!>
!>          If UPLO = 'U',
!>          ( in factorization order, k decreases from N to 1 ):
!>            a) A single positive entry IPIV(k) > 0 means:
!>               D(k,k) is a 1-by-1 diagonal block.
!>               If IPIV(k) != k, rows and columns k and IPIV(k) were
!>               interchanged in the matrix A(1:N,1:N);
!>               If IPIV(k) = k, no interchange occurred.
!>
!>            b) A pair of consecutive negative entries
!>               IPIV(k) < 0 and IPIV(k-1) < 0 means:
!>               D(k-1:k,k-1:k) is a 2-by-2 diagonal block.
!>               (NOTE: negative entries in IPIV appear ONLY in pairs).
!>               1) If -IPIV(k) != k, rows and columns
!>                  k and -IPIV(k) were interchanged
!>                  in the matrix A(1:N,1:N).
!>                  If -IPIV(k) = k, no interchange occurred.
!>               2) If -IPIV(k-1) != k-1, rows and columns
!>                  k-1 and -IPIV(k-1) were interchanged
!>                  in the matrix A(1:N,1:N).
!>                  If -IPIV(k-1) = k-1, no interchange occurred.
!>
!>            c) In both cases a) and b), always ABS( IPIV(k) ) <= k.
!>
!>            d) NOTE: Any entry IPIV(k) is always NONZERO on output.
!>
!>          If UPLO = 'L',
!>          ( in factorization order, k increases from 1 to N ):
!>            a) A single positive entry IPIV(k) > 0 means:
!>               D(k,k) is a 1-by-1 diagonal block.
!>               If IPIV(k) != k, rows and columns k and IPIV(k) were
!>               interchanged in the matrix A(1:N,1:N).
!>               If IPIV(k) = k, no interchange occurred.
!>
!>            b) A pair of consecutive negative entries
!>               IPIV(k) < 0 and IPIV(k+1) < 0 means:
!>               D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
!>               (NOTE: negative entries in IPIV appear ONLY in pairs).
!>               1) If -IPIV(k) != k, rows and columns
!>                  k and -IPIV(k) were interchanged
!>                  in the matrix A(1:N,1:N).
!>                  If -IPIV(k) = k, no interchange occurred.
!>               2) If -IPIV(k+1) != k+1, rows and columns
!>                  k-1 and -IPIV(k-1) were interchanged
!>                  in the matrix A(1:N,1:N).
!>                  If -IPIV(k+1) = k+1, no interchange occurred.
!>
!>            c) In both cases a) and b), always ABS( IPIV(k) ) >= k.
!>
!>            d) NOTE: Any entry IPIV(k) is always NONZERO on output.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension ( MAX(1,LWORK) ).
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The length of WORK.  LWORK >=1.  For best performance
!>          LWORK >= N*NB, where NB is the block size returned
!>          by ILAENV.
!>
!>          If LWORK = -1, then a workspace query is assumed;
!>          the routine only calculates the optimal size of the WORK
!>          array, returns this value as the first entry of the WORK
!>          array, and no error message related to LWORK is issued
!>          by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>
!>          < 0: If INFO = -k, the k-th argument had an illegal value
!>
!>          > 0: If INFO = k, the matrix A is singular, because:
!>                 If UPLO = 'U': column k in the upper
!>                 triangular part of A contains all zeros.
!>                 If UPLO = 'L': column k in the lower
!>                 triangular part of A contains all zeros.
!>
!>               Therefore D(k,k) is exactly zero, and superdiagonal
!>               elements of column k of U (or subdiagonal elements of
!>               column k of L ) are all zeros. The factorization has
!>               been completed, but the block diagonal matrix D is
!>               exactly singular, and division by zero will occur if
!>               it is used to solve a system of equations.
!>
!>               NOTE: INFO only stores the first occurrence of
!>               a singularity, any subsequent occurrence of singularity
!>               is not stored in INFO even though the factorization
!>               always completes.
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
!> \ingroup singleSYcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!> TODO: put correct description
!> \endverbatim
!
!> \par Contributors:
!  ==================
!>
!> \verbatim
!>
!>  December 2016,  Igor Kozachenko,
!>                  Computer Science Division,
!>                  University of California, Berkeley
!>
!>  September 2007, Sven Hammarling, Nicholas J. Higham, Craig Lucas,
!>                  School of Mathematics,
!>                  University of Manchester
!>
!> \endverbatim
!
!  =====================================================================
      SUBROUTINE SSYTRF_RK(Uplo,N,A,Lda,E,Ipiv,Work,Lwork,Info)
      USE S_ILAENV
      USE S_LSAME
      USE S_SLASYF_RK
      USE S_SSWAP
      USE S_SSYTF2_RK
      USE S_XERBLA
      IMPLICIT NONE
!*--SSYTRF_RK268
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: E
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , iinfo , ip , iws , k , kb , ldwork , lwkopt , nb , &
     &           nbmin
      LOGICAL :: lquery , upper
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
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
      lquery = (Lwork==-1)
      IF ( .NOT.upper .AND. .NOT.LSAME(Uplo,'L') ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -4
      ELSEIF ( Lwork<1 .AND. .NOT.lquery ) THEN
         Info = -8
      ENDIF
!
      IF ( Info==0 ) THEN
!
!        Determine the block size
!
         nb = ILAENV(1,'SSYTRF_RK',Uplo,N,-1,-1,-1)
         lwkopt = N*nb
         Work(1) = lwkopt
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('SSYTRF_RK',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
      nbmin = 2
      ldwork = N
      IF ( nb>1 .AND. nb<N ) THEN
         iws = ldwork*nb
         IF ( Lwork<iws ) THEN
            nb = MAX(Lwork/ldwork,1)
            nbmin = MAX(2,ILAENV(2,'SSYTRF_RK',Uplo,N,-1,-1,-1))
         ENDIF
      ELSE
         iws = 1
      ENDIF
      IF ( nb<nbmin ) nb = N
!
      IF ( upper ) THEN
!
!        Factorize A as U*D*U**T using the upper triangle of A
!
!        K is the main loop index, decreasing from N to 1 in steps of
!        KB, where KB is the number of columns factorized by SLASYF_RK;
!        KB is either NB or NB-1, or K for the last block
!
         k = N
!
!        If K < 1, exit from loop
!
         DO WHILE ( k>=1 )
!
            IF ( k>nb ) THEN
!
!           Factorize columns k-kb+1:k of A and use blocked code to
!           update columns 1:k-kb
!
               CALL SLASYF_RK(Uplo,k,nb,kb,A,Lda,E,Ipiv,Work,ldwork,    &
     &                        iinfo)
            ELSE
!
!           Use unblocked code to factorize columns 1:k of A
!
               CALL SSYTF2_RK(Uplo,k,A,Lda,E,Ipiv,iinfo)
               kb = k
            ENDIF
!
!        Set INFO on the first occurrence of a zero pivot
!
            IF ( Info==0 .AND. iinfo>0 ) Info = iinfo
!
!        No need to adjust IPIV
!
!
!        Apply permutations to the leading panel 1:k-1
!
!        Read IPIV from the last block factored, i.e.
!        indices  k-kb+1:k and apply row permutations to the
!        last k+1 colunms k+1:N after that block
!        (We can do the simple loop over IPIV with decrement -1,
!        since the ABS value of IPIV( I ) represents the row index
!        of the interchange with row i in both 1x1 and 2x2 pivot cases)
!
            IF ( k<N ) THEN
               DO i = k , (k-kb+1) , -1
                  ip = ABS(Ipiv(i))
                  IF ( ip/=i ) CALL SSWAP(N-k,A(i,k+1),Lda,A(ip,k+1),   &
     &                 Lda)
               ENDDO
            ENDIF
!
!        Decrease K and return to the start of the main loop
!
            k = k - kb
         ENDDO
!
         Work(1) = lwkopt
!
!        This label is the exit from main loop over K decreasing
!        from N to 1 in steps of KB
!
!
      ELSE
!
!        Factorize A as L*D*L**T using the lower triangle of A
!
!        K is the main loop index, increasing from 1 to N in steps of
!        KB, where KB is the number of columns factorized by SLASYF_RK;
!        KB is either NB or NB-1, or N-K+1 for the last block
!
         k = 1
!
!        If K > N, exit from loop
!
         DO WHILE ( k<=N )
!
            IF ( k<=N-nb ) THEN
!
!           Factorize columns k:k+kb-1 of A and use blocked code to
!           update columns k+kb:n
!
               CALL SLASYF_RK(Uplo,N-k+1,nb,kb,A(k,k),Lda,E(k),Ipiv(k), &
     &                        Work,ldwork,iinfo)
 
 
            ELSE
!
!           Use unblocked code to factorize columns k:n of A
!
               CALL SSYTF2_RK(Uplo,N-k+1,A(k,k),Lda,E(k),Ipiv(k),iinfo)
               kb = N - k + 1
!
            ENDIF
!
!        Set INFO on the first occurrence of a zero pivot
!
            IF ( Info==0 .AND. iinfo>0 ) Info = iinfo + k - 1
!
!        Adjust IPIV
!
            DO i = k , k + kb - 1
               IF ( Ipiv(i)>0 ) THEN
                  Ipiv(i) = Ipiv(i) + k - 1
               ELSE
                  Ipiv(i) = Ipiv(i) - k + 1
               ENDIF
            ENDDO
!
!        Apply permutations to the leading panel 1:k-1
!
!        Read IPIV from the last block factored, i.e.
!        indices  k:k+kb-1 and apply row permutations to the
!        first k-1 colunms 1:k-1 before that block
!        (We can do the simple loop over IPIV with increment 1,
!        since the ABS value of IPIV( I ) represents the row index
!        of the interchange with row i in both 1x1 and 2x2 pivot cases)
!
            IF ( k>1 ) THEN
               DO i = k , (k+kb-1) , 1
                  ip = ABS(Ipiv(i))
                  IF ( ip/=i ) CALL SSWAP(k-1,A(i,1),Lda,A(ip,1),Lda)
               ENDDO
            ENDIF
!
!        Increase K and return to the start of the main loop
!
            k = k + kb
         ENDDO
         Work(1) = lwkopt
!
!        This label is the exit from main loop over K increasing
!        from 1 to N in steps of KB
!
!
!     End Lower
!
      ENDIF
!
!     End of SSYTRF_RK
!
      END SUBROUTINE SSYTRF_RK
