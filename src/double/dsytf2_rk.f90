!*==dsytf2_rk.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b DSYTF2_RK computes the factorization of a real symmetric indefinite matrix using the bounded Bunch-Kaufman (rook) diagonal pivoting method (BLAS2 unblocked algorithm).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DSYTF2_RK + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsytf2_rk.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsytf2_rk.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsytf2_rk.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSYTF2_RK( UPLO, N, A, LDA, E, IPIV, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDA, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       DOUBLE PRECISION   A( LDA, * ), E ( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!> DSYTF2_RK computes the factorization of a real symmetric matrix A
!> using the bounded Bunch-Kaufman (rook) diagonal pivoting method:
!>
!>    A = P*U*D*(U**T)*(P**T) or A = P*L*D*(L**T)*(P**T),
!>
!> where U (or L) is unit upper (or lower) triangular matrix,
!> U**T (or L**T) is the transpose of U (or L), P is a permutation
!> matrix, P**T is the transpose of P, and D is symmetric and block
!> diagonal with 1-by-1 and 2-by-2 diagonal blocks.
!>
!> This is the unblocked version of the algorithm, calling Level 2 BLAS.
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
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
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
!>          E is DOUBLE PRECISION array, dimension (N)
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
!> \ingroup doubleSYcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!> TODO: put further details
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
!>  01-01-96 - Based on modifications by
!>    J. Lewis, Boeing Computer Services Company
!>    A. Petitet, Computer Science Dept.,
!>                Univ. of Tenn., Knoxville abd , USA
!> \endverbatim
!
!  =====================================================================
      SUBROUTINE DSYTF2_RK(Uplo,N,A,Lda,E,Ipiv,Info)
      IMPLICIT NONE
!*--DSYTF2_RK245
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER Info , Lda , N
!     ..
!     .. Array Arguments ..
      INTEGER Ipiv(*)
      DOUBLE PRECISION A(Lda,*) , E(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
      DOUBLE PRECISION EIGHT , SEVTEN
      PARAMETER (EIGHT=8.0D+0,SEVTEN=17.0D+0)
!     ..
!     .. Local Scalars ..
      LOGICAL upper , done
      INTEGER i , imax , j , jmax , itemp , k , kk , kp , kstep , p , ii
      DOUBLE PRECISION absakk , alpha , colmax , d11 , d12 , d21 , d22 ,&
     &                 rowmax , dtemp , t , wk , wkm1 , wkp1 , sfmin
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER IDAMAX
      DOUBLE PRECISION DLAMCH
      EXTERNAL LSAME , IDAMAX , DLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL DSCAL , DSWAP , DSYR , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , SQRT
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
      IF ( Info/=0 ) THEN
         CALL XERBLA('DSYTF2_RK',-Info)
         RETURN
      ENDIF
!
!     Initialize ALPHA for use in choosing pivot block size.
!
      alpha = (ONE+SQRT(SEVTEN))/EIGHT
!
!     Compute machine safe minimum
!
      sfmin = DLAMCH('S')
!
      IF ( upper ) THEN
!
!        Factorize A as U*D*U**T using the upper triangle of A
!
!        Initialize the first entry of array E, where superdiagonal
!        elements of D are stored
!
         E(1) = ZERO
!
!        K is the main loop index, decreasing from N to 1 in steps of
!        1 or 2
!
         k = N
!
!        If K < 1, exit from loop
!
         DO WHILE ( k>=1 )
            kstep = 1
            p = k
!
!        Determine rows and columns to be interchanged and whether
!        a 1-by-1 or 2-by-2 pivot block will be used
!
            absakk = ABS(A(k,k))
!
!        IMAX is the row-index of the largest off-diagonal element in
!        column K, and COLMAX is its absolute value.
!        Determine both COLMAX and IMAX.
!
            IF ( k>1 ) THEN
               imax = IDAMAX(k-1,A(1,k),1)
               colmax = ABS(A(imax,k))
            ELSE
               colmax = ZERO
            ENDIF
!
            IF ( (MAX(absakk,colmax)==ZERO) ) THEN
!
!           Column K is zero or underflow: set INFO and continue
!
               IF ( Info==0 ) Info = k
               kp = k
!
!           Set E( K ) to zero
!
               IF ( k>1 ) E(k) = ZERO
!
            ELSE
!
!           Test for interchange
!
!           Equivalent to testing for (used to handle NaN and Inf)
!           ABSAKK.GE.ALPHA*COLMAX
!
               IF ( absakk>=alpha*colmax ) THEN
!
!              no interchange,
!              use 1-by-1 pivot block
!
                  kp = k
               ELSE
!
                  done = .FALSE.
                  DO
!
!              Loop until pivot found
!
!
!                 Begin pivot search loop body
!
!                 JMAX is the column-index of the largest off-diagonal
!                 element in row IMAX, and ROWMAX is its absolute value.
!                 Determine both ROWMAX and JMAX.
!
                     IF ( imax/=k ) THEN
                        jmax = imax + IDAMAX(k-imax,A(imax,imax+1),Lda)
                        rowmax = ABS(A(imax,jmax))
                     ELSE
                        rowmax = ZERO
                     ENDIF
!
                     IF ( imax>1 ) THEN
                        itemp = IDAMAX(imax-1,A(1,imax),1)
                        dtemp = ABS(A(itemp,imax))
                        IF ( dtemp>rowmax ) THEN
                           rowmax = dtemp
                           jmax = itemp
                        ENDIF
                     ENDIF
!
!                 Equivalent to testing for (used to handle NaN and Inf)
!                 ABS( A( IMAX, IMAX ) ).GE.ALPHA*ROWMAX
!
                     IF ( ABS(A(imax,imax))>=alpha*rowmax ) THEN
!
!                    interchange rows and columns K and IMAX,
!                    use 1-by-1 pivot block
!
                        kp = imax
                        done = .TRUE.
!
!                 Equivalent to testing for ROWMAX .EQ. COLMAX,
!                 used to handle NaN and Inf
!
                     ELSEIF ( (p==jmax) .OR. (rowmax<=colmax) ) THEN
!
!                    interchange rows and columns K+1 and IMAX,
!                    use 2-by-2 pivot block
!
                        kp = imax
                        kstep = 2
                        done = .TRUE.
                     ELSE
!
!                    Pivot NOT found, set variables and repeat
!
                        p = imax
                        colmax = rowmax
                        imax = jmax
                     ENDIF
!
!                 End pivot search loop body
!
                     IF ( done ) EXIT
                  ENDDO
!
               ENDIF
!
!           Swap TWO rows and TWO columns
!
!           First swap
!
               IF ( (kstep==2) .AND. (p/=k) ) THEN
!
!              Interchange rows and column K and P in the leading
!              submatrix A(1:k,1:k) if we have a 2-by-2 pivot
!
                  IF ( p>1 ) CALL DSWAP(p-1,A(1,k),1,A(1,p),1)
                  IF ( p<(k-1) ) CALL DSWAP(k-p-1,A(p+1,k),1,A(p,p+1),  &
     &                 Lda)
                  t = A(k,k)
                  A(k,k) = A(p,p)
                  A(p,p) = t
!
!              Convert upper triangle of A into U form by applying
!              the interchanges in columns k+1:N.
!
                  IF ( k<N ) CALL DSWAP(N-k,A(k,k+1),Lda,A(p,k+1),Lda)
!
               ENDIF
!
!           Second swap
!
               kk = k - kstep + 1
               IF ( kp/=kk ) THEN
!
!              Interchange rows and columns KK and KP in the leading
!              submatrix A(1:k,1:k)
!
                  IF ( kp>1 ) CALL DSWAP(kp-1,A(1,kk),1,A(1,kp),1)
                  IF ( (kk>1) .AND. (kp<(kk-1)) )                       &
     &                 CALL DSWAP(kk-kp-1,A(kp+1,kk),1,A(kp,kp+1),Lda)
                  t = A(kk,kk)
                  A(kk,kk) = A(kp,kp)
                  A(kp,kp) = t
                  IF ( kstep==2 ) THEN
                     t = A(k-1,k)
                     A(k-1,k) = A(kp,k)
                     A(kp,k) = t
                  ENDIF
!
!              Convert upper triangle of A into U form by applying
!              the interchanges in columns k+1:N.
!
                  IF ( k<N ) CALL DSWAP(N-k,A(kk,k+1),Lda,A(kp,k+1),Lda)
!
               ENDIF
!
!           Update the leading submatrix
!
               IF ( kstep/=1 ) THEN
!
!              2-by-2 pivot block D(k): columns k and k-1 now hold
!
!              ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)
!
!              where U(k) and U(k-1) are the k-th and (k-1)-th columns
!              of U
!
!              Perform a rank-2 update of A(1:k-2,1:k-2) as
!
!              A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )**T
!                 = A - ( ( A(k-1)A(k) )*inv(D(k)) ) * ( A(k-1)A(k) )**T
!
!              and store L(k) and L(k+1) in columns k and k+1
!
                  IF ( k>2 ) THEN
!
                     d12 = A(k-1,k)
                     d22 = A(k-1,k-1)/d12
                     d11 = A(k,k)/d12
                     t = ONE/(d11*d22-ONE)
!
                     DO j = k - 2 , 1 , -1
!
                        wkm1 = t*(d11*A(j,k-1)-A(j,k))
                        wk = t*(d22*A(j,k)-A(j,k-1))
!
                        DO i = j , 1 , -1
                           A(i,j) = A(i,j) - (A(i,k)/d12)               &
     &                              *wk - (A(i,k-1)/d12)*wkm1
                        ENDDO
!
!                    Store U(k) and U(k-1) in cols k and k-1 for row J
!
                        A(j,k) = wk/d12
                        A(j,k-1) = wkm1/d12
!
                     ENDDO
!
                  ENDIF
!
!              Copy superdiagonal elements of D(K) to E(K) and
!              ZERO out superdiagonal entry of A
!
                  E(k) = A(k-1,k)
                  E(k-1) = ZERO
                  A(k-1,k) = ZERO
!
!              1-by-1 pivot block D(k): column k now holds
!
!              W(k) = U(k)*D(k)
!
!              where U(k) is the k-th column of U
!
               ELSEIF ( k>1 ) THEN
!
!                 Perform a rank-1 update of A(1:k-1,1:k-1) and
!                 store U(k) in column k
!
                  IF ( ABS(A(k,k))>=sfmin ) THEN
!
!                    Perform a rank-1 update of A(1:k-1,1:k-1) as
!                    A := A - U(k)*D(k)*U(k)**T
!                       = A - W(k)*1/D(k)*W(k)**T
!
                     d11 = ONE/A(k,k)
                     CALL DSYR(Uplo,k-1,-d11,A(1,k),1,A,Lda)
!
!                    Store U(k) in column k
!
                     CALL DSCAL(k-1,d11,A(1,k),1)
                  ELSE
!
!                    Store L(k) in column K
!
                     d11 = A(k,k)
                     DO ii = 1 , k - 1
                        A(ii,k) = A(ii,k)/d11
                     ENDDO
!
!                    Perform a rank-1 update of A(k+1:n,k+1:n) as
!                    A := A - U(k)*D(k)*U(k)**T
!                       = A - W(k)*(1/D(k))*W(k)**T
!                       = A - (W(k)/D(k))*(D(k))*(W(k)/D(K))**T
!
                     CALL DSYR(Uplo,k-1,-d11,A(1,k),1,A,Lda)
                  ENDIF
!
!                 Store the superdiagonal element of D in array E
!
                  E(k) = ZERO
!
               ENDIF
!
!           End column K is nonsingular
!
            ENDIF
!
!        Store details of the interchanges in IPIV
!
            IF ( kstep==1 ) THEN
!
!
               Ipiv(k) = kp
            ELSE
               Ipiv(k) = -p
               Ipiv(k-1) = -kp
            ENDIF
!
!        Decrease K and return to the start of the main loop
!
            k = k - kstep
         ENDDO
!
!
      ELSE
!
!        Factorize A as L*D*L**T using the lower triangle of A
!
!        Initialize the unused last entry of the subdiagonal array E.
!
         E(N) = ZERO
!
!        K is the main loop index, increasing from 1 to N in steps of
!        1 or 2
!
         k = 1
!
!        If K > N, exit from loop
!
         DO WHILE ( k<=N )
            kstep = 1
            p = k
!
!        Determine rows and columns to be interchanged and whether
!        a 1-by-1 or 2-by-2 pivot block will be used
!
            absakk = ABS(A(k,k))
!
!        IMAX is the row-index of the largest off-diagonal element in
!        column K, and COLMAX is its absolute value.
!        Determine both COLMAX and IMAX.
!
            IF ( k<N ) THEN
               imax = k + IDAMAX(N-k,A(k+1,k),1)
               colmax = ABS(A(imax,k))
            ELSE
               colmax = ZERO
            ENDIF
!
            IF ( (MAX(absakk,colmax)==ZERO) ) THEN
!
!           Column K is zero or underflow: set INFO and continue
!
               IF ( Info==0 ) Info = k
               kp = k
!
!           Set E( K ) to zero
!
               IF ( k<N ) E(k) = ZERO
!
            ELSE
!
!           Test for interchange
!
!           Equivalent to testing for (used to handle NaN and Inf)
!           ABSAKK.GE.ALPHA*COLMAX
!
               IF ( absakk>=alpha*colmax ) THEN
!
!              no interchange, use 1-by-1 pivot block
!
                  kp = k
!
               ELSE
!
                  done = .FALSE.
                  DO
!
!              Loop until pivot found
!
!
!                 Begin pivot search loop body
!
!                 JMAX is the column-index of the largest off-diagonal
!                 element in row IMAX, and ROWMAX is its absolute value.
!                 Determine both ROWMAX and JMAX.
!
                     IF ( imax/=k ) THEN
                        jmax = k - 1 + IDAMAX(imax-k,A(imax,k),Lda)
                        rowmax = ABS(A(imax,jmax))
                     ELSE
                        rowmax = ZERO
                     ENDIF
!
                     IF ( imax<N ) THEN
                        itemp = imax + IDAMAX(N-imax,A(imax+1,imax),1)
                        dtemp = ABS(A(itemp,imax))
                        IF ( dtemp>rowmax ) THEN
                           rowmax = dtemp
                           jmax = itemp
                        ENDIF
                     ENDIF
!
!                 Equivalent to testing for (used to handle NaN and Inf)
!                 ABS( A( IMAX, IMAX ) ).GE.ALPHA*ROWMAX
!
                     IF ( ABS(A(imax,imax))>=alpha*rowmax ) THEN
!
!                    interchange rows and columns K and IMAX,
!                    use 1-by-1 pivot block
!
                        kp = imax
                        done = .TRUE.
!
!                 Equivalent to testing for ROWMAX .EQ. COLMAX,
!                 used to handle NaN and Inf
!
                     ELSEIF ( (p==jmax) .OR. (rowmax<=colmax) ) THEN
!
!                    interchange rows and columns K+1 and IMAX,
!                    use 2-by-2 pivot block
!
                        kp = imax
                        kstep = 2
                        done = .TRUE.
                     ELSE
!
!                    Pivot NOT found, set variables and repeat
!
                        p = imax
                        colmax = rowmax
                        imax = jmax
                     ENDIF
!
!                 End pivot search loop body
!
                     IF ( done ) EXIT
                  ENDDO
!
               ENDIF
!
!           Swap TWO rows and TWO columns
!
!           First swap
!
               IF ( (kstep==2) .AND. (p/=k) ) THEN
!
!              Interchange rows and column K and P in the trailing
!              submatrix A(k:n,k:n) if we have a 2-by-2 pivot
!
                  IF ( p<N ) CALL DSWAP(N-p,A(p+1,k),1,A(p+1,p),1)
                  IF ( p>(k+1) ) CALL DSWAP(p-k-1,A(k+1,k),1,A(p,k+1),  &
     &                 Lda)
                  t = A(k,k)
                  A(k,k) = A(p,p)
                  A(p,p) = t
!
!              Convert lower triangle of A into L form by applying
!              the interchanges in columns 1:k-1.
!
                  IF ( k>1 ) CALL DSWAP(k-1,A(k,1),Lda,A(p,1),Lda)
!
               ENDIF
!
!           Second swap
!
               kk = k + kstep - 1
               IF ( kp/=kk ) THEN
!
!              Interchange rows and columns KK and KP in the trailing
!              submatrix A(k:n,k:n)
!
                  IF ( kp<N ) CALL DSWAP(N-kp,A(kp+1,kk),1,A(kp+1,kp),1)
                  IF ( (kk<N) .AND. (kp>(kk+1)) )                       &
     &                 CALL DSWAP(kp-kk-1,A(kk+1,kk),1,A(kp,kk+1),Lda)
                  t = A(kk,kk)
                  A(kk,kk) = A(kp,kp)
                  A(kp,kp) = t
                  IF ( kstep==2 ) THEN
                     t = A(k+1,k)
                     A(k+1,k) = A(kp,k)
                     A(kp,k) = t
                  ENDIF
!
!              Convert lower triangle of A into L form by applying
!              the interchanges in columns 1:k-1.
!
                  IF ( k>1 ) CALL DSWAP(k-1,A(kk,1),Lda,A(kp,1),Lda)
!
               ENDIF
!
!           Update the trailing submatrix
!
               IF ( kstep/=1 ) THEN
!
!              2-by-2 pivot block D(k): columns k and k+1 now hold
!
!              ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k)
!
!              where L(k) and L(k+1) are the k-th and (k+1)-th columns
!              of L
!
!
!              Perform a rank-2 update of A(k+2:n,k+2:n) as
!
!              A := A - ( L(k) L(k+1) ) * D(k) * ( L(k) L(k+1) )**T
!                 = A - ( ( A(k)A(k+1) )*inv(D(k) ) * ( A(k)A(k+1) )**T
!
!              and store L(k) and L(k+1) in columns k and k+1
!
                  IF ( k<N-1 ) THEN
!
                     d21 = A(k+1,k)
                     d11 = A(k+1,k+1)/d21
                     d22 = A(k,k)/d21
                     t = ONE/(d11*d22-ONE)
!
                     DO j = k + 2 , N
!
!                    Compute  D21 * ( W(k)W(k+1) ) * inv(D(k)) for row J
!
                        wk = t*(d11*A(j,k)-A(j,k+1))
                        wkp1 = t*(d22*A(j,k+1)-A(j,k))
!
!                    Perform a rank-2 update of A(k+2:n,k+2:n)
!
                        DO i = j , N
                           A(i,j) = A(i,j) - (A(i,k)/d21)               &
     &                              *wk - (A(i,k+1)/d21)*wkp1
                        ENDDO
!
!                    Store L(k) and L(k+1) in cols k and k+1 for row J
!
                        A(j,k) = wk/d21
                        A(j,k+1) = wkp1/d21
!
                     ENDDO
!
                  ENDIF
!
!              Copy subdiagonal elements of D(K) to E(K) and
!              ZERO out subdiagonal entry of A
!
                  E(k) = A(k+1,k)
                  E(k+1) = ZERO
                  A(k+1,k) = ZERO
!
!              1-by-1 pivot block D(k): column k now holds
!
!              W(k) = L(k)*D(k)
!
!              where L(k) is the k-th column of L
!
               ELSEIF ( k<N ) THEN
!
!              Perform a rank-1 update of A(k+1:n,k+1:n) and
!              store L(k) in column k
!
                  IF ( ABS(A(k,k))>=sfmin ) THEN
!
!                    Perform a rank-1 update of A(k+1:n,k+1:n) as
!                    A := A - L(k)*D(k)*L(k)**T
!                       = A - W(k)*(1/D(k))*W(k)**T
!
                     d11 = ONE/A(k,k)
                     CALL DSYR(Uplo,N-k,-d11,A(k+1,k),1,A(k+1,k+1),Lda)
!
!                    Store L(k) in column k
!
                     CALL DSCAL(N-k,d11,A(k+1,k),1)
                  ELSE
!
!                    Store L(k) in column k
!
                     d11 = A(k,k)
                     DO ii = k + 1 , N
                        A(ii,k) = A(ii,k)/d11
                     ENDDO
!
!                    Perform a rank-1 update of A(k+1:n,k+1:n) as
!                    A := A - L(k)*D(k)*L(k)**T
!                       = A - W(k)*(1/D(k))*W(k)**T
!                       = A - (W(k)/D(k))*(D(k))*(W(k)/D(K))**T
!
                     CALL DSYR(Uplo,N-k,-d11,A(k+1,k),1,A(k+1,k+1),Lda)
                  ENDIF
!
!                 Store the subdiagonal element of D in array E
!
                  E(k) = ZERO
!
               ENDIF
!
!           End column K is nonsingular
!
            ENDIF
!
!        Store details of the interchanges in IPIV
!
            IF ( kstep==1 ) THEN
!
!
               Ipiv(k) = kp
            ELSE
               Ipiv(k) = -p
               Ipiv(k+1) = -kp
            ENDIF
!
!        Increase K and return to the start of the main loop
!
            k = k + kstep
         ENDDO
!
!
      ENDIF
!
!
!     End of DSYTF2_RK
!
      END SUBROUTINE DSYTF2_RK
