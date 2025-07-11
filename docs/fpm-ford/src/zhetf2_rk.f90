!*==zhetf2_rk.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b ZHETF2_RK computes the factorization of a complex Hermitian indefinite matrix using the bounded Bunch-Kaufman (rook) diagonal pivoting method (BLAS2 unblocked algorithm).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZHETF2_RK + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhetf2_rk.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhetf2_rk.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhetf2_rk.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZHETF2_RK( UPLO, N, A, LDA, E, IPIV, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDA, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX*16         A( LDA, * ), E ( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!> ZHETF2_RK computes the factorization of a complex Hermitian matrix A
!> using the bounded Bunch-Kaufman (rook) diagonal pivoting method:
!>
!>    A = P*U*D*(U**H)*(P**T) or A = P*L*D*(L**H)*(P**T),
!>
!> where U (or L) is unit upper (or lower) triangular matrix,
!> U**H (or L**H) is the conjugate of U (or L), P is a permutation
!> matrix, P**T is the transpose of P, and D is Hermitian and block
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
!>          Hermitian matrix A is stored:
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
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          On entry, the Hermitian matrix A.
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
!>            a) ONLY diagonal elements of the Hermitian block diagonal
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
!>          E is COMPLEX*16 array, dimension (N)
!>          On exit, contains the superdiagonal (or subdiagonal)
!>          elements of the Hermitian block diagonal matrix D
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
!>          the Hermitian block diagonal matrix D with 1-by-1 or 2-by-2
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
!> \ingroup complex16HEcomputational
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
      SUBROUTINE ZHETF2_RK(Uplo,N,A,Lda,E,Ipiv,Info)
      IMPLICIT NONE
!*--ZHETF2_RK245
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
      COMPLEX*16 A(Lda,*) , E(*)
!     ..
!
!  ======================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
      DOUBLE PRECISION EIGHT , SEVTEN
      PARAMETER (EIGHT=8.0D+0,SEVTEN=17.0D+0)
      COMPLEX*16 CZERO
      PARAMETER (CZERO=(0.0D+0,0.0D+0))
!     ..
!     .. Local Scalars ..
      LOGICAL done , upper
      INTEGER i , ii , imax , itemp , j , jmax , k , kk , kp , kstep , p
      DOUBLE PRECISION absakk , alpha , colmax , d , d11 , d22 , r1 ,   &
     &                 dtemp , rowmax , tt , sfmin
      COMPLEX*16 d12 , d21 , t , wk , wkm1 , wkp1 , z
!     ..
!     .. External Functions ..
!
      LOGICAL LSAME
      INTEGER IZAMAX
      DOUBLE PRECISION DLAMCH , DLAPY2
      EXTERNAL LSAME , IZAMAX , DLAMCH , DLAPY2
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA , ZDSCAL , ZHER , ZSWAP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , DBLE , DCMPLX , DCONJG , DIMAG , MAX , SQRT
!     ..
!     .. Statement Functions ..
      DOUBLE PRECISION CABS1
!     ..
!     .. Statement Function definitions ..
      CABS1(z) = ABS(DBLE(z)) + ABS(DIMAG(z))
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
         CALL XERBLA('ZHETF2_RK',-Info)
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
!        Factorize A as U*D*U**H using the upper triangle of A
!
!        Initialize the first entry of array E, where superdiagonal
!        elements of D are stored
!
         E(1) = CZERO
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
            absakk = ABS(DBLE(A(k,k)))
!
!        IMAX is the row-index of the largest off-diagonal element in
!        column K, and COLMAX is its absolute value.
!        Determine both COLMAX and IMAX.
!
            IF ( k>1 ) THEN
               imax = IZAMAX(k-1,A(1,k),1)
               colmax = CABS1(A(imax,k))
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
               A(k,k) = DBLE(A(k,k))
!
!           Set E( K ) to zero
!
               IF ( k>1 ) E(k) = CZERO
!
            ELSE
!
!           ============================================================
!
!           BEGIN pivot search
!
!           Case(1)
!           Equivalent to testing for ABSAKK.GE.ALPHA*COLMAX
!           (used to handle NaN and Inf)
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
!                 BEGIN pivot search loop body
!
!
!                 JMAX is the column-index of the largest off-diagonal
!                 element in row IMAX, and ROWMAX is its absolute value.
!                 Determine both ROWMAX and JMAX.
!
                     IF ( imax/=k ) THEN
                        jmax = imax + IZAMAX(k-imax,A(imax,imax+1),Lda)
                        rowmax = CABS1(A(imax,jmax))
                     ELSE
                        rowmax = ZERO
                     ENDIF
!
                     IF ( imax>1 ) THEN
                        itemp = IZAMAX(imax-1,A(1,imax),1)
                        dtemp = CABS1(A(itemp,imax))
                        IF ( dtemp>rowmax ) THEN
                           rowmax = dtemp
                           jmax = itemp
                        ENDIF
                     ENDIF
!
!                 Case(2)
!                 Equivalent to testing for
!                 ABS( REAL( W( IMAX,KW-1 ) ) ).GE.ALPHA*ROWMAX
!                 (used to handle NaN and Inf)
!
                     IF ( ABS(DBLE(A(imax,imax)))>=alpha*rowmax ) THEN
!
!                    interchange rows and columns K and IMAX,
!                    use 1-by-1 pivot block
!
                        kp = imax
                        done = .TRUE.
!
!                 Case(3)
!                 Equivalent to testing for ROWMAX.EQ.COLMAX,
!                 (used to handle NaN and Inf)
!
                     ELSEIF ( (p==jmax) .OR. (rowmax<=colmax) ) THEN
!
!                    interchange rows and columns K-1 and IMAX,
!                    use 2-by-2 pivot block
!
                        kp = imax
                        kstep = 2
                        done = .TRUE.
!
!                 Case(4)
                     ELSE
!
!                    Pivot not found: set params and repeat
!
                        p = imax
                        colmax = rowmax
                        imax = jmax
                     ENDIF
!
!                 END pivot search loop body
!
                     IF ( done ) EXIT
                  ENDDO
!
               ENDIF
!
!           END pivot search
!
!           ============================================================
!
!           KK is the column of A where pivoting step stopped
!
               kk = k - kstep + 1
!
!           For only a 2x2 pivot, interchange rows and columns K and P
!           in the leading submatrix A(1:k,1:k)
!
               IF ( (kstep==2) .AND. (p/=k) ) THEN
!              (1) Swap columnar parts
                  IF ( p>1 ) CALL ZSWAP(p-1,A(1,k),1,A(1,p),1)
!              (2) Swap and conjugate middle parts
                  DO j = p + 1 , k - 1
                     t = DCONJG(A(j,k))
                     A(j,k) = DCONJG(A(p,j))
                     A(p,j) = t
                  ENDDO
!              (3) Swap and conjugate corner elements at row-col interserction
                  A(p,k) = DCONJG(A(p,k))
!              (4) Swap diagonal elements at row-col intersection
                  r1 = DBLE(A(k,k))
                  A(k,k) = DBLE(A(p,p))
                  A(p,p) = r1
!
!              Convert upper triangle of A into U form by applying
!              the interchanges in columns k+1:N.
!
                  IF ( k<N ) CALL ZSWAP(N-k,A(k,k+1),Lda,A(p,k+1),Lda)
!
               ENDIF
!
!           For both 1x1 and 2x2 pivots, interchange rows and
!           columns KK and KP in the leading submatrix A(1:k,1:k)
!
               IF ( kp/=kk ) THEN
!              (1) Swap columnar parts
                  IF ( kp>1 ) CALL ZSWAP(kp-1,A(1,kk),1,A(1,kp),1)
!              (2) Swap and conjugate middle parts
                  DO j = kp + 1 , kk - 1
                     t = DCONJG(A(j,kk))
                     A(j,kk) = DCONJG(A(kp,j))
                     A(kp,j) = t
                  ENDDO
!              (3) Swap and conjugate corner elements at row-col interserction
                  A(kp,kk) = DCONJG(A(kp,kk))
!              (4) Swap diagonal elements at row-col intersection
                  r1 = DBLE(A(kk,kk))
                  A(kk,kk) = DBLE(A(kp,kp))
                  A(kp,kp) = r1
!
                  IF ( kstep==2 ) THEN
!                 (*) Make sure that diagonal element of pivot is real
                     A(k,k) = DBLE(A(k,k))
!                 (5) Swap row elements
                     t = A(k-1,k)
                     A(k-1,k) = A(kp,k)
                     A(kp,k) = t
                  ENDIF
!
!              Convert upper triangle of A into U form by applying
!              the interchanges in columns k+1:N.
!
                  IF ( k<N ) CALL ZSWAP(N-k,A(kk,k+1),Lda,A(kp,k+1),Lda)
!
               ELSE
!              (*) Make sure that diagonal element of pivot is real
                  A(k,k) = DBLE(A(k,k))
                  IF ( kstep==2 ) A(k-1,k-1) = DBLE(A(k-1,k-1))
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
!                 D = |A12|
                     d = DLAPY2(DBLE(A(k-1,k)),DIMAG(A(k-1,k)))
                     d11 = A(k,k)/d
                     d22 = A(k-1,k-1)/d
                     d12 = A(k-1,k)/d
                     tt = ONE/(d11*d22-ONE)
!
                     DO j = k - 2 , 1 , -1
!
!                    Compute  D21 * ( W(k)W(k+1) ) * inv(D(k)) for row J
!
                        wkm1 = tt*(d11*A(j,k-1)-DCONJG(d12)*A(j,k))
                        wk = tt*(d22*A(j,k)-d12*A(j,k-1))
!
!                    Perform a rank-2 update of A(1:k-2,1:k-2)
!
                        DO i = j , 1 , -1
                           A(i,j) = A(i,j) - (A(i,k)/d)*DCONJG(wk)      &
     &                              - (A(i,k-1)/d)*DCONJG(wkm1)
                        ENDDO
!
!                    Store U(k) and U(k-1) in cols k and k-1 for row J
!
                        A(j,k) = wk/d
                        A(j,k-1) = wkm1/d
!                    (*) Make sure that diagonal element of pivot is real
                        A(j,j) = DCMPLX(DBLE(A(j,j)),ZERO)
!
                     ENDDO
!
                  ENDIF
!
!              Copy superdiagonal elements of D(K) to E(K) and
!              ZERO out superdiagonal entry of A
!
                  E(k) = A(k-1,k)
                  E(k-1) = CZERO
                  A(k-1,k) = CZERO
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
                  IF ( ABS(DBLE(A(k,k)))>=sfmin ) THEN
!
!                    Perform a rank-1 update of A(1:k-1,1:k-1) as
!                    A := A - U(k)*D(k)*U(k)**T
!                       = A - W(k)*1/D(k)*W(k)**T
!
                     d11 = ONE/DBLE(A(k,k))
                     CALL ZHER(Uplo,k-1,-d11,A(1,k),1,A,Lda)
!
!                    Store U(k) in column k
!
                     CALL ZDSCAL(k-1,d11,A(1,k),1)
                  ELSE
!
!                    Store L(k) in column K
!
                     d11 = DBLE(A(k,k))
                     DO ii = 1 , k - 1
                        A(ii,k) = A(ii,k)/d11
                     ENDDO
!
!                    Perform a rank-1 update of A(k+1:n,k+1:n) as
!                    A := A - U(k)*D(k)*U(k)**T
!                       = A - W(k)*(1/D(k))*W(k)**T
!                       = A - (W(k)/D(k))*(D(k))*(W(k)/D(K))**T
!
                     CALL ZHER(Uplo,k-1,-d11,A(1,k),1,A,Lda)
                  ENDIF
!
!                 Store the superdiagonal element of D in array E
!
                  E(k) = CZERO
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
!        Factorize A as L*D*L**H using the lower triangle of A
!
!        Initialize the unused last entry of the subdiagonal array E.
!
         E(N) = CZERO
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
            absakk = ABS(DBLE(A(k,k)))
!
!        IMAX is the row-index of the largest off-diagonal element in
!        column K, and COLMAX is its absolute value.
!        Determine both COLMAX and IMAX.
!
            IF ( k<N ) THEN
               imax = k + IZAMAX(N-k,A(k+1,k),1)
               colmax = CABS1(A(imax,k))
            ELSE
               colmax = ZERO
            ENDIF
!
            IF ( MAX(absakk,colmax)==ZERO ) THEN
!
!           Column K is zero or underflow: set INFO and continue
!
               IF ( Info==0 ) Info = k
               kp = k
               A(k,k) = DBLE(A(k,k))
!
!           Set E( K ) to zero
!
               IF ( k<N ) E(k) = CZERO
!
            ELSE
!
!           ============================================================
!
!           BEGIN pivot search
!
!           Case(1)
!           Equivalent to testing for ABSAKK.GE.ALPHA*COLMAX
!           (used to handle NaN and Inf)
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
!                 BEGIN pivot search loop body
!
!
!                 JMAX is the column-index of the largest off-diagonal
!                 element in row IMAX, and ROWMAX is its absolute value.
!                 Determine both ROWMAX and JMAX.
!
                     IF ( imax/=k ) THEN
                        jmax = k - 1 + IZAMAX(imax-k,A(imax,k),Lda)
                        rowmax = CABS1(A(imax,jmax))
                     ELSE
                        rowmax = ZERO
                     ENDIF
!
                     IF ( imax<N ) THEN
                        itemp = imax + IZAMAX(N-imax,A(imax+1,imax),1)
                        dtemp = CABS1(A(itemp,imax))
                        IF ( dtemp>rowmax ) THEN
                           rowmax = dtemp
                           jmax = itemp
                        ENDIF
                     ENDIF
!
!                 Case(2)
!                 Equivalent to testing for
!                 ABS( REAL( W( IMAX,KW-1 ) ) ).GE.ALPHA*ROWMAX
!                 (used to handle NaN and Inf)
!
                     IF ( ABS(DBLE(A(imax,imax)))>=alpha*rowmax ) THEN
!
!                    interchange rows and columns K and IMAX,
!                    use 1-by-1 pivot block
!
                        kp = imax
                        done = .TRUE.
!
!                 Case(3)
!                 Equivalent to testing for ROWMAX.EQ.COLMAX,
!                 (used to handle NaN and Inf)
!
                     ELSEIF ( (p==jmax) .OR. (rowmax<=colmax) ) THEN
!
!                    interchange rows and columns K+1 and IMAX,
!                    use 2-by-2 pivot block
!
                        kp = imax
                        kstep = 2
                        done = .TRUE.
!
!                 Case(4)
                     ELSE
!
!                    Pivot not found: set params and repeat
!
                        p = imax
                        colmax = rowmax
                        imax = jmax
                     ENDIF
!
!
!                 END pivot search loop body
!
                     IF ( done ) EXIT
                  ENDDO
!
               ENDIF
!
!           END pivot search
!
!           ============================================================
!
!           KK is the column of A where pivoting step stopped
!
               kk = k + kstep - 1
!
!           For only a 2x2 pivot, interchange rows and columns K and P
!           in the trailing submatrix A(k:n,k:n)
!
               IF ( (kstep==2) .AND. (p/=k) ) THEN
!              (1) Swap columnar parts
                  IF ( p<N ) CALL ZSWAP(N-p,A(p+1,k),1,A(p+1,p),1)
!              (2) Swap and conjugate middle parts
                  DO j = k + 1 , p - 1
                     t = DCONJG(A(j,k))
                     A(j,k) = DCONJG(A(p,j))
                     A(p,j) = t
                  ENDDO
!              (3) Swap and conjugate corner elements at row-col interserction
                  A(p,k) = DCONJG(A(p,k))
!              (4) Swap diagonal elements at row-col intersection
                  r1 = DBLE(A(k,k))
                  A(k,k) = DBLE(A(p,p))
                  A(p,p) = r1
!
!              Convert lower triangle of A into L form by applying
!              the interchanges in columns 1:k-1.
!
                  IF ( k>1 ) CALL ZSWAP(k-1,A(k,1),Lda,A(p,1),Lda)
!
               ENDIF
!
!           For both 1x1 and 2x2 pivots, interchange rows and
!           columns KK and KP in the trailing submatrix A(k:n,k:n)
!
               IF ( kp/=kk ) THEN
!              (1) Swap columnar parts
                  IF ( kp<N ) CALL ZSWAP(N-kp,A(kp+1,kk),1,A(kp+1,kp),1)
!              (2) Swap and conjugate middle parts
                  DO j = kk + 1 , kp - 1
                     t = DCONJG(A(j,kk))
                     A(j,kk) = DCONJG(A(kp,j))
                     A(kp,j) = t
                  ENDDO
!              (3) Swap and conjugate corner elements at row-col interserction
                  A(kp,kk) = DCONJG(A(kp,kk))
!              (4) Swap diagonal elements at row-col intersection
                  r1 = DBLE(A(kk,kk))
                  A(kk,kk) = DBLE(A(kp,kp))
                  A(kp,kp) = r1
!
                  IF ( kstep==2 ) THEN
!                 (*) Make sure that diagonal element of pivot is real
                     A(k,k) = DBLE(A(k,k))
!                 (5) Swap row elements
                     t = A(k+1,k)
                     A(k+1,k) = A(kp,k)
                     A(kp,k) = t
                  ENDIF
!
!              Convert lower triangle of A into L form by applying
!              the interchanges in columns 1:k-1.
!
                  IF ( k>1 ) CALL ZSWAP(k-1,A(kk,1),Lda,A(kp,1),Lda)
!
               ELSE
!              (*) Make sure that diagonal element of pivot is real
                  A(k,k) = DBLE(A(k,k))
                  IF ( kstep==2 ) A(k+1,k+1) = DBLE(A(k+1,k+1))
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
!                 D = |A21|
                     d = DLAPY2(DBLE(A(k+1,k)),DIMAG(A(k+1,k)))
                     d11 = DBLE(A(k+1,k+1))/d
                     d22 = DBLE(A(k,k))/d
                     d21 = A(k+1,k)/d
                     tt = ONE/(d11*d22-ONE)
!
                     DO j = k + 2 , N
!
!                    Compute  D21 * ( W(k)W(k+1) ) * inv(D(k)) for row J
!
                        wk = tt*(d11*A(j,k)-d21*A(j,k+1))
                        wkp1 = tt*(d22*A(j,k+1)-DCONJG(d21)*A(j,k))
!
!                    Perform a rank-2 update of A(k+2:n,k+2:n)
!
                        DO i = j , N
                           A(i,j) = A(i,j) - (A(i,k)/d)*DCONJG(wk)      &
     &                              - (A(i,k+1)/d)*DCONJG(wkp1)
                        ENDDO
!
!                    Store L(k) and L(k+1) in cols k and k+1 for row J
!
                        A(j,k) = wk/d
                        A(j,k+1) = wkp1/d
!                    (*) Make sure that diagonal element of pivot is real
                        A(j,j) = DCMPLX(DBLE(A(j,j)),ZERO)
!
                     ENDDO
!
                  ENDIF
!
!              Copy subdiagonal elements of D(K) to E(K) and
!              ZERO out subdiagonal entry of A
!
                  E(k) = A(k+1,k)
                  E(k+1) = CZERO
                  A(k+1,k) = CZERO
!
!              1-by-1 pivot block D(k): column k of A now holds
!
!              W(k) = L(k)*D(k),
!
!              where L(k) is the k-th column of L
!
               ELSEIF ( k<N ) THEN
!
!                 Perform a rank-1 update of A(k+1:n,k+1:n) and
!                 store L(k) in column k
!
!                 Handle division by a small number
!
                  IF ( ABS(DBLE(A(k,k)))>=sfmin ) THEN
!
!                    Perform a rank-1 update of A(k+1:n,k+1:n) as
!                    A := A - L(k)*D(k)*L(k)**T
!                       = A - W(k)*(1/D(k))*W(k)**T
!
                     d11 = ONE/DBLE(A(k,k))
                     CALL ZHER(Uplo,N-k,-d11,A(k+1,k),1,A(k+1,k+1),Lda)
!
!                    Store L(k) in column k
!
                     CALL ZDSCAL(N-k,d11,A(k+1,k),1)
                  ELSE
!
!                    Store L(k) in column k
!
                     d11 = DBLE(A(k,k))
                     DO ii = k + 1 , N
                        A(ii,k) = A(ii,k)/d11
                     ENDDO
!
!                    Perform a rank-1 update of A(k+1:n,k+1:n) as
!                    A := A - L(k)*D(k)*L(k)**T
!                       = A - W(k)*(1/D(k))*W(k)**T
!                       = A - (W(k)/D(k))*(D(k))*(W(k)/D(K))**T
!
                     CALL ZHER(Uplo,N-k,-d11,A(k+1,k),1,A(k+1,k+1),Lda)
                  ENDIF
!
!                 Store the subdiagonal element of D in array E
!
                  E(k) = CZERO
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
!     End of ZHETF2_RK
!
      END SUBROUTINE ZHETF2_RK
