!*==ssytf2_rook.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SSYTF2_ROOK computes the factorization of a real symmetric indefinite matrix using the bounded Bunch-Kaufman ("rook") diagonal pivoting method (unblocked algorithm).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SSYTF2_ROOK + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssytf2_rook.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssytf2_rook.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssytf2_rook.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SSYTF2_ROOK( UPLO, N, A, LDA, IPIV, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDA, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       REAL               A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SSYTF2_ROOK computes the factorization of a real symmetric matrix A
!> using the bounded Bunch-Kaufman ("rook") diagonal pivoting method:
!>
!>    A = U*D*U**T  or  A = L*D*L**T
!>
!> where U (or L) is a product of permutation and unit upper (lower)
!> triangular matrices, U**T is the transpose of U, and D is symmetric and
!> block diagonal with 1-by-1 and 2-by-2 diagonal blocks.
!>
!> This is the unblocked version of the algorithm, calling Level 2 BLAS.
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
!>          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!>          n-by-n upper triangular part of A contains the upper
!>          triangular part of the matrix A, and the strictly lower
!>          triangular part of A is not referenced.  If UPLO = 'L', the
!>          leading n-by-n lower triangular part of A contains the lower
!>          triangular part of the matrix A, and the strictly upper
!>          triangular part of A is not referenced.
!>
!>          On exit, the block diagonal matrix D and the multipliers used
!>          to obtain the factor U or L (see below for further details).
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
!>          Details of the interchanges and the block structure of D.
!>
!>          If UPLO = 'U':
!>             If IPIV(k) > 0, then rows and columns k and IPIV(k)
!>             were interchanged and D(k,k) is a 1-by-1 diagonal block.
!>
!>             If IPIV(k) < 0 and IPIV(k-1) < 0, then rows and
!>             columns k and -IPIV(k) were interchanged and rows and
!>             columns k-1 and -IPIV(k-1) were inerchaged,
!>             D(k-1:k,k-1:k) is a 2-by-2 diagonal block.
!>
!>          If UPLO = 'L':
!>             If IPIV(k) > 0, then rows and columns k and IPIV(k)
!>             were interchanged and D(k,k) is a 1-by-1 diagonal block.
!>
!>             If IPIV(k) < 0 and IPIV(k+1) < 0, then rows and
!>             columns k and -IPIV(k) were interchanged and rows and
!>             columns k+1 and -IPIV(k+1) were inerchaged,
!>             D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -k, the k-th argument had an illegal value
!>          > 0: if INFO = k, D(k,k) is exactly zero.  The factorization
!>               has been completed, but the block diagonal matrix D is
!>               exactly singular, and division by zero will occur if it
!>               is used to solve a system of equations.
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
!> \date November 2013
!
!> \ingroup realSYcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  If UPLO = 'U', then A = U*D*U**T, where
!>     U = P(n)*U(n)* ... *P(k)U(k)* ...,
!>  i.e., U is a product of terms P(k)*U(k), where k decreases from n to
!>  1 in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
!>  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
!>  defined by IPIV(k), and U(k) is a unit upper triangular matrix, such
!>  that if the diagonal block D(k) is of order s (s = 1 or 2), then
!>
!>             (   I    v    0   )   k-s
!>     U(k) =  (   0    I    0   )   s
!>             (   0    0    I   )   n-k
!>                k-s   s   n-k
!>
!>  If s = 1, D(k) overwrites A(k,k), and v overwrites A(1:k-1,k).
!>  If s = 2, the upper triangle of D(k) overwrites A(k-1,k-1), A(k-1,k),
!>  and A(k,k), and v overwrites A(1:k-2,k-1:k).
!>
!>  If UPLO = 'L', then A = L*D*L**T, where
!>     L = P(1)*L(1)* ... *P(k)*L(k)* ...,
!>  i.e., L is a product of terms P(k)*L(k), where k increases from 1 to
!>  n in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
!>  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
!>  defined by IPIV(k), and L(k) is a unit lower triangular matrix, such
!>  that if the diagonal block D(k) is of order s (s = 1 or 2), then
!>
!>             (   I    0     0   )  k-1
!>     L(k) =  (   0    I     0   )  s
!>             (   0    v     I   )  n-k-s+1
!>                k-1   s  n-k-s+1
!>
!>  If s = 1, D(k) overwrites A(k,k), and v overwrites A(k+1:n,k).
!>  If s = 2, the lower triangle of D(k) overwrites A(k,k), A(k+1,k),
!>  and A(k+1,k+1), and v overwrites A(k+2:n,k:k+1).
!> \endverbatim
!
!> \par Contributors:
!  ==================
!>
!> \verbatim
!>
!>  November 2013,     Igor Kozachenko,
!>                  Computer Science Division,
!>                  University of California, Berkeley
!>
!>  September 2007, Sven Hammarling, Nicholas J. Higham, Craig Lucas,
!>                  School of Mathematics,
!>                  University of Manchester
!>
!>  01-01-96 - Based on modifications by
!>    J. Lewis, Boeing Computer Services Company
!>    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville abd , USA
!> \endverbatim
!
!  =====================================================================
      SUBROUTINE SSYTF2_ROOK(Uplo,N,A,Lda,Ipiv,Info)
      IMPLICIT NONE
!*--SSYTF2_ROOK198
!
!  -- LAPACK computational routine (version 3.5.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2013
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER Info , Lda , N
!     ..
!     .. Array Arguments ..
      INTEGER Ipiv(*)
      REAL A(Lda,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
      REAL EIGHT , SEVTEN
      PARAMETER (EIGHT=8.0E+0,SEVTEN=17.0E+0)
!     ..
!     .. Local Scalars ..
      LOGICAL upper , done
      INTEGER i , imax , j , jmax , itemp , k , kk , kp , kstep , p , ii
      REAL absakk , alpha , colmax , d11 , d12 , d21 , d22 , rowmax ,   &
     &     stemp , t , wk , wkm1 , wkp1 , sfmin
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER ISAMAX
      REAL SLAMCH
      EXTERNAL LSAME , ISAMAX , SLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL SSCAL , SSWAP , SSYR , XERBLA
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
         CALL XERBLA('SSYTF2_ROOK',-Info)
         RETURN
      ENDIF
!
!     Initialize ALPHA for use in choosing pivot block size.
!
      alpha = (ONE+SQRT(SEVTEN))/EIGHT
!
!     Compute machine safe minimum
!
      sfmin = SLAMCH('S')
!
      IF ( upper ) THEN
!
!        Factorize A as U*D*U**T using the upper triangle of A
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
               imax = ISAMAX(k-1,A(1,k),1)
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
                        jmax = imax + ISAMAX(k-imax,A(imax,imax+1),Lda)
                        rowmax = ABS(A(imax,jmax))
                     ELSE
                        rowmax = ZERO
                     ENDIF
!
                     IF ( imax>1 ) THEN
                        itemp = ISAMAX(imax-1,A(1,imax),1)
                        stemp = ABS(A(itemp,imax))
                        IF ( stemp>rowmax ) THEN
                           rowmax = stemp
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
                  IF ( p>1 ) CALL SSWAP(p-1,A(1,k),1,A(1,p),1)
                  IF ( p<(k-1) ) CALL SSWAP(k-p-1,A(p+1,k),1,A(p,p+1),  &
     &                 Lda)
                  t = A(k,k)
                  A(k,k) = A(p,p)
                  A(p,p) = t
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
                  IF ( kp>1 ) CALL SSWAP(kp-1,A(1,kk),1,A(1,kp),1)
                  IF ( (kk>1) .AND. (kp<(kk-1)) )                       &
     &                 CALL SSWAP(kk-kp-1,A(kp+1,kk),1,A(kp,kp+1),Lda)
                  t = A(kk,kk)
                  A(kk,kk) = A(kp,kp)
                  A(kp,kp) = t
                  IF ( kstep==2 ) THEN
                     t = A(k-1,k)
                     A(k-1,k) = A(kp,k)
                     A(kp,k) = t
                  ENDIF
               ENDIF
!
!           Update the leading submatrix
!
               IF ( kstep==1 ) THEN
!
!              1-by-1 pivot block D(k): column k now holds
!
!              W(k) = U(k)*D(k)
!
!              where U(k) is the k-th column of U
!
                  IF ( k>1 ) THEN
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
                        CALL SSYR(Uplo,k-1,-d11,A(1,k),1,A,Lda)
!
!                    Store U(k) in column k
!
                        CALL SSCAL(k-1,d11,A(1,k),1)
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
                        CALL SSYR(Uplo,k-1,-d11,A(1,k),1,A,Lda)
                     ENDIF
                  ENDIF
!
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
               ELSEIF ( k>2 ) THEN
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
                        A(i,j) = A(i,j) - (A(i,k)/d12)                  &
     &                           *wk - (A(i,k-1)/d12)*wkm1
                     ENDDO
!
!                    Store U(k) and U(k-1) in cols k and k-1 for row J
!
                     A(j,k) = wk/d12
                     A(j,k-1) = wkm1/d12
!
                  ENDDO
!
!
               ENDIF
            ENDIF
!
!        Store details of the interchanges in IPIV
!
            IF ( kstep==1 ) THEN
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
      ELSE
!
!        Factorize A as L*D*L**T using the lower triangle of A
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
               imax = k + ISAMAX(N-k,A(k+1,k),1)
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
                        jmax = k - 1 + ISAMAX(imax-k,A(imax,k),Lda)
                        rowmax = ABS(A(imax,jmax))
                     ELSE
                        rowmax = ZERO
                     ENDIF
!
                     IF ( imax<N ) THEN
                        itemp = imax + ISAMAX(N-imax,A(imax+1,imax),1)
                        stemp = ABS(A(itemp,imax))
                        IF ( stemp>rowmax ) THEN
                           rowmax = stemp
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
                  IF ( p<N ) CALL SSWAP(N-p,A(p+1,k),1,A(p+1,p),1)
                  IF ( p>(k+1) ) CALL SSWAP(p-k-1,A(k+1,k),1,A(p,k+1),  &
     &                 Lda)
                  t = A(k,k)
                  A(k,k) = A(p,p)
                  A(p,p) = t
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
                  IF ( kp<N ) CALL SSWAP(N-kp,A(kp+1,kk),1,A(kp+1,kp),1)
                  IF ( (kk<N) .AND. (kp>(kk+1)) )                       &
     &                 CALL SSWAP(kp-kk-1,A(kk+1,kk),1,A(kp,kk+1),Lda)
                  t = A(kk,kk)
                  A(kk,kk) = A(kp,kp)
                  A(kp,kp) = t
                  IF ( kstep==2 ) THEN
                     t = A(k+1,k)
                     A(k+1,k) = A(kp,k)
                     A(kp,k) = t
                  ENDIF
               ENDIF
!
!           Update the trailing submatrix
!
               IF ( kstep==1 ) THEN
!
!              1-by-1 pivot block D(k): column k now holds
!
!              W(k) = L(k)*D(k)
!
!              where L(k) is the k-th column of L
!
                  IF ( k<N ) THEN
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
                        CALL SSYR(Uplo,N-k,-d11,A(k+1,k),1,A(k+1,k+1),  &
     &                            Lda)
!
!                    Store L(k) in column k
!
                        CALL SSCAL(N-k,d11,A(k+1,k),1)
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
                        CALL SSYR(Uplo,N-k,-d11,A(k+1,k),1,A(k+1,k+1),  &
     &                            Lda)
                     ENDIF
                  ENDIF
!
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
               ELSEIF ( k<N-1 ) THEN
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
                        A(i,j) = A(i,j) - (A(i,k)/d21)                  &
     &                           *wk - (A(i,k+1)/d21)*wkp1
                     ENDDO
!
!                    Store L(k) and L(k+1) in cols k and k+1 for row J
!
                     A(j,k) = wk/d21
                     A(j,k+1) = wkp1/d21
!
                  ENDDO
!
!
               ENDIF
            ENDIF
!
!        Store details of the interchanges in IPIV
!
            IF ( kstep==1 ) THEN
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
      ENDIF
!
!
!
!     End of SSYTF2_ROOK
!
      END SUBROUTINE SSYTF2_ROOK
