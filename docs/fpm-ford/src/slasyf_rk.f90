!*==slasyf_rk.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SLASYF_RK computes a partial factorization of a real symmetric indefinite matrix using bounded Bunch-Kaufman (rook) diagonal pivoting method.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLASYF_RK + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasyf_rk.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasyf_rk.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasyf_rk.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLASYF_RK( UPLO, N, NB, KB, A, LDA, E, IPIV, W, LDW,
!                             INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, KB, LDA, LDW, N, NB
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       REAL               A( LDA, * ), E( * ), W( LDW, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!> SLASYF_RK computes a partial factorization of a real symmetric
!> matrix A using the bounded Bunch-Kaufman (rook) diagonal
!> pivoting method. The partial factorization has the form:
!>
!> A  =  ( I  U12 ) ( A11  0  ) (  I       0    )  if UPLO = 'U', or:
!>       ( 0  U22 ) (  0   D  ) ( U12**T U22**T )
!>
!> A  =  ( L11  0 ) (  D   0  ) ( L11**T L21**T )  if UPLO = 'L',
!>       ( L21  I ) (  0  A22 ) (  0       I    )
!>
!> where the order of D is at most NB. The actual order is returned in
!> the argument KB, and is either NB or NB-1, or N if N <= NB.
!>
!> SLASYF_RK is an auxiliary routine called by SSYTRF_RK. It uses
!> blocked code (calling Level 3 BLAS) to update the submatrix
!> A11 (if UPLO = 'U') or A22 (if UPLO = 'L').
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
!> \param[in] NB
!> \verbatim
!>          NB is INTEGER
!>          The maximum number of columns of the matrix A that should be
!>          factored.  NB should be at least 2 to allow for 2-by-2 pivot
!>          blocks.
!> \endverbatim
!>
!> \param[out] KB
!> \verbatim
!>          KB is INTEGER
!>          The number of columns of A that were actually factored.
!>          KB is either NB-1 or NB, or N if N <= NB.
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
!>          at each factorization step.
!>
!>          If UPLO = 'U',
!>          ( in factorization order, k decreases from N to 1 ):
!>            a) A single positive entry IPIV(k) > 0 means:
!>               D(k,k) is a 1-by-1 diagonal block.
!>               If IPIV(k) != k, rows and columns k and IPIV(k) were
!>               interchanged in the submatrix A(1:N,N-KB+1:N);
!>               If IPIV(k) = k, no interchange occurred.
!>
!>
!>            b) A pair of consecutive negative entries
!>               IPIV(k) < 0 and IPIV(k-1) < 0 means:
!>               D(k-1:k,k-1:k) is a 2-by-2 diagonal block.
!>               (NOTE: negative entries in IPIV appear ONLY in pairs).
!>               1) If -IPIV(k) != k, rows and columns
!>                  k and -IPIV(k) were interchanged
!>                  in the matrix A(1:N,N-KB+1:N).
!>                  If -IPIV(k) = k, no interchange occurred.
!>               2) If -IPIV(k-1) != k-1, rows and columns
!>                  k-1 and -IPIV(k-1) were interchanged
!>                  in the submatrix A(1:N,N-KB+1:N).
!>                  If -IPIV(k-1) = k-1, no interchange occurred.
!>
!>            c) In both cases a) and b) is always ABS( IPIV(k) ) <= k.
!>
!>            d) NOTE: Any entry IPIV(k) is always NONZERO on output.
!>
!>          If UPLO = 'L',
!>          ( in factorization order, k increases from 1 to N ):
!>            a) A single positive entry IPIV(k) > 0 means:
!>               D(k,k) is a 1-by-1 diagonal block.
!>               If IPIV(k) != k, rows and columns k and IPIV(k) were
!>               interchanged in the submatrix A(1:N,1:KB).
!>               If IPIV(k) = k, no interchange occurred.
!>
!>            b) A pair of consecutive negative entries
!>               IPIV(k) < 0 and IPIV(k+1) < 0 means:
!>               D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
!>               (NOTE: negative entries in IPIV appear ONLY in pairs).
!>               1) If -IPIV(k) != k, rows and columns
!>                  k and -IPIV(k) were interchanged
!>                  in the submatrix A(1:N,1:KB).
!>                  If -IPIV(k) = k, no interchange occurred.
!>               2) If -IPIV(k+1) != k+1, rows and columns
!>                  k-1 and -IPIV(k-1) were interchanged
!>                  in the submatrix A(1:N,1:KB).
!>                  If -IPIV(k+1) = k+1, no interchange occurred.
!>
!>            c) In both cases a) and b) is always ABS( IPIV(k) ) >= k.
!>
!>            d) NOTE: Any entry IPIV(k) is always NONZERO on output.
!> \endverbatim
!>
!> \param[out] W
!> \verbatim
!>          W is REAL array, dimension (LDW,NB)
!> \endverbatim
!>
!> \param[in] LDW
!> \verbatim
!>          LDW is INTEGER
!>          The leading dimension of the array W.  LDW >= max(1,N).
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
      SUBROUTINE SLASYF_RK(Uplo,N,Nb,Kb,A,Lda,E,Ipiv,W,Ldw,Info)
      IMPLICIT NONE
!*--SLASYF_RK265
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER Info , Kb , Lda , Ldw , N , Nb
!     ..
!     .. Array Arguments ..
      INTEGER Ipiv(*)
      REAL A(Lda,*) , E(*) , W(Ldw,*)
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
      LOGICAL done
      INTEGER imax , itemp , j , jb , jj , jmax , k , kk , kw , kkw ,   &
     &        kp , kstep , p , ii
      REAL absakk , alpha , colmax , d11 , d12 , d21 , d22 , stemp ,    &
     &     r1 , rowmax , t , sfmin
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER ISAMAX
      REAL SLAMCH
      EXTERNAL LSAME , ISAMAX , SLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL SCOPY , SGEMM , SGEMV , SSCAL , SSWAP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , MIN , SQRT
!     ..
!     .. Executable Statements ..
!
      Info = 0
!
!     Initialize ALPHA for use in choosing pivot block size.
!
      alpha = (ONE+SQRT(SEVTEN))/EIGHT
!
!     Compute machine safe minimum
!
      sfmin = SLAMCH('S')
!
      IF ( LSAME(Uplo,'U') ) THEN
!
!        Factorize the trailing columns of A using the upper triangle
!        of A and working backwards, and compute the matrix W = U12*D
!        for use in updating A11
!
!        Initialize the first entry of array E, where superdiagonal
!        elements of D are stored
!
         E(1) = ZERO
!
!        K is the main loop index, decreasing from N in steps of 1 or 2
!
         k = N
         DO
!
!        KW is the column of W which corresponds to column K of A
!
            kw = Nb + k - N
!
!        Exit from loop
!
            IF ( (k<=N-Nb+1 .AND. Nb<N) .OR. k<1 ) THEN
!
!
!        Update the upper triangle of A11 (= A(1:k,1:k)) as
!
!        A11 := A11 - U12*D*U12**T = A11 - U12*W**T
!
!        computing blocks of NB columns at a time
!
               DO j = ((k-1)/Nb)*Nb + 1 , 1 , -Nb
                  jb = MIN(Nb,k-j+1)
!
!           Update the upper triangle of the diagonal block
!
                  DO jj = j , j + jb - 1
                     CALL SGEMV('No transpose',jj-j+1,N-k,-ONE,A(j,k+1),&
     &                          Lda,W(jj,kw+1),Ldw,ONE,A(j,jj),1)
                  ENDDO
!
!           Update the rectangular superdiagonal block
!
                  IF ( j>=2 ) CALL SGEMM('No transpose','Transpose',j-1,&
     &                 jb,N-k,-ONE,A(1,k+1),Lda,W(j,kw+1),Ldw,ONE,A(1,j)&
     &                 ,Lda)
               ENDDO
!
!        Set KB to the number of columns factorized
!
               Kb = N - k
               EXIT
            ELSE
!
               kstep = 1
               p = k
!
!        Copy column K of A to column KW of W and update it
!
               CALL SCOPY(k,A(1,k),1,W(1,kw),1)
               IF ( k<N ) CALL SGEMV('No transpose',k,N-k,-ONE,A(1,k+1),&
     &                               Lda,W(k,kw+1),Ldw,ONE,W(1,kw),1)
!
!        Determine rows and columns to be interchanged and whether
!        a 1-by-1 or 2-by-2 pivot block will be used
!
               absakk = ABS(W(k,kw))
!
!        IMAX is the row-index of the largest off-diagonal element in
!        column K, and COLMAX is its absolute value.
!        Determine both COLMAX and IMAX.
!
               IF ( k>1 ) THEN
                  imax = ISAMAX(k-1,W(1,kw),1)
                  colmax = ABS(W(imax,kw))
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
                  CALL SCOPY(k,W(1,kw),1,A(1,k),1)
!
!           Set E( K ) to zero
!
                  IF ( k>1 ) E(k) = ZERO
!
               ELSE
!
!           ============================================================
!
!           Test for interchange
!
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
!                 Begin pivot search loop body
!
!
!                 Copy column IMAX to column KW-1 of W and update it
!
                        CALL SCOPY(imax,A(1,imax),1,W(1,kw-1),1)
                        CALL SCOPY(k-imax,A(imax,imax+1),Lda,           &
     &                             W(imax+1,kw-1),1)
!
                        IF ( k<N ) CALL SGEMV('No transpose',k,N-k,-ONE,&
     &                       A(1,k+1),Lda,W(imax,kw+1),Ldw,ONE,W(1,kw-1)&
     &                       ,1)
!
!                 JMAX is the column-index of the largest off-diagonal
!                 element in row IMAX, and ROWMAX is its absolute value.
!                 Determine both ROWMAX and JMAX.
!
                        IF ( imax/=k ) THEN
                           jmax = imax + ISAMAX(k-imax,W(imax+1,kw-1),1)
                           rowmax = ABS(W(jmax,kw-1))
                        ELSE
                           rowmax = ZERO
                        ENDIF
!
                        IF ( imax>1 ) THEN
                           itemp = ISAMAX(imax-1,W(1,kw-1),1)
                           stemp = ABS(W(itemp,kw-1))
                           IF ( stemp>rowmax ) THEN
                              rowmax = stemp
                              jmax = itemp
                           ENDIF
                        ENDIF
!
!                 Equivalent to testing for
!                 ABS( W( IMAX, KW-1 ) ).GE.ALPHA*ROWMAX
!                 (used to handle NaN and Inf)
!
                        IF ( ABS(W(imax,kw-1))>=alpha*rowmax ) THEN
!
!                    interchange rows and columns K and IMAX,
!                    use 1-by-1 pivot block
!
                           kp = imax
!
!                    copy column KW-1 of W to column KW of W
!
                           CALL SCOPY(k,W(1,kw-1),1,W(1,kw),1)
!
                           done = .TRUE.
!
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
                        ELSE
!
!                    Pivot not found: set params and repeat
!
                           p = imax
                           colmax = rowmax
                           imax = jmax
!
!                    Copy updated JMAXth (next IMAXth) column to Kth of W
!
                           CALL SCOPY(k,W(1,kw-1),1,W(1,kw),1)
!
                        ENDIF
!
!                 End pivot search loop body
!
                        IF ( done ) EXIT
                     ENDDO
!
                  ENDIF
!
!           ============================================================
!
                  kk = k - kstep + 1
!
!           KKW is the column of W which corresponds to column KK of A
!
                  kkw = Nb + kk - N
!
                  IF ( (kstep==2) .AND. (p/=k) ) THEN
!
!              Copy non-updated column K to column P
!
                     CALL SCOPY(k-p,A(p+1,k),1,A(p,p+1),Lda)
                     CALL SCOPY(p,A(1,k),1,A(1,p),1)
!
!              Interchange rows K and P in last N-K+1 columns of A
!              and last N-K+2 columns of W
!
                     CALL SSWAP(N-k+1,A(k,k),Lda,A(p,k),Lda)
                     CALL SSWAP(N-kk+1,W(k,kkw),Ldw,W(p,kkw),Ldw)
                  ENDIF
!
!           Updated column KP is already stored in column KKW of W
!
                  IF ( kp/=kk ) THEN
!
!              Copy non-updated column KK to column KP
!
                     A(kp,k) = A(kk,k)
                     CALL SCOPY(k-1-kp,A(kp+1,kk),1,A(kp,kp+1),Lda)
                     CALL SCOPY(kp,A(1,kk),1,A(1,kp),1)
!
!              Interchange rows KK and KP in last N-KK+1 columns
!              of A and W
!
                     CALL SSWAP(N-kk+1,A(kk,kk),Lda,A(kp,kk),Lda)
                     CALL SSWAP(N-kk+1,W(kk,kkw),Ldw,W(kp,kkw),Ldw)
                  ENDIF
!
                  IF ( kstep==1 ) THEN
!
!              1-by-1 pivot block D(k): column KW of W now holds
!
!              W(k) = U(k)*D(k)
!
!              where U(k) is the k-th column of U
!
!              Store U(k) in column k of A
!
                     CALL SCOPY(k,W(1,kw),1,A(1,k),1)
                     IF ( k>1 ) THEN
                        IF ( ABS(A(k,k))>=sfmin ) THEN
                           r1 = ONE/A(k,k)
                           CALL SSCAL(k-1,r1,A(1,k),1)
                        ELSEIF ( A(k,k)/=ZERO ) THEN
                           DO ii = 1 , k - 1
                              A(ii,k) = A(ii,k)/A(k,k)
                           ENDDO
                        ENDIF
!
!                 Store the superdiagonal element of D in array E
!
                        E(k) = ZERO
!
                     ENDIF
!
                  ELSE
!
!              2-by-2 pivot block D(k): columns KW and KW-1 of W now
!              hold
!
!              ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)
!
!              where U(k) and U(k-1) are the k-th and (k-1)-th columns
!              of U
!
                     IF ( k>2 ) THEN
!
!                 Store U(k) and U(k-1) in columns k and k-1 of A
!
                        d12 = W(k-1,kw)
                        d11 = W(k,kw)/d12
                        d22 = W(k-1,kw-1)/d12
                        t = ONE/(d11*d22-ONE)
                        DO j = 1 , k - 2
                           A(j,k-1) = t*((d11*W(j,kw-1)-W(j,kw))/d12)
                           A(j,k) = t*((d22*W(j,kw)-W(j,kw-1))/d12)
                        ENDDO
                     ENDIF
!
!              Copy diagonal elements of D(K) to A,
!              copy superdiagonal element of D(K) to E(K) and
!              ZERO out superdiagonal entry of A
!
                     A(k-1,k-1) = W(k-1,kw-1)
                     A(k-1,k) = ZERO
                     A(k,k) = W(k,kw)
                     E(k) = W(k-1,kw)
                     E(k-1) = ZERO
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
                  Ipiv(k) = kp
               ELSE
                  Ipiv(k) = -p
                  Ipiv(k-1) = -kp
               ENDIF
!
!        Decrease K and return to the start of the main loop
!
               k = k - kstep
            ENDIF
         ENDDO
!
      ELSE
!
!        Factorize the leading columns of A using the lower triangle
!        of A and working forwards, and compute the matrix W = L21*D
!        for use in updating A22
!
!        Initialize the unused last entry of the subdiagonal array E.
!
         E(N) = ZERO
!
!        K is the main loop index, increasing from 1 in steps of 1 or 2
!
         k = 1
!
!        Exit from loop
!
         DO WHILE ( .NOT.((k>=Nb .AND. Nb<N) .OR. k>N) )
!
            kstep = 1
            p = k
!
!        Copy column K of A to column K of W and update it
!
            CALL SCOPY(N-k+1,A(k,k),1,W(k,k),1)
            IF ( k>1 ) CALL SGEMV('No transpose',N-k+1,k-1,-ONE,A(k,1), &
     &                            Lda,W(k,1),Ldw,ONE,W(k,k),1)
!
!        Determine rows and columns to be interchanged and whether
!        a 1-by-1 or 2-by-2 pivot block will be used
!
            absakk = ABS(W(k,k))
!
!        IMAX is the row-index of the largest off-diagonal element in
!        column K, and COLMAX is its absolute value.
!        Determine both COLMAX and IMAX.
!
            IF ( k<N ) THEN
               imax = k + ISAMAX(N-k,W(k+1,k),1)
               colmax = ABS(W(imax,k))
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
               CALL SCOPY(N-k+1,W(k,k),1,A(k,k),1)
!
!           Set E( K ) to zero
!
               IF ( k<N ) E(k) = ZERO
!
            ELSE
!
!           ============================================================
!
!           Test for interchange
!
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
!                 Begin pivot search loop body
!
!
!                 Copy column IMAX to column K+1 of W and update it
!
                     CALL SCOPY(imax-k,A(imax,k),Lda,W(k,k+1),1)
                     CALL SCOPY(N-imax+1,A(imax,imax),1,W(imax,k+1),1)
                     IF ( k>1 ) CALL SGEMV('No transpose',N-k+1,k-1,    &
     &                    -ONE,A(k,1),Lda,W(imax,1),Ldw,ONE,W(k,k+1),1)
!
!                 JMAX is the column-index of the largest off-diagonal
!                 element in row IMAX, and ROWMAX is its absolute value.
!                 Determine both ROWMAX and JMAX.
!
                     IF ( imax/=k ) THEN
                        jmax = k - 1 + ISAMAX(imax-k,W(k,k+1),1)
                        rowmax = ABS(W(jmax,k+1))
                     ELSE
                        rowmax = ZERO
                     ENDIF
!
                     IF ( imax<N ) THEN
                        itemp = imax + ISAMAX(N-imax,W(imax+1,k+1),1)
                        stemp = ABS(W(itemp,k+1))
                        IF ( stemp>rowmax ) THEN
                           rowmax = stemp
                           jmax = itemp
                        ENDIF
                     ENDIF
!
!                 Equivalent to testing for
!                 ABS( W( IMAX, K+1 ) ).GE.ALPHA*ROWMAX
!                 (used to handle NaN and Inf)
!
                     IF ( ABS(W(imax,k+1))>=alpha*rowmax ) THEN
!
!                    interchange rows and columns K and IMAX,
!                    use 1-by-1 pivot block
!
                        kp = imax
!
!                    copy column K+1 of W to column K of W
!
                        CALL SCOPY(N-k+1,W(k,k+1),1,W(k,k),1)
!
                        done = .TRUE.
!
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
                     ELSE
!
!                    Pivot not found: set params and repeat
!
                        p = imax
                        colmax = rowmax
                        imax = jmax
!
!                    Copy updated JMAXth (next IMAXth) column to Kth of W
!
                        CALL SCOPY(N-k+1,W(k,k+1),1,W(k,k),1)
!
                     ENDIF
!
!                 End pivot search loop body
!
                     IF ( done ) EXIT
                  ENDDO
!
               ENDIF
!
!           ============================================================
!
               kk = k + kstep - 1
!
               IF ( (kstep==2) .AND. (p/=k) ) THEN
!
!              Copy non-updated column K to column P
!
                  CALL SCOPY(p-k,A(k,k),1,A(p,k),Lda)
                  CALL SCOPY(N-p+1,A(p,k),1,A(p,p),1)
!
!              Interchange rows K and P in first K columns of A
!              and first K+1 columns of W
!
                  CALL SSWAP(k,A(k,1),Lda,A(p,1),Lda)
                  CALL SSWAP(kk,W(k,1),Ldw,W(p,1),Ldw)
               ENDIF
!
!           Updated column KP is already stored in column KK of W
!
               IF ( kp/=kk ) THEN
!
!              Copy non-updated column KK to column KP
!
                  A(kp,k) = A(kk,k)
                  CALL SCOPY(kp-k-1,A(k+1,kk),1,A(kp,k+1),Lda)
                  CALL SCOPY(N-kp+1,A(kp,kk),1,A(kp,kp),1)
!
!              Interchange rows KK and KP in first KK columns of A and W
!
                  CALL SSWAP(kk,A(kk,1),Lda,A(kp,1),Lda)
                  CALL SSWAP(kk,W(kk,1),Ldw,W(kp,1),Ldw)
               ENDIF
!
               IF ( kstep==1 ) THEN
!
!              1-by-1 pivot block D(k): column k of W now holds
!
!              W(k) = L(k)*D(k)
!
!              where L(k) is the k-th column of L
!
!              Store L(k) in column k of A
!
                  CALL SCOPY(N-k+1,W(k,k),1,A(k,k),1)
                  IF ( k<N ) THEN
                     IF ( ABS(A(k,k))>=sfmin ) THEN
                        r1 = ONE/A(k,k)
                        CALL SSCAL(N-k,r1,A(k+1,k),1)
                     ELSEIF ( A(k,k)/=ZERO ) THEN
                        DO ii = k + 1 , N
                           A(ii,k) = A(ii,k)/A(k,k)
                        ENDDO
                     ENDIF
!
!                 Store the subdiagonal element of D in array E
!
                     E(k) = ZERO
!
                  ENDIF
!
               ELSE
!
!              2-by-2 pivot block D(k): columns k and k+1 of W now hold
!
!              ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k)
!
!              where L(k) and L(k+1) are the k-th and (k+1)-th columns
!              of L
!
                  IF ( k<N-1 ) THEN
!
!                 Store L(k) and L(k+1) in columns k and k+1 of A
!
                     d21 = W(k+1,k)
                     d11 = W(k+1,k+1)/d21
                     d22 = W(k,k)/d21
                     t = ONE/(d11*d22-ONE)
                     DO j = k + 2 , N
                        A(j,k) = t*((d11*W(j,k)-W(j,k+1))/d21)
                        A(j,k+1) = t*((d22*W(j,k+1)-W(j,k))/d21)
                     ENDDO
                  ENDIF
!
!              Copy diagonal elements of D(K) to A,
!              copy subdiagonal element of D(K) to E(K) and
!              ZERO out subdiagonal entry of A
!
                  A(k,k) = W(k,k)
                  A(k+1,k) = ZERO
                  A(k+1,k+1) = W(k+1,k+1)
                  E(k) = W(k+1,k)
                  E(k+1) = ZERO
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
!        Update the lower triangle of A22 (= A(k:n,k:n)) as
!
!        A22 := A22 - L21*D*L21**T = A22 - L21*W**T
!
!        computing blocks of NB columns at a time
!
         DO j = k , N , Nb
            jb = MIN(Nb,N-j+1)
!
!           Update the lower triangle of the diagonal block
!
            DO jj = j , j + jb - 1
               CALL SGEMV('No transpose',j+jb-jj,k-1,-ONE,A(jj,1),Lda,  &
     &                    W(jj,1),Ldw,ONE,A(jj,jj),1)
            ENDDO
!
!           Update the rectangular subdiagonal block
!
            IF ( j+jb<=N ) CALL SGEMM('No transpose','Transpose',       &
     &                                N-j-jb+1,jb,k-1,-ONE,A(j+jb,1),   &
     &                                Lda,W(j,1),Ldw,ONE,A(j+jb,j),Lda)
         ENDDO
!
!        Set KB to the number of columns factorized
!
         Kb = k - 1
!
      ENDIF
!
!
!     End of SLASYF_RK
!
      END SUBROUTINE SLASYF_RK
