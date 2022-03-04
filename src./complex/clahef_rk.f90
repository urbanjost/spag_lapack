!*==clahef_rk.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CLAHEF_RK computes a partial factorization of a complex Hermitian indefinite matrix using bounded Bunch-Kaufman (rook) diagonal pivoting method.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLAHEF_RK + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clahef_rk.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clahef_rk.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clahef_rk.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLAHEF_RK( UPLO, N, NB, KB, A, LDA, E, IPIV, W, LDW,
!                             INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, KB, LDA, LDW, N, NB
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX            A( LDA, * ), E( * ), W( LDW, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!> CLAHEF_RK computes a partial factorization of a complex Hermitian
!> matrix A using the bounded Bunch-Kaufman (rook) diagonal
!> pivoting method. The partial factorization has the form:
!>
!> A  =  ( I  U12 ) ( A11  0  ) (  I       0    )  if UPLO = 'U', or:
!>       ( 0  U22 ) (  0   D  ) ( U12**H U22**H )
!>
!> A  =  ( L11  0 ) (  D   0  ) ( L11**H L21**H )  if UPLO = 'L',
!>       ( L21  I ) (  0  A22 ) (  0       I    )
!>
!> where the order of D is at most NB. The actual order is returned in
!> the argument KB, and is either NB or NB-1, or N if N <= NB.
!>
!> CLAHEF_RK is an auxiliary routine called by CHETRF_RK. It uses
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
!>          A is COMPLEX array, dimension (LDA,N)
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
!>          E is COMPLEX array, dimension (N)
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
!>          W is COMPLEX array, dimension (LDW,NB)
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
!> \ingroup complexHEcomputational
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
      SUBROUTINE CLAHEF_RK(Uplo,N,Nb,Kb,A,Lda,E,Ipiv,W,Ldw,Info)
      USE S_CCOPY
      USE S_CGEMM
      USE S_CGEMV
      USE S_CLACGV
      USE S_CSSCAL
      USE S_CSWAP
      USE S_ICAMAX
      USE S_LSAME
      USE S_SLAMCH
      IMPLICIT NONE
!*--CLAHEF_RK274
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 ,              &
     &                      EIGHT = 8.0E+0 , SEVTEN = 17.0E+0
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0) ,                 &
     &                         CZERO = (0.0E+0,0.0E+0)
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Nb
      INTEGER , INTENT(OUT) :: Kb
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , INTENT(OUT) , DIMENSION(*) :: E
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldw,*) :: W
      INTEGER :: Ldw
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL :: absakk , alpha , colmax , r1 , rowmax , sfmin , stemp , t
      REAL :: CABS1
      COMPLEX :: d11 , d21 , d22 , z
      LOGICAL :: done
      INTEGER :: ii , imax , itemp , j , jb , jj , jmax , k , kk , kkw ,&
     &           kp , kstep , kw , p
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
!     .. Statement Functions ..
!     ..
!     .. Statement Function definitions ..
      CABS1(z) = ABS(REAL(z)) + ABS(AIMAG(z))
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
!        for use in updating A11 (note that conjg(W) is actually stored)
!
!        Initialize the first entry of array E, where superdiagonal
!        elements of D are stored
!
         E(1) = CZERO
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
!        A11 := A11 - U12*D*U12**H = A11 - U12*W**H
!
!        computing blocks of NB columns at a time (note that conjg(W) is
!        actually stored)
!
               DO j = ((k-1)/Nb)*Nb + 1 , 1 , -Nb
                  jb = MIN(Nb,k-j+1)
!
!           Update the upper triangle of the diagonal block
!
                  DO jj = j , j + jb - 1
                     A(jj,jj) = REAL(A(jj,jj))
                     CALL CGEMV('No transpose',jj-j+1,N-k,-CONE,A(j,k+1)&
     &                          ,Lda,W(jj,kw+1),Ldw,CONE,A(j,jj),1)
                     A(jj,jj) = REAL(A(jj,jj))
                  ENDDO
!
!           Update the rectangular superdiagonal block
!
                  IF ( j>=2 ) CALL CGEMM('No transpose','Transpose',j-1,&
     &                 jb,N-k,-CONE,A(1,k+1),Lda,W(j,kw+1),Ldw,CONE,    &
     &                 A(1,j),Lda)
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
               IF ( k>1 ) CALL CCOPY(k-1,A(1,k),1,W(1,kw),1)
               W(k,kw) = REAL(A(k,k))
               IF ( k<N ) THEN
                  CALL CGEMV('No transpose',k,N-k,-CONE,A(1,k+1),Lda,   &
     &                       W(k,kw+1),Ldw,CONE,W(1,kw),1)
                  W(k,kw) = REAL(W(k,kw))
               ENDIF
!
!        Determine rows and columns to be interchanged and whether
!        a 1-by-1 or 2-by-2 pivot block will be used
!
               absakk = ABS(REAL(W(k,kw)))
!
!        IMAX is the row-index of the largest off-diagonal element in
!        column K, and COLMAX is its absolute value.
!        Determine both COLMAX and IMAX.
!
               IF ( k>1 ) THEN
                  imax = ICAMAX(k-1,W(1,kw),1)
                  colmax = CABS1(W(imax,kw))
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
                  A(k,k) = REAL(W(k,kw))
                  IF ( k>1 ) CALL CCOPY(k-1,W(1,kw),1,A(1,k),1)
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
                  IF ( absakk>=alpha*colmax ) THEN
!
!              no interchange, use 1-by-1 pivot block
!
                     kp = k
!
                  ELSE
!
!              Lop until pivot found
!
                     done = .FALSE.
                     DO
!
!
!                 BEGIN pivot search loop body
!
!
!                 Copy column IMAX to column KW-1 of W and update it
!
                        IF ( imax>1 )                                   &
     &                       CALL CCOPY(imax-1,A(1,imax),1,W(1,kw-1),1)
                        W(imax,kw-1) = REAL(A(imax,imax))
!
                        CALL CCOPY(k-imax,A(imax,imax+1),Lda,           &
     &                             W(imax+1,kw-1),1)
                        CALL CLACGV(k-imax,W(imax+1,kw-1),1)
!
                        IF ( k<N ) THEN
                           CALL CGEMV('No transpose',k,N-k,-CONE,       &
     &                                A(1,k+1),Lda,W(imax,kw+1),Ldw,    &
     &                                CONE,W(1,kw-1),1)
                           W(imax,kw-1) = REAL(W(imax,kw-1))
                        ENDIF
!
!                 JMAX is the column-index of the largest off-diagonal
!                 element in row IMAX, and ROWMAX is its absolute value.
!                 Determine both ROWMAX and JMAX.
!
                        IF ( imax/=k ) THEN
                           jmax = imax + ICAMAX(k-imax,W(imax+1,kw-1),1)
                           rowmax = CABS1(W(jmax,kw-1))
                        ELSE
                           rowmax = ZERO
                        ENDIF
!
                        IF ( imax>1 ) THEN
                           itemp = ICAMAX(imax-1,W(1,kw-1),1)
                           stemp = CABS1(W(itemp,kw-1))
                           IF ( stemp>rowmax ) THEN
                              rowmax = stemp
                              jmax = itemp
                           ENDIF
                        ENDIF
!
!                 Case(2)
!                 Equivalent to testing for
!                 ABS( REAL( W( IMAX,KW-1 ) ) ).GE.ALPHA*ROWMAX
!                 (used to handle NaN and Inf)
!
                        IF ( ABS(REAL(W(imax,kw-1)))>=alpha*rowmax )    &
     &                       THEN
!
!                    interchange rows and columns K and IMAX,
!                    use 1-by-1 pivot block
!
                           kp = imax
!
!                    copy column KW-1 of W to column KW of W
!
                           CALL CCOPY(k,W(1,kw-1),1,W(1,kw),1)
!
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
!
!                    Copy updated JMAXth (next IMAXth) column to Kth of W
!
                           CALL CCOPY(k,W(1,kw-1),1,W(1,kw),1)
!
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
                  kk = k - kstep + 1
!
!           KKW is the column of W which corresponds to column KK of A
!
                  kkw = Nb + kk - N
!
!           Interchange rows and columns P and K.
!           Updated column P is already stored in column KW of W.
!
                  IF ( (kstep==2) .AND. (p/=k) ) THEN
!
!              Copy non-updated column K to column P of submatrix A
!              at step K. No need to copy element into columns
!              K and K-1 of A for 2-by-2 pivot, since these columns
!              will be later overwritten.
!
                     A(p,p) = REAL(A(k,k))
                     CALL CCOPY(k-1-p,A(p+1,k),1,A(p,p+1),Lda)
                     CALL CLACGV(k-1-p,A(p,p+1),Lda)
                     IF ( p>1 ) CALL CCOPY(p-1,A(1,k),1,A(1,p),1)
!
!              Interchange rows K and P in the last K+1 to N columns of A
!              (columns K and K-1 of A for 2-by-2 pivot will be
!              later overwritten). Interchange rows K and P
!              in last KKW to NB columns of W.
!
                     IF ( k<N ) CALL CSWAP(N-k,A(k,k+1),Lda,A(p,k+1),   &
     &                    Lda)
                     CALL CSWAP(N-kk+1,W(k,kkw),Ldw,W(p,kkw),Ldw)
                  ENDIF
!
!           Interchange rows and columns KP and KK.
!           Updated column KP is already stored in column KKW of W.
!
                  IF ( kp/=kk ) THEN
!
!              Copy non-updated column KK to column KP of submatrix A
!              at step K. No need to copy element into column K
!              (or K and K-1 for 2-by-2 pivot) of A, since these columns
!              will be later overwritten.
!
                     A(kp,kp) = REAL(A(kk,kk))
                     CALL CCOPY(kk-1-kp,A(kp+1,kk),1,A(kp,kp+1),Lda)
                     CALL CLACGV(kk-1-kp,A(kp,kp+1),Lda)
                     IF ( kp>1 ) CALL CCOPY(kp-1,A(1,kk),1,A(1,kp),1)
!
!              Interchange rows KK and KP in last K+1 to N columns of A
!              (columns K (or K and K-1 for 2-by-2 pivot) of A will be
!              later overwritten). Interchange rows KK and KP
!              in last KKW to NB columns of W.
!
                     IF ( k<N ) CALL CSWAP(N-k,A(kk,k+1),Lda,A(kp,k+1), &
     &                    Lda)
                     CALL CSWAP(N-kk+1,W(kk,kkw),Ldw,W(kp,kkw),Ldw)
                  ENDIF
!
                  IF ( kstep==1 ) THEN
!
!              1-by-1 pivot block D(k): column kw of W now holds
!
!              W(kw) = U(k)*D(k),
!
!              where U(k) is the k-th column of U
!
!              (1) Store subdiag. elements of column U(k)
!              and 1-by-1 block D(k) in column k of A.
!              (NOTE: Diagonal element U(k,k) is a UNIT element
!              and not stored)
!                 A(k,k) := D(k,k) = W(k,kw)
!                 A(1:k-1,k) := U(1:k-1,k) = W(1:k-1,kw)/D(k,k)
!
!              (NOTE: No need to use for Hermitian matrix
!              A( K, K ) = REAL( W( K, K) ) to separately copy diagonal
!              element D(k,k) from W (potentially saves only one load))
                     CALL CCOPY(k,W(1,kw),1,A(1,k),1)
                     IF ( k>1 ) THEN
!
!                 (NOTE: No need to check if A(k,k) is NOT ZERO,
!                  since that was ensured earlier in pivot search:
!                  case A(k,k) = 0 falls into 2x2 pivot case(3))
!
!                 Handle division by a small number
!
                        t = REAL(A(k,k))
                        IF ( ABS(t)>=sfmin ) THEN
                           r1 = ONE/t
                           CALL CSSCAL(k-1,r1,A(1,k),1)
                        ELSE
                           DO ii = 1 , k - 1
                              A(ii,k) = A(ii,k)/t
                           ENDDO
                        ENDIF
!
!                 (2) Conjugate column W(kw)
!
                        CALL CLACGV(k-1,W(1,kw),1)
!
!                 Store the superdiagonal element of D in array E
!
                        E(k) = CZERO
!
                     ENDIF
!
                  ELSE
!
!              2-by-2 pivot block D(k): columns kw and kw-1 of W now hold
!
!              ( W(kw-1) W(kw) ) = ( U(k-1) U(k) )*D(k)
!
!              where U(k) and U(k-1) are the k-th and (k-1)-th columns
!              of U
!
!              (1) Store U(1:k-2,k-1) and U(1:k-2,k) and 2-by-2
!              block D(k-1:k,k-1:k) in columns k-1 and k of A.
!              (NOTE: 2-by-2 diagonal block U(k-1:k,k-1:k) is a UNIT
!              block and not stored)
!                 A(k-1:k,k-1:k) := D(k-1:k,k-1:k) = W(k-1:k,kw-1:kw)
!                 A(1:k-2,k-1:k) := U(1:k-2,k:k-1:k) =
!                 = W(1:k-2,kw-1:kw) * ( D(k-1:k,k-1:k)**(-1) )
!
                     IF ( k>2 ) THEN
!
!                 Factor out the columns of the inverse of 2-by-2 pivot
!                 block D, so that each column contains 1, to reduce the
!                 number of FLOPS when we multiply panel
!                 ( W(kw-1) W(kw) ) by this inverse, i.e. by D**(-1).
!
!                 D**(-1) = ( d11 cj(d21) )**(-1) =
!                           ( d21    d22 )
!
!                 = 1/(d11*d22-|d21|**2) * ( ( d22) (-cj(d21) ) ) =
!                                          ( (-d21) (     d11 ) )
!
!                 = 1/(|d21|**2) * 1/((d11/cj(d21))*(d22/d21)-1) *
!
!                   * ( d21*( d22/d21 ) conj(d21)*(           - 1 ) ) =
!                     (     (      -1 )           ( d11/conj(d21) ) )
!
!                 = 1/(|d21|**2) * 1/(D22*D11-1) *
!
!                   * ( d21*( D11 ) conj(d21)*(  -1 ) ) =
!                     (     (  -1 )           ( D22 ) )
!
!                 = (1/|d21|**2) * T * ( d21*( D11 ) conj(d21)*(  -1 ) ) =
!                                      (     (  -1 )           ( D22 ) )
!
!                 = ( (T/conj(d21))*( D11 ) (T/d21)*(  -1 ) ) =
!                   (               (  -1 )         ( D22 ) )
!
!                 Handle division by a small number. (NOTE: order of
!                 operations is important)
!
!                 = ( T*(( D11 )/conj(D21)) T*((  -1 )/D21 ) )
!                   (   ((  -1 )          )   (( D22 )     ) ),
!
!                 where D11 = d22/d21,
!                       D22 = d11/conj(d21),
!                       D21 = d21,
!                       T = 1/(D22*D11-1).
!
!                 (NOTE: No need to check for division by ZERO,
!                  since that was ensured earlier in pivot search:
!                  (a) d21 != 0 in 2x2 pivot case(4),
!                      since |d21| should be larger than |d11| and |d22|;
!                  (b) (D22*D11 - 1) != 0, since from (a),
!                      both |D11| < 1, |D22| < 1, hence |D22*D11| << 1.)
!
                        d21 = W(k-1,kw)
                        d11 = W(k,kw)/CONJG(d21)
                        d22 = W(k-1,kw-1)/d21
                        t = ONE/(REAL(d11*d22)-ONE)
!
!                 Update elements in columns A(k-1) and A(k) as
!                 dot products of rows of ( W(kw-1) W(kw) ) and columns
!                 of D**(-1)
!
                        DO j = 1 , k - 2
                           A(j,k-1) = t*((d11*W(j,kw-1)-W(j,kw))/d21)
                           A(j,k) = t*((d22*W(j,kw)-W(j,kw-1))/CONJG(d21&
     &                              ))
                        ENDDO
                     ENDIF
!
!              Copy diagonal elements of D(K) to A,
!              copy superdiagonal element of D(K) to E(K) and
!              ZERO out superdiagonal entry of A
!
                     A(k-1,k-1) = W(k-1,kw-1)
                     A(k-1,k) = CZERO
                     A(k,k) = W(k,kw)
                     E(k) = W(k-1,kw)
                     E(k-1) = CZERO
!
!              (2) Conjugate columns W(kw) and W(kw-1)
!
                     CALL CLACGV(k-1,W(1,kw),1)
                     CALL CLACGV(k-2,W(1,kw-1),1)
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
!        for use in updating A22 (note that conjg(W) is actually stored)
!
!        Initialize the unused last entry of the subdiagonal array E.
!
         E(N) = CZERO
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
!        Copy column K of A to column K of W and update column K of W
!
            W(k,k) = REAL(A(k,k))
            IF ( k<N ) CALL CCOPY(N-k,A(k+1,k),1,W(k+1,k),1)
            IF ( k>1 ) THEN
               CALL CGEMV('No transpose',N-k+1,k-1,-CONE,A(k,1),Lda,    &
     &                    W(k,1),Ldw,CONE,W(k,k),1)
               W(k,k) = REAL(W(k,k))
            ENDIF
!
!        Determine rows and columns to be interchanged and whether
!        a 1-by-1 or 2-by-2 pivot block will be used
!
            absakk = ABS(REAL(W(k,k)))
!
!        IMAX is the row-index of the largest off-diagonal element in
!        column K, and COLMAX is its absolute value.
!        Determine both COLMAX and IMAX.
!
            IF ( k<N ) THEN
               imax = k + ICAMAX(N-k,W(k+1,k),1)
               colmax = CABS1(W(imax,k))
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
               A(k,k) = REAL(W(k,k))
               IF ( k<N ) CALL CCOPY(N-k,W(k+1,k),1,A(k+1,k),1)
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
!                 Copy column IMAX to column k+1 of W and update it
!
                     CALL CCOPY(imax-k,A(imax,k),Lda,W(k,k+1),1)
                     CALL CLACGV(imax-k,W(k,k+1),1)
                     W(imax,k+1) = REAL(A(imax,imax))
!
                     IF ( imax<N ) CALL CCOPY(N-imax,A(imax+1,imax),1,  &
     &                    W(imax+1,k+1),1)
!
                     IF ( k>1 ) THEN
                        CALL CGEMV('No transpose',N-k+1,k-1,-CONE,A(k,1)&
     &                             ,Lda,W(imax,1),Ldw,CONE,W(k,k+1),1)
                        W(imax,k+1) = REAL(W(imax,k+1))
                     ENDIF
!
!                 JMAX is the column-index of the largest off-diagonal
!                 element in row IMAX, and ROWMAX is its absolute value.
!                 Determine both ROWMAX and JMAX.
!
                     IF ( imax/=k ) THEN
                        jmax = k - 1 + ICAMAX(imax-k,W(k,k+1),1)
                        rowmax = CABS1(W(jmax,k+1))
                     ELSE
                        rowmax = ZERO
                     ENDIF
!
                     IF ( imax<N ) THEN
                        itemp = imax + ICAMAX(N-imax,W(imax+1,k+1),1)
                        stemp = CABS1(W(itemp,k+1))
                        IF ( stemp>rowmax ) THEN
                           rowmax = stemp
                           jmax = itemp
                        ENDIF
                     ENDIF
!
!                 Case(2)
!                 Equivalent to testing for
!                 ABS( REAL( W( IMAX,K+1 ) ) ).GE.ALPHA*ROWMAX
!                 (used to handle NaN and Inf)
!
                     IF ( ABS(REAL(W(imax,k+1)))>=alpha*rowmax ) THEN
!
!                    interchange rows and columns K and IMAX,
!                    use 1-by-1 pivot block
!
                        kp = imax
!
!                    copy column K+1 of W to column K of W
!
                        CALL CCOPY(N-k+1,W(k,k+1),1,W(k,k),1)
!
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
!
!                    Copy updated JMAXth (next IMAXth) column to Kth of W
!
                        CALL CCOPY(N-k+1,W(k,k+1),1,W(k,k),1)
!
                     ENDIF
!
!
!                 End pivot search loop body
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
!           Interchange rows and columns P and K (only for 2-by-2 pivot).
!           Updated column P is already stored in column K of W.
!
               IF ( (kstep==2) .AND. (p/=k) ) THEN
!
!              Copy non-updated column KK-1 to column P of submatrix A
!              at step K. No need to copy element into columns
!              K and K+1 of A for 2-by-2 pivot, since these columns
!              will be later overwritten.
!
                  A(p,p) = REAL(A(k,k))
                  CALL CCOPY(p-k-1,A(k+1,k),1,A(p,k+1),Lda)
                  CALL CLACGV(p-k-1,A(p,k+1),Lda)
                  IF ( p<N ) CALL CCOPY(N-p,A(p+1,k),1,A(p+1,p),1)
!
!              Interchange rows K and P in first K-1 columns of A
!              (columns K and K+1 of A for 2-by-2 pivot will be
!              later overwritten). Interchange rows K and P
!              in first KK columns of W.
!
                  IF ( k>1 ) CALL CSWAP(k-1,A(k,1),Lda,A(p,1),Lda)
                  CALL CSWAP(kk,W(k,1),Ldw,W(p,1),Ldw)
               ENDIF
!
!           Interchange rows and columns KP and KK.
!           Updated column KP is already stored in column KK of W.
!
               IF ( kp/=kk ) THEN
!
!              Copy non-updated column KK to column KP of submatrix A
!              at step K. No need to copy element into column K
!              (or K and K+1 for 2-by-2 pivot) of A, since these columns
!              will be later overwritten.
!
                  A(kp,kp) = REAL(A(kk,kk))
                  CALL CCOPY(kp-kk-1,A(kk+1,kk),1,A(kp,kk+1),Lda)
                  CALL CLACGV(kp-kk-1,A(kp,kk+1),Lda)
                  IF ( kp<N ) CALL CCOPY(N-kp,A(kp+1,kk),1,A(kp+1,kp),1)
!
!              Interchange rows KK and KP in first K-1 columns of A
!              (column K (or K and K+1 for 2-by-2 pivot) of A will be
!              later overwritten). Interchange rows KK and KP
!              in first KK columns of W.
!
                  IF ( k>1 ) CALL CSWAP(k-1,A(kk,1),Lda,A(kp,1),Lda)
                  CALL CSWAP(kk,W(kk,1),Ldw,W(kp,1),Ldw)
               ENDIF
!
               IF ( kstep==1 ) THEN
!
!              1-by-1 pivot block D(k): column k of W now holds
!
!              W(k) = L(k)*D(k),
!
!              where L(k) is the k-th column of L
!
!              (1) Store subdiag. elements of column L(k)
!              and 1-by-1 block D(k) in column k of A.
!              (NOTE: Diagonal element L(k,k) is a UNIT element
!              and not stored)
!                 A(k,k) := D(k,k) = W(k,k)
!                 A(k+1:N,k) := L(k+1:N,k) = W(k+1:N,k)/D(k,k)
!
!              (NOTE: No need to use for Hermitian matrix
!              A( K, K ) = REAL( W( K, K) ) to separately copy diagonal
!              element D(k,k) from W (potentially saves only one load))
                  CALL CCOPY(N-k+1,W(k,k),1,A(k,k),1)
                  IF ( k<N ) THEN
!
!                 (NOTE: No need to check if A(k,k) is NOT ZERO,
!                  since that was ensured earlier in pivot search:
!                  case A(k,k) = 0 falls into 2x2 pivot case(3))
!
!                 Handle division by a small number
!
                     t = REAL(A(k,k))
                     IF ( ABS(t)>=sfmin ) THEN
                        r1 = ONE/t
                        CALL CSSCAL(N-k,r1,A(k+1,k),1)
                     ELSE
                        DO ii = k + 1 , N
                           A(ii,k) = A(ii,k)/t
                        ENDDO
                     ENDIF
!
!                 (2) Conjugate column W(k)
!
                     CALL CLACGV(N-k,W(k+1,k),1)
!
!                 Store the subdiagonal element of D in array E
!
                     E(k) = CZERO
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
!              (1) Store L(k+2:N,k) and L(k+2:N,k+1) and 2-by-2
!              block D(k:k+1,k:k+1) in columns k and k+1 of A.
!              NOTE: 2-by-2 diagonal block L(k:k+1,k:k+1) is a UNIT
!              block and not stored.
!                 A(k:k+1,k:k+1) := D(k:k+1,k:k+1) = W(k:k+1,k:k+1)
!                 A(k+2:N,k:k+1) := L(k+2:N,k:k+1) =
!                 = W(k+2:N,k:k+1) * ( D(k:k+1,k:k+1)**(-1) )
!
                  IF ( k<N-1 ) THEN
!
!                 Factor out the columns of the inverse of 2-by-2 pivot
!                 block D, so that each column contains 1, to reduce the
!                 number of FLOPS when we multiply panel
!                 ( W(kw-1) W(kw) ) by this inverse, i.e. by D**(-1).
!
!                 D**(-1) = ( d11 cj(d21) )**(-1) =
!                           ( d21    d22 )
!
!                 = 1/(d11*d22-|d21|**2) * ( ( d22) (-cj(d21) ) ) =
!                                          ( (-d21) (     d11 ) )
!
!                 = 1/(|d21|**2) * 1/((d11/cj(d21))*(d22/d21)-1) *
!
!                   * ( d21*( d22/d21 ) conj(d21)*(           - 1 ) ) =
!                     (     (      -1 )           ( d11/conj(d21) ) )
!
!                 = 1/(|d21|**2) * 1/(D22*D11-1) *
!
!                   * ( d21*( D11 ) conj(d21)*(  -1 ) ) =
!                     (     (  -1 )           ( D22 ) )
!
!                 = (1/|d21|**2) * T * ( d21*( D11 ) conj(d21)*(  -1 ) ) =
!                                      (     (  -1 )           ( D22 ) )
!
!                 = ( (T/conj(d21))*( D11 ) (T/d21)*(  -1 ) ) =
!                   (               (  -1 )         ( D22 ) )
!
!                 Handle division by a small number. (NOTE: order of
!                 operations is important)
!
!                 = ( T*(( D11 )/conj(D21)) T*((  -1 )/D21 ) )
!                   (   ((  -1 )          )   (( D22 )     ) ),
!
!                 where D11 = d22/d21,
!                       D22 = d11/conj(d21),
!                       D21 = d21,
!                       T = 1/(D22*D11-1).
!
!                 (NOTE: No need to check for division by ZERO,
!                  since that was ensured earlier in pivot search:
!                  (a) d21 != 0 in 2x2 pivot case(4),
!                      since |d21| should be larger than |d11| and |d22|;
!                  (b) (D22*D11 - 1) != 0, since from (a),
!                      both |D11| < 1, |D22| < 1, hence |D22*D11| << 1.)
!
                     d21 = W(k+1,k)
                     d11 = W(k+1,k+1)/d21
                     d22 = W(k,k)/CONJG(d21)
                     t = ONE/(REAL(d11*d22)-ONE)
!
!                 Update elements in columns A(k) and A(k+1) as
!                 dot products of rows of ( W(k) W(k+1) ) and columns
!                 of D**(-1)
!
                     DO j = k + 2 , N
                        A(j,k) = t*((d11*W(j,k)-W(j,k+1))/CONJG(d21))
                        A(j,k+1) = t*((d22*W(j,k+1)-W(j,k))/d21)
                     ENDDO
                  ENDIF
!
!              Copy diagonal elements of D(K) to A,
!              copy subdiagonal element of D(K) to E(K) and
!              ZERO out subdiagonal entry of A
!
                  A(k,k) = W(k,k)
                  A(k+1,k) = CZERO
                  A(k+1,k+1) = W(k+1,k+1)
                  E(k) = W(k+1,k)
                  E(k+1) = CZERO
!
!              (2) Conjugate columns W(k) and W(k+1)
!
                  CALL CLACGV(N-k,W(k+1,k),1)
                  CALL CLACGV(N-k-1,W(k+2,k+1),1)
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
!        A22 := A22 - L21*D*L21**H = A22 - L21*W**H
!
!        computing blocks of NB columns at a time (note that conjg(W) is
!        actually stored)
!
         DO j = k , N , Nb
            jb = MIN(Nb,N-j+1)
!
!           Update the lower triangle of the diagonal block
!
            DO jj = j , j + jb - 1
               A(jj,jj) = REAL(A(jj,jj))
               CALL CGEMV('No transpose',j+jb-jj,k-1,-CONE,A(jj,1),Lda, &
     &                    W(jj,1),Ldw,CONE,A(jj,jj),1)
               A(jj,jj) = REAL(A(jj,jj))
            ENDDO
!
!           Update the rectangular subdiagonal block
!
            IF ( j+jb<=N ) CALL CGEMM('No transpose','Transpose',       &
     &                                N-j-jb+1,jb,k-1,-CONE,A(j+jb,1),  &
     &                                Lda,W(j,1),Ldw,CONE,A(j+jb,j),Lda)
         ENDDO
!
!        Set KB to the number of columns factorized
!
         Kb = k - 1
!
      ENDIF
!
!     End of CLAHEF_RK
!
      END SUBROUTINE CLAHEF_RK
