!*==dlasyf_rook.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b DLASYF_ROOK *> DLASYF_ROOK computes a partial factorization of a real symmetric matrix using the bounded Bunch-Kaufman ("rook") diagonal pivoting method.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLASYF_ROOK + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasyf_rook.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlasyf_rook.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlasyf_rook.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLASYF_ROOK( UPLO, N, NB, KB, A, LDA, IPIV, W, LDW, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, KB, LDA, LDW, N, NB
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       DOUBLE PRECISION   A( LDA, * ), W( LDW, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLASYF_ROOK computes a partial factorization of a real symmetric
!> matrix A using the bounded Bunch-Kaufman ("rook") diagonal
!> pivoting method. The partial factorization has the form:
!>
!> A  =  ( I  U12 ) ( A11  0  ) (  I       0    )  if UPLO = 'U', or:
!>       ( 0  U22 ) (  0   D  ) ( U12**T U22**T )
!>
!> A  =  ( L11  0 ) (  D   0  ) ( L11**T L21**T )  if UPLO = 'L'
!>       ( L21  I ) (  0  A22 ) (  0       I    )
!>
!> where the order of D is at most NB. The actual order is returned in
!> the argument KB, and is either NB or NB-1, or N if N <= NB.
!>
!> DLASYF_ROOK is an auxiliary routine called by DSYTRF_ROOK. It uses
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
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!>          n-by-n upper triangular part of A contains the upper
!>          triangular part of the matrix A, and the strictly lower
!>          triangular part of A is not referenced.  If UPLO = 'L', the
!>          leading n-by-n lower triangular part of A contains the lower
!>          triangular part of the matrix A, and the strictly upper
!>          triangular part of A is not referenced.
!>          On exit, A contains details of the partial factorization.
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
!>             Only the last KB elements of IPIV are set.
!>
!>             If IPIV(k) > 0, then rows and columns k and IPIV(k) were
!>             interchanged and D(k,k) is a 1-by-1 diagonal block.
!>
!>             If IPIV(k) < 0 and IPIV(k-1) < 0, then rows and
!>             columns k and -IPIV(k) were interchanged and rows and
!>             columns k-1 and -IPIV(k-1) were inerchaged,
!>             D(k-1:k,k-1:k) is a 2-by-2 diagonal block.
!>
!>          If UPLO = 'L':
!>             Only the first KB elements of IPIV are set.
!>
!>             If IPIV(k) > 0, then rows and columns k and IPIV(k)
!>             were interchanged and D(k,k) is a 1-by-1 diagonal block.
!>
!>             If IPIV(k) < 0 and IPIV(k+1) < 0, then rows and
!>             columns k and -IPIV(k) were interchanged and rows and
!>             columns k+1 and -IPIV(k+1) were inerchaged,
!>             D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
!> \endverbatim
!>
!> \param[out] W
!> \verbatim
!>          W is DOUBLE PRECISION array, dimension (LDW,NB)
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
!>          > 0: if INFO = k, D(k,k) is exactly zero.  The factorization
!>               has been completed, but the block diagonal matrix D is
!>               exactly singular.
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
!> \ingroup doubleSYcomputational
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
!> \endverbatim
!
!  =====================================================================
      SUBROUTINE DLASYF_ROOK(Uplo,N,Nb,Kb,A,Lda,Ipiv,W,Ldw,Info)
      IMPLICIT NONE
!*--DLASYF_ROOK187
!
!  -- LAPACK computational routine (version 3.5.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2013
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER Info , Kb , Lda , Ldw , N , Nb
!     ..
!     .. Array Arguments ..
      INTEGER Ipiv(*)
      DOUBLE PRECISION A(Lda,*) , W(Ldw,*)
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
      LOGICAL done
      INTEGER imax , itemp , j , jb , jj , jmax , jp1 , jp2 , k , kk ,  &
     &        kw , kkw , kp , kstep , p , ii
 
      DOUBLE PRECISION absakk , alpha , colmax , d11 , d12 , d21 , d22 ,&
     &                 dtemp , r1 , rowmax , t , sfmin
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER IDAMAX
      DOUBLE PRECISION DLAMCH
      EXTERNAL LSAME , IDAMAX , DLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL DCOPY , DGEMM , DGEMV , DSCAL , DSWAP
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
      sfmin = DLAMCH('S')
!
      IF ( LSAME(Uplo,'U') ) THEN
!
!        Factorize the trailing columns of A using the upper triangle
!        of A and working backwards, and compute the matrix W = U12*D
!        for use in updating A11
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
                     CALL DGEMV('No transpose',jj-j+1,N-k,-ONE,A(j,k+1),&
     &                          Lda,W(jj,kw+1),Ldw,ONE,A(j,jj),1)
                  ENDDO
!
!           Update the rectangular superdiagonal block
!
                  IF ( j>=2 ) CALL DGEMM('No transpose','Transpose',j-1,&
     &                 jb,N-k,-ONE,A(1,k+1),Lda,W(j,kw+1),Ldw,ONE,A(1,j)&
     &                 ,Lda)
               ENDDO
!
!        Put U12 in standard form by partially undoing the interchanges
!        in columns k+1:n
!
               j = k + 1
               DO
!
                  kstep = 1
                  jp1 = 1
                  jj = j
                  jp2 = Ipiv(j)
                  IF ( jp2<0 ) THEN
                     jp2 = -jp2
                     j = j + 1
                     jp1 = -Ipiv(j)
                     kstep = 2
                  ENDIF
!
                  j = j + 1
                  IF ( jp2/=jj .AND. j<=N )                             &
     &                 CALL DSWAP(N-j+1,A(jp2,j),Lda,A(jj,j),Lda)
                  jj = j - 1
                  IF ( jp1/=jj .AND. kstep==2 )                         &
     &                 CALL DSWAP(N-j+1,A(jp1,j),Lda,A(jj,j),Lda)
                  IF ( j>N ) THEN
!
!        Set KB to the number of columns factorized
!
                     Kb = N - k
                     GOTO 99999
                  ENDIF
               ENDDO
            ELSE
!
               kstep = 1
               p = k
!
!        Copy column K of A to column KW of W and update it
!
               CALL DCOPY(k,A(1,k),1,W(1,kw),1)
               IF ( k<N ) CALL DGEMV('No transpose',k,N-k,-ONE,A(1,k+1),&
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
                  imax = IDAMAX(k-1,W(1,kw),1)
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
                  CALL DCOPY(k,W(1,kw),1,A(1,k),1)
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
                        CALL DCOPY(imax,A(1,imax),1,W(1,kw-1),1)
                        CALL DCOPY(k-imax,A(imax,imax+1),Lda,           &
     &                             W(imax+1,kw-1),1)
!
                        IF ( k<N ) CALL DGEMV('No transpose',k,N-k,-ONE,&
     &                       A(1,k+1),Lda,W(imax,kw+1),Ldw,ONE,W(1,kw-1)&
     &                       ,1)
!
!                 JMAX is the column-index of the largest off-diagonal
!                 element in row IMAX, and ROWMAX is its absolute value.
!                 Determine both ROWMAX and JMAX.
!
                        IF ( imax/=k ) THEN
                           jmax = imax + IDAMAX(k-imax,W(imax+1,kw-1),1)
                           rowmax = ABS(W(jmax,kw-1))
                        ELSE
                           rowmax = ZERO
                        ENDIF
!
                        IF ( imax>1 ) THEN
                           itemp = IDAMAX(imax-1,W(1,kw-1),1)
                           dtemp = ABS(W(itemp,kw-1))
                           IF ( dtemp>rowmax ) THEN
                              rowmax = dtemp
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
                           CALL DCOPY(k,W(1,kw-1),1,W(1,kw),1)
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
                           CALL DCOPY(k,W(1,kw-1),1,W(1,kw),1)
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
                     CALL DCOPY(k-p,A(p+1,k),1,A(p,p+1),Lda)
                     CALL DCOPY(p,A(1,k),1,A(1,p),1)
!
!              Interchange rows K and P in last N-K+1 columns of A
!              and last N-K+2 columns of W
!
                     CALL DSWAP(N-k+1,A(k,k),Lda,A(p,k),Lda)
                     CALL DSWAP(N-kk+1,W(k,kkw),Ldw,W(p,kkw),Ldw)
                  ENDIF
!
!           Updated column KP is already stored in column KKW of W
!
                  IF ( kp/=kk ) THEN
!
!              Copy non-updated column KK to column KP
!
                     A(kp,k) = A(kk,k)
                     CALL DCOPY(k-1-kp,A(kp+1,kk),1,A(kp,kp+1),Lda)
                     CALL DCOPY(kp,A(1,kk),1,A(1,kp),1)
!
!              Interchange rows KK and KP in last N-KK+1 columns
!              of A and W
!
                     CALL DSWAP(N-kk+1,A(kk,kk),Lda,A(kp,kk),Lda)
                     CALL DSWAP(N-kk+1,W(kk,kkw),Ldw,W(kp,kkw),Ldw)
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
                     CALL DCOPY(k,W(1,kw),1,A(1,k),1)
                     IF ( k>1 ) THEN
                        IF ( ABS(A(k,k))>=sfmin ) THEN
                           r1 = ONE/A(k,k)
                           CALL DSCAL(k-1,r1,A(1,k),1)
                        ELSEIF ( A(k,k)/=ZERO ) THEN
                           DO ii = 1 , k - 1
                              A(ii,k) = A(ii,k)/A(k,k)
                           ENDDO
                        ENDIF
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
!              Copy D(k) to A
!
                     A(k-1,k-1) = W(k-1,kw-1)
                     A(k-1,k) = W(k-1,kw)
                     A(k,k) = W(k,kw)
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
            ENDIF
         ENDDO
!
      ELSE
!
!        Factorize the leading columns of A using the lower triangle
!        of A and working forwards, and compute the matrix W = L21*D
!        for use in updating A22
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
            CALL DCOPY(N-k+1,A(k,k),1,W(k,k),1)
            IF ( k>1 ) CALL DGEMV('No transpose',N-k+1,k-1,-ONE,A(k,1), &
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
               imax = k + IDAMAX(N-k,W(k+1,k),1)
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
               CALL DCOPY(N-k+1,W(k,k),1,A(k,k),1)
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
                     CALL DCOPY(imax-k,A(imax,k),Lda,W(k,k+1),1)
                     CALL DCOPY(N-imax+1,A(imax,imax),1,W(imax,k+1),1)
                     IF ( k>1 ) CALL DGEMV('No transpose',N-k+1,k-1,    &
     &                    -ONE,A(k,1),Lda,W(imax,1),Ldw,ONE,W(k,k+1),1)
!
!                 JMAX is the column-index of the largest off-diagonal
!                 element in row IMAX, and ROWMAX is its absolute value.
!                 Determine both ROWMAX and JMAX.
!
                     IF ( imax/=k ) THEN
                        jmax = k - 1 + IDAMAX(imax-k,W(k,k+1),1)
                        rowmax = ABS(W(jmax,k+1))
                     ELSE
                        rowmax = ZERO
                     ENDIF
!
                     IF ( imax<N ) THEN
                        itemp = imax + IDAMAX(N-imax,W(imax+1,k+1),1)
                        dtemp = ABS(W(itemp,k+1))
                        IF ( dtemp>rowmax ) THEN
                           rowmax = dtemp
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
                        CALL DCOPY(N-k+1,W(k,k+1),1,W(k,k),1)
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
                        CALL DCOPY(N-k+1,W(k,k+1),1,W(k,k),1)
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
                  CALL DCOPY(p-k,A(k,k),1,A(p,k),Lda)
                  CALL DCOPY(N-p+1,A(p,k),1,A(p,p),1)
!
!              Interchange rows K and P in first K columns of A
!              and first K+1 columns of W
!
                  CALL DSWAP(k,A(k,1),Lda,A(p,1),Lda)
                  CALL DSWAP(kk,W(k,1),Ldw,W(p,1),Ldw)
               ENDIF
!
!           Updated column KP is already stored in column KK of W
!
               IF ( kp/=kk ) THEN
!
!              Copy non-updated column KK to column KP
!
                  A(kp,k) = A(kk,k)
                  CALL DCOPY(kp-k-1,A(k+1,kk),1,A(kp,k+1),Lda)
                  CALL DCOPY(N-kp+1,A(kp,kk),1,A(kp,kp),1)
!
!              Interchange rows KK and KP in first KK columns of A and W
!
                  CALL DSWAP(kk,A(kk,1),Lda,A(kp,1),Lda)
                  CALL DSWAP(kk,W(kk,1),Ldw,W(kp,1),Ldw)
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
                  CALL DCOPY(N-k+1,W(k,k),1,A(k,k),1)
                  IF ( k<N ) THEN
                     IF ( ABS(A(k,k))>=sfmin ) THEN
                        r1 = ONE/A(k,k)
                        CALL DSCAL(N-k,r1,A(k+1,k),1)
                     ELSEIF ( A(k,k)/=ZERO ) THEN
                        DO ii = k + 1 , N
                           A(ii,k) = A(ii,k)/A(k,k)
                        ENDDO
                     ENDIF
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
!              Copy D(k) to A
!
                  A(k,k) = W(k,k)
                  A(k+1,k) = W(k+1,k)
                  A(k+1,k+1) = W(k+1,k+1)
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
               CALL DGEMV('No transpose',j+jb-jj,k-1,-ONE,A(jj,1),Lda,  &
     &                    W(jj,1),Ldw,ONE,A(jj,jj),1)
            ENDDO
!
!           Update the rectangular subdiagonal block
!
            IF ( j+jb<=N ) CALL DGEMM('No transpose','Transpose',       &
     &                                N-j-jb+1,jb,k-1,-ONE,A(j+jb,1),   &
     &                                Lda,W(j,1),Ldw,ONE,A(j+jb,j),Lda)
         ENDDO
!
!        Put L21 in standard form by partially undoing the interchanges
!        in columns 1:k-1
!
         j = k - 1
         DO
!
            kstep = 1
            jp1 = 1
            jj = j
            jp2 = Ipiv(j)
            IF ( jp2<0 ) THEN
               jp2 = -jp2
               j = j - 1
               jp1 = -Ipiv(j)
               kstep = 2
            ENDIF
!
            j = j - 1
            IF ( jp2/=jj .AND. j>=1 )                                   &
     &           CALL DSWAP(j,A(jp2,1),Lda,A(jj,1),Lda)
            jj = j + 1
            IF ( jp1/=jj .AND. kstep==2 )                               &
     &           CALL DSWAP(j,A(jp1,1),Lda,A(jj,1),Lda)
            IF ( j<1 ) THEN
!
!        Set KB to the number of columns factorized
!
               Kb = k - 1
               EXIT
            ENDIF
         ENDDO
!
      ENDIF
!
!     End of DLASYF_ROOK
!
99999 END SUBROUTINE DLASYF_ROOK
