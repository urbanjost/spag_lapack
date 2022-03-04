!*==slasyf.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SLASYF computes a partial factorization of a real symmetric matrix using the Bunch-Kaufman diagonal pivoting method.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLASYF + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasyf.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasyf.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasyf.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLASYF( UPLO, N, NB, KB, A, LDA, IPIV, W, LDW, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, KB, LDA, LDW, N, NB
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       REAL               A( LDA, * ), W( LDW, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLASYF computes a partial factorization of a real symmetric matrix A
!> using the Bunch-Kaufman diagonal pivoting method. The partial
!> factorization has the form:
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
!> SLASYF is an auxiliary routine called by SSYTRF. It uses blocked code
!> (calling Level 3 BLAS) to update the submatrix A11 (if UPLO = 'U') or
!> A22 (if UPLO = 'L').
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
!>             If IPIV(k) = IPIV(k-1) < 0, then rows and columns
!>             k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
!>             is a 2-by-2 diagonal block.
!>
!>          If UPLO = 'L':
!>             Only the first KB elements of IPIV are set.
!>
!>             If IPIV(k) > 0, then rows and columns k and IPIV(k) were
!>             interchanged and D(k,k) is a 1-by-1 diagonal block.
!>
!>             If IPIV(k) = IPIV(k+1) < 0, then rows and columns
!>             k+1 and -IPIV(k) were interchanged and D(k:k+1,k:k+1)
!>             is a 2-by-2 diagonal block.
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
!> \ingroup realSYcomputational
!
!> \par Contributors:
!  ==================
!>
!> \verbatim
!>
!>  November 2013,  Igor Kozachenko,
!>                  Computer Science Division,
!>                  University of California, Berkeley
!> \endverbatim
!
!  =====================================================================
      SUBROUTINE SLASYF(Uplo,N,Nb,Kb,A,Lda,Ipiv,W,Ldw,Info)
      IMPLICIT NONE
!*--SLASYF180
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
      REAL A(Lda,*) , W(Ldw,*)
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
      INTEGER imax , j , jb , jj , jmax , jp , k , kk , kkw , kp ,      &
     &        kstep , kw
      REAL absakk , alpha , colmax , d11 , d21 , d22 , r1 , rowmax , t
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER ISAMAX
      EXTERNAL LSAME , ISAMAX
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
      IF ( LSAME(Uplo,'U') ) THEN
!
!        Factorize the trailing columns of A using the upper triangle
!        of A and working backwards, and compute the matrix W = U12*D
!        for use in updating A11
!
!        K is the main loop index, decreasing from N in steps of 1 or 2
!
!        KW is the column of W which corresponds to column K of A
!
         k = N
         DO
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
                  CALL SGEMM('No transpose','Transpose',j-1,jb,N-k,-ONE,&
     &                       A(1,k+1),Lda,W(j,kw+1),Ldw,ONE,A(1,j),Lda)
               ENDDO
!
!        Put U12 in standard form by partially undoing the interchanges
!        in columns k+1:n looping backwards from k+1 to n
!
               j = k + 1
               DO
!
!           Undo the interchanges (if any) of rows JJ and JP at each
!           step J
!
!           (Here, J is a diagonal index)
                  jj = j
                  jp = Ipiv(j)
                  IF ( jp<0 ) THEN
                     jp = -jp
!              (Here, J is a diagonal index)
                     j = j + 1
                  ENDIF
!           (NOTE: Here, J is used to determine row length. Length N-J+1
!           of the rows to swap back doesn't include diagonal element)
                  j = j + 1
                  IF ( jp/=jj .AND. j<=N )                              &
     &                 CALL SSWAP(N-j+1,A(jp,j),Lda,A(jj,j),Lda)
                  IF ( j>=N ) THEN
!
!        Set KB to the number of columns factorized
!
                     Kb = N - k
                     GOTO 99999
                  ENDIF
               ENDDO
            ELSE
!
!        Copy column K of A to column KW of W and update it
!
               CALL SCOPY(k,A(1,k),1,W(1,kw),1)
               IF ( k<N ) CALL SGEMV('No transpose',k,N-k,-ONE,A(1,k+1),&
     &                               Lda,W(k,kw+1),Ldw,ONE,W(1,kw),1)
!
               kstep = 1
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
               ELSE
                  IF ( absakk>=alpha*colmax ) THEN
!
!              no interchange, use 1-by-1 pivot block
!
                     kp = k
                  ELSE
!
!              Copy column IMAX to column KW-1 of W and update it
!
                     CALL SCOPY(imax,A(1,imax),1,W(1,kw-1),1)
                     CALL SCOPY(k-imax,A(imax,imax+1),Lda,W(imax+1,kw-1)&
     &                          ,1)
                     IF ( k<N ) CALL SGEMV('No transpose',k,N-k,-ONE,   &
     &                    A(1,k+1),Lda,W(imax,kw+1),Ldw,ONE,W(1,kw-1),1)
!
!              JMAX is the column-index of the largest off-diagonal
!              element in row IMAX, and ROWMAX is its absolute value
!
                     jmax = imax + ISAMAX(k-imax,W(imax+1,kw-1),1)
                     rowmax = ABS(W(jmax,kw-1))
                     IF ( imax>1 ) THEN
                        jmax = ISAMAX(imax-1,W(1,kw-1),1)
                        rowmax = MAX(rowmax,ABS(W(jmax,kw-1)))
                     ENDIF
!
                     IF ( absakk>=alpha*colmax*(colmax/rowmax) ) THEN
!
!                 no interchange, use 1-by-1 pivot block
!
                        kp = k
                     ELSEIF ( ABS(W(imax,kw-1))>=alpha*rowmax ) THEN
!
!                 interchange rows and columns K and IMAX, use 1-by-1
!                 pivot block
!
                        kp = imax
!
!                 copy column KW-1 of W to column KW of W
!
                        CALL SCOPY(k,W(1,kw-1),1,W(1,kw),1)
                     ELSE
!
!                 interchange rows and columns K-1 and IMAX, use 2-by-2
!                 pivot block
!
                        kp = imax
                        kstep = 2
                     ENDIF
                  ENDIF
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
                     A(kp,kp) = A(kk,kk)
                     CALL SCOPY(kk-1-kp,A(kp+1,kk),1,A(kp,kp+1),Lda)
                     IF ( kp>1 ) CALL SCOPY(kp-1,A(1,kk),1,A(1,kp),1)
!
!              Interchange rows KK and KP in last K+1 to N columns of A
!              (columns K (or K and K-1 for 2-by-2 pivot) of A will be
!              later overwritten). Interchange rows KK and KP
!              in last KKW to NB columns of W.
!
                     IF ( k<N ) CALL SSWAP(N-k,A(kk,k+1),Lda,A(kp,k+1), &
     &                    Lda)
                     CALL SSWAP(N-kk+1,W(kk,kkw),Ldw,W(kp,kkw),Ldw)
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
!              Store subdiag. elements of column U(k)
!              and 1-by-1 block D(k) in column k of A.
!              NOTE: Diagonal element U(k,k) is a UNIT element
!              and not stored.
!                 A(k,k) := D(k,k) = W(k,kw)
!                 A(1:k-1,k) := U(1:k-1,k) = W(1:k-1,kw)/D(k,k)
!
                     CALL SCOPY(k,W(1,kw),1,A(1,k),1)
                     r1 = ONE/A(k,k)
                     CALL SSCAL(k-1,r1,A(1,k),1)
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
!              Store U(1:k-2,k-1) and U(1:k-2,k) and 2-by-2
!              block D(k-1:k,k-1:k) in columns k-1 and k of A.
!              NOTE: 2-by-2 diagonal block U(k-1:k,k-1:k) is a UNIT
!              block and not stored.
!                 A(k-1:k,k-1:k) := D(k-1:k,k-1:k) = W(k-1:k,kw-1:kw)
!                 A(1:k-2,k-1:k) := U(1:k-2,k:k-1:k) =
!                 = W(1:k-2,kw-1:kw) * ( D(k-1:k,k-1:k)**(-1) )
!
                     IF ( k>2 ) THEN
!
!                 Compose the columns of the inverse of 2-by-2 pivot
!                 block D in the following way to reduce the number
!                 of FLOPS when we myltiply panel ( W(kw-1) W(kw) ) by
!                 this inverse
!
!                 D**(-1) = ( d11 d21 )**(-1) =
!                           ( d21 d22 )
!
!                 = 1/(d11*d22-d21**2) * ( ( d22 ) (-d21 ) ) =
!                                        ( (-d21 ) ( d11 ) )
!
!                 = 1/d21 * 1/((d11/d21)*(d22/d21)-1) *
!
!                   * ( ( d22/d21 ) (      -1 ) ) =
!                     ( (      -1 ) ( d11/d21 ) )
!
!                 = 1/d21 * 1/(D22*D11-1) * ( ( D11 ) (  -1 ) ) =
!                                           ( ( -1  ) ( D22 ) )
!
!                 = 1/d21 * T * ( ( D11 ) (  -1 ) )
!                               ( (  -1 ) ( D22 ) )
!
!                 = D21 * ( ( D11 ) (  -1 ) )
!                         ( (  -1 ) ( D22 ) )
!
                        d21 = W(k-1,kw)
                        d11 = W(k,kw)/d21
                        d22 = W(k-1,kw-1)/d21
                        t = ONE/(d11*d22-ONE)
                        d21 = t/d21
!
!                 Update elements in columns A(k-1) and A(k) as
!                 dot products of rows of ( W(kw-1) W(kw) ) and columns
!                 of D**(-1)
!
                        DO j = 1 , k - 2
                           A(j,k-1) = d21*(d11*W(j,kw-1)-W(j,kw))
                           A(j,k) = d21*(d22*W(j,kw)-W(j,kw-1))
                        ENDDO
                     ENDIF
!
!              Copy D(k) to A
!
                     A(k-1,k-1) = W(k-1,kw-1)
                     A(k-1,k) = W(k-1,kw)
                     A(k,k) = W(k,kw)
!
                  ENDIF
!
               ENDIF
!
!        Store details of the interchanges in IPIV
!
               IF ( kstep==1 ) THEN
                  Ipiv(k) = kp
               ELSE
                  Ipiv(k) = -kp
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
!        Copy column K of A to column K of W and update it
!
            CALL SCOPY(N-k+1,A(k,k),1,W(k,k),1)
            CALL SGEMV('No transpose',N-k+1,k-1,-ONE,A(k,1),Lda,W(k,1), &
     &                 Ldw,ONE,W(k,k),1)
!
            kstep = 1
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
            ELSE
               IF ( absakk>=alpha*colmax ) THEN
!
!              no interchange, use 1-by-1 pivot block
!
                  kp = k
               ELSE
!
!              Copy column IMAX to column K+1 of W and update it
!
                  CALL SCOPY(imax-k,A(imax,k),Lda,W(k,k+1),1)
                  CALL SCOPY(N-imax+1,A(imax,imax),1,W(imax,k+1),1)
                  CALL SGEMV('No transpose',N-k+1,k-1,-ONE,A(k,1),Lda,  &
     &                       W(imax,1),Ldw,ONE,W(k,k+1),1)
!
!              JMAX is the column-index of the largest off-diagonal
!              element in row IMAX, and ROWMAX is its absolute value
!
                  jmax = k - 1 + ISAMAX(imax-k,W(k,k+1),1)
                  rowmax = ABS(W(jmax,k+1))
                  IF ( imax<N ) THEN
                     jmax = imax + ISAMAX(N-imax,W(imax+1,k+1),1)
                     rowmax = MAX(rowmax,ABS(W(jmax,k+1)))
                  ENDIF
!
                  IF ( absakk>=alpha*colmax*(colmax/rowmax) ) THEN
!
!                 no interchange, use 1-by-1 pivot block
!
                     kp = k
                  ELSEIF ( ABS(W(imax,k+1))>=alpha*rowmax ) THEN
!
!                 interchange rows and columns K and IMAX, use 1-by-1
!                 pivot block
!
                     kp = imax
!
!                 copy column K+1 of W to column K of W
!
                     CALL SCOPY(N-k+1,W(k,k+1),1,W(k,k),1)
                  ELSE
!
!                 interchange rows and columns K+1 and IMAX, use 2-by-2
!                 pivot block
!
                     kp = imax
                     kstep = 2
                  ENDIF
               ENDIF
!
!           ============================================================
!
!           KK is the column of A where pivoting step stopped
!
               kk = k + kstep - 1
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
                  A(kp,kp) = A(kk,kk)
                  CALL SCOPY(kp-kk-1,A(kk+1,kk),1,A(kp,kk+1),Lda)
                  IF ( kp<N ) CALL SCOPY(N-kp,A(kp+1,kk),1,A(kp+1,kp),1)
!
!              Interchange rows KK and KP in first K-1 columns of A
!              (columns K (or K and K+1 for 2-by-2 pivot) of A will be
!              later overwritten). Interchange rows KK and KP
!              in first KK columns of W.
!
                  IF ( k>1 ) CALL SSWAP(k-1,A(kk,1),Lda,A(kp,1),Lda)
                  CALL SSWAP(kk,W(kk,1),Ldw,W(kp,1),Ldw)
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
!              Store subdiag. elements of column L(k)
!              and 1-by-1 block D(k) in column k of A.
!              (NOTE: Diagonal element L(k,k) is a UNIT element
!              and not stored)
!                 A(k,k) := D(k,k) = W(k,k)
!                 A(k+1:N,k) := L(k+1:N,k) = W(k+1:N,k)/D(k,k)
!
                  CALL SCOPY(N-k+1,W(k,k),1,A(k,k),1)
                  IF ( k<N ) THEN
                     r1 = ONE/A(k,k)
                     CALL SSCAL(N-k,r1,A(k+1,k),1)
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
!              Store L(k+2:N,k) and L(k+2:N,k+1) and 2-by-2
!              block D(k:k+1,k:k+1) in columns k and k+1 of A.
!              (NOTE: 2-by-2 diagonal block L(k:k+1,k:k+1) is a UNIT
!              block and not stored)
!                 A(k:k+1,k:k+1) := D(k:k+1,k:k+1) = W(k:k+1,k:k+1)
!                 A(k+2:N,k:k+1) := L(k+2:N,k:k+1) =
!                 = W(k+2:N,k:k+1) * ( D(k:k+1,k:k+1)**(-1) )
!
                  IF ( k<N-1 ) THEN
!
!                 Compose the columns of the inverse of 2-by-2 pivot
!                 block D in the following way to reduce the number
!                 of FLOPS when we myltiply panel ( W(k) W(k+1) ) by
!                 this inverse
!
!                 D**(-1) = ( d11 d21 )**(-1) =
!                           ( d21 d22 )
!
!                 = 1/(d11*d22-d21**2) * ( ( d22 ) (-d21 ) ) =
!                                        ( (-d21 ) ( d11 ) )
!
!                 = 1/d21 * 1/((d11/d21)*(d22/d21)-1) *
!
!                   * ( ( d22/d21 ) (      -1 ) ) =
!                     ( (      -1 ) ( d11/d21 ) )
!
!                 = 1/d21 * 1/(D22*D11-1) * ( ( D11 ) (  -1 ) ) =
!                                           ( ( -1  ) ( D22 ) )
!
!                 = 1/d21 * T * ( ( D11 ) (  -1 ) )
!                               ( (  -1 ) ( D22 ) )
!
!                 = D21 * ( ( D11 ) (  -1 ) )
!                         ( (  -1 ) ( D22 ) )
!
                     d21 = W(k+1,k)
                     d11 = W(k+1,k+1)/d21
                     d22 = W(k,k)/d21
                     t = ONE/(d11*d22-ONE)
                     d21 = t/d21
!
!                 Update elements in columns A(k) and A(k+1) as
!                 dot products of rows of ( W(k) W(k+1) ) and columns
!                 of D**(-1)
!
                     DO j = k + 2 , N
                        A(j,k) = d21*(d11*W(j,k)-W(j,k+1))
                        A(j,k+1) = d21*(d22*W(j,k+1)-W(j,k))
                     ENDDO
                  ENDIF
!
!              Copy D(k) to A
!
                  A(k,k) = W(k,k)
                  A(k+1,k) = W(k+1,k)
                  A(k+1,k+1) = W(k+1,k+1)
!
               ENDIF
!
            ENDIF
!
!        Store details of the interchanges in IPIV
!
            IF ( kstep==1 ) THEN
               Ipiv(k) = kp
            ELSE
               Ipiv(k) = -kp
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
!        Put L21 in standard form by partially undoing the interchanges
!        of rows in columns 1:k-1 looping backwards from k-1 to 1
!
         j = k - 1
         DO
!
!           Undo the interchanges (if any) of rows JJ and JP at each
!           step J
!
!           (Here, J is a diagonal index)
            jj = j
            jp = Ipiv(j)
            IF ( jp<0 ) THEN
               jp = -jp
!              (Here, J is a diagonal index)
               j = j - 1
            ENDIF
!           (NOTE: Here, J is used to determine row length. Length J
!           of the rows to swap back doesn't include diagonal element)
            j = j - 1
            IF ( jp/=jj .AND. j>=1 )                                    &
     &           CALL SSWAP(j,A(jp,1),Lda,A(jj,1),Lda)
            IF ( j<=1 ) THEN
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
!     End of SLASYF
!
99999 END SUBROUTINE SLASYF
