!*==csytf2.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
 
!> \brief \b CSYTF2 computes the factorization of a real symmetric indefinite matrix, using the diagonal pivoting method (unblocked algorithm).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CSYTF2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csytf2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csytf2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csytf2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CSYTF2( UPLO, N, A, LDA, IPIV, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDA, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX            A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CSYTF2 computes the factorization of a complex symmetric matrix A
!> using the Bunch-Kaufman diagonal pivoting method:
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
!>          A is COMPLEX array, dimension (LDA,N)
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
!>             If IPIV(k) > 0, then rows and columns k and IPIV(k) were
!>             interchanged and D(k,k) is a 1-by-1 diagonal block.
!>
!>             If IPIV(k) = IPIV(k-1) < 0, then rows and columns
!>             k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
!>             is a 2-by-2 diagonal block.
!>
!>          If UPLO = 'L':
!>             If IPIV(k) > 0, then rows and columns k and IPIV(k) were
!>             interchanged and D(k,k) is a 1-by-1 diagonal block.
!>
!>             If IPIV(k) = IPIV(k+1) < 0, then rows and columns
!>             k+1 and -IPIV(k) were interchanged and D(k:k+1,k:k+1)
!>             is a 2-by-2 diagonal block.
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
!> \date December 2016
!
!> \ingroup complexSYcomputational
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
!>  09-29-06 - patch from
!>    Bobby Cheng, MathWorks
!>
!>    Replace l.209 and l.377
!>         IF( MAX( ABSAKK, COLMAX ).EQ.ZERO ) THEN
!>    by
!>         IF( (MAX( ABSAKK, COLMAX ).EQ.ZERO) .OR. SISNAN(ABSAKK) ) THEN
!>
!>  1-96 - Based on modifications by J. Lewis, Boeing Computer Services
!>         Company
!> \endverbatim
!
!  =====================================================================
      SUBROUTINE CSYTF2(Uplo,N,A,Lda,Ipiv,Info)
      IMPLICIT NONE
!*--CSYTF2196
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
      COMPLEX A(Lda,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
      REAL EIGHT , SEVTEN
      PARAMETER (EIGHT=8.0E+0,SEVTEN=17.0E+0)
      COMPLEX CONE
      PARAMETER (CONE=(1.0E+0,0.0E+0))
!     ..
!     .. Local Scalars ..
      LOGICAL upper
      INTEGER i , imax , j , jmax , k , kk , kp , kstep
      REAL absakk , alpha , colmax , rowmax
      COMPLEX d11 , d12 , d21 , d22 , r1 , t , wk , wkm1 , wkp1 , z
!     ..
!     .. External Functions ..
      LOGICAL LSAME , SISNAN
      INTEGER ICAMAX
      EXTERNAL LSAME , ICAMAX , SISNAN
!     ..
!     .. External Subroutines ..
      EXTERNAL CSCAL , CSWAP , CSYR , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , AIMAG , MAX , REAL , SQRT
!     ..
!     .. Statement Functions ..
      REAL CABS1
!     ..
!     .. Statement Function definitions ..
      CABS1(z) = ABS(REAL(z)) + ABS(AIMAG(z))
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
         CALL XERBLA('CSYTF2',-Info)
         RETURN
      ENDIF
!
!     Initialize ALPHA for use in choosing pivot block size.
!
      alpha = (ONE+SQRT(SEVTEN))/EIGHT
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
!
!        Determine rows and columns to be interchanged and whether
!        a 1-by-1 or 2-by-2 pivot block will be used
!
            absakk = CABS1(A(k,k))
!
!        IMAX is the row-index of the largest off-diagonal element in
!        column K, and COLMAX is its absolute value.
!        Determine both COLMAX and IMAX.
!
            IF ( k>1 ) THEN
               imax = ICAMAX(k-1,A(1,k),1)
               colmax = CABS1(A(imax,k))
            ELSE
               colmax = ZERO
            ENDIF
!
            IF ( MAX(absakk,colmax)==ZERO .OR. SISNAN(absakk) ) THEN
!
!           Column K is zero or underflow, or contains a NaN:
!           set INFO and continue
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
!              JMAX is the column-index of the largest off-diagonal
!              element in row IMAX, and ROWMAX is its absolute value
!
                  jmax = imax + ICAMAX(k-imax,A(imax,imax+1),Lda)
                  rowmax = CABS1(A(imax,jmax))
                  IF ( imax>1 ) THEN
                     jmax = ICAMAX(imax-1,A(1,imax),1)
                     rowmax = MAX(rowmax,CABS1(A(jmax,imax)))
                  ENDIF
!
                  IF ( absakk>=alpha*colmax*(colmax/rowmax) ) THEN
!
!                 no interchange, use 1-by-1 pivot block
!
                     kp = k
                  ELSEIF ( CABS1(A(imax,imax))>=alpha*rowmax ) THEN
!
!                 interchange rows and columns K and IMAX, use 1-by-1
!                 pivot block
!
                     kp = imax
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
               kk = k - kstep + 1
               IF ( kp/=kk ) THEN
!
!              Interchange rows and columns KK and KP in the leading
!              submatrix A(1:k,1:k)
!
                  CALL CSWAP(kp-1,A(1,kk),1,A(1,kp),1)
                  CALL CSWAP(kk-kp-1,A(kp+1,kk),1,A(kp,kp+1),Lda)
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
!              Perform a rank-1 update of A(1:k-1,1:k-1) as
!
!              A := A - U(k)*D(k)*U(k)**T = A - W(k)*1/D(k)*W(k)**T
!
                  r1 = CONE/A(k,k)
                  CALL CSYR(Uplo,k-1,-r1,A(1,k),1,A,Lda)
!
!              Store U(k) in column k
!
                  CALL CSCAL(k-1,r1,A(1,k),1)
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
!                 = A - ( W(k-1) W(k) )*inv(D(k))*( W(k-1) W(k) )**T
!
               ELSEIF ( k>2 ) THEN
!
                  d12 = A(k-1,k)
                  d22 = A(k-1,k-1)/d12
                  d11 = A(k,k)/d12
                  t = CONE/(d11*d22-CONE)
                  d12 = t/d12
!
                  DO j = k - 2 , 1 , -1
                     wkm1 = d12*(d11*A(j,k-1)-A(j,k))
                     wk = d12*(d22*A(j,k)-A(j,k-1))
                     DO i = j , 1 , -1
                        A(i,j) = A(i,j) - A(i,k)*wk - A(i,k-1)*wkm1
                     ENDDO
                     A(j,k) = wk
                     A(j,k-1) = wkm1
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
               Ipiv(k) = -kp
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
!
!        Determine rows and columns to be interchanged and whether
!        a 1-by-1 or 2-by-2 pivot block will be used
!
            absakk = CABS1(A(k,k))
!
!        IMAX is the row-index of the largest off-diagonal element in
!        column K, and COLMAX is its absolute value.
!        Determine both COLMAX and IMAX.
!
            IF ( k<N ) THEN
               imax = k + ICAMAX(N-k,A(k+1,k),1)
               colmax = CABS1(A(imax,k))
            ELSE
               colmax = ZERO
            ENDIF
!
            IF ( MAX(absakk,colmax)==ZERO .OR. SISNAN(absakk) ) THEN
!
!           Column K is zero or underflow, or contains a NaN:
!           set INFO and continue
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
!              JMAX is the column-index of the largest off-diagonal
!              element in row IMAX, and ROWMAX is its absolute value
!
                  jmax = k - 1 + ICAMAX(imax-k,A(imax,k),Lda)
                  rowmax = CABS1(A(imax,jmax))
                  IF ( imax<N ) THEN
                     jmax = imax + ICAMAX(N-imax,A(imax+1,imax),1)
                     rowmax = MAX(rowmax,CABS1(A(jmax,imax)))
                  ENDIF
!
                  IF ( absakk>=alpha*colmax*(colmax/rowmax) ) THEN
!
!                 no interchange, use 1-by-1 pivot block
!
                     kp = k
                  ELSEIF ( CABS1(A(imax,imax))>=alpha*rowmax ) THEN
!
!                 interchange rows and columns K and IMAX, use 1-by-1
!                 pivot block
!
                     kp = imax
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
               kk = k + kstep - 1
               IF ( kp/=kk ) THEN
!
!              Interchange rows and columns KK and KP in the trailing
!              submatrix A(k:n,k:n)
!
                  IF ( kp<N ) CALL CSWAP(N-kp,A(kp+1,kk),1,A(kp+1,kp),1)
                  CALL CSWAP(kp-kk-1,A(kk+1,kk),1,A(kp,kk+1),Lda)
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
!                 Perform a rank-1 update of A(k+1:n,k+1:n) as
!
!                 A := A - L(k)*D(k)*L(k)**T = A - W(k)*(1/D(k))*W(k)**T
!
                     r1 = CONE/A(k,k)
                     CALL CSYR(Uplo,N-k,-r1,A(k+1,k),1,A(k+1,k+1),Lda)
!
!                 Store L(k) in column K
!
                     CALL CSCAL(N-k,r1,A(k+1,k),1)
                  ENDIF
!
!              2-by-2 pivot block D(k)
!
               ELSEIF ( k<N-1 ) THEN
!
!                 Perform a rank-2 update of A(k+2:n,k+2:n) as
!
!                 A := A - ( L(k) L(k+1) )*D(k)*( L(k) L(k+1) )**T
!                    = A - ( W(k) W(k+1) )*inv(D(k))*( W(k) W(k+1) )**T
!
!                 where L(k) and L(k+1) are the k-th and (k+1)-th
!                 columns of L
!
                  d21 = A(k+1,k)
                  d11 = A(k+1,k+1)/d21
                  d22 = A(k,k)/d21
                  t = CONE/(d11*d22-CONE)
                  d21 = t/d21
!
                  DO j = k + 2 , N
                     wk = d21*(d11*A(j,k)-A(j,k+1))
                     wkp1 = d21*(d22*A(j,k+1)-A(j,k))
                     DO i = j , N
                        A(i,j) = A(i,j) - A(i,k)*wk - A(i,k+1)*wkp1
                     ENDDO
                     A(j,k) = wk
                     A(j,k+1) = wkp1
                  ENDDO
               ENDIF
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
      ENDIF
!
!
!     End of CSYTF2
!
      END SUBROUTINE CSYTF2
