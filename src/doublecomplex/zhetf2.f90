!*==zhetf2.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b ZHETF2 computes the factorization of a complex Hermitian matrix, using the diagonal pivoting method (unblocked algorithm, calling Level 2 BLAS).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZHETF2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhetf2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhetf2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhetf2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZHETF2( UPLO, N, A, LDA, IPIV, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDA, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX*16         A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZHETF2 computes the factorization of a complex Hermitian matrix A
!> using the Bunch-Kaufman diagonal pivoting method:
!>
!>    A = U*D*U**H  or  A = L*D*L**H
!>
!> where U (or L) is a product of permutation and unit upper (lower)
!> triangular matrices, U**H is the conjugate transpose of U, and D is
!> Hermitian and block diagonal with 1-by-1 and 2-by-2 diagonal blocks.
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
!>          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading
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
!> \date November 2013
!
!> \ingroup complex16HEcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  If UPLO = 'U', then A = U*D*U**H, where
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
!>  If UPLO = 'L', then A = L*D*L**H, where
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
!>  09-29-06 - patch from
!>    Bobby Cheng, MathWorks
!>
!>    Replace l.210 and l.393
!>         IF( MAX( ABSAKK, COLMAX ).EQ.ZERO ) THEN
!>    by
!>         IF( (MAX( ABSAKK, COLMAX ).EQ.ZERO) .OR. DISNAN(ABSAKK) ) THEN
!>
!>  01-01-96 - Based on modifications by
!>    J. Lewis, Boeing Computer Services Company
!>    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA
!> \endverbatim
!
!  =====================================================================
      SUBROUTINE ZHETF2(Uplo,N,A,Lda,Ipiv,Info)
      IMPLICIT NONE
!*--ZHETF2195
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
      COMPLEX*16 A(Lda,*)
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
      LOGICAL upper
      INTEGER i , imax , j , jmax , k , kk , kp , kstep
      DOUBLE PRECISION absakk , alpha , colmax , d , d11 , d22 , r1 ,   &
     &                 rowmax , tt
      COMPLEX*16 d12 , d21 , t , wk , wkm1 , wkp1 , zdum
!     ..
!     .. External Functions ..
      LOGICAL LSAME , DISNAN
      INTEGER IZAMAX
      DOUBLE PRECISION DLAPY2
      EXTERNAL LSAME , IZAMAX , DLAPY2 , DISNAN
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
      CABS1(zdum) = ABS(DBLE(zdum)) + ABS(DIMAG(zdum))
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
         CALL XERBLA('ZHETF2',-Info)
         RETURN
      ENDIF
!
!     Initialize ALPHA for use in choosing pivot block size.
!
      alpha = (ONE+SQRT(SEVTEN))/EIGHT
!
      IF ( upper ) THEN
!
!        Factorize A as U*D*U**H using the upper triangle of A
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
            IF ( (MAX(absakk,colmax)==ZERO) .OR. DISNAN(absakk) ) THEN
!
!           Column K is zero or underflow, or contains a NaN:
!           set INFO and continue
!
               IF ( Info==0 ) Info = k
               kp = k
               A(k,k) = DBLE(A(k,k))
            ELSE
!
!           ============================================================
!
!           Test for interchange
!
               IF ( absakk>=alpha*colmax ) THEN
!
!              no interchange, use 1-by-1 pivot block
!
                  kp = k
               ELSE
!
!              JMAX is the column-index of the largest off-diagonal
!              element in row IMAX, and ROWMAX is its absolute value.
!              Determine only ROWMAX.
!
                  jmax = imax + IZAMAX(k-imax,A(imax,imax+1),Lda)
                  rowmax = CABS1(A(imax,jmax))
                  IF ( imax>1 ) THEN
                     jmax = IZAMAX(imax-1,A(1,imax),1)
                     rowmax = MAX(rowmax,CABS1(A(jmax,imax)))
                  ENDIF
!
                  IF ( absakk>=alpha*colmax*(colmax/rowmax) ) THEN
!
!                 no interchange, use 1-by-1 pivot block
!
                     kp = k
!
                  ELSEIF ( ABS(DBLE(A(imax,imax)))>=alpha*rowmax ) THEN
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
!
               ENDIF
!
!           ============================================================
!
               kk = k - kstep + 1
               IF ( kp/=kk ) THEN
!
!              Interchange rows and columns KK and KP in the leading
!              submatrix A(1:k,1:k)
!
                  CALL ZSWAP(kp-1,A(1,kk),1,A(1,kp),1)
                  DO j = kp + 1 , kk - 1
                     t = DCONJG(A(j,kk))
                     A(j,kk) = DCONJG(A(kp,j))
                     A(kp,j) = t
                  ENDDO
                  A(kp,kk) = DCONJG(A(kp,kk))
                  r1 = DBLE(A(kk,kk))
                  A(kk,kk) = DBLE(A(kp,kp))
                  A(kp,kp) = r1
                  IF ( kstep==2 ) THEN
                     A(k,k) = DBLE(A(k,k))
                     t = A(k-1,k)
                     A(k-1,k) = A(kp,k)
                     A(kp,k) = t
                  ENDIF
               ELSE
                  A(k,k) = DBLE(A(k,k))
                  IF ( kstep==2 ) A(k-1,k-1) = DBLE(A(k-1,k-1))
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
!              A := A - U(k)*D(k)*U(k)**H = A - W(k)*1/D(k)*W(k)**H
!
                  r1 = ONE/DBLE(A(k,k))
                  CALL ZHER(Uplo,k-1,-r1,A(1,k),1,A,Lda)
!
!              Store U(k) in column k
!
                  CALL ZDSCAL(k-1,r1,A(1,k),1)
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
!              A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )**H
!                 = A - ( W(k-1) W(k) )*inv(D(k))*( W(k-1) W(k) )**H
!
               ELSEIF ( k>2 ) THEN
!
                  d = DLAPY2(DBLE(A(k-1,k)),DIMAG(A(k-1,k)))
                  d22 = DBLE(A(k-1,k-1))/d
                  d11 = DBLE(A(k,k))/d
                  tt = ONE/(d11*d22-ONE)
                  d12 = A(k-1,k)/d
                  d = tt/d
!
                  DO j = k - 2 , 1 , -1
                     wkm1 = d*(d11*A(j,k-1)-DCONJG(d12)*A(j,k))
                     wk = d*(d22*A(j,k)-d12*A(j,k-1))
                     DO i = j , 1 , -1
                        A(i,j) = A(i,j) - A(i,k)*DCONJG(wk) - A(i,k-1)  &
     &                           *DCONJG(wkm1)
                     ENDDO
                     A(j,k) = wk
                     A(j,k-1) = wkm1
                     A(j,j) = DCMPLX(DBLE(A(j,j)),0.0D+0)
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
!        Factorize A as L*D*L**H using the lower triangle of A
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
            IF ( (MAX(absakk,colmax)==ZERO) .OR. DISNAN(absakk) ) THEN
!
!           Column K is zero or underflow, or contains a NaN:
!           set INFO and continue
!
               IF ( Info==0 ) Info = k
               kp = k
               A(k,k) = DBLE(A(k,k))
            ELSE
!
!           ============================================================
!
!           Test for interchange
!
               IF ( absakk>=alpha*colmax ) THEN
!
!              no interchange, use 1-by-1 pivot block
!
                  kp = k
               ELSE
!
!              JMAX is the column-index of the largest off-diagonal
!              element in row IMAX, and ROWMAX is its absolute value.
!              Determine only ROWMAX.
!
                  jmax = k - 1 + IZAMAX(imax-k,A(imax,k),Lda)
                  rowmax = CABS1(A(imax,jmax))
                  IF ( imax<N ) THEN
                     jmax = imax + IZAMAX(N-imax,A(imax+1,imax),1)
                     rowmax = MAX(rowmax,CABS1(A(jmax,imax)))
                  ENDIF
!
                  IF ( absakk>=alpha*colmax*(colmax/rowmax) ) THEN
!
!                 no interchange, use 1-by-1 pivot block
!
                     kp = k
!
                  ELSEIF ( ABS(DBLE(A(imax,imax)))>=alpha*rowmax ) THEN
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
!
               ENDIF
!
!           ============================================================
!
               kk = k + kstep - 1
               IF ( kp/=kk ) THEN
!
!              Interchange rows and columns KK and KP in the trailing
!              submatrix A(k:n,k:n)
!
                  IF ( kp<N ) CALL ZSWAP(N-kp,A(kp+1,kk),1,A(kp+1,kp),1)
                  DO j = kk + 1 , kp - 1
                     t = DCONJG(A(j,kk))
                     A(j,kk) = DCONJG(A(kp,j))
                     A(kp,j) = t
                  ENDDO
                  A(kp,kk) = DCONJG(A(kp,kk))
                  r1 = DBLE(A(kk,kk))
                  A(kk,kk) = DBLE(A(kp,kp))
                  A(kp,kp) = r1
                  IF ( kstep==2 ) THEN
                     A(k,k) = DBLE(A(k,k))
                     t = A(k+1,k)
                     A(k+1,k) = A(kp,k)
                     A(kp,k) = t
                  ENDIF
               ELSE
                  A(k,k) = DBLE(A(k,k))
                  IF ( kstep==2 ) A(k+1,k+1) = DBLE(A(k+1,k+1))
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
!                 A := A - L(k)*D(k)*L(k)**H = A - W(k)*(1/D(k))*W(k)**H
!
                     r1 = ONE/DBLE(A(k,k))
                     CALL ZHER(Uplo,N-k,-r1,A(k+1,k),1,A(k+1,k+1),Lda)
!
!                 Store L(k) in column K
!
                     CALL ZDSCAL(N-k,r1,A(k+1,k),1)
                  ENDIF
!
!              2-by-2 pivot block D(k)
!
               ELSEIF ( k<N-1 ) THEN
!
!                 Perform a rank-2 update of A(k+2:n,k+2:n) as
!
!                 A := A - ( L(k) L(k+1) )*D(k)*( L(k) L(k+1) )**H
!                    = A - ( W(k) W(k+1) )*inv(D(k))*( W(k) W(k+1) )**H
!
!                 where L(k) and L(k+1) are the k-th and (k+1)-th
!                 columns of L
!
                  d = DLAPY2(DBLE(A(k+1,k)),DIMAG(A(k+1,k)))
                  d11 = DBLE(A(k+1,k+1))/d
                  d22 = DBLE(A(k,k))/d
                  tt = ONE/(d11*d22-ONE)
                  d21 = A(k+1,k)/d
                  d = tt/d
!
                  DO j = k + 2 , N
                     wk = d*(d11*A(j,k)-d21*A(j,k+1))
                     wkp1 = d*(d22*A(j,k+1)-DCONJG(d21)*A(j,k))
                     DO i = j , N
                        A(i,j) = A(i,j) - A(i,k)*DCONJG(wk) - A(i,k+1)  &
     &                           *DCONJG(wkp1)
                     ENDDO
                     A(j,k) = wk
                     A(j,k+1) = wkp1
                     A(j,j) = DCMPLX(DBLE(A(j,j)),0.0D+0)
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
!     End of ZHETF2
!
      END SUBROUTINE ZHETF2
