!*==zsptrf.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZSPTRF
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZSPTRF + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsptrf.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsptrf.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsptrf.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZSPTRF( UPLO, N, AP, IPIV, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX*16         AP( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZSPTRF computes the factorization of a complex symmetric matrix A
!> stored in packed format using the Bunch-Kaufman diagonal pivoting
!> method:
!>
!>    A = U*D*U**T  or  A = L*D*L**T
!>
!> where U (or L) is a product of permutation and unit upper (lower)
!> triangular matrices, and D is symmetric and block diagonal with
!> 1-by-1 and 2-by-2 diagonal blocks.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  Upper triangle of A is stored;
!>          = 'L':  Lower triangle of A is stored.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] AP
!> \verbatim
!>          AP is COMPLEX*16 array, dimension (N*(N+1)/2)
!>          On entry, the upper or lower triangle of the symmetric matrix
!>          A, packed columnwise in a linear array.  The j-th column of A
!>          is stored in the array AP as follows:
!>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
!>          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
!>
!>          On exit, the block diagonal matrix D and the multipliers used
!>          to obtain the factor U or L, stored as a packed triangular
!>          matrix overwriting A (see below for further details).
!> \endverbatim
!>
!> \param[out] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          Details of the interchanges and the block structure of D.
!>          If IPIV(k) > 0, then rows and columns k and IPIV(k) were
!>          interchanged and D(k,k) is a 1-by-1 diagonal block.
!>          If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and
!>          columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
!>          is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =
!>          IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were
!>          interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
!>          > 0: if INFO = i, D(i,i) is exactly zero.  The factorization
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
!> \ingroup complex16OTHERcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  5-96 - Based on modifications by J. Lewis, Boeing Computer Services
!>         Company
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
!>
!  =====================================================================
      SUBROUTINE ZSPTRF(Uplo,N,Ap,Ipiv,Info)
      USE F77KINDS                        
      USE S_IZAMAX
      USE S_LSAME
      USE S_XERBLA
      USE S_ZSCAL
      USE S_ZSPR
      USE S_ZSWAP
      IMPLICIT NONE
!*--ZSPTRF169
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              EIGHT = 8.0D+0 , SEVTEN = 17.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0)
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Ap
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: absakk , alpha , colmax , rowmax
      REAL(R8KIND) :: CABS1
      COMPLEX(CX16KIND) :: d11 , d12 , d21 , d22 , r1 , t , wk , wkm1 , &
     &                     wkp1 , zdum
      INTEGER :: i , imax , j , jmax , k , kc , kk , knc , kp , kpc ,   &
     &           kstep , kx , npp
      LOGICAL :: upper
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
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZSPTRF',-Info)
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
         kc = (N-1)*N/2 + 1
         DO
            knc = kc
!
!        If K < 1, exit from loop
!
            IF ( k<1 ) EXIT
            kstep = 1
!
!        Determine rows and columns to be interchanged and whether
!        a 1-by-1 or 2-by-2 pivot block will be used
!
            absakk = CABS1(Ap(kc+k-1))
!
!        IMAX is the row-index of the largest off-diagonal element in
!        column K, and COLMAX is its absolute value
!
            IF ( k>1 ) THEN
               imax = IZAMAX(k-1,Ap(kc),1)
               colmax = CABS1(Ap(kc+imax-1))
            ELSE
               colmax = ZERO
            ENDIF
!
            IF ( MAX(absakk,colmax)==ZERO ) THEN
!
!           Column K is zero: set INFO and continue
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
                  rowmax = ZERO
                  jmax = imax
                  kx = imax*(imax+1)/2 + imax
                  DO j = imax + 1 , k
                     IF ( CABS1(Ap(kx))>rowmax ) THEN
                        rowmax = CABS1(Ap(kx))
                        jmax = j
                     ENDIF
                     kx = kx + j
                  ENDDO
                  kpc = (imax-1)*imax/2 + 1
                  IF ( imax>1 ) THEN
                     jmax = IZAMAX(imax-1,Ap(kpc),1)
                     rowmax = MAX(rowmax,CABS1(Ap(kpc+jmax-1)))
                  ENDIF
!
                  IF ( absakk>=alpha*colmax*(colmax/rowmax) ) THEN
!
!                 no interchange, use 1-by-1 pivot block
!
                     kp = k
                  ELSEIF ( CABS1(Ap(kpc+imax-1))>=alpha*rowmax ) THEN
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
               IF ( kstep==2 ) knc = knc - k + 1
               IF ( kp/=kk ) THEN
!
!              Interchange rows and columns KK and KP in the leading
!              submatrix A(1:k,1:k)
!
                  CALL ZSWAP(kp-1,Ap(knc),1,Ap(kpc),1)
                  kx = kpc + kp - 1
                  DO j = kp + 1 , kk - 1
                     kx = kx + j - 1
                     t = Ap(knc+j-1)
                     Ap(knc+j-1) = Ap(kx)
                     Ap(kx) = t
                  ENDDO
                  t = Ap(knc+kk-1)
                  Ap(knc+kk-1) = Ap(kpc+kp-1)
                  Ap(kpc+kp-1) = t
                  IF ( kstep==2 ) THEN
                     t = Ap(kc+k-2)
                     Ap(kc+k-2) = Ap(kc+kp-1)
                     Ap(kc+kp-1) = t
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
                  r1 = CONE/Ap(kc+k-1)
                  CALL ZSPR(Uplo,k-1,-r1,Ap(kc),1,Ap)
!
!              Store U(k) in column k
!
                  CALL ZSCAL(k-1,r1,Ap(kc),1)
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
                  d12 = Ap(k-1+(k-1)*k/2)
                  d22 = Ap(k-1+(k-2)*(k-1)/2)/d12
                  d11 = Ap(k+(k-1)*k/2)/d12
                  t = CONE/(d11*d22-CONE)
                  d12 = t/d12
!
                  DO j = k - 2 , 1 , -1
                     wkm1 = d12*(d11*Ap(j+(k-2)*(k-1)/2)-Ap(j+(k-1)*k/2)&
     &                      )
                     wk = d12*(d22*Ap(j+(k-1)*k/2)-Ap(j+(k-2)*(k-1)/2))
                     DO i = j , 1 , -1
                        Ap(i+(j-1)*j/2) = Ap(i+(j-1)*j/2)               &
     &                     - Ap(i+(k-1)*k/2)*wk - Ap(i+(k-2)*(k-1)/2)   &
     &                     *wkm1
                     ENDDO
                     Ap(j+(k-1)*k/2) = wk
                     Ap(j+(k-2)*(k-1)/2) = wkm1
                  ENDDO
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
            kc = knc - k
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
         kc = 1
         npp = N*(N+1)/2
         DO
            knc = kc
!
!        If K > N, exit from loop
!
            IF ( k>N ) EXIT
            kstep = 1
!
!        Determine rows and columns to be interchanged and whether
!        a 1-by-1 or 2-by-2 pivot block will be used
!
            absakk = CABS1(Ap(kc))
!
!        IMAX is the row-index of the largest off-diagonal element in
!        column K, and COLMAX is its absolute value
!
            IF ( k<N ) THEN
               imax = k + IZAMAX(N-k,Ap(kc+1),1)
               colmax = CABS1(Ap(kc+imax-k))
            ELSE
               colmax = ZERO
            ENDIF
!
            IF ( MAX(absakk,colmax)==ZERO ) THEN
!
!           Column K is zero: set INFO and continue
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
                  rowmax = ZERO
                  kx = kc + imax - k
                  DO j = k , imax - 1
                     IF ( CABS1(Ap(kx))>rowmax ) THEN
                        rowmax = CABS1(Ap(kx))
                        jmax = j
                     ENDIF
                     kx = kx + N - j
                  ENDDO
                  kpc = npp - (N-imax+1)*(N-imax+2)/2 + 1
                  IF ( imax<N ) THEN
                     jmax = imax + IZAMAX(N-imax,Ap(kpc+1),1)
                     rowmax = MAX(rowmax,CABS1(Ap(kpc+jmax-imax)))
                  ENDIF
!
                  IF ( absakk>=alpha*colmax*(colmax/rowmax) ) THEN
!
!                 no interchange, use 1-by-1 pivot block
!
                     kp = k
                  ELSEIF ( CABS1(Ap(kpc))>=alpha*rowmax ) THEN
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
               IF ( kstep==2 ) knc = knc + N - k + 1
               IF ( kp/=kk ) THEN
!
!              Interchange rows and columns KK and KP in the trailing
!              submatrix A(k:n,k:n)
!
                  IF ( kp<N )                                           &
     &                 CALL ZSWAP(N-kp,Ap(knc+kp-kk+1),1,Ap(kpc+1),1)
                  kx = knc + kp - kk
                  DO j = kk + 1 , kp - 1
                     kx = kx + N - j + 1
                     t = Ap(knc+j-kk)
                     Ap(knc+j-kk) = Ap(kx)
                     Ap(kx) = t
                  ENDDO
                  t = Ap(knc)
                  Ap(knc) = Ap(kpc)
                  Ap(kpc) = t
                  IF ( kstep==2 ) THEN
                     t = Ap(kc+1)
                     Ap(kc+1) = Ap(kc+kp-k)
                     Ap(kc+kp-k) = t
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
                     r1 = CONE/Ap(kc)
                     CALL ZSPR(Uplo,N-k,-r1,Ap(kc+1),1,Ap(kc+N-k+1))
!
!                 Store L(k) in column K
!
                     CALL ZSCAL(N-k,r1,Ap(kc+1),1)
                  ENDIF
!
!              2-by-2 pivot block D(k): columns K and K+1 now hold
!
!              ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k)
!
!              where L(k) and L(k+1) are the k-th and (k+1)-th columns
!              of L
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
                  d21 = Ap(k+1+(k-1)*(2*N-k)/2)
                  d11 = Ap(k+1+k*(2*N-k-1)/2)/d21
                  d22 = Ap(k+(k-1)*(2*N-k)/2)/d21
                  t = CONE/(d11*d22-CONE)
                  d21 = t/d21
!
                  DO j = k + 2 , N
                     wk = d21*(d11*Ap(j+(k-1)*(2*N-k)/2)                &
     &                    -Ap(j+k*(2*N-k-1)/2))
                     wkp1 = d21*(d22*Ap(j+k*(2*N-k-1)/2)                &
     &                      -Ap(j+(k-1)*(2*N-k)/2))
                     DO i = j , N
                        Ap(i+(j-1)*(2*N-j)/2) = Ap(i+(j-1)*(2*N-j)/2)   &
     &                     - Ap(i+(k-1)*(2*N-k)/2)                      &
     &                     *wk - Ap(i+k*(2*N-k-1)/2)*wkp1
                     ENDDO
                     Ap(j+(k-1)*(2*N-k)/2) = wk
                     Ap(j+k*(2*N-k-1)/2) = wkp1
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
            kc = knc + N - k + 2
         ENDDO
!
      ENDIF
!
!
!     End of ZSPTRF
!
      END SUBROUTINE ZSPTRF
