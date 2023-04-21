!*==ssptrd.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SSPTRD
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SSPTRD + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssptrd.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssptrd.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssptrd.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SSPTRD( UPLO, N, AP, D, E, TAU, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, N
!       ..
!       .. Array Arguments ..
!       REAL               AP( * ), D( * ), E( * ), TAU( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SSPTRD reduces a real symmetric matrix A stored in packed form to
!> symmetric tridiagonal form T by an orthogonal similarity
!> transformation: Q**T * A * Q = T.
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
!>          AP is REAL array, dimension (N*(N+1)/2)
!>          On entry, the upper or lower triangle of the symmetric matrix
!>          A, packed columnwise in a linear array.  The j-th column of A
!>          is stored in the array AP as follows:
!>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
!>          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.
!>          On exit, if UPLO = 'U', the diagonal and first superdiagonal
!>          of A are overwritten by the corresponding elements of the
!>          tridiagonal matrix T, and the elements above the first
!>          superdiagonal, with the array TAU, represent the orthogonal
!>          matrix Q as a product of elementary reflectors; if UPLO
!>          = 'L', the diagonal and first subdiagonal of A are over-
!>          written by the corresponding elements of the tridiagonal
!>          matrix T, and the elements below the first subdiagonal, with
!>          the array TAU, represent the orthogonal matrix Q as a product
!>          of elementary reflectors. See Further Details.
!> \endverbatim
!>
!> \param[out] D
!> \verbatim
!>          D is REAL array, dimension (N)
!>          The diagonal elements of the tridiagonal matrix T:
!>          D(i) = A(i,i).
!> \endverbatim
!>
!> \param[out] E
!> \verbatim
!>          E is REAL array, dimension (N-1)
!>          The off-diagonal elements of the tridiagonal matrix T:
!>          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is REAL array, dimension (N-1)
!>          The scalar factors of the elementary reflectors (see Further
!>          Details).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
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
!> \ingroup realOTHERcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  If UPLO = 'U', the matrix Q is represented as a product of elementary
!>  reflectors
!>
!>     Q = H(n-1) . . . H(2) H(1).
!>
!>  Each H(i) has the form
!>
!>     H(i) = I - tau * v * v**T
!>
!>  where tau is a real scalar, and v is a real vector with
!>  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in AP,
!>  overwriting A(1:i-1,i+1), and tau is stored in TAU(i).
!>
!>  If UPLO = 'L', the matrix Q is represented as a product of elementary
!>  reflectors
!>
!>     Q = H(1) H(2) . . . H(n-1).
!>
!>  Each H(i) has the form
!>
!>     H(i) = I - tau * v * v**T
!>
!>  where tau is a real scalar, and v is a real vector with
!>  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in AP,
!>  overwriting A(i+2:n,i), and tau is stored in TAU(i).
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE SSPTRD(Uplo,N,Ap,D,E,Tau,Info)
      IMPLICIT NONE
!*--SSPTRD154
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER Info , N
!     ..
!     .. Array Arguments ..
      REAL Ap(*) , D(*) , E(*) , Tau(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ONE , ZERO , HALF
      PARAMETER (ONE=1.0,ZERO=0.0,HALF=1.0/2.0)
!     ..
!     .. Local Scalars ..
      LOGICAL upper
      INTEGER i , i1 , i1i1 , ii
      REAL alpha , taui
!     ..
!     .. External Subroutines ..
      EXTERNAL SAXPY , SLARFG , SSPMV , SSPR2 , XERBLA
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      REAL SDOT
      EXTERNAL LSAME , SDOT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      Info = 0
      upper = LSAME(Uplo,'U')
      IF ( .NOT.upper .AND. .NOT.LSAME(Uplo,'L') ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('SSPTRD',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N<=0 ) RETURN
!
      IF ( upper ) THEN
!
!        Reduce the upper triangle of A.
!        I1 is the index in AP of A(1,I+1).
!
         i1 = N*(N-1)/2 + 1
         DO i = N - 1 , 1 , -1
!
!           Generate elementary reflector H(i) = I - tau * v * v**T
!           to annihilate A(1:i-1,i+1)
!
            CALL SLARFG(i,Ap(i1+i-1),Ap(i1),1,taui)
            E(i) = Ap(i1+i-1)
!
            IF ( taui/=ZERO ) THEN
!
!              Apply H(i) from both sides to A(1:i,1:i)
!
               Ap(i1+i-1) = ONE
!
!              Compute  y := tau * A * v  storing y in TAU(1:i)
!
               CALL SSPMV(Uplo,i,taui,Ap,Ap(i1),1,ZERO,Tau,1)
!
!              Compute  w := y - 1/2 * tau * (y**T *v) * v
!
               alpha = -HALF*taui*SDOT(i,Tau,1,Ap(i1),1)
               CALL SAXPY(i,alpha,Ap(i1),1,Tau,1)
!
!              Apply the transformation as a rank-2 update:
!                 A := A - v * w**T - w * v**T
!
               CALL SSPR2(Uplo,i,-ONE,Ap(i1),1,Tau,1,Ap)
!
               Ap(i1+i-1) = E(i)
            ENDIF
            D(i+1) = Ap(i1+i)
            Tau(i) = taui
            i1 = i1 - i
         ENDDO
         D(1) = Ap(1)
      ELSE
!
!        Reduce the lower triangle of A. II is the index in AP of
!        A(i,i) and I1I1 is the index of A(i+1,i+1).
!
         ii = 1
         DO i = 1 , N - 1
            i1i1 = ii + N - i + 1
!
!           Generate elementary reflector H(i) = I - tau * v * v**T
!           to annihilate A(i+2:n,i)
!
            CALL SLARFG(N-i,Ap(ii+1),Ap(ii+2),1,taui)
            E(i) = Ap(ii+1)
!
            IF ( taui/=ZERO ) THEN
!
!              Apply H(i) from both sides to A(i+1:n,i+1:n)
!
               Ap(ii+1) = ONE
!
!              Compute  y := tau * A * v  storing y in TAU(i:n-1)
!
               CALL SSPMV(Uplo,N-i,taui,Ap(i1i1),Ap(ii+1),1,ZERO,Tau(i),&
     &                    1)
!
!              Compute  w := y - 1/2 * tau * (y**T *v) * v
!
               alpha = -HALF*taui*SDOT(N-i,Tau(i),1,Ap(ii+1),1)
               CALL SAXPY(N-i,alpha,Ap(ii+1),1,Tau(i),1)
!
!              Apply the transformation as a rank-2 update:
!                 A := A - v * w**T - w * v**T
!
               CALL SSPR2(Uplo,N-i,-ONE,Ap(ii+1),1,Tau(i),1,Ap(i1i1))
!
               Ap(ii+1) = E(i)
            ENDIF
            D(i) = Ap(ii)
            Tau(i) = taui
            ii = i1i1
         ENDDO
         D(N) = Ap(ii)
      ENDIF
!
!
!     End of SSPTRD
!
      END SUBROUTINE SSPTRD
