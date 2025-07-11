!*==zhetd2.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
 
!> \brief \b ZHETD2 reduces a Hermitian matrix to real symmetric tridiagonal form by an unitary similarity transformation (unblocked algorithm).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZHETD2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhetd2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhetd2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhetd2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZHETD2( UPLO, N, A, LDA, D, E, TAU, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDA, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   D( * ), E( * )
!       COMPLEX*16         A( LDA, * ), TAU( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZHETD2 reduces a complex Hermitian matrix A to real symmetric
!> tridiagonal form T by a unitary similarity transformation:
!> Q**H * A * Q = T.
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
!>          On exit, if UPLO = 'U', the diagonal and first superdiagonal
!>          of A are overwritten by the corresponding elements of the
!>          tridiagonal matrix T, and the elements above the first
!>          superdiagonal, with the array TAU, represent the unitary
!>          matrix Q as a product of elementary reflectors; if UPLO
!>          = 'L', the diagonal and first subdiagonal of A are over-
!>          written by the corresponding elements of the tridiagonal
!>          matrix T, and the elements below the first subdiagonal, with
!>          the array TAU, represent the unitary matrix Q as a product
!>          of elementary reflectors. See Further Details.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>          The diagonal elements of the tridiagonal matrix T:
!>          D(i) = A(i,i).
!> \endverbatim
!>
!> \param[out] E
!> \verbatim
!>          E is DOUBLE PRECISION array, dimension (N-1)
!>          The off-diagonal elements of the tridiagonal matrix T:
!>          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is COMPLEX*16 array, dimension (N-1)
!>          The scalar factors of the elementary reflectors (see Further
!>          Details).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
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
!>
!>  If UPLO = 'U', the matrix Q is represented as a product of elementary
!>  reflectors
!>
!>     Q = H(n-1) . . . H(2) H(1).
!>
!>  Each H(i) has the form
!>
!>     H(i) = I - tau * v * v**H
!>
!>  where tau is a complex scalar, and v is a complex vector with
!>  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
!>  A(1:i-1,i+1), and tau in TAU(i).
!>
!>  If UPLO = 'L', the matrix Q is represented as a product of elementary
!>  reflectors
!>
!>     Q = H(1) H(2) . . . H(n-1).
!>
!>  Each H(i) has the form
!>
!>     H(i) = I - tau * v * v**H
!>
!>  where tau is a complex scalar, and v is a complex vector with
!>  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
!>  and tau in TAU(i).
!>
!>  The contents of A on exit are illustrated by the following examples
!>  with n = 5:
!>
!>  if UPLO = 'U':                       if UPLO = 'L':
!>
!>    (  d   e   v2  v3  v4 )              (  d                  )
!>    (      d   e   v3  v4 )              (  e   d              )
!>    (          d   e   v4 )              (  v1  e   d          )
!>    (              d   e  )              (  v1  v2  e   d      )
!>    (                  d  )              (  v1  v2  v3  e   d  )
!>
!>  where d and e denote diagonal and off-diagonal elements of T, and vi
!>  denotes an element of the vector defining H(i).
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE ZHETD2(Uplo,N,A,Lda,D,E,Tau,Info)
      IMPLICIT NONE
!*--ZHETD2180
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
      DOUBLE PRECISION D(*) , E(*)
      COMPLEX*16 A(Lda,*) , Tau(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16 ONE , ZERO , HALF
      PARAMETER (ONE=(1.0D+0,0.0D+0),ZERO=(0.0D+0,0.0D+0),              &
     &           HALF=(0.5D+0,0.0D+0))
!     ..
!     .. Local Scalars ..
      LOGICAL upper
      INTEGER i
      COMPLEX*16 alpha , taui
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA , ZAXPY , ZHEMV , ZHER2 , ZLARFG
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      COMPLEX*16 ZDOTC
      EXTERNAL LSAME , ZDOTC
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DBLE , MAX , MIN
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
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -4
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZHETD2',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N<=0 ) RETURN
!
      IF ( upper ) THEN
!
!        Reduce the upper triangle of A
!
         A(N,N) = DBLE(A(N,N))
         DO i = N - 1 , 1 , -1
!
!           Generate elementary reflector H(i) = I - tau * v * v**H
!           to annihilate A(1:i-1,i+1)
!
            alpha = A(i,i+1)
            CALL ZLARFG(i,alpha,A(1,i+1),1,taui)
            E(i) = alpha
!
            IF ( taui/=ZERO ) THEN
!
!              Apply H(i) from both sides to A(1:i,1:i)
!
               A(i,i+1) = ONE
!
!              Compute  x := tau * A * v  storing x in TAU(1:i)
!
               CALL ZHEMV(Uplo,i,taui,A,Lda,A(1,i+1),1,ZERO,Tau,1)
!
!              Compute  w := x - 1/2 * tau * (x**H * v) * v
!
               alpha = -HALF*taui*ZDOTC(i,Tau,1,A(1,i+1),1)
               CALL ZAXPY(i,alpha,A(1,i+1),1,Tau,1)
!
!              Apply the transformation as a rank-2 update:
!                 A := A - v * w**H - w * v**H
!
               CALL ZHER2(Uplo,i,-ONE,A(1,i+1),1,Tau,1,A,Lda)
!
            ELSE
               A(i,i) = DBLE(A(i,i))
            ENDIF
            A(i,i+1) = E(i)
            D(i+1) = A(i+1,i+1)
            Tau(i) = taui
         ENDDO
         D(1) = A(1,1)
      ELSE
!
!        Reduce the lower triangle of A
!
         A(1,1) = DBLE(A(1,1))
         DO i = 1 , N - 1
!
!           Generate elementary reflector H(i) = I - tau * v * v**H
!           to annihilate A(i+2:n,i)
!
            alpha = A(i+1,i)
            CALL ZLARFG(N-i,alpha,A(MIN(i+2,N),i),1,taui)
            E(i) = alpha
!
            IF ( taui/=ZERO ) THEN
!
!              Apply H(i) from both sides to A(i+1:n,i+1:n)
!
               A(i+1,i) = ONE
!
!              Compute  x := tau * A * v  storing y in TAU(i:n-1)
!
               CALL ZHEMV(Uplo,N-i,taui,A(i+1,i+1),Lda,A(i+1,i),1,ZERO, &
     &                    Tau(i),1)
!
!              Compute  w := x - 1/2 * tau * (x**H * v) * v
!
               alpha = -HALF*taui*ZDOTC(N-i,Tau(i),1,A(i+1,i),1)
               CALL ZAXPY(N-i,alpha,A(i+1,i),1,Tau(i),1)
!
!              Apply the transformation as a rank-2 update:
!                 A := A - v * w**H - w * v**H
!
               CALL ZHER2(Uplo,N-i,-ONE,A(i+1,i),1,Tau(i),1,A(i+1,i+1), &
     &                    Lda)
!
            ELSE
               A(i+1,i+1) = DBLE(A(i+1,i+1))
            ENDIF
            A(i+1,i) = E(i)
            D(i) = A(i,i)
            Tau(i) = taui
         ENDDO
         D(N) = A(N,N)
      ENDIF
!
!
!     End of ZHETD2
!
      END SUBROUTINE ZHETD2
