!*==ssytd2.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
 
!> \brief \b SSYTD2 reduces a symmetric matrix to real symmetric tridiagonal form by an orthogonal similarity transformation (unblocked algorithm).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SSYTD2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssytd2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssytd2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssytd2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SSYTD2( UPLO, N, A, LDA, D, E, TAU, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDA, N
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), D( * ), E( * ), TAU( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SSYTD2 reduces a real symmetric matrix A to symmetric tridiagonal
!> form T by an orthogonal similarity transformation: Q**T * A * Q = T.
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
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
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
!> \ingroup realSYcomputational
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
!>     H(i) = I - tau * v * v**T
!>
!>  where tau is a real scalar, and v is a real vector with
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
      SUBROUTINE SSYTD2(Uplo,N,A,Lda,D,E,Tau,Info)
      USE S_LSAME
      USE S_SAXPY
      USE S_SDOT
      USE S_SLARFG
      USE S_SSYMV
      USE S_SSYR2
      USE S_XERBLA
      IMPLICIT NONE
!*--SSYTD2185
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ONE = 1.0 , ZERO = 0.0 , HALF = 1.0/2.0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(OUT) , DIMENSION(*) :: D
      REAL , INTENT(INOUT) , DIMENSION(*) :: E
      REAL , DIMENSION(*) :: Tau
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL :: alpha , taui
      INTEGER :: i
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
!     .. External Subroutines ..
!     ..
!     .. External Functions ..
!     ..
!     .. Intrinsic Functions ..
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
         CALL XERBLA('SSYTD2',-Info)
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
         DO i = N - 1 , 1 , -1
!
!           Generate elementary reflector H(i) = I - tau * v * v**T
!           to annihilate A(1:i-1,i+1)
!
            CALL SLARFG(i,A(i,i+1),A(1,i+1),1,taui)
            E(i) = A(i,i+1)
!
            IF ( taui/=ZERO ) THEN
!
!              Apply H(i) from both sides to A(1:i,1:i)
!
               A(i,i+1) = ONE
!
!              Compute  x := tau * A * v  storing x in TAU(1:i)
!
               CALL SSYMV(Uplo,i,taui,A,Lda,A(1,i+1),1,ZERO,Tau,1)
!
!              Compute  w := x - 1/2 * tau * (x**T * v) * v
!
               alpha = -HALF*taui*SDOT(i,Tau,1,A(1,i+1),1)
               CALL SAXPY(i,alpha,A(1,i+1),1,Tau,1)
!
!              Apply the transformation as a rank-2 update:
!                 A := A - v * w**T - w * v**T
!
               CALL SSYR2(Uplo,i,-ONE,A(1,i+1),1,Tau,1,A,Lda)
!
               A(i,i+1) = E(i)
            ENDIF
            D(i+1) = A(i+1,i+1)
            Tau(i) = taui
         ENDDO
         D(1) = A(1,1)
      ELSE
!
!        Reduce the lower triangle of A
!
         DO i = 1 , N - 1
!
!           Generate elementary reflector H(i) = I - tau * v * v**T
!           to annihilate A(i+2:n,i)
!
            CALL SLARFG(N-i,A(i+1,i),A(MIN(i+2,N),i),1,taui)
            E(i) = A(i+1,i)
!
            IF ( taui/=ZERO ) THEN
!
!              Apply H(i) from both sides to A(i+1:n,i+1:n)
!
               A(i+1,i) = ONE
!
!              Compute  x := tau * A * v  storing y in TAU(i:n-1)
!
               CALL SSYMV(Uplo,N-i,taui,A(i+1,i+1),Lda,A(i+1,i),1,ZERO, &
     &                    Tau(i),1)
!
!              Compute  w := x - 1/2 * tau * (x**T * v) * v
!
               alpha = -HALF*taui*SDOT(N-i,Tau(i),1,A(i+1,i),1)
               CALL SAXPY(N-i,alpha,A(i+1,i),1,Tau(i),1)
!
!              Apply the transformation as a rank-2 update:
!                 A := A - v * w**T - w * v**T
!
               CALL SSYR2(Uplo,N-i,-ONE,A(i+1,i),1,Tau(i),1,A(i+1,i+1), &
     &                    Lda)
!
               A(i+1,i) = E(i)
            ENDIF
            D(i) = A(i,i)
            Tau(i) = taui
         ENDDO
         D(N) = A(N,N)
      ENDIF
!
!
!     End of SSYTD2
!
      END SUBROUTINE SSYTD2
