!*==cgebd2.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CGEBD2 reduces a general matrix to bidiagonal form using an unblocked algorithm.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CGEBD2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgebd2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgebd2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgebd2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGEBD2( M, N, A, LDA, D, E, TAUQ, TAUP, WORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, M, N
!       ..
!       .. Array Arguments ..
!       REAL               D( * ), E( * )
!       COMPLEX            A( LDA, * ), TAUP( * ), TAUQ( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CGEBD2 reduces a complex general m by n matrix A to upper or lower
!> real bidiagonal form B by a unitary transformation: Q**H * A * P = B.
!>
!> If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows in the matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns in the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          On entry, the m by n general matrix to be reduced.
!>          On exit,
!>          if m >= n, the diagonal and the first superdiagonal are
!>            overwritten with the upper bidiagonal matrix B; the
!>            elements below the diagonal, with the array TAUQ, represent
!>            the unitary matrix Q as a product of elementary
!>            reflectors, and the elements above the first superdiagonal,
!>            with the array TAUP, represent the unitary matrix P as
!>            a product of elementary reflectors;
!>          if m < n, the diagonal and the first subdiagonal are
!>            overwritten with the lower bidiagonal matrix B; the
!>            elements below the first subdiagonal, with the array TAUQ,
!>            represent the unitary matrix Q as a product of
!>            elementary reflectors, and the elements above the diagonal,
!>            with the array TAUP, represent the unitary matrix P as
!>            a product of elementary reflectors.
!>          See Further Details.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[out] D
!> \verbatim
!>          D is REAL array, dimension (min(M,N))
!>          The diagonal elements of the bidiagonal matrix B:
!>          D(i) = A(i,i).
!> \endverbatim
!>
!> \param[out] E
!> \verbatim
!>          E is REAL array, dimension (min(M,N)-1)
!>          The off-diagonal elements of the bidiagonal matrix B:
!>          if m >= n, E(i) = A(i,i+1) for i = 1,2,...,n-1;
!>          if m < n, E(i) = A(i+1,i) for i = 1,2,...,m-1.
!> \endverbatim
!>
!> \param[out] TAUQ
!> \verbatim
!>          TAUQ is COMPLEX array, dimension (min(M,N))
!>          The scalar factors of the elementary reflectors which
!>          represent the unitary matrix Q. See Further Details.
!> \endverbatim
!>
!> \param[out] TAUP
!> \verbatim
!>          TAUP is COMPLEX array, dimension (min(M,N))
!>          The scalar factors of the elementary reflectors which
!>          represent the unitary matrix P. See Further Details.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (max(M,N))
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value.
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
!> \date June 2017
!
!> \ingroup complexGEcomputational
! @precisions normal c -> s d z
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The matrices Q and P are represented as products of elementary
!>  reflectors:
!>
!>  If m >= n,
!>
!>     Q = H(1) H(2) . . . H(n)  and  P = G(1) G(2) . . . G(n-1)
!>
!>  Each H(i) and G(i) has the form:
!>
!>     H(i) = I - tauq * v * v**H  and G(i) = I - taup * u * u**H
!>
!>  where tauq and taup are complex scalars, and v and u are complex
!>  vectors; v(1:i-1) = 0, v(i) = 1, and v(i+1:m) is stored on exit in
!>  A(i+1:m,i); u(1:i) = 0, u(i+1) = 1, and u(i+2:n) is stored on exit in
!>  A(i,i+2:n); tauq is stored in TAUQ(i) and taup in TAUP(i).
!>
!>  If m < n,
!>
!>     Q = H(1) H(2) . . . H(m-1)  and  P = G(1) G(2) . . . G(m)
!>
!>  Each H(i) and G(i) has the form:
!>
!>     H(i) = I - tauq * v * v**H  and G(i) = I - taup * u * u**H
!>
!>  where tauq and taup are complex scalars, v and u are complex vectors;
!>  v(1:i) = 0, v(i+1) = 1, and v(i+2:m) is stored on exit in A(i+2:m,i);
!>  u(1:i-1) = 0, u(i) = 1, and u(i+1:n) is stored on exit in A(i,i+1:n);
!>  tauq is stored in TAUQ(i) and taup in TAUP(i).
!>
!>  The contents of A on exit are illustrated by the following examples:
!>
!>  m = 6 and n = 5 (m > n):          m = 5 and n = 6 (m < n):
!>
!>    (  d   e   u1  u1  u1 )           (  d   u1  u1  u1  u1  u1 )
!>    (  v1  d   e   u2  u2 )           (  e   d   u2  u2  u2  u2 )
!>    (  v1  v2  d   e   u3 )           (  v1  e   d   u3  u3  u3 )
!>    (  v1  v2  v3  d   e  )           (  v1  v2  e   d   u4  u4 )
!>    (  v1  v2  v3  v4  d  )           (  v1  v2  v3  e   d   u5 )
!>    (  v1  v2  v3  v4  v5 )
!>
!>  where d and e denote diagonal and off-diagonal elements of B, vi
!>  denotes an element of the vector defining H(i), and ui an element of
!>  the vector defining G(i).
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE CGEBD2(M,N,A,Lda,D,E,Tauq,Taup,Work,Info)
      USE S_CLACGV
      USE S_CLARF
      USE S_CLARFG
      USE S_XERBLA
      IMPLICIT NONE
!*--CGEBD2198
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0) ,                 &
     &                         ONE = (1.0E+0,0.0E+0)
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(INOUT) , DIMENSION(*) :: D
      REAL , INTENT(INOUT) , DIMENSION(*) :: E
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Tauq
      COMPLEX , DIMENSION(*) :: Taup
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      COMPLEX :: alpha
      INTEGER :: i
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
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      Info = 0
      IF ( M<0 ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Lda<MAX(1,M) ) THEN
         Info = -4
      ENDIF
      IF ( Info<0 ) THEN
         CALL XERBLA('CGEBD2',-Info)
         RETURN
      ENDIF
!
      IF ( M>=N ) THEN
!
!        Reduce to upper bidiagonal form
!
         DO i = 1 , N
!
!           Generate elementary reflector H(i) to annihilate A(i+1:m,i)
!
            alpha = A(i,i)
            CALL CLARFG(M-i+1,alpha,A(MIN(i+1,M),i),1,Tauq(i))
            D(i) = alpha
            A(i,i) = ONE
!
!           Apply H(i)**H to A(i:m,i+1:n) from the left
!
            IF ( i<N ) CALL CLARF('Left',M-i+1,N-i,A(i,i),1,            &
     &                            CONJG(Tauq(i)),A(i,i+1),Lda,Work)
            A(i,i) = D(i)
!
            IF ( i<N ) THEN
!
!              Generate elementary reflector G(i) to annihilate
!              A(i,i+2:n)
!
               CALL CLACGV(N-i,A(i,i+1),Lda)
               alpha = A(i,i+1)
               CALL CLARFG(N-i,alpha,A(i,MIN(i+2,N)),Lda,Taup(i))
               E(i) = alpha
               A(i,i+1) = ONE
!
!              Apply G(i) to A(i+1:m,i+1:n) from the right
!
               CALL CLARF('Right',M-i,N-i,A(i,i+1),Lda,Taup(i),         &
     &                    A(i+1,i+1),Lda,Work)
               CALL CLACGV(N-i,A(i,i+1),Lda)
               A(i,i+1) = E(i)
            ELSE
               Taup(i) = ZERO
            ENDIF
         ENDDO
      ELSE
!
!        Reduce to lower bidiagonal form
!
         DO i = 1 , M
!
!           Generate elementary reflector G(i) to annihilate A(i,i+1:n)
!
            CALL CLACGV(N-i+1,A(i,i),Lda)
            alpha = A(i,i)
            CALL CLARFG(N-i+1,alpha,A(i,MIN(i+1,N)),Lda,Taup(i))
            D(i) = alpha
            A(i,i) = ONE
!
!           Apply G(i) to A(i+1:m,i:n) from the right
!
            IF ( i<M ) CALL CLARF('Right',M-i,N-i+1,A(i,i),Lda,Taup(i), &
     &                            A(i+1,i),Lda,Work)
            CALL CLACGV(N-i+1,A(i,i),Lda)
            A(i,i) = D(i)
!
            IF ( i<M ) THEN
!
!              Generate elementary reflector H(i) to annihilate
!              A(i+2:m,i)
!
               alpha = A(i+1,i)
               CALL CLARFG(M-i,alpha,A(MIN(i+2,M),i),1,Tauq(i))
               E(i) = alpha
               A(i+1,i) = ONE
!
!              Apply H(i)**H to A(i+1:m,i+1:n) from the left
!
               CALL CLARF('Left',M-i,N-i,A(i+1,i),1,CONJG(Tauq(i)),     &
     &                    A(i+1,i+1),Lda,Work)
               A(i+1,i) = E(i)
            ELSE
               Tauq(i) = ZERO
            ENDIF
         ENDDO
      ENDIF
!
!     End of CGEBD2
!
      END SUBROUTINE CGEBD2
