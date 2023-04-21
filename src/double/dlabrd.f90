!*==dlabrd.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b DLABRD reduces the first nb rows and columns of a general matrix to a bidiagonal form.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLABRD + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlabrd.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlabrd.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlabrd.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLABRD( M, N, NB, A, LDA, D, E, TAUQ, TAUP, X, LDX, Y,
!                          LDY )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LDX, LDY, M, N, NB
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), D( * ), E( * ), TAUP( * ),
!      $                   TAUQ( * ), X( LDX, * ), Y( LDY, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLABRD reduces the first NB rows and columns of a real general
!> m by n matrix A to upper or lower bidiagonal form by an orthogonal
!> transformation Q**T * A * P, and returns the matrices X and Y which
!> are needed to apply the transformation to the unreduced part of A.
!>
!> If m >= n, A is reduced to upper bidiagonal form; if m < n, to lower
!> bidiagonal form.
!>
!> This is an auxiliary routine called by DGEBRD
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows in the matrix A.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns in the matrix A.
!> \endverbatim
!>
!> \param[in] NB
!> \verbatim
!>          NB is INTEGER
!>          The number of leading rows and columns of A to be reduced.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, the m by n general matrix to be reduced.
!>          On exit, the first NB rows and columns of the matrix are
!>          overwritten; the rest of the array is unchanged.
!>          If m >= n, elements on and below the diagonal in the first NB
!>            columns, with the array TAUQ, represent the orthogonal
!>            matrix Q as a product of elementary reflectors; and
!>            elements above the diagonal in the first NB rows, with the
!>            array TAUP, represent the orthogonal matrix P as a product
!>            of elementary reflectors.
!>          If m < n, elements below the diagonal in the first NB
!>            columns, with the array TAUQ, represent the orthogonal
!>            matrix Q as a product of elementary reflectors, and
!>            elements on and above the diagonal in the first NB rows,
!>            with the array TAUP, represent the orthogonal matrix P as
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
!>          D is DOUBLE PRECISION array, dimension (NB)
!>          The diagonal elements of the first NB rows and columns of
!>          the reduced matrix.  D(i) = A(i,i).
!> \endverbatim
!>
!> \param[out] E
!> \verbatim
!>          E is DOUBLE PRECISION array, dimension (NB)
!>          The off-diagonal elements of the first NB rows and columns of
!>          the reduced matrix.
!> \endverbatim
!>
!> \param[out] TAUQ
!> \verbatim
!>          TAUQ is DOUBLE PRECISION array, dimension (NB)
!>          The scalar factors of the elementary reflectors which
!>          represent the orthogonal matrix Q. See Further Details.
!> \endverbatim
!>
!> \param[out] TAUP
!> \verbatim
!>          TAUP is DOUBLE PRECISION array, dimension (NB)
!>          The scalar factors of the elementary reflectors which
!>          represent the orthogonal matrix P. See Further Details.
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension (LDX,NB)
!>          The m-by-nb matrix X required to update the unreduced part
!>          of A.
!> \endverbatim
!>
!> \param[in] LDX
!> \verbatim
!>          LDX is INTEGER
!>          The leading dimension of the array X. LDX >= max(1,M).
!> \endverbatim
!>
!> \param[out] Y
!> \verbatim
!>          Y is DOUBLE PRECISION array, dimension (LDY,NB)
!>          The n-by-nb matrix Y required to update the unreduced part
!>          of A.
!> \endverbatim
!>
!> \param[in] LDY
!> \verbatim
!>          LDY is INTEGER
!>          The leading dimension of the array Y. LDY >= max(1,N).
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
!> \ingroup doubleOTHERauxiliary
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The matrices Q and P are represented as products of elementary
!>  reflectors:
!>
!>     Q = H(1) H(2) . . . H(nb)  and  P = G(1) G(2) . . . G(nb)
!>
!>  Each H(i) and G(i) has the form:
!>
!>     H(i) = I - tauq * v * v**T  and G(i) = I - taup * u * u**T
!>
!>  where tauq and taup are real scalars, and v and u are real vectors.
!>
!>  If m >= n, v(1:i-1) = 0, v(i) = 1, and v(i:m) is stored on exit in
!>  A(i:m,i); u(1:i) = 0, u(i+1) = 1, and u(i+1:n) is stored on exit in
!>  A(i,i+1:n); tauq is stored in TAUQ(i) and taup in TAUP(i).
!>
!>  If m < n, v(1:i) = 0, v(i+1) = 1, and v(i+1:m) is stored on exit in
!>  A(i+2:m,i); u(1:i-1) = 0, u(i) = 1, and u(i:n) is stored on exit in
!>  A(i,i+1:n); tauq is stored in TAUQ(i) and taup in TAUP(i).
!>
!>  The elements of the vectors v and u together form the m-by-nb matrix
!>  V and the nb-by-n matrix U**T which are needed, with X and Y, to apply
!>  the transformation to the unreduced part of the matrix, using a block
!>  update of the form:  A := A - V*Y**T - X*U**T.
!>
!>  The contents of A on exit are illustrated by the following examples
!>  with nb = 2:
!>
!>  m = 6 and n = 5 (m > n):          m = 5 and n = 6 (m < n):
!>
!>    (  1   1   u1  u1  u1 )           (  1   u1  u1  u1  u1  u1 )
!>    (  v1  1   1   u2  u2 )           (  1   1   u2  u2  u2  u2 )
!>    (  v1  v2  a   a   a  )           (  v1  1   a   a   a   a  )
!>    (  v1  v2  a   a   a  )           (  v1  v2  a   a   a   a  )
!>    (  v1  v2  a   a   a  )           (  v1  v2  a   a   a   a  )
!>    (  v1  v2  a   a   a  )
!>
!>  where a denotes an element of the original matrix which is unchanged,
!>  vi denotes an element of the vector defining H(i), and ui an element
!>  of the vector defining G(i).
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE DLABRD(M,N,Nb,A,Lda,D,E,Tauq,Taup,X,Ldx,Y,Ldy)
      IMPLICIT NONE
!*--DLABRD213
!
!  -- LAPACK auxiliary routine (version 3.7.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2017
!
!     .. Scalar Arguments ..
      INTEGER Lda , Ldx , Ldy , M , N , Nb
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(Lda,*) , D(*) , E(*) , Taup(*) , Tauq(*) ,     &
     &                 X(Ldx,*) , Y(Ldy,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
!     ..
!     .. Local Scalars ..
      INTEGER i
!     ..
!     .. External Subroutines ..
      EXTERNAL DGEMV , DLARFG , DSCAL
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MIN
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF ( M<=0 .OR. N<=0 ) RETURN
!
      IF ( M>=N ) THEN
!
!        Reduce to upper bidiagonal form
!
         DO i = 1 , Nb
!
!           Update A(i:m,i)
!
            CALL DGEMV('No transpose',M-i+1,i-1,-ONE,A(i,1),Lda,Y(i,1), &
     &                 Ldy,ONE,A(i,i),1)
            CALL DGEMV('No transpose',M-i+1,i-1,-ONE,X(i,1),Ldx,A(1,i), &
     &                 1,ONE,A(i,i),1)
!
!           Generate reflection Q(i) to annihilate A(i+1:m,i)
!
            CALL DLARFG(M-i+1,A(i,i),A(MIN(i+1,M),i),1,Tauq(i))
            D(i) = A(i,i)
            IF ( i<N ) THEN
               A(i,i) = ONE
!
!              Compute Y(i+1:n,i)
!
               CALL DGEMV('Transpose',M-i+1,N-i,ONE,A(i,i+1),Lda,A(i,i),&
     &                    1,ZERO,Y(i+1,i),1)
               CALL DGEMV('Transpose',M-i+1,i-1,ONE,A(i,1),Lda,A(i,i),1,&
     &                    ZERO,Y(1,i),1)
               CALL DGEMV('No transpose',N-i,i-1,-ONE,Y(i+1,1),Ldy,     &
     &                    Y(1,i),1,ONE,Y(i+1,i),1)
               CALL DGEMV('Transpose',M-i+1,i-1,ONE,X(i,1),Ldx,A(i,i),1,&
     &                    ZERO,Y(1,i),1)
               CALL DGEMV('Transpose',i-1,N-i,-ONE,A(1,i+1),Lda,Y(1,i), &
     &                    1,ONE,Y(i+1,i),1)
               CALL DSCAL(N-i,Tauq(i),Y(i+1,i),1)
!
!              Update A(i,i+1:n)
!
               CALL DGEMV('No transpose',N-i,i,-ONE,Y(i+1,1),Ldy,A(i,1),&
     &                    Lda,ONE,A(i,i+1),Lda)
               CALL DGEMV('Transpose',i-1,N-i,-ONE,A(1,i+1),Lda,X(i,1), &
     &                    Ldx,ONE,A(i,i+1),Lda)
!
!              Generate reflection P(i) to annihilate A(i,i+2:n)
!
               CALL DLARFG(N-i,A(i,i+1),A(i,MIN(i+2,N)),Lda,Taup(i))
               E(i) = A(i,i+1)
               A(i,i+1) = ONE
!
!              Compute X(i+1:m,i)
!
               CALL DGEMV('No transpose',M-i,N-i,ONE,A(i+1,i+1),Lda,    &
     &                    A(i,i+1),Lda,ZERO,X(i+1,i),1)
               CALL DGEMV('Transpose',N-i,i,ONE,Y(i+1,1),Ldy,A(i,i+1),  &
     &                    Lda,ZERO,X(1,i),1)
               CALL DGEMV('No transpose',M-i,i,-ONE,A(i+1,1),Lda,X(1,i),&
     &                    1,ONE,X(i+1,i),1)
               CALL DGEMV('No transpose',i-1,N-i,ONE,A(1,i+1),Lda,      &
     &                    A(i,i+1),Lda,ZERO,X(1,i),1)
               CALL DGEMV('No transpose',M-i,i-1,-ONE,X(i+1,1),Ldx,     &
     &                    X(1,i),1,ONE,X(i+1,i),1)
               CALL DSCAL(M-i,Taup(i),X(i+1,i),1)
            ENDIF
         ENDDO
      ELSE
!
!        Reduce to lower bidiagonal form
!
         DO i = 1 , Nb
!
!           Update A(i,i:n)
!
            CALL DGEMV('No transpose',N-i+1,i-1,-ONE,Y(i,1),Ldy,A(i,1), &
     &                 Lda,ONE,A(i,i),Lda)
            CALL DGEMV('Transpose',i-1,N-i+1,-ONE,A(1,i),Lda,X(i,1),Ldx,&
     &                 ONE,A(i,i),Lda)
!
!           Generate reflection P(i) to annihilate A(i,i+1:n)
!
            CALL DLARFG(N-i+1,A(i,i),A(i,MIN(i+1,N)),Lda,Taup(i))
            D(i) = A(i,i)
            IF ( i<M ) THEN
               A(i,i) = ONE
!
!              Compute X(i+1:m,i)
!
               CALL DGEMV('No transpose',M-i,N-i+1,ONE,A(i+1,i),Lda,    &
     &                    A(i,i),Lda,ZERO,X(i+1,i),1)
               CALL DGEMV('Transpose',N-i+1,i-1,ONE,Y(i,1),Ldy,A(i,i),  &
     &                    Lda,ZERO,X(1,i),1)
               CALL DGEMV('No transpose',M-i,i-1,-ONE,A(i+1,1),Lda,     &
     &                    X(1,i),1,ONE,X(i+1,i),1)
               CALL DGEMV('No transpose',i-1,N-i+1,ONE,A(1,i),Lda,A(i,i)&
     &                    ,Lda,ZERO,X(1,i),1)
               CALL DGEMV('No transpose',M-i,i-1,-ONE,X(i+1,1),Ldx,     &
     &                    X(1,i),1,ONE,X(i+1,i),1)
               CALL DSCAL(M-i,Taup(i),X(i+1,i),1)
!
!              Update A(i+1:m,i)
!
               CALL DGEMV('No transpose',M-i,i-1,-ONE,A(i+1,1),Lda,     &
     &                    Y(i,1),Ldy,ONE,A(i+1,i),1)
               CALL DGEMV('No transpose',M-i,i,-ONE,X(i+1,1),Ldx,A(1,i),&
     &                    1,ONE,A(i+1,i),1)
!
!              Generate reflection Q(i) to annihilate A(i+2:m,i)
!
               CALL DLARFG(M-i,A(i+1,i),A(MIN(i+2,M),i),1,Tauq(i))
               E(i) = A(i+1,i)
               A(i+1,i) = ONE
!
!              Compute Y(i+1:n,i)
!
               CALL DGEMV('Transpose',M-i,N-i,ONE,A(i+1,i+1),Lda,       &
     &                    A(i+1,i),1,ZERO,Y(i+1,i),1)
               CALL DGEMV('Transpose',M-i,i-1,ONE,A(i+1,1),Lda,A(i+1,i),&
     &                    1,ZERO,Y(1,i),1)
               CALL DGEMV('No transpose',N-i,i-1,-ONE,Y(i+1,1),Ldy,     &
     &                    Y(1,i),1,ONE,Y(i+1,i),1)
               CALL DGEMV('Transpose',M-i,i,ONE,X(i+1,1),Ldx,A(i+1,i),1,&
     &                    ZERO,Y(1,i),1)
               CALL DGEMV('Transpose',i,N-i,-ONE,A(1,i+1),Lda,Y(1,i),1, &
     &                    ONE,Y(i+1,i),1)
               CALL DSCAL(N-i,Tauq(i),Y(i+1,i),1)
            ENDIF
         ENDDO
      ENDIF
!
!     End of DLABRD
!
      END SUBROUTINE DLABRD
