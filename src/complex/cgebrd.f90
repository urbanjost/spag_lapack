!*==cgebrd.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b CGEBRD
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CGEBRD + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgebrd.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgebrd.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgebrd.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGEBRD( M, N, A, LDA, D, E, TAUQ, TAUP, WORK, LWORK,
!                          INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       REAL               D( * ), E( * )
!       COMPLEX            A( LDA, * ), TAUP( * ), TAUQ( * ),
!      $                   WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CGEBRD reduces a general complex M-by-N matrix A to upper or lower
!> bidiagonal form B by a unitary transformation: Q**H * A * P = B.
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
!>          On entry, the M-by-N general matrix to be reduced.
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
!>          WORK is COMPLEX array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The length of the array WORK.  LWORK >= max(1,M,N).
!>          For optimum performance LWORK >= (M+N)*NB, where NB
!>          is the optimal blocksize.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
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
!> \date November 2017
!
!> \ingroup complexGEcomputational
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
!>  where tauq and taup are complex scalars, and v and u are complex
!>  vectors; v(1:i) = 0, v(i+1) = 1, and v(i+2:m) is stored on exit in
!>  A(i+2:m,i); u(1:i-1) = 0, u(i) = 1, and u(i+1:n) is stored on exit in
!>  A(i,i+1:n); tauq is stored in TAUQ(i) and taup in TAUP(i).
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
      SUBROUTINE CGEBRD(M,N,A,Lda,D,E,Tauq,Taup,Work,Lwork,Info)
      IMPLICIT NONE
!*--CGEBRD209
!
!  -- LAPACK computational routine (version 3.8.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2017
!
!     .. Scalar Arguments ..
      INTEGER Info , Lda , Lwork , M , N
!     ..
!     .. Array Arguments ..
      REAL D(*) , E(*)
      COMPLEX A(Lda,*) , Taup(*) , Tauq(*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX ONE
      PARAMETER (ONE=(1.0E+0,0.0E+0))
!     ..
!     .. Local Scalars ..
      LOGICAL lquery
      INTEGER i , iinfo , j , ldwrkx , ldwrky , lwkopt , minmn , nb ,   &
     &        nbmin , nx , ws
!     ..
!     .. External Subroutines ..
      EXTERNAL CGEBD2 , CGEMM , CLABRD , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN , REAL
!     ..
!     .. External Functions ..
      INTEGER ILAENV
      EXTERNAL ILAENV
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      Info = 0
      nb = MAX(1,ILAENV(1,'CGEBRD',' ',M,N,-1,-1))
      lwkopt = (M+N)*nb
      Work(1) = REAL(lwkopt)
      lquery = (Lwork==-1)
      IF ( M<0 ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Lda<MAX(1,M) ) THEN
         Info = -4
      ELSEIF ( Lwork<MAX(1,M,N) .AND. .NOT.lquery ) THEN
         Info = -10
      ENDIF
      IF ( Info<0 ) THEN
         CALL XERBLA('CGEBRD',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible
!
      minmn = MIN(M,N)
      IF ( minmn==0 ) THEN
         Work(1) = 1
         RETURN
      ENDIF
!
      ws = MAX(M,N)
      ldwrkx = M
      ldwrky = N
!
      IF ( nb>1 .AND. nb<minmn ) THEN
!
!        Set the crossover point NX.
!
         nx = MAX(nb,ILAENV(3,'CGEBRD',' ',M,N,-1,-1))
!
!        Determine when to switch from blocked to unblocked code.
!
         IF ( nx<minmn ) THEN
            ws = (M+N)*nb
            IF ( Lwork<ws ) THEN
!
!              Not enough work space for the optimal NB, consider using
!              a smaller block size.
!
               nbmin = ILAENV(2,'CGEBRD',' ',M,N,-1,-1)
               IF ( Lwork>=(M+N)*nbmin ) THEN
                  nb = Lwork/(M+N)
               ELSE
                  nb = 1
                  nx = minmn
               ENDIF
            ENDIF
         ENDIF
      ELSE
         nx = minmn
      ENDIF
!
      DO i = 1 , minmn - nx , nb
!
!        Reduce rows and columns i:i+ib-1 to bidiagonal form and return
!        the matrices X and Y which are needed to update the unreduced
!        part of the matrix
!
         CALL CLABRD(M-i+1,N-i+1,nb,A(i,i),Lda,D(i),E(i),Tauq(i),Taup(i)&
     &               ,Work,ldwrkx,Work(ldwrkx*nb+1),ldwrky)
!
!        Update the trailing submatrix A(i+ib:m,i+ib:n), using
!        an update of the form  A := A - V*Y**H - X*U**H
!
         CALL CGEMM('No transpose','Conjugate transpose',M-i-nb+1,      &
     &              N-i-nb+1,nb,-ONE,A(i+nb,i),Lda,Work(ldwrkx*nb+nb+1),&
     &              ldwrky,ONE,A(i+nb,i+nb),Lda)
         CALL CGEMM('No transpose','No transpose',M-i-nb+1,N-i-nb+1,nb, &
     &              -ONE,Work(nb+1),ldwrkx,A(i,i+nb),Lda,ONE,           &
     &              A(i+nb,i+nb),Lda)
!
!        Copy diagonal and off-diagonal elements of B back into A
!
         IF ( M>=N ) THEN
            DO j = i , i + nb - 1
               A(j,j) = D(j)
               A(j,j+1) = E(j)
            ENDDO
         ELSE
            DO j = i , i + nb - 1
               A(j,j) = D(j)
               A(j+1,j) = E(j)
            ENDDO
         ENDIF
      ENDDO
!
!     Use unblocked code to reduce the remainder of the matrix
!
      CALL CGEBD2(M-i+1,N-i+1,A(i,i),Lda,D(i),E(i),Tauq(i),Taup(i),Work,&
     &            iinfo)
      Work(1) = ws
!
!     End of CGEBRD
!
      END SUBROUTINE CGEBRD
