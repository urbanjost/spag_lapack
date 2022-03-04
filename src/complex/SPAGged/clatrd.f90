!*==clatrd.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CLATRD reduces the first nb rows and columns of a symmetric/Hermitian matrix A to real tridiagonal form by an unitary similarity transformation.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLATRD + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clatrd.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clatrd.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clatrd.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLATRD( UPLO, N, NB, A, LDA, E, TAU, W, LDW )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            LDA, LDW, N, NB
!       ..
!       .. Array Arguments ..
!       REAL               E( * )
!       COMPLEX            A( LDA, * ), TAU( * ), W( LDW, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLATRD reduces NB rows and columns of a complex Hermitian matrix A to
!> Hermitian tridiagonal form by a unitary similarity
!> transformation Q**H * A * Q, and returns the matrices V and W which are
!> needed to apply the transformation to the unreduced part of A.
!>
!> If UPLO = 'U', CLATRD reduces the last NB rows and columns of a
!> matrix, of which the upper triangle is supplied;
!> if UPLO = 'L', CLATRD reduces the first NB rows and columns of a
!> matrix, of which the lower triangle is supplied.
!>
!> This is an auxiliary routine called by CHETRD.
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
!>          = 'U': Upper triangular
!>          = 'L': Lower triangular
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.
!> \endverbatim
!>
!> \param[in] NB
!> \verbatim
!>          NB is INTEGER
!>          The number of rows and columns to be reduced.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading
!>          n-by-n upper triangular part of A contains the upper
!>          triangular part of the matrix A, and the strictly lower
!>          triangular part of A is not referenced.  If UPLO = 'L', the
!>          leading n-by-n lower triangular part of A contains the lower
!>          triangular part of the matrix A, and the strictly upper
!>          triangular part of A is not referenced.
!>          On exit:
!>          if UPLO = 'U', the last NB columns have been reduced to
!>            tridiagonal form, with the diagonal elements overwriting
!>            the diagonal elements of A; the elements above the diagonal
!>            with the array TAU, represent the unitary matrix Q as a
!>            product of elementary reflectors;
!>          if UPLO = 'L', the first NB columns have been reduced to
!>            tridiagonal form, with the diagonal elements overwriting
!>            the diagonal elements of A; the elements below the diagonal
!>            with the array TAU, represent the  unitary matrix Q as a
!>            product of elementary reflectors.
!>          See Further Details.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] E
!> \verbatim
!>          E is REAL array, dimension (N-1)
!>          If UPLO = 'U', E(n-nb:n-1) contains the superdiagonal
!>          elements of the last NB columns of the reduced matrix;
!>          if UPLO = 'L', E(1:nb) contains the subdiagonal elements of
!>          the first NB columns of the reduced matrix.
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is COMPLEX array, dimension (N-1)
!>          The scalar factors of the elementary reflectors, stored in
!>          TAU(n-nb:n-1) if UPLO = 'U', and in TAU(1:nb) if UPLO = 'L'.
!>          See Further Details.
!> \endverbatim
!>
!> \param[out] W
!> \verbatim
!>          W is COMPLEX array, dimension (LDW,NB)
!>          The n-by-nb matrix W required to update the unreduced part
!>          of A.
!> \endverbatim
!>
!> \param[in] LDW
!> \verbatim
!>          LDW is INTEGER
!>          The leading dimension of the array W. LDW >= max(1,N).
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
!> \ingroup complexOTHERauxiliary
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  If UPLO = 'U', the matrix Q is represented as a product of elementary
!>  reflectors
!>
!>     Q = H(n) H(n-1) . . . H(n-nb+1).
!>
!>  Each H(i) has the form
!>
!>     H(i) = I - tau * v * v**H
!>
!>  where tau is a complex scalar, and v is a complex vector with
!>  v(i:n) = 0 and v(i-1) = 1; v(1:i-1) is stored on exit in A(1:i-1,i),
!>  and tau in TAU(i-1).
!>
!>  If UPLO = 'L', the matrix Q is represented as a product of elementary
!>  reflectors
!>
!>     Q = H(1) H(2) . . . H(nb).
!>
!>  Each H(i) has the form
!>
!>     H(i) = I - tau * v * v**H
!>
!>  where tau is a complex scalar, and v is a complex vector with
!>  v(1:i) = 0 and v(i+1) = 1; v(i+1:n) is stored on exit in A(i+1:n,i),
!>  and tau in TAU(i).
!>
!>  The elements of the vectors v together form the n-by-nb matrix V
!>  which is needed, with W, to apply the transformation to the unreduced
!>  part of the matrix, using a Hermitian rank-2k update of the form:
!>  A := A - V*W**H - W*V**H.
!>
!>  The contents of A on exit are illustrated by the following examples
!>  with n = 5 and nb = 2:
!>
!>  if UPLO = 'U':                       if UPLO = 'L':
!>
!>    (  a   a   a   v4  v5 )              (  d                  )
!>    (      a   a   v4  v5 )              (  1   d              )
!>    (          a   1   v5 )              (  v1  1   a          )
!>    (              d   1  )              (  v1  v2  a   a      )
!>    (                  d  )              (  v1  v2  a   a   a  )
!>
!>  where d denotes a diagonal element of the reduced matrix, a denotes
!>  an element of the original matrix that is unchanged, and vi denotes
!>  an element of the vector defining H(i).
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE CLATRD(Uplo,N,Nb,A,Lda,E,Tau,W,Ldw)
      USE S_CAXPY
      USE S_CDOTC
      USE S_CGEMV
      USE S_CHEMV
      USE S_CLACGV
      USE S_CLARFG
      USE S_CSCAL
      USE S_LSAME
      IMPLICIT NONE
!*--CLATRD211
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0) ,                 &
     &                         ONE = (1.0E+0,0.0E+0) ,                  &
     &                         HALF = (0.5E+0,0.0E+0)
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Nb
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(OUT) , DIMENSION(*) :: E
      COMPLEX , DIMENSION(*) :: Tau
      COMPLEX , DIMENSION(Ldw,*) :: W
      INTEGER :: Ldw
!
! Local variable declarations rewritten by SPAG
!
      COMPLEX :: alpha
      INTEGER :: i , iw
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
!     Quick return if possible
!
      IF ( N<=0 ) RETURN
!
      IF ( LSAME(Uplo,'U') ) THEN
!
!        Reduce last NB columns of upper triangle
!
         DO i = N , N - Nb + 1 , -1
            iw = i - N + Nb
            IF ( i<N ) THEN
!
!              Update A(1:i,i)
!
               A(i,i) = REAL(A(i,i))
               CALL CLACGV(N-i,W(i,iw+1),Ldw)
               CALL CGEMV('No transpose',i,N-i,-ONE,A(1,i+1),Lda,       &
     &                    W(i,iw+1),Ldw,ONE,A(1,i),1)
               CALL CLACGV(N-i,W(i,iw+1),Ldw)
               CALL CLACGV(N-i,A(i,i+1),Lda)
               CALL CGEMV('No transpose',i,N-i,-ONE,W(1,iw+1),Ldw,      &
     &                    A(i,i+1),Lda,ONE,A(1,i),1)
               CALL CLACGV(N-i,A(i,i+1),Lda)
               A(i,i) = REAL(A(i,i))
            ENDIF
            IF ( i>1 ) THEN
!
!              Generate elementary reflector H(i) to annihilate
!              A(1:i-2,i)
!
               alpha = A(i-1,i)
               CALL CLARFG(i-1,alpha,A(1,i),1,Tau(i-1))
               E(i-1) = alpha
               A(i-1,i) = ONE
!
!              Compute W(1:i-1,i)
!
               CALL CHEMV('Upper',i-1,ONE,A,Lda,A(1,i),1,ZERO,W(1,iw),1)
               IF ( i<N ) THEN
                  CALL CGEMV('Conjugate transpose',i-1,N-i,ONE,W(1,iw+1)&
     &                       ,Ldw,A(1,i),1,ZERO,W(i+1,iw),1)
                  CALL CGEMV('No transpose',i-1,N-i,-ONE,A(1,i+1),Lda,  &
     &                       W(i+1,iw),1,ONE,W(1,iw),1)
                  CALL CGEMV('Conjugate transpose',i-1,N-i,ONE,A(1,i+1),&
     &                       Lda,A(1,i),1,ZERO,W(i+1,iw),1)
                  CALL CGEMV('No transpose',i-1,N-i,-ONE,W(1,iw+1),Ldw, &
     &                       W(i+1,iw),1,ONE,W(1,iw),1)
               ENDIF
               CALL CSCAL(i-1,Tau(i-1),W(1,iw),1)
               alpha = -HALF*Tau(i-1)*CDOTC(i-1,W(1,iw),1,A(1,i),1)
               CALL CAXPY(i-1,alpha,A(1,i),1,W(1,iw),1)
            ENDIF
!
         ENDDO
      ELSE
!
!        Reduce first NB columns of lower triangle
!
         DO i = 1 , Nb
!
!           Update A(i:n,i)
!
            A(i,i) = REAL(A(i,i))
            CALL CLACGV(i-1,W(i,1),Ldw)
            CALL CGEMV('No transpose',N-i+1,i-1,-ONE,A(i,1),Lda,W(i,1), &
     &                 Ldw,ONE,A(i,i),1)
            CALL CLACGV(i-1,W(i,1),Ldw)
            CALL CLACGV(i-1,A(i,1),Lda)
            CALL CGEMV('No transpose',N-i+1,i-1,-ONE,W(i,1),Ldw,A(i,1), &
     &                 Lda,ONE,A(i,i),1)
            CALL CLACGV(i-1,A(i,1),Lda)
            A(i,i) = REAL(A(i,i))
            IF ( i<N ) THEN
!
!              Generate elementary reflector H(i) to annihilate
!              A(i+2:n,i)
!
               alpha = A(i+1,i)
               CALL CLARFG(N-i,alpha,A(MIN(i+2,N),i),1,Tau(i))
               E(i) = alpha
               A(i+1,i) = ONE
!
!              Compute W(i+1:n,i)
!
               CALL CHEMV('Lower',N-i,ONE,A(i+1,i+1),Lda,A(i+1,i),1,    &
     &                    ZERO,W(i+1,i),1)
               CALL CGEMV('Conjugate transpose',N-i,i-1,ONE,W(i+1,1),   &
     &                    Ldw,A(i+1,i),1,ZERO,W(1,i),1)
               CALL CGEMV('No transpose',N-i,i-1,-ONE,A(i+1,1),Lda,     &
     &                    W(1,i),1,ONE,W(i+1,i),1)
               CALL CGEMV('Conjugate transpose',N-i,i-1,ONE,A(i+1,1),   &
     &                    Lda,A(i+1,i),1,ZERO,W(1,i),1)
               CALL CGEMV('No transpose',N-i,i-1,-ONE,W(i+1,1),Ldw,     &
     &                    W(1,i),1,ONE,W(i+1,i),1)
               CALL CSCAL(N-i,Tau(i),W(i+1,i),1)
               alpha = -HALF*Tau(i)*CDOTC(N-i,W(i+1,i),1,A(i+1,i),1)
               CALL CAXPY(N-i,alpha,A(i+1,i),1,W(i+1,i),1)
            ENDIF
!
         ENDDO
      ENDIF
!
!
!     End of CLATRD
!
      END SUBROUTINE CLATRD
