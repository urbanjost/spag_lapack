!*==zgeqrt2.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b ZGEQRT2 computes a QR factorization of a general real or complex matrix using the compact WY representation of Q.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZGEQRT2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgeqrt2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgeqrt2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgeqrt2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZGEQRT2( M, N, A, LDA, T, LDT, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER   INFO, LDA, LDT, M, N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16   A( LDA, * ), T( LDT, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZGEQRT2 computes a QR factorization of a complex M-by-N matrix A,
!> using the compact WY representation of Q.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= N.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          On entry, the complex M-by-N matrix A.  On exit, the elements on and
!>          above the diagonal contain the N-by-N upper triangular matrix R; the
!>          elements below the diagonal are the columns of V.  See below for
!>          further details.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[out] T
!> \verbatim
!>          T is COMPLEX*16 array, dimension (LDT,N)
!>          The N-by-N upper triangular factor of the block reflector.
!>          The elements on and above the diagonal contain the block
!>          reflector T; the elements below the diagonal are not used.
!>          See below for further details.
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is INTEGER
!>          The leading dimension of the array T.  LDT >= max(1,N).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
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
!> \ingroup complex16GEcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The matrix V stores the elementary reflectors H(i) in the i-th column
!>  below the diagonal. For example, if M=5 and N=3, the matrix V is
!>
!>               V = (  1       )
!>                   ( v1  1    )
!>                   ( v1 v2  1 )
!>                   ( v1 v2 v3 )
!>                   ( v1 v2 v3 )
!>
!>  where the vi's represent the vectors which define H(i), which are returned
!>  in the matrix A.  The 1's along the diagonal of V are not stored in A.  The
!>  block reflector H is then given by
!>
!>               H = I - V * T * V**H
!>
!>  where V**H is the conjugate transpose of V.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE ZGEQRT2(M,N,A,Lda,T,Ldt,Info)
      IMPLICIT NONE
!*--ZGEQRT2131
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Info , Lda , Ldt , M , N
!     ..
!     .. Array Arguments ..
      COMPLEX*16 A(Lda,*) , T(Ldt,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16 ONE , ZERO
      PARAMETER (ONE=(1.0D+00,0.0D+00),ZERO=(0.0D+00,0.0D+00))
!     ..
!     .. Local Scalars ..
      INTEGER i , k
      COMPLEX*16 aii , alpha
!     ..
!     .. External Subroutines ..
      EXTERNAL ZLARFG , ZGEMV , ZGERC , ZTRMV , XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      Info = 0
      IF ( M<0 ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Lda<MAX(1,M) ) THEN
         Info = -4
      ELSEIF ( Ldt<MAX(1,N) ) THEN
         Info = -6
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZGEQRT2',-Info)
         RETURN
      ENDIF
!
      k = MIN(M,N)
!
      DO i = 1 , k
!
!        Generate elem. refl. H(i) to annihilate A(i+1:m,i), tau(I) -> T(I,1)
!
         CALL ZLARFG(M-i+1,A(i,i),A(MIN(i+1,M),i),1,T(i,1))
         IF ( i<N ) THEN
!
!           Apply H(i) to A(I:M,I+1:N) from the left
!
            aii = A(i,i)
            A(i,i) = ONE
!
!           W(1:N-I) := A(I:M,I+1:N)^H * A(I:M,I) [W = T(:,N)]
!
            CALL ZGEMV('C',M-i+1,N-i,ONE,A(i,i+1),Lda,A(i,i),1,ZERO,    &
     &                 T(1,N),1)
!
!           A(I:M,I+1:N) = A(I:m,I+1:N) + alpha*A(I:M,I)*W(1:N-1)^H
!
            alpha = -CONJG(T(i,1))
            CALL ZGERC(M-i+1,N-i,alpha,A(i,i),1,T(1,N),1,A(i,i+1),Lda)
            A(i,i) = aii
         ENDIF
      ENDDO
!
      DO i = 2 , N
         aii = A(i,i)
         A(i,i) = ONE
!
!        T(1:I-1,I) := alpha * A(I:M,1:I-1)**H * A(I:M,I)
!
         alpha = -T(i,1)
         CALL ZGEMV('C',M-i+1,i-1,alpha,A(i,1),Lda,A(i,i),1,ZERO,T(1,i),&
     &              1)
         A(i,i) = aii
!
!        T(1:I-1,I) := T(1:I-1,1:I-1) * T(1:I-1,I)
!
         CALL ZTRMV('U','N','N',i-1,T,Ldt,T(1,i),1)
!
!           T(I,I) = tau(I)
!
         T(i,i) = T(i,1)
         T(i,1) = ZERO
      ENDDO
 
!
!     End of ZGEQRT2
!
      END SUBROUTINE ZGEQRT2
