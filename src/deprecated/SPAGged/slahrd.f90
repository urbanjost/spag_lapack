!*==slahrd.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SLAHRD reduces the first nb columns of a general rectangular matrix A so that elements below the k-th subdiagonal are zero, and returns auxiliary matrices which are needed to apply the transformation to the unreduced part of A.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLAHRD + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slahrd.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slahrd.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slahrd.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLAHRD( N, K, NB, A, LDA, TAU, T, LDT, Y, LDY )
!
!       .. Scalar Arguments ..
!       INTEGER            K, LDA, LDT, LDY, N, NB
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), T( LDT, NB ), TAU( NB ),
!      $                   Y( LDY, NB )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> This routine is deprecated and has been replaced by routine SLAHR2.
!>
!> SLAHRD reduces the first NB columns of a real general n-by-(n-k+1)
!> matrix A so that elements below the k-th subdiagonal are zero. The
!> reduction is performed by an orthogonal similarity transformation
!> Q**T * A * Q. The routine returns the matrices V and T which determine
!> Q as a block reflector I - V*T*V**T, and also the matrix Y = A * V * T.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The offset for the reduction. Elements below the k-th
!>          subdiagonal in the first NB columns are reduced to zero.
!> \endverbatim
!>
!> \param[in] NB
!> \verbatim
!>          NB is INTEGER
!>          The number of columns to be reduced.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N-K+1)
!>          On entry, the n-by-(n-k+1) general matrix A.
!>          On exit, the elements on and above the k-th subdiagonal in
!>          the first NB columns are overwritten with the corresponding
!>          elements of the reduced matrix; the elements below the k-th
!>          subdiagonal, with the array TAU, represent the matrix Q as a
!>          product of elementary reflectors. The other columns of A are
!>          unchanged. See Further Details.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is REAL array, dimension (NB)
!>          The scalar factors of the elementary reflectors. See Further
!>          Details.
!> \endverbatim
!>
!> \param[out] T
!> \verbatim
!>          T is REAL array, dimension (LDT,NB)
!>          The upper triangular matrix T.
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is INTEGER
!>          The leading dimension of the array T.  LDT >= NB.
!> \endverbatim
!>
!> \param[out] Y
!> \verbatim
!>          Y is REAL array, dimension (LDY,NB)
!>          The n-by-nb matrix Y.
!> \endverbatim
!>
!> \param[in] LDY
!> \verbatim
!>          LDY is INTEGER
!>          The leading dimension of the array Y. LDY >= N.
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
!> \ingroup realOTHERauxiliary
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The matrix Q is represented as a product of nb elementary reflectors
!>
!>     Q = H(1) H(2) . . . H(nb).
!>
!>  Each H(i) has the form
!>
!>     H(i) = I - tau * v * v**T
!>
!>  where tau is a real scalar, and v is a real vector with
!>  v(1:i+k-1) = 0, v(i+k) = 1; v(i+k+1:n) is stored on exit in
!>  A(i+k+1:n,i), and tau in TAU(i).
!>
!>  The elements of the vectors v together form the (n-k+1)-by-nb matrix
!>  V which is needed, with T and Y, to apply the transformation to the
!>  unreduced part of the matrix, using an update of the form:
!>  A := (I - V*T*V**T) * (A - Y*V**T).
!>
!>  The contents of A on exit are illustrated by the following example
!>  with n = 7, k = 3 and nb = 2:
!>
!>     ( a   h   a   a   a )
!>     ( a   h   a   a   a )
!>     ( a   h   a   a   a )
!>     ( h   h   a   a   a )
!>     ( v1  h   a   a   a )
!>     ( v1  v2  a   a   a )
!>     ( v1  v2  a   a   a )
!>
!>  where a denotes an element of the original matrix A, h denotes a
!>  modified element of the upper Hessenberg matrix H, and vi denotes an
!>  element of the vector defining H(i).
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE SLAHRD(N,K,Nb,A,Lda,Tau,T,Ldt,Y,Ldy)
      USE S_SAXPY
      USE S_SCOPY
      USE S_SGEMV
      USE S_SLARFG
      USE S_SSCAL
      USE S_STRMV
      IMPLICIT NONE
!*--SLAHRD177
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER :: N
      INTEGER , INTENT(IN) :: K
      INTEGER , INTENT(IN) :: Nb
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Nb) :: Tau
      REAL , DIMENSION(Ldt,Nb) :: T
      INTEGER :: Ldt
      REAL , DIMENSION(Ldy,Nb) :: Y
      INTEGER :: Ldy
!
! Local variable declarations rewritten by SPAG
!
      REAL :: ei
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
!     Quick return if possible
!
      IF ( N<=1 ) RETURN
!
      DO i = 1 , Nb
         IF ( i>1 ) THEN
!
!           Update A(1:n,i)
!
!           Compute i-th column of A - Y * V**T
!
            CALL SGEMV('No transpose',N,i-1,-ONE,Y,Ldy,A(K+i-1,1),Lda,  &
     &                 ONE,A(1,i),1)
!
!           Apply I - V * T**T * V**T to this column (call it b) from the
!           left, using the last column of T as workspace
!
!           Let  V = ( V1 )   and   b = ( b1 )   (first I-1 rows)
!                    ( V2 )             ( b2 )
!
!           where V1 is unit lower triangular
!
!           w := V1**T * b1
!
            CALL SCOPY(i-1,A(K+1,i),1,T(1,Nb),1)
            CALL STRMV('Lower','Transpose','Unit',i-1,A(K+1,1),Lda,     &
     &                 T(1,Nb),1)
!
!           w := w + V2**T *b2
!
            CALL SGEMV('Transpose',N-K-i+1,i-1,ONE,A(K+i,1),Lda,A(K+i,i)&
     &                 ,1,ONE,T(1,Nb),1)
!
!           w := T**T *w
!
            CALL STRMV('Upper','Transpose','Non-unit',i-1,T,Ldt,T(1,Nb),&
     &                 1)
!
!           b2 := b2 - V2*w
!
            CALL SGEMV('No transpose',N-K-i+1,i-1,-ONE,A(K+i,1),Lda,    &
     &                 T(1,Nb),1,ONE,A(K+i,i),1)
!
!           b1 := b1 - V1*w
!
            CALL STRMV('Lower','No transpose','Unit',i-1,A(K+1,1),Lda,  &
     &                 T(1,Nb),1)
            CALL SAXPY(i-1,-ONE,T(1,Nb),1,A(K+1,i),1)
!
            A(K+i-1,i-1) = ei
         ENDIF
!
!        Generate the elementary reflector H(i) to annihilate
!        A(k+i+1:n,i)
!
         CALL SLARFG(N-K-i+1,A(K+i,i),A(MIN(K+i+1,N),i),1,Tau(i))
         ei = A(K+i,i)
         A(K+i,i) = ONE
!
!        Compute  Y(1:n,i)
!
         CALL SGEMV('No transpose',N,N-K-i+1,ONE,A(1,i+1),Lda,A(K+i,i), &
     &              1,ZERO,Y(1,i),1)
         CALL SGEMV('Transpose',N-K-i+1,i-1,ONE,A(K+i,1),Lda,A(K+i,i),1,&
     &              ZERO,T(1,i),1)
         CALL SGEMV('No transpose',N,i-1,-ONE,Y,Ldy,T(1,i),1,ONE,Y(1,i),&
     &              1)
         CALL SSCAL(N,Tau(i),Y(1,i),1)
!
!        Compute T(1:i,i)
!
         CALL SSCAL(i-1,-Tau(i),T(1,i),1)
         CALL STRMV('Upper','No transpose','Non-unit',i-1,T,Ldt,T(1,i), &
     &              1)
         T(i,i) = Tau(i)
!
      ENDDO
      A(K+Nb,Nb) = ei
!
!
!     End of SLAHRD
!
      END SUBROUTINE SLAHRD
