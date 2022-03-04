!*==zlahr2.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b ZLAHR2 reduces the specified number of first columns of a general rectangular matrix A so that elements below the specified subdiagonal are zero, and returns auxiliary matrices which are needed to apply the transformation to the unreduced part of A.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLAHR2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlahr2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlahr2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlahr2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLAHR2( N, K, NB, A, LDA, TAU, T, LDT, Y, LDY )
!
!       .. Scalar Arguments ..
!       INTEGER            K, LDA, LDT, LDY, N, NB
!       ..
!       .. Array Arguments ..
!       COMPLEX*16        A( LDA, * ), T( LDT, NB ), TAU( NB ),
!      $                   Y( LDY, NB )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLAHR2 reduces the first NB columns of A complex general n-BY-(n-k+1)
!> matrix A so that elements below the k-th subdiagonal are zero. The
!> reduction is performed by an unitary similarity transformation
!> Q**H * A * Q. The routine returns the matrices V and T which determine
!> Q as a block reflector I - V*T*V**H, and also the matrix Y = A * V * T.
!>
!> This is an auxiliary routine called by ZGEHRD.
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
!>          K < N.
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
!>          A is COMPLEX*16 array, dimension (LDA,N-K+1)
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
!>          TAU is COMPLEX*16 array, dimension (NB)
!>          The scalar factors of the elementary reflectors. See Further
!>          Details.
!> \endverbatim
!>
!> \param[out] T
!> \verbatim
!>          T is COMPLEX*16 array, dimension (LDT,NB)
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
!>          Y is COMPLEX*16 array, dimension (LDY,NB)
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
!> \ingroup complex16OTHERauxiliary
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
!>     H(i) = I - tau * v * v**H
!>
!>  where tau is a complex scalar, and v is a complex vector with
!>  v(1:i+k-1) = 0, v(i+k) = 1; v(i+k+1:n) is stored on exit in
!>  A(i+k+1:n,i), and tau in TAU(i).
!>
!>  The elements of the vectors v together form the (n-k+1)-by-nb matrix
!>  V which is needed, with T and Y, to apply the transformation to the
!>  unreduced part of the matrix, using an update of the form:
!>  A := (I - V*T*V**H) * (A - Y*V**H).
!>
!>  The contents of A on exit are illustrated by the following example
!>  with n = 7, k = 3 and nb = 2:
!>
!>     ( a   a   a   a   a )
!>     ( a   a   a   a   a )
!>     ( a   a   a   a   a )
!>     ( h   h   a   a   a )
!>     ( v1  h   a   a   a )
!>     ( v1  v2  a   a   a )
!>     ( v1  v2  a   a   a )
!>
!>  where a denotes an element of the original matrix A, h denotes a
!>  modified element of the upper Hessenberg matrix H, and vi denotes an
!>  element of the vector defining H(i).
!>
!>  This subroutine is a slight modification of LAPACK-3.0's DLAHRD
!>  incorporating improvements proposed by Quintana-Orti and Van de
!>  Gejin. Note that the entries of A(1:K,2:NB) differ from those
!>  returned by the original LAPACK-3.0's DLAHRD routine. (This
!>  subroutine is not backward compatible with LAPACK-3.0's DLAHRD.)
!> \endverbatim
!
!> \par References:
!  ================
!>
!>  Gregorio Quintana-Orti and Robert van de Geijn, "Improving the
!>  performance of reduction to Hessenberg form," ACM Transactions on
!>  Mathematical Software, 32(2):180-194, June 2006.
!>
!  =====================================================================
      SUBROUTINE ZLAHR2(N,K,Nb,A,Lda,Tau,T,Ldt,Y,Ldy)
      IMPLICIT NONE
!*--ZLAHR2185
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER K , Lda , Ldt , Ldy , N , Nb
!     ..
!     .. Array Arguments ..
      COMPLEX*16 A(Lda,*) , T(Ldt,Nb) , Tau(Nb) , Y(Ldy,Nb)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16 ZERO , ONE
      PARAMETER (ZERO=(0.0D+0,0.0D+0),ONE=(1.0D+0,0.0D+0))
!     ..
!     .. Local Scalars ..
      INTEGER i
      COMPLEX*16 ei
!     ..
!     .. External Subroutines ..
      EXTERNAL ZAXPY , ZCOPY , ZGEMM , ZGEMV , ZLACPY , ZLARFG , ZSCAL ,&
     &         ZTRMM , ZTRMV , ZLACGV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MIN
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
!           Update A(K+1:N,I)
!
!           Update I-th column of A - Y * V**H
!
            CALL ZLACGV(i-1,A(K+i-1,1),Lda)
            CALL ZGEMV('NO TRANSPOSE',N-K,i-1,-ONE,Y(K+1,1),Ldy,        &
     &                 A(K+i-1,1),Lda,ONE,A(K+1,i),1)
            CALL ZLACGV(i-1,A(K+i-1,1),Lda)
!
!           Apply I - V * T**H * V**H to this column (call it b) from the
!           left, using the last column of T as workspace
!
!           Let  V = ( V1 )   and   b = ( b1 )   (first I-1 rows)
!                    ( V2 )             ( b2 )
!
!           where V1 is unit lower triangular
!
!           w := V1**H * b1
!
            CALL ZCOPY(i-1,A(K+1,i),1,T(1,Nb),1)
            CALL ZTRMV('Lower','Conjugate transpose','UNIT',i-1,A(K+1,1)&
     &                 ,Lda,T(1,Nb),1)
!
!           w := w + V2**H * b2
!
            CALL ZGEMV('Conjugate transpose',N-K-i+1,i-1,ONE,A(K+i,1),  &
     &                 Lda,A(K+i,i),1,ONE,T(1,Nb),1)
!
!           w := T**H * w
!
            CALL ZTRMV('Upper','Conjugate transpose','NON-UNIT',i-1,T,  &
     &                 Ldt,T(1,Nb),1)
!
!           b2 := b2 - V2*w
!
            CALL ZGEMV('NO TRANSPOSE',N-K-i+1,i-1,-ONE,A(K+i,1),Lda,    &
     &                 T(1,Nb),1,ONE,A(K+i,i),1)
!
!           b1 := b1 - V1*w
!
            CALL ZTRMV('Lower','NO TRANSPOSE','UNIT',i-1,A(K+1,1),Lda,  &
     &                 T(1,Nb),1)
            CALL ZAXPY(i-1,-ONE,T(1,Nb),1,A(K+1,i),1)
!
            A(K+i-1,i-1) = ei
         ENDIF
!
!        Generate the elementary reflector H(I) to annihilate
!        A(K+I+1:N,I)
!
         CALL ZLARFG(N-K-i+1,A(K+i,i),A(MIN(K+i+1,N),i),1,Tau(i))
         ei = A(K+i,i)
         A(K+i,i) = ONE
!
!        Compute  Y(K+1:N,I)
!
         CALL ZGEMV('NO TRANSPOSE',N-K,N-K-i+1,ONE,A(K+1,i+1),Lda,      &
     &              A(K+i,i),1,ZERO,Y(K+1,i),1)
         CALL ZGEMV('Conjugate transpose',N-K-i+1,i-1,ONE,A(K+i,1),Lda, &
     &              A(K+i,i),1,ZERO,T(1,i),1)
         CALL ZGEMV('NO TRANSPOSE',N-K,i-1,-ONE,Y(K+1,1),Ldy,T(1,i),1,  &
     &              ONE,Y(K+1,i),1)
         CALL ZSCAL(N-K,Tau(i),Y(K+1,i),1)
!
!        Compute T(1:I,I)
!
         CALL ZSCAL(i-1,-Tau(i),T(1,i),1)
         CALL ZTRMV('Upper','No Transpose','NON-UNIT',i-1,T,Ldt,T(1,i), &
     &              1)
         T(i,i) = Tau(i)
!
      ENDDO
      A(K+Nb,Nb) = ei
!
!     Compute Y(1:K,1:NB)
!
      CALL ZLACPY('ALL',K,Nb,A(1,2),Lda,Y,Ldy)
      CALL ZTRMM('RIGHT','Lower','NO TRANSPOSE','UNIT',K,Nb,ONE,A(K+1,1)&
     &           ,Lda,Y,Ldy)
      IF ( N>K+Nb ) CALL ZGEMM('NO TRANSPOSE','NO TRANSPOSE',K,Nb,      &
     &                         N-K-Nb,ONE,A(1,2+Nb),Lda,A(K+1+Nb,1),Lda,&
     &                         ONE,Y,Ldy)
      CALL ZTRMM('RIGHT','Upper','NO TRANSPOSE','NON-UNIT',K,Nb,ONE,T,  &
     &           Ldt,Y,Ldy)
!
!
!     End of ZLAHR2
!
      END SUBROUTINE ZLAHR2
