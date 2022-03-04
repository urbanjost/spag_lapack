!*==zgelqt.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b ZGELQT
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DGEQRT + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgelqt.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgelqt.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgelqt.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZGELQT( M, N, MB, A, LDA, T, LDT, WORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER INFO, LDA, LDT, M, N, MB
!       ..
!       .. Array Arguments ..
!       COMPLEX*16 A( LDA, * ), T( LDT, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZGELQT computes a blocked LQ factorization of a complex M-by-N matrix A
!> using the compact WY representation of Q.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] MB
!> \verbatim
!>          MB is INTEGER
!>          The block size to be used in the blocked QR.  MIN(M,N) >= MB >= 1.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          On entry, the M-by-N matrix A.
!>          On exit, the elements on and below the diagonal of the array
!>          contain the M-by-MIN(M,N) lower trapezoidal matrix L (L is
!>          lower triangular if M <= N); the elements above the diagonal
!>          are the rows of V.
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
!>          T is COMPLEX*16 array, dimension (LDT,MIN(M,N))
!>          The upper triangular block reflectors stored in compact form
!>          as a sequence of upper triangular blocks.  See below
!>          for further details.
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is INTEGER
!>          The leading dimension of the array T.  LDT >= MB.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (MB*N)
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
!> \date June 2017
!
!> \ingroup doubleGEcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The matrix V stores the elementary reflectors H(i) in the i-th row
!>  above the diagonal. For example, if M=5 and N=3, the matrix V is
!>
!>               V = (  1  v1 v1 v1 v1 )
!>                   (     1  v2 v2 v2 )
!>                   (         1 v3 v3 )
!>
!>
!>  where the vi's represent the vectors which define H(i), which are returned
!>  in the matrix A.  The 1's along the diagonal of V are not stored in A.
!>  Let K=MIN(M,N).  The number of blocks is B = ceiling(K/MB), where each
!>  block is of order MB except for the last block, which is of order
!>  IB = K - (B-1)*MB.  For each of the B blocks, a upper triangular block
!>  reflector factor is computed: T1, T2, ..., TB.  The MB-by-MB (and IB-by-IB
!>  for the last block) T's are stored in the MB-by-K matrix T as
!>
!>               T = (T1 T2 ... TB).
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE ZGELQT(M,N,Mb,A,Lda,T,Ldt,Work,Info)
      IMPLICIT NONE
!*--ZGELQT143
!
!  -- LAPACK computational routine (version 3.7.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2017
!
!     .. Scalar Arguments ..
      INTEGER Info , Lda , Ldt , M , N , Mb
!     ..
!     .. Array Arguments ..
      COMPLEX*16 A(Lda,*) , T(Ldt,*) , Work(*)
!     ..
!
! =====================================================================
!
!     ..
!     .. Local Scalars ..
      INTEGER i , ib , iinfo , k
!     ..
!     .. External Subroutines ..
      EXTERNAL ZGELQT3 , ZLARFB , XERBLA
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
      ELSEIF ( Mb<1 .OR. (Mb>MIN(M,N) .AND. MIN(M,N)>0) ) THEN
         Info = -3
      ELSEIF ( Lda<MAX(1,M) ) THEN
         Info = -5
      ELSEIF ( Ldt<Mb ) THEN
         Info = -7
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZGELQT',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      k = MIN(M,N)
      IF ( k==0 ) RETURN
!
!     Blocked loop of length K
!
      DO i = 1 , k , Mb
         ib = MIN(k-i+1,Mb)
!
!     Compute the LQ factorization of the current block A(I:M,I:I+IB-1)
!
         CALL ZGELQT3(ib,N-i+1,A(i,i),Lda,T(1,i),Ldt,iinfo)
!
!     Update by applying H**T to A(I:M,I+IB:N) from the right
!
         IF ( i+ib<=M ) CALL ZLARFB('R','N','F','R',M-i-ib+1,N-i+1,ib,  &
     &                              A(i,i),Lda,T(1,i),Ldt,A(i+ib,i),Lda,&
     &                              Work,M-i-ib+1)
      ENDDO
!
!     End of ZGELQT
!
      END SUBROUTINE ZGELQT
