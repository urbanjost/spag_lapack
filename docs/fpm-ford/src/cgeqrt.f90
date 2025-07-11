!*==cgeqrt.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b CGEQRT
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CGEQRT + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgeqrt.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgeqrt.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgeqrt.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGEQRT( M, N, NB, A, LDA, T, LDT, WORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER INFO, LDA, LDT, M, N, NB
!       ..
!       .. Array Arguments ..
!       COMPLEX A( LDA, * ), T( LDT, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CGEQRT computes a blocked QR factorization of a complex M-by-N matrix A
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
!> \param[in] NB
!> \verbatim
!>          NB is INTEGER
!>          The block size to be used in the blocked QR.  MIN(M,N) >= NB >= 1.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          On entry, the M-by-N matrix A.
!>          On exit, the elements on and above the diagonal of the array
!>          contain the min(M,N)-by-N upper trapezoidal matrix R (R is
!>          upper triangular if M >= N); the elements below the diagonal
!>          are the columns of V.
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
!>          T is COMPLEX array, dimension (LDT,MIN(M,N))
!>          The upper triangular block reflectors stored in compact form
!>          as a sequence of upper triangular blocks.  See below
!>          for further details.
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is INTEGER
!>          The leading dimension of the array T.  LDT >= NB.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (NB*N)
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
!> \ingroup complexGEcomputational
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
!>  in the matrix A.  The 1's along the diagonal of V are not stored in A.
!>
!>  Let K=MIN(M,N).  The number of blocks is B = ceiling(K/NB), where each
!>  block is of order NB except for the last block, which is of order
!>  IB = K - (B-1)*NB.  For each of the B blocks, a upper triangular block
!>  reflector factor is computed: T1, T2, ..., TB.  The NB-by-NB (and IB-by-IB
!>  for the last block) T's are stored in the NB-by-K matrix T as
!>
!>               T = (T1 T2 ... TB).
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE CGEQRT(M,N,Nb,A,Lda,T,Ldt,Work,Info)
      IMPLICIT NONE
!*--CGEQRT145
!
!  -- LAPACK computational routine (version 3.7.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2017
!
!     .. Scalar Arguments ..
      INTEGER Info , Lda , Ldt , M , N , Nb
!     ..
!     .. Array Arguments ..
      COMPLEX A(Lda,*) , T(Ldt,*) , Work(*)
!     ..
!
! =====================================================================
!
!     ..
!     .. Local Scalars ..
      INTEGER i , ib , iinfo , k
      LOGICAL USE_RECURSIVE_QR
      PARAMETER (USE_RECURSIVE_QR=.TRUE.)
!     ..
!     .. External Subroutines ..
      EXTERNAL CGEQRT2 , CGEQRT3 , CLARFB , XERBLA
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
      ELSEIF ( Nb<1 .OR. (Nb>MIN(M,N) .AND. MIN(M,N)>0) ) THEN
         Info = -3
      ELSEIF ( Lda<MAX(1,M) ) THEN
         Info = -5
      ELSEIF ( Ldt<Nb ) THEN
         Info = -7
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('CGEQRT',-Info)
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
      DO i = 1 , k , Nb
         ib = MIN(k-i+1,Nb)
!
!     Compute the QR factorization of the current block A(I:M,I:I+IB-1)
!
         IF ( USE_RECURSIVE_QR ) THEN
            CALL CGEQRT3(M-i+1,ib,A(i,i),Lda,T(1,i),Ldt,iinfo)
         ELSE
            CALL CGEQRT2(M-i+1,ib,A(i,i),Lda,T(1,i),Ldt,iinfo)
         ENDIF
!
!     Update by applying H**H to A(I:M,I+IB:N) from the left
!
         IF ( i+ib<=N ) CALL CLARFB('L','C','F','C',M-i+1,N-i-ib+1,ib,  &
     &                              A(i,i),Lda,T(1,i),Ldt,A(i,i+ib),Lda,&
     &                              Work,N-i-ib+1)
      ENDDO
!
!     End of CGEQRT
!
      END SUBROUTINE CGEQRT
