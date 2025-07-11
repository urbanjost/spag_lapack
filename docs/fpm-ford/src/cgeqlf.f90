!*==cgeqlf.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b CGEQLF
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CGEQLF + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgeqlf.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgeqlf.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgeqlf.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGEQLF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       COMPLEX            A( LDA, * ), TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CGEQLF computes a QL factorization of a complex M-by-N matrix A:
!> A = Q * L.
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
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          On entry, the M-by-N matrix A.
!>          On exit,
!>          if m >= n, the lower triangle of the subarray
!>          A(m-n+1:m,1:n) contains the N-by-N lower triangular matrix L;
!>          if m <= n, the elements on and below the (n-m)-th
!>          superdiagonal contain the M-by-N lower trapezoidal matrix L;
!>          the remaining elements, with the array TAU, represent the
!>          unitary matrix Q as a product of elementary reflectors
!>          (see Further Details).
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is COMPLEX array, dimension (min(M,N))
!>          The scalar factors of the elementary reflectors (see Further
!>          Details).
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
!>          The dimension of the array WORK.  LWORK >= max(1,N).
!>          For optimum performance LWORK >= N*NB, where NB is
!>          the optimal blocksize.
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
!> \date December 2016
!
!> \ingroup complexGEcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The matrix Q is represented as a product of elementary reflectors
!>
!>     Q = H(k) . . . H(2) H(1), where k = min(m,n).
!>
!>  Each H(i) has the form
!>
!>     H(i) = I - tau * v * v**H
!>
!>  where tau is a complex scalar, and v is a complex vector with
!>  v(m-k+i+1:m) = 0 and v(m-k+i) = 1; v(1:m-k+i-1) is stored on exit in
!>  A(1:m-k+i-1,n-k+i), and tau in TAU(i).
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE CGEQLF(M,N,A,Lda,Tau,Work,Lwork,Info)
      IMPLICIT NONE
!*--CGEQLF142
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Info , Lda , Lwork , M , N
!     ..
!     .. Array Arguments ..
      COMPLEX A(Lda,*) , Tau(*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL lquery
      INTEGER i , ib , iinfo , iws , k , ki , kk , ldwork , lwkopt ,    &
     &        mu , nb , nbmin , nu , nx
!     ..
!     .. External Subroutines ..
      EXTERNAL CGEQL2 , CLARFB , CLARFT , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN
!     ..
!     .. External Functions ..
      INTEGER ILAENV
      EXTERNAL ILAENV
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      Info = 0
      lquery = (Lwork==-1)
      IF ( M<0 ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Lda<MAX(1,M) ) THEN
         Info = -4
      ENDIF
!
      IF ( Info==0 ) THEN
         k = MIN(M,N)
         IF ( k==0 ) THEN
            lwkopt = 1
         ELSE
            nb = ILAENV(1,'CGEQLF',' ',M,N,-1,-1)
            lwkopt = N*nb
         ENDIF
         Work(1) = lwkopt
!
         IF ( Lwork<MAX(1,N) .AND. .NOT.lquery ) Info = -7
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('CGEQLF',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( k==0 ) RETURN
!
      nbmin = 2
      nx = 1
      iws = N
      IF ( nb>1 .AND. nb<k ) THEN
!
!        Determine when to cross over from blocked to unblocked code.
!
         nx = MAX(0,ILAENV(3,'CGEQLF',' ',M,N,-1,-1))
         IF ( nx<k ) THEN
!
!           Determine if workspace is large enough for blocked code.
!
            ldwork = N
            iws = ldwork*nb
            IF ( Lwork<iws ) THEN
!
!              Not enough workspace to use optimal NB:  reduce NB and
!              determine the minimum value of NB.
!
               nb = Lwork/ldwork
               nbmin = MAX(2,ILAENV(2,'CGEQLF',' ',M,N,-1,-1))
            ENDIF
         ENDIF
      ENDIF
!
      IF ( nb>=nbmin .AND. nb<k .AND. nx<k ) THEN
!
!        Use blocked code initially.
!        The last kk columns are handled by the block method.
!
         ki = ((k-nx-1)/nb)*nb
         kk = MIN(k,ki+nb)
!
         DO i = k - kk + ki + 1 , k - kk + 1 , -nb
            ib = MIN(k-i+1,nb)
!
!           Compute the QL factorization of the current block
!           A(1:m-k+i+ib-1,n-k+i:n-k+i+ib-1)
!
            CALL CGEQL2(M-k+i+ib-1,ib,A(1,N-k+i),Lda,Tau(i),Work,iinfo)
            IF ( N-k+i>1 ) THEN
!
!              Form the triangular factor of the block reflector
!              H = H(i+ib-1) . . . H(i+1) H(i)
!
               CALL CLARFT('Backward','Columnwise',M-k+i+ib-1,ib,       &
     &                     A(1,N-k+i),Lda,Tau(i),Work,ldwork)
!
!              Apply H**H to A(1:m-k+i+ib-1,1:n-k+i-1) from the left
!
               CALL CLARFB('Left','Conjugate transpose','Backward',     &
     &                     'Columnwise',M-k+i+ib-1,N-k+i-1,ib,A(1,N-k+i)&
     &                     ,Lda,Work,ldwork,A,Lda,Work(ib+1),ldwork)
            ENDIF
         ENDDO
         mu = M - k + i + nb - 1
         nu = N - k + i + nb - 1
      ELSE
         mu = M
         nu = N
      ENDIF
!
!     Use unblocked code to factor the last or only block
!
      IF ( mu>0 .AND. nu>0 ) CALL CGEQL2(mu,nu,A,Lda,Tau,Work,iinfo)
!
      Work(1) = iws
!
!     End of CGEQLF
!
      END SUBROUTINE CGEQLF
