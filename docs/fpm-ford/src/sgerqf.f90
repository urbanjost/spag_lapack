!*==sgerqf.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SGERQF
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SGERQF + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgerqf.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgerqf.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgerqf.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGERQF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SGERQF computes an RQ factorization of a real M-by-N matrix A:
!> A = R * Q.
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
!>          A is REAL array, dimension (LDA,N)
!>          On entry, the M-by-N matrix A.
!>          On exit,
!>          if m <= n, the upper triangle of the subarray
!>          A(1:m,n-m+1:n) contains the M-by-M upper triangular matrix R;
!>          if m >= n, the elements on and above the (m-n)-th subdiagonal
!>          contain the M-by-N upper trapezoidal matrix R;
!>          the remaining elements, with the array TAU, represent the
!>          orthogonal matrix Q as a product of min(m,n) elementary
!>          reflectors (see Further Details).
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
!>          TAU is REAL array, dimension (min(M,N))
!>          The scalar factors of the elementary reflectors (see Further
!>          Details).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.  LWORK >= max(1,M).
!>          For optimum performance LWORK >= M*NB, where NB is
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
!> \ingroup realGEcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The matrix Q is represented as a product of elementary reflectors
!>
!>     Q = H(1) H(2) . . . H(k), where k = min(m,n).
!>
!>  Each H(i) has the form
!>
!>     H(i) = I - tau * v * v**T
!>
!>  where tau is a real scalar, and v is a real vector with
!>  v(n-k+i+1:n) = 0 and v(n-k+i) = 1; v(1:n-k+i-1) is stored on exit in
!>  A(m-k+i,1:n-k+i-1), and tau in TAU(i).
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE SGERQF(M,N,A,Lda,Tau,Work,Lwork,Info)
      IMPLICIT NONE
!*--SGERQF142
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
      REAL A(Lda,*) , Tau(*) , Work(*)
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
      EXTERNAL SGERQ2 , SLARFB , SLARFT , XERBLA
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
      ELSEIF ( Lwork<MAX(1,M) .AND. .NOT.lquery ) THEN
         Info = -7
      ENDIF
!
      IF ( Info==0 ) THEN
         k = MIN(M,N)
         IF ( k==0 ) THEN
            lwkopt = 1
         ELSE
            nb = ILAENV(1,'SGERQF',' ',M,N,-1,-1)
            lwkopt = M*nb
            Work(1) = lwkopt
         ENDIF
         Work(1) = lwkopt
!
         IF ( Lwork<MAX(1,M) .AND. .NOT.lquery ) Info = -7
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('SGERQF',-Info)
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
      iws = M
      IF ( nb>1 .AND. nb<k ) THEN
!
!        Determine when to cross over from blocked to unblocked code.
!
         nx = MAX(0,ILAENV(3,'SGERQF',' ',M,N,-1,-1))
         IF ( nx<k ) THEN
!
!           Determine if workspace is large enough for blocked code.
!
            ldwork = M
            iws = ldwork*nb
            IF ( Lwork<iws ) THEN
!
!              Not enough workspace to use optimal NB:  reduce NB and
!              determine the minimum value of NB.
!
               nb = Lwork/ldwork
               nbmin = MAX(2,ILAENV(2,'SGERQF',' ',M,N,-1,-1))
            ENDIF
         ENDIF
      ENDIF
!
      IF ( nb>=nbmin .AND. nb<k .AND. nx<k ) THEN
!
!        Use blocked code initially.
!        The last kk rows are handled by the block method.
!
         ki = ((k-nx-1)/nb)*nb
         kk = MIN(k,ki+nb)
!
         DO i = k - kk + ki + 1 , k - kk + 1 , -nb
            ib = MIN(k-i+1,nb)
!
!           Compute the RQ factorization of the current block
!           A(m-k+i:m-k+i+ib-1,1:n-k+i+ib-1)
!
            CALL SGERQ2(ib,N-k+i+ib-1,A(M-k+i,1),Lda,Tau(i),Work,iinfo)
            IF ( M-k+i>1 ) THEN
!
!              Form the triangular factor of the block reflector
!              H = H(i+ib-1) . . . H(i+1) H(i)
!
               CALL SLARFT('Backward','Rowwise',N-k+i+ib-1,ib,A(M-k+i,1)&
     &                     ,Lda,Tau(i),Work,ldwork)
!
!              Apply H to A(1:m-k+i-1,1:n-k+i+ib-1) from the right
!
               CALL SLARFB('Right','No transpose','Backward','Rowwise', &
     &                     M-k+i-1,N-k+i+ib-1,ib,A(M-k+i,1),Lda,Work,   &
     &                     ldwork,A,Lda,Work(ib+1),ldwork)
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
      IF ( mu>0 .AND. nu>0 ) CALL SGERQ2(mu,nu,A,Lda,Tau,Work,iinfo)
!
      Work(1) = iws
!
!     End of SGERQF
!
      END SUBROUTINE SGERQF
