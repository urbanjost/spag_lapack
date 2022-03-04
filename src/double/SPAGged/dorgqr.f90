!*==dorgqr.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DORGQR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DORGQR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorgqr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorgqr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorgqr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, K, LDA, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DORGQR generates an M-by-N real matrix Q with orthonormal columns,
!> which is defined as the first N columns of a product of K elementary
!> reflectors of order M
!>
!>       Q  =  H(1) H(2) . . . H(k)
!>
!> as returned by DGEQRF.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix Q. M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix Q. M >= N >= 0.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The number of elementary reflectors whose product defines the
!>          matrix Q. N >= K >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, the i-th column must contain the vector which
!>          defines the elementary reflector H(i), for i = 1,2,...,k, as
!>          returned by DGEQRF in the first k columns of its array
!>          argument A.
!>          On exit, the M-by-N matrix Q.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The first dimension of the array A. LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is DOUBLE PRECISION array, dimension (K)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by DGEQRF.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK. LWORK >= max(1,N).
!>          For optimum performance LWORK >= N*NB, where NB is the
!>          optimal blocksize.
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
!>          < 0:  if INFO = -i, the i-th argument has an illegal value
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
!> \ingroup doubleOTHERcomputational
!
!  =====================================================================
      SUBROUTINE DORGQR(M,N,K,A,Lda,Tau,Work,Lwork,Info)
      USE F77KINDS                        
      USE S_DLARFB
      USE S_DLARFT
      USE S_DORG2R
      USE S_ILAENV
      USE S_XERBLA
      IMPLICIT NONE
!*--DORGQR138
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , ib , iinfo , iws , j , ki , kk , l , ldwork ,      &
     &           lwkopt , nb , nbmin , nx
      LOGICAL :: lquery
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
!     .. External Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      Info = 0
      nb = ILAENV(1,'DORGQR',' ',M,N,K,-1)
      lwkopt = MAX(1,N)*nb
      Work(1) = lwkopt
      lquery = (Lwork==-1)
      IF ( M<0 ) THEN
         Info = -1
      ELSEIF ( N<0 .OR. N>M ) THEN
         Info = -2
      ELSEIF ( K<0 .OR. K>N ) THEN
         Info = -3
      ELSEIF ( Lda<MAX(1,M) ) THEN
         Info = -5
      ELSEIF ( Lwork<MAX(1,N) .AND. .NOT.lquery ) THEN
         Info = -8
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('DORGQR',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N<=0 ) THEN
         Work(1) = 1
         RETURN
      ENDIF
!
      nbmin = 2
      nx = 0
      iws = N
      IF ( nb>1 .AND. nb<K ) THEN
!
!        Determine when to cross over from blocked to unblocked code.
!
         nx = MAX(0,ILAENV(3,'DORGQR',' ',M,N,K,-1))
         IF ( nx<K ) THEN
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
               nbmin = MAX(2,ILAENV(2,'DORGQR',' ',M,N,K,-1))
            ENDIF
         ENDIF
      ENDIF
!
      IF ( nb>=nbmin .AND. nb<K .AND. nx<K ) THEN
!
!        Use blocked code after the last block.
!        The first kk columns are handled by the block method.
!
         ki = ((K-nx-1)/nb)*nb
         kk = MIN(K,ki+nb)
!
!        Set A(1:kk,kk+1:n) to zero.
!
         DO j = kk + 1 , N
            DO i = 1 , kk
               A(i,j) = ZERO
            ENDDO
         ENDDO
      ELSE
         kk = 0
      ENDIF
!
!     Use unblocked code for the last or only block.
!
      IF ( kk<N ) CALL DORG2R(M-kk,N-kk,K-kk,A(kk+1,kk+1),Lda,Tau(kk+1),&
     &                        Work,iinfo)
!
      IF ( kk>0 ) THEN
!
!        Use blocked code
!
         DO i = ki + 1 , 1 , -nb
            ib = MIN(nb,K-i+1)
            IF ( i+ib<=N ) THEN
!
!              Form the triangular factor of the block reflector
!              H = H(i) H(i+1) . . . H(i+ib-1)
!
               CALL DLARFT('Forward','Columnwise',M-i+1,ib,A(i,i),Lda,  &
     &                     Tau(i),Work,ldwork)
!
!              Apply H to A(i:m,i+ib:n) from the left
!
               CALL DLARFB('Left','No transpose','Forward','Columnwise',&
     &                     M-i+1,N-i-ib+1,ib,A(i,i),Lda,Work,ldwork,    &
     &                     A(i,i+ib),Lda,Work(ib+1),ldwork)
            ENDIF
!
!           Apply H to rows i:m of current block
!
            CALL DORG2R(M-i+1,ib,ib,A(i,i),Lda,Tau(i),Work,iinfo)
!
!           Set rows 1:i-1 of current block to zero
!
            DO j = i , i + ib - 1
               DO l = 1 , i - 1
                  A(l,j) = ZERO
               ENDDO
            ENDDO
         ENDDO
      ENDIF
!
      Work(1) = iws
!
!     End of DORGQR
!
      END SUBROUTINE DORGQR
