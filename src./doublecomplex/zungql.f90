!*==zungql.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZUNGQL
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZUNGQL + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zungql.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zungql.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zungql.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZUNGQL( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, K, LDA, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZUNGQL generates an M-by-N complex matrix Q with orthonormal columns,
!> which is defined as the last N columns of a product of K elementary
!> reflectors of order M
!>
!>       Q  =  H(k) . . . H(2) H(1)
!>
!> as returned by ZGEQLF.
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
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          On entry, the (n-k+i)-th column must contain the vector which
!>          defines the elementary reflector H(i), for i = 1,2,...,k, as
!>          returned by ZGEQLF in the last k columns of its array
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
!>          TAU is COMPLEX*16 array, dimension (K)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by ZGEQLF.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
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
!> \ingroup complex16OTHERcomputational
!
!  =====================================================================
      SUBROUTINE ZUNGQL(M,N,K,A,Lda,Tau,Work,Lwork,Info)
      USE F77KINDS                        
      USE S_ILAENV
      USE S_XERBLA
      USE S_ZLARFB
      USE S_ZLARFT
      USE S_ZUNG2L
      IMPLICIT NONE
!*--ZUNGQL138
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , ib , iinfo , iws , j , kk , l , ldwork , lwkopt ,  &
     &           nb , nbmin , nx
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
      lquery = (Lwork==-1)
      IF ( M<0 ) THEN
         Info = -1
      ELSEIF ( N<0 .OR. N>M ) THEN
         Info = -2
      ELSEIF ( K<0 .OR. K>N ) THEN
         Info = -3
      ELSEIF ( Lda<MAX(1,M) ) THEN
         Info = -5
      ENDIF
!
      IF ( Info==0 ) THEN
         IF ( N==0 ) THEN
            lwkopt = 1
         ELSE
            nb = ILAENV(1,'ZUNGQL',' ',M,N,K,-1)
            lwkopt = N*nb
         ENDIF
         Work(1) = lwkopt
!
         IF ( Lwork<MAX(1,N) .AND. .NOT.lquery ) Info = -8
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZUNGQL',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N<=0 ) RETURN
!
      nbmin = 2
      nx = 0
      iws = N
      IF ( nb>1 .AND. nb<K ) THEN
!
!        Determine when to cross over from blocked to unblocked code.
!
         nx = MAX(0,ILAENV(3,'ZUNGQL',' ',M,N,K,-1))
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
               nbmin = MAX(2,ILAENV(2,'ZUNGQL',' ',M,N,K,-1))
            ENDIF
         ENDIF
      ENDIF
!
      IF ( nb>=nbmin .AND. nb<K .AND. nx<K ) THEN
!
!        Use blocked code after the first block.
!        The last kk columns are handled by the block method.
!
         kk = MIN(K,((K-nx+nb-1)/nb)*nb)
!
!        Set A(m-kk+1:m,1:n-kk) to zero.
!
         DO j = 1 , N - kk
            DO i = M - kk + 1 , M
               A(i,j) = ZERO
            ENDDO
         ENDDO
      ELSE
         kk = 0
      ENDIF
!
!     Use unblocked code for the first or only block.
!
      CALL ZUNG2L(M-kk,N-kk,K-kk,A,Lda,Tau,Work,iinfo)
!
      IF ( kk>0 ) THEN
!
!        Use blocked code
!
         DO i = K - kk + 1 , K , nb
            ib = MIN(nb,K-i+1)
            IF ( N-K+i>1 ) THEN
!
!              Form the triangular factor of the block reflector
!              H = H(i+ib-1) . . . H(i+1) H(i)
!
               CALL ZLARFT('Backward','Columnwise',M-K+i+ib-1,ib,       &
     &                     A(1,N-K+i),Lda,Tau(i),Work,ldwork)
!
!              Apply H to A(1:m-k+i+ib-1,1:n-k+i-1) from the left
!
               CALL ZLARFB('Left','No transpose','Backward',            &
     &                     'Columnwise',M-K+i+ib-1,N-K+i-1,ib,A(1,N-K+i)&
     &                     ,Lda,Work,ldwork,A,Lda,Work(ib+1),ldwork)
            ENDIF
!
!           Apply H to rows 1:m-k+i+ib-1 of current block
!
            CALL ZUNG2L(M-K+i+ib-1,ib,ib,A(1,N-K+i),Lda,Tau(i),Work,    &
     &                  iinfo)
!
!           Set rows m-k+i+ib:m of current block to zero
!
            DO j = N - K + i , N - K + i + ib - 1
               DO l = M - K + i + ib , M
                  A(l,j) = ZERO
               ENDDO
            ENDDO
         ENDDO
      ENDIF
!
      Work(1) = iws
!
!     End of ZUNGQL
!
      END SUBROUTINE ZUNGQL
