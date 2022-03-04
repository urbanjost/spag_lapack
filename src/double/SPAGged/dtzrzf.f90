!*==dtzrzf.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DTZRZF
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DTZRZF + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtzrzf.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtzrzf.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtzrzf.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DTZRZF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LWORK, M, N
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
!> DTZRZF reduces the M-by-N ( M<=N ) real upper trapezoidal matrix A
!> to upper triangular form by means of orthogonal transformations.
!>
!> The upper trapezoidal matrix A is factored as
!>
!>    A = ( R  0 ) * Z,
!>
!> where Z is an N-by-N orthogonal matrix and R is an M-by-M upper
!> triangular matrix.
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
!>          The number of columns of the matrix A.  N >= M.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, the leading M-by-N upper trapezoidal part of the
!>          array A must contain the matrix to be factorized.
!>          On exit, the leading M-by-M upper triangular part of A
!>          contains the upper triangular matrix R, and elements M+1 to
!>          N of the first M rows of A, with the array TAU, represent the
!>          orthogonal matrix Z as a product of M elementary reflectors.
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
!>          TAU is DOUBLE PRECISION array, dimension (M)
!>          The scalar factors of the elementary reflectors.
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
!> \date April 2012
!
!> \ingroup doubleOTHERcomputational
!
!> \par Contributors:
!  ==================
!>
!>    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The N-by-N matrix Z can be computed by
!>
!>     Z =  Z(1)*Z(2)* ... *Z(M)
!>
!>  where each N-by-N Z(k) is given by
!>
!>     Z(k) = I - tau(k)*v(k)*v(k)**T
!>
!>  with v(k) is the kth row vector of the M-by-N matrix
!>
!>     V = ( I   A(:,M+1:N) )
!>
!>  I is the M-by-M identity matrix, A(:,M+1:N)
!>  is the output stored in A on exit from DTZRZF,
!>  and tau(k) is the kth element of the array TAU.
!>
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE DTZRZF(M,N,A,Lda,Tau,Work,Lwork,Info)
      USE F77KINDS                        
      USE S_DLARZB
      USE S_DLARZT
      USE S_DLATRZ
      USE S_ILAENV
      USE S_XERBLA
      IMPLICIT NONE
!*--DTZRZF161
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER :: M
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , ib , iws , ki , kk , ldwork , lwkmin , lwkopt ,    &
     &           m1 , mu , nb , nbmin , nx
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
      ELSEIF ( N<M ) THEN
         Info = -2
      ELSEIF ( Lda<MAX(1,M) ) THEN
         Info = -4
      ENDIF
!
      IF ( Info==0 ) THEN
         IF ( M==0 .OR. M==N ) THEN
            lwkopt = 1
            lwkmin = 1
         ELSE
!
!           Determine the block size.
!
            nb = ILAENV(1,'DGERQF',' ',M,N,-1,-1)
            lwkopt = M*nb
            lwkmin = MAX(1,M)
         ENDIF
         Work(1) = lwkopt
!
         IF ( Lwork<lwkmin .AND. .NOT.lquery ) Info = -7
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('DTZRZF',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( M==0 ) THEN
         RETURN
      ELSEIF ( M==N ) THEN
         DO i = 1 , N
            Tau(i) = ZERO
         ENDDO
         RETURN
      ENDIF
!
      nbmin = 2
      nx = 1
      iws = M
      IF ( nb>1 .AND. nb<M ) THEN
!
!        Determine when to cross over from blocked to unblocked code.
!
         nx = MAX(0,ILAENV(3,'DGERQF',' ',M,N,-1,-1))
         IF ( nx<M ) THEN
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
               nbmin = MAX(2,ILAENV(2,'DGERQF',' ',M,N,-1,-1))
            ENDIF
         ENDIF
      ENDIF
!
      IF ( nb>=nbmin .AND. nb<M .AND. nx<M ) THEN
!
!        Use blocked code initially.
!        The last kk rows are handled by the block method.
!
         m1 = MIN(M+1,N)
         ki = ((M-nx-1)/nb)*nb
         kk = MIN(M,ki+nb)
!
         DO i = M - kk + ki + 1 , M - kk + 1 , -nb
            ib = MIN(M-i+1,nb)
!
!           Compute the TZ factorization of the current block
!           A(i:i+ib-1,i:n)
!
            CALL DLATRZ(ib,N-i+1,N-M,A(i,i),Lda,Tau(i),Work)
            IF ( i>1 ) THEN
!
!              Form the triangular factor of the block reflector
!              H = H(i+ib-1) . . . H(i+1) H(i)
!
               CALL DLARZT('Backward','Rowwise',N-M,ib,A(i,m1),Lda,     &
     &                     Tau(i),Work,ldwork)
!
!              Apply H to A(1:i-1,i:n) from the right
!
               CALL DLARZB('Right','No transpose','Backward','Rowwise', &
     &                     i-1,N-i+1,ib,N-M,A(i,m1),Lda,Work,ldwork,    &
     &                     A(1,i),Lda,Work(ib+1),ldwork)
            ENDIF
         ENDDO
         mu = i + nb - 1
      ELSE
         mu = M
      ENDIF
!
!     Use unblocked code to factor the last or only block
!
      IF ( mu>0 ) CALL DLATRZ(mu,N,N-M,A,Lda,Tau,Work)
!
      Work(1) = lwkopt
!
!
!     End of DTZRZF
!
      END SUBROUTINE DTZRZF
