!*==cgetri.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CGETRI
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CGETRI + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgetri.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgetri.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgetri.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LWORK, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX            A( LDA, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CGETRI computes the inverse of a matrix using the LU factorization
!> computed by CGETRF.
!>
!> This method inverts U and then computes inv(A) by solving the system
!> inv(A)*L = inv(U) for inv(A).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          On entry, the factors L and U from the factorization
!>          A = P*L*U as computed by CGETRF.
!>          On exit, if INFO = 0, the inverse of the original matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          The pivot indices from CGETRF; for 1<=i<=N, row i of the
!>          matrix was interchanged with row IPIV(i).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (MAX(1,LWORK))
!>          On exit, if INFO=0, then WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.  LWORK >= max(1,N).
!>          For optimal performance LWORK >= N*NB, where NB is
!>          the optimal blocksize returned by ILAENV.
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
!>          > 0:  if INFO = i, U(i,i) is exactly zero; the matrix is
!>                singular and its inverse could not be computed.
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
!  =====================================================================
      SUBROUTINE CGETRI(N,A,Lda,Ipiv,Work,Lwork,Info)
      USE S_CGEMM
      USE S_CGEMV
      USE S_CSWAP
      USE S_CTRSM
      USE S_CTRTRI
      USE S_ILAENV
      USE S_XERBLA
      IMPLICIT NONE
!*--CGETRI125
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0) ,                 &
     &                         ONE = (1.0E+0,0.0E+0)
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , iws , j , jb , jj , jp , ldwork , lwkopt , nb ,    &
     &           nbmin , nn
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
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
      nb = ILAENV(1,'CGETRI',' ',N,-1,-1,-1)
      lwkopt = N*nb
      Work(1) = lwkopt
      lquery = (Lwork==-1)
      IF ( N<0 ) THEN
         Info = -1
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -3
      ELSEIF ( Lwork<MAX(1,N) .AND. .NOT.lquery ) THEN
         Info = -6
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('CGETRI',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
!     Form inv(U).  If INFO > 0 from CTRTRI, then U is singular,
!     and the inverse is not computed.
!
      CALL CTRTRI('Upper','Non-unit',N,A,Lda,Info)
      IF ( Info>0 ) RETURN
!
      nbmin = 2
      ldwork = N
      IF ( nb>1 .AND. nb<N ) THEN
         iws = MAX(ldwork*nb,1)
         IF ( Lwork<iws ) THEN
            nb = Lwork/ldwork
            nbmin = MAX(2,ILAENV(2,'CGETRI',' ',N,-1,-1,-1))
         ENDIF
      ELSE
         iws = N
      ENDIF
!
!     Solve the equation inv(A)*L = inv(U) for inv(A).
!
      IF ( nb<nbmin .OR. nb>=N ) THEN
!
!        Use unblocked code.
!
         DO j = N , 1 , -1
!
!           Copy current column of L to WORK and replace with zeros.
!
            DO i = j + 1 , N
               Work(i) = A(i,j)
               A(i,j) = ZERO
            ENDDO
!
!           Compute current column of inv(A).
!
            IF ( j<N ) CALL CGEMV('No transpose',N,N-j,-ONE,A(1,j+1),   &
     &                            Lda,Work(j+1),1,ONE,A(1,j),1)
         ENDDO
      ELSE
!
!        Use blocked code.
!
         nn = ((N-1)/nb)*nb + 1
         DO j = nn , 1 , -nb
            jb = MIN(nb,N-j+1)
!
!           Copy current block column of L to WORK and replace with
!           zeros.
!
            DO jj = j , j + jb - 1
               DO i = jj + 1 , N
                  Work(i+(jj-j)*ldwork) = A(i,jj)
                  A(i,jj) = ZERO
               ENDDO
            ENDDO
!
!           Compute current block column of inv(A).
!
            IF ( j+jb<=N ) CALL CGEMM('No transpose','No transpose',N,  &
     &                                jb,N-j-jb+1,-ONE,A(1,j+jb),Lda,   &
     &                                Work(j+jb),ldwork,ONE,A(1,j),Lda)
            CALL CTRSM('Right','Lower','No transpose','Unit',N,jb,ONE,  &
     &                 Work(j),ldwork,A(1,j),Lda)
         ENDDO
      ENDIF
!
!     Apply column interchanges.
!
      DO j = N - 1 , 1 , -1
         jp = Ipiv(j)
         IF ( jp/=j ) CALL CSWAP(N,A(1,j),1,A(1,jp),1)
      ENDDO
!
      Work(1) = iws
!
!     End of CGETRI
!
      END SUBROUTINE CGETRI
