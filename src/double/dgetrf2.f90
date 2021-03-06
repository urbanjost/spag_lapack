!*==dgetrf2.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b DGETRF2
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       RECURSIVE SUBROUTINE DGETRF2( M, N, A, LDA, IPIV, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, M, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       DOUBLE PRECISION   A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGETRF2 computes an LU factorization of a general M-by-N matrix A
!> using partial pivoting with row interchanges.
!>
!> The factorization has the form
!>    A = P * L * U
!> where P is a permutation matrix, L is lower triangular with unit
!> diagonal elements (lower trapezoidal if m > n), and U is upper
!> triangular (upper trapezoidal if m < n).
!>
!> This is the recursive version of the algorithm. It divides
!> the matrix into four submatrices:
!>
!>        [  A11 | A12  ]  where A11 is n1 by n1 and A22 is n2 by n2
!>    A = [ -----|----- ]  with n1 = min(m,n)/2
!>        [  A21 | A22  ]       n2 = n-n1
!>
!>                                       [ A11 ]
!> The subroutine calls itself to factor [ --- ],
!>                                       [ A12 ]
!>                 [ A12 ]
!> do the swaps on [ --- ], solve A12, update A22,
!>                 [ A22 ]
!>
!> then calls itself to factor A22 and do the swaps on A21.
!>
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
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, the M-by-N matrix to be factored.
!>          On exit, the factors L and U from the factorization
!>          A = P*L*U; the unit diagonal elements of L are not stored.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[out] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (min(M,N))
!>          The pivot indices; for 1 <= i <= min(M,N), row i of the
!>          matrix was interchanged with row IPIV(i).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
!>                has been completed, but the factor U is exactly
!>                singular, and division by zero will occur if it is used
!>                to solve a system of equations.
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
!> \date June 2016
!
!> \ingroup doubleGEcomputational
!
!  =====================================================================
      RECURSIVE SUBROUTINE DGETRF2(M,N,A,Lda,Ipiv,Info)
      IMPLICIT NONE
!*--DGETRF2117
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!
!     .. Scalar Arguments ..
      INTEGER Info , Lda , M , N
!     ..
!     .. Array Arguments ..
      INTEGER Ipiv(*)
      DOUBLE PRECISION A(Lda,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ONE , ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION sfmin , temp
      INTEGER i , iinfo , n1 , n2
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH
      INTEGER IDAMAX
      EXTERNAL DLAMCH , IDAMAX
!     ..
!     .. External Subroutines ..
      EXTERNAL DGEMM , DSCAL , DLASWP , DTRSM , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      Info = 0
      IF ( M<0 ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Lda<MAX(1,M) ) THEN
         Info = -4
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('DGETRF2',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( M==0 .OR. N==0 ) RETURN
 
      IF ( M==1 ) THEN
!
!        Use unblocked code for one row case
!        Just need to handle IPIV and INFO
!
         Ipiv(1) = 1
         IF ( A(1,1)==ZERO ) Info = 1
!
      ELSEIF ( N==1 ) THEN
!
!        Use unblocked code for one column case
!
!
!        Compute machine safe minimum
!
         sfmin = DLAMCH('S')
!
!        Find pivot and test for singularity
!
         i = IDAMAX(M,A(1,1),1)
         Ipiv(1) = i
         IF ( A(i,1)/=ZERO ) THEN
!
!           Apply the interchange
!
            IF ( i/=1 ) THEN
               temp = A(1,1)
               A(1,1) = A(i,1)
               A(i,1) = temp
            ENDIF
!
!           Compute elements 2:M of the column
!
            IF ( ABS(A(1,1))>=sfmin ) THEN
               CALL DSCAL(M-1,ONE/A(1,1),A(2,1),1)
            ELSE
               DO i = 1 , M - 1
                  A(1+i,1) = A(1+i,1)/A(1,1)
               ENDDO
            ENDIF
!
         ELSE
            Info = 1
         ENDIF
!
      ELSE
!
!        Use recursive code
!
         n1 = MIN(M,N)/2
         n2 = N - n1
!
!               [ A11 ]
!        Factor [ --- ]
!               [ A21 ]
!
         CALL DGETRF2(M,n1,A,Lda,Ipiv,iinfo)
 
         IF ( Info==0 .AND. iinfo>0 ) Info = iinfo
!
!                              [ A12 ]
!        Apply interchanges to [ --- ]
!                              [ A22 ]
!
         CALL DLASWP(n2,A(1,n1+1),Lda,1,n1,Ipiv,1)
!
!        Solve A12
!
         CALL DTRSM('L','L','N','U',n1,n2,ONE,A,Lda,A(1,n1+1),Lda)
!
!        Update A22
!
         CALL DGEMM('N','N',M-n1,n2,n1,-ONE,A(n1+1,1),Lda,A(1,n1+1),Lda,&
     &              ONE,A(n1+1,n1+1),Lda)
!
!        Factor A22
!
         CALL DGETRF2(M-n1,n2,A(n1+1,n1+1),Lda,Ipiv(n1+1),iinfo)
!
!        Adjust INFO and the pivot indices
!
         IF ( Info==0 .AND. iinfo>0 ) Info = iinfo + n1
         DO i = n1 + 1 , MIN(M,N)
            Ipiv(i) = Ipiv(i) + n1
         ENDDO
!
!        Apply interchanges to A21
!
         CALL DLASWP(n1,A(1,1),Lda,n1+1,MIN(M,N),Ipiv,1)
!
      ENDIF
!
!     End of DGETRF2
!
      END SUBROUTINE DGETRF2
