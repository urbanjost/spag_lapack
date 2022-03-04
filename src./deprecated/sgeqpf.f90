!*==sgeqpf.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SGEQPF
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SGEQPF + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgeqpf.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgeqpf.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgeqpf.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGEQPF( M, N, A, LDA, JPVT, TAU, WORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, M, N
!       ..
!       .. Array Arguments ..
!       INTEGER            JPVT( * )
!       REAL               A( LDA, * ), TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> This routine is deprecated and has been replaced by routine SGEQP3.
!>
!> SGEQPF computes a QR factorization with column pivoting of a
!> real M-by-N matrix A: A*P = Q*R.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A. M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A. N >= 0
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>          On entry, the M-by-N matrix A.
!>          On exit, the upper triangle of the array contains the
!>          min(M,N)-by-N upper triangular matrix R; the elements
!>          below the diagonal, together with the array TAU,
!>          represent the orthogonal matrix Q as a product of
!>          min(m,n) elementary reflectors.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A. LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in,out] JPVT
!> \verbatim
!>          JPVT is INTEGER array, dimension (N)
!>          On entry, if JPVT(i) .ne. 0, the i-th column of A is permuted
!>          to the front of A*P (a leading column); if JPVT(i) = 0,
!>          the i-th column of A is a free column.
!>          On exit, if JPVT(i) = k, then the i-th column of A*P
!>          was the k-th column of A.
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is REAL array, dimension (min(M,N))
!>          The scalar factors of the elementary reflectors.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (3*N)
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
!>     Q = H(1) H(2) . . . H(n)
!>
!>  Each H(i) has the form
!>
!>     H = I - tau * v * v**T
!>
!>  where tau is a real scalar, and v is a real vector with
!>  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i).
!>
!>  The matrix P is represented in jpvt as follows: If
!>     jpvt(j) = i
!>  then the jth column of P is the ith canonical unit vector.
!>
!>  Partial column norm updating strategy modified by
!>    Z. Drmac and Z. Bujanovic, Dept. of Mathematics,
!>    University of Zagreb, Croatia.
!>  -- April 2011                                                      --
!>  For more details see LAPACK Working Note 176.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE SGEQPF(M,N,A,Lda,Jpvt,Tau,Work,Info)
      USE S_ISAMAX
      USE S_SGEQR2
      USE S_SLAMCH
      USE S_SLARF
      USE S_SLARFG
      USE S_SNRM2
      USE S_SORM2R
      USE S_SSWAP
      USE S_XERBLA
      IMPLICIT NONE
!*--SGEQPF155
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER :: M
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Jpvt
      REAL , DIMENSION(*) :: Tau
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL :: aii , temp , temp2 , tol3z
      INTEGER :: i , itemp , j , ma , mn , pvt
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
      IF ( M<0 ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Lda<MAX(1,M) ) THEN
         Info = -4
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('SGEQPF',-Info)
         RETURN
      ENDIF
!
      mn = MIN(M,N)
      tol3z = SQRT(SLAMCH('Epsilon'))
!
!     Move initial columns up front
!
      itemp = 1
      DO i = 1 , N
         IF ( Jpvt(i)/=0 ) THEN
            IF ( i/=itemp ) THEN
               CALL SSWAP(M,A(1,i),1,A(1,itemp),1)
               Jpvt(i) = Jpvt(itemp)
               Jpvt(itemp) = i
            ELSE
               Jpvt(i) = i
            ENDIF
            itemp = itemp + 1
         ELSE
            Jpvt(i) = i
         ENDIF
      ENDDO
      itemp = itemp - 1
!
!     Compute the QR factorization and update remaining columns
!
      IF ( itemp>0 ) THEN
         ma = MIN(itemp,M)
         CALL SGEQR2(M,ma,A,Lda,Tau,Work,Info)
         IF ( ma<N ) CALL SORM2R('Left','Transpose',M,N-ma,ma,A,Lda,Tau,&
     &                           A(1,ma+1),Lda,Work,Info)
      ENDIF
!
      IF ( itemp<mn ) THEN
!
!        Initialize partial column norms. The first n elements of
!        work store the exact column norms.
!
         DO i = itemp + 1 , N
            Work(i) = SNRM2(M-itemp,A(itemp+1,i),1)
            Work(N+i) = Work(i)
         ENDDO
!
!        Compute factorization
!
         DO i = itemp + 1 , mn
!
!           Determine ith pivot column and swap if necessary
!
            pvt = (i-1) + ISAMAX(N-i+1,Work(i),1)
!
            IF ( pvt/=i ) THEN
               CALL SSWAP(M,A(1,pvt),1,A(1,i),1)
               itemp = Jpvt(pvt)
               Jpvt(pvt) = Jpvt(i)
               Jpvt(i) = itemp
               Work(pvt) = Work(i)
               Work(N+pvt) = Work(N+i)
            ENDIF
!
!           Generate elementary reflector H(i)
!
            IF ( i<M ) THEN
               CALL SLARFG(M-i+1,A(i,i),A(i+1,i),1,Tau(i))
            ELSE
               CALL SLARFG(1,A(M,M),A(M,M),1,Tau(M))
            ENDIF
!
            IF ( i<N ) THEN
!
!              Apply H(i) to A(i:m,i+1:n) from the left
!
               aii = A(i,i)
               A(i,i) = ONE
               CALL SLARF('LEFT',M-i+1,N-i,A(i,i),1,Tau(i),A(i,i+1),Lda,&
     &                    Work(2*N+1))
               A(i,i) = aii
            ENDIF
!
!           Update partial column norms
!
            DO j = i + 1 , N
               IF ( Work(j)/=ZERO ) THEN
!
!                 NOTE: The following 4 lines follow from the analysis in
!                 Lapack Working Note 176.
!
                  temp = ABS(A(i,j))/Work(j)
                  temp = MAX(ZERO,(ONE+temp)*(ONE-temp))
                  temp2 = temp*(Work(j)/Work(N+j))**2
                  IF ( temp2>tol3z ) THEN
                     Work(j) = Work(j)*SQRT(temp)
                  ELSEIF ( M>i ) THEN
                     Work(j) = SNRM2(M-i,A(i+1,j),1)
                     Work(N+j) = Work(j)
                  ELSE
                     Work(j) = ZERO
                     Work(N+j) = ZERO
                  ENDIF
               ENDIF
            ENDDO
!
         ENDDO
      ENDIF
!
!     End of SGEQPF
!
      END SUBROUTINE SGEQPF
