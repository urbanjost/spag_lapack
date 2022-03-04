!*==cgeqpf.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b CGEQPF
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CGEQPF + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgeqpf.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgeqpf.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgeqpf.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGEQPF( M, N, A, LDA, JPVT, TAU, WORK, RWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, M, N
!       ..
!       .. Array Arguments ..
!       INTEGER            JPVT( * )
!       REAL               RWORK( * )
!       COMPLEX            A( LDA, * ), TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> This routine is deprecated and has been replaced by routine CGEQP3.
!>
!> CGEQPF computes a QR factorization with column pivoting of a
!> complex M-by-N matrix A: A*P = Q*R.
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
!>          A is COMPLEX array, dimension (LDA,N)
!>          On entry, the M-by-N matrix A.
!>          On exit, the upper triangle of the array contains the
!>          min(M,N)-by-N upper triangular matrix R; the elements
!>          below the diagonal, together with the array TAU,
!>          represent the unitary matrix Q as a product of
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
!>          TAU is COMPLEX array, dimension (min(M,N))
!>          The scalar factors of the elementary reflectors.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (N)
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (2*N)
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
!>     Q = H(1) H(2) . . . H(n)
!>
!>  Each H(i) has the form
!>
!>     H = I - tau * v * v**H
!>
!>  where tau is a complex scalar, and v is a complex vector with
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
      SUBROUTINE CGEQPF(M,N,A,Lda,Jpvt,Tau,Work,Rwork,Info)
      IMPLICIT NONE
!*--CGEQPF152
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Info , Lda , M , N
!     ..
!     .. Array Arguments ..
      INTEGER Jpvt(*)
      REAL Rwork(*)
      COMPLEX A(Lda,*) , Tau(*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
!     ..
!     .. Local Scalars ..
      INTEGER i , itemp , j , ma , mn , pvt
      REAL temp , temp2 , tol3z
      COMPLEX aii
!     ..
!     .. External Subroutines ..
      EXTERNAL CGEQR2 , CLARF , CLARFG , CSWAP , CUNM2R , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , CMPLX , CONJG , MAX , MIN , SQRT
!     ..
!     .. External Functions ..
      INTEGER ISAMAX
      REAL SCNRM2 , SLAMCH
      EXTERNAL ISAMAX , SCNRM2 , SLAMCH
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
         CALL XERBLA('CGEQPF',-Info)
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
               CALL CSWAP(M,A(1,i),1,A(1,itemp),1)
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
         CALL CGEQR2(M,ma,A,Lda,Tau,Work,Info)
         IF ( ma<N ) CALL CUNM2R('Left','Conjugate transpose',M,N-ma,ma,&
     &                           A,Lda,Tau,A(1,ma+1),Lda,Work,Info)
      ENDIF
!
      IF ( itemp<mn ) THEN
!
!        Initialize partial column norms. The first n elements of
!        work store the exact column norms.
!
         DO i = itemp + 1 , N
            Rwork(i) = SCNRM2(M-itemp,A(itemp+1,i),1)
            Rwork(N+i) = Rwork(i)
         ENDDO
!
!        Compute factorization
!
         DO i = itemp + 1 , mn
!
!           Determine ith pivot column and swap if necessary
!
            pvt = (i-1) + ISAMAX(N-i+1,Rwork(i),1)
!
            IF ( pvt/=i ) THEN
               CALL CSWAP(M,A(1,pvt),1,A(1,i),1)
               itemp = Jpvt(pvt)
               Jpvt(pvt) = Jpvt(i)
               Jpvt(i) = itemp
               Rwork(pvt) = Rwork(i)
               Rwork(N+pvt) = Rwork(N+i)
            ENDIF
!
!           Generate elementary reflector H(i)
!
            aii = A(i,i)
            CALL CLARFG(M-i+1,aii,A(MIN(i+1,M),i),1,Tau(i))
            A(i,i) = aii
!
            IF ( i<N ) THEN
!
!              Apply H(i) to A(i:m,i+1:n) from the left
!
               aii = A(i,i)
               A(i,i) = CMPLX(ONE)
               CALL CLARF('Left',M-i+1,N-i,A(i,i),1,CONJG(Tau(i)),      &
     &                    A(i,i+1),Lda,Work)
               A(i,i) = aii
            ENDIF
!
!           Update partial column norms
!
            DO j = i + 1 , N
               IF ( Rwork(j)/=ZERO ) THEN
!
!                 NOTE: The following 4 lines follow from the analysis in
!                 Lapack Working Note 176.
!
                  temp = ABS(A(i,j))/Rwork(j)
                  temp = MAX(ZERO,(ONE+temp)*(ONE-temp))
                  temp2 = temp*(Rwork(j)/Rwork(N+j))**2
                  IF ( temp2>tol3z ) THEN
                     Rwork(j) = Rwork(j)*SQRT(temp)
                  ELSEIF ( M>i ) THEN
                     Rwork(j) = SCNRM2(M-i,A(i+1,j),1)
                     Rwork(N+j) = Rwork(j)
                  ELSE
                     Rwork(j) = ZERO
                     Rwork(N+j) = ZERO
                  ENDIF
               ENDIF
            ENDDO
!
         ENDDO
      ENDIF
!
!     End of CGEQPF
!
      END SUBROUTINE CGEQPF
