!*==clatrz.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b CLATRZ factors an upper trapezoidal matrix by means of unitary transformations.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLATRZ + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clatrz.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clatrz.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clatrz.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLATRZ( M, N, L, A, LDA, TAU, WORK )
!
!       .. Scalar Arguments ..
!       INTEGER            L, LDA, M, N
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
!> CLATRZ factors the M-by-(M+L) complex upper trapezoidal matrix
!> [ A1 A2 ] = [ A(1:M,1:M) A(1:M,N-L+1:N) ] as ( R  0 ) * Z by means
!> of unitary transformations, where  Z is an (M+L)-by-(M+L) unitary
!> matrix and, R and A1 are M-by-M upper triangular matrices.
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
!> \param[in] L
!> \verbatim
!>          L is INTEGER
!>          The number of columns of the matrix A containing the
!>          meaningful part of the Householder vectors. N-M >= L >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          On entry, the leading M-by-N upper trapezoidal part of the
!>          array A must contain the matrix to be factorized.
!>          On exit, the leading M-by-M upper triangular part of A
!>          contains the upper triangular matrix R, and elements N-L+1 to
!>          N of the first M rows of A, with the array TAU, represent the
!>          unitary matrix Z as a product of M elementary reflectors.
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
!>          TAU is COMPLEX array, dimension (M)
!>          The scalar factors of the elementary reflectors.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (M)
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
!> \ingroup complexOTHERcomputational
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
!>  The factorization is obtained by Householder's method.  The kth
!>  transformation matrix, Z( k ), which is used to introduce zeros into
!>  the ( m - k + 1 )th row of A, is given in the form
!>
!>     Z( k ) = ( I     0   ),
!>              ( 0  T( k ) )
!>
!>  where
!>
!>     T( k ) = I - tau*u( k )*u( k )**H,   u( k ) = (   1    ),
!>                                                 (   0    )
!>                                                 ( z( k ) )
!>
!>  tau is a scalar and z( k ) is an l element vector. tau and z( k )
!>  are chosen to annihilate the elements of the kth row of A2.
!>
!>  The scalar tau is returned in the kth element of TAU and the vector
!>  u( k ) in the kth row of A2, such that the elements of z( k ) are
!>  in  a( k, l + 1 ), ..., a( k, n ). The elements of R are returned in
!>  the upper triangular part of A1.
!>
!>  Z is given by
!>
!>     Z =  Z( 1 ) * Z( 2 ) * ... * Z( m ).
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE CLATRZ(M,N,L,A,Lda,Tau,Work)
      IMPLICIT NONE
!*--CLATRZ144
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER L , Lda , M , N
!     ..
!     .. Array Arguments ..
      COMPLEX A(Lda,*) , Tau(*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX ZERO
      PARAMETER (ZERO=(0.0E+0,0.0E+0))
!     ..
!     .. Local Scalars ..
      INTEGER i
      COMPLEX alpha
!     ..
!     .. External Subroutines ..
      EXTERNAL CLACGV , CLARFG , CLARZ
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC CONJG
!     ..
!     .. Executable Statements ..
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
      DO i = M , 1 , -1
!
!        Generate elementary reflector H(i) to annihilate
!        [ A(i,i) A(i,n-l+1:n) ]
!
         CALL CLACGV(L,A(i,N-L+1),Lda)
         alpha = CONJG(A(i,i))
         CALL CLARFG(L+1,alpha,A(i,N-L+1),Lda,Tau(i))
         Tau(i) = CONJG(Tau(i))
!
!        Apply H(i) to A(1:i-1,i:n) from the right
!
         CALL CLARZ('Right',i-1,N-i+1,L,A(i,N-L+1),Lda,CONJG(Tau(i)),   &
     &              A(1,i),Lda,Work)
         A(i,i) = CONJG(alpha)
!
      ENDDO
!
!
!     End of CLATRZ
!
      END SUBROUTINE CLATRZ
