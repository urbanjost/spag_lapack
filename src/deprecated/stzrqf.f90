!*==stzrqf.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b STZRQF
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download STZRQF + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stzrqf.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stzrqf.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stzrqf.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE STZRQF( M, N, A, LDA, TAU, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, M, N
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), TAU( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> This routine is deprecated and has been replaced by routine STZRZF.
!>
!> STZRQF reduces the M-by-N ( M<=N ) real upper trapezoidal matrix A
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
!>          A is REAL array, dimension (LDA,N)
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
!>          TAU is REAL array, dimension (M)
!>          The scalar factors of the elementary reflectors.
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
!> \ingroup realOTHERcomputational
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
!>     T( k ) = I - tau*u( k )*u( k )**T,   u( k ) = (   1    ),
!>                                                   (   0    )
!>                                                   ( z( k ) )
!>
!>  tau is a scalar and z( k ) is an ( n - m ) element vector.
!>  tau and z( k ) are chosen to annihilate the elements of the kth row
!>  of X.
!>
!>  The scalar tau is returned in the kth element of TAU and the vector
!>  u( k ) in the kth row of A, such that the elements of z( k ) are
!>  in  a( k, m + 1 ), ..., a( k, n ). The elements of R are returned in
!>  the upper triangular part of A.
!>
!>  Z is given by
!>
!>     Z =  Z( 1 ) * Z( 2 ) * ... * Z( m ).
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE STZRQF(M,N,A,Lda,Tau,Info)
      IMPLICIT NONE
!*--STZRQF142
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
      REAL A(Lda,*) , Tau(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ONE , ZERO
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
!     ..
!     .. Local Scalars ..
      INTEGER i , k , m1
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN
!     ..
!     .. External Subroutines ..
      EXTERNAL SAXPY , SCOPY , SGEMV , SGER , SLARFG , XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
      IF ( M<0 ) THEN
         Info = -1
      ELSEIF ( N<M ) THEN
         Info = -2
      ELSEIF ( Lda<MAX(1,M) ) THEN
         Info = -4
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('STZRQF',-Info)
         RETURN
      ENDIF
!
!     Perform the factorization.
!
      IF ( M==0 ) RETURN
      IF ( M==N ) THEN
         DO i = 1 , N
            Tau(i) = ZERO
         ENDDO
      ELSE
         m1 = MIN(M+1,N)
         DO k = M , 1 , -1
!
!           Use a Householder reflection to zero the kth row of A.
!           First set up the reflection.
!
            CALL SLARFG(N-M+1,A(k,k),A(k,m1),Lda,Tau(k))
!
            IF ( (Tau(k)/=ZERO) .AND. (k>1) ) THEN
!
!              We now perform the operation  A := A*P( k ).
!
!              Use the first ( k - 1 ) elements of TAU to store  a( k ),
!              where  a( k ) consists of the first ( k - 1 ) elements of
!              the  kth column  of  A.  Also  let  B  denote  the  first
!              ( k - 1 ) rows of the last ( n - m ) columns of A.
!
               CALL SCOPY(k-1,A(1,k),1,Tau,1)
!
!              Form   w = a( k ) + B*z( k )  in TAU.
!
               CALL SGEMV('No transpose',k-1,N-M,ONE,A(1,m1),Lda,A(k,m1)&
     &                    ,Lda,ONE,Tau,1)
!
!              Now form  a( k ) := a( k ) - tau*w
!              and       B      := B      - tau*w*z( k )**T.
!
               CALL SAXPY(k-1,-Tau(k),Tau,1,A(1,k),1)
               CALL SGER(k-1,N-M,-Tau(k),Tau,1,A(k,m1),Lda,A(1,m1),Lda)
            ENDIF
         ENDDO
      ENDIF
!
!
!     End of STZRQF
!
      END SUBROUTINE STZRQF
