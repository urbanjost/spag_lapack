!*==ztzrqf.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZTZRQF
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZTZRQF + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztzrqf.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztzrqf.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztzrqf.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZTZRQF( M, N, A, LDA, TAU, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, M, N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * ), TAU( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> This routine is deprecated and has been replaced by routine ZTZRZF.
!>
!> ZTZRQF reduces the M-by-N ( M<=N ) complex upper trapezoidal matrix A
!> to upper triangular form by means of unitary transformations.
!>
!> The upper trapezoidal matrix A is factored as
!>
!>    A = ( R  0 ) * Z,
!>
!> where Z is an N-by-N unitary matrix and R is an M-by-M upper
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
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          On entry, the leading M-by-N upper trapezoidal part of the
!>          array A must contain the matrix to be factorized.
!>          On exit, the leading M-by-M upper triangular part of A
!>          contains the upper triangular matrix R, and elements M+1 to
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
!>          TAU is COMPLEX*16 array, dimension (M)
!>          The scalar factors of the elementary reflectors.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
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
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The  factorization is obtained by Householder's method.  The kth
!>  transformation matrix, Z( k ), whose conjugate transpose is used to
!>  introduce zeros into the (m - k + 1)th row of A, is given in the form
!>
!>     Z( k ) = ( I     0   ),
!>              ( 0  T( k ) )
!>
!>  where
!>
!>     T( k ) = I - tau*u( k )*u( k )**H,   u( k ) = (   1    ),
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
      SUBROUTINE ZTZRQF(M,N,A,Lda,Tau,Info)
      USE F77KINDS                        
      USE S_XERBLA
      USE S_ZAXPY
      USE S_ZCOPY
      USE S_ZGEMV
      USE S_ZGERC
      USE S_ZLACGV
      USE S_ZLARFG
      IMPLICIT NONE
!*--ZTZRQF150
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0) ,       &
     &                 CZERO = (0.0D+0,0.0D+0)
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Tau
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      COMPLEX(CX16KIND) :: alpha
      INTEGER :: i , k , m1
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
! =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. External Subroutines ..
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
         CALL XERBLA('ZTZRQF',-Info)
         RETURN
      ENDIF
!
!     Perform the factorization.
!
      IF ( M==0 ) RETURN
      IF ( M==N ) THEN
         DO i = 1 , N
            Tau(i) = CZERO
         ENDDO
      ELSE
         m1 = MIN(M+1,N)
         DO k = M , 1 , -1
!
!           Use a Householder reflection to zero the kth row of A.
!           First set up the reflection.
!
            A(k,k) = DCONJG(A(k,k))
            CALL ZLACGV(N-M,A(k,m1),Lda)
            alpha = A(k,k)
            CALL ZLARFG(N-M+1,alpha,A(k,m1),Lda,Tau(k))
            A(k,k) = alpha
            Tau(k) = DCONJG(Tau(k))
!
            IF ( Tau(k)/=CZERO .AND. k>1 ) THEN
!
!              We now perform the operation  A := A*P( k )**H.
!
!              Use the first ( k - 1 ) elements of TAU to store  a( k ),
!              where  a( k ) consists of the first ( k - 1 ) elements of
!              the  kth column  of  A.  Also  let  B  denote  the  first
!              ( k - 1 ) rows of the last ( n - m ) columns of A.
!
               CALL ZCOPY(k-1,A(1,k),1,Tau,1)
!
!              Form   w = a( k ) + B*z( k )  in TAU.
!
               CALL ZGEMV('No transpose',k-1,N-M,CONE,A(1,m1),Lda,      &
     &                    A(k,m1),Lda,CONE,Tau,1)
!
!              Now form  a( k ) := a( k ) - conjg(tau)*w
!              and       B      := B      - conjg(tau)*w*z( k )**H.
!
               CALL ZAXPY(k-1,-DCONJG(Tau(k)),Tau,1,A(1,k),1)
               CALL ZGERC(k-1,N-M,-DCONJG(Tau(k)),Tau,1,A(k,m1),Lda,    &
     &                    A(1,m1),Lda)
            ENDIF
         ENDDO
      ENDIF
!
!
!     End of ZTZRQF
!
      END SUBROUTINE ZTZRQF
