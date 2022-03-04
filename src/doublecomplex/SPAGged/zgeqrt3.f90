!*==zgeqrt3.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZGEQRT3 recursively computes a QR factorization of a general real or complex matrix using the compact WY representation of Q.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZGEQRT3 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgeqrt3.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgeqrt3.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgeqrt3.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       RECURSIVE SUBROUTINE ZGEQRT3( M, N, A, LDA, T, LDT, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER   INFO, LDA, M, N, LDT
!       ..
!       .. Array Arguments ..
!       COMPLEX*16   A( LDA, * ), T( LDT, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZGEQRT3 recursively computes a QR factorization of a complex M-by-N
!> matrix A, using the compact WY representation of Q.
!>
!> Based on the algorithm of Elmroth and Gustavson,
!> IBM J. Res. Develop. Vol 44 No. 4 July 2000.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= N.
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
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          On entry, the complex M-by-N matrix A.  On exit, the elements on
!>          and above the diagonal contain the N-by-N upper triangular matrix R;
!>          the elements below the diagonal are the columns of V.  See below for
!>          further details.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[out] T
!> \verbatim
!>          T is COMPLEX*16 array, dimension (LDT,N)
!>          The N-by-N upper triangular factor of the block reflector.
!>          The elements on and above the diagonal contain the block
!>          reflector T; the elements below the diagonal are not used.
!>          See below for further details.
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is INTEGER
!>          The leading dimension of the array T.  LDT >= max(1,N).
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
!> \date June 2016
!
!> \ingroup complex16GEcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The matrix V stores the elementary reflectors H(i) in the i-th column
!>  below the diagonal. For example, if M=5 and N=3, the matrix V is
!>
!>               V = (  1       )
!>                   ( v1  1    )
!>                   ( v1 v2  1 )
!>                   ( v1 v2 v3 )
!>                   ( v1 v2 v3 )
!>
!>  where the vi's represent the vectors which define H(i), which are returned
!>  in the matrix A.  The 1's along the diagonal of V are not stored in A.  The
!>  block reflector H is then given by
!>
!>               H = I - V * T * V**H
!>
!>  where V**H is the conjugate transpose of V.
!>
!>  For details of the algorithm, see Elmroth and Gustavson (cited above).
!> \endverbatim
!>
!  =====================================================================
      RECURSIVE SUBROUTINE ZGEQRT3(M,N,A,Lda,T,Ldt,Info)
      USE F77KINDS                        
      USE S_XERBLA
      USE S_ZGEMM
      USE S_ZLARFG
      USE S_ZTRMM
      IMPLICIT NONE
!*--ZGEQRT3141
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+00,0.0D+00)
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , i1 , iinfo , j , j1 , n1 , n2
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
!     .. Executable Statements ..
!
      Info = 0
      IF ( N<0 ) THEN
         Info = -2
      ELSEIF ( M<N ) THEN
         Info = -1
      ELSEIF ( Lda<MAX(1,M) ) THEN
         Info = -4
      ELSEIF ( Ldt<MAX(1,N) ) THEN
         Info = -6
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZGEQRT3',-Info)
         RETURN
      ENDIF
!
      IF ( N==1 ) THEN
!
!        Compute Householder transform when N=1
!
         CALL ZLARFG(M,A(1,1),A(MIN(2,M),1),1,T(1,1))
!
      ELSE
!
!        Otherwise, split A into blocks...
!
         n1 = N/2
         n2 = N - n1
         j1 = MIN(n1+1,N)
         i1 = MIN(N+1,M)
!
!        Compute A(1:M,1:N1) <- (Y1,R1,T1), where Q1 = I - Y1 T1 Y1^H
!
         CALL ZGEQRT3(M,n1,A,Lda,T,Ldt,iinfo)
!
!        Compute A(1:M,J1:N) = Q1^H A(1:M,J1:N) [workspace: T(1:N1,J1:N)]
!
         DO j = 1 , n2
            DO i = 1 , n1
               T(i,j+n1) = A(i,j+n1)
            ENDDO
         ENDDO
         CALL ZTRMM('L','L','C','U',n1,n2,ONE,A,Lda,T(1,j1),Ldt)
!
         CALL ZGEMM('C','N',n1,n2,M-n1,ONE,A(j1,1),Lda,A(j1,j1),Lda,ONE,&
     &              T(1,j1),Ldt)
!
         CALL ZTRMM('L','U','C','N',n1,n2,ONE,T,Ldt,T(1,j1),Ldt)
!
         CALL ZGEMM('N','N',M-n1,n2,n1,-ONE,A(j1,1),Lda,T(1,j1),Ldt,ONE,&
     &              A(j1,j1),Lda)
!
         CALL ZTRMM('L','L','N','U',n1,n2,ONE,A,Lda,T(1,j1),Ldt)
!
         DO j = 1 , n2
            DO i = 1 , n1
               A(i,j+n1) = A(i,j+n1) - T(i,j+n1)
            ENDDO
         ENDDO
!
!        Compute A(J1:M,J1:N) <- (Y2,R2,T2) where Q2 = I - Y2 T2 Y2^H
!
         CALL ZGEQRT3(M-n1,n2,A(j1,j1),Lda,T(j1,j1),Ldt,iinfo)
!
!        Compute T3 = T(1:N1,J1:N) = -T1 Y1^H Y2 T2
!
         DO i = 1 , n1
            DO j = 1 , n2
               T(i,j+n1) = CONJG(A(j+n1,i))
            ENDDO
         ENDDO
!
         CALL ZTRMM('R','L','N','U',n1,n2,ONE,A(j1,j1),Lda,T(1,j1),Ldt)
!
         CALL ZGEMM('C','N',n1,n2,M-N,ONE,A(i1,1),Lda,A(i1,j1),Lda,ONE, &
     &              T(1,j1),Ldt)
!
         CALL ZTRMM('L','U','N','N',n1,n2,-ONE,T,Ldt,T(1,j1),Ldt)
!
         CALL ZTRMM('R','U','N','N',n1,n2,ONE,T(j1,j1),Ldt,T(1,j1),Ldt)
!
!        Y = (Y1,Y2); R = [ R1  A(1:N1,J1:N) ];  T = [T1 T3]
!                         [  0        R2     ]       [ 0 T2]
!
      ENDIF
!
!
!     End of ZGEQRT3
!
      END SUBROUTINE ZGEQRT3
