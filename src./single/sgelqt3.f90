!*==sgelqt3.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SGELQT3
!
!  Definition:
!  ===========
!
!       RECURSIVE SUBROUTINE SGELQT3( M, N, A, LDA, T, LDT, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER   INFO, LDA, M, N, LDT
!       ..
!       .. Array Arguments ..
!       REAL   A( LDA, * ), T( LDT, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGELQT3 recursively computes a LQ factorization of a real M-by-N
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
!>          The number of rows of the matrix A.  M =< N.
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
!>          A is REAL array, dimension (LDA,N)
!>          On entry, the real M-by-N matrix A.  On exit, the elements on and
!>          below the diagonal contain the N-by-N lower triangular matrix L; the
!>          elements above the diagonal are the rows of V.  See below for
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
!>          T is REAL array, dimension (LDT,N)
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
!> \date November 2017
!
!> \ingroup doubleGEcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The matrix V stores the elementary reflectors H(i) in the i-th row
!>  above the diagonal. For example, if M=5 and N=3, the matrix V is
!>
!>               V = (  1  v1 v1 v1 v1 )
!>                   (     1  v2 v2 v2 )
!>                   (     1  v3 v3 v3 )
!>
!>
!>  where the vi's represent the vectors which define H(i), which are returned
!>  in the matrix A.  The 1's along the diagonal of V are not stored in A.  The
!>  block reflector H is then given by
!>
!>               H = I - V * T * V**T
!>
!>  where V**T is the transpose of V.
!>
!>  For details of the algorithm, see Elmroth and Gustavson (cited above).
!> \endverbatim
!>
!  =====================================================================
      RECURSIVE SUBROUTINE SGELQT3(M,N,A,Lda,T,Ldt,Info)
      USE S_SGEMM
      USE S_SLARFG
      USE S_STRMM
      USE S_XERBLA
      IMPLICIT NONE
!*--SGELQT3124
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ONE = 1.0E+00
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: M
      INTEGER :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(INOUT) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , i1 , iinfo , j , j1 , m1 , m2
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
      IF ( M<0 ) THEN
         Info = -1
      ELSEIF ( N<M ) THEN
         Info = -2
      ELSEIF ( Lda<MAX(1,M) ) THEN
         Info = -4
      ELSEIF ( Ldt<MAX(1,M) ) THEN
         Info = -6
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('SGELQT3',-Info)
         RETURN
      ENDIF
!
      IF ( M==1 ) THEN
!
!        Compute Householder transform when N=1
!
         CALL SLARFG(N,A,A(1,MIN(2,N)),Lda,T)
!
      ELSE
!
!        Otherwise, split A into blocks...
!
         m1 = M/2
         m2 = M - m1
         i1 = MIN(m1+1,M)
         j1 = MIN(M+1,N)
!
!        Compute A(1:M1,1:N) <- (Y1,R1,T1), where Q1 = I - Y1 T1 Y1^H
!
         CALL SGELQT3(m1,N,A,Lda,T,Ldt,iinfo)
!
!        Compute A(J1:M,1:N) = Q1^H A(J1:M,1:N) [workspace: T(1:N1,J1:N)]
!
         DO i = 1 , m2
            DO j = 1 , m1
               T(i+m1,j) = A(i+m1,j)
            ENDDO
         ENDDO
         CALL STRMM('R','U','T','U',m2,m1,ONE,A,Lda,T(i1,1),Ldt)
!
         CALL SGEMM('N','T',m2,m1,N-m1,ONE,A(i1,i1),Lda,A(1,i1),Lda,ONE,&
     &              T(i1,1),Ldt)
!
         CALL STRMM('R','U','N','N',m2,m1,ONE,T,Ldt,T(i1,1),Ldt)
!
         CALL SGEMM('N','N',m2,N-m1,m1,-ONE,T(i1,1),Ldt,A(1,i1),Lda,ONE,&
     &              A(i1,i1),Lda)
!
         CALL STRMM('R','U','N','U',m2,m1,ONE,A,Lda,T(i1,1),Ldt)
!
         DO i = 1 , m2
            DO j = 1 , m1
               A(i+m1,j) = A(i+m1,j) - T(i+m1,j)
               T(i+m1,j) = 0
            ENDDO
         ENDDO
!
!        Compute A(J1:M,J1:N) <- (Y2,R2,T2) where Q2 = I - Y2 T2 Y2^H
!
         CALL SGELQT3(m2,N-m1,A(i1,i1),Lda,T(i1,i1),Ldt,iinfo)
!
!        Compute T3 = T(J1:N1,1:N) = -T1 Y1^H Y2 T2
!
         DO i = 1 , m2
            DO j = 1 , m1
               T(j,i+m1) = (A(j,i+m1))
            ENDDO
         ENDDO
!
         CALL STRMM('R','U','T','U',m1,m2,ONE,A(i1,i1),Lda,T(1,i1),Ldt)
!
         CALL SGEMM('N','T',m1,m2,N-M,ONE,A(1,j1),Lda,A(i1,j1),Lda,ONE, &
     &              T(1,i1),Ldt)
!
         CALL STRMM('L','U','N','N',m1,m2,-ONE,T,Ldt,T(1,i1),Ldt)
!
         CALL STRMM('R','U','N','N',m1,m2,ONE,T(i1,i1),Ldt,T(1,i1),Ldt)
!
!
!
!        Y = (Y1,Y2); L = [ L1            0  ];  T = [T1 T3]
!                         [ A(1:N1,J1:N)  L2 ]       [ 0 T2]
!
      ENDIF
!
!
!     End of SGELQT3
!
      END SUBROUTINE SGELQT3
