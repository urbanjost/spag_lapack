!*==dbdt01.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b DBDT01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DBDT01( M, N, KD, A, LDA, Q, LDQ, D, E, PT, LDPT, WORK,
!                          RESID )
!
!       .. Scalar Arguments ..
!       INTEGER            KD, LDA, LDPT, LDQ, M, N
!       DOUBLE PRECISION   RESID
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), D( * ), E( * ), PT( LDPT, * ),
!      $                   Q( LDQ, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DBDT01 reconstructs a general matrix A from its bidiagonal form
!>    A = Q * B * P'
!> where Q (m by min(m,n)) and P' (min(m,n) by n) are orthogonal
!> matrices and B is bidiagonal.
!>
!> The test ratio to test the reduction is
!>    RESID = norm( A - Q * B * PT ) / ( n * norm(A) * EPS )
!> where PT = P' and EPS is the machine precision.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrices A and Q.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrices A and P'.
!> \endverbatim
!>
!> \param[in] KD
!> \verbatim
!>          KD is INTEGER
!>          If KD = 0, B is diagonal and the array E is not referenced.
!>          If KD = 1, the reduction was performed by xGEBRD; B is upper
!>          bidiagonal if M >= N, and lower bidiagonal if M < N.
!>          If KD = -1, the reduction was performed by xGBBRD; B is
!>          always upper bidiagonal.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          The m by n matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in] Q
!> \verbatim
!>          Q is DOUBLE PRECISION array, dimension (LDQ,N)
!>          The m by min(m,n) orthogonal matrix Q in the reduction
!>          A = Q * B * P'.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>          The leading dimension of the array Q.  LDQ >= max(1,M).
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (min(M,N))
!>          The diagonal elements of the bidiagonal matrix B.
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is DOUBLE PRECISION array, dimension (min(M,N)-1)
!>          The superdiagonal elements of the bidiagonal matrix B if
!>          m >= n, or the subdiagonal elements of B if m < n.
!> \endverbatim
!>
!> \param[in] PT
!> \verbatim
!>          PT is DOUBLE PRECISION array, dimension (LDPT,N)
!>          The min(m,n) by n orthogonal matrix P' in the reduction
!>          A = Q * B * P'.
!> \endverbatim
!>
!> \param[in] LDPT
!> \verbatim
!>          LDPT is INTEGER
!>          The leading dimension of the array PT.
!>          LDPT >= max(1,min(M,N)).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (M+N)
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is DOUBLE PRECISION
!>          The test ratio:  norm(A - Q * B * P') / ( n * norm(A) * EPS )
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
!> \ingroup double_eig
!
!  =====================================================================
      SUBROUTINE DBDT01(M,N,Kd,A,Lda,Q,Ldq,D,E,Pt,Ldpt,Work,Resid)
      IMPLICIT NONE
!*--DBDT01143
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Kd , Lda , Ldpt , Ldq , M , N
      DOUBLE PRECISION Resid
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(Lda,*) , D(*) , E(*) , Pt(Ldpt,*) , Q(Ldq,*) , &
     &                 Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
!     ..
!     .. Local Scalars ..
      INTEGER i , j
      DOUBLE PRECISION anorm , eps
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DASUM , DLAMCH , DLANGE
      EXTERNAL DASUM , DLAMCH , DLANGE
!     ..
!     .. External Subroutines ..
      EXTERNAL DCOPY , DGEMV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DBLE , MAX , MIN
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF ( M<=0 .OR. N<=0 ) THEN
         Resid = ZERO
         RETURN
      ENDIF
!
!     Compute A - Q * B * P' one column at a time.
!
      Resid = ZERO
      IF ( Kd/=0 ) THEN
!
!        B is bidiagonal.
!
         IF ( Kd/=0 .AND. M>=N ) THEN
!
!           B is upper bidiagonal and M >= N.
!
            DO j = 1 , N
               CALL DCOPY(M,A(1,j),1,Work,1)
               DO i = 1 , N - 1
                  Work(M+i) = D(i)*Pt(i,j) + E(i)*Pt(i+1,j)
               ENDDO
               Work(M+N) = D(N)*Pt(N,j)
               CALL DGEMV('No transpose',M,N,-ONE,Q,Ldq,Work(M+1),1,ONE,&
     &                    Work,1)
               Resid = MAX(Resid,DASUM(M,Work,1))
            ENDDO
         ELSEIF ( Kd<0 ) THEN
!
!           B is upper bidiagonal and M < N.
!
            DO j = 1 , N
               CALL DCOPY(M,A(1,j),1,Work,1)
               DO i = 1 , M - 1
                  Work(M+i) = D(i)*Pt(i,j) + E(i)*Pt(i+1,j)
               ENDDO
               Work(M+M) = D(M)*Pt(M,j)
               CALL DGEMV('No transpose',M,M,-ONE,Q,Ldq,Work(M+1),1,ONE,&
     &                    Work,1)
               Resid = MAX(Resid,DASUM(M,Work,1))
            ENDDO
         ELSE
!
!           B is lower bidiagonal.
!
            DO j = 1 , N
               CALL DCOPY(M,A(1,j),1,Work,1)
               Work(M+1) = D(1)*Pt(1,j)
               DO i = 2 , M
                  Work(M+i) = E(i-1)*Pt(i-1,j) + D(i)*Pt(i,j)
               ENDDO
               CALL DGEMV('No transpose',M,M,-ONE,Q,Ldq,Work(M+1),1,ONE,&
     &                    Work,1)
               Resid = MAX(Resid,DASUM(M,Work,1))
            ENDDO
         ENDIF
!
!        B is diagonal.
!
      ELSEIF ( M>=N ) THEN
         DO j = 1 , N
            CALL DCOPY(M,A(1,j),1,Work,1)
            DO i = 1 , N
               Work(M+i) = D(i)*Pt(i,j)
            ENDDO
            CALL DGEMV('No transpose',M,N,-ONE,Q,Ldq,Work(M+1),1,ONE,   &
     &                 Work,1)
            Resid = MAX(Resid,DASUM(M,Work,1))
         ENDDO
      ELSE
         DO j = 1 , N
            CALL DCOPY(M,A(1,j),1,Work,1)
            DO i = 1 , M
               Work(M+i) = D(i)*Pt(i,j)
            ENDDO
            CALL DGEMV('No transpose',M,M,-ONE,Q,Ldq,Work(M+1),1,ONE,   &
     &                 Work,1)
            Resid = MAX(Resid,DASUM(M,Work,1))
         ENDDO
      ENDIF
!
!     Compute norm(A - Q * B * P') / ( n * norm(A) * EPS )
!
      anorm = DLANGE('1',M,N,A,Lda,Work)
      eps = DLAMCH('Precision')
!
      IF ( anorm<=ZERO ) THEN
         IF ( Resid/=ZERO ) Resid = ONE/eps
      ELSEIF ( anorm>=Resid ) THEN
         Resid = (Resid/anorm)/(DBLE(N)*eps)
      ELSEIF ( anorm<ONE ) THEN
         Resid = (MIN(Resid,DBLE(N)*anorm)/anorm)/(DBLE(N)*eps)
      ELSE
         Resid = MIN(Resid/anorm,DBLE(N))/(DBLE(N)*eps)
      ENDIF
!
!
!     End of DBDT01
!
      END SUBROUTINE DBDT01
