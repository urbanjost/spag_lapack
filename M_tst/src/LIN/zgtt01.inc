!*==zgtt01.f90  processed by SPAG 7.51RB at 17:35 on  4 Mar 2022
!> \brief \b zgtt01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE zgtt01( N, DL, D, DU, DLF, DF, DUF, DU2, IPIV, WORK,
!                          LDWORK, RWORK, RESID )
!
!       .. Scalar Arguments ..
!       INTEGER            LDWORK, N
!       DOUBLE PRECISION   RESID
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       DOUBLE PRECISION   RWORK( * )
!       COMPLEX*16         D( * ), DF( * ), DL( * ), DLF( * ), DU( * ),
!      $                   DU2( * ), DUF( * ), WORK( LDWORK, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> zgtt01 reconstructs a tridiagonal matrix A from its LU factorization
!> and computes the residual
!>    norm(L*U - A) / ( norm(A) * EPS ),
!> where EPS is the machine epsilon.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGTER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] DL
!> \verbatim
!>          DL is COMPLEX*16 array, dimension (N-1)
!>          The (n-1) sub-diagonal elements of A.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is COMPLEX*16 array, dimension (N)
!>          The diagonal elements of A.
!> \endverbatim
!>
!> \param[in] DU
!> \verbatim
!>          DU is COMPLEX*16 array, dimension (N-1)
!>          The (n-1) super-diagonal elements of A.
!> \endverbatim
!>
!> \param[in] DLF
!> \verbatim
!>          DLF is COMPLEX*16 array, dimension (N-1)
!>          The (n-1) multipliers that define the matrix L from the
!>          LU factorization of A.
!> \endverbatim
!>
!> \param[in] DF
!> \verbatim
!>          DF is COMPLEX*16 array, dimension (N)
!>          The n diagonal elements of the upper triangular matrix U from
!>          the LU factorization of A.
!> \endverbatim
!>
!> \param[in] DUF
!> \verbatim
!>          DUF is COMPLEX*16 array, dimension (N-1)
!>          The (n-1) elements of the first super-diagonal of U.
!> \endverbatim
!>
!> \param[in] DU2
!> \verbatim
!>          DU2 is COMPLEX*16 array, dimension (N-2)
!>          The (n-2) elements of the second super-diagonal of U.
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          The pivot indices; for 1 <= i <= n, row i of the matrix was
!>          interchanged with row IPIV(i).  IPIV(i) will always be either
!>          i or i+1; IPIV(i) = i indicates a row interchange was not
!>          required.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (LDWORK,N)
!> \endverbatim
!>
!> \param[in] LDWORK
!> \verbatim
!>          LDWORK is INTEGER
!>          The leading dimension of the array WORK.  LDWORK >= max(1,N).
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is DOUBLE PRECISION
!>          The scaled residual:  norm(L*U - A) / (norm(A) * EPS)
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
!> \ingroup complex16_lin
!
!  =====================================================================
      SUBROUTINE zgtt01(N,Dl,D,Du,Dlf,Df,Duf,Du2,Ipiv,Work,Ldwork,Rwork,&
     &                  Resid)
      IMPLICIT NONE
!*--ZGTT01138
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Ldwork , N
      DOUBLE PRECISION Resid
!     ..
!     .. Array Arguments ..
      INTEGER Ipiv(*)
      DOUBLE PRECISION Rwork(*)
      COMPLEX*16 D(*) , Df(*) , Dl(*) , Dlf(*) , Du(*) , Du2(*) , Duf(*)&
     &           , Work(Ldwork,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ONE , ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
      INTEGER i , ip , j , lastj
      DOUBLE PRECISION anorm , eps
      COMPLEX*16 li
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH , ZLANGT , ZLANHS
EXTERNAL    DLAMCH,   ZLANGT,   ZLANHS
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MIN
!     ..
!     .. External Subroutines ..
EXTERNAL    ZAXPY,   ZSWAP
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF ( N<=0 ) THEN
         Resid = ZERO
         RETURN
      ENDIF
!
      eps = DLAMCH('Epsilon')
!
!     Copy the matrix U to WORK.
!
      DO j = 1 , N
         DO i = 1 , N
            Work(i,j) = ZERO
         ENDDO
      ENDDO
      DO i = 1 , N
         IF ( i==1 ) THEN
            Work(i,i) = Df(i)
            IF ( N>=2 ) Work(i,i+1) = Duf(i)
            IF ( N>=3 ) Work(i,i+2) = Du2(i)
         ELSEIF ( i==N ) THEN
            Work(i,i) = Df(i)
         ELSE
            Work(i,i) = Df(i)
            Work(i,i+1) = Duf(i)
            IF ( i<N-1 ) Work(i,i+2) = Du2(i)
         ENDIF
      ENDDO
!
!     Multiply on the left by L.
!
      lastj = N
      DO i = N - 1 , 1 , -1
         li = Dlf(i)
         CALL ZAXPY(lastj-i+1,li,Work(i,i),Ldwork,Work(i+1,i),Ldwork)
         ip = Ipiv(i)
         IF ( ip==i ) THEN
            lastj = MIN(i+2,N)
         ELSE
            CALL ZSWAP(lastj-i+1,Work(i,i),Ldwork,Work(i+1,i),Ldwork)
         ENDIF
      ENDDO
!
!     Subtract the matrix A.
!
      Work(1,1) = Work(1,1) - D(1)
      IF ( N>1 ) THEN
         Work(1,2) = Work(1,2) - Du(1)
         Work(N,N-1) = Work(N,N-1) - Dl(N-1)
         Work(N,N) = Work(N,N) - D(N)
         DO i = 2 , N - 1
            Work(i,i-1) = Work(i,i-1) - Dl(i-1)
            Work(i,i) = Work(i,i) - D(i)
            Work(i,i+1) = Work(i,i+1) - Du(i)
         ENDDO
      ENDIF
!
!     Compute the 1-norm of the tridiagonal matrix A.
!
      anorm = ZLANGT('1',N,Dl,D,Du)
!
!     Compute the 1-norm of WORK, which is only guaranteed to be
!     upper Hessenberg.
!
      Resid = ZLANHS('1',N,Work,Ldwork,Rwork)
!
!     Compute norm(L*U - A) / (norm(A) * EPS)
!
      IF ( anorm<=ZERO ) THEN
         IF ( Resid/=ZERO ) Resid = ONE/eps
      ELSE
         Resid = (Resid/anorm)/eps
      ENDIF
!
!
!     End of zgtt01
!
      END SUBROUTINE zgtt01
