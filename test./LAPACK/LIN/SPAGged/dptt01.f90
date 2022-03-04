!*==dptt01.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b DPTT01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DPTT01( N, D, E, DF, EF, WORK, RESID )
!
!       .. Scalar Arguments ..
!       INTEGER            N
!       DOUBLE PRECISION   RESID
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   D( * ), DF( * ), E( * ), EF( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DPTT01 reconstructs a tridiagonal matrix A from its L*D*L'
!> factorization and computes the residual
!>    norm(L*D*L' - A) / ( n * norm(A) * EPS ),
!> where EPS is the machine epsilon.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGTER
!>          The order of the matrix A.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>          The n diagonal elements of the tridiagonal matrix A.
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is DOUBLE PRECISION array, dimension (N-1)
!>          The (n-1) subdiagonal elements of the tridiagonal matrix A.
!> \endverbatim
!>
!> \param[in] DF
!> \verbatim
!>          DF is DOUBLE PRECISION array, dimension (N)
!>          The n diagonal elements of the factor L from the L*D*L'
!>          factorization of A.
!> \endverbatim
!>
!> \param[in] EF
!> \verbatim
!>          EF is DOUBLE PRECISION array, dimension (N-1)
!>          The (n-1) subdiagonal elements of the factor L from the
!>          L*D*L' factorization of A.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (2*N)
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is DOUBLE PRECISION
!>          norm(L*D*L' - A) / (n * norm(A) * EPS)
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
!> \ingroup double_lin
!
!  =====================================================================
      SUBROUTINE DPTT01(N,D,E,Df,Ef,Work,Resid)
      IMPLICIT NONE
!*--DPTT0195
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER N
      DOUBLE PRECISION Resid
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION D(*) , Df(*) , E(*) , Ef(*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ONE , ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
      INTEGER i
      DOUBLE PRECISION anorm , de , eps
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH
      EXTERNAL DLAMCH
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , DBLE , MAX
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
!     Construct the difference L*D*L' - A.
!
      Work(1) = Df(1) - D(1)
      DO i = 1 , N - 1
         de = Df(i)*Ef(i)
         Work(N+i) = de - E(i)
         Work(1+i) = de*Ef(i) + Df(i+1) - D(i+1)
      ENDDO
!
!     Compute the 1-norms of the tridiagonal matrices A and WORK.
!
      IF ( N==1 ) THEN
         anorm = D(1)
         Resid = ABS(Work(1))
      ELSE
         anorm = MAX(D(1)+ABS(E(1)),D(N)+ABS(E(N-1)))
         Resid = MAX(ABS(Work(1))+ABS(Work(N+1)),ABS(Work(N))           &
     &           +ABS(Work(2*N-1)))
         DO i = 2 , N - 1
            anorm = MAX(anorm,D(i)+ABS(E(i))+ABS(E(i-1)))
            Resid = MAX(Resid,ABS(Work(i))+ABS(Work(N+i-1))             &
     &              +ABS(Work(N+i)))
         ENDDO
      ENDIF
!
!     Compute norm(L*D*L' - A) / (n * norm(A) * EPS)
!
      IF ( anorm<=ZERO ) THEN
         IF ( Resid/=ZERO ) Resid = ONE/eps
      ELSE
         Resid = ((Resid/DBLE(N))/anorm)/eps
      ENDIF
!
!
!     End of DPTT01
!
      END SUBROUTINE DPTT01
