!*==cptt01.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b cptt01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE cptt01( N, D, E, DF, EF, WORK, RESID )
!
!       .. Scalar Arguments ..
!       INTEGER            N
!       REAL               RESID
!       ..
!       .. Array Arguments ..
!       REAL               D( * ), DF( * )
!       COMPLEX            E( * ), EF( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> cptt01 reconstructs a tridiagonal matrix A from its L*D*L'
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
!>          D is REAL array, dimension (N)
!>          The n diagonal elements of the tridiagonal matrix A.
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is COMPLEX array, dimension (N-1)
!>          The (n-1) subdiagonal elements of the tridiagonal matrix A.
!> \endverbatim
!>
!> \param[in] DF
!> \verbatim
!>          DF is REAL array, dimension (N)
!>          The n diagonal elements of the factor L from the L*D*L'
!>          factorization of A.
!> \endverbatim
!>
!> \param[in] EF
!> \verbatim
!>          EF is COMPLEX array, dimension (N-1)
!>          The (n-1) subdiagonal elements of the factor L from the
!>          L*D*L' factorization of A.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (2*N)
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is REAL
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
!> \ingroup complex_lin
!
!  =====================================================================
      SUBROUTINE cptt01(N,D,E,Df,Ef,Work,Resid)
      IMPLICIT NONE
!*--CPTT0196
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER N
      REAL Resid
!     ..
!     .. Array Arguments ..
      REAL D(*) , Df(*)
      COMPLEX E(*) , Ef(*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ONE , ZERO
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
!     ..
!     .. Local Scalars ..
      INTEGER i
      REAL anorm , eps
      COMPLEX de
!     ..
!     .. External Functions ..
      REAL SLAMCH
EXTERNAL    SLAMCH
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , CONJG , MAX , REAL
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
      eps = SLAMCH('Epsilon')
!
!     Construct the difference L*D*L' - A.
!
      Work(1) = Df(1) - D(1)
      DO i = 1 , N - 1
         de = Df(i)*Ef(i)
         Work(N+i) = de - E(i)
         Work(1+i) = de*CONJG(Ef(i)) + Df(i+1) - D(i+1)
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
         Resid = ((Resid/REAL(N))/anorm)/eps
      ENDIF
!
!
!     End of cptt01
!
      END SUBROUTINE cptt01
