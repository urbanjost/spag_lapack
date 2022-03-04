!*==slagts.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SLAGTS solves the system of equations (T-λI)x = y or (T-λI)Tx = y,where T is a general tridiagonal matrix and λ a scalar, using the LU factorization computed by slagtf.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLAGTS + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slagts.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slagts.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slagts.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLAGTS( JOB, N, A, B, C, D, IN, Y, TOL, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, JOB, N
!       REAL               TOL
!       ..
!       .. Array Arguments ..
!       INTEGER            IN( * )
!       REAL               A( * ), B( * ), C( * ), D( * ), Y( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLAGTS may be used to solve one of the systems of equations
!>
!>    (T - lambda*I)*x = y   or   (T - lambda*I)**T*x = y,
!>
!> where T is an n by n tridiagonal matrix, for x, following the
!> factorization of (T - lambda*I) as
!>
!>    (T - lambda*I) = P*L*U ,
!>
!> by routine SLAGTF. The choice of equation to be solved is
!> controlled by the argument JOB, and in each case there is an option
!> to perturb zero or very small diagonal elements of U, this option
!> being intended for use in applications such as inverse iteration.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOB
!> \verbatim
!>          JOB is INTEGER
!>          Specifies the job to be performed by SLAGTS as follows:
!>          =  1: The equations  (T - lambda*I)x = y  are to be solved,
!>                but diagonal elements of U are not to be perturbed.
!>          = -1: The equations  (T - lambda*I)x = y  are to be solved
!>                and, if overflow would otherwise occur, the diagonal
!>                elements of U are to be perturbed. See argument TOL
!>                below.
!>          =  2: The equations  (T - lambda*I)**Tx = y  are to be solved,
!>                but diagonal elements of U are not to be perturbed.
!>          = -2: The equations  (T - lambda*I)**Tx = y  are to be solved
!>                and, if overflow would otherwise occur, the diagonal
!>                elements of U are to be perturbed. See argument TOL
!>                below.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix T.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is REAL array, dimension (N)
!>          On entry, A must contain the diagonal elements of U as
!>          returned from SLAGTF.
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is REAL array, dimension (N-1)
!>          On entry, B must contain the first super-diagonal elements of
!>          U as returned from SLAGTF.
!> \endverbatim
!>
!> \param[in] C
!> \verbatim
!>          C is REAL array, dimension (N-1)
!>          On entry, C must contain the sub-diagonal elements of L as
!>          returned from SLAGTF.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is REAL array, dimension (N-2)
!>          On entry, D must contain the second super-diagonal elements
!>          of U as returned from SLAGTF.
!> \endverbatim
!>
!> \param[in] IN
!> \verbatim
!>          IN is INTEGER array, dimension (N)
!>          On entry, IN must contain details of the matrix P as returned
!>          from SLAGTF.
!> \endverbatim
!>
!> \param[in,out] Y
!> \verbatim
!>          Y is REAL array, dimension (N)
!>          On entry, the right hand side vector y.
!>          On exit, Y is overwritten by the solution vector x.
!> \endverbatim
!>
!> \param[in,out] TOL
!> \verbatim
!>          TOL is REAL
!>          On entry, with  JOB < 0, TOL should be the minimum
!>          perturbation to be made to very small diagonal elements of U.
!>          TOL should normally be chosen as about eps*norm(U), where eps
!>          is the relative machine precision, but if TOL is supplied as
!>          non-positive, then it is reset to eps*max( abs( u(i,j) ) ).
!>          If  JOB > 0  then TOL is not referenced.
!>
!>          On exit, TOL is changed as described above, only if TOL is
!>          non-positive on entry. Otherwise TOL is unchanged.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
!>          > 0: overflow would occur when computing the INFO(th)
!>               element of the solution vector x. This can only occur
!>               when JOB is supplied as positive and either means
!>               that a diagonal element of U is very small, or that
!>               the elements of the right-hand side vector y are very
!>               large.
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
!> \ingroup OTHERauxiliary
!
!  =====================================================================
      SUBROUTINE SLAGTS(Job,N,A,B,C,D,In,Y,Tol,Info)
      USE S_SLAMCH
      USE S_XERBLA
      IMPLICIT NONE
!*--SLAGTS167
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: Job
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) , DIMENSION(*) :: A
      REAL , INTENT(IN) , DIMENSION(*) :: B
      REAL , INTENT(IN) , DIMENSION(*) :: C
      REAL , INTENT(IN) , DIMENSION(*) :: D
      INTEGER , INTENT(IN) , DIMENSION(*) :: In
      REAL , INTENT(INOUT) , DIMENSION(*) :: Y
      REAL , INTENT(INOUT) :: Tol
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL :: absak , ak , bignum , eps , pert , sfmin , temp
      INTEGER :: k
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
!     .. Intrinsic Functions ..
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Executable Statements ..
!
      Info = 0
      IF ( (ABS(Job)>2) .OR. (Job==0) ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('SLAGTS',-Info)
         RETURN
      ENDIF
!
      IF ( N==0 ) RETURN
!
      eps = SLAMCH('Epsilon')
      sfmin = SLAMCH('Safe minimum')
      bignum = ONE/sfmin
!
      IF ( Job<0 ) THEN
         IF ( Tol<=ZERO ) THEN
            Tol = ABS(A(1))
            IF ( N>1 ) Tol = MAX(Tol,ABS(A(2)),ABS(B(1)))
            DO k = 3 , N
               Tol = MAX(Tol,ABS(A(k)),ABS(B(k-1)),ABS(D(k-2)))
            ENDDO
            Tol = Tol*eps
            IF ( Tol==ZERO ) Tol = eps
         ENDIF
      ENDIF
!
      IF ( ABS(Job)==1 ) THEN
         DO k = 2 , N
            IF ( In(k-1)==0 ) THEN
               Y(k) = Y(k) - C(k-1)*Y(k-1)
            ELSE
               temp = Y(k-1)
               Y(k-1) = Y(k)
               Y(k) = temp - C(k-1)*Y(k)
            ENDIF
         ENDDO
         IF ( Job==1 ) THEN
            DO k = N , 1 , -1
               IF ( k<=N-2 ) THEN
                  temp = Y(k) - B(k)*Y(k+1) - D(k)*Y(k+2)
               ELSEIF ( k==N-1 ) THEN
                  temp = Y(k) - B(k)*Y(k+1)
               ELSE
                  temp = Y(k)
               ENDIF
               ak = A(k)
               absak = ABS(ak)
               IF ( absak<ONE ) THEN
                  IF ( absak<sfmin ) THEN
                     IF ( absak==ZERO .OR. ABS(temp)*sfmin>absak ) THEN
                        Info = k
                        RETURN
                     ELSE
                        temp = temp*bignum
                        ak = ak*bignum
                     ENDIF
                  ELSEIF ( ABS(temp)>absak*bignum ) THEN
                     Info = k
                     RETURN
                  ENDIF
               ENDIF
               Y(k) = temp/ak
            ENDDO
         ELSE
            DO k = N , 1 , -1
               IF ( k<=N-2 ) THEN
                  temp = Y(k) - B(k)*Y(k+1) - D(k)*Y(k+2)
               ELSEIF ( k==N-1 ) THEN
                  temp = Y(k) - B(k)*Y(k+1)
               ELSE
                  temp = Y(k)
               ENDIF
               ak = A(k)
               pert = SIGN(Tol,ak)
               DO
                  absak = ABS(ak)
                  IF ( absak<ONE ) THEN
                     IF ( absak<sfmin ) THEN
                        IF ( absak==ZERO .OR. ABS(temp)*sfmin>absak )   &
     &                       THEN
                           ak = ak + pert
                           pert = 2*pert
                           CYCLE
                        ELSE
                           temp = temp*bignum
                           ak = ak*bignum
                        ENDIF
                     ELSEIF ( ABS(temp)>absak*bignum ) THEN
                        ak = ak + pert
                        pert = 2*pert
                        CYCLE
                     ENDIF
                  ENDIF
                  Y(k) = temp/ak
                  EXIT
               ENDDO
            ENDDO
         ENDIF
      ELSE
!
!        Come to here if  JOB = 2 or -2
!
         IF ( Job==2 ) THEN
            DO k = 1 , N
               IF ( k>=3 ) THEN
                  temp = Y(k) - B(k-1)*Y(k-1) - D(k-2)*Y(k-2)
               ELSEIF ( k==2 ) THEN
                  temp = Y(k) - B(k-1)*Y(k-1)
               ELSE
                  temp = Y(k)
               ENDIF
               ak = A(k)
               absak = ABS(ak)
               IF ( absak<ONE ) THEN
                  IF ( absak<sfmin ) THEN
                     IF ( absak==ZERO .OR. ABS(temp)*sfmin>absak ) THEN
                        Info = k
                        RETURN
                     ELSE
                        temp = temp*bignum
                        ak = ak*bignum
                     ENDIF
                  ELSEIF ( ABS(temp)>absak*bignum ) THEN
                     Info = k
                     RETURN
                  ENDIF
               ENDIF
               Y(k) = temp/ak
            ENDDO
         ELSE
            DO k = 1 , N
               IF ( k>=3 ) THEN
                  temp = Y(k) - B(k-1)*Y(k-1) - D(k-2)*Y(k-2)
               ELSEIF ( k==2 ) THEN
                  temp = Y(k) - B(k-1)*Y(k-1)
               ELSE
                  temp = Y(k)
               ENDIF
               ak = A(k)
               pert = SIGN(Tol,ak)
               DO
                  absak = ABS(ak)
                  IF ( absak<ONE ) THEN
                     IF ( absak<sfmin ) THEN
                        IF ( absak==ZERO .OR. ABS(temp)*sfmin>absak )   &
     &                       THEN
                           ak = ak + pert
                           pert = 2*pert
                           CYCLE
                        ELSE
                           temp = temp*bignum
                           ak = ak*bignum
                        ENDIF
                     ELSEIF ( ABS(temp)>absak*bignum ) THEN
                        ak = ak + pert
                        pert = 2*pert
                        CYCLE
                     ENDIF
                  ENDIF
                  Y(k) = temp/ak
                  EXIT
               ENDDO
            ENDDO
         ENDIF
!
         DO k = N , 2 , -1
            IF ( In(k-1)==0 ) THEN
               Y(k-1) = Y(k-1) - C(k-1)*Y(k)
            ELSE
               temp = Y(k-1)
               Y(k-1) = Y(k)
               Y(k) = temp - C(k-1)*Y(k)
            ENDIF
         ENDDO
      ENDIF
!
!     End of SLAGTS
!
      END SUBROUTINE SLAGTS
