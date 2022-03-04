!*==dlacon.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DLACON estimates the 1-norm of a square matrix, using reverse communication for evaluating matrix-vector products.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLACON + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlacon.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlacon.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlacon.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLACON( N, V, X, ISGN, EST, KASE )
!
!       .. Scalar Arguments ..
!       INTEGER            KASE, N
!       DOUBLE PRECISION   EST
!       ..
!       .. Array Arguments ..
!       INTEGER            ISGN( * )
!       DOUBLE PRECISION   V( * ), X( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLACON estimates the 1-norm of a square, real matrix A.
!> Reverse communication is used for evaluating matrix-vector products.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>         The order of the matrix.  N >= 1.
!> \endverbatim
!>
!> \param[out] V
!> \verbatim
!>          V is DOUBLE PRECISION array, dimension (N)
!>         On the final return, V = A*W,  where  EST = norm(V)/norm(W)
!>         (W is not returned).
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension (N)
!>         On an intermediate return, X should be overwritten by
!>               A * X,   if KASE=1,
!>               A**T * X,  if KASE=2,
!>         and DLACON must be re-called with all the other parameters
!>         unchanged.
!> \endverbatim
!>
!> \param[out] ISGN
!> \verbatim
!>          ISGN is INTEGER array, dimension (N)
!> \endverbatim
!>
!> \param[in,out] EST
!> \verbatim
!>          EST is DOUBLE PRECISION
!>         On entry with KASE = 1 or 2 and JUMP = 3, EST should be
!>         unchanged from the previous call to DLACON.
!>         On exit, EST is an estimate (a lower bound) for norm(A).
!> \endverbatim
!>
!> \param[in,out] KASE
!> \verbatim
!>          KASE is INTEGER
!>         On the initial call to DLACON, KASE should be 0.
!>         On an intermediate return, KASE will be 1 or 2, indicating
!>         whether X should be overwritten by A * X  or A**T * X.
!>         On the final return from DLACON, KASE will again be 0.
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
!> \ingroup doubleOTHERauxiliary
!
!> \par Contributors:
!  ==================
!>
!>  Nick Higham, University of Manchester. \n
!>  Originally named SONEST, dated March 16, 1988.
!
!> \par References:
!  ================
!>
!>  N.J. Higham, "FORTRAN codes for estimating the one-norm of
!>  a real or complex matrix, with applications to condition estimation",
!>  ACM Trans. Math. Soft., vol. 14, no. 4, pp. 381-396, December 1988.
!>
!  =====================================================================
      SUBROUTINE DLACON(N,V,X,Isgn,Est,Kase)
      USE F77KINDS                        
      USE S_DASUM
      USE S_DCOPY
      USE S_IDAMAX
      IMPLICIT NONE
!*--DLACON123
!
! PARAMETER definitions rewritten by SPAG
!
      INTEGER , PARAMETER  ::  ITMAX = 5
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              TWO = 2.0D+0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: V
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: X
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Isgn
      REAL(R8KIND) , INTENT(INOUT) :: Est
      INTEGER , INTENT(INOUT) :: Kase
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) , SAVE :: altsgn , estold , temp
      INTEGER , SAVE :: i , iter , j , jlast , jump
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
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Save statement ..
!     ..
!     .. Executable Statements ..
!
      IF ( Kase==0 ) THEN
         DO i = 1 , N
            X(i) = ONE/DBLE(N)
         ENDDO
         Kase = 1
         jump = 1
         RETURN
      ENDIF
!
      IF ( jump==2 ) THEN
!
!     ................ ENTRY   (JUMP = 2)
!     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.
!
         j = IDAMAX(N,X,1)
         iter = 2
      ELSEIF ( jump==3 ) THEN
!
!     ................ ENTRY   (JUMP = 3)
!     X HAS BEEN OVERWRITTEN BY A*X.
!
         CALL DCOPY(N,X,1,V,1)
         estold = Est
         Est = DASUM(N,V,1)
         DO i = 1 , N
            IF ( NINT(SIGN(ONE,X(i)))/=Isgn(i) ) GOTO 200
         ENDDO
!     REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED.
         GOTO 300
      ELSEIF ( jump==4 ) THEN
!
!     ................ ENTRY   (JUMP = 4)
!     X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.
!
         jlast = j
         j = IDAMAX(N,X,1)
         IF ( (X(jlast)/=ABS(X(j))) .AND. (iter<ITMAX) ) THEN
            iter = iter + 1
            GOTO 100
         ENDIF
         GOTO 300
      ELSEIF ( jump==5 ) THEN
!
!     ................ ENTRY   (JUMP = 5)
!     X HAS BEEN OVERWRITTEN BY A*X.
!
         temp = TWO*(DASUM(N,X,1)/DBLE(3*N))
         IF ( temp>Est ) THEN
            CALL DCOPY(N,X,1,V,1)
            Est = temp
         ENDIF
!
         Kase = 0
         GOTO 99999
      ELSE
!
!     ................ ENTRY   (JUMP = 1)
!     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.
!
         IF ( N==1 ) THEN
            V(1) = X(1)
            Est = ABS(V(1))
!        ... QUIT
            Kase = 0
            GOTO 99999
         ENDIF
         Est = DASUM(N,X,1)
!
         DO i = 1 , N
            X(i) = SIGN(ONE,X(i))
            Isgn(i) = NINT(X(i))
         ENDDO
         Kase = 2
         jump = 2
         RETURN
      ENDIF
!
!     MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
!
 100  DO i = 1 , N
         X(i) = ZERO
      ENDDO
      X(j) = ONE
      Kase = 1
      jump = 3
      RETURN
!
!     TEST FOR CYCLING.
 200  IF ( Est>estold ) THEN
!
         DO i = 1 , N
            X(i) = SIGN(ONE,X(i))
            Isgn(i) = NINT(X(i))
         ENDDO
         Kase = 2
         jump = 4
         RETURN
      ENDIF
!
!     ITERATION COMPLETE.  FINAL STAGE.
!
 300  altsgn = ONE
      DO i = 1 , N
         X(i) = altsgn*(ONE+DBLE(i-1)/DBLE(N-1))
         altsgn = -altsgn
      ENDDO
      Kase = 1
      jump = 5
!
!     End of DLACON
!
99999 END SUBROUTINE DLACON
