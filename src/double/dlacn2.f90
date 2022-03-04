!*==dlacn2.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b DLACN2 estimates the 1-norm of a square matrix, using reverse communication for evaluating matrix-vector products.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLACN2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlacn2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlacn2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlacn2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLACN2( N, V, X, ISGN, EST, KASE, ISAVE )
!
!       .. Scalar Arguments ..
!       INTEGER            KASE, N
!       DOUBLE PRECISION   EST
!       ..
!       .. Array Arguments ..
!       INTEGER            ISGN( * ), ISAVE( 3 )
!       DOUBLE PRECISION   V( * ), X( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLACN2 estimates the 1-norm of a square, real matrix A.
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
!>         and DLACN2 must be re-called with all the other parameters
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
!>         On entry with KASE = 1 or 2 and ISAVE(1) = 3, EST should be
!>         unchanged from the previous call to DLACN2.
!>         On exit, EST is an estimate (a lower bound) for norm(A).
!> \endverbatim
!>
!> \param[in,out] KASE
!> \verbatim
!>          KASE is INTEGER
!>         On the initial call to DLACN2, KASE should be 0.
!>         On an intermediate return, KASE will be 1 or 2, indicating
!>         whether X should be overwritten by A * X  or A**T * X.
!>         On the final return from DLACN2, KASE will again be 0.
!> \endverbatim
!>
!> \param[in,out] ISAVE
!> \verbatim
!>          ISAVE is INTEGER array, dimension (3)
!>         ISAVE is used to save variables between calls to DLACN2
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
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Originally named SONEST, dated March 16, 1988.
!>
!>  This is a thread safe version of DLACON, which uses the array ISAVE
!>  in place of a SAVE statement, as follows:
!>
!>     DLACON     DLACN2
!>      JUMP     ISAVE(1)
!>      J        ISAVE(2)
!>      ITER     ISAVE(3)
!> \endverbatim
!
!> \par Contributors:
!  ==================
!>
!>     Nick Higham, University of Manchester
!
!> \par References:
!  ================
!>
!>  N.J. Higham, "FORTRAN codes for estimating the one-norm of
!>  a real or complex matrix, with applications to condition estimation",
!>  ACM Trans. Math. Soft., vol. 14, no. 4, pp. 381-396, December 1988.
!>
!  =====================================================================
      SUBROUTINE DLACN2(N,V,X,Isgn,Est,Kase,Isave)
      IMPLICIT NONE
!*--DLACN2140
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Kase , N
      DOUBLE PRECISION Est
!     ..
!     .. Array Arguments ..
      INTEGER Isgn(*) , Isave(3)
      DOUBLE PRECISION V(*) , X(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER ITMAX
      PARAMETER (ITMAX=5)
      DOUBLE PRECISION ZERO , ONE , TWO
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0,TWO=2.0D+0)
!     ..
!     .. Local Scalars ..
      INTEGER i , jlast
      DOUBLE PRECISION altsgn , estold , temp , xs
!     ..
!     .. External Functions ..
      INTEGER IDAMAX
      DOUBLE PRECISION DASUM
      EXTERNAL IDAMAX , DASUM
!     ..
!     .. External Subroutines ..
      EXTERNAL DCOPY
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , DBLE , NINT
!     ..
!     .. Executable Statements ..
!
      IF ( Kase==0 ) THEN
         DO i = 1 , N
            X(i) = ONE/DBLE(N)
         ENDDO
         Kase = 1
         Isave(1) = 1
         RETURN
      ENDIF
!
      IF ( Isave(1)==2 ) THEN
!
!     ................ ENTRY   (ISAVE( 1 ) = 2)
!     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.
!
         Isave(2) = IDAMAX(N,X,1)
         Isave(3) = 2
      ELSEIF ( Isave(1)==3 ) THEN
!
!     ................ ENTRY   (ISAVE( 1 ) = 3)
!     X HAS BEEN OVERWRITTEN BY A*X.
!
         CALL DCOPY(N,X,1,V,1)
         estold = Est
         Est = DASUM(N,V,1)
         DO i = 1 , N
            IF ( X(i)>=ZERO ) THEN
               xs = ONE
            ELSE
               xs = -ONE
            ENDIF
            IF ( NINT(xs)/=Isgn(i) ) GOTO 200
         ENDDO
!     REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED.
         GOTO 300
      ELSEIF ( Isave(1)==4 ) THEN
!
!     ................ ENTRY   (ISAVE( 1 ) = 4)
!     X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.
!
         jlast = Isave(2)
         Isave(2) = IDAMAX(N,X,1)
         IF ( (X(jlast)/=ABS(X(Isave(2)))) .AND. (Isave(3)<ITMAX) ) THEN
            Isave(3) = Isave(3) + 1
            GOTO 100
         ENDIF
         GOTO 300
      ELSEIF ( Isave(1)==5 ) THEN
!
!     ................ ENTRY   (ISAVE( 1 ) = 5)
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
!     ................ ENTRY   (ISAVE( 1 ) = 1)
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
            IF ( X(i)>=ZERO ) THEN
               X(i) = ONE
            ELSE
               X(i) = -ONE
            ENDIF
            Isgn(i) = NINT(X(i))
         ENDDO
         Kase = 2
         Isave(1) = 2
         RETURN
      ENDIF
!
!     MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
!
 100  DO i = 1 , N
         X(i) = ZERO
      ENDDO
      X(Isave(2)) = ONE
      Kase = 1
      Isave(1) = 3
      RETURN
!
!     TEST FOR CYCLING.
 200  IF ( Est>estold ) THEN
!
         DO i = 1 , N
            IF ( X(i)>=ZERO ) THEN
               X(i) = ONE
            ELSE
               X(i) = -ONE
            ENDIF
            Isgn(i) = NINT(X(i))
         ENDDO
         Kase = 2
         Isave(1) = 4
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
      Isave(1) = 5
      RETURN
!
!     End of DLACN2
!
99999 END SUBROUTINE DLACN2
