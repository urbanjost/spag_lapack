!*==zlacon.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b ZLACON estimates the 1-norm of a square matrix, using reverse communication for evaluating matrix-vector products.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLACON + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlacon.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlacon.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlacon.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLACON( N, V, X, EST, KASE )
!
!       .. Scalar Arguments ..
!       INTEGER            KASE, N
!       DOUBLE PRECISION   EST
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         V( N ), X( N )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLACON estimates the 1-norm of a square, complex matrix A.
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
!>          V is COMPLEX*16 array, dimension (N)
!>         On the final return, V = A*W,  where  EST = norm(V)/norm(W)
!>         (W is not returned).
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is COMPLEX*16 array, dimension (N)
!>         On an intermediate return, X should be overwritten by
!>               A * X,   if KASE=1,
!>               A**H * X,  if KASE=2,
!>         where A**H is the conjugate transpose of A, and ZLACON must be
!>         re-called with all the other parameters unchanged.
!> \endverbatim
!>
!> \param[in,out] EST
!> \verbatim
!>          EST is DOUBLE PRECISION
!>         On entry with KASE = 1 or 2 and JUMP = 3, EST should be
!>         unchanged from the previous call to ZLACON.
!>         On exit, EST is an estimate (a lower bound) for norm(A).
!> \endverbatim
!>
!> \param[in,out] KASE
!> \verbatim
!>          KASE is INTEGER
!>         On the initial call to ZLACON, KASE should be 0.
!>         On an intermediate return, KASE will be 1 or 2, indicating
!>         whether X should be overwritten by A * X  or A**H * X.
!>         On the final return from ZLACON, KASE will again be 0.
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
!> \ingroup complex16OTHERauxiliary
!
!> \par Further Details:
!  =====================
!>
!>  Originally named CONEST, dated March 16, 1988. \n
!>  Last modified:  April, 1999
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
      SUBROUTINE ZLACON(N,V,X,Est,Kase)
      IMPLICIT NONE
!*--ZLACON118
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
      COMPLEX*16 V(N) , X(N)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER ITMAX
      PARAMETER (ITMAX=5)
      DOUBLE PRECISION ONE , TWO
      PARAMETER (ONE=1.0D0,TWO=2.0D0)
      COMPLEX*16 CZERO , CONE
      PARAMETER (CZERO=(0.0D0,0.0D0),CONE=(1.0D0,0.0D0))
!     ..
!     .. Local Scalars ..
      INTEGER i , iter , j , jlast , jump
      DOUBLE PRECISION absxi , altsgn , estold , safmin , temp
!     ..
!     .. External Functions ..
      INTEGER IZMAX1
      DOUBLE PRECISION DLAMCH , DZSUM1
      EXTERNAL IZMAX1 , DLAMCH , DZSUM1
!     ..
!     .. External Subroutines ..
      EXTERNAL ZCOPY
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , DBLE , DCMPLX , DIMAG
!     ..
!     .. Save statement ..
      SAVE 
!     ..
!     .. Executable Statements ..
!
      safmin = DLAMCH('Safe minimum')
      IF ( Kase==0 ) THEN
         DO i = 1 , N
            X(i) = DCMPLX(ONE/DBLE(N))
         ENDDO
         Kase = 1
         jump = 1
         RETURN
      ENDIF
!
      IF ( jump==2 ) THEN
!
!     ................ ENTRY   (JUMP = 2)
!     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY CTRANS(A)*X.
!
         j = IZMAX1(N,X,1)
         iter = 2
      ELSEIF ( jump==3 ) THEN
!
!     ................ ENTRY   (JUMP = 3)
!     X HAS BEEN OVERWRITTEN BY A*X.
!
         CALL ZCOPY(N,X,1,V,1)
         estold = Est
         Est = DZSUM1(N,V,1)
!
!     TEST FOR CYCLING.
         IF ( Est<=estold ) GOTO 200
!
         DO i = 1 , N
            absxi = ABS(X(i))
            IF ( absxi>safmin ) THEN
               X(i) = DCMPLX(DBLE(X(i))/absxi,DIMAG(X(i))/absxi)
            ELSE
               X(i) = CONE
            ENDIF
         ENDDO
         Kase = 2
         jump = 4
         RETURN
      ELSEIF ( jump==4 ) THEN
!
!     ................ ENTRY   (JUMP = 4)
!     X HAS BEEN OVERWRITTEN BY CTRANS(A)*X.
!
         jlast = j
         j = IZMAX1(N,X,1)
         IF ( (ABS(X(jlast))/=ABS(X(j))) .AND. (iter<ITMAX) ) THEN
            iter = iter + 1
            GOTO 100
         ENDIF
         GOTO 200
      ELSEIF ( jump==5 ) THEN
!
!     ................ ENTRY   (JUMP = 5)
!     X HAS BEEN OVERWRITTEN BY A*X.
!
         temp = TWO*(DZSUM1(N,X,1)/DBLE(3*N))
         IF ( temp>Est ) THEN
            CALL ZCOPY(N,X,1,V,1)
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
         Est = DZSUM1(N,X,1)
!
         DO i = 1 , N
            absxi = ABS(X(i))
            IF ( absxi>safmin ) THEN
               X(i) = DCMPLX(DBLE(X(i))/absxi,DIMAG(X(i))/absxi)
            ELSE
               X(i) = CONE
            ENDIF
         ENDDO
         Kase = 2
         jump = 2
         RETURN
      ENDIF
!
!     MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
!
 100  DO i = 1 , N
         X(i) = CZERO
      ENDDO
      X(j) = CONE
      Kase = 1
      jump = 3
      RETURN
!
!     ITERATION COMPLETE.  FINAL STAGE.
!
 200  altsgn = ONE
      DO i = 1 , N
         X(i) = DCMPLX(altsgn*(ONE+DBLE(i-1)/DBLE(N-1)))
         altsgn = -altsgn
      ENDDO
      Kase = 1
      jump = 5
      RETURN
!
!     End of ZLACON
!
99999 END SUBROUTINE ZLACON
