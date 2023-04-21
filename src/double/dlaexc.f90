!*==dlaexc.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b DLAEXC swaps adjacent diagonal blocks of a real upper quasi-triangular matrix in Schur canonical form, by an orthogonal similarity transformation.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLAEXC + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaexc.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaexc.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaexc.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAEXC( WANTQ, N, T, LDT, Q, LDQ, J1, N1, N2, WORK,
!                          INFO )
!
!       .. Scalar Arguments ..
!       LOGICAL            WANTQ
!       INTEGER            INFO, J1, LDQ, LDT, N, N1, N2
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   Q( LDQ, * ), T( LDT, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAEXC swaps adjacent diagonal blocks T11 and T22 of order 1 or 2 in
!> an upper quasi-triangular matrix T by an orthogonal similarity
!> transformation.
!>
!> T must be in Schur canonical form, that is, block upper triangular
!> with 1-by-1 and 2-by-2 diagonal blocks; each 2-by-2 diagonal block
!> has its diagonal elements equal and its off-diagonal elements of
!> opposite sign.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] WANTQ
!> \verbatim
!>          WANTQ is LOGICAL
!>          = .TRUE. : accumulate the transformation in the matrix Q;
!>          = .FALSE.: do not accumulate the transformation.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix T. N >= 0.
!> \endverbatim
!>
!> \param[in,out] T
!> \verbatim
!>          T is DOUBLE PRECISION array, dimension (LDT,N)
!>          On entry, the upper quasi-triangular matrix T, in Schur
!>          canonical form.
!>          On exit, the updated matrix T, again in Schur canonical form.
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is INTEGER
!>          The leading dimension of the array T. LDT >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] Q
!> \verbatim
!>          Q is DOUBLE PRECISION array, dimension (LDQ,N)
!>          On entry, if WANTQ is .TRUE., the orthogonal matrix Q.
!>          On exit, if WANTQ is .TRUE., the updated matrix Q.
!>          If WANTQ is .FALSE., Q is not referenced.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>          The leading dimension of the array Q.
!>          LDQ >= 1; and if WANTQ is .TRUE., LDQ >= N.
!> \endverbatim
!>
!> \param[in] J1
!> \verbatim
!>          J1 is INTEGER
!>          The index of the first row of the first block T11.
!> \endverbatim
!>
!> \param[in] N1
!> \verbatim
!>          N1 is INTEGER
!>          The order of the first block T11. N1 = 0, 1 or 2.
!> \endverbatim
!>
!> \param[in] N2
!> \verbatim
!>          N2 is INTEGER
!>          The order of the second block T22. N2 = 0, 1 or 2.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          = 1: the transformed matrix T would be too far from Schur
!>               form; the blocks are not swapped and T and Q are
!>               unchanged.
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
!  =====================================================================
      SUBROUTINE DLAEXC(Wantq,N,T,Ldt,Q,Ldq,J1,N1,N2,Work,Info)
      IMPLICIT NONE
!*--DLAEXC141
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      LOGICAL Wantq
      INTEGER Info , J1 , Ldq , Ldt , N , N1 , N2
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION Q(Ldq,*) , T(Ldt,*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
      DOUBLE PRECISION TEN
      PARAMETER (TEN=1.0D+1)
      INTEGER LDD , LDX
      PARAMETER (LDD=4,LDX=2)
!     ..
!     .. Local Scalars ..
      INTEGER ierr , j2 , j3 , j4 , k , nd
      DOUBLE PRECISION cs , dnorm , eps , scale , smlnum , sn , t11 ,   &
     &                 t22 , t33 , tau , tau1 , tau2 , temp , thresh ,  &
     &                 wi1 , wi2 , wr1 , wr2 , xnorm
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION d(LDD,4) , u(3) , u1(3) , u2(3) , x(LDX,2)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH , DLANGE
      EXTERNAL DLAMCH , DLANGE
!     ..
!     .. External Subroutines ..
      EXTERNAL DLACPY , DLANV2 , DLARFG , DLARFX , DLARTG , DLASY2 ,    &
     &         DROT
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX
!     ..
!     .. Executable Statements ..
!
      Info = 0
!
!     Quick return if possible
!
      IF ( N==0 .OR. N1==0 .OR. N2==0 ) RETURN
      IF ( J1+N1>N ) RETURN
!
      j2 = J1 + 1
      j3 = J1 + 2
      j4 = J1 + 3
!
      IF ( N1==1 .AND. N2==1 ) THEN
!
!        Swap two 1-by-1 blocks.
!
         t11 = T(J1,J1)
         t22 = T(j2,j2)
!
!        Determine the transformation to perform the interchange.
!
         CALL DLARTG(T(J1,j2),t22-t11,cs,sn,temp)
!
!        Apply transformation to the matrix T.
!
         IF ( j3<=N ) CALL DROT(N-J1-1,T(J1,j3),Ldt,T(j2,j3),Ldt,cs,sn)
         CALL DROT(J1-1,T(1,J1),1,T(1,j2),1,cs,sn)
!
         T(J1,J1) = t22
         T(j2,j2) = t11
!
!
!           Accumulate transformation in the matrix Q.
!
         IF ( Wantq ) CALL DROT(N,Q(1,J1),1,Q(1,j2),1,cs,sn)
!
      ELSE
!
!        Swapping involves at least one 2-by-2 block.
!
!        Copy the diagonal block of order N1+N2 to the local array D
!        and compute its norm.
!
         nd = N1 + N2
         CALL DLACPY('Full',nd,nd,T(J1,J1),Ldt,d,LDD)
         dnorm = DLANGE('Max',nd,nd,d,LDD,Work)
!
!        Compute machine-dependent threshold for test for accepting
!        swap.
!
         eps = DLAMCH('P')
         smlnum = DLAMCH('S')/eps
         thresh = MAX(TEN*eps*dnorm,smlnum)
!
!        Solve T11*X - X*T22 = scale*T12 for X.
!
         CALL DLASY2(.FALSE.,.FALSE.,-1,N1,N2,d,LDD,d(N1+1,N1+1),LDD,   &
     &               d(1,N1+1),LDD,scale,x,LDX,xnorm,ierr)
!
!        Swap the adjacent diagonal blocks.
!
         k = N1 + N1 + N2 - 3
         IF ( k==2 ) THEN
!
!
!        N1 = 2, N2 = 1: generate elementary reflector H so that:
!
!        H (  -X11 ) = ( * )
!          (  -X21 ) = ( 0 )
!          ( scale ) = ( 0 )
!
            u(1) = -x(1,1)
            u(2) = -x(2,1)
            u(3) = scale
            CALL DLARFG(3,u(1),u(2),1,tau)
            u(1) = ONE
            t33 = T(j3,j3)
!
!        Perform swap provisionally on diagonal block in D.
!
            CALL DLARFX('L',3,3,u,tau,d,LDD,Work)
            CALL DLARFX('R',3,3,u,tau,d,LDD,Work)
!
!        Test whether to reject swap.
!
            IF ( MAX(ABS(d(2,1)),ABS(d(3,1)),ABS(d(1,1)-t33))>thresh )  &
     &           THEN
!
!     Exit with INFO = 1 if swap was rejected.
!
               Info = 1
               GOTO 99999
            ELSE
!
!        Accept swap: apply transformation to the entire matrix T.
!
               CALL DLARFX('R',j3,3,u,tau,T(1,J1),Ldt,Work)
               CALL DLARFX('L',3,N-J1,u,tau,T(J1,j2),Ldt,Work)
!
               T(J1,J1) = t33
               T(j2,J1) = ZERO
               T(j3,J1) = ZERO
!
!
!           Accumulate transformation in the matrix Q.
!
               IF ( Wantq ) CALL DLARFX('R',N,3,u,tau,Q(1,J1),Ldq,Work)
            ENDIF
         ELSEIF ( k==3 ) THEN
!
!
!        N1 = 2, N2 = 2: generate elementary reflectors H(1) and H(2) so
!        that:
!
!        H(2) H(1) (  -X11  -X12 ) = (  *  * )
!                  (  -X21  -X22 )   (  0  * )
!                  ( scale    0  )   (  0  0 )
!                  (    0  scale )   (  0  0 )
!
            u1(1) = -x(1,1)
            u1(2) = -x(2,1)
            u1(3) = scale
            CALL DLARFG(3,u1(1),u1(2),1,tau1)
            u1(1) = ONE
!
            temp = -tau1*(x(1,2)+u1(2)*x(2,2))
            u2(1) = -temp*u1(2) - x(2,2)
            u2(2) = -temp*u1(3)
            u2(3) = scale
            CALL DLARFG(3,u2(1),u2(2),1,tau2)
            u2(1) = ONE
!
!        Perform swap provisionally on diagonal block in D.
!
            CALL DLARFX('L',3,4,u1,tau1,d,LDD,Work)
            CALL DLARFX('R',4,3,u1,tau1,d,LDD,Work)
            CALL DLARFX('L',3,4,u2,tau2,d(2,1),LDD,Work)
            CALL DLARFX('R',4,3,u2,tau2,d(1,2),LDD,Work)
!
!        Test whether to reject swap.
!
            IF ( MAX(ABS(d(3,1)),ABS(d(3,2)),ABS(d(4,1)),ABS(d(4,2)))   &
     &           >thresh ) THEN
               Info = 1
               GOTO 99999
            ELSE
!
!        Accept swap: apply transformation to the entire matrix T.
!
               CALL DLARFX('L',3,N-J1+1,u1,tau1,T(J1,J1),Ldt,Work)
               CALL DLARFX('R',j4,3,u1,tau1,T(1,J1),Ldt,Work)
               CALL DLARFX('L',3,N-J1+1,u2,tau2,T(j2,J1),Ldt,Work)
               CALL DLARFX('R',j4,3,u2,tau2,T(1,j2),Ldt,Work)
!
               T(j3,J1) = ZERO
               T(j3,j2) = ZERO
               T(j4,J1) = ZERO
               T(j4,j2) = ZERO
!
               IF ( Wantq ) THEN
!
!           Accumulate transformation in the matrix Q.
!
                  CALL DLARFX('R',N,3,u1,tau1,Q(1,J1),Ldq,Work)
                  CALL DLARFX('R',N,3,u2,tau2,Q(1,j2),Ldq,Work)
               ENDIF
            ENDIF
         ELSE
!
!
!        N1 = 1, N2 = 2: generate elementary reflector H so that:
!
!        ( scale, X11, X12 ) H = ( 0, 0, * )
!
            u(1) = scale
            u(2) = x(1,1)
            u(3) = x(1,2)
            CALL DLARFG(3,u(3),u,1,tau)
            u(3) = ONE
            t11 = T(J1,J1)
!
!        Perform swap provisionally on diagonal block in D.
!
            CALL DLARFX('L',3,3,u,tau,d,LDD,Work)
            CALL DLARFX('R',3,3,u,tau,d,LDD,Work)
!
!        Test whether to reject swap.
!
            IF ( MAX(ABS(d(3,1)),ABS(d(3,2)),ABS(d(3,3)-t11))>thresh )  &
     &           THEN
               Info = 1
               GOTO 99999
            ELSE
!
!        Accept swap: apply transformation to the entire matrix T.
!
               CALL DLARFX('L',3,N-J1+1,u,tau,T(J1,J1),Ldt,Work)
               CALL DLARFX('R',j2,3,u,tau,T(1,J1),Ldt,Work)
!
               T(j3,J1) = ZERO
               T(j3,j2) = ZERO
               T(j3,j3) = t11
!
!
!           Accumulate transformation in the matrix Q.
!
               IF ( Wantq ) CALL DLARFX('R',N,3,u,tau,Q(1,J1),Ldq,Work)
            ENDIF
         ENDIF
!
!
         IF ( N2==2 ) THEN
!
!           Standardize new 2-by-2 block T11
!
            CALL DLANV2(T(J1,J1),T(J1,j2),T(j2,J1),T(j2,j2),wr1,wi1,wr2,&
     &                  wi2,cs,sn)
            CALL DROT(N-J1-1,T(J1,J1+2),Ldt,T(j2,J1+2),Ldt,cs,sn)
            CALL DROT(J1-1,T(1,J1),1,T(1,j2),1,cs,sn)
            IF ( Wantq ) CALL DROT(N,Q(1,J1),1,Q(1,j2),1,cs,sn)
         ENDIF
!
         IF ( N1==2 ) THEN
!
!           Standardize new 2-by-2 block T22
!
            j3 = J1 + N2
            j4 = j3 + 1
            CALL DLANV2(T(j3,j3),T(j3,j4),T(j4,j3),T(j4,j4),wr1,wi1,wr2,&
     &                  wi2,cs,sn)
            IF ( j3+2<=N ) CALL DROT(N-j3-1,T(j3,j3+2),Ldt,T(j4,j3+2),  &
     &                               Ldt,cs,sn)
            CALL DROT(j3-1,T(1,j3),1,T(1,j4),1,cs,sn)
            IF ( Wantq ) CALL DROT(N,Q(1,j3),1,Q(1,j4),1,cs,sn)
         ENDIF
!
      ENDIF
      RETURN
!
!     End of DLAEXC
!
99999 END SUBROUTINE DLAEXC
